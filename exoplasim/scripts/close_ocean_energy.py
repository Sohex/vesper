#!/usr/bin/env python3
"""Close the surface energy budget against the OCEAN's and the ICE module's own books.

`close_state_energy.py` compares the atmosphere's surface flux with a heat content
this project REBUILDS from `ts`, `mld` and the rest. That comparison leaves a
residual on the ice-free ocean, and it cannot say which side carries it, because
both sides are read off the same output stream.

PlaSim writes two further streams that this project has never postprocessed, and
they are the other side. `oceanmod.f90:oceanout` writes `ocean_output` and
`icemod.f90:iceout` writes `ice_output`: raw SERVICE-format files, one 8-integer
header record plus one `NLON*NLAT` float32 record per field, with their own
accumulators and their own output interval (`nout` timesteps, `naccuout` samples).
They carry what the ocean RECEIVED and what its state DID, independently of the
atmospheric diagnostics. The codes used here:

    ice_output    701 xheata   flux from the atmosphere, as the ice module got it
                  702 xofluxa  flux from the ocean into the ice
                  703 xtsfluxa heat spent changing the surface temperature
                  704 xsmelta  heat spent melting snow
                  705 ximelta  heat spent melting ice
                  706 xcfluxa  flux passed on to the ocean
                  713 xcfluxra residual from the ice-thickness limiter
                  710/711 ice compactness and thickness, 741 snow, 772 land mask
    ocean_output  901 yheata   flux as the slab got it     (== 706 by construction)
                  902 yifluxa  flux diverted into sea ice
                  903/904/905/906 flux correction, vertical and horizontal
                              diffusion, deep ocean -- all identically zero in
                              this configuration, and checked to be
                  939 ysst     the slab temperature itself, INSTANTANEOUS
                  972 yls      land mask

THE INSTRUMENT ONLY EXISTS FOR THE LAST SEGMENT OF A RUN. `oceanmod.f90:331`
and its counterpart in `icemod.f90` open these files with a bare
`open(unit,file=...,form='unformatted')`, which truncates. Every model call
therefore discards the previous call's stream and writes only its own. This
script REFUSES a window that is not the run's final segment rather than
returning numbers for orbits whose records were overwritten, in the same spirit
as `close_term_energy.py` refusing a low-I/O window: a measurement that cannot
be made should fail loudly, not quietly return the wrong orbit.

WHAT IS TESTED. Six identities, each with a right answer of zero, on ocean cells
that carry neither ice nor snow at any point in the window -- a mask taken from
the ice stream at its own resolution, which is far stricter than the same mask
taken from twelve binned output records.

    D1  xheat  = rss + rls + hfss + hfls          the handoff
    D2  yheat  = xcflux                           ice module to ocean module
    D3  xheat - xcflux = xsmelt                   the only term withheld is the
                                                  fusion of snow falling into
                                                  open water
    D4  CRHOS*CPS*mld*d(SST)/dt = yheat           the slab integration itself
    D5  hfns - (rss+rls+hfss+hfls) = -ALF*rho*snm pyburn's own definition of hfns
    D6  hfns = CRHOS*CPS*mld*d(SST)/dt            the surface residual, end to end

D2, D3, D4 and D5 need only one stream each and are independent of how the two
are aligned in time. D1 and D6 are cross-stream and need the alignment below.

THE ALIGNMENT is bookkeeping and not the measurement. The ocean and ice streams
write every `nout` steps and the regular stream every `nafter`, and `nout`
divides `nafter`, so a whole number of stream records tiles each output bin.
This script finds that tiling by scanning the offset, and the evidence that the
two streams describe the same orbit at all is that ONE offset makes eleven
independent bins of two separately accumulated arrays agree to a fraction of a
percent of the signal while the others do not. The scan is reported in full.

THE FIRST BIN IS NOT LIKE THE OTHERS. A model call leaves a partial output
interval behind, the next call restores it, and the first record of the next
file therefore covers more timesteps than the rest -- `first-output-bin.md` for
the mechanism and its other consequences. pyburn averages the twelve bins with
EQUAL weight regardless, so an annual mean taken from a 12-bin file is a
weighted mean with the wrong weights. This script measures the size of that:
the tiling gives the span of bins 1 to 11 exactly, so the first bin's excess
span follows, and the error it puts into the annual mean is

    (1/12 - n0/(n0 + 11*n)) * (bin 0 - mean of bins 1 to 11)

which is reported as `first_bin_weight`. It is not a correction applied to
anything here.

Usage:

    python exoplasim/scripts/close_ocean_energy.py <run_dir>

with the window taken from the run's own manifest. Evidence and interpretation:
exoplasim/notes/water-and-energy-closure.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import struct
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS

# PlaSim's own constants, from plasim/src. Properties of the compiled model, not
# of this planet, so they come from the model rather than from planet.yaml.
CRHOS = 1030.0        # oceanmod.f90: density of sea water, kg/m3
CPS = 4180.0          # oceanmod.f90: specific heat of sea water, J/kg/K
ALV = 2.5008e6        # plasimmod.f90: latent heat of vaporisation, J/kg
ALS = 2.8345e6        # plasimmod.f90: latent heat of sublimation, J/kg
ALF = ALS - ALV       # fusion
RHO_WATER = 1000.0

OCEAN_CODES = {901: "yheat", 902: "yiflux", 903: "yfsst", 904: "ydsst",
               905: "yqhd", 906: "yfldo", 910: "yicec", 939: "ysst",
               972: "yls", 990: "yclsst"}
ICE_CODES = {701: "xheat", 702: "xoflux", 703: "xtsflux", 704: "xsmelt",
             705: "ximelt", 706: "xcflux", 708: "xqmelt", 709: "xflxice",
             710: "xicec", 711: "xiced", 712: "xscflx", 713: "xcfluxr",
             714: "xcfluxn", 739: "xts", 741: "xsnow", 769: "xsst", 772: "xls"}
# Zero in this configuration, and the run is not the one described if they are not.
MUST_BE_ZERO = {"yfsst": "flux correction, nfluko = 0",
                "ydsst": "ocean vertical diffusion, NLEV_OCE = 1",
                "yqhd": "ocean horizontal diffusion, nhdiff = 0",
                "yfldo": "deep-ocean flux, NLSG = 0"}


def read_service(path: Path, keep: dict[int, str]) -> dict[str, np.ndarray]:
    """Every record of a PlaSim service-format stream, as (time, lat, lon).

    Each field is two Fortran unformatted records: an 8-integer header whose
    first entry is the field code, then NLON*NLAT float32 values. Read one
    record at a time and keep only the requested codes, so the working set is
    the fields asked for rather than the whole file.
    """
    size = os.path.getsize(path)
    out: dict[int, list[np.ndarray]] = {}
    shape = None
    with open(path, "rb") as handle:
        while handle.tell() < size:
            length = struct.unpack("i", handle.read(4))[0]
            head = np.frombuffer(handle.read(length), dtype=np.int32)
            handle.read(4)
            length = struct.unpack("i", handle.read(4))[0]
            raw = handle.read(length)
            handle.read(4)
            nlon, nlat = int(head[4]), int(head[5])
            if shape is None:
                shape = (nlat, nlon)
            elif shape != (nlat, nlon):
                raise SystemExit(f"{path} changes grid mid-file")
            code = int(head[0])
            if code in keep:
                out.setdefault(code, []).append(
                    np.frombuffer(raw, dtype=np.float32).reshape(shape))
    missing = [name for code, name in keep.items() if code not in out]
    if missing:
        raise SystemExit(f"{path} does not carry {missing}")
    counts = {len(v) for v in out.values()}
    if len(counts) != 1:
        raise SystemExit(f"{path} has unequal record counts per field: {counts}")
    return {keep[code]: np.stack(v).astype(float) for code, v in out.items()}


def final_orbit(manifest: dict) -> int:
    """The last orbit of the run, which is the only one the streams still hold.

    The truncating `open` discards the stream at every model call, and the
    wrapper calls the model once per orbit, so what survives on disk is the
    final orbit and nothing else -- not the final segment, which may be several.
    """
    segments = manifest.get("segments", [])
    if not segments:
        raise SystemExit("the run manifest records no segments, so which orbit "
                         "the ocean stream covers cannot be established")
    return int(segments[-1]["end_year_index"])


def align(stream: np.ndarray, binned: np.ndarray, weights: np.ndarray,
          mask: np.ndarray) -> dict:
    """Find the offset at which stream records tile the binned output.

    `nbin` bins of `per` stream records each; the offset is how many records
    precede the first bin. Returns the scan so the discrimination is visible.
    """
    nbin = binned.shape[0]
    per = stream.shape[0] // nbin
    scan = []
    for offset in range(stream.shape[0] - per * nbin + 1):
        tiled = stream[offset:offset + per * nbin]
        tiled = tiled.reshape(nbin, per, *stream.shape[1:]).mean(1)
        diff = np.array([_mean(tiled[k] - binned[k], weights, mask)
                         for k in range(nbin)])
        scan.append({"offset": offset,
                     "per_bin_w_m2": [float(x) for x in diff],
                     # bin 0 straddles the restart and is excluded from the score
                     "score_w_m2": float(np.abs(diff[1:]).mean())})
    best = min(scan, key=lambda s: s["score_w_m2"])
    return {"records_per_bin": per, "chosen_offset": best["offset"],
            "score_w_m2": best["score_w_m2"], "scan": scan}


def _mean(field: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> float:
    """Gaussian-weighted mean over the mask, per unit area of the mask."""
    area = float((mask * weights).sum())
    if area <= 0.0:
        return float("nan")
    return float((field * weights * mask).sum()) / area


def close_ocean(run_dir: Path) -> dict:
    manifest_path = run_dir / "run_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest_path.is_file() else {})
    first = last = final_orbit(manifest)

    ocean = read_service(run_dir / "ocean_output", OCEAN_CODES)
    ice = read_service(run_dir / "ice_output", ICE_CODES)
    nrec = ocean["ysst"].shape[0]
    if ice["xheat"].shape[0] != nrec:
        raise SystemExit("the ocean and ice streams hold different record counts; "
                         "they are written by the same `nout` and should not")

    nlat, nlon = ocean["ysst"].shape[1:]
    weights = np.broadcast_to(leggauss(nlat)[1][::-1][:, None],
                              (nlat, nlon)) / (2.0 * nlon)

    orbits = list(range(first, last + 1))
    binned: dict[str, np.ndarray] = {}
    fields = ["ts", "hfns", "rss", "rls", "hfss", "hfls", "snm", "sic", "snd",
              "mld", "lsm"]
    for name in fields:
        stack = []
        for index in orbits:
            path = run_dir / f"MOST.{index:05d}.nc"
            if not path.is_file():
                raise SystemExit(f"missing annual output {path}")
            with Dataset(path) as nc:
                stack.append(np.asarray(nc[name][:], dtype=float))
        binned[name] = np.concatenate(stack)
    lsm = binned["lsm"][0]

    # The grid is shared, not reconstructed: the streams carry the land mask and
    # it must be the same array the regular output carries, index for index.
    # CLAUDE.md rule 3 -- if these disagree, nothing below means anything.
    for name, field in (("ice", ice["xls"][0]), ("ocean", ocean["yls"][0])):
        if not np.array_equal(field, lsm):
            raise SystemExit(f"the {name} stream's land mask is not the regular "
                             "stream's; the two are not on the same grid")

    for name, why in MUST_BE_ZERO.items():
        peak = float(np.abs(ocean[name]).max())
        if peak != 0.0:
            raise SystemExit(f"{name} is not identically zero (peak {peak}), but "
                             f"this closure assumes {why}")

    is_ocean = lsm < 0.5
    # Strict: no ice and no snow in ANY stream record. The same mask taken from
    # the binned output admits cells that froze between two output records.
    strict = (is_ocean & (ice["xicec"].max(0) <= 0) & (ice["xiced"].max(0) <= 0)
              & (ice["xsnow"].max(0) <= 0))
    loose = (is_ocean & (binned["sic"].max(0) <= 0) & (binned["snd"].max(0) <= 0))

    atm = binned["rss"] + binned["rls"] + binned["hfss"] + binned["hfls"]
    alignment = align(ice["xheat"], atm, weights, strict)
    per = alignment["records_per_bin"]
    offset = alignment["chosen_offset"]
    nbin = atm.shape[0]
    # The alignment is established, or it is not and the cross-stream identities
    # are not computed. The criteria, stated before the scan is looked at: the
    # chosen offset must beat the worst by a factor of five, and must agree to
    # better than one percent of the seasonal swing. Two accumulations of the
    # same physical quantity that do not tile are not two views of one orbit.
    swing = max(abs(_mean(binned["hfns"][k], weights, strict))
                for k in range(1, nbin))
    worst = max(s["score_w_m2"] for s in alignment["scan"])
    alignment["seasonal_scale_w_m2"] = float(swing)
    alignment["established"] = bool(
        alignment["score_w_m2"] * 5.0 <= worst
        and alignment["score_w_m2"] <= 0.01 * swing)
    if not alignment["established"]:
        raise SystemExit(
            "the ice stream does not tile the regular output for this window: "
            f"best offset scores {alignment['score_w_m2']:.3f} W/m2 against a "
            f"worst of {worst:.3f} and a seasonal scale of {swing:.3f}. The "
            "streams are truncated at every model call, so this usually means "
            "the surviving stream is not the orbit whose MOST file was read.")

    # The window is bins 1 to nbin-1: bin 0 covers a different number of model
    # timesteps from the rest (see the module docstring), so it is the one bin
    # for which the two streams cannot describe the same span.
    span = slice(offset + per, offset + per * nbin)
    nspan = per * (nbin - 1)
    mld = float(np.median(binned["mld"][:, is_ocean]))
    derived = manifest.get("derived_parameters", {})
    seconds = derived.get("orbital_year_seconds")
    steps_per_orbit = derived.get("runsteps_per_orbit")
    if not seconds or not steps_per_orbit:
        raise SystemExit("the manifest does not give orbital_year_seconds and "
                         "runsteps_per_orbit, so the stream's record interval "
                         "cannot be established")
    step_seconds = seconds / steps_per_orbit
    total_steps = steps_per_orbit * len(orbits)
    # The stream writes every `nout` timesteps. Its value is not in the manifest,
    # so take it from the record count and check it against the tiling: `per`
    # stream records tile one output bin, so `per * nout` is the bin's span in
    # timesteps and eleven of those must fit inside the orbit with room to spare
    # for the first bin, which is longer.
    nout = int(round(total_steps / nrec))
    ordinary_bin_steps = per * nout
    first_bin_steps = total_steps - (nbin - 1) * ordinary_bin_steps
    if not 0 < ordinary_bin_steps <= first_bin_steps:
        raise SystemExit(
            f"the tiling implies {ordinary_bin_steps} timesteps per bin and "
            f"{first_bin_steps} in the first, which is not a configuration this "
            "closure understands")
    record_seconds = nout * step_seconds
    bin_seconds = ordinary_bin_steps * step_seconds

    storage = (CRHOS * CPS * mld
               * (ocean["ysst"][offset + per * nbin - 1] - ocean["ysst"][offset + per - 1])
               / (nspan * record_seconds))

    result = {}
    for name, mask in (("strict", strict), ("loose", loose),
                       ("ocean", is_ocean)):
        atm_m = atm[1:].mean(0)
        hfns_m = binned["hfns"][1:].mean(0)
        snm_m = ALF * RHO_WATER * binned["snm"][1:].mean(0)
        xheat = ice["xheat"][span].mean(0)
        xcflux = ice["xcflux"][span].mean(0)
        yheat = ocean["yheat"][span].mean(0)
        xsmelt = ice["xsmelt"][span].mean(0)
        result[name] = {
            "area_fraction": float((mask * weights).sum()),
            "terms_w_m2": {
                "atmosphere_four_flux": _mean(atm_m, weights, mask),
                "atmosphere_hfns": _mean(hfns_m, weights, mask),
                "atmosphere_snow_fusion": _mean(snm_m, weights, mask),
                "ice_stream_xheat": _mean(xheat, weights, mask),
                "ice_stream_xsmelt": _mean(xsmelt, weights, mask),
                "ice_stream_xcflux": _mean(xcflux, weights, mask),
                "ocean_stream_yheat": _mean(yheat, weights, mask),
                "slab_storage": _mean(storage, weights, mask),
            },
            "identities_w_m2": {
                "D1_handoff": _mean(xheat - atm_m, weights, mask),
                "D2_ice_to_ocean": _mean(yheat - xcflux, weights, mask),
                "D3_snow_into_water": _mean((xheat - xcflux) - xsmelt, weights, mask),
                "D4_slab_integration": _mean(storage - yheat, weights, mask),
                "D5_hfns_definition": _mean((atm_m - hfns_m) - snm_m, weights, mask),
                "D6_surface_residual": _mean(hfns_m - storage, weights, mask),
            },
            "seasonal_swing_w_m2": {
                "min": float(min(_mean(binned["hfns"][k], weights, mask)
                                 for k in range(1, nbin))),
                "max": float(max(_mean(binned["hfns"][k], weights, mask)
                                 for k in range(1, nbin))),
            },
        }

    # The first bin's span, from the tiling, and what its equal weight costs.
    b = np.array([_mean(binned["hfns"][k], weights, strict) for k in range(nbin)])
    coefficient = 1.0 / nbin - first_bin_steps / total_steps
    first_bin = {
        "stream_interval_steps": nout,
        "ordinary_bin_steps": ordinary_bin_steps,
        "first_bin_steps": first_bin_steps,
        "first_bin_excess": first_bin_steps / ordinary_bin_steps,
        "coefficient": coefficient,
        "bin0_minus_rest_w_m2": float(b[0] - b[1:].mean()),
        "error_in_annual_mean_w_m2": float(coefficient * (b[0] - b[1:].mean())),
        "statement": "pyburn weights the 12 bins equally; the first covers more "
                     "model timesteps than the rest, so an annual mean from a "
                     "binned file carries this error",
    }

    return {
        "manifest": manifest,
        "window": {"first_orbit": first, "last_orbit": last,
                   "orbits": len(orbits), "stream_records": nrec,
                   "record_seconds": record_seconds, "bin_seconds": bin_seconds,
                   "mixed_layer_depth_m": mld},
        "alignment": alignment,
        "by_mask": result,
        "first_bin_weight": first_bin,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output", type=Path, default=ANALYSIS / "ocean_energy")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    manifest_path = run_dir / "run_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest_path.is_file() else {})
    for name in ("ocean_output", "ice_output"):
        if not (run_dir / name).is_file():
            raise SystemExit(
                f"{run_dir/name} does not exist. The ocean and ice streams are "
                "written whenever NOCEAN and NICE are on and noutput > 0; a run "
                "without them cannot be closed this way.")
    print(f"# ocean and ice streams cover only orbit {final_orbit(manifest)}, "
          "the run's last; every earlier call's stream was truncated by the "
          "next call's open()")

    result = close_ocean(run_dir)
    report = {
        "schema_version": 1,
        "generator": "exoplasim/scripts/close_ocean_energy.py",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "run_dir": str(run_dir),
        "run_id": manifest.get("run_id"),
        "source_build": manifest.get("source_build"),
        "config_sha256": manifest.get("config_sha256"),
        "geography": manifest.get("physical", {}).get("geography"),
        "window": result["window"],
        "alignment": result["alignment"],
        "by_mask": result["by_mask"],
        "first_bin_weight": result["first_bin_weight"],
        "executable_sha256": manifest.get("executable", {}).get("sha256"),
        "software": manifest.get("software"),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    # Named by the run AND the window, as the other two closures are: the window
    # is a parameter of the measurement, and here it is not even a choice.
    path = (args.output /
            f"{run_dir.name}_ocean_energy_{result['window']['first_orbit']}"
            f"-{result['window']['last_orbit']}.json")
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"alignment": {k: v for k, v in report["alignment"].items()
                                    if k != "scan"},
                      "strict": report["by_mask"]["strict"],
                      "first_bin_weight": report["first_bin_weight"],
                      "report": str(path)}, indent=2))


if __name__ == "__main__":
    main()
