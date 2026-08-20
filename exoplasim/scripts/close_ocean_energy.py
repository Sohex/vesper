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

WHICH ORBITS THE INSTRUMENT COVERS depends on when the run was made, and this
script reads what is there rather than assuming either. `oceanmod.f90:331` and
`icemod.f90:429` open these files with a bare
`open(unit,file=...,form='unformatted')`, which truncates, so every model call
discarded the previous call's stream and wrote only its own. CLIM-12 moved each
call's stream aside to `MOST_OCEAN.NNNNN` and `MOST_ICE.NNNNN` before the next
call starts, the way `MOST.NNNNN` was always handled -- 200 MB an orbit and no
model time -- so a run made after it carries the whole block and the closure is
an exact endpoint difference rather than a trend fitted through one orbit.

A run made BEFORE it has the single in-place file holding its final orbit alone.
That is read, not refused, because one orbit is what that run has and the
identities still close on it; the coverage is printed and recorded in the report
so no result is quoted as a block when it is a single orbit. What must never
happen is the third case -- reading a single file as though it were the block --
and it cannot, because the orbit list comes from the filenames.

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

THE BINS DO NOT ALL HOLD THE SAME NUMBER OF RECORDS, and this reports what that
costs an annual mean. pyburn reduces the raw stream with
`np.linspace(0, ntimes, nbin+1).astype(int)`, which spreads the remainder of
`ntimes / nbin` across the bins by truncation, and then averages the bins with
equal weight. At `NLOWIO = 0` there are 182 records per orbit, so two bins in
twelve hold sixteen and the other ten hold fifteen. `lib/climatology.py` derives
the weights and this reports them as `bin_weights`. It is a diagnostic, not a
correction applied to anything here.

An earlier version of this block claimed the FIRST bin was long -- 570 timesteps
against 480 -- and explained it by a model call restoring the previous call's
partial output interval. Both were wrong, and the error was in this file: the
span was derived as `total_steps - 11 * ordinary_bin_steps`, which assumes
eleven equal bins and therefore forces every discrepancy onto the twelfth. The
bin centres say bin 0 holds the FEWEST records, and the restore mechanism it
invoked belongs to `NLOWIO = 1`, which is the regime whose 36 records per orbit
bin exactly. `TASKS.md` CLIM-13 has the measurement.

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
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS  # noqa: F401  (also puts lib/ on sys.path)

import climatology  # noqa: E402  from lib/, via _paths

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


def read_service(paths: Sequence[Path], keep: dict[int, str]) -> dict[str, np.ndarray]:
    """Every record of one or more PlaSim service-format streams, (time, lat, lon).

    Each field is two Fortran unformatted records: an 8-integer header whose
    first entry is the field code, then NLON*NLAT float32 values. Read one
    record at a time and keep only the requested codes, so the working set is
    the fields asked for rather than the whole file.

    `paths` is a sequence because a run's stream is now one file per orbit and
    they concatenate in orbit order. Pass them already sorted; this does not
    re-sort, because the caller knows the orbit numbering and this does not.
    """
    out: dict[int, list[np.ndarray]] = {}
    shape = None
    for path in paths:
        size = os.path.getsize(path)
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
        raise SystemExit(f"{[str(p) for p in paths]} does not carry {missing}")
    counts = {len(v) for v in out.values()}
    if len(counts) != 1:
        raise SystemExit(f"unequal record counts per field across "
                         f"{[str(p) for p in paths]}: {counts}")
    return {keep[code]: np.stack(v).astype(float) for code, v in out.items()}


# The per-call names the wrapper moves the streams to, and the in-place names
# the model writes them under while a call is running.
STREAM_NAMES = {"ocean": ("MOST_OCEAN", "ocean_output"),
                "ice": ("MOST_ICE", "ice_output")}


def stream_files(run_dir: Path, which: str) -> tuple[list[Path], list[int]]:
    """The stream's files in orbit order, and the orbits they cover.

    Two layouts exist and the difference is CLIM-12. `oceanmod.f90:331` and
    `icemod.f90:429` open these files without `position='append'`, so the model
    truncates them at every call; the wrapper now moves each call's stream to
    `MOST_OCEAN.NNNNN` / `MOST_ICE.NNNNN` before the next call starts, exactly
    as it already did for `MOST.NNNNN`. A run made before that patch has only
    the single in-place file, holding its FINAL ORBIT ALONE, and is read that
    way rather than refused -- but the orbit list says so, and the caller states
    the coverage rather than assuming a block.
    """
    prefix, legacy = STREAM_NAMES[which]
    per_call = sorted(run_dir.glob(f"{prefix}.[0-9][0-9][0-9][0-9][0-9]"))
    if per_call:
        return per_call, [int(p.suffix[1:]) for p in per_call]
    single = run_dir / legacy
    if single.is_file():
        return [single], []
    raise SystemExit(
        f"{run_dir} carries neither {prefix}.NNNNN nor {legacy}. The ocean and "
        "ice streams are written whenever NOCEAN and NICE are on and "
        "noutput > 0; a run without them cannot be closed this way.")


def final_orbit(manifest: dict) -> int:
    """The run's last orbit, from the manifest.

    Before CLIM-12 this was the only orbit the streams held, because the
    truncating `open` discarded the previous call's at every call. It is now
    just the end of the range; `stream_files` reports what is actually present.
    """
    segments = manifest.get("segments", [])
    if not segments:
        raise SystemExit("the run manifest records no segments, so which orbit "
                         "the ocean stream covers cannot be established")
    return int(segments[-1]["end_year_index"])


def _bin_mean(field: np.ndarray, records_per_bin: np.ndarray) -> np.ndarray:
    """Time mean over bins 1..n-1, weighted by how many records each bin holds.

    Bin 0 is dropped because it straddles the restart, which is the same
    exclusion `align` makes. The remaining bins are weighted by their record
    counts rather than equally: pyburn splits an orbit with
    `linspace(...).astype(int)`, so the bins are not the same length, and
    treating them as if they were is CLIM-13.
    """
    w = np.asarray(records_per_bin, dtype=float)[1:]
    w = w / w.sum()
    return np.tensordot(w, field[1:], axes=(1 if field.ndim == 1 else 0, 0))


def align(stream: np.ndarray, binned: np.ndarray, weights: np.ndarray,
          mask: np.ndarray, bin_edges: np.ndarray | None = None) -> dict:
    """Find the offset at which stream records tile the binned output.

    `nbin` bins of `per` stream records each; the offset is how many records
    precede the first bin. Returns the scan so the discrimination is visible.
    """
    nbin = binned.shape[0]
    nrec = stream.shape[0]
    # FRACTIONAL, and that is the whole correction. This was
    # `per = nrec // nbin` with a reshape, which requires a whole number of
    # stream records per output bin. Vesper's year is 182.801 days, so a daily
    # stream writes 183 records in most orbits and 182 about every fifth, and
    # 183/12 is 15.25. Flooring to 15 leaves each bin a quarter-record short of
    # where the next one starts, the deficit accumulates across the orbit, and
    # the residual grows LINEARLY with the window -- measured at about
    # 0.42 W/m2 per orbit, which is how a ten-orbit window scored 4.2 against a
    # criterion of 0.26. It went unseen because the closure was validated on a
    # single pre-CLIM-12 orbit, where one orbit's worth of drift is small.
    #
    # Nothing about the model is wrong here and nothing about the streams is:
    # an orbit that is not a whole number of days simply cannot be cut into
    # twelve whole numbers of records.
    per = nrec / nbin
    # BIN WIDTHS COME FROM PYBURN, not from dividing by twelve. Its bins hold
    # different numbers of raw records -- [15,15,15,15,15,16,...] for a
    # 182-record orbit -- which is CLIM-13's finding, and `lib/climatology.py`
    # exists to carry exactly that arithmetic. This alignment never asked it,
    # and equal-width bins put the boundaries in the wrong place: 0.353 W/m2
    # against 0.191 with the real widths, on the same window.
    base = (bin_edges if bin_edges is not None
            else per * np.arange(nbin + 1, dtype=float))
    scan = []
    for offset in range(max(1, nrec - int(base[-1] - base[0]) + 1)):
        # Each bin covers a half-open interval of record index, and a record
        # straddling a boundary is split between the two bins by overlap. With
        # an integer `per` this reduces exactly to the old reshape-and-mean.
        edges = offset + base
        w = np.clip(np.minimum(edges[1:, None], np.arange(nrec)[None, :] + 1.0)
                    - np.maximum(edges[:-1, None], np.arange(nrec)[None, :]),
                    0.0, None)
        rowsum = w.sum(1, keepdims=True)
        if np.any(rowsum <= 0):
            continue
        w = w / rowsum
        tiled = np.tensordot(w, stream, axes=(1, 0))
        diff = np.array([_mean(tiled[k] - binned[k], weights, mask)
                         for k in range(nbin)])
        scan.append({"offset": offset,
                     "per_bin_w_m2": [float(x) for x in diff],
                     # bin 0 straddles the restart and is excluded from the score
                     "score_w_m2": float(np.abs(diff[1:]).mean())})
    best = min(scan, key=lambda s: s["score_w_m2"])
    return {"records_per_bin": float(per), "chosen_offset": best["offset"],
            "score_w_m2": best["score_w_m2"], "scan": scan}


def _mean(field: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> float:
    """Gaussian-weighted mean over the mask, per unit area of the mask."""
    area = float((mask * weights).sum())
    if area <= 0.0:
        return float("nan")
    return float((field * weights * mask).sum()) / area


def close_ocean(run_dir: Path, first_orbit=None, last_orbit=None) -> dict:
    manifest_path = run_dir / "run_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest_path.is_file() else {})
    ocean_paths, ocean_orbits = stream_files(run_dir, "ocean")
    ice_paths, ice_orbits = stream_files(run_dir, "ice")
    if ocean_orbits != ice_orbits:
        raise SystemExit(f"the ocean stream covers orbits {ocean_orbits} and the "
                         f"ice stream {ice_orbits}; they are written by the same "
                         "model call and must cover the same block")
    if ocean_orbits:
        # A WINDOW, because "all the orbits there are" is not always the block
        # worth closing. A run seeded with --restart-from opens on a state taken
        # under different physics, and a run still relaxing is not in the steady
        # state these identities assume. Both put records into the scan that
        # describe a different world from the rest, and the alignment is the
        # first thing to fail when they do. Defaults to everything, so the
        # unwindowed call is unchanged.
        if first_orbit is not None or last_orbit is not None:
            lo = ocean_orbits[0] if first_orbit is None else int(first_orbit)
            hi = ocean_orbits[-1] if last_orbit is None else int(last_orbit)
            keep = [o for o in ocean_orbits if lo <= o <= hi]
            if not keep:
                raise SystemExit(f"no stream orbits in [{lo}, {hi}]; the run has "
                                 f"{ocean_orbits[0]}-{ocean_orbits[-1]}")
            ocean_paths = [p for p, o in zip(ocean_paths, ocean_orbits) if lo <= o <= hi]
            ice_paths = [p for p, o in zip(ice_paths, ice_orbits) if lo <= o <= hi]
            ocean_orbits = keep
        first, last = ocean_orbits[0], ocean_orbits[-1]
        if ocean_orbits != list(range(first, last + 1)):
            raise SystemExit(f"the stream orbits {ocean_orbits} have a gap; this "
                             "closure integrates a contiguous block")
    else:
        # Pre-CLIM-12 run: one in-place file holding the final orbit alone.
        first = last = final_orbit(manifest)

    ocean = read_service(ocean_paths, OCEAN_CODES)
    ice = read_service(ice_paths, ICE_CODES)
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
    with Dataset(run_dir / f"MOST.{orbits[0]:05d}.nc") as nc:
        binned_time = np.asarray(nc["time"][:], dtype=float)
    # NOT bin 0. `align` below already excludes it because it straddles the
    # restart, and the mask is where that bites hardest: under NLOWIO = 1 the
    # first record of a RUN carries a partial accumulation window, so a field
    # that cannot vary comes back scaled -- the land mask reads 0.953 instead of
    # 1 in orbit 0 bin 0 and exactly 1 everywhere else. Taking the maximum over
    # bins is exact for a binarised mask (oceanmod.f90 hard-binarises yls, so
    # the only values are 0 and 1) and does not care which bin is the partial
    # one. Reading bin 0 here while distrusting it eight lines down was the
    # script disagreeing with itself, and it made low-I/O output unreadable.
    lsm = binned["lsm"].max(0)

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
    # One orbit's bin counts, repeated per orbit and scaled to the stream's own
    # record count. Scaling absorbs the calendar: the year is 182.801 days, so
    # the stream writes 183 records in most orbits and 182 in some, while
    # pyburn's counts describe the orbit it binned.
    per_orbit = climatology.counts_for(
        climatology.infer_ntimes(binned_time), len(binned_time))
    # THE TWO ARE WRITTEN ON DIFFERENT CADENCES AND ONLY COINCIDE AT NWPD = 1.
    # The regular output writes every `nafter = mtspd / nwpd` steps; the ocean
    # and ice streams follow their own counters and ignore nwpd entirely. At
    # nwpd = 8 the regular output holds 1462 raw records an orbit and the
    # streams still hold 182, so the tiling below stops describing anything.
    # Reported rather than silently rescaled: a closure whose two sides sample
    # at different rates is measuring the rates.
    raw_per_orbit = int(per_orbit.sum())
    stream_per_orbit = ice["xheat"].shape[0] / max(len(orbits), 1)
    # A RATIO, not a difference. The two legitimately disagree by up to one
    # record an orbit -- the year is not a whole number of write intervals, so
    # the stream holds 183 in most orbits and 182 in about every fifth, and a
    # mean over several orbits lands between them. What this is for is the
    # eightfold gap a changed NWPD opens, so 10% separates the two cleanly.
    if not 0.9 < stream_per_orbit / max(raw_per_orbit, 1) < 1.1:
        raise SystemExit(
            f"the regular output holds {raw_per_orbit} raw records an orbit and "
            f"the streams hold {stream_per_orbit:.1f}; they are written on "
            "different cadences, which happens whenever NWPD is not 1, and the "
            "stream cannot tile bins it does not share a clock with. Re-run at "
            "NWPD = 1, or compare the annual means directly instead.")
    # Records per bin across the WHOLE window, which the identity means below
    # weight by. One orbit's counts repeated per orbit.
    bin_records = np.tile(per_orbit, len(orbits)).astype(float)
    bin_edges = np.concatenate([[0.0], np.cumsum(np.tile(per_orbit, len(orbits)))])
    bin_edges = bin_edges * (ice["xheat"].shape[0] / bin_edges[-1])
    alignment = align(ice["xheat"], atm, weights, strict, bin_edges)
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
    # TWO CRITERIA, and one of them is not always evaluable. The five-fold test
    # discriminates AMONG candidate offsets, and it existed because integer
    # tiling left spare records so several offsets were possible and one had to
    # win clearly. With the bin edges taken from pyburn the bins consume every
    # record, there is exactly one candidate, and `best` and `worst` are the same
    # number -- so the test is INAPPLICABLE rather than failed, and reporting it
    # as failed would mean the closure can never establish anything, which is
    # not a criterion but a wall.
    #
    # It is recorded rather than quietly dropped: `discrimination` says whether
    # it was evaluated, so a reader cannot mistake one criterion for two. The
    # hazard it guarded -- a stream file that is not the orbit whose MOST file
    # was read -- is now structurally impossible: each orbit has its own
    # numbered stream and `stream_files` refuses a gap.
    candidates = len(alignment["scan"])
    discriminated = candidates > 1 and alignment["score_w_m2"] * 5.0 <= worst
    alignment["candidate_offsets"] = candidates
    alignment["discrimination"] = (
        "not applicable: one candidate offset, so there is nothing to "
        "discriminate among" if candidates == 1
        else f"chosen offset beats the worst by "
             f"{worst / max(alignment['score_w_m2'], 1e-12):.1f}x, needs 5x")
    alignment["established"] = bool(
        (discriminated or candidates == 1)
        and alignment["score_w_m2"] <= 0.01 * swing)
    if not alignment["established"]:
        raise SystemExit(
            "the ice stream does not tile the regular output for this window: "
            f"best offset scores {alignment['score_w_m2']:.3f} W/m2 against a "
            f"worst of {worst:.3f} and a seasonal scale of {swing:.3f}. The "
            "streams are truncated at every model call, so this usually means "
            "the surviving stream is not the orbit whose MOST file was read.")

    # The window is bins 1 to nbin-1. Bin 0 is dropped as conservatism rather
    # than because it is anomalous -- the measured spans put it at the SHORTEST,
    # not the longest -- but it is the bin whose leading edge the uniform tiling
    # is least able to place, so the cross-stream identities skip it.
    # The records of bins 1..nbin-1, i.e. everything but the first bin, which
    # straddles the restart. Taken from the BIN EDGES because `per` is
    # fractional once the widths are pyburn's: `offset + per*k` is no longer a
    # record index, and slicing with it silently meant something else.
    span = slice(int(round(offset + bin_edges[1])),
                 int(round(offset + bin_edges[-1])))
    nspan = float(bin_edges[-1] - bin_edges[1])
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
    ordinary_bin_steps = float(per) * nout
    if ordinary_bin_steps <= 0:
        raise SystemExit(f"the tiling implies {ordinary_bin_steps} timesteps per "
                         "bin, which is not a configuration this closure understands")
    record_seconds = nout * step_seconds
    bin_seconds = ordinary_bin_steps * step_seconds

    # Endpoint records from the BIN EDGES, not from `offset + per*k`. `per` is
    # fractional now -- the bins are pyburn's and have different widths -- so
    # arithmetic on it is not a record index. These are the last record of the
    # final bin and the last record of the first, which is the span the storage
    # is differenced over and is unchanged in meaning.
    last_record = int(round(offset + bin_edges[-1])) - 1
    end_of_first_bin = int(round(offset + bin_edges[1])) - 1
    nrec_total = ocean["ysst"].shape[0]
    last_record = min(max(last_record, 0), nrec_total - 1)
    end_of_first_bin = min(max(end_of_first_bin, 0), nrec_total - 1)
    storage = (CRHOS * CPS * mld
               * (ocean["ysst"][last_record] - ocean["ysst"][end_of_first_bin])
               / (nspan * record_seconds))

    result = {}
    for name, mask in (("strict", strict), ("loose", loose),
                       ("ocean", is_ocean)):
        # WEIGHTED BY RECORDS PER BIN, not `.mean(0)`. The stream terms below
        # are record means and are correctly weighted by construction; the
        # atmosphere terms are BIN means, and pyburn's bins hold different
        # numbers of records, so equal-weighting them is the CLIM-13 error.
        # This file diagnoses that error two hundred lines down -- it reports
        # hfns at -0.4195 equal-weighted against -0.2819 true-weighted, a
        # difference of +0.1376 W/m2 -- and then computed these identities with
        # the equal-weighted mean anyway. Most of D1 was that.
        atm_m = _bin_mean(atm, bin_records)
        hfns_m = _bin_mean(binned["hfns"], bin_records)
        snm_m = ALF * RHO_WATER * _bin_mean(binned["snm"], bin_records)
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

    # What pyburn's equal weighting of the bins costs, MEASURED from the bin
    # centres the binned file carries rather than inferred from the tiling.
    #
    # This block used to derive the first bin's span as the residual
    # `total_steps - 11 * ordinary_bin_steps` and report a 1.19x excess on bin 0.
    # That was an artefact of its own construction: assuming eleven equal bins
    # forces every mismatch onto the twelfth. The centres say otherwise -- the
    # spans run 480 to 496 steps with bin 0 at the MINIMUM -- and pyburn's
    # binning cannot produce a long first bin anyway, because it splits the raw
    # records with `np.linspace(0, ntimes, nbin+1).astype(int)` and divides each
    # bin by its own count. See TASKS.md CLIM-13 for the correction.
    # ONE ORBIT's worth of bins, because `counts` and `w` below describe one
    # orbit and `nbin` is now the whole window -- 36 bins over three orbits, not
    # 12. This block reports the equal-weight error pyburn's binning produces,
    # which is a per-orbit property, so it is computed on the first orbit rather
    # than on the concatenation. It broadcast fine while every window was a
    # single orbit and raised the moment one was not.
    # PER ORBIT throughout. `binned_time` holds one orbit's bin centres, so the
    # bin count here is len(binned_time) and NOT `nbin`, which is the whole
    # window -- 36 bins over three orbits. The equal-weight error this block
    # reports is a property of how pyburn splits ONE orbit, so mixing the two
    # counts asks counts_for for a split that was never made.
    nbin_orbit = len(binned_time)
    b = np.array([_mean(binned["hfns"][k], weights, strict)
                  for k in range(nbin_orbit)])
    ntimes = climatology.infer_ntimes(binned_time)
    counts = climatology.counts_for(ntimes, nbin_orbit)
    w = climatology.bin_weights(binned_time)
    equal = 1.0 / nbin_orbit
    bin_weights = {
        "method": "records per bin from pyburn's own binning arithmetic, with "
                  "the record count recovered from the bin centres and checked "
                  "against them. lib/climatology.py.",
        "raw_records_per_orbit": int(ntimes),
        "records_per_bin": [int(x) for x in counts],
        "weights": [float(x) for x in w],
        "equal_weight": equal,
        "max_deviation_fraction": float(np.abs(w - equal).max() / equal),
        "bins_per_orbit": int(nbin_orbit),
        "longest_bin": int(np.argmax(counts)),
        "shortest_bin": int(np.argmin(counts)),
        "hfns_equal_weight_w_m2": float(b.mean()),
        "hfns_true_weight_w_m2": float((b * w).sum()),
        "error_in_annual_mean_w_m2": float((b * w).sum() - b.mean()),
        "statement": "pyburn weights the bins equally and they hold different "
                     "numbers of raw records, because linspace().astype(int) "
                     "spreads the remainder by truncation. The error is the "
                     "difference above. It is NOT in the first bin, which holds "
                     "the fewest records.",
    }

    return {
        "manifest": manifest,
        "window": {"first_orbit": first, "last_orbit": last,
                   "orbits": len(orbits), "stream_records": nrec,
                   "record_seconds": record_seconds, "bin_seconds": bin_seconds,
                   "mixed_layer_depth_m": mld,
                   # Which layout this was read from, so a single-orbit result
                   # can never be quoted as a block. CLIM-12.
                   "stream_layout": "per_call" if ocean_orbits else "final_orbit_only",
                   "stream_orbits": ocean_orbits or [last]},
        "alignment": alignment,
        "by_mask": result,
        "bin_weights": bin_weights,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output", type=Path, default=ANALYSIS / "ocean_energy")
    parser.add_argument("--first", type=int, default=None,
                        help="first orbit of the window. Default: the earliest "
                             "the streams cover. Use it to drop the opening "
                             "orbits of a seeded run, which follow a restart "
                             "taken under different physics")
    parser.add_argument("--last", type=int, default=None,
                        help="last orbit of the window; default the latest")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    manifest_path = run_dir / "run_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest_path.is_file() else {})
    _, orbits = stream_files(run_dir, "ocean")     # raises if neither layout is present
    stream_files(run_dir, "ice")
    if orbits:
        print(f"# ocean and ice streams cover orbits {orbits[0]}-{orbits[-1]}, "
              f"{len(orbits)} of them, one file per model call")
    else:
        print(f"# ocean and ice streams cover only orbit {final_orbit(manifest)}, "
              "the run's last: this run predates CLIM-12 and every earlier "
              "call's stream was truncated by the next call's open()")

    result = close_ocean(run_dir, args.first, args.last)
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
        "bin_weights": result["bin_weights"],
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
                      "bin_weights": report["bin_weights"],
                      "report": str(path)}, indent=2))


if __name__ == "__main__":
    main()
