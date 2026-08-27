#!/usr/bin/env python3
"""Turn the offline dust burden into the boundary field the model reads. DUST-11.

    python exoplasim/scripts/build_surface_dust.py
    python exoplasim/scripts/build_surface_dust.py --variant low   # shelter bracket

`exoplasim-3.4.2-prescribed-dust.patch` adds surface code 1811, `ddustcol`, and
reads it in `radini` as the BAND 1 (0.34-0.75 um) column extinction optical
depth. This writes that field from `aeolian/analysis/dust_baseline.nc`, and
writes beside it the two constants the run's namelist needs and nothing else
knows: the thermal-infrared absorption ratio `dustqlw` and the scale height
`dusthsc`.

## Why band 1 and not "the" optical depth

The offline chain reports a single optical depth at a flux-weighted mass
extinction efficiency across both of ExoPlaSim's shortwave bands, and the model
needs one per band. Converting once here, rather than in the Fortran, keeps the
band split with the optics that produced it:

    tau_band = column_mass * MEE_band,   column_mass = aod_chain / MEE_chain

so `tau_1 = aod_chain * MEE_1 / MEE_chain`, and the model derives band 2 from
the aerofile's own `qex2/qex1`. That ratio is the same quantity as `MEE_2/MEE_1`
because `apart` and `rhop` cancel out of it, which is checked here rather than
asserted: if the aerofile and the optics ever stop agreeing, this refuses to
write.

## The longwave ratio

`dustqlw` is the Planck-weighted thermal-infrared mass ABSORPTION efficiency
divided by the band-1 mass extinction efficiency, so multiplying the prescribed
band-1 optical depth by it gives the grey absorption optical depth the patched
longwave solver wants.

Planck-weighted rather than flat, because dust absorption peaks in the silicate
band near 9-10 um, which is where a 290 K Planck function also peaks, and their
product is what a flux cares about. `analysis/dust_forcing.json` reports an
UNWEIGHTED mean over 4-40 um for a different purpose; the two are not the same
number and should not be read across.

## What this does NOT do

Nothing here responds to the climate. The field is a fixed annual mean, so the
run cannot show a dust-wind feedback and is not asked to. That is the whole
premise of DUST-11 being cheap: the question is what a given burden does.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, INPUTS, PROJECT_ROOT  # noqa: E402  (puts lib/ on sys.path)
import climatology  # noqa: E402  from lib/, put on sys.path by _paths
from paths import best_available_climatology, rel, require_configured_grid  # noqa: E402
from provenance import input_stamp
from sra import write_sra

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dust_indices import check_coverage, indices, selection  # noqa: E402
from mie_dust import lognormal_integrate  # noqa: E402

DUST_CODE = 1811

# The thermal band the absorption ratio is integrated over. The upper limit is
# where the longwave indices stop, not where the physics does; above it k is held
# at its last value, which carries a few percent of a 290 K body's emission.
# WHICH indices those are is declared in aeolian/config/dust.yaml under
# optics.indices.longwave and read through dust_indices, exactly as the two
# shortwave bands are. DUST-12.
LW_LO_UM, LW_HI_UM = 4.0, 40.0
LW_GRID = 160

AEROFILE = PROJECT_ROOT / "exoplasim" / "data" / "dust" / "vesper_dust_aerosol.dat"
OPTICS = PROJECT_ROOT / "analysis" / "dust_optics.json"
DUST_CFG = PROJECT_ROOT / "aeolian" / "config" / "dust.yaml"
DUST_NC = PROJECT_ROOT / "aeolian" / "analysis" / "dust_baseline.nc"
DUST_JSON = PROJECT_ROOT / "aeolian" / "analysis" / "dust_baseline.json"


def planck(lam_um, temperature_k):
    """Spectral radiant exitance, W/m2 per metre of wavelength.

    Same constant as `dust_forcing.py`: c1 is 2*pi*h*c^2, so the pi is already
    in it and this is a flux out of a surface rather than a radiance.
    """
    lam = np.asarray(lam_um, dtype=float) * 1e-6
    c1, c2 = 3.7418e-16, 1.4388e-2
    return c1 / (lam ** 5 * (np.exp(c2 / (lam * temperature_k)) - 1.0))


def check_planck(temperature_k: float, tol: float = 0.02) -> float:
    """Integral of `planck` against sigma T^4. Raises if the normalisation slips."""
    lam = np.exp(np.linspace(np.log(0.1), np.log(500.0), 4000))
    total = float(np.trapezoid(planck(lam, temperature_k), lam * 1e-6))
    ratio = total / (5.670374e-8 * temperature_k ** 4)
    if abs(ratio - 1.0) > tol:
        raise SystemExit(
            f"planck() integrates to {ratio:.4f} of Stefan-Boltzmann at "
            f"{temperature_k} K. The normalisation is wrong; nothing below is usable.")
    return ratio


def lw_mass_absorption(dataset: str, r_mod_um: float, sigma_g: float,
                       rho_g_cm3: float, temperature_k: float,
                       n_grid: int) -> float:
    """Planck-weighted thermal-infrared mass absorption efficiency, m2/g."""
    check_coverage(dataset, LW_LO_UM, LW_HI_UM)
    n_of, k_of = indices(dataset)
    grid = np.linspace(LW_LO_UM, LW_HI_UM, n_grid)
    b = planck(grid, temperature_k)
    kabs = np.empty(n_grid)
    for i, lam in enumerate(grid):
        _, ssa, _, mee = lognormal_integrate(
            lam, n_of(lam), k_of(lam), r_mod_um, sigma_g, 0.01, 25.0,
            rho_g_cm3, n_r=200)
        kabs[i] = mee * (1.0 - ssa)
    return float(np.trapezoid(b * kabs, grid) / np.trapezoid(b, grid))


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--dust", type=Path, default=DUST_NC)
    ap.add_argument("--dust-report", type=Path, default=DUST_JSON)
    ap.add_argument("--variant", default="central", choices=("low", "central", "high"),
                    help="shelter bracket end; `central` is the one the run uses")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="supplies the grid and the mean surface temperature; "
                         "defaults to the best available climatology, which is "
                         "the baseline once config/planet.yaml names one")
    ap.add_argument("--optics", type=Path, default=OPTICS)
    ap.add_argument("--aerofile", type=Path, default=AEROFILE)
    ap.add_argument("--dust-config", type=Path, default=DUST_CFG)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    # THE BEST AVAILABLE, which is the baseline once one is named and the
    # bootstrap before that. The grid and the land mask this takes are the same
    # at either stage, but the longwave ratio is not: it Planck-weights the
    # mass absorption at the global-mean surface temperature, which is the
    # climate STATE. That this writes a staged `.sra` an ExoPlaSim run reads is
    # ITERATION rather than circularity -- the dust the model sees is the dust
    # the previous climate produced -- and the `needs: bootstrap_climatology`
    # edge in config/pipeline.yaml still states what must EXIST for a first
    # pass to run.
    args.climatology, clim_stage = best_available_climatology(args.climatology)

    import netCDF4 as nc

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    # nlat is not needed here: require_configured_grid below checks both.
    nlon = int(model["longitudes"])
    resolution = str(model["resolution"]).upper()

    dust_cfg = yaml.safe_load(args.dust_config.read_text(encoding="utf-8"))
    optics = json.loads(args.optics.read_text(encoding="utf-8"))
    report = json.loads(args.dust_report.read_text(encoding="utf-8"))

    # -- the two band efficiencies, from the indices the chain actually chose --
    sel = selection(dust_cfg)
    band1 = next(r for r in optics["results"]
                 if r["indices"] == sel["band1"] and r["band"] == "band 1")
    band2 = next(r for r in optics["results"]
                 if r["indices"] == sel["band2"] and r["band"] == "band 2")
    mee1 = band1["mass_extinction_efficiency_m2_g"] * 1000.0     # m2/g -> m2/kg
    mee2 = band2["mass_extinction_efficiency_m2_g"] * 1000.0
    f1 = float(optics["stellar_flux_fraction_band1"])
    mee_chain = float(report["mass_extinction_efficiency_m2_kg"])

    # CHECK 1, an identity: the chain's own efficiency is the flux-weighted mean
    # of the two bands. If it is not, one of the three files has moved and the
    # conversion below would silently rescale every optical depth.
    #
    # The tolerance is 2e-3 and not tighter ON PURPOSE. The two band
    # efficiencies differ by 1.6%, so the band share is a weak lever here: a 2%
    # change in it moves the reconstruction by 0.03%. That share is canonical in
    # `lib/stellar.py`, and this check has to survive it being corrected there
    # before `analysis/dust_optics.json` is regenerated, without becoming so
    # loose that a real disagreement between the files passes.
    rebuilt = f1 * mee1 + (1.0 - f1) * mee2
    if abs(rebuilt / mee_chain - 1.0) > 2e-3:
        raise SystemExit(
            f"the chain reports a mass extinction efficiency of {mee_chain:.2f} "
            f"m2/kg but the band efficiencies rebuild to {rebuilt:.2f}. One of "
            f"{rel(args.optics)}, {rel(args.dust_config)} and "
            f"{rel(args.dust_report)} has moved without the others.")

    # CHECK 2, another identity: the model derives band 2 from the aerofile's
    # ratio of extinction efficiencies, and apart and rhop cancel out of it, so
    # it must equal MEE2/MEE1. This is the one thing linking the boundary field
    # to the optics file the model reads, and it can fail.
    qs = [float(x) for x in args.aerofile.read_text().splitlines()[1:9]]
    qratio = qs[4] / qs[0]
    if abs(qratio / (mee2 / mee1) - 1.0) > 1e-5:
        raise SystemExit(
            f"the aerofile's qex2/qex1 is {qratio:.6f} but the optics give "
            f"{mee2 / mee1:.6f}. Regenerate the aerofile with "
            f"exoplasim/scripts/dust_aerofile.py before writing this field.")

    # -- the grid, taken from the climatology and asserted, never matched ------
    #
    # Both files were built on the same grid, so the check is equality of the
    # coordinate arrays and not a nearest-neighbour match. Matching by longitude
    # across a boundary is what silently paired zero cells on three scripts.
    with nc.Dataset(args.climatology) as ds:
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
        lsm = np.asarray(ds["lsm"][0], dtype=float) > 0.5
        ts_mean = float(np.average(
            climatology.annual_mean_of(ds, "ts"),
            weights=np.broadcast_to(np.cos(np.deg2rad(lat))[:, None],
                                    (len(lat), len(lon)))))
    # The same guard the resolver applies, called explicitly because
    # `--climatology` can hand this an arbitrary file that never went through
    # the resolver. One expression, in lib/paths.py.
    require_configured_grid(args.climatology, config)

    with nc.Dataset(args.dust) as ds:
        dlat = np.asarray(ds["lat"][:], dtype=float)
        dlon = np.asarray(ds["lon"][:], dtype=float)
        aod = np.asarray(ds[f"aod_{args.variant}"][:], dtype=float)
        dust_terrain = getattr(ds, "terrain_hash", None)
    if not (np.allclose(dlat, lat, atol=1e-9) and np.allclose(dlon, lon, atol=1e-9)):
        raise SystemExit(
            "the dust field and the climatology are on different grids. They "
            "should be the same grid, index for index; do not reconcile them.")

    sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    import builds as _b
    want_terrain = _b.terrain_hash(config)
    if dust_terrain and dust_terrain != want_terrain:
        raise SystemExit(
            f"the dust field is on terrain {dust_terrain[:16]} and the active "
            f"build is {want_terrain[:16]}. Re-run aeolian/scripts/build_dust.py.")

    if not np.isfinite(aod).all() or aod.min() < 0.0:
        raise SystemExit("the dust optical depth field has negatives or NaNs")

    # -- the conversion -------------------------------------------------------
    tau1 = aod * (mee1 / mee_chain)

    # -- the longwave ratio ---------------------------------------------------
    dist = optics["size_distribution"]
    check_planck(ts_mean)
    kabs = lw_mass_absorption(sel["longwave"],
                              float(dist["number_median_radius_um"]),
                              float(dist["sigma_g"]),
                              float(dist["density_g_cm3"]),
                              ts_mean, LW_GRID) * 1000.0          # m2/g -> m2/kg
    kabs_half = lw_mass_absorption(sel["longwave"],
                                   float(dist["number_median_radius_um"]),
                                   float(dist["sigma_g"]),
                                   float(dist["density_g_cm3"]),
                                   ts_mean, LW_GRID // 2) * 1000.0
    convergence = abs(kabs / kabs_half - 1.0)
    if convergence > 0.01:
        raise SystemExit(
            f"the thermal-infrared integral moves {convergence:.1%} between "
            f"{LW_GRID // 2} and {LW_GRID} points; it is not converged")
    dustqlw = kabs / mee1

    scale_height = float(dust_cfg["transport"]["dust_scale_height_m"])

    # -- write ----------------------------------------------------------------
    output = args.output or (INPUTS / resolution.lower()
                             / f"orogen_{resolution}_surf_{DUST_CODE:04d}.sra")
    output.parent.mkdir(parents=True, exist_ok=True)
    write_sra(output, DUST_CODE, tau1)

    weights = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nlon))
    land_mean = float(np.average(tau1[lsm], weights=weights[lsm]))
    global_mean = float(np.average(tau1, weights=weights))

    provenance = {
        "note": "DUST-11 prescribed dust. Surface code 1811 is the BAND 1 column "
                "extinction optical depth; the model derives band 2 from the "
                "aerofile and the thermal infrared from dustqlw below. "
                "Generated by exoplasim/scripts/build_surface_dust.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code": DUST_CODE,
        "field": "ddustcol, band-1 column extinction optical depth, dimensionless",
        "variant": args.variant,
        "source_build": str(config["source_build"]),
        "terrain_hash": want_terrain,
        "dust_field": rel(args.dust),
        "dust_field_sha256": sha256_of(args.dust),
        "dust_report": rel(args.dust_report),
        "dust_report_sha256": sha256_of(args.dust_report),
        "climatology": rel(args.climatology),
        # WHICH STAGE the Planck weighting was taken at. lib/paths.py.
        "climatology_stage": clim_stage,
        "optics": rel(args.optics),
        "optics_sha256": sha256_of(args.optics),
        "optics_indices": sel,
        "aerofile": rel(args.aerofile),
        "aerofile_sha256": sha256_of(args.aerofile),
        "config_sha256": sha256_of(args.config),
        "conversion": {
            "mass_extinction_efficiency_chain_m2_kg": mee_chain,
            "mass_extinction_efficiency_band1_m2_kg": mee1,
            "mass_extinction_efficiency_band2_m2_kg": mee2,
            "stellar_flux_fraction_band1": f1,
            "stellar_flux_fraction_band1_source":
                "analysis/dust_optics.json. The canonical value is lib/stellar.py; "
                "it is used here only in the consistency check above, which is "
                "insensitive to it, and NOT in the conversion, which is a ratio "
                "of two band efficiencies.",
            "band1_over_chain": mee1 / mee_chain,
            "aerofile_qex2_over_qex1": qratio,
            "optics_mee2_over_mee1": mee2 / mee1,
        },
        "namelist_values": {
            "NDUSTRAD": 1,
            "DUSTSC": 1.0,
            "DUSTHSC": scale_height,
            "DUSTQLW": dustqlw,
            # AEROQLW is the SAME RATIO for the INTERACTIVE aerosol, and on this
            # world it is the same number. radmod keeps the two keys apart
            # because in general the transported tracer is whatever aero_nl's
            # apart and rhop describe and need not be the prescribed field's
            # dust; here it is. `dust_aerofile.py` builds the aerofile the
            # interactive path reads by rescaling these same Mie efficiencies to
            # `apart` while holding the mass extinction efficiency invariant --
            # it refuses to write the file if the rescale moves it by more than
            # 1e-9 -- so the band-1 extinction optical depth per unit dust mass
            # is the same on both paths, and this ratio, which is a thermal
            # absorption per unit band-1 extinction, is the same with it.
            # `aerofile_qex2_over_qex1` above is the per-band form of that
            # agreement, checked against the optics to 1e-5.
            #
            # It is written here rather than derived a second time because two
            # derivations of one number is how they drift. world-24v: nothing in
            # the project wrote AEROQLW at all, and radini ABORTS on an
            # interactive aerosol with it left at zero, so the emitted-dust path
            # had no longwave term structurally available.
            "AEROQLW": dustqlw,
        },
        "longwave": {
            "planck_weighted_mass_absorption_m2_kg": kabs,
            "weighting_temperature_k": ts_mean,
            "band_um": [LW_LO_UM, LW_HI_UM],
            "indices": sel["longwave"],
            "grid_points": LW_GRID,
            "grid_convergence": convergence,
            "note": "Planck-weighted, so it is NOT the unweighted 4-40 um mean "
                    "that analysis/dust_forcing.json reports for its own "
                    "purpose. Dust absorption peaks in the silicate band near "
                    "9-10 um, where the Planck function does too, so the "
                    "weighted value is the larger of the two.",
            "k_held_flat_above_um": 40.0,
        },
        "field_statistics": {
            "land_mean_band1_optical_depth": land_mean,
            "global_mean_band1_optical_depth": global_mean,
            "max_band1_optical_depth": float(tau1.max()),
            "land_cells": int(lsm.sum()),
        },
        "requires": "exoplasim/patches/exoplasim-3.4.2-prescribed-dust.patch, "
                    "applied and every binary rebuilt. A binary without it cannot "
                    "parse NDUSTRAD and aborts in radini_, which is the loud "
                    "failure and the one to want.",
        "output": rel(output),
        "output_sha256": sha256_of(output),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    # The derived files this generator read, hashed on the one footing
    # `check_consistency.py` compares. The bespoke `*_sha256` keys above are the
    # same facts in a shape nothing could check; this step is absent from
    # `INERT_CONFIG_KEYS` -- its config trace has not been done -- and the input
    # check does not depend on that set, so it is covered either way.
    provenance.update(input_stamp([
        args.dust, args.dust_report, args.optics, args.aerofile,
        args.dust_config, args.climatology]))

    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"variant            {args.variant}")
    print(f"chain MEE          {mee_chain:.2f} m2/kg, band 1 {mee1:.2f}, "
          f"band 2 {mee2:.2f}")
    print(f"band ratio         aerofile {qratio:.6f} vs optics "
          f"{mee2 / mee1:.6f}  (agree)")
    print(f"land-mean tau_1    {land_mean:.5f}   global {global_mean:.5f}   "
          f"max {tau1.max():.4f}")
    print(f"DUSTQLW            {dustqlw:.6f}  "
          f"(k_abs {kabs:.2f} m2/kg, Planck-weighted at {ts_mean:.2f} K)")
    print(f"DUSTHSC            {scale_height:.0f} m")
    print(f"\nwrote {rel(output)}")
    print(f"      {report_path.name}")


if __name__ == "__main__":
    main()
