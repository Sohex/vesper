#!/usr/bin/env python3
"""Decide which basins overflow, from an ExoPlaSim climatology.

This closes the loop the fork left open. World Orogen measures basin geometry
and refuses to say which basins are endorheic, because that is a water balance.
Here is the water balance.

A basin overflows, and over geological time therefore incises its outlet and
stops being a basin, exactly when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

The right-hand side is pure geometry and ships in `basins.nc` as
`critical_aridity_index`. The left-hand side is climate, integrated over each
basin's catchment through the sparse coupling matrix.

The evaporation problem, stated rather than hidden
--------------------------------------------------

`E` in that expression is evaporation from *open water*. The model does not have
a lake there, so it does not report one. Land evaporation is moisture-limited:
ExoPlaSim scales it by a wetness factor that reaches 1 only when soil water
exceeds 40% of field capacity. Using land evaporation therefore understates E,
understates (E-P)/runoff, and carves too many basins.

Rather than pick one estimate, this brackets it, the same way the albedo
endmembers were bracketed:

  wet    E = land evaporation as reported. A lower bound on open-water
         evaporation, biased toward carving.
  potential  E = land evaporation divided by the model's own wetness factor,
         reconstructed from soil moisture as min(1, mrso / (0.4 * wsmax)).
         This is what the same cell would evaporate if it were saturated, which
         is what a lake is. Biased toward preserving.

Basins that carve under both bounds are robust. Basins that differ between them
are genuinely undecided by this climate and should be flagged, not resolved by
picking whichever bound gives a tidier answer.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, DATA
from lake_balance import BasinSet, carve_verdict, solve

DRHSFULL = 0.4          # landmod.f90: wetness reaches 1 above this fraction
WSMAX_EARTH = 0.5       # landmod.f90 default field capacity, metres


def annual_mean(ds: Dataset, name: str) -> np.ndarray:
    """Annual mean of a (time, lat, lon) field over equal-length bins."""
    return np.asarray(ds[name][:]).mean(axis=0)


def basin_means(coupling: Path, fields: dict[str, np.ndarray], n_basins: int):
    """Catchment-area-weighted mean of each field, per basin."""
    with Dataset(coupling) as ds:
        basin = np.asarray(ds["basin"][:]).astype(np.int64)
        cell = np.asarray(ds["cell"][:]).astype(np.int64)
        area = np.asarray(ds["area_km2"][:])
        nlon = int(ds.n_lon)
    row, col = np.divmod(cell, nlon)
    weight = np.zeros(n_basins)
    np.add.at(weight, basin, area)
    out = {}
    for name, grid in fields.items():
        acc = np.zeros(n_basins)
        np.add.at(acc, basin, area * grid[row, col])
        out[name] = np.where(weight > 0, acc / np.maximum(weight, 1e-30), np.nan)
    return out, weight


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--climatology", type=Path,
                    default=Path("exoplasim/analysis/climatology_s096/"
                                 "baseline_regular_climatology.nc"))
    ap.add_argument("--coupling", type=Path,
                    default=DATA / "coupling_exoplasim-T42.nc")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--output", type=Path, default=ANALYSIS / "carve_verdict.json")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    basins = BasinSet()
    n = basins.n

    with Dataset(args.climatology) as ds:
        pr = annual_mean(ds, "pr")            # m/s
        evap = -annual_mean(ds, "evap")       # code 182 is negative upward
        mrro = annual_mean(ds, "mrro")        # m/s
        mrso = annual_mean(ds, "mrso")        # m
        lsm = annual_mean(ds, "lsm")

    # ExoPlaSim's own wetness factor, reconstructed. Where soil is wet this is 1
    # and land evaporation is already the potential rate.
    wetness = np.clip(mrso / (DRHSFULL * WSMAX_EARTH), 0.0, 1.0)
    potential = np.where(wetness > 1e-3, evap / np.maximum(wetness, 1e-3), evap)
    potential = np.maximum(potential, evap)

    fields = {"pr": pr, "evap": evap, "potential": potential,
              "mrro": mrro, "lsm": lsm}
    means, catch_area = basin_means(args.coupling, fields, n)

    # Units cancel in the aridity index, but keep them physical for the solver.
    year_s = 189.6145 * 86400.0
    to_km_per_year = year_s / 1000.0
    runoff = means["mrro"] * to_km_per_year
    precip = means["pr"] * to_km_per_year

    results = {}
    for label, e_field in (("wet", "evap"), ("potential", "potential")):
        evapo = means[e_field] * to_km_per_year
        with np.errstate(divide="ignore", invalid="ignore"):
            index = np.where(runoff > 0, (evapo - precip) / runoff, np.inf)
        carve, endorheic_km2 = carve_verdict(basins, index, basins.catchment_km2)
        lakes = solve(basins, runoff, evapo, precip)
        # A carved basin has drained; only survivors hold water. Summing lake
        # area over everything counts lakes in depressions the verdict has just
        # said should not exist.
        survives = ~carve
        area = np.where(survives, lakes["area_km2"], 0.0)
        results[label] = {
            "aridity_index": index,
            "carve": carve,
            "endorheic_land_km2": endorheic_km2,
            "lake_area_km2": area,
            "dry_survivors": int((survives & (runoff <= 0)).sum()),
            "wet_survivors": int((survives & (area > 1.0)).sum()),
            "converged": lakes["converged"],
        }

    crit = basins.catchment_km2 / np.maximum(basins.area_at_spill_km2, 1e-9) - 1.0
    both = results["wet"]["carve"] & results["potential"]["carve"]
    neither = ~results["wet"]["carve"] & ~results["potential"]["carve"]
    disputed = ~(both | neither)
    planet = 734_492_839.55

    print(f"{'bound':>10} {'carve':>7} {'survive':>8} {'dry':>6} {'with lake':>10} "
          f"{'lake % planet':>14}")
    for label in ("wet", "potential"):
        r = results[label]
        print(f"{label:>10} {int(r['carve'].sum()):7d} {int((~r['carve']).sum()):8d} "
              f"{r['dry_survivors']:6d} {r['wet_survivors']:10d} "
              f"{r['lake_area_km2'].sum()/planet*100:13.3f}%")
    print(f"\n  carve under both bounds : {int(both.sum()):5d}")
    print(f"  survive under both      : {int(neither.sum()):5d}")
    print(f"  disputed between bounds : {int(disputed.sum()):5d} "
          f"({disputed.sum()/n*100:.1f}%)")
    print(f"\n  median critical index (geometry): {np.median(crit):.2f}")

    payload = {
        "climatology": str(args.climatology),
        "coupling": str(args.coupling),
        "basins": n,
        "note": ("Two evaporation bounds. 'wet' uses land evaporation as "
                 "reported, a lower bound on open-water evaporation biased "
                 "toward carving. 'potential' divides it by the model's own "
                 "wetness factor, biased toward preserving. Basins disputed "
                 "between them are undecided by this climate."),
        "bounds": {
            label: {
                "basins_carved": int(r["carve"].sum()),
                "endorheic_land_km2": r["endorheic_land_km2"],
                "lake_area_km2": float(r["lake_area_km2"].sum()),
                "dry_survivors": r["dry_survivors"],
                "wet_survivors": r["wet_survivors"],
                "basins_without_runoff": int((runoff <= 0).sum()),
                "lake_area_fraction_of_planet": float(r["lake_area_km2"].sum() / planet),
                "solver_converged": bool(r["converged"]),
                "median_aridity_index": float(np.median(
                    r["aridity_index"][np.isfinite(r["aridity_index"])])),
            } for label, r in results.items()
        },
        "agreement": {
            "carve_under_both": int(both.sum()),
            "survive_under_both": int(neither.sum()),
            "disputed": int(disputed.sum()),
        },
        "carve_list_robust": [basins_id for basins_id, flag in
                              zip(_basin_ids(), both) if flag],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {args.output}")


def _basin_ids() -> list[str]:
    with Dataset(DATA / "basins.nc") as ds:
        return [str(x) for x in ds["basin_id"][:]]


if __name__ == "__main__":
    main()
