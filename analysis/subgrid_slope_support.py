#!/usr/bin/env python3
"""What does `build_soil.py:subgrid_slope` carry, the terrain or the mesh? WORLD-UYFB.

    python analysis/subgrid_slope_support.py

Worldbuilding. Vesper is an invented super-Earth; every quantity here is a
modelled field of that planet, measured on two World Orogen exports of the same
planet at the same seed and different region counts. Nothing is run: no climate
model, no soil solve, no biosphere.

## The question

`pedology/scripts/build_soil.py:subgrid_slope` divides the within-cell
elevation spread by the MESH spacing and hands the quotient to the catena term
as a gradient. `notes/audits/nonlinear-spatial-reductions.md` section 6 reports
that the quotient doubles between this project's two builds while the spread
moves 1.6 per cent. A statistic that shifts with the region count has one of
three mechanisms, and the repair differs for each:

  A  A LENGTH IS IN THE DEFINITION. The mesh spacing appears explicitly, so the
     shift is the spacing ratio and nothing else, it is the same at every
     quantile, and it is predictable before the second build is measured.
     Repaired by declaring the length instead of reading it from the mesh.
     `hydrography/notes/subgrid-water-table.md` finds this for the compound
     topographic index, whose `a` carries a length: 0.86, about ln 2.

  B  THE STATISTIC IS BIASED BY THE BALL POPULATION. An extreme-value reduction
     -- a minimum or a maximum over a neighbourhood -- walks with how many
     regions the neighbourhood holds even where the terrain is identical. Not
     repairable by any normalisation; the estimator has to change.
     `notes/audits/orogen-resolution.md` finds this for the minimum-over-30-km
     relief form and rejects it on those grounds.

  C  THE TERRAIN IS SELF-AFFINE AND THE RUN IS THE MESH'S. A ratio of two
     lengths carries no length, yet a gradient sampled over one mesh edge still
     steepens on a finer mesh because the terrain has spectral content at that
     scale. The shift is empirical, is NOT one number, and drifts with the
     quantile. `notes/audits/orogen-resolution.md` measures 1.408 at p50 rising
     to 1.695 at p99 for `computeScarpPotential`'s one-edge gradient. Repaired
     only by fixing the run at a declared physical length.

## The three tests, and the bar, declared before anything was run

The bar is the one `notes/audits/orogen-resolution.md` fixed for the same class
and is not re-chosen here: a quantity TRANSPORTS when every land quantile from
p50 to p99 agrees between the 2.5M and the 10M build within 1.15x.

  test A   Is the shift in `subgrid_slope` the mesh-spacing ratio? Mechanism A
           predicts, before the run, that
               ratio(subgrid_slope) = ratio(spacing) x ratio(spread)
           holds at EVERY quantile, with no residual. A residual above 1.15x
           anywhere means something other than the divisor is moving.

  test B   Is the numerator biased by the per-cell region population? The fine
           build holds about four times as many regions in a T42 cell as the
           coarse one. Drawing a seeded random one-in-four subsample of the fine
           build's land regions and recomputing the spread holds the terrain
           fixed and changes only the population. Mechanism B predicts a shift
           above the bar; mechanisms A and C predict none.

  test C   Does the numerator itself carry the support? The spread is a
           within-cell statistic of a fixed cell, so the two builds sample the
           SAME ground at two spacings. Mechanism C predicts a spread ratio that
           departs from one and drifts with the quantile the way the one-edge
           gradient does. Mechanism A predicts a flat ratio inside the bar.

## What is reported, and in whose units

The consuming step is `soil`, and its instrument is `regolith.minimum_depth_m`,
the thinnest profile the pedogenesis model distinguishes. So the consequence of
the divisor is reported as the change in the catena divisor
`1 + slope_transport * tan(beta)` and in metres of regolith at the land mean,
against that bar, and not as a percentage of the gradient.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))

from gridding import cell_moments, region_cells                # noqa: E402
from orogen import Export, LAND                                # noqa: E402

OUTPUT = PROJECT_ROOT / "analysis" / "subgrid_slope_support.json"

BUILDS = ("precarve-craton", "precarve-craton-10m")
RUNG = "T42"

# Declared before the run: a cell needs this many land regions before its
# elevation spread is a statistic rather than a pair of samples. Thirty is one
# tenth of what a T42 cell holds on the COARSE build, so the cut is set by the
# thinner support and is the same cut on both.
MIN_LAND_REGIONS = 30

QUANTILES = (0.50, 0.75, 0.90, 0.95, 0.99)
TRANSPORT_BAR = 1.15
SUBSAMPLE_KEEP = 0.25
SUBSAMPLE_SEED = 16236323


def spread_by_cell(export: Export, grid_dir: Path, population: np.ndarray | None = None):
    """Area-weighted within-cell standard deviation of land elevation, in metres.

    `lib/gridding.py:cell_moments` is the operator; the two-pass variance it
    takes is the one that does not cancel on a field whose spread is tens of
    metres about a mean of thousands.
    """
    cell, nlat, nlon = region_cells(export, grid_dir)
    ncell = nlat * nlon
    area = export.cell_area.astype(np.float64)
    elevation_m = export.elevation_km.astype(np.float64) * 1000.0
    land = export.surface_class == LAND
    if population is not None:
        land = land & population
    _, var, count, _ = cell_moments(cell, ncell, area, elevation_m, population=land)
    return np.sqrt(var), count


def quantiles(values: np.ndarray) -> dict[str, float]:
    return {f"p{int(q * 100)}": float(np.quantile(values, q)) for q in QUANTILES}


def ratio(fine: dict[str, float], coarse: dict[str, float]) -> dict[str, float]:
    return {k: float(fine[k] / coarse[k]) for k in fine}


def main() -> int:
    pedo = yaml.safe_load(
        (PROJECT_ROOT / "pedology" / "config" / "pedogenesis.yaml").read_text())
    slope_transport = float(pedo["catena"]["slope_transport"])
    baseline_km = float(pedo["catena"]["gradient_baseline_km"])

    per_build = {}
    keep = {}
    for name in BUILDS:
        root = PROJECT_ROOT / "source" / name / f"exoplasim-{RUNG}"
        export = Export(root)
        spacing_km = float(export.manifest["basins"]["resolution"]["avgEdgeKm"])
        spread, count = spread_by_cell(export, root)
        usable = count >= MIN_LAND_REGIONS
        keep[name] = usable
        per_build[name] = {
            "regions": export.n_regions,
            "terrain_hash": export.terrain_hash[:8],
            "mesh_spacing_km": spacing_km,
            "cells_used": int(usable.sum()),
            "land_regions_per_cell_median": float(np.median(count[usable])),
            "spread_m": spread,
        }

    # Marginal quantiles, on each build's own usable cells, and the paired
    # comparison on the cells both builds call usable. The marginal form is what
    # the 1.15x bar was fixed on; the paired form is the tighter instrument and
    # cannot be confused by a change in which cells qualify.
    common = keep[BUILDS[0]] & keep[BUILDS[1]]
    coarse, fine = (per_build[b] for b in BUILDS)
    for build in (coarse, fine):
        u = build["spread_m"][common]
        build["spread_quantiles_m"] = quantiles(u)
        build["tan_beta_mesh"] = quantiles(u / (build["mesh_spacing_km"] * 1000.0))
        build["tan_beta_baseline"] = quantiles(u / (baseline_km * 1000.0))
        build["land_mean_tan_beta_mesh"] = float(
            np.mean(u / (build["mesh_spacing_km"] * 1000.0)))
        build["land_mean_tan_beta_baseline"] = float(np.mean(u / (baseline_km * 1000.0)))

    spacing_ratio = coarse["mesh_spacing_km"] / fine["mesh_spacing_km"]
    spread_ratio = ratio(fine["spread_quantiles_m"], coarse["spread_quantiles_m"])
    mesh_form_ratio = ratio(fine["tan_beta_mesh"], coarse["tan_beta_mesh"])
    baseline_form_ratio = ratio(fine["tan_beta_baseline"], coarse["tan_beta_baseline"])
    # test A: the residual after the divisor and the numerator are accounted for.
    residual = {k: mesh_form_ratio[k] / (spacing_ratio * spread_ratio[k])
                for k in spread_ratio}

    # test B: population bias, terrain held fixed.
    rng = np.random.default_rng(SUBSAMPLE_SEED)
    fine_root = PROJECT_ROOT / "source" / BUILDS[1] / f"exoplasim-{RUNG}"
    fine_export = Export(fine_root)
    draw = rng.random(fine_export.n_regions) < SUBSAMPLE_KEEP
    sub_spread, sub_count = spread_by_cell(fine_export, fine_root, population=draw)
    sub_usable = common & (sub_count >= MIN_LAND_REGIONS)
    sub_q = quantiles(sub_spread[sub_usable])
    full_q = quantiles(fine["spread_m"][sub_usable])
    population_ratio = {k: float(sub_q[k] / full_q[k]) for k in sub_q}

    # The consequence, in the consuming step's own units.
    depth_scale = {}
    for label, form in (("mesh", "tan_beta_mesh"), ("baseline", "tan_beta_baseline")):
        divisor = {b["regions"]: 1.0 + slope_transport * b[f"land_mean_{form}"]
                   for b in (coarse, fine)}
        depth_scale[label] = divisor

    verdict = {
        "test_A_residual_within_bar": all(
            abs(np.log(v)) <= np.log(TRANSPORT_BAR) for v in residual.values()),
        "test_B_population_within_bar": all(
            abs(np.log(v)) <= np.log(TRANSPORT_BAR) for v in population_ratio.values()),
        "test_C_numerator_within_bar": all(
            abs(np.log(v)) <= np.log(TRANSPORT_BAR) for v in spread_ratio.values()),
        "mesh_form_transports": all(
            abs(np.log(v)) <= np.log(TRANSPORT_BAR) for v in mesh_form_ratio.values()),
        "baseline_form_transports": all(
            abs(np.log(v)) <= np.log(TRANSPORT_BAR) for v in baseline_form_ratio.values()),
    }

    out = {
        "measured": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "rung": RUNG,
        "bar": TRANSPORT_BAR,
        "min_land_regions_per_cell": MIN_LAND_REGIONS,
        "cells_compared": int(common.sum()),
        "mesh_spacing_ratio": spacing_ratio,
        "gradient_baseline_km": baseline_km,
        "slope_transport": slope_transport,
        "builds": {
            name: {k: v for k, v in per_build[name].items() if k != "spread_m"}
            for name in BUILDS},
        "spread_ratio": spread_ratio,
        "mesh_form_ratio": mesh_form_ratio,
        "baseline_form_ratio": baseline_form_ratio,
        "test_A_residual": residual,
        "test_B_population_ratio": population_ratio,
        "catena_divisor": depth_scale,
        "verdict": verdict,
    }
    OUTPUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    def row(label, d):
        return f"  {label:<34}" + "".join(f"{d[f'p{int(q*100)}']:>9.3f}" for q in QUANTILES)

    print(f"cells compared: {out['cells_compared']}   spacing ratio "
          f"{spacing_ratio:.4f}   bar {TRANSPORT_BAR}x")
    print("  " + " " * 34 + "".join(f"{'p'+str(int(q*100)):>9}" for q in QUANTILES))
    print(row("spread, fine/coarse (test C)", spread_ratio))
    print(row("spread, 1-in-4 population (test B)", population_ratio))
    print(row("spread / mesh spacing, fine/coarse", mesh_form_ratio))
    print(row("residual after divisor (test A)", residual))
    print(row("spread / declared baseline", baseline_form_ratio))
    for label, d in depth_scale.items():
        print(f"  catena divisor, {label:<10}" +
              "".join(f"  {n}: {v:.4f}" for n, v in d.items()))
    for k, v in verdict.items():
        print(f"  {k}: {v}")
    print(f"wrote {OUTPUT.relative_to(PROJECT_ROOT)}")
    return 0



if __name__ == "__main__":
    raise SystemExit(main())
