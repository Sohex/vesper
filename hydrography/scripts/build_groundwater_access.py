#!/usr/bin/env python3
"""Groundwater access on the climate grid: an AREA, never a mean depth. PLHY-5.

    python hydrography/scripts/build_groundwater_access.py
    python hydrography/scripts/build_groundwater_access.py --selftest

Worldbuilding. Vesper is an invented super-Earth; everything here is about a
steady-state water table solved on its native mesh and about what survives when
that mesh is reduced to the grid the climate model and the biosphere run on.

## Why this artifact exists

`build_groundwater.py` solves the water table PER REGION, on a 7.60 km mesh.
LPJ is driven per climate-grid cell. Between them somebody has to reduce ten
million regions to a few thousand cells, and the reduction that suggests itself
-- average the depth -- destroys exactly the thing groundwater is in the model
for. Root access and discharge vegetation depend on the FRACTION of a cell with
a shallow water table, not on its mean: averaging a narrow riparian strip into a
dry cell erases the strip, and averaging a broad shallow table with one deep
corner manufactures access the cell does not have. That is
`biosphere/notes/plant-hydraulics-groundwater-audit.md` finding 6.

So the criterion is evaluated ON THE NATIVE REGIONS, where the depth actually
lives, and what crosses to the grid is an area. **No hypsometry is needed and
none is used**: a hypsometric reconstruction of the within-cell depth
distribution would be an inference where the distribution itself is in hand.
GRID-2's hypsometry is for criteria that are terrain-conditioned and have no
native evaluation; this one is not.

A fraction is a fraction of a POPULATION and the only population this project
holds is the mesh regions inside a climate-grid cell, so every fraction here
lives at the grid and there is no per-region access fraction. That is the same
constraint `build_topographic_index.py` states for `f_sat_max` and it comes from
the same place, GW-6.

## What is written, and at which scale each thing lives

  per GRID CELL   accessible area and its share of the cell's land, one per
                  declared access depth; the conserved seepage total; the share
                  of land pinned at the surface; and the share of land in each
                  `sink_fraction` regime, including the two kinds of no-answer.
  per GRID CELL   the land-mean depth and its within-cell spread, carried with
                  a long_name that says what they are NOT for.

  data/<build>/groundwater_access_<grid>.nc
  analysis/groundwater_access_report.json

## What is NOT decided here

The access depths are DECLARED AND SWEPT, never chosen. One of them is sourced:
1.5 m is `SOILDEPTH_UPPER + SOILDEPTH_LOWER` in the vendored LPJ-GUESS's
`framework/guess.h`, so it is the depth below which the biosphere as it stands
has no roots at all whatever a real phreatophyte would do. The deeper arms are
model-form brackets for an access a phreatophytic extension would reach, and
they are declared rather than sourced. A consumer picks with its own argument
and the artifact carries all of them.

The water table itself is a bracket: permeability carries 1.5 to 2.5 orders of
Gleeson spread and `build_groundwater.py --sigma` is how that is explored. Pass
`--water-table` once per member and every fraction is written along a `member`
axis labelled by each file's own `sigma`, so the uncertainty regime crosses the
reduction instead of being averaged away inside it.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

from _paths import ANALYSIS, CONFIG, DATA          # noqa: F401  puts lib/ on the path
import builds
import gridding
from gridding import cell_fraction, cell_moments, cell_sum, transfer_ledger
from orogen import Export, LAND

CFG_PATH = Path(__file__).resolve().parents[1] / "config" / "groundwater.yaml"

# The access depths, in metres below the land surface. DECLARED HERE and not in
# the config because they are the axis of the artifact rather than a setting: an
# artifact whose axis moved between two builds would not be comparable with
# itself. 1.5 m is the vendored LPJ-GUESS soil column, SOILDEPTH_UPPER +
# SOILDEPTH_LOWER in `framework/guess.h`; 5 and 20 m are declared model-form
# brackets for phreatophytic access that model cannot currently represent.
ACCESS_DEPTHS_M = (1.5, 5.0, 20.0)

# `sink_fraction` regimes. The two cuts are the ones `build_groundwater.py`
# already prints and `water_table.nc` already documents: above 0.9 the
# evapotranspiration sink took the recharge locally and the depth is a recharge
# map in a water table's units; below 0.5 the lateral solve set it. The band
# between is neither, and it is reported as its own share rather than pushed to
# whichever side is convenient.
SINK_DOMINATED = 0.9
FLOW_DOMINATED = 0.5


def reduce_access(cell, ncell: int, area_km2, land, depth_m, seepage_m3_s,
                  at_surface, sink_fraction, excluded,
                  depths_m=ACCESS_DEPTHS_M) -> dict:
    """The whole reduction, on arrays, so it can be checked without a build.

    Every field is per region and the population is the cell's LAND, taken from
    `surface_class` by the caller -- `CLAUDE.md` rule 1, and this world's dry
    closed-basin floors below sea level are precisely where a shallow table
    matters most, so the mask that floods them would delete the answer.

    The operators are `lib/gridding.py`'s and the semantics are theirs:
    `cell_sum` for the seepage, which is a flux and must close globally;
    `cell_fraction` for every access and regime share, which are categorical
    retentions; `cell_moments` for the depth, which is the only intensive here
    and is not what the artifact is for.
    """
    land = np.asarray(land, dtype=bool)
    depth = np.asarray(depth_m, dtype=np.float64)
    out: dict[str, np.ndarray] = {}

    land_area = cell_sum(cell, ncell, np.asarray(area_km2, np.float64), land)
    covered = land_area > 0
    out["land_area_km2"] = land_area
    out["covered"] = covered

    # ACCESS. Evaluated per region and reduced as an area, which is the whole
    # point: `mean(depth) <= d` and `fraction(depth <= d)` are different
    # questions and only the second one is the habitat.
    for d in depths_m:
        frac, _ = cell_fraction(cell, ncell, area_km2, land & (depth <= d), land)
        out[f"access_fraction_{d:g}m"] = frac
        out[f"access_area_km2_{d:g}m"] = frac * land_area

    # SEEPAGE. Extensive: the cell total, and the global total is preserved.
    out["seepage_m3_s"] = cell_sum(cell, ncell,
                                   np.asarray(seepage_m3_s, np.float64), land)

    # THE REGIMES, each as a retained fraction of the cell's land.
    at_surface = np.asarray(at_surface, dtype=bool)
    excluded = np.asarray(excluded, dtype=bool)
    sink = np.asarray(sink_fraction, dtype=np.float64)
    undefined = ~np.isfinite(sink)
    out["at_surface_fraction"], _ = cell_fraction(cell, ncell, area_km2,
                                                  land & at_surface, land)
    out["excluded_fraction"], _ = cell_fraction(cell, ncell, area_km2,
                                                land & excluded, land)
    out["sink_undefined_fraction"], _ = cell_fraction(cell, ncell, area_km2,
                                                      land & undefined, land)
    defined = land & ~undefined
    for name, sel in (
            ("sink_dominated_fraction", defined & (sink >= SINK_DOMINATED)),
            ("sink_mixed_fraction",
             defined & (sink < SINK_DOMINATED) & (sink > FLOW_DOMINATED)),
            ("flow_dominated_fraction", defined & (sink <= FLOW_DOMINATED))):
        out[name], _ = cell_fraction(cell, ncell, area_km2, sel, land)

    # THE DEPTH, carried because SIMTOP's saturated-fraction closure needs a
    # cell-mean depth and for no other reason. It is NOT the access area.
    mean, var, count, _ = cell_moments(cell, ncell, area_km2, depth, land)
    out["depth_land_mean_m"] = mean
    out["depth_within_cell_stdev_m"] = np.sqrt(var)
    out["regions_per_cell"] = count.astype(np.int64)
    return out


def _read_water_table(path: Path, export: Export):
    with Dataset(path) as ds:
        if getattr(ds, "terrain_hash", None) not in (None, export.terrain_hash):
            raise SystemExit(
                f"{path} was solved on terrain {ds.terrain_hash[:16]} and the "
                f"mesh is {export.terrain_hash[:16]}; they are not the same "
                "world and no mapping between them is defined")
        n = ds.dimensions["region"].size
        if n != export.n_regions:
            raise SystemExit(
                f"{path} carries {n} regions and the mesh has "
                f"{export.n_regions}; not the same support")
        return {
            "depth_m": np.asarray(ds["depth_m"][:], dtype=np.float64),
            "seepage_m3_s": np.asarray(ds["seepage_m3_s"][:], dtype=np.float64),
            "at_surface": np.asarray(ds["at_surface"][:]).astype(bool),
            "sink_fraction": np.asarray(ds["sink_fraction"][:], dtype=np.float64),
            "excluded": np.asarray(ds["excluded"][:]).astype(bool),
            "sigma": float(getattr(ds, "sigma", 0.0)),
        }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--water-table", type=Path, action="append", default=None,
                    help="water_table.nc; repeat once per permeability member "
                         "and every fraction is written along a member axis")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="run the reduction against its invariants on "
                         "synthetic input; needs no build and no solve")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    export = Export(builds.mesh_export(config))
    resolution = str(config["model"]["resolution"]).upper()
    grid_dir = builds.grid_export(config, resolution)
    data = builds.component_data("hydrography", config, strict=True)

    members = args.water_table or [data / "water_table.nc"]
    for path in members:
        if not path.is_file():
            raise SystemExit(
                f"{path} does not exist. Run "
                "hydrography/scripts/build_groundwater.py first; this step "
                "reduces its output and cannot stand in for it.")

    cell, nlat, nlon = gridding.region_cells(export, grid_dir)
    ncell = nlat * nlon
    area_km2 = export.cell_area.astype(np.float64)
    land = export.surface_class == LAND

    # THE INDEX GUARD, and it is a real one. The coupling matrix for this rung
    # was built on the same column expression this reduction uses, so if the two
    # are not index-aligned almost none of a catchment's area would land on the
    # cells the model calls land -- `CLAUDE.md` rule 3, which has silently
    # matched zero cells on three scripts. `require_index_alignment` refuses.
    coupling = gridding.coupling_path(data, config)
    alignment = None
    if coupling.is_file():
        land_share, _ = cell_fraction(cell, ncell, area_km2, land)
        alignment = gridding.require_index_alignment(
            coupling, (land_share.reshape(nlat, nlon) >= 0.5))

    ledger = transfer_ledger(cell, ncell, area_km2, land)
    fields: dict[str, list[np.ndarray]] = {}
    sigmas = []
    for path in members:
        wt = _read_water_table(path, export)
        sigmas.append(wt["sigma"])
        reduced = reduce_access(cell, ncell, area_km2, land, wt["depth_m"],
                                wt["seepage_m3_s"], wt["at_surface"],
                                wt["sink_fraction"], wt["excluded"])
        for name, value in reduced.items():
            fields.setdefault(name, []).append(value)

    out = args.output or (data / f"groundwater_access_exoplasim-{resolution}.nc")
    if out.exists():
        raise SystemExit(
            f"{out} exists. This step writes a new artifact and never "
            "overwrites one; move it aside or pass --output.")
    covered = fields.pop("covered")[0]
    with Dataset(out, "w") as ds:
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        ds.createDimension("member", len(members))
        ds.title = "Groundwater access on the climate grid, as an area"
        ds.terrain_hash = export.terrain_hash
        ds.setncattr("vesper_source_build", builds.build_root(config).name)
        ds.resolution = resolution
        ds.access_depths_m = list(ACCESS_DEPTHS_M)
        ds.sigma = sigmas
        ds.sources = [str(p) for p in members]
        ds.caveat = (
            "ACCESS IS AN AREA. depth_land_mean_m is carried for SIMTOP's "
            "saturated-fraction closure and is NOT an oasis or riparian area "
            "proxy: the fraction of a cell with a shallow table is what root "
            "access depends on, and a mean depth erases a narrow discharge "
            "zone and manufactures a broad one. Every fraction is a fraction "
            "of the cell's own LAND population, taken from surface_class; "
            "there is no per-region access fraction and GW-6 is why."
        )
        v = ds.createVariable("covered", "i1", ("lat", "lon"), zlib=True)
        v.long_name = ("this cell holds mesh land; every fraction below is "
                       "undefined where it is zero and must be masked, not read")
        v[:] = covered.reshape(nlat, nlon).astype(np.int8)
        for name, members_of in fields.items():
            stacked = np.stack([m.reshape(nlat, nlon) for m in members_of])
            dtype = "i4" if name == "regions_per_cell" else "f4"
            v = ds.createVariable(name, dtype, ("member", "lat", "lon"),
                                  zlib=True, complevel=4)
            v[:] = stacked
            if name.startswith("access_fraction"):
                v.units = "1"
                v.long_name = ("share of the cell's land whose water table is "
                               "within this depth, evaluated per mesh region "
                               "and reduced as an area")
            elif name.startswith("access_area"):
                v.units = "km2"
            elif name == "seepage_m3_s":
                v.units = "m3 s-1"
                v.long_name = ("groundwater returning to the surface, summed "
                               "over the cell's land; conserved exactly")
            elif name.endswith("_fraction"):
                v.units = "1"
            elif name.startswith("depth"):
                v.units = "m"
                v.long_name = ("NOT an access-area proxy; see the caveat "
                               "attribute")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "issue": "PLHY-5",
        "build": builds.build_root(config).name,
        "terrain_hash": export.terrain_hash,
        "resolution": resolution,
        "access_depths_m": list(ACCESS_DEPTHS_M),
        "members": [{"path": str(p), "sigma": s} for p, s in zip(members, sigmas)],
        "transfer_ledger": ledger,
        "coupling_ocean_fraction": alignment,
        "land_totals": {
            f"access_area_km2_{d:g}m": float(fields[f"access_area_km2_{d:g}m"][0].sum())
            for d in ACCESS_DEPTHS_M},
        "seepage_total_m3_s": float(fields["seepage_m3_s"][0].sum()),
        "file": str(out),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    rp = ANALYSIS / "groundwater_access_report.json"
    rp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}\nwrote {rp}")
    return 0


# --- the checks, on synthetic input, with right answers ---------------------

def _selftest() -> int:
    """Six invariants, each with an answer that can be wrong.

    The important one is the third: a cell whose MEAN depth is far below any
    access threshold but which holds a narrow shallow strip must return that
    strip's area share and not zero. That is the failure this artifact exists
    to prevent, and a reduction that averaged the depth would fail it while
    every conservation check above still passed.
    """
    problems: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    rng = np.random.default_rng(20260824)
    n, ncell = 6000, 40
    cell = rng.integers(0, ncell, n)
    area = rng.uniform(0.5, 5.0, n)
    land = rng.random(n) < 0.7
    depth = rng.uniform(0.0, 60.0, n)
    seepage = rng.uniform(0.0, 2.0, n)
    at_surface = depth < 0.01
    sink = rng.uniform(0.0, 1.0, n)
    excluded = rng.random(n) < 0.05

    r = reduce_access(cell, ncell, area, land, depth, seepage, at_surface,
                      sink, excluded)
    covered = r["covered"]

    # 1. The flux closes: nothing is lost and nothing is created.
    got, want = float(r["seepage_m3_s"].sum()), float(seepage[land].sum())
    check("seepage is conserved over the reduction",
          abs(got - want) <= 1e-9 * abs(want), f"{got} against {want}")

    # 2. The land area closes against the population it was taken from.
    got, want = float(r["land_area_km2"].sum()), float(area[land].sum())
    check("land area closes", abs(got - want) <= 1e-9 * abs(want),
          f"{got} against {want}")

    # 3. THE NEGATIVE CONTROL. One cell, all of it deep except a narrow strip.
    cell2 = np.zeros(1000, dtype=np.int64)
    area2 = np.ones(1000)
    land2 = np.ones(1000, dtype=bool)
    depth2 = np.full(1000, 50.0)
    depth2[:30] = 0.4                      # 3% of the cell, within every depth
    r2 = reduce_access(cell2, 1, area2, land2, depth2, np.zeros(1000),
                       depth2 < 0.01, np.full(1000, 0.5),
                       np.zeros(1000, dtype=bool))
    frac = float(r2["access_fraction_1.5m"][0])
    mean_depth = float(r2["depth_land_mean_m"][0])
    check("a narrow shallow strip survives a deep cell's mean",
          abs(frac - 0.03) <= 1e-12 and mean_depth > 40.0,
          f"fraction {frac}, mean depth {mean_depth}")

    # 4. Access is monotone in the threshold: a deeper root reaches at least as
    #    much. A reduction that averaged first could violate this by clipping.
    a15 = r["access_fraction_1.5m"][covered]
    a5 = r["access_fraction_5m"][covered]
    a20 = r["access_fraction_20m"][covered]
    check("access fraction is monotone in the access depth",
          bool((a5 >= a15 - 1e-12).all() and (a20 >= a5 - 1e-12).all()),
          "a deeper threshold reached less land somewhere")

    # 5. The sink regimes partition the defined land exactly.
    total = (r["sink_dominated_fraction"] + r["sink_mixed_fraction"]
             + r["flow_dominated_fraction"] + r["sink_undefined_fraction"])
    check("the sink regimes partition the cell's land",
          bool(np.allclose(total[covered], 1.0, rtol=0, atol=1e-12)),
          f"max deviation {np.abs(total[covered] - 1.0).max():.3g}")

    # 6. An all-accessible population returns exactly one, and an inaccessible
    #    one exactly zero. Constant preservation at both ends of the criterion.
    r3 = reduce_access(cell, ncell, area, land, np.zeros(n), seepage,
                       at_surface, sink, excluded)
    r4 = reduce_access(cell, ncell, area, land, np.full(n, 1e6), seepage,
                       at_surface, sink, excluded)
    check("a wholly shallow cell reads one and a wholly deep one reads zero",
          bool(np.allclose(r3["access_fraction_1.5m"][covered], 1.0, atol=0)
               and np.all(r4["access_fraction_20m"][covered] == 0.0)),
          "the criterion did not saturate")

    print(f"\n6 checks, {len(problems)} failed")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
