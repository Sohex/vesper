#!/usr/bin/env python3
"""The compound topographic index, and the saturated-area closure on it. GW-26.

Worldbuilding. Vesper is an invented planet and everything here is about a
terrain generator's output and a statistical closure applied to it.

WHY THIS EXISTS. GW-6 established that the valley-to-ridge water table gradient
is below this mesh and stays below it, on five separate measurements: the
solver's depth field is very nearly a function of its recharge forcing where the
real one is not, and neither a finer generation nor a finer cell recovers the
difference. TOPMODEL does not resolve that gradient either. It carries the
terrain control STATISTICALLY, and what it returns is a saturated FRACTION
rather than a depth -- which is the quantity WET-2's wetness classification,
SURF-7's discharge mask and LSHY-6's mosaic actually consume.

WHAT IS COMPUTED, and at which scale each thing lives:

  per REGION      the index itself, `ln(a / tan beta)`, from the export's own
                  upslope contributing area and a local slope. Terrain only:
                  no climate, no water table, no solve.
  per GRID CELL   `f_sat_max`, the share of the cell's land AREA whose index
                  exceeds the cell's area-weighted mean index. A saturated
                  fraction is a fraction of a cell's AREA, so the population it
                  is a rank statistic over is that area and not the region
                  count: this mesh's regions are not equal-area, and a count
                  weights a sliver and a full cell alike in both the share and
                  the mean the share is taken about. `lib/gridding.py` owns the
                  operators; `cell_moments` and `cell_fraction` are the two used
                  and nothing is reimplemented here.
                  A region has no sub-population and GW-6 refuses to invent one,
                  so there is no per-region saturated fraction here and there
                  cannot be. world-d9u4 decides what the three rows waiting on
                  that get instead; `hydrography/notes/subgrid-water-table.md`
                  section 5 carries it.

  data/<build>/topographic_index_<grid>.nc
  analysis/topographic_index_report.json

WHAT IS NOT COMPUTED. `f_sat` itself, because it needs a cell-mean water table
depth and because no consumer may take it before the score declared in
`config/topographic_index.yaml` has been run. `saturated_fraction` below is the
closure, for the scorer and for whichever consumer is licensed first; it is a
function and not an artifact on purpose.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, DATA          # noqa: F401  puts lib/ on the path
import builds
import gridding
from orogen import Export, LAND

CFG_PATH = Path(__file__).resolve().parents[1] / "config" / "topographic_index.yaml"
SCORE_PATH = ANALYSIS / "topographic_index_score.json"


def license_state() -> dict:
    """Whether a consumer may take `f_sat`, read from the score and never set here.

    The score is run by the Earth harness the config names, on a real planet
    with real bores, and it is the only thing that can license this closure. So
    the state is READ rather than declared: a build of this artifact cannot
    license itself, and one run before any score has been taken says so with the
    same words a failed score does.
    """
    if not SCORE_PATH.exists():
        return {"scored": False, "consumers_licensed": False,
                "reason": f"no score at {SCORE_PATH.name}; the harness in "
                          "config/topographic_index.yaml has not been run"}
    sc = json.loads(SCORE_PATH.read_text())
    runs = sc.get("runs", {})
    return {
        "scored": bool(runs),
        "consumers_licensed": bool(sc.get("consumers_licensed", False)),
        "observation_sets_required": sc.get("observation_sets_required"),
        "scored_sets": sc.get("scored_sets"),
        "per_set": {k: {"passes": v.get("passes"), "verdict": v.get("verdict")}
                    for k, v in runs.items()},
    }


def saturated_fraction(f_sat_max, water_table_depth_m, f_grad_per_m: float):
    """`f_sat = min(f_sat_max * exp(-f_grad * z_wt), 1)`. SIMTOP, Niu et al. (2005).

    THE CAP IS EXPLICIT because the product is not otherwise bounded: ClimaLand
    writes `min(..., 1)` and so does this. `water_table_depth_m` is positive
    downward, so a table at the surface returns `f_sat_max` and a deep one
    returns nothing.

    THE CONVENTION IS A TRAP, and it is why `f_grad_per_m` has no default here.
    CLIMBER-X writes `exp(-f_wtab * w_table)` and ships `f_wtab = 2.5` per metre;
    ClimaLand writes `exp(-f_over/2 * z_wt)`, so the same numeral means half as
    much per metre there. This function takes the decay in the FIRST convention,
    per metre of depth, and the config brackets exactly that factor of two.
    """
    f = np.asarray(f_sat_max, dtype=np.float64)
    z = np.asarray(water_table_depth_m, dtype=np.float64)
    return np.minimum(f * np.exp(-float(f_grad_per_m) * z), 1.0)


def receiver_slope(export: Export):
    """tan(beta) as the drop to a region's own drainage receiver over the
    distance to it, and the distance.

    The steepest-descent estimate a grid TOPMODEL uses, computed from
    `drain_to`, which is the export's own descent on the depression-filled
    surface. It is ZERO inside a filled pit by construction -- a filled surface
    is flat there -- which is one of the two reasons the index needs a floor.

    The unit-sphere coordinates are `x, y, z` in the generator's own y-up frame,
    which is safe here because a great-circle distance is rotation invariant;
    `source/README.md` says why mixing them with lat/lon is not.
    """
    n = export.n_regions
    pos = np.stack([export.field("x"), export.field("y"), export.field("z")],
                   axis=1).astype(np.float64)
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)
    recv = export.field("drain_to").astype(np.int64)
    ok = (recv >= 0) & (recv < n) & (recv != np.arange(n))
    to = np.where(ok, recv, 0)
    dot = np.clip(np.einsum("ij,ij->i", pos, pos[to]), -1.0, 1.0)
    d_m = export.radius_km * 1000.0 * np.arccos(dot)
    elev_m = export.field("elevation_km").astype(np.float64) * 1000.0
    drop = elev_m - elev_m[to]
    return np.where(ok & (d_m > 0.0), drop / np.maximum(d_m, 1.0), 0.0), d_m


def compound_index(export: Export, cfg: dict):
    """`CTI = ln(a / tan beta)` per region, both slope arms, and the diagnostics.

    `a` is the export's `flow_accumulation` -- upstream contributing area on the
    filled surface, in real area rather than cell counts -- over a contour width
    taken as the square root of the cell's own area. There is no exact contour
    width on a Voronoi mesh; this is the characteristic width of the cell the
    water crosses, and it is declared in the config rather than derived.
    """
    icfg = cfg["index"]
    if str(icfg["contour_width"]) != "sqrt_cell_area":
        raise SystemExit(
            f"index.contour_width {icfg['contour_width']!r} is not one this "
            "script has; it is 'sqrt_cell_area'")
    floor = float(icfg["min_tan_slope"])

    area_m2 = export.cell_area.astype(np.float64) * 1e6
    accum_m2 = export.field("flow_accumulation").astype(np.float64) * 1e6
    a_m = accum_m2 / np.sqrt(area_m2)

    tan_plane = np.tan(np.deg2rad(export.local_slope_deg.astype(np.float64)))
    tan_recv, _ = receiver_slope(export)

    arms = {}
    for name, tanb in (("plane_fit", tan_plane), ("receiver_drop", tan_recv)):
        arms[name] = np.log(a_m / np.maximum(tanb, floor))
    primary = str(icfg["slope"])
    if primary not in arms:
        raise SystemExit(f"index.slope {primary!r} is not one this script has")

    land = export.surface_class == LAND
    diag = {
        "min_tan_slope": floor,
        "at_slope_floor_land_fraction": {
            name: float((np.maximum(t, 0.0)[land] <= floor).mean())
            for name, t in (("plane_fit", tan_plane), ("receiver_drop", tan_recv))},
        "land_percentiles": {
            name: {str(p): float(np.percentile(v[land], p))
                   for p in (1, 5, 25, 50, 75, 95, 99)}
            for name, v in arms.items()},
        "land_mean": {name: float(v[land].mean()) for name, v in arms.items()},
    }
    return arms, primary, land, diag


def per_cell(cti, area_km2, land, cell, ncell, min_regions: int):
    """Cell-mean index and `f_sat_max`, the share of the cell's AREA above it.

    CLIMBER-X computes `1 - cti_cdf(cti_mean)` from a CDF tabulated on integer
    index bins; the same quantity is available directly from the population,
    without a lookup table and without committing to a bin range this world's
    index does not fit inside.

    IT IS A SHARE OF AREA AND NOT OF REGIONS, which is the whole of world-d9u4's
    answer expressed in code. `f_sat_max` multiplies into a saturated FRACTION
    OF A CELL, so the population it is a rank statistic over is the cell's land
    area; CLIMBER-X's CDF is over equal-area DEM pixels, where the two
    coincide, and the regions of this mesh are not equal-area -- the land
    population spans a factor of 5.2 from its 5th to its 95th percentile. A
    region count therefore weights a sliver and a full cell alike, in both the
    share AND the cell mean that share is taken about, and the mean is the class
    boundary. `lib/gridding.py` owns the operators and names which is right for
    which field semantics: `cell_moments` for the INTENSIVE mean and the spread,
    `cell_fraction` for the CATEGORICAL share of a population in one class.
    Neither is reimplemented here.

    THE SPREAD IS THE TWO-PASS VARIANCE for the reason `cell_moments` exists:
    the difference-of-moments form this used cancels, and does not return zero
    on a constant field.

    A cell with too few regions is left unset rather than given a fraction
    computed from a handful of them: a share above a mean is meaningless on
    three samples, and returning one would be a number where there is no
    estimate. That guard stays a REGION COUNT, because it is about sample size
    and not about area.
    """
    mean, var, count, covered = gridding.cell_moments(
        cell, ncell, area_km2, cti, land)
    f_sat_max, _ = gridding.cell_fraction(
        cell, ncell, area_km2, cti > mean[cell], land)
    have = covered & (count >= int(min_regions))
    return (count.astype(np.int64), np.where(have, mean, np.nan),
            np.sqrt(var), np.where(have, f_sat_max, np.nan), have)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="GW-26: the compound topographic index on the native mesh, "
                    "and the saturated-fraction closure's terrain half.")
    ap.add_argument("--resolution", default=None,
                    help="grid to aggregate onto; defaults to config's model "
                         "resolution")
    ap.add_argument("--build", default=None,
                    help="build directory under source/; defaults to the active one")
    ap.add_argument("--min-regions", type=int, default=10,
                    help="fewest regions a grid cell may hold and still be given "
                         "an f_sat_max")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    cfg = yaml.safe_load(CFG_PATH.read_text(encoding="utf-8"))
    if str(cfg["closure"]["f_sat_max"]) != "area_share_above_cell_mean":
        raise SystemExit(
            f"closure.f_sat_max {cfg['closure']['f_sat_max']!r} is not the "
            "quantity this script computes. It computes the share of a cell's "
            "land AREA above that cell's area-weighted mean index; a region "
            "COUNT share is a different quantity on this unequal-area mesh and "
            "is not available here.")
    if str(cfg["closure"]["absolute_thresholds"]) != "refused":
        raise SystemExit(
            "closure.absolute_thresholds is not 'refused'. This index is not on "
            "the scale any published absolute threshold is calibrated against, "
            "and the config carries the measurement; there is no product here "
            "that keys on one.")

    export = Export(builds.mesh_export(config) if args.build is None
                    else builds.mesh_export_of(builds.SOURCE / args.build))
    grid_dir = builds.grid_export(config, args.resolution, build=args.build)
    grid_name = grid_dir.name
    print(f"mesh {export.n_regions:,} regions from {export.root.name}, "
          f"aggregating onto {grid_name}")

    arms, primary, land, diag = compound_index(export, cfg)
    cti = arms[primary]
    print(f"  index on {int(land.sum()):,} land regions, slope arm {primary!r}: "
          f"land mean {diag['land_mean'][primary]:.2f}, 5th to 95th "
          f"{diag['land_percentiles'][primary]['5']:.2f} to "
          f"{diag['land_percentiles'][primary]['95']:.2f}")
    print(f"  the other arm {('receiver_drop' if primary == 'plane_fit' else 'plane_fit')!r} "
          f"gives a land mean of "
          f"{diag['land_mean']['receiver_drop' if primary == 'plane_fit' else 'plane_fit']:.2f}; "
          f"the two are the bracket and both are shipped")

    cell, nlat, nlon = gridding.region_cells(export, grid_dir)
    ncell = nlat * nlon
    area_km2 = export.cell_area.astype(np.float64)
    count, cmean, csd, f_sat_max, have = per_cell(
        cti, area_km2, land, cell, ncell, args.min_regions)
    ledger = gridding.transfer_ledger(cell, ncell, area_km2, land)
    print(f"  {int(have.sum()):,} of {ncell:,} grid cells hold at least "
          f"{args.min_regions} land regions")
    print(f"  cell-mean index 5th/50th/95th: "
          f"{np.nanpercentile(cmean[have], 5):.2f} / "
          f"{np.nanpercentile(cmean[have], 50):.2f} / "
          f"{np.nanpercentile(cmean[have], 95):.2f}")
    print(f"  f_sat_max 5th/50th/95th: "
          f"{np.nanpercentile(f_sat_max[have], 5):.3f} / "
          f"{np.nanpercentile(f_sat_max[have], 50):.3f} / "
          f"{np.nanpercentile(f_sat_max[have], 95):.3f}")

    # THE SCALE STATEMENT, printed rather than buried, because it is what
    # decides which half of this artifact transports. Every published absolute
    # threshold on this index is calibrated against one computed at about a
    # kilometre on Earth; CLIMBER-X tabulates its CDF on integer bins 1 to 15
    # and zeroes the maximum wetland fraction above a cell mean of 14.
    over = float((cmean[have] > 14.0).mean())
    print(f"  and {over:.1%} of those cells have a mean index above 14, which "
          f"is the top of the range\n  CLIMBER-X's tabulated CDF covers. The "
          f"RANK statistic f_sat_max survives that offset;\n  an absolute "
          f"threshold does not, and none is written.")

    lic = license_state()
    print(f"\n  consumers_licensed: {lic['consumers_licensed']}"
          + ("" if lic["scored"] else f"   ({lic['reason']})"))

    data_dir = (builds.component_data("hydrography", config)
                if args.build is None
                else DATA / args.build)
    data_dir.mkdir(parents=True, exist_ok=True)
    out = args.output or (data_dir / f"topographic_index_{grid_name}.nc")

    with Dataset(out, "w", format="NETCDF4") as ds:
        ds.title = "Compound topographic index and saturated-fraction closure"
        ds.summary = (
            "GW-26. Per region the index ln(a / tan beta); per grid cell the "
            "share of that cell's regions whose index exceeds the cell mean, "
            "which is the saturated fraction's terrain half. A saturated "
            "fraction is a fraction of a population and the only population "
            "here is the land AREA of the regions inside a grid cell, so there "
            "is no per-region fraction and GW-6 says why there cannot be. "
            "world-d9u4 records what the three rows that wanted a native-mesh "
            "one get instead.")
        ds.consumers_licensed = "yes" if lic["consumers_licensed"] else "no"
        ds.uncertified = (
            "NO CONSUMER MAY TAKE f_sat FROM THIS UNLESS consumers_licensed "
            "SAYS YES. The license is the score declared in "
            "hydrography/config/topographic_index.yaml, run by the Earth "
            "harness and recorded in hydrography/analysis/"
            "topographic_index_score.json; this file reports it and cannot "
            "grant it. The depth field this closure multiplies failed its own "
            "bar; a fraction derived from it is not licensed by the "
            "derivation.")
        ds.f_grad_bracket_per_m = np.asarray(
            cfg["closure"]["f_grad_bracket_per_m"], dtype=np.float64)
        ds.f_grad_convention = (
            "f_sat = min(f_sat_max * exp(-f_grad * z_wt), 1), z_wt positive "
            "downward, f_grad per metre in CLIMBER-X's convention. ClimaLand "
            "writes exp(-f_over/2 * z_wt), where the same numeral is half this "
            "one; the bracket IS that factor of two.")
        ds.absolute_thresholds = (
            "REFUSED. This index is not on the scale published absolute "
            "thresholds are calibrated against, and no rescaling supplies one.")
        ds.source_build = export.root.parent.name
        ds.grid = grid_name
        ds.slope_arm = primary
        ds.created = datetime.now(timezone.utc).isoformat()

        ds.createDimension("region", export.n_regions)
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)

        def var(name, values, dtype, dims, units, desc):
            v = ds.createVariable(name, dtype, dims, zlib=True, complevel=4)
            v.units = units
            v.description = desc
            v[:] = values
            return v

        var("cti", cti, "f4", ("region",), "1",
            f"ln(a / tan beta) with the {primary} slope; a is the export's "
            "flow_accumulation over the square root of the cell area")
        other = "receiver_drop" if primary == "plane_fit" else "plane_fit"
        var("cti_alt", arms[other], "f4", ("region",), "1",
            f"the same index with the {other} slope. The bracket, not a second "
            "answer: a slope on this mesh has two defensible estimators and "
            "neither is a hillslope gradient")
        var("cell_index", cell, "i4", ("region",), "1",
            "flat grid-cell index this region falls in, row * nlon + col, from "
            "lib/gridding.py. CLAUDE.md rule 3: index, never longitude")

        shape = (nlat, nlon)
        var("cti_mean", cmean.reshape(shape), "f4", ("lat", "lon"), "1",
            "AREA-WEIGHTED mean index over the cell's land, lib/gridding.py's "
            "INTENSIVE reduction; NaN where the cell holds too few regions to "
            "estimate one. It is the boundary f_sat_max is a share above, so "
            "the weighting is part of the definition and not a refinement")
        var("cti_sd", csd.reshape(shape), "f4", ("lat", "lon"), "1",
            "spread of the index WITHIN the cell, which is the sub-grid "
            "population this closure exists to use. Area-weighted and taken in "
            "two passes about the cell mean, lib/gridding.py's cell_moments")
        var("f_sat_max", f_sat_max.reshape(shape), "f4", ("lat", "lon"), "1",
            "share of the cell's land AREA whose index exceeds the cell mean, "
            "lib/gridding.py's CATEGORICAL reduction. It multiplies into a "
            "fraction of a cell, so it is a share of area and not of regions, "
            "which on this unequal-area mesh are different quantities. A RANK "
            "statistic, so it survives an additive shift of the index; that is "
            "what makes it transportable where the absolute scale is not")
        var("land_regions", count.reshape(shape), "i4", ("lat", "lon"), "1",
            "land regions in the cell, the sample f_sat_max is computed from")

    report = {
        "script": "hydrography/scripts/build_topographic_index.py",
        "issue": "GW-26",
        "source_build": export.root.parent.name,
        "grid": grid_name,
        "regions": int(export.n_regions),
        "land_regions": int(land.sum()),
        "index": diag,
        "slope_arm": primary,
        "grid_cells_with_estimate": int(have.sum()),
        "min_regions_per_cell": int(args.min_regions),
        # WHAT THE REDUCTION IS AND WHAT IT DROPPED. `f_sat_max` is an area
        # share and the cell mean it is taken about is an area mean, so the
        # ledger is part of the result rather than a diagnostic beside it:
        # a cell holding none of the land population has no fraction, and the
        # closure residual says whether every region that was offered landed.
        "reduction": {
            "cell_mean_index": "gridding.cell_moments, INTENSIVE, area-weighted",
            "within_cell_index_sd": "gridding.cell_moments, two-pass about the cell mean",
            "f_sat_max": "gridding.cell_fraction, CATEGORICAL, share of the cell's land AREA above the cell mean",
            "population": "surface_class == LAND",
            "ledger": ledger,
        },
        "cell_mean_index": {
            str(p): float(np.nanpercentile(cmean[have], p))
            for p in (5, 25, 50, 75, 95)},
        "within_cell_index_sd": {
            str(p): float(np.nanpercentile(csd[have], p)) for p in (5, 50, 95)},
        "between_cell_index_sd": float(np.nanstd(cmean[have])),
        "f_sat_max": {
            str(p): float(np.nanpercentile(f_sat_max[have], p))
            for p in (5, 25, 50, 75, 95)},
        "cell_mean_above_climberx_cdf_top": over,
        "absolute_thresholds": "refused",
        "score": cfg["score"],
        "license": lic,
        "scored": lic["scored"],
        "consumers_licensed": lic["consumers_licensed"],
        "created": datetime.now(timezone.utc).isoformat(),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    rp = ANALYSIS / "topographic_index_report.json"
    rp.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {out}\nwrote {rp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
