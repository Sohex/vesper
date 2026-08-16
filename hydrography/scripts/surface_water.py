"""Where the standing and running water is, under a real climate.

`lake_balance.py` has always been able to solve for lake extent; what it lacked
was forcing. This runs it against the ExoPlaSim baseline climatology through
`coupling_*.nc`, and accumulates the same runoff down the drainage network to
get river discharge. Two products, one water balance:

    data/surface_water.nc   per region: lake, lake depth, river discharge
                            per basin: area, level, volume, overflow

This is the first thing in the project to decide `surface_class == 2`, which
World Orogen deliberately leaves empty. It is a result, not a picture, and it
carries the caveats of its forcing: the baseline run is T42 and predates the
carve, and the lakes it implies are not fed back into it. A world with this much
open water would evaporate more and be cloudier, so the climate that produced
these lakes is not the climate that would exist with them.

    python hydrography/scripts/surface_water.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ANALYSIS, DATA, PROJECT_ROOT  # noqa: E402

import carve_verdict as cv  # noqa: E402
import lake_balance as lb  # noqa: E402
from orogen import LAND, Export  # noqa: E402

import builds  # noqa: E402

CLIMATOLOGY = PROJECT_ROOT / "exoplasim/analysis/climatology_s096"
SECONDS_PER_DAY = 86400.0

# Set by main() from --data. Hydrography products are per build, because drainage
# is a property of the terrain, and this script previously read the flat
# hydrography/data/ whatever build was configured. That silently pairs one
# terrain's basins with another's coupling matrix.
_DATA = DATA


def data_dir() -> Path:
    return _DATA


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def climate_fields(config):
    """Runoff, lake precipitation and open-water evaporation, m/s, on the grid.

    Runoff is P-E, not the model's `mrro`. At steady state the two are the same
    thing: whatever falls on land and does not evaporate has to leave, and the
    baseline run's global budget closes to 0.03% of the mean, with the land's
    surplus matching the sea's deficit to three decimals. Its `mrro` accounts
    for only 15% of that surplus, so the diagnostic is the unreliable one here
    and P-E is what the water balance can be built on. Both are reported.

    Lake evaporation is `carve_verdict.penman_open_water`, the same estimate the
    carve verdict is decided on, so the lakes and the terrain they sit in are
    judged by one rule. It is the one that can be checked: applied to ocean
    cells, which already are open water, Penman reproduces the model's own
    evaporation to 2.5%. A Priestley-Taylor estimate stood here first and came
    out 7.8% high on that same test, which is three times the error for no
    gain, and it duplicated a validated function in this directory.
    """
    clim = globals().get("_CLIM_FILE") or (
        CLIMATOLOGY / "baseline_regular_climatology.nc")
    with Dataset(clim) as ds:
        pr = cv.annual_mean(ds, "pr")
        evap = -cv.annual_mean(ds, "evap")     # code 182 is negative upward
        mrro = cv.annual_mean(ds, "mrro")
        rss, rls = cv.annual_mean(ds, "rss"), cv.annual_mean(ds, "rls")
        ts, tas = cv.annual_mean(ds, "ts"), cv.annual_mean(ds, "tas")
        ps_pa = cv.annual_mean(ds, "ps") * 100.0
        q_air = np.asarray(ds["hus"][:]).mean(axis=0)[-1]
        wind = np.asarray(ds["spd"][:]).mean(axis=0)[-1]
        lat = np.asarray(ds["lat"][:])
        lon = np.asarray(ds["lon"][:])

    runoff = np.clip(pr - evap, 0.0, None)

    resolution = str(config["model"]["resolution"]).upper()
    land_albedo = cv.read_sra_field(
        PROJECT_ROOT / "exoplasim" / "inputs" / resolution.lower()
        / f"orogen_{resolution}_surf_0174.sra", *ps_pa.shape)
    evaporation = cv.penman_open_water(
        ts, tas, q_air, wind, ps_pa, rss, rls, land_albedo,
        float(config["planet"]["gravity_m_s2"]))
    # Same floor the verdict uses: Penman linearises around air temperature and
    # can fall below the model's own evaporation where the ground runs hotter,
    # which is impossible for a saturated surface under the same forcing.
    evaporation = np.maximum(evaporation, evap)
    return lat, lon, runoff, np.clip(pr, 0.0, None), evaporation, np.clip(mrro, 0.0, None)


def per_basin_forcing(n_basins, lat, lon, runoff, precip, evaporation, sinks, export):
    """Catchment-mean runoff, and lake fluxes taken at each basin's sink.

    Runoff has to be integrated over the catchment, which is what the coupling
    matrix is for. The lake fluxes do not: they act on the water surface, which
    sits at the sink, so taking them from the catchment mean would charge a lake
    in a desert with the evaporation of the mountains that feed it.
    """
    # cv.basin_means rather than a second aggregation here: it is the same
    # sum, and it is where the coupling grid's longitude convention is
    # reconciled with the climatology's. Doing it twice is how they came to
    # disagree in the first place.
    means, _ = cv.basin_means(data_dir() / "coupling_exoplasim-T42.nc",
                              {"runoff": runoff}, n_basins, field_lon=lon)
    catchment_runoff = np.nan_to_num(means["runoff"])

    nlat, nlon = runoff.shape
    sink_lat = export.field("lat")[sinks]
    sink_lon = export.field("lon")[sinks]
    row = np.abs(lat[None, :] - sink_lat[:, None]).argmin(axis=1)
    col = np.mod(np.round((sink_lon - lon[0]) / (360.0 / nlon)).astype(int), nlon)
    return catchment_runoff, precip[row, col], evaporation[row, col]


def region_grid_cells(export, field_lon, field_lat):
    """Per-region climate cell, binned exactly as the coupling matrix bins.

    The coupling matrix assigns a region to the cell containing its centre on
    the Orogen grid, whose columns are offset half a cell from an ExoPlaSim
    climatology's and run -180 to 180 against its 0 to 360. Sampling per region
    by nearest centre instead put the mesh and the coupling on different cells,
    so the same basin got a different runoff depending on which route was taken:
    a 36% spread at the ninetieth percentile. One binning for both.
    """
    grid_dir = builds.grid_export()
    gm = json.loads((grid_dir / "manifest.json").read_text())["grid"]
    glat = np.fromfile(grid_dir / gm["coords"]["lat"]["path"], dtype="float64")
    glon = np.fromfile(grid_dir / gm["coords"]["lon"]["path"], dtype="float64")
    nlat, nlon = glat.size, glon.size

    edges = np.empty(nlat + 1)
    edges[1:-1] = 0.5 * (glat[:-1] + glat[1:])
    edges[0], edges[-1] = 90.0, -90.0
    row = np.clip(np.searchsorted(-edges, -export.lat, side="right") - 1, 0, nlat - 1)
    col = np.clip(((export.lon + 180.0) / 360.0 * nlon).astype(np.int64), 0, nlon - 1)

    def wrap(a):
        return (np.asarray(a) + 180.0) % 360.0 - 180.0

    remap = np.abs(wrap(field_lon)[None, :] - wrap(glon)[:, None]).argmin(axis=1)
    rows = np.abs(np.asarray(field_lat)[None, :] - glat[:, None]).argmin(axis=1)
    return rows[row], remap[col]


def paint_lakes(terminal, filled_km, area_km2, solved_area_km2):
    """Which regions lie under water, filling each basin to its solved area.

    Thresholding on the solved level instead would paint whole regions, and a
    region is 230 km2. Basins whose lake is smaller than that would each get a
    15 km blob, so a wet landscape of many small basins comes out speckled with
    lakes that are mostly rounding error. Filling by area in ascending order of
    the flooded surface is the same curve the hypsometry was built from, and it
    stops a basin at what it actually holds: a lake too small to reach a whole
    region simply is not drawn.
    """
    wet = np.zeros(terminal.size, bool)
    sel = np.flatnonzero(terminal >= 0)
    order = sel[np.lexsort((filled_km[sel], terminal[sel]))]
    groups = terminal[order]
    bounds = np.searchsorted(groups, np.arange(solved_area_km2.size + 1))

    cumulative = np.cumsum(area_km2[order])
    start = np.zeros(order.size)
    for b in range(solved_area_km2.size):
        lo, hi = bounds[b], bounds[b + 1]
        if hi > lo:
            start[lo:hi] = cumulative[lo] - area_km2[order[lo]]
    within = cumulative - start
    budget = np.zeros(order.size)
    for b in range(solved_area_km2.size):
        budget[bounds[b]:bounds[b + 1]] = solved_area_km2[b]
    wet[order] = within <= budget
    return wet


def route_overflow(discharge, receiver, spill_exit, overflow_m3_s, land):
    """Carry each overflowing basin's outflow from its saddle to the next sink.

    A basin pinned at its spill passes water on, and until the saddle was
    recorded there was nowhere on the mesh to start that river: the level was
    known and the place was not. So the largest flows on the planet were absent
    from the network entirely, over a million m3/s of it.

    No double counting: the solver has already resolved the cascade, so a
    basin's overflow is its final equilibrium value with everything upstream
    included, and the paths are disjoint segments. Basin A's water runs from
    A's saddle to B's sink, and B's from B's saddle onward.
    """
    rec = receiver.tolist()
    flow = discharge.tolist()
    routed = 0
    for b in np.flatnonzero(overflow_m3_s > 0):
        exit_region = int(spill_exit[b])
        if exit_region < 0 or not land[exit_region]:
            continue          # the saddle opens straight onto the sea
        q = float(overflow_m3_s[b])
        r = exit_region
        seen = 0
        while r >= 0:
            flow[r] += q
            r = rec[r]
            seen += 1
            if seen > receiver.size:
                raise SystemExit("overflow routing looped; the tree is not a tree")
        routed += 1
    return np.array(flow), routed


def river_discharge(export, receiver, runoff_per_region):
    """Accumulate runoff down the drainage tree.

    Every region drains to `receiver`, a forest rooted at the ocean margin and
    at the basin sinks, so this is one pass from the leaves down. The order is
    found by peeling leaves rather than by sorting on elevation: inside a filled
    depression the surface is flat, so elevation does not order a region against
    the one it drains into, and sorting there would drop whole tributaries.

    A region's own total must be complete before it is passed on, so this is
    sequential by construction and `np.add.at` cannot do it. Lists rather than
    arrays, because the loop is two and a half million iterations long.
    """
    land = export.surface_class == LAND
    discharge = np.where(land, runoff_per_region, 0.0)

    children = np.zeros(export.n_regions, dtype=np.int64)
    has_receiver = receiver >= 0
    np.add.at(children, receiver[has_receiver], 1)

    rec = receiver.tolist()
    kids = children.tolist()
    flow = discharge.tolist()
    stack = np.flatnonzero(has_receiver & (children == 0)).tolist()
    processed = 0
    while stack:
        r = stack.pop()
        d = rec[r]
        flow[d] += flow[r]
        processed += 1
        kids[d] -= 1
        if kids[d] == 0 and rec[d] >= 0:
            stack.append(d)
    if processed != int(has_receiver.sum()):
        raise SystemExit("the drainage tree has a cycle; accumulation is invalid")
    return np.array(flow)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=None,
                    help="hydrography products for this build; defaults to "
                         "data/<source_build>/ when it exists")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="baseline_regular_climatology.nc to force the lakes with")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    global _DATA, CLIMATOLOGY
    import yaml as _yaml
    _cfg = _yaml.safe_load((PROJECT_ROOT / "config/planet.yaml").read_text())
    if args.data is not None:
        _DATA = args.data
    else:
        # Strict: no fallback to the flat directory. Falling back when the
        # per-build file is missing is the same trap one level down -- it turns
        # a missing input into a silent read of a terrain nobody chose.
        _DATA = builds.component_data("hydrography", _cfg, strict=True)
    if args.climatology is not None:
        # Resolve before storing: the provenance write takes relative_to
        # PROJECT_ROOT, which raises on a path given relative to the cwd.
        clim_path = args.climatology.resolve()
        CLIMATOLOGY = clim_path.parent
        globals()["_CLIM_FILE"] = clim_path

    build = builds.build_root()
    export = Export(builds.mesh_export())
    n = export.n_regions

    with Dataset(data_dir() / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:])
        filled_km = np.asarray(ds["filled_km"][:])
        if "receiver" not in ds.variables:
            raise SystemExit(
                "regions.nc has no receiver field. Re-run build_hydrography.py."
            )
        receiver = np.asarray(ds["receiver"][:])

    basins = lb.BasinSet(data_dir() / "basins.nc")
    if basins.terrain_hash != export.terrain_hash:
        raise SystemExit("basins.nc was built from a different terrain")

    import orbit
    import yaml

    config = yaml.safe_load((PROJECT_ROOT / "config/planet.yaml").read_text())
    print("reading the climatology")
    lat, lon, runoff, precip, evaporation, model_runoff = climate_fields(config)
    sinks = np.array([b.sink for b in export.basins])
    catchment_runoff, lake_precip, lake_evap = per_basin_forcing(
        basins.n, lat, lon, runoff, precip, evaporation, sinks, export
    )

    # The solver works in km/year; the year is this world's, from lib/orbit.
    year_days = orbit.orbital_year_days(config)
    per_year = year_days * SECONDS_PER_DAY
    to_km_per_year = per_year / 1000.0

    print(f"solving lake levels ({year_days:.1f}-day year)")
    solution = lb.solve(
        basins,
        catchment_runoff * to_km_per_year,
        lake_evap * to_km_per_year,
        lake_precip * to_km_per_year,
    )
    if not solution["converged"]:
        raise SystemExit("the overflow cascade did not converge; do not use this")

    level = solution["level_km"]
    area = export.cell_area.astype(np.float64)
    wet = paint_lakes(terminal, filled_km, area, solution["area_km2"])
    lake_depth = np.where(wet, level[np.maximum(terminal, 0)] - filled_km, 0.0)

    planet_km2 = float(area.sum())
    lake_km2 = float(area[wet].sum())
    print(f"  {int(solution['area_km2'].size - solution['dry'].sum())} basins hold water, "
          f"{int(solution['fills_to_spill'].sum())} to their spill")
    print(f"  lakes cover {lake_km2 / planet_km2 * 100:.2f}% of the planet, "
          f"painted area {lake_km2:,.0f} km2 against a solved "
          f"{solution['area_km2'].sum():,.0f} km2")

    print("accumulating rivers")
    row, col = region_grid_cells(export, lon, lat)
    runoff_per_region = runoff[row, col] * area * 1e6  # m3/s, area km2 to m2
    discharge = river_discharge(export, receiver, runoff_per_region)
    before = discharge.max()

    with Dataset(data_dir() / "basins.nc") as ds:
        spill_exit = np.asarray(ds["spill_exit_region"][:])
    overflow_m3_s = solution["overflow_km3_per_year"] * 1e9 / per_year
    discharge, routed = route_overflow(
        discharge, receiver, spill_exit, overflow_m3_s, export.surface_class == LAND)
    print(f"  routed the outflow of {routed} overflowing basins from their saddles")
    print(f"  largest river {discharge.max():,.0f} m3/s, against {before:,.0f} "
          f"before the overflow was routed; "
          f"{int((discharge > 1000).sum()):,} regions above 1000 m3/s")

    out = args.output or (data_dir() / "surface_water.nc")
    with Dataset(out, "w") as ds:
        ds.createDimension("region", n)
        ds.createDimension("basin", basins.n)
        ds.title = "Lakes and rivers under the baseline climatology"
        ds.terrain_hash = export.terrain_hash
        ds.forcing = str(CLIMATOLOGY.relative_to(PROJECT_ROOT))
        ds.caveat = (
            "The forcing is a T42 run on the pre-carve terrain and these lakes "
            "are not fed back into it. Open-water evaporation is the Penman "
            "combination, the same estimate the carve verdict uses."
        )
        for name, data, dtype, dim, units, note in [
            ("lake", wet, "i1", "region", "1", "region lies under a lake"),
            ("lake_depth_km", lake_depth, "f4", "region", "km",
             "lake surface less the filled terrain"),
            ("discharge_m3_s", discharge, "f4", "region", "m3 s-1",
             "runoff accumulated down the drainage tree"),
            ("basin_area_km2", solution["area_km2"], "f8", "basin", "km2",
             "solved equilibrium lake area"),
            ("basin_level_km", level, "f8", "basin", "km", "solved lake surface"),
            ("basin_volume_km3", solution["volume_km3"], "f8", "basin", "km3", ""),
            ("basin_overflow_km3_yr", solution["overflow_km3_per_year"], "f8", "basin",
             "km3 yr-1", "passed downstream once pinned at spill"),
            ("basin_fills_to_spill", solution["fills_to_spill"], "i1", "basin", "1", ""),
        ]:
            v = ds.createVariable(name, dtype, (dim,), zlib=True)
            v.units = units
            if note:
                v.long_name = note
            v[:] = data

    with Dataset(CLIMATOLOGY / "baseline_regular_climatology.nc") as ds:
        land_mask = cv.annual_mean(ds, "lsm") > 0.5
        model_evap = -cv.annual_mean(ds, "evap")
    area_weight_all = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, runoff.shape[1]))
    area_weight_sea = area_weight_all * ~land_mask
    cell_weight = area_weight_all * land_mask
    report = {
        "source_build": build.name,
        "terrain_hash": export.terrain_hash,
        "forcing": str(CLIMATOLOGY.relative_to(PROJECT_ROOT)),
        "forcing_sha256": sha256(CLIMATOLOGY / "baseline_regular_climatology.nc"),
        "orbital_year_days": year_days,
        "runoff_source": {
            "field": "P-E from the baseline climatology, not mrro",
            "why": ("the run's global water budget closes to 0.03% of the mean "
                    "and the land surplus matches the sea deficit, but mrro "
                    "accounts for only 15% of that surplus"),
            "land_mean_p_minus_e_mm_per_day": float(
                np.average(runoff, weights=cell_weight) * SECONDS_PER_DAY * 1000),
            "land_mean_mrro_mm_per_day": float(
                np.average(model_runoff, weights=cell_weight) * SECONDS_PER_DAY * 1000),
        },
        "open_water_evaporation": {
            "method": ("Penman combination, water albedo and roughness, floored at "
                       "the model's land rate; shared with carve_verdict.py"),
            "global_mean_mm_per_day": float(
                np.average(evaporation, weights=area_weight_all) * SECONDS_PER_DAY * 1000),
            "ocean_mean_mm_per_day": float(
                np.average(evaporation, weights=area_weight_sea) * SECONDS_PER_DAY * 1000),
            "model_ocean_mean_mm_per_day": float(
                np.average(model_evap, weights=area_weight_sea) * SECONDS_PER_DAY * 1000),
            "validation": ("ratio over ocean cells, which already are open water; "
                           "the only place the estimate can be checked"),
        },
        "lakes": {
            "basins_holding_water": int(basins.n - solution["dry"].sum()),
            "basins_at_spill": int(solution["fills_to_spill"].sum()),
            "solved_area_km2": float(solution["area_km2"].sum()),
            "painted_area_km2": lake_km2,
            "planet_fraction": lake_km2 / planet_km2,
            "total_volume_km3": float(solution["volume_km3"].sum()),
        },
        "rivers": {
            "max_discharge_m3_s": float(discharge.max()),
            "basins_whose_outflow_was_routed": int(routed),
            "overflow_routed_m3_s": float(overflow_m3_s.sum()),
            "regions_above_1000_m3_s": int((discharge > 1000).sum()),
            "land_runoff_km3_per_year": float(
                runoff_per_region[export.surface_class == LAND].sum() * per_year / 1e9),
            "land_runoff_m3_s": float(
                runoff_per_region[export.surface_class == LAND].sum()),
        },
        "git_commit": subprocess.run(
            ["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True).stdout.strip(),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    (ANALYSIS / "surface_water_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
