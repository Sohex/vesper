"""Where the standing and running water is, under a real climate.

`lake_balance.py` has always been able to solve for lake extent; what it lacked
was forcing. This runs it against whichever ExoPlaSim climatology
`config/planet.yaml` names, through `coupling_*.nc`, and accumulates the same
runoff down the drainage network to get river discharge. Two products, one water
balance:

    data/<build>/surface_water.nc   per region: lake, lake depth, river discharge,
                            and the periodic lake cycle painted onto regions
                            per basin: area, level, volume, overflow, and the
                            periodic steady state through the climatology's bins

This is the first thing in the project to decide `surface_class == 2`, which
World Orogen deliberately leaves empty. It is a result, not a picture, and it
carries the caveats of its forcing, and the forcing is named on the artifact
rather than assumed: a bootstrap climatology gives bootstrap lakes, and its
numbers are not the baseline. The lakes are not fed back into the run that
produced them either. A world with this much open water would evaporate more and
be cloudier, so the climate that produced these lakes is not the climate that
would exist with them.

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
from paths import best_available_climatology, rel  # noqa: E402
from provenance import staged_surface_field  # noqa: E402
from write_door import refuse_a_write_through_a_symlink  # noqa: E402

import builds  # noqa: E402
import climatology  # noqa: E402
import gridding  # noqa: E402
import nc_geometry  # noqa: E402

# Set by main(), from lib/paths.py:best_available_climatology or
# --climatology. There is no module-level default on purpose; see main().
_CLIM_FILE = None
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


def climate_fields(config, bin_index=None):
    """Runoff, lake precipitation and open-water evaporation, m/s, on the grid.

    `bin_index` evaluates the LAKE terms on one of the climatology's time bins
    instead of on the annual mean, which is what the periodic lake balance is
    integrated against. Everything acting on the water surface -- open-water
    evaporation and the precipitation falling on it -- can be taken per bin,
    because a surface flux has no storage between bins. Catchment runoff cannot;
    see the note where it is computed. With `bin_index` absent this returns
    exactly what it always did.

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
    evaporation to within a few percent, and the run's own figure is written
    into `carve_verdict.json` rather than quoted here. A Priestley-Taylor
    estimate stood here first and came out several times worse on that same
    test for no gain, and it duplicated a validated function in this directory.
    """
    clim = _CLIM_FILE
    if clim is None:
        raise RuntimeError(
            "no climatology resolved; main() sets it from "
            "the best available climatology or --climatology")
    with Dataset(clim) as ds:
        if bin_index is None:
            def take(name):
                return cv.annual_mean(ds, name)
        else:
            k = int(bin_index)

            def take(name):
                return np.asarray(ds[name][k])
        pr = take("pr")
        evap = -take("evap")                   # code 182 is negative upward
        mrro = take("mrro")
        rss, rls = take("rss"), take("rls")
        lat = np.asarray(ds["lat"][:])
        lon = np.asarray(ds["lon"][:])
        lsm = take("lsm")
        diurnal = take("maxt") - take("mint")
        annual_pe = cv.annual_mean(ds, "pr") + cv.annual_mean(ds, "evap")
    t_air, q_air, wind, p_air = cv.reference_level_air(clim, bin_index)

    # RUNOFF STAYS ANNUAL EVEN IN A BIN, and it is not an oversight. Over one
    # cycle at steady state a land cell's storage returns to where it started,
    # so annual(P - E) is exactly the runoff that cell generated; clamping the
    # negative bins away instead counts the wet season's supply twice, once as
    # runoff and once as the water that refilled the soil the dry season
    # emptied. `carve_verdict.seasonal_rectification` measures what that would
    # add. So the value here is the same annual number in every bin, and the
    # seasonal PHASE of catchment delivery is a declared absence in
    # `config/land_water_ledger.yaml` rather than a term computed from a
    # quantity that cannot carry it.
    runoff = np.clip(annual_pe, 0.0, None)

    # Through the one door. The path this used to build is keyed by the RUNG
    # alone while `surface_albedo` rewrites it per BUILD, so it could not say
    # which build's field it was reading. No cross-build read here: surface
    # water is computed on the build the config names. world-z7bu, rule 5.
    staged_albedo = staged_surface_field(174, config)
    land_albedo = cv.read_sra_field(PROJECT_ROOT / staged_albedo["path"],
                                    *p_air.shape)
    evaporation = cv.penman_open_water(
        t_air, q_air, wind, p_air, rss, rls, land_albedo,
        float(config["planet"]["gravity_m_s2"]), diurnal_range=diurnal,
        cfg=config)
    # No floor at the land rate, for the reason `carve_verdict.py` sets out at
    # length: land here is 6.4x aerodynamically rougher than water, so a lake
    # evaporating less than the wet ground around it is physical rather than a
    # failure of the estimate. Flooring it made every lake in a wet catchment
    # evaporate at the land rate.
    return (lat, lon, runoff, np.clip(pr, 0.0, None), evaporation,
            np.clip(mrro, 0.0, None), lsm, staged_albedo)


def per_basin_forcing(n_basins, lat, lon, runoff, precip, evaporation, sinks,
                      export, lsm, config):
    """Catchment-mean runoff, and lake fluxes taken at each basin's sink.

    Runoff has to be integrated over the catchment, which is what the coupling
    matrix is for. The lake fluxes do not: they act on the water surface, which
    sits at the sink, so taking them from the catchment mean would charge a lake
    in a desert with the evaporation of the mountains that feed it.
    """
    # cv.basin_means rather than a second aggregation here: it is the same sum,
    # and it is the one place that owns how a coupling column maps to a
    # climatology column. Doing it twice is how they came to disagree in the
    # first place, and the mapping is the identity -- see basin_means.
    coupling = cv.coupling_path(data_dir(), config)
    cv.require_index_alignment(coupling, lsm)
    means, _ = cv.basin_means(coupling, {"runoff": runoff}, n_basins)
    catchment_runoff = np.nan_to_num(means["runoff"])

    # The sink's own cell, through the module that owns the convention.
    #
    # This line measured an Orogen longitude from the CLIMATOLOGY's first label,
    # `col = mod(round((sink_lon - lon[0]) / dlon), nlon)`, which is the same
    # 0..360-against--180..180 label match that shifted `basin_means` and
    # `region_grid_cells` by half the grid. It survived the sweep that fixed
    # those two because it did not look like a remap. Every lake on this planet
    # was taking its precipitation and its open-water evaporation from its
    # antipode. The rows are checked rather than joined for the same reason as
    # in `region_grid_cells`.
    nlon = runoff.shape[1]
    gridding.require_same_rows(gridding.grid_geometry(builds.grid_export())[0], lat,
                               "the grid export and the climatology")
    row, col = gridding.cells(export.field("lat")[sinks],
                              export.field("lon")[sinks], lat, nlon)
    return catchment_runoff, precip[row, col], evaporation[row, col]


def region_grid_cells(export, field_lat):
    """Per-region climate cell, binned exactly as the coupling matrix bins.

    Delegates to `gridding.climatology_cells`, which is the one copy of this
    convention. It used to live here, and the reason it moved is that a second
    consumer appeared: the derived-surface classifier needs the same join, and a
    grid convention this project has got wrong twice is the last thing to keep
    two versions of. `lib/gridding.py` owns the convention and the invariant that
    proves it held.
    """
    from gridding import climatology_cells
    return climatology_cells(export, builds.grid_export(), field_lat)


def _fill_order(terminal, filled_km, area_km2, n_basins):
    """The ascending-flooded-surface order each basin is filled in, once.

    Depends on the terrain and the basin membership and NOT on how much water
    is in any basin, so it is the same order in every time bin of a cycle. It
    is hoisted out because the sort is over every region in a basin and
    `paint_lake_cycle` would otherwise repeat it once per bin.

    Returns `(order, within, bounds)`: the region indices grouped by basin and
    ascending in filled surface, the area cumulated within each basin along
    that order, and the group boundaries.
    """
    sel = np.flatnonzero(terminal >= 0)
    order = sel[np.lexsort((filled_km[sel], terminal[sel]))]
    bounds = np.searchsorted(terminal[order], np.arange(n_basins + 1))
    cumulative = np.cumsum(area_km2[order])
    start = np.zeros(order.size)
    for b in range(n_basins):
        lo, hi = bounds[b], bounds[b + 1]
        if hi > lo:
            start[lo:hi] = cumulative[lo] - area_km2[order[lo]]
    return order, cumulative - start, bounds


def _paint(order, within, bounds, n_regions, budget_km2):
    """Regions under water once each basin is filled to its own area budget."""
    wet = np.zeros(n_regions, bool)
    budget = np.zeros(order.size)
    for b in range(budget_km2.size):
        budget[bounds[b]:bounds[b + 1]] = budget_km2[b]
    wet[order] = within <= budget
    return wet


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
    order, within, bounds = _fill_order(terminal, filled_km, area_km2,
                                        solved_area_km2.size)
    return _paint(order, within, bounds, terminal.size, solved_area_km2)


def paint_lake_cycle(terminal, filled_km, area_km2, area_by_bin_km2,
                     bin_weight, cycle_closed):
    """The periodic lake cycle, painted onto regions bin by bin.

    WORLD-9IY5, and it is one operation rather than a new one: the periodic
    solve gives each closed basin a lake AREA in every time bin, and painting
    an area onto regions is exactly what `paint_lakes` does for the annual
    equilibrium. This runs the same fill, once per bin, and reduces the bins to
    per-region quantities that a classification can partition on.

    **THIS IS THE CLOSED-BASIN THIRD OF THE SEASONAL QUANTITY AND NOTHING
    ELSE.** `hydrography/notes/land-water-ledger.md` splits seasonally
    inundated land in three. A FLOODPLAIN's inundated area needs a
    height-above-nearest-drainage distribution and a routing model, neither of
    which exists in this project, and a seasonally saturated SOIL is a water
    content rather than an area. Merging either into what this returns would
    put three quantities under one name, and the merge would be invisible in
    the file.

    **A BASIN WHOSE YEAR DID NOT CLOSE HAS NO CYCLE**, so it is not painted at
    all and `decided` is false over every region that drains to it. That is the
    same refusal `solve_periodic` makes, carried across the crossing rather
    than dropped by it: a reader must be able to tell a region that is dry all
    year from one whose basin was refused, and a zero cannot say both.

    Returns a dict of per-region arrays:
      `cycle_fraction`  share of ONE CYCLE the region spends under water,
                        weighted by the bins' own lengths rather than counted,
                        because the bins need not be equal.
      `bins_wet`        how many bins it is under water, which is the exact
                        integer the permanent/seasonal split is taken on.
      `decided`         its basin closed its year, so the two above mean
                        something.
    """
    n_regions = terminal.size
    nbin, n_basins = np.asarray(area_by_bin_km2).shape
    bin_weight = np.asarray(bin_weight, dtype=float)
    if bin_weight.shape != (nbin,):
        raise ValueError(f"bin_weight is {bin_weight.shape}, expected ({nbin},)")
    closed = np.asarray(cycle_closed, dtype=bool)
    order, within, bounds = _fill_order(terminal, filled_km, area_km2, n_basins)

    cycle_fraction = np.zeros(n_regions)
    bins_wet = np.zeros(n_regions, np.int16)
    for k in range(nbin):
        # A refused basin gets a zero budget, so nothing of it is painted in
        # any bin and `decided` below is what says so.
        budget = np.where(closed, np.asarray(area_by_bin_km2[k], dtype=float), 0.0)
        wet_k = _paint(order, within, bounds, n_regions, budget)
        cycle_fraction += bin_weight[k] * wet_k
        bins_wet += wet_k
    decided = np.zeros(n_regions, bool)
    in_basin = terminal >= 0
    decided[in_basin] = closed[terminal[in_basin]]
    return {"cycle_fraction": cycle_fraction, "bins_wet": bins_wet,
            "decided": decided, "bins": nbin}


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


def _selftest() -> int:
    """The crossing, against identities the painting must satisfy.

    Every check can fail and each has a right answer fixed by construction
    rather than by a result. The first is the one that makes the crossing the
    same operation as the annual one rather than a second implementation of it.
    """
    problems: list[str] = []
    n_checks = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal n_checks
        n_checks += 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    rng = np.random.default_rng(20260825)
    n_regions, n_basins, nbin = 600, 5, 6
    terminal = rng.integers(-1, n_basins, n_regions)
    filled_km = rng.uniform(0.0, 1.0, n_regions)
    area_km2 = rng.uniform(50.0, 400.0, n_regions)
    # Basin capacities well inside what the members can hold, so the fill stops
    # partway and the seasonal ring is not empty by construction.
    per_basin = np.array([area_km2[terminal == b].sum() for b in range(n_basins)])
    annual = 0.5 * per_basin
    closed = np.array([True, True, True, False, True])
    weights = rng.dirichlet(np.ones(nbin))

    # ONE BIN AT FULL WEIGHT IS THE ANNUAL PAINT. Same operation, so it must
    # give the same regions, not merely a similar count.
    one = paint_lake_cycle(terminal, filled_km, area_km2, annual[None, :],
                           np.array([1.0]), np.ones(n_basins, bool))
    direct = paint_lakes(terminal, filled_km, area_km2, annual)
    check("one bin at full weight reproduces the annual paint region for region",
          bool(((one["bins_wet"] > 0) == direct).all()),
          f"{int(((one['bins_wet'] > 0) != direct).sum())} regions differ")

    # A FORCING WITH NO SEASON MUST PAINT NO SEASON. A per-bin painter that
    # rounded, or that mixed the order between bins, would leak regions here.
    steady = np.tile(annual, (nbin, 1))
    st = paint_lake_cycle(terminal, filled_km, area_km2, steady, weights,
                          np.ones(n_basins, bool))
    wet_all = st["bins_wet"] == nbin
    check("a steady area gives every painted region the whole cycle",
          bool(((st["bins_wet"] == 0) | wet_all).all())
          and float(np.abs(st["cycle_fraction"][wet_all] - 1.0).max()) < 1e-12,
          f"{int(((st['bins_wet'] != 0) & ~wet_all).sum())} regions part-wet")

    # THE DURATION IS IN CYCLE UNITS, NOT BIN COUNTS. A painter that averaged
    # over bins instead of weighting by their lengths passes every check above
    # and fails this one, and unequal bins are the normal case.
    lop = np.zeros((nbin, n_basins))
    lop[0] = annual
    lo = paint_lake_cycle(terminal, filled_km, area_km2, lop, weights,
                          np.ones(n_basins, bool))
    touched = lo["bins_wet"] > 0
    check("the duration is weighted by the bins' own lengths",
          float(np.abs(lo["cycle_fraction"][touched] - weights[0]).max()) < 1e-12,
          f"got {lo['cycle_fraction'][touched][:3]} for a bin of {weights[0]:.4f}")

    # NESTED BY CONSTRUCTION: the fill order does not depend on the budget, so
    # a bigger lake covers everything a smaller one did. bins_wet must then be
    # exactly how many bins reach the region, which is an independent count.
    varying = annual[None, :] * rng.uniform(0.2, 1.0, (nbin, 1))
    va = paint_lake_cycle(terminal, filled_km, area_km2, varying, weights,
                          np.ones(n_basins, bool))
    tally = np.zeros(n_regions, np.int16)
    for k in range(nbin):
        tally += paint_lakes(terminal, filled_km, area_km2, varying[k])
    check("bins_wet equals the bins that reach the region, counted separately",
          bool((va["bins_wet"] == tally).all()),
          f"{int((va['bins_wet'] != tally).sum())} regions disagree")
    check("a varying forcing does produce a seasonal ring",
          bool(((va["bins_wet"] > 0) & (va["bins_wet"] < nbin)).any()),
          "the fixture has no season to detect")

    # A REFUSED BASIN IS NOT PAINTED AND SAYS SO. A zero that meant both dry
    # and unsolved is the failure this separates.
    ref = paint_lake_cycle(terminal, filled_km, area_km2, steady, weights, closed)
    members = terminal == 3
    check("a basin whose year did not close is refused, not painted dry",
          not bool(ref["bins_wet"][members].any())
          and not bool(ref["decided"][members].any())
          and bool(ref["decided"][terminal == 0].all()),
          "the refusal did not cross")
    check("a region in no basin is never decided",
          not bool(ref["decided"][terminal < 0].any()), "off-basin regions decided")

    print(f"\n{n_checks} checks, {len(problems)} failed")
    return 1 if problems else 0


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=None,
                    help="hydrography products for this build; defaults to "
                         "data/<source_build>/ when it exists")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="regular climatology to force the lakes with; defaults "
                         "to the best available, which is the configured "
                         "baseline_climatology once one is named and the "
                         "bootstrap before that")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="run the basin-to-region crossing's identity checks; "
                         "needs no build and no climatology")
    args = ap.parse_args()

    if args.selftest:
        raise SystemExit(_selftest())

    global _DATA, _CLIM_FILE
    import yaml as _yaml
    _cfg = _yaml.safe_load((PROJECT_ROOT / "config/planet.yaml").read_text())
    if args.data is not None:
        _DATA = args.data
    else:
        # Strict: no fallback to the flat directory. Falling back when the
        # per-build file is missing is the same trap one level down -- it turns
        # a missing input into a silent read of a terrain nobody chose.
        _DATA = builds.component_data("hydrography", _cfg, strict=True)
    # THE BEST AVAILABLE, which is the baseline once one is named and the
    # bootstrap before that. A lake extent is a function of the climate STATE --
    # it is solved from P-E, routed runoff and open-water evaporation and from
    # nothing else -- and `lib/paths.py` makes that the test for which resolver
    # a step takes. The `needs: bootstrap_climatology` edge in
    # config/pipeline.yaml is unchanged and states what must EXIST for a first
    # pass to run.
    #
    # THE CIRCULARITY IS THE LOOP AND NOT A DEFECT. The lakes are one of the
    # surface fields the baseline run is run ON, so on a later pass they are
    # forced by a climate their own extent helped produce. That is true of every
    # derived surface field, and `build_dust.py`, `build_sea_salt.py` and
    # `build_volcanic_sulfate.py` already read the best available for exactly
    # this reason: a step pinned to the bootstrap forever is right on the first
    # pass and holds the loop back on every pass after. This file was the last
    # of the group still pinned, and the pin put it in disagreement with
    # `build_groundwater.py` about which world it was describing -- which its
    # own coupling guard then reported, correctly, as two climates in one water
    # balance.
    #
    # Resolved through `lib/paths.py` rather than out of the config here. This
    # file kept a private copy of the resolution, which is how it came to still
    # name `climatology_s096` -- pre-carve terrain under the superseded k2
    # spectrum -- months after pedology and biosphere had moved. The shared
    # resolver also carries `require_configured_grid`, which a hand-rolled read
    # of the config key does not.
    clim_path, clim_stage = best_available_climatology(
        args.climatology, root=PROJECT_ROOT)
    # NAMED BY WHAT WAS DECLARED, READ BY WHERE THE BYTES ARE. A worktree links
    # the climatology payload in from the main checkout, so `resolve()` follows
    # the link out of the project root and `rel()` can only hand back an
    # absolute path -- which would be stamped on `surface_water.nc` as its
    # `forcing` and would name a file that exists in no other tree.
    _CLIM_NAME = rel(args.climatology if args.climatology is not None
                     else clim_path)
    clim_path = clim_path.resolve()
    if not clim_path.is_file():
        raise SystemExit(
            f"config/planet.yaml names {_CLIM_NAME} as the {clim_stage} "
            "climatology and it does not exist")
    _CLIM_FILE = clim_path

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
    (lat, lon, runoff, precip, evaporation, model_runoff, lsm,
     staged_albedo) = climate_fields(config)
    sinks = np.array([b.sink for b in export.basins])
    catchment_runoff, lake_precip, lake_evap = per_basin_forcing(
        basins.n, lat, lon, runoff, precip, evaporation, sinks, export, lsm,
        config
    )

    # The solver works in km/year; the year is this world's, from lib/orbit.
    year_days = orbit.orbital_year_days(config)
    per_year = year_days * SECONDS_PER_DAY
    to_km_per_year = per_year / 1000.0

    area = export.cell_area.astype(np.float64)

    # --- the lake forcing, per climatology bin ----------------------------
    # THE LAKE TERMS ARE READ PER BIN AND THE ANNUAL SOLVE READS THEIR MEAN, so
    # the two solves below answer two questions about ONE forcing rather than
    # answering them on two. `land_water_ledger.yaml` puts
    # `open_water_evaporation` at `interval_floor: climatology_bin` and the
    # reason it gives is Jensen: Penman is nonlinear in the air it reads, which
    # is why it already integrates the DIURNAL cycle instead of reading a daily
    # mean, and the seasonal cycle is the same argument at a longer period.
    # Reading the annual mean of the air and evaluating Penman once on it
    # understates open-water evaporation over this catalogue by 17%, and lake
    # extent is what that lands on. The bin-weighted mean of the per-bin
    # evaporation is reported beside the single annual evaluation so the size of
    # that is on the artifact rather than only in a note.
    #
    # The CATCHMENT term stays annual and is not averaged from bins, for the
    # separate reason `climate_fields` gives where it computes it.
    print("reading the lake terms per climatology bin")
    with Dataset(_CLIM_FILE) as _ds:
        bin_weights = climatology.bin_weights(np.asarray(_ds["time"][:]))
    nbin = bin_weights.size
    season_evap = np.zeros((nbin, basins.n))
    season_precip = np.zeros((nbin, basins.n))
    for k in range(nbin):
        # Eight values, and the staged albedo is the one this path discards:
        # the lake balance reads no albedo. Unpack the whole tuple rather than a
        # prefix of it, so a ninth return raises here instead of silently
        # shifting every name one place along.
        (f_lat, f_lon, f_run, f_pr, f_ev, _, f_lsm,
         _) = climate_fields(config, bin_index=k)
        _, p_k, e_k = per_basin_forcing(basins.n, f_lat, f_lon, f_run, f_pr, f_ev,
                                        sinks, export, f_lsm, config)
        season_precip[k] = p_k * to_km_per_year
        season_evap[k] = e_k * to_km_per_year
    bin_share = bin_weights / bin_weights.sum()
    mean_evap = (bin_share[:, None] * season_evap).sum(0)
    mean_precip = (bin_share[:, None] * season_precip).sum(0)
    jensen = float(mean_evap.sum() / max((lake_evap * to_km_per_year).sum(), 1e-30))
    print(f"  open-water evaporation over the basin sinks is {jensen:.3f}x what "
          "one evaluation on the annual-mean air gives")

    print(f"solving lake levels ({year_days:.1f}-day year)")
    solution = lb.solve(
        basins,
        catchment_runoff * to_km_per_year,
        mean_evap,
        mean_precip,
    )
    if not solution["converged"]:
        raise SystemExit("the overflow cascade did not converge; do not use this")

    # THE OTHER END OF THE INTERVAL BRACKET, solved and reported rather than
    # carried as a caveat. `cv._INTERVAL_BRACKET` has the argument: the bin mean
    # is the limit for a lake with no heat storage and the annual evaluation the
    # limit for one deep enough to hold its temperature through the year, and
    # what sits between is the lake's own depth, which this project has no term
    # for. The carve verdict carves the intersection over both ends because its
    # product is an instruction that changes the terrain; lake extent is output,
    # so what is owed here is the SIZE of the spread on the artifact.
    solution_annual = lb.solve(
        basins,
        catchment_runoff * to_km_per_year,
        lake_evap * to_km_per_year,
        lake_precip * to_km_per_year,
    )
    interval_bracket = {
        "solved_area_km2_bin_mean": float(solution["area_km2"].sum()),
        "solved_area_km2_annual_evaluation": float(solution_annual["area_km2"].sum()),
        "annual_evaporation_converged": bool(solution_annual["converged"]),
        "what_it_is": cv._INTERVAL_BRACKET["what_is_bracketed"],
        "missing_term": cv._INTERVAL_BRACKET["missing_term"],
        "which_end_is_used_here": ("the bin mean, matching the carve verdict's "
                                   "primary bound. The spread is the number to "
                                   "read, not a second answer"),
    }
    interval_bracket["spread_as_share_of_bin_mean"] = float(
        abs(interval_bracket["solved_area_km2_annual_evaluation"]
            - interval_bracket["solved_area_km2_bin_mean"])
        / max(interval_bracket["solved_area_km2_bin_mean"], 1e-30))
    print(f"  the evaporation interval brackets lake area at "
          f"{interval_bracket['solved_area_km2_bin_mean']:,.0f} to "
          f"{interval_bracket['solved_area_km2_annual_evaluation']:,.0f} km2 "
          f"({interval_bracket['spread_as_share_of_bin_mean'] * 100:.0f}% spread)")

    level = solution["level_km"]

    # --- the seasonal cycle, as a periodic steady state -------------------
    # WORLD-15I0. The annual solve above answers where a basin settles if the
    # forcing never changes; this integrates the same balance through the
    # climatology's own bins and asks for the YEAR to close on itself.
    #
    # Bin lengths in the SAME year the fluxes are per, which here is Vesper's,
    # so the normalised weights are already that. See solve_periodic.
    print("solving the periodic steady state")
    periodic = lb.solve_periodic(
        basins,
        np.tile(catchment_runoff * to_km_per_year, (nbin, 1)),
        season_evap,
        season_precip,
        bin_weights,
        area_resolution_km2=float(area.mean()),
    )
    n_refused = int((~periodic["closed"] & basins.has_impoundment).sum())
    print(f"  closed in {periodic['cycles']} cycles; "
          f"{int(periodic['seasonal'].sum())} basins carry a publishable season, "
          f"{n_refused} refused for not closing")

    # DOES THE WATER CLOSE, which is not the same question as whether the cycle
    # repeats. A basin's storage can return to where it started every year while
    # the integration quietly creates or destroys water inside the year, and
    # `closed` cannot see that. The residual is read against the throughput
    # rather than reported bare: a km3 means nothing without the km3 that passed
    # through the basin it came off.
    live_basin = basins.has_impoundment
    throughput = periodic["annual_throughput_km3"]
    scale = np.maximum(np.maximum(basins.capacity_km3, np.abs(throughput)), 1e-30)
    # Named for what they are and NOT `per_bin` and `per_year`: `per_year` is
    # already the seconds in this world's year in this function, and shadowing
    # it here divided every basin's overflow by a residual.
    bin_residual_rel = np.abs(periodic["bin_balance_residual_km3"]).max(axis=0) / scale
    year_residual_rel = np.abs(periodic["annual_water_residual_km3"]) / scale
    water_balance = {
        "declared_tolerance_relative": lb.BALANCE_RELATIVE,
        "worst_bin_residual_relative": float(bin_residual_rel[live_basin].max()),
        "worst_annual_residual_relative": float(year_residual_rel[live_basin].max()),
        "set_residual_relative": float(
            np.abs(periodic["annual_water_residual_km3"]).sum()
            / max(throughput.sum(), 1e-30)),
        "annual_throughput_km3": float(throughput.sum()),
        "read_against": ("the larger of the basin's own capacity and what passed "
                         "through it in the year"),
    }
    print(f"  water closes to {water_balance['worst_annual_residual_relative']:.2g} "
          f"per basin and {water_balance['set_residual_relative']:.2g} over the set, "
          f"against a declared {lb.BALANCE_RELATIVE:g}")
    if max(water_balance["worst_bin_residual_relative"],
           water_balance["worst_annual_residual_relative"]) > lb.BALANCE_RELATIVE:
        raise SystemExit(
            "the lake balance does not conserve water to the tolerance declared "
            "in lake_balance.BALANCE_RELATIVE; do not use this")

    wet = paint_lakes(terminal, filled_km, area, solution["area_km2"])
    lake_depth = np.where(wet, level[np.maximum(terminal, 0)] - filled_km, 0.0)

    # THE CROSSING FROM BASIN TO REGION. The periodic solve is per basin and
    # every consumer of a seasonal wetness is per region, so the cycle is
    # painted the same way the annual solve is, once per bin. world-9iy5.
    cycle = paint_lake_cycle(terminal, filled_km, area,
                             periodic["area_km2"], bin_weights,
                             periodic["closed"])
    seasonal_region = (cycle["decided"] & (cycle["bins_wet"] > 0)
                       & (cycle["bins_wet"] < nbin))
    permanent_region = cycle["decided"] & (cycle["bins_wet"] == nbin)
    print(f"  painted the cycle onto regions: {int(permanent_region.sum()):,} "
          f"under water in every bin, {int(seasonal_region.sum()):,} in some "
          f"bins and not all, over {int(cycle['decided'].sum()):,} regions "
          "whose basin closed its year")

    planet_km2 = float(area.sum())
    lake_km2 = float(area[wet].sum())
    print(f"  {int(solution['area_km2'].size - solution['dry'].sum())} basins hold water, "
          f"{int(solution['fills_to_spill'].sum())} to their spill")
    print(f"  lakes cover {lake_km2 / planet_km2 * 100:.2f}% of the planet, "
          f"painted area {lake_km2:,.0f} km2 against a solved "
          f"{solution['area_km2'].sum():,.0f} km2")

    print("accumulating rivers")
    row, col = region_grid_cells(export, lat)
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
    refuse_a_write_through_a_symlink(
        out, what="the solved lake and river field the carve verdict reads",
        instead=("Pass --output to a path inside this worktree, or run the "
                 "generator in the main checkout."))
    with Dataset(out, "w") as ds:
        ds.createDimension("region", n)
        ds.createDimension("basin", basins.n)
        ds.createDimension("time_bin", nbin)
        # Named from the forcing rather than asserted. A bootstrap run's numbers
        # are NOT the baseline, `--climatology` can point this anywhere, and
        # stamping "baseline" on a bootstrap-forced lake set is how a limit
        # comes to be quoted as a state.
        ds.title = f"Lakes and rivers under {_CLIM_FILE.stem}"
        ds.terrain_hash = export.terrain_hash
        ds.setncattr("vesper_source_build", build.name)
        ds.forcing = _CLIM_NAME
        # WHICH STAGE OF THE WORLD THIS IS, recorded rather than inferred from
        # the filename. `lib/paths.py` requires it of anything built through
        # `best_available_climatology`: the choice between stages is defensible
        # because a reader can tell a first-pass artifact from a later one
        # without re-deriving it.
        ds.climatology_stage = clim_stage
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
            # The periodic cycle, crossed to regions. world-9iy5. CLOSED-BASIN
            # LAKE AREA ONLY: not a floodplain's inundation, which needs a
            # height-above-nearest-drainage distribution and a routing model
            # this project does not have, and not a seasonally saturated soil,
            # which is a water content rather than an area.
            ("lake_cycle_decided", cycle["decided"], "i1", "region", "1",
             "the basin this region drains to closed its year, so the two "
             "fields below mean something. READ THIS FIRST: a zero in them is "
             "otherwise ambiguous between dry all year and never solved"),
            ("lake_cycle_fraction", cycle["cycle_fraction"], "f4", "region", "1",
             "share of one cycle the region spends under the periodic lake, "
             "weighted by the time bins' own lengths"),
            ("lake_cycle_bins_wet", cycle["bins_wet"], "i2", "region", "1",
             "time bins the region is under the periodic lake, out of the "
             "climatology's own. The exact integer a permanent/seasonal split "
             "is taken on; the fraction above is the duration"),
            ("basin_area_km2", solution["area_km2"], "f8", "basin", "km2",
             "solved equilibrium lake area"),
            ("basin_level_km", level, "f8", "basin", "km", "solved lake surface"),
            ("basin_volume_km3", solution["volume_km3"], "f8", "basin", "km3", ""),
            ("basin_overflow_km3_yr", solution["overflow_km3_per_year"], "f8", "basin",
             "km3 yr-1", "passed downstream once pinned at spill"),
            ("basin_fills_to_spill", solution["fills_to_spill"], "i1", "basin", "1", ""),
            # The periodic steady state. WORLD-15I0.
            ("basin_area_swing_km2", periodic["area_swing_km2"], "f8", "basin", "km2",
             "peak-to-trough lake area over the cycle"),
            ("basin_seasonal_amplitude", periodic["seasonal_amplitude"], "f8", "basin",
             "1", "that swing as a share of the cycle-mean area"),
            ("basin_residence_time_years", periodic["residence_time_years"], "f8",
             "basin", "yr", "cycle-mean volume over annual inflow; what sets the "
             "amplitude, and why the product is per basin rather than a headline"),
            ("basin_has_season", periodic["seasonal"], "i1", "basin", "1",
             "the cycle is published for this basin: the swing is at least a tenth "
             "of its mean area AND at least one mean mesh cell. A zero is not a "
             "claim that the basin is steady, it is a refusal to publish a swing "
             "the hypsometric curve cannot express"),
            ("basin_cycle_closed", periodic["closed"], "i1", "basin", "1",
             "the year closed on itself to the declared tolerance. READ THIS "
             "BEFORE the per-bin fields: a zero means the arrays below are a "
             "transient, not a cycle, and the basin is refused"),
        ]:
            v = ds.createVariable(name, dtype, (dim,), zlib=True)
            v.units = units
            if note:
                v.long_name = note
            v[:] = data
        for name, data, units, note in [
            ("bin_weight", bin_weights, "1",
             "share of one cycle each time bin spans, from lib/climatology.py"),
        ]:
            v = ds.createVariable(name, "f8", ("time_bin",), zlib=True)
            v.units, v.long_name = units, note
            v[:] = data
        for name, data, units, note in [
            ("basin_area_by_bin_km2", periodic["area_km2"], "km2",
             "lake area in each time bin of the periodic steady state"),
            ("basin_level_by_bin_km", periodic["level_km"], "km",
             "lake surface in each time bin"),
            ("basin_volume_by_bin_km3", periodic["volume_km3"], "km3",
             "lake volume in each time bin"),
        ]:
            v = ds.createVariable(name, "f8", ("time_bin", "basin"), zlib=True)
            v.units, v.long_name = units, note
            v[:] = data
        # LAST in the block. THE SUPPORT IS THE MESH AND NOT A GRID: the
        # regions are unequal in area, so nothing derived from this file by
        # a reduction over the region axis is an area quantity unless it
        # carries the export's region areas. The declaration says so and
        # names where they live. lib/nc_geometry.py.
        nc_geometry.declare_region_mesh(
            ds, n_regions=n, terrain_hash=export.terrain_hash)

    with Dataset(_CLIM_FILE) as ds:
        land_mask = cv.annual_mean(ds, "lsm") > 0.5
        model_evap = -cv.annual_mean(ds, "evap")
    area_weight_all = gridding.gaussian_area_weights(lat, runoff.shape[1],
                                                     what=str(_CLIM_FILE))
    area_weight_sea = area_weight_all * ~land_mask
    cell_weight = area_weight_all * land_mask
    report = {
        "source_build": build.name,
        "terrain_hash": export.terrain_hash,
        "forcing": _CLIM_NAME,
        "climatology_stage": clim_stage,
        "forcing_sha256": sha256(_CLIM_FILE),
        "orbital_year_days": year_days,
        "runoff_source": {
            "field": "P-E from the forcing climatology, not mrro",
            "why": ("the run's global water budget closes to 0.03% of the mean "
                    "and the land surplus matches the sea deficit, but mrro "
                    "accounts for only 15% of that surplus"),
            "land_mean_p_minus_e_mm_per_day": float(
                np.average(runoff, weights=cell_weight) * SECONDS_PER_DAY * 1000),
            "land_mean_mrro_mm_per_day": float(
                np.average(model_runoff, weights=cell_weight) * SECONDS_PER_DAY * 1000),
        },
        "periodic_steady_state": {
            "time_bins": int(nbin),
            "cycles_to_close": int(periodic["cycles"]),
            "closure_relative_to_capacity": lb.PERIODIC_CLOSURE_RELATIVE,
            "closure_floor_km3": lb.PERIODIC_CLOSURE_FLOOR_KM3,
            "amplitude_publication_threshold": lb.SEASONAL_AMPLITUDE_MIN_FRACTION,
            "amplitude_area_floor_km2": float(area.mean()),
            "basins_closed": int(periodic["closed"].sum()),
            "basins_refused_not_closed": n_refused,
            "basins_refused_no_impoundment": int(periodic["refused_no_impoundment"].sum()),
            "basins_with_published_season": int(periodic["seasonal"].sum()),
            "seasonal_area_km2_at_peak": float(periodic["area_km2"].max(axis=0).sum()),
            "seasonal_area_km2_at_trough": float(periodic["area_km2"].min(axis=0).sum()),
            # THE WATER BALANCE, and it is a measurement rather than a claim.
            # `closed` says the cycle repeats; this says the cycle conserves,
            # which is a different question and the one an integrator can fail
            # silently. Read against what passed through the basins, because a
            # residual is only small relative to something.
            "water_balance": water_balance,
            # THE CROSSING, per region rather than per basin. world-9iy5. The
            # areas above are what the solve holds; these are what a
            # classification can partition, and they are smaller because a
            # lake narrower than a mesh region is not painted.
            "regions_with_a_decided_cycle": int(cycle["decided"].sum()),
            "regions_under_water_in_every_bin": int(permanent_region.sum()),
            "regions_under_water_in_some_bins": int(seasonal_region.sum()),
            "painted_seasonal_area_km2": float(area[seasonal_region].sum()),
            "painted_permanent_area_km2": float(area[permanent_region].sum()),
            "crossing_note": ("closed-basin lake area only. A floodplain's "
                              "inundated area needs a height-above-nearest-"
                              "drainage distribution and a routing model that "
                              "do not exist here, and a seasonally saturated "
                              "soil is a water content rather than an area. "
                              "Neither is merged into this."),
            "what_carries_a_season": ("residence time, and the product is per basin "
                                      "because of it: a deep terminal lake holds "
                                      "years of supply and its surface barely moves "
                                      "inside one, while a shallow playa's area is "
                                      "almost all seasonal"),
            "omission": ("the catchment's delivery is flat through the cycle. Only "
                         "the lake's own surface fluxes carry a season here; "
                         "land_water_ledger.yaml declares that as "
                         "seasonal_phase_of_catchment_delivery"),
        },
        "open_water_evaporation": {
            "method": ("Penman combination, water albedo and roughness, NOT floored "
                       "at the model's land rate; shared with carve_verdict.py"),
            "global_mean_mm_per_day": float(
                np.average(evaporation, weights=area_weight_all) * SECONDS_PER_DAY * 1000),
            "ocean_mean_mm_per_day": float(
                np.average(evaporation, weights=area_weight_sea) * SECONDS_PER_DAY * 1000),
            "model_ocean_mean_mm_per_day": float(
                np.average(model_evap, weights=area_weight_sea) * SECONDS_PER_DAY * 1000),
            "validation": ("ratio over ocean cells, which already are open water; "
                           "the only place the estimate can be checked"),
            # THE SEASONAL JENSEN TERM, and it is the reason the solves below
            # are forced with the bin mean rather than with one evaluation on
            # annual-mean air. Penman is nonlinear in the air it reads, so a
            # single evaluation on a mean over a cycle it varies within is not
            # the mean of the evaluations; that is why it already integrates the
            # diurnal cycle, and the seasonal cycle is the same argument at a
            # longer period. land_water_ledger.yaml holds open_water_evaporation
            # at interval_floor: climatology_bin for it.
            "bin_mean_over_annual_evaluation": jensen,
            "what_the_lake_solves_are_forced_with": (
                "the bin-weighted mean of the per-bin evaluation, over the "
                "basin sinks; the annual-mean evaluation is not used"),
        },
        "evaporation_interval_bracket": interval_bracket,
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
        "staged_background_albedo": staged_albedo,
        "git_commit": subprocess.run(
            ["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True).stdout.strip(),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    _report = ANALYSIS / "surface_water_report.json"
    refuse_a_write_through_a_symlink(
        _report, what="this component's record of the solved lake field",
        instead="Run the generator in the main checkout.")
    _report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
