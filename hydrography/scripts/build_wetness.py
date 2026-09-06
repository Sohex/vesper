#!/usr/bin/env python3
"""One mutually exclusive wetness classification, and its crossing to the grid.

    python hydrography/scripts/build_wetness.py
    python hydrography/scripts/build_wetness.py --selftest

WET-2 and LSHY-6. Worldbuilding: Vesper is an invented super-Earth, and
everything here is about a solved lake surface, a drainage network and a
depression catalogue computed on its native mesh, and about what survives when
that mesh is reduced to the grid the climate model and the biosphere run on.

## Why this artifact exists

One `PEATLAND` fraction collapses states whose carbon and gas physics differ,
and `biosphere/notes/wetlands-peat-methane-audit.md` finding 2 asks instead for
mutually exclusive fractions of persistent peat-forming land, saturated but
non-inundated mineral soil, seasonally flooded land, open lake, river and playa
water, and dry mineral soil. From the other side,
`biosphere/notes/soil-land-surface-hydraulic-consistency-audit.md` finding 7
asks that the crossing to the climate grid stop averaging a dry rootable upland
and a small groundwater-fed playa into one medium-wet surface. Those are the
same artifact: a partition of a cell's land area whose pieces are named.

## Mutual exclusivity is the product, not a property of it

Overlapping wetness classes double-count area, and every ledger downstream
inherits the double count with no way to see it. So a region gets ONE label
rather than a set of masks that happen not to overlap, the shares are required
to partition the cell, and `--selftest` runs fixtures built to violate each of
those in a named way that the checks are required to catch.

## Two supports, which is the cost of a decision already taken

`hydrography/notes/subgrid-water-table.md` section 5: a saturated fraction is a
fraction of a POPULATION, the only population this project holds is the mesh
regions inside a climate-grid cell, and manufacturing a per-region one needs a
sub-region hypsometry GW-6 says is not there. So the resolved classes stay per
REGION and the saturated non-inundated mineral class can only arrive as a
climate-grid AREA share. Exclusivity is then stated at two supports: the
resolved classes partition a cell's land area between them, and the saturated
share is taken OUT of what they leave rather than added beside it. Adding it
beside them is the double count, and it is one of the fixtures.

## What is written, and at which scale each thing lives

  per REGION      the class label, one per region, from the periodic lake
                  cycle, the drainage network's basin membership and the
                  depression catalogue's spill levels. Terrain and the lake
                  solve only: no latitude, no threshold on a cell mean.
  per GRID CELL   the AREA share of each class over the cell's land, which is
                  `lib/gridding.py`'s CATEGORICAL reduction. A class is not a
                  quantity with a mean, so a majority label is not a coarse
                  version of it: a majority discards every minority surface
                  however much leverage it has, and a narrow riparian strip in
                  a dry cell is the surface these rows exist to keep.

  data/<build>/wetness_<grid>.nc
  analysis/wetness_report.json

## What is declared and REFUSED rather than approximated

`config/wetness.yaml` carries each class with the artifact it comes from, and
for the ones with no source it carries the reason. Three of the audit's five
have no source here. The saturated mineral class needed `f_sat`, and the
saturated-area closure is WITHDRAWN in `config/topographic_index.yaml`: the
score that was the only thing able to license it was run and missed, and no
narrower support can resolve the gain it would have had to show. So that class
is absent permanently rather than pending, and this script forms no saturated
share from any input. Peat needs persistence of SATURATION over a cycle, and
the cycle this component carries is one of INUNDATION by a solved lake: the
groundwater solve is a steady state and every term reaching those stores in
`config/land_water_ledger.yaml` has an annual interval floor. The
CLOSED-BASIN third of seasonal inundation IS formed, and the partition is cut
against the periodic cycle rather than the annual equilibrium:
`surface_water.py` solves the lake balance through the climatology's own time
bins and paints that cycle onto regions, so a region wet in every bin is
permanent open water, one wet in some bins is seasonally inundated, and one
never wet inside the depression footprint is playa. The other two thirds of
that quantity are not this component's and are not pending here. A wetness
fraction that is a guess is indistinguishable in the file from one that is a
measurement, so an unavailable class is absent and says why rather than being
estimated.
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
import nc_geometry
from orogen import Export, LAND
from write_door import refuse_a_write_through_a_symlink  # noqa: E402

CFG_PATH = Path(__file__).resolve().parents[1] / "config" / "wetness.yaml"
TI_CFG = Path(__file__).resolve().parents[1] / "config" / "topographic_index.yaml"
SCORE_PATH = ANALYSIS / "topographic_index_score.json"

# The resolved classes, in the order a region is offered to them. ORDER IS NOT
# PRECEDENCE: the checks below refuse a region that two of them claim, so a
# rule that overlaps another is an error rather than something the order
# quietly resolves. The residual is last and is the only class allowed to be
# defined as what is left.
RESOLVED = ("open_water", "seasonal_inundation", "playa", "dry_mineral")


def assign_exclusive(masks: dict, population, residual: str | None = None):
    """One label per region, and a refusal where the masks do not partition.

    THIS IS THE CHECK THE ROW ASKED FOR. `masks` is class name to boolean array;
    every member of `population` must be claimed by exactly one. A region two
    classes claim is an AREA COUNTED TWICE, and it does not announce itself
    downstream -- the shares still look like shares, they simply sum above one,
    and a ledger that renormalises them hides it completely. A region no class
    claims is land that has left the partition.

    `residual` names the one class permitted to absorb the unclaimed remainder.
    It is optional so that a caller which believes its masks are exhaustive can
    say so and be checked on it; passing None and being wrong is the second
    fixture.

    Returns `(labels, names, counts)` with `labels` an int8 index into `names`
    and -1 off the population.
    """
    names = list(masks)
    if residual is not None and residual not in names:
        raise ValueError(f"residual class {residual!r} is not one of {names}")
    population = np.asarray(population, dtype=bool)
    claimed = np.zeros(population.shape, dtype=np.int16)
    for m in masks.values():
        claimed += (np.asarray(m, dtype=bool) & population).astype(np.int16)

    over = population & (claimed > 1)
    if over.any():
        pairs = [n for n in names
                 if (np.asarray(masks[n], dtype=bool) & over).any()]
        raise SystemExit(
            f"{int(over.sum())} regions are claimed by more than one wetness "
            f"class, among {pairs}. That is area counted twice and it does not "
            "announce itself downstream: the shares still look like shares and "
            "simply sum above one. The classes are not exclusive and this "
            "artifact will not be written.")

    under = population & (claimed == 0)
    if under.any():
        if residual is None:
            raise SystemExit(
                f"{int(under.sum())} regions of the population are claimed by "
                "no wetness class and no residual class was named. Land that "
                "has left the partition is not the same as land in a class "
                "called other, and this artifact will not invent the second "
                "from the first.")
        masks = dict(masks)
        masks[residual] = np.asarray(masks[residual], dtype=bool) | under

    labels = np.full(population.shape, -1, dtype=np.int8)
    for i, n in enumerate(names):
        labels[np.asarray(masks[n], dtype=bool) & population] = i
    counts = {n: int((labels == i).sum()) for i, n in enumerate(names)}
    return labels, names, counts


def classify_regions(export: Export, lake, terminal, filled_km, spill_km,
                     cycle_bins_wet, cycle_decided, n_bins):
    """The resolved classes, from the artifacts that already decided them.

    THE CUT IS TAKEN FROM THE PERIODIC CYCLE, NOT FROM THE ANNUAL EQUILIBRIUM,
    and WORLD-T8I5 is the row that measured the two against each other before
    choosing. `surface_water.py` solves the lake balance twice: once as an
    annual equilibrium, which paints `lake`, and once as a periodic steady state
    through the climatology's own time bins, which paints `lake_cycle_bins_wet`.
    The annual area sits between the cycle's trough and its peak, so the two
    paints cut the same ground differently and the strandline is where they
    disagree. Taking the cut from the cycle is what makes `seasonal_inundation`
    a class rather than a remainder against a different solve, which is the
    double count `assign_exclusive` exists to refuse. The measurement that
    licenses it, including how much area the two paints disagree over and where
    that area sits, is `notes/audits/wetness-partition-cut.md`, and this script
    re-measures it on every run and writes it into the report.

    `open_water` is a region under the periodic lake in EVERY time bin: the
    always-wet set, which is the permanent lake surface the cycle expresses.

    `seasonal_inundation` is a region under the lake in some bins and not all.
    IT IS THE CLOSED-BASIN THIRD OF THE SEASONAL QUANTITY AND NOTHING ELSE.
    `hydrography/notes/land-water-ledger.md` splits seasonally inundated land in
    three; a floodplain's inundated area needs a height-above-nearest-drainage
    distribution and a routing model, and a seasonally saturated SOIL is a water
    content rather than an area. Merging either into this class would put three
    quantities under one name.

    `playa` is land inside a closed depression -- its depression-filled surface
    at or below the basin's spill -- that neither of the two wet classes takes:
    the floor and strandline that stay exposed through the whole cycle, which is
    what `pedology/scripts/build_surface_classes.py` already means by the word.
    A basin that fills to its spill in every bin contributes none, which falls
    out of the rule rather than needing a second one. Note that this is the
    depression FOOTPRINT and not the endorheic catchment: the catchment is the
    part of the land that drains inward, and calling that playa would overstate
    the class by more than an order of magnitude.

    A BASIN WHOSE YEAR DID NOT CLOSE HAS NO CYCLE, so its regions carry a
    `bins_wet` of zero that is ambiguous between dry all year and never solved.
    For those regions the annual equilibrium paint is the only statement that
    exists, and it decides: an undecided region the annual solve calls lake is
    `open_water`. The count is written into the report rather than absorbed,
    because it is the one place two solves meet in one partition.

    LAND COMES FROM `surface_class`, never from `land_mask`. CLAUDE.md rule 1:
    the two disagree over dry closed-basin floor below sea level, which is
    exactly the terrain this class is about.
    """
    land = export.surface_class == LAND
    lake = np.asarray(lake, dtype=bool) & land
    terminal = np.asarray(terminal, dtype=np.int64)
    filled_km = np.asarray(filled_km, dtype=np.float64)
    bins_wet = np.asarray(cycle_bins_wet, dtype=np.int64)
    decided = np.asarray(cycle_decided, dtype=bool)
    n_bins = int(n_bins)

    in_basin = land & (terminal >= 0)
    spill = np.zeros(land.shape, dtype=np.float64)
    spill[in_basin] = np.asarray(spill_km, dtype=np.float64)[terminal[in_basin]]
    footprint = in_basin & (filled_km <= spill)

    undecided_lake = land & ~decided & lake
    open_water = (land & decided & (bins_wet >= n_bins)) | undecided_lake
    seasonal = land & decided & (bins_wet > 0) & (bins_wet < n_bins)
    playa = footprint & ~open_water & ~seasonal

    masks = {"open_water": open_water, "seasonal_inundation": seasonal,
             "playa": playa, "dry_mineral": np.zeros(land.shape, dtype=bool)}
    labels, names, counts = assign_exclusive(masks, land, residual="dry_mineral")
    paints = {
        "annual_lake_regions": int(lake.sum()),
        "cycle_always_wet_regions": int((land & decided & (bins_wet >= n_bins)).sum()),
        "cycle_some_bins_wet_regions": int(seasonal.sum()),
        "cycle_undecided_regions": int((land & ~decided).sum()),
        "undecided_and_annual_lake_regions": int(undecided_lake.sum()),
        "annual_lake_not_always_wet_regions": int((lake & ~open_water).sum()),
        "always_wet_not_annual_lake_regions": int((open_water & ~lake).sum()),
        "disagreement_inside_footprint_regions":
            int(((lake ^ open_water) & footprint).sum()),
        "disagreement_outside_footprint_regions":
            int(((lake ^ open_water) & ~footprint).sum()),
    }
    return land, labels, names, counts, (in_basin, footprint), paints, lake


def check_footprint(area_km2, footprint, area_at_spill_km2, tol: float) -> float:
    """The depression footprint against a quantity the catalogue already knows.

    THIS IS A TEST AND NOT A COMPARISON, which is the standard `CLAUDE.md` sets:
    the rule above selects the land a basin would flood at its spill, and
    `basins.nc` carries `area_at_spill_km2` for every basin, computed by
    `build_hydrography.py` from its own hypsometry. Those are the same area
    arrived at two ways, so their totals are an identity with a right answer,
    and a rule that quietly picked the endorheic CATCHMENT instead of the
    depression would fail it by a factor of several -- which is the mistake the
    class invites, because the catchment is the majority of this world's land
    and the footprint is a fifth of it.

    Returns the relative residual; raises above `tol`.
    """
    got = float(np.asarray(area_km2)[np.asarray(footprint, dtype=bool)].sum())
    want = float(np.asarray(area_at_spill_km2, dtype=np.float64).sum())
    rel = abs(got - want) / want if want > 0 else 0.0
    if rel > tol:
        raise SystemExit(
            f"the depression footprint covers {got:,.0f} km2 and the basin "
            f"catalogue's own area at spill totals {want:,.0f} km2, a relative "
            f"difference of {rel:.3g} against a declared {tol:.0e}. They are "
            "the same area arrived at two ways, so this is the playa rule "
            "selecting a different set of regions from the one the hypsometry "
            "was measured over, and the classification is not written on it.")
    return rel


def cross_to_grid(cell, ncell: int, area_km2, land, labels, names,
                  min_land_regions: int):
    """The AREA share of each class per cell, and the account of what was dropped.

    `lib/gridding.py`'s `cell_fraction` is the CATEGORICAL operator and is what
    a class gets: there is no mean of a class, and what survives a reduction is
    the retained fraction of each. The population is the cell's land and is
    passed explicitly, because a fraction is a fraction OF something and a
    reduction that assumes a population is one that can pick the wrong one.

    `covered` is every cell holding any land at all; `estimated` is the subset
    holding enough land regions for the partition to mean anything. A thin cell
    is MARKED and not dropped: its area is real, and a consumer that silently
    loses it has lost land rather than noise.
    """
    shares, areas = {}, {}
    covered = None
    for i, n in enumerate(names):
        f, covered = gridding.cell_fraction(cell, ncell, area_km2,
                                            labels == i, land)
        shares[n] = f
        # THE AREA AS WELL AS THE SHARE, and it comes from the EXTENSIVE
        # operator rather than from share times land area. Two reasons. A
        # consumer converting a depth to a mass needs the DENOMINATOR the
        # ledger names -- lake evaporation is a depth over the open-water area
        # and catchment runoff is a depth over the land cell, and averaging
        # those into a cell mean is the finding this row comes from. And
        # computing it independently makes the two operators check each other:
        # a categorical share and an extensive sum over the same class must
        # agree, and `check_area_closure` requires it.
        areas[n] = gridding.cell_sum(cell, ncell, area_km2, labels == i)
    land_area = gridding.cell_sum(cell, ncell, area_km2, land)
    count = np.bincount(np.asarray(cell)[land], minlength=ncell)
    estimated = covered & (count >= int(min_land_regions))
    ledger = gridding.transfer_ledger(cell, ncell, area_km2, land)
    return shares, areas, land_area, covered, estimated, count.astype(np.int64), ledger


def check_area_closure(shares: dict, areas: dict, land_area, covered,
                       tol: float) -> float:
    """The categorical share and the extensive area must be the same statement.

    `share * land_area` and the class's own summed area are two reductions of
    one population by two operators, so their agreement is an identity with a
    right answer rather than a plausibility check. It is here because the pair
    is exactly where a population mismatch hides: a share taken over the land
    and an area summed over everything would both look ordinary and would
    disagree by the cell's ocean.

    Returns the worst relative departure; raises above `tol`.
    """
    land_area = np.asarray(land_area, dtype=np.float64)
    worst = 0.0
    for n in shares:
        want = np.asarray(shares[n], dtype=np.float64) * land_area
        got = np.asarray(areas[n], dtype=np.float64)
        d = np.abs(got - want)[covered]
        scale = np.maximum(land_area[covered], 1e-30)
        worst = max(worst, float((d / scale).max()) if d.size else 0.0)
    if worst > tol:
        raise SystemExit(
            f"the class shares and the class areas disagree by {worst:.3g} of "
            f"the cell's land area against a declared {tol:.0e}. They are the "
            "same population reduced by two operators, so this is a population "
            "mismatch and not a rounding difference.")
    return worst


def check_partition(shares: dict, covered, tol: float, what: str) -> None:
    """Every share is a share, and together they are all of the cell's land.

    Non-negative, at most one, and summing to one in every covered cell. A
    partition that sums above one is area counted twice; below one is area that
    has left. Both are refused here rather than renormalised, because
    renormalising turns a defect into a plausible number.
    """
    total = np.zeros_like(next(iter(shares.values())), dtype=np.float64)
    for n, f in shares.items():
        f = np.asarray(f, dtype=np.float64)
        bad = covered & ((f < -tol) | (f > 1.0 + tol))
        if bad.any():
            raise SystemExit(
                f"{what}: class {n!r} is not a share in {int(bad.sum())} cells "
                f"(range {f[covered].min():.6g} to {f[covered].max():.6g})")
        total += f
    dev = np.abs(total[covered] - 1.0)
    if dev.size and dev.max() > tol:
        worst = float(dev.max())
        raise SystemExit(
            f"{what}: the class shares do not partition the cell's land. The "
            f"worst cell is off by {worst:.3g}, which is area counted twice if "
            "it is above one and area that has left the partition if it is "
            "below. It is refused rather than renormalised.")


def closure_state() -> dict:
    """The saturated-area closure's standing, read and never decided here.

    `hydrography/config/topographic_index.yaml` carries the decision, because it
    is a statement about what that component publishes. This script reads it and
    records the score beside it as the evidence, so the artifact says why the
    saturated class is absent without a reader having to follow two files.
    """
    ti = yaml.safe_load(TI_CFG.read_text(encoding="utf-8"))["closure"]
    status = str(ti.get("status", "active"))
    # `withdrawn_on` is a YAML 1.1 DATE, which json refuses. Rendered here
    # rather than at the point of writing: the report and the netCDF attribute
    # both read this dict, and a value that is a date in one and a string in the
    # other is the same defect twice.
    withdrawn_on = ti.get("withdrawn_on")
    out = {"closure_status": status,
           "withdrawn": status == "withdrawn",
           "withdrawn_on": None if withdrawn_on is None else str(withdrawn_on),
           "withdrawn_because": ti.get("withdrawn_because"),
           "reason": "hydrography/config/topographic_index.yaml closure.status"}
    if SCORE_PATH.exists():
        sc = json.loads(SCORE_PATH.read_text())
        out["score"] = {
            "observation_sets_required": sc.get("observation_sets_required"),
            "per_set": {k: v.get("passes") for k, v in sc.get("runs", {}).items()}}
    return out


def read_config() -> dict:
    cfg = yaml.safe_load(CFG_PATH.read_text(encoding="utf-8"))
    r = cfg["refusals"]
    # THE REFUSALS ARE ENFORCED, not documented. Each of these is a route that
    # returns an ordinary-looking number and a wrong one, and each was named as
    # a defect somewhere before it was refused here.
    if str(r["regime_selector"]) != "none":
        raise SystemExit(
            f"refusals.regime_selector is {r['regime_selector']!r}. No latitude "
            "and no other geographic selector may pick a wetness class: the "
            "vendored decomposition already switches on one, and a "
            "classification that switched too would put the physics in two "
            "places. This script has no such selector to enable.")
    if str(r["area_proxy_from_cell_mean"]) != "refused":
        raise SystemExit(
            "refusals.area_proxy_from_cell_mean is not 'refused'. A cell-mean "
            "water table depth or available water capacity may not stand in "
            "for an area fraction; runoff generation, groundwater access, lake "
            "evaporation and root access happen on different fractions and no "
            "choice of mean fixes it. The quantity here is an AREA.")
    if str(r["saturated_fraction_support"]) != "climate_grid":
        raise SystemExit(
            f"refusals.saturated_fraction_support is "
            f"{r['saturated_fraction_support']!r}. A saturated fraction is a "
            "fraction of a population and the only population this project "
            "holds is the mesh regions inside a climate-grid cell. There is no "
            "per-region one and manufacturing one needs a sub-region "
            "hypsometry GW-6 says is not there.")
    if str(r["absolute_index_threshold"]) != "forbidden":
        raise SystemExit(
            "refusals.absolute_index_threshold is not 'forbidden'. The "
            "topographic index carries a length, so its whole distribution "
            "shifts with the mesh and only its rank statistics transport.")
    return cfg


def _terrain_check(path: Path, export: Export, regions_expected: int) -> None:
    with Dataset(path) as ds:
        th = getattr(ds, "terrain_hash", None)
        if th is not None and th != export.terrain_hash:
            raise SystemExit(
                f"{path} was built on terrain {th[:16]} and the mesh is "
                f"{export.terrain_hash[:16]}; they are not the same world and "
                "no mapping between them is defined")
        if "region" in ds.dimensions:
            n = ds.dimensions["region"].size
            if n != regions_expected:
                raise SystemExit(
                    f"{path} carries {n} regions and the mesh has "
                    f"{regions_expected}; not the same support")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="WET-2 and LSHY-6: one mutually exclusive wetness "
                    "classification on the native mesh, and its area shares on "
                    "the climate grid.")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--resolution", default=None,
                    help="grid to cross onto; defaults to config's resolution")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="run the exclusivity and partition checks against "
                         "fixtures built to violate them; needs no build")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()

    cfg = read_config()
    tol = float(cfg["closure"]["partition_tolerance"])
    min_regions = int(cfg["crossing"]["min_land_regions_per_cell"])
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))

    export = Export(builds.mesh_export(config))
    resolution = (args.resolution or str(config["model"]["resolution"])).upper()
    grid_dir = builds.grid_export(config, resolution)
    grid_name = grid_dir.name
    data = builds.component_data("hydrography", config, strict=True)

    for name in ("surface_water.nc", "regions.nc", "basins.nc"):
        if not (data / name).is_file():
            raise SystemExit(
                f"{data / name} does not exist. Run "
                "hydrography/scripts/build_hydrography.py and "
                "hydrography/scripts/surface_water.py first; this step "
                "classifies their output and cannot stand in for either.")
    _terrain_check(data / "surface_water.nc", export, export.n_regions)
    _terrain_check(data / "basins.nc", export, export.n_regions)

    with Dataset(data / "surface_water.nc") as ds:
        lake = np.asarray(ds["lake"][:]).astype(bool)
        cycle_bins_wet = np.asarray(ds["lake_cycle_bins_wet"][:], dtype=np.int64)
        cycle_decided = np.asarray(ds["lake_cycle_decided"][:]).astype(bool)
        n_bins = ds.dimensions["time_bin"].size
    with Dataset(data / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:], dtype=np.int64)
        filled_km = np.asarray(ds["filled_km"][:], dtype=np.float64)
    with Dataset(data / "basins.nc") as ds:
        spill_km = np.asarray(ds["spill_km"][:], dtype=np.float64)

    print(f"mesh {export.n_regions:,} regions from {export.root.name}, "
          f"crossing onto {grid_name}")
    land, labels, names, counts, (in_basin, footprint), paints, annual_lake = (
        classify_regions(export, lake, terminal, filled_km, spill_km,
                         cycle_bins_wet, cycle_decided, n_bins))
    area_km2 = export.cell_area.astype(np.float64)
    with Dataset(data / "basins.nc") as ds:
        area_at_spill = np.asarray(ds["area_at_spill_km2"][:], dtype=np.float64)
    footprint_residual = check_footprint(area_km2, footprint, area_at_spill,
                                         float(cfg["closure"]["footprint_tolerance"]))
    land_area = float(area_km2[land].sum())
    print(f"  depression footprint {float(area_km2[footprint].sum()) / land_area:6.2%} "
          f"of land, against the catalogue's own area at spill to "
          f"{footprint_residual:.1e} relative")
    print(f"  the endorheic CATCHMENT is "
          f"{float(area_km2[in_basin].sum()) / land_area:6.2%} of land and is a "
          f"different quantity; the class is the floor, not the catchment")
    for i, n in enumerate(names):
        a = float(area_km2[labels == i].sum())
        print(f"  {n:19s} {counts[n]:>10,} regions   "
              f"{a / area_km2[land].sum():6.2%} of land area")

    # THE TWO PAINTS, RE-MEASURED ON EVERY RUN RATHER THAN CITED. WORLD-T8I5
    # chose the cycle cut on this measurement, and a measurement that only ever
    # ran once is a number the artifact would carry without checking. The
    # quantity that licenses the choice is where the disagreement SITS: the
    # annual equilibrium area lies between the cycle's trough and its peak, so
    # the two paints may differ only in the strandline, and strandline is
    # inside the depression footprint by construction. Area outside it is the
    # two solves disagreeing about something else.
    open_water = labels == names.index("open_water")
    disagreement = (annual_lake ^ open_water) & land
    paints["open_water_land_area_fraction"] = float(
        area_km2[open_water].sum() / land_area)
    paints["annual_lake_land_area_fraction"] = float(
        area_km2[annual_lake].sum() / land_area)
    paints["disagreement_land_area_fraction"] = float(
        area_km2[disagreement].sum() / land_area)
    paints["disagreement_share_of_annual_lake"] = float(
        area_km2[disagreement].sum()
        / max(float(area_km2[annual_lake].sum()), 1e-30))
    paints["disagreement_outside_footprint_land_area_fraction"] = float(
        area_km2[disagreement & ~footprint].sum() / land_area)
    print(f"  the two lake paints disagree over "
          f"{paints['disagreement_land_area_fraction']:.4%} of land, "
          f"{paints['disagreement_share_of_annual_lake']:.3f} of the annual "
          f"lake area, of which "
          f"{paints['disagreement_outside_footprint_land_area_fraction']:.4%} "
          f"of land lies outside the depression footprint")
    if paints["undecided_and_annual_lake_regions"]:
        print(f"  {paints['undecided_and_annual_lake_regions']:,} regions whose "
              "basin did not close its year are open water on the annual paint "
              "alone")

    cell, nlat, nlon = gridding.region_cells(export, grid_dir)
    ncell = nlat * nlon
    shares, areas, land_area_cell, covered, estimated, count, ledger = cross_to_grid(
        cell, ncell, area_km2, land, labels, names, min_regions)
    check_partition(shares, covered, tol, "the resolved classes")
    area_residual = check_area_closure(shares, areas, land_area_cell, covered, tol)
    print(f"  {int(covered.sum()):,} of {ncell:,} cells hold land, "
          f"{int(estimated.sum()):,} hold at least {min_regions} land regions")

    lic = closure_state()
    unresolved = {k: v["unavailable_reason"]
                  for k, v in cfg["classes"].items()
                  if "unavailable_reason" in v}
    if not lic["withdrawn"]:
        raise SystemExit(
            "the saturated-area closure's status in "
            "hydrography/config/topographic_index.yaml is "
            f"{lic['closure_status']!r} and this script has no route that forms "
            "a saturated class. It was withdrawn there on the declared score's "
            "verdict; reviving it is a new closure with a new criterion, and "
            "that criterion decides what a bracketed class means here before "
            "any share is written. A class formed by the first plausible "
            "reading of a revived closure is worse than an absent one.")
    unresolved["saturated_mineral"] = (
        cfg["classes"]["saturated_mineral"]["unavailable_reason"]
        + f" Read here from {lic['reason']}: {lic['closure_status']} on "
        + f"{lic['withdrawn_on']}, {lic['withdrawn_because']}.")
    print("  classes with no source, declared and left absent:")
    for k in unresolved:
        print(f"    {k}")

    out = args.output or (data / f"wetness_{grid_name}.nc")
    refuse_a_write_through_a_symlink(
        out, what="the wetness classification a consumer of this build reads",
        instead=("Pass --output to a path inside this worktree, or run the "
                 "generator in the main checkout."))
    with Dataset(out, "w", format="NETCDF4") as ds:
        ds.title = "Mutually exclusive wetness classification and its area shares"
        ds.summary = (
            "WET-2 and LSHY-6. Per region one wetness class, from the solved "
            "lake surface, basin membership and the depression catalogue's "
            "spill levels; per grid cell the AREA share of each class over the "
            "cell's land, which is lib/gridding.py's CATEGORICAL reduction. "
            "The shares partition the cell's land exactly and the partition is "
            "checked rather than assumed, because overlapping wetness classes "
            "double-count area and every ledger downstream inherits it.")
        ds.classes_absent = json.dumps(unresolved)
        ds.exclusivity = (
            "One label per region, refused where two classes claim one or "
            "where none does. The saturated non-inundated mineral class is not "
            "here and is not pending: the saturated-area closure it would have "
            "come from is withdrawn in hydrography/config/topographic_index."
            "yaml on its declared score's verdict, so dry_mineral is the whole "
            "of the land the resolved classes leave. hydrography/notes/"
            "subgrid-water-table.md sections 5 to 7 carry the decision.")
        ds.saturated_area_closure = lic["closure_status"]
        ds.source_build = export.root.parent.name
        ds.terrain_hash = export.terrain_hash
        ds.grid = grid_name
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

        wc = var("wetness_class", labels, "i1", ("region",), "1",
                 "one wetness class per region, -1 off the land population. "
                 "Land is surface_class == LAND and never land_mask")
        wc.flag_values = np.arange(len(names), dtype=np.int8)
        wc.flag_meanings = " ".join(names)
        var("cell_index", cell, "i4", ("region",), "1",
            "flat grid-cell index this region falls in, row * nlon + col, from "
            "lib/gridding.py. CLAUDE.md rule 3: index, never longitude")

        shape = (nlat, nlon)
        for n in names:
            var(f"f_{n}", shares[n].reshape(shape), "f4", ("lat", "lon"), "1",
                f"share of the cell's land AREA in class {n}, "
                "lib/gridding.py's CATEGORICAL reduction. The shares partition "
                "the cell's land and are checked to; a majority label is not a "
                "coarse version of them")
        var("land_area_km2", land_area_cell.reshape(shape), "f8",
            ("lat", "lon"), "km2",
            "the cell's land area, lib/gridding.py's EXTENSIVE reduction. The "
            "denominator a share above is a share OF; a depth converted to a "
            "mass on the wrong one of these is finding 7 of the "
            "hydraulic-consistency audit")
        for n in names:
            var(f"area_{n}_km2", areas[n].reshape(shape), "f8",
                ("lat", "lon"), "km2",
                f"area of class {n} in the cell, lib/gridding.py's EXTENSIVE "
                "reduction. Emitted beside the share because a lake "
                "evaporation is a depth over the open-water area and a "
                "catchment runoff is a depth over the land cell, and the two "
                "denominators are different")
        var("land_regions", count.reshape(shape), "i4", ("lat", "lon"), "1",
            "land regions in the cell; below the configured minimum the shares "
            "are still area and are marked rather than dropped")
        var("land_partition_estimated", estimated.reshape(shape).astype(np.int8),
            "i1", ("lat", "lon"), "1",
            "1 where the cell holds enough land regions for the partition to "
            "mean anything, 0 where it holds land but too little of it")

        # LAST in the block. The axis is constructed on the export's centres,
        # the convention gridding.region_cells binned these cells on, and is
        # named as such so it cannot be matched against a model-labelled one.
        nc_geometry.declare_grid(ds, convention=nc_geometry.EXPORT_CENTRES,
                                 what="the wetness classification")

    report = {
        "script": "hydrography/scripts/build_wetness.py",
        "issues": ["WET-2", "LSHY-6"],
        "source_build": export.root.parent.name,
        "grid": grid_name,
        "regions": int(export.n_regions),
        "land_regions": int(land.sum()),
        "depression_footprint_land_fraction": float(area_km2[footprint].sum() / land_area),
        "endorheic_catchment_land_fraction": float(area_km2[in_basin].sum() / land_area),
        "footprint_against_catalogue_area_at_spill_relative": footprint_residual,
        "classes_resolved": {
            n: {"regions": counts[n],
                "land_area_fraction": float(area_km2[labels == i].sum()
                                            / area_km2[land].sum())}
            for i, n in enumerate(names)},
        "classes_absent": unresolved,
        "two_lake_paints": paints,
        "f_sat_license": lic,
        "reduction": {
            "class_share": "gridding.cell_fraction, CATEGORICAL, share of the cell's land AREA",
            "class_area_km2": "gridding.cell_sum, EXTENSIVE, the area itself",
            "land_area_km2": "gridding.cell_sum, EXTENSIVE",
            "population": "surface_class == LAND",
            "share_against_area_worst_relative": area_residual,
            "ledger": ledger,
        },
        "cells_with_land": int(covered.sum()),
        "cells_with_partition_estimated": int(estimated.sum()),
        "min_land_regions_per_cell": min_regions,
        "class_share_percentiles": {
            n: {str(p): float(np.percentile(shares[n][estimated], p))
                for p in (5, 50, 95)} for n in names},
        "partition_tolerance": tol,
        "created": datetime.now(timezone.utc).isoformat(),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    rp = ANALYSIS / "wetness_report.json"
    refuse_a_write_through_a_symlink(
        rp, what="this component's record of the wetness classification",
        instead="Run the generator in the main checkout.")
    rp.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {out}\nwrote {rp}")
    return 0


def _selftest() -> int:
    """Nine checks, and four of them are fixtures built to be WRONG in a named way.

    A check that has only ever seen correct input has not been tested, so every
    rule this script enforces is handed a case that violates it and is required
    to refuse. The negative control at the end is the one the artifact exists
    for: a cell that is almost all dry upland with a small lake in it must
    return the lake's area share, which a majority label returns as zero.
    """
    problems: list[str] = []
    n_checks = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal n_checks
        n_checks += 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    def refuses(fn) -> tuple[bool, str]:
        try:
            fn()
        except SystemExit as exc:
            return True, str(exc)
        return False, "accepted"

    # The config the rest of this script is nothing without. It went
    # unparseable once and stayed that way, because --selftest returned before
    # it was ever read and no other check loads it: a class description gained
    # a bare "three: " inside a plain scalar and every run of the real script
    # died on the first line. So the first check is that the file loads and
    # that its refusals still say what the enforcement below expects.
    try:
        _cfg = read_config()
        cfg_ok, cfg_detail = True, ""
    except Exception as exc:                                # noqa: BLE001
        _cfg, cfg_ok, cfg_detail = None, False, f"{type(exc).__name__}: {exc}"
    check("config/wetness.yaml parses and its refusals hold", cfg_ok, cfg_detail)
    if _cfg is not None:
        # Every class either names a source or says why it has none. A class
        # with neither is a share a consumer would form from nothing.
        stated = [k for k, v in _cfg["classes"].items()
                  if not (v.get("source") or v.get("unavailable_reason"))]
        check("every class names a source or says why it has none",
              not stated, f"neither on {stated}")

    rng = np.random.default_rng(20260825)
    n, ncell = 4000, 25
    cell = rng.integers(0, ncell, n)
    area = rng.uniform(0.5, 5.0, n)
    land = rng.random(n) < 0.7
    lake = land & (rng.random(n) < 0.12)
    playa = land & ~lake & (rng.random(n) < 0.08)

    # 1. A well-formed case partitions the cell exactly.
    labels, names, _ = assign_exclusive(
        {"open_water": lake, "playa": playa,
         "dry_mineral": np.zeros(n, bool)}, land, residual="dry_mineral")
    shares, areas, land_area, covered, _, _, ledger = cross_to_grid(
        cell, ncell, area, land, labels, names, 1)
    ok = True
    try:
        check_partition(shares, covered, 1e-12, "fixture")
    except SystemExit as exc:
        ok, why = False, str(exc)
    check("well-formed classes partition the cell's land", ok,
          "" if ok else why)

    # 2. FIXTURE, overlapping classes. A playa rule that forgot to exclude the
    #    lake is the real version of this bug: every lake region inside a
    #    non-filling basin's footprint would be claimed twice.
    caught, msg = refuses(lambda: assign_exclusive(
        {"open_water": lake, "playa": playa | lake,
         "dry_mineral": np.zeros(n, bool)}, land, residual="dry_mineral"))
    check("a region claimed by two classes is refused", caught, msg)

    # 3. FIXTURE, land in no class and no residual offered.
    caught, msg = refuses(lambda: assign_exclusive(
        {"open_water": lake, "playa": playa}, land, residual=None))
    check("land claimed by no class is refused when no residual is named",
          caught, msg)

    # 4. FIXTURE, an unresolved climate-grid share ADDED beside the resolved
    #    classes instead of taken out of them. That is the double count this
    #    artifact exists to prevent, and the partition check is what has to
    #    catch it: the shares still look like shares and simply sum above one.
    #    No such class is written -- the one that would have been, the saturated
    #    mineral share, comes from a closure that is withdrawn -- so this is a
    #    standing guard on whatever unresolved class arrives next rather than a
    #    test of a live path.
    added = dict(shares)
    added["an_unresolved_grid_share"] = np.full(ncell, 0.2) * covered
    caught, msg = refuses(
        lambda: check_partition(added, covered, 1e-12, "fixture"))
    check("an unresolved grid share added beside the resolved classes is refused",
          caught, msg)

    # 5. THE CUT IS TAKEN FROM THE CYCLE, and this is the fixture that says so
    #    with a case that can fail. A region the ANNUAL paint calls lake, whose
    #    basin's cycle leaves it dry in some bins, must be seasonal_inundation
    #    and not open_water; a region wet in every bin must be open_water; and a
    #    region whose basin never closed its year has no cycle, so the annual
    #    paint is the only statement about it and it decides. Getting any of the
    #    three wrong returns an ordinary-looking partition, which is why each is
    #    asserted by name rather than by the totals agreeing.
    class _FakeExport:
        def __init__(self, surface_class):
            self.surface_class = surface_class

    n5 = 6
    sc = np.full(n5, LAND, dtype=np.int8)
    fake = _FakeExport(sc)
    #  region: 0 always wet, 1 wet in some bins, 2 dry floor, 3 upland,
    #          4 undecided and painted lake by the annual solve,
    #          5 undecided and not painted
    terminal5 = np.array([0, 0, 0, -1, 1, 1], dtype=np.int64)
    filled5 = np.array([1.0, 1.0, 1.0, 9.0, 1.0, 1.0])
    spill5 = np.array([2.0, 2.0])
    lake5 = np.array([1, 1, 0, 0, 1, 0], dtype=bool)
    bins5 = np.array([12, 5, 0, 0, 0, 0], dtype=np.int64)
    dec5 = np.array([1, 1, 1, 0, 0, 0], dtype=bool)
    _land, lab5, nm5, cnt5, (_ib, fp5), paints5, _al = classify_regions(
        fake, lake5, terminal5, filled5, spill5, bins5, dec5, 12)
    got = {nm5[int(lab5[i])] for i in range(n5)}
    check("the cycle decides open water and the annual paint does not",
          nm5[int(lab5[1])] == "seasonal_inundation"
          and nm5[int(lab5[0])] == "open_water",
          f"region 0 is {nm5[int(lab5[0])]!r}, region 1 is "
          f"{nm5[int(lab5[1])]!r} while the annual paint calls both lake")
    check("a region whose basin never closed its year takes the annual paint",
          nm5[int(lab5[4])] == "open_water"
          and paints5["undecided_and_annual_lake_regions"] == 1,
          f"region 4 is {nm5[int(lab5[4])]!r} and the report counts "
          f"{paints5['undecided_and_annual_lake_regions']}")
    check("a seasonally inundated region inside the footprint is not playa",
          nm5[int(lab5[2])] == "playa" and "playa" not in
          {nm5[int(lab5[0])], nm5[int(lab5[1])]},
          f"labels {[nm5[int(x)] for x in lab5]}")
    check("land outside every depression is the residual",
          nm5[int(lab5[3])] == "dry_mineral" and sum(cnt5.values()) == n5,
          f"region 3 is {nm5[int(lab5[3])]!r}, counts {cnt5}")

    # 6. THE NEGATIVE CONTROL. One cell, 97% dry upland and 3% lake. The share
    #    must be the lake's area, which a majority label returns as zero.
    c2 = np.zeros(1000, dtype=np.int64)
    a2 = np.ones(1000)
    l2 = np.ones(1000, dtype=bool)
    lk = np.zeros(1000, dtype=bool)
    lk[:30] = True
    lab2, nm2, _ = assign_exclusive(
        {"open_water": lk, "playa": np.zeros(1000, bool),
         "dry_mineral": np.zeros(1000, bool)}, l2, residual="dry_mineral")
    s2, _, _, cov2, _, _, _ = cross_to_grid(c2, 1, a2, l2, lab2, nm2, 1)
    got = float(s2["open_water"][0])
    check("a small lake survives a cell that is almost all upland",
          abs(got - 0.03) <= 1e-12, f"open-water share {got}")

    # 7. FIXTURE, the playa rule selecting the endorheic CATCHMENT instead of
    #    the depression floor. That is the mistake this class invites, because
    #    on this world the catchment is most of the land and the floor is a
    #    fifth of it, and the identity against the catalogue's own area at spill
    #    is what catches it.
    fa = np.ones(500)
    floor = np.zeros(500, dtype=bool)
    floor[:100] = True
    catchment = np.zeros(500, dtype=bool)
    catchment[:400] = True
    got = check_footprint(fa, floor, np.array([60.0, 40.0]), 1e-9)
    caught, msg = refuses(
        lambda: check_footprint(fa, catchment, np.array([60.0, 40.0]), 1e-9))
    check("the floor matches the catalogue's area at spill and the catchment "
          "is refused",
          got == 0.0 and caught, msg if not caught else f"residual {got}")

    # 8. The categorical share and the extensive area are the same statement,
    #    and a share taken over the land against an area summed over everything
    #    -- the population mismatch, which is what actually goes wrong -- is
    #    refused.
    got = check_area_closure(shares, areas, land_area, covered, 1e-12)
    wrong = {n: gridding.cell_sum(cell, ncell, area, labels == i)
             for i, n in enumerate(names)}
    wrong["dry_mineral"] = wrong["dry_mineral"] + gridding.cell_sum(
        cell, ncell, area, ~land)
    caught, msg = refuses(
        lambda: check_area_closure(shares, wrong, land_area, covered, 1e-12))
    check("the class share and the class area agree, and a population mismatch "
          "is refused",
          got <= 1e-12 and caught, msg if not caught else f"residual {got}")

    print(f"\n{n_checks} checks, {len(problems)} failed")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
