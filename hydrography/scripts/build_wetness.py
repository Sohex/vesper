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

  per REGION      the class label, one per region, from the solved lake
                  surface, the drainage network's basin membership and the
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
have no source here. The saturated mineral class needs `f_sat`, which is
unlicensed until the score in `config/topographic_index.yaml` passes, and which
needs a water table as well. Seasonal inundation and peat need a season, and
nothing in this component carries one: the lake solve is an annual equilibrium,
the groundwater solve is a steady state, and every term reaching those stores in
`config/land_water_ledger.yaml` has an annual interval floor. A wetness fraction
that is a guess is indistinguishable in the file from one that is a measurement,
so an unavailable class is absent and says why rather than being estimated.
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

CFG_PATH = Path(__file__).resolve().parents[1] / "config" / "wetness.yaml"
SCORE_PATH = ANALYSIS / "topographic_index_score.json"

# The resolved classes, in the order a region is offered to them. ORDER IS NOT
# PRECEDENCE: the checks below refuse a region that two of them claim, so a
# rule that overlaps another is an error rather than something the order
# quietly resolves. The residual is last and is the only class allowed to be
# defined as what is left.
RESOLVED = ("open_water", "playa", "dry_mineral")


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


def classify_regions(export: Export, lake, terminal, filled_km, spill_km):
    """The resolved classes, from the artifacts that already decided them.

    `open_water` is `surface_water.py`'s solved equilibrium lake surface,
    painted by area up each basin's own hypsometry. It is the only per-region
    open-water statement in the project and it is a result rather than a level
    threshold.

    `playa` is land inside a closed depression -- its depression-filled surface
    at or below the basin's spill -- that the solved lake does not cover: the
    exposed floor and strandline of a basin that does not fill, which is what
    `pedology/scripts/build_surface_classes.py` already means by the word. A
    basin that DOES fill to its spill contributes none, and that falls out of
    the rule rather than needing a second one, because its lake covers the
    depression by definition. Note that this is the depression FOOTPRINT and
    not the endorheic catchment: the catchment is the fifth of the land that
    drains inward, and calling that playa would overstate the class by more
    than an order of magnitude.

    LAND COMES FROM `surface_class`, never from `land_mask`. CLAUDE.md rule 1:
    the two disagree over dry closed-basin floor below sea level, which is
    exactly the terrain this class is about.
    """
    land = export.surface_class == LAND
    lake = np.asarray(lake, dtype=bool) & land
    terminal = np.asarray(terminal, dtype=np.int64)
    filled_km = np.asarray(filled_km, dtype=np.float64)

    in_basin = land & (terminal >= 0)
    spill = np.zeros(land.shape, dtype=np.float64)
    spill[in_basin] = np.asarray(spill_km, dtype=np.float64)[terminal[in_basin]]
    playa = in_basin & ~lake & (filled_km <= spill)

    masks = {"open_water": lake, "playa": playa,
             "dry_mineral": np.zeros(land.shape, dtype=bool)}
    labels, names, counts = assign_exclusive(masks, land, residual="dry_mineral")
    return land, labels, names, counts, (in_basin, in_basin & (filled_km <= spill))


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


def take_saturated_out(shares: dict, f_sat, covered, tol: float):
    """The unresolved saturated share, taken OUT of the mineral class.

    `f_sat` is a share of the CELL, from a rank statistic over that cell's whole
    land area, and the resolved classes have already claimed part of that area.
    So the saturated non-inundated mineral class is `min(f_sat, dry_mineral)`
    and the mineral class keeps the rest. ADDING IT BESIDE THE RESOLVED CLASSES
    IS THE DOUBLE COUNT this artifact exists to prevent, and it is a fixture.

    THE CLAMP IS A DECISION AND IS REPORTED RATHER THAN ABSORBED. Where `f_sat`
    exceeds what the resolved classes leave, the closure has predicted saturated
    area that is already realised as open water or exposed basin floor -- the
    wettest part of the cell, which is where a lake is. The resolved classes
    win, because they are a solved result and the closure is a statistic; the
    share of cells where that binds, and the area it removes, go in the report,
    because a clamp that binds everywhere is a disagreement and not a detail.
    """
    out = dict(shares)
    mineral = np.asarray(shares["dry_mineral"], dtype=np.float64)
    f_sat = np.asarray(f_sat, dtype=np.float64)
    sat = np.minimum(f_sat, mineral)
    out["saturated_mineral"] = sat
    out["dry_mineral"] = mineral - sat
    binds = covered & (f_sat > mineral + tol)
    return out, binds


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


def license_state() -> dict:
    """Whether `f_sat` may be taken, read from the score and never decided here."""
    if not SCORE_PATH.exists():
        return {"consumers_licensed": False,
                "reason": f"no score at {SCORE_PATH.name}"}
    sc = json.loads(SCORE_PATH.read_text())
    return {"consumers_licensed": bool(sc.get("consumers_licensed", False)),
            "observation_sets_required": sc.get("observation_sets_required"),
            "per_set": {k: v.get("passes") for k, v in
                        sc.get("runs", {}).items()},
            "reason": "the declared score in "
                      "hydrography/config/topographic_index.yaml"}


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
    with Dataset(data / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:], dtype=np.int64)
        filled_km = np.asarray(ds["filled_km"][:], dtype=np.float64)
    with Dataset(data / "basins.nc") as ds:
        spill_km = np.asarray(ds["spill_km"][:], dtype=np.float64)

    print(f"mesh {export.n_regions:,} regions from {export.root.name}, "
          f"crossing onto {grid_name}")
    land, labels, names, counts, (in_basin, footprint) = classify_regions(
        export, lake, terminal, filled_km, spill_km)
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
        print(f"  {n:14s} {counts[n]:>10,} regions   "
              f"{a / area_km2[land].sum():6.2%} of land area")

    cell, nlat, nlon = gridding.region_cells(export, grid_dir)
    ncell = nlat * nlon
    shares, areas, land_area_cell, covered, estimated, count, ledger = cross_to_grid(
        cell, ncell, area_km2, land, labels, names, min_regions)
    check_partition(shares, covered, tol, "the resolved classes")
    area_residual = check_area_closure(shares, areas, land_area_cell, covered, tol)
    print(f"  {int(covered.sum()):,} of {ncell:,} cells hold land, "
          f"{int(estimated.sum()):,} hold at least {min_regions} land regions")

    lic = license_state()
    unresolved = {k: v["unavailable_reason"]
                  for k, v in cfg["classes"].items()
                  if "unavailable_reason" in v}
    ti = data / f"topographic_index_{grid_name}.nc"
    wt = data / "water_table.nc"
    saturated = None
    if lic["consumers_licensed"] and ti.is_file() and wt.is_file():
        raise SystemExit(
            "f_sat is licensed and both inputs are present, and this script "
            "does not yet form the saturated class from them. take_saturated_"
            "out() is the rule and is under --selftest; wiring it needs the "
            "f_grad arm to be carried along a member axis the way "
            "build_groundwater_access.py carries the permeability bracket, "
            "which is its own row. Refusing is deliberate: a class formed by "
            "the first plausible reading of a bracketed closure is worse than "
            "an absent one.")
    reason = ("f_sat is not licensed: " + str(lic.get("reason"))
              if not lic["consumers_licensed"] else
              "f_sat is licensed but " + (
                  f"{ti.name} is absent" if not ti.is_file()
                  else f"{wt.name} is absent"))
    unresolved["saturated_mineral"] = (
        cfg["classes"]["saturated_mineral"]["unavailable_reason"]
        + " Measured here: " + reason)
    print("  classes with no source, declared and left absent:")
    for k in unresolved:
        print(f"    {k}")

    out = args.output or (data / f"wetness_{grid_name}.nc")
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
            "where none does. The saturated non-inundated mineral class is "
            "unresolved on this mesh and can only arrive as a climate-grid "
            "area share taken OUT of dry_mineral, never added beside the "
            "resolved classes; hydrography/notes/subgrid-water-table.md "
            "section 5 is the decision and this file carries its cost.")
        ds.consumers_licensed_f_sat = "yes" if lic["consumers_licensed"] else "no"
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
    rp.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {out}\nwrote {rp}")
    return 0


def _selftest() -> int:
    """Six checks, and four of them are fixtures built to be WRONG in a named way.

    A check that has only ever seen correct input has not been tested, so every
    rule this script enforces is handed a case that violates it and is required
    to refuse. The negative control at the end is the one the artifact exists
    for: a cell that is almost all dry upland with a small lake in it must
    return the lake's area share, which a majority label returns as zero.
    """
    problems: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    def refuses(fn) -> tuple[bool, str]:
        try:
            fn()
        except SystemExit as exc:
            return True, str(exc)
        return False, "accepted"

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

    # 4. FIXTURE, the saturated share ADDED beside the resolved classes instead
    #    of taken out of them. This is the double count the row is about, and
    #    the partition check is what has to catch it: the shares still look like
    #    shares and simply sum above one.
    f_sat = np.full(ncell, 0.2)
    added = dict(shares)
    added["saturated_mineral"] = f_sat * covered
    caught, msg = refuses(
        lambda: check_partition(added, covered, 1e-12, "fixture"))
    check("a saturated share added beside the resolved classes is refused",
          caught, msg)

    # 5. Taken OUT of the mineral class instead, the partition still closes and
    #    the clamp is reported where the closure asks for more area than the
    #    resolved classes leave.
    taken, binds = take_saturated_out(shares, np.full(ncell, 0.9), covered, 1e-12)
    ok = True
    try:
        check_partition(taken, covered, 1e-12, "fixture")
    except SystemExit as exc:
        ok, why = False, str(exc)
    check("taken out of the mineral class, the partition still closes and the "
          "clamp is reported",
          ok and bool(binds[covered].any()),
          "" if ok else why)

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

    print(f"\n8 checks, {len(problems)} failed")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
