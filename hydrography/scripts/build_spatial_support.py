#!/usr/bin/env python3
"""The 10M native mesh under every ladder cell: area, hypsometry and class shares.

    python hydrography/scripts/build_spatial_support.py                    # configured rung
    python hydrography/scripts/build_spatial_support.py --rungs T21 T42 T85 T127 T170
    python hydrography/scripts/build_spatial_support.py --selftest         # needs no build

SPAT-3, which extends GRID-2. Worldbuilding: Vesper is an invented super-Earth.
Everything below is a statement about a terrain a generator made on a 10M-region
mesh, and about what survives when that mesh is reduced to a Gaussian grid the
climate model and the biosphere run on. Nothing is run: no climate model, no
ocean, no biosphere.

## Why this artifact exists

GRID-2: the distribution of MESH elevations inside a model grid cell is computed
three times in this tree and persisted nowhere, so three implementations of one
criterion sit on a field no artifact carries. It is not a small quantity. A cell
is evaluated at its MEAN elevation and its high ground stands well above that,
so the model runs its own peaks warm by the peak excess times the lapse rate and
the ice that would sit on them never forms.

The spatial-support audit's finding 2 is the same absence from the other side:
several consumers independently rebuild counts, moments or class shares from the
mesh, and the support artifact they should all read does not exist. This is that
artifact.

## What it is NOT, and each of these is a row's own reservation

  NOT a second land or bathymetry product. The land MASK the model reads is
  `exoplasim/scripts/build_boundary_conditions.py`'s and the ocean bathymetry
  contract is OCN-11's. What is here is the SUPPORT underneath both: exact
  areas, region counts, the distribution and the class shares, from which a mask
  or a bed depth is a decision someone else takes. The CLIMBER-X precedent hands
  a coarse ocean cell the 90th-percentile bed elevation of its ocean part, which
  is a choice this artifact SUPPLIES the quantile for and does not make.

  NOT the below-mesh terrain. GW-6 reserves the relief inside a mesh region, and
  `notes/audits/orogen-resolution.md` finds the generator designs nothing below
  about 20 km, so there is nothing under the mesh to persist. Every number here
  is a statistic of the regions a cell actually holds: `cell_quantiles` returns
  values the mesh has and interpolates between none of them.

  NOT a coarse grid treated as truth. The reference support is the same 10M mesh
  on every rung, whatever the rung, and where a rung is finer than the mesh the
  cells that hold no region are REPORTED rather than filled. This artifact
  applies no nearest-region fallback at all. The mesh's own mean edge is
  measured into the report rather than written here.

## The support, and why it is the mesh rather than the export's own grid fields

The export ships `orog_mean/std/min/max/count` per cell and they are not this.
They are taken over EVERY region in the cell, land and seabed together, so a
coastal cell's mean is dragged toward the seabed and a peak excess measured
against it is inflated by the ocean floor rather than by the terrain. Every
statistic here names its POPULATION, because `lib/gridding.py`'s operators
require one and `CLAUDE.md` rule 1 is why: land is `surface_class` and never
`land_mask`, or the dry closed-basin floors below sea level drop out of the very
distribution this exists to hold. The export's own fields are compared against
this artifact's ALL-region moments and the difference is reported, which is what
says the two binnings are the same binning: both place every region of one mesh
in exactly one cell by its centre, so the totals are an identity and are checked
as one, and a per-cell count may differ only where a region sits on a boundary.

## Two partitions of the sphere, and the artifact carries both

`cell_area_km2` and `mesh_area_km2` are not two estimates of one thing. The
first is the QUADRATURE partition -- the Gauss-Legendre weights laid end to end
in the sine of latitude, which is the weight the spectral model's own global
budget is taken over. The second is what the BINNING delivered, and the binning
cuts at the midpoints between row centres. Those are different partitions and
`lib/gridding.py:export_grid` documents the gap: it is small away from the poles
and does not shrink with the rung, because it is a property of the two
constructions rather than of the mesh. Every share written here is a share of
the mesh area or of the land inside it, so nothing in the file mixes the two,
and the report carries the ratio per rung because a consumer turning a share
into a model-cell quantity crosses exactly that gap. The report separates it
from the other thing the same ratio shows -- the per-cell tails, which ARE the
mesh being coarse against the cell and which widen with the rung until a polar
cell holds no region at all.

## What is written

  data/<build>/support_exoplasim-<rung>.nc   one per requested rung
  analysis/spatial_support_report.json       the report, tracked
  data/<build>/support_checkpoint.json       what is done, for a resumed run

## Batch, chunk, checkpoint, report

`docs/src/reference/large-data.md` binds here: the input is a 10M-region mesh
across five rungs. THE CHUNK IS THE RUNG and it is declared rather than
inherited -- a rung is where the binning changes, so it is the only boundary at
which partial output is a complete artifact rather than half a file. Each rung
is written as it finishes and recorded in the checkpoint against the hashes of
the inputs it was built from, so re-running the same command resumes; kill it at
80 per cent and running it again costs the last rung. Per chunk the script
prints which of how many, elapsed, an estimate of what is left, and resident
memory.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import pathlib
import sys
import time

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, PROJECT_ROOT   # noqa: F401  puts lib/ on the path
import builds
import gridding
import rungs
from orogen import Export, LAND

sys.path.insert(0, str(PROJECT_ROOT))
from lib import provenance                          # noqa: E402

REPORT = ANALYSIS / "spatial_support_report.json"
CHECKPOINT = "support_checkpoint.json"

# THE QUANTILE VECTOR, DECLARED HERE AND WRITTEN INTO EVERY FILE. It is not a
# resolution knob: it states which part of the distribution a consumer needs
# resolved, and on this world that is the TAILS. Permanent ice sits on the top
# per cent of a cell's land and the abyssal floor on the bottom of its ocean, so
# a vector uniform in probability would resolve the middle nobody asked about
# and blur the ends every consumer named. The centre keeps a 0.05 step because a
# median and an interquartile range are what a roughness or a drag scheme takes.
#
# THE STEP IS THE RESOLUTION, and it is the one number a consumer has to carry.
# `area_fraction_above` interpolates the cumulative curve between tabulated
# points, so a share it returns is accurate to the probability STEP bracketing
# the threshold: 0.0025 in the top and bottom per cent, 0.05 across the middle.
# The tails are stepped at 0.0025 because that is where the answers are -- land
# below freezing in the warmest month is a fraction of a per cent of this
# world's land, and a table stepped at 0.05 there would return it as noise.
PROBS = np.array([
    0.0, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02,
    0.03, 0.04, 0.05, 0.07, 0.10,
    0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60,
    0.65, 0.70, 0.75, 0.80, 0.85,
    0.90, 0.93, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.9925, 0.995,
    0.997, 0.998, 0.999, 0.9995, 1.0])
TAIL_STEP = 0.0025      # the coarsest step inside the top and bottom per cent

# The bars, fixed before anything was measured. Every one is an IDENTITY rather
# than an agreement between two estimates, so what it bounds is floating-point
# order and the bar is set by the arithmetic instead of by the result.
CLOSURE_TOL = 1e-12        # partitions and area closures: exact but for summation order
EXTENSIVE_TOL = 1e-9       # a share times its denominator against its own extensive sum
MESH_CLOSURE_TOL = 1e-6    # the mesh's region areas against the manifest's own sphere

# The surface classes, in export code order. `inland_water` is empty by design
# in the export and is carried anyway, because a partition with a member left
# out is not a partition and downstream lake water arrives through exactly it.
SURFACE = ("ocean", "land", "inland_water")

# The drainage classes, over the cell's LAND. They come from `regions.nc`'s
# `terminal`, which is the network this component solved, and they partition the
# land because every land region drains somewhere.
HYDROLOGIC = ("exorheic", "endorheic")


def rss_mb() -> float:
    """Resident memory in MB. One line, and it is what separates a slow job from
    one that will not finish."""
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return float(line.split()[1]) / 1024.0
    except OSError:
        pass
    return float("nan")


def sha256(path) -> str | None:
    path = Path(path)
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def rock_classes(mesh_root: Path) -> list[dict]:
    """The lithology legend, from the export rather than from a copy of it."""
    lit = json.loads((mesh_root / "manifest.json").read_text(encoding="utf-8"))["lithology"]
    return [{"id": int(r["id"]), "code": str(r["code"])} for r in lit["rockClasses"]]


def barren_selector(mesh: Export, codes) -> tuple[np.ndarray, list[str], list[str]]:
    """The barren substrate mask, and which declared classes actually formed it.

    `config/planet.yaml`'s `model.barren_rock_classes` is the declaration, and it
    is the SAME key `exoplasim/scripts/build_surface_albedo.py` reads to mask
    vegetation off barren ground. This artifact persists that SUBSTRATE half as
    `f_nonbarren`; it is not BIO-11's rootable fraction, which additionally
    subtracts solved open water and lives in biosphere/data/<build>/.

    A class the export does not carry is recorded as absent rather than quietly
    skipped: an empty barren mask and a misspelled class code otherwise produce
    identical files.
    """
    legend = {r["code"]: r["id"] for r in rock_classes(mesh.root)}
    rock = np.asarray(mesh.substrate_class, dtype=np.int64)
    mask = np.zeros(rock.shape, dtype=bool)
    applied, absent = [], []
    for code in codes:
        if code in legend:
            mask |= rock == legend[code]
            applied.append(str(code))
        else:
            absent.append(str(code))
    return mask, applied, absent


def hydrologic_labels(terminal: np.ndarray, land: np.ndarray) -> np.ndarray:
    """One drainage class per land region, refusing a land-definition mismatch.

    `terminal` is -2 off the land population, -1 at the world ocean, and the
    preserved basin index otherwise. The two land definitions have to be the
    same land: `build_hydrography.py` takes it from `surface_class` per
    `CLAUDE.md` rule 1 and so does this, so a disagreement is not two
    conventions rounding differently, it is one of the two having read
    `land_mask`. It is refused rather than reconciled.
    """
    hydro_land = terminal != -2
    if not np.array_equal(hydro_land, land):
        raise SystemExit(
            f"regions.nc calls {int(hydro_land.sum()):,} regions land and "
            f"surface_class calls {int(land.sum()):,} land, differing on "
            f"{int((hydro_land ^ land).sum()):,}. Both are supposed to be "
            "surface_class == LAND. One of them has read land_mask, which "
            "floods the dry closed-basin floors this world is built around: "
            "CLAUDE.md rule 1 and source/README.md.")
    labels = np.full(terminal.shape, -1, dtype=np.int8)
    labels[land & (terminal == -1)] = HYDROLOGIC.index("exorheic")
    labels[land & (terminal >= 0)] = HYDROLOGIC.index("endorheic")
    unclaimed = int((land & (labels < 0)).sum())
    if unclaimed:
        raise SystemExit(
            f"{unclaimed:,} land regions carry a terminal that is neither the "
            "world ocean nor a preserved basin. The drainage classes are meant "
            "to partition the land, and this is land that has left the "
            "partition rather than a class to add.")
    return labels


def check_partition(shares, covered, what: str, tol: float = CLOSURE_TOL) -> float:
    """Shares of a partition are non-negative, at most one, and sum to one.

    Refused rather than renormalised: above one is area counted twice and below
    one is area that has left the partition, and renormalising turns either into
    a plausible number. Returns the worst departure.
    """
    covered = np.asarray(covered, dtype=bool)
    total = np.zeros(np.asarray(shares[0]).shape, dtype=np.float64)
    for f in shares:
        f = np.asarray(f, dtype=np.float64)
        bad = covered & ((f < -tol) | (f > 1.0 + tol))
        if bad.any():
            raise SystemExit(
                f"{what}: {int(bad.sum())} cells hold a share outside [0, 1] "
                f"(range {f[covered].min():.6g} to {f[covered].max():.6g})")
        total += f
    dev = np.abs(total[covered] - 1.0)
    worst = float(dev.max()) if dev.size else 0.0
    if worst > tol:
        raise SystemExit(
            f"{what}: the shares do not partition their population; the worst "
            f"cell is off by {worst:.3g} against a declared {tol:.0e}. Above "
            "one is area counted twice, below one is area that has left.")
    return worst


def check_share_against_area(share, denominator, area, covered, what: str) -> float:
    """A categorical share times its denominator IS the extensive area.

    Two operators over one population, so their agreement is an identity with a
    right answer rather than a plausibility check -- and the pair is exactly
    where a population mismatch hides, because a share taken over the land and
    an area summed over everything both look ordinary and differ by the cell's
    ocean. Returns the worst relative departure.
    """
    want = np.asarray(share, dtype=np.float64) * np.asarray(denominator, dtype=np.float64)
    got = np.asarray(area, dtype=np.float64)
    scale = np.maximum(np.asarray(denominator, dtype=np.float64), 1e-30)
    d = (np.abs(got - want) / scale)[np.asarray(covered, dtype=bool)]
    worst = float(d.max()) if d.size else 0.0
    if worst > EXTENSIVE_TOL:
        raise SystemExit(
            f"{what}: the share and the area disagree by {worst:.3g} of the "
            f"denominator against a declared {EXTENSIVE_TOL:.0e}. They are one "
            "population reduced by two operators, so this is a population "
            "mismatch and not a rounding difference.")
    return worst


def check_quantiles(q, covered, lo, hi, what: str) -> None:
    """A quantile table's ends are the population's true extremes, and it does
    not decrease. Both are exact identities, so both are asserted rather than
    toleranced -- an interpolating convention fails the first by construction,
    which is how ground the mesh does not have would get in."""
    c = np.asarray(covered, dtype=bool)
    if not np.array_equal(q[c, 0], lo[c]):
        n = int((q[c, 0] != lo[c]).sum())
        raise SystemExit(f"{what}: Q(0) is not the minimum in {n} cells")
    if not np.array_equal(q[c, -1], hi[c]):
        n = int((q[c, -1] != hi[c]).sum())
        raise SystemExit(f"{what}: Q(1) is not the maximum in {n} cells")
    if not (np.diff(q[c], axis=1) >= 0).all():
        n = int((np.diff(q[c], axis=1) < 0).sum())
        raise SystemExit(f"{what}: the table decreases at {n} points")


def check_inverse_against_mesh(cell, ncell: int, area, values, population,
                               denominator, block, rung: str) -> dict:
    """The whole point of the artifact, tested against the mesh it came from.

    A consumer does not read the hypsometry, it reads a SHARE out of it: how
    much of this cell's land stands above a height. So the test is that read,
    against the exact answer computed on the 10M mesh, which is a right answer
    rather than a second opinion. Everything else here checks that a reduction
    closed; this checks that the reduction can be inverted.

    THE THRESHOLD IS THE CELL'S OWN mean plus one standard deviation. It needs
    no climatology, so this check runs wherever the artifact does, and it lands
    in the upper part of the distribution, which is where the question that
    motivates the artifact lives.

    THE BAR IS A THEOREM AND NOT A MEASUREMENT. `area_fraction_above`
    interpolates the cumulative curve linearly in probability, and the true
    curve is monotone, so both lie between the two probabilities bracketing the
    threshold and the error cannot exceed the widest step in `PROBS`. A run
    that exceeds it has a defect in the operator or in the table, not a coarse
    vector, and it is refused.

    SINGLE-VALUED CELLS ARE COUNTED, NOT AVERAGED IN. A cell whose land is one
    repeated elevation has no distribution: the threshold lands on the table's
    own `Q(0)`, where the convention returns all of the cell rather than none.
    `area_fraction_above` carries the argument; here they are separated so the
    error statistic is about interpolation and the convention is a count.
    """
    values = np.asarray(values, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    pop = np.asarray(population, dtype=bool)
    threshold = block["mean"] + block["sd"]
    sel = pop & (values > threshold[np.asarray(cell)])
    exact_area = gridding.cell_sum(cell, ncell, area, sel)
    live = np.asarray(denominator, dtype=np.float64) > 0
    exact = exact_area[live] / np.asarray(denominator, dtype=np.float64)[live]
    read = gridding.area_fraction_above(block["quantiles"][live], PROBS,
                                        threshold[live])
    d = np.abs(read - exact)
    flat = threshold[live] <= block["min"][live]
    real = ~flat
    w = np.asarray(denominator, dtype=np.float64)[live]
    bar = float(np.diff(PROBS).max())
    worst = float(d[real].max()) if real.any() else 0.0
    if worst > bar + 1e-9:
        raise SystemExit(
            f"{rung}: reading a share back out of the hypsometry is off the "
            f"mesh's own answer by {worst:.4g} in the worst cell, against the "
            f"widest step in PROBS at {bar:.4g}. Linear interpolation of a "
            "monotone curve cannot exceed that step, so this is a defect in the "
            "operator or in the table rather than a coarse quantile vector.")
    return {
        "threshold": "the cell's own land elevation mean plus one sd",
        "bar": bar,
        "bar_is": ("the widest step in PROBS. The interpolated curve and the "
                   "true one both lie between the two bracketing "
                   "probabilities, so this is a bound rather than a tolerance"),
        "cells": int(live.sum()),
        "single_valued_cells": int(flat.sum()),
        "single_valued_land_area_fraction":
            float(w[flat].sum() / w.sum()) if w.sum() > 0 else 0.0,
        "mean_abs_error": (float(np.average(d[real], weights=w[real]))
                           if real.any() else 0.0),
        "p99_abs_error": float(np.percentile(d[real], 99)) if real.any() else 0.0,
        "worst_abs_error": worst,
    }


def extremes(cell, ncell: int, values, population):
    """Per-cell minimum and maximum of a population, for the quantile identity.

    Computed independently of `cell_quantiles` on purpose: checking a table's
    ends against extremes the same sort produced would check nothing.
    """
    sel = np.asarray(population, dtype=bool)
    lo = np.full(ncell, np.inf)
    hi = np.full(ncell, -np.inf)
    np.minimum.at(lo, np.asarray(cell)[sel], np.asarray(values)[sel])
    np.maximum.at(hi, np.asarray(cell)[sel], np.asarray(values)[sel])
    return lo, hi


def population_block(cell, ncell, area, values, population, probs):
    """Moments, extremes, hypsometry and the extensive integral for one population.

    One function because land elevation and ocean depth are the same reduction
    of two populations, and writing it twice is how the two drift apart.
    """
    mean, var, count, covered = gridding.cell_moments(cell, ncell, area, values,
                                                      population)
    q, q_covered = gridding.cell_quantiles(cell, ncell, area, values, probs,
                                           population)
    if not np.array_equal(covered, q_covered):
        raise SystemExit(
            "the moments and the quantile table disagree about which cells hold "
            "the population; they are the same operator's binning and cannot")
    lo, hi = extremes(cell, ncell, values, population)
    integral = gridding.cell_sum(cell, ncell,
                                 np.asarray(area, dtype=np.float64)
                                 * np.asarray(values, dtype=np.float64),
                                 population)
    return {"mean": mean, "sd": np.sqrt(var), "count": count.astype(np.int64),
            "covered": covered, "quantiles": q, "min": lo, "max": hi,
            "integral": integral}


def export_orography(grid_dir: Path):
    """The export's own per-cell orography statistics, or None if it ships none.

    Read for the DIAGNOSTIC only. `vendor/orogen/tools/README.md` is the
    authority on them: they are physical kilometres on the same vertical scale
    as `elevation_km`, with the hypsometric curve and the 1/g relief scaling both
    applied, for every build generated after 2026-08-16.
    """
    with Dataset(grid_dir / "planet.nc") as ds:
        if "orog_mean" not in ds.variables:
            return None
        out = {}
        for name in ("orog_mean", "orog_std", "orog_min", "orog_max", "orog_count"):
            if name in ds.variables:
                out[name] = np.asarray(ds[name][:], dtype=np.float64).ravel()
    return out or None


def build_rung(mesh: Export, grid_dir: Path, rung: str, ctx: dict) -> tuple[dict, dict]:
    """One rung. Returns (arrays to write, the rung's part of the report)."""
    nlat, nlon, _ = rungs.geometry(rung)
    spec = gridding.export_grid(grid_dir, name=f"{rung}-gaussian")
    if spec.shape != (nlat, nlon):
        raise SystemExit(
            f"{grid_dir} is {spec.shape} and lib/rungs.py says {rung} is "
            f"{(nlat, nlon)}; they are not the same grid")
    ncell = nlat * nlon

    # The binning, and the ONLY place a region meets a cell. `CLAUDE.md` rule 3:
    # index, never longitude, and never a second derivation of the column.
    cell, _, _ = gridding.region_cells(mesh, grid_dir)

    area = ctx["area"]
    land = ctx["land"]
    ocean = ctx["ocean"]
    inland = ctx["inland"]
    elev = ctx["elev"]

    # -- area and count. THESE ARE TWO DIFFERENT PARTITIONS OF THE SPHERE and
    # the artifact carries both on purpose. `cell_area_km2` is the QUADRATURE
    # partition, the Gauss-Legendre weights laid end to end in the sine of
    # latitude, which is the weight the spectral model's own global budget is
    # taken over. `mesh_area_km2` is what the BINNING delivered, and the binning
    # is `row()`, which cuts at the midpoints between row centres. The gap is
    # the one `export_grid` documents, it is systematic in the polar row, and it
    # does not shrink with the rung because it belongs to the two constructions
    # rather than to the mesh; the report measures it per rung. A share here
    # is a share of the MESH area or of the land area inside it and never of
    # `cell_area_km2`, so nothing in this file mixes the two; a consumer
    # converting a share into a model-cell quantity needs the model's weight and
    # that is what `cell_area_km2` is for.
    cell_area = spec.cell_area(ctx["radius_m"]).ravel() * 1e-6
    mesh_area = gridding.cell_sum(cell, ncell, area)
    land_area = gridding.cell_sum(cell, ncell, area, land)
    ocean_area = gridding.cell_sum(cell, ncell, area, ocean)
    inland_area = gridding.cell_sum(cell, ncell, area, inland)
    counts = np.bincount(cell, minlength=ncell).astype(np.int64)
    land_count = np.bincount(cell[land], minlength=ncell).astype(np.int64)
    ocean_count = np.bincount(cell[ocean], minlength=ncell).astype(np.int64)
    covered = mesh_area > 0
    empty = ~covered

    # -- the surface-class partition, over the cell's MESH area.
    surface_shares = []
    for name, sel in zip(SURFACE, (ocean, land, inland)):
        f, _ = gridding.cell_fraction(cell, ncell, area, sel)
        surface_shares.append(f)
    surface_worst = check_partition(surface_shares, covered,
                                    f"{rung}: the surface-class shares")

    # -- the partial-surface fractions, over the cell's LAND area.
    land_covered = land_area > 0
    f_barren, _ = gridding.cell_fraction(cell, ncell, area, ctx["barren"], land)
    f_nonbarren = 1.0 - f_barren
    f_lake, _ = gridding.cell_fraction(cell, ncell, area, ctx["lake"], land)
    lake_absent = ctx["lake_absent"]
    barren_area = gridding.cell_sum(cell, ncell, area, ctx["barren"] & land)
    lake_area = gridding.cell_sum(cell, ncell, area, ctx["lake"] & land)
    nonbarren_worst = check_partition([f_barren, f_nonbarren], land_covered,
                                      f"{rung}: nonbarren against barren")
    barren_extensive = check_share_against_area(f_barren, land_area, barren_area,
                                                land_covered, f"{rung}: barren")
    lake_extensive = check_share_against_area(f_lake, land_area, lake_area,
                                              land_covered, f"{rung}: solved lake")
    if lake_absent:
        # The lake solve has not run on this build, so the share is ABSENT and
        # not zero. A zero share and an unsolved balance are the same number and
        # different facts, and this component's standing rule is that an
        # unavailable class is absent and says why rather than being estimated.
        # The schema is kept so a consumer reads one shape on every build.
        f_lake = np.full(ncell, np.nan)
        lake_area = np.full(ncell, np.nan)
        lake_extensive = None

    # -- the class shares, over the cell's LAND area. CATEGORICAL: what survives
    # a reduction of a class field is the retained share of each class, and a
    # majority label discards every minority surface however much leverage it
    # has.
    rock_share = np.zeros((ncell, len(ctx["rocks"])))
    for i, r in enumerate(ctx["rocks"]):
        rock_share[:, i], _ = gridding.cell_fraction(
            cell, ncell, area, ctx["rock"] == r["id"], land)
    rock_worst = check_partition([rock_share[:, i] for i in range(len(ctx["rocks"]))],
                                 land_covered, f"{rung}: the lithology shares")

    hydro_share = np.zeros((ncell, len(HYDROLOGIC)))
    for i in range(len(HYDROLOGIC)):
        hydro_share[:, i], _ = gridding.cell_fraction(
            cell, ncell, area, ctx["hydro"] == i, land)
    hydro_worst = check_partition([hydro_share[:, i] for i in range(len(HYDROLOGIC))],
                                  land_covered, f"{rung}: the drainage shares")

    # -- the distributions. Land elevation and ocean depth are the same
    # reduction of two populations; depth is positive downward so a bathymetric
    # quantile reads the way a consumer asks for it.
    land_block = population_block(cell, ncell, area, elev, land, PROBS)
    depth = -elev
    ocean_block = population_block(cell, ncell, area, depth, ocean, PROBS)
    check_quantiles(land_block["quantiles"], land_block["covered"],
                    land_block["min"], land_block["max"],
                    f"{rung}: the land hypsometry")
    check_quantiles(ocean_block["quantiles"], ocean_block["covered"],
                    ocean_block["min"], ocean_block["max"],
                    f"{rung}: the ocean hypsometry")
    inverse = check_inverse_against_mesh(cell, ncell, area, elev, land,
                                         land_area, land_block, rung)

    # -- the closures. Each has a right answer.
    sphere_km2 = 4.0 * np.pi * (ctx["radius_m"] * 1e-3) ** 2
    cell_area_closure = abs(float(cell_area.sum()) - sphere_km2) / sphere_km2
    if cell_area_closure > CLOSURE_TOL:
        raise SystemExit(
            f"{rung}: the constructed cell areas sum to {cell_area.sum():.6e} km2 "
            f"and the sphere is {sphere_km2:.6e}, off by {cell_area_closure:.3g}. "
            "The grid is not a partition of this planet.")
    mesh_closure = abs(float(mesh_area.sum()) - ctx["mesh_area_total"]) / ctx["mesh_area_total"]
    if mesh_closure > MESH_CLOSURE_TOL:
        raise SystemExit(
            f"{rung}: the mesh area that landed is off the mesh area offered by "
            f"{mesh_closure:.3g}; regions have been lost in the binning")
    land_closure = abs(float(land_area.sum()) - ctx["land_area_total"]) / ctx["land_area_total"]
    if land_closure > CLOSURE_TOL:
        raise SystemExit(
            f"{rung}: the binned land area is off the mesh's own land area by "
            f"{land_closure:.3g}")
    want_volume = float((area[land] * elev[land]).sum())
    got_volume = float(land_block["integral"].sum())
    volume_closure = abs(got_volume - want_volume) / max(abs(want_volume), 1e-30)

    ledger_all = gridding.transfer_ledger(cell, ncell, area)
    ledger_land = gridding.transfer_ledger(cell, ncell, area, land)
    ledger_ocean = gridding.transfer_ledger(cell, ncell, area, ocean)

    # -- the DIAGNOSTIC against the export's own orography. Declared as a
    # diagnostic and not a gate before it was first run: the exporter bins by
    # its own rule and this module bins by `column()` and `row()`, so the two
    # are independently defined and this artifact's correctness rests on the
    # closure identities above. What it can show is that they are the same
    # binning, which is the check nothing else in the tree makes.
    diag = None
    orog = export_orography(grid_dir)
    if orog is not None:
        all_mean, all_var, _, all_cov = gridding.cell_moments(cell, ncell, area, elev)
        lo_all, hi_all = extremes(cell, ncell, elev, np.ones(cell.shape, bool))
        good = all_cov & (orog["orog_count"] > 0)
        # THE BINNING CHECK, and it is the one with a right answer: the export
        # counts the regions whose CENTRE falls in each cell and so does
        # `region_cells`, so the two counts must sum to the same mesh and may
        # differ per cell only where a region sits on a boundary. A convention
        # difference would show as a systematic offset, not as a handful.
        count_diff = counts - orog["orog_count"].astype(np.int64)
        if int(counts.sum()) != int(orog["orog_count"].sum()):
            raise SystemExit(
                f"{rung}: this binning placed {int(counts.sum()):,} regions and "
                f"the export's own placed {int(orog['orog_count'].sum()):,}. "
                "Both are supposed to be every region of one mesh assigned to "
                "exactly one cell, so one of them is not a partition.")
        diag = {
            "cells_compared": int(good.sum()),
            "note": ("the export's orog_* are taken over EVERY region in the "
                     "cell, so they are compared against this artifact's "
                     "ALL-region moments and not against its land ones"),
            "peak_excess_note": (
                "the two peak excesses below are NOT two estimates of one "
                "number and their difference is not an error bar. They differ "
                "in three ways at once: the population inside a cell (every "
                "region against the cell's land), the set of cells (those whose "
                "CENTRE region is land against every cell holding any land), "
                "and the weighting (one vote per cell against one vote per "
                "square kilometre of land). The land one is what a consumer "
                "asking how much colder a cell's high ground is should read, "
                "because each of those three is the correct choice for that "
                "question and none of the export's three is"),
            "region_count_total_here": int(counts.sum()),
            "region_count_total_export": int(orog["orog_count"].sum()),
            "cells_whose_region_count_differs": int((count_diff != 0).sum()),
            "worst_region_count_difference": int(np.abs(count_diff).max()),
            "median_regions_per_cell": float(np.median(counts[good])) if good.any() else 0.0,
            "max_abs_mean_km": float(np.abs(all_mean[good] - orog["orog_mean"][good]).max()),
            "max_abs_sd_km": float(np.abs(np.sqrt(all_var[good]) - orog["orog_std"][good]).max()),
            "max_abs_min_km": float(np.abs(lo_all[good] - orog["orog_min"][good]).max()),
            "max_abs_max_km": float(np.abs(hi_all[good] - orog["orog_max"][good]).max()),
            "export_cells_with_no_region": int((orog["orog_count"] <= 0).sum()),
        }
        # What the export's ALL-region peak excess says against the land-only
        # one. This is the number GRID-2 and PHYS-13 quote, and the difference
        # between the two is the seabed inside a coastal cell.
        gl = orog["orog_count"] > 0
        with Dataset(grid_dir / "planet.nc") as ds:
            gsc = np.asarray(ds["surface_class"][:], dtype=np.int64).ravel()
        gland = gl & (gsc == LAND)
        if gland.any():
            exc = (orog["orog_max"] - orog["orog_mean"])[gland]
            diag["export_all_region_peak_excess_km"] = {
                "mean": float(exc.mean()), "p90": float(np.percentile(exc, 90)),
                "max": float(exc.max())}
        lc = land_block["covered"]
        if lc.any():
            w = land_area[lc]
            exl = (land_block["max"] - land_block["mean"])[lc]
            diag["land_population_peak_excess_km"] = {
                "mean": float(np.average(exl, weights=w)),
                "p90": float(np.percentile(exl, 90)),
                "max": float(exl.max()),
                "weighting": "land area, over cells that hold land"}

    shape = (nlat, nlon)
    arrays = {
        "spec": spec, "nlat": nlat, "nlon": nlon,
        "cell_area_km2": cell_area.reshape(shape),
        "mesh_area_km2": mesh_area.reshape(shape),
        "land_area_km2": land_area.reshape(shape),
        "ocean_area_km2": ocean_area.reshape(shape),
        "inland_water_area_km2": inland_area.reshape(shape),
        "region_count": counts.reshape(shape),
        "land_region_count": land_count.reshape(shape),
        "ocean_region_count": ocean_count.reshape(shape),
        "empty_cell": empty.reshape(shape).astype(np.int8),
        "f_surface": [f.reshape(shape) for f in surface_shares],
        "f_barren": f_barren.reshape(shape),
        "f_nonbarren": f_nonbarren.reshape(shape),
        "f_solved_lake": f_lake.reshape(shape),
        "barren_area_km2": barren_area.reshape(shape),
        "solved_lake_area_km2": lake_area.reshape(shape),
        "substrate_share": rock_share.reshape(nlat, nlon, -1),
        "hydrologic_share": hydro_share.reshape(nlat, nlon, -1),
        "land": land_block, "ocean": ocean_block, "shape": shape,
    }
    report = {
        "grid": grid_dir.name,
        "latitudes": nlat, "longitudes": nlon, "cells": ncell,
        "grid_spec_source": spec.source,
        "cells_with_any_region": int(covered.sum()),
        "cells_with_no_region": int(empty.sum()),
        "cells_with_no_region_area_fraction_of_sphere":
            float(cell_area[empty].sum() / sphere_km2),
        "cells_with_land": int(land_covered.sum()),
        "cells_with_ocean": int((ocean_area > 0).sum()),
        "median_regions_per_cell": float(np.median(counts[covered])) if covered.any() else 0.0,
        "min_regions_per_covered_cell": int(counts[covered].min()) if covered.any() else 0,
        "closure": {
            "cell_area_against_sphere_relative": cell_area_closure,
            "mesh_area_against_export_relative": mesh_closure,
            "land_area_against_mesh_relative": land_closure,
            "land_volume_against_direct_relative": volume_closure,
            "surface_shares_worst_partition_departure": surface_worst,
            "lithology_shares_worst_partition_departure": rock_worst,
            "drainage_shares_worst_partition_departure": hydro_worst,
            "nonbarren_against_barren_worst_departure": nonbarren_worst,
            "barren_share_against_area_worst_relative": barren_extensive,
            "solved_lake_share_against_area_worst_relative": lake_extensive,
            "tolerances": {"partition_and_area": CLOSURE_TOL,
                           "share_against_extensive": EXTENSIVE_TOL,
                           "mesh_against_manifest": MESH_CLOSURE_TOL},
        },
        "two_partitions": {
            "note": ("cell_area_km2 is the QUADRATURE partition the model "
                     "budgets over and mesh_area_km2 is what the midpoint "
                     "binning delivered. They are not the same partition, and "
                     "the ratio below separates two effects that look alike. "
                     "The POLAR ROW total is systematic and is the partitions "
                     "themselves disagreeing, the gap "
                     "lib/gridding.py:export_grid documents; it does not fall "
                     "with the rung. The per-cell TAILS are the mesh being "
                     "coarse against the cell, and they widen with the rung "
                     "until a polar cell can hold no region at all, which is "
                     "what cells_with_no_region counts"),
            "mesh_over_cell_area": {
                str(p): float(np.percentile(
                    (mesh_area[covered] / cell_area[covered]), p))
                for p in (0, 1, 50, 99, 100)},
            "polar_row_mesh_over_cell_area": [
                float(mesh_area.reshape((nlat, nlon))[i].sum()
                      / cell_area.reshape((nlat, nlon))[i].sum())
                for i in (0, -1)],
        },
        "hypsometry_inverse_against_the_mesh": inverse,
        "ledger": {"all_regions": ledger_all, "land": ledger_land, "ocean": ledger_ocean},
        "fallback": {
            "applied": "none",
            "note": ("this artifact substitutes nothing for an empty cell. "
                     "lib/gridding.py:land_weighted falls back to the nearest "
                     "region centre and would have substituted a value in every "
                     "cell counted in cells_with_no_region; the count is "
                     "reported so a consumer taking that door knows its size"),
        },
        "export_orography_diagnostic": diag,
        "absent": ({"solved_lake": lake_absent} if lake_absent else {}),
    }
    return arrays, report


def write_rung(path: Path, a: dict, mesh: Export, grid_dir: Path, rung: str,
               ctx: dict, rung_report: dict) -> None:
    nlat, nlon = a["nlat"], a["nlon"]
    lat, lon, _ = gridding.grid_geometry(grid_dir)
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.title = ("The 10M native mesh under every cell of one ladder rung: "
                    "area, hypsometry and class shares")
        ds.summary = (
            "SPAT-3, extending GRID-2. Per grid cell, statistics of the World "
            "Orogen mesh regions inside it, each named with the POPULATION it "
            "was taken over. The land population is surface_class == LAND and "
            "never land_mask (CLAUDE.md rule 1). This is a SUPPORT artifact: it "
            "is not a land mask, not a bathymetry, and carries nothing below "
            "the mesh.")
        ds.reduction = (
            "lib/gridding.py's operators by field semantics. Areas and volumes "
            "are EXTENSIVE (cell_sum), means are INTENSIVE (cell_mean / "
            "cell_moments, area-weighted, two-pass variance), shares are "
            "CATEGORICAL (cell_fraction) and never a majority label, and the "
            "hypsometry is the DISTRIBUTION operator (cell_quantiles).")
        ds.hypsometry_convention = (
            "Q(p) is the value of the first region whose cumulative AREA share "
            "in the cell reaches p, with no interpolation between regions, so "
            "every tabulated number is ground the mesh has. Q(0) is the "
            "minimum and Q(1) the maximum exactly. Read a share back out with "
            "lib/gridding.py:area_fraction_above.")
        ds.empty_cells = (
            "region_count == 0 marks a cell finer than the mesh. Every "
            "population statistic is zero there and covered by empty_cell; NO "
            "fallback is applied, so a consumer chooses one and records it.")
        ds.below_mesh = (
            "Nothing here describes relief inside a mesh region. That is GW-6's "
            "and notes/audits/orogen-resolution.md finds the generator designs "
            "nothing below about 20 km, so there is nothing to persist.")
        ds.grid_convention = (
            "lon is Orogen's -180..180 labelling of the same columns ExoPlaSim "
            "labels 0..360. The correspondence is BY INDEX and there is no "
            "mapping between the labels: CLAUDE.md rule 3.")
        ds.source_build = mesh.root.parent.name
        ds.vesper_source_build = mesh.root.parent.name
        ds.terrain_hash = mesh.terrain_hash
        ds.grid = grid_dir.name
        ds.rung = rung
        ds.grid_spec_source = a["spec"].source
        ds.barren_rock_classes = json.dumps(ctx["barren_applied"])
        ds.barren_rock_classes_absent = json.dumps(ctx["barren_absent"])
        ds.absent = json.dumps(
            {"solved_lake": ctx["lake_absent"]} if ctx["lake_absent"] else {})
        ds.provenance = json.dumps(ctx["stamp"])
        ds.closure = json.dumps(rung_report["closure"])
        ds.precision = (
            "Shares, the quantile tables and anything that is a ratio are "
            "stored float32; areas, volumes and moments are float64. Every "
            "closure in the `closure` attribute was taken in float64 BEFORE "
            "the cast, so a consumer re-summing the stored shares should "
            "expect float32 rounding and not the tolerance those numbers "
            "report.")
        ds.created = datetime.now(timezone.utc).isoformat()

        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        ds.createDimension("prob", PROBS.size)
        ds.createDimension("rock_class", len(ctx["rocks"]))
        ds.createDimension("surface_class", len(SURFACE))
        ds.createDimension("hydrologic_class", len(HYDROLOGIC))

        def var(name, values, dtype, dims, units, desc, **attrs):
            v = ds.createVariable(name, dtype, dims, zlib=True, complevel=4)
            v.units = units
            v.description = desc
            for k, x in attrs.items():
                setattr(v, k, x)
            v[:] = values
            return v

        var("lat", lat, "f8", ("lat",), "degrees_north",
            "the export's own Gaussian latitudes, north to south. Checked "
            "against lib/gridding.py's constructed rows rather than trusted")
        var("lon", lon, "f8", ("lon",), "degrees_east",
            "Orogen's labelling of the columns; the correspondence to a model "
            "field is by INDEX and never by this label")
        var("prob", PROBS, "f8", ("prob",), "1",
            "cumulative AREA share the quantile tables are taken at")

        # -- area and count
        var("cell_area_km2", a["cell_area_km2"], "f8", ("lat", "lon"), "km2",
            "the cell's area under the QUADRATURE partition: the Gauss-Legendre "
            "weights laid end to end in the sine of latitude. Sums to 4 pi R^2 "
            "and is the weight the spectral model's own global budget is taken "
            "over. Use it to turn a share here into a model-cell quantity")
        var("mesh_area_km2", a["mesh_area_km2"], "f8", ("lat", "lon"), "km2",
            "mesh region area that landed in the cell, EXTENSIVE, and the "
            "denominator every share here is a share of. It is a DIFFERENT "
            "partition from cell_area_km2: the binning cuts at the midpoints "
            "between row centres and the quadrature does not, and the gap is a "
            "property of the two constructions rather than of the mesh, so it "
            "does not shrink with the rung. The report's two_partitions block "
            "carries their ratio and separates that from the other thing it "
            "shows, which is the mesh being coarse against the cell; where a "
            "rung is finer than the mesh the binning delivers nothing at all "
            "and empty_cell says so")
        var("land_area_km2", a["land_area_km2"], "f8", ("lat", "lon"), "km2",
            "land area in the cell, EXTENSIVE, over surface_class == LAND")
        var("ocean_area_km2", a["ocean_area_km2"], "f8", ("lat", "lon"), "km2",
            "ocean area in the cell, EXTENSIVE")
        var("inland_water_area_km2", a["inland_water_area_km2"], "f8",
            ("lat", "lon"), "km2",
            "inland-water area in the cell, EXTENSIVE. Empty in the export by "
            "design; water levels are solved downstream and arrive here")
        var("region_count", a["region_count"], "i4", ("lat", "lon"), "1",
            "native mesh regions in the cell. 0 is a cell finer than the mesh "
            "and every statistic in it is undefined rather than filled")
        var("land_region_count", a["land_region_count"], "i4", ("lat", "lon"), "1",
            "land regions in the cell; what says whether a spread means anything")
        var("ocean_region_count", a["ocean_region_count"], "i4", ("lat", "lon"), "1",
            "ocean regions in the cell")
        var("empty_cell", a["empty_cell"], "i1", ("lat", "lon"), "1",
            "1 where the cell holds no mesh region. No fallback is applied here")

        # -- the surface-class partition
        f_surface = np.stack(a["f_surface"], axis=-1)
        v = var("surface_share", f_surface, "f4", ("lat", "lon", "surface_class"),
                "1",
                "share of the cell's MESH area in each surface class, "
                "CATEGORICAL. Partitions the cell and is checked to")
        v.flag_meanings = " ".join(SURFACE)
        for i, name in enumerate(SURFACE):
            var(f"f_{name}", f_surface[..., i], "f4", ("lat", "lon"), "1",
                f"share of the cell's mesh area that is {name}; the "
                f"{name} column of surface_share, beside it because a partial "
                "land fraction is what a coastline threshold rounds")

        # -- the partial-surface fractions, over the cell's land
        var("f_barren", a["f_barren"], "f4", ("lat", "lon"), "1",
            "share of the cell's LAND area on barren substrate, from "
            "config/planet.yaml model.barren_rock_classes -- the same key "
            "build_surface_albedo.py masks vegetation with")
        var("f_nonbarren", a["f_nonbarren"], "f4", ("lat", "lon"), "1",
            "1 - f_barren over the cell's land. This is only a substrate "
            "partition and is deliberately NOT called rootable: BIO-11's "
            "canonical rootable_fraction also subtracts solved open water")
        v = var("f_solved_lake", a["f_solved_lake"], "f4", ("lat", "lon"), "1",
                "share of the cell's LAND area under the solved lake surface, "
                "from hydrography's surface_water.nc. The lake is a solved "
                "water balance and not a terrain class, so it is the one field "
                "here that needs a climatology; where none has been solved this "
                "is NOT A NUMBER rather than a zero, because a zero share and "
                "an unsolved balance are the same number and different facts")
        if ctx["lake_absent"]:
            v.absent_reason = ctx["lake_absent"]
        var("barren_area_km2", a["barren_area_km2"], "f8", ("lat", "lon"), "km2",
            "barren area, EXTENSIVE. Beside the share because a flux over an "
            "area and a share of a cell take different denominators")
        v = var("solved_lake_area_km2", a["solved_lake_area_km2"], "f8",
                ("lat", "lon"), "km2",
                "solved lake area, EXTENSIVE. A lake evaporation is a depth "
                "over this and never over the land cell")
        if ctx["lake_absent"]:
            v.absent_reason = ctx["lake_absent"]

        # -- class shares
        v = var("substrate_share", a["substrate_share"], "f4",
                ("lat", "lon", "rock_class"), "1",
                "share of the cell's LAND area in each substrate class, "
                "CATEGORICAL. Partitions the land and is checked to; a majority "
                "label would discard every minority lithology")
        v.flag_meanings = " ".join(r["code"] for r in ctx["rocks"])
        v.flag_values = np.array([r["id"] for r in ctx["rocks"]], dtype=np.int32)
        v = var("hydrologic_share", a["hydrologic_share"], "f4",
                ("lat", "lon", "hydrologic_class"), "1",
                "share of the cell's LAND area draining to the world ocean and "
                "to a preserved closed basin, CATEGORICAL, from the drainage "
                "network in regions.nc")
        v.flag_meanings = " ".join(HYDROLOGIC)

        # -- the distributions
        for pop, block, sign in (("land_elevation", a["land"], "elevation"),
                                 ("ocean_depth", a["ocean"], "depth")):
            var(f"{pop}_mean_km", block["mean"].reshape(a["shape"]), "f8",
                ("lat", "lon"), "km",
                f"area-weighted mean {sign} over the cell's "
                f"{'land' if sign == 'elevation' else 'ocean'}, INTENSIVE. "
                "Zero where the cell holds none of that population")
            var(f"{pop}_sd_km", block["sd"].reshape(a["shape"]), "f8",
                ("lat", "lon"), "km",
                "area-weighted standard deviation about the cell's own mean, "
                "taken in two passes so a flat cell returns exactly zero")
            var(f"{pop}_min_km", block["min"].reshape(a["shape"]), "f8",
                ("lat", "lon"), "km", f"least {sign} in the cell")
            var(f"{pop}_max_km", block["max"].reshape(a["shape"]), "f8",
                ("lat", "lon"), "km", f"greatest {sign} in the cell")
            var(f"{pop}_km", block["quantiles"].reshape(*a["shape"], PROBS.size),
                "f4", ("lat", "lon", "prob"), "km",
                f"the cell's {sign} hypsometry: the {sign} at each cumulative "
                "area share in prob. No interpolation between regions, so every "
                "value is ground the mesh has")
        var("land_volume_km3", a["land"]["integral"].reshape(a["shape"]), "f8",
            ("lat", "lon"), "km3",
            "integral of land elevation over land area, EXTENSIVE and signed. "
            "Negative where the cell's land is dry floor below the datum, which "
            "is terrain a land_mask would have flooded")
        var("ocean_volume_km3", a["ocean"]["integral"].reshape(a["shape"]), "f8",
            ("lat", "lon"), "km3",
            "integral of ocean depth over ocean area, EXTENSIVE. The water a "
            "cell holds above its bed, not a circulation quantity")


def load_checkpoint(path: Path, identity: dict) -> dict:
    """The checkpoint, VALIDATED against the inputs before it is trusted.

    A state file that has lost track of which input it belongs to is worse than
    none, because it silently skips work that was never done for this build.
    """
    if not path.is_file():
        return {"identity": identity, "rungs": {}}
    try:
        rec = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {"identity": identity, "rungs": {}}
    if rec.get("identity") != identity:
        print("  checkpoint belongs to different inputs; starting from nothing")
        return {"identity": identity, "rungs": {}}
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(
        description="SPAT-3: the 10M native mesh under every cell of every "
                    "requested ladder rung.")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--rungs", nargs="+", default=None,
                    help="ladder rungs; defaults to the configured resolution")
    ap.add_argument("--force", action="store_true",
                    help="rebuild rungs the checkpoint says are done")
    ap.add_argument("--selftest", action="store_true",
                    help="run the closure checks against fixtures built to "
                         "violate them; needs no build on disk")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    mesh = Export(builds.mesh_export(config))
    data = builds.component_data("hydrography", config, strict=True)
    wanted = [str(r).upper() for r in
              (args.rungs or [str(config["model"]["resolution"])])]
    for r in wanted:
        rungs.geometry(r)                      # refuses anything off the ladder

    inputs = [data / "regions.nc", data / "surface_water.nc"]
    if not inputs[0].is_file():
        raise SystemExit(
            f"{inputs[0]} does not exist. Run hydrography/scripts/"
            "build_hydrography.py first: the drainage network is its product "
            "and this step cannot stand in for it.")
    # THE LAKE IS THE ONE THING HERE THAT NEEDS A CLIMATOLOGY, because a lake
    # extent is a precipitation-against-evaporation balance and not a terrain
    # class. Everything else is terrain and the drainage network, so the
    # artifact is reachable on a build whose lake balance has not been solved,
    # and the share is then ABSENT and says so rather than being written as a
    # zero -- a zero share and an unsolved balance are the same number and
    # different facts.
    lake_absent = None
    if not inputs[1].is_file():
        lake_absent = (
            f"{inputs[1].name} does not exist on this build; the lake extent is "
            "a solved water balance and hydrography/scripts/surface_water.py "
            "needs a climatology to solve it. Run that step and rebuild with "
            "--force to fill the share.")
        print(f"  solved lake share ABSENT: {lake_absent}")

    stamp = provenance.config_stamp(
        config, "hydrography/scripts/build_spatial_support.py", inputs=inputs)
    identity = {
        "terrain_hash": mesh.terrain_hash,
        "regions": int(mesh.n_regions),
        "probs": PROBS.tolist(),
        "inputs": stamp[provenance.INPUT_STAMP_KEY],
        "config_sha256": stamp["config_sha256"],
        "barren_rock_classes": list(config["model"].get("barren_rock_classes") or []),
    }
    ckpt_path = data / CHECKPOINT
    ckpt = load_checkpoint(ckpt_path, identity)

    todo = [r for r in wanted
            if args.force or r not in ckpt["rungs"]
            or not (data / f"support_exoplasim-{r}.nc").is_file()]
    done_already = [r for r in wanted if r not in todo]
    if done_already:
        print(f"already built and unchanged: {', '.join(done_already)}")
    if not todo:
        print("nothing to do")
        return 0

    # THE FOOTPRINT, ESTIMATED BEFORE ANYTHING STARTS rather than after. The
    # per-region cost is the fields held over the whole run plus what the
    # largest rung's sort needs; the arithmetic is short and it is the
    # difference between a scheduling decision and an OOM.
    n = mesh.n_regions
    per_region_b = 8 * 6 + 1 * 3 + 4 * 2       # the held fields, roughly
    sort_b = 8 * 5                              # order, cs, vs, csum, key
    peak_gb = n * (per_region_b + sort_b) / 2**30
    biggest = max(rungs.geometry(r)[0] * rungs.geometry(r)[1] for r in todo)
    print(f"mesh {n:,} regions from {mesh.root.parent.name}/{mesh.root.name}")
    print(f"rungs to build: {', '.join(todo)}  (largest {biggest:,} cells)")
    print(f"estimated peak footprint about {peak_gb:.1f} GB; "
          f"resident now {rss_mb():.0f} MB", flush=True)

    area = mesh.cell_area.astype(np.float64)
    surface = np.asarray(mesh.surface_class, dtype=np.int64)
    land = surface == LAND
    ocean = surface == SURFACE.index("ocean")
    inland = surface == SURFACE.index("inland_water")
    if not (land | ocean | inland).all():
        raise SystemExit(
            "the export carries a surface class outside "
            f"{SURFACE}; the partition here is not this export's")
    elev = np.asarray(mesh.elevation_km, dtype=np.float64)
    rock = np.asarray(mesh.substrate_class, dtype=np.int64)
    barren, barren_applied, barren_absent = barren_selector(
        mesh, config["model"].get("barren_rock_classes") or [])
    if barren_absent:
        print(f"  barren classes not in this export, recorded absent: {barren_absent}")

    with Dataset(inputs[0]) as ds:
        terminal = np.asarray(ds["terminal"][:], dtype=np.int64)
    if lake_absent:
        lake = np.zeros(n, dtype=bool)
    else:
        with Dataset(inputs[1]) as ds:
            lake = np.asarray(ds["lake"][:]).astype(bool)
            sw_hash = getattr(ds, "terrain_hash", None)
        if sw_hash is not None and sw_hash != mesh.terrain_hash:
            raise SystemExit(
                f"surface_water.nc was solved on terrain {sw_hash[:16]} and the "
                f"mesh is {mesh.terrain_hash[:16]}; not the same world")
        if lake.size != n:
            raise SystemExit(
                f"surface_water.nc has {lake.size:,} regions against the mesh's "
                f"{n:,}; not the same support")
    if terminal.size != n:
        raise SystemExit(
            f"regions.nc has {terminal.size:,} regions against the mesh's "
            f"{n:,}; not the same support")
    hydro = hydrologic_labels(terminal, land)

    ctx = {
        "area": area, "land": land, "ocean": ocean, "inland": inland,
        "elev": elev, "rock": rock, "barren": barren, "lake": lake,
        "hydro": hydro, "rocks": rock_classes(mesh.root),
        "lake_absent": lake_absent,
        "barren_applied": barren_applied, "barren_absent": barren_absent,
        "radius_m": float(mesh.radius_km) * 1e3,
        "mesh_area_total": float(area.sum()),
        "land_area_total": float(area[land].sum()),
        "stamp": stamp,
    }
    print(f"  land {ctx['land_area_total'] / ctx['mesh_area_total']:6.2%} of the "
          f"mesh by area, {int(land.sum()):,} regions; barren classes applied "
          f"{barren_applied}")

    report = (json.loads(REPORT.read_text(encoding="utf-8"))
              if REPORT.is_file() else {})
    report.update({
        "script": "hydrography/scripts/build_spatial_support.py",
        "issues": ["SPAT-3", "GRID-2"],
        "source_build": mesh.root.parent.name,
        "terrain_hash": mesh.terrain_hash,
        "regions": int(n),
        "mesh_mean_edge_km": float(np.pi * mesh.radius_km / np.sqrt(n)),
        "probs": PROBS.tolist(),
        "barren_rock_classes_applied": barren_applied,
        "barren_rock_classes_absent": barren_absent,
        "hydrologic_classes": list(HYDROLOGIC),
        "solved_lake_absent": lake_absent,
        "surface_classes": list(SURFACE),
        "chunk": "one ladder rung; see the module docstring",
        "provenance": stamp,
        "created": datetime.now(timezone.utc).isoformat(),
    })
    report.setdefault("rungs", {})

    t0 = time.time()
    for i, rung in enumerate(todo, start=1):
        t1 = time.time()
        grid_dir = builds.grid_export(config, rung)
        arrays, rung_report = build_rung(mesh, grid_dir, rung, ctx)
        out = data / f"support_exoplasim-{rung}.nc"
        write_rung(out, arrays, mesh, grid_dir, rung, ctx, rung_report)
        del arrays
        rung_report["seconds"] = time.time() - t1
        rung_report["output"] = str(out.relative_to(PROJECT_ROOT.parent)) \
            if PROJECT_ROOT.parent in out.parents else str(out)
        rung_report["output_bytes"] = out.stat().st_size
        report["rungs"][rung] = rung_report
        ckpt["rungs"][rung] = {"sha256": sha256(out),
                               "bytes": out.stat().st_size,
                               "seconds": rung_report["seconds"]}
        ckpt_path.write_text(json.dumps(ckpt, indent=2) + "\n")
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(json.dumps(report, indent=2) + "\n")
        elapsed = time.time() - t0
        left = elapsed / i * (len(todo) - i)
        print(f"[{i}/{len(todo)}] {rung:>5} {rung_report['cells']:>7,} cells  "
              f"{rung_report['seconds']:6.1f} s  "
              f"empty {rung_report['cells_with_no_region']:>5}  "
              f"rss {rss_mb():6.0f} MB  elapsed {elapsed:6.1f} s  "
              f"about {left:6.1f} s left", flush=True)
        print(f"        wrote {out} ({out.stat().st_size / 2**20:.1f} MB)",
              flush=True)

    print(f"\nwrote {REPORT}")
    print(f"wrote {ckpt_path}")
    return 0


def _selftest() -> int:
    """Twelve checks, and six are fixtures built to be WRONG in a named way.

    Every rule this script enforces is handed a case that violates it and is
    required to refuse, because a check that has only ever seen correct input
    has not been tested. They run on synthetic arrays and need no build.
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

    rng = np.random.default_rng(20260826)
    ncell = 40
    covered = np.ones(ncell, dtype=bool)
    a = rng.random(ncell)
    b = 1.0 - a
    check("a true partition is accepted",
          check_partition([a, b], covered, "control") <= CLOSURE_TOL)

    ok, msg = refuses(lambda: check_partition([a, b * 1.05], covered, "double count"))
    check("shares summing above one are refused", ok, msg)
    ok, msg = refuses(lambda: check_partition([a * 0.9, b * 0.9], covered, "leaked"))
    check("shares summing below one are refused", ok, msg)
    ok, msg = refuses(lambda: check_partition([a + 1.5, b - 1.5], covered, "outside"))
    check("a share outside [0, 1] is refused", ok, msg)

    denom = rng.uniform(1.0, 100.0, ncell)
    check("a share against its own extensive area is accepted",
          check_share_against_area(a, denom, a * denom, covered, "control")
          <= EXTENSIVE_TOL)
    ok, msg = refuses(lambda: check_share_against_area(
        a, denom, a * denom * 1.01, covered, "population mismatch"))
    check("a share against the WRONG population's area is refused", ok, msg)

    # The land-definition refusal. `terminal` and `surface_class` are both meant
    # to be the same land, and the fixture is one region the elevation-sign test
    # would have flooded -- a dry closed-basin floor below the datum, which is
    # the exact case CLAUDE.md rule 1 exists for.
    land = rng.random(500) < 0.4
    terminal = np.where(land, -1, -2).astype(np.int64)
    check("matching land definitions are accepted",
          bool((hydrologic_labels(terminal, land)[land] == 0).all()))
    flooded = terminal.copy()
    flooded[np.flatnonzero(land)[0]] = -2
    ok, msg = refuses(lambda: hydrologic_labels(flooded, land))
    check("a land definition that floods a dry basin floor is refused", ok, msg)

    # The quantile identity, on a table deliberately broken at one end. Without
    # this the write path would accept a hypsometry whose ends are not the
    # cell's own extremes, which is how an interpolated elevation gets in.
    q = np.sort(rng.random((30, PROBS.size)), axis=1)
    lo, hi = q[:, 0].copy(), q[:, -1].copy()
    cov = np.ones(30, dtype=bool)
    check_quantiles(q, cov, lo, hi, "control")
    bad = q.copy()
    bad[3, 0] -= 0.1
    ok, msg = refuses(lambda: check_quantiles(bad, cov, lo, hi, "not the minimum"))
    check("a table whose end is not the true extreme is refused", ok, msg)

    # THE CHECKPOINT IS VALIDATED, not trusted. A state file that has lost track
    # of which inputs it belongs to is worse than none: it silently skips work
    # that was never done for the build in hand. large-data.md.
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        ck = pathlib.Path(d) / CHECKPOINT
        mine = {"terrain_hash": "abc", "probs": [0.0, 1.0]}
        theirs = {"terrain_hash": "def", "probs": [0.0, 1.0]}
        ck.write_text(json.dumps({"identity": mine, "rungs": {"T42": {}}}))
        check("a checkpoint for these inputs is resumed",
              load_checkpoint(ck, mine)["rungs"] == {"T42": {}})
        check("a checkpoint for other inputs is discarded",
              load_checkpoint(ck, theirs)["rungs"] == {})

    # THE NEGATIVE CONTROL the artifact exists for, and it is GRID-2's own
    # criterion: a cell that is mostly low ground with one per cent of its area
    # standing 4 km up must report that one per cent. A cell mean returns 0.238
    # km and no ice anywhere, which is why `glac` is identically zero over every
    # cell of the run while the manifest records glaciers as on.
    #
    # The bar is the TABLE'S OWN RESOLUTION and is declared above rather than
    # chosen here: `area_fraction_above` interpolates the cumulative curve, so
    # it recovers a share to the probability step bracketing the threshold, and
    # in the top per cent that step is TAIL_STEP. Asking for more than the table
    # resolves would be asking the artifact to invent the distribution between
    # its own points.
    cell = np.zeros(1000, dtype=np.int64)
    elev = np.concatenate([np.full(990, 0.2), np.full(10, 4.0)])
    area = np.ones(1000)
    tbl, _ = gridding.cell_quantiles(cell, 1, area, elev, PROBS)
    above = float(gridding.area_fraction_above(tbl, PROBS, 2.0)[0])
    mean, _, _ = gridding.cell_mean(cell, 1, area, elev)
    flat = np.full((1, PROBS.size), float(mean[0]))
    lost = float(gridding.area_fraction_above(flat, PROBS, 2.0)[0])
    check("a small high fraction survives the reduction a mean loses",
          abs(above - 0.01) <= TAIL_STEP and lost == 0.0 and float(mean[0]) < 2.0,
          f"share above 2 km {above:.4g} against 0.01 (table step "
          f"{TAIL_STEP}); the cell mean is {float(mean[0]):.4g} km and loses "
          f"{lost:.4g}")

    print(f"\n{n_checks} checks, {len(problems)} failed")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
