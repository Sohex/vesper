"""Bin World Orogen mesh regions onto a model grid.

The export already resamples every field onto each grid it emits, but for
boundary conditions that resampling is the wrong tool twice over. Categorical
fields, `surface_class` among them, take the value of the region containing the
cell centre, which at coarse truncation throws away the coastline: a T21 cell
holds roughly 1200 mesh regions. And continuous fields are averaged over *all*
regions in a cell, so a coastal cell's `rock_albedo` is dragged toward open
water's 0.06 and its elevation toward the seabed.

Everything here therefore integrates from the native mesh, which lives in the
T42 export's `raw/`, onto whichever grid is being targeted. The mesh is the same
for every export, so mesh and grid can come from different directories.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from orogen import Export, LAND

# --- the convention, and the only place it is written down ------------------
#
# THIS MODULE OWNS THE COLUMN. Nothing else in the project may derive one, and
# `scripts/smoke_test.py` lints for it, because the same defect has now been
# introduced three times: once as the original bug that superseded the
# `carved-zoned` build, once in the shape of its own fix, and once in a sink
# lookup that the fix's sweep did not reach.
#
# There is nothing to translate. A World Orogen export and the ExoPlaSim grid
# it was integrated onto are the same columns in the same order by
# construction, because `column()` below is the expression that put every mesh
# region in a column in the first place -- for the boundary conditions the
# model reads, for the coupling matrix the catchments are integrated over, and
# for anything sampled per region afterwards. The model then labels that axis
# 0..360 because it has never heard of Orogen, whose own axis reads -180..180.
# A label is not a coordinate correspondence: `CLAUDE.md` rule 3.
#
# A translation layer would be worse than the bug, because it implies there is
# something to translate. So there is no such function here, and
# `display_longitude` -- the only longitude arithmetic left in the project
# outside this block -- is for the eye and says so.


# The ladder registry lives in `lib/rungs.py`, which carries no dependencies so
# that `build_model.py` and the shell probes can reach it without pulling in
# numpy and the export reader. It is re-exported here because this module is
# where a reader looks for anything about the grid, and because one definition
# with two doors is the point: SPAT-2 exists because the ladder had six copies.
from rungs import RUNGS, fft_module, geometry, model_grid, rung_of_latitudes  # noqa: F401,E402


def column(lon, nlon: int) -> np.ndarray:
    """Grid column for a World Orogen longitude. THE column expression.

    Cell `k` spans `[-180 + k*dlon, -180 + (k+1)*dlon)`, which is where the
    export's own column centres sit and where `build_boundary_conditions.py`
    puts the land the model reads back. Index for index the mask the model
    returns is bit-identical to the one it was handed; shifted by half the grid
    it agrees on 0.5955 of cells. See notes/audits/grid-convention-and-runoff.md.
    """
    return np.clip(((np.asarray(lon, dtype=np.float64) + 180.0) / 360.0
                    * nlon).astype(np.int64), 0, nlon - 1)


def column_fraction(lon, nlon: int) -> np.ndarray:
    """Continuous column position for a World Orogen longitude, for resampling.

    `column()` rounds down to the cell a point falls in; this is the same
    expression left continuous and measured from the CENTRE of column 0, which
    is what a bilinear resample onto a finer raster needs. Column `k` runs from
    `k - 0.5` to `k + 0.5` here and the result lies in `[-0.5, nlon - 0.5)`, so
    a caller padding a wrap column at each end indexes at `result + 1`.

    It exists so that resampling a model field is not a reason to write the
    convention out a second time. It was: `maps/build_basemap.py` interpolated
    against the climatology's own 0..360 LABELS and put every climate-derived
    layer on the basemap 180 degrees from the terrain under it.
    """
    return ((np.asarray(lon, dtype=np.float64) + 180.0) / 360.0 * nlon - 0.5)


def row(lat, row_centres) -> np.ndarray:
    """Grid row for a latitude, binned on the midpoints between row centres.

    Rows run north to south and sit at Gauss-Legendre latitudes on a spectral
    grid, so the spacing is not uniform and nearest-centre is not the same
    thing as the cell a point falls in.
    """
    centres = np.asarray(row_centres, dtype=np.float64)
    nlat = centres.size
    edges = np.empty(nlat + 1)
    edges[1:-1] = 0.5 * (centres[:-1] + centres[1:])
    edges[0], edges[-1] = 90.0, -90.0
    return np.clip(np.searchsorted(-edges, -np.asarray(lat, dtype=np.float64),
                                   side="right") - 1, 0, nlat - 1)


def cells(lat, lon, row_centres, nlon: int):
    """(row, column) on a model grid for points in World Orogen coordinates.

    The one entry point for "which cell is this mesh region in". Everything
    that needs a column goes through here or through `coupling_cells`, and
    nothing derives one of its own.
    """
    return row(lat, row_centres), column(lon, nlon)


def coupling_cells(coupling: Path):
    """Decode a coupling matrix into (basin, row, col, area_km2, nlat, nlon).

    The coupling's `cell` is `row * n_lon + col` on the columns `column()`
    assigned, so this decode is the same convention read back rather than a
    second opinion about it.
    """
    from netCDF4 import Dataset
    with Dataset(coupling) as ds:
        basin = np.asarray(ds["basin"][:]).astype(np.int64)
        cell = np.asarray(ds["cell"][:]).astype(np.int64)
        area = np.asarray(ds["area_km2"][:]).astype(np.float64)
        nlat, nlon = int(ds.n_lat), int(ds.n_lon)
    r, c = np.divmod(cell, nlon)
    return basin, r, c, area, nlat, nlon


def require_same_rows(a, b, what: str = "the two grids") -> None:
    """Refuse two latitude axes that are not the same Gaussian latitudes.

    An export's grid and the climatology integrated onto it carry the SAME
    rows, to netCDF's float32 rounding. Asserting that is better than matching
    them by nearest centre, because nearest-centre quietly succeeds on axes
    that are genuinely different and there would be nothing to notice.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape or not np.allclose(a, b, atol=1e-3):
        raise SystemExit(
            f"{what}: latitude axes differ ({a.size} rows against {b.size}"
            + (f", max |da| = {np.abs(a - b).max():.4g} deg" if a.shape == b.shape else "")
            + "). They are not the same grid and no mapping between them is defined.")


def display_longitude(lon, field=None):
    """Rotate a model field onto the -180..180 axis `maps/` draws on. FOR THE EYE.

    ExoPlaSim labels its longitude axis 0..360 and Orogen labels the same
    columns -180..180, so a figure drawn on the model's own labels sits half a
    world away from the same feature in `maps/`.

    **Nothing computed may go through here, and nothing does.** The index is
    canonical and the labels are decoration. This function lives in this module
    rather than beside the plotting code so that the project contains exactly
    one place where a longitude is arithmetic on, and it returns a permutation
    of a field rather than a mapping between grids so that it cannot be
    mistaken for one. Its existence was once read as evidence that a coordinate
    transform was needed between the export and the model, which is the false
    premise that put half of every catchment integral over open ocean, twice.
    """
    lon = np.asarray(lon, dtype=np.float64)
    display = (lon + 180.0) % 360.0 - 180.0
    order = np.argsort(display)
    if field is None:
        return display[order], order
    return display[order], np.asarray(field)[..., order]


def grid_geometry(grid_dir: Path):
    """Latitude and longitude cell centres of an export's grid."""
    gm = json.loads((grid_dir / "manifest.json").read_text(encoding="utf-8"))["grid"]
    lat = np.fromfile(grid_dir / gm["coords"]["lat"]["path"], dtype="float64")
    lon = np.fromfile(grid_dir / gm["coords"]["lon"]["path"], dtype="float64")
    return lat, lon, gm


def region_cells(export: Export, grid_dir: Path):
    """Flat grid-cell index for every mesh region, plus the grid shape."""
    lat, lon, _ = grid_geometry(grid_dir)
    nlat, nlon = lat.size, lon.size
    r, c = cells(export.lat, export.lon, lat, nlon)
    return r * nlon + c, nlat, nlon


def climatology_cells(export: Export, grid_dir: Path, clim_lat):
    """Per-region (row, col) into a climatology field, as (row, col) arrays.

    The mesh-to-grid binning is `region_cells`'s, so a region lands on the same
    cell here as it does in the coupling matrix; only the ROW is then remapped
    onto the climatology's own latitude axis. That asymmetry is the whole point
    and it is deliberate on both sides.

    Columns are NOT remapped. The grid export and the climatology are the same
    columns in the same order and differ only in how they LABEL them, -180..180
    against 0..360. Matching those labels shifts the field by half the planet:
    it put 50.71% of land mesh area onto cells the model calls ocean, against
    7.32% index for index. `CLAUDE.md` rule 3, and `coupling_ocean_fraction`
    below is the invariant that proves the convention held.

    Latitude is CHECKED, not joined. The two axes were once believed to differ,
    and they do not: measured 2026-08-17 they agree to 3.6e-06 degrees, which is
    netCDF float32 rounding of the same Gauss-Legendre rows. So the rows are the
    identity like the columns, and `require_same_rows` asserts it. A
    nearest-centre join returns the same answer here and a wrong one silently on
    axes that are genuinely different, and the bijection guard it carried does
    not catch that, because two different axes of the same length still join
    bijectively.

    One copy, because four modules kept their own version of a path resolver and
    three of them went stale; the same argument applies to a grid convention with
    more at stake. See `lib/paths.py`.
    """
    cell, nlat, nlon = region_cells(export, grid_dir)
    row, col = np.divmod(cell, nlon)
    glat, _, _ = grid_geometry(grid_dir)
    clim_lat = np.asarray(clim_lat, dtype=float)
    if clim_lat.size != glat.size:
        raise ValueError(
            f"climatology has {clim_lat.size} latitude rows and the grid export "
            f"has {glat.size}; they are not the same resolution and no mapping "
            "between them is defined here")
    require_same_rows(glat, clim_lat, "the grid export and the climatology")
    return row, col


# --- field-specific reduction operators. SPAT-4 -----------------------------
#
# A reduction from the mesh to a grid is not one operation. What is correct
# depends on what the field MEANS, and getting that wrong is silent:
#
#   EXTENSIVE   a store or a flux -- cubic metres, kilograms, watts. The cell
#               value is the SUM over the regions it contains, and the global
#               total is preserved exactly. Averaging one destroys the total.
#   INTENSIVE   a density, a rate per area, a temperature, a fraction. The cell
#               value is the AREA-WEIGHTED MEAN over the population it belongs
#               to, and the population is part of the definition: a land field
#               averaged over a coastal cell's ocean is dragged toward the sea.
#   CATEGORICAL a class. There is no mean; what survives is the RETAINED
#               FRACTION of each class, and picking the majority discards every
#               minority surface however much leverage it has.
#   NONLINEAR   a quantity a law is applied to. `f(mean(x))` is not `mean(f(x))`
#               and the gap is real, so the operator is the EXPECTATION:
#               evaluate per region, then reduce. `cell_moments` is the cheaper
#               half of the same correction where the law is smooth and only the
#               variance is wanted.
#
# The four operators below take a binning rather than an export, so they can be
# checked against invariants on synthetic input -- `python lib/gridding.py
# --selftest` -- with no build on disk. The export-facing wrappers underneath
# are thin, and `land_weighted` and `land_fraction_of_class`, which predate
# this block and have consumers, are now expressed in terms of it rather than
# kept as a second implementation.
#
# EVERY ONE OF THEM TAKES THE POPULATION EXPLICITLY. There is no default, and
# `CLAUDE.md` rule 1 is why: land comes from `surface_class` and never from
# `land_mask`, so a reduction that helpfully assumes a land population is a
# reduction that can pick the wrong one. The wrappers pass
# `export.surface_class == LAND` and say so.


def cell_sum(cell, ncell: int, values, population=None) -> np.ndarray:
    """EXTENSIVE reduction: the total of `values` in each cell.

    `values` is already a per-region total -- a volume, a mass, a flux -- so it
    is summed and not weighted by area again. A per-area density becomes
    extensive by multiplying it by `export.cell_area` first, which the caller
    does explicitly because doing it here would make the units of the argument
    ambiguous.

    The invariant is exact: `cell_sum(...).sum()` equals `values[population].sum()`
    to floating-point associativity. `--selftest` asserts it.
    """
    cell = np.asarray(cell)
    values = np.asarray(values, dtype=np.float64)
    sel = slice(None) if population is None else np.asarray(population, dtype=bool)
    return np.bincount(cell[sel], weights=values[sel], minlength=ncell)


def cell_mean(cell, ncell: int, area, values, population=None):
    """INTENSIVE reduction: the area-weighted mean of `values` over a population.

    Returns `(mean, population_area, covered)`. `covered` is the cells that hold
    any of the population at all; the mean is zero elsewhere and a caller must
    mask on `covered` rather than read it, because a cell with no population has
    no mean and a zero is a value.

    Constant preservation is the invariant that makes this a test rather than a
    comparison: a field that is `c` everywhere reduces to `c` in every covered
    cell, exactly, at any binning and any area weighting.
    """
    cell = np.asarray(cell)
    area = np.asarray(area, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    sel = (np.ones(cell.shape, dtype=bool) if population is None
           else np.asarray(population, dtype=bool))
    w = np.bincount(cell[sel], weights=area[sel], minlength=ncell)
    m = np.bincount(cell[sel], weights=area[sel] * values[sel], minlength=ncell)
    covered = w > 0
    mean = np.zeros(ncell)
    np.divide(m, w, out=mean, where=covered)
    return mean, w, covered


def cell_fraction(cell, ncell: int, area, selector, population=None):
    """CATEGORICAL reduction: the share of a cell's population in one class.

    Returns `(fraction, covered)`. Fractions over a partition of the population
    sum to one in every covered cell, which is the invariant `--selftest`
    asserts and which a majority-class rule fails by construction.
    """
    cell = np.asarray(cell)
    area = np.asarray(area, dtype=np.float64)
    pop = (np.ones(cell.shape, dtype=bool) if population is None
           else np.asarray(population, dtype=bool))
    sel = np.asarray(selector, dtype=bool) & pop
    total = np.bincount(cell[pop], weights=area[pop], minlength=ncell)
    hit = np.bincount(cell[sel], weights=area[sel], minlength=ncell)
    covered = total > 0
    frac = np.zeros(ncell)
    np.divide(hit, total, out=frac, where=covered)
    return frac, covered


def cell_moments(cell, ncell: int, area, values, population=None):
    """The first two area-weighted moments of `values` inside each cell.

    Returns `(mean, variance, count, covered)`. This is the statistic a
    cell-scale parameter is allowed to consume -- `hydrography/notes/
    subgrid-water-table.md` section 2, second clause -- and it is here because
    four consumers had each built their own copy of the same three bincounts.

    The variance is the AREA-weighted one and is taken in TWO PASSES, about the
    cell's own mean, rather than as the difference of two moments. The
    difference form is what every consumer had written and it cancels: on a
    field whose spread is far below its mean it returns round-off, and on this
    world's orography that is the common case, since a cell's elevation spread
    is tens of metres about a mean of thousands. The second pass costs one more
    bincount and makes a constant field return exactly zero, which is the
    identity `--selftest` asserts. `count` is the plain region count, which is
    what says whether a cell holds enough of the population for a spread to
    mean anything.
    """
    cell = np.asarray(cell)
    area = np.asarray(area, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    sel = (np.ones(cell.shape, dtype=bool) if population is None
           else np.asarray(population, dtype=bool))
    w = np.bincount(cell[sel], weights=area[sel], minlength=ncell)
    m1 = np.bincount(cell[sel], weights=area[sel] * values[sel], minlength=ncell)
    count = np.bincount(cell[sel], minlength=ncell)
    covered = w > 0
    mean = np.zeros(ncell)
    np.divide(m1, w, out=mean, where=covered)
    residual = values[sel] - mean[cell[sel]]
    m2 = np.bincount(cell[sel], weights=area[sel] * residual ** 2, minlength=ncell)
    var = np.zeros(ncell)
    np.divide(m2, w, out=var, where=covered)
    return mean, np.maximum(var, 0.0), count, covered


def cell_expectation(cell, ncell: int, area, law, values, population=None):
    """NONLINEAR reduction: evaluate the law per region, THEN area-average it.

    `law` is a callable applied to the per-region `values`. Taking the callable
    rather than a pre-evaluated array is the whole point of the operator: a
    caller who has already reduced cannot pass the reduced field in by
    accident, so the order is a property of the call and not of the discipline
    of whoever wrote it.

    Returns `(expectation, mean_of_input, covered)` so the Jensen gap
    `expectation - law(mean_of_input)` is available at the call site without
    computing the binning twice. That gap is what `analysis/
    spatial_reduction_gap.py` sizes per field; use this operator where it sized
    the gap above the consuming step's own instrument, and `cell_mean` where it
    did not, because an expectation is not free and a reduction that does not
    need one should not carry one.
    """
    values = np.asarray(values, dtype=np.float64)
    mean, _, covered = cell_mean(cell, ncell, area, values, population)
    expectation, _, _ = cell_mean(cell, ncell, area, np.asarray(law(values),
                                                               dtype=np.float64),
                                  population)
    return expectation, mean, covered


def transfer_ledger(cell, ncell: int, area, population=None) -> dict:
    """What a reduction onto this binning carries, and what it does not.

    A reduction that returns a finite array everywhere has said nothing about
    what it dropped. This is the ledger the audit's finding 9 asks every
    conversion to emit: the area offered, the area that landed, the cells that
    hold none of the population, and the closure residual. It is cheap and it
    belongs in the provenance of any artifact a reduction wrote.
    """
    cell = np.asarray(cell)
    area = np.asarray(area, dtype=np.float64)
    sel = (np.ones(cell.shape, dtype=bool) if population is None
           else np.asarray(population, dtype=bool))
    binned = np.bincount(cell[sel], weights=area[sel], minlength=ncell)
    offered = float(area[sel].sum())
    landed = float(binned.sum())
    return {
        "regions_offered": int(sel.sum()),
        "cells": int(ncell),
        "cells_with_no_population": int((binned <= 0).sum()),
        "area_offered": offered,
        "area_landed": landed,
        "closure_residual_relative": (abs(landed - offered) / offered
                                      if offered > 0 else 0.0),
    }


def land_weighted(export: Export, grid_dir: Path, values: np.ndarray):
    """Per-cell land fraction and the land-area-weighted mean of `values`.

    Returns (land_fraction, mean_over_land), both shaped (nlat, nlon). Cells with
    no land get a mean of zero, which callers should mask rather than use.
    """
    cell, nlat, nlon = region_cells(export, grid_dir)
    ncell = nlat * nlon
    area = export.cell_area.astype(np.float64)
    is_land = export.surface_class == LAND

    # `cell_mean` and `cell_fraction` above are the operators; this is the land
    # case of them. One implementation, so a fix to the reduction reaches every
    # caller of either door.
    mean, land_area, _ = cell_mean(cell, ncell, area, values, is_land)
    total = cell_sum(cell, ncell, area)
    empty = total <= 0
    fraction = np.zeros(ncell)
    np.divide(land_area, total, out=fraction, where=~empty)

    if empty.any():
        # Cells finer than the mesh, normally narrow polar Gaussian cells. The
        # old 2.5M build had 35 at T85, but that is not a portable count: SPAT-3
        # measures and reports it separately for every grid against the 10M
        # fine-support reference. Fall back to the nearest region centre, which
        # is what the exporter's own categorical rule does; callers must retain
        # `empty` so the substitution cannot disappear from provenance.
        nearest = _nearest_region(export, grid_dir, np.flatnonzero(empty), nlat, nlon)
        fraction[empty] = (export.surface_class[nearest] == LAND).astype(float)
        mean[empty] = np.where(export.surface_class[nearest] == LAND,
                               values[nearest], 0.0)
    return fraction.reshape(nlat, nlon), mean.reshape(nlat, nlon), int(empty.sum())


def _nearest_region(export: Export, grid_dir: Path, cells, nlat: int, nlon: int):
    """Region index nearest each listed flat cell index, on the unit sphere."""
    from scipy.spatial import cKDTree
    lat, lon, _ = grid_geometry(grid_dir)
    rows, cols = np.divmod(cells, nlon)
    clat = np.deg2rad(lat[rows])
    clon = np.deg2rad(lon[cols])
    pts = np.stack([np.cos(clat) * np.sin(clon), np.sin(clat),
                    np.cos(clat) * np.cos(clon)], axis=1)
    tree = cKDTree(np.stack([export.x, export.y, export.z], axis=1))
    return tree.query(pts)[1]


def land_fraction_of_class(export: Export, grid_dir: Path, mask: np.ndarray):
    """Fraction of each cell's *land* area for which `mask` is true.

    The land case of `cell_fraction`, which is the categorical operator: what
    survives a reduction of a class field is the retained fraction, and a
    caller wanting several classes should take one fraction per class rather
    than a majority code, because the fractions of a partition sum to one and a
    majority does not.
    """
    cell, nlat, nlon = region_cells(export, grid_dir)
    frac, _ = cell_fraction(cell, nlat * nlon, export.cell_area.astype(np.float64),
                            mask, export.surface_class == LAND)
    return frac.reshape(nlat, nlon)

# --- the grid convention, and the one check that proves it ------------------

def coupling_path(data_dir, config) -> Path:
    """The coupling matrix for the CONFIGURED rung, not for T42.

    Three scripts wrote `coupling_exoplasim-T42.nc` as a literal, which is
    SPAT-2's finding -- the artifact paths still carry resolution literals -- and
    it is not a cosmetic one: `require_index_alignment` compares the coupling's
    column count against the land mask's and raises, so a T21 run does not read
    a T42 coupling silently. It refuses, which is the good failure. But it
    refuses at every caller, and the literal has to move once rather than three
    times, so it lives here beside the code that owns the mapping.
    """
    resolution = str(config["model"]["resolution"]).upper()
    return Path(data_dir) / f"coupling_exoplasim-{resolution}.nc"


COUPLING_OCEAN_FRACTION_LIMIT = 0.10


def coupling_ocean_fraction(coupling, lsm) -> float:
    """Share of coupling catchment area landing on cells the model calls ocean.

    Lives here because this module owns the index convention: `region_cells`
    places a mesh region at `col = (lon + 180)/360 * nlon`, and everything
    downstream inherits that column. So this is the invariant that proves the
    convention held, and it belongs beside the thing it is about rather than
    being reimplemented per consumer.

    Endorheic catchments are inland, so almost none of their area may fall on a
    cell the model calls ocean. Index for index it is 1.01%; matching longitude
    LABELS instead gives 49.77%, because an ExoPlaSim climatology numbers its
    axis 0..360 and Orogen numbers the same columns -180..180. A label is not a
    coordinate correspondence -- `CLAUDE.md` rule 3.

    This is the check that was missing while the mapping was wrong twice. What
    existed asserted that a `cell_lon` variable EXISTED, which no wrong mapping
    would ever have failed. See `docs/src/practice/failure-modes.md` class 17.
    """
    _, r, c, area, _, nlon = coupling_cells(coupling)
    if lsm.shape[1] != nlon:
        raise ValueError(
            f"coupling has {nlon} columns and the mask has {lsm.shape[1]}; "
            "they are not the same grid and no mapping between them is defined")
    return float(area[lsm[r, c] < 0.5].sum() / area.sum())


def require_index_alignment(coupling, lsm) -> float:
    """Refuse to integrate over a coupling that is not index-aligned."""
    frac = coupling_ocean_fraction(coupling, lsm)
    if frac > COUPLING_OCEAN_FRACTION_LIMIT:
        raise SystemExit(
            f"{frac:.1%} of coupling catchment area lands on cells the model "
            f"calls ocean, against a limit of {COUPLING_OCEAN_FRACTION_LIMIT:.0%}. "
            "The coupling and the climatology are not index-aligned. Do NOT "
            "'fix' this by matching longitude labels; that is what produced it "
            "twice. See notes/audits/grid-convention-and-runoff.md.")
    return frac


# --- the operators' own checks. SPAT-4 --------------------------------------
#
# `python lib/gridding.py --selftest`
#
# Every one of these has a RIGHT ANSWER rather than a plausible one: an
# identity, a conservation law, or a control that must fail. That is the bar
# `docs/src/practice/failure-modes.md` class 17 sets, and it is why these run on
# synthetic input and need no build on disk -- a check that can only be run
# where the answer is unknown is not a check.


_CHECKS = 10


def _selftest() -> int:
    rng = np.random.default_rng(20260824)
    n, ncell = 20000, 97
    cell = rng.integers(0, ncell, n)
    area = rng.uniform(0.1, 10.0, n)
    population = rng.random(n) < 0.6
    values = rng.normal(3.0, 2.0, n)
    problems: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    # EXTENSIVE: the global total is preserved exactly.
    store = values * area
    total = cell_sum(cell, ncell, store, population).sum()
    want = store[population].sum()
    check("extensive sum preserves the global total",
          abs(total - want) <= 1e-9 * abs(want), f"{total} against {want}")

    # INTENSIVE: a constant field reduces to that constant in every covered cell.
    mean, _, covered = cell_mean(cell, ncell, area, np.full(n, 4.25), population)
    check("intensive mean preserves a constant",
          bool(np.allclose(mean[covered], 4.25, rtol=0, atol=1e-12)),
          f"max deviation {np.abs(mean[covered] - 4.25).max():.3g}")

    # INTENSIVE: the area-weighted global mean is recovered from the cell means.
    mean, w, covered = cell_mean(cell, ncell, area, values, population)
    got = float((mean[covered] * w[covered]).sum() / w[covered].sum())
    want = float((values[population] * area[population]).sum() / area[population].sum())
    check("intensive mean closes over the population",
          abs(got - want) <= 1e-12 * abs(want), f"{got} against {want}")

    # CATEGORICAL: the fractions of a partition sum to one, everywhere covered.
    klass = rng.integers(0, 4, n)
    fracs = [cell_fraction(cell, ncell, area, klass == k, population)[0]
             for k in range(4)]
    covered = cell_fraction(cell, ncell, area, klass == 0, population)[1]
    total = sum(fracs)
    check("categorical fractions of a partition sum to one",
          bool(np.allclose(total[covered], 1.0, rtol=0, atol=1e-12)),
          f"max deviation {np.abs(total[covered] - 1.0).max():.3g}")

    # MOMENTS: a constant field has exactly zero variance; and the moments agree
    # with a direct per-cell computation on a cell picked at random.
    # The bar is `(1e-12 * mean)^2`, and it discriminates: the two-pass form
    # lands near `(eps * mean)^2` because the cell mean of a constant is not
    # bitwise the constant, while the difference-of-moments form this replaced
    # returns about `(1e-7 * mean)^2` and misses the bar by nine orders.
    const = 1.5
    _, var, _, covered = cell_moments(cell, ncell, area, np.full(n, const), population)
    check("moments give a constant field zero variance to round-off",
          float(var[covered].max()) <= (1e-12 * const) ** 2,
          f"max {var[covered].max():.3g}")
    m, var, count, covered = cell_moments(cell, ncell, area, values, population)
    k = int(np.flatnonzero(covered & (np.bincount(cell[population], minlength=ncell) > 50))[0])
    sel = population & (cell == k)
    direct_m = float(np.average(values[sel], weights=area[sel]))
    direct_v = float(np.average((values[sel] - direct_m) ** 2, weights=area[sel]))
    check("moments agree with a direct per-cell computation",
          abs(m[k] - direct_m) <= 1e-10 * abs(direct_m)
          and abs(var[k] - direct_v) <= 1e-8 * abs(direct_v),
          f"{m[k]}/{var[k]} against {direct_m}/{direct_v}")

    # NONLINEAR: the expectation of a LINEAR law is the law of the mean, exactly
    # -- there is no Jensen term to find where there is no curvature.
    exp, mean, covered = cell_expectation(cell, ncell, area,
                                          lambda x: 2.5 * x - 1.0, values, population)
    check("expectation of a linear law equals the law of the mean",
          bool(np.allclose(exp[covered], 2.5 * mean[covered] - 1.0, rtol=0, atol=1e-10)),
          f"max deviation {np.abs(exp[covered] - (2.5 * mean[covered] - 1.0)).max():.3g}")

    # THE NEGATIVE CONTROL, and it is the point of the whole block: a CONVEX law
    # must show a strictly positive gap wherever the cell has any spread at all.
    # If this passes, the operator is not measuring anything and the linear check
    # above would pass on a broken implementation too.
    exp, mean, covered = cell_expectation(cell, ncell, area,
                                          lambda x: np.exp(x), values, population)
    _, var, _, _ = cell_moments(cell, ncell, area, values, population)
    spread = covered & (var > 1e-6)
    gap = exp[spread] - np.exp(mean[spread])
    check("convex control: the expectation exceeds the law of the mean",
          bool(spread.any()) and bool((gap > 0).all()),
          f"{int((gap <= 0).sum())} of {gap.size} cells did not")

    # THE LEDGER closes: every region offered lands in exactly one cell.
    ledger = transfer_ledger(cell, ncell, area, population)
    check("the transfer ledger closes on area",
          ledger["closure_residual_relative"] <= 1e-12,
          f"{ledger['closure_residual_relative']:.3g}")

    # PERIODICITY: the column expression wraps rather than folding. A longitude
    # just inside -180 and one just inside +180 are adjacent columns, not the
    # same one, and neither escapes the grid.
    nlon = 128
    edge = column(np.array([-179.999, 179.999, -180.0, 180.0]), nlon)
    check("the column expression is periodic and in range",
          bool((edge >= 0).all() and (edge < nlon).all()
               and edge[0] == 0 and edge[1] == nlon - 1),
          f"{edge}")

    print(f"\n{_CHECKS} checks, {len(problems)} failed")
    return 1 if problems else 0


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true",
                    help="run the reduction operators against their invariants "
                         "on synthetic input; needs no build on disk")
    args = ap.parse_args()
    if not args.selftest:
        ap.error("this module is a library; --selftest is the only thing to run")
    return _selftest()


if __name__ == "__main__":
    raise SystemExit(main())
