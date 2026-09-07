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
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from orogen import Export, LAND

# One full turn of longitude, in degrees. Named because this module is the only
# one allowed to do arithmetic with it and a bare literal reads as a magic number.
_FULL_TURN = 360.0

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
# something to translate. So there is no function here that maps one axis onto
# the other, and `display_longitude` is a permutation for the eye and says so.
#
# What there IS, because artifacts do carry the model's labels, is one statement
# of that label axis and one inverse from a label back to the INDEX it names:
# `model_longitude_labels`, `model_label_column`, `model_label_cells` and
# `require_model_labels`, below. A per-cell product written from a climatology
# -- `pedology/data/<build>/land_column_states_<res>.txt` is one -- is keyed by
# those labels, and reading it needs the index they stand for. The inverse
# REFUSES a label that is not on the axis rather than snapping it to the nearest
# column, which is what makes an axis handed in on the export's own centres an
# error instead of a field rotated by half the planet.


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


def export_grid(grid_dir: Path, name: str | None = None) -> "GridSpec":
    """The cell boundaries of an export's gridded output, CONSTRUCTED not read.

    An export ships three area quantities and only one of them is the weight for
    a gridded field: `raw/cell_area.bin` is the mesh region's dual area,
    `grid/gauss_weights.bin` is the quadrature the spectral model's global budget
    is taken over, and `grid/grid_cell_area.bin` is meant to be the second of
    those in km2. Exports written before the generator was corrected instead
    carry the NEAREST-ROW partition there: a valid partition of the sphere that
    closes to 4 pi R^2 and reports a global mean the model does not take, wider
    than the quadrature interval by 22 per cent in the polar row at every
    truncation and not converging. `notes/audits/ocean-grid-crossing.md` has the
    measurement.

    So this constructs the partition from the manifest's `gridType` and shape
    rather than reading a field, which makes it right on both sides of that
    change, and then checks the rows it built against the axis the export
    actually ships so a disagreement is an error instead of a quiet wrong
    weight. Take `spec.cell_area(radius_m)` for areas and
    `spec.cell_area_fraction()` where the radius would cancel.
    """
    lat, lon, gm = grid_geometry(grid_dir)
    nlat, nlon = lat.size, lon.size
    kind = str(gm.get("gridType", "uniform"))
    if kind == "gaussian":
        spec = gaussian_grid(nlat, nlon, name=name or f"{grid_dir.name}-gaussian")
        require_gaussian_rows(spec, lat, what=f"the export grid at {grid_dir}")
    elif kind == "uniform":
        lat_edges = 90.0 - np.arange(nlat + 1, dtype=np.float64) * (180.0 / nlat)
        sin_edges = np.sin(np.deg2rad(lat_edges))
        sin_edges[0], sin_edges[-1] = 1.0, -1.0
        lon_edges = -180.0 + np.arange(nlon + 1, dtype=np.float64) * (360.0 / nlon)
        spec = GridSpec(name=name or f"{grid_dir.name}-uniform",
                        lon_edges=lon_edges, sin_edges=sin_edges,
                        source="lib/gridding.py:export_grid, equally spaced in latitude")
        want = 90.0 - (np.arange(nlat) + 0.5) * (180.0 / nlat)
        if float(np.abs(want - lat).max()) > 1e-6:
            raise ValueError(
                f"{grid_dir}: manifest says gridType=uniform but grid/lat.bin is not "
                "equally spaced; the cell boundaries this would build are not the "
                "export's own")
    else:
        raise ValueError(
            f"{grid_dir}: gridType {kind!r} is not one this module builds boundaries "
            "for. Add the arm rather than defaulting one, because the wrong "
            "partition closes to the sphere just as exactly as the right one.")
    return spec


def region_cells(export: Export, grid_dir: Path):
    """Flat grid-cell index for every mesh region, plus the grid shape."""
    lat, lon, _ = grid_geometry(grid_dir)
    nlat, nlon = lat.size, lon.size
    r, c = cells(export.lat, export.lon, lat, nlon)
    return r * nlon + c, nlat, nlon


def spec_cells(spec: "GridSpec", lat, lon) -> np.ndarray:
    """Flat cell index on ANY constructed grid for points at `lat`, `lon`.

    THE MESH-TO-GRID DOOR FOR A GRID THE EXPORTER DOES NOT WRITE, and only for
    those. `region_cells` is the door for a grid the export ships: it reads that
    grid's own axis off the manifest and bins to the row the export itself used,
    so a region lands where the export put it. GOLDSTEIN's rows and columns
    exist only as a constructor in this module, and a mesh region can then be
    placed only against `GridSpec`'s cell BOUNDARIES.

    DO NOT REACH FOR IT ON THE EXPORT'S OWN GAUSSIAN GRID. The two doors
    disagree there by construction and neither is a rounding of the other: a
    Gaussian node is a quadrature abscissa and not the centre of its cell, so
    binning to the nearest node and binning between the quadrature edges put a
    band of mesh either side of every row boundary in different rows.
    `ocean/scripts/build_ocean_grid.py` measures the band on the active build.

    RULE 3, AND THE FORM IN WHICH THIS IS NOT IT. What rule 3 forbids is
    matching one grid's longitude LABEL against another's: two conventions
    naming the same columns of ONE grid, where the mapping is the index and a
    label match rotates the field by half the planet. A mesh and a grid share no
    index, so a placement has to go through coordinates, and the question is
    which FRAME those coordinates are in. `lat` and `lon` are the generator's
    own, reproducing from `x, y, z` to within float32, and they are what defines
    this world's prime meridian. `spec.lon_edges` is `phi0` plus whole columns,
    and `phi0` is a rotation of the ocean grid ON that planet rather than a
    second naming of it -- which is what `analysis/ocean_remap.py`'s origin
    check asserts, by rolling the answer one column when `phi0` moves one
    column. One frame, two partitions: `lib/remap.py:_periodic_overlap` already
    compares the same two constructors on the same argument, and folding the
    source into the destination's window by a WHOLE number of turns is the same
    shift it applies.

    THE CONDITION THAT KEEPS THAT TRUE is that the ocean grid's own topography
    is written through this placement, so the frame the ocean model runs in is
    the generator's by construction rather than by agreement. A `.k1` built any
    other way carries a frame of its own and the comparison is then between two
    of them.

    Rows are binned in the SINE of latitude, because that is what a `GridSpec`
    carries and because binning in the angle puts a region in the wrong row of
    an equal-area grid. Either row order is accepted -- this project's grids run
    north to south and GOLDSTEIN's `j` runs south to north -- and the index
    comes back in the spec's OWN order, which is the order the model that reads
    it wants.
    """
    if not isinstance(spec, GridSpec):
        raise TypeError(
            "spec_cells bins against cell BOUNDARIES, so it takes a GridSpec. "
            "Centre arrays read off a file are labels, not a partition.")
    s = np.sin(np.deg2rad(np.asarray(lat, dtype=np.float64)))
    edges = spec.sin_edges
    if edges[0] < edges[-1]:                       # south to north
        r = np.searchsorted(edges, s, side="right") - 1
    else:                                          # north to south
        r = np.searchsorted(-edges, -s, side="right") - 1
    np.clip(r, 0, spec.nlat - 1, out=r)
    turn = float(spec.lon_edges[-1] - spec.lon_edges[0])
    x = np.asarray(lon, dtype=np.float64)
    x = x - turn * np.floor((x - spec.lon_edges[0]) / turn)
    c = np.searchsorted(spec.lon_edges, x, side="right") - 1
    np.clip(c, 0, spec.nlon - 1, out=c)
    return r * spec.nlon + c


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
#   DISTRIBUTION a field whose consumers ask what SHARE of a cell lies past a
#               threshold rather than what its average is. No moment answers
#               that on a skewed cell, so what survives is the area-weighted
#               QUANTILE TABLE -- `cell_quantiles`, with `area_fraction_above`
#               reading a share back out of it. This is the hypsometry GRID-2
#               found three implementations of and no artifact carrying.
#
# The operators below take a binning rather than an export, so they can be
# checked against invariants on synthetic input -- `python lib/gridding.py
# --selftest` -- with no build on disk. The export-facing wrappers underneath
# are thin, and `land_weighted` and `land_fraction_of_class`, which predate
# this block and have consumers, are now expressed in terms of it rather than
# kept as a second implementation.
#
# EVERY REDUCTION IN THIS BLOCK IS NAMED `cell_<what it reduces to>`, AND THE
# PREFIX IS LOAD-BEARING RATHER THAN A HABIT. `config/spatial_support.yaml`'s
# `aggregation_operator_reductions` says which vocabulary term each one is, and
# `scripts/spatial_support_gate.py` requires that every module-level `cell_*`
# callable here is claimed by at least one of them, and that every name a term
# claims exists here. The relation is many-to-one on purpose -- `cell_fraction`
# is both the vocabulary's `area_weighted_fraction` and its
# `categorical_histogram`, because one class share and a partition of them are
# the same computation under two declared semantics -- so what the gate checks
# is COVERAGE, which is the direction that drifted. They drifted once, when this
# module gained the moment, expectation and
# distribution operators and the contract vocabulary could still only name a
# mean and a fraction, so an artifact carrying a hypsometry had to declare
# itself a cell mean or fail validation. A new reduction added here fails that
# gate until the vocabulary can name what it is; a helper that is not a
# reduction -- `area_fraction_above` reads an answer back out of one,
# `transfer_ledger` says what a binning dropped -- takes a name without the
# prefix and is outside the comparison.
#
# DURATION is not here. A mean over time is `lib/climatology.py`'s time-bin
# weights, and a field that is both -- an area mean of a time mean -- takes them
# in that order and from that module. Two places, because a spatial support and
# a temporal one are different facts about an artifact and collapsing them into
# one operator would hide which of the two a caller got wrong.
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



def cell_quantiles(cell, ncell: int, area, values, probs, population=None):
    """DISTRIBUTION reduction: the area-weighted quantiles of `values` per cell.

    Returns `(quantiles, covered)` with `quantiles` shaped `(ncell, len(probs))`
    and zero where the cell holds none of the population -- which a caller masks
    on `covered` rather than reads, for the reason `cell_mean` gives.

    This is the operator a mesh-under-cell hypsometry IS, and it is here rather
    than in a consumer because the same criterion was computed three separate
    ways off a field no artifact carried (GRID-2). A moment is not a substitute
    for it: a cell's land is not normal about its mean, and what a consumer asks
    is what SHARE of the cell lies above a height, which only the distribution
    answers.

    THE DEFINITION, which a caller has to know because quantile conventions
    differ and the difference is not small on a skewed cell. Order the
    population in a cell by value and let `F_k` be the cumulative AREA share up
    to and including region `k`. Then `Q(p)` is the value of the first region
    with `F_k >= p`. It is a sample quantile of an area-weighted distribution,
    with NO INTERPOLATION BETWEEN REGIONS: an interpolated one would report an
    elevation no ground in the cell has, which is the below-mesh terrain GW-6
    reserves. `Q(0)` is the minimum and `Q(1)` the maximum, exactly, and that is
    the identity `--selftest` asserts.

    `probs` is the caller's and is declared beside the result rather than
    defaulted here: a quantile vector states which part of the distribution a
    consumer needs resolved, and on this world the tails are where the ice and
    the abyssal floor are.
    """
    probs = np.asarray(probs, dtype=np.float64)
    if probs.ndim != 1 or probs.size == 0:
        raise ValueError("probs must be a non-empty one-dimensional array")
    if not np.all(np.diff(probs) > 0) or probs[0] < 0.0 or probs[-1] > 1.0:
        raise ValueError(
            "probs must strictly increase inside [0, 1]; an unordered vector "
            "returns an unordered quantile table and nothing says so")
    cell = np.asarray(cell)
    area = np.asarray(area, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    sel = (np.ones(cell.shape, dtype=bool) if population is None
           else np.asarray(population, dtype=bool))
    total = np.bincount(cell[sel], weights=area[sel], minlength=ncell)
    covered = total > 0
    out = np.zeros((ncell, probs.size), dtype=np.float64)
    if not covered.any():
        return out, covered

    # Cell major, value minor. `lexsort` takes its keys least-significant first.
    order = np.lexsort((values[sel], cell[sel]))
    cs = cell[sel][order]
    vs = values[sel][order]
    csum = np.cumsum(area[sel][order])

    # The cumulative area share WITHIN each cell, in (0, 1]. `starts` is where a
    # cell's block begins in the sorted array, so the running total before the
    # block is what comes off.
    counts = np.bincount(cs, minlength=ncell)
    ends = np.cumsum(counts)
    base = np.concatenate([[0.0], csum])[ends - counts]
    frac = (csum - base[cs]) / total[cs]
    # The last region of a cell is its maximum by construction, so its share is
    # exactly one rather than whatever `csum / total` rounded to. Without this
    # the `Q(1) == max` identity fails by a rounding at a handful of cells.
    frac[ends[counts > 0] - 1] = 1.0

    # ONE strictly increasing key over the whole sorted array: the cell index
    # plus the within-cell share. A query for (cell j, probability p) is then a
    # single search for `j + p`, which is what makes this one `searchsorted`
    # over `covered * len(probs)` queries rather than a loop over cells. The key
    # is strictly increasing because every region carries a positive area, so a
    # cell's last key is exactly `j + 1` and the next cell's first is above it.
    key = cs.astype(np.float64) + frac
    live = np.flatnonzero(covered)
    query = (live[:, None].astype(np.float64) + probs[None, :]).ravel()
    hit = np.searchsorted(key, query, side="left").reshape(live.size, probs.size)
    # CLAMP INTO THE CELL'S OWN BLOCK, which is the invariant and not a guard: a
    # cell's quantile is a value that cell holds. It also repairs the one
    # boundary the composite key cannot separate -- a cell's last key is exactly
    # `j + 1`, which is the same number as the NEXT cell's `p = 0` query, so
    # without it `Q(0)` returns the previous cell's maximum on every cell that
    # follows a populated one.
    starts = ends - counts
    np.clip(hit, starts[live][:, None], (ends[live] - 1)[:, None], out=hit)
    out[live] = vs[hit]
    return out, covered


def area_fraction_above(quantiles, probs, threshold):
    """Share of a cell's population above `threshold`, read off its quantiles.

    The inverse of `cell_quantiles`, and it lives beside it so that reading an
    answer out of a hypsometry table is one expression in one place rather than
    the fourth implementation of the criterion GRID-2 found three of.

    `quantiles` is `(ncell, nq)`, `probs` the vector it was taken at, and
    `threshold` a scalar or one value per cell. Per cell is the case that
    matters: a freezing height is a property of the cell's own climate.

    The cumulative curve is interpolated LINEARLY IN PROBABILITY between
    tabulated quantiles. That interpolates the CURVE and not the terrain -- the
    result is a share and never an elevation -- so nothing here manufactures
    ground the mesh does not have. What comes back is accurate to the
    probability STEP bracketing the threshold, which is a property of the
    caller's `probs` and is the bound a check against the mesh actually reaches.

    THE ENDPOINTS ARE A CONVENTION AND ONE OF THEM CAN BE READ WRONG. A
    threshold at or below the cell's minimum returns 1 and one at or above its
    maximum returns `1 - probs[-1]`, which is 0 for a vector reaching one. The
    top is exactly the share strictly above. The bottom is the share strictly
    above only where the minimum is UNIQUE, because the table records `Q(0)` as
    standing at probability zero and cannot say how much area sits at that one
    value. So a cell whose population is a single repeated value -- one mesh
    region of land in a coastal cell -- reports ALL of itself above a threshold
    equal to that value, where the exact answer is none of it. That is not an
    interpolation error and no denser `probs` removes it: such a cell has no
    distribution to read. Carry the population's region count beside the share
    and treat those cells as what they are.
    """
    q = np.asarray(quantiles, dtype=np.float64)
    p = np.asarray(probs, dtype=np.float64)
    if q.ndim != 2 or q.shape[1] != p.size:
        raise ValueError(
            f"quantiles is {q.shape} and probs has {p.size} entries; they are "
            "not the same table")
    z = np.asarray(threshold, dtype=np.float64)
    z = np.full(q.shape[0], float(z)) if z.ndim == 0 else z.reshape(-1)
    if z.size != q.shape[0]:
        raise ValueError(
            f"{z.size} thresholds for {q.shape[0]} cells; pass one per cell or "
            "a scalar")

    n = p.size
    rows = np.arange(q.shape[0])
    k = (q < z[:, None]).sum(axis=1)             # first tabulated point >= z
    j = np.clip(k, 1, n - 1)
    lo_q, hi_q = q[rows, j - 1], q[rows, j]
    lo_p, hi_p = p[j - 1], p[j]
    width = hi_q - lo_q
    step = np.where(width > 0, (z - lo_q) / np.where(width > 0, width, 1.0), 1.0)
    f = lo_p + np.clip(step, 0.0, 1.0) * (hi_p - lo_p)
    f = np.where(k == 0, 0.0, f)                 # at or below the minimum
    f = np.where(k >= n, p[-1], f)               # above the maximum
    return np.clip(1.0 - f, 0.0, 1.0)


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



# --- grids as CELL BOUNDARIES, for the crossings that are not the identity ---
#
# Everything above bins a mesh onto ONE grid, where the answer is an index. A
# crossing between two DIFFERENT grids has no index to appeal to, so it needs
# each grid's cell boundaries and their areas, and `lib/remap.py` is the
# operator that consumes them. The boundaries are constructed here because this
# module owns the column expression and a second file deriving one is the defect
# `check_one_grid_convention` exists to catch.
#
# Two facts make the boundaries constructible rather than negotiable:
#
#   A Gaussian grid HAS NO CELL BOUNDARIES. It has nodes and quadrature
#   weights, and the only construction whose cell areas reproduce those weights
#   is the one that lays the weights end to end in the sine of latitude. Taking
#   midpoints between Gaussian latitudes instead gives a partition of the
#   sphere that is perfectly self-consistent and is NOT the one the model's own
#   budget is taken over, so a crossing built on it would conserve against a
#   grid the model does not use and would do it quietly.
#
#   GOLDSTEIN's `igrid = 0` is uniform in the sine of latitude by construction:
#   `initialise_goldstein.F` sets `sv(j) = s0 + j*(s1-s0)/jmax` and
#   `asurf(j) = rsc*rsc*ds(j)*dphi`, so every ocean cell has exactly the same
#   area. `igrid = 1` is uniform in latitude and `igrid = 2` takes the
#   atmosphere's own rows, which is what the shipped PlaSim coupling
#   (`genie-main/configs/pl_go_gs_GMD.xml`) selects.
#
# Neither grid's row boundaries depend on longitude and neither grid's column
# boundaries depend on latitude, so an overlap area factorises into a longitude
# overlap times a sine-of-latitude overlap. That is why a `GridSpec` carries two
# one-dimensional edge arrays and no polygon.
#
# THE PLANETARY RADIUS IS NOT HERE. Every remap weight is a ratio of areas, so
# the radius cancels; `cell_area_fraction` returns the share of the sphere and
# `cell_area` multiplies it by a radius the CALLER reads from
# `config/planet.yaml`. `CLAUDE.md` rule 2.


@dataclass(frozen=True)
class GridSpec:
    """A grid as its cell boundaries: what an area weight needs and a centre is not.

    `lon_edges` is `nlon + 1` degrees, strictly increasing, spanning exactly
    360. `sin_edges` is `nlat + 1` values of the SINE of latitude in the grid's
    OWN row order, so row `i` is bounded by `sin_edges[i]` and `sin_edges[i+1]`
    and the array decreases for a north-to-south grid and increases for a
    south-to-north one. Both orders are real: this project's grids run north to
    south (`row()` above) and GOLDSTEIN's `j` runs south to north, and a
    remapped field comes out in the DESTINATION grid's order because that is the
    order the model that reads it wants.

    `source` records what constructed the edges. It is not decoration: the
    crossing refuses a spec that was assembled from centre labels read off a
    file, because reconstructing a coordinate from labels is what `CLAUDE.md`
    rule 3 forbids and it has already matched zero cells three times.
    """

    name: str
    lon_edges: np.ndarray
    sin_edges: np.ndarray
    source: str

    def __post_init__(self):
        lon = np.asarray(self.lon_edges, dtype=np.float64)
        sin = np.asarray(self.sin_edges, dtype=np.float64)
        if lon.ndim != 1 or lon.size < 2 or sin.ndim != 1 or sin.size < 2:
            raise ValueError(f"{self.name}: edges must be one-dimensional and hold a cell")
        if not np.all(np.diff(lon) > 0):
            raise ValueError(f"{self.name}: longitude edges must strictly increase")
        span = float(lon[-1] - lon[0])
        if abs(span - _FULL_TURN) > 1e-9:
            raise ValueError(f"{self.name}: longitude edges span {span} degrees, not a full turn")
        step = np.diff(sin)
        if not (np.all(step > 0) or np.all(step < 0)):
            raise ValueError(f"{self.name}: sine-of-latitude edges must be monotonic")
        if abs(abs(sin[0]) - 1.0) > 1e-12 or abs(abs(sin[-1]) - 1.0) > 1e-12:
            raise ValueError(f"{self.name}: sine-of-latitude edges must run pole to pole")
        object.__setattr__(self, "lon_edges", lon)
        object.__setattr__(self, "sin_edges", sin)

    @property
    def nlat(self) -> int:
        return self.sin_edges.size - 1

    @property
    def nlon(self) -> int:
        return self.lon_edges.size - 1

    @property
    def ncell(self) -> int:
        return self.nlat * self.nlon

    @property
    def shape(self) -> tuple[int, int]:
        return self.nlat, self.nlon

    @property
    def dsin(self) -> np.ndarray:
        """Sine-of-latitude extent of each row. For a Gaussian grid these ARE
        the Gauss-Legendre quadrature weights, which is the property that makes
        a crossing built on this spec conserve against the model's own budget."""
        return np.abs(np.diff(self.sin_edges))

    @property
    def dlon(self) -> np.ndarray:
        """Longitude extent of each column, in degrees."""
        return np.diff(self.lon_edges)

    def cell_area_fraction(self) -> np.ndarray:
        """Share of the sphere in each cell, shaped (nlat, nlon). Sums to one.

        A cell between two sines of latitude and two longitudes covers
        `dsin * dlon / (2 * full turn)` of the sphere, exactly, at any latitude.
        No radius, because every remap weight is a ratio of areas and the radius
        cancels out of all of them.
        """
        return np.outer(self.dsin, self.dlon) / (2.0 * _FULL_TURN)

    def cell_area(self, radius_m: float) -> np.ndarray:
        """Cell areas in square metres, for a radius the CALLER read from
        `config/planet.yaml`. Nothing here knows the planet."""
        return self.cell_area_fraction() * (4.0 * np.pi * float(radius_m) ** 2)

    def cell_centres(self):
        """(lat, lon) cell centres in degrees, for distance work and for the eye.

        The latitude is the centre of the cell in the SINE, which is the cell's
        centroid on the sphere, and for a Gaussian grid it is NOT the Gaussian
        node -- the node is a quadrature abscissa, not a centre. Nothing in a
        remap weight goes through here; `_FULL_TURN` arithmetic on a centre is
        how a coordinate gets reconstructed.
        """
        mid = 0.5 * (self.sin_edges[:-1] + self.sin_edges[1:])
        lat = np.rad2deg(np.arcsin(np.clip(mid, -1.0, 1.0)))
        lon = 0.5 * (self.lon_edges[:-1] + self.lon_edges[1:])
        return lat, lon

    def cell_centroids(self) -> np.ndarray:
        """The AREA CENTROID of each cell as a 3-vector, shaped (ncell, 3).

        The mean of the position vector over the cell, which lies INSIDE the
        sphere and is shorter than one. It is not `cell_centres` lifted onto the
        sphere and the difference is the point: a field that is linear in
        position has this vector as its exact cell mean, and no point on the
        surface has that property. A rigid rotation is such a field, which is
        what makes `lib/remap.py:apply_vector`'s acceptance test an identity
        rather than a comparison.

        Both integrals are analytic on a cell bounded by two sines of latitude
        and two longitudes, so nothing here is quadrature. The polar axis is z,
        matching `cell_centres`; `source/README.md`'s y-up warning is about the
        EXPORT's mesh coordinates and does not reach a grid constructed here.
        """
        lam = np.deg2rad(self.lon_edges)
        s = self.sin_edges
        dlam = np.diff(lam)
        ds = np.diff(s)
        area = np.outer(ds, dlam)
        root = np.sqrt(np.clip(1.0 - s * s, 0.0, None))
        prim = 0.5 * (s * root + np.arcsin(np.clip(s, -1.0, 1.0)))
        band = np.diff(prim)                       # integral of the cosine over the row
        ix = np.outer(band, np.diff(np.sin(lam)))
        iy = np.outer(band, -np.diff(np.cos(lam)))
        iz = np.outer(0.5 * np.diff(s * s), dlam)
        return np.stack([(ix / area).ravel(), (iy / area).ravel(),
                         (iz / area).ravel()], axis=1)

    def cell_frames(self):
        """Local (east, north, up) unit vectors per cell, each shaped (ncell, 3).

        The frame is taken at the cell's own area centroid projected onto the
        sphere, so it is the one reference point a cell has that is defined by
        the cell rather than by a convention. This is what a VECTOR field's
        components are components IN, and it is why `lib/remap.py:apply_vector`
        can move a vector between two grids without ever asking one grid what
        the other calls east: it converts to the sphere's own three cartesian
        components and back.
        """
        up = self.cell_centroids()
        up = up / np.linalg.norm(up, axis=1, keepdims=True)
        polar = np.array([0.0, 0.0, 1.0])
        east = np.cross(polar, up)
        n = np.linalg.norm(east, axis=1, keepdims=True)
        # At a pole the cross product vanishes and east is undefined. No cell
        # centroid reaches a pole -- a polar cell's centroid sits inside its own
        # sine band -- so this is a guard and not a case.
        east = np.divide(east, n, out=np.zeros_like(east), where=(n > 0))
        north = np.cross(up, east)
        return east, north, up


def gaussian_grid(nlat: int, nlon: int | None = None, name: str | None = None) -> GridSpec:
    """The spectral grid ExoPlaSim and the export share, as cell boundaries.

    Rows run north to south, matching `row()` and the export's own
    `rowOrder`. The row boundaries are the Gauss-Legendre weights laid end to
    end from the north pole, so `dsin` reproduces the weights and the cell areas
    reproduce the quadrature the model's global budget is taken over. Columns
    are `column()`'s: cell `k` spans `[-180 + k*dlon, -180 + (k+1)*dlon)`.

    `nlon` defaults to `2 * nlat`, which is the ladder's shape; passing it
    explicitly is for a grid that is Gaussian but not on the ladder.
    """
    nlat = int(nlat)
    nlon = 2 * nlat if nlon is None else int(nlon)
    _, weights = np.polynomial.legendre.leggauss(nlat)
    sin_edges = 1.0 - np.concatenate([[0.0], np.cumsum(weights[::-1])])
    sin_edges[0], sin_edges[-1] = 1.0, -1.0
    lon_edges = -180.0 + np.arange(nlon + 1, dtype=np.float64) * (360.0 / nlon)
    return GridSpec(name=name or f"gaussian-{nlat}x{nlon}",
                    lon_edges=lon_edges, sin_edges=sin_edges,
                    source="lib/gridding.py:gaussian_grid, Gauss-Legendre weights")


def gaussian_latitudes(nlat: int) -> np.ndarray:
    """The Gaussian NODES, north to south, in degrees. Not cell centres.

    This is the axis the export writes and the model transforms on. It is here
    so that a caller checking a spec against an artifact has one expression to
    check against rather than a second construction of its own.
    """
    nodes, _ = np.polynomial.legendre.leggauss(int(nlat))
    return np.rad2deg(np.arcsin(nodes))[::-1]


def require_gaussian_rows(spec: GridSpec, lat_axis, what: str = "the grid") -> None:
    """Refuse a Gaussian spec whose rows are not the axis the artifact carries.

    FOR AN AXIS STORED AS FLOAT64, which is an export's. A CLIMATOLOGY's axis is
    netCDF float32 and agrees with the constructed nodes to a few parts in a
    million and no closer, so this bar refuses a correct axis there;
    `require_same_rows` is the door for that comparison and carries the looser
    tolerance for exactly this reason.

    The spec is CONSTRUCTED and the artifact's axis is READ, and this is the one
    place they are compared. If they disagree the spec is describing a different
    grid from the field it is about to remap, and every weight is then a weight
    between two grids nobody is using. Measured against the export's own axis
    the two agree to about 1e-14 degrees, which is the same expression evaluated
    twice; the bar below is loose against that and tight against any real
    disagreement.
    """
    axis = np.asarray(lat_axis, dtype=np.float64)
    want = gaussian_latitudes(spec.nlat)
    if axis.size != want.size:
        raise SystemExit(
            f"{what}: the spec has {want.size} rows and the axis has {axis.size}. "
            "They are not the same grid and no mapping between them is defined.")
    off = float(np.abs(axis - want).max())
    if off > 1e-6:
        raise SystemExit(
            f"{what}: the constructed Gaussian rows and the axis on disk differ "
            f"by up to {off:.3g} degrees. Do NOT reconcile them by matching "
            "nearest centres; construct both from one source. CLAUDE.md rule 3.")


# --- THE AREA WEIGHT, and the one thing it is not ---------------------------
#
# A weighted mean over a Gaussian grid takes the QUADRATURE WEIGHTS, which are
# the sine-of-latitude extent of each row and are what `GridSpec.dsin` is. The
# two accessors below are that weight with the axis checked, so that a caller
# needing a global or land mean has one door and does not reach for the cosine.
#
# `cos(lat)` IS NOT AN AREA WEIGHT HERE. It is the metric factor of a band that
# is equally spaced in latitude, and a Gaussian row is not one: the rows crowd
# toward the equator and the nodes are quadrature abscissae rather than cell
# centres. Against the quadrature weights the cosine is wide in the polar row by
# 1.767 per cent at T21, 1.796 at T42, 1.805 at T85 and 1.808 at T170, which
# converges to about 1.81 rather than to zero; the interior converges normally,
# 0.0292 per cent median at T21 down to 0.0009 at T170. So a cosine-weighted
# mean is a mean over a grid nobody is running, and it is quiet, because a
# weight that is wrong by a part in ten thousand returns an ordinary-looking
# number everywhere. `notes/audits/ocean-grid-crossing.md` section 3 has the
# measurement and `source/README.md` names the cosine as a fourth thing that is
# none of the three area quantities an export ships.
#
# A cosine that is a METRIC term is a different creature and stays: the
# east-west cell width an advection scheme divides by, the polar-convergence
# blur a projection applies, the unit vector a KD-tree is built on. What is
# replaced here is the cosine used as a WEIGHT.


def _require_gaussian_axis(lat_axis, what: str) -> int:
    """The row count of an axis that IS this many Gaussian nodes, or a refusal.

    The bar is `require_same_rows`'s rather than `require_gaussian_rows`'s,
    because a climatology's latitude axis is netCDF float32 and agrees with the
    constructed nodes to a few parts in a million and no closer, and because a
    per-cell text product prints its rows to a handful of decimals. It is still
    two orders tighter than the gap between a Gaussian axis and an equally
    spaced one at any rung on the ladder, which is what it exists to catch.
    """
    axis = np.asarray(lat_axis, dtype=np.float64)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError(
            f"{what}: a latitude axis is one-dimensional and holds a row; this "
            f"has shape {axis.shape}")
    off = float(np.abs(axis - gaussian_latitudes(axis.size)).max())
    if off > 1e-3:
        raise SystemExit(
            f"{what}: the axis differs from the {axis.size}-row Gauss-Legendre "
            f"nodes by up to {off:.4g} degrees, so it is not the grid these "
            "weights are the quadrature of. Weight the field on the partition "
            "its own grid has rather than on one it does not.")
    return int(axis.size)


def gaussian_row_weights(lat_axis, what: str = "the area weight") -> np.ndarray:
    """Area weight of each ROW of the Gaussian grid this axis names.

    The Gauss-Legendre quadrature weights north to south, summing to two. A
    cell's share of the sphere is one of these divided by `2 * nlon`, so on a
    grid whose columns are all one width a mean over rows needs nothing more,
    and the constant cancels out of any ratio of weighted sums.

    Take this where the weight is broadcast against an axis the caller already
    holds -- a spectrum over wavenumber, a column over levels -- and
    `gaussian_area_weights` where a two-dimensional field is being averaged.
    """
    return gaussian_grid(_require_gaussian_axis(lat_axis, what)).dsin


def gaussian_area_weights(lat_axis, nlon: int,
                          what: str = "the area weight") -> np.ndarray:
    """Share of the sphere in every cell, shaped `(nlat, nlon)` and summing to one.

    So `(field * w).sum()` is a global mean, `(field * w)[mask].sum() /
    w[mask].sum()` is a mean over a subset, and multiplying by `4 pi R^2` for a
    radius the CALLER reads from `config/planet.yaml` turns a weighted sum into
    an integral. `CLAUDE.md` rule 2: nothing here knows the planet.
    """
    return gaussian_grid(_require_gaussian_axis(lat_axis, what),
                         int(nlon)).cell_area_fraction()


# ---------------------------------------------------------------------------
# THE MODEL'S OWN LABELS. `column()` above places an Orogen longitude, which
# runs -180..180; ExoPlaSim labels the SAME columns 0..360 from index 0 and has
# never heard of Orogen. Both label column k of one grid, and the mapping
# between them is the identity ON THE INDEX and a half-grid rotation on the
# label -- which is CLAUDE.md rule 3 and what
# notes/audits/grid-convention-and-runoff.md measured at 0.5955 of cells
# agreeing when the labels are matched instead of the indices.
#
# These exist because artifacts DO carry the model's labels: a climatology's
# `lon` variable, and every per-cell text product written from one --
# `pedology/data/<build>/land_column_states_<res>.txt` keys its rows that way.
# Reading such a product needs the label axis, and a consumer that builds one of
# its own is the third occurrence of the defect. So the axis is constructed once
# here, the inverse is a function rather than a subtraction at the point of use,
# and a label that is not ON the axis is refused rather than snapped to the
# nearest column.
# ---------------------------------------------------------------------------

def model_longitude_labels(nlon: int) -> np.ndarray:
    """The longitude ExoPlaSim WRITES for grid column i. A LABEL, not a coordinate.

    Degrees east from index 0, which is the column `column()` places Orogen
    longitude -180 in. Nothing computed may be done with these against Orogen
    coordinates; they exist to be matched against an artifact that carries them.
    """
    nlon = int(nlon)
    return np.arange(nlon, dtype=np.float64) * (360.0 / nlon)


def model_label_column(lon_label, nlon: int, what: str = "a longitude label",
                       tolerance_cells: float = 1.0e-3) -> np.ndarray:
    """Grid column for a longitude on ExoPlaSim's own label axis.

    NOT `column()`, and the two are not interchangeable: `column()` reads an
    Orogen longitude and this reads a model label, and on the same grid they
    differ by half the columns.

    A label is a name for a column, so this REFUSES one that does not sit on a
    column rather than rounding it to the nearest. That is the check that
    catches an axis handed in on the wrong convention: an Orogen centre offered
    here lands exactly half a column from every model label and fails by 500
    times the bar, where a nearest-column read would have accepted it and
    rotated the field by half a world.
    """
    nlon = int(nlon)
    exact = np.asarray(lon_label, dtype=np.float64) % 360.0 * nlon / 360.0
    col = np.rint(exact)
    off = float(np.abs(exact - col).max()) if exact.size else 0.0
    if off > tolerance_cells:
        raise SystemExit(
            f"{what}: a label is {off:.4g} of a column away from any column of "
            f"the {nlon}-column model label axis. It is not on this axis. Half "
            "a column means it is on the export's -180..180 centres instead, "
            "which is a DIFFERENT NAME FOR THE SAME COLUMNS and not a "
            "coordinate transform. CLAUDE.md rule 3.")
    return col.astype(np.int64) % nlon


def model_label_cells(lat_label, lon_label, nlat: int, nlon: int,
                      what: str = "a labelled point"):
    """(row, column) on the model grid for points carrying the MODEL's labels.

    The one entry point for "which cell is this row of a per-cell model-labelled
    product". Rows go through `row()` against the constructed Gaussian nodes,
    which bins on the midpoints and so is indifferent to the few parts in a
    million a float32 axis printed to four decimals costs; columns go through
    `model_label_column`, which refuses a label off the axis.
    """
    rows = row(lat_label, gaussian_latitudes(int(nlat)))
    cols = model_label_column(lon_label, nlon, what=what)
    return rows, cols


def label_row_weights(lat_label, nlat: int, what: str = "a labelled point",
                      tolerance_rows: float = 0.05) -> np.ndarray:
    """Row area weight for each POINT of a per-cell product keyed by latitude.

    For an artifact whose rows are addressed by the latitude they carry rather
    than by an index -- an LPJ-GUESS `.out` table, `pedology/data/<build>/
    land_column_states_<res>.txt` -- where the rows present are whatever subset
    was simulated, so there is no axis on the file to hand `gaussian_row_weights`.
    Each point gets the quadrature weight of the row its latitude names, which
    is proportional to that cell's area because every column is one width, and
    the constant cancels out of a weighted mean.

    THE ROW IS LOOKED UP AND THEN CHECKED. `row()` bins on the midpoints and so
    always returns a row; what makes this refuse rather than snap is that the
    latitude is then compared against the node of the row it landed in, which is
    how a table written on a different rung from the one asked for here is
    caught instead of being weighted by whichever rows its latitudes happen to
    fall in.

    THE BAR IS IN ROWS AND NOT IN DEGREES, because these products print their
    coordinates to a couple of decimals: an LPJ-GUESS table on the T21 rows is
    up to 0.005 degrees off the nodes it was written from, which is a thousandth
    of a row and would fail any absolute bar tight enough to be worth having.
    A table on another rung's rows misses by a real fraction of a row, so the
    ratio is the quantity that separates the two.
    """
    nodes = gaussian_latitudes(int(nlat))
    gaps = np.abs(np.diff(nodes))
    spacing = np.empty(nodes.size)
    spacing[0], spacing[-1] = gaps[0], gaps[-1]
    spacing[1:-1] = np.minimum(gaps[:-1], gaps[1:])
    labels = np.asarray(lat_label, dtype=np.float64)
    rows = row(labels, nodes)
    off = (float((np.abs(labels - nodes[rows]) / spacing[rows]).max())
           if labels.size else 0.0)
    if off > tolerance_rows:
        raise SystemExit(
            f"{what}: a latitude is {off:.4g} of a row from the nearest row of "
            f"the {int(nlat)}-row Gaussian grid. It is not on this grid, and "
            "the weight it would be given belongs to a row it does not sit in.")
    return gaussian_grid(int(nlat)).dsin[rows]


def require_model_labels(lat_axis, lon_axis, nlat: int, nlon: int,
                         what: str = "the axes") -> None:
    """Refuse an axis pair that is not the labels ExoPlaSim writes for this grid.

    For an artifact the model wrote: a climatology, or a product derived from
    one. `require_same_rows` carries the latitude comparison because such an
    axis is netCDF float32 and agrees with the constructed nodes to a few parts
    in a million and no closer.

    This is what says a file handed to a builder is on the grid the builder
    thinks it is, in the one direction that matters: a file on another rung has
    the wrong number of rows, and a file on the export's own longitude labels is
    half a grid out.
    """
    require_same_rows(np.asarray(lat_axis, dtype=np.float64),
                      gaussian_latitudes(int(nlat)), what=what)
    axis = np.asarray(lon_axis, dtype=np.float64)
    want = model_longitude_labels(nlon)
    if axis.size != want.size:
        raise SystemExit(
            f"{what}: the longitude axis has {axis.size} labels and this grid "
            f"has {want.size} columns. They are not the same grid.")
    off = float(np.abs(axis - want).max())
    if off > 1.0e-3:
        raise SystemExit(
            f"{what}: the longitude axis differs from the labels ExoPlaSim "
            f"writes for a {nlon}-column grid by up to {off:.4g} degrees. Do "
            "NOT reconcile the two by matching labels; they name the same "
            "columns and the mapping is the index. CLAUDE.md rule 3.")


GOLDSTEIN_EQUAL_AREA = 0        # igrid = 0: uniform in the sine of latitude
GOLDSTEIN_EQUAL_ANGLE = 1       # igrid = 1: uniform in latitude
GOLDSTEIN_ATMOSPHERE_ROWS = 2   # igrid = 2: the atmosphere's own rows
GOLDSTEIN_LON_ORIGIN = -260.0   # phi0 in initialise_goldstein.F, degrees


def goldstein_grid(nlon: int, nlat: int, igrid: int = GOLDSTEIN_EQUAL_AREA,
                   lon_origin_deg: float = GOLDSTEIN_LON_ORIGIN,
                   atmosphere_rows: np.ndarray | None = None,
                   name: str | None = None) -> GridSpec:
    """The GOLDSTEIN ocean grid as cell boundaries, read off the model's own setup.

    `nlon` is `imax` and `nlat` is `jmax`, and rows run SOUTH TO NORTH because
    `j` does. The three `igrid` arms are the three the model implements and the
    arm is a required fact about the ocean configuration rather than a default
    to be inherited quietly: `igrid = 0` makes every ocean cell exactly equal in
    area, and it is what the shipped 36 x 36 configurations use, while the
    shipped PlaSim coupling uses `igrid = 2`.

    `lon_origin_deg` is `phi0`, whose shipped value puts the first column edge
    at 260 degrees west. It is a parameter and not a label: the crossing takes
    both grids' columns from their constructors and never compares longitude
    numbers between them, which is the form `CLAUDE.md` rule 3 takes here.
    """
    nlon, nlat, igrid = int(nlon), int(nlat), int(igrid)
    if igrid == GOLDSTEIN_EQUAL_AREA:
        sin_edges = -1.0 + np.arange(nlat + 1, dtype=np.float64) * (2.0 / nlat)
        sin_edges[-1] = 1.0
        arm = "igrid=0, uniform in the sine of latitude"
    elif igrid == GOLDSTEIN_EQUAL_ANGLE:
        lat_edges = -90.0 + np.arange(nlat + 1, dtype=np.float64) * (180.0 / nlat)
        sin_edges = np.sin(np.deg2rad(lat_edges))
        sin_edges[0], sin_edges[-1] = -1.0, 1.0
        arm = "igrid=1, uniform in latitude"
    elif igrid == GOLDSTEIN_ATMOSPHERE_ROWS:
        if atmosphere_rows is None:
            raise ValueError(
                "igrid=2 takes the atmosphere's own row boundaries; pass "
                "atmosphere_rows=gaussian_grid(nlat).sin_edges rather than "
                "letting this module guess which atmosphere is meant")
        sin_edges = np.sort(np.asarray(atmosphere_rows, dtype=np.float64))
        if sin_edges.size != nlat + 1:
            raise ValueError(
                f"igrid=2 with jmax={nlat} needs {nlat + 1} row boundaries and "
                f"was given {sin_edges.size}")
        arm = "igrid=2, the atmosphere's rows"
    else:
        raise ValueError(
            f"igrid={igrid} is not one GOLDSTEIN implements. "
            "initialise_goldstein.F has 0, 1 and 2 and nothing else.")
    lon_edges = float(lon_origin_deg) + np.arange(nlon + 1, dtype=np.float64) * (360.0 / nlon)
    return GridSpec(name=name or f"goldstein-{nlon}x{nlat}",
                    lon_edges=lon_edges, sin_edges=sin_edges,
                    source=f"lib/gridding.py:goldstein_grid, {arm}, phi0={lon_origin_deg}")


# --- the operators' own checks. SPAT-4 --------------------------------------
#
# `python lib/gridding.py --selftest`
#
# Every one of these has a RIGHT ANSWER rather than a plausible one: an
# identity, a conservation law, or a control that must fail. That is the bar
# `docs/src/practice/failure-modes.md` class 17 sets, and it is why these run on
# synthetic input and need no build on disk -- a check that can only be run
# where the answer is unknown is not a check.


_CHECKS = 44


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

    # DISTRIBUTION. The ends of a quantile table are the true extremes, exactly,
    # at any area weighting -- that is the identity a convention that
    # interpolates between regions fails, and interpolating is how ground the
    # mesh does not have gets manufactured.
    probs = np.array([0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0])
    q, covered = cell_quantiles(cell, ncell, area, values, probs, population)
    lo = np.full(ncell, np.inf)
    hi = np.full(ncell, -np.inf)
    np.minimum.at(lo, cell[population], values[population])
    np.maximum.at(hi, cell[population], values[population])
    check("a quantile table's ends are the cell's true extremes",
          bool(np.array_equal(q[covered, 0], lo[covered])
               and np.array_equal(q[covered, -1], hi[covered])),
          f"{int((q[covered, 0] != lo[covered]).sum())} minima and "
          f"{int((q[covered, -1] != hi[covered]).sum())} maxima differ")
    check("a quantile table is non-decreasing in probability",
          bool((np.diff(q[covered], axis=1) >= 0).all()),
          f"{int((np.diff(q[covered], axis=1) < 0).sum())} inversions")

    # Every tabulated value is a value some region in that cell actually has.
    k = int(np.flatnonzero(covered & (np.bincount(cell[population], minlength=ncell) > 50))[0])
    have = set(values[population & (cell == k)].tolist())
    check("every tabulated quantile is a value the cell holds",
          all(v in have for v in q[k].tolist()),
          f"{sum(v not in have for v in q[k].tolist())} of {probs.size} are not")

    # The median of an EQUALLY weighted cell is the ordinary sample median under
    # this convention, which is what says the area weighting is the only thing
    # separating the two.
    flat = np.ones(n)
    qf, cf = cell_quantiles(cell, ncell, flat, values, np.array([0.0, 0.5, 1.0]),
                            population)
    sel = population & (cell == k)
    v = np.sort(values[sel])
    want = v[int(np.ceil(0.5 * v.size)) - 1]
    check("with equal weights the median is the ordinary sample median",
          qf[k, 1] == want, f"{qf[k, 1]} against {want}")

    # THE INVERSE, and it is exact where it has to be: a threshold at the cell
    # minimum leaves all of the population above it and one at the maximum
    # leaves none.
    at_min = area_fraction_above(q[covered], probs, lo[covered])
    at_max = area_fraction_above(q[covered], probs, hi[covered])
    check("the inverse returns all above the minimum and none above the maximum",
          bool(np.allclose(at_min, 1.0, atol=0) and np.allclose(at_max, 0.0, atol=0)),
          f"min {at_min.min():.3g}..{at_min.max():.3g}, "
          f"max {at_max.min():.3g}..{at_max.max():.3g}")

    # THE NEGATIVE CONTROL for the pair. Read the share above the cell's own
    # TABULATED 0.75 quantile back out of the table: it has to come out at 0.25.
    # A table read with the wrong convention, or an inverse that interpolated in
    # the wrong variable, returns an ordinary-looking number here and not this
    # one.
    # Restricted to cells whose table has no repeated value, because a repeated
    # one is genuinely ambiguous about which probability it stands at and that
    # is a property of the definition rather than a defect in the inverse.
    strict = np.zeros(ncell, dtype=bool)
    strict[covered] = (np.diff(q[covered], axis=1) > 0).all(axis=1)
    got = area_fraction_above(q[strict], probs, q[strict, probs.searchsorted(0.75)])
    worst = float(np.abs(got - 0.25).max()) if strict.any() else float("nan")
    check("the inverse recovers the probability a tabulated quantile stands at",
          bool(strict.any()) and worst <= 1e-12,
          f"{int(strict.sum())} tie-free cells, worst {worst:.3g}")

    # THE ENDPOINT CONVENTION, pinned rather than only documented. A cell whose
    # population is one repeated value has no distribution, and the table cannot
    # say how much area stands AT its own Q(0), so the inverse returns all of it
    # above a threshold equal to that value where the exact answer is none. A
    # consumer that does not know this reads a coastal cell holding one mesh
    # region as entirely above every height.
    one = np.zeros(12, dtype=np.int64)
    flat_v = np.full(12, 1.75)
    tbl, _ = cell_quantiles(one, 1, np.ones(12), flat_v, probs)
    check("a single-valued cell reports all of itself above its own value",
          float(area_fraction_above(tbl, probs, 1.75)[0]) == 1.0
          and float(area_fraction_above(tbl, probs, 1.7500001)[0]) == 0.0,
          f"at the value {area_fraction_above(tbl, probs, 1.75)[0]}, just above "
          f"{area_fraction_above(tbl, probs, 1.7500001)[0]}")

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

    # GRID SPECS. A crossing between two grids is only as good as the two
    # partitions of the sphere it is built from, and each of these has an
    # answer rather than a range. `analysis/ocean_remap.py` carries the rest of
    # the suite, including the negative control, because it needs the export.
    gauss = gaussian_grid(64)
    check("a Gaussian spec's cells partition the sphere",
          abs(gauss.cell_area_fraction().sum() - 1.0) <= 1e-13,
          f"{gauss.cell_area_fraction().sum() - 1.0:.3g}")
    _, gl = np.polynomial.legendre.leggauss(64)
    check("a Gaussian spec's rows ARE the quadrature weights",
          float(np.abs(gauss.dsin - gl[::-1]).max()) <= 1e-13,
          f"{np.abs(gauss.dsin - gl[::-1]).max():.3g}")
    ocean = goldstein_grid(36, 36)
    spread = float(np.ptp(ocean.cell_area_fraction()))
    check("GOLDSTEIN igrid=0 makes every cell exactly equal in area",
          spread <= 1e-16, f"spread {spread:.3g}")
    equal_angle = goldstein_grid(36, 36, igrid=GOLDSTEIN_EQUAL_ANGLE)
    lat_edges = np.rad2deg(np.arcsin(equal_angle.sin_edges))
    check("GOLDSTEIN igrid=1 makes every row exactly equal in latitude",
          float(np.ptp(np.diff(lat_edges))) <= 1e-12
          and abs(equal_angle.cell_area_fraction().sum() - 1.0) <= 1e-13,
          f"row spread {np.ptp(np.diff(lat_edges)):.3g}")

    # CENTROIDS AND FRAMES, which is what a VECTOR field's components are
    # components in. The whole sphere as ONE cell has its centroid exactly at
    # the origin, analytically and by symmetry, and it is the one case where the
    # answer is a number rather than a limit.
    whole = GridSpec(name="the whole sphere", lon_edges=np.array([-180.0, 180.0]),
                     sin_edges=np.array([1.0, -1.0]),
                     source="the selftest's one-cell grid")
    check("the sphere as one cell has its centroid at the origin",
          float(np.abs(whole.cell_centroids()).max()) <= 1e-15,
          f"{np.abs(whole.cell_centroids()).max():.3g}")
    coarse = float(np.linalg.norm(gaussian_grid(32).cell_centroids(), axis=1).min())
    fine = float(np.linalg.norm(gaussian_grid(256).cell_centroids(), axis=1).min())
    check("a centroid lies inside the sphere and approaches it as cells shrink",
          coarse < fine < 1.0, f"{coarse:.6f} at T21 against {fine:.6f} at T170")
    fe, fn, fu = gaussian_grid(64).cell_frames()
    ortho = max(float(np.abs((fe * fn).sum(axis=1)).max()),
                float(np.abs((fe * fu).sum(axis=1)).max()),
                float(np.abs((fn * fu).sum(axis=1)).max()),
                float(np.abs(np.linalg.norm(fe, axis=1) - 1.0).max()),
                float(np.abs(np.linalg.norm(fn, axis=1) - 1.0).max()))
    check("the local frame is orthonormal in every cell",
          ortho <= 1e-14 and bool((fn[:, 2] > 0).all()),
          f"worst departure {ortho:.3g}")

    # BINNING ONTO A CONSTRUCTED GRID. `spec_cells` is the mesh-to-grid door for
    # a grid the exporter does not write, and the identity it has to satisfy is
    # that a cell's own CENTRE lands back in that cell -- on both grid kinds and
    # in both row orders, one running north to south and the other south to
    # north -- and that the same centre offered a whole turn either way lands
    # there too, because a longitude is periodic and the fold is exact. It fails
    # if the sine is binned as an angle, if the row order is assumed, or if the
    # longitude origin is read as a label.
    for spec in (gaussian_grid(64), goldstein_grid(36, 36),
                 goldstein_grid(72, 72, igrid=GOLDSTEIN_EQUAL_ANGLE)):
        clat, clon = spec.cell_centres()
        own = np.arange(spec.ncell)
        la = np.repeat(clat, spec.nlon)
        lo = np.tile(clon, spec.nlat)
        wrong = max(int((spec_cells(spec, la, lo + turn) != own).sum())
                    for turn in (-360.0, 0.0, 360.0))
        check(f"{spec.name}: every cell centre bins back into its own cell",
              wrong == 0, f"{wrong} of {spec.ncell} misplaced")

    # CONTAINMENT, against the constructor's OWN edge arrays rather than against
    # a restatement of them: whatever cell a point is placed in, the point lies
    # between that cell's boundaries. This is the property `region_cells` gets
    # from the export's axis and a constructed grid has to be given.
    ocn = goldstein_grid(36, 36)
    m = 50000
    rand_lat = np.rad2deg(np.arcsin(rng.uniform(-1.0, 1.0, m)))
    rand_lon = rng.uniform(-180.0, 180.0, m)
    binned = spec_cells(ocn, rand_lat, rand_lon)
    br, bc = np.divmod(binned, ocn.nlon)
    sn = np.sin(np.deg2rad(rand_lat))
    turn = float(ocn.lon_edges[-1] - ocn.lon_edges[0])
    folded = rand_lon - turn * np.floor((rand_lon - ocn.lon_edges[0]) / turn)
    lo_s = np.minimum(ocn.sin_edges[br], ocn.sin_edges[br + 1])
    hi_s = np.maximum(ocn.sin_edges[br], ocn.sin_edges[br + 1])
    outside = int(((sn < lo_s - 1e-12) | (sn > hi_s + 1e-12)
                   | (folded < ocn.lon_edges[bc] - 1e-12)
                   | (folded > ocn.lon_edges[bc + 1] + 1e-12)).sum())
    check("a placed point lies between the boundaries of the cell it landed in",
          outside == 0 and int(binned.min()) >= 0 and int(binned.max()) < ocn.ncell,
          f"{outside} of {m} outside, range {binned.min()}..{binned.max()}")

    # THE CONTROL, and it is the defect the sine binning exists to stop. On an
    # equal-area grid the rows are uniform in the SINE, so a restatement that
    # bins the same points uniformly in the ANGLE is a real placement onto a
    # real grid and it is not this one. It has to disagree on a large share of
    # the sphere; if it did not, the sine would not be load-bearing.
    angle_row = np.clip(((rand_lat + 90.0) / (180.0 / ocn.nlat)).astype(int),
                        0, ocn.nlat - 1)
    disagree = float((angle_row != br).mean())
    check("binning in the angle misplaces points on an equal-area grid",
          disagree >= 0.2, f"only {disagree:.3%} of points move")

    # THE LONGITUDE ORIGIN IS A PARAMETER AND NOT A LABEL. GOLDSTEIN's first
    # column edge sits at phi0, whose shipped value is 260 degrees west, so the
    # SAME physical longitude is column 0 on the ocean's grid and something else
    # on the atmosphere's. Reading either grid's column number as if it were the
    # other's is `CLAUDE.md` rule 3 at this boundary.
    atm_spec = gaussian_grid(32)
    just_in = float(ocn.lon_edges[0]) + 1e-9
    just_out = float(ocn.lon_edges[0]) - 1e-9
    first = int(spec_cells(ocn, np.array([0.0]), np.array([just_in]))[0]) % ocn.nlon
    last = int(spec_cells(ocn, np.array([0.0]), np.array([just_out]))[0]) % ocn.nlon
    same = int(spec_cells(atm_spec, np.array([0.0]), np.array([just_in]))[0]) % atm_spec.nlon
    check("the longitude origin is honoured and the two grids disagree on the column",
          first == 0 and last == ocn.nlon - 1 and same != 0,
          f"phi0 column {first}, just outside {last}, atmosphere column {same}")

    # THE LEDGER CLOSES on a constructed grid too: every point offered lands in
    # exactly one cell of it, and nothing falls off the edge of the partition.
    led = transfer_ledger(binned, ocn.ncell, rng.uniform(0.1, 10.0, m))
    check("binning onto a constructed grid loses no area",
          led["closure_residual_relative"] <= 1e-12,
          f"{led['closure_residual_relative']:.3g}")

    # AND THE REFUSAL: a centre array read off a file is a label and not a
    # partition, so it is refused rather than turned into edges here.
    try:
        spec_cells(ocn.cell_centres()[0], rand_lat, rand_lon)
    except TypeError:
        refused_centres = True
    else:
        refused_centres = False
    check("a centre array offered as a grid is refused", refused_centres,
          "spec_cells accepted something that is not a GridSpec")

    # THE CONTROL: a spec that is not a partition of the sphere must be refused
    # at construction. Without this the four checks above would pass on an
    # object that silently describes a different world from the one it names.
    refused = 0
    for bad in (dict(lon_edges=np.linspace(-180.0, 90.0, 5),
                     sin_edges=np.linspace(1.0, -1.0, 5)),
                dict(lon_edges=np.linspace(-180.0, 180.0, 5),
                     sin_edges=np.array([1.0, 0.5, 0.7, -1.0]))):
        try:
            GridSpec(name="control", source="the selftest's control", **bad)
        except ValueError:
            refused += 1
    check("a spec that is not a partition of the sphere is refused",
          refused == 2, f"{refused} of 2 refused")

    # THE MODEL'S LABEL AXIS: an identity. Every label this constructs names the
    # column it was constructed from, so the round trip is the index it started
    # at and nothing about the export enters.
    for nlon_lab in (64, 128, 170):
        got = model_label_column(model_longitude_labels(nlon_lab), nlon_lab)
        check(f"model labels round-trip to their own columns at nlon={nlon_lab}",
              bool(np.array_equal(got, np.arange(nlon_lab))),
              f"{int((got != np.arange(nlon_lab)).sum())} columns wrong")

    # THE CONTROL, and it is the defect this pair exists to stop. The export's
    # own column centres are a valid longitude axis for the SAME grid, half a
    # column from every model label. Offered as model labels they must be
    # refused, not snapped to the nearest column and silently rotated.
    _, export_lon = gaussian_grid(32, 64).cell_centres()
    refused_labels = 0
    for bad_lat, bad_lon in ((gaussian_latitudes(32), export_lon),
                             (gaussian_latitudes(64), model_longitude_labels(64))):
        try:
            require_model_labels(bad_lat, bad_lon, 32, 64, what="the control")
        except SystemExit:
            refused_labels += 1
    try:
        model_label_column(export_lon, 64, what="the control")
    except SystemExit:
        refused_labels += 1
    check("an axis on the export's labels is refused as a model label axis",
          refused_labels == 3, f"{refused_labels} of 3 refused")

    # And the axis the model actually writes is accepted, so the control above
    # is refusing the convention rather than refusing everything.
    require_model_labels(gaussian_latitudes(32), model_longitude_labels(64),
                         32, 64, what="the model's own axes")
    check("the model's own axes are accepted", True)

    # THE AREA WEIGHT: a Gauss-Legendre quadrature is an IDENTITY rather than an
    # approximation, so this is asserted exactly and not compared loosely. The
    # weights integrate every polynomial in the sine of latitude up to degree
    # 2n-1 over the sphere exactly, and the moments below are three such
    # integrals with a known closed form: 2, 2/3 and 2/5. Nothing about the new
    # accessors can be right if these are not.
    # The rungs come from the ladder registry rather than being spelled here;
    # `lib/rungs.py` is the one table and a second copy of it is a lint failure.
    for nlat_w in sorted(set(RUNGS.values()))[:3]:
        nodes = gaussian_latitudes(nlat_w)
        w = gaussian_row_weights(nodes, what="the selftest's axis")
        mu = np.sin(np.deg2rad(nodes))
        moments = [float((w * mu ** k).sum()) for k in (0, 2, 4)]
        exact = [2.0, 2.0 / 3.0, 2.0 / 5.0]
        worst = max(abs(g - e) / e for g, e in zip(moments, exact))
        # 1e-12 is the round trip through degrees and back, not slack: over
        # these three rungs the quadrature misses by at most 1.6e-14 and the
        # cosine below misses by at least 3.4e-04, so the bar discriminates by
        # ten orders.
        check(f"the row weights integrate the sphere's moments exactly at "
              f"nlat={nlat_w}", worst <= 1e-12, f"worst relative miss {worst:.3g}")

    # THE CONTROL, and it is the defect this row exists to remove: the same
    # three moments taken with cos(lat) as the weight must MISS. If they did
    # not, the cosine would be the quadrature and there would be nothing here to
    # fix. The zeroth moment is the sphere itself, so the miss is reported on
    # the weight sum a caller normalises by.
    cos_w = np.cos(np.deg2rad(gaussian_latitudes(64)))
    cos_w = cos_w * (2.0 / cos_w.sum())
    mu = np.sin(np.deg2rad(gaussian_latitudes(64)))
    cos_miss = abs(float((cos_w * mu ** 2).sum()) / (2.0 / 3.0) - 1.0)
    check("cos(lat) as a weight misses the same identity",
          cos_miss > 1e-6, f"missed by only {cos_miss:.3g}")

    # And the area weight is the row weight spread over the columns: it sums to
    # one, so a weighted sum is a mean, and it agrees with the spec it is built
    # from rather than being a second construction.
    aw = gaussian_area_weights(gaussian_latitudes(48), 96,
                               what="the selftest's axis")
    check("the area weights partition the sphere",
          abs(float(aw.sum()) - 1.0) <= 1e-14
          and float(np.abs(aw - gaussian_grid(48, 96).cell_area_fraction()).max()) <= 0.0,
          f"sum {aw.sum() - 1.0:.3g}")

    # THE CONTROL for the axis check: an equally spaced latitude axis of the
    # same length is a real grid and is not this one, so it must be refused
    # rather than weighted as if it were Gaussian. Same for a latitude label
    # that sits between two rows.
    refused_axes = 0
    uniform = 90.0 - (np.arange(64) + 0.5) * (180.0 / 64)
    for bad in (lambda: gaussian_row_weights(uniform, what="the control"),
                lambda: gaussian_area_weights(uniform, 128, what="the control"),
                lambda: label_row_weights(gaussian_latitudes(48), 32,
                                          what="the control")):
        try:
            bad()
        except SystemExit:
            refused_axes += 1
    check("an axis that is not Gaussian is refused as an area weight",
          refused_axes == 3, f"{refused_axes} of 3 refused")

    # And the labels the model does write are accepted and land on their own
    # rows, so the control above refuses the convention rather than everything.
    labels = gaussian_latitudes(32)[[0, 5, 16, 31]]
    got = label_row_weights(labels, 32, what="the model's own rows")
    check("a label on the axis takes its own row's weight",
          bool(np.array_equal(got, gaussian_grid(32).dsin[[0, 5, 16, 31]])),
          f"{got}")

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
