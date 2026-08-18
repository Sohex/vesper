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


def land_weighted(export: Export, grid_dir: Path, values: np.ndarray):
    """Per-cell land fraction and the land-area-weighted mean of `values`.

    Returns (land_fraction, mean_over_land), both shaped (nlat, nlon). Cells with
    no land get a mean of zero, which callers should mask rather than use.
    """
    cell, nlat, nlon = region_cells(export, grid_dir)
    ncell = nlat * nlon
    area = export.cell_area.astype(np.float64)
    is_land = export.surface_class == LAND

    total = np.zeros(ncell)
    land_area = np.zeros(ncell)
    weighted = np.zeros(ncell)
    np.add.at(total, cell, area)
    np.add.at(land_area, cell[is_land], area[is_land])
    np.add.at(weighted, cell[is_land], area[is_land] * values[is_land])

    empty = total <= 0
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(land_area > 0, weighted / np.maximum(land_area, 1e-30), 0.0)
    fraction = np.where(empty, 0.0, land_area / np.maximum(total, 1e-30))

    if empty.any():
        # Cells finer than the mesh. At T85 this is 35 of 32,768, all polar,
        # because the 2.5M-region mesh averages 17 km spacing against a polar
        # cell far narrower than that. Fall back to the nearest region centre,
        # which is what the exporter's own categorical rule does. Refining the
        # mesh is not a free fix: region count is a generation parameter, so a
        # finer mesh is a different planet.
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
    """Fraction of each cell's *land* area for which `mask` is true."""
    cell, nlat, nlon = region_cells(export, grid_dir)
    ncell = nlat * nlon
    area = export.cell_area.astype(np.float64)
    is_land = export.surface_class == LAND
    sel = is_land & mask

    land_area = np.zeros(ncell)
    hit_area = np.zeros(ncell)
    np.add.at(land_area, cell[is_land], area[is_land])
    np.add.at(hit_area, cell[sel], area[sel])
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(land_area > 0, hit_area / np.maximum(land_area, 1e-30), 0.0)
    return out.reshape(nlat, nlon)

# --- the grid convention, and the one check that proves it ------------------

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
    would ever have failed. See `notes/failure-modes.md` class 17.
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
