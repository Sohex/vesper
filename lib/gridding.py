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


def grid_geometry(grid_dir: Path):
    """Latitude and longitude cell centres of an export's grid."""
    gm = json.loads((grid_dir / "manifest.json").read_text(encoding="utf-8"))["grid"]
    lat = np.fromfile(grid_dir / gm["coords"]["lat"]["path"], dtype="float64")
    lon = np.fromfile(grid_dir / gm["coords"]["lon"]["path"], dtype="float64")
    return lat, lon, gm


def region_cells(export: Export, grid_dir: Path):
    """Flat grid-cell index for every mesh region, plus the grid shape.

    Rows run north to south and sit at Gauss-Legendre latitudes on a spectral
    grid, so bin on the midpoints between row centres rather than assuming
    uniform spacing.
    """
    lat, lon, _ = grid_geometry(grid_dir)
    nlat, nlon = lat.size, lon.size
    edges = np.empty(nlat + 1)
    edges[1:-1] = 0.5 * (lat[:-1] + lat[1:])
    edges[0], edges[-1] = 90.0, -90.0
    row = np.clip(np.searchsorted(-edges, -export.lat, side="right") - 1, 0, nlat - 1)
    col = np.clip(((export.lon + 180.0) / 360.0 * nlon).astype(np.int64), 0, nlon - 1)
    return row * nlon + col, nlat, nlon


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

    Latitude IS matched by nearest centre, because the export grid and the
    climatology are genuinely different Gaussian axes there -- both descending
    Gauss-Legendre rows, but resolved and labelled independently -- so index
    identity is not available and nearest centre is the right join.

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
    rows = np.abs(clim_lat[None, :] - glat[:, None]).argmin(axis=1)
    if np.unique(rows).size != rows.size:
        raise RuntimeError(
            "the nearest-centre latitude join is not a bijection: two export "
            "rows chose the same climatology row, so the two axes are not the "
            "same Gaussian grid")
    return rows[row], col


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
    import numpy as np
    from netCDF4 import Dataset
    with Dataset(coupling) as ds:
        cell = np.asarray(ds["cell"][:]).astype(np.int64)
        area = np.asarray(ds["area_km2"][:])
        nlon = int(ds.n_lon)
    if lsm.shape[1] != nlon:
        raise ValueError(
            f"coupling has {nlon} columns and the mask has {lsm.shape[1]}; "
            "they are not the same grid and no mapping between them is defined")
    row, col = np.divmod(cell, nlon)
    return float(area[lsm[row, col] < 0.5].sum() / area.sum())


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
