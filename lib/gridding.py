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
