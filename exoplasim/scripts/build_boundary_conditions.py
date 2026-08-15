#!/usr/bin/env python3
"""Write ExoPlaSim land mask and topography from the World Orogen mesh.

Replaces `convert_orogen.py`, which read the equirectangular map PNGs and did its
own conservative remap onto a Gaussian grid. That path is obsolete twice over:
the fork now emits Gaussian grids directly off the mesh, and the PNG land mask
uses the `elevation > 0` convention, which floods every dry closed-basin floor.

Two decisions worth stating, because they are where this differs from just
reading `planet.nc`:

**Land comes from `surface_class`, not `land_mask`.** The two disagree by 1.9% of
the planet, all of it dry basin floor below sea level, which is the terrain the
fork exists to preserve. `land_mask` would put it under water.

**Both fields are integrated from the native mesh, not sampled from the gridded
export.** The export resamples categorical fields, `surface_class` among them, by
taking the value of the region containing the cell centre. At T42 a cell holds
roughly 300 mesh regions, so a point sample throws away the coastline. Here the
land fraction of each cell is the area-weighted fraction of its regions that are
land, thresholded at `model.geography_land_threshold`, and topography is the
land-area-weighted mean elevation over the land regions only, so ocean depths
never drag a coastal cell's elevation down.

Elevation is written as geopotential, which is what SRA code 129 expects, using
this planet's gravity. Dry basin floors keep their negative elevation: a floor
562 m below sea level is 562 m below sea level, and flattening it to zero would
undo the preservation the whole pipeline is built around.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, INPUTS, SOURCE
from convert_orogen import write_sra
from orogen import Export, LAND

LAND_MASK_CODE = 172
TOPOGRAPHY_CODE = 129


def bin_to_grid(export: Export, grid_dir: Path):
    """Assign every mesh region to the grid cell containing its centre."""
    gm = json.loads((grid_dir / "manifest.json").read_text(encoding="utf-8"))["grid"]
    lat = np.fromfile(grid_dir / gm["coords"]["lat"]["path"], dtype="float64")
    lon = np.fromfile(grid_dir / gm["coords"]["lon"]["path"], dtype="float64")
    nlat, nlon = lat.size, lon.size

    # Rows run north to south and sit at Gauss-Legendre latitudes, so bin on the
    # midpoints between row centres rather than assuming uniform spacing.
    edges = np.empty(nlat + 1)
    edges[1:-1] = 0.5 * (lat[:-1] + lat[1:])
    edges[0], edges[-1] = 90.0, -90.0
    row = np.clip(np.searchsorted(-edges, -export.lat, side="right") - 1, 0, nlat - 1)
    col = np.clip(((export.lon + 180.0) / 360.0 * nlon).astype(np.int64), 0, nlon - 1)
    return row * nlon + col, nlat, nlon, lat, lon


def build(export: Export, grid_dir: Path, threshold: float, gravity: float):
    cell, nlat, nlon, lat, lon = bin_to_grid(export, grid_dir)
    ncell = nlat * nlon
    area = export.cell_area.astype(np.float64)
    elev = export.elevation_km.astype(np.float64) * 1000.0
    is_land = export.surface_class == LAND

    total = np.zeros(ncell)
    land_area = np.zeros(ncell)
    land_elev = np.zeros(ncell)
    np.add.at(total, cell, area)
    np.add.at(land_area, cell[is_land], area[is_land])
    np.add.at(land_elev, cell[is_land], area[is_land] * elev[is_land])

    if (total <= 0).any():
        raise RuntimeError(
            f"{int((total <= 0).sum())} grid cells contain no mesh region; "
            "the grid is finer than the mesh and the mask would be undefined"
        )

    fraction = land_area / total
    mask = (fraction >= threshold).astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_elev = np.where(land_area > 0, land_elev / np.maximum(land_area, 1e-30), 0.0)
    # Elevation is only meaningful where the cell is land; elsewhere sea level.
    mean_elev = np.where(mask > 0, mean_elev, 0.0)

    return {
        "land_mask": mask.reshape(nlat, nlon),
        "land_fraction": fraction.reshape(nlat, nlon),
        "geopotential": (mean_elev * gravity).reshape(nlat, nlon),
        "elevation_m": mean_elev.reshape(nlat, nlon),
        "lat": lat, "lon": lon,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--export", type=Path, default=SOURCE / "exoplasim-T42")
    ap.add_argument("--output", type=Path, default=INPUTS / "t42")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model, planet = config["model"], config["planet"]
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    threshold = float(model["geography_land_threshold"])
    gravity = float(planet["gravity_m_s2"])

    ex = Export(args.export)
    out = build(ex, args.export, threshold, gravity)
    if out["land_mask"].shape != (nlat, nlon):
        raise RuntimeError(
            f"export grid is {out['land_mask'].shape}, config asks for {(nlat, nlon)}"
        )

    args.output.mkdir(parents=True, exist_ok=True)
    write_sra(args.output / f"orogen_T42_surf_{LAND_MASK_CODE:04d}.sra",
              LAND_MASK_CODE, out["land_mask"])
    write_sra(args.output / f"orogen_T42_surf_{TOPOGRAPHY_CODE:04d}.sra",
              TOPOGRAPHY_CODE, out["geopotential"])

    gw = np.fromfile(args.export / "grid" / "gauss_weights.bin", dtype="float64")
    m, e = out["land_mask"], out["elevation_m"]
    land_cells = m > 0
    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "source": ex.provenance(),
        "method": ("Area-weighted land fraction per cell from the native mesh, "
                   "thresholded; topography is the land-area-weighted mean "
                   "elevation over land regions only."),
        "land_definition": "surface_class == land (subaerial)",
        "land_threshold": threshold,
        "gravity_m_s2": gravity,
        "land_fraction_gauss_weighted": float(
            (m.mean(axis=1) * gw).sum() / gw.sum()),
        "mesh_land_fraction": float(
            ex.cell_area[ex.surface_class == LAND].sum() / ex.cell_area.sum()),
        "land_cells": int(land_cells.sum()),
        "elevation_m": {
            "min": float(e[land_cells].min()),
            "max": float(e[land_cells].max()),
            "mean": float(e[land_cells].mean()),
            "cells_below_sea_level": int((e[land_cells] < 0).sum()),
        },
        "codes": [LAND_MASK_CODE, TOPOGRAPHY_CODE],
    }
    (args.output / "boundary_conditions_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
