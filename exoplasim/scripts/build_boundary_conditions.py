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
from gridding import land_weighted
from orogen import Export, LAND

LAND_MASK_CODE = 172
TOPOGRAPHY_CODE = 129


def build(mesh: Export, grid_dir: Path, threshold: float, gravity: float):
    """Land mask and geopotential on `grid_dir`'s grid, integrated from the mesh.

    Elevation is metres, so the land-weighted mean is the mean elevation of the
    land in each cell. Ocean depths never enter it.
    """
    elev_m = mesh.elevation_km.astype(np.float64) * 1000.0
    fraction, mean_elev, empty = land_weighted(mesh, grid_dir, elev_m)
    mask = (fraction >= threshold).astype(np.float64)
    mean_elev = np.where(mask > 0, mean_elev, 0.0)
    return {
        "land_mask": mask,
        "cells_below_mesh_resolution": empty,
        "land_fraction": fraction,
        "geopotential": mean_elev * gravity,
        "elevation_m": mean_elev,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--mesh", type=Path, default=SOURCE / "exoplasim-T42",
                    help="export carrying raw/; the mesh is the same for all")
    ap.add_argument("--grid", type=Path, default=None,
                    help="export whose grid to target; defaults to the config resolution")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model, planet = config["model"], config["planet"]
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    threshold = float(model["geography_land_threshold"])
    gravity = float(planet["gravity_m_s2"])

    resolution = str(model["resolution"]).upper()
    grid_dir = args.grid or (SOURCE / f"exoplasim-{resolution}")
    output = args.output or (INPUTS / resolution.lower())
    ex = Export(args.mesh)
    out = build(ex, grid_dir, threshold, gravity)
    if out["land_mask"].shape != (nlat, nlon):
        raise RuntimeError(
            f"export grid is {out['land_mask'].shape}, config asks for {(nlat, nlon)}"
        )

    output.mkdir(parents=True, exist_ok=True)
    write_sra(output / f"orogen_{resolution}_surf_{LAND_MASK_CODE:04d}.sra",
              LAND_MASK_CODE, out["land_mask"])
    write_sra(output / f"orogen_{resolution}_surf_{TOPOGRAPHY_CODE:04d}.sra",
              TOPOGRAPHY_CODE, out["geopotential"])

    gw = np.fromfile(grid_dir / "grid" / "gauss_weights.bin", dtype="float64")
    m, e = out["land_mask"], out["elevation_m"]
    land_cells = m > 0
    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "mesh": ex.provenance(),
        "grid": str(grid_dir),
        "resolution": resolution,
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
        "cells_below_mesh_resolution": out["cells_below_mesh_resolution"],
        "elevation_m": {
            "min": float(e[land_cells].min()),
            "max": float(e[land_cells].max()),
            "mean": float(e[land_cells].mean()),
            "cells_below_sea_level": int((e[land_cells] < 0).sum()),
        },
        "codes": [LAND_MASK_CODE, TOPOGRAPHY_CODE],
    }
    (output / "boundary_conditions_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
