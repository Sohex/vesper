#!/usr/bin/env python3
"""SUPERSEDED by build_boundary_conditions.py. Kept only to reproduce the runs
that were built with it, and because `write_sra` still lives here.

It reads the equirectangular map PNGs and does its own conservative remap onto a
Gaussian grid. Both halves are now obsolete: the fork emits Gaussian grids
directly off the mesh, so the remap is a lossy intermediate, and the PNG land
mask uses the `elevation > 0` convention, which floods every dry closed-basin
floor. Do not use it to prepare new runs.

Conservatively remap World Orogen exports to an ExoPlaSim Gaussian grid.

The source exports are never modified. The explicit Orogen land mask controls
land/sea; the land-only heightmap supplies physical elevation on land. Output
topography is geopotential (m2 s-2), as required by PlaSim/ExoPlaSim SRA code
129. SRA arrays are north-to-south and 0..360 degrees eastward.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import warnings

_MPL_CACHE = Path("/tmp/world-matplotlib-cache")
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from numpy.polynomial.legendre import leggauss
from PIL import Image
import yaml

from _paths import ANALYSIS, CONFIG, INPUTS, SOURCE

Image.MAX_IMAGE_PIXELS = None
# netCDF4 1.7.4 emits this harmless warning from its NumPy 2.5 assignment shim.
warnings.filterwarnings(
    "ignore", message="Setting the shape on a NumPy array has been deprecated"
)

EARTH_STANDARD_GRAVITY = 9.80665
OROGEN_FULL_MIN_M = -5000.0
OROGEN_FULL_MAX_M = 6000.0
OROGEN_LAND_MAX_M = 6000.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_one(source: Path, pattern: str) -> Path:
    matches = sorted(source.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one {pattern!r} in {source}, found {matches}")
    return matches[0]


def inventory_source(source: Path) -> dict[str, dict]:
    """Record every canonical export without interpreting unused color maps."""
    inventory = {}
    for path in sorted(item for item in source.iterdir() if item.is_file()):
        record = {"sha256": sha256(path), "bytes": path.stat().st_size}
        if path.suffix.lower() == ".png":
            with Image.open(path) as image:
                record.update(
                    {
                        "dimensions": list(image.size),
                        "mode": image.mode,
                        "bands": list(image.getbands()),
                        "bit_depth": 16 if image.mode == "I;16" else 8,
                    }
                )
        inventory[str(path)] = record
    return inventory


def gaussian_grid(nlat: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return north-to-south latitude centers/bounds and quadrature weights."""
    nodes, weights = leggauss(nlat)
    nodes = nodes[::-1]
    weights = weights[::-1]
    lat = np.degrees(np.arcsin(nodes))
    mu_bounds = np.empty(nlat + 1)
    mu_bounds[0] = 1.0
    mu_bounds[1:] = 1.0 - np.cumsum(weights)
    mu_bounds[-1] = -1.0
    lat_bounds = np.degrees(np.arcsin(np.clip(mu_bounds, -1.0, 1.0)))
    return lat, lat_bounds, weights


def latitude_overlap_matrix(source_height: int, target_mu_bounds: np.ndarray) -> np.ndarray:
    """Area weights mapping equirectangular source rows to target latitude bands."""
    source_lat_bounds = np.linspace(90.0, -90.0, source_height + 1)
    source_mu_bounds = np.sin(np.deg2rad(source_lat_bounds))
    target_north = target_mu_bounds[:-1, None]
    target_south = target_mu_bounds[1:, None]
    source_north = source_mu_bounds[:-1][None, :]
    source_south = source_mu_bounds[1:][None, :]
    overlap = np.maximum(
        0.0,
        np.minimum(target_north, source_north) - np.maximum(target_south, source_south),
    )
    overlap /= target_north - target_south
    if not np.allclose(overlap.sum(axis=1), 1.0, atol=2e-13):
        raise RuntimeError("Latitude overlap weights do not sum to one")
    return overlap


def longitude_block_average(field: np.ndarray, nlon: int) -> np.ndarray:
    """Conservatively average source pixels into 0..360 target longitude cells."""
    height, width = field.shape
    if width % nlon:
        raise ValueError(f"Source width {width} is not divisible by target nlon {nlon}")
    block = width // nlon
    if block % 2:
        raise ValueError("An even number of source pixels per target longitude is required")
    # Orogen spans [-180, 180] with pixel edges at the endpoints. ExoPlaSim's
    # first gridpoint is 0E, so its first cell straddles the prime meridian.
    start = width // 2 - block // 2
    shifted = np.roll(field, -start, axis=1)
    return shifted.reshape(height, nlon, block).mean(axis=2)


def write_sra(path: Path, code: int, field: np.ndarray) -> None:
    nlat, nlon = field.shape
    flat = np.asarray(field, dtype=np.float64).ravel(order="C")
    if flat.size % 8:
        raise ValueError("SRA field size must be divisible by 8")
    header = [code, 0, 20260811, 0, nlon, nlat, 0, 0]
    with path.open("w", encoding="ascii") as handle:
        handle.write("".join(f" {value:11d}" for value in header) + "\n")
        for row in flat.reshape(-1, 8):
            handle.write("".join(f" {value:12.5f}" for value in row) + "\n")


def save_netcdf(
    path: Path,
    lat: np.ndarray,
    lat_bounds: np.ndarray,
    lon: np.ndarray,
    land_fraction: np.ndarray,
    land_mask: np.ndarray,
    elevation_m: np.ndarray,
    geopotential: np.ndarray,
) -> None:
    with Dataset(path, "w") as nc:
        nc.createDimension("lat", lat.size)
        nc.createDimension("lon", lon.size)
        nc.createDimension("bounds", 2)
        latv = nc.createVariable("lat", "f8", ("lat",))
        lonv = nc.createVariable("lon", "f8", ("lon",))
        latbv = nc.createVariable("lat_bounds", "f8", ("lat", "bounds"))
        latv[:] = lat
        lonv[:] = lon
        latbv[:, 0] = lat_bounds[:-1]
        latbv[:, 1] = lat_bounds[1:]
        latv.units = "degrees_north"
        lonv.units = "degrees_east"
        latv.bounds = "lat_bounds"
        lf = nc.createVariable("source_land_fraction", "f4", ("lat", "lon"), zlib=True)
        lm = nc.createVariable("land_mask", "i1", ("lat", "lon"), zlib=True)
        el = nc.createVariable("surface_elevation", "f4", ("lat", "lon"), zlib=True)
        gp = nc.createVariable("surface_geopotential", "f4", ("lat", "lon"), zlib=True)
        lf[:] = land_fraction
        lm[:] = land_mask.astype(np.int8)
        el[:] = elevation_m
        gp[:] = geopotential
        lf.units = "1"
        lm.units = "1"
        el.units = "m"
        gp.units = "m2 s-2"
        nc.title = "World Orogen geography conservatively remapped to ExoPlaSim T42"
        nc.source = "Explicit Orogen land mask plus 0--6 km land-only heightmap"


def plot_diagnostic(
    path: Path,
    preview_mask: np.ndarray,
    preview_elevation: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    land_fraction: np.ndarray,
    mask: np.ndarray,
    elevation: np.ndarray,
    source_area: float,
    target_area: float,
) -> None:
    lon_plot = ((lon + 180.0) % 360.0) - 180.0
    order = np.argsort(lon_plot)
    extent = [-180, 180, -90, 90]
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.5), constrained_layout=True)
    axes[0, 0].imshow(preview_mask, origin="upper", extent=extent, cmap="gray", vmin=0, vmax=1)
    axes[0, 0].set_title(f"Canonical Orogen land mask\narea-weighted land = {source_area:.2%}")
    p = axes[0, 1].imshow(
        preview_elevation,
        origin="upper",
        extent=extent,
        cmap="terrain",
        vmin=0,
        vmax=6000,
    )
    axes[0, 1].set_title("Orogen land-only elevation (preview)")
    fig.colorbar(p, ax=axes[0, 1], label="m")
    q = axes[1, 0].pcolormesh(
        lon_plot[order], lat, land_fraction[:, order], shading="nearest", cmap="Blues", vmin=0, vmax=1
    )
    axes[1, 0].contour(lon_plot[order], lat, mask[:, order], levels=[0.5], colors="black", linewidths=0.45)
    axes[1, 0].set_title("T42 source land fraction; black = retained coastline")
    fig.colorbar(q, ax=axes[1, 0], label="land fraction")
    shown = np.where(mask[:, order], elevation[:, order], np.nan)
    r = axes[1, 1].pcolormesh(
        lon_plot[order], lat, shown, shading="nearest", cmap="terrain", vmin=0, vmax=6000
    )
    axes[1, 1].set_title(f"T42 land elevation\nmodel land = {target_area:.2%}")
    fig.colorbar(r, ax=axes[1, 1], label="m")
    for ax in axes.ravel():
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
    fig.suptitle("World Orogen → ExoPlaSim T42 geography sanity check", fontsize=15)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=INPUTS / "t42")
    parser.add_argument("--analysis", type=Path, default=ANALYSIS / "geography")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    threshold = float(model["geography_land_threshold"])
    if (nlon != 2 * nlat):
        raise ValueError("ExoPlaSim Gaussian grid requires nlon=2*nlat")

    mask_path = find_one(args.source, "orogen-landmask-*.png")
    land_height_path = find_one(args.source, "orogen-land-heightmap-*.png")
    full_height_path = find_one(args.source, "orogen-heightmap-*.png")

    with Image.open(mask_path) as image:
        if image.mode != "RGBA":
            raise ValueError(f"Expected RGBA land mask, got {image.mode}")
        mask_rgba = np.asarray(image)
    channel_values = [np.unique(mask_rgba[..., channel]).tolist() for channel in range(4)]
    if channel_values != [[0, 255], [0, 255], [0, 255], [255]]:
        raise ValueError(f"Unexpected land-mask channel encoding: {channel_values}")
    if not (np.array_equal(mask_rgba[..., 0], mask_rgba[..., 1]) and np.array_equal(mask_rgba[..., 0], mask_rgba[..., 2])):
        raise ValueError("Land-mask RGB channels differ")
    source_mask = mask_rgba[..., 0] == 255
    del mask_rgba

    with Image.open(land_height_path) as image:
        if image.mode != "I;16":
            raise ValueError(f"Expected 16-bit grayscale land heightmap, got {image.mode}")
        land_raw = np.asarray(image)
    if land_raw.shape != source_mask.shape:
        raise ValueError("Mask and heightmap dimensions differ")
    with Image.open(full_height_path) as image:
        if image.mode != "I;16":
            raise ValueError(f"Expected 16-bit grayscale full heightmap, got {image.mode}")
        full_raw = np.asarray(image)
    if full_raw.shape != source_mask.shape:
        raise ValueError("Mask and full heightmap dimensions differ")
    full_raw_min = int(full_raw.min())
    full_raw_max = int(full_raw.max())
    full_saturated_low = int(np.count_nonzero(full_raw == 0))
    full_saturated_high = int(np.count_nonzero(full_raw == 65535))
    full_sea_level_raw = 65535.0 * (-OROGEN_FULL_MIN_M) / (OROGEN_FULL_MAX_M - OROGEN_FULL_MIN_M)
    full_positive = full_raw > full_sea_level_raw
    full_land_nonpositive = int(np.count_nonzero(source_mask & ~full_positive))
    full_ocean_positive = int(np.count_nonzero(~source_mask & full_positive))
    del full_raw, full_positive
    source_elevation = land_raw.astype(np.float64) * (OROGEN_LAND_MAX_M / 65535.0)
    source_elevation[~source_mask] = 0.0
    land_raw_min = int(land_raw.min())
    land_raw_max = int(land_raw.max())
    land_raw_zero = int(np.count_nonzero(land_raw == 0))
    land_raw_saturated_high = int(np.count_nonzero(land_raw == 65535))

    height, width = source_mask.shape
    lat, lat_bounds, gaussian_weights = gaussian_grid(nlat)
    target_mu_bounds = np.sin(np.deg2rad(lat_bounds))
    lat_overlap = latitude_overlap_matrix(height, target_mu_bounds)
    lon = np.arange(nlon, dtype=np.float64) * (360.0 / nlon)

    lon_mask = longitude_block_average(source_mask.astype(np.float64), nlon)
    land_fraction = lat_overlap @ lon_mask
    target_mask = land_fraction >= threshold

    # Average physical elevation over only the land portion of each target cell.
    lon_elev_land = longitude_block_average(source_elevation * source_mask, nlon)
    land_elevation_numerator = lat_overlap @ lon_elev_land
    target_elevation = np.divide(
        land_elevation_numerator,
        land_fraction,
        out=np.zeros_like(land_elevation_numerator),
        where=land_fraction > 0,
    )
    target_elevation[~target_mask] = 0.0

    gravity = EARTH_STANDARD_GRAVITY * float(config["planet"]["mass_earth"]) / float(config["planet"]["radius_earth"]) ** 2
    geopotential = target_elevation * gravity

    source_lat_centers = 90.0 - (np.arange(height) + 0.5) * 180.0 / height
    source_area = float(np.average(source_mask.mean(axis=1), weights=np.cos(np.deg2rad(source_lat_centers))))
    target_area = float(np.sum(gaussian_weights[:, None] * target_mask) / (2.0 * nlon))
    mixed = (land_fraction > 0.0) & (land_fraction < 1.0)

    args.output.mkdir(parents=True, exist_ok=True)
    args.analysis.mkdir(parents=True, exist_ok=True)
    land_sra = args.output / "orogen_T42_surf_0172.sra"
    topo_sra = args.output / "orogen_T42_surf_0129.sra"
    nc_path = args.output / "orogen_T42_geography.nc"
    write_sra(land_sra, 172, target_mask.astype(float))
    write_sra(topo_sra, 129, geopotential)
    save_netcdf(nc_path, lat, lat_bounds, lon, land_fraction, target_mask, target_elevation, geopotential)

    preview_size = (2048, 1024)
    with Image.open(mask_path) as image:
        preview_mask = np.asarray(image.getchannel("R").resize(preview_size, Image.Resampling.BOX), dtype=float) / 255.0
    with Image.open(land_height_path) as image:
        preview_elevation = np.asarray(image.resize(preview_size, Image.Resampling.BOX), dtype=float) * (6000.0 / 65535.0)
    preview_elevation[preview_mask < 0.5] = np.nan
    plot_path = args.analysis / "t42_geography_sanity.png"
    plot_diagnostic(
        plot_path,
        preview_mask,
        preview_elevation,
        lat,
        lon,
        land_fraction,
        target_mask,
        target_elevation,
        source_area,
        target_area,
    )

    report = {
        "source": {
            "dimensions": [width, height],
            "projection": "equirectangular; -180..180 longitude left-to-right; 90..-90 latitude top-to-bottom",
            "land_mask": "RGBA8; RGB exactly 0 ocean or 255 land; alpha 255",
            "full_heightmap": "grayscale16 linear; 0=-5000 m, 65535=6000 m (reference only)",
            "land_heightmap": "grayscale16 linear; 0=0 m/ocean, 65535=6000 m",
            "pixel_land_fraction": float(source_mask.mean()),
            "land_heightmap_validation": {
                "raw_min": land_raw_min,
                "raw_max": land_raw_max,
                "zero_pixels_including_ocean": land_raw_zero,
                "saturated_max_pixels": land_raw_saturated_high,
            },
            "full_heightmap_validation": {
                "raw_min": full_raw_min,
                "raw_max": full_raw_max,
                "saturated_min_pixels": full_saturated_low,
                "saturated_max_pixels": full_saturated_high,
                "mask_land_but_full_height_nonpositive_pixels": full_land_nonpositive,
                "mask_ocean_but_full_height_positive_pixels": full_ocean_positive,
                "coastal_disagreement_explanation": "Orogen flat-shades its mask but interpolates heightmap vertices",
            },
            "area_weighted_land_fraction": source_area,
            "files": inventory_source(args.source),
        },
        "target": {
            "grid": f"T42 {nlat}x{nlon} Gaussian",
            "orientation": "north-to-south rows; 0..360 eastward columns",
            "land_threshold": threshold,
            "area_weighted_land_fraction": target_area,
            "land_fraction_change_percentage_points": 100.0 * (target_area - source_area),
            "land_cells": int(target_mask.sum()),
            "ocean_cells": int((~target_mask).sum()),
            "mixed_source_cells": int(mixed.sum()),
            "mixed_source_cells_percent": float(100.0 * mixed.mean()),
            "land_elevation_m": {
                "min": float(target_elevation[target_mask].min()),
                "mean": float(np.average(target_elevation, weights=gaussian_weights[:, None] * target_mask)),
                "max": float(target_elevation.max()),
            },
            "gravity_m_s2": gravity,
            "outputs": [str(land_sra), str(topo_sra), str(nc_path), str(plot_path)],
        },
    }
    report_path = args.analysis / "geography_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
