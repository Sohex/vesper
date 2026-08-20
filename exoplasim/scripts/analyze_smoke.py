#!/usr/bin/env python3
"""Audit and plot a completed one-orbit ExoPlaSim smoke run."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

_MPL_CACHE = Path("/tmp/world-matplotlib-cache")
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS, INPUTS

import climatology  # noqa: E402  from lib/, put on sys.path by _paths
import gridding


def weighted_mean(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Global mean over the last two (latitude, longitude) dimensions."""
    return np.sum(field * weights[None, :, None], axis=(-2, -1)) / (2.0 * field.shape[-1])


def periodic_plot_order(lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Display order only; `lib/gridding.py` owns every longitude expression."""
    return gridding.display_longitude(lon)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output", type=Path, default=ANALYSIS / "smoke")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    args.output.mkdir(parents=True, exist_ok=True)

    regular_path = run_dir / "MOST.00000.nc"
    with Dataset(regular_path) as nc:
        lat = np.asarray(nc["lat"][:], dtype=float)
        lon = np.asarray(nc["lon"][:], dtype=float)
        lev = np.asarray(nc["lev"][:], dtype=float)
        centres = np.asarray(nc["time"][:], dtype=float)
        data = {
            name: np.asarray(nc[name][:], dtype=float)
            for name in [
                "ts", "tas", "pr", "ua", "va", "snd", "sic", "sit",
                "mrso", "mrro", "evap", "ntr", "hfns", "lsm", "sg",
            ]
        }
    weights = leggauss(lat.size)[1][::-1]
    lsm = data["lsm"][0]
    ocean = 1.0 - lsm
    ocean_area = np.sum(ocean * weights[:, None])
    land_area = np.sum(lsm * weights[:, None])

    ts = data["ts"]
    pr = data["pr"] * 86400.0 * 1000.0
    sic = data["sic"]
    summary = {
        "warning": "One-orbit smoke output; not an equilibrium climate.",
        "run_dir": str(run_dir),
        "time_bins": int(ts.shape[0]),
        "lowest_model_level_sigma": float(lev[-1]),
        "global_orbit_means": {
            "surface_temperature_k": float(weighted_mean(ts, weights).mean()),
            "air_temperature_2m_k": float(weighted_mean(data["tas"], weights).mean()),
            "precipitation_mm_day": float(weighted_mean(pr, weights).mean()),
            "toa_net_radiation_w_m2": float(weighted_mean(data["ntr"], weights).mean()),
            "surface_downward_heat_flux_w_m2": float(weighted_mean(data["hfns"], weights).mean()),
            "planetary_sea_ice_fraction": float(weighted_mean(sic, weights).mean()),
            "ocean_sea_ice_fraction": float(
                np.mean(np.sum(sic * ocean[None, :, :] * weights[None, :, None], axis=(-2, -1)) / ocean_area)
            ),
        },
        "global_monthly_bins": {
            "surface_temperature_k": weighted_mean(ts, weights).tolist(),
            "toa_net_radiation_w_m2": weighted_mean(data["ntr"], weights).tolist(),
            "surface_downward_heat_flux_w_m2": weighted_mean(data["hfns"], weights).tolist(),
            "precipitation_mm_day": weighted_mean(pr, weights).tolist(),
        },
    }

    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    gravity = float(manifest["derived_parameters"]["gravity_m_s2"])
    # The supplied topography is surf code 0129 (geopotential); divide by g.
    import sra as _sra
    supplied_elevation = _sra.read_sra(
        INPUTS / "t42" / "orogen_T42_surf_0129.sra",
        lat.size, lon.size) / gravity
    model_elevation = data["sg"][0] / gravity
    difference = model_elevation - supplied_elevation
    land_bool = lsm > 0.5
    summary["spectral_topography"] = {
        "supplied_land_min_mean_max_m": [
            float(supplied_elevation[land_bool].min()),
            float(np.sum(supplied_elevation * lsm * weights[:, None]) / land_area),
            float(supplied_elevation[land_bool].max()),
        ],
        "model_land_min_mean_max_m": [
            float(model_elevation[land_bool].min()),
            float(np.sum(model_elevation * lsm * weights[:, None]) / land_area),
            float(model_elevation[land_bool].max()),
        ],
        "land_rmse_m": float(np.sqrt(np.mean(difference[land_bool] ** 2))),
        "land_correlation": float(np.corrcoef(supplied_elevation[land_bool], model_elevation[land_bool])[0, 1]),
        "negative_model_land_cells": int(np.count_nonzero(model_elevation[land_bool] < 0)),
        "negative_model_ocean_cells": int(np.count_nonzero(model_elevation[~land_bool] < 0)),
    }

    lon_plot, order = periodic_plot_order(lon)
    extent = [-180, 180, -90, 90]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    fields = [supplied_elevation, model_elevation, difference]
    titles = ["Supplied T42 elevation", "ExoPlaSim spectral elevation", "Model minus supplied"]
    cmaps = ["terrain", "terrain", "RdBu_r"]
    limits = [(0, 6000), (-1000, 6000), (-3000, 3000)]
    for ax, field, title, cmap, (vmin, vmax) in zip(axes, fields, titles, cmaps, limits):
        shown = field[:, order]
        image = ax.imshow(shown, origin="upper", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.contour(lon_plot, lat, lsm[:, order], levels=[0.5], colors="black", linewidths=0.35)
        ax.set_title(title)
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        fig.colorbar(image, ax=ax, label="m", shrink=0.83)
    fig.suptitle("Custom topography after ExoPlaSim T42 spectral truncation")
    topo_plot = args.output / "model_used_topography.png"
    fig.savefig(topo_plot, dpi=180)
    plt.close(fig)

    # Annual means are weighted by records per bin, which are unequal because
    # pyburn's linspace().astype(int) split spreads the remainder. CLIM-13.

    def yearly(field):
        return climatology.annual_mean(np.asarray(field), centres)

    annual_ts = yearly(ts)
    seasonal_range = ts.max(axis=0) - ts.min(axis=0)
    annual_pr = yearly(pr)
    annual_sic = yearly(sic)
    u_low = yearly(data["ua"])[-1]
    v_low = yearly(data["va"])[-1]
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.5), constrained_layout=True)
    panels = [
        (annual_ts, "First-orbit mean surface temperature", "coolwarm", 220, 310, "K"),
        (seasonal_range, "First-orbit surface-temperature range", "magma", 0, 80, "K"),
        (annual_pr, "First-orbit mean precipitation", "YlGnBu", 0, 8, "mm day$^{-1}$"),
        (annual_sic, "First-orbit mean sea-ice cover", "Blues", 0, 1, "fraction"),
    ]
    for ax, (field, title, cmap, vmin, vmax, label) in zip(axes.ravel(), panels):
        image = ax.imshow(field[:, order], origin="upper", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.contour(lon_plot, lat, lsm[:, order], levels=[0.5], colors="black", linewidths=0.35)
        ax.set_title(title)
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        fig.colorbar(image, ax=ax, label=label, shrink=0.83)
    skip = (slice(2, -2, 4), slice(None, None, 4))
    lon_mesh, lat_mesh = np.meshgrid(lon_plot, lat)
    axes[0, 0].quiver(
        lon_mesh[skip], lat_mesh[skip], u_low[:, order][skip], v_low[:, order][skip],
        color="black", alpha=0.55, scale=450, width=0.0016,
    )
    axes[0, 0].text(
        0.01, 0.02, f"Arrows: mean wind at lowest model level (sigma={lev[-1]:.3f})",
        transform=axes[0, 0].transAxes, fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none"},
    )
    fig.suptitle("Runtime sanity check — first orbit only, NOT equilibrated")
    climate_plot = args.output / "first_orbit_climate_sanity.png"
    fig.savefig(climate_plot, dpi=180)
    plt.close(fig)

    summary["outputs"] = [str(topo_plot), str(climate_plot)]
    summary_path = args.output / "first_orbit_audit.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
