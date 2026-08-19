#!/usr/bin/env python3
"""Create climate diagnostics, maps, and a rate-normalized Koppen interpretation."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import warnings

_MPL_CACHE = Path("/tmp/world-matplotlib-cache")
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS  # noqa: F401  (puts lib/ on the path)
import climatology
import gridding


EARTH_YEAR_DAYS = 365.2425
SEASONS = [
    ("Northern spring", 0.0, 90.0),
    ("Northern summer", 90.0, 180.0),
    ("Northern autumn", 180.0, 270.0),
    ("Northern winter", 270.0, 360.0),
]

KOPPEN_COLORS = {
    "Af": "#006400", "Am": "#007f32", "Aw": "#b5c900", "As": "#d6d500",
    "BWh": "#ff3300", "BWk": "#ff8c00", "BSh": "#f5a742", "BSk": "#ffd27f",
    "Csa": "#ffff00", "Csb": "#c8c800", "Csc": "#969600",
    "Cwa": "#96ff96", "Cwb": "#64c864", "Cwc": "#329632",
    "Cfa": "#00ffff", "Cfb": "#64ffff", "Cfc": "#c8ffff",
    "Dsa": "#ff00ff", "Dsb": "#c800c8", "Dsc": "#963296", "Dsd": "#966496",
    "Dwa": "#a00000", "Dwb": "#c80000", "Dwc": "#dc6464", "Dwd": "#c86464",
    "Dfa": "#0000ff", "Dfb": "#3232c8", "Dfc": "#6464ff", "Dfd": "#9696ff",
    "ET": "#b2b2b2", "EF": "#eeeeee",
}


def load(path: Path, names: list[str]) -> tuple[dict[str, np.ndarray], dict]:
    with Dataset(path) as nc:
        values = {name: np.asarray(nc[name][:], dtype=float) for name in names}
        attrs = {name: nc.getncattr(name) for name in nc.ncattrs()}
    return values, attrs


def shifted(lon: np.ndarray, field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reorder a field for DISPLAY on a -180..180 axis. Presentation only.

    The rotation itself is `lib/gridding.py:display_longitude`, and it lives
    there rather than here because that module owns every piece of longitude
    arithmetic in this project. This is a two-line alias so the plotting code
    reads the way it always has.

    **Nothing computed may go through here, and nothing does.** The index is
    canonical and the labels are decoration. See the note beside
    `display_longitude` and `notes/audits/grid-convention-and-runoff.md`.
    """
    return gridding.display_longitude(lon, field)


def decorate(ax: plt.Axes, lat: np.ndarray, lon: np.ndarray, land: np.ndarray) -> None:
    plon, pland = shifted(lon, land)
    ax.contour(plon, lat, pland, levels=[0.5], colors="black", linewidths=0.32, alpha=0.7)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-60, 61, 30))
    ax.set_xlabel("longitude (-180..180, the maps/ convention; the model's own axis is 0..360)")
    ax.set_ylabel("latitude")


def panel(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray, field: np.ndarray,
          land: np.ndarray, title: str, cmap: str | mcolors.Colormap,
          vmin: float | None = None, vmax: float | None = None,
          norm: mcolors.Normalize | None = None):
    plon, pdata = shifted(lon, field)
    mesh = ax.pcolormesh(plon, lat, pdata, shading="auto", cmap=cmap,
                         vmin=vmin, vmax=vmax, norm=norm, rasterized=True)
    decorate(ax, lat, lon, land)
    ax.set_title(title)
    return mesh


def save_map(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def season_masks(lambda_deg: np.ndarray) -> list[np.ndarray]:
    phase = np.mod(lambda_deg, 360.0)
    return [(phase >= lower) & (phase < upper) for _, lower, upper in SEASONS]


def reconstruct_regular_phase(
    regular_time: np.ndarray,
    snapshot_time: np.ndarray,
    snapshot_lambda: np.ndarray,
) -> np.ndarray:
    unwrapped = np.rad2deg(np.unwrap(np.deg2rad(snapshot_lambda)))
    slope = np.median(np.diff(unwrapped) / np.diff(snapshot_time))
    phase = np.interp(regular_time, snapshot_time, unwrapped)
    before = regular_time < snapshot_time[0]
    after = regular_time > snapshot_time[-1]
    phase[before] = unwrapped[0] + slope * (regular_time[before] - snapshot_time[0])
    phase[after] = unwrapped[-1] + slope * (regular_time[after] - snapshot_time[-1])
    return np.mod(phase, 360.0)


def koppen(
    temperature_k: np.ndarray,
    precipitation_m_s: np.ndarray,
    phase_deg: np.ndarray,
    lat: np.ndarray,
    land: np.ndarray,
    bin_centres: np.ndarray,
) -> np.ndarray:
    """Koppen-Geiger classes using Earth-year-normalized precipitation rates."""
    temp = temperature_k - 273.15
    # Each bin is charged with the share of the year its own record count
    # covers, not a flat twelfth: the bins are not equal. CLIM-13.
    share = climatology.bin_weights(bin_centres)
    pmonth = (np.maximum(precipitation_m_s, 0.0) * 86400.0
              * (EARTH_YEAR_DAYS * share[:, None, None]) * 1000.0)
    pann = pmonth.sum(axis=0)
    tmean = climatology.annual_mean(temp, bin_centres)
    tmin = temp.min(axis=0)
    tmax = temp.max(axis=0)
    result = np.full(land.shape, "Ocean", dtype="U8")

    for iy, latitude in enumerate(lat):
        north_summer = (phase_deg >= 0.0) & (phase_deg < 180.0)
        summer = north_summer if latitude >= 0.0 else ~north_summer
        winter = ~summer
        for ix in np.flatnonzero(land[iy] >= 0.5):
            p = pmonth[:, iy, ix]
            psum_s = p[summer].sum()
            psum_w = p[winter].sum()
            wet_fraction_s = psum_s / max(pann[iy, ix], 1e-12)
            if wet_fraction_s >= 0.70:
                aridity = 20.0 * tmean[iy, ix] + 280.0
            elif psum_w / max(pann[iy, ix], 1e-12) >= 0.70:
                aridity = 20.0 * tmean[iy, ix]
            else:
                aridity = 20.0 * tmean[iy, ix] + 140.0
            aridity = max(aridity, 0.0)

            if pann[iy, ix] < 0.5 * aridity:
                result[iy, ix] = "BWh" if tmean[iy, ix] >= 18.0 else "BWk"
                continue
            if pann[iy, ix] < aridity:
                result[iy, ix] = "BSh" if tmean[iy, ix] >= 18.0 else "BSk"
                continue
            if tmin[iy, ix] >= 18.0:
                pmin = p.min()
                if pmin >= 60.0:
                    result[iy, ix] = "Af"
                elif pmin >= 100.0 - pann[iy, ix] / 25.0:
                    result[iy, ix] = "Am"
                else:
                    result[iy, ix] = "As" if p[summer].min() < p[winter].min() else "Aw"
                continue
            if tmax[iy, ix] < 10.0:
                result[iy, ix] = "EF" if tmax[iy, ix] < 0.0 else "ET"
                continue

            group = "C" if tmin[iy, ix] > 0.0 else "D"
            dry_s = p[summer].min() < 40.0 and p[summer].min() < p[winter].max() / 3.0
            dry_w = p[winter].min() < p[summer].max() / 10.0
            second = "s" if dry_s else ("w" if dry_w else "f")
            warm_count = np.count_nonzero(temp[:, iy, ix] > 10.0)
            if tmax[iy, ix] >= 22.0:
                third = "a"
            elif warm_count >= temperature_k.shape[0] / 3.0:
                third = "b"
            else:
                third = "d" if group == "D" and tmin[iy, ix] <= -38.0 else "c"
            result[iy, ix] = group + second + third
    return result


def biome_for(code: str) -> str:
    if code == "Ocean": return "Ocean"
    if code == "Af": return "Everwet tropical forest"
    if code == "Am": return "Monsoon tropical forest"
    if code in {"Aw", "As"}: return "Seasonal tropical woodland/savanna"
    if code.startswith("BW"): return "Desert"
    if code.startswith("BS"): return "Steppe/semidesert"
    if code.startswith("Cs"): return "Mediterranean woodland/shrubland"
    if code.startswith("C") and code.endswith(("a", "b")): return "Temperate/subtropical forest"
    if code.startswith("C"): return "Cool oceanic forest/heath"
    if code.startswith("D") and code.endswith(("a", "b")): return "Continental temperate/mixed forest"
    if code.startswith("D"): return "Boreal forest/woodland"
    if code == "ET": return "Tundra/alpine"
    if code == "EF": return "Permanent ice"
    return "Unclassified"


def area_fractions(labels: np.ndarray, lat: np.ndarray, land: np.ndarray) -> dict[str, float]:
    weights = leggauss(len(lat))[1][::-1, None] * np.ones((1, land.shape[1]))
    weights = np.where(land >= 0.5, weights, 0.0)
    total = weights.sum()
    return {
        str(label): float(weights[labels == label].sum() / total)
        for label in sorted(set(labels[land >= 0.5].ravel()))
    }


def categorical_map(
    labels: np.ndarray,
    colors: dict[str, str],
    lon: np.ndarray,
    lat: np.ndarray,
    land: np.ndarray,
    title: str,
    output: Path,
) -> tuple[np.ndarray, list[str]]:
    categories = [name for name in colors if np.any(labels == name)]
    indices = np.zeros(labels.shape, dtype=np.int16)
    for index, name in enumerate(categories, start=1):
        indices[labels == name] = index
    cmap = mcolors.ListedColormap(["#9bc9e2"] + [colors[name] for name in categories])
    norm = mcolors.BoundaryNorm(np.arange(-0.5, len(categories) + 1.5), cmap.N)
    fig, ax = plt.subplots(figsize=(14, 6.5), constrained_layout=True)
    mesh = panel(ax, lon, lat, indices, land, title, cmap, norm=norm)
    cbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.10, fraction=0.08)
    cbar.set_ticks(np.arange(1, len(categories) + 1))
    cbar.set_ticklabels(categories)
    cbar.ax.tick_params(labelsize=8, rotation=45)
    save_map(fig, output)
    return indices, categories


def main() -> None:
    parser = argparse.ArgumentParser()
    # --label is required and the climatology paths are DERIVED from it.
    #
    # These used to default to a fixed `baseline_regular_climatology.nc`, so
    # running with no arguments silently analysed whichever file had been given
    # that name -- which was the superseded first-era climatology, reporting
    # 280.9 K for a world whose baseline is 293.8 K. build_climatology.py names
    # its products by label; this reads them back by the same label, so the two
    # cannot disagree.
    parser.add_argument("--label", required=True,
                        help="climatology label, as passed to build_climatology.py")
    parser.add_argument("--regular", type=Path, default=None,
                        help="override; normally derived from --label")
    parser.add_argument("--snapshots", type=Path, default=None,
                        help="override; normally derived from --label")
    parser.add_argument("--output", type=Path, default=ANALYSIS / "climatology")
    args = parser.parse_args()
    clim = ANALYSIS / "climatology"
    if args.regular is None:
        args.regular = clim / f"{args.label}_regular_climatology.nc"
    if args.snapshots is None:
        args.snapshots = clim / f"{args.label}_snapshot_climatology.nc"
    for pth in (args.regular, args.snapshots):
        if not pth.is_file():
            raise SystemExit(f"no climatology at {pth}; build it first with "
                             f"build_climatology.py --label {args.label}")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    regular, regular_attrs = load(args.regular, [
        "time", "lat", "lon", "lev", "ts", "tas", "pr", "ua", "va", "snd",
        "sic", "sit", "mrso", "mrro", "evap", "ntr", "hfns", "lsm", "sg",
    ])
    snap, snapshot_attrs = load(args.snapshots, [
        "time", "lambda", "zdec", "ts", "tas", "snd", "sic", "sit", "lsm",
    ])
    lat, lon = regular["lat"], regular["lon"]
    land = regular["lsm"][0]
    gravity = 10.215260416666666
    elevation = regular["sg"][0] / gravity
    # Bins hold unequal numbers of raw records; weight by them. CLIM-13.
    rtime = regular["time"]

    def rmean(field):
        return climatology.annual_mean(field, rtime)

    annual_ts = rmean(regular["ts"])
    annual_tas = rmean(regular["tas"])
    annual_pr = rmean(np.maximum(regular["pr"], 0.0)) * 86400.0 * 1000.0
    display_label = args.label.replace("_", " ").replace("-", " ").title()
    phase_regular = reconstruct_regular_phase(regular["time"], snap["time"], snap["lambda"])

    fig, ax = plt.subplots(figsize=(14, 6), constrained_layout=True)
    mesh = panel(ax, lon, lat, annual_ts - 273.15, land,
                 f"{display_label} annual mean surface temperature (5-orbit climatology)",
                 "coolwarm", -35, 35)
    fig.colorbar(mesh, ax=ax, label="°C", orientation="horizontal", pad=0.09)
    save_map(fig, output / f"{args.label}_annual_surface_temperature.png")

    masks = season_masks(snap["lambda"])
    # snap is the 32-sample snapshot stream on its OWN evenly spaced axis,
    # not the 12 regular bins, so no bin weighting applies here.
    seasonal_ts = [snap["ts"][mask].mean(axis=0) - 273.15 for mask in masks]
    fig, axes = plt.subplots(2, 2, figsize=(15, 8), constrained_layout=True)
    for ax, (name, _, _), field in zip(axes.ravel(), SEASONS, seasonal_ts):
        mesh = panel(ax, lon, lat, field, land, name, "coolwarm", -45, 40)
    fig.colorbar(mesh, ax=axes.ravel().tolist(), label="surface temperature (°C)", shrink=0.75)
    fig.suptitle("Seasonal surface temperature by simulated solar longitude")
    save_map(fig, output / f"{args.label}_seasonal_surface_temperature.png")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), constrained_layout=True)
    mesh = panel(axes[0], lon, lat, annual_pr, land, "Annual mean precipitation rate", "YlGnBu", 0, 8)
    fig.colorbar(mesh, ax=axes[0], label="mm day⁻¹")
    # PlaSim code 182 is negative for upward evaporation, so P-E is pr+evap.
    p_minus_e = rmean(np.maximum(regular["pr"], 0.0) + regular["evap"]) * 86400.0 * 1000.0
    mesh = panel(axes[1], lon, lat, p_minus_e, land, "Annual mean precipitation minus evaporation", "BrBG", -4, 4)
    fig.colorbar(mesh, ax=axes[1], label="mm day⁻¹")
    save_map(fig, output / f"{args.label}_precipitation_and_water_balance.png")

    regular_season_masks = season_masks(phase_regular)
    fig, axes = plt.subplots(2, 2, figsize=(15, 8), constrained_layout=True)
    for ax, (name, _, _), mask in zip(axes.ravel(), SEASONS, regular_season_masks):
        field = climatology.masked_mean(np.maximum(regular["pr"], 0.0), rtime, mask) \
            * 86400.0 * 1000.0
        mesh = panel(ax, lon, lat, field, land, name, "YlGnBu", 0, 10)
    fig.colorbar(mesh, ax=axes.ravel().tolist(), label="precipitation rate (mm day⁻¹)", shrink=0.75)
    fig.suptitle("Seasonal precipitation by simulated solar longitude")
    save_map(fig, output / f"{args.label}_seasonal_precipitation.png")

    u = rmean(regular["ua"][:, -1])
    v = rmean(regular["va"][:, -1])
    speed = np.hypot(u, v)
    fig, ax = plt.subplots(figsize=(14, 6), constrained_layout=True)
    mesh = panel(ax, lon, lat, speed, land,
                 f"Annual circulation at lowest model level (σ={regular['lev'][-1]:.3f})",
                 "viridis", 0, 12)
    plon, pu = shifted(lon, u)
    _, pv = shifted(lon, v)
    skip = (slice(None, None, 4), slice(None, None, 4))
    ax.quiver(plon[::4], lat[::4], pu[skip], pv[skip], color="white", scale=230, width=0.002)
    fig.colorbar(mesh, ax=ax, label="wind speed (m s⁻¹)", orientation="horizontal", pad=0.09)
    save_map(fig, output / f"{args.label}_annual_lowest_level_winds.png")

    fig, axes = plt.subplots(2, 2, figsize=(15, 8), constrained_layout=True)
    fields = [
        (np.where(land >= 0.5, rmean(regular["snd"]) * 100.0, np.nan), "Mean land snow thickness", "Blues", 0, 100, "cm"),
        (np.where(land >= 0.5, snap["snd"].max(axis=0) * 100.0, np.nan), "Seasonal maximum snow thickness", "Blues", 0, 200, "cm"),
        (np.where(land < 0.5, rmean(regular["sic"]), np.nan), "Mean sea-ice concentration", "PuBu", 0, 1, "fraction"),
        (np.where(land < 0.5, snap["sic"].max(axis=0), np.nan), "Seasonal maximum sea-ice concentration", "PuBu", 0, 1, "fraction"),
    ]
    for ax, (field, title, cmap, low, high, units) in zip(axes.ravel(), fields):
        mesh = panel(ax, lon, lat, field, land, title, cmap, low, high)
        fig.colorbar(mesh, ax=ax, label=units)
    fig.suptitle("Snow and sea ice")
    save_map(fig, output / f"{args.label}_snow_and_sea_ice.png")

    fig, axes = plt.subplots(3, 1, figsize=(14, 14), constrained_layout=True)
    hydro = [
        (np.where(land >= 0.5, rmean(regular["mrso"]), np.nan), "Mean soil-water equivalent", "YlGnBu", 0, 0.4, "m"),
        (np.where(land >= 0.5, np.maximum(rmean(regular["mrro"]), 0.0) * 86400.0 * 1000.0, np.nan), "Mean river-routed net water flux", "Blues", 0, 5, "mm day⁻¹"),
        (np.where(land >= 0.5, -rmean(regular["evap"]) * 86400.0 * 1000.0, np.nan), "Mean upward land evaporation", "YlGn", 0, 5, "mm day⁻¹"),
    ]
    for ax, (field, title, cmap, low, high, units) in zip(axes, hydro):
        mesh = panel(ax, lon, lat, field, land, title, cmap, low, high)
        fig.colorbar(mesh, ax=ax, label=units)
    save_map(fig, output / f"{args.label}_hydrology.png")

    classes = koppen(regular["tas"], regular["pr"], phase_regular, lat, land, rtime)
    class_indices, class_names = categorical_map(
        classes, KOPPEN_COLORS, lon, lat, land,
        "Rate-normalized Köppen–Geiger climate interpretation", output / f"{args.label}_koppen_geiger.png",
    )
    biomes = np.vectorize(biome_for)(classes)
    biome_palette = {
        "Everwet tropical forest": "#075b21", "Monsoon tropical forest": "#168c35",
        "Seasonal tropical woodland/savanna": "#a6b735", "Desert": "#e8c681",
        "Steppe/semidesert": "#c7aa62", "Mediterranean woodland/shrubland": "#879b39",
        "Temperate/subtropical forest": "#45a05a", "Cool oceanic forest/heath": "#72a786",
        "Continental temperate/mixed forest": "#367b59", "Boreal forest/woodland": "#356354",
        "Tundra/alpine": "#9fa99b", "Permanent ice": "#eeeeee",
    }
    biome_indices, biome_names = categorical_map(
        biomes, biome_palette, lon, lat, land,
        "Broad ecological interpretation from simulated climate", output / f"{args.label}_biome_interpretation.png",
    )

    # The label already names the product; "baseline" was hardcoded here from
    # when there was only ever one, and a run labelled `baseline` produced
    # `baseline_baseline_classification.nc`.
    classification_nc = output / f"{args.label}_classification.nc"
    # netCDF4 1.7.4 emits a harmless NumPy 2.5 deprecation from assignment internals.
    with warnings.catch_warnings(), Dataset(classification_nc, "w", format="NETCDF4") as nc:
        warnings.filterwarnings(
            "ignore", message="Setting the shape on a NumPy array has been deprecated"
        )
        nc.createDimension("lat", len(lat)); nc.createDimension("lon", len(lon))
        nc.createVariable("lat", "f8", ("lat",))[:] = lat
        nc.createVariable("lon", "f8", ("lon",))[:] = lon
        koppen_var = nc.createVariable("koppen_class_index", "i2", ("lat", "lon"), zlib=True)
        koppen_var[:] = class_indices
        biome_var = nc.createVariable("biome_index", "i2", ("lat", "lon"), zlib=True)
        biome_var[:] = biome_indices
        temperature_var = nc.createVariable("annual_surface_temperature", "f4", ("lat", "lon"), zlib=True)
        temperature_var[:] = annual_ts
        temperature_var.units = "K"
        precipitation_var = nc.createVariable("annual_precipitation_rate", "f4", ("lat", "lon"), zlib=True)
        precipitation_var[:] = annual_pr
        precipitation_var.units = "mm day-1"
        elevation_var = nc.createVariable("surface_elevation", "f4", ("lat", "lon"), zlib=True)
        elevation_var[:] = elevation
        elevation_var.units = "m"
        nc.koppen_class_mapping = json.dumps({i + 1: name for i, name in enumerate(class_names)})
        nc.biome_mapping = json.dumps({i + 1: name for i, name in enumerate(biome_names)})
        nc.precipitation_normalization = "Koppen thresholds use rates annualized to 365.2425 days"
        nc.temperature_boundary = "0 degC coldest-period boundary between C and D"

    weights = leggauss(len(lat))[1][::-1, None] * np.ones_like(land)
    global_mean = lambda field: float(np.sum(field * weights) / np.sum(weights))
    ocean = land < 0.5
    report = {
        "label": args.label,
        "source_regular": str(args.regular.resolve()),
        "source_snapshots": str(args.snapshots.resolve()),
        "climatology_orbits": int(regular_attrs["climatology_orbit_count"]),
        "year_indices": [int(regular_attrs["climatology_start_year_index"]), int(regular_attrs["climatology_end_year_index"])],
        "global_metrics": {
            "surface_temperature_k": global_mean(annual_ts),
            "air_temperature_2m_k": global_mean(annual_tas),
            "precipitation_mm_day": global_mean(annual_pr),
            "toa_net_radiation_w_m2": global_mean(rmean(regular["ntr"])),
            "surface_downward_heat_flux_w_m2": global_mean(rmean(regular["hfns"])),
            "planetary_sea_ice_fraction": global_mean(rmean(regular["sic"])),
            "ocean_mean_sea_ice_concentration": float(np.sum(rmean(regular["sic"])[ocean] * weights[ocean]) / np.sum(weights[ocean])),
            "land_area_fraction": float(np.sum(weights[land >= 0.5]) / np.sum(weights)),
        },
        "koppen_land_area_fractions": area_fractions(classes, lat, land),
        "biome_land_area_fractions": area_fractions(biomes, lat, land),
        "classification_method": {
            "precipitation": "mean rates annualized to 365.2425 days before applying empirical Earth Koppen thresholds",
            "seasonality": "twelve model time-bin means assigned by reconstructed simulated solar longitude",
            "temperature_boundary_c_d_celsius": 0.0,
            "interpretation": "heuristic worldbuilding diagnostic, not a dynamic vegetation model",
        },
        "hydrology_sign_convention": "ExoPlaSim evap is negative upward; plotted evaporation is -evap and P-E is pr+evap",
        "season_definitions_solar_longitude_degrees": {
            name: [lower, upper] for name, lower, upper in SEASONS
        },
        "lowest_wind_sigma": float(regular["lev"][-1]),
        "snapshot_product": snapshot_attrs.get("climatology_product"),
    }
    # Named by label, like every other product of this pipeline. A fixed
    # filename meant each run silently overwrote the previous one's report, so
    # two climatologies could never coexist and the file's name told you nothing
    # about which world it described.
    (output / f"{args.label}_climate_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["global_metrics"], indent=2))


if __name__ == "__main__":
    main()
