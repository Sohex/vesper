#!/usr/bin/env python3
"""Assess ExoPlaSim spin-up convergence from validated annual NetCDF files."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

_MPL_CACHE = Path("/tmp/world-matplotlib-cache")
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS


def output_files(run_dir: Path) -> list[Path]:
    outputs = []
    for path in run_dir.glob("MOST.*.nc"):
        match = re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)
        if match:
            outputs.append((int(match.group(1)), path))
    outputs.sort()
    if [year for year, _ in outputs] != list(range(len(outputs))):
        raise RuntimeError("Annual outputs are not contiguous from year zero")
    return [path for _, path in outputs]


def global_mean(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.sum(field * weights[None, :, None], axis=(-2, -1)) / (2.0 * field.shape[-1])


def slope(series: np.ndarray, window: int) -> float:
    y = series[-window:]
    x = np.arange(window, dtype=float)
    return float(np.polyfit(x, y, 1)[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--output", type=Path, default=ANALYSIS / "convergence")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    files = output_files(run_dir)
    if len(files) < args.window:
        raise RuntimeError("Not enough annual outputs for requested window")

    records = []
    for year, path in enumerate(files):
        with Dataset(path) as nc:
            weights = leggauss(len(nc.dimensions["lat"]))[1][::-1]
            record = {"year_index": year}
            for name in ["ts", "ntr", "hfns", "sic", "pr"]:
                field = np.asarray(nc[name][:], dtype=float)
                value = float(global_mean(field, weights).mean())
                if name == "pr":
                    value *= 86400.0 * 1000.0
                record[name] = value
            records.append(record)

    arrays = {key: np.array([record[key] for record in records]) for key in ["ts", "ntr", "hfns", "sic", "pr"]}
    w = args.window
    metrics = {
        "temperature_slope_k_per_orbit": slope(arrays["ts"], w),
        "toa_balance_slope_w_m2_per_orbit": slope(arrays["ntr"], w),
        "surface_balance_slope_w_m2_per_orbit": slope(arrays["hfns"], w),
        "sea_ice_slope_fraction_per_orbit": slope(arrays["sic"], w),
        "mean_toa_balance_w_m2": float(arrays["ntr"][-w:].mean()),
        "mean_surface_balance_w_m2": float(arrays["hfns"][-w:].mean()),
        "temperature_mean_k": float(arrays["ts"][-w:].mean()),
        "temperature_last_k": float(arrays["ts"][-1]),
        "sea_ice_mean_fraction": float(arrays["sic"][-w:].mean()),
    }
    criteria = {
        "abs_temperature_slope_lt_0.05_k_per_orbit": abs(metrics["temperature_slope_k_per_orbit"]) < 0.05,
        "abs_toa_slope_lt_0.05_w_m2_per_orbit": abs(metrics["toa_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_surface_slope_lt_0.05_w_m2_per_orbit": abs(metrics["surface_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_sea_ice_slope_lt_0.001_per_orbit": abs(metrics["sea_ice_slope_fraction_per_orbit"]) < 0.001,
        "abs_mean_toa_lt_0.5_w_m2": abs(metrics["mean_toa_balance_w_m2"]) < 0.5,
        "abs_mean_surface_lt_0.5_w_m2": abs(metrics["mean_surface_balance_w_m2"]) < 0.5,
    }
    report = {
        "run_dir": str(run_dir),
        "completed_orbits": len(files),
        "window_orbits": w,
        "metrics": metrics,
        "criteria": criteria,
        "sufficiently_equilibrated_for_worldbuilding": bool(all(criteria.values())),
        "annual_records": records,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    report_path = args.output / "baseline_convergence.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if report["sufficiently_equilibrated_for_worldbuilding"]:
        manifest_path = run_dir / "run_manifest.json"
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["status"] = "equilibrated_for_worldbuilding"
            manifest["equilibrium_cutoff_year_index"] = len(files) - 1
            manifest["convergence_assessment"] = {
                "report": str(report_path.resolve()),
                "window_orbits": w,
                "metrics": metrics,
                "criteria": criteria,
            }
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )

    years = np.arange(len(files))
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    panels = [
        ("ts", "Global surface temperature", "K"),
        ("ntr", "TOA net radiation", "W m$^{-2}$"),
        ("hfns", "Surface downward heat flux", "W m$^{-2}$"),
        ("sic", "Planetary sea-ice fraction", "fraction"),
    ]
    for ax, (key, title, units) in zip(axes.ravel(), panels):
        ax.plot(years, arrays[key], marker="o", markersize=2.5, linewidth=1.1)
        ax.axvspan(years[-w], years[-1], color="tab:orange", alpha=0.12, label=f"{w}-orbit test window")
        if key in {"ntr", "hfns"}:
            ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(title)
        ax.set_xlabel("orbit index")
        ax.set_ylabel(units)
        ax.grid(alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(
        f"Baseline spin-up convergence: {'PASS' if all(criteria.values()) else 'NOT YET'} "
        f"({len(files)} orbits, {w}-orbit window)"
    )
    plot_path = args.output / "baseline_convergence.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    payload = {"metrics": metrics, "criteria": criteria,
               "pass": all(criteria.values()), "orbits": len(files),
               "window_orbits": w}

    # Record the verdict with the run as well as in the analysis directory.
    # Provenance travels with the artifact everywhere else in this project, and a
    # convergence result that lives only in an output folder cannot be found from
    # the run it describes.
    manifest_path = args.run_dir / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["convergence_assessment"] = payload
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
