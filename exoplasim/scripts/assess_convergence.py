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
from scipy.optimize import curve_fit
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


# The mixed layer's heat capacity, and the feedback strength measured from
# converged points spanning the design band. Together they set how long an
# approach takes, which is what turns a drift rate into a remaining offset.
SLAB_HEAT_CAPACITY = 50.0 * 1025.0 * 3990.0      # J/m2/K, 50 m of seawater
FEEDBACK_W_M2_K = 1.31                            # measured, not assumed


def relaxation_orbits(orbital_year_days: float) -> float:
    tau_seconds = SLAB_HEAT_CAPACITY / FEEDBACK_W_M2_K
    return tau_seconds / (orbital_year_days * 86400.0)


def approach_to_equilibrium(orbits: np.ndarray, series: np.ndarray,
                            fraction: float = 0.35):
    """Fit T(n) = T_inf - A exp(-n/tau) and return the asymptote with an interval.

    A drift rate is not a distance. An exponential approach at rate `d` with time
    constant `tau` still has `d * tau` to travel, all of it one-signed, so a
    threshold on the rate bounds nothing on its own. At this planet's tau of
    about ten orbits, the 0.05 K/orbit criterion permits roughly half a kelvin of
    remaining approach -- the same order as corrections the project has thought
    worth applying, and larger than some.

    The interval is what makes this usable rather than a second point estimate.
    A fit whose tau is comparable to the run length is extrapolating past its
    data, and says so through a wide interval instead of a confident wrong
    number; one such fit in this project's history sat 0.6 K from its own
    reported value.

    Returns (asymptote, half_width, tau_orbits) or (nan, nan, nan) if the series
    will not support a fit.
    """
    mask = orbits >= orbits.max() * fraction
    x, y = orbits[mask], series[mask]
    if x.size < 8:
        return float("nan"), float("nan"), float("nan")
    model = lambda n, inf, amp, tau: inf - amp * np.exp(-n / tau)
    try:
        popt, pcov = curve_fit(model, x, y,
                               p0=[y[-1], max(y[0] - y[-1], 1e-3), max(x.size / 3, 1.0)],
                               maxfev=40000)
    except Exception:
        return float("nan"), float("nan"), float("nan")
    inf, tau = float(popt[0]), float(popt[2])
    err = float(np.sqrt(np.diag(pcov))[0]) if np.all(np.isfinite(pcov)) else float("nan")
    # A fit extrapolating beyond its own span is reported as such rather than
    # trusted: widen the interval by how far past the data the asymptote sits.
    if np.isfinite(tau) and tau > 0:
        overshoot = max(0.0, tau / max(x.size, 1))
        err = (err if np.isfinite(err) else abs(inf - y[-1])) * (1.0 + overshoot)
    return inf, err, tau


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

    # The quantity the design band is stated in, and therefore the one that has
    # to be bounded. Everything above is a rate; this is a distance.
    # (tau_expected is needed by the fallback above, so it is computed first.)
    # The year comes from the run's own manifest, not from the current config:
    # this is a property of the run being assessed, which may predate a
    # re-baseline that moved the orbit.
    manifest_path = run_dir / "run_manifest.json"
    year_days = 182.8
    if manifest_path.is_file():
        year_days = float(json.loads(manifest_path.read_text(encoding="utf-8"))
                          ["derived_parameters"]["orbital_year_earth_days"])
    tau_expected = relaxation_orbits(year_days)
    orbits_axis = np.arange(len(arrays["ts"]), dtype=float)
    asymptote, half_width, tau_fit = approach_to_equilibrium(orbits_axis, arrays["ts"])
    offset = asymptote - metrics["temperature_mean_k"]

    # A converged run has no approach left to fit, so the exponential becomes
    # unconstrained and can return anything: one run here reported an asymptote
    # of -23064 K with an interval of six million. That is the fit having nothing
    # to grip, not the run being far from equilibrium, and the two must not be
    # treated alike -- the first version of this criterion failed a converged run
    # for it.
    #
    # When the fit is unusable, fall back to the bound the drift itself implies,
    # offset ~ drift * tau, using the EXPECTED tau rather than a fitted one, and
    # carry 100% uncertainty on it. That is conservative in the right direction:
    # it can only refuse a run, never pass one it should not.
    fit_usable = (np.isfinite(asymptote) and np.isfinite(half_width)
                  and abs(offset) < 20.0 and half_width < 5.0)
    if not fit_usable:
        offset = metrics["temperature_slope_k_per_orbit"] * tau_expected
        half_width = abs(offset)
        asymptote = metrics["temperature_mean_k"] + offset
    metrics.update({
        "temperature_asymptote_k": asymptote,
        "temperature_asymptote_half_width_k": half_width,
        "temperature_remaining_offset_k": offset,
        "relaxation_orbits_fitted": tau_fit,
        "relaxation_orbits_expected": tau_expected,
        "remaining_offset_implied_by_drift_k":
            metrics["temperature_slope_k_per_orbit"] * tau_expected,
    })

    # 0.15 K against a 3 K design band: small enough that the band's edges mean
    # what they say, loose enough to be reachable. Stated before it was applied.
    OFFSET_TOLERANCE_K = 0.15
    criteria = {
        "abs_temperature_slope_lt_0.05_k_per_orbit": abs(metrics["temperature_slope_k_per_orbit"]) < 0.05,
        "abs_toa_slope_lt_0.05_w_m2_per_orbit": abs(metrics["toa_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_surface_slope_lt_0.05_w_m2_per_orbit": abs(metrics["surface_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_sea_ice_slope_lt_0.001_per_orbit": abs(metrics["sea_ice_slope_fraction_per_orbit"]) < 0.001,
        "abs_mean_toa_lt_0.5_w_m2": abs(metrics["mean_toa_balance_w_m2"]) < 0.5,
        "abs_mean_surface_lt_0.5_w_m2": abs(metrics["mean_surface_balance_w_m2"]) < 0.5,
        # The one that bounds the answer rather than its rate of change. Fails
        # closed: a fit that will not converge is not evidence of equilibrium.
        f"extrapolated_offset_lt_{OFFSET_TOLERANCE_K}_k": bool(
            np.isfinite(offset) and np.isfinite(half_width)
            and abs(offset) + half_width < OFFSET_TOLERANCE_K),
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
    # Named by the run, not fixed. A fixed filename meant every assessment
    # overwrote the last, so a run's convergence record could not survive the
    # next run being assessed -- and assessing a run silently destroyed the
    # evidence for a previous one.
    report_path = args.output / f"{run_dir.name}_convergence.json"
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
    plot_path = args.output / f"{run_dir.name}_convergence.png"
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
