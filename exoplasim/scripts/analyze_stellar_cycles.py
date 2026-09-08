#!/usr/bin/env python3
"""Analyze periodic response in the paired ExoPlaSim stellar-cycle runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

from _paths import ANALYSIS, RUNS
import gridding  # noqa: E402  from lib/, via _paths: the one Gaussian quadrature


OUTDIR = ANALYSIS / "stellar_cycles"
FIELDS = ["ts", "ntr", "hfns", "sic", "pr", "snd"]


def global_mean(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.sum(field * weights[None, :, None], axis=(-2, -1)) / (
        2.0 * field.shape[-1]
    )


def load_run(run_dir: Path, fold_on: str | None = None) -> dict:
    """Load one cycle run, phase-folded on a single named component.

    The forcing carries every component, but a phase fold is defined against one
    period at a time -- with non-commensurate periods there is no single phase
    that describes both. `fold_on` names which; the default is the longest,
    since that is the one a composite needs the most orbits to resolve.
    """
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    if manifest["status"] != "simulation_complete":
        raise RuntimeError(f"Run is not complete: {run_dir}")
    cycle = manifest["cycle"]
    components = cycle["components"]
    if fold_on is None:
        fold_on = max(components, key=lambda n: components[n]["period_earth_years"])
    if fold_on not in components:
        raise RuntimeError(
            f"{run_dir.name} has no cycle component {fold_on!r}; "
            f"it has {sorted(components)}"
        )
    orbit_steps = int(manifest["fixed_orbit_parameters"]["runsteps_per_orbit"])
    period_steps = float(components[fold_on]["period_model_steps"])
    files = sorted(run_dir.glob("MOST.[0-9][0-9][0-9][0-9][0-9].nc"))
    if len(files) != int(manifest["completed_orbits"]):
        raise RuntimeError(f"Output count disagrees with manifest in {run_dir}")

    annual = {name: [] for name in FIELDS}
    series = {name: [] for name in FIELDS}
    series_age = []
    composite_sums = {
        state: {name: None for name in ["ts", "sic", "pr"]}
        for state in ["bright", "dim"]
    }
    composite_counts = {"bright": 0, "dim": 0}
    lon = lat = weights = None
    for path in files:
        year_index = int(path.name.split(".")[1])
        with Dataset(path) as nc:
            if weights is None:
                lat = np.asarray(nc["lat"][:], dtype=float)
                lon = np.asarray(nc["lon"][:], dtype=float)
                weights = gridding.gaussian_row_weights(
                    lat, what=f"{path}'s latitude axis")
            for name in FIELDS:
                values = np.asarray(nc[name][:], dtype=float)
                gm = global_mean(values, weights)
                if name == "pr":
                    gm *= 86400.0 * 1000.0
                series[name].extend(gm.tolist())
                annual[name].append(float(gm.mean()))

            times = np.asarray(nc["time"][:], dtype=float)
            ages = (year_index * orbit_steps + times) / period_steps
            series_age.extend(ages.tolist())
            phases = ages % 1.0
            final = (ages >= 2.0) & (ages < 4.0)
            for state, center in [("bright", 0.25), ("dim", 0.75)]:
                selected = final & (circular_phase_distance(phases, center) <= 0.125)
                if not selected.any():
                    continue
                composite_counts[state] += int(selected.sum())
                for name in composite_sums[state]:
                    values = np.asarray(nc[name][:], dtype=float)[selected]
                    if name == "pr":
                        values *= 86400.0 * 1000.0
                    subtotal = values.sum(axis=0)
                    if composite_sums[state][name] is None:
                        composite_sums[state][name] = subtotal
                    else:
                        composite_sums[state][name] += subtotal

    for name in annual:
        annual[name] = np.asarray(annual[name])
        series[name] = np.asarray(series[name])
    series_age = np.asarray(series_age)
    series_phase = series_age % 1.0
    start_step = float(cycle["start_model_step"])
    series_orbit_phase = (
        (start_step + series_age * period_steps) / orbit_steps
    ) % 1.0
    composites = {
        state: {
            name: composite_sums[state][name] / composite_counts[state]
            for name in composite_sums[state]
        }
        for state in composite_sums
    }
    orbit_number = np.arange(len(files), dtype=float) + 0.5
    age_steps = orbit_number * orbit_steps
    cycle_age = age_steps / period_steps
    phase = cycle_age % 1.0
    # Every component contributes at its own period, exactly as radmod sums them:
    # flux = mean + sum_i amp_i * sin(2pi * (age/period_i + phase_i)). Folding on
    # one component does not remove the others from the forcing, so reconstruct
    # from all of them rather than from the folded one alone.
    forcing = np.full(age_steps.shape, float(cycle["mean_flux_earth"]))
    for spec in components.values():
        forcing += float(spec["semi_amplitude_flux_earth"]) * np.sin(
            2.0
            * np.pi
            * (
                age_steps / float(spec["period_model_steps"])
                + float(spec.get("phase_cycles", 0.0))
            )
        )
    return {
        "manifest": manifest,
        "fold_on": fold_on,
        "cycle_age": cycle_age,
        "phase": phase,
        "forcing": forcing,
        "annual": annual,
        "series_age": series_age,
        "series_phase": series_phase,
        "series_orbit_phase": series_orbit_phase,
        "series": series,
        "composites": composites,
        "composite_counts": composite_counts,
        "lat": lat,
        "lon": lon,
        "weights": weights,
    }


def phase_curve(data: dict, cycle_index: int, name: str, phases: np.ndarray) -> np.ndarray:
    x = data["cycle_age"]
    y = data["annual"][name]
    sample = cycle_index + phases
    return np.interp(sample, x, y)


def circular_phase_distance(a: np.ndarray, b: float) -> np.ndarray:
    return np.abs((a - b + 0.5) % 1.0 - 0.5)


def run_report(label: str, data: dict) -> dict:
    phases = np.linspace(0.01, 0.99, 99)
    comparisons = {}
    for name in ["ts", "sic", "pr", "ntr", "hfns"]:
        c3 = phase_curve(data, 2, name, phases)
        c4 = phase_curve(data, 3, name, phases)
        comparisons[name] = {
            "cycle_3_mean": float(c3.mean()),
            "cycle_4_mean": float(c4.mean()),
            "cycle_4_minus_3_mean": float((c4 - c3).mean()),
            "cycle_4_vs_3_rmse": float(np.sqrt(np.mean((c4 - c3) ** 2))),
        }

    final = (data["series_age"] >= 2.0) & (data["series_age"] < 4.0)
    last_cycle = (data["series_age"] >= 3.0) & (data["series_age"] < 4.0)
    x = data["series_age"][last_cycle]
    temperature = data["series"]["ts"][last_cycle]
    warm_i = int(np.argmax(temperature))
    warm_phase = float(x[warm_i] % 1.0)
    harmonic_x = data["series_age"][final]
    orbital_x = data["series_orbit_phase"][final]
    columns = [np.ones(harmonic_x.size)]
    for harmonic in range(1, 4):
        columns.extend([
            np.sin(2.0 * np.pi * harmonic * harmonic_x),
            np.cos(2.0 * np.pi * harmonic * harmonic_x),
        ])
    for harmonic in range(1, 5):
        columns.extend([
            np.sin(2.0 * np.pi * harmonic * orbital_x),
            np.cos(2.0 * np.pi * harmonic * orbital_x),
        ])
    design = np.column_stack(columns)
    coefficients = np.linalg.lstsq(
        design, data["series"]["ts"][final], rcond=None
    )[0]
    phase_grid = np.linspace(0.0, 1.0, 2001)
    star_component = np.full(phase_grid.shape, coefficients[0])
    coefficient_index = 1
    for harmonic in range(1, 4):
        star_component += (
            coefficients[coefficient_index]
            * np.sin(2.0 * np.pi * harmonic * phase_grid)
            + coefficients[coefficient_index + 1]
            * np.cos(2.0 * np.pi * harmonic * phase_grid)
        )
        coefficient_index += 2
    harmonic_peak_phase = float(phase_grid[np.argmax(star_component)])
    lag_phase = float((harmonic_peak_phase - 0.25 + 0.5) % 1.0 - 0.5)
    _comps = data["manifest"]["cycle"]["components"]
    period_years = float(_comps[data["fold_on"]]["period_earth_years"])
    result = {
        "label": label,
        "fold_on": data["fold_on"],
        "cycles_completed": float(
            data["manifest"]["completed_cycles"][data["fold_on"]]),
        "periodic_comparison_cycles_3_and_4": comparisons,
        "final_two_cycles_binned": {},
        "warmest_raw_bin_phase": warm_phase,
        "harmonic_temperature_peak_phase": harmonic_peak_phase,
        "seasonally_adjusted_temperature_peak_to_peak_k": float(
            star_component.max() - star_component.min()
        ),
        "temperature_lag_after_flux_max_phase": lag_phase,
        "temperature_lag_earth_years": lag_phase * period_years,
        "bright_composite_samples": data["composite_counts"]["bright"],
        "dim_composite_samples": data["composite_counts"]["dim"],
        "bright_minus_dim_composites": {},
        "final_two_cycles_orbit_means": {},
    }
    for name in ["ts", "sic", "pr"]:
        difference = (
            data["composites"]["bright"][name]
            - data["composites"]["dim"][name]
        )
        result["bright_minus_dim_composites"][name] = {
            "area_weighted_mean": float(
                global_mean(difference[None, :, :], data["weights"])[0]
            ),
            "minimum": float(difference.min()),
            "maximum": float(difference.max()),
            "fifth_percentile": float(np.percentile(difference, 5.0)),
            "ninety_fifth_percentile": float(np.percentile(difference, 95.0)),
        }
    for name in ["ts", "sic", "pr", "ntr", "hfns"]:
        values = data["series"][name][final]
        result["final_two_cycles_binned"][name] = {
            "mean": float(values.mean()),
            "minimum": float(values.min()),
            "maximum": float(values.max()),
            "peak_to_peak": float(values.max() - values.min()),
        }
        orbit_final = (data["cycle_age"] >= 2.0) & (data["cycle_age"] < 4.0)
        orbit_values = data["annual"][name][orbit_final]
        result["final_two_cycles_orbit_means"][name] = {
            "mean": float(orbit_values.mean()),
            "minimum": float(orbit_values.min()),
            "maximum": float(orbit_values.max()),
            "peak_to_peak": float(orbit_values.max() - orbit_values.min()),
        }
    return result


def plot_timeseries(all_data: dict[str, dict]) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True, constrained_layout=True)
    colors = ["#9b2226", "#005f73"]
    for (label, data), color in zip(all_data.items(), colors, strict=True):
        x = data["series_age"]
        # The multi-component forcing load_run already reconstructed, on orbit
        # centres.
        axes[0].plot(data["cycle_age"], data["forcing"], color=color,
                     label=label)
        axes[1].plot(x, data["series"]["ts"], color=color, lw=0.8)
        axes[2].plot(x, data["series"]["sic"], color=color, lw=0.8)
        axes[3].plot(x, data["series"]["pr"], color=color, lw=0.8)
    axes[0].set_ylabel("Stellar flux\n($S_\\oplus$)")
    axes[1].set_ylabel("Surface T\n(K)")
    axes[2].set_ylabel("Sea-ice\nfraction")
    axes[3].set_ylabel("Precipitation\n(mm day$^{-1}$)")
    axes[3].set_xlabel("Elapsed cycles of the folded component")
    axes[0].legend(loc="best")
    for ax in axes:
        ax.grid(alpha=0.25)
        for boundary in range(1, 5):
            ax.axvline(boundary, color="0.6", lw=0.7, ls="--")
    fig.suptitle("Climate response to the K-dwarf starspot cycle")
    fig.savefig(OUTDIR / "stellar_cycle_timeseries.png", dpi=180)
    plt.close(fig)


def plot_phase_response(all_data: dict[str, dict]) -> None:
    phases = np.linspace(0.0, 1.0, 101)
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True, constrained_layout=True)
    colors = ["#9b2226", "#005f73"]
    for (label, data), color in zip(all_data.items(), colors, strict=True):
        for cycle_index, alpha, ls in [(2, 0.45, "--"), (3, 1.0, "-")]:
            suffix = "cycle 3" if cycle_index == 2 else "cycle 4"
            axes[0].plot(
                phases, phase_curve(data, cycle_index, "ts", phases),
                color=color, alpha=alpha, ls=ls, label=f"{label}, {suffix}",
            )
            axes[1].plot(
                phases, phase_curve(data, cycle_index, "sic", phases),
                color=color, alpha=alpha, ls=ls,
            )
            axes[2].plot(
                phases, phase_curve(data, cycle_index, "pr", phases),
                color=color, alpha=alpha, ls=ls,
            )
    axes[0].set_ylabel("Surface T (K)")
    axes[1].set_ylabel("Sea-ice fraction")
    axes[2].set_ylabel("Precipitation\n(mm day$^{-1}$)")
    axes[2].set_xlabel("Stellar-cycle phase (flux maximum at 0.25)")
    axes[0].legend(fontsize=8, ncol=2)
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.axvline(0.25, color="#e9c46a", lw=1.0)
        ax.axvline(0.75, color="#264653", lw=1.0)
    fig.savefig(OUTDIR / "stellar_cycle_phase_response.png", dpi=180)
    plt.close(fig)


def plot_composite_maps(all_data: dict[str, dict]) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), constrained_layout=True)
    for row, (label, data) in enumerate(all_data.items()):
        delta_t = data["composites"]["bright"]["ts"] - data["composites"]["dim"]["ts"]
        ice_delta = data["composites"]["bright"]["sic"] - data["composites"]["dim"]["sic"]
        precip_delta = data["composites"]["bright"]["pr"] - data["composites"]["dim"]["pr"]
        fields = [delta_t, ice_delta, precip_delta]
        titles = ["Temperature (K)", "Sea-ice fraction", "Precip. (mm/day)"]
        ranges = [(-12, 12), (-0.8, 0.8), (-3, 3)]
        for col, (field, title, limits) in enumerate(zip(fields, titles, ranges, strict=True)):
            image = axes[row, col].pcolormesh(
                data["lon"], data["lat"], field,
                shading="auto", cmap="RdBu_r", vmin=limits[0], vmax=limits[1],
            )
            axes[row, col].set_title(f"{label}\nbright − dim: {title}")
            axes[row, col].set_xlabel("Longitude")
            axes[row, col].set_ylabel("Latitude")
            fig.colorbar(image, ax=axes[row, col], shrink=0.82)
    fig.savefig(OUTDIR / "stellar_cycle_bright_minus_dim_maps.png", dpi=180)
    plt.close(fig)


def label_for(run_dir: Path, data: dict) -> str:
    """Name a run by what it physically is, since run ids are UUIDs."""
    cycle = data["manifest"]["cycle"]
    mean = float(cycle["mean_flux_earth"])
    ptp = float(cycle["total_amplitude_flux_peak_to_peak"])
    return f"{run_dir.name} ({mean - ptp / 2:.3f}-{mean + ptp / 2:.3f} S-Earth)"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Phase-fold and compare completed stellar-cycle runs.",
    )
    parser.add_argument(
        "run",
        nargs="+",
        help="run ids under exoplasim/runs/ (ids are UUIDs -- see index_runs.py)",
    )
    parser.add_argument(
        "--fold-on",
        default=None,
        help="cycle component to phase-fold on; default is the longest period",
    )
    args = parser.parse_args()

    run_dirs = []
    for run_id in args.run:
        run_dir = RUNS / run_id
        if not (run_dir / "run_manifest.json").is_file():
            raise SystemExit(
                f"no run manifest at {run_dir}. Run ids are UUIDs and carry no "
                f"meaning; list what exists with index_runs.py."
            )
        run_dirs.append(run_dir)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    loaded = [(d, load_run(d, args.fold_on)) for d in run_dirs]
    all_data = {label_for(d, data): data for d, data in loaded}
    report = {label: run_report(label, data) for label, data in all_data.items()}
    plot_timeseries(all_data)
    plot_phase_response(all_data)
    plot_composite_maps(all_data)
    (OUTDIR / "stellar_cycle_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
