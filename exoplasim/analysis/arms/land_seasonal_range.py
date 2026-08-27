#!/usr/bin/env python
"""The land seasonal surface-temperature range of an arm, against its control.

    python exoplasim/analysis/arms/land_seasonal_range.py \\
        --control exoplasim/runs/run_xxxxxxxxxxxx --out sizing.json
    python exoplasim/analysis/arms/land_seasonal_range.py \\
        --control CONTROL --arm ARM --restart MOST_REST.000nn \\
        --window-start 15 --out result.json

WORLDBUILDING CONTEXT: Vesper is a fictional super-Earth. Everything below is a
diagnostic of the climate model that simulates it, over its simulated land.

WHY THIS EXISTS BESIDE `compare_forcing_arms.py`. That one reads an arm on
GLOBAL MEANS, which is the right instrument for a term whose entry in the bundle
sum is a top-of-atmosphere flux. The soil heat solver's thermal pair has no
top-of-atmosphere forcing at all -- a closed column's annual-mean heat flux is
zero whatever its capacity and conductivity are -- so its registered prediction
is an AMPLITUDE: the seasonal peak-to-trough range of the modelled surface
temperature over land, which a global mean of anything cannot see.

WHAT IS MEASURED, and each of these answers one registered "what would mean
wrong" in `exoplasim/notes/forcing-bundle-predictions.md`:

- The per-orbit, per-cell seasonal range of `ts`, max minus min over the orbit's
  own output bins, area-weighted over the model's LAND cells on the Gaussian
  latitudes. Its arm-minus-control difference carries the paired difference's
  own autocorrelation-corrected standard error, from `lib/autocorrelation.py`.
- The share of land cells whose range FALLS. Lowering a column's thermal inertia
  cannot damp its surface more, so a cell that damps more is an implementation
  fault and not a small result.
- The correlation between a cell's fractional range change and its soil water
  DEFICIT. The relation the solver runs is a function of the column's own water,
  so a change uncorrelated with wetness says the water is not being read.
  `dwmax` comes from the restart both arms branched from, because a restarted
  run takes the soil water capacity from the restart and not from the staged
  field.
- The seasonal amplitude of each soil TEMPERATURE layer the postprocessor emits,
  and the bin its maximum falls in. A lower diffusivity shortens the seasonal
  penetration depth, so the deep layers' amplitude ratio and phase lag are where
  that shows.

WITH NO `--arm` this reports the INSTRUMENT alone: the control's own orbit-to-
orbit scatter in the land-mean range, which is what a resolution bar has to be
fixed against before an arm is bought (`docs/src/practice/failure-modes.md`
class 34).
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

REPO = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "exoplasim" / "scripts"))
import lib.autocorrelation as ac  # noqa: E402


def orbit_files(run_dir: pathlib.Path) -> list[pathlib.Path]:
    found = []
    for path in run_dir.glob("MOST.*.nc"):
        m = re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)
        if m:
            found.append((int(m.group(1)), path))
    found.sort()
    return [p for _, p in found]


def gaussian_weights(nlat: int) -> np.ndarray:
    """Area weights for the model's own Gaussian latitudes, north to south."""
    nodes, weights = np.polynomial.legendre.leggauss(nlat)
    return weights[np.argsort(-nodes)]


SOIL_LEVELS = ("tsod", "tso2", "tso3", "tso4")


def read_run(run_dir: pathlib.Path) -> dict:
    """Per-orbit fields this comparison is taken on.

    `range_cell` is (orbits, nlat, nlon) and is the seasonal peak-to-trough of
    `ts` inside each orbit. The land mask is read from the model's own `lsm`
    rather than reconstructed, because what the soil solver applies to is the
    model's land and not the export's.
    """
    files = orbit_files(run_dir)
    if not files:
        raise SystemExit(f"{run_dir} holds no MOST.*.nc")
    ranges, wetness, soil_amp, soil_phase = [], [], [], []
    land = None
    for path in files:
        with Dataset(path) as nc:
            ts = np.asarray(nc.variables["ts"][:], dtype=float)
            if land is None:
                lsm = np.asarray(nc.variables["lsm"][:], dtype=float)
                land = (lsm[0] if lsm.ndim == 3 else lsm) > 0.5
            ranges.append(ts.max(axis=0) - ts.min(axis=0))
            wetness.append(np.asarray(nc.variables["mrso"][:], dtype=float).mean(axis=0))
            amps, phases = {}, {}
            for name in SOIL_LEVELS:
                if name not in nc.variables:
                    continue
                field = np.asarray(nc.variables[name][:], dtype=float)
                amps[name] = field.max(axis=0) - field.min(axis=0)
                phases[name] = np.argmax(field, axis=0).astype(float)
            soil_amp.append(amps)
            soil_phase.append(phases)
    return {"run": run_dir.name, "orbits": len(files), "land": land,
            "range_cell": np.stack(ranges), "mrso_cell": np.stack(wetness),
            "soil_amp": soil_amp, "soil_phase": soil_phase,
            "weights": gaussian_weights(land.shape[0])}


def land_mean(field: np.ndarray, land: np.ndarray, weights: np.ndarray) -> float:
    """Area-weighted mean over land. Longitude is uniform; latitude is not."""
    w = np.broadcast_to(weights[:, None], field.shape) * land
    return float((field * w).sum() / w.sum())


def series_land_mean(run: dict) -> np.ndarray:
    return np.array([land_mean(r, run["land"], run["weights"])
                     for r in run["range_cell"]])


def window_stats(diff: np.ndarray) -> dict:
    """Mean of a paired difference with the standard error its memory allows."""
    if diff.size < 4:
        return {"orbits": int(diff.size), "mean": float(diff.mean()) if diff.size else None,
                "tau_orbits": None, "standard_error": None}
    tau = float(ac.integrated_time(diff).get("tau", float("nan")))
    sem = (float(ac.mean_standard_error(diff, tau))
           if np.isfinite(tau) and tau > 0 else None)
    return {"orbits": int(diff.size), "mean": float(diff.mean()),
            "scatter": float(diff.std(ddof=1)), "tau_orbits": tau,
            "standard_error": sem}


def restart_dwmax(path: pathlib.Path, shape: tuple[int, int]) -> np.ndarray:
    import restart_format
    raw = restart_format.payloads(path)["dwmax"]
    return np.frombuffer(raw, dtype="<f8").reshape(shape)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--control", type=pathlib.Path, required=True)
    ap.add_argument("--arm", type=pathlib.Path)
    ap.add_argument("--restart", type=pathlib.Path,
                    help="the restart both arms branched from; `dwmax` is read "
                         "from it, because a restarted run takes the soil water "
                         "capacity from the restart and not from the staged field")
    ap.add_argument("--window-start", type=int, default=0,
                    help="first orbit of the measurement window; everything "
                         "before it is the settling block")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()

    control = read_run(args.control.resolve())
    c_series = series_land_mean(control)
    w0 = args.window_start
    payload = {
        "provenance": {
            "generator": "exoplasim/analysis/arms/land_seasonal_range.py",
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                         capture_output=True, text=True).stdout.strip(),
            "python": sys.version.split()[0],
        },
        "control": {
            "run": control["run"], "orbits": control["orbits"],
            "land_cells": int(control["land"].sum()),
            "land_mean_range_k_per_orbit": [float(x) for x in c_series],
            "window_starts_at_orbit": w0,
            "window_land_mean_range_k": float(c_series[w0:].mean()),
            # The instrument's own scatter, which is what a resolution bar on an
            # unpaired quantity has to be fixed against.
            "window_scatter_k": (float(c_series[w0:].std(ddof=1))
                                 if c_series[w0:].size > 1 else None),
            **{f"window_{k}": v for k, v in window_stats(c_series[w0:]).items()},
        },
    }

    if args.arm is not None:
        arm = read_run(args.arm.resolve())
        a_series = series_land_mean(arm)
        n = min(arm["orbits"], control["orbits"])
        diff = a_series[:n] - c_series[:n]
        stats = window_stats(diff[w0:n])
        sem = stats.get("standard_error")
        # The bar registered with the prediction, and the same one `vdiff_lamm`
        # was reported under: two arms each carrying their own error, so the
        # separation has to clear sqrt(2) of the larger.
        bar = 2.0 * np.sqrt(2.0) * sem if sem else None
        c_mean = control["range_cell"][w0:n].mean(axis=0)
        a_mean = arm["range_cell"][w0:n].mean(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.where(c_mean > 0, (a_mean - c_mean) / c_mean, np.nan)
        land = control["land"]
        falls = int(((frac < 0) & land).sum())

        correlation = None
        if args.restart is not None:
            dwmax = restart_dwmax(args.restart.resolve(), land.shape)
            with np.errstate(invalid="ignore", divide="ignore"):
                zf = np.clip(np.where(dwmax > 0,
                                      control["mrso_cell"][w0:n].mean(axis=0) / dwmax,
                                      np.nan), 0.0, 1.0)
            good = land & np.isfinite(zf) & np.isfinite(frac)
            deficit = 1.0 - zf[good]
            correlation = {
                "cells": int(good.sum()),
                "pearson_r_frac_change_vs_deficit":
                    float(np.corrcoef(deficit, frac[good])[0, 1]),
                "spearman_r_frac_change_vs_deficit":
                    float(np.corrcoef(
                        np.argsort(np.argsort(deficit)).astype(float),
                        np.argsort(np.argsort(frac[good])).astype(float))[0, 1]),
                "degree_of_saturation": {
                    "min": float(zf[good].min()), "median": float(np.median(zf[good])),
                    "max": float(zf[good].max())},
            }

        soil = {}
        for name in SOIL_LEVELS:
            if name not in control["soil_amp"][0] or name not in arm["soil_amp"][0]:
                continue
            c_amp = np.mean([a[name] for a in control["soil_amp"][w0:n]], axis=0)
            a_amp = np.mean([a[name] for a in arm["soil_amp"][w0:n]], axis=0)
            c_ph = np.mean([p[name] for p in control["soil_phase"][w0:n]], axis=0)
            a_ph = np.mean([p[name] for p in arm["soil_phase"][w0:n]], axis=0)
            soil[name] = {
                "control_land_mean_amplitude_k": land_mean(c_amp, land, control["weights"]),
                "arm_land_mean_amplitude_k": land_mean(a_amp, land, control["weights"]),
                "land_mean_phase_bin_shift":
                    land_mean(a_ph - c_ph, land, control["weights"]),
            }

        payload["arm"] = {
            "run": arm["run"], "orbits": arm["orbits"],
            "land_mean_range_k_per_orbit": [float(x) for x in a_series],
            "window_land_mean_range_k": float(a_series[w0:n].mean()),
            "difference_k": stats,
            "resolution_bar_k": (float(bar) if bar else None),
            "resolved": bool(bar and abs(stats["mean"]) > bar),
            "fractional_change_of_land_mean_range":
                float(a_series[w0:n].mean() / c_series[w0:n].mean() - 1.0),
            "land_cells_whose_range_falls": falls,
            "land_cells": int(land.sum()),
            "per_cell_fractional_change": {
                "p05": float(np.nanpercentile(frac[land], 5)),
                "median": float(np.nanmedian(frac[land])),
                "p95": float(np.nanpercentile(frac[land], 95)),
                "min": float(np.nanmin(frac[land])),
                "max": float(np.nanmax(frac[land])),
            },
            "wetness_correlation": correlation,
            "soil_temperature_levels": soil,
        }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in payload.items() if k != "provenance"},
                     indent=2, default=str)[:4000])


if __name__ == "__main__":
    main()
