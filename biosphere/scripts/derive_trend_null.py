"""Measure the equilibrium-window trend test's false-positive rate on this model.

BIO-12's per-cell trend test flags a gridcell when its cycle-mean series drifts
more than `relative_end_to_end_limit` end to end with a slope over
`slope_standard_errors` standard errors, and the gate refuses a field when the
share of its occupied cells that are flagged exceeds a limit. That limit is only
meaningful above the rate the test flags cells when there is NO trend to find,
and this script measures that rate from the model itself.

`lib/lpj_output.py` now measures the same null inside the reducer, from the record
before the acceptance window of the run being judged. This script is the offline
instrument that established what that null looks like across fields, run lengths
and patch counts, and it is how a claim about the null is re-checked.

The forcing cycle is one simulation year and every year replays byte-identical
forcing, so a long run at fixed forcing carries no forced interannual signal:
all interannual variation in it is internal stochasticity, which is exactly the
null the limit has to clear. A per-cell linear fit over the whole run is removed
so any genuine secular approach to equilibrium cannot enter the null, the record
is cut into disjoint windows of the contract's own length, and the contract's own
test is applied to each. Because the gate refuses a run when ANY field exceeds
the limit, the null distribution of a refusal is the per-window MAXIMUM flagged
fraction over fields.

Power is measured on the same surrogate windows by adding a coherent linear trend
of a stated multiple of the contract's end-to-end limit to every cell.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml

import _paths  # noqa: F401  (puts lib/ on the path)

from lpj_output import (POLICY_PATH, _cell_fraction, _manifest_for,  # noqa: E402
                        _read_rows, forcing_cycle_years, read_policy, sha256)

PROJECT_ROOT = _paths.PROJECT_ROOT
ACCEPTANCE_PATH = PROJECT_ROOT / "biosphere/config/lpj_acceptance.yaml"


def cell_flag_fraction(windows: np.ndarray, policy: dict) -> np.ndarray:
    """Flagged-cell fraction per window per field, from the contract's own test."""
    return np.asarray([_cell_fraction(window, policy) for window in windows])


def surrogate(cube: np.ndarray, cycles: int,
              cycle_years: int) -> tuple[np.ndarray, np.ndarray]:
    """Detrend each cell's whole record, keep its mean, and cut disjoint windows.

    A window is `cycles` COMPLETE forcing cycles, so its length in retained years
    is that times the cycle's own length. Reading the record as though a cycle
    were a year would judge a multi-year forcing cycle on the wrong support.
    """
    x = np.arange(cube.shape[0], dtype=float)
    x -= x.mean()
    slope = np.einsum("t,tcf->cf", x, cube) / float(np.sum(x * x))
    mean = cube.mean(axis=0)
    flat = cube - x[:, None, None] * slope[None]
    window_years = cycles * cycle_years
    nwindow = cube.shape[0] // window_years
    used = flat[: nwindow * window_years]
    return used.reshape(nwindow, cycles, cycle_years,
                        cube.shape[1], cube.shape[2]).mean(axis=2), mean


def inject(windows: np.ndarray, mean: np.ndarray, policy: dict,
           multiple: float) -> np.ndarray:
    """Add a coherent end-to-end drift of `multiple` times the contract's limit."""
    cycles = windows.shape[1]
    x = np.arange(cycles, dtype=float)
    x -= x.mean()
    limit = float(policy["trend"]["relative_end_to_end_limit"])
    slope = multiple * limit * np.abs(mean) / (cycles - 1)
    return windows + x[None, :, None, None] * slope[None, None]


def assess(run_dir: Path, outputs: list[str], policy: dict, powers: list[float],
           cycle_years: int) -> dict:
    cycles = int(policy["complete_forcing_cycles"])
    per_field: dict[str, dict] = {}
    null_stack, power_stacks = [], {m: [] for m in powers}
    for output in outputs:
        path = run_dir / output
        names, cells, years, cube = _read_rows(path)
        if not np.isfinite(cube).all():
            raise ValueError(f"{path} has missing cell-year rows")
        windows, mean = surrogate(cube, cycles, cycle_years)
        null = cell_flag_fraction(windows, policy)
        null_stack.append(null)
        for multiple in powers:
            power_stacks[multiple].append(
                cell_flag_fraction(inject(windows, mean, policy, multiple), policy))
        occupancy = (np.abs(mean) > float(policy["trend"]["absolute_scale_floor"])
                     ).mean(axis=0)
        per_field[output] = {
            "cells": len(cells), "windows": int(null.shape[0]),
            "first_year": years[0], "last_year": years[-1],
            "fields": {name: {
                "occupied_cell_fraction": float(occupancy[i]),
                "null_median": float(np.median(null[:, i])),
                "null_p95": float(np.quantile(null[:, i], 0.95)),
                "null_max": float(null[:, i].max()),
                **{f"power_x{multiple:g}_median":
                   float(np.median(power_stacks[multiple][-1][:, i]))
                   for multiple in powers},
            } for i, name in enumerate(names)},
        }
    null_all = np.concatenate(null_stack, axis=1)
    per_window_max = null_all.max(axis=1)
    report = {
        "null": {
            "windows": int(per_window_max.size),
            "per_window_max_median": float(np.median(per_window_max)),
            "per_window_max_p90": float(np.quantile(per_window_max, 0.90)),
            "per_window_max_p95": float(np.quantile(per_window_max, 0.95)),
            "per_window_max_p99": float(np.quantile(per_window_max, 0.99)),
            "per_window_max_maximum": float(per_window_max.max()),
        },
        "power": {},
        "by_output": per_field,
    }
    for multiple in powers:
        stacked = np.concatenate(power_stacks[multiple], axis=1).max(axis=1)
        report["power"][f"x{multiple:g}"] = {
            "per_window_max_median": float(np.median(stacked)),
            "per_window_max_p05": float(np.quantile(stacked, 0.05)),
        }
    half = per_window_max.size // 2
    report["null"]["drift_contamination_check"] = {
        "first_half_median": float(np.median(per_window_max[:half])),
        "second_half_median": float(np.median(per_window_max[half:])),
    }
    return report


def render(report: dict, powers: list[float]) -> None:
    """Report to a reader. This instrument generates no artifact.

    The null it measures is measured inside the reducer per run, so there is no
    standing product for it to write; the finding and its numbers are recorded in
    `biosphere/notes/equilibrium-trend-null.md`.
    """
    identity, null = report["identity"], report["null"]
    print(f"{identity['run_id']}  npatch {identity['npatch']}  "
          f"nyear {identity['nyear']}  forcing cycle {identity['forcing_cycle_years']}y")
    print(f"policy {identity['policy']} ({identity['policy_sha256'][:12]})")
    print()
    print(f"per-window maximum over fields, {null['windows']} detrended windows:")
    print("  median {per_window_max_median:.4f}  p90 {per_window_max_p90:.4f}  "
          "p95 {per_window_max_p95:.4f}  p99 {per_window_max_p99:.4f}  "
          "max {per_window_max_maximum:.4f}".format(**null))
    drift = null["drift_contamination_check"]
    print(f"  residual-drift check: first-half median "
          f"{drift['first_half_median']:.4f}, second-half "
          f"{drift['second_half_median']:.4f}")
    for multiple in powers:
        row = report["power"][f"x{multiple:g}"]
        print(f"  power at {multiple:g}x the contract's drift limit: median "
              f"{row['per_window_max_median']:.4f}  p05 {row['per_window_max_p05']:.4f}")
    print()
    head = f"{'output':<16}{'field':<9}{'occupied':>9}{'null med':>10}{'null p95':>10}"
    head += "".join(f"{'power x' + f'{m:g}':>10}" for m in powers)
    print(head)
    for output, block in report["by_output"].items():
        for name, field in block["fields"].items():
            row = (f"{output:<16}{name:<9}{field['occupied_cell_fraction']:>9.2f}"
                   f"{field['null_median']:>10.4f}{field['null_p95']:>10.4f}")
            row += "".join(f"{field[f'power_x{m:g}_median']:>10.4f}" for m in powers)
            print(row)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--policy", type=Path, default=POLICY_PATH)
    parser.add_argument("--power-multiple", type=float, action="append")
    args = parser.parse_args()
    policy = read_policy(args.policy)
    outputs = yaml.safe_load(ACCEPTANCE_PATH.read_text())["stability_outputs"]
    powers = args.power_multiple or [1.0, 2.0]
    manifest = json.loads((args.run / "run_manifest.json").read_text())
    cycle_years, _ = forcing_cycle_years(
        args.run / outputs[0], _manifest_for(args.run / outputs[0])[1])
    report = assess(args.run, outputs, policy, powers, cycle_years)
    report["identity"] = {
        "run_id": manifest.get("run_id", args.run.name),
        "npatch": manifest.get("physical", {}).get("npatch"),
        "nyear": manifest.get("physical", {}).get("nyear"),
        "root_seed": manifest.get("stochastic_randomness", {}).get("root_seed"),
        "forcing_cycle_years": manifest.get("forcing", {}).get("cycle_years"),
        "source_build": manifest.get("source_build"),
        "policy": str(args.policy), "policy_sha256": sha256(args.policy),
        "stability_outputs": outputs,
    }
    render(report, powers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
