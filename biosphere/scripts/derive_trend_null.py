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
from scipy import stats
from scipy.optimize import curve_fit
import yaml

import _paths  # noqa: F401  (puts lib/ on the path)

from autocorrelation import (integrated_time,  # noqa: E402
                             stationary_enough)

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


SIGMA_SWEEP = (1.0, 1.5, 2.0, 2.306, 2.5, 3.0, 4.0)


def _leave_one_out_limits(fractions: np.ndarray) -> np.ndarray:
    """Window i's limit: the per-field maximum over every OTHER null window.

    The contract builds a field's limit from the windows beside the one it judges,
    so the null rate of a refusal is 1/(N+1) by exchangeability. Measuring that
    rate needs the same construction applied to a null window, which is this.
    """
    order = np.argsort(fractions, axis=0)
    largest, second = order[-1], order[-2]
    top = np.take_along_axis(fractions, largest[None], axis=0)[0]
    runner = np.take_along_axis(fractions, second[None], axis=0)[0]
    index = np.arange(fractions.shape[0])[:, None]
    return np.where(index == largest[None], runner[None], top[None])


def significance_rate(windows: np.ndarray, policy: dict) -> np.ndarray:
    """The bare slope test's flag rate, without the relative-drift condition.

    This is the quantity the row was originally about: what share of a field's
    occupied cells clear `slope_standard_errors` when nothing is drifting. Its
    nominal value under independent samples is the two-sided probability of the
    t distribution at the window's own degrees of freedom.
    """
    floor = float(policy["trend"]["absolute_scale_floor"])
    sigma = float(policy["trend"]["slope_standard_errors"])
    rates = []
    for window in windows:
        ncycle = window.shape[0]
        x = np.arange(ncycle, dtype=float)
        x -= x.mean()
        den = float(np.sum(x * x))
        slope = np.einsum("t,tcf->cf", x, window) / den
        intercept = window.mean(axis=0)
        residual = window - (intercept[None] + x[:, None, None] * slope[None])
        se = np.sqrt(np.sum(residual * residual, axis=0) / (ncycle - 2) / den)
        with np.errstate(divide="ignore", invalid="ignore"):
            significance = np.where(se > 0, np.abs(slope) / se,
                                    np.where(np.abs(slope) > 0, np.inf, 0.0))
        occupied = np.abs(intercept) > floor
        counts = occupied.sum(axis=0)
        rates.append(np.where(counts > 0,
                              ((significance > sigma) & occupied).sum(axis=0)
                              / np.maximum(counts, 1), 0.0))
    return np.asarray(rates)


def sigma_sweep(run_dir: Path, outputs: list[str], policy: dict,
                cycle_years: int, multiple: float = 2.0) -> dict:
    """Size and power of the contract's cell half against `slope_standard_errors`.

    Under contract 2 the cell half's size is fixed by the leave-one-out
    construction rather than by this number, so the number is a shape parameter of
    the statistic and its only effect is on power. Both are measured here rather
    than assumed: size that moves with sigma would mean the exchangeability
    argument the contract rests on is wrong.
    """
    cycles = int(policy["complete_forcing_cycles"])
    null_by_sigma = {sigma: [] for sigma in SIGMA_SWEEP}
    drift_by_sigma = {sigma: [] for sigma in SIGMA_SWEEP}
    bare, names_by_output = {}, {}
    for output in outputs:
        names, cells, years, cube = _read_rows(run_dir / output)
        windows, mean = surrogate(cube, cycles, cycle_years)
        drifted = inject(windows, mean, policy, multiple)
        names_by_output[output] = names
        for sigma in SIGMA_SWEEP:
            tuned = {**policy, "trend": {**policy["trend"],
                                         "slope_standard_errors": sigma}}
            null_by_sigma[sigma].append(cell_flag_fraction(windows, tuned))
            drift_by_sigma[sigma].append(cell_flag_fraction(drifted, tuned))
        bare[output] = significance_rate(windows, policy).mean(axis=0)
    rows = {}
    reference = None
    for sigma in SIGMA_SWEEP:
        null = np.concatenate(null_by_sigma[sigma], axis=1)
        drift = np.concatenate(drift_by_sigma[sigma], axis=1)
        limits = _leave_one_out_limits(null)
        refused_null = (null > limits).any(axis=1)
        refused_drift = (drift > limits).any(axis=1)
        if sigma == float(policy["trend"]["slope_standard_errors"]):
            reference = refused_drift
        rows[sigma] = {"size": float(refused_null.mean()),
                       "power": float(refused_drift.mean()),
                       "refused_drift": refused_drift}
    for sigma, row in rows.items():
        flips = row.pop("refused_drift")
        if reference is None:
            row["paired_standard_error"] = None
            continue
        discordant = int((flips != reference).sum())
        row["paired_standard_error"] = float(
            np.sqrt(discordant) / len(flips)) if discordant else 0.0
    return {"windows": int(len(reference)), "drift_multiple": multiple,
            "declared_sigma": float(policy["trend"]["slope_standard_errors"]),
            "expected_size": 1.0 / len(reference),
            "by_sigma": rows,
            "bare_significance_rate": {
                output: dict(zip(names_by_output[output],
                                 (float(v) for v in values)))
                for output, values in bare.items()}}


def render_sigma_sweep(sweep: dict, nominal: float) -> None:
    print(f"SIGMA SWEEP, {sweep['windows']} leave-one-out null windows, drift at "
          f"{sweep['drift_multiple']:g}x the contract's limit")
    print(f"  size is expected to be {sweep['expected_size']:.4f} at every sigma "
          "if the contract's exchangeability argument holds")
    print(f"\n{'sigma':>8}{'size':>10}{'power':>10}{'paired SE':>12}")
    for sigma, row in sweep["by_sigma"].items():
        mark = "   <-- declared" if sigma == sweep["declared_sigma"] else ""
        se = "-" if row["paired_standard_error"] is None else \
            f"{row['paired_standard_error']:.4f}"
        print(f"{sigma:>8.3f}{row['size']:>10.4f}{row['power']:>10.4f}"
              f"{se:>12}{mark}")
    print(f"\nbare slope-test flag rate at sigma {sweep['declared_sigma']:g}, "
          f"nominal {nominal:.4f} under independent samples:")
    for output, fields in sweep["bare_significance_rate"].items():
        shown = {k: v for k, v in fields.items() if v > 0}
        if not shown:
            continue
        worst = max(shown.items(), key=lambda kv: kv[1])
        print(f"  {output:<16} median over fields {np.median(list(shown.values())):.4f}"
              f"   worst {worst[0]} {worst[1]:.4f}")


def timescales(run_dir: Path, outputs: list[str], policy: dict,
               cycle_years: int, cells_sampled: int = 400) -> dict:
    """The biosphere's memory and relaxation times, at fixed forcing.

    `lib/run_lengths.py` keeps these two apart and never interchanges them, and
    the same split decides two different things here: the MEMORY time bounds the
    window and the retained record, because a window shorter than the memory time
    is one excursion of a process rather than a sample of it; the RELAXATION time
    bounds the SPIN-UP, because a run still approaching equilibrium carries a
    drift its own null construction would detrend away from its reference.

    Both are read through `lib/autocorrelation.py`, which owns the estimator and
    the two verdicts that decide whether its answer means anything: `reliable`
    for a span that can support the tau it produced, and `stationary_enough` for
    a span flat enough for a tau to describe its variability rather than its
    approach.
    """
    cycles = int(policy["complete_forcing_cycles"])
    floor = float(policy["trend"]["absolute_scale_floor"])
    rng = np.random.default_rng(20260830)
    out = {}
    for output in outputs:
        names, cells, years, cube = _read_rows(run_dir / output)
        record = cube[: (len(years) // cycle_years) * cycle_years].reshape(
            -1, cycle_years, len(cells), len(names)).mean(axis=1)
        span = record.shape[0]
        x = np.arange(span, dtype=float)
        pick = rng.choice(len(cells), min(cells_sampled, len(cells)), replace=False)
        mean = record.mean(axis=0)
        fields = {}
        for f, name in enumerate(names):
            spatial = record[:, :, f].mean(axis=1)
            if not np.isfinite(spatial).all() or spatial.std() == 0:
                continue
            flat = spatial - np.polyval(np.polyfit(x, spatial, 1), x) + spatial.mean()
            memory = integrated_time(flat)
            steady = stationary_enough(spatial)
            taus = []
            for c in pick:
                if abs(mean[c, f]) <= floor:
                    continue
                series = record[:, c, f]
                if series.std() == 0:
                    continue
                detrended = (series - np.polyval(np.polyfit(x, series, 1), x)
                             + series.mean())
                taus.append(integrated_time(detrended)["tau"])
            fields[name] = {
                "spatial_mean_tau": memory["tau"],
                "spatial_mean_tau_reliable": memory["reliable"],
                "spatial_mean_stationary": steady["stationary"],
                "cell_tau_median": float(np.median(taus)) if taus else None,
                "cell_tau_p90": float(np.quantile(taus, 0.90)) if taus else None,
                "cell_tau_over_window": (float(np.mean(np.asarray(taus) > cycles))
                                         if taus else None),
                "relaxation": _relaxation(spatial, span),
            }
        out[output] = {"span_cycles": span, "window_cycles": cycles,
                       "cells_sampled": int(len(pick)), "fields": fields}
    return out


def _relaxation(series: np.ndarray, span: int) -> dict:
    """The e-folding time of an exponential approach to the record's asymptote.

    ADMISSIBLE, declared before it was applied: the fitted asymptote must lie
    inside the record's own range and the e-folding time must be shorter than the
    record. A fit whose timescale exceeds the span that produced it has measured
    the span and not the timescale, which is failure-modes class 34, and it is
    reported as a lower bound rather than as a value.
    """
    x = np.arange(span, dtype=float)
    lo, hi = float(series.min()), float(series.max())
    try:
        (a, b, rate), _ = curve_fit(
            lambda t, a, b, r: a + b * np.exp(-r * t), x, series,
            p0=[float(series[-1]), float(series[0] - series[-1]), 1.0 / max(span / 5, 1)],
            maxfev=20000)
    except (RuntimeError, ValueError):
        return {"admissible": False, "reason": "no exponential fit converged"}
    if rate <= 0:
        return {"admissible": False, "reason": "fitted approach diverges",
                "lower_bound_cycles": float(span)}
    tau = 1.0 / rate
    if not lo <= a <= hi:
        return {"admissible": False, "reason": "asymptote outside the record",
                "e_folding_cycles": float(tau)}
    if tau >= span:
        return {"admissible": False, "reason": "timescale exceeds the record",
                "lower_bound_cycles": float(span)}
    return {"admissible": True, "e_folding_cycles": float(tau),
            "asymptote": float(a), "remaining_at_record_end":
                float(abs(b) * np.exp(-rate * (span - 1)))}


def render_timescales(scales: dict) -> None:
    print(f"{'output':<16}{'field':<9}{'sm tau':>8}{'rel':>5}{'stat':>6}"
          f"{'cell tau':>10}{'p90':>8}{'>window':>9}   relaxation")
    for output, block in scales.items():
        for name, field in block["fields"].items():
            if field["cell_tau_median"] is None:
                continue
            relax = field["relaxation"]
            note = (f"{relax['e_folding_cycles']:.1f} cycles"
                    if relax["admissible"] else relax["reason"])
            print(f"{output:<16}{name:<9}{field['spatial_mean_tau']:>8.1f}"
                  f"{'y' if field['spatial_mean_tau_reliable'] else 'n':>5}"
                  f"{'y' if field['spatial_mean_stationary'] else 'n':>6}"
                  f"{field['cell_tau_median']:>10.1f}{field['cell_tau_p90']:>8.1f}"
                  f"{field['cell_tau_over_window']:>9.2f}   {note}")


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
    parser.add_argument("--timescales", action="store_true",
                        help="the memory and relaxation times this record can "
                             "support (world-3erv)")
    parser.add_argument("--sigma-sweep", action="store_true",
                        help="size and power of the cell half against "
                             "trend.slope_standard_errors (world-ioxr)")
    args = parser.parse_args()
    policy = read_policy(args.policy)
    outputs = yaml.safe_load(ACCEPTANCE_PATH.read_text())["stability_outputs"]
    powers = args.power_multiple or [1.0, 2.0]
    manifest = json.loads((args.run / "run_manifest.json").read_text())
    cycle_years, _ = forcing_cycle_years(
        args.run / outputs[0], _manifest_for(args.run / outputs[0])[1])
    if args.timescales:
        render_timescales(timescales(args.run, outputs, policy, cycle_years))
        return 0
    if args.sigma_sweep:
        cycles = int(policy["complete_forcing_cycles"])
        nominal = float(2.0 * stats.t.sf(
            policy["trend"]["slope_standard_errors"], cycles - 2))
        render_sigma_sweep(
            sigma_sweep(args.run, outputs, policy, cycle_years), nominal)
        return 0
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
