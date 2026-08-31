"""Validate the equilibrium contract's drift bound against a known answer.

Vesper is a simulated world. This instrument validates the statistic that
decides whether a LPJ-GUESS run of its simulated biosphere has settled. It runs
no model: every series here is either a synthetic AR(1) whose memory time is
known by construction, or a series read from a completed run's own tables.

WHY A VALIDATION AND NOT A CALIBRATION. `lib/lpj_output.py:drift_bound` bounds
a field's end-to-end relative drift from above and the contract passes the
field when that bound is inside its tolerance. The quantity the contract
declares is the rate at which a field that truly drifts AT the tolerance is
nonetheless passed, and that rate is asserted analytically. An assertion of that
shape has a right answer, so it is testable: generate series whose drift is
known and count. Nothing here is fitted and nothing here is tuned; a failure
changes the statistic rather than a number in it.

THE FOUR ARMS.

  size      A true end-to-end relative drift of exactly the contract's own
            limit is injected and the share of trials the field PASSES is
            counted. This is the declared rate and the only barred arm.
  cost      Zero drift, and the share of trials REFUSED. It is not barred; it
            is what sizes the retained record and it is reported per memory
            time.
  known     At a memory time of one the construction must reproduce the
            textbook two-sample t, and the estimator must recover the AR(1)
            memory time it was given. A harness that cannot reproduce an answer
            it already knows is not a validation.
  model     size and cost again, at the per-field scatter and memory time
            measured on a real run's own record, so the reported cost is this
            model's rather than a synthetic's.

The row that asked for this proposed a SIGNIFICANCE form instead -- reject the
null of no drift, with the per-field level derived from a declared family rate
by a Sidak correction over the counted family size. It is measured here on the
same trials, for size and for power at one and two times the limit, so the two
forms are compared on one set of numbers rather than on two.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.signal import lfilter
import yaml

import _paths  # noqa: F401  (puts lib/ on the path)

from autocorrelation import integrated_time  # noqa: E402

from lpj_output import (ACCEPTANCE_PATH, POLICY_PATH, _manifest_for,  # noqa: E402
                        _read_rows, drift_bound, forcing_cycle_years,
                        read_policy, record_cycles_for_bound)

PROJECT_ROOT = _paths.PROJECT_ROOT

# The measured memory times of the fields this statistic judges span 11.1 to
# 212.7 complete forcing cycles on the longest record this project has. The
# sweep covers that range and exceeds it, and it starts at 1 because an
# independent series is the case whose answer is already known.
TAU_SWEEP = (1.0, 10.0, 20.0, 40.0, 80.0, 125.0, 175.0, 213.0, 300.0, 400.0)
TRIALS = 2000
POWER_MULTIPLES = (1.0, 2.0)


def ar1_paths(rng, tau: float, cycles: int, trials: int) -> np.ndarray:
    """`trials` AR(1) series of unit marginal variance and memory time `tau`.

    An AR(1) with lag-1 correlation r has an integrated autocorrelation time of
    (1 + r) / (1 - r), which inverts to r = (tau - 1) / (tau + 1). The series is
    started from its own stationary distribution rather than from zero, so no
    part of it is a burn-in the statistic would read as a trend.
    """
    r = (tau - 1.0) / (tau + 1.0)
    if r <= 0.0:
        return rng.standard_normal((trials, cycles))
    innovation = rng.standard_normal((trials, cycles)) * np.sqrt(1.0 - r * r)
    innovation[:, 0] = rng.standard_normal(trials)
    return lfilter([1.0], [1.0, -r], innovation, axis=1)


def sidak_level(family_rate: float, family_size: int) -> float:
    """The per-field two-sided level a declared FAMILY rate implies.

    The significance form's multiplicity correction, kept here because it is the
    candidate being compared against and not because the shipped form needs one.
    """
    return 1.0 - (1.0 - family_rate) ** (1.0 / family_size)


def significance_verdict(bound: dict, level: float, limit: float) -> bool:
    """The row's candidate: reject "no drift" at a Sidak-corrected level.

    Same doubled half-mean difference and same memory-corrected standard error,
    read off the same `drift_bound` call; what differs is the direction. A field
    is REFUSED when its drift is both over the limit and resolved, so a record
    that cannot resolve anything passes everything -- which is the property that
    makes it the wrong instrument for an acceptance gate, and the reason it is
    measured here rather than argued about.
    """
    if bound["relative_standard_error"] > 0:
        ratio = bound["relative_drift"] / bound["relative_standard_error"]
    else:
        ratio = np.inf
    critical = float(stats.t.isf(level / 2.0, bound["degrees_of_freedom"]))
    return bool(bound["relative_drift"] > limit and ratio > critical)


def arm(rng, tau: float, cycles: int, trials: int, level: float,
        alpha: float, limit: float, floor: float, scatter: float,
        drift_multiple: float) -> dict:
    """One (memory time, injected drift) cell, both forms measured together.

    `scatter` is the series' standard deviation as a fraction of its mean, which
    is what sets how hard the cell is: the drift is stated relative to the same
    mean, so the whole problem is the ratio of an end-to-end change to the
    scatter it has to be seen through.
    """
    paths = ar1_paths(rng, tau, cycles, trials)
    x = np.arange(cycles, dtype=float)
    x = (x - x.mean()) / (cycles - 1)
    series = 1.0 + scatter * paths + drift_multiple * limit * x[None, :]
    accepted, refused, significance = 0, 0, 0
    bounds, errors = [], []
    for row in series:
        bound = drift_bound(row, alpha, floor)
        bounds.append(bound["upper_bound"])
        errors.append(bound["relative_standard_error"])
        if bound["upper_bound"] <= limit:
            accepted += 1
        else:
            refused += 1
        if significance_verdict(bound, level, limit):
            significance += 1
    return {"tau": tau, "drift_multiple": drift_multiple,
            "trials": trials,
            "equivalence_accept_rate": accepted / trials,
            "equivalence_refuse_rate": refused / trials,
            "significance_refuse_rate": significance / trials,
            "median_upper_bound": float(np.median(bounds)),
            "median_relative_standard_error": float(np.median(errors))}


def known_answer(rng, cycles: int, trials: int, alpha: float,
                 limit: float, floor: float, scatter: float) -> dict:
    """The two answers this harness must reproduce before any of it is believed.

    At a memory time of one the doubled half-mean difference is an ordinary
    two-sample comparison of independent normals, so the share of trials whose
    bound excludes a drift injected exactly at the limit must be `alpha`. And
    the estimator must return the memory time it was handed, inside the regime
    `lib/autocorrelation.py` itself calls reliable.
    """
    paths = ar1_paths(rng, 1.0, cycles, trials)
    x = np.arange(cycles, dtype=float)
    x = (x - x.mean()) / (cycles - 1)
    series = 1.0 + scatter * paths + limit * x[None, :]
    accepted = sum(1 for row in series
                   if drift_bound(row, alpha, floor)["upper_bound"] <= limit)
    recovered = {}
    for tau in (1.0, 10.0, 40.0, 125.0, 213.0):
        sample = ar1_paths(rng, tau, cycles, 200)
        estimates = [integrated_time(row) for row in sample]
        recovered[tau] = {
            "median": float(np.median([e["tau"] for e in estimates])),
            "median_upper": float(np.median([e["upper"] for e in estimates])),
            # How often the upper end actually covers the memory time the series
            # was built with. The point estimate is biased low, so this is the
            # number that says whether the conservatism is enough.
            "upper_covers_truth": float(np.mean(
                [e["upper"] >= tau for e in estimates])),
            "reliable_bar_cycles": 10.0 * tau,
            "record_cycles": cycles}
    return {"independent_accept_rate": accepted / trials,
            "declared_alpha": alpha,
            "recovered_tau": recovered}


def model_scales(run_dir: Path, outputs: list[str], policy: dict,
                 cycle_years: int, alpha: float, limit: float) -> dict:
    """Per field, the scatter and memory time this model actually has.

    Read from the run's own retained record so the cost arm is quoted at this
    model's numbers. The scatter is taken about the record's own linear fit,
    because the arm that uses it injects a drift of its own and would otherwise
    charge the field twice for one.
    """
    scales = {}
    for output in outputs:
        names, cells, years, cube = _read_rows(run_dir / output)
        usable = (len(years) // cycle_years) * cycle_years
        record = cube[:usable].reshape(usable // cycle_years, cycle_years,
                                       len(cells), len(names)).mean(axis=1)
        span = record.shape[0]
        x = np.arange(span, dtype=float)
        floor = float(policy["trend"]["absolute_scale_floor"])
        for index, name in enumerate(names):
            spatial = record[:, :, index].mean(axis=1)
            level = float(abs(spatial.mean()))
            if level <= floor or not np.isfinite(spatial).all() or spatial.std() == 0:
                continue
            flat = spatial - np.polyval(np.polyfit(x, spatial, 1), x)
            scales[f"{output} {name}"] = {
                "record_cycles": int(span),
                "relative_scatter": float(np.std(flat, ddof=1) / level),
                "tau_cycles": float(integrated_time(
                    flat + spatial.mean())["tau"]),
                # Through the same helper the reducer uses, so the length this
                # instrument quotes and the length a refusal quotes are one
                # number rather than two that agree by inspection.
                "cycles_for_bound": float(record_cycles_for_bound(
                    spatial, span, alpha, limit, floor)),
            }
    return scales


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path,
                        help="a completed run directory, for the model arm")
    parser.add_argument("--policy", type=Path, default=POLICY_PATH)
    parser.add_argument("--cycles", type=int, default=1253,
                        help="synthetic record length; the default is the "
                             "longest record this project has")
    parser.add_argument("--trials", type=int, default=TRIALS)
    parser.add_argument("--scatter", type=float, default=0.02,
                        help="synthetic relative scatter of the series")
    parser.add_argument("--family-rate", type=float, default=0.05,
                        help="the declared FAMILY rate the significance form "
                             "corrects to, for comparison only")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    policy = read_policy(args.policy)
    limit = float(policy["trend"]["relative_end_to_end_limit"])
    floor = float(policy["trend"]["absolute_scale_floor"])
    alpha = float(policy["trend"].get("maximum_false_acceptance_rate", 0.05))
    outputs = yaml.safe_load(ACCEPTANCE_PATH.read_text())["stability_outputs"]
    cycle_years, _ = forcing_cycle_years(
        args.run / outputs[0], _manifest_for(args.run / outputs[0])[1])
    scales = model_scales(args.run, outputs, policy, cycle_years, alpha, limit)
    family_size = len(scales)
    level = sidak_level(args.family_rate, family_size)

    rng = np.random.default_rng(20260831)
    known = known_answer(rng, args.cycles, args.trials, alpha, limit, floor,
                         args.scatter)
    rows = []
    for tau in TAU_SWEEP:
        for multiple in (0.0,) + POWER_MULTIPLES:
            rows.append(arm(rng, tau, args.cycles, args.trials, level, alpha,
                            limit, floor, args.scatter, multiple))

    model_rows = [{"field": name, **scale} for name, scale in sorted(scales.items())]

    report = {
        "contract_version": policy["contract_version"],
        "declared_alpha": alpha,
        "relative_end_to_end_limit": limit,
        "synthetic": {"cycles": args.cycles, "trials": args.trials,
                      "relative_scatter": args.scatter,
                      "family_size": family_size,
                      "significance_family_rate": args.family_rate,
                      "significance_per_field_level": level},
        "known_answer": known,
        "sweep": rows,
        "model_fields": model_rows,
    }

    print(f"KNOWN ANSWER, {args.trials} trials at memory time 1, drift injected "
          f"at the contract's own limit {limit:g}")
    print(f"  share accepted {known['independent_accept_rate']:.4f} against a "
          f"declared {alpha:g}; Monte-Carlo standard error "
          f"{np.sqrt(alpha * (1 - alpha) / args.trials):.4f}")
    for tau, row in known["recovered_tau"].items():
        print(f"  AR(1) tau {tau:>6.1f} recovered {row['median']:>8.2f}, upper "
              f"{row['median_upper']:>8.2f}, upper covers the truth "
              f"{row['upper_covers_truth']:.3f} of the time over "
              f"{row['record_cycles']} cycles")
    print(f"\nSWEEP, {args.cycles} cycles, relative scatter {args.scatter:g}, "
          f"significance form at a family rate of {args.family_rate:g} over "
          f"{family_size} fields (per-field level {level:.2e})")
    print(f"\n{'tau':>7}{'drift':>7}{'accept':>9}{'refuse':>9}"
          f"{'sig refuse':>12}{'median U':>10}{'median seR':>12}")
    for row in rows:
        print(f"{row['tau']:>7.0f}{row['drift_multiple']:>7.1f}"
              f"{row['equivalence_accept_rate']:>9.4f}"
              f"{row['equivalence_refuse_rate']:>9.4f}"
              f"{row['significance_refuse_rate']:>12.4f}"
              f"{row['median_upper_bound']:>10.4f}"
              f"{row['median_relative_standard_error']:>12.4f}")
    print(f"\nTHIS MODEL'S FIELDS, from {args.run.name}")
    print(f"\n{'field':<26}{'record':>8}{'scatter':>10}{'tau':>8}{'needs':>10}")
    for row in sorted(model_rows, key=lambda item: -item["cycles_for_bound"]):
        needed = ("no finite record" if not np.isfinite(row["cycles_for_bound"])
                  else f"{row['cycles_for_bound']:.0f}")
        print(f"{row['field']:<26}{row['record_cycles']:>8}"
              f"{row['relative_scatter']:>10.5f}{row['tau_cycles']:>8.1f}"
              f"{needed:>10}")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n")
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
