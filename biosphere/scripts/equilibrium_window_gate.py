#!/usr/bin/env python3
"""No-simulation fixtures for BIO-12's shared LPJ equilibrium reducer."""

from __future__ import annotations

import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT.parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))

from lpj_output import (EquilibriumWindowError, cycles_for_bound, drift_bound,
                        read_policy, reduce_table, relaxation_time)
import run_lengths


def _relaxation_of(series: np.ndarray) -> dict:
    """`relaxation_time` at the memory time of the series it is handed.

    The estimator takes the memory time as an argument because the block means
    it differences carry it, and a fixture that passed a wrong one would be
    testing something other than the estimator. The memory time is taken about
    the series' own linear fit, exactly as `lpj_output.py:timescale_report`
    takes it: an approach left in inflates every lag correlation toward one and
    the tau that comes out describes the approach rather than the variability.
    """
    from autocorrelation import integrated_time
    x = np.arange(series.size, dtype=float)
    flat = series - np.polyval(np.polyfit(x, series, 1), x) + series.mean()
    return relaxation_time(series, integrated_time(flat)["tau"])

REPORT = PROJECT_ROOT / "biosphere/generated/equilibrium_window_report.json"


# The contract measures its own trending-cell null on the record BEFORE the
# acceptance window, so a fixture has to retain one: at a declared false-refusal
# rate of 0.05 that is 20 null windows, and 21 windows of two-year cycles is 420
# years. A fixture shorter than that is the insufficient-record case, not the
# stationary one.
FIXTURE_CYCLES = 210


def write_run(root: Path, name: str, *, seed: int = 1, npatch: int = 5,
              cycle_years: int = 2, cycles: int = FIXTURE_CYCLES,
              trend: float = 0.0, offset: float = 0.0, dipole: bool = False,
              jitter: float = 0.0, omit: tuple[int, int] | None = None) -> Path:
    run = root / name
    run.mkdir()
    manifest = {
        "run_id": name,
        "source_build": "fixture",
        "physical": {"nyear": cycles * cycle_years, "npatch": npatch,
                     "root_seed": seed, "nfix_a": 0.234, "nfix_b": -0.172},
        "stochastic_randomness": {"root_seed": seed},
        "forcing": {"format": "VESPDRV8", "cycle_years": cycle_years},
        "inputs": {"driver": {"sha256": "driver"},
                   "soilmap": {"sha256": "soil"},
                   "pfts": {"sha256": "pfts"}},
    }
    (run / "run_manifest.json").write_text(json.dumps(manifest) + "\n")
    # Cycle-to-cycle jitter, deterministic so a fixture is reproducible. Without
    # it a fixture's cycle means are exactly smooth, the record's integrated
    # autocorrelation time is a large fraction of the record, and the memory guard
    # refuses before the criterion under test is ever reached. Real output is
    # noisy at this scale; a fixture that is not cannot exercise the trend test.
    noise = np.random.default_rng(20260830).normal(0.0, jitter, (2, cycles)) \
        if jitter else np.zeros((2, cycles))
    lines = ["Lon Lat Year A B\n"]
    for cell, (lon, lat) in enumerate(((-10.0, 20.0), (30.0, -40.0))):
        for year in range(cycles * cycle_years):
            if omit == (cell, year):
                continue
            phase = -2.0 if year % cycle_years == 0 else 2.0
            cycle = year // cycle_years
            # A dipole gives the two cells equal means and opposite drift, so the
            # spatial mean the global half reads is flat to machine precision and
            # only the per-cell half can see it.
            sign = (1.0 if cell == 0 else -1.0) if dipole else 1.0
            base = 0.0 if dipole else float(cell)
            # The drift starts at the acceptance window, so the record the
            # contract measures its null on is stationary and the window it
            # judges is not. A drift spread over the whole record instead is
            # what the global half is for, and it is bounded there: an additive
            # slope over 210 cycles can never move the last ten by 5 per cent of
            # a mean it has already raised.
            drift = trend * max(0, cycle - (cycles - 10))
            wobble = noise[cell, cycle]
            a = 10.0 + base + phase + sign * drift + offset + wobble
            b = (100.0 + 2.0 * base + 0.5 * phase + sign * 3.0 * drift + offset
                 + 3.0 * wobble)
            lines.append(f"{lon} {lat} {year} {a:.8f} {b:.8f}\n")
    table = run / "fixture.out"
    table.write_text("".join(lines))
    return table


def main() -> None:
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    def refuses(name: str, function, phrase: str) -> None:
        try:
            function()
        except EquilibriumWindowError as exc:
            check(name, phrase in str(exc), str(exc))
        else:
            check(name, False, "fixture was accepted")

    with tempfile.TemporaryDirectory(prefix="vesper-equilibrium-") as directory:
        root = Path(directory)
        baseline = write_run(root, "baseline")
        reduced = reduce_table(baseline)
        check("complete-cycle mean", np.allclose(
            reduced.values[(-10.0, 20.0)], [10.0, 100.0]),
            f"mean {reduced.values[(-10.0, 20.0)].tolist()}")
        # The reported spread is over the SAME span as the reported value, which
        # is the whole retained record and not the per-cell half's window.
        expected_std_a = np.std(np.tile([-2.0, 2.0], FIXTURE_CYCLES), ddof=1)
        check("temporal spread is over the reported span", np.isclose(
            reduced.temporal_std[(-10.0, 20.0)][0], expected_std_a),
            f"sample std {reduced.temporal_std[(-10.0, 20.0)][0]:.8f} over "
            f"{reduced.report['reported']['annual_values']} annual values")
        check("the reported span is the record and the window is the cell half's",
              reduced.report["reported"]["complete_forcing_cycles"]
              == FIXTURE_CYCLES
              and reduced.report["window"]["complete_forcing_cycles"] == 10
              and reduced.report["window"]["annual_values"] == 20,
              f"reported {reduced.report['reported']} against cell-half window "
              f"{reduced.report['window']}")
        check("stationary periodic phase accepted",
              reduced.report["trend"]["verdict"] == "PASS",
              "a repeating two-year phase has constant cycle means")
        check("missing stochastic ensemble is explicit",
              reduced.report["uncertainty"]["seed"]["status"] == "not_measured"
              and reduced.report["uncertainty"]["patch_count"]["status"]
              == "not_measured",
              "one run is never reported as zero seed/patch uncertainty")

        check("trending-cell limit is measured, not declared",
              reduced.report["trend"]["cell_fraction_null"]["windows"] == 20
              and reduced.report["trend"]["cell_fraction_null"][
                  "false_refusal_rate"] <= reduced.report["trend"]["rule"][
                      "cell_fraction"]["per_field_false_refusal_rate"],
              str(reduced.report["trend"]["cell_fraction_null"]))
        check("every assessed field carries its own limit",
              all("trending_cell_fraction_limit" in item
                  for item in reduced.report["trend"]["fields"]),
              "one declared fraction cannot serve fields whose null spans "
              "three orders of magnitude")

        # A record whose spatial mean drifts is refused by the GLOBAL half, which
        # bounds the drift over the whole retained record. Both fixtures below
        # drift; the noiseless one shows the refusal is the drift and not the
        # scatter, and the noisy one shows the bound still catches it once the
        # series has memory to widen the error.
        memoryful = write_run(root, "memoryful", jitter=0.0, trend=1.0)
        refuses("a drifting record is refused by its own drift bound",
                lambda: reduce_table(memoryful), "does not bound its drift")
        trending = write_run(root, "trending", trend=1.0, jitter=2.0)
        refuses("a drifting record with memory is refused too",
                lambda: reduce_table(trending), "does not bound its drift")
        # The one case a spatial mean cannot see, and the reason the per-cell
        # half exists: the two cells drift in opposite directions, so the global
        # half has nothing to bound and hands the record on.
        dipole = write_run(root, "dipole", trend=1.0, dipole=True, jitter=0.5)
        refuses("cancelling regional drift refused, though the mean is flat",
                lambda: reduce_table(dipole), "still trending")
        starved = write_run(root, "starved", cycle_years=1, cycles=100)
        refuses("a record too short to measure the null is refused, not passed",
                lambda: reduce_table(starved), "retain at least 200 years")
        short = write_run(root, "short", cycle_years=1, cycles=9)
        refuses("one-year forcing still needs a multi-year window",
                lambda: reduce_table(short), "needs 10 consecutive end years")
        incomplete = write_run(root, "incomplete",
                               omit=(1, FIXTURE_CYCLES * 2 - 1))
        refuses("rank/cell truncation refused", lambda: reduce_table(incomplete),
                "ended at different years")

        seed_peer = write_run(root, "seed-peer", seed=2, offset=0.5)
        patch_peer = write_run(root, "patch-peer", seed=1, npatch=10,
                               offset=0.2)
        ensemble = reduce_table(baseline, [seed_peer, patch_peer])
        check("seed uncertainty measured",
              ensemble.report["uncertainty"]["seed"]["status"] == "measured",
              "two roots at npatch=5 produce a sample standard deviation")
        check("patch-count uncertainty measured",
              ensemble.report["uncertainty"]["patch_count"]["status"] == "measured",
              "npatch 5 and 10 produce a separately labelled range")

    # THE DRIFT BOUND'S OWN IDENTITIES. Each has a right answer that is known
    # before the code runs, which is what makes them tests rather than reports.
    policy = read_policy()
    limit = float(policy["trend"]["relative_end_to_end_limit"])
    alpha = float(policy["trend"]["maximum_false_acceptance_rate"])
    floor = float(policy["trend"]["absolute_scale_floor"])
    span = 1000
    ramp = np.linspace(-0.5, 0.5, span)
    # A half-to-half mean difference is HALF the end-to-end change of a steady
    # drift. Undoubled it silently doubles the tolerance it is judged against,
    # so the estimator has to return the change itself: a ramp of exactly the
    # limit must read as the limit and not as half of it.
    exact = drift_bound(1.0 + limit * ramp, alpha, floor)
    check("the drift estimate is the END-TO-END change, not the half difference",
          abs(exact["relative_drift"] - limit) < 0.002 * limit,
          f"a ramp of exactly {limit:g} reads "
          f"{exact['relative_drift']:.6f}; half of it would be {limit / 2:g}")
    # A flat series has no drift and the estimator must say so exactly, not
    # approximately: there is nothing for the doubling or the scale to act on.
    flat = drift_bound(np.full(span, 3.0) + np.tile([-1.0, 1.0], span // 2),
                       alpha, floor)
    check("an exactly alternating series carries no drift",
          flat["relative_drift"] < 1e-12,
          f"drift {flat['relative_drift']:.3e}")
    # The bound is one over the root of the record, so quadrupling the record
    # halves the standard error. That is an identity of the construction and it
    # is what makes the record floor converge.
    quiet = np.random.default_rng(20260831).normal(0.0, 0.01, 4 * span)
    short_error = drift_bound(1.0 + quiet[:span], alpha, floor)[
        "relative_standard_error"]
    long_error = drift_bound(1.0 + quiet, alpha, floor)[
        "relative_standard_error"]
    check("the standard error falls as one over the root of the record",
          0.4 < long_error / short_error < 0.6,
          f"{short_error:.5f} over {span} cycles against {long_error:.5f} over "
          f"{4 * span}, a ratio of {long_error / short_error:.3f} against 0.5")
    # And the inversion agrees with the test: at the length `cycles_for_bound`
    # returns, a settled field of that scatter bounds inside the limit. The
    # scatter here is deliberately wide enough that `span` cycles do NOT resolve
    # the limit, so the check exercises the inversion rather than its identity.
    wide = drift_bound(1.0 + 30.0 * quiet[:span], alpha, floor)
    needed = cycles_for_bound(wide["relative_standard_error"], span,
                              wide["tau_cycles"], alpha, limit)
    scaled = wide["relative_standard_error"] * np.sqrt(span / needed)
    check("the record floor is the length at which the bound closes",
          needed > span and 0.8 * limit
          <= (np.sqrt(2.0 / np.pi) + 1.65) * scaled <= 1.05 * limit,
          f"a standard error of {wide['relative_standard_error']:.4f} over "
          f"{span} cycles needs {needed:.0f}, where it falls to {scaled:.4f} "
          f"and the bound closes on {limit:g}")

    # THE RELAXATION ESTIMATOR AGAINST A TIMESCALE IT ALREADY KNOWS. `tau =
    # -Q / ln(ratio)` diverges as the contraction approaches one, so the two
    # answers here are a real exponential, whose e-folding time must come back,
    # and a straight line, which has no e-folding time and must be declined
    # rather than read out as a very long one.
    relax_noise = np.random.default_rng(20260901).normal(0.0, 0.002, 1200)
    steps = np.arange(1200, dtype=float)
    known_tau = 200.0
    approach = 5.0 - 2.0 * np.exp(-steps / known_tau) + relax_noise
    measured = _relaxation_of(approach)
    check("an exponential approach returns the e-folding time it was built from",
          measured["admissible"]
          and abs(measured["tau_cycles"] - known_tau) < 0.05 * known_tau,
          f"{measured.get('tau_cycles', float('nan')):.1f} cycles against "
          f"{known_tau:.0f}")
    straight = _relaxation_of(5.0 + 0.001 * steps + relax_noise)
    check("a straight line has no e-folding time and is declined",
          not straight["admissible"] and "carries no curvature" in
          straight.get("reason", ""),
          straight.get("reason", "admitted"))

    # THE ECOLOGICAL RUN LENGTHS, checked against something that can fail rather
    # than reported. The derivation is an identity: the spin-up it returns must
    # actually satisfy the inequality it was derived from, and a spin-up one per
    # cent shorter must not. That is a right answer, so it is a test.
    try:
        derived = run_lengths.ecological_run_cycles(PROJECT_ROOT)
    except RuntimeError as exc:
        check("ecological run lengths are derived from measured timescales",
              False, str(exc))
    else:
        import math as _math
        tolerance = run_lengths.ecological_drift_tolerance()
        record, spinup = derived["record_cycles"], derived["spinup_cycles"]

        # THE SPIN-UP IS ALWAYS DERIVED, and the check is that it satisfies the
        # inequality it came from AT THE RELAXATION TIME IT CLAIMS TO COVER --
        # the measured one where there is one, and EVERY one where there is not.
        # A number one per cent shorter must fail, or the derivation is loose
        # rather than a floor. A default fails in either branch.
        def residual_drift(span: float, tau: float) -> float:
            return _math.exp(-span / tau) * (1.0 - _math.exp(-record / tau))

        if derived["spinup_is_minimax"]:
            multiple, worst = run_lengths.ecological_spinup_multiple(tolerance)
            # The whole claim is "at every relaxation time", so the sweep is what
            # tests it: decades either side of the worst case, which no single
            # point could distinguish from a lucky choice.
            taus = [record / (worst * factor)
                    for factor in (0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0)]
            covered = max(residual_drift(spinup, tau) for tau in taus)
            at_worst = record / worst
            check("the minimax spin-up covers every relaxation time",
                  covered <= tolerance * 1.000001
                  and residual_drift(spinup * 0.99, at_worst) > tolerance,
                  f"{spinup:.0f} cycles leaves at most {covered:.4f} across "
                  f"{record:.0f} retained over relaxation times "
                  f"{min(taus):.0f} to {max(taus):.0f}, against a tolerance of "
                  f"{tolerance:g}; one per cent shorter leaves "
                  f"{residual_drift(spinup * 0.99, at_worst):.4f} at the worst "
                  f"case of {at_worst:.0f} cycles")
            check("the minimax spin-up is a multiple of the retained record",
                  abs(spinup - multiple * record) <= 1.0e-6 * spinup
                  and derived["spinup_basis"].endswith(
                      "ecological_spinup_multiple"),
                  f"{spinup:.1f} against {multiple:.5f} x {record:.0f}")
        else:
            tau = derived["brackets"]["relaxation_cycles_bracket"][1]
            check("the derived spin-up satisfies the drift it was derived from",
                  residual_drift(spinup, tau) <= tolerance * 1.000001
                  and residual_drift(spinup * 0.99, tau) > tolerance,
                  f"{spinup:.0f} cycles leaves "
                  f"{residual_drift(spinup, tau):.4f} "
                  f"across {record:.0f} retained, against a tolerance of "
                  f"{tolerance:g}; one per cent shorter leaves "
                  f"{residual_drift(spinup * 0.99, tau):.4f}")
        check("the spin-up names the derivation that produced it",
              bool(derived["spinup_basis"]) and bool(derived["spinup_reason"])
              and derived["total_cycles"] == spinup + record,
              f"{derived['spinup_basis']}")
        # The record floor is an inversion, so it has a right answer too: the
        # field that set it must resolve the tolerance AT that length and must
        # not resolve it one per cent short. `cycles_for_bound` is monotone in
        # the standard error, so re-asking it at the derived length with the
        # standard error that length implies has to return the length itself.
        check("the retained record resolves the drift the contract refuses",
              record >= derived["brackets"]["resolving_cycles_bracket"][1]
              and record > derived["brackets"]["resolving_cycles_bracket"][0],
              f"{record:.0f} cycles against a resolving bracket of "
              f"{derived['brackets']['resolving_cycles_bracket'][0]:.0f} to "
              f"{derived['brackets']['resolving_cycles_bracket'][1]:.0f}, and "
              f"{derived['brackets']['fields_with_no_finite_record']} fields "
              "for which no finite record resolves it")
        check("the ecological lengths are a floor and say why",
              derived["is_a_floor"] and bool(derived["floor_because"]),
              derived["floor_because"])
        # The two pairs are in different units and the module must not let one
        # be read as the other. A cycle is one orbit only while the driver's
        # cycle is one year, so the reading carries the cycle length with it.
        check("the ecological pair carries its own unit",
              derived["brackets"]["forcing_cycle_years"] >= 1
              and all(name.endswith("_cycles") for name in
                      ("spinup_cycles", "record_cycles", "total_cycles")),
              f"forcing cycle {derived['brackets']['forcing_cycle_years']} "
              "simulation years, and every ecological length is named in cycles "
              "while every climate length is named in orbits")

    consumers = {
        "feedback": PROJECT_ROOT / "exoplasim/scripts/build_surface_albedo.py",
        "scoring": PROJECT_ROOT / "biosphere/scripts/score_prediction.py",
        "soil carbon": PROJECT_ROOT / "pedology/scripts/build_soil.py",
    }
    for name, path in consumers.items():
        source = path.read_text()
        check(f"{name} uses shared reducer", "reduce_table(" in source,
              str(path.relative_to(PROJECT_ROOT)))
        check(f"{name} requires BIO-14 acceptance",
              "require_lpj_acceptance(" in source,
              str(path.relative_to(PROJECT_ROOT)))
        check(f"{name} has no latest-year selector",
              "latest" not in source[source.find("def read_"):
                                     source.find("def main")],
              "consumer reader contains no independent maximum-year rule")

    failed = [item for item in checks if not item["passed"]]
    report = {
        "contract_version": policy["contract_version"],
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "checks": checks,
        "summary": {"passed": len(checks) - len(failed), "failed": len(failed)},
        "verdict": "PASS" if not failed else "FAIL",
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    for item in checks:
        print(f"{'PASS' if item['passed'] else 'FAIL'}  {item['name']}: "
              f"{item['detail']}")
    print(f"\n{report['verdict']}: {len(checks) - len(failed)} passed, "
          f"{len(failed)} failed; {REPORT.relative_to(PROJECT_ROOT)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
