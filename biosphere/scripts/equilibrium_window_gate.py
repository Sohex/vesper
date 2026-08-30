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

from lpj_output import EquilibriumWindowError, reduce_table
import run_lengths

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
        expected_std_a = np.std(np.tile([-2.0, 2.0], 10), ddof=1)
        check("temporal spread", np.isclose(
            reduced.temporal_std[(-10.0, 20.0)][0], expected_std_a),
            f"sample std {reduced.temporal_std[(-10.0, 20.0)][0]:.8f}")
        check("fixed ten-cycle window",
              reduced.report["window"]["complete_forcing_cycles"] == 10
              and reduced.report["window"]["annual_values"] == 20,
              str(reduced.report["window"]))
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

        memoryful = write_run(root, "memoryful", jitter=0.0, trend=1.0)
        refuses("a window inside one memory time is refused, not judged",
                lambda: reduce_table(memoryful), "inside one memory time")

        # A record that drifts at its end is refused -- but by the memory guard,
        # not by the trend test, because a series with any memory at all cannot
        # support a slope fitted over ten cycles: the guard needs the window to be
        # ten times the memory time, so it admits only a memoryless series. That
        # the trend test is unreachable this way IS the contract's current state
        # and is asserted here rather than reached by contorting a fixture. The
        # trend test's own refusal path is exercised by the dipole below, whose
        # spatial mean is flat by construction so the guard has nothing to
        # establish and hands the record on.
        trending = write_run(root, "trending", trend=1.0, jitter=2.0)
        refuses("a drifting record is refused, and by the earlier guard",
                lambda: reduce_table(trending), "inside one memory time")
        dipole = write_run(root, "dipole", trend=1.0, dipole=True, jitter=2.0)
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
        tau = derived["brackets"]["relaxation_cycles_bracket"][1]
        tolerance = run_lengths.ecological_drift_tolerance()
        record, spinup = derived["record_cycles"], derived["spinup_cycles"]

        def residual_drift(span: float) -> float:
            return (_math.exp(-span / tau)
                    * (1.0 - _math.exp(-record / tau)))

        check("the derived spin-up satisfies the drift it was derived from",
              residual_drift(spinup) <= tolerance * 1.000001
              and residual_drift(spinup * 0.99) > tolerance,
              f"{spinup:.0f} cycles leaves {residual_drift(spinup):.4f} across "
              f"{record:.0f} retained, against a tolerance of {tolerance:g}; "
              f"one per cent shorter leaves {residual_drift(spinup*0.99):.4f}")
        check("the retained record can establish its own memory time",
              record >= run_lengths.RELIABLE_SPAN_MULTIPLE
              * derived["brackets"]["memory_cycles_bracket"][1],
              f"{record:.0f} cycles against a memory time of "
              f"{derived['brackets']['memory_cycles_bracket'][1]:.1f}")
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
        "contract_version": "vesper-lpj-equilibrium-window/2",
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
