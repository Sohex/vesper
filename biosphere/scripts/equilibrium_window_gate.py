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

REPORT = PROJECT_ROOT / "biosphere/generated/equilibrium_window_report.json"


def write_run(root: Path, name: str, *, seed: int = 1, npatch: int = 5,
              cycle_years: int = 2, cycles: int = 10, trend: float = 0.0,
              offset: float = 0.0, omit: tuple[int, int] | None = None) -> Path:
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
    lines = ["Lon Lat Year A B\n"]
    for cell, (lon, lat) in enumerate(((-10.0, 20.0), (30.0, -40.0))):
        for year in range(cycles * cycle_years):
            if omit == (cell, year):
                continue
            phase = -2.0 if year % cycle_years == 0 else 2.0
            cycle = year // cycle_years
            a = 10.0 + cell + phase + trend * cycle + offset
            b = 100.0 + 2.0 * cell + 0.5 * phase + 3.0 * trend * cycle + offset
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

        trending = write_run(root, "trending", trend=1.0)
        refuses("trending end window refused", lambda: reduce_table(trending),
                "still trending")
        short = write_run(root, "short", cycle_years=1, cycles=9)
        refuses("one-year forcing still needs a multi-year window",
                lambda: reduce_table(short), "needs 10 consecutive end years")
        incomplete = write_run(root, "incomplete", omit=(1, 19))
        refuses("rank/cell truncation refused", lambda: reduce_table(incomplete),
                "ended at different years")

        seed_peer = write_run(root, "seed-peer", seed=2, offset=0.5)
        patch_peer = write_run(root, "patch-peer", seed=1, npatch=10, offset=0.2)
        ensemble = reduce_table(baseline, [seed_peer, patch_peer])
        check("seed uncertainty measured",
              ensemble.report["uncertainty"]["seed"]["status"] == "measured",
              "two roots at npatch=5 produce a sample standard deviation")
        check("patch-count uncertainty measured",
              ensemble.report["uncertainty"]["patch_count"]["status"] == "measured",
              "npatch 5 and 10 produce a separately labelled range")

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
        "contract_version": "vesper-lpj-equilibrium-window/1",
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
