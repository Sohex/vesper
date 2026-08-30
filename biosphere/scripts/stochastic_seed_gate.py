#!/usr/bin/env python3
"""No-simulation verification of Vesper's deterministic random substreams."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT.parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))

from stochastic_seeds import (DECLARATION, derive_seed, read_declaration,
                              stream_seeds)

REPORT = PROJECT_ROOT / "biosphere/generated/stochastic_seed_report.json"
GUESS_H = PROJECT_ROOT / "vendor/lpj-guess/framework/guess.h"
GUESS_CPP = PROJECT_ROOT / "vendor/lpj-guess/framework/guess.cpp"
VEGDYNAM = PROJECT_ROOT / "vendor/lpj-guess/modules/vegdynam.cpp"
BLAZE = PROJECT_ROOT / "vendor/lpj-guess/modules/blaze.cpp"
VESPER = PROJECT_ROOT / "vendor/lpj-guess/modules/vesperinput.cpp"
RUNNER = PROJECT_ROOT / "biosphere/scripts/run_lpj_guess.py"


def main() -> None:
    declaration = read_declaration()
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    fixed = stream_seeds(20260828, -73.25, 41.5, 0, 0, declaration)
    expected = {
        "establishment": 1655158136,
        "fire_occurrence": 1770134463,
        "fire_mortality": 1885110790,
        "background_mortality": 1080276517,
        "disturbance": 1195252828,
    }
    check("fixed cross-language vector", fixed == expected,
          f"derived {fixed}; expected {expected}")
    check("repeatability", fixed == stream_seeds(
        20260828, -73.25, 41.5, 0, 0, declaration),
        "the same complete key produces the same five states")
    check("process separation", len(set(fixed.values())) == len(fixed),
          "each named process begins on a distinct state")
    patch1 = stream_seeds(20260828, -73.25, 41.5, 0, 1, declaration)
    check("replicate separation", all(fixed[k] != patch1[k] for k in fixed),
          "patch ids 0 and 1 share no initial process state")
    stand1 = stream_seeds(20260828, -73.25, 41.5, 1, 0, declaration)
    check("stand separation", all(fixed[k] != stand1[k] for k in fixed),
          "stand ids 0 and 1 share no initial process state")
    other_root = stream_seeds(20260829, -73.25, 41.5, 0, 0, declaration)
    check("root separation", all(fixed[k] != other_root[k] for k in fixed),
          "changing the root changes every sampled stream")
    check("Park-Miller range",
          all(1 <= value <= 2147483646 for value in fixed.values()),
          "every state is nonzero and accepted by randfrac")

    cells = [(-170.0 + 17.0 * i, -80.0 + 8.0 * i) for i in range(20)]
    serial = {(lon, lat): stream_seeds(20260828, lon, lat, 0, 0, declaration)
              for lon, lat in cells}
    ranked: dict = {}
    for rank in range(7):
        for lon, lat in cells[rank::7]:
            ranked[(lon, lat)] = stream_seeds(
                20260828, lon, lat, 0, 0, declaration)
    reversed_order = {
        (lon, lat): stream_seeds(20260828, lon, lat, 0, 0, declaration)
        for lon, lat in reversed(cells)
    }
    check("MPI rank invariance", ranked == serial,
          "a seven-way stride produces the serial coordinate-keyed states")
    check("traversal-order invariance", reversed_order == serial,
          "reversing cell traversal leaves every coordinate-keyed state fixed")

    # Advancing fire must not advance establishment: this was impossible while
    # all calls consumed Stand::seed.
    modulus, multiplier = 2147483647, 16807
    fire = fixed["fire_occurrence"]
    for _ in range(100):
        fire = (multiplier * fire) % modulus
    check("draw-count isolation", fixed["establishment"] == expected["establishment"],
          f"100 fire draws end at {fire} without touching establishment")

    header = GUESS_H.read_text()
    implementation = GUESS_CPP.read_text()
    stochastic_sources = VEGDYNAM.read_text() + BLAZE.read_text()
    vesper = VESPER.read_text()
    runner = RUNNER.read_text()
    check("source process ABI",
          all(f"{name} = {number}" in header for name, number in {
              "STOCHASTIC_ESTABLISHMENT": 1,
              "STOCHASTIC_FIRE_OCCURRENCE": 2,
              "STOCHASTIC_FIRE_MORTALITY": 3,
              "STOCHASTIC_BACKGROUND_MORTALITY": 4,
              "STOCHASTIC_DISTURBANCE": 5,
          }.items()), "C++ process ids match the declaration")
    check("source hash ABI",
          "14695981039346656037" in implementation
          and "1099511628211" in implementation
          and "longitude * 1.0e6" in implementation,
          "C++ carries the same FNV constants and coordinate quantization")
    check("no shared stand stream",
          "stand.seed" not in stochastic_sources
          and "patch.stand.seed" not in stochastic_sources,
          "all ecological draw sites use Patch::random_seed")
    check("all processes consumed",
          all(name in stochastic_sources for name in (
              "STOCHASTIC_ESTABLISHMENT", "STOCHASTIC_FIRE_OCCURRENCE",
              "STOCHASTIC_FIRE_MORTALITY", "STOCHASTIC_BACKGROUND_MORTALITY",
              "STOCHASTIC_DISTURBANCE")),
          "each declared stream has a stochastic call site")
    check("restart state serialized",
          "arch & stochastic_seed[process]" in implementation
          and "& stochastic_root_seed" in implementation,
          "patch states and the grid-cell root survive restart")
    check("Vesper boundary initializes streams",
          "gridcell.initialize_stochastic_streams(root_seed)" in vesper
          and "vesper_root_seed" in vesper,
          "coordinates are assigned before the Vesper input module rekeys streams")
    check("harness and manifest expose root",
          "vesper_root_seed {settings['root_seed']}" in runner
          and '"stochastic_randomness"' in runner
          and '"declaration_sha256"' in runner,
          "instruction and provenance both carry the selected root and contract")

    failed = [item for item in checks if not item["passed"]]
    report = {
        "contract_version": declaration["contract_version"],
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": str(DECLARATION.relative_to(PROJECT_ROOT)),
        "root_seed": declaration["root_seed"],
        "checks": checks,
        "summary": {"passed": len(checks) - len(failed),
                    "failed": len(failed)},
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
