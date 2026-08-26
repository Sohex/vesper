#!/usr/bin/env python3
"""Is the DERIVED relaxation time a ceiling on the FITTED one?

Worldbuilding frame: every number here is a diagnostic of the Vesper climate
model. Nothing in this script is about the real world.

WHY THIS EXISTS. `assess_convergence.py`'s offset criterion, when the
exponential fit has nothing to grip, falls back to `slope * tau_expected` and
holds `tau_expected` FIXED at what the modelled mixed layer's heat capacity and
the run's own radiative damping imply. Holding it fixed is what makes the test
cheap, and it is conservative only while the derived time is a CEILING on the
true relaxation: too long a tau inflates the remaining offset, which can only
refuse a run. Three reasons were given for expecting a ceiling -- the damping is
an equilibrium response and a transient sees a stronger one, the heat capacity
is a full ocean mixed-layer column applied to a planet with land, and the
reservoirs that would lengthen it hold under a tenth of that column -- and all
three are arguments. Every convergence report already records
`relaxation_orbits_fitted` beside `relaxation_orbits_expected`, so the
comparison accumulates; this takes it.

THE CHECK, stated as a rule that can fail.

  1. An artifact is EVIDENCE only if all three hold:
     a. the exponential fit was USED rather than the drift fallback. The
        fallback is identifiable from the artifact alone: it sets
        `temperature_remaining_offset_k` exactly equal to
        `remaining_offset_implied_by_drift_k`. A tau from a fit the assessment
        itself discarded is not a measurement of anything.
     b. the assessment reports the fitted tau as IDENTIFIABLE. Over a span
        short compared with tau the exponential is indistinguishable from a
        straight line, so tau is whatever the optimiser drifted to: 44018
        orbits on one run here and -337429 on another. Reports written since
        `relaxation_fit_identifiable` existed carry the verdict; older ones are
        judged by the same rule applied here, which is that the fitted tau must
        be finite, positive, and no longer than the span it was fitted over.
  2. On an evidence artifact the ceiling HOLDS when tau_fitted <= tau_expected,
     and is FALSIFIED only when tau_fitted exceeds tau_expected by more than the
     fit's own standard error on tau. One fitted value above a derived one is
     otherwise indistinguishable from the fit's noise. Where the artifact does
     not record that error the comparison is point against point and the row is
     labelled UNBRACKETED: it can raise a doubt and it cannot settle one.
  3. The ceiling SURVIVES when no evidence artifact falsifies it. It is
     UNTESTED, not confirmed, when no artifact is evidence.

WHAT A FALSIFICATION WOULD COST is in `exoplasim/notes/convergence-lengths.md`:
the criterion would need a bracket over the relaxation time, and the record
length it needs grows with the top of that bracket.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from _paths import ANALYSIS

# `approach_to_equilibrium` fits the last 65 per cent of the series. Restated
# here rather than imported because this script reads ARTIFACTS, including ones
# written before the span was recorded in them, and has to reconstruct it the
# same way the assessment did.
FIT_TAIL_FRACTION = 0.35


def fit_span_orbits(report: dict) -> int:
    """How many orbits the exponential fit saw, recorded or reconstructed."""
    recorded = report.get("metrics", {}).get("relaxation_fit_span_orbits")
    if recorded is not None:
        return int(recorded)
    end = report.get("window_end_year_index")
    n = int(end) + 1 if end is not None else int(report.get("completed_orbits", 0))
    if n <= 0:
        return 0
    # orbits are 0 .. n-1 and the mask is `orbit >= (n-1) * 0.35`.
    cut = (n - 1) * FIT_TAIL_FRACTION
    return sum(1 for orbit in range(n) if orbit >= cut)


def assess(path: Path) -> dict:
    report = json.loads(path.read_text(encoding="utf-8"))
    metrics = report.get("metrics", {})
    fitted = metrics.get("relaxation_orbits_fitted")
    expected = metrics.get("relaxation_orbits_expected")
    offset = metrics.get("temperature_remaining_offset_k")
    drift_offset = metrics.get("remaining_offset_implied_by_drift_k")
    standard_error = metrics.get("relaxation_orbits_fitted_standard_error")
    identifiable = metrics.get("relaxation_fit_identifiable")
    span = fit_span_orbits(report)

    row = {
        "artifact": path.name,
        "run": report.get("run_dir", "").rsplit("/", 1)[-1],
        "completed_orbits": report.get("completed_orbits"),
        "window_orbits": report.get("window_orbits"),
        "relaxation_orbits_fitted": fitted,
        "relaxation_orbits_fitted_standard_error": standard_error,
        "relaxation_orbits_expected": expected,
        "fit_span_orbits": span,
        "relaxation_fit_identifiable": identifiable,
    }

    if fitted is None and metrics.get("relaxation_fit_raw_tau_orbits") is not None:
        # The assessment refused to report a tau. That refusal IS the answer for
        # this artifact, and reading the raw number past it would be reading the
        # field the refusal exists to withhold.
        row.update(evidence=False,
                   why=metrics.get("relaxation_fit_verdict",
                                   "the assessment withheld the fitted tau"),
                   verdict="not evidence")
        return row
    if fitted is None or expected is None:
        row.update(evidence=False, why="the artifact records no relaxation pair",
                   verdict="not evidence")
        return row
    fell_back = (offset is not None and drift_offset is not None
                 and offset == drift_offset)
    if fell_back:
        row.update(evidence=False, why="the assessment took the drift fallback, "
                                       "so this tau is from a fit it discarded",
                   verdict="not evidence")
        return row
    if identifiable is False:
        row.update(evidence=False,
                   why=metrics.get("relaxation_fit_verdict",
                                   "the assessment reports the fitted tau as "
                                   "unidentifiable from its own series"),
                   verdict="not evidence")
        return row
    if not math.isfinite(fitted) or fitted <= 0.0:
        row.update(evidence=False, why="the fitted tau is not a positive finite "
                                       "number; the exponential had nothing to grip",
                   verdict="not evidence")
        return row
    if fitted > span:
        row.update(evidence=False,
                   why=f"the fitted tau of {fitted:.3g} orbits is longer than the "
                       f"{span} orbits it was fitted over, so it extrapolates "
                       "past the record",
                   verdict="not evidence")
        return row

    bracketed = standard_error is not None and math.isfinite(standard_error)
    row["bracketed"] = bool(bracketed)
    excess = fitted - expected
    row["fitted_minus_expected_orbits"] = excess
    if excess <= 0.0:
        verdict = "ceiling holds"
        why = (f"the fitted {fitted:.3f} is below the derived {expected:.3f} "
               f"by {-excess:.3f} orbits")
    elif bracketed and excess <= standard_error:
        verdict = "ceiling holds"
        why = (f"the fitted {fitted:.3f} is above the derived {expected:.3f} by "
               f"{excess:.3f} orbits, within the fit's own {standard_error:.3f}")
    else:
        verdict = "ceiling falsified"
        why = (f"the fitted {fitted:.3f} exceeds the derived {expected:.3f} by "
               f"{excess:.3f} orbits"
               + (f", more than the fit's own {standard_error:.3f}" if bracketed
                  else ", and the artifact records no error on the fit"))
    row.update(evidence=True, verdict=verdict, why=why)
    # HOW MUCH THE RUN COULD SEE. A fitted tau whose own error is comparable to
    # the gap it is being compared against cannot settle the question either
    # way, and saying "holds" on such a row would be reading the instrument
    # below its own scatter.
    if bracketed:
        row["discriminates"] = bool(abs(excess) > standard_error)
    else:
        row["discriminates"] = None
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--convergence-dir", type=Path,
                        default=ANALYSIS / "convergence",
                        help="directory of *_convergence.json reports")
    parser.add_argument("--output", type=Path, default=None,
                        help="where to write the result; defaults to "
                             "relaxation_ceiling.json beside the reports")
    args = parser.parse_args()
    reports = sorted(p for p in args.convergence_dir.glob("*_convergence*.json"))
    rows = [assess(path) for path in reports]
    evidence = [row for row in rows if row.get("evidence")]
    falsified = [row for row in evidence if row["verdict"] == "ceiling falsified"]
    discriminating = [row for row in evidence if row.get("discriminates")]

    if not evidence:
        verdict = "untested"
    elif falsified:
        verdict = "falsified"
    else:
        verdict = "survives"

    result = {
        "generator": "exoplasim/scripts/check_relaxation_ceiling.py",
        "question": "is relaxation_orbits_expected a ceiling on "
                    "relaxation_orbits_fitted",
        "rule": {
            "evidence": "the fit was used rather than the drift fallback, and "
                        "the fitted tau is identifiable from its own series: "
                        "finite, positive, and no longer than the span fitted",
            "falsified": "an evidence artifact whose fitted tau exceeds the "
                         "derived one by more than the fit's own standard error, "
                         "or by any amount where no such error is recorded",
            "survives": "no evidence artifact falsifies it",
            "untested": "no artifact is evidence",
        },
        "artifacts_read": len(rows),
        "evidence_artifacts": len(evidence),
        "artifacts_that_discriminate": len(discriminating),
        "verdict": verdict,
        "rows": rows,
    }
    output = args.output or (args.convergence_dir / "relaxation_ceiling.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=2))
    for row in rows:
        print(f"  {row['artifact']:<52} {row['verdict']:<18} {row.get('why', '')}")
    print(f"\nwritten to {output}")


if __name__ == "__main__":
    main()
