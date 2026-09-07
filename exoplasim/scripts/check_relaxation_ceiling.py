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

  1. An artifact is EVIDENCE only if the assessment reports the fitted tau as
     IDENTIFIABLE. Over a span short compared with tau the exponential is
     indistinguishable from a straight line, so tau is whatever the optimiser
     drifted to: 44018 orbits on one run here and -337429 on another. Reports
     written since `relaxation_fit_identifiable` existed carry the verdict.

     THE QUESTION IS ABOUT TAU AND NOT ABOUT WHICH ESTIMATOR THE OFFSET
     CRITERION TOOK, and those were one test here until they came apart. The
     rule read "the exponential fit was USED rather than the drift fallback",
     on the argument that a tau from a fit the assessment discarded is not a
     measurement of anything. That argument held only while the assessment
     discarded a fit for being degenerate. It now also declines a fit whose
     ASYMPTOTE cannot resolve the offset criterion's 0.15 K threshold, which is
     a statement about the asymptote and says nothing at all about whether the
     series determined tau -- so reading the estimator choice here would throw
     away every fit in the tree and empty the bracket rule 4 emits. world-jejw.

     AN ARTIFACT THAT PREDATES THE VERDICT is judged by the fallback signature
     instead, which is what the old rule was reaching for: the fallback is
     identifiable from such an artifact alone, because it sets
     `temperature_remaining_offset_k` exactly equal to
     `remaining_offset_implied_by_drift_k`. Those artifacts are also held to the
     identifiability rule applied here, which is that the fitted tau must be
     finite, positive, and no longer than the span it was fitted over.
  2. On an evidence artifact the ceiling HOLDS when tau_fitted <= tau_expected,
     and is FALSIFIED only when tau_fitted exceeds tau_expected by more than the
     fit's own standard error on tau. One fitted value above a derived one is
     otherwise indistinguishable from the fit's noise. Where the artifact does
     not record that error the comparison is point against point and the row is
     labelled UNBRACKETED: it can raise a doubt and it cannot settle one.
  3. The ceiling SURVIVES when no evidence artifact falsifies it. It is
     UNTESTED, not confirmed, when no artifact is evidence.

  4. The EVIDENCE ROWS ARE ALSO THE BRACKET. The same rule that decides what
     may test the ceiling decides what may bound the relaxation time, because
     it is one question -- which fits are measurements of it -- and answering it
     twice is how the two answers drift apart. `fitted_bracket_orbits` is the
     range of the evidence rows' fitted taus, and `lib/run_lengths.py` reads it
     rather than declaring a bracket of its own.

     A fit whose own standard error is comparable to the value does NOT drop out
     of the bracket. It cannot settle whether the ceiling holds, which is what
     `discriminates` says, and it is still this project's best reading of how
     fast this model returns: dropping it would narrow the bracket on the
     strength of the reading being uncertain, which is backwards. The bracket is
     what a settling block is bought in, so a wide one costs orbits and a
     falsely narrow one costs a transient in a climatology.

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
# THE FIT'S TAIL FRACTION IS STATED ONCE, in the assessment that owns the fit,
# and the span is RECONSTRUCTED BY THE ASSESSMENT'S OWN FUNCTION rather than by
# a second implementation of it here. This script reads artifacts written before
# the span was recorded in them, so it has to reconstruct the span the way the
# assessment did; calling the assessment IS that, and a copy of the number was
# only ever an approximation of it.
from assess_convergence import FIT_TAIL_FRACTION, fit_span_orbits


def recorded_or_reconstructed_span(report: dict) -> int:
    """How many orbits the exponential fit saw, recorded or reconstructed."""
    recorded = report.get("metrics", {}).get("relaxation_fit_span_orbits")
    if recorded is not None:
        return int(recorded)
    end = report.get("window_end_year_index")
    n = int(end) + 1 if end is not None else int(report.get("completed_orbits", 0))
    return fit_span_orbits(n)


def assess(path: Path) -> dict:
    report = json.loads(path.read_text(encoding="utf-8"))
    metrics = report.get("metrics", {})
    fitted = metrics.get("relaxation_orbits_fitted")
    expected = metrics.get("relaxation_orbits_expected")
    offset = metrics.get("temperature_remaining_offset_k")
    drift_offset = metrics.get("remaining_offset_implied_by_drift_k")
    standard_error = metrics.get("relaxation_orbits_fitted_standard_error")
    identifiable = metrics.get("relaxation_fit_identifiable")
    span = recorded_or_reconstructed_span(report)

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
    # ONLY WHERE THE ARTIFACT CARRIES NO IDENTIFIABILITY VERDICT OF ITS OWN.
    # See rule 1: on a report that has one, which estimator the offset criterion
    # took is a fact about the asymptote and not about tau.
    fell_back = (identifiable is None
                 and offset is not None and drift_offset is not None
                 and offset == drift_offset)
    if fell_back:
        row.update(evidence=False,
                   why="this artifact predates the identifiability verdict and "
                       "its assessment took the drift fallback, so the tau in "
                       "it is from a fit that was discarded as degenerate",
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

    # THE BRACKET OVER THE RELAXATION TIME, emitted from the same evidence rows
    # rather than from a second reading of them. `lib/run_lengths.py` reads this
    # and declares nothing; where there is no evidence there is no bracket, and
    # a null here makes a settling length refuse rather than quietly take a
    # default nobody measured.
    fits = sorted(row["relaxation_orbits_fitted"] for row in evidence)
    bracket = [fits[0], fits[-1]] if fits else None

    result = {
        "generator": "exoplasim/scripts/check_relaxation_ceiling.py",
        "question": "is relaxation_orbits_expected a ceiling on "
                    "relaxation_orbits_fitted",
        "rule": {
            "evidence": "the fitted tau is identifiable from its own series: "
                        "finite, positive, and no longer than the span fitted. "
                        "Which estimator the offset criterion took is a fact "
                        "about the fitted ASYMPTOTE and is not read here; an "
                        "artifact predating the identifiability verdict is "
                        "judged by the fallback signature instead",
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
        "fitted_bracket_orbits": bracket,
        "fitted_bracket_rule": "the range of the fitted relaxation times over "
                               "the evidence artifacts, which are the rows this "
                               "file's own evidence rule admits. Read by "
                               "lib/run_lengths.py:tau_relaxation_orbits_bracket",
        "fitted_bracket_from": [row["run"] for row in evidence],
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
