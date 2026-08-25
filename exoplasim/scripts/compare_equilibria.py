#!/usr/bin/env python3
"""Are two runs at the same equilibrium, against their own internal variability?

Worldbuilding frame: both runs are integrations of the Vesper climate model.
Nothing here is about the real world.

    python exoplasim/scripts/compare_equilibria.py RUN_A RUN_B --window 10

WHY A TOLERANCE CANNOT BE PICKED. The model is chaotic, so two runs from
different initial conditions never agree pointwise and their climatologies
never agree exactly. "Approximately the same state" is meaningless until it is
said what a difference would have to exceed to mean anything, and the only
honest scale is the run's OWN year-to-year scatter: a ten-orbit mean has a
standard error, and two independent ten-orbit means of the SAME climate differ
by about sqrt(2) times it.

THE CRITERION, fixed here before either arm was compared:

    |mean_A - mean_B|  <  2 * sqrt(2) * max(SEM_A, SEM_B)

for every metric. Two sigma on a difference of two means. A metric that misses
is reported with how far it missed, and the verdict is all-or-nothing.

THE SEM IS NOT THE SCATTER OVER THE ROOT OF THE COUNT. That form assumes the
orbits in a window are independent samples and they are not: this model has
variability at the window's own timescale, so a mean over n consecutive orbits
carries far fewer than n independent samples. Taking the standard error from
the raw count understated it by about a factor of two on the one pair where it
has been measured, and a test that understates its own error rejects pairs
differing by less than the model's own drift -- with a sigma figure that looks
decisive. It called a pair 7.24 sigma and failed four of six metrics on a
whole-run difference of four tenths of a sigma. world-yj9o.

What replaces it is `lib/autocorrelation.py`: sigma * sqrt(tau / n), with tau
the integrated autocorrelation time. Two things follow and both are visible in
the report.

  THE AUTOCORRELATION IS MEASURED ON A LONGER SPAN THAN THE WINDOW. tau cannot
  be recovered from a span comparable to itself -- estimated from inside a
  twenty-orbit window it comes back at about half its true value, so a
  correction taken there recovers half the missing error and still looks
  rigorous. Each run's tau is taken over the LONGEST TAIL of its production
  series that is stationary, which is as much data as the run has to offer.

  AN INSTRUMENT THAT CANNOT BOUND ITS OWN ERROR REFUSES. Where no tail is long
  enough to support a tau, the metric is reported indeterminate and the verdict
  is not "the same equilibrium". That is the same stance
  `assess_convergence.py` already takes on a window it cannot close against the
  state, and it is the only honest one: two runs cannot be declared to agree
  within an error nobody has measured.

WHAT THIS IS FOR. Deciding whether a state converted up the resolution ladder
settles where a run started cold at that resolution settles. If it does, the
ladder can be climbed by conversion and the wall clock that buys is the whole
point; if it does not, conversion is an optimisation that changes the answer.
SPAT-8.

WHAT IT DOES NOT DO is tell you the two runs are the same run. They are not,
and at equilibrium they should not be: the criterion is about the CLIMATE, and
a pattern comparison is reported alongside precisely so a metric-by-metric pass
cannot be mistaken for a field-by-field one.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import _paths  # noqa: F401  -- puts lib and this directory on the path
from _paths import ANALYSIS
from segments import production_window, non_production_orbits
# CLAUDE.md names lib/autocorrelation.py as the one place a standard error is
# taken over a series with memory. The defect this script carried was a second,
# wrong copy of that arithmetic.
import autocorrelation as ac

# Two sigma on the difference of two independent means of the same climate.
SIGMA = 2.0

# Unit conversions applied AFTER the area weighting, so the bound scales with
# the metric and the printed table is legible. The model writes precipitation
# as a rate in m/s; the outgoing radiation terms are signed downward-positive,
# so an upward flux reads negative and is flipped to be read as a magnitude.
SCALE = {"precipitation_mm_day": 1000.0 * 86400.0,
         "toa_shortwave_up_w_m2": -1.0,
         "toa_longwave_up_w_m2": -1.0}

# What is compared, and why each earns its place. Global means unless the name
# says otherwise; every one is area weighted.
METRICS = {
    "surface_temperature_k": ("ts", "the quantity the design flux is set by"),
    "air_temperature_2m_k": ("tas", "what the biosphere is driven by"),
    "precipitation_mm_day": ("pr", "the hydrological cycle's amplitude"),
    # pyburn writes pr in m/s; a global mean of 3e-8 has no readable scale and
    # its standard error rounds to zero in the printed table. Converted here so
    # the number and its bound are both legible, and named for what it is.
    "sea_ice_fraction": ("sic", "the albedo feedback's state variable"),
    "toa_shortwave_up_w_m2": ("rsut", "reflected, where a cloud or ice change shows"),
    "toa_longwave_up_w_m2": ("rlut", "emitted, the other half of the balance"),
}


# THE LONGEST STATIONARY TAIL, and the floor below which there is no point
# asking. `stationary_enough` refuses a span whose fitted trend moves it by
# more than its own scatter, so the search walks the start forward until the
# guard passes -- the longest tail of the production series that is varying
# rather than still approaching. The rule is mechanical and is declared here,
# before it was run on anything: no span is chosen by looking at which one
# gives a convenient tau.
MIN_TAU_SPAN = 12


def stationary_tail(series) -> tuple:
    """(start index, verdict) for the longest tail that passes the guard.

    Returns a start of None when no tail of at least `MIN_TAU_SPAN` orbits is
    stationary, which is a statement about the run and not a failure here: a
    run still approaching equilibrium has no variability to measure yet, only
    an approach, and a tau taken across it describes the approach.
    """
    n = len(series)
    for start in range(0, n - MIN_TAU_SPAN + 1):
        verdict = ac.stationary_enough(series[start:])
        if verdict["stationary"]:
            return start, verdict
    return None, ac.stationary_enough(series[-MIN_TAU_SPAN:]) if n >= 3 else {
        "stationary": False, "reason": "too few orbits to fit a trend"}


def variability(series, window_orbits: int) -> dict:
    """What one `window_orbits`-long mean of this series is worth.

    sigma and tau both come from the longest stationary tail, which is the most
    data the run has; n is the window's, because that is what the mean is taken
    over. A tau the span cannot support is reported and NOT used -- the caller
    refuses rather than correcting by an unknown factor.
    """
    start, flat = stationary_tail(series)
    if start is None:
        return {"resolved": False,
                "reason": "no tail of this run is stationary: " + flat.get("reason", ""),
                "tau_span_orbits": 0}
    tail = np.asarray(series[start:], dtype=float)
    t = ac.integrated_time(tail)
    out = {"tau_span_start_orbit": int(start), "tau_span_orbits": int(tail.size),
           "tau_orbits": t["tau"], "lag1_autocorrelation": t["lag1"],
           "effective_samples_in_window": window_orbits / t["tau"],
           "orbit_scatter": float(np.std(tail, ddof=1)),
           "tail_drift_per_orbit": flat["drift_per_sample"]}
    if not t["reliable"]:
        out.update({"resolved": False,
                    "reason": (f"tau {t['tau']:.2f} was estimated over "
                               f"{tail.size} orbits, under the "
                               f"{ac.RELIABLE_SPAN_MULTIPLE:g}x span it needs "
                               "to mean anything")})
        return out
    out["resolved"] = True
    out["sem"] = float(out["orbit_scatter"] * np.sqrt(t["tau"] / window_orbits))
    return out


def annual_series(run_dir: Path, window: int) -> dict:
    """Per-orbit global means: the whole production series, and the window in it.

    THE WHOLE SERIES, not the window alone. The window is what the means are
    taken over; the rest is what the autocorrelation is measured on, and
    reading only the window is what left the standard error with nothing longer
    than itself to be calibrated against.
    """
    files = sorted(run_dir.glob("MOST.[0-9]*.nc"))
    if len(files) < window:
        raise SystemExit(f"{run_dir} holds {len(files)} orbits, fewer than the "
                         f"{window}-orbit window")
    start, end = production_window(run_dir, len(files), window)
    excluded = set(non_production_orbits(run_dir, range(len(files))))
    # Every production orbit up to the end of the window. A diagnostic orbit is
    # not the run's trajectory and is no more admissible in a variability
    # estimate than in the window itself.
    span = [y for y in range(end + 1) if y not in excluded]
    out: dict = {k: [] for k in METRICS}
    fields = {}
    for year in span:
        with Dataset(files[year]) as nc:
            w = np.cos(np.deg2rad(nc.variables["lat"][:]))[:, None]
            for key, (var, _) in METRICS.items():
                if var not in nc.variables:
                    out[key] = None
                    continue
                x = np.asarray(nc.variables[var][:], dtype=float)
                if x.ndim == 4:                    # a level axis: take the surface
                    x = x[:, -1]
                # Mean over the orbit's output records first, then area
                # weighted over the globe. The other order weights a record
                # by nothing and gives the same answer only for equal-length
                # records, which is not something to rely on.
                m = x.mean(axis=0)
                out[key].append(SCALE.get(key, 1.0)
                                * float((m * w).sum() / (w.sum() * m.shape[-1])))
            if year == end:
                v = np.asarray(nc.variables["ts"][:], dtype=float).mean(axis=0)
                fields["ts"] = v
                fields["weights"] = np.repeat(w, v.shape[-1], axis=1)
    return {"window": (start, end), "window_orbits": window,
            "production_orbits": span, "series": out, "fields": fields}


def compare(a: dict, b: dict) -> dict:
    rows = []
    window = a["window_orbits"]
    for key in METRICS:
        sa, sb = a["series"].get(key), b["series"].get(key)
        if not sa or not sb:
            rows.append({"metric": key, "verdict": "absent from one run"})
            continue
        # The window is the last `window` production orbits of each series;
        # `production_window` refuses a window with a diagnostic hole in it, so
        # these are contiguous.
        wa, wb = np.asarray(sa[-window:]), np.asarray(sb[-window:])
        ma, mb = float(wa.mean()), float(wb.mean())
        va, vb = variability(sa, window), variability(sb, window)
        row = {"metric": key, "a": ma, "b": mb, "difference": mb - ma,
               "variability_a": va, "variability_b": vb}
        if not (va["resolved"] and vb["resolved"]):
            row["within"] = None
            reasons = dict.fromkeys(v["reason"] for v in (va, vb)
                                    if not v["resolved"])
            row["verdict"] = "indeterminate: " + "; ".join(reasons)
            rows.append(row)
            continue
        sem = max(va["sem"], vb["sem"])
        bound = SIGMA * np.sqrt(2.0) * sem
        row.update({"sem": sem, "bound": bound,
                    "sigmas": abs(mb - ma) / (np.sqrt(2.0) * sem) if sem else None,
                    "within": bool(abs(mb - ma) <= bound),
                    # What this comparison could have seen. A difference below
                    # it is not a null result, it is an unasked question.
                    "smallest_resolvable_difference": bound})
        rows.append(row)
    # FAILS CLOSED on an indeterminate metric as well as on a miss. "Not shown
    # to differ" is not "shown to agree", and the two were the same verdict
    # while the error bar was taken from the wrong count.
    return {"metrics": rows,
            "same_equilibrium": all(r.get("within") is True for r in rows),
            "indeterminate_metrics": sorted(r["metric"] for r in rows
                                            if r.get("within") is None)}


def pattern(a: dict, b: dict) -> dict | None:
    """Reported, never thresholded: two chaotic runs differ in pattern."""
    fa, fb = a["fields"].get("ts"), b["fields"].get("ts")
    if fa is None or fb is None or fa.shape != fb.shape:
        return None
    w = a["fields"]["weights"]
    d = fa - fb
    rms = float(np.sqrt((d * d * w).sum() / w.sum()))
    aa, bb = fa - (fa * w).sum() / w.sum(), fb - (fb * w).sum() / w.sum()
    corr = float((aa * bb * w).sum()
                 / np.sqrt((aa * aa * w).sum() * (bb * bb * w).sum()))
    return {"surface_temperature_rms_k": rms, "spatial_correlation": corr,
            "note": "reported for context; two runs of a chaotic model differ "
                    "in pattern at equilibrium and no bound is applied here"}


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_a", type=Path)
    ap.add_argument("run_b", type=Path)
    ap.add_argument("--window", type=int, default=10,
                    help="orbits in the equilibrium window of each run")
    ap.add_argument("--output", type=Path, default=ANALYSIS / "ladder")
    args = ap.parse_args()

    a = annual_series(args.run_a.resolve(), args.window)
    b = annual_series(args.run_b.resolve(), args.window)
    verdict = compare(a, b)
    report = {
        "schema_version": 2,
        "criterion": (f"|mean_A - mean_B| <= {SIGMA} * sqrt(2) * max(SEM_A, SEM_B), "
                      "declared before either arm was compared"),
        "standard_error": (
            "sigma * sqrt(tau / n_window), with sigma and the integrated "
            "autocorrelation time tau taken over the longest stationary tail "
            "of each run's production series. Orbits within a window are not "
            "independent samples and the scatter over the root of their count "
            "understated this by about a factor of two. world-yj9o."),
        "window_orbits": args.window,
        "run_a": {"path": str(args.run_a), "window": a["window"],
                  "production_orbits_read": len(a["production_orbits"])},
        "run_b": {"path": str(args.run_b), "window": b["window"],
                  "production_orbits_read": len(b["production_orbits"])},
        "pattern": pattern(a, b),
        **verdict,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    name = f"{args.run_a.name}__vs__{args.run_b.name}.json"
    (args.output / name).write_text(json.dumps(report, indent=2) + "\n")

    print(f"{'metric':24s} {'A':>12s} {'B':>12s} {'B-A':>11s} {'bound':>11s} "
          f"{'sigma':>7s} {'tau':>6s} {'n_eff':>6s}")
    for r in report["metrics"]:
        if "a" not in r:
            print(f"{r['metric']:24s}  {r['verdict']}")
            continue
        if r["within"] is None:
            print(f"{r['metric']:24s} {r['a']:12.5f} {r['b']:12.5f} "
                  f"{r['difference']:11.5f}  {r['verdict']}")
            continue
        va, vb = r["variability_a"], r["variability_b"]
        tau = max(va["tau_orbits"], vb["tau_orbits"])
        mark = "ok " if r["within"] else "OVER"
        print(f"{r['metric']:24s} {r['a']:12.5f} {r['b']:12.5f} {r['difference']:11.5f} "
              f"{r['bound']:11.5f} {r['sigmas']:6.2f} {tau:6.2f} "
              f"{args.window / tau:6.2f}  {mark}")
    p = report["pattern"]
    if p:
        print(f"\npattern: surface temperature RMS {p['surface_temperature_rms_k']:.4f} K, "
              f"spatial correlation {p['spatial_correlation']:.6f}")
    if report["indeterminate_metrics"]:
        # NOT A PASS AND NOT A FAIL OF THE ARMS. The runs have not been shown
        # to disagree; they have not been shown to agree either, because
        # neither carries a stationary span long enough to say what a
        # window mean of it is worth.
        print("\nindeterminate: " + ", ".join(report["indeterminate_metrics"]))
        print("The comparison cannot bound its own error on those metrics, so "
              "it makes no claim about them. Longer production spans at "
              "equilibrium are what would settle it.")
    print(f"\nsame equilibrium: {report['same_equilibrium']}")
    print(f"wrote {args.output / name}")
    return 0 if report["same_equilibrium"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
