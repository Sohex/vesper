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

for every metric, where SEM is the standard error of that metric's mean over
the window -- the scatter of the annual means divided by the root of their
count. Two sigma on a difference of two means. A metric that misses is
reported with how far it missed, and the verdict is all-or-nothing.

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
from segments import production_window

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


def annual_series(run_dir: Path, window: int) -> dict:
    """Per-orbit area-weighted global means over the production window."""
    files = sorted(run_dir.glob("MOST.[0-9]*.nc"))
    if len(files) < window:
        raise SystemExit(f"{run_dir} holds {len(files)} orbits, fewer than the "
                         f"{window}-orbit window")
    start, end = production_window(run_dir, len(files), window)
    out: dict = {k: [] for k in METRICS}
    fields = {}
    for year in range(start, end + 1):
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
    return {"window": (start, end), "series": out, "fields": fields}


def compare(a: dict, b: dict) -> dict:
    rows = []
    for key in METRICS:
        sa, sb = a["series"].get(key), b["series"].get(key)
        if not sa or not sb:
            rows.append({"metric": key, "verdict": "absent from one run"})
            continue
        ma, mb = float(np.mean(sa)), float(np.mean(sb))
        # Standard error of each window's mean, from its own year-to-year
        # scatter. ddof=1: these are a sample, not the population.
        sem = max(float(np.std(sa, ddof=1) / np.sqrt(len(sa))),
                  float(np.std(sb, ddof=1) / np.sqrt(len(sb))))
        bound = SIGMA * np.sqrt(2.0) * sem
        rows.append({"metric": key, "a": ma, "b": mb, "difference": mb - ma,
                     "sem": sem, "bound": bound,
                     "sigmas": abs(mb - ma) / (np.sqrt(2.0) * sem) if sem else None,
                     "within": bool(abs(mb - ma) <= bound)})
    return {"metrics": rows,
            "same_equilibrium": all(r.get("within", False) for r in rows)}


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
        "schema_version": 1,
        "criterion": (f"|mean_A - mean_B| <= {SIGMA} * sqrt(2) * max(SEM_A, SEM_B), "
                      "declared before either arm was compared"),
        "window_orbits": args.window,
        "run_a": {"path": str(args.run_a), "window": a["window"]},
        "run_b": {"path": str(args.run_b), "window": b["window"]},
        "pattern": pattern(a, b),
        **verdict,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    name = f"{args.run_a.name}__vs__{args.run_b.name}.json"
    (args.output / name).write_text(json.dumps(report, indent=2) + "\n")

    print(f"{'metric':24s} {'A':>12s} {'B':>12s} {'B-A':>11s} {'bound':>11s} {'sigma':>7s}")
    for r in report["metrics"]:
        if "a" not in r:
            print(f"{r['metric']:24s}  {r['verdict']}")
            continue
        mark = "ok " if r["within"] else "OVER"
        print(f"{r['metric']:24s} {r['a']:12.5f} {r['b']:12.5f} {r['difference']:11.5f} "
              f"{r['bound']:11.5f} {r['sigmas']:6.2f}  {mark}")
    p = report["pattern"]
    if p:
        print(f"\npattern: surface temperature RMS {p['surface_temperature_rms_k']:.4f} K, "
              f"spatial correlation {p['spatial_correlation']:.6f}")
    print(f"\nsame equilibrium: {report['same_equilibrium']}")
    print(f"wrote {args.output / name}")
    return 0 if report["same_equilibrium"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
