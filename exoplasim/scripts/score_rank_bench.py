#!/usr/bin/env python3
"""Score the rank-layout benchmark against the rules declared before it ran.

    python exoplasim/scripts/score_rank_bench.py <run_dir> [<run_dir> ...]

Reads per-orbit wall time as the gap between consecutive `MOST_REST.NNNNN`,
which excludes preparation and is the same instrument
`notes/audits/nlowio-collective-deadlock.md` and `docs/src/pipeline/costs.md` use.

The thresholds are `exoplasim/notes/rank-layout-benchmark.md`'s and are NOT
arguments: under 5% is no difference, and arm 4 has to beat arm 3 by 10% on
throughput. They are here rather than on the command line so the answer cannot
be fitted by moving them.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import statistics as st
from pathlib import Path

NO_DIFFERENCE_PCT = 5.0
THROUGHPUT_ADOPT_PCT = 10.0


def orbit_gaps(run_dir: Path):
    """Seconds per orbit, from restart mtimes. First orbit has no predecessor."""
    stamps = []
    for f in run_dir.glob("MOST_REST.*"):
        m = re.fullmatch(r"MOST_REST\.(\d{5})", f.name)
        if m:
            stamps.append((int(m.group(1)), os.path.getmtime(f)))
    stamps.sort()
    return [(b[0], b[1] - a[1]) for a, b in zip(stamps, stamps[1:])]


def describe(run_dir: Path) -> dict:
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    gaps = [g for _, g in orbit_gaps(run_dir)]
    return {
        "run": run_dir.name,
        "ranks": manifest["physical"]["ranks"],
        "n": len(gaps),
        "gaps": [round(g, 1) for g in gaps],
        "median": st.median(gaps) if gaps else float("nan"),
        "spread_pct": (100 * (max(gaps) - min(gaps)) / st.median(gaps)) if len(gaps) > 1 else float("nan"),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dirs", type=Path, nargs="+")
    args = ap.parse_args()

    rows = [describe(d.resolve()) for d in args.run_dirs]
    print(f"{'run':22} {'ranks':>5} {'n':>2} {'median s':>9} {'spread':>7}  gaps")
    for r in rows:
        print(f"{r['run']:22} {r['ranks']:5d} {r['n']:2d} {r['median']:9.1f} "
              f"{r['spread_pct']:6.1f}%  {r['gaps']}")

    # Any arm whose own orbits scatter more than the recorded 3% is suspect:
    # something else was running, and the note says to distrust it rather than
    # average it in.
    for r in rows:
        if r["n"] > 1 and r["spread_pct"] > NO_DIFFERENCE_PCT:
            print(f"\nSUSPECT: {r['run']} scatters {r['spread_pct']:.1f}% within itself, "
                  f"above the {NO_DIFFERENCE_PCT}% floor. Do not score it; find the load.")

    print("\nPairwise, on the declared floor:")
    for i, a in enumerate(rows):
        for b in rows[i + 1:]:
            d = 100 * (b["median"] - a["median"]) / a["median"]
            verdict = ("no difference" if abs(d) < NO_DIFFERENCE_PCT
                       else f"{a['run']} faster" if d > 0 else f"{b['run']} faster")
            print(f"  {a['run']} vs {b['run']}: {d:+.1f}%  -> {verdict}")

    print(f"\nThroughput note: two concurrent 8-rank runs beat one 16-rank run when a "
          f"concurrent 8-rank orbit is under TWICE a 16-rank orbit, and are ADOPTED "
          f"only past {THROUGHPUT_ADOPT_PCT}%. Compute that from the arm-4 pair and arm 3.")


if __name__ == "__main__":
    main()
