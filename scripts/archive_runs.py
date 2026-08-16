#!/usr/bin/env python3
"""Extract what a dead run is worth, then delete the rest of it.

    python scripts/archive_runs.py                 # dry run, shows everything
    python scripts/archive_runs.py --execute       # do it

A finished run on a superseded terrain is 2-8 GB of NetCDF describing a world we
no longer model. What it is actually worth is a few hundred kilobytes: the
manifest, the convergence assessment, the climate series, and the derived report.
Those are what anything downstream ever cites.

So "archive" here does not mean move the directory. Moving 50 GB into a folder
called `archive/` leaves 50 GB. It means extract the derived products into
`archive/runs/<directory>/`, which is tracked, and delete the raw output, which
is not.

**Live is defined by the active build**, not by recency. A run is live if its
`source_build` matches `config/planet.yaml`. Everything else is dead by
definition, and the definition is worth more than a hand-kept list: a build being
superseded is exactly what makes its runs uninteresting.

Nothing is decided by directory name or sort order. The run set comes from
`exoplasim/runs/INDEX.json`.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "exoplasim" / "runs"
INDEX = RUNS / "INDEX.json"
ARCHIVE = ROOT / "archive" / "runs"

# Small enough to keep, and the only things anything cites.
KEEP_NAMES = ["run_manifest.json", "baseline_convergence.json",
              "baseline_climate_report.json"]
KEEP_GLOBS = ["*_climate_series.json", "*.cfg", "*_namelist", "*.nl"]

# Everything not in the keep list goes. A denylist of output patterns missed
# `ice_output` and `ocean_output` -- 220 MB per run in extensionless files -- and
# would have missed the next one too. For a dead run the safe default is delete,
# because the terrain that produced it is superseded and the keep list is the
# complete set of things anything cites.


def load_index() -> list:
    if not INDEX.is_file():
        raise SystemExit("run `python exoplasim/scripts/index_runs.py` first")
    return json.loads(INDEX.read_text(encoding="utf-8"))["runs"]


def size_of(paths) -> int:
    return sum(p.stat().st_size for p in paths if p.is_file())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true",
                    help="actually extract and delete; default is a dry run")
    args = ap.parse_args()

    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    active = cfg.get("source_build")

    rows = load_index()
    live = [r for r in rows if r.get("source_build") == active]
    dead = [r for r in rows if r.get("source_build") != active]

    print(f"active build: {active}\n")
    print(f"LIVE, untouched ({len(live)} runs, "
          f"{sum(r['size_gb'] for r in live):.1f} GB)")
    for r in live:
        print(f"  {r['directory'][:64]:66}{r['size_gb']:>7.2f} GB")

    print(f"\nDEAD, extract then delete ({len(dead)} runs, "
          f"{sum(r['size_gb'] for r in dead):.1f} GB)")

    freed = 0
    for r in dead:
        run = RUNS / r["directory"]
        if not run.is_dir():
            continue
        keep = [run / n for n in KEEP_NAMES if (run / n).is_file()]
        for g in KEEP_GLOBS:
            keep += [p for p in run.glob(g) if p.is_file()]
        drop = {p for p in run.rglob("*") if p.is_file()} - set(keep)
        freed += size_of(drop)
        print(f"  {r['directory'][:60]:62}"
              f"keep {len(keep):>3} files  drop {len(drop):>4} files"
              f"  {size_of(drop)/1e9:>6.2f} GB")

        if args.execute:
            dest = ARCHIVE / r["directory"]
            dest.mkdir(parents=True, exist_ok=True)
            for p in keep:
                shutil.copy2(p, dest / p.name)
            (dest / "INDEX_ENTRY.json").write_text(
                json.dumps(r, indent=2) + "\n", encoding="utf-8")
            for p in drop:
                p.unlink()
            # Directories that held only dropped files.
            for d in sorted((q for q in run.rglob("*") if q.is_dir()),
                            key=lambda q: -len(q.parts)):
                if not any(d.iterdir()):
                    d.rmdir()
            if not any(run.iterdir()):
                run.rmdir()

    print(f"\n{'freed' if args.execute else 'would free'} {freed/1e9:.1f} GB")
    if not args.execute:
        print("dry run; pass --execute to do it")
        return

    (ARCHIVE / "README.md").write_text(
        "# Archived runs\n\n"
        "Derived products from runs on superseded terrains. The raw NetCDF is\n"
        "gone: it described a world this project no longer models, and it was\n"
        "50 GB. What remains is what anything ever cited -- the manifest, the\n"
        "convergence assessment, the climate series, and the namelists.\n\n"
        "`INDEX_ENTRY.json` in each directory is that run's row from\n"
        "`exoplasim/runs/INDEX.json` at the time it was archived, so a result\n"
        "computed from one of these stays readable and datable without the\n"
        "output being present.\n\n"
        f"Archived {datetime.now(timezone.utc).date().isoformat()}.\n",
        encoding="utf-8")
    print(f"wrote {ARCHIVE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
