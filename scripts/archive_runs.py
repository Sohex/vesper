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
`archive/runs/<directory>/`, which is tracked, and delete the run, which is not.

The whole run directory goes, kept files included, because they are in the
archive by then. An earlier version deleted only the raw output and left the
small files in place, which meant every dead run still appeared in INDEX.json as
though it were present, and the same file existed in two places with only one of
them tracked.

**Nothing is deleted until the archive has been checked.** `verify()` confirms
each extracted file is present at the same size and that the JSON among them
parses, and raises otherwise. A half-written manifest is silent until someone
reads it, by which time the original is long gone.

**Live is defined by the active build**, not by recency. A run is live if its
`source_build` matches `config/planet.yaml`. Everything else is dead by
definition, and the definition is worth more than a hand-kept list: a build being
superseded is exactly what makes its runs uninteresting.

That rule knows about one of the two ways a run dies. The other is a
RE-COMMISSIONING: the terrain stands, but something that is not the terrain has
invalidated the climatology, and by CLAUDE.md rule 7 everything below a
climatology is worthless rather than stale once the climatology is. Those runs
are on the ACTIVE build and the live rule protects them, correctly, because it
cannot see the invalidation -- nothing in a run's manifest records that the model
or the config moved underneath it. `--include-live` is that judgment made by the
caller: it says every run in the index is dead, whatever build it names, and
archives the lot. It is DECLARED rather than inferred for the same reason a
segment's purpose is.

Nothing is decided by directory name or sort order. The run set comes from
`exoplasim/runs/INDEX.json`, and from the part of it that still has a payload:
that file is a LEDGER of every run that has existed, so most of its rows are
runs this script has already dealt with and nothing about them is a work item.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "exoplasim" / "runs"
INDEX = RUNS / "INDEX.json"
ARCHIVE = ROOT / "archive" / "runs"

# Small enough to keep, and the only things anything cites.
KEEP_NAMES = ["run_manifest.json", "baseline_convergence.json",
              "baseline_climate_report.json", "climate_report.json"]
KEEP_GLOBS = ["*_climate_series.json", "*.cfg", "*_namelist", "*.nl"]

# Everything not in the keep list goes. A denylist of output patterns missed
# `ice_output` and `ocean_output` -- 220 MB per run in extensionless files -- and
# would have missed the next one too. For a dead run the safe default is delete,
# because the terrain that produced it is superseded and the keep list is the
# complete set of things anything cites.


def verify(dest: Path, kept: list[Path]) -> None:
    """Confirm the archive is complete and readable before anything is deleted.

    Extract-then-delete is only safe if the extract worked. This checks the
    copies exist at the same size, and that the JSON among them parses -- a
    truncated or half-written manifest is exactly the failure that would make
    the archive worthless, and it is silent until someone tries to read it
    months later with the original long gone.

    Raises rather than returning a flag: there is no sensible way to continue.
    """
    for src in kept:
        copy = dest / src.name
        if not copy.is_file():
            raise RuntimeError(f"archive incomplete: {copy} missing; nothing deleted")
        if copy.stat().st_size != src.stat().st_size:
            raise RuntimeError(
                f"archive corrupt: {copy} is {copy.stat().st_size} bytes against "
                f"{src.stat().st_size} at source; nothing deleted")
    for js in sorted(dest.glob("*.json")):
        try:
            json.loads(js.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"archive corrupt: {js} does not parse ({exc}); "
                               f"nothing deleted") from exc


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
    # `--include-live --keep <id>` is the triage case: a re-commissioning where
    # every run on the active build is spent EXCEPT the one still being worked.
    # Without it the caller either archives the live run too or moves
    # directories out of the way by hand, and moving a run directory to protect
    # it from a tool that deletes run directories is not a safe manoeuvre.
    ap.add_argument("--keep", action="append", default=[], metavar="RUN_ID",
                    help="spare this run whatever the build rule says. Repeat "
                         "for several. An id that matches nothing is an error, "
                         "not a no-op: a typo would otherwise delete the run it "
                         "was meant to protect")
    ap.add_argument("--include-live", action="store_true",
                    help="treat every run as dead, including runs on the active "
                         "build. For a re-commissioning, where the terrain "
                         "stands but the climatology below it does not")
    args = ap.parse_args()

    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    active = cfg.get("source_build")

    # THE LEDGER OUTLIVES THE PAYLOAD, so most of it is not a work list. Since
    # world-ww6z, `exoplasim/runs/INDEX.json` keeps the row of every run that
    # has existed and marks the ones whose directory has gone, which is what
    # makes a run identifiable after deletion. A row with no payload has nothing
    # left to extract and nothing left to delete, so it is not a candidate here:
    # counting it would report the same archived run as dead on every pass.
    rows = [r for r in load_index() if r.get("payload_present", True)]
    if args.include_live:
        live, dead = [], list(rows)
    else:
        live = [r for r in rows if r.get("source_build") == active]
        dead = [r for r in rows if r.get("source_build") != active]
    if args.keep:
        names = {r["directory"] for r in rows}
        missing = [k for k in args.keep if k not in names]
        if missing:
            raise SystemExit(
                f"--keep names {missing} which the index does not have. "
                "Refusing, because a typo here spares nothing and deletes the "
                f"run it was meant to spare. Known: {sorted(names)}")
        spared = [r for r in dead if r["directory"] in set(args.keep)]
        dead = [r for r in dead if r["directory"] not in set(args.keep)]
        live = live + spared

    print(f"active build: {active}")
    print("every run treated as dead (--include-live)\n" if args.include_live
          else "")
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
            # Check before deleting anything at all, including the raw output.
            verify(dest, keep)
            # The WHOLE directory goes, kept files included. They live in the
            # tracked archive now, and leaving copies behind under runs/ left
            # every dead run still registering in INDEX.json as though it were
            # present -- the same file in two places, one of them the record and
            # one of them not.
            shutil.rmtree(run)

    print(f"\n{'freed' if args.execute else 'would free'} {freed/1e9:.1f} GB")
    if not args.execute:
        print("dry run; pass --execute to do it")
        return

    # A row's payload flags come from the run directories, so archiving moves
    # them the moment a directory goes. Reindex here so the ledger says which
    # runs still have a payload; it KEEPS the row of every run whose directory
    # this pass deleted and marks it, so nothing is lost by running it. The
    # extracted products live in archive/runs/<dir>/, which is tracked.
    reindex = ROOT / "exoplasim" / "scripts" / "index_runs.py"
    if reindex.is_file():
        subprocess.run([sys.executable, str(reindex)], check=True,
                       stdout=subprocess.DEVNULL)
        print("reindexed exoplasim/runs/INDEX.json")

    # `archive/runs/README.md` IS NOT WRITTEN HERE. It is a tracked document
    # that says what the three record shapes under this directory are -- the
    # INDEX_ENTRY.json stub, the RECONSTRUCTED.json for a run that was never
    # registered, and the RECORDLESS.json enumeration -- and only one of the
    # three is anything this script knows about. Rewriting it from here reverted
    # every hand edit on the next archive pass, and stamped an "Archived <date>"
    # line whose date named no run and identified nothing.
    print(f"archive is {ARCHIVE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
