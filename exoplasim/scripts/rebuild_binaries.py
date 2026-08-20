#!/usr/bin/env python3
"""Rebuild every ExoPlaSim executable and record which patches each contains.

    python exoplasim/scripts/rebuild_binaries.py            # rebuild all
    python exoplasim/scripts/rebuild_binaries.py --verify   # check, build nothing

ExoPlaSim compiles a separate executable per (resolution, layers, ranks) triple.
Patching the model source and running rebuilds **only the configuration you are
running**; every other binary keeps the old code, silently, until something asks
for it. That is `docs/src/practice/failure-modes.md` class 11 and it fired three times in one
day: the T42 16-rank binary was three days stale, the star-cycle tree was four
days stale and had no 16-rank binary at all, and of six executables on disk
exactly one was newer than the patched source.

That failure is survivable only by luck. The ozone patch added a *namelist key*,
which an old binary cannot parse, so it aborts loudly in `radini_`. A patch that
changed the value or meaning of an existing quantity would have run to
convergence on unpatched physics and nothing would have said so.

The model source is the `vendor/exoplasim` subtree, so there is no patch stack
to police any more: what the binaries are built from is whatever is in that
directory. This script records the sha of every `plasim/src` file that goes into
a build, and `--verify` reports any executable whose recorded sources no longer
match what is on disk.

So: **after any change under `vendor/exoplasim`, rebuild everything.** This
script is that operation, and it writes `exoplasim/binary_manifest.json` mapping
every executable's sha256 to the sources it was built from and the subtree
commit it came from. `scripts/check_consistency.py` reads that manifest and
fails when a binary on disk is absent from it or was built from a source that
has since changed.

Hashing the files rather than trusting the subtree commit is deliberate: a
commit says what was committed, the shas say what is on disk, and an
uncommitted edit under `vendor/` is caught rather than waved through.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "vendor" / "exoplasim" / "exoplasim"
SRC = PKG / "plasim" / "src"
RUN = PKG / "plasim" / "run"
BIN = PKG / "plasim" / "bin"
MANIFEST = ROOT / "exoplasim" / "binary_manifest.json"

# Precision is declared, never inferred: compile.sh takes it in BYTES and
# silently falls back to 4 for anything it does not recognise.
PRECISION = int(yaml.safe_load(
    (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
)["model"]["precision_bytes"])

# The matrix. NLAT must divide by ranks: 32 at T21, 64 at T42, 128 at T85.
MATRIX = [("T21", 10, 8), ("T21", 10, 16),
          ("T42", 10, 8), ("T42", 10, 16),
          ("T85", 10, 16)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model_sources() -> dict[str, str]:
    """Every file a binary is compiled FROM, with its current sha.

    This is what makes a stale executable detectable. `docs/src/practice/failure-modes.md`
    class 11 is a binary built from source that has since moved, and it has
    fired repeatedly: a T42 binary three days stale, a cycle tree four days
    stale, and of six executables on disk exactly one newer than its source.
    Vendoring the model as a subtree does not fix that -- a stale binary is
    still a stale binary -- so the check survives the migration even though the
    patch stack it used to police does not.

    Hashing the files rather than reading the subtree's commit is deliberate.
    A commit says what was COMMITTED; these shas say what is on disk right now,
    so an uncommitted edit under vendor/ is caught rather than waved through.
    """
    files = (sorted(SRC.glob("*.f90")) + sorted(SRC.glob("*.c"))
             + sorted(SRC.glob("make_*")))
    return {str(f.relative_to(PKG)): sha256(f) for f in files if f.is_file()}


def expected(res: str, lev: int, ranks: int) -> str:
    return f"most_plasim_{res.lower()}_l{lev}_p{ranks}.x"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true",
                    help="report staleness and exit without building")
    args = ap.parse_args()

    if not SRC.is_dir():
        raise SystemExit(f"no ExoPlaSim source at {SRC}; is the vendor/exoplasim subtree present?")

    # `configure.sh` writes most_compiler_mpi and friends, and ExoPlaSim only runs
    # it from `sysconfigure()`. A pip install used to trigger that; a subtree
    # checkout has never been configured, so do it here when the output is
    # absent. It is idempotent.
    if not (PKG / "most_compiler_mpi").is_file():
        print("configuring the vendored source (most_compiler_mpi absent) ...")
        import exoplasim
        exoplasim.sysconfigure()

    rev = subprocess.run(["git", "-C", str(ROOT), "log", "-1", "--format=%h %s",
                          "--", "vendor/exoplasim"],
                         capture_output=True, text=True)
    print(f"model source: vendor/exoplasim @ {rev.stdout.strip() or 'unknown'}")
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain",
                            "--", "vendor/exoplasim"],
                           capture_output=True, text=True).stdout.strip()
    if dirty:
        print(f"  UNCOMMITTED changes under vendor/exoplasim ({len(dirty.splitlines())} paths). "
              "The build will use them; the manifest records their shas.")

    sources = model_sources()

    if args.verify:
        stale, absent = [], []
        prior = (json.loads(MANIFEST.read_text(encoding="utf-8"))
                 if MANIFEST.is_file() else {"binaries": {}})
        for exe in sorted(RUN.glob("most_plasim_*.x")):
            rec = prior["binaries"].get(exe.name)
            if rec is None or rec.get("sha256") != sha256(exe):
                absent.append(exe.name)
            elif rec.get("sources") != sources:
                stale.append(exe.name)
        if absent:
            print("NOT IN MANIFEST (provenance unknown):", ", ".join(absent))
        if stale:
            print("BUILT FROM OLDER SOURCE:", ", ".join(stale))
        if not absent and not stale:
            print(f"all {len(list(RUN.glob('most_plasim_*.x')))} binaries current")
        raise SystemExit(1 if (absent or stale) else 0)

    removed = 0
    for d in (RUN, BIN):
        for exe in sorted(d.glob("most_plasim_*.x")):
            exe.unlink()
            removed += 1
    print(f"removed {removed} existing executables\n")

    built = {}
    for res, lev, ranks in MATRIX:
        name = expected(res, lev, ranks)
        print(f"building {name} ...", flush=True)
        # `-n` is ncpus and `-p` is PRECISION IN BYTES. This passed the rank
        # count to both, so an 8-rank build got -p 8 and was double precision by
        # accident, while a 16-rank build got -p 16, matched no case in
        # compile.sh, and fell back to its default of 4. Every 16-rank binary
        # was therefore single precision while config/planet.yaml declared 8,
        # and ExoPlaSim quietly recompiled at -r8 on first use -- so the
        # manifest described binaries that never ran a single orbit.
        r = subprocess.run(["./compile.sh", "-n", str(ranks), "-p", str(PRECISION),
                            "-r", res, "-v", str(lev)],
                           cwd=PKG, capture_output=True, text=True)
        target = RUN / name
        if r.returncode != 0 or not target.is_file():
            print(r.stdout[-1500:])
            print(r.stderr[-1500:])
            raise SystemExit(f"build failed for {name}")
        built[name] = {"sha256": sha256(target),
                       "resolution": res, "layers": lev, "ranks": ranks,
                       "sources": sources}
        print(f"  {built[name]['sha256'][:16]}")

    payload = {
        "note": "Which model source each executable was compiled from. Generated "
                "by exoplasim/scripts/rebuild_binaries.py; do not edit. A binary "
                "absent from here has unknown provenance and must not be "
                "trusted -- rebuild rather than reason about it.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "exoplasim_version": "3.4.2",
        "model_source": subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:vendor/exoplasim"],
            capture_output=True, text=True).stdout.strip() or None,
        "binaries": built,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {MANIFEST.relative_to(ROOT)} ({len(built)} executables)")
    print("\nThe stellar cycle is a namelist switch (nsolcycle) on these same "
          "binaries; there is no separate cycle executable to build.")


if __name__ == "__main__":
    main()
