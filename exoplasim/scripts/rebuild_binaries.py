#!/usr/bin/env python3
"""Rebuild every ExoPlaSim executable and record which patches each contains.

    python exoplasim/scripts/rebuild_binaries.py            # rebuild all
    python exoplasim/scripts/rebuild_binaries.py --verify   # check, build nothing

ExoPlaSim compiles a separate executable per (resolution, layers, ranks) triple.
Patching the model source and running rebuilds **only the configuration you are
running**; every other binary keeps the old code, silently, until something asks
for it. That is `notes/failure-modes.md` class 11 and it fired three times in one
day: the T42 16-rank binary was three days stale, the star-cycle tree was four
days stale and had no 16-rank binary at all, and of six executables on disk
exactly one was newer than the patched source.

That failure is survivable only by luck. The ozone patch added a *namelist key*,
which an old binary cannot parse, so it aborts loudly in `radini_`. A patch that
changed the value or meaning of an existing quantity would have run to
convergence on unpatched physics and nothing would have said so.

So: **after any patch, rebuild everything.** This script is that operation, and
it writes `exoplasim/patches/binary_manifest.json` mapping every executable's
sha256 to the patch stack it contains and the source files it was built from.
`scripts/check_consistency.py` reads that manifest and fails when a binary on
disk is absent from it or was built from a source that has since changed.

## The other half: .venv is untracked

`.venv/` is not in git and is expected to be reinstallable. A reinstall restores
pristine ExoPlaSim and silently discards every applied patch. Nothing warns. So a
rebuild is not only for new patches -- it is mandatory after any reinstall, and
`--verify` is the cheap check that tells you whether you need one.
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

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / ".venv" / "lib" / "python3.12" / "site-packages" / "exoplasim"
SRC = PKG / "plasim" / "src"
RUN = PKG / "plasim" / "run"
BIN = PKG / "plasim" / "bin"
PATCHES = ROOT / "exoplasim" / "patches"
MANIFEST = PATCHES / "binary_manifest.json"

# Patches that must be resident in the source for a build to be correct. The
# star-cycle patch is deliberately NOT here: it is applied and reversed around
# its own build by build_star_cycle_exoplasim.sh, because a cycle binary and a
# steady binary are different things and only one tree can hold it at a time.
RESIDENT_PATCHES = ["exoplasim-3.4.2-ozone-band-weights.patch"]

# The matrix. NLAT must divide by ranks: 32 at T21, 64 at T42, 128 at T85.
MATRIX = [("T21", 10, 8), ("T21", 10, 16),
          ("T42", 10, 8), ("T42", 10, 16),
          ("T85", 10, 16)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def patched_sources() -> dict[str, str]:
    """Every source file any resident patch touches, with its current sha."""
    touched: set[str] = set()
    for name in RESIDENT_PATCHES:
        for line in (PATCHES / name).read_text(encoding="utf-8").splitlines():
            if line.startswith("+++ "):
                touched.add(Path(line.split()[1]).name)
    return {n: sha256(SRC / n) for n in sorted(touched) if (SRC / n).is_file()}


def resident_ok() -> tuple[bool, list[str]]:
    """Is every resident patch actually applied to the source right now?

    Tested by asking `patch` to apply it in reverse as a dry run: that succeeds
    only if the change is already present.
    """
    missing = []
    for name in RESIDENT_PATCHES:
        r = subprocess.run(
            ["patch", "--dry-run", "--reverse", "--force", "-p1",
             "--directory", str(SRC)],
            stdin=(PATCHES / name).open("rb"),
            capture_output=True)
        if r.returncode != 0:
            missing.append(name)
    return (not missing), missing


def expected(res: str, lev: int, ranks: int) -> str:
    return f"most_plasim_{res.lower()}_l{lev}_p{ranks}.x"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true",
                    help="report staleness and exit without building")
    args = ap.parse_args()

    if not SRC.is_dir():
        raise SystemExit(f"no ExoPlaSim source at {SRC}; is .venv installed?")

    ok, missing = resident_ok()
    if not ok:
        print("RESIDENT PATCHES ARE NOT APPLIED:", ", ".join(missing))
        print("\nThis is what a .venv reinstall looks like. Apply them first:")
        for m in missing:
            print(f"  patch -p1 -d {SRC} < {PATCHES / m}")
        raise SystemExit(1)
    print(f"resident patches applied: {', '.join(RESIDENT_PATCHES)}")

    sources = patched_sources()

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
        r = subprocess.run(["./compile.sh", "-n", str(ranks), "-p", str(ranks),
                            "-r", res, "-v", str(lev)],
                           cwd=PKG, capture_output=True, text=True)
        target = RUN / name
        if r.returncode != 0 or not target.is_file():
            print(r.stdout[-1500:])
            print(r.stderr[-1500:])
            raise SystemExit(f"build failed for {name}")
        built[name] = {"sha256": sha256(target),
                       "resolution": res, "layers": lev, "ranks": ranks,
                       "patches": list(RESIDENT_PATCHES),
                       "sources": sources}
        print(f"  {built[name]['sha256'][:16]}")

    payload = {
        "note": "Which patches are compiled into which executable. Generated by "
                "exoplasim/scripts/rebuild_binaries.py; do not edit. A binary "
                "absent from here has unknown provenance and must not be "
                "trusted -- rebuild rather than reason about it.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "exoplasim_version": "3.4.2",
        "resident_patches": {n: sha256(PATCHES / n) for n in RESIDENT_PATCHES},
        "binaries": built,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {MANIFEST.relative_to(ROOT)} ({len(built)} executables)")
    print("\nThe star-cycle tree is separate and is NOT rebuilt here. Run "
          "exoplasim/scripts/build_star_cycle_exoplasim.sh after this.")


if __name__ == "__main__":
    main()
