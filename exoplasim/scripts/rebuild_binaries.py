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

Two lists live here. `RESIDENT_PATCHES` is what is applied to the vendored source
RIGHT NOW; `PENDING_PATCHES` is what has been authored and verified but not
applied, which is where a patch waits between being written and being merged.
A patch moves between them in the commit that applies it and rebuilds.

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
#
# WHAT IS APPLIED TO THE VENDORED SOURCE RIGHT NOW, not what ought to be. A patch
# joins this list when it is applied and the binaries are rebuilt: `resident_ok`
# reports anything here that is missing as a problem, so a patch listed before it
# is applied reports a problem that is not one, and that is how a check stops
# being read.
#
# Each entry carries its own ROOT, because they do not all strip to the same
# depth: the patches written against the plasim source use `a/radmod.f90` and
# apply in `plasim/src`, while the ones written against the package use
# `a/plasim/src/radmod.f90` or `a/makestellarspec.py` and apply one level up. A
# single root silently reported the last two as unapplied, which is why the list
# carried one entry for so long while three patches were resident.
#
# Several of these are no-ops until something turns them on, and that is
# deliberate rather than incidental. `h2osww` defaults to 1.0, `ndustrad` to 0,
# and `co2sww` to 0.0, so a rebuilt binary reproduces the physics the current
# runs already have. What residency buys is that the run which turns one on and
# the run it is compared against come from the SAME executable, which is the only
# way the comparison measures the physics rather than a rebuild.
#
# `co2sww` is the one whose default reads backwards, so it is worth saying here.
# Zero means the CO2 shortwave term is ABSENT, because upstream has no shortwave
# CO2 at all and zero is what reproduces upstream; 1.0 is the solar-weighted
# term, not the off switch.
RESIDENT_PATCHES = [
    ("exoplasim-3.4.2-ozone-band-weights.patch", "src"),
    ("exoplasim-3.4.2-energy-diagnostics.patch", "src"),
    ("exoplasim-3.4.2-denergy-accumulator.patch", "src"),
    ("exoplasim-3.4.2-prescribed-dust.patch", "src"),
    ("exoplasim-3.4.2-rayleigh-reference-grid.patch", "src"),
    ("exoplasim-3.4.2-h2o-shortwave-weight.patch", "pkg"),
    ("exoplasim-3.4.2-co2-shortwave.patch", "pkg"),
    ("exoplasim-3.4.2-aerocore-defects.patch", "src"),
    ("exoplasim-3.4.2-aerosol-deposition.patch", "src"),
    ("exoplasim-3.4.2-dust-emission.patch", "src"),
    ("exoplasim-3.4.2-aerosol-apart.patch", "src"),
    ("exoplasim-3.4.2-aerosol-longwave.patch", "src"),
    ("exoplasim-3.4.2-makestellarspec.patch", "pkg"),
    ("exoplasim-3.4.2-nlowio-broadcast.patch", "src"),
    ("exoplasim-3.4.2-arasc-output.patch", "src"),
    ("exoplasim-3.4.2-lowio-first-record.patch", "src"),
]

# Patches that are AUTHORED and verified but NOT applied to the vendored source
# yet. They exist because a patch has to be written, reviewed and merged before
# anyone applies it and rebuilds, and RESIDENT_PATCHES above means "applied right
# now" -- listing one there early reports a problem that is not one, which is how
# a check stops being read.
#
# A patch moves from here to RESIDENT_PATCHES in the SAME commit that applies it
# and rebuilds every binary. Nothing here affects any executable, so `--verify`
# reports this list as information and never as a failure. The point of it is
# findability: a patch file sitting in exoplasim/patches/ and named by nothing is
# a patch that gets reimplemented beside itself.
#
# Verify one without touching .venv the way resident_ok() does: copy the files it
# touches to a scratch directory, apply with `patch -p1`, and reverse.
PENDING_PATCHES: list[tuple[str, str]] = [
    # Empty: everything authored has been applied and rebuilt. An entry here means
    # a patch exists in exoplasim/patches/ and is NOT in the vendored source, which
    # is the state between authoring a patch and the rebuild that lands it.
]

ROOTS = {"src": SRC, "pkg": PKG}


# The matrix. NLAT must divide by ranks: 32 at T21, 64 at T42, 128 at T85.
MATRIX = [("T21", 10, 8), ("T21", 10, 16),
          ("T42", 10, 8), ("T42", 10, 16),
          ("T85", 10, 16)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def patched_sources() -> dict[str, str]:
    """Every source file any resident patch touches, with its current sha.

    Keyed by the path relative to the package, not by basename: two patches with
    different roots can name the same file and a basename would collide them.
    """
    touched: set[Path] = set()
    for name, root in RESIDENT_PATCHES:
        for line in (PATCHES / name).read_text(encoding="utf-8").splitlines():
            if line.startswith("+++ "):
                rel = Path(line.split()[1])
                rel = Path(*rel.parts[1:])          # strip the leading b/
                touched.add((ROOTS[root] / rel).resolve())
    return {str(f.relative_to(PKG)): sha256(f)
            for f in sorted(touched) if f.is_file()}


def resident_ok() -> tuple[bool, list[str]]:
    """Is every resident patch actually applied to the source right now?

    The whole STACK is unwound in a scratch mirror, newest first, and each patch
    must reverse cleanly in turn. Testing each patch against the live source
    independently does not work once there is more than one: a later patch that
    touches adjacent lines changes the CONTEXT an earlier patch's reverse needs,
    so the earlier one reports as missing while being perfectly present. That
    happened the moment this list went from one entry to five, and taking the
    report at face value would have meant re-applying a patch already in place.

    Unwinding in order is the exact test, and it is the only one that can fail
    for the right reason. What it is for is a `.venv` reinstall, which restores
    pristine sources and discards the lot without a word.
    """
    import tempfile

    touched: set[Path] = set()
    for name, root in RESIDENT_PATCHES:
        for line in (PATCHES / name).read_text(encoding="utf-8").splitlines():
            if line.startswith("+++ "):
                rel = Path(*Path(line.split()[1]).parts[1:])
                touched.add((ROOTS[root] / rel).resolve().relative_to(PKG))

    with tempfile.TemporaryDirectory() as td:
        mirror = Path(td)
        for rel in touched:
            src = PKG / rel
            if not src.is_file():
                continue
            (mirror / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, mirror / rel)
        roots = {"src": mirror / "plasim" / "src", "pkg": mirror}
        missing = []
        for name, root in reversed(RESIDENT_PATCHES):
            r = subprocess.run(
                ["patch", "--dry-run" if missing else "--forward",
                 "--reverse", "--force", "-p1",
                 "--directory", str(roots[root])],
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
        print("\nEither a .venv reinstall discarded them or one is new to the "
              "list. Apply them first:")
        for m in missing:
            root = ROOTS[dict(RESIDENT_PATCHES)[m]]
            print(f"  patch -p1 -d {root} < {PATCHES / m}")
        raise SystemExit(1)
    print("resident patches applied: "
          + ", ".join(n for n, _ in RESIDENT_PATCHES))
    if PENDING_PATCHES:
        print("authored but NOT applied (PENDING_PATCHES): "
              + ", ".join(n for n, _ in PENDING_PATCHES))

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
                       "patches": [n for n, _ in RESIDENT_PATCHES],
                       "sources": sources}
        print(f"  {built[name]['sha256'][:16]}")

    payload = {
        "note": "Which patches are compiled into which executable. Generated by "
                "exoplasim/scripts/rebuild_binaries.py; do not edit. A binary "
                "absent from here has unknown provenance and must not be "
                "trusted -- rebuild rather than reason about it.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "exoplasim_version": "3.4.2",
        "resident_patches": {n: sha256(PATCHES / n)
                             for n, _ in RESIDENT_PATCHES},
        "binaries": built,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {MANIFEST.relative_to(ROOT)} ({len(built)} executables)")
    print("\nThe star-cycle tree is separate and is NOT rebuilt here. Run "
          "exoplasim/scripts/build_star_cycle_exoplasim.sh after this.")


if __name__ == "__main__":
    main()
