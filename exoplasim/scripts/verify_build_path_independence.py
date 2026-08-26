#!/usr/bin/env python3
"""Is the executable a function of the model source, or of where it is checked out?

    python exoplasim/scripts/verify_build_path_independence.py
    python exoplasim/scripts/verify_build_path_independence.py --res T21 --ranks 16 --jobs 4

Worldbuilding frame: this checks a property of the build of the Vesper project's
climate model, a hard fork of ExoPlaSim. Nothing here concerns the simulated
planet.

THE CLAIM UNDER TEST. `notes/audits/aocl-and-model-build-flags.md` records that
the model build is bit-reproducible: same sources, same flags, same compiler give
a byte-identical executable, which is what lets `exoplasim/binary_manifest.json`
identify a binary by its sha at all. That claim was measured twice from ONE
directory, so it could not see the term this gate adds: the absolute path of the
checkout. `-g` is not what carries it. gfortran stores the source path of every
I/O statement in .rodata for its runtime diagnostics, that string survives with
no debug information at all, and gcc(1) states that `-ffile-prefix-map` does not
reach the cpp line markers CMake's preprocess-and-scan puts it in.

WHAT IT COST, and why the gate is worth two builds. A binary compiled in a git
worktree differed from the manifest's, so `rebuild_binaries.py --verify` reported
it absent and every consumer check refused it. An agent could compile the model,
which is the one way to verify a Fortran change, and could not then run what it
had built through the registry. world-ynpx.

HOW IT CAN FAIL, which is the whole point. It stages the driver, the declaration
and the model source into TWO roots that differ in both depth and length, checks
that the two stagings are byte-identical so that a difference in the executables
can only be the path, builds each through `build_model.py` itself -- the real
driver, anchored on its own file location, so each copy builds its own root --
and compares the two shas. Different shas is the failure, and the report then
names the embedded strings that differ so the next path leak is attributable
rather than merely present.

IT NEVER PUBLISHES AND NEVER TOUCHES THE REGISTRY. Both arms run with
`--no-publish`, and each arm's `MODEL_RUN` is inside its own staging directory in
any case, so neither `exoplasim/binary_manifest.json` nor an installed executable
is reachable from here.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

from _paths import PROJECT_ROOT  # noqa: E402

# What a checkout has to hold for `build_model.py` to build the model out of it:
# the driver and its path anchor, the shared library it puts on the path, and the
# declaration the flag line comes from. Copied rather than symlinked, because a
# symlink would give both arms one path and the gate would then be unable to
# fail.
STAGED = ("exoplasim/scripts", "lib", "config",
          "vendor/exoplasim/exoplasim/plasim")

# SHTns is linked, not compiled, and its prefix reaches the compiler as an
# include path -- so it is a path the executable could record, and each arm gets
# its own name for it.
SHTNS = "vendor/shtns-install"

# TWO ROOTS THAT DIFFER IN DEPTH AND IN LENGTH. Length alone would miss a leak
# that survives as a relative path with a fixed number of `..` segments, and
# depth alone would miss one that pads to a fixed width. The real case this was
# written for -- a git worktree under `.claude/worktrees/` against the checkout
# it was made from -- differs in both.
ARMS = {"one": "one",
        "two": "a-second-checkout/nested/one/level/deeper"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage(root: Path) -> None:
    if root.exists():
        shutil.rmtree(root)
    for rel in STAGED:
        shutil.copytree(PROJECT_ROOT / rel, root / rel, symlinks=False)
    link = root / SHTNS
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to((PROJECT_ROOT / SHTNS).resolve())


def staged_tree_hash(root: Path) -> str:
    """One hash over everything staged, so the two arms can be shown equal.

    Without this the gate has no right answer: two shas that differ would be
    read as a path leak when they could equally be a stale copy. Every regular
    file under the staged directories goes in, by relative path and content.
    """
    h = hashlib.sha256()
    for rel in STAGED:
        for path in sorted((root / rel).rglob("*")):
            if path.is_file() and not path.is_symlink():
                h.update(str(path.relative_to(root)).encode())
                h.update(path.read_bytes())
    return h.hexdigest()


def build(root: Path, res: str, levels: int, ranks: int,
          profile: str | None, jobs: int) -> Path:
    cmd = [sys.executable, str(root / "exoplasim/scripts/build_model.py"),
           "--res", res, "--levels", str(levels), "--ranks", str(ranks),
           "--jobs", str(jobs), "--no-publish", "--print-path"]
    if profile:
        cmd += ["--profile", profile]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout or "")
        sys.stderr.write(r.stderr or "")
        raise SystemExit(f"the arm under {root} did not build")
    return Path(r.stdout.strip().splitlines()[-1])


def embedded_paths(exe: Path, root: Path) -> list[str]:
    """Strings in the executable that name this arm's root, if any.

    Read with `strings` when it is there and with a scan of the bytes when it is
    not, so the report does not depend on binutils being installed.
    """
    needle = str(root).encode()
    data = exe.read_bytes()
    out, start = [], 0
    while True:
        i = data.find(needle, start)
        if i < 0:
            break
        lo = i
        while lo > 0 and 32 <= data[lo - 1] < 127:
            lo -= 1
        hi = i
        while hi < len(data) and 32 <= data[hi] < 127:
            hi += 1
        out.append(data[lo:hi].decode("ascii", "replace"))
        start = hi
    return sorted(set(out))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build the same model source at two different absolute "
                    "paths and compare the executable shas.")
    ap.add_argument("--res", default="T21",
                    help="rung to build the arms at (default T21, the cheapest)")
    ap.add_argument("--levels", type=int, default=10)
    ap.add_argument("--ranks", type=int, default=16)
    ap.add_argument("--profile", default=None,
                    help="build profile (default: model.compile_flags.profile)")
    ap.add_argument("--jobs", type=int, default=4,
                    help="parallel compile jobs per arm (default 4). The arms "
                         "build one after the other, so this is the whole load.")
    ap.add_argument("--work", type=Path,
                    default=PROJECT_ROOT / "vendor/exoplasim/build/path-independence",
                    help="where the two stagings go (default under the build "
                         "root, which git does not track)")
    ap.add_argument("--keep", action="store_true",
                    help="leave the stagings behind for inspection")
    a = ap.parse_args()

    roots = {name: a.work / rel for name, rel in ARMS.items()}
    for root in roots.values():
        stage(root)

    trees = {name: staged_tree_hash(root) for name, root in roots.items()}
    if len(set(trees.values())) != 1:
        for name, h in trees.items():
            print(f"  {name}: staged tree {h[:16]}")
        raise SystemExit("the two stagings differ, so the comparison below "
                         "could not attribute a difference to the path")
    print(f"staged tree {next(iter(trees.values()))[:16]} in both arms")

    exes, shas = {}, {}
    for name, root in roots.items():
        exes[name] = build(root, a.res, a.levels, a.ranks, a.profile, a.jobs)
        shas[name] = sha256(exes[name])
        print(f"  {name}: {shas[name][:16]}  {exes[name]}")

    ok = len(set(shas.values())) == 1
    if ok:
        print(f"\n{a.res} l{a.levels} p{a.ranks} is byte-identical from both "
              f"paths: the sha is a function of the source, not the checkout")
    else:
        print("\nTHE EXECUTABLE DEPENDS ON WHERE IT WAS BUILT. The two arms "
              "hold the same source and differ only in path, so one manifest "
              "cannot describe both checkouts. What each binary records of its "
              "own root:")
        for name, root in roots.items():
            found = embedded_paths(exes[name], root)
            print(f"  {name}: {len(found)} strings name {root}")
            for line in found[:5]:
                print(f"    {line}")

    if not a.keep:
        shutil.rmtree(a.work, ignore_errors=True)
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
