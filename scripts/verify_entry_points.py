#!/usr/bin/env python3
"""Every script with an entry point starts: it imports, and its parser builds.

    python scripts/verify_entry_points.py
    python scripts/verify_entry_points.py --jobs 1          # serial
    python scripts/verify_entry_points.py exoplasim/scripts # a subtree only

WHAT IT ASSERTS. Every `.py` under the script directories that has a `__main__`
block answers `--help` with exit status 0, in an interpreter of its own. That is
two properties in one act: the module's import graph resolves on THIS
interpreter, and the argparse parser is constructed without raising. The second
is not hypothetical -- a strict per-build default raises while the parser is
being built when the build directory is missing, and that is the intended
failure.

WHY IT IS NOT IN `smoke_test.py` ANY MORE. `smoke_test.py` is the gate every
session runs before every commit, so its cost is multiplied by every commit in
the project; six agents paying it at once is a measurable share of this host's
load. This pass is about 136 fresh interpreters, each importing numpy, netCDF4
and whatever else the module reaches for at import time, and its own opt-out
flag used to say it "dominates runtime". A smoke test answers "is this tree
coherent" from static reads in seconds. This is a real check with real value and
it is a PRE-BUILD and PRE-PUSH check, not a per-commit one, so it has its own
entry point and its own place in `CLAUDE.md` rule 8.

WHY NOT A STATIC PASS INSTEAD. Because a static pass answers a DIFFERENT
question, and saying otherwise would be the cheap gate pretending to be the
expensive one. `smoke_test.py` parses every module, runs pyflakes F821 over it,
and asks CPython's own symbol table which loaded names nothing binds; all three
are about NAME BINDING within a source text. None of them can say whether the
import graph resolves against the packages installed in this venv, whether a
module-level statement raises, or whether an argparse default blows up while
being evaluated. Importing is a dynamic act and only an import answers it.

WHY A FRESH INTERPRETER PER SCRIPT, when one process importing all of them would
be far cheaper. Isolation, and it is load-bearing rather than tidy. Many scripts
here begin by inserting `lib/` or their own directory on `sys.path`; run them in
one process and the second script's `import rungs` succeeds because the first
one edited the path, not because the second one is correct. Module basenames
also repeat across the script directories, so `sys.modules` would hand the
second one the first one's module. A single process would report a pass this
tree has not earned, which is worse than the cost it saves.

PARALLELISM IS DELIBERATELY MODEST. This host runs several agents at once, so
this takes a quarter of the logical cores and no more; `--jobs` overrides it.
The work is subprocesses, so threads carry it and the GIL is not in the way.
Anything heavier than this takes /tmp/world.lock first -- see CLAUDE.md. This
does not, because a quarter of the cores for a few seconds is not the kind of
load that lock exists to serialise.

IT IS NOT THE COMPILE GATE'S SIBLING BY ACCIDENT.
`exoplasim/scripts/verify_model_compiles.py` is the same shape for the model
Fortran -- a process per translation unit, too expensive per commit, run before
paying for a build. Rule 8 names them together.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess
import sys

# ONE CORE PER PROBE. `--help` is not a numerical workload, but most of these
# scripts import numpy on the way to their parser, and numpy's bundled OpenBLAS
# opens a spinning thread pool sized to the logical core count. Multiplied by
# the worker count that is the whole machine, for nothing. Sized when the
# library loads, so it has to be in the environment the child starts with.
PROBE_ENV = dict(os.environ)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    PROBE_ENV.setdefault(_v, "1")
PROBE_ENV.setdefault("OMP_WAIT_POLICY", "passive")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from smoke_test import SCRIPT_DIRS  # noqa: E402  -- one declaration of the tree

ROOT = Path(__file__).resolve().parents[1]

# Entry points that must not be started: they are slow, or they mutate the tree,
# and `--help` on them is not free of side effects.
SKIP = {"rebuild_binaries.py"}


def entry_points(roots: list[Path]) -> list[Path]:
    """Every script under `roots` that has an entry point worth starting.

    `_`-prefixed files are path helpers imported by their neighbours rather than
    run, and a file with no `__main__` has no entry point to check.
    """
    out = set()
    for d in roots:
        if not d.is_dir():
            continue
        for f in d.glob("*.py"):
            if f.name in SKIP or f.name.startswith("_"):
                continue
            if "__main__" not in f.read_text(encoding="utf-8", errors="replace"):
                continue
            out.add(f)
    return sorted(out)


def start(path: Path) -> str:
    """`--help` on one script, in an interpreter of its own. "" is a pass."""
    r = subprocess.run([sys.executable, str(path), "--help"],
                       capture_output=True, text=True, timeout=120, cwd=ROOT,
                       env=PROBE_ENV)
    if r.returncode == 0:
        return ""
    tail = (r.stderr.strip().splitlines() or ["(no output)"])[-1]
    return f"{path.relative_to(ROOT)}: {tail}"


def verify(roots: list[Path], jobs: int, verbose: bool) -> list[str]:
    files = entry_points(roots)
    print(f"{len(files)} entry points, {jobs} at a time")
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results = list(pool.map(start, files))
    if verbose:
        for f, err in zip(files, results):
            print(f"  {'FAIL' if err else ' ok '}  {f.relative_to(ROOT)}")
    return [e for e in results if e]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("roots", nargs="*", type=Path,
                    help="directories to check (default: every script "
                         "directory smoke_test.py declares)")
    ap.add_argument("--jobs", type=int, default=None,
                    help="workers (default: a quarter of the logical cores, "
                         "so a shared host keeps most of itself)")
    ap.add_argument("--verbose", action="store_true",
                    help="print each entry point with its verdict")
    a = ap.parse_args()
    roots = [ROOT / r for r in a.roots] if a.roots else SCRIPT_DIRS
    jobs = a.jobs if a.jobs else max(1, (os.cpu_count() or 4) // 4)
    problems = verify(roots, jobs, a.verbose)
    for p in problems:
        print(p, file=sys.stderr)
    print(f"{len(problems)} entry points do not start" if problems
          else "every entry point starts")
    raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()
