#!/usr/bin/env python3
"""Build the pyfft extensions pyburn needs, for the interpreter that will use them.

    python exoplasim/scripts/build_pyfft.py            # build both, verify, report
    python exoplasim/scripts/build_pyfft.py --check    # report only, build nothing

Worldbuilding frame: a build dependency of the Vesper climate model's
postprocessor. Nothing here is about the simulated planet.

WHY THIS EXISTS. `pyburn.readfile` imports `exoplasim.pyfft` to build the
Gaussian grid, so without these two extensions no raw model output can be read
and every run fails at postprocessing. They are f2py artifacts tagged with a
CPython ABI, and git tracks only `pyfft.f90` and `pyfft991.f90` -- so an
interpreter change invalidates them and leaves the sources looking fine.

That happened. The venv moved to Python 3.14 and left `cpython-312` extensions
behind; the next run reported it as `ExoPlaSim has crashed or begun producing
garbage`, because `__init__.py` turns any postprocessing exception into
`_crash()`. The model had integrated the orbit perfectly. It cost two
diagnostic runs to find, which is why the rebuild is a registered step with a
`--check` arm rather than a thing to remember.

THIS IS THE DEFAULT PATH, and the only one. The upstream self-healing route --
`sysconfigure()` calling `compile_pyfft()` behind a `firstrun` marker -- is gone
from this fork: constructing a Model created that marker, which disabled the
rebuild permanently whether or not anything had been built, and
`compile_pyfft()` had its failure `raise` commented out so a build that produced
nothing returned normally. `compile_pyfft()` is rewritten to raise and to verify,
and this script is its entry point.

Needs meson, ninja and gfortran on PATH: numpy refuses the distutils backend on
Python >= 3.12, so `f2py -c` shells out to meson.
"""
from __future__ import annotations

import argparse
import importlib
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401  -- anchors paths from the file, not the cwd

MODULES = ("exoplasim.pyfft", "exoplasim.pyfft991")
TOOLS = ("meson", "ninja", "gfortran")


def loadable() -> dict[str, str | None]:
    """Which extensions this interpreter can actually import.

    Imports rather than looking for a file of about the right name, because the
    question is whether the running interpreter can load them.
    """
    out = {}
    for name in MODULES:
        try:
            importlib.invalidate_caches()
            importlib.import_module(name)
            out[name] = None
        except Exception as exc:
            out[name] = f"{type(exc).__name__}: {exc}"
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="report whether the extensions load and exit non-zero "
                         "if they do not; build nothing")
    args = ap.parse_args()

    print(f"interpreter: {sys.version.split()[0]} at {sys.executable}")

    if args.check:
        status = loadable()
        for name, problem in status.items():
            print(f"  {name}: {'ok' if problem is None else problem}")
        raise SystemExit(1 if any(status.values()) else 0)

    missing = [t for t in TOOLS if shutil.which(t) is None]
    if missing:
        # Named rather than worked around: a build tool absent from the host is
        # the host's to fix, and a fallback that silently produces something
        # else is how the stale extensions survived in the first place.
        raise SystemExit(
            f"missing build tools: {', '.join(missing)}. f2py needs meson and "
            f"ninja on Python >= 3.12, and gfortran to compile the sources.")

    from exoplasim import compile_pyfft
    built = compile_pyfft(verify=True)
    print(f"built and verified: {', '.join(built)}")
    for name, problem in loadable().items():
        print(f"  {name}: {'ok' if problem is None else problem}")


if __name__ == "__main__":
    main()
