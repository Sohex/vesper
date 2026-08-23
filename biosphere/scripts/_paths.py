"""Canonical locations for the biosphere component, resolved from this file.

Mirrors `exoplasim/scripts/_paths.py`: scripts here may be invoked from anywhere,
so nothing depends on the caller's working directory. `config/` and `source/` are
project-level and shared; `lib/` holds the readers every component reuses.
"""

from __future__ import annotations

from pathlib import Path
import sys

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = COMPONENT_ROOT.parent

CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
SOURCE = PROJECT_ROOT / "source"

GENERATED = COMPONENT_ROOT / "generated"
RUNS = COMPONENT_ROOT / "runs"

# The CNP fork is a git subtree, on the same terms as ExoPlaSim and Orogen: the
# source read here is the source compiled and run, and its commit is repository
# provenance rather than unrecorded state in an external checkout.
GUESS_SOURCE = PROJECT_ROOT / "vendor" / "lpj-guess"
GUESS_BUILD = GUESS_SOURCE / "build"
GUESS_BINARY = GUESS_BUILD / "guess"

LIB = PROJECT_ROOT / "lib"
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))


def climatology_path(name: str | None = None) -> Path:
    """The climatology this component is driven from.

    Delegates to `lib/paths.py`, which is the one copy. This module
    kept its own, and so did biosphere, `surface_water.py` and two
    exoplasim builders; they did not stay in step, and three of them
    were still naming the superseded `climatology_s096` when a
    baseline re-run ran them for the first time in months.
    """
    import sys as _sys
    if str(PROJECT_ROOT / "lib") not in _sys.path:
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    from paths import climatology_path as _resolve
    return _resolve(name, root=PROJECT_ROOT)
