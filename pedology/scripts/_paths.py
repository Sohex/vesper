"""Canonical locations for the pedology component, resolved from this file.

Mirrors the other components' `_paths.py`: scripts may be invoked from anywhere,
so nothing depends on the caller's working directory.
"""

from __future__ import annotations

from pathlib import Path
import sys

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = COMPONENT_ROOT.parent

CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
SOURCE = PROJECT_ROOT / "source"

PEDOGENESIS = COMPONENT_ROOT / "config" / "pedogenesis.yaml"
DATA = COMPONENT_ROOT / "data"
ANALYSIS = COMPONENT_ROOT / "analysis"

LIB = PROJECT_ROOT / "lib"
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))


def climatology_path(name: str | None = None) -> Path:
    """The BASELINE climatology: the run on the full surface fields.

    Delegates to `lib/paths.py`, which is the one copy. This module
    kept its own, and so did biosphere, `surface_water.py` and two
    exoplasim builders; they did not stay in step, and three of them
    were still naming the superseded `climatology_s096` when a
    baseline re-run ran them for the first time in months.

    NOT every step in this component wants it. `config/pipeline.yaml`
    names one of the two climatologies per step; `build_soil.py` is a
    `bootstrap_climatology` step and takes the other resolver below.
    """
    import sys as _sys
    if str(PROJECT_ROOT / "lib") not in _sys.path:
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    from paths import climatology_path as _resolve
    return _resolve(name, root=PROJECT_ROOT)


def bootstrap_climatology_path() -> Path:
    """The BOOTSTRAP climatology: the run on terrain-only surface fields.

    The soil map is built FROM it, so it cannot be built from the baseline:
    the baseline is the run on the surface fields the soil map is one of.
    Delegates to `lib/paths.py` on the same terms as the resolver above.
    """
    import sys as _sys
    if str(PROJECT_ROOT / "lib") not in _sys.path:
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    from paths import bootstrap_climatology_path as _resolve
    return _resolve(root=PROJECT_ROOT)


WEATHERING_SCHEMES = COMPONENT_ROOT / "config" / "weathering_schemes.yaml"
OUTGASSING = COMPONENT_ROOT / "config" / "outgassing.yaml"
