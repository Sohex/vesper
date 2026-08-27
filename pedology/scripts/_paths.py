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
    names one of the two climatologies per step, and a step that needs the
    bootstrap to EXIST reads through `best_available_climatology` below.
    """
    import sys as _sys
    if str(PROJECT_ROOT / "lib") not in _sys.path:
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    from paths import climatology_path as _resolve
    return _resolve(name, root=PROJECT_ROOT)


def bootstrap_climatology_path() -> Path:
    """The BOOTSTRAP climatology: the run on terrain-only surface fields.

    For a quantity that does not depend on the climate state. Where the answer
    does depend on it, `best_available_climatology` below is the resolver.
    Delegates to `lib/paths.py` on the same terms as the resolver above.
    """
    import sys as _sys
    if str(PROJECT_ROOT / "lib") not in _sys.path:
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    from paths import bootstrap_climatology_path as _resolve
    return _resolve(root=PROJECT_ROOT)


def best_available_climatology(override: Path | None = None):
    """The most determined climatology that exists, WITH its stage.

    Soil texture and the weathering fluxes are functions of temperature and
    runoff, so on a pass where a baseline exists the bootstrap is no longer
    the best available, only the earliest. Returns `(path, stage)`; the stage
    goes into the product's provenance record. `lib/paths.py` carries the
    argument, including why this is not the fallback the no-fallback rule
    forbids.
    """
    import sys as _sys
    if str(PROJECT_ROOT / "lib") not in _sys.path:
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    from paths import best_available_climatology as _resolve
    return _resolve(override, root=PROJECT_ROOT)


WEATHERING_SCHEMES = COMPONENT_ROOT / "config" / "weathering_schemes.yaml"
OUTGASSING = COMPONENT_ROOT / "config" / "outgassing.yaml"
