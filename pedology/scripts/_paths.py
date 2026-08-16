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
    """The regular (binned) climatology downstream components are driven from.

    Read from `config/planet.yaml`'s `baseline_climatology` rather than
    hardcoded. It used to default to `climatology_s096`, which is on pre-carve
    terrain under the superseded `k2` spectrum and is the surface the antipodal
    carve verdict was taken from -- and six scripts across two components took
    that default silently. A stale default is worse than a missing one: it
    produces a plausible number instead of an error.

    `name` still accepts a directory under exoplasim/analysis/ for the old
    layout, so existing callers that pass one keep working.
    """
    if name is not None:
        return (PROJECT_ROOT / "exoplasim" / "analysis" / name
                / "baseline_regular_climatology.nc")
    import yaml
    config = yaml.safe_load(
        (PROJECT_ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    declared = config.get("baseline_climatology")
    if not declared:
        raise SystemExit(
            "config/planet.yaml has no `baseline_climatology`. Name one there "
            "or pass --climatology; there is deliberately no fallback.")
    return PROJECT_ROOT / declared
