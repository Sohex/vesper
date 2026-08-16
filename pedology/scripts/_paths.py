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


def climatology_path(name: str = "climatology_s096") -> Path:
    return (PROJECT_ROOT / "exoplasim" / "analysis" / name
            / "baseline_regular_climatology.nc")
