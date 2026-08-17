"""Canonical locations, resolved from this file rather than the working directory.

`config/` and `source/` are project-level and shared with other components;
`config/dust.yaml`, `data/` and `analysis/` belong to the aeolian component.
"""

from __future__ import annotations

from pathlib import Path
import sys

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = COMPONENT_ROOT.parent

CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
SOURCE = PROJECT_ROOT / "source"

DUST_CONFIG = COMPONENT_ROOT / "config" / "dust.yaml"
DATA = COMPONENT_ROOT / "data"
ANALYSIS = COMPONENT_ROOT / "analysis"

# Project-level shared modules: the export reader, the mesh-to-grid integration,
# the climatology resolver and the provenance check.
LIB = PROJECT_ROOT / "lib"
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))
