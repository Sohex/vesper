"""Canonical ocean-component paths, independent of the caller's cwd."""

from __future__ import annotations

from pathlib import Path
import sys

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = COMPONENT_ROOT.parent

PLANET_CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
PIPELINE_CONFIG = PROJECT_ROOT / "config" / "pipeline.yaml"
LOOP_CONFIG = COMPONENT_ROOT / "config" / "transport_loop.yaml"
CARBON_CONFIG = COMPONENT_ROOT / "config" / "carbon_feedback.yaml"
DATA = COMPONENT_ROOT / "data"
ANALYSIS = COMPONENT_ROOT / "analysis"
RUNS = COMPONENT_ROOT / "runs"
CGENIE = PROJECT_ROOT / "vendor" / "cgenie"

LIB = PROJECT_ROOT / "lib"
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))
