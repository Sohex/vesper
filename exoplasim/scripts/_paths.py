"""Canonical locations, resolved from this file rather than the working directory.

Scripts under `exoplasim/scripts/` may be invoked from anywhere; every default
path they expose is anchored here so behaviour does not depend on the caller's
cwd. `config/` and `source/` are project-level and shared with other components;
everything else belongs to the ExoPlaSim component.
"""

from __future__ import annotations

from pathlib import Path
import sys

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = COMPONENT_ROOT.parent

CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
SOURCE = PROJECT_ROOT / "source"

INPUTS = COMPONENT_ROOT / "inputs"
RUNS = COMPONENT_ROOT / "runs"
ANALYSIS = COMPONENT_ROOT / "analysis"
PATCHES = COMPONENT_ROOT / "patches"
# The vendored model, and where its compiled executables land.
MODEL_SRC = PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim"
MODEL_RUN = MODEL_SRC / "plasim" / "run"

# Project-level shared modules, notably the World Orogen export reader.
LIB = PROJECT_ROOT / "lib"
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))
