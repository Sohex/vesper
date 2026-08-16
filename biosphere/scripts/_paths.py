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
PATCHES = COMPONENT_ROOT / "patches"
SRC = COMPONENT_ROOT / "src"
RUNS = COMPONENT_ROOT / "runs"

# LPJ-GUESS is third-party source this project modifies but does not own, so it
# lives beside ExoPlaSim and the Orogen fork rather than in this repository.
GUESS_ROOT = Path("/home/cfutro/git/lpj-guess")
GUESS_SOURCE = GUESS_ROOT / "guess_4.1"
GUESS_BUILD = GUESS_ROOT / "build"
GUESS_BINARY = GUESS_BUILD / "guess"

LIB = PROJECT_ROOT / "lib"
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))


def climatology_path(name: str = "climatology_s096") -> Path:
    """The regular (binned) climatology LPJ-GUESS is driven from.

    Defaults to the only described climatology that exists. It is on
    `precarve-unzoned` under the superseded `k2` spectrum, so anything derived
    from it is provisional; pass an explicit path once a current one is built.
    """
    return (PROJECT_ROOT / "exoplasim" / "analysis" / name
            / "baseline_regular_climatology.nc")
