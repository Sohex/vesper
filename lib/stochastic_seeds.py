"""Stable LPJ-GUESS stochastic substreams for the Vesper run arm.

The implementation mirrors the deliberately small C++ function in
``vendor/lpj-guess/framework/guess.cpp``.  Inputs are encoded as fixed-width
little-endian integers before FNV-1a hashing.  Consequently a stream depends on
the physical cell and ecological replicate, never its MPI rank or traversal
position, while adding draws to one named process cannot advance another.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = PROJECT_ROOT / "biosphere/config/stochastic_seeds.yaml"

_FNV_OFFSET = 14695981039346656037
_FNV_PRIME = 1099511628211
_MASK64 = (1 << 64) - 1
_PARK_MILLER_MAX_STATE = 2147483646


class SeedContractError(ValueError):
    """The declared seed contract or a requested key is invalid."""


def read_declaration(path: Path = DECLARATION) -> dict:
    declaration = yaml.safe_load(path.read_text())
    expected = {
        "establishment": 1,
        "fire_occurrence": 2,
        "fire_mortality": 3,
        "background_mortality": 4,
        "disturbance": 5,
    }
    if declaration.get("contract_version") != "vesper-stochastic-streams/1":
        raise SeedContractError("unsupported stochastic stream contract")
    if declaration.get("algorithm") != "fnv1a64-components-le-v1":
        raise SeedContractError("unsupported stochastic stream algorithm")
    if declaration.get("processes") != expected:
        raise SeedContractError("stochastic process ids changed or are incomplete")
    root = declaration.get("root_seed")
    if not isinstance(root, int) or not 1 <= root <= _PARK_MILLER_MAX_STATE:
        raise SeedContractError("root_seed must be in the Park-Miller state range")
    if declaration.get("park_miller_state_range") != [1, _PARK_MILLER_MAX_STATE]:
        raise SeedContractError("declared Park-Miller state range is inconsistent")
    return declaration


def coordinate_microdegrees(value: float) -> int:
    """Round a coordinate exactly as C++ ``llround(value * 1e6)`` does."""
    scaled = float(value) * 1_000_000.0
    # Geographic coordinates in this project are not half-microdegree ties,
    # but spell out C++'s away-from-zero rule rather than inherit Python's
    # bankers' rounding.
    return int(scaled + 0.5) if scaled >= 0.0 else int(scaled - 0.5)


def _fnv_components(components: Iterable[int]) -> int:
    value = _FNV_OFFSET
    for component in components:
        unsigned = int(component) & _MASK64
        for shift in range(0, 64, 8):
            value ^= (unsigned >> shift) & 0xFF
            value = (value * _FNV_PRIME) & _MASK64
    return value


def derive_seed(root_seed: int, lon: float, lat: float, stand_id: int,
                patch_id: int, process_id: int) -> int:
    """Return one nonzero Park-Miller state for a named ecological stream."""
    values = (root_seed, coordinate_microdegrees(lon),
              coordinate_microdegrees(lat), stand_id, patch_id, process_id)
    if not 1 <= int(root_seed) <= _PARK_MILLER_MAX_STATE:
        raise SeedContractError("root_seed is outside the Park-Miller state range")
    if stand_id < 0 or patch_id < 0 or process_id < 1:
        raise SeedContractError("stand, patch and process ids must be nonnegative")
    return (_fnv_components(values) % _PARK_MILLER_MAX_STATE) + 1


def stream_seeds(root_seed: int, lon: float, lat: float, stand_id: int,
                 patch_id: int, declaration: dict | None = None) -> dict[str, int]:
    declaration = declaration or read_declaration()
    return {
        name: derive_seed(root_seed, lon, lat, stand_id, patch_id, process_id)
        for name, process_id in declaration["processes"].items()
    }
