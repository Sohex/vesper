"""Sea water as the MODEL has it: the four numbers salinity reaches it through.

WORLDBUILDING CONTEXT: Vesper is a fictional planet and this module reads
constants out of the climate model that simulates it. Sea water below names a
modelled substance.

`icemod.f90` states the rule this module exists to keep: the ocean's salinity
reaches the model through FOUR numbers, all of them `icemod_nl` keys, and
`icemod` passes three of them to `oceanini` so the ice and ocean modules cannot
hold different sea water.

    TFREEZE  the freezing point, which sets where sea ice forms at all
    CRHOS    density; the mixed-layer heat capacity with CPS, and the snow-ice
             flooding threshold as the difference CRHOS - CRHOI
    CPS      specific heat
    CLFI     the heat of fusion of sea ice, depressed below pure ice's by brine

WHY THEY ARE NOT COPIED. Four analysis scripts carried `CRHOS = 1030.0` and
`CPS = 4180.0` as literals attributed to `oceanmod.f90`, which does not own
them. When the model moved CPS to sea water's value at S = 34.7 and its freezing
point, the copies stayed at fresh water's 4180 at about 25 C, and every slab
heat capacity built from them was 4.75 per cent too large. Nothing failed:
`assess_convergence.py` had already replaced an uncited 3990 with the copied
4180 on the argument that the model's own value is the right one, so the
correction moved it away from the model.

A run may set any of the four, so a caller that has the run directory passes it
and gets what that run used; a caller that has only the model gets the compiled
declaration. Either way the number is read, never written down.
"""

from __future__ import annotations

import re
from pathlib import Path

# The model source is the declaration. `lib/` sits beside `vendor/`.
ICEMOD_SOURCE = (Path(__file__).resolve().parents[1] / "vendor" / "exoplasim"
                 / "exoplasim" / "plasim" / "src" / "icemod.f90")

NAMES = ("TFREEZE", "CRHOS", "CPS", "CLFI")


def _number(text: str, pattern: str) -> float | None:
    m = re.search(pattern, text, re.IGNORECASE | re.MULTILINE)
    if not m:
        return None
    return float(m.group(1).replace("d", "e").replace("D", "e"))


def constants(run_dir: Path | None = None) -> dict:
    """The four, plus `source` saying where each value came from.

    Raises if the model no longer declares one: a wrong slab heat capacity does
    not fail, it rescales a closure residual, so an absent declaration must not
    fall back to anything.
    """
    source = ICEMOD_SOURCE.read_text(encoding="utf-8", errors="replace")
    out: dict[str, float | str] = {}
    for name in NAMES:
        value = _number(source, rf"^\s*real\s*::\s*{name}\s*=\s*([-+0-9.eEdD]+)")
        if value is None:
            raise SystemExit(
                f"{ICEMOD_SOURCE} no longer declares {name}. It is an icemod_nl "
                "key and the value must be read from the model, never copied "
                "into the script that uses it.")
        out[name] = value
    out["source"] = f"{ICEMOD_SOURCE.name} compiled defaults"

    namelist = (run_dir / "icemod_namelist") if run_dir is not None else None
    if namelist is not None and namelist.is_file():
        text = namelist.read_text(encoding="utf-8", errors="replace")
        overridden = []
        for name in NAMES:
            value = _number(text, rf"^\s*{name}\s*=\s*([-+0-9.eEdD]+)")
            if value is not None:
                out[name] = value
                overridden.append(name)
        if overridden:
            out["source"] = (f"{namelist.name} sets {', '.join(overridden)}; "
                             f"the rest are {ICEMOD_SOURCE.name} defaults")
    return out


def slab_heat_capacity(mixed_layer_depth_m: float,
                       run_dir: Path | None = None) -> float:
    """CRHOS * CPS * mld, J/m2/K. The mixed layer's heat capacity per unit area."""
    c = constants(run_dir)
    return float(c["CRHOS"]) * float(c["CPS"]) * float(mixed_layer_depth_m)
