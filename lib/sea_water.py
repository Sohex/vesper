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

`melting_point()` at the foot of this file is a FIFTH number and is NOT one of
the four: it is the melting point of the modelled water, pumamod's `tmelt` and a
`planet_nl` key, not the sea-water freezing point. It lives here because the two
are adjacent and were confused -- a script asked for the melting point and cited
TFREEZE's owner for it -- and a reader who comes here for one now meets both with
the difference stated.
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


# ---------------------------------------------------------------------------
# NOT ONE OF THE FOUR, and here so that the two are never confused again.
#
# TFREEZE above is the SEA WATER FREEZING POINT at the declared salinity: where
# sea ice forms at all, an `icemod_nl` key, depressed below pure water's by the
# brine. `tmelt` below is the MELTING POINT OF THE MODELLED WATER: what every
# soil, snow, sea and ice routine tests a skin temperature against, and what the
# condensation schemes switch phase at. They are different quantities and they
# have different owners, and a script that wanted the second and cited the first
# is why this accessor exists rather than a literal.
#
# THE OWNER IS pumamod, THROUGH planet_nl. `exoplasim/config/ocean_tier.yaml`
# declares it as a single-declaration handoff: `p_earth.f90` assigns it in
# `planet_ini` before reading `planet_nl`, so a run can set it, and `icemod`
# takes it through `iceini` rather than holding a compile-time copy. That is why
# the run's `planet_namelist` is consulted first and `p_earth.f90` second: the
# namelist is what the run integrated with, and the assignment is what it
# integrated with when the namelist is silent.
PLANET_SOURCE = (Path(__file__).resolve().parents[1] / "vendor" / "exoplasim"
                 / "exoplasim" / "plasim" / "src" / "p_earth.f90")


def melting_point(run_dir: Path | None = None) -> dict:
    """The modelled water's melting point in K, plus where the value came from.

    NOT the sea-water freezing point; that is `constants()["TFREEZE"]`.

    Raises if neither the run nor the planet module states it. A melting point
    that silently reverts to Earth's rescales the sea-ice part of a rebuilt heat
    content without failing, which is the same shape of defect the four above
    are read for.
    """
    namelist = (run_dir / "planet_namelist") if run_dir is not None else None
    if namelist is not None and namelist.is_file():
        value = _number(namelist.read_text(encoding="utf-8", errors="replace"),
                        r"^\s*TMELT\b\s*=\s*([-+0-9.eEdD]+)")
        if value is not None:
            return {"TMELT": value, "source": f"{namelist.name} sets TMELT"}

    source = PLANET_SOURCE.read_text(encoding="utf-8", errors="replace")
    value = _number(source, r"^\s*tmelt\b\s*=\s*([-+0-9.eEdD]+)")
    if value is None:
        raise SystemExit(
            f"{PLANET_SOURCE} no longer assigns tmelt in planet_ini. It is a "
            "planet_nl key and the melting point must be read from the run or "
            "the model, never copied into the script that uses it.")
    return {"TMELT": value, "source": f"{PLANET_SOURCE.name} planet_ini default"}
