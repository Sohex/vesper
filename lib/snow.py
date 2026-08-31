"""The thermal conductivity of the modelled snow, declared once for every column.

Worldbuilding. Vesper is an invented planet and this module is about the models
that simulate it. The snow below is a modelled substance: the layer the climate
model's `landmod` puts between the simulated atmosphere and the simulated soil,
and the layer the vegetation model's `soil.cpp` puts between the same two. Snow
is a porous arrangement of ice, and a published relation for its conductivity is
a laboratory result about a MATERIAL, which is what lets it transfer to a planet
that was never measured.

## Why this module exists

Two components modelled the same snow with two different conductivity relations.
The climate column derived its `snowdiff` from Fourteau et al. (2021) Eq. (18);
the ecology column computed its `Ksnow` from Sturm et al. (1997). At the density
the climate column declares, the two differ by close to a factor of two, so one
snowfall insulated the ecology column's soil about twice as well as the climate
column's. That is one material property with two values inside one project.

They cannot be reconciled by matching NUMBERS. The vegetation model's snow
density is prognostic across a compaction range, and the climate model's is a
single namelist key, so equal values at one density is a coincidence at one point
of two curves that diverge everywhere else. It is the RELATION that has to be
shared, and this module is where it is stated.

## Why Fourteau, on three checks that can each fail

- Sturm et al. (1997) is a needle-probe regression. Riche and Schneebeli (2013)
  ran a long-heating needle probe, a guarded heat flux plate and a direct
  numerical simulation on IDENTICAL samples and concluded the simulation is the
  most reliable of the three, with a horizontally inserted needle probe wrong by
  up to a quarter either way through the anisotropy of the pack.
- Sturm's relation sits BELOW BOTH ARMS of the published kinetics bracket at
  every density the two components span, so it is not a point inside the honest
  uncertainty; it is outside it. `analysis/ice_properties.py` evaluates all three
  relations across that range and reports the comparison.
- Fourteau computes the effective conductivity on tomographic microstructures
  INCLUDING the latent heat carried by water vapour diffusing through the pore
  space, which is a real term in a pack under a temperature gradient and which a
  conduction-only computation omits.

## What is still open, and it is one-signed

The adopted row is the FAST kinetics limit, which is the UPPER endpoint of the
bracket rather than a point inside it. Fourteau's Sect. 4.1 says which limit snow
is in is unresolved; the slow limit is Calonne et al. (2011) Eq. (12) and lies
below this line at every density. So the residual uncertainty after this
declaration has a sign: **the modelled snow may conduct LESS than this module
says, and it cannot conduct more.**

## How the one declaration reaches two compiled models

A Fortran model and a C++ model cannot import this at runtime, so each carries
the adopted row as a literal and `check_restatements()` holds both to the table
here. That is `lib/rungs.py`'s arrangement and it is the reason a restatement is
allowed at all: a copy with no gate is the defect this module was written to
remove. A generated header was the alternative and covers only one of the two,
because `landmod.f90` is committed model source that is read and edited by hand.

`notes/audits/cryosphere-material-properties.md` argues the choice and the
bracket; `biosphere/config/snow_thermal.yaml` registers the ecology column's
departure from mainline LPJ-GUESS and `biosphere/scripts/snow_thermal_gate.py`
checks both halves of it.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Fourteau, Domine, Hagenmuller (2021), The Cryosphere 15, 2739-2755,
# 10.5194/tc-15-2739-2021. Equation (18): the VERTICAL effective thermal
# conductivity under the fast kinetics hypothesis, as a second-order polynomial
# in the ice volume fraction, at each of the five temperatures the paper
# simulated. `references/pdf/fourteau_2021_impact-of-water-vapor-diffusion-and-latent-heat-on-the-effective-therm.pdf`
# carries it and `references/INDEX.md` records it as read.
#
#     k = a (rho/rhoice)**2 + b (rho/rhoice) + c        W/m/K
# ---------------------------------------------------------------------------
FOURTEAU_2021_VERTICAL = {
    223.0: (2.564, -0.059, 0.0205),
    248.0: (2.172, 0.015, 0.0252),
    263.0: (1.985, 0.073, 0.0336),
    268.0: (1.883, 0.107, 0.0386),
    273.0: (1.776, 0.147, 0.0455),
}

RHOICE_F2021 = 917.0
"""The normalising ice density the volume fraction is taken against, kg/m3.

PART OF THE FIT AND NOT A CONSTANT OF THIS WORLD. The polynomial is in
rho/rhoice and 917 is the value its coefficients were regressed against, which
Fourteau states. Changing it would not describe denser ice; it would misread the
fit. It is NOT `icemod`'s density of the modelled SEA ice and NOT `soil.h`'s
`ice_density`, and none of the three may be deduplicated into the others.
"""

ADOPTED_TEMPERATURE_K = 263.0
"""Which of the five rows both columns evaluate.

The middle of the five, and a DECLARATION rather than a tuning. Across the whole
223 to 273 K span the value at the declared snow density moves by a few per
cent, against a factor of two between published relations, so this is not what
the constant turns on. Making it a per-cell function would be a change to the
snow scheme rather than to this relation.
"""


def conductivity(density_kg_m3: float, temperature_k: float = ADOPTED_TEMPERATURE_K) -> float:
    """The modelled snow's thermal conductivity in W/m/K, from its density.

    The fast-kinetics arm, which is the bracket's upper endpoint: the modelled
    snow may conduct less than this and cannot conduct more.
    """
    a, b, c = FOURTEAU_2021_VERTICAL[temperature_k]
    vf = density_kg_m3 / RHOICE_F2021
    return float(a * vf * vf + b * vf + c)


def resistance(water_equivalent_m: float, density_kg_m3: float,
               water_density_kg_m3: float = 1000.0) -> float:
    """The conductive resistance of a snow layer, m2 K/W, from its water equivalent.

    The quantity both columns actually use the conductivity for: a layer's
    thickness over its conductivity, which is what sets the temperature drop
    across the pack per unit ground heat flux. Carried here because the
    conductivity and the thickness both move with the density and in opposite
    directions, and quoting either alone says nothing about what a change is
    worth.
    """
    thickness = water_equivalent_m * water_density_kg_m3 / density_kg_m3
    return float(thickness / conductivity(density_kg_m3))


# ---------------------------------------------------------------------------
# The restatements. Neither is removable: a Fortran model and a C++ model cannot
# import this module at runtime, so each writes the adopted row out as a literal
# and is held to this table instead.
#
# `target` is the variable the quadratic is assigned to and `vf` the volume
# fraction it is written in, so the check reads the source's own arithmetic
# rather than a comment beside it.
# ---------------------------------------------------------------------------
RESTATEMENTS = (
    {"path": "vendor/exoplasim/exoplasim/plasim/src/landmod.f90",
     "target": "snowdiff", "vf": "zsnowvf",
     "what": "landini, the climate column's snow"},
    {"path": "vendor/lpj-guess/modules/soil.cpp",
     "target": "Ksnow", "vf": "snowdens_vf",
     "what": "update_snow_properties, the ecology column's snow"},
)

# The climate column also declares a DEFAULT for `snowdiff`, which `landini`
# overwrites from the namelist density. It is a placeholder and is checked
# anyway: a default that has drifted from the relation is what a build reads if
# `landini` is ever bypassed, and it is written at four decimal places, so it is
# held to one unit in its last printed place rather than to equality.
PLACEHOLDER = {
    "path": "vendor/exoplasim/exoplasim/plasim/src/landmod.f90",
    "declared": "snowdiff", "density": "rhosnow", "tolerance": 1.0e-4,
}


def check_restatements(root) -> list[str]:
    """Every literal restatement of the relation, against this table. Empty when they agree.

    A CHECK WITH A RIGHT ANSWER: each restatement has to carry exactly the three
    coefficients of the adopted row and exactly the normalising ice density the
    fit was made against, and any disagreement is reported as both numbers.
    Takes the repository root rather than resolving one, so this module keeps
    knowing nothing but the relation.
    """
    import re
    from pathlib import Path

    want = FOURTEAU_2021_VERTICAL[ADOPTED_TEMPERATURE_K]
    problems: list[str] = []

    def number(token: str) -> float:
        return float(token.rstrip("."))

    for entry in RESTATEMENTS:
        path = Path(root) / entry["path"]
        if not path.is_file():
            problems.append(f"{entry['path']} is not there, and it restates the "
                            f"snow conductivity relation ({entry['what']})")
            continue
        text = path.read_text(encoding="utf-8")
        vf = re.escape(entry["vf"])
        pattern = (rf"{re.escape(entry['target'])}\s*=\s*([-\d.eE+]+)\s*\*\s*{vf}"
                   rf"\s*\*\s*{vf}\s*\+\s*([-\d.eE+]+)\s*\*\s*{vf}"
                   rf"\s*\+\s*([-\d.eE+]+)\s*;?")
        match = re.search(pattern, text)
        if match is None:
            problems.append(
                f"{entry['path']} no longer writes the relation as a quadratic "
                f"in {entry['vf']} assigned to {entry['target']} ({entry['what']}); "
                f"it restates lib/snow.py and cannot go unchecked")
            continue
        got = tuple(number(g) for g in match.groups())
        if got != want:
            problems.append(
                f"{entry['path']} restates Fourteau Eq. (18) as {got} and the "
                f"adopted {ADOPTED_TEMPERATURE_K:.0f} K row is {want} "
                f"({entry['what']})")
        rho = re.search(r"RHOICE_F2021\s*=\s*([\d.eE+]+)", text)
        if rho is None:
            problems.append(
                f"{entry['path']} restates the relation and declares no "
                f"RHOICE_F2021, so its volume fraction is taken against "
                f"something this check cannot see")
        elif number(rho.group(1)) != RHOICE_F2021:
            problems.append(
                f"{entry['path']} normalises the volume fraction by "
                f"{number(rho.group(1))} and the fit was made against "
                f"{RHOICE_F2021}")

    path = Path(root) / PLACEHOLDER["path"]
    if path.is_file():
        text = path.read_text(encoding="utf-8")
        declared = {}
        for name in ("declared", "density"):
            key = PLACEHOLDER[name]
            found = re.search(rf"^\s*real\s*::\s*{key}\s*=\s*([-\d.eE+]+)",
                              text, re.M | re.I)
            declared[name] = None if found is None else number(found.group(1))
        if declared["declared"] is None or declared["density"] is None:
            problems.append(
                f"{PLACEHOLDER['path']} no longer declares both "
                f"{PLACEHOLDER['declared']} and {PLACEHOLDER['density']} as "
                f"reals, so the default cannot be checked against the relation")
        else:
            want_value = conductivity(declared["density"])
            miss = abs(declared["declared"] - want_value)
            if miss > PLACEHOLDER["tolerance"]:
                problems.append(
                    f"{PLACEHOLDER['path']} defaults {PLACEHOLDER['declared']} to "
                    f"{declared['declared']} and the relation at "
                    f"{PLACEHOLDER['density']} {declared['density']} gives "
                    f"{want_value:.4f}. landini overwrites it, so this is the "
                    f"value a build that bypasses landini would run")
    return problems
