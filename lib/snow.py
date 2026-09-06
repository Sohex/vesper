"""The material properties of the modelled snow, declared once for every column.

Worldbuilding. Vesper is an invented planet and this module is about the models
that simulate it. The snow below is a modelled substance: the layer the climate
model's `landmod` puts between the simulated atmosphere and the simulated soil,
and the layer the vegetation model's `soil.cpp` puts between the same two. Snow
is a porous arrangement of ice, and a published relation for its conductivity is
a laboratory result about a MATERIAL, which is what lets it transfer to a planet
that was never measured.

## What is declared here

TWO relations, and they are different kinds of thing. The CONDUCTIVITY is a
property of a porous arrangement of ice and is a function of the pack's density;
the SPECIFIC HEAT is a property of the ice itself and is a function of its
temperature. Neither settles the other and each has its own source. The
volumetric heat capacity both columns actually solve with is the second times
the pack's own density, which is why this module carries that product too rather
than leaving a caller to form it.

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

## The same defect one constant over

The two columns then carried two SPECIFIC HEATS of the same ice. The climate
column declared a fixed 2090 J/kg/K with no citation; the ecology column
evaluated Fukusako's linear relation in absolute temperature, which this project
holds no copy of. One of the two is a function of temperature and the other is
that function evaluated somewhere: 2090 is Fukusako's line at 276.5 K and
IAPWS-06's ice at 272.2 K, both of which are AT OR ABOVE the melting point of
the modelled snow. So the fixed value was the specific heat the modelled pack has
at the moment it melts, applied at every temperature it reaches, and at 233 K it
was 15.8 per cent high.

Neither number needed to be chosen over the other, because a better source was
already in the tree. IAPWS-06 gives the specific heat of ice Ih exactly and
`analysis/ice_properties.py` implements it against the release's own check
table; `glaciermod`'s `CPGLAC` already comes from it. The relation below is that
standard, represented in the one closed form two compiled models can restate.

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


# ---------------------------------------------------------------------------
# THE SPECIFIC HEAT OF THE ICE THE MODELLED SNOW IS MADE OF.
#
# IAPWS R10-06(2009), the equation of state 2006 for H2O ice Ih, is the source
# and `analysis/ice_properties.py` is where it is implemented: a Gibbs function
# whose second temperature derivative gives `cp` exactly, checked against every
# quantity at every state of the release's own Table 6 before it reports
# anything. `references/pdf/iapws_2009_revised-release-on-the-equation-of-state-2006-for-h2o-ice-ih.pdf`
# carries the release and `references/INDEX.md` records it as read.
#
# The modelled snowpack is ice Ih plus air, and at every density either column
# reaches the air carries under a thousandth of the mass, so a kilogram of the
# modelled pack stores what a kilogram of ice stores. The specific heat of the
# modelled snow IS ice Ih's.
#
# A Gibbs function in complex arithmetic is not something a Fortran model and a
# C++ model restate, so what they restate is the quadratic below: IAPWS-06's
# `cp` least-squares represented over the temperature span the modelled snow
# reaches. That is a DERIVATION and not a fit to data -- the standard is exact
# and this is a closed form for it, with the error of the representation stated
# and gated rather than assumed small.
#
#     cp = a T**2 + b T + c        J/kg/K, T in kelvin
# ---------------------------------------------------------------------------
IAPWS_06_SPECIFIC_HEAT = (2.89232e-3, 5.85840748, 281.225695)

SPECIFIC_HEAT_DOMAIN_K = (170.0, 273.16)
"""The span the quadratic was derived over: 170 K to the triple point.

The upper end is where ice Ih stops. The lower end is below anything the
simulated air over the modelled snow has reached -- the coldest monthly bin of
the two bounding climates reaches -75.8 degC, which is 197.4 K -- and the
representation degrades gracefully rather than sharply below it: the quadratic
is within 0.9 J/kg/K of IAPWS-06 down to 140 K and only leaves 0.5 per cent at
100 K, so a simulated excursion under the domain is a widening error and never a
discontinuity. The domain is what `SPECIFIC_HEAT_MAX_RESIDUAL_J_KG_K` is
certified over, not a bound the models clamp to.
"""

SPECIFIC_HEAT_MAX_RESIDUAL_J_KG_K = 0.554
"""What the closed form costs against the standard, over the domain above.

0.554 J/kg/K is 0.034 per cent, and it is the whole price of representing
IAPWS-06 as something two compiled models can carry. Against it, the fixed
2090 this replaced in the climate column is 68 J/kg/K out at the temperature
that column states its snow at, and Fukusako's linear relation, which the
ecology column ran, is 29.7 J/kg/K low at the melting point. The representation
error is two orders below either.

A CHECK THAT CAN FAIL: `scripts/smoke_test.py` re-derives `cp` from
`analysis/ice_properties.py`'s IAPWS-06 across the domain and refuses if the
quadratic misses by more than this. Change a coefficient, the domain or this
bound and the gate fires.
"""


def specific_heat(temperature_k: float) -> float:
    """The specific heat of the ice the modelled snow is made of, J/kg/K.

    IAPWS-06's `cp` for ice Ih at normal pressure, in the closed form both
    compiled models restate. The pressure dependence is not carried: over the
    whole overburden a modelled snowpack can put on itself it is parts per
    million, which `analysis/ice_properties.py` evaluates.
    """
    a, b, c = IAPWS_06_SPECIFIC_HEAT
    t = float(temperature_k)
    return float(a * t * t + b * t + c)


def volumetric_heat_capacity(density_kg_m3: float,
                             temperature_k: float = ADOPTED_TEMPERATURE_K) -> float:
    """The modelled pack's heat capacity per unit volume, J/m3/K.

    The pack's OWN density times the specific heat of the ice in it. The density
    is the pack's and not solid ice's, and that is a dimensional statement rather
    than an empirical one: a volumetric heat capacity is the density of the
    substance occupying the volume. Snow reaches both columns as a water
    equivalent and a pack's thickness goes as one over its density, so at fixed
    water equivalent the thermal mass does not depend on the density at all --
    the two factors cancel, and they only cancel when this product is formed
    with the pack's density. Pinning it at solid ice's broke that cancellation
    in both columns and was repaired in both.
    """
    return float(density_kg_m3) * specific_heat(temperature_k)


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

# The SPECIFIC HEAT's restatements, on the same footing and for the same reason.
#
# `t` is the temperature symbol the quadratic is written in, so the check reads
# each model's own arithmetic. The two are not the same symbol and must not be:
# the climate column's `snowcap` is one compiled scalar, so it evaluates the
# relation at a DECLARED temperature -- the same one it evaluates the
# conductivity at, so its modelled snow is one material stated at one
# temperature -- while the ecology column has a simulated air temperature in
# hand every day and evaluates the relation at it. Sharing the RELATION is what
# makes those one material; sharing a value would have made them one material
# only on the days the ecology column's snow happened to sit at 263 K.
SPECIFIC_HEAT_RESTATEMENTS = (
    {"path": "vendor/exoplasim/exoplasim/plasim/src/landmod.f90",
     "target": "zsncp", "t": "TSNOWREF",
     "what": "landini, the climate column's snow"},
    {"path": "vendor/lpj-guess/modules/soil.cpp",
     "target": "Csnow", "t": "snowtemp_k",
     "what": "update_snow_properties, the ecology column's snow"},
)

# The temperature the climate column evaluates BOTH relations at, as that model
# declares it. It is `ADOPTED_TEMPERATURE_K` and the check says so, because a
# reference temperature that drifts from the conductivity row's is two
# temperatures for one material again in the one place it would not show.
CLIMATE_COLUMN_REFERENCE_TEMPERATURE = {
    "path": "vendor/exoplasim/exoplasim/plasim/src/landmod.f90",
    "declared": "TSNOWREF",
}

# The climate column also declares DEFAULTS for `snowdiff` and `snowcap`, which
# `landini` overwrites from the namelist density. They are placeholders and are
# checked anyway: a default that has drifted from the relation is what a build
# reads if `landini` is ever bypassed. Each is held to one unit in its own last
# printed place rather than to equality.
PLACEHOLDERS = (
    {"path": "vendor/exoplasim/exoplasim/plasim/src/landmod.f90",
     "declared": "snowdiff", "density": "rhosnow", "tolerance": 1.0e-4,
     "relation": "conductivity"},
    {"path": "vendor/exoplasim/exoplasim/plasim/src/landmod.f90",
     "declared": "snowcap", "density": "rhosnow", "tolerance": 0.1,
     "relation": "volumetric_heat_capacity"},
)


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

    want_cp = IAPWS_06_SPECIFIC_HEAT
    for entry in SPECIFIC_HEAT_RESTATEMENTS:
        path = Path(root) / entry["path"]
        if not path.is_file():
            problems.append(f"{entry['path']} is not there, and it restates the "
                            f"snow specific heat relation ({entry['what']})")
            continue
        text = path.read_text(encoding="utf-8")
        t = re.escape(entry["t"])
        pattern = (rf"{re.escape(entry['target'])}\s*=\s*([-\d.eE+]+)\s*\*\s*{t}"
                   rf"\s*\*\s*{t}\s*\+\s*([-\d.eE+]+)\s*\*\s*{t}"
                   rf"\s*\+\s*([-\d.eE+]+)\s*;?")
        match = re.search(pattern, text)
        if match is None:
            problems.append(
                f"{entry['path']} no longer writes the specific heat as a "
                f"quadratic in {entry['t']} assigned to {entry['target']} "
                f"({entry['what']}); it restates lib/snow.py and cannot go "
                f"unchecked")
            continue
        got = tuple(number(g) for g in match.groups())
        if got != want_cp:
            problems.append(
                f"{entry['path']} restates IAPWS-06's specific heat of ice Ih "
                f"as {got} and the declared closed form is {want_cp} "
                f"({entry['what']})")

    path = Path(root) / CLIMATE_COLUMN_REFERENCE_TEMPERATURE["path"]
    key = CLIMATE_COLUMN_REFERENCE_TEMPERATURE["declared"]
    if path.is_file():
        found = re.search(rf"^\s*real\s*,\s*parameter\s*::\s*{key}\s*=\s*"
                          rf"([-\d.eE+]+)", path.read_text(encoding="utf-8"),
                          re.M | re.I)
        if found is None:
            problems.append(
                f"{CLIMATE_COLUMN_REFERENCE_TEMPERATURE['path']} no longer "
                f"declares {key}, and it is the temperature that column states "
                f"both of its snow material properties at")
        elif number(found.group(1)) != ADOPTED_TEMPERATURE_K:
            problems.append(
                f"{CLIMATE_COLUMN_REFERENCE_TEMPERATURE['path']} declares {key} "
                f"= {number(found.group(1))} and lib/snow.py evaluates the "
                f"conductivity at {ADOPTED_TEMPERATURE_K}. One material, one "
                f"temperature: the specific heat and the conductivity of that "
                f"column's snow would be stated at two")

    relations = {"conductivity": conductivity,
                 "volumetric_heat_capacity": volumetric_heat_capacity}
    for entry in PLACEHOLDERS:
        path = Path(root) / entry["path"]
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        declared = {}
        for name in ("declared", "density"):
            key = entry[name]
            found = re.search(rf"^\s*real\s*::\s*{key}\s*=\s*([-\d.eE+]+)",
                              text, re.M | re.I)
            declared[name] = None if found is None else number(found.group(1))
        if declared["declared"] is None or declared["density"] is None:
            problems.append(
                f"{entry['path']} no longer declares both "
                f"{entry['declared']} and {entry['density']} as "
                f"reals, so the default cannot be checked against the relation")
            continue
        want_value = relations[entry["relation"]](declared["density"])
        if abs(declared["declared"] - want_value) > entry["tolerance"]:
            problems.append(
                f"{entry['path']} defaults {entry['declared']} to "
                f"{declared['declared']} and {entry['relation']} at "
                f"{entry['density']} {declared['density']} gives "
                f"{want_value:.4f}. landini overwrites it, so this is the "
                f"value a build that bypasses landini would run")
    return problems
