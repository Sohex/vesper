#!/usr/bin/env python3
"""What the modelled sea ice and its snow are made of, computed from IAPWS-06 rather than quoted.

    python analysis/ice_properties.py

WORLDBUILDING CONTEXT: Vesper is an invented planet and this file is about the
climate model that simulates it. The ice and snow below are modelled
substances. IAPWS-06 is a laboratory standard for the properties of a
MATERIAL -- water ice Ih -- and a material does not change when the star does,
which is the whole reason it transfers here.

## The question

`icemod.f90` carries six compile-time constants of the modelled sea ice and its
snow, each an Earth measurement standing where nothing said it was chosen for
this world. Most of them are properties of water ice as a material: they depend
on temperature, density and, for the ice, on brine content -- not on which star
the planet orbits and not, at the pressures involved, on how hard it pulls. So
the right outcome for most of the set is the same number WITH ITS SOURCE and
its range of validity, and this script is that source made computable instead of
copied.

## What is derivable and what is not, and the line between them

The model carries the sea ice as one number per cell: a thickness. It has no
salinity, no brine volume, no porosity and no vertical structure. So of the six:

- The SNOW constants are pure ice Ih plus air. Snow has no brine, and the air
  carries a negligible share of the mass, so the specific heat and the enthalpy
  of fusion of the modelled snow ARE ice Ih's, and IAPWS-06 gives them exactly.
  The snow's conductivity is not a property of ice at all but of a porous
  arrangement of it, and it is settled from the density elsewhere.
- The SEA ICE constants are not derivable here, and the reason is the model's
  and not the standard's. Density, specific heat and conductivity of sea ice are
  all functions of the brine volume, which is a function of the salinity and the
  temperature of the ice -- and this model carries neither. `icemod.f90` already
  says exactly this about the heat of fusion `CLFI`, which is DECLARED for that
  reason. The same sentence covers the other three, and what this script
  supplies for them is the PURE ICE Ih value as the bound the declaration sits
  against, so that a declared number is at least placed against a computed one.

## Gravity, which is where a super-Earth would enter and does not

Surface gravity reaches these constants through overburden pressure and nowhere
else. `overburden` below evaluates the isothermal compressibility against the
pressure under the thickest ice the model can carry, and the answer is parts per
million: this world's stronger gravity CANCELS out of the ice's material
properties. It does not cancel out of the SNOW's density, which is set by
compaction, and that is the one route by which the surface gravity reaches this
set at all -- through `landmod`'s `rhosnow`, and then through the conductivity
that now follows it.

## The check that can fail

IAPWS-06 publishes a table of numerical check values at three states, and the
release says double-precision code should reproduce the digits given. This
implementation is checked against all three states and every quantity in that
table before it reports anything, so a transcription error in a coefficient
cannot reach the output. `CHECK_TOLERANCE` is relative and is set at the
printed precision of the table.
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.paths import PROJECT_ROOT, rel  # noqa: E402
from lib import snow  # noqa: E402

ROOT = PROJECT_ROOT
PLANET = ROOT / "config" / "planet.yaml"
ICEMOD = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src" / "icemod.f90"
LANDMOD = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src" / "landmod.f90"
GLACIERMOD = (ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src"
              / "glaciermod.f90")
LPJ_SOIL_H = ROOT / "vendor" / "lpj-guess" / "modules" / "soil.h"
OUTPUT = ROOT / "analysis" / "ice_properties.json"

# ---------------------------------------------------------------------------
# IAPWS R10-06(2009), the Gibbs energy of ice Ih. Coefficients from Table 2 of
# the release, which `references/pdf/iapws_2009_revised-release-on-the-equation-of-state-2006-for-h2o-ice-ih.pdf`
# carries and which Appendix I of `references/pdf/TEOS-10_Manual.pdf` restates.
# `s0` is the IAPWS-95 value, which is the one the release's own check table is
# computed with and the one used wherever ice and liquid water are compared.
# ---------------------------------------------------------------------------
TT = 273.16          # triple-point temperature, K
PT = 611.657         # triple-point pressure, Pa
P0 = 101325.0        # normal pressure, Pa

G0 = (-0.632020233335886e6, 0.655022213658955,
      -0.189369929326131e-7, 0.339746123271053e-14, -0.556464869058991e-21)
S0 = -0.332733756492168e4
T1 = complex(0.368017112855051e-1, 0.510878114959572e-1)
R1 = complex(0.447050716285388e2, 0.656876847463481e2)
T2 = complex(0.337315741065416, 0.335449415919309)
R2 = (complex(-0.725974574329220e2, -0.781008427112870e2),
      complex(-0.557107698030123e-4, 0.464578634580806e-4),
      complex(0.234801409215913e-10, -0.285651142904972e-10))

CHECK_TOLERANCE = 1.0e-9
"""Relative agreement required against every entry of the release's Table 6.

FIXED BEFORE THE COMPARISON WAS RUN. The release states that the digits it
prints "can reasonably be expected to be reproduced by rounded results from
double-precision code", and it prints twelve significant figures, so a bar at
1e-9 is looser than the table's own precision and far tighter than any
transcription error could hide inside.
"""

# Table 6 of the release: three states, and every quantity it tabulates. This is
# the fixed point the implementation is checked against; nothing else in this
# file is trusted until all of it reproduces.
CHECK_TABLE = [
    {"T": TT, "p": PT, "label": "triple point",
     "g": 0.611784135, "g_p": 0.109085812737e-2, "g_T": 0.122069433940e4,
     "g_pp": -0.128495941571e-12, "g_Tp": 0.174387964700e-6,
     "g_TT": -0.767602985875e1, "h": -0.333444253966e6, "s": -0.122069433940e4,
     "cp": 0.209678431622e4, "rho": 0.916709492200e3},
    {"T": 273.152519, "p": P0, "label": "normal pressure melting point",
     "g": 0.101342740690e3, "g_p": 0.109084388214e-2, "g_T": 0.122076932550e4,
     "g_pp": -0.128485364928e-12, "g_Tp": 0.174362219972e-6,
     "g_TT": -0.767598233365e1, "h": -0.333354873637e6, "s": -0.122076932550e4,
     "cp": 0.209671391024e4, "rho": 0.916721463419e3},
    {"T": 100.0, "p": 100.0e6, "label": "100 K at 100 MPa",
     "g": -0.222296513088e6, "g_p": 0.106193389260e-2, "g_T": 0.261195122589e4,
     "g_pp": -0.941807981761e-13, "g_Tp": 0.274505162488e-7,
     "g_TT": -0.866333195517e1, "h": -0.483491635676e6, "s": -0.261195122589e4,
     "cp": 0.866333195517e3, "rho": 0.941678203297e3},
]

# ---------------------------------------------------------------------------
# Yen (1981), CRREL Report 81-10, the thermal conductivity of PURE ice as a
# function of temperature. IAPWS-06 is a Gibbs function and gives density,
# specific heat and compressibility; it says NOTHING about conductivity, which
# is a transport property. So the bound the modelled sea ice's declared
# conductivity sits against is Yen's and not the standard's.
#
# Yen regressed the polycrystalline-ice measurements of Jakob and Erk (1929),
# Powell (1958), Ratcliffe (1962), Dean and Timmerhaus (1963), Wolfe and Thieme
# (1964), Dillard and Timmerhaus (1966) and Ashworth (1972) on the form
# `lambda = a exp(b T)`, T in kelvin, and reports three arms in his Table 3
# because the data have a gap between 150 and 195 K. Eq. (33) is the whole-range
# arm, which he recommends for practical use because it has the highest
# correlation coefficient; the `> 195 K` arm is the one fitted where the
# modelled ice actually lives, and the two are carried as a bracket rather than
# collapsed, because their disagreement at the melting point is the honest width
# of "pure ice's conductivity".
# ---------------------------------------------------------------------------
YEN_1981_ICE_CONDUCTIVITY = {
    "eq_33_all_temperatures": (9.828, -0.0057, 0.9313),
    "table_3_above_195_k": (6.727, -0.0041, 0.5962),
}
"""Yen's `a`, `b` and his correlation coefficient, for `lambda = a exp(b T)`."""

YEN_1981_SEA_ICE_LAMBDA_I = 2.09
"""The pure-ice conductivity Yen's OWN sea-ice model uses over 0 to -20 C.

Carried because it is what his Figure 23's sea-ice curves were computed with,
so the band those curves span is only comparable against this number. W/m/K.

It is 1.3 per cent above what his own Eq. (33) gives at the melting point, which
is the width of his internal inconsistency about the conductivity of pure ice
and is the reason the temperature a declaration is evaluated at has to be
stated rather than inherited.
"""

YEN_1981_AIR_CONDUCTIVITY = 2.51e-2
"""The conductivity of air Yen's Eq. (70) is evaluated with, W/m/K."""

YEN_1981_RHOICE = 917.0
"""The pure-ice density Yen's air-fraction relations are written against.

Part of THOSE relations and not a property of this model. `icemod`'s CRHOI is
the density of the modelled SEA ice and `RHOICE_F2021` is the normalising
density of a snow conductivity fit; the three are the same order and must not be
deduplicated into one another. kg/m3.
"""

# ---------------------------------------------------------------------------
# The snow conductivity bracket. Snow's conductivity is not a property of ice
# but of a porous arrangement of it, and the two arms below differ in WHICH
# PROCESSES they count, not in how well they were measured.
#
# Fourteau et al. (2021) Eq. (18) computes the effective conductivity on
# tomographic microstructures WITH the latent heat carried by water vapour
# diffusing through the pore space, under the FAST kinetics limit in which the
# vapour is at saturation at every ice surface. It is a quadratic in the ice
# volume fraction at each of five temperatures, and `landmod` derives `snowdiff`
# from the 263 K row of it.
#
# Calonne et al. (2011) Eq. (12) computes the same quantity on the same kind of
# data with conduction through ice and interstitial air ONLY. Fourteau's Sect.
# 2.1 identifies that as the SLOW kinetics limit -- the case where deposition is
# too slow for latent heat to reach either the temperature field or the
# conduction -- and treats it as the lower bound of the same bracket. Calonne's
# is a quadratic in the DENSITY rather than the volume fraction, and its fit was
# made on his entire sample set at 271 K, which Fourteau states.
#
# Neither paper claims to know which limit snow is in. Fourteau's Sect. 4.1 says
# so in as many words, which is why this is reported as a bracket.
# ---------------------------------------------------------------------------
CALONNE_2011_SNOW = (2.5e-6, -1.23e-4, 0.024)
"""Eq. (12), `k = a rho**2 + b rho + c`, rho in kg/m3 and k in W/m/K.

Fitted so that k goes to the conductivity of air at zero snow density. The
correlation coefficient over his 30 samples is 0.985 and the standard deviation
of the residuals is 0.025 W/m/K, which is the scatter to compare any claimed
effect against before believing it.
"""
CALONNE_2011_RESIDUAL_SD = 0.025
CALONNE_2011_FIT_TEMPERATURE_K = 271.0

# `lib/snow.py` IS THE ONE DECLARATION OF THE FAST ARM and this file imports it
# rather than restating it. That module is what `landmod.f90` and
# `vendor/lpj-guess/modules/soil.cpp` are both held to, so the bracket reported
# here is a bracket around the relation the two models actually run rather than
# around a third copy of it. WORLD-GJOV.
FOURTEAU_2021_SNOW = snow.FOURTEAU_2021_VERTICAL
RHOICE_F2021 = snow.RHOICE_F2021


H_LIQUID_AT_TRIPLE_POINT = 0.611783e-3 * 1.0e3
"""IAPWS-95's specific enthalpy of liquid water at the triple point, J/kg.

The IAPWS-95 reference state sets the internal energy and entropy of liquid
water to zero at the triple point, so the enthalpy there is `p*v`, a fraction of
a joule per kilogram. It is carried explicitly rather than assumed zero because
the melting enthalpy below is a DIFFERENCE of two enthalpies and both have to be
on the same reference; at this size the term is irrelevant to the answer and
stating it is what makes that visible rather than lucky.
"""


def _reduced(p: float) -> float:
    """Eq. (1)'s polynomial argument, `pi - pi0`, with both reduced by the
    triple-point pressure. It is the pressure ABOVE normal pressure, not the
    absolute pressure: at `p0` the two polynomials are their constant terms."""
    return (p - P0) / PT


def _r2(p: float) -> complex:
    x = _reduced(p)
    return R2[0] + R2[1] * x + R2[2] * x ** 2


def _g0_derivatives(p: float) -> tuple[float, float, float]:
    """g0, dg0/dp and d2g0/dp2 at pressure p."""
    x = _reduced(p)
    g0 = sum(c * x ** k for k, c in enumerate(G0))
    g0_p = sum(k * c * x ** (k - 1) / PT for k, c in enumerate(G0) if k >= 1)
    g0_pp = sum(k * (k - 1) * c * x ** (k - 2) / PT ** 2
                for k, c in enumerate(G0) if k >= 2)
    return g0, g0_p, g0_pp


def gibbs(temperature_k: float, pressure_pa: float) -> dict:
    """The ice Ih Gibbs function and the derivatives Table 6 tabulates.

    Eq. (1) of the release, differentiated analytically. The complex arithmetic
    is the release's own device for compactness and has no physical content;
    only the real part is taken, and it is taken at the end.
    """
    tau = temperature_k / TT
    g0, g0_p, g0_pp = _g0_derivatives(pressure_pa)
    r2 = _r2(pressure_pa)
    r2_p = R2[1] / PT + 2.0 * R2[2] * _reduced(pressure_pa) / PT
    r2_pp = 2.0 * R2[2] / PT ** 2

    def block(t, r):
        return r * ((t - tau) * np.log(t - tau) + (t + tau) * np.log(t + tau)
                    - 2.0 * t * np.log(t) - tau * tau / t)

    def d_dtau(t, r):
        return r * (-np.log(t - tau) + np.log(t + tau) - 2.0 * tau / t)

    def d2_dtau2(t, r):
        return r * (1.0 / (t - tau) + 1.0 / (t + tau) - 2.0 / t)

    g = g0 - S0 * TT * tau + TT * (block(T1, R1) + block(T2, r2)).real
    g_T = -S0 + (d_dtau(T1, R1) + d_dtau(T2, r2)).real
    g_TT = (d2_dtau2(T1, R1) + d2_dtau2(T2, r2)).real / TT
    g_p = g0_p + TT * block(T2, r2_p).real
    g_pp = g0_pp + TT * block(T2, r2_pp).real
    g_Tp = d_dtau(T2, r2_p).real

    return {"g": g, "g_T": g_T, "g_TT": g_TT, "g_p": g_p, "g_pp": g_pp,
            "g_Tp": g_Tp,
            "s": -g_T,
            "h": g - temperature_k * g_T,
            "cp": -temperature_k * g_TT,
            "rho": 1.0 / g_p,
            "kappa_T": -g_pp / g_p}


def verify_against_release() -> tuple[float, list[str]]:
    """Every quantity at every state of Table 6. Returns the worst relative miss."""
    failures, worst = [], 0.0
    for state in CHECK_TABLE:
        got = gibbs(state["T"], state["p"])
        for key, want in state.items():
            if key in ("T", "p", "label"):
                continue
            value = got[key]
            miss = abs(value - want) / abs(want)
            worst = max(worst, miss)
            if miss > CHECK_TOLERANCE:
                failures.append(f"  {state['label']}, {key}: {value!r} "
                                f"against {want!r} ({miss:.3e} relative)")
    return worst, failures


def declared(path: Path, name: str) -> float:
    """One Fortran declaration's default, read from the source that compiles it.

    The comparisons below are only comparisons if the declared side is the
    number the model actually carries, so it is PARSED rather than restated.
    """
    pattern = re.compile(rf"^\s*real\s*::\s*{name}\s*=\s*([-\d.eE+]+)", re.M | re.I)
    match = pattern.search(path.read_text())
    if match is None:
        raise SystemExit(f"{rel(path)} no longer declares `{name}`, so the "
                         f"comparison this file makes against it is not a "
                         f"comparison. Re-locate it before trusting anything here.")
    return float(match.group(1).rstrip("."))


def declared_parameter(path: Path, name: str) -> float:
    """One Fortran `real, parameter` default, read the same way and for the same
    reason as `declared`: what this file compares against has to be the number
    the model compiles, not a restatement of it."""
    pattern = re.compile(
        rf"^\s*real,\s*parameter\s*::\s*{name}\s*=\s*([-\d.eE+]+)", re.M | re.I)
    match = pattern.search(path.read_text())
    if match is None:
        raise SystemExit(f"{rel(path)} no longer declares parameter `{name}`, "
                         f"so the comparison this file makes against it is not "
                         f"a comparison. Re-locate it before trusting anything "
                         f"here.")
    return float(match.group(1).rstrip("."))


def lpj_constant(name: str) -> float:
    """One LPJ-GUESS constant, read the same way and for the same reason."""
    pattern = re.compile(rf"^\s*const\s+double\s+{name}\s*=\s*([-\d.eE+]+)", re.M)
    match = pattern.search(LPJ_SOIL_H.read_text())
    if match is None:
        raise SystemExit(f"{rel(LPJ_SOIL_H)} no longer declares `{name}`.")
    return float(match.group(1))


def sturm_1997_snow(density_kg_m3: float) -> float:
    """Sturm et al. (1997)'s needle-probe regression, W/m/K. NOT ADOPTED.

    `notes/audits/cryosphere-material-properties.md` argues against it and
    `vendor/lpj-guess`'s `soil.cpp` no longer runs it: WORLD-GJOV made both
    columns model snow from the one relation `lib/snow.py` declares. It stays
    here because it is what the SUPERSEDED ecology column ran and because it is
    the evidence for the choice: it sits below BOTH arms of the bracket at
    every density the two components span, so it is outside the honest
    uncertainty rather than a point inside it, and that is a number in the
    output rather than an assertion in a note.
    """
    x = density_kg_m3 / 1000.0
    if density_kg_m3 > 156.0:
        return float(0.138 - 1.01 * x + 3.233 * x * x)
    return float(0.023 + 0.234 * x)


def yen_1981_snow(density_kg_m3: float) -> float:
    """Yen (1981) Eq. (34), `2.22362 (rho/rho_water)**1.885`, W/m/K.

    An APPARENT conductivity: Yen states that a measured snow conductivity
    "includes vapor diffusion", so his curve counts the latent heat term that
    separates the two arms of the bracket rather than sitting on one side of it.
    Calonne's Sect. 3.1 reports his own purely conductive data agreeing with it,
    which is the coincidence that makes carrying it worthwhile.
    """
    return float(2.22362 * (density_kg_m3 / 1000.0) ** 1.885)


def pure_ice_conductivity(temperature_k: float, arm: str) -> float:
    """Yen (1981)'s pure ice conductivity, W/m/K, on one of his two warm arms."""
    a, b, _ = YEN_1981_ICE_CONDUCTIVITY[arm]
    return float(a * np.exp(b * temperature_k))


def bubbly_ice_factor(density_kg_m3: float) -> float:
    """Yen (1981) Eq. (37): bubbly ice's conductivity over pure ice's.

    `2 rho / (3 rhoice - rho)`, Schwerdtfeger's reduction of Maxwell's
    effective-medium result for randomly distributed spherical air inclusions,
    once the conductivity of air is dropped against the ice's. Yen gives the
    unreduced form twice -- Eq. (36) for dense snow and Eq. (70) for the bubbly
    ice inside sea ice -- and `bubbly_ice_factor_full` below is that form, kept
    so the reduction is checked rather than assumed.

    `YEN_1981_RHOICE` is the pure-ice density Yen's air-fraction relations are
    written against and is part of THEM: it is neither `icemod`'s CRHOI, the
    density of the modelled sea ice, nor `RHOICE_F2021`, the normalising density
    of a snow conductivity fit.
    """
    return float(2.0 * density_kg_m3
                 / (3.0 * YEN_1981_RHOICE - density_kg_m3))


def bubbly_ice_factor_full(density_kg_m3: float) -> float:
    """Yen (1981) Eq. (70), the Maxwell form the one above reduces from."""
    va = 1.0 - density_kg_m3 / YEN_1981_RHOICE
    li, la = YEN_1981_SEA_ICE_LAMBDA_I, YEN_1981_AIR_CONDUCTIVITY
    return float((2 * li + la - 2 * va * (li - la))
                 / (2 * li + la + va * (li - la)))


def snow_conductivity_fast(density_kg_m3: float, temperature_k: float = 263.0) -> float:
    """Fourteau (2021) Eq. (18): the fast-kinetics arm, W/m/K.

    `lib/snow.py`'s, which is the relation `landmod`'s `landini` and
    `soil.cpp`'s `update_snow_properties` both restate under a check. Delegated
    rather than reproduced: a third copy would be free to agree with neither.
    """
    return snow.conductivity(density_kg_m3, temperature_k)


def snow_conductivity_slow(density_kg_m3: float) -> float:
    """Calonne (2011) Eq. (12): the slow-kinetics arm, W/m/K.

    Conduction through ice and interstitial air, with no latent heat term. One
    temperature only: the fit is his whole sample set at 271 K.
    """
    a, b, c = CALONNE_2011_SNOW
    return float(a * density_kg_m3 * density_kg_m3 + b * density_kg_m3 + c)


# ---------------------------------------------------------------------------
# SNOW COMPACTION, and what the gravity term inside it is worth. GRAV-8.
#
# This model's snow density is the constant `rhosnow` and there is no compaction
# at all, so the one snow term that carries gravity is absent. PALADYN's is the
# formulation at this project's throughput -- Willeit and Ganopolski (2016)
# Eqs. (46) to (48), self-loading after Kojima (1967) as implemented by Pitman
# et al. (1991), fresh-snow density after Anderson (1976):
#
#     d(rho)/dt = 0.5 g rho w / eta  +  P (rho_fresh - rho) / w
#     eta       = eta_0 exp[k_T (T_0 - T_sn) + k_rho rho]
#     rho_fresh = rho_min + 1.7 (T_a - T_0 + 15)**1.5, over T_0-15 < T_a < T_0+2
#
# with `w` the water equivalent as a mass per area and `P` the snowfall reaching
# the ground. GRAVITY ENTERS ONE TERM, linearly, so the compaction RATE runs
# 1.306 times Earth's here for the same load and viscosity. The density it
# produces does not: the pack also relaxes toward the fresh-snow density at a
# rate set by the snowfall, so what the rate ratio is worth on the state depends
# on how long the modelled snow survives.
#
# THE THREE VISCOSITY CONSTANTS HAVE NO STATED DERIVATION. Willeit and
# Ganopolski give eta_0, k_T and k_rho in a table with a description and no
# source column, no range of validity and no sensitivity test, and take them
# from a 1991 technical report's implementation of a 1967 conference paper. So
# adding this scheme imports three constants of exactly the class this project
# refuses, which is why GRAV-8 asks for the term to be PRICED before it is
# added rather than the other way round. The paper's own stated limitation is
# scope: metamorphism and the effect of melting on density are neglected.
# ---------------------------------------------------------------------------
PALADYN_ETA_0 = 9.0e6      # Pa s, reference snow viscosity
PALADYN_K_T = 0.06         # 1/K
PALADYN_K_RHO = 0.02       # m3/kg
PALADYN_RHO_MIN = 50.0     # kg/m3, and the paper declares no maximum
PALADYN_FRESH_A = 1.7      # kg/m3/K**1.5
PALADYN_FRESH_EXP = 1.5
PALADYN_FRESH_OFFSET = 15.0  # K
# The freezing temperature of water the viscosity and the fresh-snow relation
# are both written against. NOT `TT`, which is the TRIPLE POINT and is 0.01 K
# above it: the two are different quantities and the paper states this one.
PALADYN_T_0 = 273.15
G_EARTH = 9.80665

# The snow water equivalent below which a cell is treated as bare, kg/m2. One
# millimetre: a numerical floor on a store that is a ratio's divisor, not a
# threshold with physical content.
SNOW_PRESENT_KG_M2 = 1.0


def paladyn_fresh_density(air_temperature_k):
    """Anderson (1976)'s fresh-snow density as PALADYN states it, kg/m3.

    Clipped to the stated validity window rather than extrapolated: the paper
    gives the relation only over `T_0 - 15 < T_a < T_0 + 2` and says nothing
    about outside it, so the ends are held rather than invented.
    """
    x = np.clip(np.asarray(air_temperature_k, dtype=float) - PALADYN_T_0 + PALADYN_FRESH_OFFSET,
                0.0, PALADYN_FRESH_OFFSET + 2.0)
    return PALADYN_RHO_MIN + PALADYN_FRESH_A * x ** PALADYN_FRESH_EXP


def paladyn_density(gravity, water_equivalent, snowfall, air_temperature,
                    snow_temperature, year_seconds, substeps=240, cycles=6):
    """Integrate PALADYN's prognostic snow density over a repeating year.

    `water_equivalent` (kg/m2), `snowfall` (kg/m2/s) and the two temperatures
    are (bin, cell) arrays from a climatology, held constant within each bin.
    The water equivalent is PRESCRIBED rather than integrated, because the
    climatology already carries what this model's own snow scheme produced, and
    the question is what a density would do on that snow rather than what a
    different mass balance would.

    Cycled until periodic: the density is a state with a memory of order the
    pack's lifetime, so a single pass through the year would carry the initial
    condition into the answer.
    """
    nbin = water_equivalent.shape[0]
    dt = year_seconds / (nbin * substeps)
    rho = paladyn_fresh_density(air_temperature[0]).copy()
    out = np.zeros_like(water_equivalent)
    for _ in range(cycles):
        for b in range(nbin):
            w = water_equivalent[b]
            fresh = paladyn_fresh_density(air_temperature[b])
            # The snow layer's own temperature is what the viscosity is a
            # function of, and this model does not report it. The surface
            # temperature capped at melting stands in for it, which is the
            # warm end of the pack and therefore the SOFT end: it makes the
            # viscosity smaller and the compaction faster than the pack mean
            # would, so the compaction reported here is an upper bound.
            tsn = np.minimum(snow_temperature[b], PALADYN_T_0)
            bare = w <= SNOW_PRESENT_KG_M2
            for _step in range(substeps):
                eta = PALADYN_ETA_0 * np.exp(PALADYN_K_T * (PALADYN_T_0 - tsn)
                                             + PALADYN_K_RHO * rho)
                load = 0.5 * gravity * rho * w / eta
                source = np.where(w > 0.0,
                                  snowfall[b] * (fresh - rho)
                                  / np.maximum(w, SNOW_PRESENT_KG_M2), 0.0)
                rho = rho + dt * (load + source)
                rho = np.where(bare, fresh, rho)
                rho = np.clip(rho, PALADYN_RHO_MIN, RHOICE_F2021)
            out[b] = rho
    return out


def price_snow_compaction(path, gravity, year_seconds, rhosnow):
    """What compaction, and the gravity inside it, are worth on this world's snow.

    Offline, against an existing climatology, and BEFORE any source change --
    which is GRAV-8's own sequencing, and the reason is that adding a prognostic
    density imports three Earth-calibrated viscosity constants with no stated
    derivation. What this answers is whether that price buys anything.

    THE QUANTITY IS THE PACK'S CONDUCTIVE RESISTANCE, `lib/snow.py:resistance`,
    and not its depth. Depth and conductivity both move with the density and in
    OPPOSITE directions, so quoting either alone says nothing; the resistance is
    what sets the temperature drop across the pack per unit ground heat flux and
    is what `landmod`'s `zdiff1` is built from.

    AND IT IS NOT ALBEDO, which is the correction this pricing makes to its own
    row. `landmod` takes the snow-covered fraction as
    `dsnow/(dsnow + snowcovz)` in WATER EQUIVALENT, so the density does not
    enter it. The one path from density to albedo is the canopy burial depth
    handed to `snowcanopymask`, and that is read only when `forhgt` is positive,
    which it is not. So at the shipped configuration the density reaches the
    modelled surface albedo through nothing at all, and what a compaction term
    would be worth there is exactly zero until BIO-33 and GRAV-7 give the canopy
    a height.
    """
    import netCDF4 as nc

    with nc.Dataset(path) as data:
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        swe = np.asarray(data["snd"][:], dtype=float) * 1000.0    # m w.e. -> kg/m2
        snowfall = np.asarray(data["prsn"][:], dtype=float) * 1000.0
        air = np.asarray(data["tas"][:], dtype=float)
        surface = np.asarray(data["ts"][:], dtype=float)

    rows, cols = np.where(land)
    swe = swe[:, rows, cols]
    snowfall = snowfall[:, rows, cols]
    air = air[:, rows, cols]
    surface = surface[:, rows, cols]
    snowy = swe > SNOW_PRESENT_KG_M2
    if not snowy.any():
        return {"climatology": rel(Path(path)),
                "verdict": "this climatology carries no snow on land, so "
                           "nothing here can be priced against it"}

    rho_here = paladyn_density(gravity, swe, snowfall, air, surface, year_seconds)
    rho_earth = paladyn_density(G_EARTH, swe, snowfall, air, surface, year_seconds)

    resist = np.vectorize(snow.resistance)
    cases = {
        "constant": (np.full_like(rho_here, rhosnow), None),
        "prognostic_earth_gravity": (rho_earth, G_EARTH),
        "prognostic_this_gravity": (rho_here, gravity),
    }
    priced = {}
    for name, (rho, g) in cases.items():
        res = resist(swe / 1000.0, rho)
        priced[name] = {
            "gravity_m_s2": g,
            "density_kg_m3": _spread(rho[snowy]),
            "physical_depth_m": _spread((swe / rho)[snowy]),
            "conductive_resistance_m2_k_w": _spread(res[snowy]),
        }
        priced[name]["_res"] = res

    def ratio(a, b):
        return _spread((priced[a]["_res"] / priced[b]["_res"])[snowy])

    out = {
        "climatology": rel(Path(path)),
        "snowy_land_bins": int(snowy.sum()),
        "land_bins": int(snowy.size),
        "scheme": "Willeit and Ganopolski (2016) Eqs. (46)-(48); self-loading "
                  "after Kojima (1967) as implemented by Pitman et al. (1991), "
                  "fresh-snow density after Anderson (1976)",
        "constants": {"eta_0_pa_s": PALADYN_ETA_0, "k_T_per_k": PALADYN_K_T,
                      "k_rho_m3_per_kg": PALADYN_K_RHO,
                      "rho_min_kg_m3": PALADYN_RHO_MIN,
                      "derivation": "none stated. The paper tabulates the three "
                                    "viscosity constants with a description and "
                                    "no source column, no range of validity and "
                                    "no sensitivity test, and takes them from a "
                                    "1991 technical report's implementation of a "
                                    "1967 conference paper. Adopting the scheme "
                                    "imports them."},
        "cases": {k: {n: v for n, v in body.items() if not n.startswith("_")}
                  for k, body in priced.items()},
        "resistance_ratios": {
            "prognostic_earth_over_constant": ratio("prognostic_earth_gravity",
                                                    "constant"),
            "this_gravity_over_earth_gravity": ratio("prognostic_this_gravity",
                                                     "prognostic_earth_gravity"),
            "prognostic_this_over_constant": ratio("prognostic_this_gravity",
                                                   "constant"),
        },
        "albedo": {
            "worth": 0.0,
            "why": "landmod takes the snow-covered fraction as "
                   "dsnow/(dsnow + snowcovz) in WATER EQUIVALENT, so the "
                   "density does not enter it. The only path from density to "
                   "albedo is the canopy burial depth handed to "
                   "snowcanopymask, read only when forhgt is positive, and it "
                   "is not. So the answer is exactly zero rather than small, "
                   "and it becomes non-zero when the canopy gets a height.",
            "peer_practice": "ClimaLand multiplies snow albedo by "
                             "min(1 - beta*(rho/rho_liq - x0), 1). Its beta is "
                             "a FREE PARAMETER with no citation, 0.97 in its "
                             "calibrated set and 0 in its uncalibrated one, so "
                             "adopting that path would import a tuned constant "
                             "rather than a mechanism. PALADYN puts no density "
                             "in albedo at all and uses a snow AGE factor as "
                             "the grain-size proxy instead.",
        },
        # CHECK THE INSTRUMENT AGAINST THE SIZE OF THE EFFECT. The gravity term
        # is read through the pack's conductive resistance, and the conductivity
        # that resistance is built from already carries a BRACKET at fixed
        # density: the fast and slow limits of the vapour deposition kinetics,
        # which no paper claims to resolve. If the gravity term is smaller than
        # that bracket, a pair of runs differing only in gravity cannot separate
        # it from a choice the model has already had to make.
        "instrument_versus_effect": _instrument_check(
            priced["prognostic_this_gravity"]["density_kg_m3"]["p50"],
            priced["prognostic_earth_gravity"]["density_kg_m3"]["p50"],
            ratio("prognostic_this_gravity", "prognostic_earth_gravity")["p50"]),
        "limitations": [
            "the forcing is a climatology, so an intermittent pack appears as a "
            "persistent thin one. Compaction is linear in the load, so a "
            "time-mean water equivalent understates the compaction of a "
            "transient deep pack and the reported densities are a floor for "
            "those cells",
            "the water equivalent is PRESCRIBED from a run that used the "
            "constant density, so this prices the density on that snow rather "
            "than on the snow a prognostic density would itself produce",
            "the snow layer temperature is proxied by the surface temperature "
            "capped at melting, which is the pack's warm end and makes the "
            "viscosity soft, so the compaction here is an upper bound",
            "PALADYN itself neglects metamorphism and the effect of melting on "
            "density, and this inherits both",
        ],
    }
    return out


def _instrument_check(density_here, density_earth, gravity_resistance_ratio):
    """Is the gravity term bigger than the uncertainty already in the relation?

    The comparison is made in ONE quantity -- the pack's conductive resistance
    -- so the two sides are commensurable. The gravity side is what moving the
    density from the Earth-gravity answer to this world's does to it. The
    instrument side is what choosing the other arm of the kinetics bracket does
    to it at ONE density, which is a choice the model has already had to make
    and cannot currently justify either way.
    """
    fast = snow_conductivity_fast(density_here)
    slow = snow_conductivity_slow(density_here)
    gravity_effect = abs(1.0 / gravity_resistance_ratio - 1.0)
    bracket_effect = abs(fast / slow - 1.0)
    return {
        "quantity": "the pack's conductive resistance, z/k",
        "gravity_effect": round(gravity_effect, 4),
        "gravity_effect_is": f"the resistance at {density_here:.1f} kg/m3 "
                             f"against {density_earth:.1f}, which is what the "
                             "1.306x compaction rate is worth on the STATE",
        "kinetics_bracket_effect": round(bracket_effect, 4),
        "kinetics_bracket_is": "Fourteau's fast arm over Calonne's slow arm at "
                               "the same density, which is the width of an "
                               "open question rather than an error bar",
        "gravity_over_bracket": round(gravity_effect / bracket_effect, 3),
        "verdict": ("the gravity term is SMALLER than the bracket the "
                    "conductivity already carries, so a pair of runs differing "
                    "only in gravity would not separate it from the arm the "
                    "model happens to evaluate"
                    if gravity_effect < bracket_effect else
                    "the gravity term is LARGER than the bracket the "
                    "conductivity already carries, so it is separable"),
    }


def _spread(values):
    values = np.asarray(values, dtype=float)
    return {"p5": round(float(np.percentile(values, 5)), 4),
            "p50": round(float(np.median(values)), 4),
            "p95": round(float(np.percentile(values, 95)), 4),
            "max": round(float(values.max()), 4)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--climatology", type=Path, default=None,
                    help="an existing climatology to price snow compaction "
                         "against, GRAV-8. Named rather than resolved because "
                         "the pricing is a question about a PARTICULAR run's "
                         "snow, and because the answer is reported with the "
                         "climatology it was measured on")
    args = ap.parse_args()

    worst, failures = verify_against_release()
    if failures:
        raise SystemExit(
            "this implementation of IAPWS-06 does not reproduce the release's "
            "own check table, so nothing computed from it means anything:\n"
            + "\n".join(failures))

    config = yaml.safe_load(PLANET.read_text())
    gravity = float(config["planet"]["gravity_m_s2"])
    sys.path.insert(0, str(ROOT / "lib"))
    import orbit  # noqa: E402  from lib/
    year_seconds = orbit.orbital_year_days(config) * 86400.0

    # The melting enthalpy of the modelled snow: liquid water minus ice Ih, both
    # on the IAPWS-95 reference, at the triple point where the liquid's value is
    # stated exactly.
    ice_triple = gibbs(TT, PT)
    fusion_j_kg = H_LIQUID_AT_TRIPLE_POINT - ice_triple["h"]

    # Ice Ih over the temperature range the modelled sea ice and its snow
    # occupy. The lower end is well below any sea surface this model produces
    # and is carried so the range of validity is stated rather than assumed.
    temperatures = [233.15, 243.15, 253.15, 263.15, 268.15, TT]
    profile = []
    for t in temperatures:
        state = gibbs(t, P0)
        profile.append({"temperature_k": round(t, 4),
                        "density_kg_m3": round(state["rho"], 4),
                        "specific_heat_j_kg_k": round(state["cp"], 4)})

    # Gravity, priced. The overburden under a column of ice is rho*g*h, and the
    # density response to it is the isothermal compressibility. `hmax` is a
    # deliberately generous column: the model's own maximum thickness is
    # switched off, so this is the thickness at which the question would have to
    # be revisited rather than one the model enforces.
    hmax = 20.0
    state = gibbs(263.15, P0)
    overburden_pa = state["rho"] * gravity * hmax
    overburden = {
        "gravity_m_s2": gravity,
        "column_m": hmax,
        "overburden_pa": round(overburden_pa, 1),
        "isothermal_compressibility_pa": float(f"{state['kappa_T']:.6g}"),
        "relative_density_change": float(f"{state['kappa_T'] * overburden_pa:.3g}"),
        "earth_relative_density_change": float(
            f"{state['kappa_T'] * state['rho'] * 9.80665 * hmax:.3g}"),
        "verdict": "the surface gravity CANCELS out of the material properties "
                   "of the modelled ice. Under a column far thicker than this "
                   "model carries, the density response to overburden is parts "
                   "per million, and the specific heat and conductivity "
                   "responses are smaller still. Gravity reaches this set "
                   "through the SNOW's density, which is set by compaction, and "
                   "through nothing else.",
    }

    # The conductivity of pure ice, which IAPWS-06 does not give: a Gibbs
    # function carries no transport property. Yen (1981) is the bound the
    # modelled sea ice's declared conductivity sits against, and the two arms
    # are carried as a bracket rather than averaged.
    ckapi = declared(ICEMOD, "CKAPI")
    conductivity_profile = []
    for t in temperatures:
        arms = {arm: pure_ice_conductivity(t, arm) for arm in YEN_1981_ICE_CONDUCTIVITY}
        conductivity_profile.append({
            "temperature_k": round(t, 4),
            "eq_33_all_temperatures_w_m_k": round(arms["eq_33_all_temperatures"], 4),
            "table_3_above_195_k_w_m_k": round(arms["table_3_above_195_k"], 4),
            "declared_over_pure_ice": round(ckapi / max(arms.values()), 4),
        })
    ice_conductivity = {
        "source": "Yen (1981), CRREL Report 81-10, Eq. (33) and Table 3",
        "form": "lambda = a exp(b T), T in kelvin, lambda in W/m/K",
        "arms": {k: {"a": v[0], "b": v[1], "correlation_coefficient": v[2]}
                 for k, v in YEN_1981_ICE_CONDUCTIVITY.items()},
        "profile": conductivity_profile,
        "declared_w_m_k": ckapi,
        "verdict": "the declared conductivity of the modelled sea ice sits "
                   "BELOW pure ice's over the whole range, which is the side "
                   "brine puts it on: Yen's Eq. (71) subtracts a brine term "
                   "from bubbly pure ice, so sea ice is the less conductive of "
                   "the two. Yen's own Figure 23 computes that band and the "
                   "declared value sits at the top of it, against his "
                   "zero-salinity curves -- which is the ice this model "
                   "carries, since it has no salinity and no brine volume. So "
                   "the declaration is the brine-free limit of the bound and "
                   "not a value for salty ice.",
    }

    # The snow conductivity bracket. Two arms, differing in which processes they
    # count rather than in how well they were measured, and neither paper claims
    # to know which limit snow is in.
    rhosnow = declared(LANDMOD, "rhosnow")
    snowdiff = declared(LANDMOD, "snowdiff")
    tsnowref = declared_parameter(LANDMOD, "TSNOWREF")

    # THE CLOSED FORM AGAINST THE STANDARD, at the temperature the climate
    # column states its snow at. `lib/snow.py` carries the specific heat of ice
    # Ih as a quadratic because two compiled models restate it, and the
    # quadratic is a representation of exactly what is computed here. This is
    # the same check the glacier literals get below and it fails the same way:
    # if the closed form and the standard part, one of them has been edited.
    cp_closed = snow.specific_heat(tsnowref)
    cp_standard = gibbs(tsnowref, P0)["cp"]
    if abs(cp_closed - cp_standard) > snow.SPECIFIC_HEAT_MAX_RESIDUAL_J_KG_K:
        raise SystemExit(
            f"lib/snow.py's closed form gives {cp_closed:.4f} J/kg/K for the "
            f"specific heat of ice Ih at {tsnowref} K and IAPWS-06 gives "
            f"{cp_standard:.4f}; they are the same quantity and the module "
            f"declares them within "
            f"{snow.SPECIFIC_HEAT_MAX_RESIDUAL_J_KG_K} J/kg/K")
    densities = sorted({lpj_constant("snowdens_start"), rhosnow,
                        lpj_constant("snowdens_end")}
                       | {400.0})
    bracket = []
    for rho in densities:
        slow = snow_conductivity_slow(rho)
        fast = snow_conductivity_fast(rho, 263.0)
        fast_warm = snow_conductivity_fast(rho, 273.0)
        bracket.append({
            "density_kg_m3": rho,
            "slow_kinetics_w_m_k": round(slow, 4),
            "fast_kinetics_263k_w_m_k": round(fast, 4),
            "fast_kinetics_273k_w_m_k": round(fast_warm, 4),
            "fast_over_slow_263k": round(fast / slow, 4),
            "width_over_slow_arm_residual_sd": round(
                (fast - slow) / CALONNE_2011_RESIDUAL_SD, 3),
            "yen_1981_apparent_w_m_k": round(yen_1981_snow(rho), 4),
            "sturm_1997_as_lpj_guess_runs_it_w_m_k": round(sturm_1997_snow(rho), 4),
        })
    snow_bracket = {
        "slow_arm": "Calonne et al. (2011) Eq. (12), conduction through ice and "
                    "interstitial air only, his whole sample set at 271 K",
        "fast_arm": "Fourteau et al. (2021) Eq. (18), the same computation WITH "
                    "the latent heat carried by water vapour through the pore "
                    "space, at 263 K and at 273 K",
        "which_limit_applies": "unsettled in the literature. Fourteau's Sect. "
                               "4.1 says so, citing Krol and Lowe (2016) for "
                               "isothermal metamorphism looking like slow "
                               "kinetics and temperature-gradient metamorphism "
                               "looking like fast. That is why this is a "
                               "bracket and not a correction.",
        "adopted": "the FAST arm, at 263 K, by BOTH columns: lib/snow.py "
                   "declares it, landmod's landini and soil.cpp's "
                   "update_snow_properties each restate it under a check. So "
                   "the adopted relation is the bracket's upper endpoint "
                   "rather than a point inside it, and the slow arm below is "
                   "how far under it the other limit sits.",
        "declared_density_kg_m3": rhosnow,
        "declared_conductivity_w_m_k": snowdiff,
        # `landini` derives it as `rhosnow` times the specific heat at
        # `TSNOWREF`, so it FOLLOWS the density and is not independent of it.
        # Emitted here because `pedology/config/land_column_properties.yaml`
        # restates all four of these thermal constants and had no producer to be
        # held against.
        "declared_heat_capacity_j_m3_k": round(
            snow.volumetric_heat_capacity(rhosnow, tsnowref), 1),
        "declared_specific_heat_j_kg_k": round(snow.specific_heat(tsnowref), 4),
        "declared_specific_heat_source": (
            "IAPWS-06 for ice Ih at TSNOWREF, through lib/snow.py's closed "
            "form. It replaced a fixed 2090 J/kg/K that carried no citation "
            "and was that standard's value at 272.24 K, a degree below the "
            "melting point, applied at every temperature the modelled snow "
            "reaches. WORLD-A2LV."),
        "superseded_specific_heat_j_kg_k": 2090.0,
        "slow_arm_residual_sd_w_m_k": CALONNE_2011_RESIDUAL_SD,
        "per_density": bracket,
        "instrument_versus_effect": "the bracket is wider than the scatter of "
                                    "the fit that defines its lower arm at "
                                    "every density carried, so it is a real "
                                    "disagreement and not the noise of one "
                                    "regression.",
        "cross_component": "SETTLED. vendor/lpj-guess/modules/soil.cpp computed "
                           "its snow conductivity from Sturm et al. (1997), "
                           "which is BELOW the slow arm at every density here, "
                           "so the ecology column and the climate column "
                           "insulated their soil differently from the same "
                           "snowfall at the same density by more than the "
                           "bracket between the two published limits. Both now "
                           "restate lib/snow.py's one relation, held there by "
                           "snow.check_restatements, which scripts/"
                           "smoke_test.py and biosphere/scripts/"
                           "snow_thermal_gate.py both run. The Sturm column "
                           "below is what the superseded ecology column ran "
                           "and is the evidence for the choice.",
    }

    # ------------------------------------------------------------------
    # THE MODELLED GLACIAL ICE. WORLD-FG8W.
    #
    # Two constants of the ice `glaciermod` grows, both DERIVED from its own
    # declared density rather than declared beside it. The capacity is that
    # density times ice Ih's specific heat from the standard; the conductivity is
    # Yen's pure ice reduced for the air the density implies. Both are computed
    # here and the compiled literals are held to them, which is a check that can
    # fail on either side.
    #
    # NOT the sea-ice constants and not deduplicable with them: Yen's Eq. (71)
    # subtracts a brine term from exactly the bubbly ice computed here, so the
    # two sit on opposite sides of pure ice for different reasons.
    # ------------------------------------------------------------------
    rhoglac = declared(GLACIERMOD, "rhoglac")
    tglacref = declared_parameter(GLACIERMOD, "TGLACREF")
    cpglac = declared_parameter(GLACIERMOD, "CPGLAC")
    factor = bubbly_ice_factor(rhoglac)
    glac_state = gibbs(tglacref, P0)
    glac_cap = rhoglac * glac_state["cp"]
    glac_diff = pure_ice_conductivity(tglacref, "eq_33_all_temperatures") * factor
    glac_diff_arm = pure_ice_conductivity(tglacref, "table_3_above_195_k") * factor

    sweep = []
    for t in (233.15, 243.15, 253.15, 263.15, 268.15, TT):
        sweep.append({
            "temperature_k": round(t, 4),
            "heat_capacity_j_m3_k": round(rhoglac * gibbs(t, P0)["cp"], 1),
            "conductivity_eq33_w_m_k": round(
                pure_ice_conductivity(t, "eq_33_all_temperatures") * factor, 4),
            "conductivity_above_195k_w_m_k": round(
                pure_ice_conductivity(t, "table_3_above_195_k") * factor, 4),
        })
    caps = [row["heat_capacity_j_m3_k"] for row in sweep]
    diffs = [row["conductivity_eq33_w_m_k"] for row in sweep]

    glacier_failures = []
    for name, path, computed in (
            ("CPGLAC", GLACIERMOD, glac_state["cp"]),
            ("sicecap", LANDMOD, glac_cap),
            ("sicediff", LANDMOD, glac_diff)):
        literal = (declared_parameter(path, name) if name == "CPGLAC"
                   else declared(path, name))
        if abs(literal - computed) > 1.0e-4 * abs(computed):
            glacier_failures.append(
                f"{rel(path)} declares {name} = {literal!r} and this file "
                f"computes {computed:.6g} from the standard and the regression; "
                f"they are the same quantity and must agree")
    if glacier_failures:
        raise SystemExit("\n".join(glacier_failures))

    glacial_ice = {
        "density_kg_m3": rhoglac,
        "declared_temperature_k": tglacref,
        "why_this_temperature": "the melting point of the modelled snow less "
                                "ten kelvin, which is within 0.15 K of "
                                "landmod's TSNOWREF, so the two cryosphere "
                                "materials are stated at one temperature to "
                                "the precision either of them turns on",
        "specific_heat_j_kg_k": round(glac_state["cp"], 4),
        "heat_capacity_j_m3_k": round(glac_cap, 1),
        "heat_capacity_source": "IAPWS-06 at the declared temperature, times "
                                "the modelled glacial ice's own density",
        "superseded_heat_capacity_j_m3_k": 2.07e6,
        "what_the_superseded_value_was": "one thousand times ice's specific "
                                         "heat, which is LIQUID WATER's "
                                         "density. It followed nothing, so a "
                                         "bracket on rhoglac moved the ice "
                                         "orography and left the thermal mass "
                                         "of the ice behind.",
        "bubbly_factor_eq37": round(factor, 6),
        "bubbly_factor_eq70": round(bubbly_ice_factor_full(rhoglac), 6),
        "bubbly_forms_agree_to_relative": float(
            f"{abs(bubbly_ice_factor_full(rhoglac) - factor) / factor:.3g}"),
        "air_volume_fraction": round(1.0 - rhoglac / YEN_1981_RHOICE, 5),
        "conductivity_w_m_k": round(glac_diff, 6),
        "conductivity_source": "Yen (1981) Eq. (33) for pure ice at the "
                               "declared temperature, reduced by his Eq. (37) "
                               "for the air the density implies. Eq. (70) is "
                               "the unreduced Maxwell form and is carried "
                               "beside it so the reduction is checked.",
        "conductivity_on_yens_other_arm_w_m_k": round(glac_diff_arm, 6),
        "superseded_conductivity_w_m_k": 2.03,
        "brackets": {
            "yen_arms_at_the_declared_temperature": round(
                glac_diff_arm / glac_diff, 4),
            "temperature_233k_to_melting_conductivity": round(
                max(diffs) / min(diffs), 4),
            "temperature_233k_to_melting_heat_capacity": round(
                max(caps) / min(caps), 4),
            "which_dominates": "the declared TEMPERATURE, by a wide margin. "
                               "Over 233.15 K to the melting point the "
                               "conductivity moves about four times as far as "
                               "the difference between Yen's two regression "
                               "arms does at one temperature. A "
                               "temperature-dependent pair is a change to the "
                               "soil heat solver's material model rather than "
                               "to a constant, and is named rather than made "
                               "here.",
        },
        "sweep": sweep,
        "not_ckapi": "icemod's CKAPI carried the same number as the superseded "
                     "sicediff and is a different quantity. Yen's Eq. (71) "
                     "subtracts a brine term from exactly the bubbly ice "
                     "computed here, so the modelled sea ice conducts LESS "
                     "than the modelled glacial ice of the same density and "
                     "the two must not be deduplicated.",
        "what_this_is_worth_today": "nothing measured. dglac is identically "
                                    "zero on every cell of the climatology "
                                    "this world has, so the pair reaches no "
                                    "result until glaciermod grows ice.",
    }

    # GRAV-8. Priced only when a climatology is named, because the question is
    # what compaction is worth on a PARTICULAR run's snow and there is no
    # defensible default run to answer it against.
    compaction = (price_snow_compaction(args.climatology, gravity, year_seconds,
                                        rhosnow)
                  if args.climatology is not None
                  else {"priced": False,
                        "how": "pass --climatology; GRAV-8 asks for the term to "
                               "be priced offline against an existing "
                               "climatology before any source change"})

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/ice_properties.py",
        "standard": "IAPWS R10-06(2009), the equation of state 2006 for H2O ice Ih",
        "standard_reference": rel(
            ROOT / "references"
            / "iapws_2009_revised-release-on-the-equation-of-state-2006-for-h2o-ice-ih.pdf"),
        "check_tolerance": CHECK_TOLERANCE,
        "check_worst_relative": float(f"{worst:.3g}"),
        "check_states": [s["label"] for s in CHECK_TABLE],
        "validity": "IAPWS-06 is stated valid from 0 K to the melting curve, "
                    "and to 210 MPa. Every state evaluated here is inside it "
                    "by a wide margin.",
        "pure_ice_ih": {
            "melting_enthalpy_j_kg": round(fusion_j_kg, 2),
            "melting_enthalpy_state": "at the triple point, against IAPWS-95's "
                                      "liquid water on the same reference",
            "density_at_melting_kg_m3": round(gibbs(273.152519, P0)["rho"], 4),
            "specific_heat_at_melting_j_kg_k": round(gibbs(273.152519, P0)["cp"], 4),
            "profile_at_normal_pressure": profile,
        },
        "gravity": overburden,
        "pure_ice_conductivity": ice_conductivity,
        "glacial_ice": glacial_ice,
        "snow_conductivity_bracket": snow_bracket,
        "snow_compaction": compaction,
        "what_is_not_derivable_here": "the density, specific heat and "
                                      "conductivity of SEA ice are functions of "
                                      "brine volume, hence of the ice's own "
                                      "salinity and temperature, and this model "
                                      "carries neither as a variable. icemod.f90 "
                                      "already says so of the heat of fusion. "
                                      "The pure ice Ih values above are the "
                                      "bound those declarations sit against, "
                                      "not substitutes for them.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    print(f"IAPWS-06 reproduced Table 6 to {worst:.2e} relative, "
          f"inside {CHECK_TOLERANCE}\n")
    print(f"pure ice Ih, melting enthalpy      {fusion_j_kg:11.2f} J/kg")
    print(f"pure ice Ih, density at melting    {gibbs(273.152519, P0)['rho']:11.4f} kg/m3")
    print(f"pure ice Ih, cp at melting         {gibbs(273.152519, P0)['cp']:11.4f} J/kg/K\n")
    print(f"{'T (K)':>8}  {'rho (kg/m3)':>12}  {'cp (J/kg/K)':>12}")
    for row in profile:
        print(f"{row['temperature_k']:8.2f}  {row['density_kg_m3']:12.4f}  "
              f"{row['specific_heat_j_kg_k']:12.4f}")
    print(f"\npure ice conductivity, Yen (1981), W/m/K")
    print(f"{'T (K)':>8}  {'Eq. (33)':>9}  {'>195 K':>9}")
    for row in conductivity_profile:
        print(f"{row['temperature_k']:8.2f}  "
              f"{row['eq_33_all_temperatures_w_m_k']:9.4f}  "
              f"{row['table_3_above_195_k_w_m_k']:9.4f}")
    print(f"declared conductivity of the modelled sea ice: {ckapi:.4f}\n")

    print(f"modelled glacial ice at {rhoglac:.0f} kg/m3 and {tglacref:.2f} K")
    print(f"  bubbly factor, Yen Eq. (37)      {factor:11.6f} "
          f"(Eq. (70) gives {bubbly_ice_factor_full(rhoglac):.6f})")
    print(f"  heat capacity                    {glac_cap:11.4e} J/m3/K "
          f"against a superseded 2.07e6")
    print(f"  conductivity                     {glac_diff:11.6f} W/m/K "
          f"against a superseded 2.03")
    print(f"  over 233 K to melting it spans a factor "
          f"{max(diffs)/min(diffs):.3f} in conductivity and "
          f"{max(caps)/min(caps):.3f} in capacity\n")

    print("snow conductivity bracket, W/m/K")
    print(f"{'rho':>6}  {'slow':>7}  {'fast263':>8}  {'ratio':>6}  "
          f"{'Yen':>7}  {'Sturm':>7}")
    for row in bracket:
        print(f"{row['density_kg_m3']:6.0f}  {row['slow_kinetics_w_m_k']:7.4f}  "
              f"{row['fast_kinetics_263k_w_m_k']:8.4f}  "
              f"{row['fast_over_slow_263k']:6.3f}  "
              f"{row['yen_1981_apparent_w_m_k']:7.4f}  "
              f"{row['sturm_1997_as_lpj_guess_runs_it_w_m_k']:7.4f}")

    if compaction.get("priced") is not False:
        print("\nsnow compaction, PALADYN's Kojima self-loading, priced on "
              f"{compaction['snowy_land_bins']} snowy land bins of "
              f"{compaction['climatology']}")
        for name, body in compaction["cases"].items():
            print(f"  {name:26s} density p50 "
                  f"{body['density_kg_m3']['p50']:6.1f} kg/m3   depth p50 "
                  f"{body['physical_depth_m']['p50']:.4f} m   z/k p50 "
                  f"{body['conductive_resistance_m2_k_w']['p50']:.4f} m2K/W")
        for name, body in compaction["resistance_ratios"].items():
            print(f"  resistance {name:38s} p50 {body['p50']:.3f} "
                  f"(p5 {body['p5']:.3f} to p95 {body['p95']:.3f})")
        inst = compaction["instrument_versus_effect"]
        print(f"  the gravity term moves the resistance by "
              f"{100*inst['gravity_effect']:.1f} per cent; the kinetics "
              f"bracket at the same density moves it by "
              f"{100*inst['kinetics_bracket_effect']:.1f} per cent")
        print(f"  what it is worth on the modelled surface albedo: "
              f"{compaction['albedo']['worth']:.1f}, exactly")

    print(f"\ngravity: {overburden['relative_density_change']:.2e} relative "
          f"density change under {hmax:.0f} m of ice, against "
          f"{overburden['earth_relative_density_change']:.2e} at Earth's")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
