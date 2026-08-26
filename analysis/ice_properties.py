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

ROOT = PROJECT_ROOT
PLANET = ROOT / "config" / "planet.yaml"
ICEMOD = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src" / "icemod.f90"
LANDMOD = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src" / "landmod.f90"
LPJ_SOIL_H = ROOT / "vendor" / "lpj-guess" / "modules" / "soil.h"
OUTPUT = ROOT / "analysis" / "ice_properties.json"

# ---------------------------------------------------------------------------
# IAPWS R10-06(2009), the Gibbs energy of ice Ih. Coefficients from Table 2 of
# the release, which `references/iapws_2009_revised-release-on-the-equation-of-state-2006-for-h2o-ice-ih.pdf`
# carries and which Appendix I of `references/TEOS-10_Manual.pdf` restates.
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

FOURTEAU_2021_SNOW = {
    223.0: (2.564, -0.059, 0.0205),
    248.0: (2.172, 0.015, 0.0252),
    263.0: (1.985, 0.073, 0.0336),
    268.0: (1.883, 0.107, 0.0386),
    273.0: (1.776, 0.147, 0.0455),
}
"""Eq. (18), `k = a vf**2 + b vf + c` in the ice volume fraction, per temperature."""

RHOICE_F2021 = 917.0
"""The normalising ice density Fourteau's volume fraction is taken against.

The same constant `landmod.f90` declares as `RHOICE_F2021`, and it belongs to
the RELATION rather than to this world's ice: changing it would not describe
denser ice, it would misread the fit.
"""


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


def lpj_constant(name: str) -> float:
    """One LPJ-GUESS constant, read the same way and for the same reason."""
    pattern = re.compile(rf"^\s*const\s+double\s+{name}\s*=\s*([-\d.eE+]+)", re.M)
    match = pattern.search(LPJ_SOIL_H.read_text())
    if match is None:
        raise SystemExit(f"{rel(LPJ_SOIL_H)} no longer declares `{name}`.")
    return float(match.group(1))


def sturm_1997_snow(density_kg_m3: float) -> float:
    """Sturm et al. (1997)'s needle-probe regression, W/m/K, as LPJ-GUESS runs it.

    Carried NOT as a candidate -- `notes/audits/cryosphere-material-properties.md`
    argues against adopting it -- but because `vendor/lpj-guess`'s `soil.cpp`
    computes its snow conductivity from exactly this, so the ecology column and
    the climate column stand on two different relations for one material
    property and the gap between them is a number rather than a suspicion.
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


def snow_conductivity_fast(density_kg_m3: float, temperature_k: float = 263.0) -> float:
    """Fourteau (2021) Eq. (18): the fast-kinetics arm, W/m/K.

    This is the relation `landmod`'s `landini` evaluates to set `snowdiff`, and
    it is reproduced here rather than restated so that the two cannot drift.
    """
    a, b, c = FOURTEAU_2021_SNOW[temperature_k]
    vf = density_kg_m3 / RHOICE_F2021
    return float(a * vf * vf + b * vf + c)


def snow_conductivity_slow(density_kg_m3: float) -> float:
    """Calonne (2011) Eq. (12): the slow-kinetics arm, W/m/K.

    Conduction through ice and interstitial air, with no latent heat term. One
    temperature only: the fit is his whole sample set at 271 K.
    """
    a, b, c = CALONNE_2011_SNOW
    return float(a * density_kg_m3 * density_kg_m3 + b * density_kg_m3 + c)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()

    worst, failures = verify_against_release()
    if failures:
        raise SystemExit(
            "this implementation of IAPWS-06 does not reproduce the release's "
            "own check table, so nothing computed from it means anything:\n"
            + "\n".join(failures))

    config = yaml.safe_load(PLANET.read_text())
    gravity = float(config["planet"]["gravity_m_s2"])

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
        "adopted": "the FAST arm. `landmod`'s landini evaluates Eq. (18) at 263 "
                   "K, so the adopted relation is the bracket's upper endpoint "
                   "rather than a point inside it, and the slow arm below is "
                   "how far under it the other limit sits.",
        "declared_density_kg_m3": rhosnow,
        "declared_conductivity_w_m_k": snowdiff,
        "slow_arm_residual_sd_w_m_k": CALONNE_2011_RESIDUAL_SD,
        "per_density": bracket,
        "instrument_versus_effect": "the bracket is wider than the scatter of "
                                    "the fit that defines its lower arm at "
                                    "every density carried, so it is a real "
                                    "disagreement and not the noise of one "
                                    "regression.",
        "cross_component": "vendor/lpj-guess/modules/soil.cpp computes its snow "
                           "conductivity from Sturm et al. (1997), which is "
                           "BELOW the slow arm at every density here. So the "
                           "ecology column and the climate column insulate "
                           "their soil differently from the same snowfall at "
                           "the same density, and the difference is larger "
                           "than the bracket between the two published limits.",
    }

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
        "snow_conductivity_bracket": snow_bracket,
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

    print("snow conductivity bracket, W/m/K")
    print(f"{'rho':>6}  {'slow':>7}  {'fast263':>8}  {'ratio':>6}  "
          f"{'Yen':>7}  {'Sturm':>7}")
    for row in bracket:
        print(f"{row['density_kg_m3']:6.0f}  {row['slow_kinetics_w_m_k']:7.4f}  "
              f"{row['fast_kinetics_263k_w_m_k']:8.4f}  "
              f"{row['fast_over_slow_263k']:6.3f}  "
              f"{row['yen_1981_apparent_w_m_k']:7.4f}  "
              f"{row['sturm_1997_as_lpj_guess_runs_it_w_m_k']:7.4f}")

    print(f"\ngravity: {overburden['relative_density_change']:.2e} relative "
          f"density change under {hmax:.0f} m of ice, against "
          f"{overburden['earth_relative_density_change']:.2e} at Earth's")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
