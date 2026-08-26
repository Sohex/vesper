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
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.paths import PROJECT_ROOT, rel  # noqa: E402

ROOT = PROJECT_ROOT
PLANET = ROOT / "config" / "planet.yaml"
ICEMOD = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src" / "icemod.f90"
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
    print(f"\ngravity: {overburden['relative_density_change']:.2e} relative "
          f"density change under {hmax:.0f} m of ice, against "
          f"{overburden['earth_relative_density_change']:.2e} at Earth's")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
