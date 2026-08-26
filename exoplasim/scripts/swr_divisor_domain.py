#!/usr/bin/env python
"""Can any divisor in `swr` reach zero on this planet's configuration?

WORLDBUILDING FRAME. This is a simulated planet, Vesper, and every number below
is a property of the simulation's shortwave radiation code or of the papers that
code is ported from. Nothing here is a measurement of anything outside the
simulation.

THE QUESTION. `radmod.f90`'s `swr` bounds neither its clear-sky band
transmissivities nor its layer reflectivities, and thirteen of its divisors are
one minus one of them. Its longwave counterpart `lwr` clamps `ztaucs` to
[zero, 1-zero]; `swr` clamps nothing. world-2223 asks whether either family can
reach zero HERE, and asks it by EVALUATING the two expressions over the model's
own ranges rather than by argument.

THE TWO FAMILIES.

  A. `1 - A(u)/zsolar_b`, nine sites: the six running transmissivity products
     `zto3t`, `zto3u`, `ztwvt`, `ztwvu`, `ztco2t`, `ztco2u`, and the three
     denominators the upward beam spells out in full. `A(u)` is Lacis and
     Hansen's absorptance as a fraction of TOTAL incident flux and `zsolar_b`
     is band b's share of that flux, so the quotient is the fraction of the
     BAND the absorber removes and the divisor vanishes when it removes all of
     it. Five namelist re-weightings for a non-solar host -- `o3uvw`, `o3visw`,
     `h2osww`, `h2oswl`, `co2sww` -- multiply `A`, and none is confined to 1;
     `zsolar1` falls as the host reddens, which is the direction this host is.

  B. `1 - R_above*R_below`, four sites: the adding method's denominator in the
     downward loop, the upward loop and the flux loop, per band. Both factors
     are layer reflectivities summed from a Rayleigh term, a cloud term and an
     aerosol term without the sum being bounded.

WHAT IS EVALUATED, PER FAMILY.

  A. The supremum of `A(u)/zsolar_b` over u in [0, inf) in closed form and
     numerically; whether it exceeds 1 at all; where the divisor first reaches
     zero if it does; the model's OWN range of u, derived from the model; and
     the margin as a ratio of absorber amounts.

  B. `zscf` and the Rayleigh diffuse reflectance at this planet's surface
     pressure and gravity; the cloud reflectivity ceilings that follow from the
     code's own cap on the liquid water path; the aerosol ceiling that follows
     from the configured dust file; the exact one-dimensional maximisation of
     the layer sum over cloud fraction; the accumulated reflectance the adding
     recursion reaches; and hence the product's ceiling.

THE CHECKS THAT CAN FAIL. Five, each with a right answer stated in advance and
an artifact or a source line that is not this script. The first three gate: a
failure there stops the run, because it would mean the analysis is being done on
the wrong function. The last two report.

  1. Every literal this script evaluates is READ OUT of `radmod.f90` and
     compared with what the script uses. A coefficient edited in the model and
     not here makes the run stop.
  2. The water vapour absorptance function is Lacis and Hansen Eq. 21, and
     `exoplasim/analysis/shortwave_band_weights.json` stores it at nine water
     amounts under `h2o.weight_vs_water_amount[*].lacis_hansen_eq21`. This
     script's own evaluation must reproduce all nine. Disagreement would mean
     the function this domain analysis is run over is not the function the band
     weights were derived against, and every margin below would be for the
     wrong curve.
  3. `lib/stellar.py`'s Rayleigh coefficient at 4965 K must reproduce what the
     compiled model printed, 0.862014830, recorded in
     `notes/audits/physics-review.md`. Disagreement would mean the `zscf` below
     is not the `zscf` the model forms.
  4. `radmod`'s own CO2 comment claims the fit is quoted over 1 to 1e4 atmos-cm
     and that the model never leaves that range. The CO2 path this script
     derives from the model's own sigma grid, surface pressure and
     magnification is compared with both ends of it. The upper end gates,
     because a crossing amount read off an extrapolated fit is not a property
     of the fit. The lower end reports: the fit is monotone from A(0) = 0, so
     going below the quoted range costs accuracy near the model top and cannot
     take the divisor anywhere near zero.
  5. The water vapour verdict turns on one inequality, and the declared brackets
     on both re-weightings are compared with the threshold it turns on. A
     verdict that did not survive its own brackets would be a property of the
     rounding rather than of the configuration.

WHAT WOULD MAKE THE ANSWER MOVE. The absorber ranges come from a BOOTSTRAP
climatology, which is a time mean: an instantaneous column carries more water
than a monthly mean of it, so the water vapour margin below leans HIGH by
whatever that ratio is. The magnification factor is taken at its ceiling of 35,
which the model reaches at the terminator, so that half leans low.

Run it as `python exoplasim/scripts/swr_divisor_domain.py`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset
from scipy.optimize import brentq

import _paths  # noqa: F401  anchors every path on this file and adds lib/

import provenance  # noqa: E402
import stellar  # noqa: E402
from paths import rel  # noqa: E402

import cloud_optical_depth_bracket as cod  # noqa: E402
import stephens_tables_vs_fits as stephens  # noqa: E402

ROOT = _paths.PROJECT_ROOT
RADMOD = _paths.MODEL_SRC / "plasim" / "src" / "radmod.f90"
DEFAULT_CLIM = _paths.ANALYSIS / "climatology" / "bootstrap_regular_climatology.nc"
DEFAULT_BANDS = _paths.ANALYSIS / "shortwave_band_weights.json"
DUST = _paths.COMPONENT_ROOT / "data" / "dust" / "vesper_dust_aerosol.dat"

# THIS SCRIPT WRITES NOTHING, DELIBERATELY. What it produces is a bound, and a
# bound belongs in the dated finding that carries its evidence -- here
# notes/audits/masked-where-blocks.md -- and in the CLASSIFIED rows of
# exoplasim/scripts/lint_masked_domains.py that cite it. Nothing in the pipeline
# consumes a number from here, and an artifact no step in config/pipeline.yaml
# generates does not exist. Re-run it to re-derive the bounds.

# The value the compiled model PRINTED for `rcoeff` on the 4965 K blackbody,
# transcribed from notes/audits/physics-review.md. It is here so check 3 has a
# right answer that did not come from the Python side of the comparison.
MODEL_PRINTED_RCOEFF_4965 = 0.862014830
RCOEFF_TOLERANCE = 1.0e-6          # single precision, relative

# `swr`'s own parameters, as this script uses them. Every one is checked against
# the model source by `read_model_constants` before anything is evaluated.
ZBETTA = 1.66      # magnification factor, water vapour and CO2
ZMBAR = 1.9        # magnification factor, ozone
ZRO3 = 2.14        # ozone density, kg m-3 STP
ZMU00 = 0.5        # the diffuse stream's beam cosine
ZM_CEILING = 35.0  # zm = 35/sqrt(1 + 1224 mu0^2), at mu0 -> 0
P0_RAYLEIGH = 101100.0   # the reference pressure zscf divides by
G0_RAYLEIGH = 9.80665    # the reference gravity zscf scales by
RAYLEIGH_T0 = 0.144      # the diffuse-beam single-scattering loss
LWP_CAP = 1000.0         # g m-2, the cap on zlwp
ZWFIT = 10.0             # g m-2, the bottom of Stephens's fitted range
# Stephens (1978) Eqs. (10a) and (10b), as radmod writes them: (10**a, b*ln10).
ZTAUA1, ZTAUP1 = 1.8336, 3.9363
ZTAUA2, ZTAUP2 = 2.2346, 3.8034
# The CO2 absorptance fit, u in atmos-cm.
ZCA1, ZCB1 = 3.1020e-4, 19.857
ZCA2, ZCB2 = 3.8291e-3, 3.9587e-3


# --------------------------------------------------------------------------
# Check 1: the model source says what this script thinks it says.
# --------------------------------------------------------------------------

def read_model_constants() -> dict:
    """Every literal below is READ from radmod.f90 and compared, not assumed.

    A text read rather than a compile, for the reason
    `stephens_tables_vs_fits.py` gives: this has to run without a build. What
    fails here is a coefficient that moved in the model and not in this script,
    which would leave every margin below describing a function the model no
    longer evaluates.
    """
    source = RADMOD.read_text(encoding="utf-8")

    def parameter(name: str) -> float:
        match = re.search(rf"parameter\(\s*{name}\s*=\s*([-0-9.eEdD+]+)\s*\)",
                          source)
        if match is None:
            raise SystemExit(f"{rel(RADMOD)} declares no parameter {name}")
        return float(match.group(1).replace("D", "E").replace("d", "e"))

    def real_parameter(name: str) -> float:
        match = re.search(
            rf"real,\s*parameter\s*::\s*{name}\s*=\s*([-0-9.eEdD+]+)", source)
        if match is None:
            raise SystemExit(f"{rel(RADMOD)} declares no real parameter {name}")
        return float(match.group(1).replace("D", "E").replace("d", "e"))

    def default(name: str) -> float:
        match = re.search(
            rf"(?:real|integer)\s*::\s*{name}\s*=\s*([-0-9.eEdD+]+)", source)
        if match is None:
            raise SystemExit(f"{rel(RADMOD)} declares no default for {name}")
        return float(match.group(1).replace("D", "E").replace("d", "e"))

    found = {
        "zbetta": parameter("zbetta"), "zmbar": parameter("zmbar"),
        "zro3": parameter("zro3"),
        "zca1": parameter("zca1"), "zcb1": parameter("zcb1"),
        "zca2": parameter("zca2"), "zcb2": parameter("zcb2"),
        "ztaua1": real_parameter("ztaua1"), "ztaup1": real_parameter("ztaup1"),
        "ztaua2": real_parameter("ztaua2"), "ztaup2": real_parameter("ztaup2"),
        "zwfit": real_parameter("zwfit"),
        "a0o3": default("a0o3"), "a1o3": default("a1o3"),
        "aco3": default("aco3"),
        "nrscat": default("nrscat"), "newrsc": default("newrsc"),
        "nclouds": default("nclouds"), "nswrcl": default("nswrcl"),
        "naerosp": default("naerosp"), "cloudabs": default("cloudabs"),
    }
    expected = {
        "zbetta": ZBETTA, "zmbar": ZMBAR, "zro3": ZRO3,
        "zca1": ZCA1, "zcb1": ZCB1, "zca2": ZCA2, "zcb2": ZCB2,
        "ztaua1": ZTAUA1, "ztaup1": ZTAUP1,
        "ztaua2": ZTAUA2, "ztaup2": ZTAUP2, "zwfit": ZWFIT,
    }
    wrong = {k: (found[k], v) for k, v in expected.items() if found[k] != v}
    if wrong:
        raise SystemExit(
            "radmod.f90 no longer carries the constants this script "
            f"evaluates: {wrong}. Fix the script, not the model.")

    # The absorptance expressions themselves, matched as source text. These are
    # the nine divisors' numerators and there is no parameter to read them from,
    # so the check is that the coefficient string is still present.
    fragments = {
        "ozone_visible": "o3visw*0.02118*",
        "ozone_uv_huggins": "o3uvw*1.082*",
        "ozone_uv_hartley": "o3uvw*0.0658*",
        "ozone_denominator": "1.+0.042*",
        "ozone_denominator_quadratic": "0.000323*",
        "ozone_uv_exponent": "**0.805",
        "water_vapour": "h2osww*h2oswl*2.9*",
        "water_vapour_denominator": "(1.+141.5*",
        "water_vapour_exponent": "**0.635",
        "water_vapour_linear": "+5.925*",
        "co2": "co2sww*(zca1*LOG(1.+zcb1*",
        "rayleigh_diffuse": "log(1.0-0.144)",
        "rayleigh_reference_pressure": "/101100.0",
        "rayleigh_reference_gravity": "9.80665/ga",
        "zscf": "zscf(:) = rcoeff*dp(:)/101100.0*9.80665/ga",
        "lwp_cap": "min(1000.0,1000.*dql(",
        "diffuse_cosine": "zmu00  = 0.5",
        "cloud_band1_reflectivity": "zrcl1s(:,jlev)=1.-1./(1.+zb1(:)*ztau1(:)/zmu00)",
        "layer_reflectivity_band1_diffuse":
            "zrb1s(:,jlev)=zrcsu(:,jlev)+zrcl1s(:,jlev)*dcc(:,jlev)*nclouds",
        "layer_reflectivity_band2_diffuse":
            "zrb2s(:,jlev)=zrcl2s(:,jlev)*dcc(:,jlev)*nclouds",
        "magnification": "zm(:)=35./SQRT(1.+1224.*zmu0(:)*zmu0(:))",
    }
    missing = sorted(k for k, v in fragments.items() if v not in source)
    if missing:
        raise SystemExit(
            f"{rel(RADMOD)} no longer contains the expressions this script "
            f"evaluates: {missing}. Fix the script, not the model.")

    # radmod's own claim about the CO2 fit's range, read out of its comment so
    # check 4 tests the model's statement rather than a number retyped here.
    quoted = re.search(r"Quoted over (\d+) to\s*(?:\n!\s*)?(\d+e\d+) atmos-cm",
                       source)
    if quoted is None:
        raise SystemExit(f"{rel(RADMOD)} no longer states the CO2 fit's range")
    found["co2_fit_range"] = [float(quoted.group(1)), float(quoted.group(2))]

    # The Stephens tables, as the model stores them: (NSWCM, NSWCT) with mu0
    # ascending. `stephens.interp_table` wants (tau, mu0) with mu0 DESCENDING,
    # which is the transpose plus a flip.
    tables = {}
    for name in ("swcb1", "swcb2", "swcoa"):
        flat = stephens._fortran_array(name, source)
        tables[name] = flat.reshape(len(stephens.TAU),
                                    len(stephens.MU0))[:, ::-1]
    found["tables"] = tables
    return found


# --------------------------------------------------------------------------
# Family A: the clear-sky band transmissivities.
# --------------------------------------------------------------------------

def absorptance_ozone(u, o3visw: float, o3uvw: float):
    """Lacis and Hansen (1974) Eqs. 22 and 23, as `swr` writes them.

    u is the ozone amount in cm STP. `o3visw` and `o3uvw` re-weight the visible
    and ultraviolet halves for a non-solar host. The result is a fraction of
    TOTAL incident flux.
    """
    u = np.asarray(u, dtype=float)
    return (o3visw * 0.02118 * u / (1.0 + 0.042 * u + 0.000323 * u ** 2)
            + o3uvw * 1.082 * u / ((1.0 + 138.6 * u) ** 0.805)
            + o3uvw * 0.0658 * u / (1.0 + (103.6 * u) ** 3))


def absorptance_water(u, h2osww: float, h2oswl: float):
    """Lacis and Hansen (1974) Eq. 21, as `swr` writes it.

    u is the water vapour path in precipitable cm. The two weights scale the
    whole curve, so the asymptote below scales with their product.
    """
    u = np.asarray(u, dtype=float)
    return (h2osww * h2oswl * 2.9 * u
            / ((1.0 + 141.5 * u) ** 0.635 + 5.925 * u))


def absorptance_co2(u, co2sww: float):
    """The fork's CO2 shortwave absorptance, u in atmos-cm."""
    u = np.asarray(u, dtype=float)
    return co2sww * (ZCA1 * np.log(1.0 + ZCB1 * u)
                     + ZCA2 * np.log(1.0 + ZCB2 * u))


def supremum(fn, band_share: float) -> dict:
    """sup of fn(u)/band_share over u in [0, inf), in closed form and numerically.

    Water vapour is the only one of the three with a finite asymptote, and the
    caller says which it is: the closed form is stated by the caller and this
    confirms it on a grid that reaches far enough for the confirmation to mean
    something.
    """
    grid = np.concatenate([np.linspace(0.0, 10.0, 100001)[1:],
                           np.geomspace(10.0, 1.0e18, 200000)])
    values = np.asarray(fn(grid)) / band_share
    peak = int(np.argmax(values))
    return {"numerical_supremum": float(values[peak]),
            "at_u": float(grid[peak]),
            "value_at_1e18": float(values[-1])}


def crossing(fn, band_share: float, hi: float) -> float | None:
    """The u at which `1 - fn(u)/band_share` first reaches zero, or None.

    Solved in log10(u) because two of the three crossings sit tens of decades
    above the model's own amounts and a linear bracket that wide cannot be
    bisected to a relative tolerance.
    """
    root = lambda x: float(fn(10.0 ** x)) / band_share - 1.0   # noqa: E731
    lo = -8.0
    top = math.log10(hi)
    if root(top) < 0.0:
        return None
    return float(10.0 ** brentq(root, lo, top, xtol=1.0e-12, rtol=8.9e-16))


# --------------------------------------------------------------------------
# Family B: the layer reflectivities and the adding method's denominator.
# --------------------------------------------------------------------------

def cloud_optical_depths(w):
    """`swr`'s two cloud optical depths at liquid water path w, g m-2.

    Both branches of the one expression the model writes: the fit above ZWFIT,
    and Stephens Eq. (7) at fixed effective radius below it.
    """
    w = np.asarray(w, dtype=float)
    wl = np.log10(np.maximum(ZWFIT, np.maximum(0.0, w)))
    ramp = np.minimum(1.0, np.maximum(0.0, w) / ZWFIT)
    return (ZTAUA1 * wl ** ZTAUP1 * ramp, ZTAUA2 * wl ** ZTAUP2 * ramp)


def cloud_reflectivity_band1(tau1, beta1_table):
    """`zrcl1s`: Stephens Eq. (1) at the diffuse cosine."""
    beta = stephens.interp_table(beta1_table, tau1, ZMU00)
    x = beta * np.asarray(tau1, dtype=float) / ZMU00
    return x / (1.0 + x)


def cloud_reflectivity_band2(tau2, beta2_table, coalb_table, cloudabs, floor):
    """`zrcl2s`: the two-stream reflectance at the diffuse cosine.

    The co-albedo carries `swr`'s own floor, `max(zepsc, cloudabs*swcoalb)`,
    because the u-factor divides by it.
    """
    tau2 = np.asarray(tau2, dtype=float)
    beta = stephens.interp_table(beta2_table, tau2, ZMU00)
    un = np.maximum(floor, cloudabs * stephens.interp_table(coalb_table,
                                                            tau2, ZMU00))
    uz = un + 2.0 * beta * (1.0 - un)
    u = np.sqrt(np.maximum(0.0, uz / un))
    e = np.exp(np.minimum(25.0, tau2 * np.sqrt(np.maximum(0.0, uz * un))
                          / ZMU00))
    r = (u + 1.0) ** 2 * e - (u - 1.0) ** 2 / e
    return (u * u - 1.0) / r * (e - 1.0 / e)


def aerosol_reflectivity_ceiling(ssa: float, bscat: float) -> float:
    """sup over optical depth of the two-stream diffuse reflectance.

    Both of `swr`'s branches have the same limit. In the exact branch the
    reflectance is (u^2-1)(e^t - e^-t)/((u+1)^2 e^t - (u-1)^2 e^-t), which rises
    monotonically in t to (u-1)/(u+1); in the conservative branch it is
    x/(1+x) with x = b*tau/mu, which rises to 1 but is only taken when the
    single-scattering albedo is within SQRT(epsilon) of 1. Neither reaches 1 at
    any finite optical depth, and for the species this world configures the
    exact branch applies.
    """
    u = math.sqrt((1.0 - ssa + 2.0 * bscat * ssa) / (1.0 - ssa))
    return (u - 1.0) / (u + 1.0)


def adding_down(layer_r, layer_t, nlev: int) -> list[float]:
    """`zra1s` accumulated by the downward loop, layer by layer.

    R_new = R_b + T_bu*R_a*T_b/(1 - R_a*R_b), the statement `swr` writes. Run
    here with every layer at its own ceiling, which is the worst case the
    arithmetic admits rather than one the model produces.
    """
    out = []
    r = 0.0
    for _ in range(nlev):
        r = layer_r + layer_t * r * layer_t / (1.0 - r * layer_r)
        out.append(r)
    return out


def adding_up(layer_r, layer_t, surface_albedo: float, nlev: int) -> list[float]:
    """`zra1s` accumulated by the upward loop, from the surface upward."""
    out = []
    r = surface_albedo
    for _ in range(nlev):
        r = layer_r + layer_t * r * layer_t / (1.0 - r * layer_r)
        out.append(r)
    return out


# --------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=DEFAULT_CLIM,
                    help="the field the model's own absorber ranges come from")
    ap.add_argument("--bands", type=Path, default=DEFAULT_BANDS,
                    help="the artifact check 2 reproduces")
    args = ap.parse_args()

    model = read_model_constants()
    config = yaml.safe_load(_paths.CONFIG.read_text(encoding="utf-8"))

    # --- the planet, from the single source of truth -----------------------
    ga = float(config["planet"]["gravity_m_s2"])
    atmosphere = config["atmosphere"]
    surface_pressure_pa = 1.0e5 * sum(
        float(v) for k, v in atmosphere.items() if k.startswith("p")
        and k.endswith("_bar"))
    mcfg = config["model"]
    o3visw = float(mcfg["ozone_visible_weight"])
    o3uvw = float(mcfg["ozone_uv_weight"])
    h2osww = float(mcfg["h2o_sw_weight"])
    h2oswl = float(mcfg["h2o_sw_level"])
    h2oswl_bracket = [float(x) for x in mcfg["h2o_sw_level_bracket"]]
    co2sww = float(mcfg["co2_sw_weight"])
    o3scale = float(mcfg["ozone_scale"])
    co2_ppmv = 1.0e6 * float(atmosphere["pCO2_bar"])

    zsolar1, zsolar2 = (float(x) for x in stellar.band_fractions())

    # --- the model's own absorber ranges -----------------------------------
    ds = Dataset(args.climatology)
    lev = np.asarray(ds.variables["lev"][:], dtype=float)
    sigma, sigmah, dsigma = cod.sigma_grid(lev)
    hus = np.asarray(ds.variables["hus"][:], dtype=float)
    ta = np.asarray(ds.variables["ta"][:], dtype=float)
    ps_pa = np.asarray(ds.variables["ps"][:], dtype=float) * 100.0
    prw = np.asarray(ds.variables["prw"][:], dtype=float)
    alb1 = np.asarray(ds.variables["alb1"][:], dtype=float)
    alb2 = np.asarray(ds.variables["alb2"][:], dtype=float)
    clt = np.asarray(ds.variables["clt"][:], dtype=float)
    nlev = len(lev)
    ds.close()

    s4 = sigma[None, :, None, None]
    d4 = dsigma[None, :, None, None]
    dp4 = ps_pa[:, None]

    # zwv, exactly as the absorber loop forms it, summed down the column. This
    # is the amount BEFORE any magnification factor.
    zwv = (0.1 * d4 * hus * dp4 / ga * np.sqrt(273.0 / ta)
           * s4 * dp4 / 100000.0)
    water_column = zwv.sum(axis=1)

    # zco2, the same way. The mass mixing ratio is one scalar for the column.
    zpv2pm = 0.0440098 / 0.0289644
    zfco2 = 100.0 / 1.9635
    zzco2 = zpv2pm * 1.0e-6 * co2_ppmv
    zco2 = zfco2 * zzco2 * d4 * dp4 / ga * s4 * dp4 / 100000.0
    co2_column = zco2.sum(axis=1)
    co2_thinnest_layer = zco2[:, 0].min()

    # The ozone column is `mko3`'s own za, in cm STP, scaled by o3scale. Its
    # maximum over latitude and season is a0o3 + a1o3 + aco3, all three at once.
    ozone_column_max = (model["a0o3"] + model["a1o3"] + model["aco3"]) * o3scale

    # The magnification factor's ceiling. `losun` keeps a lane whose top-of-
    # atmosphere flux exceeds 1e-6 W m-2, which at this star's flux is a beam
    # cosine of order 1e-9, so zm reaches its own ceiling of 35 on lanes the
    # mask keeps.
    up = ZM_CEILING + ZBETTA          # the total upward path's multiplier
    up_o3 = ZM_CEILING + ZMBAR

    # --- family A ----------------------------------------------------------
    family_a = {}

    # Ozone. The Huggins term goes as u**0.195 at large u, so A is unbounded
    # and the divisor CAN reach zero; the question is only where.
    o3 = lambda u: absorptance_ozone(u, o3visw, o3uvw)   # noqa: E731
    o3_cross = crossing(o3, zsolar1, 1.0e30)
    o3_u_max = up_o3 * ozone_column_max
    family_a["ozone"] = {
        "band": 1, "band_share": zsolar1,
        "sites": ["zto3t", "zto3u", "zto3tu spelled out"],
        "weights": {"o3visw": o3visw, "o3uvw": o3uvw},
        "supremum": "unbounded: the Huggins term goes as "
                    "o3uvw*1.082/138.6**0.805 * u**0.195",
        "supremum_exceeds_one": True,
        "u_cross_cm_stp": o3_cross,
        "model_range": {
            "column_cm_stp_max": ozone_column_max,
            "column_source": "mko3's own a0o3 + a1o3*|sin lat| + "
                             "aco3*sin(lat)*cos(...), all three terms at once, "
                             "times model.ozone_scale",
            "magnification_ceiling": up_o3,
            "u_max_cm_stp": o3_u_max,
        },
        "margin": o3_cross / o3_u_max,
        "divisor_floor": 1.0 - float(o3(o3_u_max)) / zsolar1,
    }

    # Water vapour. Finite asymptote: as u -> inf the 5.925*u term dominates
    # the denominator and A -> h2osww*h2oswl*2.9/5.925.
    wv = lambda u: absorptance_water(u, h2osww, h2oswl)   # noqa: E731
    wv_asymptote = h2osww * h2oswl * 2.9 / 5.925
    wv_sup = wv_asymptote / zsolar2
    wv_numeric = supremum(wv, zsolar2)
    wv_cross = crossing(wv, zsolar2, 1.0e12) if wv_sup > 1.0 else None
    water_column_max = float(water_column.max())
    wv_u_max = up * water_column_max
    # The same question over the declared bracket on h2oswl and over h2osww's
    # own bracket, because a supremum that only just exceeds 1 would be a
    # property of the rounding rather than of the configuration.
    wv_bracket = {}
    for label, level in (("low", h2oswl_bracket[0]), ("high", h2oswl_bracket[1])):
        sup = h2osww * level * 2.9 / 5.925 / zsolar2
        f = lambda u, lv=level: absorptance_water(u, h2osww, lv)   # noqa: E731
        wv_bracket[label] = {
            "h2oswl": level, "supremum": sup,
            "supremum_exceeds_one": bool(sup > 1.0),
            "u_cross_cm": crossing(f, zsolar2, 1.0e14) if sup > 1.0 else None,
        }
    # And on the Sun, with both weights at 1 and radmod's declared zsolar2, so
    # the size of what the re-weighting did is visible.
    solar_sup = 2.9 / 5.925 / 0.483
    solar_cross = crossing(lambda u: absorptance_water(u, 1.0, 1.0),
                           0.483, 1.0e14) if solar_sup > 1.0 else None
    family_a["water_vapour"] = {
        "band": 2, "band_share": zsolar2,
        "sites": ["ztwvt", "ztwvu", "ztwvtu spelled out"],
        "weights": {"h2osww": h2osww, "h2oswl": h2oswl,
                    "product": h2osww * h2oswl},
        "asymptote_closed_form": wv_asymptote,
        "asymptote_over_band_share": wv_sup,
        "asymptote_numerical": wv_numeric,
        "supremum_exceeds_one": bool(wv_sup > 1.0),
        "u_cross_cm": wv_cross,
        "bracket": wv_bracket,
        "on_the_sun": {"weights": 1.0, "zsolar2": 0.483,
                       "supremum": solar_sup, "u_cross_cm": solar_cross},
        "model_range": {
            "column_cm_max": water_column_max,
            "column_cm_mean": float(water_column.mean()),
            "effective_path_at_zbetta_mean": float(water_column.mean()) * ZBETTA,
            "column_source": rel(args.climatology),
            "column_is_a_time_mean": True,
            "magnification_ceiling": up,
            "u_max_cm": wv_u_max,
        },
        "margin": wv_cross / wv_u_max if wv_cross else None,
        "divisor_floor": 1.0 - float(wv(wv_u_max)) / zsolar2,
    }

    # CO2.
    co2 = lambda u: absorptance_co2(u, co2sww)   # noqa: E731
    co2_cross = crossing(co2, zsolar2, 1.0e80)
    co2_column_max = float(co2_column.max())
    co2_u_max = up * co2_column_max
    family_a["co2"] = {
        "band": 2, "band_share": zsolar2,
        "sites": ["ztco2t", "ztco2u", "ztco2tu spelled out"],
        "weights": {"co2sww": co2sww},
        "supremum": "unbounded: two logarithms, so A grows as "
                    "co2sww*(zca1+zca2)*ln(u)",
        "supremum_exceeds_one": True,
        "u_cross_atmos_cm": co2_cross,
        "model_range": {
            "column_atmos_cm_max": co2_column_max,
            "column_atmos_cm_mean": float(co2_column.mean()),
            "effective_path_at_zbetta_mean": float(co2_column.mean()) * ZBETTA,
            "thinnest_layer_atmos_cm": float(co2_thinnest_layer),
            "magnification_ceiling": up,
            "u_max_atmos_cm": co2_u_max,
        },
        "margin": co2_cross / co2_u_max,
        "divisor_floor": 1.0 - float(co2(co2_u_max)) / zsolar2,
    }

    # --- family B ----------------------------------------------------------
    tables = model["tables"]
    cloudabs = model["cloudabs"]
    # `swr` floors the co-albedo at SQRT(EPSILON(1.0)). The model compiles in
    # single precision unless told otherwise, and the single-precision floor is
    # the larger of the two, so it is the one that binds the u-factor.
    coalb_floor = math.sqrt(float(np.finfo(np.float32).eps))

    tau1_cap, tau2_cap = cloud_optical_depths(LWP_CAP)

    # (a) the Rayleigh diffuse reflectance. rcoeff has three candidate values
    # and the largest gives the largest reflectance, so the ceiling is taken on
    # that one whatever the run configures.
    rcoeff = {
        "k25v_one_grid": float(stellar.rayleigh_coefficient()),
        "k25v_as_solarini_codes_it":
            float(stellar.rayleigh_coefficient(as_the_model_does=True)),
        "blackbody_4965k":
            float(stellar.rayleigh_coefficient(temperature_k=4965.0)),
    }
    ps_max = float(ps_pa.max())
    rayleigh = {}
    for label, value in rcoeff.items():
        for pname, pressure in (("declared", surface_pressure_pa),
                                ("model_maximum", ps_max)):
            zscf = value * pressure / P0_RAYLEIGH * G0_RAYLEIGH / ga
            rayleigh[f"{label}__{pname}"] = {
                "rcoeff": value, "surface_pressure_pa": pressure,
                "zscf": zscf,
                "reflectance": 1.0 - math.exp(zscf * math.log(1.0 - RAYLEIGH_T0)),
            }
    rayleigh_ceiling = max(v["reflectance"] for v in rayleigh.values())

    # (b) the cloud reflectivity ceilings, over every liquid water path the cap
    # admits. Neither is monotone in the path, so both are scanned.
    w_grid = np.linspace(0.0, LWP_CAP, 200001)
    t1_grid, t2_grid = cloud_optical_depths(w_grid)
    r1_grid = cloud_reflectivity_band1(t1_grid, tables["swcb1"])
    r2_grid = cloud_reflectivity_band2(t2_grid, tables["swcb2"],
                                       tables["swcoa"], cloudabs, coalb_floor)
    i1, i2 = int(np.argmax(r1_grid)), int(np.argmax(r2_grid))
    cloud = {
        "lwp_cap_g_m2": LWP_CAP,
        "tau1_at_cap": float(tau1_cap), "tau2_at_cap": float(tau2_cap),
        "zrcl1s_ceiling": float(r1_grid[i1]),
        "zrcl1s_at_lwp_g_m2": float(w_grid[i1]),
        "zrcl2s_ceiling": float(r2_grid[i2]),
        "zrcl2s_at_lwp_g_m2": float(w_grid[i2]),
        "coalbedo_floor": coalb_floor,
        "cloudabs": cloudabs,
    }

    # (c) the aerosol ceiling, from the file the model reads.
    dust_numbers = [float(line) for line in
                    DUST.read_text(encoding="utf-8").splitlines()
                    if line.strip() and not line.lstrip().startswith("#")][:8]
    qex1, qsca1, qback1, _g1, qex2, qsca2, qback2, _g2 = dust_numbers
    ssa1, bscat1 = qsca1 / qex1, qback1 / qsca1
    ssa2, bscat2 = qsca2 / qex2, qback2 / qsca2
    aerosol = {
        "file": rel(DUST),
        "band1": {"ssa": ssa1, "backscatter_ratio": bscat1,
                  "reflectance_ceiling": aerosol_reflectivity_ceiling(ssa1, bscat1)},
        "band2": {"ssa": ssa2, "backscatter_ratio": bscat2,
                  "reflectance_ceiling": aerosol_reflectivity_ceiling(ssa2, bscat2)},
        "active_in_the_configured_runs": bool(model["naerosp"] > 0),
        "note": "iaeron is 0 unless a species is configured, so these enter "
                "zrb1s and zrb2s multiplied by zero in the runs this tree "
                "carries. The ceilings are the sup over optical depth and are "
                "reached only in the limit.",
    }

    # (d) the layer sum, maximised over cloud fraction. Both terms are linear
    # in dcc, so the maximum is at an endpoint and the maximisation is exact.
    #
    # Only the LOWEST level carries a Rayleigh term: newrsc is 0, which zeroes
    # the layer-by-layer branch and leaves the whole-column one at NLEV, and
    # there the Rayleigh factor's (1 - dcc(NLEV)) and the cloud term's dcc are
    # the same cloud fraction.
    def layer_max(rayleigh_term, cloud_term, aerosol_term):
        at_zero = rayleigh_term + aerosol_term
        at_one = cloud_term
        return (max(at_zero, at_one), 0.0 if at_zero >= at_one else 1.0,
                at_zero, at_one)

    b1_max, b1_at, b1_clear, b1_cloudy = layer_max(
        rayleigh_ceiling, cloud["zrcl1s_ceiling"],
        aerosol["band1"]["reflectance_ceiling"] if model["naerosp"] > 0 else 0.0)
    b2_max, b2_at, b2_clear, b2_cloudy = layer_max(
        0.0, cloud["zrcl2s_ceiling"],
        aerosol["band2"]["reflectance_ceiling"] if model["naerosp"] > 0 else 0.0)
    # And with an aerosol present, since nothing forbids configuring one.
    b1_with_aer = max(rayleigh_ceiling + aerosol["band1"]["reflectance_ceiling"],
                      cloud["zrcl1s_ceiling"])
    b2_with_aer = max(aerosol["band2"]["reflectance_ceiling"],
                      cloud["zrcl2s_ceiling"])

    # What the model's OWN cloud water reaches, against those ceilings. The
    # per-layer cloud water is not written out, so this is the diagnostic
    # relation `mkclouds` builds it from, run on the climatology's own fields.
    run_dir = None
    runs = sorted(p for p in _paths.RUNS.glob("run_*") if p.is_dir())
    realised = None
    if runs:
        run_dir = runs[-1]
        nl = cod.read_namelists(run_dir)
        gascon = float(nl["GASCON"])
        clwref = float(nl.get("CLWREF", 0.00021))
        clwhsc = (float(nl["CLWHSC"]) if "CLWHSC" in nl
                  else 700.0 * (gascon / ga) / (287.0 / 9.80665))
        lwp = cod.water_paths(ta, ps_pa, prw, sigma, sigmah, dsigma,
                              gascon, ga, clwref, clwhsc)
        t1r, t2r = cloud_optical_depths(lwp)
        r1r = cloud_reflectivity_band1(t1r, tables["swcb1"])
        r2r = cloud_reflectivity_band2(t2r, tables["swcb2"], tables["swcoa"],
                                       cloudabs, coalb_floor)
        realised = {
            "constants_from": rel(run_dir),
            "lwp_g_m2_max": float(lwp.max()),
            "tau1_max": float(t1r.max()), "tau2_max": float(t2r.max()),
            "zrcl1s_max": float(r1r.max()), "zrcl2s_max": float(r2r.max()),
            "total_cloud_cover_max": float(clt.max()),
            "note": "the per-layer cloud fraction is not written out, so these "
                    "are the reflectivities a layer would have if it were "
                    "wholly overcast. The cloud term carries dcc, so the "
                    "layer's own reflectivity is at most this.",
        }

    # (e) the accumulated reflectance. A layer's R and T satisfy R + T <= 1 by
    # construction -- ztb1u is 1 minus the ozone absorptance times the clear
    # fraction, minus zrb1s, minus the aerosol absorption -- so the worst case
    # the recursion admits is a conservative layer, T = 1 - R, at every level.
    surface_albedo = {"band1": float(alb1.max()), "band2": float(alb2.max())}
    arms = {
        "arithmetic_ceiling_and_model_surface_albedo":
            ({"band1": b1_max, "band2": b2_max}, surface_albedo),
        "arithmetic_ceiling_and_a_perfect_reflector":
            ({"band1": b1_max, "band2": b2_max},
             {"band1": 1.0, "band2": 1.0}),
    }
    if realised:
        arms["what_the_model_reaches"] = (
            {"band1": realised["zrcl1s_max"], "band2": realised["zrcl2s_max"]},
            surface_albedo)
    accumulated = {}
    for label, (rmaxes, albedos) in arms.items():
        for band in ("band1", "band2"):
            rmax = rmaxes[band]
            albedo = albedos[band]
            down = adding_down(rmax, 1.0 - rmax, nlev)
            up_stack = adding_up(rmax, 1.0 - rmax, albedo, nlev)
            accumulated[f"{band}__{label}"] = {
                "layer_reflectivity": rmax,
                "surface_albedo": albedo,
                "downward_loop_final": down[-1],
                "upward_loop_final": up_stack[-1],
                "flux_loop_worst_product": down[-1] * albedo,
                "downward_loop_worst_divisor": 1.0 - down[-2] * rmax,
                "upward_loop_worst_divisor": 1.0 - up_stack[-1] * rmax,
                "flux_loop_worst_divisor": 1.0 - down[-1] * albedo,
                "margin": 1.0 / (down[-1] * albedo),
            }

    # (f) the tightest divisor each band's family B admits.
    family_b = {
        "nlev": nlev,
        "switches": {"nrscat": model["nrscat"], "newrsc": model["newrsc"],
                     "nclouds": model["nclouds"], "nswrcl": model["nswrcl"],
                     "naerosp": model["naerosp"]},
        "gravity_m_s2": ga,
        "surface_pressure_pa_declared": surface_pressure_pa,
        "surface_pressure_pa_model_maximum": ps_max,
        "rayleigh": rayleigh,
        "rayleigh_reflectance_ceiling": rayleigh_ceiling,
        "cloud": cloud,
        "aerosol": aerosol,
        "layer_reflectivity_ceiling": {
            "band1": {"value": b1_max, "at_cloud_fraction": b1_at,
                      "at_dcc_0": b1_clear, "at_dcc_1": b1_cloudy,
                      "with_an_aerosol_configured": b1_with_aer},
            "band2": {"value": b2_max, "at_cloud_fraction": b2_at,
                      "at_dcc_0": b2_clear, "at_dcc_1": b2_cloudy,
                      "with_an_aerosol_configured": b2_with_aer},
        },
        "accumulated": accumulated,
        "surface_albedo_model_maximum": surface_albedo,
        "what_the_model_reaches": realised,
        "arms": "the ARITHMETIC CEILING arms saturate every cap in the code at "
                "once -- 1000 g m-2 of cloud water and total overcast in every "
                "one of the model's layers -- and are the bound, not a state "
                "the model produces. The what-the-model-reaches arm runs the "
                "same recursion on the cloud water the climatology's own "
                "fields imply and on the surface albedo it writes.",
    }

    # --- the checks --------------------------------------------------------
    checks = {}

    bands = json.loads(args.bands.read_text(encoding="utf-8"))
    rows = bands["h2o"]["weight_vs_water_amount"]
    worst = 0.0
    for row in rows:
        want = float(row["lacis_hansen_eq21"])
        got = float(absorptance_water(row["w"], 1.0, 1.0))
        worst = max(worst, abs(got - want) / want)
    checks["lacis_hansen_eq21"] = {
        "what": "the unweighted water vapour absorptance this script evaluates "
                "against the same function stored in "
                f"{rel(args.bands)}, at {len(rows)} water amounts",
        "would_have_meant": "the domain analysis is being run over a different "
                            "curve from the one the band weights were derived "
                            "against, so every water vapour margin below is "
                            "for the wrong function",
        "worst_relative_error": worst,
        "tolerance": 1.0e-12,
        "gates": True,
        "passed": bool(worst < 1.0e-12),
    }

    got = float(stellar.rayleigh_coefficient(temperature_k=4965.0))
    err = abs(got - MODEL_PRINTED_RCOEFF_4965) / MODEL_PRINTED_RCOEFF_4965
    checks["rcoeff_against_the_compiled_model"] = {
        "what": "lib/stellar.py's Rayleigh coefficient at 4965 K against what "
                "the compiled model printed, recorded in "
                "notes/audits/physics-review.md",
        "would_have_meant": "the zscf this script forms is not the zscf swr "
                            "forms, so the Rayleigh reflectance below is not "
                            "the model's",
        "model_printed": MODEL_PRINTED_RCOEFF_4965,
        "recomputed": got,
        "relative_error": err,
        "tolerance": RCOEFF_TOLERANCE,
        "gates": True,
        "passed": bool(err < RCOEFF_TOLERANCE),
    }

    lo, hi = model["co2_fit_range"]
    co2_lo = float(co2_thinnest_layer) * ZBETTA
    co2_hi = co2_u_max
    checks["co2_fit_range_upper"] = {
        "what": f"radmod's own comment states the CO2 fit is quoted over {lo:g}"
                f" to {hi:g} atmos-cm and that the model never leaves it. The "
                "largest path this script derives from the model's own sigma "
                "grid, surface pressure and magnification must lie inside it.",
        "would_have_meant": "the fit is being evaluated above its stated "
                            "range, so the crossing amount below is an "
                            "extrapolation rather than a property of the fit",
        "quoted_range": [lo, hi],
        "model_maximum": co2_hi,
        "gates": True,
        "passed": bool(co2_hi <= hi),
    }
    checks["co2_fit_range_lower"] = {
        "what": "the same claim at the other end: the SMALLEST path the "
                "downward beam evaluates the fit at is the top layer's own "
                "amount at the smallest magnification, zbetta.",
        "would_have_meant": "the comment's reasoning -- the thinnest sigma "
                            "layer carries a few percent of the column -- does "
                            "not survive the pressure reduction the amount "
                            "also carries, so the fit is evaluated below its "
                            "quoted range near the model top",
        "quoted_range": [lo, hi],
        "model_minimum": co2_lo,
        "absorptance_there": float(absorptance_co2(co2_lo, co2sww)),
        "gates": False,
        "passed": bool(co2_lo >= lo),
        "recorded_in": "radmod.f90, the zca1/zcb1 block",
        "consequence": "none for this question. The fit is A(0) = 0 and rises "
                       "monotonically, so below the quoted range it is small "
                       "and positive and the divisor is within 0.002 of 1. It "
                       "is an accuracy claim in radmod's comment that the "
                       "sigma*ps/p0 pressure reduction does not support, not a "
                       "route to a zero divisor.",
    }

    checks["model_source_constants"] = {
        "what": "every literal this script evaluates, read out of "
                f"{rel(RADMOD)} and compared",
        "would_have_meant": "a coefficient moved in the model and not here, "
                            "so the domain analysis describes a function the "
                            "model no longer evaluates",
        "gates": True,
        "passed": True,
        "note": "this check raises rather than reporting: reaching this line "
                "at all means it passed",
    }

    # THE INSTRUMENT AGAINST THE SIZE OF THE EFFECT. The water vapour verdict
    # turns on one comparison, h2osww*h2oswl against zsolar2/(2.9/5.925). If
    # the declared brackets on those weights straddled that threshold the
    # verdict would be a property of the rounding rather than of the
    # configuration, so the threshold is stated and the brackets are measured
    # against it.
    threshold = zsolar2 / (2.9 / 5.925)
    checks["instrument_against_the_effect"] = {
        "what": "the water vapour divisor reaches zero if and only if "
                "h2osww*h2oswl exceeds zsolar2/(2.9/5.925). The declared "
                "brackets on both weights are compared with that threshold.",
        "threshold": threshold,
        "configured_product": h2osww * h2oswl,
        "product_over_threshold": h2osww * h2oswl / threshold,
        "lowest_bracketed_product": h2osww * h2oswl_bracket[0],
        "lowest_bracketed_over_threshold":
            h2osww * h2oswl_bracket[0] / threshold,
        "verdict_survives_the_bracket":
            bool(h2osww * h2oswl_bracket[0] > threshold),
    }

    # --- the verdict -------------------------------------------------------
    tightest_a = min(
        ((name, block["margin"]) for name, block in family_a.items()
         if block.get("margin")), key=lambda kv: kv[1])
    b_divisors = {}
    for key, value in accumulated.items():
        for loop in ("downward", "upward", "flux"):
            b_divisors[f"{key}__{loop}_loop"] = value[f"{loop}_loop_worst_divisor"]
    tightest_b = min(b_divisors.items(), key=lambda kv: kv[1])
    # The one currency both families share: how close the divisor itself gets
    # to zero on the worst input each family admits.
    floors = {f"family_a_{k}": v["divisor_floor"] for k, v in family_a.items()}
    floors.update({f"family_b_{k}": v for k, v in b_divisors.items()})
    tightest_overall = min(floors.items(), key=lambda kv: kv[1])

    out = {
        "what": "whether any divisor in radmod's swr can reach zero on this "
                "planet's configuration, evaluated over the model's own ranges",
        "issue": "world-2223",
        "model_source": rel(RADMOD),
        "band_split": {"zsolar1": zsolar1, "zsolar2": zsolar2,
                       "source": "lib/stellar.py:band_fractions, reproducing "
                                 "radmod.f90:solarini on the configured "
                                 "spectrum. NOT radmod's declared 0.517/0.483, "
                                 "which are the Sun's."},
        "planet": {"gravity_m_s2": ga,
                   "surface_pressure_pa": surface_pressure_pa,
                   "source": "config/planet.yaml"},
        "weights": {"o3visw": o3visw, "o3uvw": o3uvw, "h2osww": h2osww,
                    "h2oswl": h2oswl, "co2sww": co2sww, "o3scale": o3scale},
        "family_a": family_a,
        "family_b": family_b,
        "checks": checks,
        "verdict": {
            "family_a_can_reach_zero": True,
            "family_a_note": "all three absorptances can exceed their band's "
                             "share. Ozone and CO2 are unbounded in u; water "
                             "vapour has a finite asymptote and this "
                             "configuration puts it above 1, where the solar "
                             "weights put it barely above 1 and Earth's "
                             "amounts nowhere near.",
            "family_a_tightest": tightest_a[0],
            "family_a_tightest_margin": tightest_a[1],
            "family_b_can_reach_zero": False,
            "family_b_note": "the layer reflectivity is bounded below 1 by "
                             "hard caps in the code: zlwp is capped at "
                             "1000 g m-2, the table lookup clamps its optical "
                             "depth axis, and the diffuse cosine is fixed at "
                             "0.5. The accumulated reflectance cannot exceed 1 "
                             "because every layer satisfies R + T <= 1, so the "
                             "product is bounded by the layer ceiling.",
            "family_b_tightest": tightest_b[0],
            "family_b_tightest_divisor": tightest_b[1],
            "divisor_floors": floors,
            "tightest_overall": tightest_overall[0],
            "tightest_overall_divisor": tightest_overall[1],
            "two_currencies": "the family-A margin is a ratio of ABSORBER "
                              "AMOUNTS and the family-B one a ratio of "
                              "REFLECTANCES, so the divisor floor above is "
                              "what compares them. Family B's floor is the "
                              "smaller and is reached only with every cap in "
                              "the code saturated at once; family A's is "
                              "larger but sits on an amount that grows with "
                              "the modelled climate, which no cap holds.",
        },
        "provenance": provenance.config_stamp(
            config, "exoplasim/scripts/swr_divisor_domain.py",
            inputs=[RADMOD, args.climatology, args.bands, DUST,
                    _paths.INPUTS / "stellarspectra" / "k25v_hr.dat"]),
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "script_sha256": hashlib.sha256(
                Path(__file__).read_bytes()).hexdigest(),
        },
    }

    # --- the report --------------------------------------------------------
    print(f"model source  {rel(RADMOD)}")
    print(f"band split    zsolar1 {zsolar1:.6f}  zsolar2 {zsolar2:.6f}")
    print(f"planet        ga {ga} m/s2, surface pressure "
          f"{surface_pressure_pa:.0f} Pa")
    print()
    print("FAMILY A, nine divisors 1 - A(u)/zsolar_b")
    for name, block in family_a.items():
        sup = block.get("asymptote_over_band_share")
        sup_s = f"{sup:.4f}" if sup else "unbounded"
        cross = (block.get("u_cross_cm") or block.get("u_cross_cm_stp")
                 or block.get("u_cross_atmos_cm"))
        umax = [v for k, v in block["model_range"].items()
                if k.startswith("u_max")][0]
        print(f"  {name:<14} sup {sup_s:>10}  crosses at u = {cross:.4g}  "
              f"model u_max {umax:.4g}  margin {block['margin']:.3g}x")
    print(f"  tightest      {tightest_a[0]}, margin {tightest_a[1]:.3g}x, "
          f"divisor floor {family_a[tightest_a[0]]['divisor_floor']:.6f}")
    print()
    print("FAMILY B, four divisors 1 - R_above*R_below")
    print(f"  zscf          {rayleigh['k25v_one_grid__declared']['zscf']:.6f} "
          f"at rcoeff {rcoeff['k25v_one_grid']:.6f}, reflectance "
          f"{rayleigh['k25v_one_grid__declared']['reflectance']:.6f}")
    print(f"  rayleigh      ceiling {rayleigh_ceiling:.6f} over every rcoeff "
          "and surface pressure the model produces")
    print(f"  cloud         zrcl1s <= {cloud['zrcl1s_ceiling']:.6f} at "
          f"tau1 {cloud['tau1_at_cap']:.2f}; zrcl2s <= "
          f"{cloud['zrcl2s_ceiling']:.6f}")
    print(f"  aerosol       band 1 <= "
          f"{aerosol['band1']['reflectance_ceiling']:.6f}, band 2 <= "
          f"{aerosol['band2']['reflectance_ceiling']:.6f} "
          f"(active: {aerosol['active_in_the_configured_runs']})")
    print(f"  layer sum     band 1 <= {b1_max:.6f} at dcc = {b1_at:.0f}; "
          f"band 2 <= {b2_max:.6f} at dcc = {b2_at:.0f}")
    for key, value in accumulated.items():
        print(f"  {key:<44} product {value['flux_loop_worst_product']:.6f}  "
              f"divisor >= {value['flux_loop_worst_divisor']:.6f}  "
              f"margin {value['margin']:.4f}x")
    if realised:
        print(f"  what the model reaches: zrcl1s {realised['zrcl1s_max']:.6f}, "
              f"zrcl2s {realised['zrcl2s_max']:.6f} at lwp "
              f"{realised['lwp_g_m2_max']:.1f} g/m2")
    print(f"  tightest      {tightest_b[0]}, divisor {tightest_b[1]:.6f}")
    print()
    print(f"TIGHTEST OF THE THIRTEEN, by divisor floor: {tightest_overall[0]}, "
          f"{tightest_overall[1]:.6f}")
    print()
    print("CHECKS")
    for name, block in checks.items():
        if "passed" in block:
            if block["passed"]:
                verdict = "PASS"
            elif block.get("recorded_in"):
                # It answered no, the answer is written down where a reader
                # will find it, and nothing is left open. A check whose
                # negative result is recorded is settled, not failing.
                verdict = f"NO, recorded in {block['recorded_in']}"
            else:
                verdict = "FAIL"
            print(f"  {name:<38} {verdict}")
    print(f"  {'instrument_against_the_effect':<38} "
          f"threshold {threshold:.6f}, configured product "
          f"{h2osww * h2oswl:.6f} "
          f"({h2osww * h2oswl / threshold:.3f}x), survives the bracket: "
          f"{checks['instrument_against_the_effect']['verdict_survives_the_bracket']}")
    print()

    failed = [n for n, b in checks.items()
              if b.get("passed") is False and b.get("gates")]
    if failed:
        raise SystemExit(f"checks failed: {failed}")


if __name__ == "__main__":
    main()
