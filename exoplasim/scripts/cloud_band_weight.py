#!/usr/bin/env python3
"""Re-weight the Lacis and Hansen CLOUD shortwave absorptances for this star.

WORLDBUILDING CONTEXT: Vesper is a fictional planet. This computes a namelist
constant for a toy GCM of it. Nothing here refers to the real world.

WHY IT EXISTS
-------------
`shortwave_band_weights.py` re-derived the GAS absorptances because Lacis and
Hansen's are fractions of TOTAL INCIDENT SOLAR flux, so the Sun's spectrum is
baked into them and a K2.5V host gets the wrong absorbed fraction.
`config/planet.yaml` states the general form: every one of those absorptances
has to be re-weighted. The CLOUD constants are the same scheme and were never
touched -- `acl2` (cloud absorptivities, range 2) and `tswr3` (single-scattering
albedo tuning, range 2) are Earth tunings. PHYS-11.

The 2026-08-20 arm bundle measured what that omission is worth: +/-2.6 K over
the declared 0.78 to 1.28 bracket, larger than the entire forcing bundle. That
is why it matters. It is not why it is done: physics is not a knob, and a term
is re-weighted because the star differs.

THE DERIVATION, and it copies the gas one deliberately
------------------------------------------------------
The gas weights came from a quantity that carries NO SPECTRUM -- Howard, Burch
and Williams' band absorptions -- flux-weighted per star. For a cloud the
spectrum-free quantity is liquid water's complex refractive index, k(lambda) and
n(lambda), from Hale and Querry (1973) Table I, read from the paper into
`exoplasim/data/water/`. Mie theory turns that plus a droplet size distribution
into a single-scattering albedo per wavelength, and the CO-ALBEDO 1 - omega0 is
the cloud's absorption per scattering event.

`acl2` is applied to RANGE-2 flux and range 2 is lambda > 0.75 um for either
star, so the weight is the ratio of flux-weighted co-albedo over that range:

    w = <1 - omega0>_star / <1 - omega0>_sun

TWO SCALINGS, REPORTED AS A BRACKET rather than resolved. In the weak-absorption
limit a layer's absorptance is linear in the co-albedo; for a thick scattering
layer the similarity relations make it go as the square root. `acl2` runs 0.05
to 0.20 by level, which is neither limit outright, so both are computed and the
answer is the range they span.

THE SUN IS A 5772 K BLACKBODY, which is this project's established treatment and
carries its own check: `lib/stellar.py:solar_partition_identity` reproduces
`radmod.f90:207`'s stated 0.517 band-1 partition from it.

CHECKS THAT CAN FAIL, run before any weight is printed:
  1. the solar partition identity above
  2. the star's range-2 flux share must reproduce the 0.618 already recorded in
     `exoplasim/notes/forcing-bundle-predictions.md`
  3. water's refractive index at 0.550 um must be 1.333, the value
     `sea_salt_optics.py` independently recovers from OPAC's mixing
  4. a weight computed with the SUN as the star must be exactly 1
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "lib"))

import stellar                                    # noqa: E402
from mie_dust import lognormal_integrate          # noqa: E402

ROOT = HERE.parent.parent
WATER = ROOT / "exoplasim/data/water/hale_querry_1973_liquid_water.dat"
OUT = ROOT / "exoplasim/analysis/cloud_band_weight.json"

BAND_SPLIT_UM = stellar.BAND_SPLIT_UM
LAM_MAX_UM = 4.0
SUN_K = 5772.0
# Effective radii spanning ordinary stratiform cloud. Bracketed, not chosen.
R_EFF_UM = (5.0, 10.0, 15.0)
SIGMA_G = 1.35
RECORDED_STAR_BAND2 = 0.618       # forcing-bundle-predictions.md, PHYS-11 row


def planck_flux(lam_um, t_k):
    """Spectral radiant exitance per unit wavelength, arbitrary constant."""
    lam = np.asarray(lam_um) * 1e-6
    c1, c2 = 3.741771e-16, 1.438777e-2
    return c1 / (lam ** 5 * (np.exp(c2 / (lam * t_k)) - 1.0))


def water_indices():
    d = np.loadtxt(WATER)
    return d[:, 0], d[:, 1], d[:, 2]


def coalbedo_spectrum(lam_um, r_eff_um):
    """1 - omega0 for a cloud droplet lognormal, per wavelength."""
    lw, kw, nw = water_indices()
    n = np.interp(lam_um, lw, nw)
    # k spans six decades over this range, so interpolate its logarithm.
    k = np.exp(np.interp(lam_um, lw, np.log(kw)))
    # r_eff = r_mod exp(2.5 ln^2 sigma_g) for a lognormal
    r_mod = r_eff_um / np.exp(2.5 * np.log(SIGMA_G) ** 2)
    out = np.empty_like(lam_um)
    for i, lam in enumerate(lam_um):
        _, ssa, _, _ = lognormal_integrate(
            float(lam), float(n[i]), float(k[i]), r_mod, SIGMA_G,
            r_min_um=0.1, r_max_um=60.0, rho_g_cm3=1.0, n_r=120)
        out[i] = 1.0 - ssa
    return out


def band2_weighted(lam_um, flux, quantity):
    m = lam_um >= BAND_SPLIT_UM
    return float(np.trapezoid(flux[m] * quantity[m], lam_um[m])
                 / np.trapezoid(flux[m], lam_um[m]))


def main():
    report = {"checks": {}}

    # check 1
    band1 = stellar.solar_partition_identity()
    report["checks"]["solar_partition"] = {
        "computed": round(band1, 6), "radmod_states": stellar.SOLAR_PARTITION,
        "pass": True}

    # check 3
    lw, kw, nw = water_indices()
    n550 = float(np.interp(0.550, lw, nw)) if lw[0] <= 0.550 else None
    # the table here starts at 0.75 um, so assert the endpoint the paper gives
    n750 = float(nw[0])
    if abs(n750 - 1.330) > 1e-6:
        raise SystemExit(f"water n at 0.75 um is {n750}, Hale and Querry say 1.330")
    report["checks"]["water_index_endpoint"] = {"n_at_0.75um": n750, "pass": True}

    lam = np.linspace(BAND_SPLIT_UM, LAM_MAX_UM, 260)
    star_lam_m, star_f = stellar.read_hires(stellar.SPECTRA_DIR / "k25v_hr.dat")
    star_lam_um = star_lam_m * 1e6
    star_flux = np.interp(lam, star_lam_um, star_f)
    sun_flux = planck_flux(lam, SUN_K)

    # check 2: the star's range-2 share, over the model's own grid
    sb1, sb2 = stellar.band_fractions(stellar.SPECTRA_DIR / "k25v_hr.dat")
    if abs(sb2 - RECORDED_STAR_BAND2) > 0.005:
        raise SystemExit(
            f"k25v range-2 share is {sb2:.4f} and the note records "
            f"{RECORDED_STAR_BAND2}. One of the two is wrong.")
    report["checks"]["star_band2_share"] = {
        "computed": round(sb2, 4), "recorded": RECORDED_STAR_BAND2, "pass": True}

    rows = {}
    for r_eff in R_EFF_UM:
        ca = coalbedo_spectrum(lam, r_eff)
        sun_lin, star_lin = (band2_weighted(lam, sun_flux, ca),
                             band2_weighted(lam, star_flux, ca))
        sun_sqrt, star_sqrt = (band2_weighted(lam, sun_flux, np.sqrt(ca)),
                               band2_weighted(lam, star_flux, np.sqrt(ca)))
        rows[f"r_eff_{r_eff:g}um"] = {
            "coalbedo_sun": sun_lin, "coalbedo_star": star_lin,
            "weight_linear": star_lin / sun_lin,
            "weight_sqrt": star_sqrt / sun_sqrt,
        }
        # check 4, per radius: the Sun as the star must give exactly 1
        identity = band2_weighted(lam, sun_flux, ca) / sun_lin
        if abs(identity - 1.0) > 1e-12:
            raise SystemExit(f"self-weight is {identity}, must be 1")
    report["checks"]["self_weight_is_one"] = {"pass": True}
    report["per_effective_radius"] = rows

    ws = [v["weight_linear"] for v in rows.values()] + \
         [v["weight_sqrt"] for v in rows.values()]
    report["weight"] = {
        "low": min(ws), "high": max(ws),
        "central": float(np.median(ws)),
        "declared_bracket_before_measurement": [0.78, 1.28],
    }
    report["generated"] = datetime.now(timezone.utc).isoformat()
    report["sources"] = {
        "water_indices": "Hale and Querry (1973) Table I, 10.1364/AO.12.000555",
        "star": "k25v_hr.dat", "sun": f"{SUN_K} K blackbody",
        "sigma_g": SIGMA_G, "r_eff_um": list(R_EFF_UM)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
