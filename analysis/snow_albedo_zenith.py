#!/usr/bin/env python3
"""What a solar zenith term is worth to the modelled snow albedo, on this star.

    python analysis/snow_albedo_zenith.py

WORLD-SPD. `landmod.f90` sets the snow albedo as a linear ramp in SURFACE
TEMPERATURE with no zenith angle anywhere, while `radmod.f90:swr` already
overwrites the ICE-FREE OCEAN albedo with a zenith form by default. Same
predictor, one surface, in the same low-sun regime. This sizes what is missing.

The form is the asymptotic semi-infinite result, `alpha(mu0) = alpha_diffuse **
K(mu0)` with van de Hulst's escape function `K(mu0) = (3/7)(1 + 2*mu0)`. It has
NO fitted coefficient: `K` is normalised by `2*int K(mu) mu dmu = 1`, which
fixes `K(2/3) = 1`, so a model carrying no zenith angle is implicitly quoting
the `mu0 = 2/3` value and the term is a redistribution about the number the
model already has. Its one input, `alpha_diffuse`, is taken PER BAND from the
model's own star-weighted `dsnowalbmn`/`dsnowalbmx`, so nothing is inherited
from an Earth calibration.

`notes/audits/snow-albedo-zenith.md` carries the argument and the verdict on
whether the open-ocean formula is the right form to reuse. This script is the
arithmetic behind that note's tables, and it runs from anywhere.
"""
from pathlib import Path
import sys
ROOT = str(Path(__file__).resolve().parents[1]) + "/"
sys.path.insert(0, ROOT + "lib")
import numpy as np, yaml, json
import stellar, sensitivity

cfg = yaml.safe_load(open(ROOT + "config/planet.yaml"))
OBLIQ = float(cfg["planet"]["obliquity_degrees"])
S0 = float(cfg["orbit"]["earth_solar_constant_w_m2"]) * float(cfg["orbit"]["baseline_flux_earth"])
F1 = float(stellar.band_fractions()[0])
F1_SUN = float(stellar.blackbody_band_fractions(5772.0)[0])

# The model's OWN k25v-weighted band albedos, from analysis/ice_albedo.py's
# RECORDED block: what radini printed under "Finalized Albedos".
MODEL = {"dsnowalbmx": (0.98278401164203333, 0.58529460781319775),
         "dsnowalb":   (0.76382640121202905, 0.40611515565148087),
         "dsnowalbmn": (0.50782484289747332, 0.27261180703024718)}

def K(mu):
    """van de Hulst escape function; 2*int K(mu) mu dmu == 1, K(2/3) == 1."""
    return (3.0 / 7.0) * (1.0 + 2.0 * mu)

def alpha_dir(a_dif, mu):
    """Asymptotic semi-infinite direct-beam albedo. Bounded in [0,1] by form."""
    return np.power(a_dif, K(mu))

def broadband(b1, b2, mu, f1=F1):
    return f1 * alpha_dir(b1, mu) + (1.0 - f1) * alpha_dir(b2, mu)

def spans(b1, b2, f1=F1):
    return {"band1_diffuse": b1, "band2_diffuse": b2,
            "broadband_diffuse": f1 * b1 + (1 - f1) * b2,
            "band1_mu1": float(alpha_dir(b1, 1.0)), "band1_mu0": float(alpha_dir(b1, 0.0)),
            "band2_mu1": float(alpha_dir(b2, 1.0)), "band2_mu0": float(alpha_dir(b2, 0.0)),
            "broadband_mu1": float(broadband(b1, b2, 1.0, f1)),
            "broadband_mu0": float(broadband(b1, b2, 0.0, f1)),
            "broadband_span": float(broadband(b1, b2, 0.0, f1) - broadband(b1, b2, 1.0, f1))}

d = np.loadtxt(ROOT + "references/exocam/tools/spectral_albedos/snow100um.txt")
k25 = stellar.band_reflectances(d[:, 0], d[:, 1])
sun = stellar.band_reflectances(d[:, 0], d[:, 1], temperature_k=5772.0)

def day(lat_deg, decl_deg, b1, b2, n=20001):
    """Flux-weighted daily-mean albedo and the daily-mean insolation factor."""
    lat, dec = np.deg2rad(lat_deg), np.deg2rad(decl_deg)
    h = np.linspace(-np.pi, np.pi, n)
    mu = np.clip(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h), 0.0, None)
    f = np.trapezoid(mu, h)
    if f <= 0:
        return None, 0.0
    a = broadband(b1, b2, mu)
    return float(np.trapezoid(a * mu, h) / f), float(f / (2 * np.pi))

B1, B2 = MODEL["dsnowalbmx"]                      # cold snow: the polar case
A_DIF = F1 * B1 + (1 - F1) * B2

rows = []
for lat in (0, 15, 30, 45, 60, 75):
    for season, dec in (("equinox", 0.0), ("summer", OBLIQ), ("winter", -OBLIQ)):
        a, insol = day(lat, dec, B1, B2)
        if a is None:
            rows.append({"lat": lat, "season": season, "polar_night": True}); continue
        rows.append({"lat": lat, "season": season, "alpha_flux_weighted": a,
                     "delta_vs_model": a - A_DIF,
                     "daily_insolation_w_m2": insol * S0,
                     "absorbed_change_w_m2": -(a - A_DIF) * insol * S0})

def cap(lat0):
    lats, decs = np.linspace(lat0, 89.5, 120), np.linspace(-OBLIQ, OBLIQ, 61)
    num = den = 0.0
    for la in lats:
        w = np.cos(np.deg2rad(la))
        for de in decs:
            a, insol = day(la, de, B1, B2, n=4001)
            den += w
            if a is not None:
                num += w * (-(a - A_DIF)) * insol * S0
    return num / den

price = {f"poleward_of_{l}": cap(l) for l in (45, 60, 70)}
kelvin = {f"albedo_{a}": {f"snow_fraction_{f}":
          sensitivity.forcing_to_kelvin(price["poleward_of_60"] * f, a, cfg)
          for f in (0.02, 0.05, 0.10, 0.20)} for a in (0.25, 0.30, 0.35)}

report = {
    "obliquity_degrees": OBLIQ, "top_of_atmosphere_w_m2": S0,
    "band1_flux_fraction_k25v": F1, "band1_flux_fraction_sun5772": F1_SUN,
    "model_finalized": {k: spans(*v) for k, v in MODEL.items()},
    "exocam_snow_k25v": {**k25, **spans(k25["band1"], k25["band2"])},
    "exocam_snow_sun5772": {**sun, **spans(sun["band1"], sun["band2"], F1_SUN)},
    "cold_snow_diffuse_broadband": A_DIF,
    "daily": rows, "cap_price_w_m2_per_unit_snow_area": price,
    "kelvin_bracket": kelvin,
    "form": "alpha(mu0) = alpha_diffuse ** K(mu0), K = (3/7)(1 + 2 mu0); "
            "K is normalised so 2*int K(mu) mu dmu = 1, which puts K(2/3) = 1, "
            "so a scheme with no zenith angle is quoting the mu0 = 2/3 value",
    "alpha_diffuse_source": "the model's own star-weighted dsnowalb* constants, "
                            "verified by analysis/ice_albedo.py",
    "note": "notes/audits/snow-albedo-zenith.md",
}
out = Path(ROOT) / "analysis" / "snow_albedo_zenith.json"
out.write_text(json.dumps(report, indent=2, default=float) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2, default=float))
print(f"\nwrote {out}")
