"""The environmental lapse rate, and the height the lowest model level sits at.

Both are the same hypsometric `R T / g`, which is why they share a module.

Earth's 6.5 K/km was hardcoded in three places (PHYS-12), none of which said it
was Earth's, and the note that bracketed it explored 5.5 to 6.5 -- entirely on
the shallow side of any defensible value for a world whose dry adiabat, `g/cp`,
is 31% steeper than Earth's. The fix is the same one `lib/sensitivity.py`
embodies: MEASURE the quantity at call time from the artifact the config names,
so it moves when the climate does, and keep one implementation instead of three
copies.

Method. Annual (or per-bin) mean temperature on the 10 sigma levels, weighted by
records per output bin via `lib/climatology.py`; level heights above the surface
by hypsometric integration with `R` and `cp` DERIVED from the configured
composition rather than declared; a least-squares slope of T against z per land
column over a stated sigma window; then the land-area-weighted mean of the
per-column rates.

The window is sigma 0.45 to 0.90, which on this atmosphere's scale height is
about 1.1 to 3.5 km above the surface -- the altitude band every consumer
extrapolates across (the freezing-height criterion reaches to 4.6 km peaks, the
dust layer's mass-weighted height is one scale height). Below it the surface
layer carries inversions and superadiabatic skin effects; above it the fit
starts to feel the tropopause. Moving the upper edge to sigma 0.35 shifts the
annual answer by 2%, which is the window sensitivity and is smaller than the
seasonal spread.

Two seasons are exposed, because the consumers want different ones:

  - `season="annual"`: the mean-state rate, for anything that extrapolates an
    annual-mean temperature (the dust layer's radiating temperature).
  - `season="warmest"`: the rate in each cell's warmest output bin, for
    anything that extrapolates the warmest month (the freezing-height criterion
    `z* = cell_mean + (T_warmest - 273.15) / lapse`, and the basemap's ice
    shading). Warm-season columns convect closer to the adiabat, so this rate
    is systematically steeper than the annual one.

Two checks that can fail, both raised on rather than warned about: the measured
rate must be POSITIVE and BELOW the dry adiabat (a violation means the window
strayed into an inversion or the fit is reading the stratosphere), and the
land-mean rms of the linear fit must stay under 1 K (a violation means the
window is not a straight line and the single number is not a faithful summary).

What is deliberately NOT a check: a moist-adiabat floor. The measured annual
rate sits BELOW the moist adiabatic rate evaluated at the window-mean state
(stable stratification does that), so flooring at the moist rate would silently
decide the answer instead of checking it -- the same failure as the Penman
floor that once decided 73% of a carve verdict.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

import climatology
from paths import climatology_path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT_ROOT / "config" / "planet.yaml"

# Molar masses (g/mol) and specific heats (J/kg/K) of the constituents the
# config can declare. Physical constants, not tunings.
MOLAR_MASS = {"pN2_bar": 28.0134, "pO2_bar": 31.9988,
              "pAr_bar": 39.948, "pCO2_bar": 44.0098}
SPECIFIC_HEAT = {"pN2_bar": 1040.0, "pO2_bar": 918.0,
                 "pAr_bar": 520.0, "pCO2_bar": 844.0}
R_UNIVERSAL = 8314.462618          # J/kmol/K

SIGMA_WINDOW = (0.45, 0.90)
FIT_RMS_LIMIT_K = 1.0

# The model's bottom full level at NLEV = 10. PlaSim sets the half levels first
# and each full level as the midpoint of the pair around it (plasim.f90:1644),
# so this is a property of the level count and not of the planet. Passed rather
# than read from a climatology so a caller that runs before any run exists can
# still ask; a caller holding the model's `lev` axis should pass its bottom
# entry instead.
SIGMA_LOWEST = 0.9828


def reference_height_m(t_air_k: float, cfg: dict | None = None,
                       sigma_lowest: float = SIGMA_LOWEST) -> float:
    """Height of the lowest model level above the surface, metres.

    Hypsometric, with the same `R T / g` scale height `_column_rates` integrates
    with, so the height a bulk transfer coefficient is derived over and the
    heights a lapse rate is fitted against come from one expression.

    THE GRAVITY IS THIS PLANET'S. Two consumers hardcoded the Earth-gravity
    answer, 141.6 m, which is high by the gravity ratio and inflates every
    `ln(z_ref/z0)` it feeds. `ce = k^2 / ln(z_ref/z0)^2` is the quantity
    roughness acts through, so a reference height 31% too high moves `ce` at
    ExoPlaSim's uniform `dz0land` by 14% and the bare-to-vegetated CONTRAST it
    is quoted for by 9%.

    Linear in the air temperature, so a caller with no climatology to measure
    one from brackets it rather than picking a value. `t_air_k` may be a field,
    in which case the height comes back with its shape.
    """
    cfg = _config(cfg)
    r_specific, _cp = gas_properties(cfg)
    gravity = float(cfg["planet"]["gravity_m_s2"])
    height = (r_specific * np.asarray(t_air_k, dtype=float) / gravity
              ) * float(np.log(1.0 / sigma_lowest))
    return float(height) if height.ndim == 0 else height


def _config(cfg: dict | None) -> dict:
    return cfg if cfg is not None else yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


def gas_properties(cfg: dict | None = None) -> tuple[float, float]:
    """(R, cp) in J/kg/K, derived from the configured composition."""
    atm = _config(cfg)["atmosphere"]
    total = sum(atm[k] for k in MOLAR_MASS)
    x = {k: atm[k] / total for k in MOLAR_MASS}
    mbar = sum(x[k] * MOLAR_MASS[k] for k in MOLAR_MASS)
    mass = {k: x[k] * MOLAR_MASS[k] / mbar for k in MOLAR_MASS}
    r_specific = R_UNIVERSAL / mbar
    cp = sum(mass[k] * SPECIFIC_HEAT[k] for k in MOLAR_MASS)
    return r_specific, cp


def dry_adiabat_k_per_km(cfg: dict | None = None) -> float:
    """g/cp for the configured gravity and composition."""
    cfg = _config(cfg)
    r_specific, cp = gas_properties(cfg)
    return 1000.0 * float(cfg["planet"]["gravity_m_s2"]) / cp


def _column_rates(ta: np.ndarray, sigma: np.ndarray, r_specific: float,
                  gravity: float) -> tuple[np.ndarray, np.ndarray]:
    """Per-column lapse rate (K/km) and fit rms (K) over the sigma window.

    `ta` is (lev, lat, lon) ordered like `sigma`. Heights are relative to the
    surface, which is all a slope needs.
    """
    order = np.argsort(-sigma)
    s, t = sigma[order], ta[order]
    z = np.zeros_like(t)
    for k in range(1, len(s)):
        z[k] = z[k - 1] + (r_specific * 0.5 * (t[k - 1] + t[k]) / gravity
                           ) * np.log(s[k - 1] / s[k])
    lo, hi = SIGMA_WINDOW
    sel = (s >= lo) & (s <= hi)
    zs, ts = z[sel], t[sel]
    zm = zs - zs.mean(axis=0)
    tm = ts - ts.mean(axis=0)
    slope = (zm * tm).sum(axis=0) / (zm * zm).sum(axis=0)
    rms = np.sqrt(((tm - slope * zm) ** 2).mean(axis=0))
    return -slope * 1000.0, rms


def environmental_lapse_k_per_km(cfg: dict | None = None,
                                 climatology_file: Path | None = None,
                                 season: str = "annual") -> float:
    """Land-area-weighted environmental lapse rate, K/km, measured at call time.

    `season="annual"` fits the bin-weighted annual-mean profile;
    `season="warmest"` fits, in each cell, the profile of the output bin whose
    lowest-level temperature is warmest there, which is the rate the
    freezing-height extrapolation wants.
    """
    from netCDF4 import Dataset
    from numpy.polynomial.legendre import leggauss

    cfg = _config(cfg)
    r_specific, _cp = gas_properties(cfg)
    gravity = float(cfg["planet"]["gravity_m_s2"])
    path = climatology_file or climatology_path()
    with Dataset(path) as ds:
        centres = np.asarray(ds["time"][:], dtype=float)
        sigma = np.asarray(ds["lev"][:], dtype=float)
        ta_bins = np.asarray(ds["ta"][:], dtype=float)
        land = climatology.annual_mean(np.asarray(ds["lsm"][:], dtype=float),
                                       centres) >= 0.5

    nlat, nlon = land.shape
    weights = leggauss(nlat)[1][::-1][:, None] * np.ones((1, nlon)) * land
    weights = weights / weights.sum()

    if season == "annual":
        rate, rms = _column_rates(climatology.annual_mean(ta_bins, centres),
                                  sigma, r_specific, gravity)
    elif season == "warmest":
        lowest = int(np.argmax(sigma))
        warm_bin = ta_bins[:, lowest].argmax(axis=0)
        rate = np.zeros((nlat, nlon))
        rms = np.zeros((nlat, nlon))
        for b in range(ta_bins.shape[0]):
            mask = warm_bin == b
            if mask.any():
                rate_b, rms_b = _column_rates(ta_bins[b], sigma, r_specific,
                                              gravity)
                rate[mask] = rate_b[mask]
                rms[mask] = rms_b[mask]
    else:
        raise ValueError(f"season must be 'annual' or 'warmest', not {season!r}")

    mean_rate = float((rate * weights).sum())
    mean_rms = float((rms * weights).sum())
    ceiling = dry_adiabat_k_per_km(cfg)
    if not 0.0 < mean_rate < ceiling:
        raise RuntimeError(
            f"measured lapse rate {mean_rate:.3f} K/km is outside (0, "
            f"{ceiling:.2f}): the fit window is reading an inversion or the "
            "stratosphere and the number must not be used")
    if mean_rms > FIT_RMS_LIMIT_K:
        raise RuntimeError(
            f"the temperature profile is not linear over the sigma window "
            f"(land-mean fit rms {mean_rms:.2f} K > {FIT_RMS_LIMIT_K}); a "
            "single rate is not a faithful summary of it")
    return mean_rate
