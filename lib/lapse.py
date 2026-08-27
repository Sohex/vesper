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
# entry instead. `sigma_levels` re-derives it from the model's own construction
# and `check_sigma_lowest` refuses when the two disagree, so this literal is a
# declaration inside a loop rather than a copy.
SIGMA_LOWEST = 0.9828

# EARTH, NAMED. Standard-atmosphere constants used ONLY as the anchor a
# scaled-to-this-planet quantity is scaled FROM, never as a value for this
# world. Each is a published standard that nothing in this tree computes, so
# none of them can drift here; what would drift is a Vesper number derived from
# one of them and then written down, which is why the derivations below return
# the number instead of declaring it.
EARTH_GRAVITY_M_S2 = 9.80665                  # ISO 2533 standard gravity
EARTH_SURFACE_TEMPERATURE_K = 288.0           # plasimmod.f90's own `tgr`
EARTH_TROPOSPHERIC_LAPSE_K_PER_M = 6.5e-3     # ISO 2533 standard lapse rate
# Dry air, by volume, the mixture the Earth constants above describe. Fed
# through `gas_properties` so Earth's R and cp come out of the SAME expression
# this planet's do; a separate pair of literals would be two derivations of one
# thing.
EARTH_ATMOSPHERE = {"atmosphere": {"pN2_bar": 0.78084, "pO2_bar": 0.20946,
                                   "pAr_bar": 0.00934, "pCO2_bar": 0.000415}}

# `setzt`'s smoothing width at the tropopause, `plasimmod.f90`'s `dttrp`,
# in kelvin: the model divides it by `ct` at read and multiplies it back inside
# `setzt`, so the number the profile sees is this one.
TROPOPAUSE_SMOOTHING_K = 2.0


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


def scale_height_m(temperature_k: float, cfg: dict | None = None) -> float:
    """R T / g, the pressure scale height, for the configured atmosphere."""
    cfg = _config(cfg)
    r_specific, _cp = gas_properties(cfg)
    return r_specific * float(temperature_k) / float(cfg["planet"]["gravity_m_s2"])


def earth_scale_height_m(temperature_k: float = EARTH_SURFACE_TEMPERATURE_K) -> float:
    """The same expression on Earth's own R and g, as the anchor to scale from."""
    r_earth, _cp = gas_properties(EARTH_ATMOSPHERE)
    return r_earth * float(temperature_k) / EARTH_GRAVITY_M_S2


def pressure_length_ratio_to_earth(cfg: dict | None = None) -> float:
    """(R/g) here over (R/g) on Earth: a geometric length as a share of Earth's.

    An ISOTHERMAL ratio, which is what a length fitted to Earth's atmosphere
    should be scaled by when what it locates is a PRESSURE. `radmod.f90`'s ozone
    profile is the case: `mko3` builds its height coordinate hypsometrically
    from the model's own gascon and ga, correctly, and then places the profile
    against `bo3` and `co3`, which upstream carries as geometric metres fitted
    to Earth. The photochemical maximum is set by pressure-like conditions --
    ultraviolet optical depth and three-body recombination density -- not by
    geometric altitude, so pressure is the invariant to hold, and holding it
    reproduces Earth's own layer-by-layer share of the ozone column on this
    world's sigma levels. `exoplasim/notes/ozone.md` carries that check.

    No temperature enters, because the two lengths being compared are compared
    against each other at the same point in each profile.
    """
    cfg = _config(cfg)
    r_specific, _cp = gas_properties(cfg)
    r_earth, _cpe = gas_properties(EARTH_ATMOSPHERE)
    return ((r_specific / float(cfg["planet"]["gravity_m_s2"]))
            / (r_earth / EARTH_GRAVITY_M_S2))


def cold_start_lapse_rate_k_per_m(cfg: dict | None = None) -> float:
    """`ALR`: the cold start's tropospheric lapse rate, K/m.

    THE FRACTION OF NEUTRALITY IS WHAT TRANSFERS, not the rate. Earth's
    standard 6.5 K/km is a fixed share of Earth's own dry adiabat, and a profile
    that kept the Earth NUMBER on this atmosphere would keep neither the
    stratification nor the share of neutrality the number was chosen at: g/cp is
    12.75 K/km here against 9.76 on Earth. So the share is measured on Earth's
    constants and applied to this planet's.

    Not `environmental_lapse_k_per_km`: that reads a climatology, and a cold
    start is what runs before there is one.
    """
    cfg = _config(cfg)
    _r, cp_earth = gas_properties(EARTH_ATMOSPHERE)
    fraction = EARTH_TROPOSPHERIC_LAPSE_K_PER_M / (EARTH_GRAVITY_M_S2 / cp_earth)
    _r2, cp = gas_properties(cfg)
    return fraction * float(cfg["planet"]["gravity_m_s2"]) / cp


def cold_start_tropopause_height_m(surface_temperature_k: float,
                                   cfg: dict | None = None,
                                   earth_tropopause_m: float = 12000.0) -> float:
    """`DTROP`: where the cold start's profile turns isothermal, metres.

    THE NUMBER OF PRESSURE SCALE HEIGHTS IS WHAT TRANSFERS. Upstream's 12 km is
    a fixed count of scale heights at Earth's own R, T and g; a fixed geometric
    height on this atmosphere would put the turn at a different pressure, which
    is the wrong invariant for the place a temperature profile stops following
    the adiabat. `earth_tropopause_m` defaults to upstream's `dtrop` and is an
    argument so a caller holding the model source's own default can pass it
    rather than trusting this one.

    Linear in the surface temperature, so it moves when the cold start's
    surface temperature does -- which is the whole reason it is derived here
    instead of written down beside it.
    """
    return (float(earth_tropopause_m)
            * scale_height_m(surface_temperature_k, cfg) / earth_scale_height_m())


def sigma_levels(cfg: dict | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(sigma, sigmah, dsigma) as `plasim.f90` builds them, from config alone.

    The model's own construction at `model.vertical_grid`: the quartic
    half-level distribution, shifted so the top half-level lands on
    `model.model_top_hpa` and normalised so the last one is the surface, with
    each full level the midpoint of the pair around it. It is a property of the
    layer count and the model top and of nothing else, so it is answerable
    before any run exists -- which is what a cold-start quantity needs.

    A climatology's `levp` is NOT this: the file writes the midpoints between
    full levels, so its differences do not sum to one. Read `lev` and invert the
    midpoint recursion if a run's own grid is wanted.
    """
    cfg = _config(cfg)
    model = cfg["model"]
    layers = int(model["layers"])
    grid = int(model.get("vertical_grid", 0))
    if grid not in (0, 4):
        raise ValueError(
            f"model.vertical_grid is {grid}; this reproduces plasim.f90's "
            "quartic branches 0 and 4 only, and a construction it does not "
            "know would return a plausible grid that is not the model's.")
    k = np.arange(1, layers + 1, dtype=float) / layers
    sigmah = 0.75 * k + 1.75 * k ** 3 - 1.5 * k ** 4
    if grid == 4:
        # The whole of what branch 4 adds: shift the top half-level to zero,
        # normalise, and map affinely onto [ptop/psurf, 1].
        total_bar = sum(float(cfg["atmosphere"][key]) for key in MOLAR_MASS)
        top = float(model["model_top_hpa"]) * 100.0 / (total_bar * 1.0e5)
        sigmah = (sigmah - sigmah[0]) / (sigmah[-1] - sigmah[0])
        sigmah = sigmah * (1.0 - top) + top
    dsigma = np.diff(np.concatenate([[0.0], sigmah]))
    sigma = np.concatenate([[0.5 * sigmah[0]], 0.5 * (sigmah[:-1] + sigmah[1:])])
    return sigma, sigmah, dsigma


def check_sigma_lowest(cfg: dict | None = None,
                       tolerance: float = 5.0e-5) -> list[str]:
    """`SIGMA_LOWEST` against the construction that produces it. Empty means agreed."""
    sigma, _sigmah, _dsigma = sigma_levels(cfg)
    got = float(sigma[-1])
    if abs(got - SIGMA_LOWEST) > tolerance:
        return [f"lib/lapse.py declares SIGMA_LOWEST = {SIGMA_LOWEST} and the "
                f"model's own level construction at this layer count and model "
                f"top gives {got:.6f}"]
    return []


def setzt_profile(surface_temperature_k: float, lapse_rate_k_per_m: float,
                  tropopause_height_m: float,
                  cfg: dict | None = None) -> np.ndarray:
    """`plasim.f90:setzt`'s restoration temperature, per full sigma level, in K.

    Reproduced rather than approximated, because what rests on it is the
    dsigma-weighted mean below and an approximation would move that by more than
    the mean's own precision. The one-sided iteration, the two-pass height
    predictor and the hyperbolic smoothing at the tropopause are the model's.
    """
    cfg = _config(cfg)
    r_specific, _cp = gas_properties(cfg)
    gravity = float(cfg["planet"]["gravity_m_s2"])
    sigma, _sigmah, _dsigma = sigma_levels(cfg)
    tgr = float(surface_temperature_k)
    alr = float(lapse_rate_k_per_m)
    dtrop = float(tropopause_height_m)
    out = np.zeros_like(sigma)
    sigprev, tprev, zprev = 1.0, tgr, 0.0

    def smoothed(height: float) -> float:
        offset = 0.5 * alr * (height - dtrop)
        return (tgr - dtrop * alr
                + float(np.hypot(offset, TROPOPAUSE_SMOOTHING_K)) - offset)

    for index in range(len(sigma) - 1, -1, -1):
        span = float(np.log(sigprev / sigma[index]))
        first = smoothed(zprev + (r_specific * tprev / gravity) * span)
        mid = 0.5 * (tprev + first)
        second = smoothed(zprev + (r_specific * mid / gravity) * span)
        out[index] = second
        zprev += (0.5 * (second + tprev) * r_specific / gravity) * span
        tprev = second
        sigprev = float(sigma[index])
    return out


def semi_implicit_reference_temperature_k(cfg: dict | None = None,
                                          profile: dict | None = None) -> float:
    """`T0`: the dsigma-weighted mean of the cold start's `setzt` profile.

    The semi-implicit scheme's reference temperature is isothermal on every
    level and is what the reference geopotential is built out of, so the
    adiabatic conversion defect the core carries is proportional to it. The rule
    is the model's own: mass-weight the profile `setzt` builds from the cold
    start's surface temperature, lapse rate and tropopause height.

    `profile` defaults to the configured `model.cold_start_profile`, with the
    lapse rate and the tropopause height DERIVED rather than read, because that
    is where they come from.
    """
    cfg = _config(cfg)
    block = profile if profile is not None else cfg["model"]["cold_start_profile"]
    surface = float(block["surface_temperature_k"])
    ztrs = setzt_profile(surface, cold_start_lapse_rate_k_per_m(cfg),
                         cold_start_tropopause_height_m(surface, cfg), cfg)
    _sigma, _sigmah, dsigma = sigma_levels(cfg)
    return float((ztrs * dsigma).sum() / dsigma.sum())


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
