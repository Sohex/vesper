#!/usr/bin/env python3
"""The wet endmember of every land surface class, per band, under this star.

    python analysis/soil_albedo_wetting.py

A worldbuilding project. Everything below is about the simulated planet Vesper:
its land surface classes, the reflectance spectra that stand in for them, and
the albedo boundary condition `exoplasim/scripts/build_surface_albedo.py` stages
into the ExoPlaSim climate column.

`analysis/rock_albedo_bands.py` gives every class a DRY band pair. This is the
other end of the same axis: the pair a fully wetted surface of the same material
would have, in the same two bands, so that the modelled soil albedo can move
between them as the modelled land column wets and dries. PHYS-15 is the row.

## Why this is derivable at all, and what changed

Before the sources below were held, this project had one wetting measurement --
Penndorf (1956) Table 1, luminous reflectance over 0.38-0.77 um -- which brackets
band 1 and says nothing about band 2, where liquid water absorbs and where the
larger share of this star's flux arrives. A wet endmember for the band carrying
most of the flux cannot be guessed and must not be carried over from an Earth
calibration.

The two papers close it WITHOUT a new measurement, because both give the wet
reflectance as a closed-form function of the DRY reflectance AT THE SAME
WAVELENGTH. Applied wavelength by wavelength to the dry spectra this project
already integrates, and then band-integrated against this star, they produce a
band-2 wet endmember from band-2 dry data. The bands come out different because
the dry spectra are different in the two bands, not because a band-2 ratio was
asserted.

  Lekner, Dorf (1988). The Angstrom mechanism, extended. A water film over a
  rough surface reflects light back down onto it by total internal reflection at
  the film's upper face, so the surface gets more chances to absorb. Closed
  form, no fitted constant: the probability of internal reflection p follows
  from the film's refractive index, and the change in the surface's own
  absorptance follows from the drop in relative index at the substrate.

  Twomey, Bohren, Mergenthaler (1986). The other mechanism. In a finely divided
  medium, replacing interstitial air with water lowers the particle-to-medium
  index ratio, which makes single scattering far more forward-peaked; a photon
  then needs more scatterings to escape and is more likely to be absorbed on the
  way. Similarity scaling turns that into a reduced effective single-scattering
  albedo, and Chandrasekhar's semi-infinite isotropic result turns that back
  into a reflectance.

Lekner and Dorf say plainly which medium each is for: theirs "would seem to
apply best to rough solid surfaces", TBM's "to finely divided media, such as
sand". This world's land is both -- outcrop and regolith -- and neither paper
claims to cover the other's case, so the two are carried as the two arms of a
BRACKET rather than one being chosen. That is also what makes the bracket
testable: it has to contain the wet/dry pairs this project already holds.

## The trap this script exists to avoid: the map is not linear in the level

Both maps are strongly nonlinear in the dry reflectance, and both fix the
endpoints -- a perfect absorber and a perfect reflector are unchanged by
wetting, and the darkening is largest in between. So the wet/dry RATIO is not a
property of the material alone: it depends on where the material's dry
reflectance sits. A laboratory powder at 0.6 and the outcrop it came from at
0.15 have different wetting ratios, and `analysis/rock_albedo.py` measures that
powder-to-slab offset at 2.2 to 5.1 on the same rock.

Ratios transfer and levels do not, which is the rule `rock_albedo_bands.py` and
`playa_albedo.py` are both built on -- but here the ratio being derived depends
on the level, so the level has to be right before the map is applied. Each
spectrum is therefore rescaled to the class's own broadband albedo before
wetting: the SHAPE from the library, the LEVEL from `lithology.js` and
`config/planet.yaml`'s overrides, exactly the division `rock_albedo_bands.py`
makes. Applying the map to a raw laboratory powder spectrum would report the
wetting of a powder.

## What is NOT derived here

**The liquid's own absorption.** Both mechanisms treat the water as
non-absorbing and describe only what it does to the geometry of scattering.
Liquid water absorbs weakly through band 1 and strongly in parts of band 2 --
`exoplasim/data/water/hale_querry_1973_liquid_water.dat` carries the constants --
so both arms are UPPER bounds on the band-2 wet endmember and the bracket is
open at its dark end there. `band2_liquid_absorption` reports how far, as the
band-2 flux-weighted absorptance of the skin's own equivalent water depth, which
is a magnitude and not a correction: turning it into one needs an absolute
scattering coefficient for the dry soil, which nothing here carries.

**The shape between the endmembers.** That is Sadeghi, Jones and Philpot (2015),
which is linear in the Kubelka-Munk transformed reflectance rather than in the
reflectance itself, and it is applied where the mixing happens rather than here.
This script derives the two ENDS.

**The level of the dry endmember**, which stays exactly what `lithology.js` and
the config overrides say, and **the sign on salt crust**, which is disputed
between the in-situ pairs and a twenty-year MODIS series and is carried as an
explicit refusal by the consumer rather than resolved here.

## The checks that can fail, all fixed before the first run

1. Both implementations reproduce their own papers' printed numbers: Lekner and
   Dorf's average interface reflectance, internal-reflection probability and
   wet-to-dry absorptance ratios at three substrate indices; TBM's worked
   example carried end to end.
2. The two arms must BRACKET every wet/dry pair this project holds. Seven pairs
   are held and one exception is expected and named.
3. The band-1 arms must bracket Penndorf's three band-1 ratios, which is the one
   band where a measurement exists to check against.
4. Both maps must DARKEN every class in both bands, neither may return a
   negative albedo, and the film arm may not fall below the reflectance of its
   own upper face. The last is a floor Lekner and Dorf's algebra guarantees, so
   violating it is an implementation error rather than a surprise about rock.

There is no check that the wet endmember stays above open water's albedo, and
the reason is worth stating because the dry side does carry one. The two arms
describe two different surfaces. A continuous film has a specular upper face and
cannot reflect less than that face does; interstitial water in a granular medium
has no continuous interface and legitimately goes darker, which is what the
interstitial arm does on this world's dark basalt classes.
`interstitial_arm_below_open_water` reports where, so the difference is visible.
"""

from __future__ import annotations

import argparse
import datetime
import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "analysis"))
from lib.paths import rel  # noqa: E402
from lib import stellar  # noqa: E402
import rock_albedo_bands as rab  # noqa: E402

OUTPUT = ROOT / "analysis" / "soil_albedo_wetting.json"
WATER_OPTICS = ROOT / "exoplasim" / "data" / "water" / "hale_querry_1973_liquid_water.dat"

# The real refractive index of the substrate, which is the only material input
# the Lekner and Dorf map takes besides the dry reflectance itself. 1.55 is the
# rock-forming silicate range -- quartz 1.544, alkali feldspar 1.52 to 1.53,
# plagioclase 1.53 to 1.59 -- and halite, the mineral `evaporite` is named for,
# is 1.544. The bracket is the range Lekner and Dorf state their own figures
# over, "most minerals have refractive indices within this range". The map's
# sensitivity to it is REPORTED rather than assumed small: see
# `substrate_index_sensitivity` in the output.
SUBSTRATE_INDEX = 1.55
SUBSTRATE_INDEX_BRACKET = (1.5, 2.5)

# Liquid water's real refractive index below 0.75 um. The Hale and Querry
# extract this project holds starts at 0.75 um, and over 0.4 to 0.75 um the real
# index runs 1.339 to 1.330 in that table's parent, so a single value carries
# band 1 to better than the map can resolve. `nl_flat_sensitivity` measures what
# using this value in band 2 as well would have changed, so the wavelength
# dependence that IS carried has a size.
WATER_INDEX_BAND1 = 1.333

# The reflectance of the water film's own upper face at normal incidence, which
# is the floor the film arm cannot go below by its own algebra. Evaluated at the
# band-1 index because water's real index falls with wavelength and this is the
# larger of the two, so it is the binding floor over the whole working range.
FILM_FACE_REFLECTANCE = ((WATER_INDEX_BAND1 - 1.0) / (WATER_INDEX_BAND1 + 1.0)) ** 2

# Twomey, Bohren and Mergenthaler's scaling ratio, r = (1 - g_wet)/(1 - g_dry),
# the single number their similarity scaling needs. TAKEN FROM THEIR OWN WORKED
# EXAMPLE in section III rather than read off their Figure 1: for a particle of
# real refractive index 1.5, "typical of sand and many common minerals at
# visible wavelengths", the asymmetry parameter goes from about 0.83 in air to
# 0.96 in water and "giving for the ratio r in Eq. (5) the value 0.23".
#
# 0.235 is (1 - 0.96)/(1 - 0.83) to the precision the two asymmetry parameters
# are stated at, and reproduces their 0.23. It is NOT wavelength dependent: the
# asymmetry parameter of a particle much larger than the wavelength is set by
# the index ratio and not by the size, which is the paper's own argument for
# why the effective-size explanation fails.
TBM_SCALING_RATIO = 0.235
TBM_G_DRY = 0.83
TBM_G_WET = 0.96

# The illumination TBM's Figures 3 and 4 are drawn at, cos(41.4 degrees). The
# map is evaluated here at the same one so the reproduction check means
# something; `illumination_sensitivity` reports what the end-to-end map does at
# the two extremes instead.
TBM_MU0 = 0.75

# Every dry-against-wet pair this project holds, and the bracket has to contain
# them. Levels are Earth surfaces and do not transfer to this world -- what is
# being tested is the MAP, which takes a dry level to a wet one, so the test is
# performed at each source's own dry level. references/INDEX.md carries each row.
HELD_PAIRS = (
    ("Penndorf 1956 clay soil", 0.150, 0.075, "band1", True),
    ("Penndorf 1956 sand", 0.310, 0.180, "band1", True),
    ("Penndorf 1956 bare rich soil", 0.072, 0.055, "band1", True),
    ("Angstrom 1925 sand, via Lekner Table III", 0.182, 0.091, "broadband", True),
    ("Angstrom 1925 black mold, via Lekner Table III", 0.141, 0.084, "broadband", True),
    ("Craft, Horel 2019 Bonneville halite crust", 0.450, 0.220, "broadband", True),
    # The named exception, and it is expected to fail rather than tolerated.
    # Malek's surface is a thin halite crust over SHALLOW BRINE and its wet
    # value is the crust with brine at the surface, which is inundation and not
    # a wetted skin; its dry value is stated only as "above 0.75", so the pair
    # is a bound in one coordinate and a different process in the other.
    ("Malek et al. 1990 Pilot Valley playa", 0.750, 0.240, "broadband", False),
)

# Open water's albedo in the export's rock table, carried here only to REPORT
# which classes wet to below it. It is not a floor and the script does not
# refuse on it, because the two arms describe two different surfaces: Lekner and
# Dorf's is a continuous water FILM, which has a specular upper face and cannot
# reflect less than that face does, and TBM's is INTERSTITIAL water in a
# granular medium, which has no continuous interface at all and legitimately
# goes darker. Open water's 0.06 is the reflectance of a smooth interface over
# an effectively infinite absorber, so it bounds the film arm and says nothing
# about the interstitial one. The consumer's own ceiling -- substrate albedo
# minus open water's -- is about full INUNDATION and is a different quantity.
OPEN_WATER_ALBEDO = 0.06

# Gauss-Legendre order for the Chandrasekhar H-function quadrature. Raised until
# the reproduction check stopped moving in its fifth figure.
NMU = 128


# --------------------------------------------------------------------------
# Fresnel averages. Lekner and Dorf's Eq (7) is defined as an integral, and it
# is evaluated here as one rather than through their Eq (8)'s closed form: the
# closed form is a long rational expression that is easy to transcribe wrongly
# and impossible to check by eye, and the integral is the definition. The
# reproduction check is what says the two agree.
# --------------------------------------------------------------------------

def fresnel_unpolarised(theta_i: np.ndarray, n: float) -> np.ndarray:
    """Reflectance at a plane interface, medium 1 to medium 2, n = n2/n1."""
    sin2_t = (np.sin(theta_i) / n) ** 2
    total = sin2_t >= 1.0
    cos_t = np.sqrt(np.clip(1.0 - sin2_t, 0.0, None))
    cos_i = np.cos(theta_i)
    r_s = ((cos_i - n * cos_t) / (cos_i + n * cos_t)) ** 2
    r_p = ((cos_t - n * cos_i) / (cos_t + n * cos_i)) ** 2
    return np.where(total, 1.0, 0.5 * (r_s + r_p))


def mean_reflectance_exact(n: float, samples: int = 200001) -> float:
    """Lekner and Dorf Eq (7): reflectance of an isotropically illuminated face."""
    theta = np.linspace(0.0, np.pi / 2, samples)
    weight = np.sin(theta) * np.cos(theta)
    return float(np.trapezoid(fresnel_unpolarised(theta, n) * weight, theta)
                 / np.trapezoid(weight, theta))


# Eq (7) is needed at one index per WAVELENGTH once the water's dispersion is
# carried, and the integral above is too slow to evaluate thousands of times.
# It is smooth and monotone in n, so it is tabulated once over the whole range
# either map can ask for and interpolated. `mean_reflectance_table_error` in the
# output is the interpolation's own residual against the exact integral, so the
# shortcut has a measured cost rather than an assumed one.
_RBAR_N = np.linspace(0.6, 3.2, 1301)
_RBAR_V = np.array([mean_reflectance_exact(float(n), 20001) for n in _RBAR_N])


def mean_reflectance(n: np.ndarray | float) -> np.ndarray | float:
    """Eq (7), tabulated. Exact within `mean_reflectance_table_error`."""
    return np.interp(n, _RBAR_N, _RBAR_V)


def internal_reflection_probability(n_liquid):
    """Lekner and Dorf Eq (9). Angstrom's Eq (2) is 1 - 1/n**2, which is lower."""
    return 1.0 - (1.0 - mean_reflectance(n_liquid)) / np.asarray(n_liquid) ** 2


def lekner_wet(alpha_dry: np.ndarray, n_liquid: np.ndarray | float,
               n_substrate: float = SUBSTRATE_INDEX,
               normal_illumination: bool = True) -> np.ndarray:
    """Wet albedo from dry albedo, Lekner and Dorf (1988) Eqs (1), (9) and (11).

    `alpha_dry` is 1 - a_d in the paper's notation; `n_liquid` may be per
    wavelength, which is the whole reason this closes band 2.

    R1, the reflectance at the film's upper face, is the normal-incidence value
    the paper's own figures are drawn at, and the departure from it is offered
    rather than taken. Eq (7)'s isotropic average would be the consistent choice
    for a flux albedo -- it is the reason the other arm is evaluated on TBM's
    Eq (4) rather than its Eq (3) -- but it puts the film's own face at 0.066
    and therefore BRIGHTENS every substrate darker than that, which is a
    branch neither paper carries and which this model's skin cannot reach: a
    continuous surface film is Sadeghi's oversaturation regime, and the
    modelled skin caps at field capacity. `illumination_sensitivity` reports
    both.
    """
    a_dry = 1.0 - np.asarray(alpha_dry, dtype=float)
    n_l = np.broadcast_to(np.asarray(n_liquid, dtype=float), a_dry.shape)
    # Eqs (10) and (11): the wet surface's own absorptance, interpolated between
    # the weak- and strong-absorption limits with weights 1 - a_d and a_d.
    small = ((1.0 - mean_reflectance(n_substrate / n_l))
             / (1.0 - mean_reflectance(n_substrate)))
    a_wet = a_dry * ((1.0 - a_dry) * small + a_dry)
    p = internal_reflection_probability(n_l)
    r1 = (((n_l - 1.0) / (n_l + 1.0)) ** 2 if normal_illumination
          else mean_reflectance(n_l))
    absorbed = (1.0 - r1) * a_wet / (1.0 - p * (1.0 - a_wet))
    return 1.0 - absorbed


# --------------------------------------------------------------------------
# Twomey, Bohren and Mergenthaler. Chandrasekhar's semi-infinite isotropic
# layer both ways: albedo from scaled single-scattering albedo, and back.
# --------------------------------------------------------------------------

_MU, _WMU = np.polynomial.legendre.leggauss(NMU)
_MU = 0.5 * (_MU + 1.0)
_WMU = 0.5 * _WMU


def h_function(omega: float, mu0: float, iterations: int = 500) -> float:
    """Chandrasekhar's H at `mu0` for isotropic scattering of albedo `omega`."""
    h = np.ones(NMU)
    for _ in range(iterations):
        kernel = ((omega / 2.0) * h * _WMU)[None, :] / (_MU[:, None] + _MU[None, :])
        nxt = 1.0 / (1.0 - _MU * kernel.sum(axis=1))
        if np.max(np.abs(nxt - h)) < 1.0e-14:
            h = nxt
            break
        h = nxt
    return float(np.interp(mu0, _MU, h))


def albedo_from_omega(omega: np.ndarray, mu0: float = TBM_MU0) -> np.ndarray:
    """TBM Eq (4). Vectorised over `omega` by evaluating H per value."""
    vals = np.atleast_1d(np.asarray(omega, dtype=float))
    out = np.array([1.0 - np.sqrt(max(0.0, 1.0 - w)) * h_function(float(w), mu0)
                    for w in vals])
    return out.reshape(np.shape(omega)) if np.shape(omega) else float(out[0])


def _omega_table(mu0: float, points: int = 257):
    """A monotone (omega', albedo) table, so the inversion is an interpolation."""
    omega = 1.0 - np.geomspace(1.0e-9, 1.0, points)[::-1]
    omega = np.clip(omega, 0.0, 1.0 - 1.0e-12)
    alb = np.array([1.0 - np.sqrt(1.0 - w) * h_function(float(w), mu0)
                    for w in omega])
    order = np.argsort(alb)
    return alb[order], omega[order]


def tbm_wet(alpha_dry: np.ndarray, ratio: float = TBM_SCALING_RATIO,
            mu0: float = TBM_MU0, table=None) -> np.ndarray:
    """Wet albedo from dry albedo, TBM (1986) Eqs (2), (4) and (5)."""
    alb, omega = table if table is not None else _omega_table(mu0)
    dry = np.asarray(alpha_dry, dtype=float)
    w_dry = np.interp(np.clip(dry, alb[0], alb[-1]), alb, omega)
    w_wet = ratio * w_dry / (ratio * w_dry + 1.0 - w_dry)
    return np.interp(w_wet, omega, alb)


# --------------------------------------------------------------------------
# Spectra
# --------------------------------------------------------------------------

def water_indices(wavelength_um: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(n, k) of liquid water on `wavelength_um`, Hale and Querry (1973).

    The extract this project holds starts at 0.75 um, so band 1 takes the flat
    band-1 value. Both are returned because n drives the two wetting maps and k
    is what neither of them carries.
    """
    table = np.loadtxt(WATER_OPTICS)
    lam, k_w, n_w = table[:, 0], table[:, 1], table[:, 2]
    n = np.where(wavelength_um < lam[0], WATER_INDEX_BAND1,
                 np.interp(wavelength_um, lam, n_w))
    k = np.where(wavelength_um < lam[0], np.interp(lam[0], lam, k_w),
                 np.interp(wavelength_um, lam, k_w))
    return n, k


def band_means(grid: np.ndarray, refl: np.ndarray, flux: np.ndarray):
    """(band1, band2, broadband) flux-weighted means of `refl` on `grid`."""
    lo = grid <= rab.BAND_SPLIT_UM

    def mean(sel):
        return float(np.trapezoid(refl[sel] * flux[sel], grid[sel])
                     / np.trapezoid(flux[sel], grid[sel]))

    return mean(lo), mean(~lo), mean(np.ones_like(grid, dtype=bool))


def qualifying_spectra(code: str, spec: dict):
    """Every (name, wavelength, reflectance) this class's selectors admit.

    The selection is `rock_albedo_bands.py`'s, imported rather than restated:
    the same taxa, the same sample-level mineralogy rules, the same
    exclusions. A second copy of that table would be a second opinion about
    which rock a class is.
    """
    out = []
    for pattern in spec.get("eco", ()):
        for path in sorted(glob.glob(str(rab.ECOSTRESS / (pattern + "spectrum.txt")))):
            head, arr = rab.read_ecostress(Path(path))
            if arr is None or rab.sample_verdict(code, head) is not None:
                continue
            out.append((head.get("Name", Path(path).name), arr[:, 0], arr[:, 1] / 100.0))
    for pattern, mineral in spec.get("eco_named", ()):
        for path in sorted(glob.glob(str(rab.ECOSTRESS / (pattern + "spectrum.txt")))):
            head, arr = rab.read_ecostress(Path(path))
            if arr is None:
                continue
            if not head.get("Name", "").lower().startswith(mineral.lower()):
                continue
            if rab.sample_verdict(code, head) is not None:
                continue
            out.append((head.get("Name", Path(path).name), arr[:, 0], arr[:, 1] / 100.0))
    for name in spec.get("poseidon", ()):
        path = rab.POSEIDON / name
        if path.is_file():
            table = np.loadtxt(path)
            order = np.argsort(table[:, 0])
            out.append((name, table[order, 0], table[order, 1]))
    return out


def paper_reproduction() -> dict:
    """Check 1: both implementations against their own papers' printed numbers."""
    n_l = 4.0 / 3.0
    got_rbar = mean_reflectance(n_l)
    got_p = internal_reflection_probability(n_l)
    small = {nr: (1.0 - mean_reflectance(nr / n_l)) / (1.0 - mean_reflectance(nr))
             for nr in (1.5, 2.0, 2.5)}
    table = _omega_table(TBM_MU0)
    tbm_example = float(tbm_wet(np.array([0.30]), table=table)[0])
    return {
        "lekner_mean_reflectance_water": {
            "computed": round(got_rbar, 5),
            "paper": 0.0667,
            "how": "implied by the paper's own p = 0.475 through its Eq (9)",
            "passes": bool(abs(got_rbar - 0.0667) < 5.0e-4),
        },
        "lekner_internal_reflection_probability": {
            "computed": round(got_p, 4), "paper": 0.475,
            "passes": bool(abs(got_p - 0.475) < 1.0e-3),
        },
        "angstrom_internal_reflection_probability": {
            "computed": round(1.0 - 1.0 / n_l ** 2, 4), "paper": 0.437,
            "passes": bool(abs((1.0 - 1.0 / n_l ** 2) - 0.437) < 1.0e-3),
        },
        "lekner_absorptance_ratio_small_absorption": {
            "computed": {str(k): round(v, 3) for k, v in small.items()},
            "paper": {"1.5": 1.07, "2.0": 1.08, "2.5": 1.10},
            "passes": all(abs(small[nr] - p) < 5.0e-3
                          for nr, p in ((1.5, 1.07), (2.0, 1.08), (2.5, 1.10))),
        },
        "mean_reflectance_table_error": {
            "computed": round(float(np.max([
                abs(float(mean_reflectance(n)) - mean_reflectance_exact(n))
                for n in (1.14, 1.20, 4.0 / 3.0, 1.45, 1.55, 2.0, 2.5)])), 8),
            "bound": 1.0e-5,
            "passes": bool(max(
                abs(float(mean_reflectance(n)) - mean_reflectance_exact(n))
                for n in (1.14, 1.20, 4.0 / 3.0, 1.45, 1.55, 2.0, 2.5)) < 1.0e-5),
            "note": "the tabulated Eq (7) against the exact integral, at the "
                    "indices the two maps actually ask for. The bound is one "
                    "part in 1e5, which is below the .sra write quantum the "
                    "staged field reaches the model through.",
        },
        "tbm_worked_example": {
            "dry": 0.30, "computed_wet": round(tbm_example, 3), "paper": 0.12,
            "passes": bool(abs(tbm_example - 0.12) < 0.01),
            "note": "the paper carries its example through Eq (3), the zenith "
                    "REFLECTANCE, and reads omega' = 0.825 off Figure 3's "
                    "reflectance curve. This implementation uses Eq (4), the "
                    "flux albedo, on both legs, because a flux albedo is what "
                    "the model reads; the end-to-end map differs from the "
                    "paper's printed answer by less than the width of its "
                    "figure. The two are not the same quantity and are not "
                    "required to agree exactly.",
        },
    }


def bracket_check(table) -> dict:
    """Check 2 and 3: the two arms against every held wet/dry pair."""
    rows = []
    for name, dry, wet, band, must in HELD_PAIRS:
        lek = float(lekner_wet(np.array([dry]), WATER_INDEX_BAND1)[0])
        tbm = float(tbm_wet(np.array([dry]), table=table)[0])
        lo, hi = min(lek, tbm), max(lek, tbm)
        inside = bool(lo <= wet <= hi)
        rows.append({
            "source": name, "quantity": band, "dry": dry, "wet_measured": wet,
            "lekner": round(lek, 4), "twomey": round(tbm, 4),
            "inside_bracket": inside, "required": must,
            "passes": bool(inside == must),
        })
    return {
        "pairs": rows,
        "passes": all(r["passes"] for r in rows),
        "band1_pairs_pass": all(r["passes"] for r in rows if r["quantity"] == "band1"),
        "what_it_tests": "the MAP, at each source's own dry level. The levels "
                         "are Earth surfaces and do not transfer; a map that "
                         "takes a dry level to a wet one is exactly what has to.",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--library", type=Path, default=rab.ECOSTRESS)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()
    if not args.library.is_dir():
        raise SystemExit(f"{args.library} is absent. It is bulk reference data, "
                         "excluded from git; unzip references/ecospeclib_*.zip "
                         "beside this path.")

    levels, overridden = rab.table_albedos()
    _low, hires = stellar.spectrum_paths()
    data = np.loadtxt(hires, skiprows=1)
    star_wl, star_flux = data[:, 0], data[:, 1]
    z1, z2 = (float(v) for v in stellar.band_fractions(hires))

    grid = star_wl[(star_wl >= rab.WORKING[0]) & (star_wl <= rab.WORKING[1])]
    flux = np.interp(grid, star_wl, star_flux)
    n_water, k_water = water_indices(grid)
    omega_table = _omega_table(TBM_MU0)

    classes: dict[str, dict] = {}
    for code, spec in rab.PROXIES.items():
        if code == "water":
            continue  # open water is not a surface that wets
        level = levels[code]
        rows = []
        for name, wl, refl in qualifying_spectra(code, spec):
            if wl.min() > rab.QUALIFY[0] or wl.max() < rab.QUALIFY[1]:
                continue
            r = np.interp(grid, wl, refl)
            _b1, _b2, broad = band_means(grid, r, flux)
            if broad <= 0.0:
                continue
            # THE LEVEL MOVE. Shape from the library, level from the table,
            # applied BEFORE the map because the map is not linear in the level.
            dry = np.clip(r * (level / broad), 0.0, 1.0)
            d1, d2, dbb = band_means(grid, dry, flux)
            lek = lekner_wet(dry, n_water)
            tbm = tbm_wet(dry, table=omega_table)
            l1, l2, lbb = band_means(grid, lek, flux)
            t1, t2, tbb = band_means(grid, tbm, flux)
            rows.append((d1, d2, dbb, l1, l2, lbb, t1, t2, tbb))
        if not rows:
            raise SystemExit(
                f"no qualifying spectrum for {code!r} after the mineralogy "
                "selection. The selector, not this script, is what to fix; a "
                "class with no wet endmember would be staged as never wetting.")
        arr = np.array(rows)
        med = np.median(arr, axis=0)
        d1, d2, dbb, l1, l2, lbb, t1, t2, tbb = (float(v) for v in med)
        # The two arms, per band, as a RATIO to the dry pair of the same
        # spectra. A ratio because the consumer holds the level and this
        # script must not move it: `build_surface_albedo.py` writes
        # level * shape * wetting.
        classes[code] = {
            "why": spec["why"],
            "spectra": len(rows),
            "dry_broadband_level": level,
            "level_source": ("config model.lithology_albedo_overrides"
                             if code in overridden
                             else "vendor/orogen/js/lithology.js ROCK_CLASSES"),
            "band1_wet_over_dry": [round(t1 / d1, 4), round(l1 / d1, 4)],
            "band2_wet_over_dry": [round(t2 / d2, 4), round(l2 / d2, 4)],
            "broadband_wet_over_dry": [round(tbb / dbb, 4), round(lbb / dbb, 4)],
            "band1_wet_albedo": [round(t1, 4), round(l1, 4)],
            "band2_wet_albedo": [round(t2, 4), round(l2, 4)],
            "band1_dry_albedo": round(d1, 4),
            "band2_dry_albedo": round(d2, 4),
            # The whole point of the exercise: the two bands wet by DIFFERENT
            # amounts, and they do so because the dry spectra differ between
            # the bands, not because a band-2 ratio was asserted.
            "band2_wets_less_than_band1": bool((t2 / d2) > (t1 / d1)),
            # THE CHECKS THAT CAN FAIL, per class. Both maps must DARKEN in
            # both bands and neither may return a negative albedo; and the
            # film arm may not fall below the reflectance of its own upper
            # face, which is a floor its algebra guarantees and therefore a
            # genuine implementation check rather than a physical assumption.
            "darkens": bool(t1 <= d1 and t2 <= d2 and l1 <= d1 and l2 <= d2),
            "positive": bool(min(t1, t2) > 0.0),
            "film_arm_above_its_own_interface": bool(
                min(l1, l2) >= FILM_FACE_REFLECTANCE - 1e-6),
            "interstitial_arm_below_open_water": bool(
                min(t1, t2) < OPEN_WATER_ALBEDO),
        }

    broken = [c for c, r in classes.items()
              if not (r["darkens"] and r["positive"]
                      and r["film_arm_above_its_own_interface"])]
    if broken:
        raise SystemExit(
            f"the wetting map failed its own sign checks on {sorted(broken)}. "
            "A map that brightens a surface on wetting, returns a negative "
            "albedo, or puts a water film below the reflectance of its own "
            "upper face is being evaluated outside the range it is derived on. "
            "Fix the map or the level before staging a field that asserts it.")

    # The band-2 term neither map carries: liquid water's own absorption. Stated
    # as the single-pass absorptance of the equivalent water depth a saturated
    # 0.02 m skin holds, flux-weighted over band 2 -- a MAGNITUDE, not a
    # correction. Turning it into one needs an absolute scattering coefficient
    # for the dry soil, which nothing in this tree carries; the measurement that
    # would sidestep it is a wet-and-dry spectrum pair on the same sample, and
    # references/INDEX.md records which one is wanted.
    hi = grid > rab.BAND_SPLIT_UM
    skin_water_m = 0.02 * 0.45          # skin thickness times a full pore volume
    alpha_w = 4.0 * np.pi * k_water / (grid * 1.0e-6)
    absorptance = 1.0 - np.exp(-alpha_w * skin_water_m)
    band2_absorptance = float(np.trapezoid(absorptance[hi] * flux[hi], grid[hi])
                              / np.trapezoid(flux[hi], grid[hi]))
    band1_absorptance = float(np.trapezoid(absorptance[~hi] * flux[~hi], grid[~hi])
                              / np.trapezoid(flux[~hi], grid[~hi]))

    # What the two declared inputs are worth, measured rather than asserted.
    probe = np.array([0.10, 0.20, 0.30, 0.50])
    idx_spread = {
        f"{nr}": [round(float(v), 4)
                  for v in lekner_wet(probe, WATER_INDEX_BAND1, nr)]
        for nr in (SUBSTRATE_INDEX, *SUBSTRATE_INDEX_BRACKET)}
    mu_spread = {
        f"{mu}": [round(float(v), 4)
                  for v in tbm_wet(probe, mu0=mu, table=_omega_table(mu))]
        for mu in (0.5, TBM_MU0, 1.0)}
    lekner_isotropic = [round(float(v), 4)
                        for v in lekner_wet(probe, WATER_INDEX_BAND1,
                                            normal_illumination=False)]

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": rel(Path(__file__)),
        "spectrum": rel(hires),
        "library": rel(rab.ECOSTRESS),
        "water_optics": rel(WATER_OPTICS),
        "band_split_um": rab.BAND_SPLIT_UM,
        "model_band_flux_fractions": [round(z1, 6), round(z2, 6)],
        "working_range_um": list(rab.WORKING),
        "arms": ["twomey", "lekner"],
        "arms_note":
            "the two ARMS of one bracket and never an average. Lekner and Dorf "
            "give the water-film mechanism and say it fits rough solid "
            "surfaces; Twomey, Bohren and Mergenthaler give the interstitial "
            "forward-scattering mechanism and it fits finely divided media. "
            "This world's land is both. Twomey is the darker arm everywhere "
            "the bracket was evaluated, so it is listed first.",
        "substrate_refractive_index": SUBSTRATE_INDEX,
        "substrate_refractive_index_bracket": list(SUBSTRATE_INDEX_BRACKET),
        "twomey_scaling_ratio": TBM_SCALING_RATIO,
        "twomey_asymmetry_pair": [TBM_G_DRY, TBM_G_WET],
        "twomey_mu0": TBM_MU0,
        "level_move":
            "each spectrum is rescaled to the class's own broadband albedo "
            "before the map is applied, because both maps are nonlinear in the "
            "dry level and a laboratory powder sits 2.2 to 5.1 times above the "
            "outcrop it came from",
        "paper_reproduction": paper_reproduction(),
        "held_pair_bracket": bracket_check(omega_table),
        "substrate_index_sensitivity": {
            "dry_probe": [float(v) for v in probe],
            "wet_by_index": idx_spread,
            "note": "Lekner and Dorf's map over the whole index range they "
                    "state their figures for. The spread is the reason the "
                    "index is declared once rather than sourced per class.",
        },
        "illumination_sensitivity": {
            "dry_probe": [float(v) for v in probe],
            "wet_by_mu0": mu_spread,
            "lekner_film_face_normal_incidence": [
                round(float(v), 4)
                for v in lekner_wet(probe, WATER_INDEX_BAND1)],
            "lekner_film_face_isotropic": lekner_isotropic,
            "film_face_reflectance_normal": round(FILM_FACE_REFLECTANCE, 4),
            "film_face_reflectance_isotropic": round(
                float(mean_reflectance(WATER_INDEX_BAND1)), 4),
            "note": "the TBM arm at three illuminations, and the film arm with "
                    "its upper face at normal incidence, which is adopted, "
                    "against Eq (7)'s isotropic average, which is not. The "
                    "isotropic face is 0.066 and therefore brightens any "
                    "substrate darker than that, which on this world is every "
                    "basalt class; a continuous film over a dark absorber is "
                    "Sadeghi's oversaturation regime and not a wetted skin.",
        },
        "interstitial_arm_below_open_water": sorted(
            c for c, r in classes.items()
            if r["interstitial_arm_below_open_water"]),
        "band2_liquid_absorption": {
            "skin_equivalent_water_m": skin_water_m,
            "band1_single_pass_absorptance": round(band1_absorptance, 4),
            "band2_single_pass_absorptance": round(band2_absorptance, 4),
            "direction": "both arms are UPPER bounds on the wet endmember in "
                         "band 2: neither carries the liquid's own absorption, "
                         "and a wetted surface can only be darker for it",
            "what_would_close_it": "an absolute scattering coefficient for the "
                                   "dry soil, or a wet-and-dry reflectance "
                                   "spectrum pair measured on one sample",
        },
        "classes": classes,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{'class':22}{'dry b1':>8}{'dry b2':>8}"
          f"{'wet/dry b1':>22}{'wet/dry b2':>22}")
    for code, row in classes.items():
        b1 = row["band1_wet_over_dry"]
        b2 = row["band2_wet_over_dry"]
        print(f"{code[:21]:22}{row['band1_dry_albedo']:8.3f}"
              f"{row['band2_dry_albedo']:8.3f}"
              f"{f'{b1[0]:.3f} to {b1[1]:.3f}':>22}"
              f"{f'{b2[0]:.3f} to {b2[1]:.3f}':>22}")
    rep = report["paper_reproduction"]
    print("\npaper reproduction:",
          "all pass" if all(v["passes"] for v in rep.values()) else "FAILED")
    bc = report["held_pair_bracket"]
    print("held wet/dry pairs:",
          "as expected" if bc["passes"] else "NOT as expected")
    for r in bc["pairs"]:
        mark = "in" if r["inside_bracket"] else "OUT"
        print(f"  {r['source'][:44]:46}{r['dry']:6.3f}{r['wet_measured']:8.3f}"
              f"  [{r['twomey']:.3f}, {r['lekner']:.3f}] {mark}")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
