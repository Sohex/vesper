"""Re-weight the Lacis and Hansen shortwave gas absorptances for this star.

WHAT THE SCHEME ACTUALLY SAYS
-----------------------------
`radmod.f90`'s shortwave clear-sky code is Lacis and Hansen (1974). Their water
vapour absorptance, their Eq. 21,

    A_wv(y) = 2.9 y / ((1 + 141.5 y)^0.635 + 5.925 y)

is a fit to Yamamoto (1962), and Yamamoto states the definition outright: "the
definition of absorptivity is given by the ratio to the solar constant of the
energy absorbed by the entire vertical air column for normal incidence". It is a
fraction of TOTAL INCIDENT FLUX, and the Sun's spectrum is inside it -- Yamamoto
built it by weighting laboratory band absorptivities with the solar flux and
summing. Lacis and Hansen say the same thing a second way in their Section 5a:
"approximately 35% of the solar flux is contained in the regions of significant
water vapor absorption", which is 1 - p(k_1) = 1 - 0.6470 = 0.3530 from their own
Table 1, and appears in their Eq. 39 as the literal constant 0.353.

`radmod.f90` divides A_wv by `zsolar2`, the star's own flux share above 0.75 um,
and then applies it to band-2 flux. The two cancel exactly: the absorbed flux is
A_wv times the TOTAL incident flux whatever the star is. So the model currently
gives a K dwarf the Sun's absorbed fraction, which is the defect. Ozone is
already corrected this way through `o3uvw` and `o3visw`; this is the same
correction in the larger term.

HOW THE WEIGHT IS DERIVED
-------------------------
Yamamoto's construction is reversible. Howard, Burch and Williams (1956) measured
the band absorption of each near-infrared H2O and CO2 band as a total absorption
(an equivalent width) in cm-1, which is a property of the molecule and carries no
spectrum at all, and fitted it to

    weak band     int A_nu d_nu = c w^(1/2) (P + p)^k
    strong band   int A_nu d_nu = C + D log10 w + K log10 (P + p)

Howard's own Eq. 11 then defines the band-average fractional absorption
Abar_i = (int A_nu d_nu) / (nu_2 - nu_1), and says in as many words that "if the
spectral distribution of the radiation from a given source is known, the fraction
of the total radiation absorbed ... can be computed". So

    A_total(w) = sum_i f_i Abar_i(w)

with f_i the fraction of the incident flux falling in band i. Put the Sun in and
this must reproduce Lacis and Hansen Eq. 21; put this star in and it gives the
absorptance this star should have. The weight is the ratio, and it is the number
the patch multiplies the 2.9 by.

THE TEST THAT CAN FAIL
----------------------
Reconstructing A_total(w) from Howard's bands and a solar spectrum has a right
answer that this project did not choose. It is NOT Eq. 21's own stated accuracy:
"fits Yamamoto's absorption curve within ~1% for 10^-2 ~< y ~< 10 cm" is a
residual between two curves, and this reconstruction is a fourth construction of
the quantity rather than a fifth fit to Eq. 21. The right answer is the ENVELOPE
of the three published determinations Lacis and Hansen plot together in their
Fig. 11 -- Yamamoto as their Eq. 21, Fowle as Eq. 22, Korb as Eq. 23 -- because
of what they say about those curves on p. 127: "the uncertainty in their absolute
value is as great as the differences among the three curves." A reconstruction
inside that envelope agrees with every published construction of the quantity;
one outside it disagrees with all three. The envelope is evaluated at each water
amount rather than as one scalar, because the spread is 1.91 at y = 0.01 and 1.04
near y = 5, and it is widened by the +/-3% Howard state for the band absorptions
the reconstruction is built from.

THAT COMPARISON IS RECORDED AND DOES NOT GATE, and the reason is a measurement
rather than a preference. HITRAN2020 correlated-k, computing the same defined
quantity on the same path, misses the same widened envelope at every water amount
from 0.1 to 10 cm and by the same 6 to 10 per cent, and that is a floor because
it carries no water vapour continuum either. A bar a modern line list fails is
not a bar on this reconstruction: what the envelope has established is that all
three published curves are LOW. It is reported rather than widened.

THE GATE IS THAT CORRELATED-K ANSWER ITSELF. `CORRK_RATIO_TO_EQ21`, from
`exoplasim/notes/corrk-cross-check.md`, is this file's own `ratio_to_eq21`
measured from a different absorption dataset -- 76 correlated-k bands on
HITRAN2020 against Howard's nine -- at the same homogeneous 760 mm Hg path. Same
numerator definition, same denominator, so the two must agree to within Howard's
+/-3%, which is the only stated accuracy either side carries, and the run raises
if they do not. That keeps two measurements of one quantity from drifting apart
in silence, which is what `h2o_sw_level` rests on.

THE TWO PART BELOW THE AMOUNTS THE MODEL EVALUATES, by 14.4% at 0.01
precipitable cm against 0.26% at the operating path, and the reason is on this
side: Howard's weak-band fit is a square-root law, which is the strong-line
regime, and every weak band leaves that regime as the path dries out. The
correlated-k band mean over Howard's own intervals is 1.02 of this
reconstruction's at 6.3 um and 0.02 at 0.81 um at 0.01 cm, monotone in how weak
the band is. Dropping the 0.72 and 0.81 um bands entirely removes less than a
fifth of the gap, so it is the extrapolation and not those two bands.
`--blue` re-measures all of it; the dry end is reported and does not gate,
because a T42 column spans roughly 0.3 to 5 cm.

The solar reference is a BT-Settl 5772 K model built through the same blend as
`build_stellar_spectrum.py` uses for the star, so grid, converter and any model
systematic divide out of the ratio. Two further checks that can fail: it must
give 0.353 above 0.9 um, which is Lacis and Hansen's Section 5a value, and 0.517
below 0.75 um, which is `radmod.f90`'s own solar reference partitioning.

AND THE SAME MACHINERY FOR CO2, WHICH IS A NEW ABSORBER RATHER THAN A WEIGHT
----------------------------------------------------------------------------
`swr` has ozone in band 1 and water vapour in band 2 and nothing else; CO2
appears only in `lwr`, from Sasamori (1968). Lacis and Hansen did not
parameterise shortwave CO2 either, so the port is faithful and the absorber is
simply missing. Howard's Table II covers CO2 as well as H2O, so the same
integration gives what the missing term is worth, and this script also FITS it,
in the Lacis and Hansen manner, to the closed form
`exoplasim/patches/exoplasim-3.4.2-co2-shortwave.patch` codes. The fit is
solar-weighted, so the patch carries the Sun's CO2 absorptance and `co2sww`
re-weights it for the host exactly as `h2osww` re-weights Eq. 21.

THE FIT MADE HERE IS NOT THE FIT THE MODEL RUNS, AND HAS NOT BEEN SINCE PHYS-10.
Howard's band set is the only absorption data this file has, so the fit it makes
is a Howard fit and always will be. `corrk_cross_check.py --fit` refitted the
same closed form to HITRAN2020 through the correlated-k tables and that is what
`radmod.f90` now carries, 7.2% weaker at this planet's CO2 path. Both are
reported: `closed_form_fit` is this file's, which is what the patch header codes
and what the cross-check compares against, and `closed_form_fit_in_radmod` is
read out of the model source. Every PREDICTION here is priced on the model's,
because a prediction about the other one is a prediction about no code.

Two things separate the CO2 half from the H2O half and both are in
`exoplasim/notes/shortwave-co2.md`. Every CO2 band shares its interval with
water vapour, so CO2 is charged only with what water vapour leaves it; and the
CO2 column follows the planet's gravity, which is the term that makes this
world's column smaller than Earth's at the same mixing ratio.

THE CO2 TEST THAT CAN FAIL
--------------------------
Earth's near-infrared CO2 solar absorption is measured at 1.5 to 2.5 W/m2, a
quantity this derivation did not choose. Running the same integration at EARTH's
column, EARTH's gravity and EARTH's mean insolation has to land inside it, and
`--verify` reports where it lands. That check is what licenses the number for
this star; without it the ratio would be two unvalidated integrals.

    python exoplasim/scripts/shortwave_band_weights.py
    python exoplasim/scripts/shortwave_band_weights.py --verify
    python exoplasim/scripts/shortwave_band_weights.py --blue

Writes `analysis/shortwave_band_weights.json`. Downloads are cached in the same
place `build_stellar_spectrum.py` caches its own; pass --refresh to refetch.
`--blue` writes nothing and is the only mode that opens the correlated-k bundle
`corrk_cross_check.py` reads.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, MODEL_SRC  # also puts lib/ on sys.path
from paths import climatology_path
import sensitivity  # noqa: E402  from lib/

SSAP = "http://svo2.cab.inta-csic.es/theory/newov2/ssap.php"
MODEL = "bt-settl"
LOGG = 4.5
METALLICITY = 0.0

# SVO record identifiers for the BT-Settl CIFIST2011 grid points, log g = 4.5,
# [M/H] = 0. The 4900 and 5000 K entries are the same fids
# `build_stellar_spectrum.py` pins, so the star built here is the star in
# `inputs/stellarspectra/k25v.dat` and not a second interpolation of it. Pinned
# so a server-side renumbering becomes a header mismatch rather than a different
# star.
FIDS = {4900: 3697, 5000: 3850, 5700: 4824, 5800: 4965}

# The solar reference. 5772 K is the effective temperature `radmod.f90:207`
# normalises its own Rayleigh cross-section to, and the temperature at which the
# scheme's zsolar1 = 0.517 partitioning is stated to hold.
SOLAR_TEFF = 5772.0

CACHE = Path(tempfile.gettempdir()) / "vesper-btsettl-cache"

# Howard, Burch and Williams (1956), Table II of papers II (CO2) and III (H2O).
# Band limits in cm-1; w in precipitable cm for H2O and atmos-cm for CO2; P and p
# in mm Hg; the fits return a total band absorption in cm-1. `transition` is the
# total absorption at which their weak fit gives way to their strong fit.
#
# The 6.3 um H2O band is in the list because Yamamoto's sum includes it and
# Lacis and Hansen fit Yamamoto; it carries almost no solar flux and rather more
# of a K dwarf's, which is exactly the kind of term this exercise exists to move.
H2O_BANDS = {
    #  um    lo_cm1  hi_cm1     c     k   C      D     K    transition
    "6.3": dict(lo=1150.0, hi=2050.0, c=356.0, k=0.30, C=302.0, D=218.0, K=157.0, transition=160.0),
    "3.2": dict(lo=2800.0, hi=3340.0, c=40.2, k=0.30, C=None, D=None, K=None, transition=500.0),
    "2.7": dict(lo=3340.0, hi=4400.0, c=316.0, k=0.32, C=337.0, D=246.0, K=150.0, transition=200.0),
    "1.87": dict(lo=4800.0, hi=5900.0, c=152.0, k=0.30, C=127.0, D=232.0, K=144.0, transition=275.0),
    "1.38": dict(lo=6500.0, hi=8000.0, c=163.0, k=0.30, C=202.0, D=460.0, K=198.0, transition=350.0),
    "1.1": dict(lo=8300.0, hi=9300.0, c=31.0, k=0.26, C=None, D=None, K=None, transition=200.0),
    "0.94": dict(lo=10100.0, hi=11500.0, c=38.0, k=0.27, C=None, D=None, K=None, transition=200.0),
}

# Howard et al. did not measure the 0.72 and 0.81 um bands. Yamamoto estimated
# them from Fowle's data over 13514-14286 and 11905-12658 cm-1 and says they
# "cannot be neglected ... because of the large solar energy in this region".
# They are carried here at the 0.94 um band's shape with a strength scaled down
# by the factor `WEAK_BLUE_SCALE`, and `--verify` reports what dropping them
# does. They matter to the weight only in the direction of making it SMALLER,
# because they sit where a K dwarf's flux boost is least, so leaving them out
# would flatter the correction.
H2O_BLUE_BANDS = {
    "0.81": dict(lo=11905.0, hi=12658.0, c=38.0, k=0.27, C=None, D=None, K=None, transition=200.0),
    "0.72": dict(lo=13514.0, hi=14286.0, c=38.0, k=0.27, C=None, D=None, K=None, transition=200.0),
}
WEAK_BLUE_SCALE = {"0.81": 0.30, "0.72": 0.10}

# THOSE TWO FACTORS ARE DECLARED WITH A MEASURED BRACKET, and the bracket is the
# correlated-k band-mean absorptance over Howard's own two intervals divided by
# what the 0.94 um shape puts there. Measured on 2026-08-26 against HITRAN2020
# through the LMD Generic PCM tables, the same absorption data BAR 4 of
# `exoplasim/notes/corrk-cross-check.md` gates on, at each water amount:
#
#     w, cm     0.01    0.03    0.1     0.3     1.0    2.7891   5.0    10.0
#     0.81 um  0.0206  0.0353  0.0624  0.0992  0.1453  0.1779  0.1919  0.2032
#     0.72 um  0.0109  0.0188  0.0336  0.0549  0.0850  0.1093  0.1206  0.1302
#
# so `WEAK_BLUE_MEASURED_RANGE` below is that over the 0.3 to 5 cm a T42 column
# spans. NO CONSTANT IS RIGHT FOR EITHER BAND: the factor runs by a factor of ten
# across the table because the 0.94 um band's SHAPE is what is wrong, not only
# its size. Howard's weak fit is a square-root law, which is the strong-line
# regime, and these two bands are unsaturated enough to go as w^0.83 and w^0.86.
#
# THE DECLARED VALUES STAY AS THEY ARE, and that is a statement about the gate
# rather than about the numbers. This reconstruction is the independent side of
# BAR 4: its evidential value is that it reaches the same defined quantity from
# Howard's laboratory data with no correlated-k input anywhere in it, and two of
# nine bands taking their strength from those tables would make the gate partly a
# comparison of the tables with themselves. The substitution is priced instead,
# by `--blue`, and it moves BOTH aggregates away from the correlated-k answer --
# `h2osww` 1.3456 to 1.3502 against its 1.3272, `ratio_to_eq21` +0.20% to -1.35%
# -- because the operating-path agreement is a cancellation of band-level
# disagreements running 0.70 to 1.24 and a partial substitution breaks the
# cancellation without touching the bands supplying the other half of it. The
# whole question is worth 0.0046 in `h2osww`, 0.034 K, and lands inside the
# bracket this file already reports.
WEAK_BLUE_MEASURED = {"0.81": 0.1779, "0.72": 0.1093}
WEAK_BLUE_MEASURED_RANGE = {"0.81": (0.0992, 0.1919), "0.72": (0.0549, 0.1206)}

CO2_BANDS = {
    "15": dict(lo=550.0, hi=800.0, c=3.16, k=0.44, C=-68.0, D=55.0, K=47.0, transition=50.0),
    "5.2": dict(lo=1870.0, hi=1980.0, c=0.024, k=0.40, C=None, D=None, K=None, transition=30.0),
    "4.8": dict(lo=1980.0, hi=2160.0, c=0.12, k=0.37, C=None, D=None, K=None, transition=60.0),
    "4.3": dict(lo=2160.0, hi=2500.0, c=None, k=None, C=27.5, D=34.0, K=31.5, transition=50.0),
    "2.7": dict(lo=3480.0, hi=3800.0, c=3.15, k=0.43, C=-137.0, D=77.0, K=68.0, transition=50.0),
    "2.0": dict(lo=4750.0, hi=5200.0, c=0.492, k=0.39, C=-536.0, D=138.0, K=114.0, transition=80.0),
    "1.6": dict(lo=6000.0, hi=6550.0, c=0.063, k=0.38, C=None, D=None, K=None, transition=80.0),
    "1.4": dict(lo=6650.0, hi=7250.0, c=0.048, k=0.41, C=None, D=None, K=None, transition=80.0),
}

# `lwr`'s own zmmair/zmmco2/zrco2 parameters in radmod.f90, reused rather than
# restated so the shortwave CO2 column and the longwave one are the same
# quantity. Cited by symbol: the line numbers moved by 37 under the fork's
# multi-species aerosol work and the old ones now land in the cloud branch.
ZMMAIR = 0.0289644          # molecular weight of air, kg/mol
ZMMCO2 = 0.0440098          # molecular weight of CO2, kg/mol
ZRCO2 = 1.9635              # CO2 density at STP, kg/m3

# The 2.7 um CO2 band lies under the strong 2.7 um H2O band, and Yamamoto (1962)
# drops it outright "because of overlapping by the strong 2.7 um H2O band". It is
# KEPT here and charged with only the fraction of its interval water vapour has
# left, which is the same correction every other CO2 band gets and is available
# because the overlap is computed band by band rather than argued about. Setting
# this False reproduces Yamamoto's choice, and the difference between the two is
# reported as the bracket rather than hidden inside one number.
CO2_KEEP_27UM = True

# The interval the closed-form CO2 fit is quoted over, in atmos-cm. The model
# never evaluates it below this: the thinnest sigma layer carries a few percent
# of the column and the smallest magnification is zbetta.
CO2_FIT_RANGE = (1.0, 1.0e4)

# Earth, for the check that can fail. Present-day gravity and surface pressure;
# the mixing ratio is deliberately the config's, so the comparison isolates the
# gravity and pressure terms rather than mixing in a different CO2 abundance.
EARTH_GRAVITY = 9.80665
EARTH_SURFACE_PRESSURE_PA = 101325.0
EARTH_MEAN_INSOLATION = 1361.0 / 4.0
# Earth's near-infrared CO2 solar absorption, the range the literature measures.
EARTH_CO2_SHORTWAVE_W_M2 = (1.5, 2.5)

# Standard pressure, in the mm Hg the Howard fits are written in. Lacis and
# Hansen Eq. 21 is stated for P0 = 1013 mb, T0 = 273 K, and the model reaches it
# by scaling the water amount rather than the pressure, so the reconstruction is
# evaluated here at standard pressure throughout.
P_STANDARD_MMHG = 760.0

# Lacis and Hansen (1974) Section 5a and Table 1: the fraction of the solar flux
# in the regions of significant water vapour absorption. Reached two ways in the
# paper, from the k-distribution as 1 - 0.6470 and from Joseph's (1971) cut at
# 0.9 um, and used as the literal constant in their Eq. 39.
LH74_SOLAR_ACTIVE_FRACTION = 0.353

# radmod.f90:125-126 and :207. The scheme's own solar partitioning at 5772 K.
RADMOD_ZSOLAR1 = 0.517


def lacis_hansen_h2o(y: np.ndarray | float) -> np.ndarray | float:
    """Lacis and Hansen (1974) Eq. 21, verbatim, as radmod.f90 codes it."""
    return 2.9 * y / ((1.0 + 141.5 * y) ** 0.635 + 5.925 * y)


def fowle_h2o(y: np.ndarray | float) -> np.ndarray | float:
    """Lacis and Hansen (1974) Eq. 22: Fowle (1915), the dotted curve of Fig. 11.

    Fowle's determination with a roughly 10% modification by Manabe and Moller
    (1961) for the 0.7 and 0.8 um bands. That 10% is a PHYSICAL adjustment for
    two bands and is not a fit residual, so nothing here reads it as one.
    """
    return 0.0946 * y ** 0.303


def korb_h2o(y: np.ndarray | float) -> np.ndarray | float:
    """Lacis and Hansen (1974) Eq. 23: Korb et al. (1956), the dashed curve.

    Curtis-Godson on Howard and Korb data, WITHOUT the weak near-infrared bands.
    The paper writes the left-hand side as log10 of TWICE the absorptance and
    the factor of 2 is inside the bracket, so the solved form carries the 0.5.
    Dropping it puts Korb at exactly double and still looks plausible on a plot:
    at y = 1 it would give 0.18 against the 0.09 the three curves converge on.
    The exponents are powers OF the logarithm, not logarithms of powers.
    """
    lg = np.log10(y)
    return 0.5 * 10.0 ** (-0.74 + 0.347 * lg - 0.056 * lg ** 2 - 0.006 * lg ** 3)


# THE BOUND ON THE RECONSTRUCTION, and what it is made of.
#
# The reconstruction here is a fourth construction of the same physical
# quantity: Howard's band absorptions, weighted by a solar spectrum and summed,
# which is Yamamoto's own construction run again. Eq. 21's stated accuracy
# CANNOT be the bar for it -- "fits Yamamoto's absorption curve within ~1% for
# 10^-2 ~< y ~< 10 cm" is a curve-reading residual between two of the curves,
# not a statement about the quantity. What Lacis and Hansen say about the
# quantity is on their p. 127, of the three curves they plot together: "Although
# the three curves in Fig. 11 are qualitatively similar, their differences are
# significant, particularly for small water vapor amounts. Moreover, the
# uncertainty in their absolute value is as great as the differences among the
# three curves."
#
# So the envelope of Eqs. 21, 22 and 23 is the bar, at each water amount rather
# than as one scalar over the decade: the spread is 1.91 at y = 0.01 and 1.04
# near y = 5, and a single number would be the dry end and nothing else. A
# reconstruction inside the envelope agrees with every published determination
# of this quantity; one outside it disagrees with all three.
#
# IT IS DECLARED AS AN INTER-FORMULA SPREAD AND NEVER AS AN ERROR BAR, which is
# the same sentence's other half: the authors put the absolute uncertainty at
# AS GREAT AS the spread, so the envelope is a floor on the expected
# disagreement and the bar here is conservative by construction. The three are
# also not independent -- Eqs. 21 and 23 both trace to Howard, Burch and
# Williams (1956) -- so the dry end of the spread is largely Fowle 1915 against
# Howard 1956.
#
# The envelope alone would be tighter at y = 5 than the data the reconstruction
# is built from, so it is widened by the accuracy Howard state for the band
# absorptions themselves, p. 244: "The expressions for the other bands give
# results which on the average agree with observed values of total absorption to
# within +/-3%." That is an average and it excludes the 0.94 um band, which is
# fitted on twelve runs and called approximate, so 3% is a floor here too.
#
# BOUNDING A HEATING RATE WOULD NEED A DIFFERENT BAND and this is not one. The
# ordering of dA/dy is not the ordering of A and swaps twice across the decade,
# and the derivative spread REOPENS to 1.51 at y = 10 where the absorptance
# spread has closed to 1.08. What this file produces is a ratio of absorptances,
# so the absorptance envelope is the right bar for it; anything derived from
# this file that consumes a heating rate needs its own.
HOWARD_BAND_ABSORPTION_ACCURACY = 0.03
# The only interval Lacis and Hansen state for Eq. 21, and the plotted range of
# all three curves. No range is printed for Eq. 22 or Eq. 23.
LH74_FIT_RANGE_CM = (0.01, 10.0)

# AND THE ENVELOPE IS NOT A GATE, because a modern line list misses it too.
#
# `exoplasim/notes/corrk-cross-check.md` computes the same defined quantity from
# HITRAN2020 correlated-k tables instead of Howard's bands, and it lands OUTSIDE
# the widened envelope at every water amount from 0.1 to 10 cm, by 6 to 10 per
# cent -- inside only at 0.01 cm, which is exactly where the reconstruction is
# inside as well. That is a floor, because the correlated-k side carries no water
# vapour continuum either and the continuum would push it further out.
#
# A bar that rejects a modern line-by-line calculation is not a bar on this
# reconstruction. What has been established is that all three published curves
# are LOW, not that this file is high. So the envelope comparison is kept and
# reported -- it is the record of how far every construction of this quantity,
# old and new, sits from what was published -- and it does not raise. It is NOT
# widened to admit the miss; that would be a criterion chosen after the run.
#
# THE GATE THAT REPLACES IT compares this file's reconstruction against that
# correlated-k answer at the path it was measured on, which is the same quantity
# by the same definition with only the absorption data differing. The two must
# not drift apart in silence: a change to the band set, to the spectra or to that
# note's table has to show up somewhere, and this is where.
CORRK_PATH_CM = 2.7891
# Six digits because `H2O_SW_LEVEL` below MULTIPLIES this and then rounds to the
# three decimals `config/planet.yaml` writes, so the fourth digit reaches the
# third of the key. `exoplasim/notes/corrk-cross-check.md` quotes it as 1.127
# because that is all the comparison it appears in can resolve; this is
# `corrk_cross_check.py`'s own absorptance at CORRK_PATH_CM over Eq. 21 at the
# same amount, and running that script reproduces it.
CORRK_RATIO_TO_EQ21 = 1.127592
# Howard's own +/-3% is the tolerance because it is the only stated accuracy
# either side of the comparison carries. Fixed before the comparison was made.
CORRK_AGREEMENT = HOWARD_BAND_ABSORPTION_ACCURACY

# THE DRY END, WHERE THE TWO PART, and which half of this reconstruction does it.
# The same correlated-k ratio at the two smallest amounts, declared from the same
# note. The two determinations agree to 1.5% from 0.1 to 10 cm and are 14.4%
# apart at 0.01, and the divergence is Howard's weak-band fit evaluated below the
# water amounts he measured: dropping the 0.72 and 0.81 um bands entirely removes
# 0.183 of the gap at 0.01 cm and 0.307 at 0.03, leaving 11.8% and 6.1%. The
# per-band evidence is monotone in band strength -- at 0.01 cm the correlated-k
# band mean over Howard's own intervals is 1.02 of the reconstruction's at 6.3 um
# and 0.02 at 0.81 um -- which is a wrongly extrapolated weak-band form and not a
# wrong scale on two bands.
#
# REPORTED AND NOT GATED, for the reason `world-njlb` states: a T42 column spans
# roughly 0.3 to 5 cm so the model never evaluates there, and `h2o_sw_weight` is
# a ratio in which a dry-end level error largely divides out. The gate is at
# CORRK_PATH_CM, which is where the number that is used comes from.
CORRK_RATIO_TO_EQ21_DRY = {0.01: 1.1529, 0.03: 1.0957}
# The share of the dry-end gap the blue pair carries, at those two amounts, as
# `--blue` measures it. Recorded so a change to the band set moves a reported
# number rather than nothing. One digit each: the correlated-k side's own
# temperature scan moves the 0.01 cm share over 0.15 to 0.23.
BLUE_SHARE_OF_DRY_GAP = {0.01: 0.18, 0.03: 0.31}

# THE WATER VAPOUR CONTINUUM NEITHER SIDE OF THAT RATIO CARRIES, as a fraction of
# the correlated-k absorptance. The correlated-k tables hold line centres to
# +/-25 cm-1 with the plinth removed, which is the MT_CKD definition, and Eq. 21
# is a fit to Yamamoto's sum over Howard's laboratory BAND absorptions, so
# neither has the window continuum. It absorbs BETWEEN the bands, so it adds to
# the line-by-line side and makes Eq. 21's deficit larger rather than smaller.
#
# DECLARED, from published numbers this project did not measure:
# `shine2012` p. 548 sizes the CAVIAR laboratory continuum against MT_CKD and
# p. 536 gives water vapour's share of clear-sky shortwave absorption;
# `mlawer2012` p. 2551 gives the laboratory-over-MT_CKD factor window by window,
# which is what turns an increment into a total. The bracket ends carry the
# transfer from a global mean to this homogeneous path as well.
# `exoplasim/notes/corrk-cross-check.md` does the arithmetic and cites each.
H2O_CONTINUUM_FRACTION = 0.031
H2O_CONTINUUM_FRACTION_BRACKET = (0.0015, 0.070)

# The level correction `config/planet.yaml` carries as `h2o_sw_level`. DERIVED
# here, from the two declared quantities above and nothing else -- no spectrum,
# no climatology, no network -- which is why `--level` can write its artifact on
# a tree that has neither. The bracket ends are arms to run rather than an error
# bar on a settled number.
H2O_SW_LEVEL = CORRK_RATIO_TO_EQ21 * (1.0 + H2O_CONTINUUM_FRACTION)
H2O_SW_LEVEL_BRACKET = tuple(
    CORRK_RATIO_TO_EQ21 * (1.0 + f) for f in H2O_CONTINUUM_FRACTION_BRACKET)


def h2o_sw_level_report() -> dict:
    """The `h2o_sw_level` derivation, its inputs and where each came from.

    Its own artifact rather than a node of the full report, because the full
    report cannot be written without a baseline climatology -- `predict()` opens
    one and there is deliberately no Earth fallback -- while this half needs
    nothing but the two declared constants above. Tying `config/planet.yaml`'s
    retyped copy to a node inside the climatology-dependent artifact would have
    made a stale file into a failing gate; this is the same tie without it.
    """
    return {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "exoplasim/scripts/shortwave_band_weights.py --level",
        "key": "model.h2o_sw_level",
        "value": H2O_SW_LEVEL,
        "bracket": list(H2O_SW_LEVEL_BRACKET),
        "value_rounded_for_config": round(H2O_SW_LEVEL, 3),
        "bracket_rounded_for_config": [round(x, 3) for x in H2O_SW_LEVEL_BRACKET],
        "inputs": {
            "correlated_k_over_eq21": {
                "value": CORRK_RATIO_TO_EQ21,
                "water_cm": CORRK_PATH_CM,
                "status": "DECLARED",
                "source": (
                    "exoplasim/notes/corrk-cross-check.md, HITRAN2020 through "
                    "the LMD Generic PCM correlated-k tables at a homogeneous "
                    "760 mm Hg path; reproduced by "
                    "exoplasim/scripts/corrk_cross_check.py, whose tables live "
                    "outside this repository. Gated inside this script against "
                    "the Howard reconstruction to +/-3%"),
            },
            "continuum_fraction": {
                "value": H2O_CONTINUUM_FRACTION,
                "bracket": list(H2O_CONTINUUM_FRACTION_BRACKET),
                "status": "DECLARED",
                "source": (
                    "shine2012 pp. 536 and 548, mlawer2012 p. 2551; the "
                    "MT_CKD-convention window continuum as a fraction of the "
                    "correlated-k absorptance, transferred to this path. "
                    "exoplasim/notes/corrk-cross-check.md"),
            },
        },
        "derivation": (
            "h2o_sw_level = correlated_k_over_eq21 * (1 + continuum_fraction), "
            "and the bracket is the same product at each end of the continuum "
            "fraction's bracket. Both factors are DECLARED, so the key is "
            "derived arithmetic on measurements this tree does not itself make; "
            "what this artifact closes is the retyping route into "
            "config/planet.yaml, not the measurements"),
        "what_it_scales": (
            "radmod.f90's h2oswl, a multiplicative level on the 2.9 of Lacis "
            "and Hansen Eq. 21 in swr (the ztwvtu and ztwv terms). The model "
            "default is 1.0, which is Eq. 21 unmodified"),
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def fetch(teff: int, refresh: bool) -> Path:
    """Download one BT-Settl grid point, or reuse the cached copy."""
    CACHE.mkdir(parents=True, exist_ok=True)
    target = CACHE / f"btsettl_{teff}_logg{LOGG}_m{METALLICITY}.txt"
    if target.is_file() and not refresh:
        return target
    url = f"{SSAP}?model={MODEL}&fid={FIDS[teff]}&format=ascii"
    print(f"fetching {teff} K from {url}")
    with urllib.request.urlopen(url, timeout=900) as response:
        target.write_bytes(response.read())
    return target


def read_btsettl(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """SVO's ASCII BT-Settl export: Angstrom against erg/cm2/s/A."""
    header: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if "=" in line:
            key, _, rest = line[1:].partition("=")
            value = rest.split("(")[0].strip().split()
            header[key.strip()] = value[0] if value else ""
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1], header


def blend(teff: float, lo_t: int, hi_t: int, refresh: bool) -> tuple[np.ndarray, np.ndarray, dict]:
    """Log-linear blend of two grid points, the same one build_stellar_spectrum uses."""
    sources = {}
    for t in (lo_t, hi_t):
        path = fetch(t, refresh)
        wave, flux, header = read_btsettl(path)
        if int(float(header.get("teff", -1))) != t:
            raise SystemExit(f"{path} declares teff={header.get('teff')}, expected {t}")
        sources[t] = (wave, flux, path)

    wave_hi, flux_hi, path_hi = sources[hi_t]
    wave_lo, flux_lo, path_lo = sources[lo_t]
    flux_lo_on_hi = np.interp(wave_hi, wave_lo, flux_lo)

    fraction = (teff - lo_t) / (hi_t - lo_t)
    floor = 1e-300
    log_blend = (1.0 - fraction) * np.log(np.maximum(flux_lo_on_hi, floor)) + fraction * np.log(
        np.maximum(flux_hi, floor)
    )
    flux = np.exp(log_blend)
    flux[(flux_lo_on_hi <= floor) | (flux_hi <= floor)] = 0.0

    meta = {
        "teff_k": teff,
        "interpolation_fraction": fraction,
        "endpoints": {
            str(lo_t): {"fid": FIDS[lo_t], "sha256": sha256(path_lo)},
            str(hi_t): {"fid": FIDS[hi_t], "sha256": sha256(path_hi)},
        },
    }
    return wave_hi, flux, meta


def planck(wave_angstrom: np.ndarray, teff: float) -> np.ndarray:
    """Planck's law per unit wavelength, on a wavelength grid in Angstrom.

    Kept as a reference: `solarini` builds a Planck curve only when no
    spectrum is staged, and `run_exoplasim.py` refuses that configuration
    (`radiation.stellar_spectrum: k25v` is the configured state; PHYS-2 and
    SPEC-2 are archived). The blackbody weight and the spectrum weight are
    different numbers, so both are computed.
    """
    h = 6.62607015e-34
    c = 2.99792458e8
    k = 1.380649e-23
    lam = np.clip(wave_angstrom, 1e-3, None) * 1e-10
    x = h * c / (lam * k * teff)
    # The absolute scale is irrelevant: everything downstream is a flux fraction.
    return np.where(x < 700.0, 1.0 / lam**5 / np.expm1(np.minimum(x, 700.0)), 0.0)


class Spectrum:
    """A stellar spectrum reduced to what a band weighting needs of it."""

    def __init__(self, wave_angstrom: np.ndarray, flux: np.ndarray, label: str, meta: dict):
        order = np.argsort(wave_angstrom)
        self.wave = wave_angstrom[order]
        self.flux = flux[order]
        self.label = label
        self.meta = meta
        self.cumulative = np.concatenate(
            [[0.0], np.cumsum(0.5 * (self.flux[1:] + self.flux[:-1]) * np.diff(self.wave))]
        )
        self.total = float(self.cumulative[-1])

    def fraction_between_um(self, lo_um: float, hi_um: float) -> float:
        lo, hi = lo_um * 1e4, hi_um * 1e4
        return float(
            (np.interp(hi, self.wave, self.cumulative) - np.interp(lo, self.wave, self.cumulative))
            / self.total
        )

    def fraction_in_band(self, lo_cm1: float, hi_cm1: float) -> float:
        """Flux fraction between two wavenumbers, given in cm-1."""
        return self.fraction_between_um(1e4 / hi_cm1, 1e4 / lo_cm1)


def band_absorption(band: dict, w: float, pressure_mmhg: float = P_STANDARD_MMHG) -> float:
    """Howard, Burch and Williams total band absorption, in cm-1.

    Their weak fit below the tabulated transition and their strong fit above it,
    which is how they say to use them. Bands with no strong fit stay on the weak
    one, as their Table II intends: the 3.2 um band is called out in the text as
    one where "the weak-band relation was satisfactory for all data obtained".
    The result is capped at the band width, since a band cannot absorb more than
    all of itself.
    """
    width = band["hi"] - band["lo"]
    weak = None
    if band["c"] is not None:
        weak = band["c"] * math.sqrt(w) * pressure_mmhg ** band["k"]
    strong = None
    if band["C"] is not None:
        strong = band["C"] + band["D"] * math.log10(w) + band["K"] * math.log10(pressure_mmhg)
    if weak is None:
        total = strong
    elif strong is None or weak < band["transition"]:
        total = weak
    else:
        total = strong
    return float(min(max(total, 0.0), width))


def absorptance(spectrum: Spectrum, bands: dict, w: float, scale: dict | None = None) -> float:
    """Fraction of the star's TOTAL flux absorbed, summed over the bands.

    Howard et al. Eq. 11: the band-average fractional absorption is the total
    band absorption divided by the band width, and the fraction of the whole
    incident flux absorbed is that weighted by the flux fraction in the band.
    """
    total = 0.0
    for name, band in bands.items():
        width = band["hi"] - band["lo"]
        mean_absorptance = band_absorption(band, w) / width
        if scale is not None and name in scale:
            mean_absorptance *= scale[name]
        total += spectrum.fraction_in_band(band["lo"], band["hi"]) * mean_absorptance
    return total


def co2_column_atmos_cm(vmr: float, surface_pressure_pa: float, gravity: float) -> float:
    """The CO2 column as Howard measures absorber amount: cm of pure gas at STP.

    An atmos-cm is the depth the gas alone would occupy at 273 K and 1 atm, so
    it is the column MASS divided by the STP density of CO2. The mass column of
    a trace species is not its partial pressure over g -- that is the mass a
    PURE atmosphere of that surface pressure would have -- but

        m = vmr * (M_CO2 / M_air) * p_surface / g

    because hydrostatic balance converts total pressure, not partial pressure,
    into mass. `radmod.f90:2395` carries exactly this factor as `zpv2pm` and
    `lwr` applies it; dropping it understates the column by 1.52, which survives
    a star-over-Sun ratio unchanged and destroys every absolute the ratio is
    quoted beside. At Earth's gravity and 450 ppmv the answer is 360 atmos-cm.
    """
    mass_column = vmr * (ZMMCO2 / ZMMAIR) * surface_pressure_pa / gravity
    return mass_column / ZRCO2 * 100.0


def co2_volume_mixing_ratio(config: dict) -> tuple[float, float]:
    """CO2 volume mixing ratio and total surface pressure, from the partials."""
    partials = {k: v for k, v in config["atmosphere"].items()
                if k.startswith("p") and k.endswith("_bar")}
    total_bar = sum(float(v) for v in partials.values())
    return float(partials["pCO2_bar"]) / total_bar, total_bar * 1.0e5


def pressure_reduction(config: dict,
                       climatology: Path | None = None) -> tuple[float, str]:
    """sum(dsigma * sigma) on the model's own grid, the amount's pressure scaling.

    `swr` does not evaluate the absorptance at the true column. It reduces each
    layer's amount to an equivalent standard-pressure amount by multiplying by
    sigma * ps / p0, and then evaluates a fit stated at standard pressure --
    `radmod.f90:1997` for water vapour and `:2532` for `lwr`'s CO2. For a
    well-mixed gas that reduction is sum(dsigma * sigma) over the column, close
    to a half, and the CO2 term has to use it or it would be the one absorber in
    the scheme evaluated on a different kind of amount from the others.

    THE SIGMA GRID IS A GRID PROPERTY and not a climate state, so any run's
    climatology answers it identically: `climatology` lets a caller that has one
    but no baseline pass it, which is what `corrk_cross_check.py` does to derive
    the amounts it quotes rather than carrying them as literals.
    """
    named = climatology or config.get("baseline_climatology")
    if named:
        path = Path(named)
        if not path.is_absolute():
            path = Path(CONFIG).resolve().parents[1] / named
        if path.is_file():
            import netCDF4

            with netCDF4.Dataset(path) as data:
                sigma = np.asarray(data.variables["lev"][:])
                sigma_half = np.asarray(data.variables["levp"][:])
            dsigma = np.diff(sigma_half)
            if len(dsigma) == len(sigma):
                return float((dsigma * sigma).sum()), f"the model's own sigma grid, from {Path(named).name}"
    return 0.5, "the well-mixed limit, because no baseline climatology is named"


def co2_absorptance(
    spectrum: Spectrum,
    column_atmos_cm: float,
    water_cm: float,
    h2o: dict,
    h2o_scale: dict,
    keep_27um: bool = CO2_KEEP_27UM,
) -> tuple[float, dict]:
    """Fraction of the star's TOTAL flux absorbed by CO2, net of the H2O overlap.

    Howard's bands for CO2, weighted by the incident flux fraction the same way
    the water vapour reconstruction is, and then multiplied by the fraction of
    the interval water vapour has left. Without that factor the two gases would
    both claim the same photons, because Lacis and Hansen Eq. 21 already carries
    everything water vapour absorbs across the whole near infrared, the CO2 band
    intervals included.
    """
    rows: dict[str, dict] = {}
    total = 0.0
    for name, band in CO2_BANDS.items():
        if name == "2.7" and not keep_27um:
            continue
        width = band["hi"] - band["lo"]
        mean_absorptance = band_absorption(band, column_atmos_cm) / width
        overlap = 0.0
        for h_name, h_band in h2o.items():
            lo, hi = max(band["lo"], h_band["lo"]), min(band["hi"], h_band["hi"])
            if hi > lo:
                overlap += (
                    (hi - lo) / width
                    * band_absorption(h_band, water_cm) / (h_band["hi"] - h_band["lo"])
                    * h2o_scale.get(h_name, 1.0)
                )
        clear = max(0.0, 1.0 - overlap)
        fraction = spectrum.fraction_in_band(band["lo"], band["hi"])
        contribution = fraction * mean_absorptance * clear
        total += contribution
        rows[name] = {
            "mean_absorptance": mean_absorptance,
            "h2o_transmission_in_band": clear,
            "flux_fraction": fraction,
            "contribution": contribution,
        }
    return total, rows


def co2_closed_form(u, a1: float, b1: float, a2: float, b2: float):
    """The form radmod.f90 codes: two logarithms, positive and monotone in u.

    Lacis and Hansen's own Eq. 21 form was tried first and fits worse over the
    range that matters while wanting a negative coefficient in its denominator,
    which is a fit that can go singular on a column this scheme has no reason to
    forbid. Two logarithms are the shape Howard's own strong-band fit has, they
    stay positive and monotone for every u, and they cost two LOGs per layer.
    """
    return a1 * np.log1p(b1 * np.asarray(u, dtype=float)) + a2 * np.log1p(b2 * np.asarray(u, dtype=float))


def fit_co2_closed_form(sun: Spectrum, water_cm: float, h2o: dict, h2o_scale: dict) -> dict:
    """Fit the closed form to the solar-weighted CO2 absorptance.

    Fitted on RELATIVE residuals, because the quantity spans a decade over the
    range and an absolute fit would spend its accuracy where the absorption is
    largest and the flux is not.
    """
    from scipy.optimize import curve_fit

    lo, hi = CO2_FIT_RANGE
    u = np.logspace(math.log10(lo), math.log10(hi), 160)
    y = np.array([co2_absorptance(sun, float(ui), water_cm, h2o, h2o_scale)[0] for ui in u])
    popt, _ = curve_fit(
        co2_closed_form, u, y, p0=[1.0e-3, 1.0, 1.0e-3, 1.0e-2], sigma=y,
        bounds=([0.0, 0.0, 0.0, 0.0], [np.inf] * 4), maxfev=400000,
    )
    coefficients = [float(f"{v:.5g}") for v in popt]      # as the patch codes them
    residual = co2_closed_form(u, *coefficients) / y - 1.0
    return {
        "form": "A(u) = a1 ln(1 + b1 u) + a2 ln(1 + b2 u), u in atmos-cm",
        "a1": coefficients[0], "b1": coefficients[1],
        "a2": coefficients[2], "b2": coefficients[3],
        "range_atmos_cm": list(CO2_FIT_RANGE),
        "max_relative_error": float(np.abs(residual).max()),
        "rms_relative_error": float(np.sqrt((residual ** 2).mean())),
        "water_path_cm": water_cm,
    }


RADMOD = MODEL_SRC / "plasim" / "src" / "radmod.f90"


def radmod_co2_fit(path: Path = RADMOD) -> dict:
    """The four coefficients `swr` ACTUALLY runs, read out of the model source.

    They are not the fit this file makes. PHYS-10 replaced the Howard-band fit
    with one to HITRAN2020 through the Generic PCM correlated-k tables, by
    `corrk_cross_check.py --fit`, and this file cannot reach that data: its own
    absorptances come from Howard's band set and always will. So the model's
    coefficients are READ rather than restated, and a divergence between the
    source and this artifact becomes impossible instead of undetectable.

    Parsed rather than hardcoded for the same reason: a later refit that edits
    radmod.f90 and forgets the artifact cannot leave a superseded fit standing
    here, and a rename or a deletion in the source raises rather than passing a
    stale number through.
    """
    text = path.read_text()
    out = {}
    for key, name in (("a1", "zca1"), ("b1", "zcb1"), ("a2", "zca2"), ("b2", "zcb2")):
        found = re.findall(rf"^\s*parameter\(\s*{name}\s*=\s*([0-9.eEdD+-]+)\s*\)",
                           text, flags=re.MULTILINE)
        if len(found) != 1:
            raise ValueError(
                f"{path}: expected exactly one parameter({name}=...), found {len(found)}"
            )
        out[key] = float(found[0].replace("D", "E").replace("d", "e"))
    out["form"] = "A(u) = a1 ln(1 + b1 u) + a2 ln(1 + b2 u), u in atmos-cm"
    out["source"] = (
        "read from vendor/exoplasim/exoplasim/plasim/src/radmod.f90; fitted to "
        "HITRAN2020 through the Generic PCM correlated-k tables by "
        "exoplasim/scripts/corrk_cross_check.py --fit (PHYS-10)"
    )
    return out


def h2o_bands(include_blue: bool) -> tuple[dict, dict]:
    bands = dict(H2O_BANDS)
    scale: dict[str, float] = {}
    if include_blue:
        bands.update(H2O_BLUE_BANDS)
        scale.update(WEAK_BLUE_SCALE)
    return bands, scale


# radmod.f90:1579. The magnification factor the scheme applies to the water path
# below cloud and on the reflected beam. The absorptance is evaluated at the
# magnified path, so the weight is quoted there too.
WATER_MAGNIFICATION = 1.66


def column_water_cm(config: dict) -> tuple[float, str]:
    """This world's own effective water path, from the baseline climatology.

    A correction's size depends on the state it acts on, so the weight is quoted
    where this planet actually sits rather than where the number was first
    measured. The path is built the way `radmod.f90:1802` builds it -- the
    column integral of specific humidity with a linear pressure scaling and a
    sqrt(273/T) temperature scaling -- and then magnified, because that is the
    argument the absorptance is evaluated at. Falls back to Earth's global mean
    only if no baseline is named, and says which it used.
    """
    named = config.get("baseline_climatology")
    if named:
        path = Path(CONFIG).resolve().parents[1] / named
        if path.is_file():
            import netCDF4

            with netCDF4.Dataset(path) as data:
                hus = np.asarray(data.variables["hus"][:])
                air_t = np.asarray(data.variables["ta"][:])
                surface_p = np.asarray(data.variables["ps"][:]) * 100.0
                sigma = np.asarray(data.variables["lev"][:])
                sigma_half = np.asarray(data.variables["levp"][:])
                lat = np.asarray(data.variables["lat"][:])
            gravity = float(config["planet"]["gravity_m_s2"])
            dsigma = np.diff(sigma_half)
            column = np.zeros_like(surface_p)
            for k in range(len(sigma)):
                column += (
                    0.1 * dsigma[k] * hus[:, k] * surface_p / gravity
                    * np.sqrt(273.0 / air_t[:, k])
                    * sigma[k] * surface_p / 1.0e5
                )
            weights = np.cos(np.deg2rad(lat))[None, :, None] * np.ones_like(column)
            mean_cm = float((column * weights).sum() / weights.sum())
            return (
                mean_cm * WATER_MAGNIFICATION,
                f"radmod's own effective path, magnified by {WATER_MAGNIFICATION}, from {named}",
            )
    return 2.5, "Earth's global mean, because no baseline climatology is named"


def predict(config: dict, weight: float, absorber: str = "h2o", co2_fit: dict | None = None) -> dict:
    """What a weight does to the baseline, priced against the baseline itself.

    A correction's size depends on how much of the surface it acts on, so this
    reimplements `radmod.f90`'s own effective absorber path against the baseline
    climatology and asks how much shortwave that gas is absorbing there, rather
    than scaling a number measured somewhere else.

    `absorber` selects which term is being priced. For `h2o` the scheme already
    carries the absorptance at weight 1.0 and the change is from 1.0 to `weight`.
    For `co2` there is no term at all upstream, so the change is from ZERO to the
    whole of it, which is why the two cannot share a default.

    Everything here is a PREDICTION about a run that has not happened, which is
    the point: it is falsifiable and the run falsifies it.
    """
    import netCDF4

    # Through the one resolver, which raises when no baseline is named. This
    # prediction needs a real climatology; there is no Earth fallback for it as
    # there is for the two absorptance paths above.
    with netCDF4.Dataset(climatology_path(root=Path(CONFIG).resolve().parents[1])) as data:
        get = lambda name: np.asarray(data.variables[name][:])
        lat, sigma, sigma_half = get("lat"), get("lev"), get("levp")
        hus, air_t = get("hus"), get("ta")
        surface_p = get("ps") * 100.0
        rst, rss, rsut, ssru = get("rst"), get("rss"), get("rsut"), get("ssru")
        precip, ts = get("pr"), get("ts")

    gravity = float(config["planet"]["gravity_m_s2"])
    dsigma = np.diff(sigma_half)
    area = np.cos(np.deg2rad(lat))[None, :, None] * np.ones_like(rst)
    mean = lambda x: float((x * area).sum() / area.sum())

    column = np.zeros_like(surface_p)
    for k in range(len(sigma)):
        column += (
            0.1 * dsigma[k] * hus[:, k] * surface_p / gravity
            * np.sqrt(273.0 / air_t[:, k]) * sigma[k] * surface_p / 1.0e5
        )
    path = WATER_MAGNIFICATION * column

    incident = rst - rsut          # rsut is signed upward-negative on this stream
    upward = np.abs(ssru)
    planetary_albedo = 1.0 - mean(rst) / mean(incident)
    surface_albedo = mean(upward) / mean(rss + upward)

    # Ozone, so the decomposition of the atmosphere's shortwave absorption closes
    # and the water vapour share is a measurement rather than a fraction assumed.
    # radmod.f90:44-49's synthetic Earth column, scaled by o3scale, magnified by
    # zmbar, through the already-re-weighted Lacis and Hansen Eqs. 8 and 9.
    model = config.get("model", {})
    o3_scale = float(model.get("ozone_scale", 1.0))
    uvw = float(model.get("ozone_uv_weight", 1.0))
    visw = float(model.get("ozone_visible_weight", 1.0))
    phase = np.linspace(0.0, 1.0, rst.shape[0], endpoint=False)
    sine = np.sin(np.deg2rad(lat))[None, :, None]
    o3 = o3_scale * (
        0.25 + 0.11 * np.abs(sine)
        + 0.08 * sine * np.cos(2 * np.pi * (phase - 0.25))[:, None, None]
    )
    x = 1.9 * o3
    a_o3 = (
        visw * 0.02118 * x / (1.0 + 0.042 * x + 0.000323 * x**2)
        + uvw * 1.082 * x / ((1.0 + 138.6 * x) ** 0.805)
        + uvw * 0.0658 * x / (1.0 + (103.6 * x) ** 3)
    )

    if absorber == "h2o":
        gas_path, gas_absorptance, base_weight = path, lacis_hansen_h2o, 1.0
    elif absorber == "co2":
        if co2_fit is None:
            raise ValueError("pricing the CO2 term needs its closed-form fit")
        vmr, _ = co2_volume_mixing_ratio(config)
        reduction, _ = pressure_reduction(config)
        # Well mixed, so the column follows surface pressure and gravity alone
        # and every cell's path differs only through ps.
        gas_path = WATER_MAGNIFICATION * reduction * co2_column_atmos_cm(
            vmr, surface_p, gravity)
        gas_absorptance = lambda u: co2_closed_form(
            u, co2_fit["a1"], co2_fit["b1"], co2_fit["a2"], co2_fit["b2"])
        base_weight = 0.0
    else:
        raise ValueError(f"unknown absorber {absorber!r}")

    a_old, a_new = base_weight * gas_absorptance(gas_path), weight * gas_absorptance(gas_path)
    a2_old = base_weight * gas_absorptance(2 * gas_path)
    a2_new = weight * gas_absorptance(2 * gas_path)
    down = mean(incident * (a_new - a_old))
    reflected_leg = mean(
        upward
        * (
            np.clip((a2_new - a_new) / np.clip(1 - a_new, 1e-6, None), 0, 1)
            - np.clip((a2_old - a_old) / np.clip(1 - a_old, 1e-6, None), 0, 1)
        )
    )
    d_atmosphere = down + reflected_leg
    d_surface = -down * (1.0 - surface_albedo)

    # The baseline decomposition is a property of the run, not of the term being
    # priced, so it is always water vapour's own absorptance at weight 1.0.
    w_old, w2_old = lacis_hansen_h2o(path), lacis_hansen_h2o(2 * path)
    water_vapour_absorption = mean(incident * w_old) + mean(
        upward * np.clip((w2_old - w_old) / np.clip(1 - w_old, 1e-6, None), 0, 1)
    )
    ozone_absorption = mean(incident * a_o3)

    # How much of what the extra absorption intercepts would otherwise have been
    # reflected back to space? If all of it sat above every reflector the answer
    # is the planetary albedo; if all of it sat below the cloud, the surface
    # albedo. Water vapour is bottom-heavy and cloud is not, so the truth is
    # inside, and the bracket is carried rather than collapsed.
    toa = {
        "below_all_cloud": down * surface_albedo + reflected_leg,
        "above_all_cloud": down * planetary_albedo + reflected_leg,
    }
    toa["central"] = 0.65 * toa["below_all_cloud"] + 0.35 * toa["above_all_cloud"]

    # The flux-to-kelvin conversion is lib/sensitivity.py and nowhere else.
    # This block used to hardcode a 196/209 segment pair measured on the
    # superseded precarve-zoned-g1281 terrain -- a fourth simultaneous
    # sensitivity, which is the thing that module exists to prevent. The
    # stellar sweep's slope is still NOT usable here: it crosses the ice
    # transition (TASKS BUDG-4, and the module docstring).
    slope_lo, slope_hi = sensitivity.SLOPE_SPREAD_K_PER_FLUX_RATIO
    k_per_w = {
        "slope_spread_low": sensitivity.kelvin_per_w_m2(
            planetary_albedo, slope=slope_lo),
        "slope_central": sensitivity.kelvin_per_w_m2(planetary_albedo),
        "slope_spread_high": sensitivity.kelvin_per_w_m2(
            planetary_albedo, slope=slope_hi),
    }
    warming = {
        name: {k: v * s for k, s in k_per_w.items()} for name, v in toa.items()
    }

    latent = 2.5008e6
    seconds_per_year = 86400.0 * 365.25
    precip_mm_yr = mean(precip) * 1000.0 * seconds_per_year
    d_precip_full = -d_atmosphere / latent * seconds_per_year

    central_warming = warming["central"]["slope_central"]
    return {
        "absorber": absorber,
        "weight": weight,
        "effective_path_priced_at": float(mean(gas_path)),
        "baseline": {
            "incident_toa_shortwave": mean(incident),
            "toa_net_shortwave_rst": mean(rst),
            "surface_net_shortwave_rss": mean(rss),
            "atmospheric_shortwave_absorption": mean(rst - rss),
            "planetary_albedo": planetary_albedo,
            "surface_albedo": surface_albedo,
            "effective_water_path_cm": mean(path),
            "water_vapour_shortwave_absorption": water_vapour_absorption,
            "ozone_shortwave_absorption": ozone_absorption,
            "cloud_and_scattering_residual": (
                mean(rst - rss) - water_vapour_absorption - ozone_absorption
            ),
            "precipitation_mm_yr": precip_mm_yr,
            "mean_surface_temperature_k": mean(ts),
        },
        "predicted_change": {
            "atmospheric_shortwave_absorption": d_atmosphere,
            "of_which_downward_beam": down,
            "of_which_reflected_leg": reflected_leg,
            "surface_net_shortwave": d_surface,
            "toa_net_shortwave": toa,
            "kelvin_per_w_m2": k_per_w,
            "mean_surface_temperature_k": warming,
            "precipitation_mm_yr_full_compensation": d_precip_full,
            # INVARIANT to the slope, and written this way because it is the
            # quantity rather than the arithmetic that is meant. `central_warming`
            # is the same forcing times `slope / absorbed`, so the slope cancels
            # and this reduces to the TOA change over absorbed flux per unit flux
            # ratio. It is the one number in this report that does NOT move when
            # `lib/sensitivity.py` is re-measured; every kelvin above does.
            "flux_ratio_to_restore_the_design_mean": -central_warming / sensitivity.SLOPE_K_PER_FLUX_RATIO,
        },
    }


def weight_curve(star: Spectrum, sun: Spectrum, bands: dict, scale: dict, amounts) -> list[dict]:
    rows = []
    for w in amounts:
        a_sun = absorptance(sun, bands, w, scale)
        a_star = absorptance(star, bands, w, scale)
        rows.append(
            {
                "w": w,
                "absorptance_solar": a_sun,
                "absorptance_star": a_star,
                "weight": a_star / a_sun,
                "lacis_hansen_eq21": float(lacis_hansen_h2o(w)),
                "fowle_eq22": float(fowle_h2o(w)),
                "korb_eq23": float(korb_h2o(w)),
            }
        )
    return rows


def blue_bands_against_correlated_k(refresh: bool = False) -> None:
    """Re-measure the dry-end attribution and the two blue scale factors.

    WRITES NOTHING. It needs the LMD Generic PCM bundle that
    `corrk_cross_check.py` reads, which lives outside this repository, so it is
    a mode rather than part of the report: everything the report carries from
    this measurement is a declared constant above, exactly as
    `CORRK_RATIO_TO_EQ21` is, and regenerating the report needs no bundle.

    The import is deferred because `corrk_cross_check` imports THIS module at
    module level, and a top-level import here would close the cycle.
    """
    import corrk_cross_check as ck  # noqa: PLC0415  deferred: see the docstring

    config = yaml.safe_load(CONFIG.read_text())
    teff = float(config["star"]["effective_temperature_k"])
    sun = Spectrum(*blend(SOLAR_TEFF, 5700, 5800, refresh)[:2], "G2V reference", {})
    star = Spectrum(*blend(teff, 4900, 5000, refresh)[:2],
                    config["star"]["spectral_type"], {})
    table = ck.CorrK(ck.TABLE_376)
    f_sun = ck.flux_fractions(table, sun)
    f_star = ck.flux_fractions(table, star)

    def corrk_absorptance(w: float, t_k: float = ck.DEFAULT_TEMPERATURE,
                          q: float = ck.DEFAULT_WATER_VMR) -> np.ndarray:
        """H2O-only absorptance per correlated-k band, the note's own isolation."""
        u = ck.air_column_for_water(w, q)
        wet = table.transmission(ck.P_STANDARD_MBAR, t_k, q, u)
        dry = table.transmission(ck.P_STANDARD_MBAR, t_k, ck.DRY, u)
        return 1.0 - wet / dry

    def interval_mean(per_band: np.ndarray, lo: float, hi: float) -> float:
        """Flux-weighted correlated-k band mean over one of Howard's intervals."""
        num = den = 0.0
        for (edge_lo, edge_hi), value in zip(table.edges, per_band):
            left, right = max(lo, edge_lo), min(hi, edge_hi)
            if right <= left:
                continue
            f = sun.fraction_in_band(left, right)
            num += f * value
            den += f
        return num / den

    bands, scale = h2o_bands(include_blue=True)
    bands_no_blue, scale_no_blue = h2o_bands(include_blue=False)
    amounts = (0.01, 0.03, 0.1, 0.3, 1.0, CORRK_PATH_CM, 5.0, 10.0)

    print("THE DRY-END ATTRIBUTION. f_blue is the share of the gap between this "
          "reconstruction\nand correlated-k that the 0.72 and 0.81 um bands carry.\n")
    print(f"{'w, cm':>8} {'R_full':>8} {'R_noblue':>9} {'R_ck':>8} {'f_blue':>8} {'noblue-ck':>10}")
    for w in amounts:
        eq21 = float(lacis_hansen_h2o(w))
        full = absorptance(sun, bands, w, scale) / eq21
        no_blue = absorptance(sun, bands_no_blue, w, scale_no_blue) / eq21
        corrk = ck.broadband(f_sun, corrk_absorptance(w)) / eq21
        # The share of a gap is only a share where there IS a gap. Inside the
        # tolerance the two determinations agree, and dividing by that residual
        # returns whatever the last digit of each happens to be.
        gap = full / corrk - 1.0
        share = (f"{(full - no_blue) / (full - corrk):8.3f}"
                 if gap > CORRK_AGREEMENT else f"{'--':>8}")
        print(f"{w:8.4f} {full:8.4f} {no_blue:9.4f} {corrk:8.4f} {share} "
              f"{no_blue / corrk - 1:9.2%}")

    print("\nPER BAND: the correlated-k mean over each Howard interval, divided by "
          "what this\nfile puts there. The dry-end error is monotone in how weak "
          "the band is.\n")
    all_bands = dict(H2O_BANDS)
    all_bands.update(H2O_BLUE_BANDS)
    print(f"{'w, cm':>8} " + " ".join(f"{name:>7}" for name in all_bands))
    for w in amounts:
        per_band = corrk_absorptance(w)
        cells = []
        for name, band in all_bands.items():
            shape = band_absorption(band, w) / (band["hi"] - band["lo"])
            cells.append(interval_mean(per_band, band["lo"], band["hi"]) / shape)
        print(f"{w:8.4f} " + " ".join(f"{c:7.4f}" for c in cells))
    print("\nThe last two columns ARE the measured WEAK_BLUE_SCALE. Declared: "
          + ", ".join(f"{k} {v:g}" for k, v in WEAK_BLUE_SCALE.items()))

    # WHY NO CONSTANT IS RIGHT: the power of w the correlated-k band mean over
    # each blue interval follows, against the 1/2 Howard's weak fit imposes. A
    # square-root law is the strong-line regime and these two bands are not in
    # it, so the scale factor has to absorb the difference between the laws and
    # cannot be one number.
    ends = (amounts[0], amounts[-1])
    span = math.log10(ends[1] / ends[0])
    per_band = [corrk_absorptance(w) for w in ends]
    for name, band in H2O_BLUE_BANDS.items():
        lo, hi = (interval_mean(p, band["lo"], band["hi"]) for p in per_band)
        print(f"  {name} um: the correlated-k band mean goes as w^"
              f"{math.log10(hi / lo) / span:.2f} over {ends[0]:g} to {ends[1]:g} cm, "
              f"against the w^0.5 Howard's weak fit imposes")

    print("\nWHAT SUBSTITUTING THE MEASURED FACTORS WOULD COST, at the operating path.")
    eq21 = float(lacis_hansen_h2o(CORRK_PATH_CM))
    arms = (("declared", bands, scale), ("measured", bands, WEAK_BLUE_MEASURED),
            ("dropped", bands_no_blue, scale_no_blue))
    for label, arm_bands, arm_scale in arms:
        a_sun = absorptance(sun, arm_bands, CORRK_PATH_CM, arm_scale)
        a_star = absorptance(star, arm_bands, CORRK_PATH_CM, arm_scale)
        print(f"  blue pair {label:>8}: ratio_to_eq21 {a_sun / eq21:.4f} "
              f"({a_sun / eq21 / CORRK_RATIO_TO_EQ21 - 1:+.2%} from correlated-k)  "
              f"h2osww {a_star / a_sun:.4f}")
    per_band = corrk_absorptance(CORRK_PATH_CM)
    a_sun = ck.broadband(f_sun, per_band)
    a_star = ck.broadband(f_star, per_band)
    print(f"  correlated-k itself: ratio_to_eq21 {a_sun / eq21:.4f}"
          f"{'':>28}h2osww {a_star / a_sun:.4f}")
    print("\nexoplasim/notes/corrk-cross-check.md carries the argument and the "
          "criterion,\nwhich was fixed before any of this was computed.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="refetch the grid points")
    parser.add_argument("--verify", action="store_true", help="print the checks in full")
    parser.add_argument(
        "--water",
        type=float,
        default=None,
        help="precipitable water in cm to quote the weight at; default is the config value",
    )
    parser.add_argument(
        "--level",
        action="store_true",
        help="derive h2o_sw_level into its own artifact and stop; needs no "
             "spectrum, no climatology and no network",
    )
    parser.add_argument(
        "--blue",
        action="store_true",
        help="re-measure the dry-end attribution and the two blue scale factors "
             "against the correlated-k tables and stop; writes nothing, and "
             "needs the bundle corrk_cross_check.py reads",
    )
    args = parser.parse_args()

    if args.blue:
        blue_bands_against_correlated_k(args.refresh)
        return

    if args.level:
        report = h2o_sw_level_report()
        ANALYSIS.mkdir(parents=True, exist_ok=True)
        out = ANALYSIS / "h2o_sw_level.json"
        out.write_text(json.dumps(report, indent=2) + "\n")
        lo, hi = report["bracket"]
        print(f"h2o_sw_level {report['value']:.6f}  bracket {lo:.6f} to {hi:.6f}")
        print(f"  = {CORRK_RATIO_TO_EQ21:.6f} x (1 + {H2O_CONTINUUM_FRACTION:g}), "
              f"the bracket over the continuum fraction's "
              f"{H2O_CONTINUUM_FRACTION_BRACKET[0]:g} to "
              f"{H2O_CONTINUUM_FRACTION_BRACKET[1]:g}")
        print(f"  config/planet.yaml carries model.h2o_sw_level to three "
              f"decimals: {report['value_rounded_for_config']:g}")
        print(f"wrote {out}")
        return

    config = yaml.safe_load(CONFIG.read_text())
    teff = float(config["star"]["effective_temperature_k"])

    star_wave, star_flux, star_meta = blend(teff, 4900, 5000, args.refresh)
    sun_wave, sun_flux, sun_meta = blend(SOLAR_TEFF, 5700, 5800, args.refresh)
    star = Spectrum(star_wave, star_flux, config["star"]["spectral_type"], star_meta)
    sun = Spectrum(sun_wave, sun_flux, "G2V reference", sun_meta)
    blackbody = Spectrum(
        star_wave,
        planck(star_wave, teff),
        f"{teff:.0f} K blackbody",
        {"note": "what solarini actually builds while NSTARFILE = 0"},
    )

    # --- checks that can fail -------------------------------------------------
    solar_above_09 = sun.fraction_between_um(0.9, 1.0e6)
    solar_below_075 = sun.fraction_between_um(0.0, 0.75)
    star_below_075 = star.fraction_between_um(0.0, 0.75)
    star_above_075 = star.fraction_between_um(0.75, 1.0e6)

    checks = {
        "solar_fraction_above_0.9um": {
            "computed": solar_above_09,
            "expected": LH74_SOLAR_ACTIVE_FRACTION,
            "source": "Lacis and Hansen (1974) Section 5a and Table 1",
            "relative_error": solar_above_09 / LH74_SOLAR_ACTIVE_FRACTION - 1.0,
        },
        "solar_fraction_below_0.75um": {
            "computed": solar_below_075,
            "expected": RADMOD_ZSOLAR1,
            "source": "radmod.f90:125 zsolar1, stated to hold at 5772 K",
            "relative_error": solar_below_075 / RADMOD_ZSOLAR1 - 1.0,
        },
    }

    bands, scale = h2o_bands(include_blue=True)
    amounts = [0.01, 0.03, 0.1, 0.3, 1.0, 2.0, 3.0, 5.0, 10.0]
    curve = weight_curve(star, sun, bands, scale, amounts)

    # THE CHECK THAT CAN FAIL. The reconstruction against the envelope of the
    # three published absorptivity determinations Lacis and Hansen plot together
    # in their Fig. 11, over the only interval they state for Eq. 21, widened by
    # the accuracy Howard state for the band absorptions the reconstruction is
    # built from. The constants and the argument are beside `korb_h2o` above.
    lo_cm, hi_cm = LH74_FIT_RANGE_CM
    inside = [row for row in curve if lo_cm <= row["w"] <= hi_cm]
    ratios = [row["absorptance_solar"] / row["lacis_hansen_eq21"] for row in inside]
    margin = HOWARD_BAND_ABSORPTION_ACCURACY
    envelope = []
    for row in inside:
        published = [row["lacis_hansen_eq21"], row["fowle_eq22"], row["korb_eq23"]]
        band_lo = min(published) * (1.0 - margin)
        band_hi = max(published) * (1.0 + margin)
        a = row["absorptance_solar"]
        envelope.append({
            "w": row["w"],
            "reconstruction": a,
            "yamamoto_eq21": row["lacis_hansen_eq21"],
            "fowle_eq22": row["fowle_eq22"],
            "korb_eq23": row["korb_eq23"],
            "band": [band_lo, band_hi],
            # How far outside the band it sits, as a factor. 1.0 is on the edge
            # and below 1.0 is inside, so the worst row is the maximum.
            "excess": max(a / band_hi, band_lo / a),
            "inside": bool(band_lo <= a <= band_hi),
        })
    worst = max(envelope, key=lambda r: r["excess"])
    checks["reconstruction_vs_published_absorptivity_envelope"] = {
        "ratio_to_eq21_min": min(ratios),
        "ratio_to_eq21_max": max(ratios),
        "ratio_to_eq21_median": float(np.median(ratios)),
        "range_cm": list(LH74_FIT_RANGE_CM),
        "howard_margin": margin,
        "per_amount": envelope,
        "worst_w": worst["w"],
        "worst_excess": worst["excess"],
        "inside": all(r["inside"] for r in envelope),
        "gates": False,
        "why_not_a_gate": (
            "HITRAN2020 correlated-k falls outside this same widened envelope "
            "at every amount from 0.1 to 10 cm, by 6 to 10 per cent, and that "
            "is a floor because it carries no water vapour continuum either. A "
            "bar a modern line list misses is not a bar on this "
            "reconstruction: exoplasim/notes/corrk-cross-check.md. Reported, "
            "not widened, and not raised on"),
        "source": (
            "Lacis and Hansen (1974) Eqs. 21, 22 and 23, the three curves of "
            "their Fig. 11, widened by Howard, Burch and Williams (1956) "
            "+/-3% on the band absorptions; declared as an inter-formula "
            "spread and not as an error bar, per LH74 p. 127"),
    }

    # THE GATE. The same quantity from a different absorption dataset, at the
    # path corrk-cross-check.md measured on. Climatology-free by construction:
    # the path is that note's recorded number, not this run's column.
    a_corrk_path = absorptance(sun, bands, CORRK_PATH_CM, scale)
    ratio_corrk_path = a_corrk_path / float(lacis_hansen_h2o(CORRK_PATH_CM))
    apart = ratio_corrk_path / CORRK_RATIO_TO_EQ21 - 1.0
    checks["reconstruction_vs_correlated_k"] = {
        "w": CORRK_PATH_CM,
        "reconstruction_over_eq21": ratio_corrk_path,
        "correlated_k_over_eq21": CORRK_RATIO_TO_EQ21,
        "apart": apart,
        "tolerance": CORRK_AGREEMENT,
        "inside": bool(abs(apart) <= CORRK_AGREEMENT),
        "gates": True,
        "source": (
            "exoplasim/notes/corrk-cross-check.md, HITRAN2020 through the "
            "LMD Generic PCM correlated-k tables at the same homogeneous "
            "760 mm Hg path; tolerance is Howard, Burch and Williams (1956) "
            "+/-3%, the only stated accuracy either side carries"),
    }
    # THE DRY END, REPORTED AND NOT GATED. The two determinations part below the
    # amounts the model evaluates, and this records which half of the
    # reconstruction does it, so a change to the band set moves a number here
    # rather than nothing. `--blue` re-measures the correlated-k side; the
    # constants it is compared against are declared, so this node needs no
    # bundle. It does not raise, for the reason `world-njlb` records: a T42
    # column spans roughly 0.3 to 5 cm, and the weight is a ratio in which a
    # dry-end level error largely divides out.
    dry_end = []
    dry_bands, dry_scale = h2o_bands(include_blue=False)
    for w, corrk in sorted(CORRK_RATIO_TO_EQ21_DRY.items()):
        eq21 = float(lacis_hansen_h2o(w))
        full = absorptance(sun, bands, w, scale) / eq21
        no_blue = absorptance(sun, dry_bands, w, dry_scale) / eq21
        dry_end.append({
            "w": w,
            "reconstruction_over_eq21": full,
            "reconstruction_without_blue_over_eq21": no_blue,
            "correlated_k_over_eq21": corrk,
            "apart": full / corrk - 1.0,
            "apart_without_blue": no_blue / corrk - 1.0,
            "blue_share_of_gap": (full - no_blue) / (full - corrk),
            "blue_share_of_gap_measured": BLUE_SHARE_OF_DRY_GAP[w],
        })
    checks["reconstruction_vs_correlated_k_dry_end"] = {
        "per_amount": dry_end,
        "attributed_to": "howard weak-band fit extrapolated below the measured amounts",
        "gates": False,
        "why_not_a_gate": (
            "the model never evaluates there: a T42 column spans roughly 0.3 "
            "to 5 cm, and h2o_sw_weight is a ratio in which a dry-end level "
            "error largely divides out. The gate is at CORRK_PATH_CM"),
        "source": (
            "exoplasim/notes/corrk-cross-check.md, the dry-end section; "
            "`--blue` re-measures the correlated-k side and the per-band "
            "attribution behind the verdict"),
    }
    # Reported here as well because it belongs beside the ratio it is built on,
    # but `exoplasim/analysis/h2o_sw_level.json` is the artifact config is
    # checked against: this report cannot be written without a climatology and
    # that one can, so the tie to `config/planet.yaml` hangs off the half that
    # is always regenerable.
    checks["h2o_sw_level"] = h2o_sw_level_report()
    checks["h2o_sw_level"]["artifact"] = "exoplasim/analysis/h2o_sw_level.json"

    water = args.water
    water_source = "given on the command line"
    if water is None:
        water, water_source = column_water_cm(config)
    a_sun = absorptance(sun, bands, water, scale)
    a_star = absorptance(star, bands, water, scale)
    h2o_weight = a_star / a_sun

    bands_no_blue, scale_no_blue = h2o_bands(include_blue=False)
    h2o_weight_no_blue = absorptance(star, bands_no_blue, water, scale_no_blue) / absorptance(
        sun, bands_no_blue, water, scale_no_blue
    )
    # The reconstruction sits above Lacis and Hansen Eq. 21 by a roughly flat
    # factor, so the level divides out of a ratio. What would NOT divide out is
    # the excess being concentrated in one part of the spectrum, so the extreme
    # attribution is tested: charge all of it to the 2.7, 3.2 and 6.3 um complex
    # and drop that complex entirely. That is Yamamoto's own curve 1, and it is
    # the low end of the bracket.
    bands_curve1 = {k: v for k, v in bands.items() if k not in ("6.3", "3.2", "2.7")}
    h2o_weight_curve1 = absorptance(star, bands_curve1, water, scale) / absorptance(
        sun, bands_curve1, water, scale
    )
    bracket = (
        min(h2o_weight, h2o_weight_no_blue, h2o_weight_curve1),
        max(h2o_weight, h2o_weight_no_blue, h2o_weight_curve1),
    )

    # The weight the model would need TODAY, while solarini is building a Planck
    # curve instead of reading k25v. Lower, because a blackbody has no line
    # blanketing pushing flux out of the blue and into the near infrared, so it
    # is less red than the star it stands for.
    h2o_weight_blackbody = absorptance(blackbody, bands, water, scale) / a_sun

    # CO2. `swr` has no term to re-weight, so this is what the missing absorber
    # is worth, what its closed form is, and the number the patch's namelist key
    # carries. Three quantities have to be kept apart and were not on the first
    # pass: the TRUE column, the pressure-REDUCED amount the scheme evaluates a
    # standard-pressure fit at, and the MAGNIFIED path it evaluates it on.
    gravity = float(config["planet"]["gravity_m_s2"])
    co2_vmr, surface_pressure_pa = co2_volume_mixing_ratio(config)
    co2_ppmv = co2_vmr * 1.0e6
    co2_column = co2_column_atmos_cm(co2_vmr, surface_pressure_pa, gravity)
    co2_column_earth = co2_column_atmos_cm(
        co2_vmr, EARTH_SURFACE_PRESSURE_PA, EARTH_GRAVITY)
    reduction, reduction_source = pressure_reduction(config)
    co2_path = co2_column * reduction * WATER_MAGNIFICATION
    co2_path_earth = co2_column_earth * reduction * WATER_MAGNIFICATION

    co2_sun, co2_rows = co2_absorptance(sun, co2_path, water, bands, scale)
    co2_star, _ = co2_absorptance(star, co2_path, water, bands, scale)
    co2_blackbody, _ = co2_absorptance(blackbody, co2_path, water, bands, scale)
    co2_weight = co2_star / co2_sun
    co2_weight_blackbody = co2_blackbody / co2_sun

    # Yamamoto's own choice, carried as the bracket: drop the 2.7 um CO2 band
    # outright "because of overlapping by the strong 2.7 um H2O band". It is kept
    # in the central number and charged with what water vapour leaves it instead,
    # so this says what the decision is worth rather than burying it.
    co2_sun_no27, _ = co2_absorptance(sun, co2_path, water, bands, scale, keep_27um=False)
    co2_star_no27, _ = co2_absorptance(star, co2_path, water, bands, scale, keep_27um=False)
    co2_sun_earth, _ = co2_absorptance(sun, co2_path_earth, water, bands, scale)
    co2_sun_earth_no27, _ = co2_absorptance(
        sun, co2_path_earth, water, bands, scale, keep_27um=False)

    # THE CHECK THAT CAN FAIL. Earth's column, Earth's gravity, Earth's mean
    # insolation, against the 1.5 to 2.5 W/m2 the literature measures for Earth's
    # near-infrared CO2 solar absorption. Nothing here was tuned to it, and it
    # has to pass with the 2.7 um band in and with it out, or the bracket would
    # be deciding the check.
    lo, hi = EARTH_CO2_SHORTWAVE_W_M2
    earth_w_m2 = co2_sun_earth * EARTH_MEAN_INSOLATION
    earth_w_m2_no27 = co2_sun_earth_no27 * EARTH_MEAN_INSOLATION
    checks["co2_earth_shortwave_absorption_w_m2"] = {
        "computed": earth_w_m2,
        "computed_dropping_the_2.7um_band": earth_w_m2_no27,
        "expected_range": [lo, hi],
        "inside": bool(lo <= earth_w_m2 <= hi),
        "inside_dropping_the_2.7um_band": bool(lo <= earth_w_m2_no27 <= hi),
        "source": "Earth's measured near-infrared CO2 solar absorption",
    }

    # Two fits, and they are not the same object. `co2_fit` is this file's own,
    # to Howard's bands, and it is what the patch header codes and what
    # `corrk_cross_check.py --fit` compares against. `model_fit` is what `swr`
    # runs today. Every PREDICTION about a model run is priced on the model's,
    # because a prediction about the other one is a prediction about no code.
    co2_fit = fit_co2_closed_form(sun, water, bands, scale)
    model_fit = radmod_co2_fit()

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "star": {
            "spectral_type": config["star"]["spectral_type"],
            "effective_temperature_k": teff,
            "flux_fraction_below_0.75um": star_below_075,
            "flux_fraction_above_0.75um": star_above_075,
            "provenance": star_meta,
        },
        "solar_reference": {
            "effective_temperature_k": SOLAR_TEFF,
            "flux_fraction_below_0.75um": solar_below_075,
            "flux_fraction_above_0.9um": solar_above_09,
            "provenance": sun_meta,
        },
        "checks": checks,
        "h2o": {
            "column_water_cm": water,
            "column_water_source": water_source,
            "absorptance_solar": a_sun,
            "absorptance_star": a_star,
            "weight": h2o_weight,
            "weight_rounded_for_namelist": round(h2o_weight, 3),
            "weight_against_the_blackbody_the_model_is_actually_using": h2o_weight_blackbody,
            "weight_without_0.72_0.81um_bands": h2o_weight_no_blue,
            "weight_without_2.7_3.2_6.3um_complex": h2o_weight_curve1,
            "bracket": list(bracket),
            "weight_vs_water_amount": curve,
            "band_flux_fractions": {
                name: {
                    "solar": sun.fraction_in_band(band["lo"], band["hi"]),
                    "star": star.fraction_in_band(band["lo"], band["hi"]),
                    "ratio": star.fraction_in_band(band["lo"], band["hi"])
                    / sun.fraction_in_band(band["lo"], band["hi"]),
                }
                for name, band in bands.items()
            },
        },
        "prediction": {
            name: predict(config, w)
            for name, w in (
                ("central", round(h2o_weight, 3)),
                ("bracket_low", round(bracket[0], 3)),
                ("bracket_high", round(bracket[1], 3)),
                ("blackbody_star", round(h2o_weight_blackbody, 3)),
            )
        },
        "co2": {
            "ppmv": co2_ppmv,
            "column_atmos_cm": co2_column,
            "column_atmos_cm_at_earth_gravity": co2_column_earth,
            "pressure_reduction": reduction,
            "pressure_reduction_source": reduction_source,
            "effective_path_atmos_cm": co2_path,
            "effective_path_atmos_cm_earth": co2_path_earth,
            "absorptance_solar": co2_sun,
            "absorptance_star": co2_star,
            "absorptance_solar_dropping_the_2.7um_band": co2_sun_no27,
            "absorptance_star_dropping_the_2.7um_band": co2_star_no27,
            "weight": co2_weight,
            "weight_rounded_for_namelist": round(co2_weight, 3),
            "weight_dropping_the_2.7um_band": co2_star_no27 / co2_sun_no27,
            "weight_against_the_blackbody_the_model_is_actually_using": co2_weight_blackbody,
            "keeps_the_2.7um_band": CO2_KEEP_27UM,
            "closed_form_fit": co2_fit,
            "closed_form_fit_in_radmod": model_fit,
            "per_band": co2_rows,
            "note": (
                "Upstream ExoPlaSim's shortwave has no CO2 absorptance at all: "
                "radmod.f90 carries CO2 only in lwr, from Sasamori (1968). The "
                "absorptance here is solar-weighted, and it is Howard's band set: "
                "it is what exoplasim/patches/exoplasim-3.4.2-co2-shortwave.patch "
                "codes and NOT what the fork runs. PHYS-10 refitted the closed "
                "form to HITRAN2020 through the correlated-k tables, 7.2% weaker "
                "at this planet's path, and closed_form_fit_in_radmod is that fit "
                "read out of the model source. `weight` is the namelist key "
                "co2sww that re-weights the absorptance for this star, and 0.0 "
                "leaves the term absent as upstream has it."
            ),
        },
        "prediction_co2": {
            "central": predict(config, round(co2_weight, 3), "co2", model_fit),
            "solar_weighted": predict(config, 1.0, "co2", model_fit),
            "absorptance_source": (
                "radmod.f90's own coefficients, closed_form_fit_in_radmod, not "
                "the Howard fit in closed_form_fit"
            ),
        },
    }

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    out = ANALYSIS / "shortwave_band_weights.json"
    out.write_text(json.dumps(report, indent=2) + "\n")

    print(f"star   {config['star']['spectral_type']} at {teff:.0f} K")
    print(f"  flux below 0.75 um  {star_below_075:.4f}   above  {star_above_075:.4f}")
    print(
        f"  the {teff:.0f} K blackbody solarini actually builds: "
        f"below {blackbody.fraction_between_um(0.0, 0.75):.4f}   "
        f"above {blackbody.fraction_between_um(0.75, 1.0e6):.4f}"
    )
    print(f"solar reference at {SOLAR_TEFF:.0f} K")
    print(f"  flux below 0.75 um  {solar_below_075:.4f} against radmod's {RADMOD_ZSOLAR1}")
    print(f"  flux above 0.90 um  {solar_above_09:.4f} against LH74's {LH74_SOLAR_ACTIVE_FRACTION}")
    rec = checks["reconstruction_vs_published_absorptivity_envelope"]
    print(
        f"reconstruction / LH74 Eq. 21 over {rec['range_cm'][0]}-"
        f"{rec['range_cm'][1]} cm: "
        f"{rec['ratio_to_eq21_min']:.3f} to {rec['ratio_to_eq21_max']:.3f}, "
        f"median {rec['ratio_to_eq21_median']:.3f}"
    )
    print(
        f"  RECORDED, not a gate: against the envelope of Eqs. 21, 22 and 23 "
        f"widened by Howard's {rec['howard_margin']*100:.0f}%: "
        f"{'inside at every water amount' if rec['inside'] else 'OUTSIDE'}"
        f"; worst at w = {rec['worst_w']} cm, "
        f"a factor of {rec['worst_excess']:.3f} past the band edge."
        " Correlated-k misses the same envelope, so it does not gate"
    )
    ck2 = checks["reconstruction_vs_correlated_k"]
    lev = checks["h2o_sw_level"]
    print(
        f"  CHECK, against HITRAN2020 correlated-k at w = {ck2['w']} cm: "
        f"{ck2['reconstruction_over_eq21']:.4f} here against "
        f"{ck2['correlated_k_over_eq21']:.4f} there, "
        f"{abs(ck2['apart'])*100:.2f}% apart, tolerance "
        f"{ck2['tolerance']*100:.0f}%: "
        f"{'agrees' if ck2['inside'] else 'DISAGREES'}"
    )
    print(
        f"  h2o_sw_level for config/planet.yaml: {lev['value']}, bracket "
        f"{lev['bracket'][0]} to {lev['bracket'][1]} -- that ratio with the "
        "water vapour continuum, which neither side of it carries"
    )
    print(f"water path {water:.3f} cm, from {water_source}")
    print(f"H2O shortwave weight at that path: {h2o_weight:.4f}")
    print(f"  without the 0.72/0.81 um bands: {h2o_weight_no_blue:.4f}")
    print(f"  without the 6.3, 3.2 and 2.7 um complex: {h2o_weight_curve1:.4f}")
    print(f"  bracket over which bands are counted: {bracket[0]:.3f} to {bracket[1]:.3f}")
    print(f"  h2o_sw_weight for config/planet.yaml: {round(h2o_weight, 3)}")
    print(f"  against the blackbody the model is running today: {h2o_weight_blackbody:.4f}")

    ck = checks["co2_earth_shortwave_absorption_w_m2"]
    print(f"\nCO2, the absorber the scheme does not have")
    print(f"  column {co2_column:.1f} atmos-cm against Earth's {co2_column_earth:.1f} "
          f"at the same mixing ratio, pressure-reduced by {reduction:.4f}, "
          f"path {co2_path:.1f}")
    print(f"  CHECK, Earth's column and Earth's insolation: {ck['computed']:.2f} W/m2 "
          f"against a measured {ck['expected_range'][0]}-{ck['expected_range'][1]} "
          f"-- {'inside' if ck['inside'] else 'OUTSIDE'}; "
          f"dropping the 2.7 um band {ck['computed_dropping_the_2.7um_band']:.2f} "
          f"-- {'inside' if ck['inside_dropping_the_2.7um_band'] else 'OUTSIDE'}")
    print(f"  absorptance solar {co2_sun:.6f}  this star {co2_star:.6f}")
    print(f"  co2_sw_weight for config/planet.yaml: {round(co2_weight, 3)} "
          f"(2.7 um band dropped: {co2_star_no27 / co2_sun_no27:.3f}; "
          f"blackbody: {co2_weight_blackbody:.3f})")
    print(f"  closed form fitted here, to Howard's bands: a1={co2_fit['a1']:.6g} "
          f"b1={co2_fit['b1']:.6g} a2={co2_fit['a2']:.6g} b2={co2_fit['b2']:.6g}")
    print(f"    max relative error {co2_fit['max_relative_error']:.4f}, "
          f"rms {co2_fit['rms_relative_error']:.4f}, over "
          f"{co2_fit['range_atmos_cm'][0]:.0f}-{co2_fit['range_atmos_cm'][1]:.0f} atmos-cm")
    print(f"  closed form radmod.f90 RUNS, fitted to the line list: "
          f"a1={model_fit['a1']:.6g} b1={model_fit['b1']:.6g} "
          f"a2={model_fit['a2']:.6g} b2={model_fit['b2']:.6g}")
    a_howard = float(co2_closed_form(co2_path, co2_fit["a1"], co2_fit["b1"],
                                     co2_fit["a2"], co2_fit["b2"]))
    a_model = float(co2_closed_form(co2_path, model_fit["a1"], model_fit["b1"],
                                    model_fit["a2"], model_fit["b2"]))
    print(f"    at this planet's path the two differ by "
          f"{100 * (a_howard / a_model - 1):+.1f}% -- the predictions below are "
          f"priced on the model's")
    pc = report["prediction_co2"]["central"]["predicted_change"]
    print(f"  atmospheric shortwave absorption {pc['atmospheric_shortwave_absorption']:+.2f} W/m2, "
          f"surface {pc['surface_net_shortwave']:+.2f}, "
          f"top of atmosphere {pc['toa_net_shortwave']['central']:+.2f} "
          f"({pc['toa_net_shortwave']['below_all_cloud']:+.2f} to "
          f"{pc['toa_net_shortwave']['above_all_cloud']:+.2f})")
    print(f"  mean surface temperature "
          f"{pc['mean_surface_temperature_k']['central']['slope_central']:+.2f} K, "
          f"SAME SIGN as the water vapour correction and about a sixth of its size")

    p = report["prediction"]["central"]
    base, change = p["baseline"], p["predicted_change"]
    print(f"\nprediction for a baseline re-run at h2osww = {p['weight']}")
    print(
        f"  of the {base['atmospheric_shortwave_absorption']:.1f} W/m2 the atmosphere takes: "
        f"water vapour {base['water_vapour_shortwave_absorption']:.1f}, "
        f"ozone {base['ozone_shortwave_absorption']:.1f}, "
        f"cloud and scattering {base['cloud_and_scattering_residual']:.1f}"
    )
    print(f"  atmospheric shortwave absorption {change['atmospheric_shortwave_absorption']:+.1f} W/m2")
    print(f"  surface net shortwave            {change['surface_net_shortwave']:+.1f} W/m2")
    print(
        f"  top-of-atmosphere net shortwave  {change['toa_net_shortwave']['central']:+.1f} W/m2 "
        f"(bracket {change['toa_net_shortwave']['below_all_cloud']:+.1f} to "
        f"{change['toa_net_shortwave']['above_all_cloud']:+.1f})"
    )
    print(
        f"  mean surface temperature         "
        f"{change['mean_surface_temperature_k']['central']['slope_central']:+.2f} K "
        f"from {base['mean_surface_temperature_k']:.2f} K"
    )
    print(
        f"  precipitation, full compensation "
        f"{change['precipitation_mm_yr_full_compensation']:+.0f} mm/yr on "
        f"{base['precipitation_mm_yr']:.0f}"
    )
    print(f"  flux ratio must move by          {change['flux_ratio_to_restore_the_design_mean']:+.4f}")

    if args.verify:
        print("\nweight against column water:")
        print(f"  {'w (cm)':>8} {'A solar':>10} {'A star':>10} {'weight':>8} {'LH74 Eq21':>10}")
        for row in curve:
            print(
                f"  {row['w']:8.2f} {row['absorptance_solar']:10.4f} "
                f"{row['absorptance_star']:10.4f} {row['weight']:8.4f} "
                f"{row['lacis_hansen_eq21']:10.4f}"
            )
        print("\nper band, flux fraction and the star/Sun ratio:")
        for name, row in report["h2o"]["band_flux_fractions"].items():
            print(
                f"  {name:>5} um  solar {row['solar']:.5f}  star {row['star']:.5f}  "
                f"ratio {row['ratio']:.3f}"
            )
        print(f"\nCO2 at {co2_path:.1f} atmos-cm of path, net of the H2O overlap:")
        print(f"  {'band':>6} {'Abar':>7} {'clear':>7} {'f_solar':>9} {'solar':>9}")
        for name, row in co2_rows.items():
            print(
                f"  {name:>6} {row['mean_absorptance']:7.4f} "
                f"{row['h2o_transmission_in_band']:7.4f} "
                f"{row['flux_fraction']:9.6f} {row['contribution']:9.6f}"
            )
        print(f"  {'total':>6} {'':>7} {'':>7} {'':>9} {co2_sun:9.6f}")
        print("\nCO2 closed form against the integration it is fitted to,")
        print("and against what radmod.f90 runs, which is fitted to the line list:")
        print(f"  {'u (atmos-cm)':>13} {'integrated':>11} {'closed form':>12} {'radmod':>12}")
        for amount in (1.0, 10.0, 100.0, co2_path, 1000.0, 10000.0):
            exact = co2_absorptance(sun, amount, water, bands, scale)[0]
            fitted = float(co2_closed_form(
                amount, co2_fit["a1"], co2_fit["b1"], co2_fit["a2"], co2_fit["b2"]))
            in_model = float(co2_closed_form(
                amount, model_fit["a1"], model_fit["b1"],
                model_fit["a2"], model_fit["b2"]))
            print(f"  {amount:13.1f} {exact:11.6f} {fitted:12.6f} {in_model:12.6f}")

    print(f"\nwrote {out}")

    # THE CHECK THAT CAN FAIL NOW DOES. `EARTH_CO2_SHORTWAVE_W_M2` is declared
    # in this file ahead of any run, the section above says the computation "has
    # to pass with the 2.7 um band in and with it out", and both verdicts were
    # computed, printed as `inside`/`OUTSIDE`, written to the report and then
    # dropped: nothing in this file raised, so every weight it prints was quoted
    # under a bar that could not stop it. Raised AFTER the report is written, so
    # the numbers that failed are on disk to read. world-60x0.
    # BOTH CHECKS THAT CAN FAIL ARE COLLECTED AND RAISED TOGETHER, so that a
    # miss on one does not hide the other's verdict. Raised AFTER the report is
    # written, so the numbers that failed are on disk to read.
    failures = []
    if not ck2["inside"]:
        failures.append(
            f"at w = {ck2['w']} cm this file's Howard reconstruction is "
            f"{ck2['reconstruction_over_eq21']:.4f} times Lacis and Hansen "
            f"Eq. 21 and the HITRAN2020 correlated-k answer for the same "
            f"defined quantity is {ck2['correlated_k_over_eq21']:.4f}, "
            f"{abs(ck2['apart'])*100:.1f}% apart against a tolerance of "
            f"{ck2['tolerance']*100:.0f}%. Two absorption datasets on one "
            "quantity: a miss means the band set, the band intervals or the "
            "spectra have moved here, or that "
            "exoplasim/notes/corrk-cross-check.md's table has, and the level "
            "correction h2o_sw_level rests on the two agreeing. "
            f"{out} carries the numbers.")

    if not (ck["inside"] and ck["inside_dropping_the_2.7um_band"]):
        failures.append(
            f"Earth's near-infrared CO2 shortwave absorption comes out at "
            f"{ck['computed']:.2f} W/m2 with the 2.7 um band and "
            f"{ck['computed_dropping_the_2.7um_band']:.2f} without it, against "
            f"the measured {lo}-{hi} W/m2. Nothing here is tuned to that range, "
            "so a miss means the band set, the path or the overlap treatment is "
            f"wrong, and every weight above is computed the same way. {out} "
            "carries the numbers.")

    if failures:
        raise SystemExit("\n\n".join(failures))


if __name__ == "__main__":
    main()
