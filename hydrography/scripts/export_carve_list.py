#!/usr/bin/env python3
"""Write the carve verdict as a World Orogen `--preserve-basins` list.

Three outcomes per basin, and the middle one is the point:

  carved     The basin overflows, and with enough discharge to cut its sill
             through within the relaxation window. Retain 0, which Orogen treats
             as not preserved, so drainage enforcement carves the outlet as it
             would have without preservation at all.
  preserved  The basin does not overflow. Retain 1, rim intact.
  marginal   The basin overflows but only trickles over. Retain between 0 and 1,
             which cuts a notch at the saddle and tapers it over the divide band,
             producing a through-flowing valley with a residual lake.

Retain is what Orogen cuts with, so it answers a geomorphic question: how much of
the rim survives the water crossing it. Two separate things decide that, and they
are computed separately here.

**Whether the basin overflows** is the balance at spill level,

    Q = runoff * (catchment - area_at_spill) - (E - P) * area_at_spill

which is the same test as `(E - P) / runoff <= catchment / area_at_spill - 1`
wherever runoff is positive, and is defined where it is not. A basin with no
catchment at all still gains water when its own lake surface receives more
precipitation than it evaporates, and it then fills until it spills. Dividing by
runoff loses that case, which is how 228 dry pans came to be carved on one pass
and pinned shut on the next.

**How deeply the outlet is cut** is Q against the depth of the basin. Incision
goes as `K Q^m S^n`, so what an overflow achieves over a relaxation window is a
LENGTH; whether that length empties the basin depends on how far it is from the
spill point down to the floor.

    cut    = coefficient * erodibility * slope**SLOPE_EXPONENT * Q**INCISION_EXPONENT
    retain = clip(1 - cut / depth_at_spill, 0, 1)

That is the fraction of the FINISHED impoundment that survives, in metres,
which is the basis stream power is written in. **It is not the basis Orogen
spends it against**: `buildBasinProtection` sets the allowance to
`(1 - retain)` times the depression's NATURAL relief in the model's
dimensionless elevation parameter, because protection is built before erosion
and the finished depth does not exist yet. `to_natural_relief_basis` converts
before the list is written, and the list and the format both declare the basis.

Three of those four terms are measured per basin and only the coefficient is
free. `erodibility` is the export's own field, a relative stream-power
multiplier normalised to 1 over land, and it is the contrast a channel network
expresses rather than the contrast between intact rock samples; see
`sill_erodibility`. `slope` is the gradient of the outflow channel below the
saddle, relative to the population's, and it is measured because it cannot be
absorbed: the outlet's gradient and the basin's depth are the same relief seen
twice, correlated at r = 0.735 on this build, so folding one into a constant
while dividing by the other counts it twice. See `outlet_gradient`.

`coefficient` is SOLVED, not declared, and `calibrate_coefficient` says why: the
relaxation window is undefined in a generator with no time axis, so it is settled
the way this project settles every other missing-time question, by an
expected-value argument over a stationary population. Earth is one randomly
chosen moment in such a population, and matching the density of its standing
through-flowing impounded basins fixes the number. It is solved on every run
rather than written down, because the count it matches depends on the climate,
so a literal goes stale the moment the verdict moves -- and it had.

**The uncertainty is a third thing, it is kept apart, and it turns out not to
decide anything.** Retain is the larger of the incision value and the margin
below, on the reasoning that either is a reason to leave a rim standing. But
taking the larger lets the margin raise a retain and never lower one, and the
margin is positive exactly when the basin does not overflow, which is exactly
when the incision term is already 1. Both read the sign of the same `Q`, so
`retain` equals `retain_incision` on every basin. The margin is still computed
and still written to the sidecar, because how close a preserved basin sits to
its threshold is worth reading; it is not a second input to the cut.

A basin balances exactly at an evaporation of `E* = P + critical * runoff`. Under
the Penman estimate it evaporates `E_penman`. The fractional margin

    margin = (E_penman - E*) / E_penman
    retain_margin = min(1, margin / TOLERANCE)

says how much open-water evaporation would have to fall before the basin starts
overflowing. A basin needing a 3% change is genuinely marginal and keeps little
of its rim; one needing 40% is comfortably closed and keeps all of it.

TOLERANCE is 0.25, set by the terms that dominate lake evaporation on this world:
the assumed biosphere at 3.7 to 7.1 K, dust at 10 to 20%, and the sub-grid dry
column at up to 9%. Penman's own method error is the smallest of them and is
computed per run rather than quoted, so it is not what sets the scale.

An earlier mapping interpolated the critical index between the two evaporation
estimates. It is retained in the sidecar as `retain_span` but is not used, because
it measures the wrong thing: the land-evaporation end is a deliberately weak lower
bound, so the span is enormous and the critical value lands near its low end
almost regardless of the basin. It put 70% of marginal basins above retain 0.9,
which is to say it called them marginal and then declined to cut them.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

import carve_verdict as cv
from _paths import CONFIG, DATA, PROJECT_ROOT  # noqa: F401
import climatology  # noqa: E402  from lib/, put on sys.path by _paths
from builds import component_data
from orbit import orbital_year_days
from paths import climatology_path, rel, require_configured_grid
from provenance import staged_surface_field
from orogen import Export, LAND
from lake_balance import BasinSet


# The stream-power exponents. `K Q^m S^n`: m is the discharge exponent, 0.5
# being the middle of the 0.4 to 0.6 range detachment-limited bedrock studies
# fit, and n is the slope exponent, 1 because that is Orogen's own -- GRAV-4
# records that n = 1 is baked into the Braun-Willett closed form the generator
# solves rather than being a parameter of it, and the post-hoc 1/g relief
# scaling this project relies on is only correct at that value. Carrying a
# different n here than the generator carves with would be two worlds.
INCISION_EXPONENT = 0.5
SLOPE_EXPONENT = 1.0
SLOPE_EXPONENT_BRACKET = (0.5, 2.0)
KM3_PER_YEAR_TO_M3_PER_S = 1e9

# The Earth calibration, which fixes the coefficient. Measured from HydroLAKES
# v1.0 joined to HydroBASINS level 5: natural lakes (`Lake_type == 1`) deeper
# than 5 m whose pour point sits in a basin with `ENDO == 0`, within 35 degrees
# of the equator, over the land those same level 5 polygons cover in that band.
# `ENDO == 0` selects
# a sill a river actually crosses; the latitude cut keeps out sills that were
# under an ice sheet 20,000 years ago and have had no time to be cut at all.
#
# THE SIZE FLOOR ON THAT SAMPLE IS A DECLARED PARAMETER AND IT IS SWEPT. It used
# to be a single number defended in this comment as "what the mesh resolves",
# and it decides more of the answer than the calibration's own error does: over
# the span the build itself derives, the solved coefficient moves by a factor of
# eight where the Poisson error on the Earth count moves it by 1.8.
# `sweep_size_floor` re-solves at every rung of the ladder below and reports the
# result in the sidecar beside the verdict; `hydrography/notes/retain-fraction.md`
# carries the sweep and the criterion it was read against.
#
# `--measure-earth-floors` regenerates the table from the two bulk sources under
# `hydrography/data/reference/`. The floor in force is the only entry the
# calibration reads; the rest exist so the choice is reported rather than
# asserted.
EARTH_SIZE_FLOOR_KM2 = 1000.0
EARTH_STANDING_BY_FLOOR = {          # measured 2026-08-26, see measure_earth_floors
    10.0: 501,
    20.0: 311,
    50.0: 164,
    100.0: 102,
    200.0: 57,
    500.0: 21,
    1000.0: 15,
    2000.0: 12,
    5000.0: 5,
    10000.0: 2,
}
EARTH_STANDING_BASINS = EARTH_STANDING_BY_FLOOR[EARTH_SIZE_FLOOR_KM2]
# THE DENOMINATOR OF THAT DENSITY, AND IT HAS TO BE THE SAME EARTH AS THE
# NUMERATOR. The count above is lakes whose pour point falls inside a
# HydroBASINS level 5 polygon, so the land it is a density over is the land
# those polygons cover -- summed `SUB_AREA` over the level 5 basins in the band.
# Any other land area makes this a ratio of two different Earths, in the same
# way that a response measured on one build over an amplification measured on
# another is a sensitivity of neither.
#
# It was 78.9 with no producer anywhere in this tree, while the numerator beside
# it had `measure_earth_floors`. `measure_earth_band_land_mkm2` is now that
# producer, `--measure-earth-floors` runs both in one pass over the same
# shapefiles, and it REFUSES when this constant disagrees with what it measures.
# The 78.9 is not what these polygons cover.
#
# Assignment is by each polygon's representative point, and the alternative --
# each polygon's `SUB_AREA` scaled by the share of it inside the band -- agrees
# to 0.06%, so the band edge is not what decides this number. Against the
# Poisson error on a count of 15, 25.8%, the 2.4% between the old figure and
# this one is a tenth of the acknowledged uncertainty; it is corrected because
# the denominator has to be the numerator's Earth, not because it is large.
EARTH_BAND_LAND_MKM2 = 77.05
# Written to two decimals, so a unit in the last place is what agreement means.
EARTH_BAND_LAND_TOLERANCE_MKM2 = 0.005
CALIBRATION_LATITUDE = 35.0
COEFFICIENT_SEARCH = (1.0, 1.0e5)   # m per (m3/s)^0.5, the bisection bracket
# The statistics guard on a rung of that ladder, fixed before the sweep ran. A
# rung is USABLE while the Earth sample keeps at least this many lakes: below
# ten the Poisson fractional error 1/sqrt(N) exceeds 0.32, worse than the 0.258
# the floor in force already carries, and a density stops being a measurement.
# Thin rungs are reported with their counts and marked, never dropped -- a
# reader has to see where the sample runs out.
EARTH_SAMPLE_MINIMUM = 10


def incision_retain(q_km3_per_year, year_s: float, depth_m, erodibility,
                    coefficient: float, slope=None,
                    slope_exponent: float = SLOPE_EXPONENT):
    """Rim surviving the overflow. 1 keeps it, 0 cuts it.

    Takes the overflow at spill level in km3 per Vesper year, the depth from the
    spill point down to the basin floor, the sill's own erodibility, and the
    outlet channel's own gradient relative to the population. A basin that does
    not overflow gets 1 by construction, since Q is then zero.

        cut     = coefficient * erodibility * slope**n * Q**m
        retain  = clip(1 - cut / depth_at_spill, 0, 1)

    **Dividing by the depth is what makes retain a fraction of the impoundment.**
    What the overflow can cut is a length, and whether that empties the basin
    depends on how deep the basin is. An older mapping compared the cut against
    a declared discharge instead and never asked about depth, so it said the
    same thing about a 14 m pan and a 1,420 m trough and made the marginal class
    a razor-thin band around one discharge rather than a property of the terrain.

    **`slope` is measured rather than absorbed, and that is HYD-15.** It used to
    be folded into the coefficient, which asserts that a sill's gradient is
    independent of the depth behind it. It is not: measured over the outflow
    channel below each saddle on this build, `log S` against `log depth` gives
    r = 0.735, with an exponent of 1.04 by ordinary regression and 1.41 by
    reduced major axis. The two are the same relief counted twice. Under
    `S ~ depth` the effective depth exponent in `retain` is `n - 1`, so at
    Orogen's own n = 1 the depth CANCELS -- a different mapping rather than a
    different coefficient, which is why this had to be settled and not
    bracketed.

    What replaces the depth is the outlet gradient itself, and it is worth
    carrying because it is nearly orthogonal to everything else here: against
    discharge it measures r = 0.072, so it is new information rather than a
    proxy for the water. Its spread is wider than depth's, so the marginal class
    stays populated. It is now decided by how steeply a basin's overflow leaves
    rather than by how deep the basin is.

    **`Q**m` IS EVALUATED ON AN ANNUAL DISCHARGE AND CANNOT BE EVALUATED ON
    ANYTHING ELSE, and that is a declared absence rather than an oversight.**
    `m` is under one, so the power is concave and `sqrt(mean Q)` exceeds
    `mean(sqrt Q)`: a basin that delivers its overflow in a season is credited
    here with more cutting than that season's discharge does. The correction is
    not computable, because the per-bin discharge does not exist. `Q` is built
    from the catchment's runoff, and `land_water_ledger.yaml` declares
    `seasonal_phase_of_catchment_delivery` unrepresentable -- annual `P - E` is
    exactly what a cell generated over a cycle, while per-bin clamping would
    count the wet season's supply twice. So there is no second arm to bracket
    against, and inventing one would be a term computed from a quantity that
    cannot carry it. Closing this is the same decision as adopting the climate
    column's `surface_runoff`, which that absence already names.

    **What the absence does NOT prevent is bounding the error, and the bound
    changes what the absence costs.** `seasonal_concavity` carries it in full;
    the two facts it turns on are these. The concavity discount is a factor on
    `cut`, `cut` is linear in the coefficient, and the coefficient is SOLVED
    from this same discharge field on every run -- so a discount common to
    every basin is met by a coefficient of exactly `coefficient / phi` and the
    verdict comes back bit-identical. The level is absorbed and only the
    basin-to-basin spread survives, which is why a bracketed discharge here
    would be a bracket of zero width rather than two carve lists. And the
    spread is bounded by the bin weights alone, without the phase: an overflow
    is never negative, so `phi` lies in `[w_min**(1 - m), 1]` and no basin's cut
    can move more than `J = 1/w_min**(1 - m)` either side of the class boundary.
    Every basin outside that band keeps its verdict under every admissible
    delivery phase; the ones inside it are flagged `seasonal_concavity_movable`
    on the sidecar, one basin at a time, because a count cannot say WHICH.

    The EVAPORATION inside `Q` is a different question and is already handled:
    `cv._INTERVAL_BRACKET` brackets the interval Penman is evaluated over, and
    both ends are run.
    """
    return np.clip(1.0 - cut_over_depth(q_km3_per_year, year_s, depth_m,
                                        erodibility, coefficient, slope=slope,
                                        slope_exponent=slope_exponent),
                   0.0, 1.0)


def cut_over_depth(q_km3_per_year, year_s: float, depth_m, erodibility,
                   coefficient: float, slope=None,
                   slope_exponent: float = SLOPE_EXPONENT):
    """The incision as a fraction of the depression, BEFORE the clip.

    `retain` is `1 - this`, clipped, and the clip is why this is its own
    function rather than a line inside `incision_retain`: a basin cut through
    twice over and one cut through exactly both come back as retain 0, and
    `seasonal_concavity` needs the difference to say how far a basin sits from
    the class boundary. One arithmetic, two readers, no copy to drift.
    """
    q = np.asarray(q_km3_per_year, dtype=float) * KM3_PER_YEAR_TO_M3_PER_S / year_s
    cut = coefficient * np.asarray(erodibility, dtype=float) * np.power(
        np.clip(q, 0.0, None), INCISION_EXPONENT)
    if slope is not None:
        cut = cut * np.power(np.asarray(slope, dtype=float), slope_exponent)
    return cut / np.maximum(np.asarray(depth_m, dtype=float), 1.0)


def to_natural_relief_basis(retain, retained_fraction):
    """Rebase a retain from the finished depression to the one Orogen spends it against.

    **The carve list's retain is a fraction of the depression's NATURAL relief,
    measured in the generator's dimensionless elevation parameter.** That is
    the basis `vendor/orogen/js/basins.js:buildBasinProtection` spends
    `(1 - retain)` against, and it is not a choice: protection is built before
    erosion runs, so the finished depth does not exist yet, and the carve
    operates on model elevations. `vendor/orogen/tools/README.md` declares it
    on the format, and this is where the verdict is put into it.

    `incision_retain` divides by the FINISHED depression in metres, because
    that is the length the overflow has to remove and metres are what stream
    power is written in. The two bases differ by exactly the catalogue's
    `retainedFraction` -- final over natural spill depth, both in model units
    -- so with `f` for the incision as a fraction of the finished depression,

        retain_finished = 1 - f
        retain_natural  = 1 - f * retained_fraction

    Two things that conversion assumes, and both are stated rather than
    absorbed. The fraction of a depression's relief is read as basis-free
    within the one depression, which is the secant of the height curve across
    that depression rather than its local slope; the curve is quartic on land,
    so this is a linearisation over the depression and not an identity. And
    `retained_fraction` is measured on THIS generation's conditioning while
    the allowance is spent on the NEXT one's, so the per-basin value is an
    estimate of a quantity that does not exist yet. It is a property of the
    conditioning rather than of the carve list, which is what makes it
    transferable; the sidecar carries the population spread as the bracket on
    it.

    Both endpoints are verdicts rather than fractions and mean the same thing
    on either basis, so they pass through exactly: retain 1 is an allowance of
    zero and is bit-identical to no allowance, and retain 0 clears `noLower`
    outright. Only the marginal class is rebased, and the marginal class is
    the one the retain machinery exists for.

    Returns (rebased retain, mask of basins the ratio could be applied to).
    A basin whose `retained_fraction` is absent or zero has no finished
    depression to have measured a ratio on; it keeps its finished-basis value
    and is counted.

    notes/audits/basin-catalogue-floor.md, "Preserved is not the same as still
    closed", measures what the mismatch was worth before this conversion:
    a median 2.1x the intended incision with a 5th-to-95th spread of 0.25 to
    8.4, so it was never a scale factor that could be divided out.
    """
    r = np.asarray(retain, dtype=float)
    ratio = np.asarray(retained_fraction, dtype=float)
    usable = np.isfinite(ratio) & (ratio > 0.0)
    rebased = np.where(usable, 1.0 - (1.0 - r) * np.where(usable, ratio, 1.0), r)
    # Endpoints are instructions, not fractions.
    rebased = np.where(r >= 1.0, 1.0, rebased)
    rebased = np.where(r <= 0.0, 0.0, rebased)
    return np.clip(rebased, 0.0, 1.0), usable


def outlet_gradient(basins_path: Path, regions_path: Path, export: Export,
                    steps: int = 10, min_path_m: float = 50e3):
    """Each basin's outflow gradient, relative to the population's own.

    Returns (normalised gradient, geometric mean in m/m, number measured). The
    normalisation is by the geometric mean over the basins where it resolves, so
    the coefficient keeps the meaning its Earth calibration gives it and only
    the spread about the median basin is new. A basin whose outflow path is too
    short to measure gets 1, the population value, rather than a guess.

    **The saddle itself is the wrong place to measure and yields nothing.** A
    saddle is a saddle: the regions either side of it sit at the same elevation
    by construction, and the cross-divide drop on this build has a median of
    0.1 m with half the basins negative. What stream power wants is the gradient
    of the channel the overflow runs DOWN, which starts at `spill_exit_region`
    and is not sub-grid at all: following the drainage ten receiver steps covers
    a median of 135 km, a dozen mesh regions.

    That corrects `hydrography/notes/retain-fraction.md`, which recorded that
    the slope term "is not in the export". The saddle's own slope is not; the
    outflow channel's is, and it is the one the law is written about.

    The estimate is stable in the length it is measured over: the same basins at
    5 steps and at 20 give r = 0.903 in the log with a median ratio of 1.00, so
    it is a property of the outlet and not of the sampling.
    """
    with Dataset(regions_path) as ds:
        receiver = np.asarray(ds["receiver"][:]).astype(np.int64)
    with Dataset(basins_path) as ds:
        start = np.asarray(ds["spill_exit_region"][:]).astype(np.int64)
    radius_m = export.radius_km * 1000.0
    elevation_m = np.asarray(export.field("elevation_km"), dtype=float) * 1000.0
    xyz = np.stack([np.asarray(export.x, dtype=float),
                    np.asarray(export.y, dtype=float),
                    np.asarray(export.z, dtype=float)], axis=1)
    xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)

    here = np.maximum(start, 0)
    top = elevation_m[here].copy()
    at = xyz[here].copy()
    travelled = np.zeros(start.size)
    alive = start >= 0
    for _ in range(steps):
        nxt = np.where(alive, receiver[np.maximum(here, 0)], -1)
        moving = alive & (nxt >= 0)
        if not moving.any():
            break
        j = np.maximum(nxt, 0)
        step_m = radius_m * np.arccos(np.clip((at * xyz[j]).sum(axis=1), -1.0, 1.0))
        travelled = np.where(moving, travelled + step_m, travelled)
        at = np.where(moving[:, None], xyz[j], at)
        here = np.where(moving, nxt, here)
        alive = moving

    drop = top - elevation_m[np.maximum(here, 0)]
    measured = (travelled >= min_path_m) & (drop > 0.0)
    if not measured.any():
        raise SystemExit("no basin outflow path was long enough to measure a "
                         "gradient; the drainage products and the basin "
                         "catalogue are probably from different builds")
    gradient = np.where(measured, drop / np.maximum(travelled, 1.0), np.nan)
    reference = float(np.exp(np.mean(np.log(gradient[measured]))))
    return (np.where(measured, gradient / reference, 1.0), reference,
            int(measured.sum()))


def _standing_solver(q_km3_per_year, year_s: float, depth_m, erodibility,
                     slope, basin_latitude, slope_exponent: float = SLOPE_EXPONENT):
    """`(standing, solve)` over one basin population, so the sweep cannot drift.

    Both the calibration and the size-floor sweep answer the same question at a
    different target, and a second copy of the bisection is how the two would
    come to disagree about what "standing" means. `standing` counts overflowing
    basins inside the calibration band that keep some rim at a given
    coefficient; `solve` inverts it.
    """
    band = np.abs(np.asarray(basin_latitude, dtype=float)) < CALIBRATION_LATITUDE
    overflowing = np.asarray(q_km3_per_year, dtype=float) > 0.0

    def standing(coefficient: float) -> int:
        retain = incision_retain(q_km3_per_year, year_s, depth_m, erodibility,
                                 coefficient, slope=slope,
                                 slope_exponent=slope_exponent)
        return int((overflowing & band & (retain > 0.0)).sum())

    def solve(target: float) -> float:
        lo, hi = COEFFICIENT_SEARCH
        for _ in range(80):
            mid = float(np.sqrt(lo * hi))
            if standing(mid) > target:
                lo = mid            # too many left standing: cut harder
            else:
                hi = mid
        return float(np.sqrt(lo * hi))

    return standing, solve


def calibrate_coefficient(q_km3_per_year, year_s: float, depth_m, erodibility,
                          slope, basin_latitude, land_band_mkm2: float,
                          slope_exponent: float = SLOPE_EXPONENT):
    """The incision coefficient, from Earth's standing-basin density.

    **Calibrated here rather than declared, because the target moves.** The
    coefficient is whatever leaves as many overflowing-but-still-standing basins
    per unit land as Earth has, and how many basins overflow is a property of
    the climate. So a literal goes stale the moment the verdict does, silently,
    and it did: the published 161 was fitted against a verdict taken before
    HYD-13, and against the verdict that replaced it that same 161 leaves 55
    standing basins in the band where the Earth density asks for 34.

    Orogen has no time axis, so the relaxation window this coefficient absorbs
    is UNDEFINED rather than unmeasured, and an expected-value argument over a
    stationary population is this project's standing answer to that. Earth is
    one randomly chosen moment in such a population and so is this terrain;
    `docs/src/reference/no-time-axis.md` carries the argument.

    Returns (coefficient, bracket, target, achieved). The bracket is the Poisson
    error on the Earth count at the floor in force, which is the dominant
    uncertainty WITHIN a floor and is reported rather than hidden. It is not the
    dominant uncertainty overall: `sweep_size_floor` measures the floor's own
    lever, which is larger, and both go into the sidecar.
    """
    standing, solve = _standing_solver(q_km3_per_year, year_s, depth_m,
                                       erodibility, slope, basin_latitude,
                                       slope_exponent)
    per_mkm2 = land_band_mkm2 / EARTH_BAND_LAND_MKM2
    target = EARTH_STANDING_BASINS * per_mkm2
    spread = np.sqrt(EARTH_STANDING_BASINS)
    bracket = tuple(sorted(solve((EARTH_STANDING_BASINS + s) * per_mkm2)
                           for s in (-spread, spread)))
    coefficient = solve(target)
    return coefficient, bracket, target, standing(coefficient)


def sweep_size_floor(q_km3_per_year, year_s: float, depth_m, erodibility,
                     slope, basin_latitude, land_band_mkm2: float,
                     land_cell_area_km2, area_at_spill_km2,
                     slope_exponent: float = SLOPE_EXPONENT,
                     table: dict | None = None) -> dict:
    """What the size floor on the Earth sample is worth, against what it is not.

    **The floor decides more of the coefficient than the calibration's own error
    does, and it used to be defended in a comment rather than reported.** The
    calibration matches a DENSITY of standing through-flowing impounded basins,
    and a density is a count over an area at a size class. Move the size class
    and the count moves by two orders of magnitude while the area does not, so
    the target moves with it and the coefficient follows. Nothing downstream can
    see that from a single solved number, which is what this puts in the sidecar.

    **The span is derived from the build rather than chosen**, which is the whole
    difference between a sweep and a defence. Both bounds are measured in the
    same run as the verdict, so they move when the mesh or the terrain does:

    - the LOWER bound is the median area of a land mesh cell. Below it an Earth
      lake counted has no representable counterpart at all, because a depression
      smaller than one cell does not exist in the generator's output.
    - the UPPER bound is the median area at spill of the overflowing basins.
      Above it more than half of the population the coefficient is being solved
      for sits outside the class the Earth sample stands for, so the Earth
      density has stopped being a density of the same object.

    The ladder is 1-2-5 per decade over the range HydroLAKES supports, and rungs
    outside the span are reported with their counts rather than dropped: the
    reader has to be able to see the whole lever, including the part the span
    excludes and why.

    **The criterion is fixed before the sweep runs, and it compares the lever
    with the noise.** `L` is the ratio of solved coefficients over the rungs that
    are both inside the span and above `EARTH_SAMPLE_MINIMUM`; `P` is the ratio
    of the Poisson bracket on the Earth count at the floor in force, solved on
    this same population. `L > P` is DECISIVE -- the floor moves the answer by
    more than the count's own error, so it has to be declared and swept. `L <= P`
    is SUBORDINATE -- it sits inside an uncertainty already reported and adds
    nothing quantitative. Two classes, no gap and no overlap. The same call is
    made a second time on the marginal-basin count, which is the landform class
    the coefficient decides the size of. `failure-modes.md` class 34 is what the
    comparison exists for: a lever smaller than the instrument's own scatter is
    not a finding.

    The counts reported per rung are on the FINISHED-depression basis, before
    `to_natural_relief_basis`. That is deliberate and it is the basis the
    coefficient is solved on, so the sweep is internally consistent; the verdict
    the list carries is the rebased one and the sidecar's own counts are that.
    """
    standing, solve = _standing_solver(q_km3_per_year, year_s, depth_m,
                                       erodibility, slope, basin_latitude,
                                       slope_exponent)
    per_mkm2 = land_band_mkm2 / EARTH_BAND_LAND_MKM2
    overflowing = np.asarray(q_km3_per_year, dtype=float) > 0.0
    lower = float(np.median(np.asarray(land_cell_area_km2, dtype=float)))
    spill = np.asarray(area_at_spill_km2, dtype=float)[overflowing]
    upper = float(np.median(spill)) if spill.size else float("inf")

    def classes(coefficient: float) -> tuple[int, int, int]:
        retain = incision_retain(q_km3_per_year, year_s, depth_m, erodibility,
                                 coefficient, slope=slope,
                                 slope_exponent=slope_exponent)
        return (int((retain <= 0.0).sum()),
                int(((retain > 0.0) & (retain < 1.0)).sum()),
                int((retain >= 1.0).sum()))

    table = EARTH_STANDING_BY_FLOOR if table is None else table
    rungs = []
    exact = {}          # unrounded, because a solved value sits ON a step
    for floor in sorted(table):
        n = table[floor]
        target = n * per_mkm2
        coefficient = solve(target)
        spread = np.sqrt(n)
        lo, hi = sorted(solve(max(n + s, 0.0) * per_mkm2) for s in (-spread, spread))
        exact[floor] = (coefficient, lo, hi)
        carve, marginal, preserve = classes(coefficient)
        rungs.append({
            "floor_km2": floor,
            "earth_lakes": n,
            "target_standing_basins": round(target, 2),
            "coefficient": round(coefficient, 3),
            "coefficient_poisson_bracket": [round(lo, 3), round(hi, 3)],
            "standing": standing(coefficient),
            "carve": carve,
            "marginal": marginal,
            "preserve": preserve,
            "inside_span": bool(lower <= floor <= upper),
            "usable": bool(n >= EARTH_SAMPLE_MINIMUM),
        })

    used = [r for r in rungs if r["inside_span"] and r["usable"]]
    # Unrounded throughout: `solve` returns a value sitting exactly on a step of
    # the standing count, so rounding it before counting classes moves a basin
    # across the boundary and the reported noise is then off by one.
    _, poisson_lo, poisson_hi = exact[EARTH_SIZE_FLOOR_KM2]
    ratio = (poisson_hi / poisson_lo) if poisson_lo > 0 else float("inf")
    marginal_poisson = sorted(classes(b)[1] for b in (poisson_lo, poisson_hi))
    verdict = {}
    if used:
        cs = [exact[r["floor_km2"]][0] for r in used]
        ms = [r["marginal"] for r in used]
        lever = max(cs) / min(cs)
        lever_marginal = (max(ms) / min(ms)) if min(ms) > 0 else float("inf")
        ratio_marginal = ((marginal_poisson[1] / marginal_poisson[0])
                          if marginal_poisson[0] > 0 else float("inf"))
        verdict = {
            "coefficient_over_span": [round(min(cs), 3), round(max(cs), 3)],
            "coefficient_lever": round(lever, 3),
            "coefficient_poisson_ratio": round(ratio, 3),
            "coefficient_call": "DECISIVE" if lever > ratio else "SUBORDINATE",
            "marginal_over_span": [min(ms), max(ms)],
            "marginal_lever": round(lever_marginal, 3),
            "marginal_poisson": marginal_poisson,
            "marginal_poisson_ratio": round(ratio_marginal, 3),
            "marginal_call": ("DECISIVE" if lever_marginal > ratio_marginal
                              else "SUBORDINATE"),
        }
    return {
        "floor_in_force_km2": EARTH_SIZE_FLOOR_KM2,
        "span_km2": [round(lower, 2), round(upper, 2)],
        "span_derivation": "lower = median land mesh cell area, below which a "
            "counted Earth lake has no representable counterpart; upper = "
            "median area at spill of the overflowing basins, above which most "
            "of the population being solved for is outside the class the Earth "
            "sample stands for. Both measured from this build in this run.",
        "sample_minimum_lakes": EARTH_SAMPLE_MINIMUM,
        "criterion": "L = ratio of solved coefficients over the rungs that are "
            "both inside the span and at or above sample_minimum_lakes; P = "
            "ratio of the Poisson bracket on the Earth count at the floor in "
            "force, solved on this population. L > P is DECISIVE, L <= P is "
            "SUBORDINATE. Fixed before the sweep was run.",
        "basis": "finished-depression retain, the basis the coefficient is "
            "solved on. The verdict counts elsewhere in this sidecar are the "
            "rebased ones and are not comparable rung for rung.",
        "verdict": verdict,
        "rungs": rungs,
    }


def seasonal_concavity(q_km3_per_year, year_s: float, depth_m, erodibility,
                       slope, coefficient: float, bin_weights,
                       coefficient_bracket, slope_exponent: float = SLOPE_EXPONENT):
    """What `Q**m` on an annual discharge is worth, bounded without the phase.

    `Q` is annual and cannot be otherwise here: `land_water_ledger.yaml`
    declares `seasonal_phase_of_catchment_delivery` an ABSENCE, so there is no
    per-bin discharge to evaluate and no second arm to run. `incision_retain`
    declares that at the site. What follows is the bound that absence still
    permits, and it needs no phase at all -- Jensen's gap for a concave
    function depends on how uneven the discharge is across bins, never on which
    bin the water arrives in.

    **THE LEVEL IS ABSORBED BY THE CALIBRATION, EXACTLY.** Write the truth as
    `mean_k sqrt(Q_k) = phi * sqrt(Q_annual)`, so `phi` is the whole of the
    concavity error on one basin. `cut` is LINEAR in the coefficient and
    linear in `Q**m`, and `calibrate_coefficient` solves the coefficient on
    every run against Earth's standing-basin density from this same discharge
    field. So a `phi` common to every basin is met by a coefficient of exactly
    `coefficient / phi` and every retain comes back bit-identical: not close,
    identical, and `_selftest` checks it. A bracket over the common level would
    therefore have ZERO WIDTH, which is why this is not two carve lists and not
    two generations. Only the basin-to-basin SPREAD of `phi` reaches the
    verdict.

    **THE SPREAD IS BOUNDED BY THE BIN WEIGHTS ALONE.** Two facts and nothing
    else: an overflow is never negative, and the cycle has the climatology's
    own bins with the weights `lib/climatology.py` gives them. A discharge
    spread evenly gives `phi = 1`; the most uneven admissible discharge puts
    the whole year's overflow in the lightest bin and gives
    `phi = w_min**(1 - m)`. So `phi` lies in `[w_min**(1 - m), 1]`, the ratio
    between two basins' `phi` lies within `J = 1/w_min**(1 - m)`, and with the
    population level pinned by the calibration a basin's cut can move by at
    most `J` either side of the class boundary.

    That makes the categorical question answerable. A basin whose annual
    `cut / depth` sits outside `(1/J, J)` KEEPS ITS VERDICT under every
    admissible delivery phase; one inside that band does not, and is reported
    per basin as `seasonal_concavity_movable`. The band is a statement about
    what this pipeline cannot see, not a claim that those basins are wrong.

    **The call uses `sweep_size_floor`'s criterion, unchanged**, because it is
    the same comparison against the same instrument: `J` against `P`, the ratio
    of the Poisson bracket on the Earth count. `J > P` is DECISIVE and has to be
    declared and carried; `J <= P` is SUBORDINATE and sits inside an uncertainty
    already reported.

    Counts are on the FINISHED-depression basis, for `sweep_size_floor`'s
    reason: it is the basis the coefficient is solved on. The class boundary is
    the same on either basis -- `to_natural_relief_basis` passes both endpoints
    through -- so a movable basin is movable on the list as written.

    What would tighten it is a per-bin lake surface term. The overflow is
    catchment delivery, whose phase is the absence, PLUS the lake's own surface
    flux, whose phase is representable and is already carried in
    `surface_water.nc`. Bounding the two separately is strictly narrower than
    bounding their sum, and needs a per-bin open-water evaporation that
    `carve_verdict.bin_mean_open_water` reduces before returning.

    Returns (report, movable mask).
    """
    w = np.asarray(bin_weights, dtype=float)
    if w.ndim != 1 or w.size < 2 or not np.isclose(w.sum(), 1.0) or (w <= 0).any():
        raise SystemExit(
            f"bin weights are {w.size} values summing to {w.sum():.6g}; the "
            "concavity bound is read off the LIGHTEST bin and means nothing "
            "unless they are the climatology's own weights and sum to one")
    # `w_min ** (1 - m)`, not `sqrt(w_min)`: the floor is a property of the
    # exponent as well as of the weighting, and hardcoding the square root would
    # go silently wrong if `m` ever moved off 0.5.
    phi_floor = float(w.min() ** (1.0 - INCISION_EXPONENT))
    factor = 1.0 / phi_floor

    x = cut_over_depth(q_km3_per_year, year_s, depth_m, erodibility,
                       coefficient, slope=slope, slope_exponent=slope_exponent)

    overflows = np.asarray(q_km3_per_year, dtype=float) > 0.0
    cuts = overflows & (x >= 1.0)
    stands = overflows & (x < 1.0)
    movable = overflows & (x > phi_floor) & (x < factor)
    poisson = float(coefficient_bracket[1]) / float(coefficient_bracket[0])
    return {
        "bins": int(w.size),
        "lightest_bin_weight": round(float(w.min()), 6),
        "phi_range": [round(phi_floor, 4), 1.0],
        "factor_J": round(factor, 4),
        "coefficient_poisson_ratio_P": round(poisson, 3),
        "call": "DECISIVE" if factor > poisson else "SUBORDINATE",
        "criterion": "J against P by sweep_size_floor's criterion, unchanged: "
            "J > P is DECISIVE and has to be declared and carried per basin, "
            "J <= P is SUBORDINATE and sits inside an uncertainty already "
            "reported. J is the widest ratio between two basins' seasonal "
            "concentration; P is the Poisson bracket on the Earth count.",
        "level_is_absorbed": "a concavity discount common to every basin is met "
            "by a coefficient of exactly coefficient/phi and returns a "
            "bit-identical verdict, because calibrate_coefficient solves the "
            "coefficient on every run from this same discharge field. Only the "
            "basin-to-basin spread reaches the verdict, and a bracket over the "
            "common level would have zero width.",
        "assumptions": ["an overflow is never negative in any bin",
                        "the cycle is the climatology's own bins at "
                        "lib/climatology.py's weights",
                        "nothing about which bin the water arrives in"],
        "basis": "finished-depression cut over depth, the basis the coefficient "
            "is solved on. The class boundary is the same on the rebased basis, "
            "because to_natural_relief_basis passes both endpoints through.",
        "overflowing": int(overflows.sum()),
        "cuts": int(cuts.sum()),
        "stands": int(stands.sum()),
        "movable": int(movable.sum()),
        "movable_of_the_cut": int((movable & cuts).sum()),
        "movable_of_the_standing": int((movable & stands).sum()),
        "held_cut": int((overflows & (x >= factor)).sum()),
        "held_standing": int((overflows & (x <= phi_floor)).sum()),
    }, movable


def measure_earth_floors(reference: Path, floors=None) -> dict:
    """Regenerate `EARTH_STANDING_BY_FLOOR` from the two bulk sources.

    The code that measured a constant belongs beside the constant, or the
    constant is unreproducible the first time anyone doubts it. This is the
    query `EARTH_STANDING_BY_FLOOR` was measured with, and running it is how a
    reader checks the table rather than trusting a date in a comment.

    Selection, declared and unchanged from the original measurement: HydroLAKES
    v1.0 natural lakes (`Lake_type == 1`) with `Depth_avg` above 5 m and
    `|Pour_lat|` under `CALIBRATION_LATITUDE`, whose pour point falls inside a
    HydroBASINS level 5 polygon with `ENDO == 0`. Counted at each floor on
    `Lake_area`, strictly above.

    geopandas is imported here rather than at module scope. It is in
    `requirements.txt` for HYD-4, but the calibration itself reads only the
    table, so the verdict must not acquire a dependency on the bulk sources
    being present -- a worktree without them still has to be able to run.
    """
    import geopandas as gpd     # noqa: PLC0415  see docstring
    import pandas as pd         # noqa: PLC0415
    import pyogrio              # noqa: PLC0415

    floors = sorted(EARTH_STANDING_BY_FLOOR) if floors is None else sorted(floors)
    lakes_path = (reference / "HydroLAKES_polys_v10_shp"
                  / "HydroLAKES_polys_v10.shp")
    basin_paths = sorted(reference.glob("hybas_*_lev05_v1c.shp"))
    if not lakes_path.exists() or not basin_paths:
        raise SystemExit(
            f"{lakes_path} or the HydroBASINS level 5 shapefiles are missing "
            f"from {reference}. This mode reads the BULK sources; the "
            "calibration itself does not, and reads EARTH_STANDING_BY_FLOOR.")

    lakes = pyogrio.read_dataframe(
        lakes_path, read_geometry=False,
        columns=["Hylak_id", "Lake_name", "Lake_type", "Lake_area",
                 "Depth_avg", "Pour_long", "Pour_lat"])
    sel = lakes[(lakes.Lake_type == 1) & (lakes.Depth_avg > 5.0)
                & (lakes.Pour_lat.abs() < CALIBRATION_LATITUDE)
                & (lakes.Lake_area > min(floors) / 10.0)].copy()
    points = gpd.GeoDataFrame(
        sel, crs="EPSG:4326",
        geometry=gpd.points_from_xy(sel.Pour_long, sel.Pour_lat))
    basins = gpd.GeoDataFrame(
        pd.concat([pyogrio.read_dataframe(p, columns=["ENDO", "SUB_AREA"])
                   for p in basin_paths], ignore_index=True),
        geometry="geometry", crs="EPSG:4326")
    joined = gpd.sjoin(points, basins[["ENDO", "geometry"]], how="left",
                       predicate="within").drop_duplicates(subset="Hylak_id")
    area = joined[joined.ENDO == 0].Lake_area.values

    # The denominator, from the SAME polygons, in the same pass. Measuring the
    # numerator without it is what left `EARTH_BAND_LAND_MKM2` with no producer.
    band = measure_earth_band_land_mkm2(basins)
    declared = EARTH_BAND_LAND_MKM2
    if abs(band["band_land_mkm2"] - declared) > EARTH_BAND_LAND_TOLERANCE_MKM2:
        raise SystemExit(
            f"EARTH_BAND_LAND_MKM2 is {declared} and these polygons cover "
            f"{band['band_land_mkm2']:.3f} Mkm2 within "
            f"{CALIBRATION_LATITUDE:g} degrees. The count above is a density "
            "over exactly this land, so the two move together or the density "
            "is a ratio of two different Earths. Update the constant.")
    return {
        "selection": "HydroLAKES v1.0 Lake_type == 1, Depth_avg > 5 m, "
            f"|Pour_lat| < {CALIBRATION_LATITUDE:g} deg, pour point inside a "
            "HydroBASINS level 5 polygon with ENDO == 0. Counted strictly above "
            "each floor on Lake_area.",
        "unmatched_pour_points": int(joined.ENDO.isna().sum()),
        "counts": {float(f): int((area > f).sum()) for f in floors},
        "band_land": band,
    }


def measure_earth_band_land_mkm2(basins) -> dict:
    """Land the HydroBASINS level 5 polygons cover within the calibration band.

    The denominator of the standing-lake density, measured from the polygons the
    numerator's lakes are assigned to, so the two are one Earth. `SUB_AREA` is
    HydroBASINS' own area for each level 5 basin, in km2.

    Assignment is by each polygon's REPRESENTATIVE POINT rather than by clipping
    it at the band edge, because clipping in a geographic CRS weights a degree
    of longitude equally at every latitude and the areas it returns are not
    areas. The two answers agree to 0.06% here, which is what says the band edge
    does not decide this number; if they ever stop agreeing, the clip is the one
    to distrust and a projected re-clip is the repair.
    """
    import numpy as np_          # noqa: PLC0415  local, like geopandas above
    if "SUB_AREA" not in basins.columns:
        raise SystemExit(
            "the HydroBASINS level 5 frame carries no SUB_AREA column, so the "
            "band land cannot be measured from the polygons the lakes were "
            "assigned to. Read the shapefiles with SUB_AREA included.")
    latitude = basins.geometry.representative_point().y.to_numpy()
    sub_area = basins.SUB_AREA.to_numpy(dtype=float)
    inside = np_.abs(latitude) < CALIBRATION_LATITUDE
    return {
        "band_land_mkm2": float(sub_area[inside].sum()) / 1.0e6,
        "global_land_mkm2": float(sub_area.sum()) / 1.0e6,
        "polygons_in_band": int(inside.sum()),
        "polygons": int(sub_area.size),
        "assignment": "each level 5 polygon's representative point, against "
                      f"|lat| < {CALIBRATION_LATITUDE:g} deg",
        "note": "HydroBASINS omits Antarctica, which is outside this band and "
                "so cannot reach the sum.",
    }


def sill_erodibility(basins_path: Path, terrain_hash: str) -> np.ndarray:
    """The rock at each basin's outlet, as a relative stream-power multiplier.

    The export's `erodibility` is exactly this quantity and says so: a relative
    stream-power multiplier from the exposed rock. So it enters as `K` does, and
    Q_full divides by it: a soft sill is cut through by less water than a hard
    one. Taken AS SHIPPED and never renormalised: Orogen normalises the land mean
    to 1 when it builds the field, which is before cover stripping, so the
    delivered mean is off 1 by a few percent in a direction that varies by
    planet. Nothing here needs it to be 1, because `calibrate_coefficient` fits
    the leading coefficient against Earth and absorbs any constant factor.

    **It is the expressed contrast, not the intact-rock one, and that is the
    number a landscape model wants.** Stock and Montgomery (1999) measure K
    across five orders of magnitude between lithologies, but Zondervan (2020)
    measures the contrast a real channel network expresses at about 4x, because
    channels adjust width and slope in response to the rock they are in. This
    field spans under 4x, which is the right order; using the intact-rock spread
    would make the rock the only term that mattered and the discharge decorative.

    The saddle is smaller than a mesh cell, so the two regions either side of it
    bracket the rock being cut rather than naming it. Their geometric mean is
    what is used, that being the right average for a multiplicative factor, and
    the two sides agree only loosely: correlation 0.34 on the current build.
    """
    export = Export()
    if export.terrain_hash != terrain_hash:
        raise SystemExit(
            f"basins.nc was built from terrain {terrain_hash[:16]} and the "
            f"configured export is {export.terrain_hash[:16]}. The sill lookup "
            "is by region index, which does not survive a terrain change."
        )
    try:
        ero = export.field("erodibility")
    except Exception:
        print("note: the export carries no erodibility field; sills are uniform")
        with Dataset(basins_path) as ds:
            return np.ones(ds.dimensions["basin"].size)

    with Dataset(basins_path) as ds:
        inner = np.asarray(ds["spill_region"][:])
        outer = np.asarray(ds["spill_exit_region"][:])
    if inner.size and max(int(inner.max()), int(outer.max())) >= ero.size:
        raise SystemExit("spill regions index past the export's mesh")

    # -1 marks a basin whose saddle was not resolved on one side or either.
    # Fall back through the side that exists to the land mean, which is 1.
    a = np.where(inner >= 0, ero[np.maximum(inner, 0)], np.nan)
    b = np.where(outer >= 0, ero[np.maximum(outer, 0)], np.nan)
    with np.errstate(invalid="ignore"):
        both = np.sqrt(a * b)
    out = np.where(np.isfinite(both), both,
                   np.where(np.isfinite(a), a, np.where(np.isfinite(b), b, 1.0)))
    return np.clip(out, 1e-3, None)



# The fractional evaporation change that would flip a basin, and the width of the
# `marginal` band the margin term produces. Module scope because `climate_terms`
# applies it and `main` records it in the sidecar, and those must be one number.
TOLERANCE = 0.25

def climate_terms(clim_path, args, config, basins, *, interval="bin_mean"):
    """Everything the carve criterion reads from ONE climate.

    Factored out of `main` because `docs/src/pipeline/loops.md` evaluates this same
    criterion at TWO bounding climates -- the warm vegetated end and the cold
    bare-rock end -- and carves only the intersection. Running one function
    twice is the point rather than an implementation detail: two hand-written
    evaluations of a criterion are two formulations that can disagree, and
    `docs/src/practice/failure-modes.md` class 17 is about exactly that. The arms differ in
    their climatology and in nothing else.

    Everything here is climate-dependent. What is NOT here, deliberately, is
    the incision coefficient: it is calibrated against Earth from the discharge
    field, so calibrating it per arm would let each arm move its own yardstick
    and the two verdicts would no longer be comparable. `main` calibrates once
    on the primary arm and applies that coefficient to both.
    """
    # `clim_path`, NOT `args.climatology`. This function exists to be run once
    # per arm of the bracket and it read the primary climatology both times, so
    # the cold arm was the warm arm and the two-climate bracket was one climate
    # compared with itself. What that produces is not an error: it is a bracket
    # reporting zero width, which reads as the arms agreeing completely. `main`
    # now refuses two arms that come out bit-identical.
    with Dataset(clim_path) as ds:
        am = cv.annual_mean
        pr, evap, mrro = am(ds, "pr"), -am(ds, "evap"), am(ds, "mrro")
        rss, rls = am(ds, "rss"), am(ds, "rls")
        diurnal = am(ds, "maxt") - am(ds, "mint")
        lsm = am(ds, "lsm")
    t_air, q_air, wind, p_air = cv.reference_level_air(clim_path)

    # Through the one door, not by rebuilding the path from `model.resolution`:
    # that path is keyed by the RUNG alone and `surface_albedo` rewrites it per
    # BUILD. `--for-build` declares the deliberate cross-build read this export
    # makes once the carved build is staged, by NAMING the build it means.
    # world-z7bu, CLAUDE.md rule 5.
    staged_albedo = staged_surface_field(174, config, for_build=args.for_build)
    land_albedo = cv.read_sra_field(PROJECT_ROOT / staged_albedo["path"],
                                    *p_air.shape)
    # ONE ESTIMATOR WITH carve_verdict.py, evaluated over the interval this arm
    # was asked for. `cv._INTERVAL_BRACKET` carries the argument for why the
    # interval is a bracket: Penman's zero clamps are a real rectification an
    # annual mean cancels away, and against that Penman carries no heat storage
    # so per bin it charges the bright season for energy the water column gave
    # back in the dark one. Two hand-written evaluations of one criterion are
    # two formulations that can disagree, which is this function's whole reason
    # for existing, so it calls what carve_verdict.py calls rather than
    # rebuilding it.
    gravity = float(config["planet"]["gravity_m_s2"])
    if interval == "bin_mean":
        penman_raw = cv.bin_mean_open_water(clim_path, land_albedo, gravity,
                                            cfg=config)
    elif interval == "annual":
        penman_raw = cv.penman_open_water(
            t_air, q_air, wind, p_air, rss, rls, land_albedo, gravity,
            diurnal_range=diurnal, cfg=config)
    else:
        raise ValueError(
            f"interval {interval!r} is neither 'bin_mean' nor 'annual'; those "
            "are the two ends of cv._INTERVAL_BRACKET and there is no third")
    # Computed here rather than quoted. The sidecar carried a hardcoded 1.017
    # for a day after `carve_verdict.py` started computing this, which is the
    # same defect one file over: a validation that cannot move is not one.
    ocean_validation = cv.validate_over_ocean(penman_raw, evap, lsm)
    penman_error_pct = abs(ocean_validation["ratio"] - 1.0) * 100.0
    # No floor at the land rate: see the long note in `carve_verdict.py`. A
    # smooth lake in a rough wet landscape evaporates less than the ground
    # around it, and clamping that away decided 73% of the overflowing set.
    penman = penman_raw

    runoff_field = mrro if args.runoff_source == "mrro" else (pr - evap)
    cv.require_index_alignment(args.coupling, lsm)
    means, _ = cv.basin_means(args.coupling,
                              {"pr": pr, "wet": evap, "pen": penman,
                               "ro": runoff_field},
                              basins.n)
    # Orbital period varies with flux, so take it from the config rather than
    # hardcoding. The literal here was 189.6145 d, the 0.90-flux year, while this
    # baseline runs at 0.96 and 180.655 d. It cancels out of the aridity index
    # and the equilibrium lake area, both ratios, but it was making the reported
    # runoff depth 5% high.
    year_s = orbital_year_days(config) * 86400.0
    # Clamped at zero for the same reason carve_verdict.py clamps it: a
    # catchment losing more to evaporation than it receives delivers nothing,
    # not a negative amount. It keeps the reported depth physical, and the
    # discharge below reads the same clamped field.
    runoff = np.maximum(means["ro"], 0.0) * year_s / 1000.0
    precip = means["pr"] * year_s / 1000.0
    crit = basins.catchment_km2 / np.maximum(basins.area_at_spill_km2, 1e-9) - 1.0

    # Reported, not decided on. The index is undefined where a basin has no
    # catchment runoff, which is exactly the set the discharge form handles.
    def index(field):
        e = means[field] * year_s / 1000.0
        return np.where(runoff > 0, (e - precip) / np.where(runoff > 0, runoff, 1.0),
                        np.inf)

    idx_wet, idx_pen = index("wet"), index("pen")

    # The verdict, as the water that has to leave at spill level. Runoff is
    # generated over the dry catchment; the lake surface itself gains P and
    # loses E, which is the term the ratio form divides away.
    dry_km2 = np.maximum(basins.catchment_km2 - basins.area_at_spill_km2, 0.0)

    def discharge(field):
        e = means[field] * year_s / 1000.0
        return runoff * dry_km2 - (e - precip) * basins.area_at_spill_km2

    q_wet, q_pen = discharge("wet"), discharge("pen")
    # THE TWO ESTIMATES DO NOT NEST, and this used to assume they did.
    #
    # `q_wet >= q_pen everywhere` held only because Penman was floored at the
    # model's land evaporation. Without that floor the ordering reverses wherever
    # the ground is wet, because land here is several times aerodynamically
    # rougher than open water and evaporates more than a lake would; see the
    # note in `carve_verdict.py`. So a basin can overflow under Penman and not
    # under the land rate, the three-way partition below is not a partition, and
    # the assertion that it summed to `basins.n` failed the moment the floor
    # went.
    #
    # What is reported instead is the primary verdict and the disagreement,
    # which is what the two estimates can honestly say. The verdict itself has
    # always followed `retain`, not this triple.
    carved = q_pen > 0.0
    overflows_wet = q_wet > 0.0
    disputed = carved != overflows_wet

    # The fractional evaporation change that would flip a basin, reported and
    # never decisive: `retain` is the LARGER of this and the incision value, and
    # the margin is positive exactly when the basin does not overflow, which is
    # exactly when the incision value is already 1.
    #
    # 0.25 stands, and HYD-12 asked whether it should now that Penman's own
    # error is measured rather than hardcoded. It should. The tolerance is set
    # by the uncertainty that dominates, and Penman's method error is not it:
    # against the model over ocean, in the configuration a lake is actually in,
    # the estimate now lands within a few percent, and the terms above it are
    # the assumed biosphere at 3.7 to 7.1 K, dust at 10 to 20% of lake
    # evaporation, and the sub-grid dry column at up to 9%. A quarter brackets
    # all three; a tolerance tied to the method's own error would have been the
    # smallest term setting the scale for the largest.
    e_pen = means["pen"] * year_s / 1000.0
    # The third of the three catchment-mean water-balance terms, in the same
    # km/year depth as `runoff_km_per_year`. All three are recorded per basin in
    # the sidecar, and the reason is that the criterion is RE-EVALUATED
    # downstream: `scripts/error_budget.py` runs the same test under uniform
    # perturbations to price a budget item in basins. A consumer that instead
    # solves for P and E by inverting the recorded aridity index and evaporation
    # margin is tied to whatever algebra this file used on the day it ran, and
    # that inversion existed and broke the first time the criterion moved. Three
    # floats per basin cost nothing and cannot go stale: if the verdict changes
    # what it means by E, it changes what it writes here.
    e_wet = means["wet"] * year_s / 1000.0
    e_balance = precip + crit * runoff              # evaporation that exactly balances
    with np.errstate(divide="ignore", invalid="ignore"):
        margin = np.where(e_pen > 0, (e_pen - e_balance) / np.where(e_pen > 0, e_pen, 1.0), 1.0)
    retain_margin = np.where(carved, 0.0, np.clip(margin / TOLERANCE, 0.0, 1.0))
    return {"staged_background_albedo": staged_albedo,
            "interval": interval,
            "mrro": mrro, "means": means, "year_s": year_s, "crit": crit, "runoff": runoff, "precip": precip, "idx_pen": idx_pen, "idx_wet": idx_wet, "q_pen": q_pen, "q_wet": q_wet, "carved": carved, "overflows_wet": overflows_wet, "disputed": disputed, "e_pen": e_pen, "e_wet": e_wet, "margin": margin, "retain_margin": retain_margin, "ocean_validation": ocean_validation, "penman_error_pct": penman_error_pct}
def _selftest() -> int:
    """The basis conversion, against identities rather than against outcomes.

    Every check has a right answer and can fail. The first is the one the
    conversion exists for: the physical incision Orogen is PERMITTED, once the
    allowance has been spent against the natural relief, must equal the
    incision hydrography INTENDED against the finished depression. Before the
    conversion that identity failed by a median factor of 2.1 on
    `precarve-craton-10m`, with a 5th-to-95th spread of 0.25 to 8.4.
    """
    problems: list[str] = []
    n_checks = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal n_checks
        n_checks += 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    rng = np.random.default_rng(20260825)
    # A ratio population spanning both sides of 1: the conditioning takes relief
    # off most basins and puts it on about a sixth.
    ratio = np.concatenate([rng.uniform(0.2, 1.0, 400),
                            rng.uniform(1.0, 2.5, 100)])
    natural_depth = rng.uniform(0.01, 1.0, ratio.size)   # model units
    final_depth = ratio * natural_depth
    intended = rng.uniform(0.02, 0.98, ratio.size)       # retain, finished basis
    sent, usable = to_natural_relief_basis(intended, ratio)

    # THE IDENTITY. The allowance is spent as (1 - retain) * natural depth; the
    # cut intended is (1 - retain_finished) * finished depth. Same length.
    permitted = (1.0 - sent) * natural_depth
    wanted = (1.0 - intended) * final_depth
    expressible = sent > 0.0
    err = np.abs(permitted - wanted)[expressible] / np.maximum(
        wanted[expressible], 1e-15)
    check("the permitted incision equals the intended incision",
          bool(usable.all()) and float(err.max()) < 1e-12,
          f"max relative error {float(err.max()):.3g}")

    # The unconverted crossing must FAIL that identity, or the check above is
    # measuring nothing. It is the defect this conversion removes.
    naive = np.abs((1.0 - intended) * natural_depth - wanted) / wanted
    check("the unconverted crossing fails it, so the check is not vacuous",
          float(np.median(naive)) > 0.1,
          f"median relative error {float(np.median(naive)):.3g}")

    # Endpoints are instructions, not fractions, and mean the same on either
    # basis. A conversion that touched them would turn every carve into a
    # partial preservation.
    for value, name in ((0.0, "carve"), (1.0, "preserve")):
        out, _ = to_natural_relief_basis(np.full(ratio.size, value), ratio)
        check(f"retain {value:.0f} ({name}) crosses exactly",
              bool((out == value).all()),
              f"{int((out != value).sum())} entries moved")

    # Order is meaning: a verdict that keeps more rim must send a larger retain.
    order = np.argsort(intended)
    one_ratio = np.full(ratio.size, 0.8)
    mono, _ = to_natural_relief_basis(intended[order], one_ratio)
    check("a larger retain stays a larger retain at one ratio",
          bool((np.diff(mono) >= -1e-15).all()), "the mapping is not monotone")

    # A basin the conditioning left with no depression has no ratio to convert
    # with, and must not be silently converted with a zero.
    out, ok = to_natural_relief_basis(np.array([0.4, 0.4, 0.4]),
                                      np.array([np.nan, 0.0, 0.5]))
    check("a basin with no measurable ratio keeps its finished-basis value",
          (not ok[0]) and (not ok[1]) and bool(ok[2])
          and out[0] == 0.4 and out[1] == 0.4 and abs(out[2] - 0.7) < 1e-12,
          f"got {out.tolist()} usable {ok.tolist()}")

    # Saturation is honest, not a clip that hides an error: an intended cut
    # deeper than the whole natural relief cannot be an allowance, and carve is
    # the closest instruction Orogen has.
    out, _ = to_natural_relief_basis(np.array([0.1]), np.array([2.0]))
    check("an intended cut deeper than the natural relief saturates to carve",
          float(out[0]) == 0.0, f"got {float(out[0])}")

    # THE EARTH SAMPLE AND ITS SIZE FLOOR. The floor is a declared parameter
    # now, so the table it indexes has to hold the properties a count over a
    # size class has, and the sweep that reports it has to be able to say
    # "no lever" when there is none.
    check("the count in force is the table entry at the floor in force",
          EARTH_STANDING_BASINS == EARTH_STANDING_BY_FLOOR[EARTH_SIZE_FLOOR_KM2],
          f"{EARTH_STANDING_BASINS} against "
          f"{EARTH_STANDING_BY_FLOOR.get(EARTH_SIZE_FLOOR_KM2)}")
    floors = sorted(EARTH_STANDING_BY_FLOOR)
    counts = [EARTH_STANDING_BY_FLOOR[f] for f in floors]
    check("the Earth count never rises with the floor",
          all(a >= b for a, b in zip(counts, counts[1:])),
          f"counts {counts} are not non-increasing")

    # A synthetic basin population, in the units the real one is in: overflow in
    # km3 per year, depth in metres, erodibility and normalised gradient about 1.
    n_syn = 3000
    q_syn = np.abs(rng.lognormal(0.0, 2.0, n_syn))
    q_syn[: n_syn // 3] = 0.0                       # basins that do not overflow
    depth_syn = rng.lognormal(4.8, 1.2, n_syn)
    ero_syn = rng.lognormal(0.0, 0.3, n_syn)
    slope_syn = rng.lognormal(0.0, 0.8, n_syn)
    lat_syn = rng.uniform(-70.0, 70.0, n_syn)
    year_syn = 3.0e7

    # THE SOLVER LANDS ON THE CROSSING. `standing` is a non-increasing STEP
    # function of the coefficient, so no coefficient reproduces an arbitrary
    # count exactly and asking for one would fail on the step rather than on a
    # defect. The identity that does hold is that the solved value sits at the
    # step: nudge it up and the count is at or below the target, nudge it down
    # and it is at or above. Without this the sweep reports differences between
    # rungs that are bisection artefacts, and the count reported beside each
    # solved coefficient would not be the count that coefficient produces.
    standing_syn, solve_syn = _standing_solver(
        q_syn, year_syn, depth_syn, ero_syn, slope_syn, lat_syn)
    straddle = []
    for c in (50.0, 200.0, 800.0):
        target = standing_syn(c)
        solved = solve_syn(target)
        straddle.append((c, target, standing_syn(solved * (1.0 + 1e-6)),
                         standing_syn(solved * (1.0 - 1e-6))))
    check("the solved coefficient sits on the step it was solved for",
          all(above <= t <= below for _, t, above, below in straddle),
          "; ".join(f"C={c:g} target {t}: above {a}, below {b}"
                    for c, t, a, b in straddle))

    # THE SWEEP CAN SAY "NO LEVER". Fed a table whose count does not move with
    # the floor, it must report a lever of exactly 1 and call it SUBORDINATE. A
    # criterion that cannot return its negative class is not a criterion.
    flat = dict.fromkeys(EARTH_STANDING_BY_FLOOR, EARTH_STANDING_BASINS)
    cells_syn = np.full(400, 50.0)
    spill_syn = np.full(n_syn, 4000.0)
    flat_verdict = sweep_size_floor(
        q_syn, year_syn, depth_syn, ero_syn, slope_syn, lat_syn,
        EARTH_BAND_LAND_MKM2, cells_syn, spill_syn, table=flat)["verdict"]
    check("a floor that moves no count is called SUBORDINATE",
          bool(flat_verdict)
          and abs(flat_verdict["coefficient_lever"] - 1.0) < 1e-9
          and flat_verdict["coefficient_call"] == "SUBORDINATE",
          f"got {flat_verdict}")

    # THE CONCAVITY IN `Q**m`, checked against identities that hold whether or
    # not the delivery phase is known. `seasonal_concavity` says what the bound
    # is for; these are the three things it asserts, each with a right answer
    # fixed in advance and each able to fail.
    w_syn = np.full(12, 1.0 / 12.0)
    w_syn[3] = 0.5 / 12.0                 # a deliberately uneven weighting
    w_syn = w_syn / w_syn.sum()

    def _phi(bins, weights):
        bins = np.asarray(bins, dtype=float)
        return float((weights * np.power(bins, INCISION_EXPONENT)).sum()
                     / float((weights * bins).sum()) ** INCISION_EXPONENT)

    # THE DIRECTION, and it does not need the phase. `m` is under one, so the
    # power is concave and the annual evaluation is at or above the mean of the
    # per-bin evaluations, with equality only when the bins are identical.
    flat_bins = np.full(w_syn.size, 7.0)
    uneven = rng.lognormal(0.0, 1.5, (400, w_syn.size))
    phis = np.array([_phi(row, w_syn) for row in uneven])
    check("identical bins evaluate to the annual value exactly",
          abs(_phi(flat_bins, w_syn) - 1.0) < 1e-12,
          f"phi = {_phi(flat_bins, w_syn):.15g}")
    check("uneven bins evaluate strictly below the annual value",
          bool((phis < 1.0 - 1e-9).all()),
          f"max phi {float(phis.max()):.12g}")

    # THE BOUND. Nothing about the phase enters it: an overflow is never
    # negative, so the worst admissible cycle puts the whole year in the
    # lightest bin, and that is exactly sqrt(w_min).
    worst = np.zeros(w_syn.size)
    worst[int(np.argmin(w_syn))] = 1.0 / float(w_syn.min())
    floor_syn = float(w_syn.min()) ** (1.0 - INCISION_EXPONENT)
    check("the bound is attained by the whole year in the lightest bin",
          abs(_phi(worst, w_syn) - floor_syn) < 1e-12,
          f"{_phi(worst, w_syn):.12g} against {floor_syn:.12g}")
    check("no admissible cycle falls below the bound",
          bool((phis >= floor_syn - 1e-12).all()),
          f"min phi {float(phis.min()):.12g}")

    # THE LEVEL IS ABSORBED BY THE CALIBRATION, which is what makes a bracket
    # over the common discount a bracket of zero width. Scaling every discharge
    # so that Q**m is multiplied by phi, and dividing the coefficient by the
    # same phi, must return the retain vector unchanged.
    base_retain = incision_retain(q_syn, year_syn, depth_syn, ero_syn, 200.0,
                                  slope=slope_syn)
    phi_common = 0.37
    scaled = q_syn * phi_common ** (1.0 / INCISION_EXPONENT)
    same = incision_retain(scaled, year_syn, depth_syn, ero_syn,
                           200.0 / phi_common, slope=slope_syn)
    moved = incision_retain(scaled, year_syn, depth_syn, ero_syn, 200.0,
                            slope=slope_syn)
    check("a common concavity discount is absorbed by the coefficient exactly",
          float(np.abs(same - base_retain).max()) < 1e-12,
          f"max retain difference {float(np.abs(same - base_retain).max()):.3g}")
    check("and is not absorbed without it, so the check is not vacuous",
          int(((moved <= 0.0) != (base_retain <= 0.0)).sum()) > 0,
          "the uncompensated discount changed no basin's class")

    # THE MOVABLE MASK IS THE WHOLE OF THE EXPOSURE. Draw an admissible phi per
    # basin, re-solve the coefficient against the same standing count, and no
    # basin outside the band may change class. The band is the claim; this is
    # the thing that can falsify it.
    concav, movable_syn = seasonal_concavity(
        q_syn, year_syn, depth_syn, ero_syn, slope_syn, 200.0, w_syn,
        (186.0, 415.0))
    standing_base, _ = _standing_solver(q_syn, year_syn, depth_syn, ero_syn,
                                        slope_syn, lat_syn)
    target_syn = standing_base(200.0)
    # The draws go to the ENDS as well as through the middle. A band that is
    # merely too narrow is invisible to interior draws, because the escape needs
    # one basin at the floor while the population that sets the coefficient sits
    # at the top; the Bernoulli rows put basins there.
    escapes = 0
    draws = ([rng.uniform(floor_syn, 1.0, n_syn) for _ in range(4)]
             + [np.where(rng.random(n_syn) < f, floor_syn, 1.0)
                for f in (0.02, 0.2, 0.5, 0.8, 0.98)])
    for phi_syn in draws:
        perturbed = q_syn * phi_syn ** (1.0 / INCISION_EXPONENT)
        _, solve_p = _standing_solver(perturbed, year_syn, depth_syn, ero_syn,
                                      slope_syn, lat_syn)
        c_p = solve_p(target_syn)
        r_p = incision_retain(perturbed, year_syn, depth_syn, ero_syn, c_p,
                              slope=slope_syn)
        changed = (r_p <= 0.0) != (base_retain <= 0.0)
        escapes += int((changed & ~movable_syn).sum())
    check("no basin outside the movable band changes class under any draw",
          escapes == 0, f"{escapes} basins changed class outside the band")
    # The band is read off the cut/depth ratio and the verdict off the retain,
    # so `retain == 1 - ratio` has to hold exactly wherever the clip is not
    # active. Otherwise the band is drawn on a different quantity from the one
    # the verdict is taken on.
    unclipped = (base_retain > 0.0) & (base_retain < 1.0)
    x_syn = cut_over_depth(q_syn, year_syn, depth_syn, ero_syn, 200.0,
                           slope=slope_syn)
    agree = float(np.abs((1.0 - x_syn) - base_retain)[unclipped].max())
    check("the band's cut/depth is the retain's, wherever the retain is unclipped",
          bool(unclipped.any()) and agree == 0.0,
          f"max difference {agree:.3g} over {int(unclipped.sum())} basins")

    check("the band is not the whole population, so it says something",
          0 < concav["movable"] < concav["overflowing"],
          f"movable {concav['movable']} of {concav['overflowing']} overflowing")

    # Weights that are not the climatology's own would silently move the bound,
    # since it is read off the lightest bin.
    try:
        seasonal_concavity(q_syn, year_syn, depth_syn, ero_syn, slope_syn,
                           200.0, np.full(12, 1.0), (186.0, 415.0))
        refused = False
    except SystemExit:
        refused = True
    check("bin weights that do not sum to one are refused", refused,
          "the bound was computed off weights that are not a partition of the cycle")

    print(f"\n{n_checks} checks, {len(problems)} failed")
    return 1 if problems else 0


def main() -> None:
    ap = argparse.ArgumentParser()
    # Every per-build path below defaults to None and is resolved AFTER
    # parse_args. Resolving them while building the parser meant --help did
    # real work and died on a missing build, which is also why nothing caught
    # that the climatology default pointed at a superseded directory.
    #
    # The paths are per-build and strict. This script's output leaves the
    # project and changes the terrain, so a default that quietly reads or writes
    # the wrong build is the most expensive one here.
    ap.add_argument("--climatology", type=Path, default=None,
                    help="regular climatology; defaults to config's "
                         "baseline_climatology")
    # `docs/src/pipeline/loops.md`: the verdict map is ANTITONE -- carving removes the
    # bright closed-basin fill, the world warms, lake evaporation rises, and
    # basins that were marginal stay closed, so a larger carve set produces a
    # SMALLER next verdict. An antitone map oscillates rather than approaching a
    # fixed point from one side, so successive verdicts bracket the answer and
    # the honest procedure is to bracket it deliberately: take the verdict at
    # both bounding climates and carve only the intersection.
    #
    # Both arms are at the SAME flux, the design flux, differing only in
    # model.land_albedo_source. Section 4 says why: placing each endmember at
    # its own habitable flux puts both worlds at the same global mean by
    # construction, which collapses the bracket to the residual difference in
    # albedo PATTERN and reports a width narrower than the real uncertainty.
    ap.add_argument("--endmember-climatology", type=Path, default=None,
                    help="the cold bare-rock arm's regular climatology. Given, "
                         "the carve set becomes the INTERSECTION: a basin is "
                         "cut only where both climates cut it, and the "
                         "disagreement is reported as the BRACKETED set and the "
                         "bracket width. Absent, this is a single-climate "
                         "verdict and the sidecar says so")
    ap.add_argument("--coupling", type=Path, default=None)
    # Hydrography is per-build now, so the basin set has to be selectable
    # alongside the coupling it was built with. Mixing a coupling matrix from one
    # terrain with a basin catalogue from another would misalign the rows in the
    # same silent way the longitude convention misaligned the columns.
    ap.add_argument("--basins", type=Path, default=None)
    # Which field stands for runoff generated over the catchment. This was
    # `mrro` and should not have been: `mrro` is not local runoff generation but
    # river-routed net divergence, because landmod.f90's roffstep calls mkradv,
    # which advects runoff downhill and modifies its argument in place. So it is
    # local generation minus river outflow plus river inflow, on ExoPlaSim's own
    # grid and its own downhill directions, which know nothing about our basins.
    #
    # P - E is the water balance the derivation actually wants: whatever falls on
    # the catchment and does not evaporate is what reaches the sink. Both of this
    # project's other consumers of catchment runoff already made this call --
    # pedology's `runoff_source: p_minus_e`, with the reasoning in
    # pedogenesis.yaml, and surface_water.py's lake solver. The carve verdict was
    # the only one left on `mrro`, and it is the one whose output changes the
    # terrain.
    ap.add_argument("--runoff-source", choices=("p_minus_e", "mrro"),
                    default="p_minus_e")
    # Iteration 2 onward. A verdict is taken on the basins a build still has,
    # but Orogen regenerates from the planet code and needs the whole catalogue,
    # including the basins an earlier pass already carved. Merging keeps the loop
    # monotone: a basin carved in a previous iteration stays carved, because the
    # terrain that justified re-examining it no longer exists. Without this the
    # list would silently re-preserve every basin the current build has already
    # lost, and the next build would undo the last one.
    ap.add_argument("--previous", type=Path, default=None,
                    help="carve_list.json from the pass that produced the "
                         "current build; its retain-0 basins are carried forward")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--out-list", type=Path, default=None)
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--for-build", default=None,
                    help="declare a DELIBERATE cross-build read of the staged "
                         "background albedo by naming the build it was staged "
                         "from. A name from lib/orogen.py's registry; any other "
                         "build is refused.")
    ap.add_argument("--selftest", action="store_true",
                    help="run the basis conversion's identity checks; needs "
                         "no build and no climatology")
    ap.add_argument("--measure-earth-floors", action="store_true",
                    help="re-measure EARTH_STANDING_BY_FLOOR from HydroLAKES "
                         "and HydroBASINS under hydrography/data/reference/ "
                         "and print the table. Needs the BULK sources and "
                         "neither a build nor a climatology; the verdict itself "
                         "reads only the table")
    args = ap.parse_args()

    if args.selftest:
        raise SystemExit(_selftest())

    if args.measure_earth_floors:
        measured = measure_earth_floors(DATA / "reference")
        print(measured["selection"])
        print(f"unmatched pour points: {measured['unmatched_pour_points']}")
        print("EARTH_STANDING_BY_FLOOR = {")
        for floor, count in measured["counts"].items():
            drift = EARTH_STANDING_BY_FLOOR.get(floor)
            note = "" if drift is None or drift == count else f"   # was {drift}"
            print(f"    {floor}: {count},{note}")
        print("}")
        band = measured["band_land"]
        print(f"EARTH_BAND_LAND_MKM2 = {band['band_land_mkm2']:.4g}   "
              f"({band['polygons_in_band']:,} of {band['polygons']:,} level 5 "
              f"polygons, {band['global_land_mkm2']:.4g} Mkm2 globally)")
        print(f"  {band['assignment']}")
        raise SystemExit(0)

    _bd = component_data("hydrography", strict=True)
    # LOADED BEFORE THE DEFAULTS: `coupling_path` reads `model.resolution`,
    # because the coupling matrix is per rung.
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if args.coupling is None:
        args.coupling = cv.coupling_path(_bd, config)
    if args.basins is None:
        args.basins = _bd / "basins.nc"
    if args.out_list is None:
        args.out_list = _bd / "carve_list.txt"
    if args.out_json is None:
        args.out_json = _bd / "carve_list.json"
    if args.climatology is None:
        args.climatology = climatology_path()
    # THE RUNG GUARD RUNS WHATEVER THE CLIMATOLOGY CAME FROM. It used to be
    # reached only through `climatology_path()`, which is called only when
    # `--climatology` is ABSENT -- so every re-take, every arm of a bracket and
    # every sensitivity, which all pass the flag, skipped it. The rung appears
    # nowhere in a climatology's name, so nothing else would have caught a T85
    # file driving a T21 configuration, and the sidecar would have recorded the
    # configured rung beside it with nothing contradicting itself.
    require_configured_grid(args.climatology, config)
    if args.endmember_climatology is not None:
        # The cold arm is a climatology like any other, and an arm on the wrong
        # grid is a bracket whose ends are not comparable.
        require_configured_grid(args.endmember_climatology, config)

    basins = BasinSet(args.basins)
    resolution = str(config["model"]["resolution"]).upper()

    # The list is written in the basis Orogen spends it against, and that
    # conversion needs the catalogue's own natural-to-final depth ratio. It is
    # not optional: without it the marginal class crosses the interface as a
    # fraction of one depression spent against another, which on this terrain
    # is a median 2.1x the intended incision.
    if basins.retained_fraction is None:
        raise SystemExit(
            f"{args.basins} carries no `retained_fraction`, so a retain "
            "computed against the finished depression cannot be put into the "
            "natural-relief basis the carve list declares. Rebuild it with "
            "hydrography/scripts/build_hydrography.py, which writes the field.")
    retained_fraction = basins.retained_fraction

    # The primary arm: the warm, vegetated end of section 4's bracket.
    primary = climate_terms(args.climatology, args, config, basins)
    (year_s, crit, runoff, precip, idx_pen, idx_wet, q_pen, carved, overflows_wet, disputed, e_pen, e_wet, margin, retain_margin, ocean_validation, penman_error_pct) = (
        primary["year_s"], primary["crit"], primary["runoff"], primary["precip"], primary["idx_pen"], primary["idx_wet"], primary["q_pen"], primary["carved"], primary["overflows_wet"], primary["disputed"], primary["e_pen"], primary["e_wet"], primary["margin"], primary["retain_margin"], primary["ocean_validation"], primary["penman_error_pct"])

    # What the water can actually cut, which is the half the margin never knew.
    # The sill's own rock and the outlet's own gradient set how much water that
    # takes, and the coefficient is solved against Earth rather than declared.
    sill_ero = sill_erodibility(args.basins, basins.terrain_hash)
    export = Export()
    if export.terrain_hash != basins.terrain_hash:
        raise SystemExit(
            f"basins.nc was built from terrain {basins.terrain_hash[:16]} and "
            f"the configured export is {export.terrain_hash[:16]}; the outlet "
            "gradient is looked up by region index, which does not survive a "
            "terrain change")
    regions_path = args.basins.parent / "regions.nc"
    slope, slope_reference, slope_measured = outlet_gradient(
        args.basins, regions_path, export)
    with Dataset(args.basins) as ds:
        sink_region = np.asarray(ds["sink"][:]).astype(np.int64)
    basin_latitude = np.asarray(export.field("lat"))[sink_region]
    band_land_mkm2 = float(
        export.cell_area[(export.surface_class == LAND)
                         & (np.abs(export.lat) < CALIBRATION_LATITUDE)].sum()) / 1e6
    coefficient, coefficient_bracket, standing_target, standing_here = (
        calibrate_coefficient(q_pen, year_s, basins.depth_at_spill_m, sill_ero,
                              slope, basin_latitude, band_land_mkm2))
    # What the size floor on the Earth sample is worth, reported BESIDE the
    # verdict rather than left in a comment. The floor is the largest single
    # lever on the coefficient and it is not the one the Poisson bracket above
    # measures; both go into the sidecar so a reader can tell them apart.
    floor_sweep = sweep_size_floor(
        q_pen, year_s, basins.depth_at_spill_m, sill_ero, slope, basin_latitude,
        band_land_mkm2,
        export.cell_area[export.surface_class == LAND],
        basins.area_at_spill_km2)
    retain_incision_finished = incision_retain(
        q_pen, year_s, basins.depth_at_spill_m, sill_ero, coefficient,
        slope=slope)
    # Into the basis the list declares. `basis_usable` is false where the
    # conditioning left no finished depression to have measured a ratio on.
    retain_incision, basis_usable = to_natural_relief_basis(
        retain_incision_finished, retained_fraction)
    print(f"outlet gradient   {slope_measured:5d} basins measured, "
          f"population {slope_reference:.5f} m/m")
    print(f"coefficient       {coefficient:7.1f} m per (m3/s)^0.5, bracket "
          f"{coefficient_bracket[0]:.1f}-{coefficient_bracket[1]:.1f}; "
          f"{standing_here} standing within {CALIBRATION_LATITUDE:g} deg "
          f"against Earth's {standing_target:.1f}")
    _span = floor_sweep["span_km2"]
    print(f"size floor        {EARTH_SIZE_FLOOR_KM2:g} km2 on the Earth sample, "
          f"swept over {_span[0]:.0f}-{_span[1]:.0f} km2 "
          f"(mesh cell to median area at spill)")
    if floor_sweep["verdict"]:
        _v = floor_sweep["verdict"]
        print(f"                  coefficient {_v['coefficient_over_span'][0]:.1f}"
              f"-{_v['coefficient_over_span'][1]:.1f} across it, a factor of "
              f"{_v['coefficient_lever']:.2f} against {_v['coefficient_poisson_ratio']:.2f} "
              f"for the Poisson error on Earth's {EARTH_STANDING_BASINS}: "
              f"{_v['coefficient_call']}")
        print(f"                  marginal {_v['marginal_over_span'][0]}"
              f"-{_v['marginal_over_span'][1]} across it, a factor of "
              f"{_v['marginal_lever']:.2f} against {_v['marginal_poisson_ratio']:.2f}: "
              f"{_v['marginal_call']}")

    # Either is a reason to leave a rim standing: that the basin may not overflow
    # at all, or that its overflow cannot cut. Taking the larger keeps both.
    retain = np.maximum(retain_incision, retain_margin)
    # The same composite on the basis the incision was COMPUTED on, kept only
    # so the sidecar can say what the conversion moved. Never written to the
    # list: Orogen cannot spend it.
    retain_finished = np.maximum(retain_incision_finished, retain_margin)

    # THE INTERVAL BRACKET, folded in by the same idiom as the arms below and
    # for the same reason. `cv._INTERVAL_BRACKET` carries the argument: the bin
    # mean is the limit for a lake with no heat storage and the annual
    # evaluation the limit for one deep enough to hold its temperature through
    # the year, and the spread between them is the water body's depth, which
    # this project has no term for. So the rim stands wherever EITHER end would
    # leave it standing.
    #
    # This is not the arm bracket. The arms vary the CLIMATE and this varies the
    # interval one climate is read over; a basin can be bracketed by one and not
    # the other, and merging them would put two uncertainties under one name.
    #
    # Same coefficient, for the reason the arm block gives: it is calibrated
    # against Earth from the discharge field, so letting each end calibrate its
    # own would move the yardstick with the estimate.
    interval_annual = climate_terms(args.climatology, args, config, basins,
                                    interval="annual")
    retain_incision_annual = incision_retain(
        interval_annual["q_pen"], year_s, basins.depth_at_spill_m, sill_ero,
        coefficient, slope=slope)
    retain_interval_annual = np.maximum(
        to_natural_relief_basis(retain_incision_annual, retained_fraction)[0],
        interval_annual["retain_margin"])
    cut_bin_mean = retain <= 0.0
    cut_annual = retain_interval_annual <= 0.0
    retain = np.maximum(retain, retain_interval_annual)
    retain_finished = np.maximum(
        retain_finished,
        np.maximum(retain_incision_annual, interval_annual["retain_margin"]))
    n_interval_both = int((cut_bin_mean & cut_annual).sum())
    n_interval_either = int((cut_bin_mean | cut_annual).sum())
    n_interval_split = n_interval_either - n_interval_both
    interval_width = 100.0 * n_interval_split / max(n_interval_either, 1)
    print(f"interval bracket  bin mean cuts {int(cut_bin_mean.sum())}, "
          f"annual evaluation cuts {int(cut_annual.sum())}")
    print(f"                  BOTH {n_interval_both}, either {n_interval_either}, "
          f"disagreement {n_interval_split} ({interval_width:.1f}% of the union)")
    print("                  the spread is the water body's heat storage, which "
          "has no term here; the rim stands where either end leaves it standing")

    # THE INTERSECTION, and it is the same idiom one line up: taking the larger
    # retain keeps the rim wherever EITHER climate would have kept it, which is
    # exactly "carve only what both carve" generalised to a fractional rim. A
    # basin cut to bare rock under the warm arm and left standing under the cold
    # one comes out standing.
    #
    # The coefficient is the primary arm's, deliberately. It is calibrated
    # against Earth from the discharge field, so letting the cold arm calibrate
    # its own would move the yardstick with the climate and the two verdicts
    # would no longer be measuring against the same thing.
    endmember = None
    retain_primary = retain
    # Per-basin bracket membership, not just its width. The width says how much
    # the two bounding climates disagree; it does not say WHICH basins, and a
    # later pass asking whether a basin that flipped had been bracketed or had
    # been in the intersection cannot answer that from a count. Recorded on
    # every basin below when both arms ran, absent when one did.
    retain_endmember = None
    if args.endmember_climatology is not None:
        endmember = climate_terms(args.endmember_climatology, args, config,
                                  basins)
        # THE ARMS MUST ACTUALLY DIFFER. `climate_terms` took a climatology
        # path and read `args.climatology` regardless, so both arms were the
        # primary and the bracket reported zero width -- which reads as the two
        # bounding climates agreeing rather than as never having been asked.
        # Two different files that produce a bit-identical discharge field are
        # not a bracket, whatever the cause.
        if np.array_equal(endmember["q_pen"], primary["q_pen"]):
            raise SystemExit(
                f"the endmember arm ({args.endmember_climatology}) produced a "
                f"discharge field bit-identical to the primary arm "
                f"({args.climatology}). Two bounding climates that agree to the "
                "last bit have not been evaluated on two climates, and a "
                "bracket of zero width would be reported as agreement")
        retain_endmember_incision_finished = incision_retain(
            endmember["q_pen"], year_s, basins.depth_at_spill_m, sill_ero,
            coefficient, slope=slope)
        retain_endmember = np.maximum(
            to_natural_relief_basis(retain_endmember_incision_finished,
                                    retained_fraction)[0],
            endmember["retain_margin"])
        retain_finished = np.maximum(
            retain_finished,
            np.maximum(retain_endmember_incision_finished,
                       endmember["retain_margin"]))
        retain = np.maximum(retain_primary, retain_endmember)
        cut_primary = retain_primary <= 0.0
        cut_endmember = retain_endmember <= 0.0
        n_both = int((cut_primary & cut_endmember).sum())
        n_either = int((cut_primary | cut_endmember).sum())
        n_split = n_either - n_both
        # The width IS the uncertainty, per section 4, and is reported as a
        # share of the union rather than of the basin count: a bracket that
        # disagrees on 200 of 220 candidates is not the same result as one that
        # disagrees on 200 of 3,600, and the denominator that carries that is
        # how many basins either arm would cut.
        width = 100.0 * n_split / max(n_either, 1)
        print(f"intersection      warm arm cuts {int(cut_primary.sum())}, "
              f"cold arm cuts {int(cut_endmember.sum())}")
        print(f"                  BOTH {n_both}, either {n_either}, "
              f"disagreement {n_split} ({width:.1f}% of the union)")
        print(f"                  the {n_split} are the BRACKETED set BY "
              f"CONSTRUCTION, not by a tolerance chosen afterwards. "
              f"`bracketed` is our uncertainty; `marginal` below is a "
              f"landform and the two are unrelated")

    # WHAT THE ANNUAL DISCHARGE UNDER THE SQUARE ROOT IS WORTH, bounded from the
    # bin weights and no delivery phase. `seasonal_concavity` carries the
    # argument; the two results it turns on are that the calibration absorbs any
    # common concavity discount exactly, so this is not a second carve list, and
    # that the residual spread cannot move a basin more than `J` either side of
    # the class boundary. Run on every arm, because the final class is a max over
    # the arms and a basin is exposed if ANY of them can be moved.
    with Dataset(args.climatology) as ds:
        concavity_weights = climatology.bin_weights(np.asarray(ds["time"][:]))
    concavity, concavity_movable = seasonal_concavity(
        q_pen, year_s, basins.depth_at_spill_m, sill_ero, slope, coefficient,
        concavity_weights, coefficient_bracket)
    concavity_arms = {"warm_vegetated_bin_mean": concavity}
    for _label, _arm in (("warm_vegetated_annual", interval_annual),
                         ("cold_bare_rock_bin_mean", endmember)):
        if _arm is None:
            continue
        _report, _mask = seasonal_concavity(
            _arm["q_pen"], year_s, basins.depth_at_spill_m, sill_ero, slope,
            coefficient, concavity_weights, coefficient_bracket)
        concavity_arms[_label] = _report
        concavity_movable = concavity_movable | _mask
    concavity["arms"] = {k: {"movable": v["movable"],
                             "movable_of_the_cut": v["movable_of_the_cut"],
                             "movable_of_the_standing": v["movable_of_the_standing"],
                             "overflowing": v["overflowing"]}
                         for k, v in concavity_arms.items()}
    concavity["movable_over_all_arms"] = int(concavity_movable.sum())
    concavity["movable_note"] = (
        "the counts above the `arms` block are the primary arm's, the arm the "
        "coefficient is solved on. `seasonal_concavity_movable` on each basin "
        "is the UNION over the arms, because the class a basin ends in is a max "
        "over them and it is exposed if any one of them can be moved.")
    print(f"concavity         Q**{INCISION_EXPONENT:g} is evaluated on an ANNUAL "
          f"discharge; the delivery phase is an absence, the level is absorbed "
          f"by the calibration")
    print(f"                  phi in {concavity['phi_range']}, so a cut moves at "
          f"most a factor {concavity['factor_J']:.3f} either side of the "
          f"boundary, against P = {concavity['coefficient_poisson_ratio_P']:.3f}: "
          f"{concavity['call']}")
    print(f"                  {concavity['movable_over_all_arms']} basins can "
          f"change class under some admissible delivery phase, of "
          f"{concavity['overflowing']} that overflow; the rest hold their "
          f"verdict under every one")

    # The superseded mapping, kept for comparison in the sidecar only.
    span = idx_pen - idx_wet
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(np.isfinite(span) & (span > 1e-9),
                     (crit - idx_wet) / np.where(span > 1e-9, span, 1.0), 0.5)
    retain_span = np.where(carved, 0.0,
                           np.where(~overflows_wet, 1.0, np.clip(1.0 - f, 0.0, 1.0)))

    # The verdict follows retain rather than the overflow test, because what
    # Orogen does to a basin is set by retain. A basin can overflow and still
    # keep most of its rim, if what crosses the sill is a trickle.
    verdict_name = np.where(retain <= 0.0, "carve",
                            np.where(retain >= 1.0, "preserve", "marginal"))
    n_carve = int((retain <= 0.0).sum())
    n_preserve = int((retain >= 1.0).sum())
    n_marginal = basins.n - n_carve - n_preserve
    # Basins the drainage conditioning already flattened. They are published as
    # preserved with retain 1 and impound nothing, so the balance above has no
    # lake to run on and depth_at_spill_m falls back to its 1 m guard, which
    # carves them at any coefficient. That is the right instruction, but it is
    # arrived at by a floor rather than by a verdict, so it is counted and
    # reported rather than left to look like a decision.
    no_impoundment = ~basins.has_impoundment
    n_no_impoundment = int(no_impoundment.sum())
    n_no_impoundment_carved = int((no_impoundment & (retain <= 0.0)).sum())

    # WHAT THE BASIS CONVERSION MOVED, and what is left over after it. The
    # criteria here are fixed by construction rather than chosen: the marginal
    # class is the one the conversion touches, the shift is reported over that
    # class, and the residual is a bracket because the ratio used is this
    # generation's estimate of the next one's.
    marginal_here = (retain > 0.0) & (retain < 1.0)
    marginal_finished = (retain_finished > 0.0) & (retain_finished < 1.0)
    n_rebased = int((marginal_here | marginal_finished).sum())
    shift = retain[marginal_here] - retain_finished[marginal_here]
    verdict_finished = np.where(retain_finished <= 0.0, "carve",
                                np.where(retain_finished >= 1.0, "preserve",
                                         "marginal"))
    n_verdict_moved = int((verdict_name != verdict_finished).sum())
    ratio_usable = retained_fraction[basis_usable]
    ratio_q = (np.percentile(ratio_usable, [5, 25, 50, 75, 95])
               if ratio_usable.size else np.full(5, np.nan))
    basis_report = {
        "basis": "fraction of the depression's NATURAL spill depth, in the "
                 "generator's dimensionless elevation parameter -- the "
                 "quantity buildBasinProtection spends (1 - retain) against. "
                 "vendor/orogen/tools/README.md declares it on the format.",
        "why_not_the_finished_depth": "protection is built before erosion "
                 "runs, so the finished depression does not exist when the "
                 "allowance is set, and the carve operates on model "
                 "elevations. Orogen cannot change basis; the conversion "
                 "belongs here.",
        "conversion": "retain = 1 - (1 - retain_finished_depth_basis) * "
                 "retained_fraction, with both endpoints passed through: "
                 "retain 1 is an allowance of zero and retain 0 clears "
                 "noLower, and both mean the same on either basis",
        "assumptions": [
            "the fraction of a depression's relief is read as basis-free "
            "within that one depression, which is the height curve's secant "
            "across it rather than its local slope. The curve is quartic on "
            "land, so this is a linearisation over the depression.",
            "retained_fraction is measured on THIS generation's conditioning "
            "and the allowance is spent on the NEXT one's. It is a property "
            "of the conditioning rather than of the carve list, which is what "
            "makes it transferable, and the residual below is its bracket.",
        ],
        "retained_fraction_quantiles_5_25_50_75_95": [round(float(x), 4) for x in ratio_q],
        "residual_ratio_bracket": (
            [round(float(ratio_q[0] / ratio_q[2]), 3),
             round(float(ratio_q[4] / ratio_q[2]), 3)]
            if ratio_usable.size and ratio_q[2] > 0 else None),
        "residual_note": "the factor between the incision Orogen spends and "
                 "the incision intended, if the next generation's "
                 "conditioning gives a basin the population's 5th or 95th "
                 "percentile ratio in place of the per-basin value used here. "
                 "A BRACKET, not a point: the per-basin conversion is exact "
                 "against this build and the next build's ratio is unknown.",
        "basins_without_ratio": int((~basis_usable).sum()),
        "basins_without_ratio_note": "no finished depression to have measured "
                 "a ratio on. They keep their finished-basis value, which for "
                 "these is the 1 m depth guard driving retain to 0, and 0 "
                 "means the same on either basis.",
        "marginal_basins_rebased": n_rebased,
        "marginal_shift_5_50_95": (
            [round(float(x), 4) for x in np.percentile(shift, [5, 50, 95])]
            if shift.size else None),
        "verdict_class_moved": n_verdict_moved,
        "verdict_class_moved_note": "basins whose class differs between the "
                 "two bases. A marginal basin whose intended incision exceeds "
                 "the whole natural relief cannot be expressed as an "
                 "allowance and saturates to carve, which is the honest "
                 "reading of an overflow that removes the depression.",
    }

    with Dataset(args.basins) as ds:
        ids = [str(x) for x in ds["basin_id"][:]]

    carried = {}
    if args.previous is not None:
        prev = json.loads(args.previous.read_text(encoding="utf-8"))
        here = set(ids)
        for entry in prev["basins"]:
            if entry["id"] not in here:
                # Absent from this build because it was carved away. Its verdict
                # is not re-decidable and is carried at 0.
                if float(entry["retain"]) > 0.0:
                    raise RuntimeError(
                        f"basin {entry['id']} is absent from this build but "
                        f"carries retain {entry['retain']}; carried-forward "
                        "entries must all be retain 0")
                carried[entry["id"]] = 0.0

    pass_number = 1 if args.previous is None else prev.get("pass_number", 1) + 1
    pass_label = ("first pass" if args.previous is None
                  else f"pass {pass_number}, {len(carried)} basins carried forward")
    n_carried = len(carried)
    n_zero = n_carve + n_carried
    n_total = basins.n + n_carried
    n_here = basins.n
    clim_name = args.climatology.name
    flux_earth = float(config["orbit"]["baseline_flux_earth"])
    with Dataset(args.climatology) as ds:
        _ts = climatology.annual_mean_of(ds, "ts")
        _lat = np.asarray(ds["lat"][:])
        _w = np.cos(np.deg2rad(_lat))[:, None] * np.ones_like(_ts)
        mean_ts = float((_w * _ts).sum() / _w.sum())
    runoff_source = ("precipitation minus evaporation over the catchment"
                     if args.runoff_source == "p_minus_e"
                     else "the model's mrro field")
    land_surface = str(config["model"].get("land_albedo_source", "uniform"))
    glaciers = ("glaciers enabled" if (config["surface"].get("glaciers") or {}).get("enabled")
                else "glaciers off")
    n_without_ratio = int((~basis_usable).sum())
    header = f"""# Vesper carve verdict, {pass_label}
#
# Produced from a converged ExoPlaSim climatology: {resolution},
# {flux_earth:g} S-Earth, {land_surface} land surface, {glaciers},
# {mean_ts:.2f} K, climatology {clim_name}.
# Terrain {basins.terrain_hash[:16]}, basin catalogue unchanged.
#
# A basin overflows when more water arrives than its lake surface can evaporate,
#     Q = runoff * (catchment - area_at_spill) - (E - P) * area_at_spill  >  0
# and retain is then what that Q can cut against how deep the basin is,
#     retain = 1 - {coefficient:.1f} * erodibility * slope^{SLOPE_EXPONENT:g} * Q^{INCISION_EXPONENT:g} / depth_at_spill
# in metres, where slope is the outlet channel's own gradient relative to the
# population's, and the coefficient is solved against the density of Earth's own
# standing through-flowing impounded basins rather than declared. Open-water
# evaporation is the Penman combination equation with water's albedo and
# roughness, validated against the model over ocean cells to within
# {penman_error_pct:.1f}%.
# Catchment runoff is {runoff_source}; see the sidecar.
#
#   retain 1.0   {n_preserve:4d} basins  closed: the lake surface evaporates all that arrives
#   retain 0<r<1 {n_marginal:4d} basins  a notch: the outlet is cut, but not to the floor
#   retain 0.0   {n_zero:4d} basins  carved
#                            {n_carried:4d} of them carried forward, {n_carve:4d} decided here
#
# {n_no_impoundment:4d} of the entries above hold no water at spill at all: the
# drainage conditioning took their rim down to their floor, and the catalogue
# publishes them as preserved anyway. {n_no_impoundment_carved:d} of them carve here, on the 1 m
# depth guard rather than on a water balance. notes/audits/basin-catalogue-floor.md.
#
# The counts above are of the whole {n_total:d}-entry catalogue. This pass could
# only decide the {n_here:d} basins the current build still has; the rest were
# carved by an earlier pass and are held at 0 to keep the loop monotone.
#
# THE BASIS. retain here is a fraction of the depression's NATURAL spill depth,
# in the generator's dimensionless elevation parameter -- the depth
# buildBasinProtection spends (1 - retain) against, and the only depth that
# exists when protection is built. The incision above is a LENGTH in metres
# against the FINISHED depression, so it is converted before it is written:
#     retain = 1 - (1 - retain_finished) * retainedFraction
# using this build's own natural-to-final ratio as the estimate of the next
# generation's. Both endpoints pass through exactly, and the sidecar carries
# the ratio's population spread as the bracket on that estimate.
# {n_without_ratio:4d} entries had no ratio to convert with and are all at retain 0.
#
# Carved basins are listed explicitly at retain 0, so this file is the
# complete verdict.
#
# id                                    retain
"""
    lines = [header]
    for i in np.argsort(-retain):
        lines.append(f"{ids[i]:<38s} {retain[i]:.4f}")
    for bid in sorted(carried):
        lines.append(f"{bid:<38s} {0.0:.4f}")
    args.out_list.parent.mkdir(parents=True, exist_ok=True)
    args.out_list.write_text("\n".join(lines) + "\n", encoding="ascii")

    sidecar = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "pass_number": pass_number,
        "terrain_hash": basins.terrain_hash,
        "climatology": str(args.climatology),
        "staged_background_albedo": primary["staged_background_albedo"],
        "climate": {
            "resolution": resolution,
            "flux_earth": flux_earth,
            "stellar_spectrum": str(config["radiation"].get("stellar_spectrum")),
            "land_surface": land_surface,
            "mean_surface_temperature_k": round(mean_ts, 2),
            "note": None if land_surface == "modelled" else
                    "The land surface is assumed rather than modelled, and the "
                    "albedo bracket puts that assumption at 3.7 to 7.1 K.",
        },
        "catalogue": {
            "basins": basins.n,
            "without_impoundment": n_no_impoundment,
            "without_impoundment_carved": n_no_impoundment_carved,
            "note": "A preserved basin is one the drainage conditioning was told "
                    "not to breach by its CARVE passes; erosion, ridge sharpening "
                    "and soil creep run over it regardless and the sink is left "
                    "unprotected. These entries came out of that with no "
                    "depression left, are published as preserved with retain 1, "
                    "and have no water balance to decide.",
        },
        "method": {
            "test": "Q = runoff*(catchment - area_at_spill) - (E - P)*area_at_spill > 0",
            "test_note": "equivalent to (E - P) / runoff <= catchment / "
                         "area_at_spill - 1 wherever catchment runoff is "
                         "positive, and defined where it is not: a basin with no "
                         "catchment still fills if its own lake surface gains "
                         "more precipitation than it evaporates",
            "open_water_evaporation": "Penman combination, water albedo and "
                                      "roughness, evaluated wholly at the "
                                      "lowest model level and NOT floored at "
                                      "the model's land rate: this world's land "
                                      "is several times rougher than open "
                                      "water, so a smooth lake in a rough wet "
                                      "landscape evaporates less than the "
                                      "ground around it",
            "runoff_source": args.runoff_source,
            "recorded_terms_note": "precipitation_km_per_year, "
                                   "lake_evaporation_km_per_year and "
                                   "land_evaporation_km_per_year are the "
                                   "catchment means the test above is built "
                                   "from. Under runoff_source p_minus_e, "
                                   "runoff_km_per_year is max(P - E_land, 0) by "
                                   "construction, which is the identity a "
                                   "consumer should check before trusting them.",
            "penman_ocean_validation_ratio": ocean_validation["ratio"],
            "retain_mapping": "max(incision, margin), the larger of what the "
                              "overflow cannot cut and what the overflow test "
                              "cannot decide, then put into the basis "
                              "carve_list_basis declares",
            "retain_incision_mapping":
                f"1 - {coefficient:.4g} * erodibility * slope**{SLOPE_EXPONENT:g} "
                f"* Q**{INCISION_EXPONENT:g} / depth_at_spill_m, clipped to [0, 1], "
                "and then rebased -- the incision is a LENGTH and depth_at_spill_m "
                "is the finished depression in metres, which is not the depth "
                "Orogen spends the allowance against",
            "carve_list_basis": basis_report,
            "retain_incision_coefficient_m_per_sqrt_q": round(coefficient, 4),
            "retain_incision_coefficient_bracket": [round(b, 4) for b in coefficient_bracket],
            "retain_incision_exponent": INCISION_EXPONENT,
            "retain_slope_exponent": SLOPE_EXPONENT,
            "retain_slope_exponent_bracket": list(SLOPE_EXPONENT_BRACKET),
            "outlet_gradient": {
                "reference_m_per_m": round(slope_reference, 6),
                "basins_measured": slope_measured,
                "range": [round(float(slope.min()), 4), round(float(slope.max()), 4)],
                "note": "gradient of the outflow channel below each saddle, ten "
                        "receiver steps from spill_exit_region, normalised by the "
                        "population's geometric mean. The saddle's OWN slope is "
                        "zero by construction and is not what stream power reads; "
                        "the channel below it spans a dozen mesh regions and is "
                        "measurable. It correlates with depth at spill at r = "
                        "0.735, which is why it can no longer be absorbed into "
                        "the coefficient: at n = 1 the depth cancels out of "
                        "retain. HYD-15",
            },
            "calibration": {
                "earth_standing_basins": EARTH_STANDING_BASINS,
                "earth_size_floor_km2": EARTH_SIZE_FLOOR_KM2,
                "earth_band_land_mkm2": EARTH_BAND_LAND_MKM2,
                "band_latitude_deg": CALIBRATION_LATITUDE,
                "vesper_band_land_mkm2": round(band_land_mkm2, 2),
                "target_standing_basins": round(standing_target, 2),
                "achieved_standing_basins": standing_here,
                "size_floor_sensitivity": floor_sweep,
                "size_floor_note": "the Earth count is a count at a SIZE CLASS, "
                    "and the floor that sets it is the largest single lever on "
                    "the coefficient -- larger than the Poisson error on the "
                    "count itself, which is what the bracket beside the "
                    "coefficient measures. size_floor_sensitivity re-solves at "
                    "every rung over a span the build derives, and calls the "
                    "lever against that noise on a criterion fixed beforehand.",
            },
            "seasonal_concavity": concavity,
            "retain_incision_rationale":
                "stream power goes as K Q^m S^n, so what the overflow achieves is "
                "a LENGTH of incision, and whether that empties the basin depends "
                "on the depth from spill point to floor. S is MEASURED per basin "
                "rather than absorbed into the coefficient, because the outlet's "
                "gradient and the basin's depth are the same relief and folding "
                "one into a constant while dividing by the other counts it twice. "
                "The coefficient absorbs K and the relaxation window, and it is "
                "SOLVED rather than declared because Orogen has no time axis: it "
                "is whatever matches the density of Earth's standing "
                "through-flowing impounded basins, at the size this mesh "
                "resolves. It is solved on every run because the count it matches "
                "depends on the climate, so a literal goes stale when the verdict "
                "does. See hydrography/notes/retain-fraction.md",
            "sill_rock": "K from the export's erodibility field, a relative "
                "stream-power multiplier mean-normalised to 1 over land, taken "
                "as the geometric mean of the two regions either side of the "
                "saddle. It is the fluvially expressed contrast rather than the "
                "intact-rock one, which is what a landscape model wants; see "
                "Zondervan (2020) in references/INDEX.md",
            "sill_erodibility_range": [round(float(sill_ero.min()), 4),
                                       round(float(sill_ero.max()), 4)],
            "retain_margin_mapping": "min(1, ((E_penman - (P + critical*runoff)) / E_penman) / 0.25)",
            "retain_tolerance": TOLERANCE,
            "retain_tolerance_rationale": "fractional change in open-water "
                "evaporation that would flip the verdict. 0.25 is set by the "
                "terms that dominate -- the assumed biosphere at 3.7 to 7.1 K, "
                "dust at 10 to 20% of lake evaporation, the sub-grid dry column "
                "at up to 9% -- and not by Penman's own method error, which is "
                "the smallest of them and is computed per run rather than "
                "quoted. It is reported and never decisive: retain is the "
                "larger of this and the incision value, and the margin is "
                "positive exactly where the incision value is already 1.",
        },
        "counts": {"carve": n_carve, "preserve": n_preserve, "marginal": n_marginal},
        # `docs/src/pipeline/loops.md`'s bracket. Recorded rather than printed, because
        # the WIDTH is the honest uncertainty on the carve and something
        # downstream will want to quote it -- error_budget.py above all, which
        # prices items in basins. `single_climate` says plainly when no bracket
        # was taken, so a reader cannot mistake one arm for two.
        "intersection": ({
            "single_climate": True,
            "note": "no --endmember-climatology given, so this is ONE arm of "
                    "section 4's bracket and carries no width. The verdict is "
                    "not robust to the vegetation question; see "
                    "docs/src/pipeline/loops.md.",
        } if endmember is None else {
            "single_climate": False,
            "endmember_climatology": rel(args.endmember_climatology),
            "cut_by_warm_vegetated_arm": int(cut_primary.sum()),
            "cut_by_cold_bare_rock_arm": int(cut_endmember.sum()),
            "cut_by_both": n_both,
            "cut_by_either": n_either,
            "bracketed": n_split,
            "width_pct_of_union": round(width, 3),
            "bracketed_is_not_marginal": "`bracketed` counts basins the two "
                    "bounding climates DISAGREE about, which is a statement "
                    "about our uncertainty. `marginal` in counts above is a "
                    "landform: retain strictly between 0 and 1, an outlet "
                    "notched but not cut to the basin floor, which Orogen "
                    "turns into a through-flowing valley with a residual lake. "
                    "A basin can be either, both or neither. They were briefly "
                    "the same word and that was a defect.",
            "note": "carved is the INTERSECTION: retain is the larger of the "
                    "two arms', so a rim survives wherever either climate "
                    "would have kept it. The disagreement is the BRACKETED set "
                    "by construction rather than by a tolerance chosen after "
                    "the fact, and width_pct_of_union is the bracket's own "
                    "width -- the honest uncertainty on this carve. Both arms "
                    "sit at the SAME flux and differ only in "
                    "model.land_albedo_source; the cold arm is a BOUND, NOT A "
                    "WORLD, and its climate must never be quoted as a "
                    "description of the planet.",
        }),
        "counts_by_overflow_test_alone": {
            "overflows_under_penman": int(carved.sum()),
            "overflows_under_land_evaporation": int(overflows_wet.sum()),
            "disputed": int(disputed.sum()),
            "note": "the two are NOT nested: without the land-rate floor, "
                    "Penman falls below the model's land evaporation wherever "
                    "the ground is wet, because land here is several times "
                    "rougher than open water. So this is agreement and "
                    "disagreement, not a bracket."},
        "basins": [
            {
                "id": ids[i],
                "verdict": str(verdict_name[i]),
                "retain": round(float(retain[i]), 4),
                # The two arms and whether they disagreed about this basin. The
                # bracket's membership, so a later pass can say whether a basin
                # that flipped had been bracketed or had been cut by both arms.
                # Null on a single-climate verdict, which is what
                # `intersection.single_climate` already says.
                "retain_warm_vegetated_arm": (
                    None if retain_endmember is None
                    else round(float(retain_primary[i]), 4)),
                "retain_cold_bare_rock_arm": (
                    None if retain_endmember is None
                    else round(float(retain_endmember[i]), 4)),
                "bracketed": (
                    None if retain_endmember is None
                    else bool((retain_primary[i] <= 0.0)
                              != (retain_endmember[i] <= 0.0))),
                "retain_incision": round(float(retain_incision[i]), 4),
                # The same quantity on the basis it was COMPUTED on: a fraction
                # of the FINISHED depression in metres. Orogen cannot spend it
                # -- see method.carve_list_basis -- and it is carried so the
                # conversion is auditable rather than implicit.
                "retain_incision_finished_depth_basis":
                    round(float(retain_incision_finished[i]), 4),
                "retained_fraction": (
                    None if not np.isfinite(retained_fraction[i])
                    else round(float(retained_fraction[i]), 4)),
                "retain_margin": round(float(retain_margin[i]), 4),
                # Whether the DELIVERY PHASE this pipeline does not carry could
                # change this basin's class. False is the strong statement: no
                # admissible seasonal concentration of the overflow moves it.
                # method.seasonal_concavity carries the bound and its argument.
                "seasonal_concavity_movable": bool(concavity_movable[i]),
                "retain_span_superseded": round(float(retain_span[i]), 4),
                "overflow_km3_per_year": round(float(q_pen[i]), 6),
                "overflow_m3_per_s": round(
                    float(q_pen[i]) * KM3_PER_YEAR_TO_M3_PER_S / year_s, 6),
                "sill_erodibility": round(float(sill_ero[i]), 4),
                "evaporation_margin": None if not np.isfinite(margin[i])
                                      else round(float(margin[i]), 4),
                "critical_aridity_index": round(float(crit[i]), 4),
                "aridity_index_penman": None if not np.isfinite(idx_pen[i])
                                        else round(float(idx_pen[i]), 4),
                "aridity_index_land_evaporation": None if not np.isfinite(idx_wet[i])
                                                  else round(float(idx_wet[i]), 4),
                "catchment_km2": round(float(basins.catchment_km2[i]), 2),
                "area_at_spill_km2": round(float(basins.area_at_spill_km2[i]), 2),
                "runoff_km_per_year": round(float(runoff[i]), 8),
                "precipitation_km_per_year": round(float(precip[i]), 8),
                "lake_evaporation_km_per_year": round(float(e_pen[i]), 8),
                "land_evaporation_km_per_year": round(float(e_wet[i]), 8),
                "no_catchment_runoff": bool(runoff[i] <= 0),
                # False where the conditioning left no depression. Every water
                # quantity on this record is then meaningless and the retain is
                # the 1 m guard, not a balance.
                "has_impoundment": bool(basins.has_impoundment[i]),
            } for i in range(basins.n)
        ] + [
            {"id": bid, "verdict": "carve", "retain": 0.0,
             "carried_from_previous_pass": True}
            for bid in sorted(carried)
        ],
    }
    args.out_json.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")

    mid = (retain > 0.0) & (retain < 1.0)
    n_trickle = int((carved & (retain > 0.0)).sum())
    n_lake_fed = int((carved & (runoff <= 0)).sum())
    if carried:
        print(f"carried   {len(carried):5d}  retain 0.0, from the previous pass")
    print(f"carve     {n_carve:5d}  retain 0.0  (this pass, of {basins.n} remaining)")
    print(f"          {n_trickle:5d}  outlet cut, but not through to the basin floor")
    print(f"          {n_lake_fed:5d}  overflow fed by the lake surface, no catchment runoff")
    if n_marginal:
        print(f"marginal  {n_marginal:5d}  retain {retain[mid].min():.3f}"
              f"-{retain[mid].max():.3f}, median {np.median(retain[mid]):.3f}")
    else:
        print(f"marginal  {n_marginal:5d}")
    print(f"preserve  {n_preserve:5d}  retain 1.0")
    print(f"basis     retain is a fraction of the NATURAL relief in model "
          f"units; {n_rebased} marginal entries rebased, "
          f"{n_verdict_moved} changed class, {n_without_ratio} had no ratio")
    print(f"\nwrote {args.out_list}\n      {args.out_json}")


if __name__ == "__main__":
    main()
