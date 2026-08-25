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

Retain is then the fraction of the impoundment that survives, which is what
Orogen cuts with and what a reader of the map sees.

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
from paths import climatology_path, rel
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

# The Earth calibration, which fixes the coefficient. Measured 2026-08-17 from
# HydroLAKES v1.0 joined to HydroBASINS level 5: natural lakes deeper than 5 m,
# above 1,000 km2, whose pour point sits in a basin with `ENDO == 0`, within 35
# degrees of the equator. Fifteen of them, over 78.9 Mkm2 of Earth land in that
# band. The size floor is what this mesh can resolve and it decides the answer:
# counting from 10 km2 instead gives a density 33 times higher and reads as
# "carve almost nothing". The latitude cut keeps out sills that were under an
# ice sheet 20,000 years ago and have had no time to be cut at all.
# `hydrography/notes/retain-fraction.md` carries the derivation.
EARTH_STANDING_BASINS = 15
EARTH_BAND_LAND_MKM2 = 78.9
CALIBRATION_LATITUDE = 35.0
COEFFICIENT_SEARCH = (1.0, 1.0e5)   # m per (m3/s)^0.5, the bisection bracket


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
    """
    q = np.asarray(q_km3_per_year, dtype=float) * KM3_PER_YEAR_TO_M3_PER_S / year_s
    cut = coefficient * np.asarray(erodibility, dtype=float) * np.power(
        np.clip(q, 0.0, None), INCISION_EXPONENT)
    if slope is not None:
        cut = cut * np.power(np.asarray(slope, dtype=float), slope_exponent)
    depth = np.maximum(np.asarray(depth_m, dtype=float), 1.0)
    return np.clip(1.0 - cut / depth, 0.0, 1.0)


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
    error on Earth's fifteen, which is the dominant uncertainty and is reported
    rather than hidden.
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

    per_mkm2 = land_band_mkm2 / EARTH_BAND_LAND_MKM2
    target = EARTH_STANDING_BASINS * per_mkm2
    spread = np.sqrt(EARTH_STANDING_BASINS)
    bracket = tuple(sorted(solve((EARTH_STANDING_BASINS + s) * per_mkm2)
                           for s in (-spread, spread)))
    coefficient = solve(target)
    return coefficient, bracket, target, standing(coefficient)


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

def climate_terms(clim_path, args, config, basins):
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
    with Dataset(args.climatology) as ds:
        am = cv.annual_mean
        pr, evap, mrro = am(ds, "pr"), -am(ds, "evap"), am(ds, "mrro")
        rss, rls = am(ds, "rss"), am(ds, "rls")
        diurnal = am(ds, "maxt") - am(ds, "mint")
        lsm = am(ds, "lsm")
    t_air, q_air, wind, p_air = cv.reference_level_air(args.climatology)

    # Through the one door, not by rebuilding the path from `model.resolution`:
    # that path is keyed by the RUNG alone and `surface_albedo` rewrites it per
    # BUILD. `--for-build` declares the deliberate cross-build read this export
    # makes once the carved build is staged, by NAMING the build it means.
    # world-z7bu, CLAUDE.md rule 5.
    staged_albedo = staged_surface_field(174, config, for_build=args.for_build)
    land_albedo = cv.read_sra_field(PROJECT_ROOT / staged_albedo["path"],
                                    *p_air.shape)
    penman_raw = cv.penman_open_water(
        t_air, q_air, wind, p_air, rss, rls, land_albedo,
        float(config["planet"]["gravity_m_s2"]), diurnal_range=diurnal)
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
            "mrro": mrro, "means": means, "year_s": year_s, "crit": crit, "runoff": runoff, "precip": precip, "idx_pen": idx_pen, "idx_wet": idx_wet, "q_pen": q_pen, "q_wet": q_wet, "carved": carved, "overflows_wet": overflows_wet, "disputed": disputed, "e_pen": e_pen, "e_wet": e_wet, "margin": margin, "retain_margin": retain_margin, "ocean_validation": ocean_validation, "penman_error_pct": penman_error_pct}
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
    args = ap.parse_args()

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

    basins = BasinSet(args.basins)
    resolution = str(config["model"]["resolution"]).upper()

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
    retain_incision = incision_retain(q_pen, year_s, basins.depth_at_spill_m,
                                      sill_ero, coefficient, slope=slope)
    print(f"outlet gradient   {slope_measured:5d} basins measured, "
          f"population {slope_reference:.5f} m/m")
    print(f"coefficient       {coefficient:7.1f} m per (m3/s)^0.5, bracket "
          f"{coefficient_bracket[0]:.1f}-{coefficient_bracket[1]:.1f}; "
          f"{standing_here} standing within {CALIBRATION_LATITUDE:g} deg "
          f"against Earth's {standing_target:.1f}")

    # Either is a reason to leave a rim standing: that the basin may not overflow
    # at all, or that its overflow cannot cut. Taking the larger keeps both.
    retain = np.maximum(retain_incision, retain_margin)

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
        retain_endmember = np.maximum(
            incision_retain(endmember["q_pen"], year_s,
                            basins.depth_at_spill_m, sill_ero, coefficient,
                            slope=slope),
            endmember["retain_margin"])
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
                              "cannot decide",
            "retain_incision_mapping":
                f"1 - {coefficient:.4g} * erodibility * slope**{SLOPE_EXPONENT:g} "
                f"* Q**{INCISION_EXPONENT:g} / depth_at_spill_m, clipped to [0, 1]",
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
                "earth_band_land_mkm2": EARTH_BAND_LAND_MKM2,
                "band_latitude_deg": CALIBRATION_LATITUDE,
                "vesper_band_land_mkm2": round(band_land_mkm2, 2),
                "target_standing_basins": round(standing_target, 2),
                "achieved_standing_basins": standing_here,
            },
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
                "retain_margin": round(float(retain_margin[i]), 4),
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
    print(f"\nwrote {args.out_list}\n      {args.out_json}")


if __name__ == "__main__":
    main()
