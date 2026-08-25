#!/usr/bin/env python3
"""Assess ExoPlaSim spin-up convergence from validated annual NetCDF files.

The verdict is taken over the last `--window` PRODUCTION orbits. Every segment
of a run declares what it was for, and orbits run to measure the model rather
than to advance the planet -- an I/O verification, a high-cadence wind sample, a
block on a differently patched binary -- are not evidence about where the run is
settling. `exoplasim/scripts/segments.py` owns that vocabulary.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

_MPL_CACHE = Path("/tmp/world-matplotlib-cache")
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy.optimize import curve_fit
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS
import close_state_energy
# lib/sea_water.py owns the four numbers salinity reaches the model through.
import sea_water
# CLAUDE.md names lib/autocorrelation.py as the one place a standard error is
# taken over a series with memory. Consecutive orbits of this model are not
# independent samples, and every uncertainty here that was scaled as though
# they were came back about half the size it should have been.
import autocorrelation as ac
# The one reader of the manifest's segment records; see
# exoplasim/scripts/segments.py for what a purpose means.
from segments import orbit_purposes, production_window
# CLAUDE.md names lib/sensitivity.py as the one flux-to-kelvin conversion, and
# the radiative damping this file relaxes at is that conversion inverted.
from sensitivity import planetary_albedo_from_fluxes, radiative_damping_w_m2_per_k


def output_files(run_dir: Path) -> list[Path]:
    outputs = []
    for path in run_dir.glob("MOST.*.nc"):
        match = re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)
        if match:
            outputs.append((int(match.group(1)), path))
    outputs.sort()
    if [year for year, _ in outputs] != list(range(len(outputs))):
        raise RuntimeError("Annual outputs are not contiguous from year zero")
    return [path for _, path in outputs]


def global_mean(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.sum(field * weights[None, :, None], axis=(-2, -1)) / (2.0 * field.shape[-1])


def slope(series: np.ndarray, window: int) -> float:
    y = series[-window:]
    x = np.arange(window, dtype=float)
    return float(np.polyfit(x, y, 1)[0])


def slope_standard_error(series: np.ndarray, window: int) -> float:
    """The standard error of that slope, with the residual memory in it.

    An ordinary least-squares slope over n points has variance
    sigma^2 / sum((x - xbar)^2) ONLY when the residuals are independent. They
    are not here, and the correction is the integrated autocorrelation time of
    the residuals: var(slope) = sigma^2 * tau / sum((x - xbar)^2).

    WHY IT MATTERS RATHER THAN BEING A REFINEMENT. The offset criterion falls
    back to `slope * tau_expected` on a converged run, and tau_expected is
    about ten orbits, so the slope's error is multiplied by ten before it is
    compared with a 0.15 K allowance. A slope error that is understated by two
    is a criterion whose statistic wanders by two, which is what made the
    verdict on this project's T21 baseline flip seven times on a state that had
    stopped moving. world-omn.
    """
    y = np.asarray(series[-window:], dtype=float)
    n = y.size
    if n < 4:
        return float("nan")
    x = np.arange(n, dtype=float)
    fit = np.polyfit(x, y, 1)
    residual = y - np.polyval(fit, x)
    # Two parameters fitted, so n - 2 degrees of freedom.
    variance = float(np.dot(residual, residual)) / (n - 2)
    tau = ac.integrated_time(residual)["tau"]
    sxx = float(np.dot(x - x.mean(), x - x.mean()))
    return float(np.sqrt(variance * tau / sxx))


def orbits_for_slope_standard_error(scatter: float, tau: float,
                                    target: float) -> float:
    """Window length at which a fitted slope's standard error reaches `target`.

    var(slope) = scatter^2 * tau * 12 / (n (n^2 - 1)), so this inverts a cubic.
    It is what prices a window: the answer to "how many orbits does this
    criterion need before it is measuring the planet rather than the model's
    own variability", in the units the run is bought in.
    """
    if target <= 0 or not all(np.isfinite([scatter, tau])) or scatter <= 0:
        return float("nan")
    need = (scatter ** 2) * max(tau, 1.0) * 12.0 / (target ** 2)
    n = max(4.0, need ** (1.0 / 3.0))
    for _ in range(64):                       # n(n^2-1) = need, by iteration
        n = (need + n) ** (1.0 / 3.0)
    return float(n)



# The mixed layer's heat capacity and the radiative damping. Together they set
# how long an approach takes, which is what turns a drift rate into a remaining
# offset -- and the offset is what `OFFSET_TOLERANCE_K` passes or fails a run
# on, so neither may be a number typed in beside the code that uses it.
#
# THE SLAB IS MODELLING THE MODEL'S OWN MIXED LAYER, so it takes the model's own
# sea water -- and it takes it FROM THE RUN BEING ASSESSED. The 1025 and 3990
# that stood here were uncited and 0.5% and 4.6% off the model.
#
# It cannot be a module constant. `sea_water.constants()` with no run directory
# returns what `icemod.f90` declares today, and the model's CPS has already
# moved once, by 4.54 per cent. Assessing a run older than that change against
# today's declaration rescales `tau_expected`, which sets the fallback remaining
# offset, which is what `OFFSET_TOLERANCE_K` fails a run on -- so a constant
# read at import turns a model change into a silent re-verdict on every run that
# predates it. Every run this project writes declares all four keys in its own
# `icemod_namelist`, so the run's own value is always available.
import yaml as _yaml
_MLD = float(_yaml.safe_load(
    (Path(__file__).resolve().parents[2] / "config" / "planet.yaml")
    .read_text(encoding="utf-8"))["surface"]["mixed_layer_depth_m"])


def slab_heat_capacity(run_dir: Path) -> tuple[float, dict]:
    """The mixed layer's heat capacity per square metre, and the sea water it used.

    `_MLD` is the configured design depth of the ocean mixed layer, which is
    what an equilibration time is a property of; `CRHOS` and `CPS` come from the
    run. Returns both so the assessment can record what it was built from.
    """
    water = sea_water.constants(run_dir)
    return _MLD * water["CRHOS"] * water["CPS"], water


# THE WINDOW, DERIVED RATHER THAN PICKED.
#
# The criterion that bounds the ANSWER rather than a rate is the extrapolated
# offset, and on a settled run the exponential fit has nothing to grip, so the
# fallback tests `slope * tau_expected` against 0.15 K. The slope's own error
# is therefore multiplied by the relaxation time before the comparison, and the
# window that criterion needs follows from
#
#     var(slope) = scatter^2 * tau * 12 / (n (n^2 - 1))
#
# with the bar that a threshold discriminates only when its statistic's
# standard error is at most a third of it (`RESOLVING_FACTOR`, declared in
# `main` before it was applied to anything).
#
# THE INPUTS ARE PRIOR MEASUREMENTS, EACH WITH A SOURCE, and the number below
# is computed from them rather than typed, so nobody can move the window
# without moving an input that has a citation. A ten-orbit window was in use
# and it does not resolve this criterion: its statistic wandered between 0.031
# and 0.255 against a 0.15 K allowance on a T21 baseline that had stopped
# moving, and the verdict flipped seven times over twenty orbits. world-omn.
#
# THE ANSWER IS A BRACKET AND THIS IS ITS ESTIMATE. At tau = 1, the value if
# consecutive orbits were independent, the window is 21 orbits; at tau = 9, the
# value for a lag-1 correlation of 0.8, it is 44. The factor tau is what the
# model's memory costs and everything else is what an independent series would
# have needed anyway. `exoplasim/notes/convergence-lengths.md` carries the
# bracket; every assessment reports the window ITS OWN run would need, which is
# what follows the resolution up the ladder.
NOMINAL_ORBIT_SCATTER_K = 0.07     # the T21 baseline's stationary spread, world-omn
NOMINAL_TAU_ORBITS = 4.2           # (1+r)/(1-r) at r = 0.615, world-yj9o
NOMINAL_RELAXATION_ORBITS = 10.0   # this planet's slab; each run reports its own
# The same 0.15 K and the same factor of three the criteria use, restated here
# only because the default has to exist before `main` runs. `main` asserts they
# agree, so the two cannot drift.
_OFFSET_TOLERANCE_K = 0.15
_RESOLVING_FACTOR = 3.0
# THE STATISTIC IS `|offset| + half_width`, NOT `offset`. Both terms are the
# same slope error multiplied by the same relaxation time, so the quantity
# compared with the tolerance has about twice the offset's scale. Sizing the
# window on the offset alone would be sizing it for half the statistic.
_OFFSET_STATISTIC_TERMS = 2.0


def window_for_offset_criterion(scatter: float, tau: float,
                                relaxation_orbits: float) -> float:
    """The window at which the offset criterion starts to discriminate."""
    return orbits_for_slope_standard_error(
        scatter, tau, _OFFSET_TOLERANCE_K
        / (_RESOLVING_FACTOR * _OFFSET_STATISTIC_TERMS * relaxation_orbits))


DEFAULT_WINDOW_ORBITS = int(np.ceil(window_for_offset_criterion(
    NOMINAL_ORBIT_SCATTER_K, NOMINAL_TAU_ORBITS, NOMINAL_RELAXATION_ORBITS)))


# HOW LONG A RUN TAKES TO GET HERE, AS OPERATIONAL EXPERIENCE. These are what
# this project has repeatedly seen, not a measured distribution, and they are
# marked as experience wherever they are quoted so nobody later reads them as
# an artifact-backed result. They belong beside the criteria because the
# criteria are what decides a run is finished, and a window that is a large
# fraction of the approach is testing the approach.
#
#   COLD START AT T21: about seventy orbits, seed- and initial-condition-
#   dependent. `exoplasim/notes/parameter-decisions.md` records one converging
#   on all six criteria after 70.
#
#   RECONVERGENCE after a timestep change or a resolution conversion: generally
#   ten to twenty orbits. The state is already at a climate; what is settling
#   is the model's response to a changed discretisation.
#
# The consequence for a window: a 10-orbit window is a seventh of a cold
# start's approach and the whole of a reconvergence, so on the second it is
# assessing orbits that are still moving by construction.
CONVERGENCE_LENGTHS = {
    "basis": "operational experience across this project's runs, not a "
             "measured distribution; every number here is bracketed by what "
             "has been seen rather than fitted",
    "cold_start_t21_orbits": 70,
    "cold_start_note": "seed- and initial-condition-dependent; "
                       "exoplasim/notes/parameter-decisions.md records a cold "
                       "start converging on all six criteria after 70 orbits",
    "reconvergence_after_timestep_or_resolution_change_orbits": [10, 20],
    "reconvergence_note": "the state is already at a climate; what settles is "
                          "the response to a changed discretisation",
    "recorded_in": "exoplasim/notes/convergence-lengths.md",
}


def relaxation_orbits(orbital_year_days: float, feedback_w_m2_k: float,
                      slab_heat_capacity_j_m2_k: float) -> float:
    """Orbits for an e-folding of the slab's approach to equilibrium.

    `feedback_w_m2_k` is the radiative damping, and it comes from
    `lib/sensitivity.py` -- the one flux-to-kelvin conversion -- evaluated at
    the planetary albedo THIS RUN reports, so it moves with the run the way
    `year_days` already does. A private constant stood here instead, 1.31
    W/m2/K commented "measured, not assumed" with nothing saying where, 11%
    above what the module's own slope implies.
    """
    return ((slab_heat_capacity_j_m2_k / feedback_w_m2_k)
            / (orbital_year_days * 86400.0))


def approach_to_equilibrium(orbits: np.ndarray, series: np.ndarray,
                            fraction: float = 0.35):
    """Fit T(n) = T_inf - A exp(-n/tau) and return the asymptote with an interval.

    A drift rate is not a distance. An exponential approach at rate `d` with time
    constant `tau` still has `d * tau` to travel, all of it one-signed, so a
    threshold on the rate bounds nothing on its own. At this planet's tau of
    about ten orbits, the 0.05 K/orbit criterion permits roughly half a kelvin of
    remaining approach -- the same order as corrections the project has thought
    worth applying, and larger than some.

    The interval is what makes this usable rather than a second point estimate.
    A fit whose tau is comparable to the run length is extrapolating past its
    data, and says so through a wide interval instead of a confident wrong
    number; one such fit in this project's history sat 0.6 K from its own
    reported value.

    THE INTERVAL IS AUTOCORRELATION-CORRECTED. `curve_fit` returns a covariance
    scaled as though the residuals were independent draws, and orbits of this
    model are not; the reported half width is widened by the root of the
    residuals' integrated autocorrelation time, which is the factor that was
    missing.

    Returns (asymptote, half_width, tau_orbits) or (nan, nan, nan) if the series
    will not support a fit.
    """
    mask = orbits >= orbits.max() * fraction
    x, y = orbits[mask], series[mask]
    if x.size < 8:
        return float("nan"), float("nan"), float("nan")
    model = lambda n, inf, amp, tau: inf - amp * np.exp(-n / tau)
    try:
        popt, pcov = curve_fit(model, x, y,
                               p0=[y[-1], max(y[0] - y[-1], 1e-3), max(x.size / 3, 1.0)],
                               maxfev=40000)
    except Exception:
        return float("nan"), float("nan"), float("nan")
    inf, tau = float(popt[0]), float(popt[2])
    err = float(np.sqrt(np.diag(pcov))[0]) if np.all(np.isfinite(pcov)) else float("nan")
    # `curve_fit` SCALES ITS COVARIANCE AS THOUGH THE RESIDUALS WERE
    # INDEPENDENT. They are not: consecutive orbits of this model carry memory,
    # so the residual variance is spread over fewer independent samples than
    # the fit counts and every diagonal of pcov comes back too small by the
    # integrated autocorrelation time. The standard error goes as its root.
    # This is the same defect `compare_equilibria.py` had in a different
    # dress -- there the count was the window's, here it is the fit's --
    # and correcting it can only widen the interval, which is the direction
    # a convergence test must err in. world-yj9o.
    residual = y - model(x, *popt)
    tau_residual = ac.integrated_time(residual)["tau"]
    if np.isfinite(err):
        err *= float(np.sqrt(tau_residual))
    # A fit extrapolating beyond its own span is reported as such rather than
    # trusted: widen the interval by how far past the data the asymptote sits.
    if np.isfinite(tau) and tau > 0:
        overshoot = max(0.0, tau / max(x.size, 1))
        err = (err if np.isfinite(err) else abs(inf - y[-1])) * (1.0 + overshoot)
    return inf, err, tau


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW_ORBITS,
                        help="orbits in the test window. Counted back from the "
                             "last PRODUCTION orbit, not the last orbit. The "
                             f"default, {DEFAULT_WINDOW_ORBITS}, is derived: it "
                             "is the shortest window at which the offset "
                             "criterion's statistic has a standard error of a "
                             "third of its threshold, at the scatter and "
                             "autocorrelation this model has been measured "
                             "with. A shorter one does not test that criterion, "
                             "it flips on it")
    parser.add_argument("--through", type=int, default=None, metavar="ORBITS",
                        help="assess the run AS IF it had stopped after this "
                             "many orbits, ignoring everything later. For "
                             "asking when a run FIRST met the criteria rather "
                             "than whether it meets them now -- which is what "
                             "a relaxation time is, and what decides whether "
                             "converting a coarse state into a finer one is "
                             "worth the wall clock. Its report carries the "
                             "truncation in its name and it writes nothing "
                             "into the run's manifest, so a sweep cannot "
                             "overwrite the run's own verdict in either "
                             "place.")
    parser.add_argument("--output", type=Path, default=ANALYSIS / "convergence")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    files = output_files(run_dir)
    if args.through is not None:
        if args.through > len(files):
            raise RuntimeError(
                f"--through {args.through} but the run has {len(files)} orbits")
        files = files[:args.through]
    if len(files) < args.window:
        raise RuntimeError("Not enough annual outputs for requested window")

    # THE WINDOW IS THE LAST `window` PRODUCTION ORBITS, NOT THE LAST `window`
    # ORBITS. A segment declares what it was for, and a diagnostic one -- an I/O
    # verification, a high-cadence wind sample, a block on a differently patched
    # binary -- is not evidence about where the run is settling. Taking the tail
    # blindly would have averaged three such orbits into this project's baseline
    # window on the day CLIM-9 was written. `production_window` drops a
    # diagnostic tail and REFUSES a window with a diagnostic hole in it, because
    # a slope across a gap is not a trend.
    window_start, window_end = production_window(run_dir, len(files), args.window)
    purposes = orbit_purposes(run_dir, range(len(files)))
    dropped = sorted(o for o in range(window_end + 1, len(files)))

    records = []
    for year, path in enumerate(files):
        with Dataset(path) as nc:
            weights = leggauss(len(nc.dimensions["lat"]))[1][::-1]
            record = {"year_index": year}
            # rst and rsut are read for the planetary albedo the radiative
            # damping is evaluated at: the run's own, not a declared one.
            for name in ["ts", "ntr", "hfns", "sic", "pr", "rst", "rsut"]:
                field = np.asarray(nc[name][:], dtype=float)
                value = float(global_mean(field, weights).mean())
                if name == "pr":
                    value *= 86400.0 * 1000.0
                record[name] = value
            records.append(record)

    # Every array is truncated at the last production orbit, so the trailing
    # window, the exponential fit and the storage term all describe the same
    # orbits. `records` keeps all of them, and the plot draws all of them, so
    # nothing is hidden -- only excluded from the verdict.
    arrays = {key: np.array([record[key] for record in records[:window_end + 1]])
              for key in ["ts", "ntr", "hfns", "sic", "pr", "rst", "rsut"]}
    w = args.window
    metrics = {
        "temperature_slope_k_per_orbit": slope(arrays["ts"], w),
        "toa_balance_slope_w_m2_per_orbit": slope(arrays["ntr"], w),
        "surface_balance_slope_w_m2_per_orbit": slope(arrays["hfns"], w),
        "sea_ice_slope_fraction_per_orbit": slope(arrays["sic"], w),
        "mean_toa_balance_w_m2": float(arrays["ntr"][-w:].mean()),
        "mean_surface_balance_w_m2": float(arrays["hfns"][-w:].mean()),
        "temperature_mean_k": float(arrays["ts"][-w:].mean()),
        "temperature_last_k": float(arrays["ts"][-1]),
        "sea_ice_mean_fraction": float(arrays["sic"][-w:].mean()),
    }

    # The storage the criterion now passes on, from the same window. Computed by
    # close_state_energy.state_energy rather than reimplemented, because the
    # criterion turning on a quantity is the strongest possible reason not to
    # have two versions of it.
    #
    # FAILS CLOSED. If the window cannot be closed against the state -- a run
    # whose namelist does not declare its own constants, a window too short for a
    # trend -- `storage` stays None and the criterion is False. A convergence
    # check that cannot measure the thing it tests must refuse, not abstain.
    storage = None
    storage_error = None
    try:
        storage = close_state_energy.state_energy(
            run_dir, window_start, window_end)
    # SystemExit is in the tuple because `close_state_energy` raises it, not
    # Exception, for the case this comment names first: a run whose namelist
    # does not declare GA, GASCON, GSOL0 and ECCEN. With only `Exception` here
    # the assessment died on that run and wrote no report at all, so the
    # criterion could not be False -- there was nothing to be False in.
    except (Exception, SystemExit) as exc:        # noqa: BLE001 - reported, not raised
        storage_error = f"{type(exc).__name__}: {exc}"

    if storage is not None:
        metrics.update({
            "state_storage_w_m2": storage["storage_w_m2_least_squares"],
            "state_storage_endpoint_w_m2": storage["storage_w_m2_endpoint"],
            "surface_storage_w_m2": storage["surface_storage_w_m2"],
            # Recorded, not tested. This is the gap CLIM-1 is open on: the
            # reported TOA net minus what the state actually stores. It is a
            # structural property of the diagnostic rather than of any run, so it
            # belongs in every assessment until it is explained.
            "reported_toa_minus_storage_w_m2":
                metrics["mean_toa_balance_w_m2"] - storage["storage_w_m2_least_squares"],
        })
    else:
        metrics["state_storage_unavailable"] = storage_error

    # The quantity the design band is stated in, and therefore the one that has
    # to be bounded. Everything above is a rate; this is a distance.
    # (tau_expected is needed by the fallback below, so it is computed first.)
    # The year comes from the run's own manifest, not from the current config:
    # this is a property of the run being assessed, which may predate a
    # baseline re-run that moved the orbit.
    manifest_path = run_dir / "run_manifest.json"
    year_days = 182.8
    if manifest_path.is_file():
        year_days = float(json.loads(manifest_path.read_text(encoding="utf-8"))
                          ["derived_parameters"]["orbital_year_earth_days"])
    # The albedo is this run's own, over the same window every other metric is
    # taken over, so the damping moves with the run exactly as `year_days` does.
    window_alpha = planetary_albedo_from_fluxes(
        float(arrays["rst"][-w:].mean()), float(arrays["rsut"][-w:].mean()))
    feedback_w_m2_k = radiative_damping_w_m2_per_k(window_alpha)
    slab_capacity, slab_water = slab_heat_capacity(run_dir)
    tau_expected = relaxation_orbits(year_days, feedback_w_m2_k, slab_capacity)
    orbits_axis = np.arange(len(arrays["ts"]), dtype=float)
    asymptote, half_width, tau_fit = approach_to_equilibrium(orbits_axis, arrays["ts"])
    offset = asymptote - metrics["temperature_mean_k"]

    # A converged run has no approach left to fit, so the exponential becomes
    # unconstrained and can return anything: one run here reported an asymptote
    # of -23064 K with an interval of six million. That is the fit having nothing
    # to grip, not the run being far from equilibrium, and the two must not be
    # treated alike -- the first version of this criterion failed a converged run
    # for it.
    #
    # When the fit is unusable, fall back to the bound the drift itself implies,
    # offset ~ drift * tau, using the EXPECTED tau rather than a fitted one, and
    # carry 100% uncertainty on it. That is conservative in the right direction:
    # it can only refuse a run, never pass one it should not.
    fit_usable = (np.isfinite(asymptote) and np.isfinite(half_width)
                  and abs(offset) < 20.0 and half_width < 5.0)
    slope_se = slope_standard_error(arrays["ts"], w)
    if not fit_usable:
        offset = metrics["temperature_slope_k_per_orbit"] * tau_expected
        # THE FALLBACK'S HALF WIDTH IS TWO UNCERTAINTIES, NOT ONE. Carrying
        # 100% of the offset covers tau_expected being wrong, and it was the
        # whole interval; but the offset is a fitted SLOPE times tau_expected,
        # and the slope over a short window has an error of its own that the
        # multiplication magnifies by tau_expected. On a settled run that
        # second term is the larger of the two, and leaving it out is what let
        # `|offset| + half_width` wander across its own allowance on a state
        # that had stopped moving. Added in quadrature: the two are
        # independent, and this reduces to the old interval exactly where the
        # slope is well determined. world-omn.
        half_width = float(np.hypot(offset, tau_expected * slope_se))
        asymptote = metrics["temperature_mean_k"] + offset
    metrics.update({
        "temperature_asymptote_k": asymptote,
        "temperature_asymptote_half_width_k": half_width,
        "temperature_remaining_offset_k": offset,
        "relaxation_orbits_fitted": tau_fit,
        "relaxation_orbits_expected": tau_expected,
        # What tau_expected was built from, so the fallback offset below can be
        # audited without re-deriving it. lib/sensitivity.py owns the damping.
        "planetary_albedo_in_window": window_alpha,
        "radiative_damping_w_m2_per_k": feedback_w_m2_k,
        "slab_heat_capacity_j_m2_k": slab_capacity,
        # Which sea water, and whether it came from the run or the
        # compiled model. An artifact that does not say cannot be
        # checked against the run it describes.
        "slab_sea_water": slab_water,
        "slab_mixed_layer_depth_m": _MLD,
        "remaining_offset_implied_by_drift_k":
            metrics["temperature_slope_k_per_orbit"] * tau_expected,
        # The slope errors, with the residual memory in them. Reported for
        # every slope the criteria threshold, so nobody has to take a
        # threshold's word for what the window can see.
        "temperature_slope_standard_error_k_per_orbit": slope_se,
        "toa_balance_slope_standard_error_w_m2_per_orbit":
            slope_standard_error(arrays["ntr"], w),
        "surface_balance_slope_standard_error_w_m2_per_orbit":
            slope_standard_error(arrays["hfns"], w),
        "sea_ice_slope_standard_error_fraction_per_orbit":
            slope_standard_error(arrays["sic"], w),
    })

    # 0.15 K against a 3 K design band: small enough that the band's edges mean
    # what they say, loose enough to be reachable. Stated before it was applied.
    OFFSET_TOLERANCE_K = 0.15

    # WHEN A THRESHOLD IS MEASURING THE PLANET AND NOT THE INSTRUMENT. A
    # criterion whose statistic has a standard error comparable to its own
    # threshold does not discriminate; it flips. The bar, fixed here before it
    # was applied to any assessment: the standard error of the statistic must
    # be at most a third of the threshold it is compared with. Three is the
    # smallest factor at which a statistic sitting on the threshold is more
    # than a chance excursion from either side of it.
    RESOLVING_FACTOR = 3.0
    # The default window was sized from these two before `main` ran. One copy
    # of each number, checked rather than trusted.
    assert (OFFSET_TOLERANCE_K, RESOLVING_FACTOR) == (_OFFSET_TOLERANCE_K,
                                                      _RESOLVING_FACTOR)

    # THE ENERGY-BALANCE CRITERION TESTS STORAGE, NOT REPORTED TOA. Changed
    # 2026-08-17, deliberately and with the reasoning recorded, because the two
    # are not the same number here and the criterion was thresholding the one
    # that is wrong. CLIM-7; the evidence is exoplasim/notes/baseline-equilibration.md.
    #
    # `ntr` is the model's own statement about whether the planet is gaining or
    # losing energy, and it carries a structural offset of -0.573 +/- 0.035 W/m2
    # measured across four runs spanning 7.8 K. Conservation gives the same
    # quantity from the prognostic state instead -- mixed layer, sea ice and snow
    # as latent heat, soil, atmospheric enthalpy, column vapour -- and THAT is
    # what a convergence criterion is trying to bound. Where the energy the
    # offset represents actually goes is CLIM-1 and is still open, so both are
    # recorded here and the difference is carried as its own metric.
    #
    # THE THRESHOLD IS DERIVED, and the 0.5 W/m2 it replaces was not: that number
    # was picked early, never revisited, and applied to two different quantities.
    # Two independent bounds, stated before this was run against any assessment:
    #
    #   What matters. A residual imbalance X leaves the world short of its
    #   asymptote by roughly X * (dT/dt per W/m2) * tau. Measured on this model
    #   the mixed layer drifts 0.128 K/orbit per W/m2 and tau is 9.9 orbits, so
    #   holding that inside the 0.15 K offset tolerance already committed to
    #   above needs X < 0.15 / (0.128 * 9.9) = 0.118 W/m2. Deriving it from the
    #   temperature tolerance rather than inventing a second number is the point:
    #   the two criteria now bound the same thing in two units.
    #
    #   What is measurable. The storage estimate on a 10-orbit block differs from
    #   the same run's 20-orbit block by 0.111, 0.115 and 0.116 W/m2 on the three
    #   runs where both exist. A threshold below about 0.11 is measuring the
    #   sampling noise of the estimator rather than the planet. That figure was
    #   taken on 10-orbit blocks and the default window is now wider, so it is a
    #   CEILING on the estimator's noise here rather than a live floor; the
    #   binding bound is the temperature-derived one either way and 0.12 stands.
    #
    # The two land within 6% of each other, which is the argument for the number
    # rather than a coincidence to note: below 0.11 is unresolvable, above 0.118
    # admits more drift than the temperature criterion already forbids. 0.12 is
    # where they meet, and it is 4x tighter than what it replaces.
    STORAGE_TOLERANCE_W_M2 = 0.12

    criteria = {
        "abs_temperature_slope_lt_0.05_k_per_orbit": abs(metrics["temperature_slope_k_per_orbit"]) < 0.05,
        "abs_toa_slope_lt_0.05_w_m2_per_orbit": abs(metrics["toa_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_surface_slope_lt_0.05_w_m2_per_orbit": abs(metrics["surface_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_sea_ice_slope_lt_0.001_per_orbit": abs(metrics["sea_ice_slope_fraction_per_orbit"]) < 0.001,
        f"abs_state_storage_lt_{STORAGE_TOLERANCE_W_M2}_w_m2": bool(
            storage is not None
            and abs(storage["storage_w_m2_least_squares"]) < STORAGE_TOLERANCE_W_M2),
        # The one that bounds the answer rather than its rate of change. Fails
        # closed: a fit that will not converge is not evidence of equilibrium.
        f"extrapolated_offset_lt_{OFFSET_TOLERANCE_K}_k": bool(
            np.isfinite(offset) and np.isfinite(half_width)
            and abs(offset) + half_width < OFFSET_TOLERANCE_K),
    }
    # WHAT THIS WINDOW CAN SEE, reported beside every verdict rather than
    # assumed. Each row is a threshold, the standard error of the statistic it
    # tests, and the window at which that error would fall to a third of it.
    #
    # THE COST OF THE MEMORY IS THE FACTOR TAU. An independent series would
    # need `orbits_for_slope_standard_error` with tau = 1; the model has
    # variability at the window's own timescale, so it needs tau times the
    # variance and the window grows as the cube root of that. The required
    # windows below are computed from THIS run's own residual scatter and its
    # own autocorrelation, so they follow the resolution rather than being a
    # number swept once at T21 and carried up the ladder.
    ts_residual = arrays["ts"][-w:] - np.polyval(
        np.polyfit(np.arange(w, dtype=float), arrays["ts"][-w:], 1),
        np.arange(w, dtype=float))
    ts_scatter = float(np.std(ts_residual, ddof=1)) if w > 2 else float("nan")
    ts_memory = (ac.integrated_time(ts_residual) if w > 3
                 else {"tau": float("nan"), "reliable": False,
                       "lag1": float("nan")})
    ts_tau = ts_memory["tau"]
    resolving = {
        "factor": RESOLVING_FACTOR,
        "criterion": ("a threshold discriminates when the standard error of "
                      "the statistic it tests is at most 1/factor of it"),
        "window_orbits": w,
        "temperature_residual_scatter_k": ts_scatter,
        "temperature_residual_tau_orbits": ts_tau,
        "temperature_residual_lag1": ts_memory["lag1"],
        # TAU MEASURED INSIDE THE WINDOW IS A LOWER BOUND ON TAU. A span
        # comparable to the correlation time cannot resolve it: on synthetic
        # series with a known answer, a twenty-sample window recovers about
        # half the true value and a ten-sample window barely more than one.
        # Every number below inherits that, so the required window is a FLOOR
        # and is labelled one until a stationary span long enough to carry the
        # estimate exists. lib/autocorrelation.py.
        "tau_estimated_from_the_window_itself": True,
        "tau_span_supports_the_estimate": bool(ts_memory["reliable"]),
        "required_window_is_a_lower_bound": not bool(ts_memory["reliable"]),
        "rows": [],
    }
    for name, statistic_se, threshold in (
            ("abs_temperature_slope", metrics[
                "temperature_slope_standard_error_k_per_orbit"], 0.05),
            ("abs_toa_slope", metrics[
                "toa_balance_slope_standard_error_w_m2_per_orbit"], 0.05),
            ("abs_surface_slope", metrics[
                "surface_balance_slope_standard_error_w_m2_per_orbit"], 0.05),
            ("abs_sea_ice_slope", metrics[
                "sea_ice_slope_standard_error_fraction_per_orbit"], 0.001),
            # The offset criterion tests a slope multiplied by tau_expected, so
            # its statistic's error is that multiple of the slope's.
            # The statistic is `|offset| + half_width` and both terms carry the
            # same slope error times tau_expected, so its scale is twice the
            # offset's. Sizing on the offset alone sizes for half the test.
            ("extrapolated_offset",
             _OFFSET_STATISTIC_TERMS * tau_expected * slope_se,
             OFFSET_TOLERANCE_K)):
        resolves = bool(np.isfinite(statistic_se)
                        and statistic_se * RESOLVING_FACTOR <= threshold)
        resolving["rows"].append({
            "criterion": name, "threshold": threshold,
            "statistic_standard_error": statistic_se,
            "resolves": resolves})
    resolving["window_orbits_for_offset_criterion"] = window_for_offset_criterion(
        ts_scatter, ts_tau, tau_expected)
    # The same window if the orbits were independent, so the price of the
    # memory is visible rather than folded into one number.
    resolving["window_orbits_for_offset_criterion_if_independent"] = \
        window_for_offset_criterion(ts_scatter, 1.0, tau_expected)
    resolving["default_window_orbits"] = DEFAULT_WINDOW_ORBITS

    report = {
        "run_dir": str(run_dir),
        "completed_orbits": len(files),
        "window_orbits": w,
        # Which orbits the verdict is actually about. `completed_orbits` above
        # is what the run contains; these are what was assessed, and the two
        # differ whenever a diagnostic segment sits at the end.
        "window_start_year_index": window_start,
        "window_end_year_index": window_end,
        "orbits_after_window_excluded_as_non_production": dropped,
        "segment_purposes": {str(k): v for k, v in sorted(purposes.items())},
        "metrics": metrics,
        "criteria": criteria,
        "criteria_provenance": {
            "energy_balance_tests": "state storage, not reported top-of-atmosphere net",
            "storage_tolerance_w_m2": STORAGE_TOLERANCE_W_M2,
            "storage_tolerance_derivation":
                "0.15 K offset tolerance / (0.128 K per orbit per W/m2 * 9.9 orbits) "
                "= 0.118 W/m2, against an estimator that resolves 0.11 W/m2 on a "
                "10-orbit block. The two meet at 0.12.",
            "supersedes": "abs_mean_toa_lt_0.5_w_m2 and abs_mean_surface_lt_0.5_w_m2",
            "why": "The reported TOA net carries a structural offset of -0.573 "
                   "+/- 0.035 W/m2 across four runs, so the criterion was "
                   "thresholding a diagnostic rather than the planet. The 0.5 "
                   "W/m2 it replaces was chosen early, never revisited, and "
                   "applied to two different quantities. CLIM-7, decided "
                   "2026-08-17; the offset itself is CLIM-1 and is still open, "
                   "which is why both numbers and their difference are recorded.",
        },
        # REPORTED, NOT THRESHOLDED, and deliberately so. The criteria already
        # fail closed on an interval too wide to clear the allowance, which is
        # what a window that cannot resolve produces; adding a sixth criterion
        # on the same fact would refuse the same runs twice. What this buys is
        # the number a reader would otherwise have to sweep for: how many
        # orbits the window needs, at THIS resolution and this run's own
        # variability.
        "resolving_power": resolving,
        "convergence_lengths": CONVERGENCE_LENGTHS,
        "sufficiently_equilibrated_for_worldbuilding": bool(all(criteria.values())),
        "failed_criteria": sorted(name for name, ok in criteria.items() if not ok),
        "annual_records": records,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    # Named by the run, not fixed. A fixed filename meant every assessment
    # overwrote the last, so a run's convergence record could not survive the
    # next run being assessed -- and assessing a run silently destroyed the
    # evidence for a previous one.
    #
    # A TRUNCATED assessment carries its length in the name and never takes the
    # run's own. Sweeping `--through` to find when a run first converged would
    # otherwise leave the verdict of whichever truncation ran last standing as
    # the run's, which is the same defect one directory up.
    suffix = "" if args.through is None else f"_through{args.through:03d}"
    report_path = args.output / f"{run_dir.name}_convergence{suffix}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    # A run that has been assessed and misses is NOT "spinup_in_progress", and
    # leaving it labelled that way is how the baseline this project built its
    # production climatology from came to be recorded as still spinning up.
    #
    # This does not round a miss into a pass. `sufficiently_equilibrated_for_
    # worldbuilding` stays strictly all-or-nothing and is what any consumer
    # should test. What is added is the third state the prose convention has
    # always used and the code never had: quasi-equilibrated, carrying WHICH
    # criteria failed and the metric each failed on, so the label can never be a
    # hand-wave. One cold case in this project is called quasi-equilibrated for
    # missing by 0.004 W/m2; that standard is preserved by recording the margin,
    # not by widening the threshold.
    failed = report["failed_criteria"]
    # A TRUNCATED ASSESSMENT NEVER TOUCHES THE RUN'S MANIFEST. `--through`
    # already keeps its report under its own filename so a sweep cannot
    # overwrite the run's verdict one directory up -- and then wrote the same
    # truncated verdict straight into the run itself: its status, its
    # equilibrium cutoff and the whole assessment payload, taken from a
    # prefix of the run. Sweeping `--through` to find when a run FIRST
    # converged left the run recorded as having converged wherever the sweep
    # happened to stop last, which is the exact defect the suffix was added to
    # prevent. The sweep is how the window is priced, so this is not a corner.
    manifest_path = run_dir / "run_manifest.json"
    if args.through is not None:
        manifest_path = None
    if manifest_path is not None and manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if report["sufficiently_equilibrated_for_worldbuilding"]:
            manifest["status"] = "equilibrated_for_worldbuilding"
        else:
            manifest["status"] = "quasi_equilibrated"
        # The last orbit this verdict is ABOUT, which is the end of the
        # production window and not the end of the run. A diagnostic tail was
        # not assessed, so it must not be counted as settled.
        manifest["equilibrium_cutoff_year_index"] = window_end
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )

    # Draw every orbit, including the ones outside the verdict, and mark them.
    # A plot that showed only the assessed orbits would make an excluded tail
    # invisible, which is the failure this whole change is about.
    years = np.arange(len(files))
    all_series = {key: np.array([record[key] for record in records])
                  for key in ["ts", "ntr", "hfns", "sic", "pr"]}
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    panels = [
        ("ts", "Global surface temperature", "K"),
        ("ntr", "TOA net radiation", "W m$^{-2}$"),
        ("hfns", "Surface downward heat flux", "W m$^{-2}$"),
        ("sic", "Planetary sea-ice fraction", "fraction"),
    ]
    for ax, (key, title, units) in zip(axes.ravel(), panels):
        ax.plot(years, all_series[key], marker="o", markersize=2.5, linewidth=1.1)
        ax.axvspan(window_start, window_end, color="tab:orange", alpha=0.12,
                   label=f"{w}-orbit test window")
        if dropped:
            ax.axvspan(dropped[0] - 0.5, dropped[-1] + 0.5, color="tab:grey",
                       alpha=0.18, label="not production; excluded")
        if key in {"ntr", "hfns"}:
            ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(title)
        ax.set_xlabel("orbit index")
        ax.set_ylabel(units)
        ax.grid(alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(
        f"Baseline spin-up convergence: {'PASS' if all(criteria.values()) else 'NOT YET'} "
        f"({len(files)} orbits, {w}-orbit window ending at {window_end})"
    )
    plot_path = args.output / f"{run_dir.name}_convergence.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    # One writer for this key. There were two, and this one silently clobbered
    # the other, losing the report path and the failure list every time.
    payload = {"metrics": metrics, "criteria": criteria,
               "pass": all(criteria.values()), "failed_criteria": failed,
               "report": str(report_path.resolve()), "orbits": len(files),
               "window_orbits": w, "window_start_year_index": window_start,
               "window_end_year_index": window_end,
               "orbits_after_window_excluded_as_non_production": dropped}

    # Record the verdict with the run as well as in the analysis directory.
    # Provenance travels with the artifact everywhere else in this project, and a
    # convergence result that lives only in an output folder cannot be found from
    # the run it describes.
    manifest_path = args.run_dir / "run_manifest.json"
    if args.through is None and manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["convergence_assessment"] = payload
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
