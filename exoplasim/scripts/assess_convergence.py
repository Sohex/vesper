#!/usr/bin/env python3
"""Assess ExoPlaSim spin-up convergence from validated annual NetCDF files.

The verdict is taken over the last `--window` PRODUCTION orbits. Every segment
of a run declares what it was for, and orbits run to measure the model rather
than to advance the planet -- an I/O verification, a high-cadence wind sample, a
block on a differently patched binary -- are not evidence about where the run is
settling. `exoplasim/scripts/segments.py` owns that vocabulary.

`--assess diagnostic` takes the same six criteria over the orbits a run
declares ARE diagnostics. That is how an A/B arm shows it had settled:
sequencing.md A3 requires an arm's segments carry that label so no tail
reaches a climatology, and the label is what hides them from the default
mode. The verdict it produces is about the experiment and not the planet, so
it carries its mode in its filename and writes nothing into the run.

EVERY CRITERION IS AN UPPER BOUND: `|statistic| + its standard error <
threshold`. The argument for that form, and why it was still honest to adopt it
after runs had been judged under the other one, is recorded in `main` beside the
thresholds. Each criterion and the `resolving_power` row that prices it are
built from one table there, so a row can never describe an estimator the verdict
did not take.
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
# MEASURED 2026-08-26 on run_432e5e46adef, the bootstrap of canonical-10m-base,
# replacing an assumed pair that had drifted a long way from this model. The
# scatter is the spread over the settled window; tau is the integrated
# autocorrelation time of the per-orbit area-weighted mean surface temperature,
# by Geyer's initial monotone positive sequence over five candidate windows,
# which returned 1.89, 2.00, 2.07, 2.11 and 2.22 and were reliable on every one.
# The largest is carried, because `stationary_enough` refused all five -- the run
# still drifts about 0.003 to 0.007 K per orbit -- and a residual trend pushes
# every lag correlation UP, so each of those is an upper bound on the truth.
#
# WHAT MOVED, and it is not a correction of the old measurement. The 4.2 was
# (1+r)/(1-r) at a lag-1 of 0.615 measured on a different model: this tree has
# since taken the Stephens cloud tables, a derived orographic roughness and a
# hyperdiffusion 1.699x shorter, and the last of those damps the model faster,
# which is the mechanism a shorter memory would come from. This run reads a
# lag-1 of 0.39 to 0.49.
#
# WHY IT MATTERS HERE: the window this file derives goes as tau, and so does
# `lib/run_lengths.py`'s production span. At 4.2 the span was 84 orbits and at
# 10.43 it was 209, which is a ladder nobody can afford. At 2.2 it is 44.
NOMINAL_ORBIT_SCATTER_K = 0.091    # settled-window spread, run_432e5e46adef
NOMINAL_TAU_ORBITS = 2.2           # measured, and an upper bound; see above
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
def assessment_window(run_dir: Path, n_orbits: int, window: int,
                      purpose: str) -> tuple[int, int]:
    """The last `window` orbits of the declared purpose, as inclusive indices.

    HOW AN A/B ARM DEMONSTRATES ITS OWN SETTLING. `production_window` is the
    default and the only one a verdict about the PLANET may be taken on: it
    drops a diagnostic tail and refuses a window with a diagnostic hole in it.
    But sequencing.md A3 requires that the segments of a short A/B be labelled
    DIAGNOSTICS, and labelling them is then exactly what hides them from the
    test that would judge whether either arm had settled -- world-u9hq recorded
    A3's fourth condition as untestable for that reason, and no arm could
    satisfy both conditions as they stood.

    It is a mode here and not a widening of `production_window` because the two
    answer different questions. `production` asks where the planet is settling
    and its answer is the run's own verdict, written into the run's manifest.
    `diagnostic` asks whether an ARM had settled far enough for its difference
    against the other arm to mean anything, and that answer is about the
    experiment rather than the planet: `main` keeps it under its own filename
    and never writes it into the run, on exactly the terms `--through` is kept
    out.

    The trailing-drop and interior-hole rules are `production_window`'s, applied
    to the declared purpose instead: a diagnostic block interrupted by a spin-up
    segment is no more a trend than the other way round.
    """
    if purpose == "production":
        return production_window(run_dir, n_orbits, window)
    purposes = orbit_purposes(run_dir, range(n_orbits))
    admissible = {orbit for orbit, seen in purposes.items() if seen == purpose}
    end = n_orbits - 1
    while end >= 0 and end not in admissible:
        end -= 1
    if end < 0:
        raise RuntimeError(
            f"{run_dir.name} has no orbit any segment declares "
            f"`{purpose}`. A run whose segments do not say what they were for "
            f"cannot be assessed as though they had: label the segments, or "
            f"assess it as production.")
    start = end - window + 1
    if start < 0:
        raise RuntimeError(
            f"{run_dir.name} has {end + 1} orbits up to the last `{purpose}` "
            f"one (index {end}), fewer than the {window}-orbit window "
            f"requested")
    interior = sorted(o for o in range(start, end + 1) if o not in admissible)
    if interior:
        raise RuntimeError(
            f"orbits {interior} are inside the {window}-orbit window ending at "
            f"{end} and are not declared `{purpose}`. A window with a hole in "
            f"it is not a trend; shorten --window, or assess the block before "
            f"them.")
    return start, end


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

    THE FITTED TAU CARRIES ITS OWN ERROR TOO, and it is returned rather than
    dropped. It is not an input to any criterion; it is the only measurement of
    the relaxation time this project has, and the criterion's cheapness rests on
    the DERIVED time being a ceiling on it. A point estimate cannot test a
    ceiling -- one fitted value above a derived one is either a falsification or
    a fit's own noise, and without the error there is no telling which.

    Returns (asymptote, half_width, tau_orbits, tau_standard_error) or four nans
    if the series will not support a fit.
    """
    mask = orbits >= orbits.max() * fraction
    x, y = orbits[mask], series[mask]
    if x.size < 8:
        return float("nan"), float("nan"), float("nan"), float("nan")
    model = lambda n, inf, amp, tau: inf - amp * np.exp(-n / tau)
    try:
        popt, pcov = curve_fit(model, x, y,
                               p0=[y[-1], max(y[0] - y[-1], 1e-3), max(x.size / 3, 1.0)],
                               maxfev=40000)
    except Exception:
        return float("nan"), float("nan"), float("nan"), float("nan")
    inf, tau = float(popt[0]), float(popt[2])
    err = float(np.sqrt(np.diag(pcov))[0]) if np.all(np.isfinite(pcov)) else float("nan")
    tau_err = float(np.sqrt(np.diag(pcov))[2]) if np.all(np.isfinite(pcov)) else float("nan")
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
    # The same correction on the same covariance, for the same reason.
    if np.isfinite(tau_err):
        tau_err *= float(np.sqrt(tau_residual))
    # A fit extrapolating beyond its own span is reported as such rather than
    # trusted: widen the interval by how far past the data the asymptote sits.
    if np.isfinite(tau) and tau > 0:
        overshoot = max(0.0, tau / max(x.size, 1))
        err = (err if np.isfinite(err) else abs(inf - y[-1])) * (1.0 + overshoot)
    # The span the fit actually saw, so a consumer can tell a tau shorter than
    # its own data from one extrapolated past the end of it.
    return inf, err, tau, tau_err


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW_ORBITS,
                        help="orbits in the test window. Counted back from the "
                             "last orbit of the ASSESSED PURPOSE, not the last "
                             "orbit of the run. The "
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
    parser.add_argument("--assess", choices=("production", "diagnostic"),
                        default="production",
                        help="which orbits the verdict is taken over. "
                             "`production` is the run's own verdict, over "
                             "orbits its segments declare are the planet's "
                             "trajectory, and it is what gets written into the "
                             "run's manifest. `diagnostic` is how an A/B arm "
                             "demonstrates its own settling: A3 requires an "
                             "arm's segments be labelled diagnostics so no "
                             "tail reaches a climatology, and that label is "
                             "what hides them from the default mode. A "
                             "diagnostic verdict is about the EXPERIMENT and "
                             "not the planet -- it carries its mode in its "
                             "filename and writes nothing into the run.")
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
    #
    # `--assess diagnostic` takes the same window over the orbits an arm
    # DECLARES are diagnostics, which is the only way an arm labelled the way
    # A3 requires can show it had settled. See `assessment_window`.
    window_start, window_end = assessment_window(
        run_dir, len(files), args.window, args.assess)
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

    # THE STORAGE ESTIMATOR'S OWN ERROR, in the units it reports.
    #
    # `storage_w_m2_least_squares` is the least-squares slope of the planetary
    # heat content across the window, divided by the orbit. So its standard
    # error is the standard error of THAT slope in the same units, with the same
    # residual memory correction every other slope here carries -- not a
    # temperature-slope proxy converted through a sensitivity. That is what puts
    # the threshold and the estimator in one unit, which is the whole of the
    # convention this row applies: work out what the effect is worth in the
    # units the instrument reports and compare it with the instrument's own
    # scatter.
    #
    # The series is the one `close_state_energy` fitted, returned rather than
    # re-summed, so there is one heat content and not two.
    storage_se = float("nan")
    storage_scatter = float("nan")
    storage_memory = {"tau": float("nan"), "reliable": False, "lag1": float("nan")}
    if storage is not None:
        total_j = np.asarray(storage["heat_content_total_j_m2"], dtype=float)
        orbit_seconds = float(storage["orbit_seconds"])
        n_storage = total_j.size
        storage_se = slope_standard_error(total_j, n_storage) / orbit_seconds
        x_storage = np.arange(n_storage, dtype=float)
        storage_residual = total_j - np.polyval(
            np.polyfit(x_storage, total_j, 1), x_storage)
        storage_scatter = (float(np.std(storage_residual, ddof=1))
                           if n_storage > 2 else float("nan"))
        if n_storage > 3:
            storage_memory = ac.integrated_time(storage_residual)

    if storage is not None:
        metrics.update({
            "state_storage_w_m2": storage["storage_w_m2_least_squares"],
            "state_storage_endpoint_w_m2": storage["storage_w_m2_endpoint"],
            "state_storage_standard_error_w_m2": storage_se,
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
    asymptote, half_width, tau_fit, tau_fit_se = approach_to_equilibrium(
        orbits_axis, arrays["ts"])
    # How many orbits the fit was taken over. `approach_to_equilibrium` keeps the
    # last 65 per cent of the series, and a fitted tau longer than the span it
    # was fitted over is extrapolation rather than a measurement, which is what
    # decides whether it is evidence about the relaxation time at all.
    fit_span_orbits = int(np.count_nonzero(
        orbits_axis >= orbits_axis.max() * 0.35))
    if not np.isfinite(tau_fit):
        relaxation_fit_identifiable = False
        relaxation_fit_verdict = "the series would not support a fit at all"
    elif tau_fit <= 0.0:
        relaxation_fit_identifiable = False
        relaxation_fit_verdict = (
            f"the fit returned tau = {tau_fit:.6g} orbits; a negative time "
            "constant is the exponential having collapsed to a line")
    elif tau_fit > fit_span_orbits:
        relaxation_fit_identifiable = False
        relaxation_fit_verdict = (
            f"the fit returned tau = {tau_fit:.6g} orbits over {fit_span_orbits} "
            "orbits of data, so the exponential is indistinguishable from a "
            "line across the fitted range and tau is not determined by it")
    else:
        relaxation_fit_identifiable = True
        relaxation_fit_verdict = (
            f"tau = {tau_fit:.6g} orbits is shorter than the {fit_span_orbits} "
            "orbits fitted, so the series carries the curvature that fixes it")
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
    #
    # THE BRANCH IS TAKEN ONCE, HERE, AND EVERYTHING THAT DESCRIBES THE OFFSET
    # CRITERION COMES OUT OF IT. The verdict tests one estimator and the
    # `resolving_power` row prices an estimator, and those were two independent
    # pieces of code: the row multiplied the temperature slope's error by
    # `tau_expected` unconditionally, which is the FALLBACK statistic's error,
    # so on every run where the fit was used the row priced a statistic the
    # verdict had not taken. A reader could not tell which case a row was, and
    # `window_orbits_for_offset_criterion` was a fallback-only number without
    # saying so. Both now come out of this branch, so a change to one case
    # cannot leave the other pricing the wrong estimator. world-sl0q.
    fit_usable = (np.isfinite(asymptote) and np.isfinite(half_width)
                  and abs(offset) < 20.0 and half_width < 5.0)
    slope_se = slope_standard_error(arrays["ts"], w)
    if fit_usable:
        offset_statistic_source = "exponential_fit"
        # THE STATISTIC IS `|offset| + half_width` AND ITS ERROR IS THE FITTED
        # ASYMPTOTE'S OWN, because `half_width` IS that error:
        # `approach_to_equilibrium` returns the asymptote's standard error
        # corrected for the residuals' memory and widened where the fit
        # extrapolates past its data. The window mean the offset is measured
        # from carries an error too, and it is left out because the asymptote
        # and that mean are taken from one series and move together, so the
        # difference's error is smaller than adding the two in quadrature
        # would give.
        offset_statistic_se = half_width
        offset_statistic_error_source = (
            "the fitted asymptote's autocorrelation-corrected standard error, "
            "widened where the fit extrapolates past its own data")
    else:
        offset_statistic_source = "drift_fallback"
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
        # THE ERROR OF THE FALLBACK STATISTIC IS NOT ITS HALF WIDTH. The half
        # width carries a 100 per cent allowance on `tau_expected`, which is a
        # stated allowance and not a sampling error; what a resolving row asks
        # is whether the WINDOW can see the threshold. Both terms of
        # `|offset| + half_width` move with the same fitted slope, and each
        # moves by `tau_expected` times that slope's error, so the statistic's
        # sampling error is at most twice it. That is the conservative end and
        # it is the same product the default window was sized on.
        offset_statistic_se = (_OFFSET_STATISTIC_TERMS * tau_expected
                               * slope_se)
        offset_statistic_error_source = (
            "the temperature slope's standard error times the expected "
            "relaxation time, counted once for the offset and once for the "
            "half width that carries the same product")
    metrics.update({
        "temperature_asymptote_k": asymptote,
        "temperature_asymptote_half_width_k": half_width,
        "temperature_remaining_offset_k": offset,
        # A FITTED TAU IS REPORTED ONLY WHERE THE SERIES DETERMINES IT.
        #
        # Over a span short compared with tau, `exp(-n/tau)` is linear in n to
        # within the fit's own noise, so the exponential and a straight line are
        # the same curve and tau is whatever the optimiser drifted to. One run
        # here returned 44018 orbits against a derived 10.10 and another
        # -337429; those are the fit finding no curvature, not relaxation times,
        # and the gate that read them failed closed on the half width exactly as
        # it should while the artifact went on carrying a plausible-looking
        # field beside the verdict. A later consumer reading it at face value
        # would conclude the derived time is wrong by four orders of magnitude.
        #
        # So the rule, and it is a property of the data rather than of any run:
        # tau is IDENTIFIABLE from a series only when it is finite, positive and
        # no longer than the span the fit saw. Otherwise the honest output is a
        # refusal plus the bound the data supports, which is the span itself:
        # the series says tau is at least of that order and says nothing more.
        #
        # NOTHING NUMERIC CHANGES. `approach_to_equilibrium` still returns the
        # raw tau and still widens the asymptote's interval by the overshoot it
        # implies, which is how a fit extrapolating past its own data refuses
        # through the criterion. This governs the REPORTED field only.
        "relaxation_orbits_fitted": (
            tau_fit if relaxation_fit_identifiable else None),
        "relaxation_fit_identifiable": relaxation_fit_identifiable,
        "relaxation_fit_verdict": relaxation_fit_verdict,
        # The number the optimiser returned, named so nobody reads it as a
        # relaxation time. Kept because a degenerate fit is itself evidence
        # about the series, and dropping it would hide that the fit ran.
        "relaxation_fit_raw_tau_orbits": tau_fit,
        "relaxation_orbits_fitted_lower_bound_orbits": (
            None if relaxation_fit_identifiable else float(fit_span_orbits)),
        # The fit's own error on tau, autocorrelation-corrected like the
        # asymptote's. Reported and thresholded by nothing: it is what makes the
        # derived time's CEILING claim testable rather than a point comparison.
        "relaxation_orbits_fitted_standard_error": (
            tau_fit_se if relaxation_fit_identifiable else None),
        "relaxation_fit_span_orbits": fit_span_orbits,
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

    # THE FORM OF A CRITERION: AN UPPER BOUND, NOT A POINT ESTIMATE.
    #
    # Every criterion here is `|statistic| + its standard error < threshold`.
    # Five of the six tested `|statistic| < threshold` instead, which is a
    # point estimate: a state sitting exactly on its threshold passes about
    # half the time, on the statistic's own scatter, and that is a coin flip
    # rather than a criterion. A convergence verdict is what licenses the claim
    # that a run is equilibrated, and every number this project publishes about
    # the modelled world stands on that claim, so the standing conventions all
    # point one way -- an estimate that cannot be verified is bracketed and the
    # bracket is reported, a convergence claim is labelled honestly when it
    # misses, and a red gate naming a real conflict beats a green one hiding
    # it. The upper-bound form is that disposition written into the test.
    #
    # THE THRESHOLDS DID NOT MOVE. What changed is the form of the comparison,
    # not its tolerance: a threshold that now looks wrong under this form is a
    # separate finding and gets its own row rather than an adjustment here.
    #
    # WHY IT WAS STILL HONEST TO CHANGE THE FORM. Adopting it flips
    # `run_ec32946bec89` from pass to fail on the storage criterion, and a
    # criterion chosen after the run it judges is not a criterion. What
    # dissolves the objection is that NOTHING CANONICAL HAS BEEN RUN: the
    # canonical climatology lineage does not exist yet, so every build, run and
    # climatology in this tree is disposable whatever has consumed it, and that
    # run is a stepping stone on a binary already superseded. The flip costs
    # nothing that is being kept, which is the window in which a criterion's
    # FORM can still be settled. After the lineage is declared this change
    # becomes impossible to make honestly, which is why it is recorded here
    # beside the criteria rather than in a commit message. world-s8n3.
    CRITERION_FORM = ("|statistic| + its standard error < threshold. A "
                      "criterion states an upper bound, so a statistic whose "
                      "interval reaches the threshold fails.")

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
    #   What is measurable. This half was stated as "the storage estimate on a
    #   10-orbit block differs from the same run's 20-orbit block by 0.111, 0.115
    #   and 0.116 W/m2, so a threshold below about 0.11 is measuring the
    #   estimator". That is not the estimator's standard error: the two blocks
    #   are nested, so their difference is neither an independent pair nor a
    #   scatter, and it was taken at a 10-orbit window that is no longer in use.
    #   The estimator's own standard error is the standard error of the slope it
    #   fits, in the units it reports, and it is now computed per run and per
    #   window as `resolving_power` row `abs_state_storage`. Measured 2026-08-25
    #   on the 35-orbit window of this project's dt-45 run: the heat content's
    #   residual scatter is 1.09e7 J/m2 with a lag-1 of 0.80, giving 0.030 W/m2
    #   at 35 orbits, 0.042 at 28 and 0.197 at 10. So the measurable floor is
    #   3 * 0.030 = 0.090 W/m2 at the window now derived, and was 0.59 W/m2 at
    #   the 10-orbit window this threshold was first written against.
    #
    # THE THRESHOLD IS UNCHANGED AT 0.12 AND THE TEMPERATURE BOUND IS WHAT HOLDS
    # IT. 0.090 is below 0.118, so the measurable half no longer binds and the
    # number rests on the offset tolerance alone, which is where it was derived
    # from and which was fixed before any of this was read. What the correction
    # changes is not the number but what a reader may conclude from it: at a
    # 10-orbit window this criterion could not see its own threshold, and the
    # `resolving_power` row is now what says so per run rather than a sentence
    # here saying it once.
    STORAGE_TOLERANCE_W_M2 = 0.12

    # ONE TABLE, TWO CONSUMERS: the verdict below and the `resolving_power` row
    # beside it. They were built by two separate pieces of code over two
    # separate lists, and a row that prices a different estimator than the
    # verdict was taken on is worse than no row -- a resolving row exists to say
    # what the window can SEE about the statistic being tested. Every criterion
    # is one entry here, so a criterion cannot gain a verdict without a row, or
    # a row without a verdict, and neither can be repriced without the other
    # moving with it. world-s8n3, world-sl0q.
    #
    # `interval` is what the upper-bound form adds to `|statistic|`, and
    # `standard_error` is what the resolving row prices. They are the same
    # number for the five direct measurements. They differ for the offset
    # criterion alone, and only in its fallback branch, where the reported half
    # width carries a stated 100 per cent allowance on `tau_expected` on top of
    # the sampling error: the verdict must carry that allowance and the
    # resolving question must not, because the window's ability to see 0.15 K
    # is not a property of how far the run still has to travel.
    criterion_table = [
        {"key": "abs_temperature_slope_lt_0.05_k_per_orbit",
         "statistic": metrics["temperature_slope_k_per_orbit"],
         "interval": metrics["temperature_slope_standard_error_k_per_orbit"],
         "standard_error": metrics["temperature_slope_standard_error_k_per_orbit"],
         "threshold": 0.05,
         "error_source": "the fitted slope's own standard error, with the "
                         "residual memory in it"},
        {"key": "abs_toa_slope_lt_0.05_w_m2_per_orbit",
         "statistic": metrics["toa_balance_slope_w_m2_per_orbit"],
         "interval": metrics["toa_balance_slope_standard_error_w_m2_per_orbit"],
         "standard_error": metrics["toa_balance_slope_standard_error_w_m2_per_orbit"],
         "threshold": 0.05,
         "error_source": "the fitted slope's own standard error, with the "
                         "residual memory in it"},
        {"key": "abs_surface_slope_lt_0.05_w_m2_per_orbit",
         "statistic": metrics["surface_balance_slope_w_m2_per_orbit"],
         "interval": metrics["surface_balance_slope_standard_error_w_m2_per_orbit"],
         "standard_error": metrics["surface_balance_slope_standard_error_w_m2_per_orbit"],
         "threshold": 0.05,
         "error_source": "the fitted slope's own standard error, with the "
                         "residual memory in it"},
        {"key": "abs_sea_ice_slope_lt_0.001_per_orbit",
         "statistic": metrics["sea_ice_slope_fraction_per_orbit"],
         "interval": metrics["sea_ice_slope_standard_error_fraction_per_orbit"],
         "standard_error": metrics["sea_ice_slope_standard_error_fraction_per_orbit"],
         "threshold": 0.001,
         "error_source": "the fitted slope's own standard error, with the "
                         "residual memory in it"},
        # FAILS CLOSED where the window could not be closed against the state:
        # `storage` is None, the statistic is nan, and a nan fails the
        # comparison below. A convergence check that cannot measure the thing
        # it tests must refuse rather than abstain.
        {"key": f"abs_state_storage_lt_{STORAGE_TOLERANCE_W_M2}_w_m2",
         "statistic": (storage["storage_w_m2_least_squares"]
                       if storage is not None else float("nan")),
         "interval": storage_se, "standard_error": storage_se,
         "threshold": STORAGE_TOLERANCE_W_M2,
         "error_source": "the standard error of the heat content's fitted "
                         "slope, in the units the estimator reports"},
        # The one that bounds the answer rather than its rate of change, and
        # the one that already had this form. Its interval and its error are
        # branch-dependent and both were settled where the branch was taken.
        {"key": f"extrapolated_offset_lt_{OFFSET_TOLERANCE_K}_k",
         "statistic": offset, "interval": half_width,
         "standard_error": offset_statistic_se,
         "threshold": OFFSET_TOLERANCE_K,
         "error_source": offset_statistic_error_source,
         "statistic_source": offset_statistic_source},
    ]

    def _number(value) -> float:
        """A criterion term as a float, with a missing one as nan so it fails."""
        return float("nan") if value is None else float(value)

    criteria = {}
    for entry in criterion_table:
        # THE THRESHOLD IN THE NAME IS THE THRESHOLD APPLIED, checked rather
        # than trusted. The key is what every consumer, every manifest and
        # every note quotes a criterion by, and it is written out in full above
        # so that grepping for one finds it; a threshold moved in the code and
        # left in the name would be a silent re-verdict on every run assessed
        # after it, and this is the one place that can catch it.
        entry["name"] = entry["key"].rsplit("_lt_", 1)[0]
        named = float(entry["key"].split("_lt_", 1)[1].split("_")[0])
        if named != entry["threshold"]:
            raise RuntimeError(
                f"criterion {entry['key']} is named for {named} and applies "
                f"{entry['threshold']}")
        statistic, interval = _number(entry["statistic"]), _number(entry["interval"])
        criteria[entry["key"]] = bool(
            np.isfinite(statistic) and np.isfinite(interval)
            and abs(statistic) + interval < entry["threshold"])

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
        # The storage estimator's own scatter, in the units it reports. The
        # temperature block above prices the four slope criteria and the offset;
        # this one prices the storage criterion, and it is a different series
        # with a different memory. Reporting the temperature figures alone was
        # what left the sixth criterion without a row at all.
        "storage_residual_scatter_j_m2": storage_scatter,
        "storage_residual_tau_orbits": storage_memory["tau"],
        "storage_residual_lag1": storage_memory["lag1"],
        "storage_tau_span_supports_the_estimate": bool(
            storage_memory.get("reliable", False)),
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
    # EVERY CRITERION HAS A ROW AND EVERY ROW PRICES ITS OWN CRITERION'S
    # STATISTIC, by construction rather than by a check afterwards: both come
    # from `criterion_table` above, in one pass. The storage criterion once had
    # no row at all and the offset row once priced an estimator its verdict had
    # not used, and both were possible because the rows were a second list.
    # Each row carries WHICH statistic its error belongs to, so a reader can
    # tell one branch of the offset criterion from the other without inferring
    # it from the metrics.
    for entry in criterion_table:
        statistic_se = _number(entry["standard_error"])
        resolves = bool(np.isfinite(statistic_se)
                        and statistic_se * RESOLVING_FACTOR <= entry["threshold"])
        row = {
            "criterion": entry["name"], "threshold": entry["threshold"],
            "statistic_standard_error": statistic_se,
            "statistic_error_source": entry["error_source"],
            # What the verdict added to `|statistic|`. The same number as the
            # standard error everywhere except the offset criterion's fallback,
            # where the verdict also carries the allowance on `tau_expected`.
            "verdict_interval": _number(entry["interval"]),
            "resolves": resolves}
        if "statistic_source" in entry:
            row["statistic_source"] = entry["statistic_source"]
        resolving["rows"].append(row)
    # WHICH ESTIMATOR THE OFFSET CRITERION WAS ON, beside the rows rather than
    # only inside one of them, because the window figures below are the
    # fallback form's arithmetic and a reader has to be able to tell whether
    # they price this run's verdict at all.
    resolving["offset_statistic_source"] = offset_statistic_source

    # THE WINDOW THE STORAGE CRITERION NEEDS TO RESOLVE ITS OWN THRESHOLD,
    # in orbits, from the heat content's own scatter and memory. The target is
    # the threshold divided by RESOLVING_FACTOR, converted from W/m2 into the
    # J/m2 per orbit the fitted slope is in, so nothing here passes through a
    # temperature sensitivity.
    resolving["window_orbits_for_storage_criterion"] = (
        orbits_for_slope_standard_error(
            storage_scatter, storage_memory["tau"],
            STORAGE_TOLERANCE_W_M2 / RESOLVING_FACTOR * float(
                storage["orbit_seconds"]))
        if storage is not None else float("nan"))
    resolving["window_orbits_for_storage_criterion_if_independent"] = (
        orbits_for_slope_standard_error(
            storage_scatter, 1.0,
            STORAGE_TOLERANCE_W_M2 / RESOLVING_FACTOR * float(
                storage["orbit_seconds"]))
        if storage is not None else float("nan"))
    # THE WINDOW THE OFFSET CRITERION NEEDS, IN ITS FALLBACK FORM. It inverts
    # `var(slope) = sigma^2 * tau * 12 / (n (n^2 - 1))` for the slope error the
    # fallback multiplies by `tau_expected`, so it prices this run's verdict
    # only where the fallback was taken. Where the fit was used the criterion's
    # error is the fit's own and no window formula inverts it; the figure is
    # still reported, as the window this run's variability would need if it
    # were on the fallback, and `resolving_power.offset_statistic_source` says
    # which case the run is. world-sl0q.
    resolving["window_orbits_for_offset_criterion"] = window_for_offset_criterion(
        ts_scatter, ts_tau, tau_expected)
    # The same window if the orbits were independent, so the price of the
    # memory is visible rather than folded into one number.
    resolving["window_orbits_for_offset_criterion_if_independent"] = \
        window_for_offset_criterion(ts_scatter, 1.0, tau_expected)
    resolving["window_orbits_for_offset_criterion_prices"] = (
        "the drift fallback form of the criterion; it prices this run's own "
        "verdict only where offset_statistic_source is drift_fallback")
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
        # WHICH ORBITS THIS VERDICT IS ABOUT, by their declared purpose. A
        # diagnostic verdict is about an A/B arm's settling and not about where
        # the planet is going, and a reader must be able to tell the two apart
        # from the report alone rather than from its filename.
        "assessed_purpose": args.assess,
        "verdict_is_about": (
            "the planet's trajectory; this is the run's own verdict"
            if args.assess == "production" else
            "the EXPERIMENT: whether the orbits an A/B arm declares as "
            "diagnostics had settled far enough for a difference against "
            "another arm to mean anything. It is not the run's verdict and is "
            "not written into the run's manifest."),
        "segment_purposes": {str(k): v for k, v in sorted(purposes.items())},
        "metrics": metrics,
        "criteria": criteria,
        "criteria_provenance": {
            "form": CRITERION_FORM,
            "form_applies_to": "all six criteria, each with its own "
                               "statistic's standard error. Five tested a "
                               "point estimate until world-s8n3; the "
                               "thresholds did not move with the form.",
            "energy_balance_tests": "state storage, not reported top-of-atmosphere net",
            "storage_tolerance_w_m2": STORAGE_TOLERANCE_W_M2,
            "storage_tolerance_derivation":
                "0.15 K offset tolerance / (0.128 K per orbit per W/m2 * 9.9 orbits) "
                "= 0.118 W/m2, rounded to 0.12. The temperature tolerance is the "
                "whole of the derivation: the measurability half it was first "
                "paired with compared two NESTED blocks rather than taking the "
                "estimator's standard error, and what that error actually is at "
                "this run's own window is the resolving_power row "
                "abs_state_storage.",
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
    #
    # A DIAGNOSTIC ASSESSMENT CARRIES ITS MODE IN THE NAME for the same reason.
    # It is a verdict about an arm's own settling, taken over orbits the run
    # declares are not its trajectory, and it must never stand where a reader
    # looks for the run's verdict.
    suffix = "" if args.through is None else f"_through{args.through:03d}"
    if args.assess != "production":
        suffix += f"_{args.assess}"
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
    #
    # A DIAGNOSTIC ASSESSMENT NEVER TOUCHES IT EITHER, and for a stronger
    # reason than the truncated one. `status` and `equilibrium_cutoff_year_
    # index` are claims about where the PLANET has got to, and a diagnostic
    # segment is by declaration not evidence about that. Writing an arm's own
    # settling into those fields would relabel a run as equilibrated on orbits
    # that were excluded from every verdict about it -- which is the failure
    # the purpose vocabulary exists to prevent, arriving through the tool that
    # reads it.
    manifest_path = run_dir / "run_manifest.json"
    if args.through is not None or args.assess != "production":
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
                       alpha=0.18,
                       label=f"not {args.assess}; excluded")
        if key in {"ntr", "hfns"}:
            ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(title)
        ax.set_xlabel("orbit index")
        ax.set_ylabel(units)
        ax.grid(alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(
        f"{args.assess.capitalize()} convergence: "
        f"{'PASS' if all(criteria.values()) else 'NOT YET'} "
        f"({len(files)} orbits, {w}-orbit window ending at {window_end})"
    )
    # The plot takes the report's suffix, so a diagnostic verdict cannot
    # overwrite the figure a reader looks at for the run's own.
    plot_path = args.output / f"{run_dir.name}_convergence{suffix}.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    # One writer for this key. There were two, and this one silently clobbered
    # the other, losing the report path and the failure list every time.
    payload = {"metrics": metrics, "criteria": criteria,
               "pass": all(criteria.values()), "failed_criteria": failed,
               "report": str(report_path.resolve()), "orbits": len(files),
               "window_orbits": w, "window_start_year_index": window_start,
               "window_end_year_index": window_end,
               "assessed_purpose": args.assess,
               "orbits_after_window_excluded_as_non_production": dropped}

    # Record the verdict with the run as well as in the analysis directory.
    # Provenance travels with the artifact everywhere else in this project, and a
    # convergence result that lives only in an output folder cannot be found from
    # the run it describes.
    #
    # THE RUN'S OWN VERDICT ONLY. `convergence_assessment` is the key every
    # consumer reads to decide whether a run may be used, so a diagnostic
    # assessment writing there would answer that question with a verdict about
    # an experiment. It stays in the analysis directory under its own name.
    manifest_path = args.run_dir / "run_manifest.json"
    if (args.through is None and args.assess == "production"
            and manifest_path.is_file()):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["convergence_assessment"] = payload
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
