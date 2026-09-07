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
import math
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
from segments import (io_regime_changes, non_production_orbits,
                      orbit_purposes, production_window,
                      refuse_a_window_spanning_an_io_regime_change)
# CLAUDE.md names lib/sensitivity.py as the one flux-to-kelvin conversion, and
# the radiative damping this file relaxes at is that conversion inverted.
from sensitivity import planetary_albedo_from_fluxes, radiative_damping_w_m2_per_k
# lib/run_lengths.py owns the two timescales a run is bought in. The memory time
# this file sizes its window on is the same quantity, so it is read from there
# rather than restated: a window derived from one number and a span derived from
# another statement of it is two ladders.
import run_lengths



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


def annual_records(files: list[Path]) -> list[dict]:
    """The global annual means every verdict in this file is taken over.

    ONE READER, because the cold-start assessment and the reconvergence
    assessment must be looking at the same numbers when both are reported for
    one run.
    """
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
    return records


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
# ALL THREE ARE UPPER BOUNDS AND THE CHECK BELOW HOLDS THEM TO IT. The window
# grows with each of them -- as `scatter^(2/3)`, as `tau^(1/3)` and as
# `relaxation^(1/3)` -- so a nominal BELOW what a run reads under-sizes the
# window that run is judged in, and the criterion then discriminates less than
# it claims. A nominal above only buys orbits. That asymmetry is the whole
# reason these are bounds rather than best estimates, and it is what makes
# "the artifact reads lower" a pass and "the artifact reads higher" a failure.
#
# `check_convergence_bounds` is what re-reads them. It is not a comparison of
# each number with the run it was taken on -- that comparison is the circle the
# audit named, since the run was assessed over the window these produced. It is
# a bound test against every report on disk, plus an ANCHOR test that fires when
# the report a bound was taken beside changes at all. The second is what tells
# "the bound still holds" from "the bound was never re-examined".

# THE SETTLED-WINDOW SPREAD, and it is a bound over EVERY report rather than a
# reading of one run. Two things have raised it, in order:
#
# 0.091 to 0.092, because the window it is read on became the clean-I/O block
# alone and the clean stream is the noisier instrument: PlaSim's low-I/O
# accumulation averages over the output interval, which suppresses per-orbit
# variance, so a window reaching back across the join read 0.0910 where the
# twelve clean orbits read 0.09103.
#
# 0.092 to 0.106, on the intersection bracket's COLD arm, which reads 0.10577
# over thirty-seven clean orbits. That arm is bare rock, three kelvin colder,
# and carries a mean sea ice fraction of 0.108 against the vegetated arm's
# 0.080 -- and sea ice is where this model's interannual variability lives, so
# a colder arm is a noisier one. Both arms are settled: the cold one passes all
# six criteria and RESOLVES all six, which the vegetated bootstrap does not.
#
# 0.106 to 0.117, when that carved baseline did converge. Its 57-orbit clean
# production window reads 0.116626 K and passes all six criteria. The carve did
# drain the world -- lakes from 4.25 per cent of the planet to 1.31, the
# endorheic share of land from 74 to 51 -- so the larger variability that first
# appeared during approach remains in the settled climate and now counts.
#
# 0.117 to 0.124, on run_0d41aa82c287, the T21 arm at dt 45 cut as the donor for
# the T42 ladder comparison. It reads 0.123493 K over a 52-orbit window and
# passes all six criteria on the same build and the same rung as the 0.117
# reading, so the difference is not the world: it is the WINDOW. That reading
# was taken over 57 orbits and this one over 108, and a longer settled window
# admits more of the low-frequency variability the short one could not see. The
# same effect is what
# `exoplasim/notes/memory-time-and-the-production-span.md` measures on the
# memory time, and it is why this is a bound over reports rather than a property
# of the process.
#
# 0.124 to 0.145, on run_9d5dbf9bd3d9, the WARM arm of the flux-to-kelvin
# bracket at flux ratio 1.000. It reads 0.144917 K over a 55-orbit clean window
# and passes all six criteria. Its cold twin, identical in structure and in
# everything but the flux, reads 0.103598 over the same window indices, so the
# 40 per cent difference is the flux and nothing else.
#
# AND IT IS THE LESS ICY ARM THAT IS NOISIER, which bounds the reading above
# that sea ice is where this model's interannual variability lives. Across that
# ALBEDO pair the colder, icier arm was the noisier one; across this FLUX pair
# the warm arm carries a mean sea-ice fraction of 0.0190 against the cold arm's
# 0.0753 and is half again as noisy. Whatever sets the scatter, ice cover does
# not order it on its own. `notes/audits/flux-slope-bracket.md`.
#
# THE FINER RUNG IS QUIETER, which is worth recording because it bounds how far
# this can run. The two T42 arms at the same step and build read 0.074925 and
# 0.079107 over the same window length, so the ladder's expensive rungs are not
# where this bound is set; T21 is.
#
# The bound therefore has to cover the noisiest world this project judges and
# not the most comfortable, per the asymmetry argued above: a nominal below
# what a run reads under-sizes the window that run is judged in, and a nominal
# above only buys orbits.
NOMINAL_ORBIT_SCATTER_K = 0.145

# THE MEMORY TIME STATES NO NUMBER HERE. `lib/run_lengths.py` declares the
# bracket, with the sweep it came from and the anchor that re-examines it, and
# the window is sized on the TOP of that bracket: the window has to be long
# enough for the longest memory the bracket admits, and sizing on the bottom
# would buy a window the model's own variability can outrun.
NOMINAL_TAU_ORBITS = max(run_lengths.TAU_MEMORY_ORBITS_BRACKET)

# THE RELAXATION TIME THE WINDOW IS SIZED ON, bounding what `relaxation_orbits`
# below returns on a real run. It is DERIVED per run, from the modelled mixed
# layer's heat capacity and that run's own radiative damping, so it moves with
# the run's albedo; the nominal exists only because the default window has to be
# a number before any run is opened.
#
# It was 10.0 and every recent report derived MORE than that, up to 11.75, so
# the default window was sized on a relaxation time no run had. The bound is
# restored here at the top of what the reports derive, to the precision that set
# supports, and `check_convergence_bounds` refuses if a report goes above it.
NOMINAL_RELAXATION_ORBITS = 11.8
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
# THE SEA-ICE SLOPE'S THRESHOLD, STATED ONCE. It was a literal inside the
# cold-start criterion table, and the reconvergence verdict below carries the
# same criterion: fixing a decay at the derived relaxation time is conservative
# only while the modelled sea ice has stopped reorganising, which is the one
# regime that time is not a ceiling in. Two criteria reading two copies of one
# threshold is two chances for one of them to move.
SEA_ICE_SLOPE_TOLERANCE_PER_ORBIT = 0.001


def window_for_offset_criterion(scatter: float, tau: float,
                                relaxation_orbits: float) -> float:
    """The window at which the offset criterion starts to discriminate."""
    return orbits_for_slope_standard_error(
        scatter, tau, _OFFSET_TOLERANCE_K
        / (_RESOLVING_FACTOR * _OFFSET_STATISTIC_TERMS * relaxation_orbits))


DEFAULT_WINDOW_ORBITS = int(np.ceil(window_for_offset_criterion(
    NOMINAL_ORBIT_SCATTER_K, NOMINAL_TAU_ORBITS, NOMINAL_RELAXATION_ORBITS)))


# ==========================================================================
# RECONVERGENCE: A DIFFERENT QUESTION, AND THEREFORE A DIFFERENT CRITERION.
# ==========================================================================
#
# Everything above asks whether a state is TRENDLESS. It is established from
# scratch, its equilibrium is unknown, and the only evidence about where the
# run is going is the run itself, so the answer has to be read out of a slope
# fitted to noise and multiplied by the relaxation time. That is why it costs
# `DEFAULT_WINDOW_ORBITS`.
#
# A RECONVERGENCE ASKS WHETHER A TRANSIENT HAS DECAYED, and there the record
# carries three things the cold-start test neither has nor uses:
#
#   1. THE TRANSIENT'S START IS STAMPED. The run began from a converted or
#      seeded restart, so orbit zero IS the conversion and the whole record is
#      the transient. A cold start does not know which part of its record is
#      approach.
#   2. THE DEPARTURE IS BOUNDED AND MEASURABLE WHERE THE SIGNAL IS LARGEST.
#      The state was AT a climate when the conversion moved it, so the record
#      spans the departure from its full size down to whatever is left, and the
#      amplitude is fitted off the steep part rather than inferred from the
#      tail.
#   3. THE RELAXATION TIME IS DERIVED FROM THE PLANET, not from the window.
#      `relaxation_orbits` is the slab's heat capacity over the run's own
#      radiative damping, both measured outside this record, so the shape of
#      the decay is known before the fit and only its amplitude and asymptote
#      are free. `lib/autocorrelation.py` records that a tau taken from inside
#      a short window comes back badly biased low; nothing here takes one.
#
# THE MODEL. T(n) = L + D exp(-n / tau_relax) over the post-conversion record,
# with tau_relax FIXED. Two linear parameters, so it is an ordinary least
# squares with a closed form -- it cannot fail to converge, which is the exact
# failure `approach_to_equilibrium` has to guard against on a settled run.
#
# THE STATISTIC IS A DISTANCE, THE SAME KIND THE OFFSET CRITERION BOUNDS:
# `remaining = D * exp(-(n-1)/tau_relax)`, how far the last orbit still sits
# from where the transient is going, compared as `|remaining| + half_width`
# against the SAME `_OFFSET_TOLERANCE_K`. Nothing is loosened; the question is
# different and so is the estimator.
#
# WHY FIXING tau_relax IS CONSERVATIVE AND NOT CONVENIENT. Fitting a decay at a
# tau larger than the truth reports MORE remaining departure than there is, so
# the direction that matters is the one where the derived tau is a CEILING on
# the true relaxation. For this model every large correction points the same
# way:
#
#   - The damping is the model's own EQUILIBRIUM response, measured between two
#     converged runs several kelvin apart (`lib/sensitivity.py`). During a
#     transient the modelled sea ice has not finished moving, so less of the
#     positive ice-albedo feedback is engaged and the damping acting on the
#     departure is stronger, which SHORTENS the relaxation.
#   - The heat capacity is a full ocean mixed-layer column applied to the whole
#     planet. Land carries a soil column of a few per cent of it, so the
#     area-weighted capacity is SMALLER, which shortens the relaxation again.
#   - The reservoirs that lengthen it -- sea ice and snow as latent heat, the
#     soil column, the atmosphere and its vapour, the ones
#     `close_state_energy.py` enumerates -- are together a small fraction of a
#     mixed-layer column, and only over the fraction of the planet that has
#     them.
#
# THE ONE REGIME WHERE THAT ARGUMENT FAILS is a state whose sea ice is
# reorganising: `lib/sensitivity.py` records a much larger response across the
# ice transition than away from it, and the relaxation time goes with it. So
# the reconvergence verdict CARRIES the sea-ice criterion rather than assuming
# it, and a run whose ice is still moving is refused instead of being assessed
# on a relaxation time its own regime does not have.
#
# WHAT IT COSTS. Two things have to be true and the record length is the larger
# of them: the modelled transient must have fallen under the tolerance, which
# takes `tau_relax * ln(D / tolerance)` orbits and is LOGARITHMIC in the
# departure; and the statistic must resolve its own threshold at
# `_RESOLVING_FACTOR`, which is `DEFAULT_RECONVERGENCE_ORBITS` and does not
# depend on the departure at all. That is the whole reason a reconvergence is
# cheaper than a cold start rather than a shorter version of it: the cold-start
# window is set by how long noise takes to average out, and this one is set by
# how long the transient takes to decay.
#
# THE LOGARITHM IS NOT RESTATED HERE. `lib/run_lengths.py:settling_orbits` is
# this project's one statement of `n = tau ln(A / residual)`, with the same
# 0.15 K residual and the argument for it, and a settling block is priced from
# it. Two statements of one relation is two chances for one of them to move.


def departure_decay(series, tau_relax: float, tau_ac: float) -> dict:
    """Fit T(n) = L + D exp(-n/tau_relax) and report what is left to travel.

    `tau_relax` is fixed, so this is a two-parameter linear least squares.
    `tau_ac` is the integrated autocorrelation time of consecutive orbits and
    must come from OUTSIDE this record: it inflates the residual variance the
    way `lib/autocorrelation.py` says a series with memory requires, and an
    estimate taken from inside a window this short is biased low, which is the
    direction that passes a run it should not.
    """
    y = np.asarray(series, dtype=float)
    n = y.size
    blank = {"level_k": float("nan"), "departure_amplitude_k": float("nan"),
             "remaining_departure_k": float("nan"),
             "half_width_k": float("nan"),
             "residual_scatter_k": float("nan"),
             "undecayed_fraction": float("nan"),
             "orbits": int(n)}
    if n < 4 or not np.isfinite(tau_relax) or tau_relax <= 0:
        return blank
    x = np.arange(n, dtype=float)
    u = np.exp(-x / tau_relax)
    design = np.column_stack([np.ones(n), u])
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - design @ beta
    variance = float(np.dot(residual, residual)) / (n - 2)
    s_uu = float(np.dot(u - u.mean(), u - u.mean()))
    if s_uu <= 0.0:
        return blank
    undecayed = float(np.exp(-(n - 1) / tau_relax))
    amplitude = float(beta[1])
    memory = ac.integrated_time(residual) if n > 3 else {
        "tau": float("nan"), "reliable": False}
    return {"level_k": float(beta[0]),
            "departure_amplitude_k": amplitude,
            "remaining_departure_k": amplitude * undecayed,
            "half_width_k": float(np.sqrt(variance * max(tau_ac, 1.0) / s_uu))
                            * undecayed,
            "residual_scatter_k": float(np.sqrt(variance)),
            "undecayed_fraction": undecayed,
            "orbits": int(n),
            # Reported, never used. It is the tau this record would give from
            # inside itself, and it is here so the bias lib/autocorrelation.py
            # measured can be seen rather than trusted.
            "residual_tau_from_this_record": memory["tau"],
            "residual_tau_span_supports_it": bool(memory.get("reliable", False))}


def orbits_for_departure_decay(scatter: float, tau_ac: float,
                               tau_relax: float, target: float) -> float:
    """Orbits at which the remaining departure's standard error reaches `target`.

    The counterpart of `orbits_for_slope_standard_error` for this criterion,
    and the thing that prices a reconvergence block. It has no closed form --
    the exponential regressor's spread and its decay both move with n -- so it
    is walked, which costs nothing at these lengths.
    """
    if target <= 0 or tau_relax <= 0 or scatter <= 0 \
            or not all(np.isfinite([scatter, tau_ac, tau_relax])):
        return float("nan")
    for n in range(4, 4001):
        x = np.arange(n, dtype=float)
        u = np.exp(-x / tau_relax)
        s_uu = float(np.dot(u - u.mean(), u - u.mean()))
        if s_uu <= 0.0:
            continue
        se = scatter * np.sqrt(max(tau_ac, 1.0) / s_uu) \
            * float(np.exp(-(n - 1) / tau_relax))
        if se <= target:
            return float(n)
    return float("nan")


def orbits_to_departure_tolerance(amplitude: float, tau_relax: float,
                                  tolerance: float) -> float:
    """When the modelled transient itself falls under the bar, as a RECORD LENGTH.

    The other half of what a reconvergence costs, and the half that carries the
    departure. It is what makes a timestep change cheaper than a resolution
    conversion without either needing its own criterion: the two differ in D,
    the fit MEASURES D, and the orbits it buys go as the logarithm of it.

    THE LOGARITHM IS `lib/run_lengths.py:settling_orbits` AND IS NOT RESTATED.
    What is added here is the offset between a decay TIME and a record LENGTH:
    the fit's last sample sits at index `n - 1`, so a record of `n` orbits has
    decayed for `n - 1` of them.
    """
    if not all(np.isfinite([amplitude, tau_relax, tolerance])) \
            or tau_relax <= 0 or tolerance <= 0:
        return float("nan")
    if abs(amplitude) <= tolerance:
        return 0.0
    return float(run_lengths.settling_orbits(abs(amplitude), tau_relax,
                                             tolerance) + 1.0)


# The block a reconvergence has to run before its criterion can decide
# anything, at the same nominal scatter, autocorrelation time and relaxation
# time the cold-start window is priced from, and against the same bar: the
# statistic's standard error at most `_RESOLVING_FACTOR` of its threshold.
# Derived here, so moving it means moving an input that has a citation.
DEFAULT_RECONVERGENCE_ORBITS = int(orbits_for_departure_decay(
    NOMINAL_ORBIT_SCATTER_K, NOMINAL_TAU_ORBITS, NOMINAL_RELAXATION_ORBITS,
    _OFFSET_TOLERANCE_K / _RESOLVING_FACTOR))


# HOW MUCH OF THE SERIES THE EXPONENTIAL FIT SEES. A DECISION and not a measured
# quantity: the head of a cold start is not on the exponential the fit is for,
# and dropping a third of it is where the fit stops being pulled by the part of
# the approach it does not model. Stated ONCE, here, because it was stated three
# times -- this default, a literal beside the span it implies, and a copy in
# `check_relaxation_ceiling.py` reconstructing the span from it -- and three
# statements of one algorithm parameter is three chances for two of them to
# agree while the third does not.
FIT_TAIL_FRACTION = 0.35


def fit_span_orbits(total_orbits: int) -> int:
    """How many orbits `approach_to_equilibrium` sees on a series this long.

    The ONE implementation of the fit's span. `check_relaxation_ceiling.py`
    reconstructs the span of artifacts written before it was recorded in them,
    and calls this rather than reimplementing the mask.
    """
    n = int(total_orbits)
    if n <= 0:
        return 0
    axis = np.arange(n, dtype=float)
    return int(np.count_nonzero(axis >= axis.max() * FIT_TAIL_FRACTION))



# THE FIT SPANS THE APPROACH BY CONSTRUCTION, SO IT CAN CROSS AN I/O JOIN.
# The verdict WINDOW refuses to span one -- `segments.py` argues why -- but the
# relaxation fit is a different quantity: it needs the approach, and on this
# project every approach was integrated in the low-I/O regime while the settled
# tail is clean. Refusing outright would make the relaxation time unavailable
# on every run that has one to measure.
#
# So the step is MEASURED against the criterion it would decide, which is
# CLAUDE.md's rule about checking the instrument against the size of the
# effect. The asymptote's whole use is the offset criterion, `asymptote -
# window mean` against `_OFFSET_TOLERANCE_K`. A step at the join enters the
# asymptote directly, so a step at or above that tolerance means the fit cannot
# decide the criterion whatever it returns, and a step well below it means the
# join is not what the answer rests on.
def io_step_at_join(series, run_dir: Path, first: int, last: int,
                    sample: int = 5) -> tuple[int | None, float]:
    """The jump in `series` across the first I/O-regime change in [first, last].

    Means over up to `sample` orbits each side, so the answer is a level change
    rather than one pair of orbits. Returns `(None, 0.0)` when the range holds
    no join.
    """
    joins = io_regime_changes(run_dir, range(first, last + 1))
    if not joins:
        return None, 0.0
    join = joins[0]
    before = series[max(first, join - sample):join]
    after = series[join:min(last + 1, join + sample)]
    if len(before) == 0 or len(after) == 0:
        return join, float("nan")
    return join, float(np.mean(after) - np.mean(before))

# THE REPORT EACH DECLARED BOUND WAS TAKEN BESIDE, recorded so that a bound
# nobody has looked at can be told from one that has been looked at and held.
# The scatter was read off this run's settled window; the relaxation time is
# derived by every report and needs no single anchor, so it has none.
#
# If the run is extended or re-assessed this observation moves, and the check
# refuses -- not because the bound is wrong but because it is now unexamined,
# which is the whole state the audit found and could not distinguish.
CONVERGENCE_BOUND_ANCHOR = {
    # THE FLUX BRACKET'S WARM ARM, 55 clean production orbits at one I/O regime.
    # It is the noisiest current settled production run. The bootstrap remains
    # the memory-time anchor in lib/run_lengths.py; the two inputs need not be
    # observations of the same run because each is independently bounded over
    # every qualifying report.
    "run": "run_9d5dbf9bd3d9",
    "artifact": "exoplasim/analysis/convergence/"
                "run_9d5dbf9bd3d9_convergence.json",
    "orbit_scatter_k": 0.14491617956981903,
}
CONVERGENCE_REPORTS = "exoplasim/analysis/convergence"


def check_convergence_bounds(root) -> list[str]:
    """The window's three nominal inputs against the reports they must bound.

    A CHECK WITH A RIGHT ANSWER in each arm, and the arms answer different
    questions:

      IDENTITY   `NOMINAL_TAU_ORBITS` is the top of `lib/run_lengths.py`'s
                 memory bracket and not a second number that happens to match.
      ANCHOR     the report the scatter was read off is on disk and still reads
                 what was recorded. Either failure means the bound is
                 unexamined, which is not the same as its holding.
      BOUND      no report reads a scatter above the declared scatter, and none
                 derives a relaxation time above the declared one. A report
                 BELOW is the expected state of an upper bound and is not a
                 finding.

    The scatter is tested only against reports that read the SETTLED
    trajectory, by `run_lengths.report_is_settled_production`: a diagnostic
    arm's scatter is the experiment's, and a window too short to fix its own tau
    is an instrument below its own scatter. The relaxation time is tested
    against every report, because it is derived from the modelled mixed layer
    and the run's own albedo rather than measured on a window, so a diagnostic
    report's value is as much a property of this model as any other's.

    Takes the repository root, and returns the empty list when the bounds hold.
    """
    base = Path(root)
    problems: list[str] = []

    top = max(run_lengths.TAU_MEMORY_ORBITS_BRACKET)
    if NOMINAL_TAU_ORBITS != top:
        problems.append(
            f"the window is sized on a memory time of {NOMINAL_TAU_ORBITS} "
            f"orbits and lib/run_lengths.py's bracket tops out at {top}. The "
            "window and the production span are bought in one number and this "
            "file states none of its own.")

    anchor = base / CONVERGENCE_BOUND_ANCHOR["artifact"]
    if not anchor.is_file():
        problems.append(
            f"the orbit scatter of {NOMINAL_ORBIT_SCATTER_K} K was read off "
            f"{CONVERGENCE_BOUND_ANCHOR['artifact']}, which is not on disk. "
            "Nothing can say whether the bound still holds, which is not the "
            "same as its holding: re-take it on a run that exists and record "
            "the new anchor.")
    else:
        report = json.loads(anchor.read_text(encoding="utf-8"))
        read = report.get("resolving_power", {}).get(
            "temperature_residual_scatter_k")
        recorded = CONVERGENCE_BOUND_ANCHOR["orbit_scatter_k"]
        if not isinstance(read, (int, float)):
            problems.append(
                f"{CONVERGENCE_BOUND_ANCHOR['artifact']} no longer reports "
                "resolving_power.temperature_residual_scatter_k, so the "
                "declared orbit scatter has nothing to be re-examined against")
        elif not math.isclose(float(read), recorded, rel_tol=1e-12):
            problems.append(
                f"the orbit scatter was anchored to {recorded} K on "
                f"{CONVERGENCE_BOUND_ANCHOR['run']} and that report now reads "
                f"{float(read)} K. The window this bound sizes was derived "
                "before that reading changed, so it is unexamined rather than "
                "wrong: re-take the scatter and record the new anchor.")

    derived_seen = False
    for path in sorted((base / CONVERGENCE_REPORTS).glob("*_convergence*.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except Exception:                                     # noqa: BLE001
            continue
        expected = report.get("metrics", {}).get("relaxation_orbits_expected")
        if isinstance(expected, (int, float)) and np.isfinite(expected):
            derived_seen = True
            if float(expected) > NOMINAL_RELAXATION_ORBITS:
                problems.append(
                    f"{path.name} derives a relaxation time of "
                    f"{float(expected):.4f} orbits and the window is sized on "
                    f"{NOMINAL_RELAXATION_ORBITS}. The window grows as the cube "
                    "root of this time, so a nominal below what a run derives "
                    "under-sizes the window that run is judged in.")
        if not run_lengths.report_is_settled_production(path.name, report):
            continue
        scatter = report.get("resolving_power", {}).get(
            "temperature_residual_scatter_k")
        if not isinstance(scatter, (int, float)) or not np.isfinite(scatter):
            continue
        if float(scatter) > NOMINAL_ORBIT_SCATTER_K:
            problems.append(
                f"{path.name} reads an orbit scatter of {float(scatter):.5f} K "
                f"on its settled production window and the window is sized on "
                f"{NOMINAL_ORBIT_SCATTER_K} K. The window grows as the two "
                "thirds power of the scatter, so this is the bound being wrong "
                "rather than merely old.")
    if not derived_seen:
        problems.append(
            f"no report under {CONVERGENCE_REPORTS} derives a relaxation time, "
            f"so the declared {NOMINAL_RELAXATION_ORBITS} orbits bounds "
            "nothing that can be read. The bound is unexamined.")
    return problems


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
#   RECONVERGENCE after a timestep change or a resolution conversion: ten to
#   twenty orbits was what this project repeatedly saw. The state is already at
#   a climate; what is settling is the model's response to a changed
#   discretisation.
#
# THAT BAND IS NO LONGER WHAT THE CRITERION DERIVES, and the experience is kept
# because the difference is informative rather than because it is still the
# answer. A reconvergence record has to be the longer of two lengths, both
# derived above: `DEFAULT_RECONVERGENCE_ORBITS`, the orbits the statistic needs
# before it resolves its own threshold, and `orbits_to_departure_tolerance`,
# the orbits the modelled transient needs to fall under the allowance. The
# first depends on the declared orbit scatter and the second does not depend on
# it at all, and the scatter has since been raised from 0.07 K -- the T21
# figure the experience was accumulated at -- to a bound over every settled
# production report. The floor moved with it, so the derivation now costs more
# orbits than the experience did. The band was experience and the floor is
# derived from inputs that each carry a citation, so the floor governs.
#
# The consequence for a window: a 10-orbit window is a seventh of a cold
# start's approach and the whole of a reconvergence, so on the second the
# COLD-START criteria are assessing orbits that are still moving by
# construction. That is an argument for asking a reconvergence a different
# question, not for charging it a cold start's window.
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
    segment is no more a trend than the other way round. So is the I/O-regime
    rule, and unmodified: an A/B arm is judged on the same instrument on both
    sides of its window or it is not judged.
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
    refuse_a_window_spanning_an_io_regime_change(run_dir, start, end, window)
    return start, end


CONVERGENCE_LENGTHS = {
    "basis": "operational experience across this project's runs, not a "
             "measured distribution; every number here is bracketed by what "
             "has been seen rather than fitted",
    "cold_start_t21_orbits": 70,
    "cold_start_note": "seed- and initial-condition-dependent; "
                       "exoplasim/notes/parameter-decisions.md records a cold "
                       "start converging on all six criteria after 70 orbits",
    "reconvergence_orbits_floor": DEFAULT_RECONVERGENCE_ORBITS,
    "reconvergence_note": "the state is already at a climate; what settles is "
                          "the response to a changed discretisation. NO BAND "
                          "IS CARRIED HERE: the record length is the larger of "
                          "reconvergence_orbits_floor, which is derived from "
                          "the same nominal scatter, memory and relaxation the "
                          "cold-start window is, and "
                          "relaxation * ln(departure / tolerance), which the "
                          "reconvergence criterion computes per run from the "
                          "departure it FITS. Ten to twenty orbits is what "
                          "this project saw when the declared orbit scatter "
                          "was 0.07 K; the floor moved when that bound did.",
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
                            fraction: float | None = None):
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
    # `None` and not the constant as a default: a default is bound once, at
    # definition, so a constant spelled there is a SNAPSHOT of it and moving the
    # constant would leave this reading the old value. That is the same defect
    # in miniature as the three copies the constant replaced.
    fraction = FIT_TAIL_FRACTION if fraction is None else float(fraction)
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


def conversion_context(run_dir: Path) -> dict | None:
    """What this run started from, or None when it cold-started.

    THE RUN SAYS WHETHER IT IS A RECONVERGENCE; nothing here infers it. A run
    seeded from another run's restart records `initial_state.restart_from_run`
    and, where the support changed, a `conversion` block naming the truncations
    it came from and went to. Both are written by `run_exoplasim.py` at prepare
    time, which is the only place that knows.

    Orbit zero of such a run IS the conversion: `convert_restart.py` produces an
    initial condition and the target model owns the equilibrium from there. That
    is the stamped start the reconvergence criterion needs and the cold-start
    criteria do not have.

    `kind` is reported and never branched on. One criterion covers both kinds
    because both are the same slab relaxing from a displaced state and they
    differ only in how far it was displaced, which the fit measures.
    """
    manifest_path = Path(run_dir) / "run_manifest.json"
    if not manifest_path.is_file():
        return None
    initial = json.loads(manifest_path.read_text(encoding="utf-8")).get(
        "initial_state") or {}
    if initial.get("cold_start", True):
        return None
    conversion = initial.get("conversion") or {}
    truncations = (conversion.get("from"), conversion.get("to"))
    kind = ("resolution_conversion"
            if all(truncations) and truncations[0] != truncations[1]
            else "timestep_change_or_seeded_restart")
    return {"kind": kind,
            "restart_from_run": initial.get("restart_from_run"),
            "from_truncation": truncations[0],
            "to_truncation": truncations[1],
            "conversion_report": conversion.get("report"),
            "conversion_orbit_index": 0}


def expected_relaxation(run_dir: Path, net_down_w_m2: float,
                        upward_w_m2: float) -> dict:
    """The slab's relaxation time, from the planet and this run's own albedo.

    ONE IMPLEMENTATION, TWO CRITERIA. The offset criterion's fallback and the
    reconvergence criterion both hold the relaxation time fixed at this value,
    and two derivations of it would be two ladders.

    The year comes from the run's own manifest rather than the current config:
    this is a property of the run being assessed, which may predate a baseline
    re-run that moved the orbit. The damping is `lib/sensitivity.py` evaluated
    at the albedo the run itself reports, so it follows the run the way the
    year does, and the heat capacity is the modelled mixed layer with the run's
    own sea water in it.
    """
    manifest_path = Path(run_dir) / "run_manifest.json"
    year_days = 182.8
    if manifest_path.is_file():
        year_days = float(json.loads(manifest_path.read_text(encoding="utf-8"))
                          ["derived_parameters"]["orbital_year_earth_days"])
    alpha = planetary_albedo_from_fluxes(net_down_w_m2, upward_w_m2)
    feedback = radiative_damping_w_m2_per_k(alpha)
    capacity, water = slab_heat_capacity(Path(run_dir))
    return {"orbital_year_earth_days": year_days,
            "planetary_albedo": alpha,
            "radiative_damping_w_m2_per_k": feedback,
            "slab_heat_capacity_j_m2_k": capacity,
            "slab_sea_water": water,
            "slab_mixed_layer_depth_m": _MLD,
            "relaxation_orbits": relaxation_orbits(year_days, feedback,
                                                   capacity)}


def reconvergence_assessment(run_dir: Path, conversion: dict,
                             temperature: np.ndarray, sea_ice: np.ndarray,
                             relaxation: dict) -> dict:
    """The verdict on whether a conversion's transient has decayed.

    `temperature` and `sea_ice` are the whole POST-CONVERSION production record,
    starting at the conversion. The window is not a tail: the early orbits are
    where the departure is largest and are what determines its amplitude, so
    trimming them would throw away the evidence this criterion exists to use.
    """
    run_dir = Path(run_dir)
    tau_relax = relaxation["relaxation_orbits"]
    # THE AUTOCORRELATION TIME IS DECLARED, NOT TAKEN FROM THIS RECORD. A tau
    # estimated inside a span of this length comes back well short of the truth
    # (lib/autocorrelation.py), and understating it understates the half width,
    # which passes runs that should be refused. The record's own value is used
    # only where it is LARGER than the declared one -- more memory than nominal
    # is a fact about this run and widening on it is conservative.
    fit = departure_decay(temperature, tau_relax, NOMINAL_TAU_ORBITS)
    record_tau = fit.get("residual_tau_from_this_record", float("nan"))
    tau_used = NOMINAL_TAU_ORBITS
    if np.isfinite(record_tau) and record_tau > NOMINAL_TAU_ORBITS:
        tau_used = float(record_tau)
        fit = departure_decay(temperature, tau_relax, tau_used)

    remaining = fit["remaining_departure_k"]
    half_width = fit["half_width_k"]
    statistic = abs(remaining) + half_width if np.isfinite(remaining) \
        and np.isfinite(half_width) else float("nan")
    n = int(np.asarray(temperature).size)
    ice_slope = slope(sea_ice, n)
    ice_slope_se = slope_standard_error(sea_ice, n)
    # THE RECORD IS THE WHOLE TRANSIENT, SO IT CAN CROSS AN I/O JOIN, and a step
    # at the join enters the fitted level and the fitted amplitude directly --
    # the same way it enters the cold-start asymptote, and refused on the same
    # terms. Measured against the criterion it would decide rather than refused
    # outright: a step at or above the allowance means the fit cannot decide
    # this criterion whatever it returns.
    io_join, io_step_k = io_step_at_join(np.asarray(temperature, dtype=float),
                                         run_dir, 0, n - 1)

    criteria = {
        # The distance still to travel, against the same allowance the
        # cold-start offset criterion is held to. Fails closed: a fit that
        # cannot be taken is not evidence that a transient has gone.
        f"remaining_departure_lt_{_OFFSET_TOLERANCE_K}_k":
            bool(np.isfinite(statistic) and statistic < _OFFSET_TOLERANCE_K),
        # A criterion whose statistic cannot resolve its own threshold does not
        # discriminate, it flips. Stated as a criterion here rather than only
        # reported, because a reconvergence block is short by design and this is
        # the case where that would matter.
        "remaining_departure_resolves_its_threshold":
            bool(np.isfinite(half_width)
                 and half_width * _RESOLVING_FACTOR <= _OFFSET_TOLERANCE_K),
        # THE PRECONDITION, CARRIED RATHER THAN ASSUMED. Fixing the decay at
        # the derived relaxation time is conservative only while that time is a
        # ceiling, and the one regime where it is not is a state whose sea ice
        # is still reorganising: the model's measured response across the ice
        # transition is much larger than away from it, and the relaxation time
        # goes with it.
        f"abs_sea_ice_slope_lt_{SEA_ICE_SLOPE_TOLERANCE_PER_ORBIT}_per_orbit":
            bool(np.isfinite(ice_slope)
                 and abs(ice_slope) < SEA_ICE_SLOPE_TOLERANCE_PER_ORBIT),
        # ONE INSTRUMENT ACROSS THE RECORD, or the step is small enough not to
        # decide the verdict. A record with no join passes this trivially.
        f"io_regime_step_lt_{_OFFSET_TOLERANCE_K}_k":
            bool(io_join is None or (np.isfinite(io_step_k)
                                     and abs(io_step_k) < _OFFSET_TOLERANCE_K)),
    }

    # What this record says the reconvergence costs, from its OWN scatter and
    # its OWN fitted departure rather than the nominal figures the default is
    # priced from. Two independent requirements; the record length is the larger.
    for_statistic = orbits_for_departure_decay(
        fit["residual_scatter_k"], tau_used, tau_relax,
        _OFFSET_TOLERANCE_K / _RESOLVING_FACTOR)
    for_transient = orbits_to_departure_tolerance(
        fit["departure_amplitude_k"], tau_relax, _OFFSET_TOLERANCE_K)
    needed = float(np.nanmax([for_statistic, for_transient])) \
        if np.any(np.isfinite([for_statistic, for_transient])) else float("nan")

    return {
        "conversion": conversion,
        "orbits_assessed": n,
        "criterion": ("the departure the conversion introduced, decaying at "
                      "the relaxation time the planet's slab and this run's "
                      "own radiative damping imply, has fallen inside the same "
                      "allowance the cold-start offset criterion is held to"),
        "tolerance_k": _OFFSET_TOLERANCE_K,
        "resolving_factor": _RESOLVING_FACTOR,
        "statistic_k": statistic,
        "fit": fit,
        "relaxation": relaxation,
        "autocorrelation_time_orbits_used": tau_used,
        "autocorrelation_time_source":
            "declared nominal (lib/autocorrelation.py: a tau taken from inside "
            "a window this short is biased low)" if tau_used == NOMINAL_TAU_ORBITS
            else "this record's residuals, which exceed the declared nominal",
        "sea_ice_slope_fraction_per_orbit": ice_slope,
        "sea_ice_slope_standard_error_fraction_per_orbit": ice_slope_se,
        "io_regime_join_orbit": io_join,
        "io_regime_step_k": None if io_join is None else io_step_k,
        "orbits_for_the_statistic_to_resolve": for_statistic,
        "orbits_for_the_transient_to_fall_under_tolerance": for_transient,
        "orbits_this_reconvergence_needs": needed,
        # WHAT THE FIT PREDICTS, NOT WHAT THE VERDICT SAYS. It can read zero
        # while the criteria still refuse: the prediction is about the modelled
        # transient and the verdict is about this record's own scatter as well.
        # When the two disagree the verdict governs and the run continues by the
        # increment the escalation declares.
        "orbits_still_needed_by_the_fitted_decay":
            (float(max(0.0, np.ceil(needed) - n))
             if np.isfinite(needed) else float("nan")),
        "default_reconvergence_orbits": DEFAULT_RECONVERGENCE_ORBITS,
        "criteria": criteria,
        "failed_criteria": sorted(name for name, ok in criteria.items()
                                  if not ok),
        "reconverged_after_conversion": bool(all(criteria.values())),
        "not_a_substitute_for_convergence":
            "This verdict says a conversion's transient has decayed. It does "
            "NOT say the run is at the climate the world's numbers may be "
            "taken from: that is the cold-start criteria, and a run whose "
            "climatology is going to be read is still held to all of them.",
    }


def post_conversion_record(run_dir: Path, n_orbits: int) -> int:
    """The last production orbit of a record that STARTS at the conversion.

    The counterpart of `production_window` for a criterion whose record is the
    whole run rather than a tail: a diagnostic tail is dropped, and a
    diagnostic orbit INSIDE the record refuses, because a decay fitted across a
    block that was not the run's trajectory is not a decay.
    """
    excluded = set(non_production_orbits(run_dir, range(n_orbits)))
    end = n_orbits - 1
    while end >= 0 and end in excluded:
        end -= 1
    if end < 0:
        raise RuntimeError(
            f"{Path(run_dir).name} has no production orbits to assess")
    interior = sorted(o for o in excluded if o <= end)
    if interior:
        raise RuntimeError(
            f"{Path(run_dir).name} declares orbits {interior} non-production "
            "inside the post-conversion record. A decay fitted across a "
            "diagnostic block is not a decay; assess the block before the "
            "interruption instead.")
    return end


def assess_reconvergence_only(run_dir: Path, files: list[Path],
                              conversion: dict, output: Path,
                              through: int | None) -> dict:
    """Report the reconvergence verdict for a run shorter than the cold-start window.

    A RECONVERGENCE IS SHORTER THAN THAT WINDOW BY DESIGN. Refusing to report
    anything until the run is long enough to answer the cold-start question is
    how the criterion that costs `DEFAULT_WINDOW_ORBITS` orbits came to be the
    only criterion a reconvergence had. The cold-start criteria are not
    evaluated here and are NOT reported as failed: they are reported as not
    evaluable, which is a different statement and the honest one.
    """
    end = post_conversion_record(run_dir, len(files))
    records = annual_records(files[:end + 1])
    arrays = {key: np.array([record[key] for record in records])
              for key in ["ts", "sic", "rst", "rsut"]}
    relaxation = expected_relaxation(run_dir, float(arrays["rst"].mean()),
                                     float(arrays["rsut"].mean()))
    block = reconvergence_assessment(run_dir, conversion, arrays["ts"],
                                     arrays["sic"], relaxation)
    report = {
        "run_dir": str(run_dir),
        "completed_orbits": len(files),
        "window_start_year_index": 0,
        "window_end_year_index": end,
        "orbits_after_window_excluded_as_non_production":
            sorted(range(end + 1, len(files))),
        "assessed_purpose": "production",
        "segment_purposes": {str(k): v for k, v in
                             sorted(orbit_purposes(run_dir,
                                                   range(len(files))).items())},
        "reconvergence": block,
        "convergence_criteria_evaluated": False,
        "convergence_criteria_not_evaluated_because":
            f"the cold-start criteria are taken over {DEFAULT_WINDOW_ORBITS} "
            f"orbits and this run has {len(files)}. They are not failed here, "
            "they are unasked.",
        "convergence_lengths": CONVERGENCE_LENGTHS,
        "annual_records": records,
    }
    output.mkdir(parents=True, exist_ok=True)
    suffix = "" if through is None else f"_through{through:03d}"
    report_path = output / f"{run_dir.name}_reconvergence{suffix}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    draw_reconvergence(records, block,
                       output / f"{run_dir.name}_reconvergence{suffix}.png")

    payload = {"reconvergence": block["reconverged_after_conversion"],
               "statistic_k": block["statistic_k"],
               "remaining_departure_k": block["fit"]["remaining_departure_k"],
               "half_width_k": block["fit"]["half_width_k"],
               "orbits_assessed": block["orbits_assessed"],
               "orbits_still_needed_by_the_fitted_decay":
                   block["orbits_still_needed_by_the_fitted_decay"],
               "failed_criteria": block["failed_criteria"],
               "report": str(report_path.resolve())}
    manifest_path = run_dir / "run_manifest.json"
    # The verdict travels with the run, beside the convergence one and never
    # instead of it: `status` is what the cold-start criteria set, and they were
    # not asked here.
    if through is None and manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["reconvergence_assessment"] = payload
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n",
                                 encoding="utf-8")
    return payload


def draw_reconvergence(records: list[dict], block: dict, path: Path) -> None:
    """The record, the fitted decay, and where the fit says it is going.

    A reader has to be able to see that the transient the criterion measured is
    the one in the record, which a pass/fail line cannot show.
    """
    ts = np.array([record["ts"] for record in records])
    orbits = np.arange(ts.size, dtype=float)
    fit = block["fit"]
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    ax.plot(orbits, ts, marker="o", markersize=3, linewidth=1.1,
            label="global surface temperature")
    tau = block["relaxation"]["relaxation_orbits"]
    if np.isfinite(fit["level_k"]) and np.isfinite(tau):
        ax.plot(orbits,
                fit["level_k"] + fit["departure_amplitude_k"]
                * np.exp(-orbits / tau),
                linewidth=1.4,
                label=f"decay fitted at the derived {tau:.1f}-orbit relaxation")
        ax.axhline(fit["level_k"], color="black", linewidth=0.7,
                   label="where the transient is going")
        ax.axhspan(fit["level_k"] - block["tolerance_k"],
                   fit["level_k"] + block["tolerance_k"],
                   color="tab:green", alpha=0.10,
                   label=f"+/- {block['tolerance_k']} K allowance")
    ax.set_xlabel("orbits since the conversion")
    ax.set_ylabel("K")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("Reconvergence: "
                 + ("DECAYED" if block["reconverged_after_conversion"]
                    else "NOT YET")
                 + f" ({block['orbits_assessed']} orbits since the conversion)")
    fig.savefig(path, dpi=180)
    plt.close(fig)

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
    # WHICH QUESTION THIS RUN IS ASKING, taken from what the run says it began
    # from rather than from a flag. A run seeded from a conversion has a stamped
    # transient start and a bounded departure; a cold start has neither, and the
    # two want different criteria.
    conversion = conversion_context(run_dir)
    if len(files) < args.window:
        if conversion is None:
            raise RuntimeError("Not enough annual outputs for requested window")
        if args.assess != "production":
            raise RuntimeError(
                f"{run_dir.name} has {len(files)} orbits, fewer than the "
                f"{args.window} the cold-start criteria are taken over. It "
                f"began from a conversion, so the reconvergence criterion "
                f"could answer for it -- but that criterion is taken over the "
                f"whole post-conversion PRODUCTION record and there is no "
                f"`{args.assess}` form of it. Assess it as production.")
        print(json.dumps(assess_reconvergence_only(
            run_dir, files, conversion, args.output, args.through), indent=2))
        return

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

    records = annual_records(files)

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
    # The albedo is this run's own, over the same window every other metric is
    # taken over, so the damping moves with the run exactly as `year_days` does.
    # `expected_relaxation` is the ONE derivation of this time, shared with the
    # reconvergence criterion, which holds the same time fixed.
    relaxation = expected_relaxation(
        run_dir, float(arrays["rst"][-w:].mean()),
        float(arrays["rsut"][-w:].mean()))
    year_days = relaxation["orbital_year_earth_days"]
    window_alpha = relaxation["planetary_albedo"]
    feedback_w_m2_k = relaxation["radiative_damping_w_m2_per_k"]
    slab_capacity = relaxation["slab_heat_capacity_j_m2_k"]
    slab_water = relaxation["slab_sea_water"]
    tau_expected = relaxation["relaxation_orbits"]
    orbits_axis = np.arange(len(arrays["ts"]), dtype=float)
    asymptote, half_width, tau_fit, tau_fit_se = approach_to_equilibrium(
        orbits_axis, arrays["ts"])
    # How many orbits the fit was taken over, from the one function that knows.
    # A fitted tau longer than the span it was fitted over is extrapolation
    # rather than a measurement, which is what decides whether it is evidence
    # about the relaxation time at all.
    fit_span = fit_span_orbits(len(orbits_axis))
    if not np.isfinite(tau_fit):
        relaxation_fit_identifiable = False
        relaxation_fit_verdict = "the series would not support a fit at all"
    elif tau_fit <= 0.0:
        relaxation_fit_identifiable = False
        relaxation_fit_verdict = (
            f"the fit returned tau = {tau_fit:.6g} orbits; a negative time "
            "constant is the exponential having collapsed to a line")
    elif tau_fit > fit_span:
        relaxation_fit_identifiable = False
        relaxation_fit_verdict = (
            f"the fit returned tau = {tau_fit:.6g} orbits over {fit_span} "
            "orbits of data, so the exponential is indistinguishable from a "
            "line across the fitted range and tau is not determined by it")
    else:
        relaxation_fit_identifiable = True
        relaxation_fit_verdict = (
            f"tau = {tau_fit:.6g} orbits is shorter than the {fit_span} "
            "orbits fitted, so the series carries the curvature that fixes it")
    # The fitted range is the tail `fit_span_orbits` selects, so the join is
    # looked for there and not over the whole series.
    fit_first = len(orbits_axis) - fit_span
    io_join, io_step_k = io_step_at_join(
        arrays["ts"], run_dir, fit_first, len(orbits_axis) - 1)
    if io_join is not None and (
            not np.isfinite(io_step_k) or abs(io_step_k) >= _OFFSET_TOLERANCE_K):
        relaxation_fit_identifiable = False
        relaxation_fit_verdict = (
            f"the fitted range crosses the I/O-regime change at orbit "
            f"{io_join}, and the surface temperature steps {io_step_k:+.4g} K "
            f"across it against the {_OFFSET_TOLERANCE_K} K the offset "
            f"criterion discriminates at. The step enters the asymptote "
            f"directly, so this fit cannot decide that criterion whatever it "
            f"returns; extend the run in one regime, or read the relaxation "
            f"time off the clean block alone")
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
    # AN ESTIMATOR IS CHOSEN BY WHETHER IT CAN DECIDE THIS CRITERION, and the
    # bar is the one this file already declares and this report already prints
    # beside every verdict: a threshold discriminates when the standard error of
    # the statistic it tests is at most `_RESOLVING_FACTOR` of it. No new
    # number enters here. The offset criterion's `resolving_power` row applies
    # exactly that test to exactly this half width, and it was reported and read
    # by nothing while the branch above it selected on four sanity bounds.
    #
    # WHAT THAT COST, AND IT IS AN INVERSION RATHER THAN A LOOSENING. The old
    # test asked only that the asymptote and its half width were finite, the
    # offset under 20 K and the half width under 5 K. A half width of 0.19 K
    # clears that by a factor of 26 and then becomes the criterion's own error
    # against a 0.15 K allowance, so the run FAILS. A fit that collapses
    # outright trips the 20 K bound instead, routes to the drift fallback, and
    # PASSES on the same data. The criterion was therefore easier to pass the
    # worse the fit was, and the flux-to-kelvin bracket's two arms -- identical
    # in build, staged surface, executable, structure and window indices, and
    # differing only in the flux they were given -- landed on opposite sides of
    # it on opposite windows, with every other criterion passing on all four
    # assessments. `notes/audits/flux-slope-bracket.md`; world-jejw.
    #
    # WHY THE HALF WIDTH IS THE RIGHT INSTRUMENT AND THE FITTED TAU IS NOT.
    # The row that found this named the fit's tau standard error -- 40.35 orbits
    # against a tau of 39.99 -- and that number is real, but it is not what this
    # criterion uses. The statistic is the ASYMPTOTE against the window mean, so
    # the question is whether the asymptote is determined, and `half_width` is
    # that answer with tau's uncertainty already inside it: it comes off the same
    # covariance and is then widened by the overshoot a tau comparable to the
    # span implies. Gating on tau instead would also drop fits whose tau is
    # loosely determined and whose asymptote is tight, and it would empty
    # `check_relaxation_ceiling.py`'s evidence set, which deliberately keeps an
    # uncertain tau in the bracket because narrowing a bracket on the strength
    # of a reading being uncertain is backwards.
    #
    # A FIT WHOSE TAU IS NOT IDENTIFIABLE STILL DOES NOT DECIDE THE CRITERION.
    # `relaxation_fit_identifiable` is computed above and was governing the
    # reported field only; a negative tau, a tau longer than the span fitted, or
    # a fitted range crossing an I/O-regime step at or above this very tolerance
    # are all cases where the fit's asymptote is not a reading of where the run
    # is going, and the last of those says so in its own verdict text.
    #
    # WHAT IT DOES TO THIS TREE, stated because a branch nobody can reach is
    # worse than no branch: NO fit on any run currently in this project clears
    # the resolving bar. The tightest asymptote half width on disk is 0.0198 K
    # -- 0.059 against the 0.05 the bar allows -- so every run takes the drift
    # form today. That is the finding rather than a hidden default: the free
    # three-parameter exponential does not resolve a 0.15 K offset on this
    # model's series, and the drift form's window is derived so that it does.
    # The branch stays because the condition is a criterion rather than a
    # verdict about the model: a longer or quieter series can clear it, and
    # where it does, the fit is the better estimator because it does not rest on
    # the derived relaxation time being a ceiling.
    # The fit's own half width, captured before the fallback overwrites it, so
    # the record of the CHOICE says what was on offer rather than what was used.
    fit_half_width_offered = half_width
    fit_resolves_the_threshold = bool(
        np.isfinite(half_width)
        and half_width * _RESOLVING_FACTOR <= _OFFSET_TOLERANCE_K)
    fit_usable = bool(np.isfinite(asymptote) and np.isfinite(half_width)
                      and abs(offset) < 20.0 and half_width < 5.0
                      and relaxation_fit_identifiable
                      and fit_resolves_the_threshold)
    # Why the fit was declined, named rather than left to be inferred from the
    # source label. `None` where the fit was taken.
    if fit_usable:
        offset_fit_declined_because = None
    elif not (np.isfinite(asymptote) and np.isfinite(half_width)):
        offset_fit_declined_because = "the series would not support a fit"
    elif not (abs(offset) < 20.0 and half_width < 5.0):
        offset_fit_declined_because = (
            f"the fit returned an offset of {offset:.4g} K with a half width "
            f"of {half_width:.4g} K, which is the exponential having nothing "
            f"to grip rather than a distance to equilibrium")
    elif not relaxation_fit_identifiable:
        offset_fit_declined_because = relaxation_fit_verdict
    else:
        offset_fit_declined_because = (
            f"the fitted asymptote's half width is {half_width:.4g} K and this "
            f"criterion discriminates at {_OFFSET_TOLERANCE_K} K, so the fit's "
            f"own error is more than {_OFFSET_TOLERANCE_K / _RESOLVING_FACTOR:g}"
            f" K and it cannot decide the threshold whatever it returns")
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
            None if relaxation_fit_identifiable else float(fit_span)),
        # The fit's own error on tau, autocorrelation-corrected like the
        # asymptote's. Reported and thresholded by nothing: it is what makes the
        # derived time's CEILING claim testable rather than a point comparison.
        "relaxation_orbits_fitted_standard_error": (
            tau_fit_se if relaxation_fit_identifiable else None),
        "relaxation_fit_span_orbits": fit_span,
        # Whether the approach and the settled tail were written by the same
        # instrument, and what the join is worth if not. Reported on every run
        # -- including `null` -- so that a fit nobody has checked for a join
        # can be told from one that has been checked and has none.
        "relaxation_fit_io_regime_join_orbit": io_join,
        "relaxation_fit_io_regime_step_k": (
            None if io_join is None else io_step_k),
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
        {"key": f"abs_sea_ice_slope_lt_{SEA_ICE_SLOPE_TOLERANCE_PER_ORBIT}"
                f"_per_orbit",
         "statistic": metrics["sea_ice_slope_fraction_per_orbit"],
         "interval": metrics["sea_ice_slope_standard_error_fraction_per_orbit"],
         "standard_error": metrics["sea_ice_slope_standard_error_fraction_per_orbit"],
         "threshold": SEA_ICE_SLOPE_TOLERANCE_PER_ORBIT,
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
    # HOW THAT ESTIMATOR WAS CHOSEN, in the same terms the rows are priced in.
    # The choice is a criterion and not a fallback of last resort, so it is
    # recorded rather than inferred from the label.
    resolving["offset_statistic_selection"] = {
        "rule": ("the exponential fit decides this criterion only where it can "
                 "resolve the criterion's own threshold at the same "
                 "`factor` every row above is priced against, and where its "
                 "fitted relaxation time is identifiable from the series. "
                 "Otherwise the drift form, whose error the verdict window is "
                 "derived to resolve."),
        "fit_asymptote_half_width_k": float(fit_half_width_offered),
        "fit_resolves_the_threshold": fit_resolves_the_threshold,
        "relaxation_fit_identifiable": relaxation_fit_identifiable,
        "fit_declined_because": offset_fit_declined_because,
    }

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
    # BOTH QUESTIONS, WHEN BOTH APPLY. A run long enough for the cold-start
    # criteria that BEGAN from a conversion still has a stamped transient in it,
    # and reporting when that transient decayed is what says whether a surprise
    # after the next conversion is attributable to the support alone. The two
    # verdicts are separate keys and neither overrides the other; `status` and
    # `sufficiently_equilibrated_for_worldbuilding` stay what the cold-start
    # criteria set.
    #
    # ONLY IN THE PRODUCTION MODE. The record this criterion fits starts at the
    # conversion and runs to the last production orbit, which is not a window
    # and has no diagnostic form.
    if conversion is not None and args.assess == "production":
        try:
            reconv_end = post_conversion_record(run_dir, len(files))
        except RuntimeError as exc:                # noqa: BLE001 - reported
            report["reconvergence_not_assessed_because"] = str(exc)
        else:
            reconv = {key: np.array([r[key] for r in records[:reconv_end + 1]])
                      for key in ["ts", "sic", "rst", "rsut"]}
            report["reconvergence"] = reconvergence_assessment(
                run_dir, conversion, reconv["ts"], reconv["sic"],
                expected_relaxation(run_dir, float(reconv["rst"].mean()),
                                    float(reconv["rsut"].mean())))
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
    if "reconvergence" in report:
        payload["reconvergence"] = \
            report["reconvergence"]["reconverged_after_conversion"]
        payload["reconvergence_failed_criteria"] = \
            report["reconvergence"]["failed_criteria"]

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
