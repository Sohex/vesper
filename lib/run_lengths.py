"""How many orbits to BUY, derived from the timescales rather than typed.

WORLDBUILDING CONTEXT: Vesper is a fictional planet and this module derives run
lengths for the climate model that simulates it. Every orbit below is a modelled
orbit.

WHY THIS EXISTS. Two instruments refuse AFTER the wall clock is spent:
`assess_convergence.py` reports the window its criteria would have needed, and
`compare_equilibria.py` returns INDETERMINATE. Nothing consumed a length when a
run was LAUNCHED, so the orbit count was typed, and a 40-orbit precision arm
that reversed sign when it was re-run at 85 is what that cost.

TWO DIFFERENT TIMES ARE BOTH CALLED TAU AND THEY ARE NOT INTERCHANGEABLE.
Conflating them is how a span derived for one gets spent on the other.

    tau_memory      the INTEGRATED AUTOCORRELATION TIME of the stationary
                    variability: how long the model's own wobble stays
                    correlated with itself. It sets how long a span must be
                    before a mean taken over it has an interval that covers.
                    `lib/autocorrelation.py` owns its arithmetic.
    tau_relaxation  the E-FOLDING TIME OF AN APPROACH to equilibrium: how fast a
                    perturbed state returns. It sets how long a settling block
                    must be before the perturbation has gone.

A commissioning span is bought in tau_memory and a settling block in
tau_relaxation, and neither number substitutes for the other.

NEITHER IS MEASURED, SO EVERY ANSWER HERE IS A BRACKET. Fixing tau_memory by
direct estimation needs of order a thousand orbits at equilibrium and is out of
reach; tau_relaxation has three usable fits on record and they do not agree with
the derived value closely enough to replace it.
`exoplasim/notes/convergence-lengths.md` carries both measurements and
`exoplasim/scripts/check_relaxation_ceiling.py` keeps the second up to date. The
honest form is the bracket, and the bottom of it is what a run buys up front:
the instruments still refuse, and a refusal is what buys the rest.

WHAT A DECLARED SPAN IS AND IS NOT. It is a FLOOR, decided before the run, so
that the length is not chosen by whoever is watching the wall clock. It is not a
stopping rule: a run that reaches it and is still refused by the convergence
criteria is not finished. The alternative of declaring nothing and running until
the instrument stops refusing was priced and rejected, and the reason is in
`docs/src/pipeline/sequencing.md`: the convergence criteria are sized for the
offset criterion's slope error, not for a batch-mean interval, so they stop
refusing about seventy orbits before the climatology's own error bar is bought.
The refusal is a weaker bar than the error bar, so it cannot stand in for it.
"""

from __future__ import annotations

import math

# THE MEMORY TIME, BRACKETED. Both ends are measured and both have a source in
# exoplasim/notes/convergence-lengths.md.
#   4.2   (1+r)/(1-r) at the lag-1 of 0.615 measured on the 85-orbit T21 pair,
#         which is the AR(1) reading of the model's variability.
#   10.43 the fast-plus-slow mixture that section stands the slower branch up
#         with, deliberately not AR(1), because a model with a slow component
#         has a longer memory than its lag-1 admits.
TAU_MEMORY_ORBITS_BRACKET = (4.2, 10.43)

# THE RELAXATION TIME, BRACKETED BY THE FITS THAT ARE EVIDENCE. The three
# convergence artifacts whose exponential fit is usable and identifiable report
# 6.27, 6.92 and 13.68 orbits; the derived value of about 10.1 sits inside that
# and is NOT a ceiling on it, which is why this is a bracket and not the derived
# number. `check_relaxation_ceiling.py` is what would move these.
TAU_RELAXATION_ORBITS_BRACKET = (6.27, 13.68)

# A SPAN OF TWENTY TAU IS WHERE AN HONEST INTERVAL STARTS, and batches of five
# tau are what the span has to be made of. Measured on synthetic series whose
# answer is known in closed form: batches of two tau miss the coverage criterion
# however many of them there are, so batch LENGTH is what has to be bought first.
PRODUCTION_SPAN_TAU_MULTIPLE = 20.0
BATCH_LENGTH_TAU_MULTIPLE = 5.0

# WHAT A SETTLING BLOCK HAS TO DECAY TO. The same 0.15 K the offset criterion
# allows, and deliberately the same number: a residual smaller than what the
# instrument that judges the NEXT state can see is a residual that state cannot
# be held responsible for. `scripts/smoke_test.py` checks the two agree.
SETTLING_RESIDUAL_K = 0.15


def production_span_orbits(tau_memory_orbits: float) -> float:
    """Orbits of production a window mean's interval needs before it covers."""
    return PRODUCTION_SPAN_TAU_MULTIPLE * float(tau_memory_orbits)


def commissioning_orbits(approach_orbits: float,
                         tau_memory_orbits: float) -> float:
    """The approach, then the production span that follows it.

    `approach_orbits` is operational experience and has no artifact behind it:
    about seventy at T21 from cold, ten to twenty on a reconvergence, and it
    moves with the rung. It is an ARGUMENT here rather than a constant so that
    the caller has to state which approach it means.
    """
    return float(approach_orbits) + production_span_orbits(tau_memory_orbits)


def commissioning_bracket(approach_orbits: float) -> tuple[float, float]:
    """The commissioning length over the whole memory-time bracket."""
    low, high = TAU_MEMORY_ORBITS_BRACKET
    return (commissioning_orbits(approach_orbits, low),
            commissioning_orbits(approach_orbits, high))


def settling_orbits(perturbation_k: float,
                    tau_relaxation_orbits: float,
                    residual_k: float = SETTLING_RESIDUAL_K) -> float:
    """Orbits for a perturbation of `perturbation_k` to decay below `residual_k`.

    An approach decays as `exp(-n / tau)`, so `n = tau * ln(A0 / residual)`.
    This is a LENGTH and not a verdict on a mean, which is the whole point: it
    makes no claim about a window mean and therefore does not have to be priced
    against a slope's standard error.

    WHAT IT DOES NOT ESTABLISH: that the state is equilibrated, that its climate
    is the rung's climate, or that any mean taken on it carries an interval.
    Those are the commissioning standard and they need
    `commissioning_orbits`. A settling block is for handing the next conversion
    a restart that is not mid-transient, and the residual it leaves is absorbed
    by the convergence that follows that conversion.
    """
    if perturbation_k <= residual_k:
        return 1.0
    return max(1.0, float(tau_relaxation_orbits)
               * math.log(float(perturbation_k) / float(residual_k)))


def settling_bracket(perturbation_k: float,
                     residual_k: float = SETTLING_RESIDUAL_K
                     ) -> tuple[float, float]:
    """The settling length over the whole relaxation-time bracket."""
    low, high = TAU_RELAXATION_ORBITS_BRACKET
    return (settling_orbits(perturbation_k, low, residual_k),
            settling_orbits(perturbation_k, high, residual_k))
