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
reach; tau_relaxation is bracketed by the fits on record, which do not agree
with the derived value closely enough to replace it.
`exoplasim/notes/convergence-lengths.md` carries the argument behind both. The
honest form is the bracket, and the bottom of it is what a run buys up front:
the instruments still refuse, and a refusal is what buys the rest.

THE TWO BRACKETS ARE HELD TO THEIR MEASUREMENTS BY DIFFERENT MEANS, and the
difference is the whole of what `notes/audits/frozen-derived-quantities.md`
tier 1 asks for.

    tau_relaxation  STATES NO NUMBER. `tau_relaxation_orbits_bracket()` reads
                    the bracket out of the artifact
                    `check_relaxation_ceiling.py` emits, so there is nothing
                    here to go stale against it.
    tau_memory      IS DECLARED, because it is an UPPER BOUND taken by a sweep
                    no artifact reproduces, and a bound may not be replaced by
                    a reading that is a LOWER bound on the same quantity --
                    that would shorten every span on the ladder on the strength
                    of an estimate the project has written down as too small.
                    So it is declared with the observation it was anchored to,
                    and `check_memory_bracket()` refuses when the anchor moves
                    or when any artifact reads above the bound.

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

import json
import math
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]

# THE MEMORY TIME, MEASURED 2026-08-26 on run_432e5e46adef and no longer
# inferred. Both ends are readings of the same series: the per-orbit
# area-weighted mean surface temperature of the bootstrap, reduced by Geyer's
# initial monotone positive sequence over five candidate windows starting at
# orbits 25, 30, 35, 40 and 45. They returned 2.22, 1.89, 2.00, 2.07 and 2.11,
# every one reliable with a span over ten tau, and the bracket is their range.
#
# BOTH ENDS ARE UPPER BOUNDS. `autocorrelation.stationary_enough` refused all
# five windows: the run still drifts 0.003 to 0.007 K per orbit, and a residual
# trend pushes every lag correlation up. A flatter series returns a smaller tau,
# never a larger one, so the production span this multiplies is bounded above.
#
# WHAT IT REPLACED, and why that pair was not wrong when it was written. It was
# 4.2 to 10.43: an AR(1) reading of a lag-1 of 0.615 on the 85-orbit T21 pair,
# and a fast-plus-slow mixture standing up a slower branch. Both were inferences
# from a model this tree no longer has -- since then the Stephens cloud tables,
# a derived orographic roughness and a hyperdiffusion 1.699x shorter have all
# landed, and the last damps the model faster, which is where a shorter memory
# would come from. The bootstrap reads a lag-1 of 0.39 to 0.49 against that
# 0.615.
#
# THE CONSEQUENCE IS THE WHOLE LADDER'S COST. The production span is twenty
# times whichever end is carried, at every commissioning rung, so the pair this
# replaced asked for about five times the orbits this one does.
# `commissioning_bracket` is what says how many; a count written here would be a
# third statement of the same measurement.
#
# `exoplasim/scripts/assess_convergence.py` sizes its default window on the TOP
# of this bracket and states no number of its own, so this is the one statement
# of the memory time in the tree.
TAU_MEMORY_ORBITS_BRACKET = (1.89, 2.45)

# WHERE THE BOUND WAS ANCHORED, so that "the bound still holds" can be told from
# "the bound was never re-examined". Those two look identical in a declared
# constant and they are not the same claim: the first is a reading, the second
# is the absence of one, and a bound nobody has looked at since the series moved
# under it is the frozen quantity this file was audited for.
#
# The bracket was taken by a five-window Geyer sweep over one run's per-orbit
# mean surface temperature. No artifact reproduces that sweep, so the anchor is
# the next best thing that is written down: what that run's own convergence
# report reads for the same quantity, recorded here at the precision the report
# carries. `check_memory_bracket` refuses on three separate conditions, and they
# are different failures:
#
#   THE ANCHOR IS GONE      the report is not on disk. Nothing can be compared,
#                           so the declaration is unexamined and says so.
#   THE ANCHOR HAS MOVED    the report is there and reads something else. The
#                           run was extended or re-assessed, the series the
#                           sweep was taken over no longer exists, and the
#                           sweep has to be re-taken. This is the condition that
#                           closes the loop: a longer run changes the reading,
#                           and the reading refuses the declaration.
#   THE BOUND IS VIOLATED   some report reads a memory time ABOVE the top of the
#                           bracket. That is the bound being wrong rather than
#                           merely old.
#
# A reading BELOW the bracket is none of those. It is what an upper bound taken
# on a drifting series is supposed to look like, and the check reports the
# headroom rather than treating it either as agreement or as a defect.
MEMORY_BRACKET_ANCHOR = {
    "run": "run_76e441e0a761",
    "artifact": "exoplasim/analysis/convergence/run_76e441e0a761_convergence.json",
    "node": ("resolving_power", "temperature_residual_tau_orbits"),
    "observation": 2.442874832379452,
    "why_below": "the anchor is the CARVED build's bootstrap, forty clean "
                 "orbits at one I/O regime, and it reads AT the top of this "
                 "bracket rather than below it. The top moved 2.22 to 2.45 to "
                 "cover it, which is what an upper bound does when a supported "
                 "reading passes it. WHAT THE BRACKET IS NOT: a property of "
                 "the process. Tau on this model grows with the window it is "
                 "measured on -- 1.00 at twenty orbits and 6.30 at a hundred "
                 "and forty on one clean block, every reading supported -- so "
                 "this is a bound on what has been READ at the lengths runs "
                 "are actually bought at, and chasing it upward with longer "
                 "windows would never terminate. That is why it no longer "
                 "prices the production span: `production_span_from_report` "
                 "does, from the run's own criteria, and this sizes only the "
                 "a-priori default window. "
                 "exoplasim/notes/memory-time-and-the-production-span.md.",
}

# WHAT COUNTS AS A READING OF THE SETTLED VARIABILITY, applied to every
# convergence report the bound is tested against. Two conditions, and each
# excludes a class of report that would otherwise falsify the bound with a
# number that is not about the settled trajectory at all:
#
#   PRODUCTION      a `--assess diagnostic` report prices an A/B arm's
#                   deliberately perturbed orbits. Its scatter and its memory
#                   are the experiment's and not the planet's.
#   RESOLVED TAU    a report whose own window is too short to fix its tau says
#                   so, through `tau_span_supports_the_estimate`. Reading a tau
#                   the report itself declines to stand behind is reading the
#                   instrument below its own scatter.
#
# Reports written before `assessed_purpose` existed carry no such key; they are
# production unless their filename says diagnostic, which is the same rule
# `check_relaxation_ceiling.py` applies to artifacts older than its own fields.
CONVERGENCE_REPORTS = "exoplasim/analysis/convergence"


def report_is_settled_production(name: str, report: dict) -> bool:
    """Is this convergence report a reading of the settled trajectory?"""
    purpose = report.get("assessed_purpose")
    if purpose is None:
        purpose = "diagnostic" if "diagnostic" in name else "production"
    if purpose != "production":
        return False
    return bool(report.get("resolving_power", {})
                .get("tau_span_supports_the_estimate", False))


def check_memory_bracket(root=None) -> list[str]:
    """`TAU_MEMORY_ORBITS_BRACKET` against the reports it was anchored to.

    A CHECK WITH A RIGHT ANSWER in each of its three arms: the anchor artifact
    is there or it is not, it reads the recorded observation or it does not, and
    every settled production report reads at or below the top of the bracket or
    one of them does not. Returns the empty list when the bound holds.

    Takes the repository root so this module keeps knowing nothing but the two
    timescales, the same arrangement `lib/rungs.py` uses for the ladder.
    """
    base = Path(root) if root is not None else _REPO_ROOT
    problems: list[str] = []
    _, high = TAU_MEMORY_ORBITS_BRACKET

    anchor = base / MEMORY_BRACKET_ANCHOR["artifact"]
    if not anchor.is_file():
        problems.append(
            f"the memory-time bracket {TAU_MEMORY_ORBITS_BRACKET} was anchored "
            f"to {MEMORY_BRACKET_ANCHOR['artifact']}, which is not on disk. "
            "Nothing can say whether the bound still holds, and that is not "
            "the same as its holding: re-take the sweep on a run that exists "
            "and record the new anchor.")
    else:
        try:
            report = json.loads(anchor.read_text(encoding="utf-8"))
        except Exception as exc:                              # noqa: BLE001
            return problems + [
                f"{MEMORY_BRACKET_ANCHOR['artifact']} cannot be read as a "
                f"convergence report: {exc}"]
        node = report
        for key in MEMORY_BRACKET_ANCHOR["node"]:
            node = node.get(key) if isinstance(node, dict) else None
        recorded = MEMORY_BRACKET_ANCHOR["observation"]
        if not isinstance(node, (int, float)):
            problems.append(
                f"{MEMORY_BRACKET_ANCHOR['artifact']} no longer reports "
                + ".".join(MEMORY_BRACKET_ANCHOR["node"])
                + ", so the memory-time bracket has nothing to be re-examined "
                  "against")
        elif not math.isclose(float(node), recorded, rel_tol=1e-12):
            problems.append(
                f"the memory-time bracket was anchored to a memory time of "
                f"{recorded} on {MEMORY_BRACKET_ANCHOR['run']} and that report "
                f"now reads {float(node)}. The series the sweep was taken over "
                "has moved, so the bracket is unexamined rather than wrong: "
                "re-take the five-window sweep and record the new anchor.")

    for path in sorted((base / CONVERGENCE_REPORTS).glob("*_convergence*.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except Exception:                                     # noqa: BLE001
            continue
        if not report_is_settled_production(path.name, report):
            continue
        tau = report.get("resolving_power", {}).get(
            "temperature_residual_tau_orbits")
        if not isinstance(tau, (int, float)) or not math.isfinite(float(tau)):
            continue
        if float(tau) > high:
            problems.append(
                f"{path.name} reads a memory time of {float(tau):.4f} orbits "
                f"on its settled production window and the declared bracket "
                f"tops out at {high}. The bracket's ends are UPPER BOUNDS, so "
                "this is the bound being wrong and not merely old: every "
                "commissioning span on the ladder is twenty times a number "
                "that no longer bounds the model.")
    return problems


# THE RELAXATION TIME STATES NO NUMBER HERE. `check_relaxation_ceiling.py` reads
# every convergence report, decides which of them are EVIDENCE about the
# relaxation time by a rule it states, and emits the bracket over the fits that
# are. This reads that bracket, so there is nothing in this file to drift
# against it, and a run whose approach the fit can grip moves the bracket the
# moment its report is written.
#
# The evidence rule is the producer's and is not restated here: a second
# statement of it would be the same defect one level up.
RELAXATION_CEILING_ARTIFACT = "exoplasim/analysis/convergence/relaxation_ceiling.json"
RELAXATION_CEILING_GENERATOR = "exoplasim/scripts/check_relaxation_ceiling.py"


def tau_relaxation_orbits_bracket(root=None) -> tuple[float, float]:
    """The relaxation-time bracket, read from the artifact that measures it.

    Raises rather than falling back on a number. A settling block bought from a
    default nobody measured is the failure this whole module exists against, and
    the artifact is regenerable from the tracked convergence reports alone.
    """
    base = Path(root) if root is not None else _REPO_ROOT
    path = base / RELAXATION_CEILING_ARTIFACT
    if not path.is_file():
        raise RuntimeError(
            f"{RELAXATION_CEILING_ARTIFACT} is not there and it is where the "
            "relaxation-time bracket lives. Run "
            f"`python {RELAXATION_CEILING_GENERATOR}`; it reads the convergence "
            "reports and nothing else.")
    result = json.loads(path.read_text(encoding="utf-8"))
    bracket = result.get("fitted_bracket_orbits")
    if (not isinstance(bracket, list) or len(bracket) != 2
            or not all(isinstance(x, (int, float)) for x in bracket)):
        raise RuntimeError(
            f"{RELAXATION_CEILING_ARTIFACT} carries no `fitted_bracket_orbits` "
            f"pair; regenerate it with `python {RELAXATION_CEILING_GENERATOR}`.")
    low, high = float(bracket[0]), float(bracket[1])
    if not 0.0 < low <= high:
        raise RuntimeError(
            f"{RELAXATION_CEILING_ARTIFACT} emits a relaxation bracket of "
            f"({low}, {high}), which is not an ordered pair of positive times")
    return (low, high)

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
    """The A-PRIORI span, for pricing a run that does not exist yet.

    IT OVERBUYS ON THIS WORLD AND THE MEASUREMENT SAYS BY HOW MUCH. The
    coverage criterion behind the multiple is sound and was taken on synthetic
    series whose answer is known in closed form; what it assumes is a tau that
    is a property of the PROCESS. This model's is a property of the WINDOW:
    measured on 140 clean orbits of `run_893e276ee029` it reads 1.00 at 20
    orbits, 1.35 at 40, 1.90 at 60 and 6.30 at 140, every one a supported
    estimate, while the lag-1 correlation stays between 0.27 and 0.52. So the
    rule does not close -- at 60 orbits it asks for 38 and at 140 for 126 --
    and each span bought raises the tau that prices the next.

    `production_span_from_report` is the operative rule wherever a run exists,
    and this stays for the case where none does.
    `exoplasim/notes/memory-time-and-the-production-span.md` has the tables.
    """
    return PRODUCTION_SPAN_TAU_MULTIPLE * float(tau_memory_orbits)


# The keys `assess_convergence.py` writes for "orbits THIS criterion needs".
# Named rather than pattern-matched over the block: `_if_independent` is the
# same number with the memory taken out and `_prices` is prose, and a rule that
# swept the block by prefix would take both.
_REQUIRED_WINDOW_KEYS = (
    "window_orbits_for_storage_criterion",
    "window_orbits_for_offset_criterion",
)


def production_span_from_report(report: dict) -> tuple[float, bool, str]:
    """The orbits a run's OWN statistics say its criteria need.

    Returns `(orbits, is_a_floor, which)`: the span, whether it is a lower
    bound rather than an answer, and the criterion that set it.

    THIS IS THE OPERATIVE RULE AND `production_span_orbits` IS NOT. A span
    priced through tau is priced through the standard error of a MEAN, and on
    this world that error has already converged: it is 0.0201 K at twenty
    orbits and 0.0255 at a hundred and forty, while the independent expectation
    halves. Seven times the orbits buy no better a mean, against criteria that
    discriminate at 0.15 K. What still tightens with the window is the SLOPE's
    error, which is what these keys are computed from, so reading them prices
    the span through the statistic that has not converged instead of the one
    that has.

    IT IS SELF-LIMITING, which is the property the tau route lacks. Each key is
    computed from the run's own scatter and its own memory time over its own
    window, so a run that has bought enough says so and a run that has not
    names the number it is short of.

    THEN ASSESS AT WHAT YOU BOUGHT, with `--window`. A block bought to this
    number can be SHORTER than `assess_convergence.py`'s default window, which
    is sized on the nominal bounds because it has to exist before any run is
    opened. The default then asks for orbits the block does not have, and
    `segments.py` refuses -- correctly, and it names the number to pass. That
    refusal is the two rules disagreeing, not a defect in either: the nominal
    is what to use when nothing is known about the run, and this is what to use
    once something is. Buying `max(this, DEFAULT_WINDOW_ORBITS)` avoids the
    friction and costs the difference in orbits, which is minutes at T21 and
    is not at T85.

    `is_a_floor` carries `required_window_is_a_lower_bound` through unchanged:
    a tau estimated inside the window it sizes is a lower bound on tau, so the
    span is a FLOOR until a span long enough to carry the estimate exists. A
    caller that treats a floor as an answer buys too few orbits, which is the
    one direction this must not fail in silently.
    """
    resolving = report.get("resolving_power") or {}
    wanted = {key: float(resolving[key]) for key in _REQUIRED_WINDOW_KEYS
              if isinstance(resolving.get(key), (int, float))}
    if not wanted:
        raise RuntimeError(
            "the report carries none of "
            f"{list(_REQUIRED_WINDOW_KEYS)} in its resolving_power block, so "
            "it cannot say how many orbits its own criteria need. It predates "
            "the block; re-run assess_convergence.py on the run.")
    which = max(wanted, key=wanted.__getitem__)
    return (wanted[which],
            bool(resolving.get("required_window_is_a_lower_bound", True)),
            which)


def commissioning_orbits(approach_orbits: float,
                         tau_memory_orbits: float) -> float:
    """The approach, then the production span that follows it.

    `approach_orbits` is operational experience and has no artifact behind it:
    about seventy at T21 from cold, ten to twenty on a reconvergence, and it
    moves with the rung. It is an ARGUMENT here rather than a constant so that
    the caller has to state which approach it means.
    """
    return float(approach_orbits) + production_span_orbits(tau_memory_orbits)


def commissioning_orbits_from_report(approach_orbits: float,
                                     report: dict) -> tuple[float, bool, str]:
    """The approach, then the span the run's own criteria ask for.

    The report-based sibling of `commissioning_orbits`, and the one to use
    wherever a run exists. Carries the floor flag out with it rather than
    resolving it here, because what a caller does about a floor is a decision:
    buy the orbits and re-read, or say the number is a lower bound.
    """
    span, is_a_floor, which = production_span_from_report(report)
    return float(approach_orbits) + span, is_a_floor, which


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
                     residual_k: float = SETTLING_RESIDUAL_K,
                     root=None) -> tuple[float, float]:
    """The settling length over the whole relaxation-time bracket.

    The bracket is READ and not declared; see `tau_relaxation_orbits_bracket`.
    """
    low, high = tau_relaxation_orbits_bracket(root)
    return (settling_orbits(perturbation_k, low, residual_k),
            settling_orbits(perturbation_k, high, residual_k))
