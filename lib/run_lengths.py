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

# THE MEMORY TIME, MEASURED and no longer inferred. The lower end is retained
# from the 2026-08-26 five-window Geyer sweep on run_432e5e46adef. The upper end
# follows run_0d41aa82c287, the T21 arm at dt 45 cut as the donor for the T42
# ladder comparison: its 52-orbit settled production window reads 3.861600
# orbits and explicitly reports enough span to resolve that estimate. Rounded
# upward, not to nearest, because this pair sizes a window and is an upper bound
# rather than a fit parameter.
#
# IT REPLACED 3.69, from run_5994d1f9624e's 46-orbit window on the same build and
# the same rung, and the two do not disagree about the model. The newer reading
# is taken over a LONGER settled window, and this model's tau grows with the
# window it is measured on -- the note below says so in as many words and the
# tables are in
# `exoplasim/notes/memory-time-and-the-production-span.md`. So the bound moved
# because a longer span was bought, which is the only direction it can move.
#
# THE TWO T42 ARMS READ WELL BELOW IT, 2.835763 and 1.177548 at the same step and
# build, so the finer rung is not what sets this. T21 is.
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
TAU_MEMORY_ORBITS_BRACKET = (1.89, 3.87)

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
    "run": "run_0d41aa82c287",
    "artifact": "exoplasim/analysis/convergence/run_0d41aa82c287_convergence.json",
    "node": ("resolving_power", "temperature_residual_tau_orbits"),
    "observation": 3.8615996294737345,
    "why_below": "the anchor is the T21 arm at dt 45 cut as the T42 ladder "
                 "comparison's donor, 108 orbits from cold on "
                 "canonical-10m-carve2, converged on all six criteria over a "
                 "52-orbit window, and the reported estimate has enough span to "
                 "resolve it. The top is the observation rounded upward to "
                 "3.87, which is what an "
                 "upper bound does when a supported reading passes it. WHAT "
                 "THE BRACKET IS NOT: a property of "
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
    """Is this convergence report a reading of the settled trajectory?

    THREE CONDITIONS, and the third was missing. The window has to be the
    run's own trajectory rather than a diagnostic, the span has to support the
    memory-time estimate it carries, and THE RUN HAS TO HAVE CONVERGED.

    The third is not a formality. A run still approaching its asymptote
    carries a residual trend, a residual trend pushes every lag correlation
    UP, and the memory time read off it is therefore an upper bound on a
    quantity the run does not yet have. The carved build's baseline is the
    case: it fails `extrapolated_offset` with 0.15 K still to go and reads tau
    3.14 against a bracket topping out at 2.45. Admitting that reading would
    raise a bound on every rung of the ladder on the strength of a run whose
    own verdict says it is not finished.

    The run's verdict is the discriminator and it already exists, so this
    reads it rather than inventing a second test. A report predating the field
    is admitted, because refusing every older report would empty the bound.
    """
    purpose = report.get("assessed_purpose")
    if purpose is None:
        purpose = "diagnostic" if "diagnostic" in name else "production"
    if purpose != "production":
        return False
    settled = report.get("sufficiently_equilibrated_for_worldbuilding")
    if settled is False:
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
    # THE OFFSET ROW CARRIES A QUALIFIER AND IT HAS TO BE READ. Its own
    # `..._prices` field says it prices the run's verdict ONLY where
    # `offset_statistic_source` is `drift_fallback`. On the exponential-fit
    # path the statistic is the FIT's half width and the number here is
    # computed from the expected relaxation time instead, so the two are
    # unrelated -- on the carved baseline they read 42.8 orbits against a
    # statistic six times its target, a factor of four and a half apart.
    #
    # Taking the number without its qualifier is how a rule that reads an
    # artifact goes wrong: the artifact was not lying, it was annotated, and
    # the annotation was ignored. Dropped rather than corrected, because the
    # window the fit path needs is not derivable from what this block carries.
    source = resolving.get("offset_statistic_source")
    if source is not None and source != "drift_fallback":
        wanted.pop("window_orbits_for_offset_criterion", None)
        unpriced = source
    if not wanted:
        if source is not None and source != "drift_fallback":
            raise RuntimeError(
                f"the offset criterion's window is priced for the drift "
                f"fallback and this run's statistic is `{source}`, so the "
                f"number does not price its verdict; nothing else in the "
                f"resolving_power block does either. The span cannot be read "
                f"off this report. What settles it is a run whose offset "
                f"statistic IS the drift fallback -- one settled enough that "
                f"the exponential has no approach left to fit.")
        raise RuntimeError(
            "the report carries none of "
            f"{list(_REQUIRED_WINDOW_KEYS)} in its resolving_power block, so "
            "it cannot say how many orbits its own criteria need. It predates "
            "the block; re-run assess_convergence.py on the run.")
    which = max(wanted, key=wanted.__getitem__)
    # DROPPING A CRITERION MAKES THE ANSWER A FLOOR, whatever the report says
    # about tau. What is left is the span the PRICEABLE criteria need, and the
    # one that was dropped may need more -- on the carved baseline it is the
    # criterion the run actually fails. A caller that took this as sufficient
    # would buy the orbits the storage criterion wants and still not have a
    # converged run.
    floor = bool(resolving.get("required_window_is_a_lower_bound", True))
    if source is not None and source != "drift_fallback":
        floor = True
    return (wanted[which], floor, which)


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


# =============================================================================
# THE ECOLOGICAL PAIR. A SECOND SET OF BOTH TIMES, FOR A DIFFERENT MODEL.
# =============================================================================
#
# READ THE TOP OF THIS FILE FIRST. Everything above is the CLIMATE model's, in
# ORBITS. Everything below is the BIOSPHERE model's, in COMPLETE FORCING CYCLES.
# They are four quantities, not two, and the ways they can be confused are worth
# naming because three of them look like arithmetic that works:
#
#   NOT THE SAME MODEL     `TAU_MEMORY_ORBITS_BRACKET` is how long ExoPlaSim's
#                          surface temperature stays correlated with itself.
#                          `ecological_timescale_brackets` is how long a
#                          gridcell's simulated vegetation does. They differ by
#                          more than an order of magnitude and neither prices
#                          the other's runs.
#   NOT THE SAME UNIT      a cycle is `forcing_cycle_years` simulation years and
#                          a simulation year is one modelled orbit, so cycles and
#                          orbits coincide NUMERICALLY only while the driver's
#                          cycle is one year. The run's own declared cycle length
#                          travels with the reading and is never assumed here.
#   NOT THE SAME SPAN      the climate side's settling block hands the next
#                          conversion a restart that is not mid-transient. The
#                          ecological spin-up precedes a RETAINED RECORD and has
#                          to leave that record clean. Same shape, different
#                          consumer, different residual.
#
# WHAT IS PRESERVED IS THE DISCIPLINE, not the names: memory sizes a span whose
# MEAN needs an interval, relaxation sizes a block for a transient to DECAY, and
# no answer here is anything but a bracket.
#
# THIS SIDE STATES NO NUMBER AT ALL. Both ecological times are read out of the
# acceptance artifacts `biosphere/scripts/assess_lpj_run.py` writes, which record
# them on a refusal as well as on a pass -- a run refused for not having settled
# is exactly the run whose timescales say how long the next one must be. That is
# the `tau_relaxation` treatment rather than the `TAU_MEMORY` one, and it is the
# stronger of the two: there is nothing here to go stale against the measurement.
ECOLOGICAL_ACCEPTANCE_DIR = "biosphere/analysis"
ECOLOGICAL_TIMESCALE_GENERATOR = "biosphere/scripts/assess_lpj_run.py"


def ecological_timescale_brackets(root=None) -> dict:
    """Both ecological times, read from every acceptance artifact that has them.

    The MEMORY bracket spans the fields, because a retained record has to
    establish the memory time of every field it will be asked to judge, so the
    top is what sizes it.

    The RELAXATION bracket is bounded BELOW by the slowest field whose approach
    is measurable and is OPEN ABOVE whenever any field's approach is not, because
    a field the estimator declines is a field that may be slower still. An open
    top is not a missing number; it is the honest shape of a timescale longer than
    the record that measured it, and it is why what comes out of it is a FLOOR.
    It may be MISSING ENTIRELY, and on the records this project has it is: no
    field's contraction is resolvably below one, so `relaxation_cycles_bracket`
    is None and `relaxation_measurable` says so. That does not touch the memory
    or resolving brackets, because the two lengths are independent -- the record
    is sized by what the acceptance contract must be able to see and the spin-up
    by how long a transient takes to decay -- and one being unmeasurable must not
    hide the other.

    A MEMORY READING THE ESTIMATOR CALLS UNRELIABLE IS STILL ADMITTED, and that is
    deliberate rather than an oversight. `reliable` false means the record is
    shorter than ten times the tau it just produced, which makes that tau a LOWER
    bound on a large number. Dropping it would shorten the record floor on the
    strength of an estimate the estimator has written down as too small, which is
    the same move `TAU_MEMORY_ORBITS_BRACKET` refuses on the climate side. Only
    the TOP of each bracket sizes anything here, so a short run pooled with a long
    one can widen the bottom and cannot shorten a floor.

    AN ARTIFACT FROM A SUPERSEDED CONTRACT IS NOT POOLED, and that is the wire
    that makes the assessed set a live declaration rather than a number frozen
    into every artifact ever written. What a `record_cycles_for_bound` means is
    "the record at which THIS quantity closes on THAT tolerance", and both halves
    are the contract's. Contract 4 assessed every column of every stability table
    at one tolerance and its artifacts ask for 33846 cycles; contract 5 assesses
    the quantities a consumer reads at the tolerance each owes, and asks for the
    record already on disk. Pooling the two would size every future run from the
    contract that no longer judges it. Superseded artifacts are counted and named
    rather than dropped in silence, and when none is current this raises with the
    command that re-takes them.
    """
    base = Path(root) if root is not None else _REPO_ROOT
    directory = base / ECOLOGICAL_ACCEPTANCE_DIR
    memory, relaxation, sources, cycle_years = [], [], [], set()
    declined = 0
    tolerance, span = None, []
    resolving, unresolvable = [], 0
    from lpj_output import read_policy      # local: this module stays light
    contract = read_policy()["contract_version"]
    superseded = []
    for path in sorted(directory.glob("lpj_*/acceptance.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        scales = report.get("timescales") or {}
        tables = scales.get("tables")
        if not tables:
            continue
        if scales.get("contract_version") != contract:
            superseded.append(
                f"{path.relative_to(base)} ({scales.get('contract_version')})")
            continue
        sources.append(str(path.relative_to(base)))
        if tolerance is None and scales.get("drift_tolerance") is not None:
            tolerance = float(scales["drift_tolerance"])
        for table in tables.values():
            cycle_years.add(int(table.get("forcing_cycle_years", 1)))
            span.append(int(table.get("record_cycles", 0)))
            for field in table.get("fields", {}).values():
                memory.append(float(field["memory"]["tau_cycles"]))
                approach = field.get("relaxation", {})
                if approach.get("admissible"):
                    relaxation.append(float(approach["tau_cycles"]))
                else:
                    declined += 1
                needed = (field.get("drift") or {}).get("record_cycles_for_bound")
                if needed is None:
                    continue
                if math.isfinite(float(needed)):
                    resolving.append(float(needed))
                else:
                    unresolvable += 1
    if not memory:
        stale = (f" {len(superseded)} artifact(s) carry a superseded contract "
                 f"and were not pooled: {', '.join(superseded)}."
                 if superseded else "")
        raise RuntimeError(
            f"no acceptance artifact under {ECOLOGICAL_ACCEPTANCE_DIR} carries "
            f"ecological timescales under {contract}, and they are where this "
            f"bracket lives.{stale} Run "
            f"`python {ECOLOGICAL_TIMESCALE_GENERATOR} <run>`; it records them "
            "whatever the verdict.")
    if len(cycle_years) != 1:
        raise RuntimeError(
            f"the acceptance artifacts declare more than one forcing cycle "
            f"length {sorted(cycle_years)}, so their cycles are not one unit and "
            "cannot be bracketed together")
    if not resolving:
        raise RuntimeError(
            "no acceptance artifact carries a per-field resolving length. The "
            "retained record floor is `lib/lpj_output.py:cycles_for_bound`'s "
            "answer for the field that needs the most, and it is written onto "
            f"the artifact by `{ECOLOGICAL_TIMESCALE_GENERATOR}`. Re-assess a "
            "run under the current contract rather than falling back to a span "
            "multiple, which is the rule this one replaced.")
    return {
        "memory_cycles_bracket": (min(memory), max(memory)),
        "relaxation_cycles_bracket": ((min(relaxation), max(relaxation))
                                      if relaxation else None),
        "relaxation_measurable": bool(relaxation),
        "relaxation_top_is_open": declined > 0,
        "fields_with_a_measurable_approach": len(relaxation),
        "fields_declined": declined,
        # The record each field needs before its drift bound falls inside the
        # contract's tolerance, as the producer computed it. The TOP sizes the
        # record because it has to serve every field it will be asked to judge,
        # and `unresolvable` counts the fields for which no finite record does --
        # a bracket with an open top, reported and not silently dropped.
        "resolving_cycles_bracket": (min(resolving), max(resolving)),
        "fields_with_no_finite_record": unresolvable,
        "forcing_cycle_years": cycle_years.pop(),
        "longest_record_cycles": max(span) if span else 0,
        "drift_tolerance_when_assessed": tolerance,
        "contract_version": contract,
        "sources": sources,
        # Named rather than dropped in silence: an artifact taken under a
        # superseded contract measured a different assessed set at a different
        # tolerance, so it cannot size a run this contract will judge.
        "superseded_sources": superseded,
    }


def ecological_drift_tolerance() -> float:
    """The residual a spin-up has to decay to, from the contract that judges it.

    `SETTLING_RESIDUAL_K` takes its value from the criterion that judges the next
    state, deliberately, and this is the same move: the drift the acceptance
    contract will absorb over a retained record is what a spin-up has to leave
    below. It is read LIVE rather than from the artifact a past run recorded,
    because the floor is for the run that has not happened yet -- tighten the
    contract and the spin-up it demands lengthens, with nothing to go stale in
    between. `lib/lpj_output.py:read_policy` is the one reader of that contract
    and is imported here rather than the YAML being parsed a second time.
    """
    from lpj_output import read_policy      # local: this module stays light
    return float(read_policy()["trend"]["relative_end_to_end_limit"])


def ecological_record_cycles(brackets: dict) -> float:
    """How much record a run must RETAIN, from what its own acceptance resolves.

    THE RECORD IS SIZED BY WHAT THE CONTRACT MUST BE ABLE TO SEE, not by a
    multiple of anything. `lib/lpj_output.py:cycles_for_bound` inverts the
    contract's own drift bound: the record at which a settled field of this
    scatter and this memory time bounds its end-to-end drift inside
    `relative_end_to_end_limit`. The memory time enters through the standard
    error of a half-record mean, which is where it belongs, and the top of the
    bracket sizes the record because a record has to serve every field it will
    be asked to judge.

    WHY IT IS NOT A SPAN MULTIPLE ANY MORE. `RELIABLE_SPAN_MULTIPLE * tau` asks
    for a record ten times a memory time read off that same record, and the
    memory time this model reports GROWS with the window it is read on: 9.3 to
    125.3 cycles on a 1000-cycle record became 11.1 to 212.7 on the 1253-cycle
    record those numbers sized. It did not converge, twice. Here the standard
    error falls as one over the root of the record while the memory time grows
    sublinearly with it, so the requirement is reached rather than chased.

    IT IS STILL A FLOOR, for the weaker reason: the memory time is held at its
    measured value while a longer record may read a larger one, which moves the
    answer rather than preventing one.

    POOLING ACROSS RUNS OF DIFFERENT LENGTHS IS SAFE HERE, and the arithmetic
    says why rather than a convention: the resolving length is proportional to a
    field's variance times its memory time and does NOT otherwise depend on the
    record it was read on, so it is a property of the field. A record too short
    to see a field's memory time underestimates it and so underestimates the
    length, which can only widen the bottom of the bracket. Only the top sizes
    anything.
    """
    return float(brackets["resolving_cycles_bracket"][1])


def ecological_spinup_cycles(tau_relaxation_cycles: float,
                             record_cycles: float,
                             drift_tolerance: float) -> float:
    """How much spin-up must precede a record for that record to be clean.

    The same shape as `settling_orbits` and for the same reason, with the residual
    taken from the criterion that judges what follows rather than typed. A
    spin-up starting from bare ground begins a full equilibrium level away, so the
    approach remaining after `S` cycles is `exp(-S / tau)` of the level, and the
    DRIFT it puts across a retained record of `L` cycles is that times
    `1 - exp(-L / tau)`. The acceptance contract refuses a record whose relative
    end-to-end change exceeds its `relative_end_to_end_limit`, so

        exp(-S / tau) * (1 - exp(-L / tau)) <= tolerance

    and this returns the smallest `S` that satisfies it. Nothing in it is typed:
    the tolerance is the contract's own, read from the artifact.

    STARTING FROM BARE GROUND IS THE CONSERVATIVE READING and is stated rather
    than hidden. A pool the CENTURY accelerator hands over part-grown begins
    closer than a full level away, so its own requirement is shorter than this.
    """
    tau = float(tau_relaxation_cycles)
    share = 1.0 - math.exp(-float(record_cycles) / tau)
    if share <= 0.0:
        return 1.0
    remaining = float(drift_tolerance) / share
    if remaining >= 1.0:
        return 1.0
    return max(1.0, tau * math.log(1.0 / remaining))


def ecological_spinup_multiple(drift_tolerance: float) -> tuple[float, float]:
    """The spin-up, in retained records, that works at EVERY relaxation time.

    THIS IS WHY THE RELAXATION TIME NEVER HAS TO BE MEASURED. The requirement
    `ecological_spinup_cycles` inverts is a function of a `tau` no record this
    project has can resolve, and the four dispositions a number without a
    derivation has do not all apply here: the requirement is BOUNDED OVER ALL
    `tau`, so the bound itself is the derivation. Substituting `u = L / tau`,

        S(tau) = L * ln((1 - exp(-u)) / tolerance) / u

    falls away at both ends -- an approach far slower than the record puts
    almost none of itself into the record, and one far faster has finished --
    so it has a finite maximum. The left-hand side of the requirement is
    decreasing in `S`, so any `S` at or above that maximum satisfies it at every
    `tau`. It is the smallest such `S`, which makes it the minimax answer rather
    than a guess at a length.

    It is exactly proportional to the retained record and depends on nothing
    else, which is what makes it a multiple:

        C(t) = max over u of ln((1 - exp(-u)) / t) / u

    Returns `(C, u)`, the multiple and the `L / tau` it is attained at, so a
    caller can check the requirement AT the worst case rather than near it. A
    MEASURED relaxation time can only ask for less than this, so it is a
    CEILING on the floor: `ecological_run_cycles` prefers a measured value where
    one is admissible and falls back here, never on a convention.
    """
    tolerance = float(drift_tolerance)
    if not 0.0 < tolerance < 1.0:
        raise ValueError("the drift tolerance must lie in (0, 1)")

    def value(u: float) -> float:
        share = 1.0 - math.exp(-u)
        if share <= tolerance:
            return -math.inf
        return math.log(share / tolerance) / u

    # Golden section on a unimodal function. The bracket has to hold the whole
    # answer: the maximum sits at u below one for every tolerance this contract
    # can carry -- 0.136 at 0.05 and 0.054 at 0.02 -- and moves toward zero as
    # the tolerance tightens, so the lower end is what has to be small.
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    low, high = 1.0e-9, 50.0
    left, right = high - ratio * (high - low), low + ratio * (high - low)
    for _ in range(400):
        if value(left) > value(right):
            high = right
        else:
            low = left
        left, right = high - ratio * (high - low), low + ratio * (high - low)
    u = 0.5 * (low + high)
    return value(u), u


def ecological_spinup_ceiling(record_cycles: float,
                              drift_tolerance: float) -> float:
    """The spin-up a record of this length needs whatever the relaxation time is."""
    multiple, _ = ecological_spinup_multiple(drift_tolerance)
    return max(1.0, multiple * float(record_cycles))


def ecological_run_cycles(root=None) -> dict:
    """The spin-up and retained record this world's ecology needs, and why.

    THE TWO ARE INDEPENDENT AND ARE RETURNED INDEPENDENTLY. The record's floor is
    what the acceptance contract has to be able to RESOLVE, inverted out of its
    own drift bound. The spin-up's is how long a TRANSIENT takes to decay.

    A SPIN-UP IS ALWAYS DERIVED HERE, and on the records this project has it is
    derived WITHOUT a relaxation time. No field's contraction is resolvably below
    one, so `relaxation_time` declines every field-record and there is no
    measured `tau` to invert the requirement at. There does not have to be:
    `ecological_spinup_multiple` shows the requirement is bounded over all `tau`,
    so its supremum satisfies it at every `tau` and is the smallest length that
    does. `spinup_basis` names which of the two the returned number is, and the
    minimax is a CEILING on the floor -- a measured relaxation time can only ask
    for less.

    BOTH ARE FLOORS, and the spin-up's floor is what a run buys UP FRONT while
    the acceptance contract's refusal buys the rest. That is the same bargain the
    climate side makes, and it is why the spin-up is not a correctness
    requirement: a residual transient enters the very statistic the contract
    gates on, so an under-spun run is REFUSED rather than silently accepted.
    """
    brackets = ecological_timescale_brackets(root)
    tolerance = ecological_drift_tolerance()
    record = ecological_record_cycles(brackets)
    if brackets["relaxation_measurable"]:
        spinup = ecological_spinup_cycles(
            brackets["relaxation_cycles_bracket"][1], record, tolerance)
        basis = "lib/run_lengths.py:ecological_spinup_cycles"
        reason = (
            "inverted at the top of the measured relaxation bracket, "
            f"{brackets['relaxation_cycles_bracket'][1]:.1f} cycles, over the "
            f"{record:.0f}-cycle record floor at a drift tolerance of "
            f"{tolerance:g}")
    else:
        multiple, worst = ecological_spinup_multiple(tolerance)
        spinup = ecological_spinup_ceiling(record, tolerance)
        basis = "lib/run_lengths.py:ecological_spinup_multiple"
        reason = (
            "no field's approach to equilibrium is measurable on the records "
            f"available: all {brackets['fields_declined']} field-records were "
            "declined by `lib/lpj_output.py:relaxation_time`, most for block "
            "differences that change sign or sit inside their own standard "
            "error, and every one that reached a contraction had an uncertainty "
            "reaching one. None is needed: the requirement is bounded over all "
            f"relaxation times, and its supremum is {multiple:.5f} times the "
            f"retained record, attained at a relaxation time of {1.0 / worst:.4f} "
            "records. That is the smallest spin-up that satisfies the "
            "requirement whatever the relaxation time is, so it is derived "
            "rather than defaulted, and a measured relaxation time could only "
            "ask for less.")
    return {
        "spinup_cycles": spinup,
        "spinup_basis": basis,
        "spinup_reason": reason,
        "spinup_is_minimax": not brackets["relaxation_measurable"],
        "record_cycles": record,
        "total_cycles": spinup + record,
        "is_a_floor": True,
        "floor_because": (
            "the memory time each resolving length was computed at is the one "
            "its own record read, and a longer record may read a larger one"
            + ("; and the relaxation bracket's top is open, because "
               f"{brackets['fields_declined']} fields' approach is not "
               "measurable on the records available and a field the estimator "
               "declines may be slower still"
               if brackets["relaxation_measurable"]
               and brackets["relaxation_top_is_open"] else "")),
        "brackets": brackets,
    }
