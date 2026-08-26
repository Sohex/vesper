"""The ExoPlaSim resolution ladder: the one place it is written down.

Worldbuilding frame: the ladder is the set of Gaussian grids the Vesper
climate model can be solved on. Nothing here is about the real world.

WHY THIS IS ITS OWN MODULE and not part of `lib/gridding.py`: the build script
and two shell probes need it, and neither should have to import numpy and the
World Orogen export reader to find out how many latitudes T85 has.
`lib/gridding.py` re-exports everything below, so there is one definition and
two doors rather than two definitions.

WHAT IT REPLACED. `config/planet.yaml` stated `resolution`, `latitudes` and
`longitudes` independently, so a valid truncation could be paired with another
grid's dimensions and nothing said so; `build_model.py`, the SHTns variant
sweep, `verify_shtns_equivalence.sh` and `run_shtns_probe.sh` each carried
their own copy of the mapping, and two of the four were missing rungs the
others had. SPAT-2.
"""
from __future__ import annotations


# THE TABLE IS SELF-CHECKING rather than merely agreed: PlaSim derives the
# truncation from the grid as `(NLON-1)/3` with `NLON = 2*NLAT`, so a rung's
# name and its latitude count are not two facts but one, and `_check_rungs()`
# below rejects any row where they disagree. That is a right answer, not a
# convention.
RUNGS = {"T21": 32, "T31": 48, "T42": 64, "T63": 96,
         "T85": 128, "T106": 160, "T127": 192, "T170": 256}


def _check_rungs() -> None:
    for rung, nlat in RUNGS.items():
        if (2 * nlat - 1) // 3 != int(rung[1:]):
            raise RuntimeError(
                f"the rung table is inconsistent: {rung} claims {nlat} "
                f"latitudes, and a {2 * nlat}-column Gaussian grid truncates "
                f"at T{(2 * nlat - 1) // 3}")


_check_rungs()


def geometry(rung: str) -> tuple[int, int, int]:
    """(latitudes, longitudes, truncation) for a ladder rung.

    The only accepted spellings are the table's. A rung outside it must be an
    ERROR rather than a default: `-r 170` once built T21 and said nothing.
    """
    key = str(rung).upper()
    if key not in RUNGS:
        raise RuntimeError(
            f"{rung!r} is not a ladder rung. The ladder is "
            f"{', '.join(RUNGS)}; anything else needs a row here and a grid "
            "in the export before it can be run.")
    nlat = RUNGS[key]
    return nlat, 2 * nlat, int(key[1:])


def rung_of_latitudes(nlat: int) -> str:
    """The rung with this many latitudes, or an error naming what is close."""
    for rung, n in RUNGS.items():
        if n == int(nlat):
            return rung
    raise RuntimeError(
        f"no ladder rung has {nlat} latitudes; the ladder is "
        + ", ".join(f"{r} ({n})" for r, n in RUNGS.items()))


def model_grid(config: dict) -> tuple[str, int, int]:
    """The configured (rung, latitudes, longitudes), refusing a stale pairing.

    `config/planet.yaml` carries all three because `read_sra` validates every
    staged surface file against `latitudes` and `longitudes`, so a resolution
    changed without them refuses its own inputs. That makes the three a single
    fact written three times, and this is where it is checked rather than
    trusted -- the same shape as gravity and mass, which `run_exoplasim.py`
    refuses when they disagree.
    """
    model = config["model"]
    rung = str(model["resolution"]).upper()
    nlat, nlon, _ = geometry(rung)
    declared = int(model["latitudes"]), int(model["longitudes"])
    if declared != (nlat, nlon):
        raise RuntimeError(
            f"config/planet.yaml sets model.resolution {rung} with "
            f"latitudes {declared[0]} and longitudes {declared[1]}, and {rung} "
            f"is {nlat} by {nlon}. The three are one fact: fix whichever is "
            "stale rather than leaving a truncation paired with another "
            "grid's dimensions.")
    return rung, nlat, nlon


def fft_module(rung: str) -> str:
    """The Fortran FFT module that can transform this rung's longitudes.

    `fftmod` is the radix 8-4-3-2 transform: `gp2fc` does one pass of radix 8,
    then radix 4 while four or more remain, then a single radix 3 or radix 2
    tail. It therefore transforms `8 * 4**k * r` with `r` in 1, 2, 3, and
    NOTHING ELSE -- its `nallowed` table is that set, not an independent fact.
    `fft991mod` is the FFT991 package, whose `set99` factorises with 8, 6, 5,
    4, 3 and 2 and covers the lengths `fftmod` cannot.

    DERIVED RATHER THAN LISTED, because a table is what went wrong: the
    postprocessor named the two longitude counts that need `fft991mod` and
    then tested them against a LATITUDE, so T63 and T106 both selected the
    module that cannot transform them and the extension's bare Fortran `stop`
    killed the interpreter with no traceback. world-i38.

    `vendor/exoplasim/exoplasim/pyburn.py:_fftmodule` states the same rule for
    the postprocessor's f2py extensions. That package stays importable on its
    own and cannot import this one, so the rule is written twice on purpose;
    each names the other.
    """
    _, nlon, _ = geometry(rung)
    return "fftmod" if _fftmod_can_transform(nlon) else "fft991mod"


def _fftmod_can_transform(nlon: int) -> bool:
    """Whether `fftmod`'s radix loop covers `nlon` columns. See `fft_module`."""
    n = int(nlon)
    if n % 8:
        return False
    r = n // 8
    while r % 4 == 0:
        r //= 4
    return r in (1, 2, 3)


# THE TWO RESTATEMENTS THAT CANNOT BE REMOVED. The model's build system and the
# vendored ExoPlaSim package each need the ladder before this module is
# reachable: CMake has no Python, and `vendor/exoplasim/exoplasim/__init__.py`
# has to stay importable as an installed package that knows nothing about this
# repository. So the ladder is written three times on purpose, and the other two
# are CHECKED against this one rather than trusted. T31 is why: CMake's list
# omitted 48 while this table declared T31, so `build_model.py --res T31`
# resolved cleanly and then died at cmake configure, and the package's dispatch
# chain accepted seven rungs where the ladder declares eight. world-ajx.
RESTATEMENTS = (
    ("vendor/exoplasim/exoplasim/plasim/CMakeLists.txt",
     r'require_one_of\(PLASIM_NLAT\s+"\$\{PLASIM_NLAT\}"\s+([0-9 ]+)\)'),
    ("vendor/exoplasim/exoplasim/__init__.py",
     r'RESOLUTIONS = (\{[^}]*\})'),
    # A COMMENT, and checked for the same reason as the two above: it is what
    # someone editing `NLAT_ATM` by hand reads, and it named itself as checked
    # by this function while sitting outside the table.
    ("vendor/exoplasim/exoplasim/plasim/src/plasimmod.f90",
     r'checked by rungs\.check_restatements\. NLAT by rung:\n!\s*([0-9, ]+)\n'),
)


def check_restatements(root) -> list[str]:
    """Every place the ladder is restated, against this table. Empty when they agree.

    A CHECK WITH A RIGHT ANSWER rather than a convention: each restatement has
    to name exactly the latitude counts in `RUNGS`, and any disagreement is
    reported as the two sets. Takes the repository root rather than resolving
    one, so this module keeps knowing nothing but the ladder.
    """
    import ast
    import re
    from pathlib import Path

    want = set(RUNGS.values())
    problems = []
    for rel, pattern in RESTATEMENTS:
        path = Path(root) / rel
        if not path.is_file():
            problems.append(f"{rel} is not there, and it restates the ladder")
            continue
        m = re.search(pattern, path.read_text(encoding="utf-8"))
        if m is None:
            problems.append(
                f"{rel} no longer carries the ladder in the shape this check "
                f"reads ({pattern!r}); it restates the ladder and cannot go "
                f"unchecked")
            continue
        body = m.group(1)
        if body.lstrip().startswith("{"):
            got = set(ast.literal_eval(body).values())
        else:
            # Whitespace or comma separated: CMake's list is spaced and the
            # Fortran comment's is a comma-separated row.
            got = {int(t) for t in body.replace(",", " ").split()}
        if got != want:
            problems.append(
                f"{rel} allows latitudes {sorted(got)} and the ladder is "
                f"{sorted(want)}: missing {sorted(want - got)}, extra "
                f"{sorted(got - want)}")
    return problems


# ---------------------------------------------------------------------------
# The timestep: THREE QUANTITIES, kept apart
# ---------------------------------------------------------------------------
#
# "The timestep at rung X" has named three different things in this tree, and
# reconciling them to one number is the wrong move: each is load-bearing
# somewhere and they answer different questions. WORLD-F997, WORLD-J37.
#
#   STABILITY_CEILING_MINUTES  the coarsest step the rung has been MEASURED to
#                              start clean at. A property of the rung and of
#                              `model.filter_kappa`, not a preference.
#   ESCALATION_ROUTE           the (rung, step) sequence the project RUNS, in
#                              order. A decision, and the only one of the three
#                              that says what a run is configured to.
#   COMMISSIONING_EVIDENCE     what a run of COMMISSIONING LENGTH has shown
#                              about one (rung, step) pair. A probe cannot say
#                              this and the ceiling table does not carry it.
#
# A fourth quantity used to live here -- every rung at ONE step, so the ladder's
# rungs would differ only by support -- and is deliberately NOT carried. SPAT-8,
# the between-rung convergence comparison, was its only consumer and closed when
# T85 was declared the operating support. `notes/audits/resolution-ladder.md`
# keeps the argument and records that it has no consumer.
#
# `config/planet.yaml` carries `model.timestep_minutes`, the ACTIVE step one run
# is configured at, and nothing else about the timestep: the ceiling table used
# to sit beside it and drifted from the grid it was copied out of.

# MEASURED. The coarsest step at which the 900-step probe starts and completes,
# at kappa 8 and `--tau-scale 1`, from the grid in
# `exoplasim/notes/physics-filter-stability.md` and the per-cell verdicts in
# `exoplasim/analysis/stability_probe.json`. `check_stability_ceilings()` below
# re-derives these from that artifact rather than trusting this copy.
#
# EVERY VALUE HERE IS FROM THE 2026-08-26 RE-TAKE (WORLD-37TN) and none is a
# tested-not-bounded figure any more: the grid now carries columns above each
# ceiling, so each of these is a real boundary with a refusing cell above it.
# T21 refuses at 150, T42 at 90, T85 at 60. T85 had never been probed on any
# earlier source at all.
#
# The re-take also invalidates every grid before it, and not only on
# attribution: NHDIFF was never written by the probe, so every cell above T21 of
# every prior grid was damped from T21's wavenumber.
#
# PROVISIONAL, in two ways that are both about what a probe can say. It is a
# FLOOR: 900 steps is about a seventh of an orbit at dt 45, so a rung marked
# clean here is qualified against REFUSAL and against nothing else, and
# `COMMISSIONING_EVIDENCE` is where endurance lives. And every cell was taken at
# `filter_power` 8 where the model runs 16, which WORLD-37TN re-measures.
#
# T31, T63 and T106 are on the ladder and have never been probed, so they have
# no entry: `stability_ceiling` raises for them rather than interpolating, which
# is what the 1/N fit would invite and what T170 is the standing warning about.
STABILITY_CEILING_MINUTES = {
    "T21": 120.0,
    "T42": 60.0,
    "T85": 45.0,
    "T127": 30.0,
    "T170": 15.0,
}

# The artifact the ceilings are read out of, and the two conditions that select
# the cells that count. A cell at another kappa or on inherited damping is a
# boundary for a different model.
#
# `tau_scale` is a multiplier on the DERIVED table in `config/planet.yaml`, and
# that table moved when `model.hyperdiffusion.eddy_wind_m_s` was re-measured:
# every timescale is now 1.699x shorter, so a cell probed at `tau_scale` 1.0
# then and one probed at 1.0 now are not the same damping. The cells stand
# anyway, and as FLOORS rather than as boundaries that need re-taking, because
# the change is MORE damping and the operator that applies it is implicit in
# the damped variable (`plasim.f90:5148`), which is unconditionally stable and
# cannot bring a refusal forward. `exoplasim/notes/resolution-tuned-parameters.md`
# carries the measurement.
STABILITY_GRID = "exoplasim/analysis/stability_probe.json"
STABILITY_GRID_KAPPA = 8.0
STABILITY_GRID_TAU_SCALE = 1.0
# The probe's own three verdicts. `refused_only_at_length` is the one that
# matters here and the reason a two-verdict reading of the grid goes wrong: the
# cell STARTED and then died inside the probe's own range, which is neither a
# refusal nor a clean run. T170 at dt 22.5 is that cell, and it was carried as
# T170's adopted step for as long as the grid was read as refuses-or-runs.
# WORLD-6QRR.
PROBE_CLEAN = "no_refusal_in_steps"

# A DECISION, and the user's: `docs/src/pipeline/sequencing.md` section D is
# authoritative for it and this is the machine-readable restatement, checked
# against that chapter by `check_timestep_restatements`.
#
# In ORDER. Each entry changes the rung or the step and never both, so exactly
# one variable moves at a time and a surprise after a conversion is
# attributable. The invariant that falls out of it -- and the one WORLD-FL9C is
# about -- is that EVERY CHANGE OF RUNG HAPPENS AT CONSTANT dt. `_check_route()`
# is where that is a right answer rather than a description.
#
# THIS ROUTE NEVER CHANGES THE STEP AT ALL, which is what makes it three entries
# rather than four. dt 45 is at or below the measured ceiling of every rung on
# it, so no conversion needs a step change in front of it, and a settling block
# exists only to make one. It is an ATTEMPT: the route the project runs first,
# to find out whether the ladder holds at one step, and the fallback if a rung
# will not take 45 is the four-entry form that reconverges T21 at 30 and runs
# both conversions there.
ESCALATION_ROUTE = (
    ("T21", 45.0),
    ("T42", 45.0),
    ("T85", 45.0),
)

# WHAT A COMMISSIONING-LENGTH RUN HAS SHOWN, per (rung, step). The probe grid
# qualifies a step against refusal and nothing else, and the two are different
# verdicts: T42 at dt 45 completes every probe arm and is in this table as a
# late blow-up. Only pairs an actual run has reached appear; a pair with no row
# has no endurance evidence, which is a different thing from passing.
#
# `orbits` is the last orbit the run reached. `verdict` is `endured` for a run
# that reached a commissioning length without failing, `blew_up` for one that
# died after starting clean.
#
# `binds` is whether the row may REFUSE a configuration on the source this tree
# has now, and it defaults to True: a blow-up bars the pair until something is
# argued about it. Setting it False takes an argument in the row, and the only
# argument that works is that the measurement cannot be checked -- the run, its
# provenance and its reproducer all gone, so there is nothing to re-read and
# nothing to re-run. A row that does not bind is still EVIDENCE and is still
# reported: `commissioning_caveats` is where it comes out, and
# `run_exoplasim.py` prints it at launch, so an attempt at that pair knows what
# it is attempting.
COMMISSIONING_EVIDENCE = {
    ("T21", 45.0): {
        "verdict": "endured",
        "orbits": 50,
        "run": "run_ec32946bec89",
        "detail": "The route's first rung. ENDURED is the verdict here and it "
                  "is not the convergence verdict: under the upper-bound "
                  "criterion form this run MISSES on state storage, 0.0961 "
                  "plus a standard error of 0.0303 against a threshold of "
                  "0.12. It met all six under the point-estimate form the "
                  "criteria carried when it was taken. WORLD-S8N3.",
    },
    ("T21", 30.0): {
        "verdict": "endured",
        "orbits": 35,
        "run": "run_14906cb7b914",
        "detail": "The route's second rung, seeded from run_ec32946bec89. "
                  "Five of six criteria met; the extrapolated offset misses. "
                  "Its equilibrium sits about 0.46 K below the dt 45 one at "
                  "the same resolution, which is a timestep-dependent shift at "
                  "fixed support. This is a SETTLING step rather than a "
                  "commissioning one and the offset criterion is a "
                  "commissioning claim nothing between rungs reads: "
                  "WORLD-NDBQ. Its fitted relaxation is degenerate and "
                  "supports nothing, WORLD-5GQY.",
    },
    ("T42", 45.0): {
        "verdict": "blew_up",
        "orbits": 46,
        "run": "run_900548ae632e",
        "binds": False,
        "detail": "SIGFPE inside the 47th orbit on a gridpoint at -12.81 K at "
                  "the second level from the top, after 46 orbits of ordinary "
                  "climate with no trend towards it. WORLD-TD3; "
                  "notes/audits/resolution-ladder.md and "
                  "exoplasim/notes/physics-filter-stability.md.",
        "does_not_bind": "The run is gone and so is everything that could "
                         "re-check it: run_900548ae632e is not in "
                         "exoplasim/runs/, not a stub under archive/runs/ and "
                         "not in any INDEX_AT_DELETION.json, so its source sha "
                         "cannot be read and the ninety-second reproducer "
                         "WORLD-TD3 rests on cannot be re-run. What it was "
                         "measured on is therefore unknown, and is certainly "
                         "not this source: the damping correction, the epilog "
                         "use-after-free and the batch-2 forcing terms all "
                         "landed after it. Refusing a step on a profile the "
                         "model no longer has is failure-modes class 34, and a "
                         "claim with no artifact to check it against cannot "
                         "refuse one. The reopen condition is WORLD-TD3's and "
                         "is unchanged: if it recurs, this row binds again on "
                         "a run that exists and the ladder falls back to dt "
                         "30.",
    },
    ("T42", 30.0): {
        "verdict": "endured",
        "orbits": 84,
        "run": "run_1d39fef9bfc2",
        "detail": "Converged on the corrected damping. Its state does not "
                  "CONVERT cleanly to T85 across a change of step, which is "
                  "WORLD-FL9C and is a property of the conversion rather than "
                  "of this run.",
    },
}


def _check_ceilings() -> None:
    """Every ceiling names a rung on the ladder, and is a positive step."""
    for rung, dt in STABILITY_CEILING_MINUTES.items():
        if rung not in RUNGS:
            raise RuntimeError(
                f"{rung} has a stability ceiling and is not a ladder rung; the "
                f"ladder is {', '.join(RUNGS)}")
        if not (dt > 0.0):
            raise RuntimeError(f"{rung}'s stability ceiling is {dt}")


def _check_route() -> None:
    """The route's own invariants, each with a right answer.

    Three, and the third is what WORLD-FL9C asks of the route:

    1. Every entry is a ladder rung with a measured ceiling.
    2. Consecutive entries change EXACTLY ONE of (rung, step). Both at once is
       a state that cannot say which moved it; neither is not a step.
    3. Every change of rung happens at CONSTANT dt. A conversion across a
       change of step reinterprets the donor's two stored leapfrog levels as
       spanning the target's step and copies `nstep` -- which is elapsed time
       divided by the step -- so the converted run's calendar and stellar phase
       move by the ratio. The route is built to never ask for one.
    """
    for rung, dt in ESCALATION_ROUTE:
        if rung not in RUNGS:
            raise RuntimeError(f"the escalation route names {rung}, which is "
                               f"not a ladder rung")
        ceiling = STABILITY_CEILING_MINUTES.get(rung)
        if ceiling is None:
            raise RuntimeError(
                f"the escalation route runs {rung} at dt {dt} and {rung} has "
                "no measured stability ceiling. A route step on an unprobed "
                "rung is a step nothing has shown the rung can take.")
        if dt > ceiling:
            raise RuntimeError(
                f"the escalation route runs {rung} at dt {dt} and the coarsest "
                f"step {rung} is measured to start clean at is {ceiling}. The "
                "ceiling bounds the route; it does not set it.")
    for (r0, d0), (r1, d1) in zip(ESCALATION_ROUTE, ESCALATION_ROUTE[1:]):
        moved = (r0 != r1) + (d0 != d1)
        if moved != 1:
            raise RuntimeError(
                f"the escalation route goes {r0} at dt {d0} to {r1} at dt {d1}, "
                + ("which moves nothing" if moved == 0 else
                   "which moves the rung and the step at once. The route "
                   "alternates so that exactly one variable moves per entry."))
        if r0 != r1 and d0 != d1:                       # unreachable via `moved`
            raise RuntimeError("a conversion at a changed step")


_check_ceilings()
_check_route()


def stability_ceiling(rung: str) -> float:
    """The coarsest step this rung is MEASURED to start clean at.

    An unprobed rung is an ERROR rather than a fit. `exoplasim/notes/physics-
    filter-stability.md` measures the 1/N rule predicting T127 exactly and
    over-predicting T170 by a full step -- the rung it would have been used to
    plan -- so interpolating here would be reading the one place the rule is
    known to break.
    """
    key = str(rung).upper()
    geometry(key)                                # refuses a non-ladder rung
    if key not in STABILITY_CEILING_MINUTES:
        raise RuntimeError(
            f"{key} is on the ladder and has never been probed, so it has no "
            f"measured stability ceiling. Probed rungs are "
            f"{', '.join(STABILITY_CEILING_MINUTES)}; run "
            "exoplasim/scripts/stability_probe.py before running it.")
    return STABILITY_CEILING_MINUTES[key]


def route_steps(rung: str) -> tuple[float, ...]:
    """Every step the escalation route runs this rung at, in route order.

    More than one is normal and is the point: the route reconverges a rung at
    the NEXT rung's step before converting, so T21 is on the route at 45 and at
    30. Empty for a ladder rung the route does not visit.
    """
    key = str(rung).upper()
    geometry(key)
    return tuple(dt for r, dt in ESCALATION_ROUTE if r == key)


def timestep_problems(rung: str, timestep_minutes: float) -> list[str]:
    """What is wrong with running `rung` at this step. Empty when nothing is.

    TWO DIFFERENT VERDICTS, and they are returned as separate lines because
    they carry different weight. Above the measured ceiling is a step the model
    is measured to refuse or blow up at. Off the route is a step that may be
    perfectly stable and is not the escalation the project declared -- which is
    what a diagnostic arm is, and `filter_timestep_matrix.py` deliberately runs
    T42 at dt 90 to find a trap boundary. A caller that judges the DECLARED
    configuration should treat both as failures; one launching a diagnostic
    should report them and go on.
    """
    key = str(rung).upper()
    problems = []
    try:
        ceiling = stability_ceiling(key)
    except RuntimeError as exc:
        problems.append(str(exc))
    else:
        if float(timestep_minutes) > ceiling:
            problems.append(
                f"dt {timestep_minutes} is coarser than the coarsest step "
                f"{key} is measured to start clean at, {ceiling} "
                f"({STABILITY_GRID}, kappa {STABILITY_GRID_KAPPA}).")
    steps = route_steps(key)
    if not steps:
        problems.append(
            f"the escalation route does not visit {key}. It is "
            + " then ".join(f"{r} at dt {d}" for r, d in ESCALATION_ROUTE)
            + " (docs/src/pipeline/sequencing.md section D).")
    elif not any(abs(float(timestep_minutes) - s) < 1e-9 for s in steps):
        problems.append(
            f"the escalation route runs {key} at "
            + " and ".join(str(s) for s in steps)
            + f", not at {timestep_minutes} "
              "(docs/src/pipeline/sequencing.md section D).")
    evidence = COMMISSIONING_EVIDENCE.get((key, float(timestep_minutes)))
    if (evidence and evidence["verdict"] == "blew_up"
            and evidence.get("binds", True)):
        problems.append(
            f"{key} at dt {timestep_minutes} has been run to commissioning "
            f"length and failed: {evidence['detail']}")
    return problems


def commissioning_caveats(rung: str, timestep_minutes: float) -> list[str]:
    """What is known against this pair that is NOT a reason to refuse it.

    A blow-up whose run no longer exists cannot bar a configuration -- there is
    no artifact to check the claim against -- but it is still the only endurance
    evidence the pair has, and dropping it would leave an attempt at that pair
    looking like an attempt at an untested one. Those are different things and
    the caller is entitled to know which one it is holding.

    Kept apart from `timestep_problems` because the two have different
    consequences: a problem refuses, a caveat is printed.
    """
    key = str(rung).upper()
    evidence = COMMISSIONING_EVIDENCE.get((key, float(timestep_minutes)))
    if not evidence or evidence.get("binds", True):
        return []
    return [f"{key} at dt {timestep_minutes} has a blow-up on record that does "
            f"not bar it. What happened: {evidence['detail']} Why it does not "
            f"bar the pair: {evidence['does_not_bind']}"]


def configured_timestep(config: dict) -> tuple[float, str]:
    """`model.timestep_minutes` and the rung it will be run at, as one fact.

    The value stays in `config/planet.yaml` because it is an OPERATIONAL choice
    that selects a binary and that a diagnostic arm deliberately moves. What
    moved here is the AUTHORITY over it: the ceiling and the route are declared
    once, above, and this is where a configuration meets them.

    Returns `(step, rung)` and raises only through `model_grid`, which refuses
    a resolution paired with another grid's dimensions. The step's own verdict
    is `timestep_problems`, kept separate because whether an off-route step is a
    defect depends on the caller. `scripts/check_consistency.py` and `scripts/smoke_test.py`
    judge the declared configuration and refuse; `run_exoplasim.py` reports.
    """
    rung, _, _ = model_grid(config)
    return float(config["model"]["timestep_minutes"]), rung


def check_stability_ceilings(root) -> list[str]:
    """`STABILITY_CEILING_MINUTES` against the probe grid it was read out of.

    A CHECK WITH A RIGHT ANSWER: the ceiling for a rung is the coarsest dt whose
    cell, at the declared kappa and tau scale, is `no_refusal_in_steps`. Every
    other outcome -- `refused`, and `refused_only_at_length` -- is not a clean
    start, and reading the grid as though it had only those first two is
    precisely how T170 came to carry dt 22.5. WORLD-6QRR.

    Takes the repository root rather than resolving one, so this module keeps
    knowing nothing but the ladder. Returns the empty list when they agree, and
    an empty list when the artifact is absent: the grid is untracked output and
    a worktree without it is not a disagreement.
    """
    import json
    from pathlib import Path

    path = Path(root) / STABILITY_GRID
    if not path.is_file():
        return []
    try:
        probes = json.loads(path.read_text(encoding="utf-8"))["probes"]
    except Exception as exc:                                # noqa: BLE001
        return [f"{STABILITY_GRID} cannot be read as a probe grid: {exc}"]

    clean: dict[str, float] = {}
    seen: set[str] = set()
    for p in probes:
        if p.get("kappa") != STABILITY_GRID_KAPPA:
            continue
        if p.get("tau_scale") != STABILITY_GRID_TAU_SCALE:
            continue
        rung, dt = str(p["rung"]).upper(), float(p["dt_minutes"])
        seen.add(rung)
        if p.get("outcome") == PROBE_CLEAN:
            clean[rung] = max(clean.get(rung, 0.0), dt)

    problems = []
    for rung in sorted(seen | set(STABILITY_CEILING_MINUTES)):
        declared = STABILITY_CEILING_MINUTES.get(rung)
        measured = clean.get(rung)
        if declared is None and measured is not None:
            problems.append(
                f"{STABILITY_GRID} measures {rung} clean at dt {measured} and "
                f"the registry declares no ceiling for it")
        elif measured is None and declared is not None and rung in seen:
            problems.append(
                f"{STABILITY_GRID} has no clean cell for {rung} at kappa "
                f"{STABILITY_GRID_KAPPA} and the registry declares {declared}")
        elif (declared is not None and measured is not None
                and abs(declared - measured) > 1e-9):
            problems.append(
                f"{rung}: the registry declares a ceiling of {declared} and "
                f"the coarsest cell {STABILITY_GRID} marks {PROBE_CLEAN!r} is "
                f"dt {measured}")
    return problems


# The escalation route, restated in the chapter that DECIDES it. The chapter is
# authoritative and this direction of the check is deliberate: the table above
# is what code reads, so it is the copy that can drift silently.
ROUTE_CHAPTER = "docs/src/pipeline/sequencing.md"
ROUTE_IN_CHAPTER = r"\bT([0-9]+),? at dt ([0-9]+(?:\.[0-9]+)?)"


def check_timestep_restatements(root) -> list[str]:
    """The route table against `docs/src/pipeline/sequencing.md` section D.

    Reads the chapter's numbered route in ORDER and requires the same sequence
    of (rung, step) pairs, because the order is the content: the route is a
    sequence in which one variable moves at a time, and a set would lose the
    only thing it asserts.
    """
    import re
    from pathlib import Path

    path = Path(root) / ROUTE_CHAPTER
    if not path.is_file():
        return [f"{ROUTE_CHAPTER} is not there, and it decides the route"]
    text = path.read_text(encoding="utf-8")
    start = text.find("The route escalates resolution and timestep")
    if start < 0:
        return [f"{ROUTE_CHAPTER} no longer carries the route in the shape this "
                "check reads; it decides the route and cannot go unchecked"]
    end = text.find("\n\n**", start)
    got = tuple((f"T{m.group(1)}", float(m.group(2)))
                for m in re.finditer(ROUTE_IN_CHAPTER,
                                     text[start:end if end > 0 else None]))
    if got != ESCALATION_ROUTE:
        return [f"{ROUTE_CHAPTER} section D runs "
                + " then ".join(f"{r} at dt {d}" for r, d in got)
                + " and lib/rungs.py declares "
                + " then ".join(f"{r} at dt {d}" for r, d in ESCALATION_ROUTE)]
    return []
