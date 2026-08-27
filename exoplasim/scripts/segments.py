#!/usr/bin/env python3
"""What a run's segments were FOR, and which orbits that makes usable.

A run is made of segments: contiguous blocks of orbits added by
`continue_exoplasim.py`, each recorded in `run_manifest.json`. A segment carries
its I/O regime and its PURPOSE, and both change what an orbit may be used for.

**The purpose is declared by the caller, never inferred from the flags.** It was
inferred once -- seasonal output plus a run past its equilibrium cutoff was
labelled a climatology segment -- and that labelled a three-orbit verification
run, made with `--low-io` on a different binary, as climatology input. The flags
say what the model was asked to write; only the person running it knows what the
orbits are for, and a rule over the flags will keep getting that wrong in a new
way each time a flag is added.

This module owns the purpose vocabulary; route new questions about what an
orbit was for through it rather than reading `segments` raw. A third ad-hoc
reader was once nearly written, which is how this project came to have four
copies of a path resolver.
"""

from __future__ import annotations

import json
from pathlib import Path

# Every value `purpose` may take. The first two are the run's own trajectory:
# orbits that carry the planet toward equilibrium, and orbits at equilibrium
# meant to be read as its climate. Both belong in a convergence window.
#
# `diagnostic` is a segment run to measure the MODEL rather than to advance the
# planet: a low-I/O verification, a high-cadence wind sample for DUST-5, a short
# block on a differently patched binary. Those orbits are real integrations and
# stay in the run, but they are not evidence about where the run is settling,
# and a tail of them must not enter a convergence window or a climatology.
SEGMENT_PURPOSES = ("spinup", "post_equilibrium_climatology", "diagnostic")
PRODUCTION_PURPOSES = frozenset({"spinup", "post_equilibrium_climatology"})


def segment_records(run_dir: Path) -> list[dict]:
    """The manifest's segment list, or empty when there is no manifest."""
    manifest = Path(run_dir) / "run_manifest.json"
    if not manifest.is_file():
        return []
    return json.loads(manifest.read_text(encoding="utf-8")).get("segments", [])


def orbit_purposes(run_dir: Path, orbits) -> dict[int, str | None]:
    """Purpose per orbit index, None where no segment covers it.

    Orbit 0 is normally uncovered: `run_exoplasim.py` writes it and only
    continuations append segments. Runs made before segments were recorded are
    uncovered throughout.
    """
    segments = segment_records(run_dir)
    out: dict[int, str | None] = {}
    for orbit in orbits:
        out[orbit] = None
        for seg in segments:
            if seg["start_year_index"] <= orbit <= seg["end_year_index"]:
                out[orbit] = seg.get("purpose")
                break
    return out


def non_production_orbits(run_dir: Path, orbits) -> list[int]:
    """Orbits in `orbits` that a segment declares are not the run's trajectory.

    **An orbit no segment covers counts as production**, which is the OPPOSITE
    default from `low_io_orbits` below, deliberately. There the unlabelled case
    is a silent corruption of a field, so the honest default is the unsafe one.
    Here the unlabelled case is every orbit of every run made before purposes
    were declared, plus orbit zero of every run ever: defaulting those out would
    not be conservative, it would refuse to assess anything. The protection here
    comes from the declaration instead -- an orbit is excluded only when its
    segment SAYS what it was for and that was not the run's trajectory.
    """
    purposes = orbit_purposes(run_dir, orbits)
    return [orbit for orbit, purpose in purposes.items()
            if purpose is not None and purpose not in PRODUCTION_PURPOSES]


def production_window(run_dir: Path, n_orbits: int, window: int) -> tuple[int, int]:
    """The last `window` production orbits, as inclusive year indices.

    Trailing non-production orbits are dropped, which is the case this exists
    for: a diagnostic tail on a different binary, silently averaged into a
    convergence window because the window was "the last ten orbits".

    A non-production orbit INSIDE the resulting window raises instead. A window
    with a hole in it is not a trend, and no rule for filling it is better than
    the caller saying which orbits they meant: shorten `--window`, or assess the
    block before the interruption.

    A window that spans a change of I/O REGIME raises on the same terms, for
    the reason `refuse_a_window_spanning_an_io_regime_change` gives: the two
    regimes are two instruments and the join between them is a step. That rule
    lives here rather than in each caller because a window is where a verdict
    is taken, and there is no verdict the step is admissible in.
    """
    refuse_orbits_no_segment_covers(run_dir, n_orbits)
    excluded = set(non_production_orbits(run_dir, range(n_orbits)))
    end = n_orbits - 1
    while end >= 0 and end in excluded:
        end -= 1
    if end < 0:
        raise RuntimeError(
            f"{Path(run_dir).name} has no production orbits: every segment "
            f"declares a purpose outside {sorted(PRODUCTION_PURPOSES)}")
    start = end - window + 1
    if start < 0:
        raise RuntimeError(
            f"{Path(run_dir).name} has {end + 1} orbits up to the last "
            f"production one (index {end}), fewer than the {window}-orbit "
            f"window requested")
    interior = sorted(o for o in excluded if start <= o <= end)
    if interior:
        raise RuntimeError(
            f"orbits {interior} are inside the {window}-orbit window ending at "
            f"{end} but are not production orbits. A convergence window with a "
            f"hole in it is not a trend; shorten --window, or assess the block "
            f"before them.")
    refuse_a_window_spanning_an_io_regime_change(run_dir, start, end, window)
    return start, end


def low_io_orbits(run_dir: Path, orbits) -> list[int]:
    """Orbits in `orbits` that were run with PlaSim's low-I/O accumulation on.

    Those hold interval accumulations rather than instantaneous samples, which
    cannot be undone; runs written before the first-record patch additionally
    carry a corrupt first record per orbit (wind 7.5x, humidity 27% low).
    Averaging them into a climatology puts that into anything downstream that
    reads a wind or a humidity -- the Penman evaporation the carve criterion
    turns on, and the gust distribution the dust emission turns on.

    **A segment with no `low_io` key is treated as low-I/O**, because every run
    made before 2026-08-17 was, and the honest default for an unlabelled orbit is
    the unsafe one.
    """
    segments = segment_records(run_dir)
    if not segments and not (Path(run_dir) / "run_manifest.json").is_file():
        return list(orbits)
    tainted = []
    for orbit in orbits:
        for seg in segments:
            if seg["start_year_index"] <= orbit <= seg["end_year_index"]:
                if seg.get("low_io", True):
                    tainted.append(orbit)
                break
        else:
            tainted.append(orbit)
    return tainted


def refuse_orbits_no_segment_covers(run_dir: Path, n_orbits: int) -> None:
    """Refuse when the disk holds orbits past the last segment the manifest has.

    A KILLED CONTINUATION LEAVES EXACTLY THIS. The model writes its orbits and
    the segment record is appended when the block finishes, so a run stopped
    partway has the orbits and no line saying what they were. `low_io_orbits`
    then reads them as low-I/O -- its documented default, and the safe one when
    "unlabelled" means a run made before segments existed -- and here
    unlabelled means UNRECORDED, where the same default is simply wrong: the
    orbits are clean, the join is invisible, and
    `refuse_a_window_spanning_an_io_regime_change` sees nothing to refuse.

    That is the guard being walked around rather than failing, which is worse
    than either. So the manifest's coverage is checked against the orbit count
    the caller is asking about, and a trailing gap refuses.

    A LEADING gap does not, and the difference is not cosmetic. Orbit 0 is
    uncovered on runs made before the runner recorded its own first block, and
    every orbit of every run made before segments existed is uncovered; those
    are legitimate and `non_production_orbits` documents why they count as
    production. What cannot be legitimate is an orbit AFTER a segment that
    claims to end earlier: something wrote it and nothing said what it was.
    """
    segments = segment_records(run_dir)
    if not segments:
        return                       # no manifest, or a run that predates them
    covered = max(int(s["end_year_index"]) for s in segments)
    if n_orbits - 1 <= covered:
        return
    raise RuntimeError(
        f"{Path(run_dir).name} has {n_orbits} orbits on disk and its manifest's "
        f"last segment ends at {covered}, so orbits {covered + 1} to "
        f"{n_orbits - 1} are recorded nowhere. Nothing says what I/O regime or "
        f"what purpose they were run at, and the unlabelled default reads them "
        f"as low-I/O production -- which is how a killed continuation's CLEAN "
        f"orbits become invisible to the regime guard.\n"
        f"  A continuation that was interrupted is the usual cause. Re-run it "
        f"so the segment is recorded, or delete the uncovered orbits; do not "
        f"hand-write the segment, because then the manifest asserts a regime "
        f"nobody observed."
    )


def io_regime_changes(run_dir: Path, orbits) -> list[int]:
    """The first orbit of each block whose I/O regime differs from the one before.

    `orbits` is taken in the order given and read as one series, so the answer
    is about that series and not about the run's segment list: a caller that
    hands over a subrange gets the changes inside its own subrange.
    """
    tainted = set(low_io_orbits(run_dir, orbits))
    changes, previous = [], None
    for orbit in orbits:
        regime = orbit in tainted
        if previous is not None and regime != previous:
            changes.append(orbit)
        previous = regime
    return changes


def refuse_a_window_spanning_an_io_regime_change(
        run_dir: Path, start: int, end: int, window: int) -> None:
    """A verdict window has to be one I/O regime, and this is why.

    PlaSim's low-I/O accumulation writes interval accumulations where the clean
    regime writes instantaneous samples, and runs written before the
    first-record patch additionally carry a corrupt first record per orbit. The
    two regimes therefore report the same planet on two instruments, and the
    join between them is a STEP in every series taken across it.

    That step is not a small bias to be tolerated. Fitted on this project's
    82-orbit bootstrap, the integrated autocorrelation time of the per-orbit
    mean surface temperature came back at 9.08 orbits; fitted on each side of
    the join separately it is 2.00 on the low-I/O block and 1.00 and stationary
    on the clean one. The whole of the difference was a +0.1616 K step at the
    boundary, read by the estimator as memory the planet does not have. A
    memory time is what sizes the window and prices the production span, so a
    tau inflated ninefold buys orbits nobody needs -- and the same step tilts
    the trend the convergence verdict is taken on.

    A REFUSAL AND NOT A CLIP, on `production_window`'s own terms: no rule for
    which side to keep is better than the caller saying which orbits they
    meant. The message names the join and the orbits available on the near side
    of it, which is what a caller needs to shorten `--window` or to decide the
    run has not yet bought enough orbits in one regime to be judged.
    """
    inside = io_regime_changes(run_dir, range(start, end + 1))
    if not inside:
        return
    join = inside[0]
    raise RuntimeError(
        f"the {window}-orbit window {start}-{end} spans a change of I/O "
        f"regime at orbit {join}: PlaSim's low-I/O accumulation and the clean "
        f"stream are two instruments, and every series across the join carries "
        f"a step that an autocorrelation reads as memory. {end - join + 1} "
        f"orbits are available on the near side; shorten --window to that or "
        f"below, or extend the run in one regime.")
