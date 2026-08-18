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

This module is the one reader of the `segments` list. There were two, and a
third was about to be written for the convergence window, which is how this
project came to have four copies of a path resolver.
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
    """
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
    return start, end


def low_io_orbits(run_dir: Path, orbits) -> list[int]:
    """Orbits in `orbits` that were run with PlaSim's low-I/O accumulation on.

    Those carry a corrupt first output record per orbit: bottom-level wind reads
    about 7.5x the other bins and humidity 27% low, while every scalar is within
    2%. Averaging them into a climatology puts that into anything downstream that
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
