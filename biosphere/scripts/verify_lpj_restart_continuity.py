#!/usr/bin/env python3
"""Is a continued LPJ-GUESS run the same experiment as the run it continues?

    python biosphere/scripts/verify_lpj_restart_continuity.py --nyear 12 --state-year 8

Worldbuilding frame: this checks the Vesper project's vegetation model against
itself. Nothing here is a claim about the simulated planet.

WHAT IT TESTS. Two runs of the same length over the same forcing: one taken
whole, one stopped at `--state-year` and resumed from the state file it wrote.
Every simulated year from `--state-year` onwards has to agree exactly, table for
table and row for row. That is a check with a WRONG answer rather than a merely
different one, which is what makes it a fixture: state the serializer does not
carry shows up as the resumed run integrating a simulated year the uninterrupted
run never had.

The defect it was built for: `Soil::serialize` carried `wtp` and `awtp` and not
`Wtot`, `wtd`, `stand_water`, `mwtp`, `Frac_ice` or `rootfrac`, so a resumed run
started its first simulated year with an empty acrotelm, a water table at the
surface and no pore ice, whatever the run it continued had ended with.
`biosphere/config/wetlands.yaml` lists the members under
`hydrology.required_restart_members` and `biosphere/scripts/wetland_gate.py`
reads the serializer statically; this is the other half, and it is the half that
can fail on behaviour.

THREE MODES, and the last two localise what the first only detects.

  default       two runs of `--nyear` years, one whole and one split at the
                year boundary before `--state-year`, compared over their OUTPUT
                TABLES from that year on. This is the property every continued
                run depends on and it is what fails first, but a lost quantity
                reaches it as a whole simulated year of divergence.
  --one-day     split at an arbitrary simulated day instead. Both arms write a
                state file for the END of the SAME day, one simulated day after
                the split, and the two files are compared byte for byte. One day
                of arithmetic separates the branch from the comparison, so what
                differs is what that one day did differently, not a year of
                accumulated consequence.
  --round-trip  restart from a state file and write the state again with NO
                simulated day in between. Nothing integrates between the two
                files, so every difference is something the write-and-read of a
                state file does not carry.

WHAT --round-trip CAN AND CANNOT SEE, stated because the answer is not the same
as it is for the climate model. A record ExoPlaSim writes and then overwrites on
the restart path shows up there; `world-8yyh` was found exactly that way. A Soil
member LPJ-GUESS omits from `Soil::serialize` entirely cannot show up here,
because it is absent from both files. So --round-trip is the check on the
serializer's own symmetry, and --one-day is the check that names lost state: one
simulated day is the shortest interval in which a member the restart did not
carry can reach a member it did.

WHY THIS WAS ONCE A YEAR BOUNDARY. LPJ-GUESS serialized exactly once, at the end
of simulation year `state_year - 1`; a restarted run resumed at day 0 of
`state_year` and could never reach that save point again, and a run could not
both restart and save. There was no second state file to compare the first
against. WORLD-FUJ4 gave `framework/framework.cpp` a save point the caller
places on any simulated day (`state_day`), a save point separable from the
restart point (`save_year`, `save_day`), and permission to do both in one run.

WHAT IT NEEDS. A compiled `guess`, a driver, a soil map and a PFT file: the same
inputs `run_lpj_guess.py` needs. None of them exists until a baseline
climatology does, and this script exits naming what is missing rather than
reporting a pass it did not earn: an absent measurement is not a pass.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml  # noqa: E402

from _paths import (CONFIG, GENERATED, GUESS_BINARY, GUESS_SOURCE,  # noqa: E402
                    PROJECT_ROOT, RUNS)

from orbit import model_year_days  # lib/orbit.py, the year length  # noqa: E402

import run_lpj_guess  # noqa: E402
import wetland_gate  # noqa: E402

# The serialization block plib takes as the later declaration. `state_path` is
# absolute for the reason every other path in a generated instruction file is:
# each rank chdirs into its own directory before reading anything.
#
# state_day and save_day name the LAST simulated day the state covers, -1 being
# the year boundary and the default this model had before WORLD-FUJ4. A run
# resumes on the day after state_year/state_day and writes its own state at the
# end of save_year/save_day.
SERIALIZATION = """
! Restart continuity fixture.
state_year {state_year}
state_day {state_day}
save_year {save_year}
save_day {save_day}
save_state {save_state}
restart {restart}
state_path "{state_path}"
save_path "{save_path}"
"""


def rel(path: Path) -> Path:
    try:
        return path.relative_to(PROJECT_ROOT)
    except ValueError:
        return path


def missing_inputs(driver: Path, soilmap: Path, pfts: Path) -> list[str]:
    return [str(p) for p in (driver, soilmap, pfts, GUESS_BINARY)
            if not Path(p).is_file()]


def build_bed(bed: Path, paths: dict, settings: dict, state_dir: Path,
              state_year: int, save_state: int, restart: int,
              state_day: int = -1, save_year: int | None = None,
              save_day: int | None = None, save_dir: Path | None = None) -> Path:
    """One runnable directory, and the instruction file it runs."""
    bed.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(paths["pfts_source"], bed / "vesper_pfts.ins")
    for extra in (GUESS_SOURCE / "data" / "ins").glob("*.ins"):
        if not (bed / extra.name).exists():
            shutil.copyfile(extra, bed / extra.name)
    instruction = bed / "run.ins"
    instruction.write_text(
        run_lpj_guess.build_instruction(paths, settings)
        + SERIALIZATION.format(
            state_year=state_year, state_day=state_day,
            save_year=state_year if save_year is None else save_year,
            save_day=state_day if save_day is None else save_day,
            save_state=save_state, restart=restart,
            state_path=state_dir.resolve(),
            save_path=(state_dir if save_dir is None else save_dir).resolve()))
    return instruction


def state_files(state_dir: Path) -> list[Path]:
    """The per-rank state files, in a stable order, excluding the metadata."""
    return sorted(p for p in state_dir.iterdir()
                  if p.is_file() and p.name != "meta.bin")


def compare_states(left: Path, right: Path, label_left: str,
                   label_right: str) -> list[str]:
    """Every way two state directories fail to hold the same bytes.

    Byte equality, not a tolerance. Both directories were written by the same
    executable from the same forcing at the same simulated instant, so anything
    but identical bytes is a different simulated state.

    The offset of the first difference is reported because it is the only handle
    this format gives on WHICH member differs: the archive is an untagged stream
    of raw object bytes in the order `serialize` writes them, with no names and
    no record boundaries. It locates a difference; it does not name one. Naming
    it is what `biosphere/config/soil_restart_state.yaml` and the daily output
    tables are for.
    """
    failures: list[str] = []
    left_files = state_files(left)
    right_files = state_files(right)
    if not left_files or not right_files:
        failures.append(
            f"no state files: {label_left} wrote {len(left_files)}, "
            f"{label_right} wrote {len(right_files)}")
        return failures
    if [p.name for p in left_files] != [p.name for p in right_files]:
        failures.append(
            f"the two arms wrote different state files: "
            f"{[p.name for p in left_files]} against "
            f"{[p.name for p in right_files]}")
        return failures

    for a, b in zip(left_files, right_files):
        first, second = a.read_bytes(), b.read_bytes()
        if first == second:
            continue
        if len(first) != len(second):
            failures.append(
                f"{a.name}: {len(first)} bytes in {label_left} against "
                f"{len(second)} in {label_right}")
            continue
        offset = next(i for i, (x, y) in enumerate(zip(first, second)) if x != y)
        differing = sum(1 for x, y in zip(first, second) if x != y)
        failures.append(
            f"{a.name}: {differing} of {len(first)} bytes differ, first at "
            f"offset {offset}")
    return failures


def run_bed(instruction: Path, ranks: int, tables: tuple[str, ...]) -> None:
    bed = instruction.parent
    for rank in range(1, ranks + 1):
        (bed / f"run{rank}").mkdir(exist_ok=True)
    command = ["mpirun", "-np", str(ranks), "--bind-to", "core",
               str(GUESS_BINARY), "-parallel", "-input", "vesper",
               str(instruction.resolve())]
    result = subprocess.run(command, cwd=bed, capture_output=True, text=True)
    (bed / "mpirun.log").write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        raise SystemExit(
            f"LPJ-GUESS exited {result.returncode} in {rel(bed)}. See "
            f"{rel(bed / 'mpirun.log')} and {rel(bed)}/run*/guess.log")
    run_lpj_guess.merge_outputs(bed, ranks, tables)


def year_column(header: list[str]) -> int | None:
    for index, name in enumerate(header):
        if name.lower() == "year":
            return index
    return None


def compare(whole: Path, resumed: Path, first_year: int,
            tables: tuple[str, ...]) -> list[str]:
    """Every way the resumed run failed to reproduce the uninterrupted one.

    Exact string equality on the emitted rows, not a tolerance. The two runs
    execute the same arithmetic in the same order over the same forcing, so
    anything but identical output is a different simulated state and not a
    rounding difference. A tolerance here would be a bar chosen to admit
    whatever the runs happened to produce.
    """
    failures: list[str] = []
    for name in tables:
        left, right = whole / name, resumed / name
        if not left.is_file() or not right.is_file():
            failures.append(f"{name}: absent from "
                            f"{'the whole run' if not left.is_file() else ''}"
                            f"{' and ' if not left.is_file() and not right.is_file() else ''}"
                            f"{'the resumed run' if not right.is_file() else ''}")
            continue
        left_lines = left.read_text().splitlines()
        right_lines = right.read_text().splitlines()
        if not left_lines or not right_lines:
            failures.append(f"{name}: empty")
            continue
        header = left_lines[0].split()
        if header != right_lines[0].split():
            failures.append(f"{name}: the two runs wrote different columns")
            continue
        column = year_column(header)
        if column is None:
            failures.append(f"{name}: no Year column, so the rows the resumed "
                            "run is responsible for cannot be selected")
            continue

        def rows(lines):
            kept = {}
            for line in lines[1:]:
                parts = line.split()
                if len(parts) <= column:
                    continue
                try:
                    year = int(float(parts[column]))
                except ValueError:
                    continue
                if year >= first_year:
                    kept[(parts[0], parts[1], year)] = line.split()
            return kept

        left_rows, right_rows = rows(left_lines), rows(right_lines)
        if set(left_rows) != set(right_rows):
            failures.append(
                f"{name}: the resumed run covers {len(right_rows)} rows from "
                f"year {first_year} against {len(left_rows)} in the "
                "uninterrupted run")
            continue
        differing = [key for key, row in left_rows.items()
                     if row != right_rows[key]]
        if differing:
            lon, lat, year = differing[0]
            failures.append(
                f"{name}: {len(differing)} of {len(left_rows)} rows differ, "
                f"first at lon {lon} lat {lat} year {year}")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nyear", type=int, default=12,
                        help="simulated years in each run")
    parser.add_argument("--state-year", type=int, default=8,
                        help="the year the split run resumes at; the first "
                             "year the two runs are compared over")
    parser.add_argument("--state-day", type=int, default=None,
                        help="the simulated day within --state-year the split "
                             "is taken at, for --one-day and --round-trip. "
                             "Default: the middle of the year, which is the one "
                             "day of it that is neither the reset nor the flush")
    parser.add_argument("--one-day", action="store_true",
                        help="instead: split at --state-day and compare the two "
                             "arms' STATE FILES one simulated day later")
    parser.add_argument("--round-trip", action="store_true",
                        help="instead: restart at --state-day and write the "
                             "state again with NO simulated day in between")
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--npatch", type=int, default=5)
    parser.add_argument("--bed", type=Path, default=None,
                        help="where to build the two runs (default: a "
                             "restart_continuity directory under runs/)")
    parser.add_argument("--driver", type=Path,
                        default=GENERATED / "vesper_driver.bin")
    parser.add_argument("--soilmap", type=Path, default=None)
    parser.add_argument("--pfts", type=Path,
                        default=GENERATED / "vesper_pfts.ins")
    args = parser.parse_args()

    if not 0 < args.state_year < args.nyear:
        raise SystemExit(
            f"--state-year {args.state_year} has to fall inside the run: a "
            "restart at year 0 continues nothing and a restart at the last "
            "year compares nothing.")

    if args.one_day and args.round_trip:
        raise SystemExit(
            "--one-day and --round-trip are two different questions: one asks "
            "what a simulated day did differently either side of a restart, the "
            "other asks what the write and read of a state file does not carry. "
            "Ask one.")

    if args.soilmap is None:
        import builds
        args.soilmap = builds.soilmap()

    absent = missing_inputs(args.driver, args.soilmap, args.pfts)
    if absent:
        raise SystemExit(
            "the restart continuity fixture has NOT been run. It needs a "
            "compiled model and a built forcing, and these are missing:\n  "
            + "\n  ".join(absent) +
            "\n\nLPJ-GUESS does not build on this tree until a baseline "
            "climatology exists. This is an absent measurement, not a pass.")

    # The retained tables. The wetland set is included only when the gate has
    # granted activation, on exactly the terms run_lpj_guess.py includes it:
    # a table the run does not write cannot be compared, and asking for one
    # would make this fixture fail for a reason that is not a restart defect.
    declaration = wetland_gate.read_declaration()
    active = wetland_gate.granted(declaration)[0]
    tables = run_lpj_guess.OUTPUTS + (
        tuple(declaration["acceptance"]["retained_outputs"]) if active else ())

    bed_root = args.bed or (RUNS / "restart_continuity")
    whole, resumed = bed_root / "whole", bed_root / "resumed"
    state_dir = bed_root / "state"
    whole_state = bed_root / "state_whole"
    for path in (whole, resumed, bed_root / "continued"):
        if path.exists():
            raise SystemExit(f"{path} exists; remove it before re-running")
    state_dir.mkdir(parents=True, exist_ok=True)

    def paths_for(bed: Path) -> dict:
        return {"driver": Path(args.driver).resolve(),
                "soilmap": Path(args.soilmap).resolve(),
                "pfts": (bed / "vesper_pfts.ins").resolve(),
                "pfts_source": Path(args.pfts).resolve()}

    settings = {"nyear": args.nyear, "npatch": args.npatch,
                "root_seed": 20260828,
                "nfix_a": 0.234, "nfix_b": -0.172, "ifbvoc": 0,
                "outputs": tables, **wetland_gate.switches(active)}

    # The default split day. The middle of the simulation year is chosen before
    # any result is seen and for a stated reason: day 0 is where every annual
    # accumulator resets and the last day is where they flush, so both are days
    # on which a member that is lost the rest of the year looks carried.
    split_day = (args.state_day if args.state_day is not None
                 else model_year_days(yaml.safe_load(CONFIG.read_text())) // 2)
    mode = ("round-trip" if args.round_trip
            else "one-day" if args.one_day else "annual")

    if mode == "annual":
        run_bed(build_bed(whole, paths_for(whole),
                          {**settings, "title": "restart_continuity_whole"},
                          state_dir, args.state_year, 1, 0), args.ranks, tables)
        run_bed(build_bed(resumed, paths_for(resumed),
                          {**settings, "title": "restart_continuity_resumed"},
                          state_dir, args.state_year, 0, 1), args.ranks, tables)
        failures = compare(whole, resumed, args.state_year, tables)
        headline = (f"{args.nyear} simulated years whole against a resume at "
                    f"year {args.state_year}; {len(tables)} tables compared")
        verdict = ("The resumed run is not the run it continues. Every "
                   "difference above is simulated state the serializer does "
                   "not carry.")

    elif mode == "round-trip":
        # One arm writes the state at the end of `split_day`. The other reads
        # that state and writes it straight back out, having simulated nothing.
        run_bed(build_bed(whole, paths_for(whole),
                          {**settings, "title": "restart_round_trip_written"},
                          state_dir, args.state_year, 1, 0,
                          state_day=split_day), args.ranks, tables)
        whole_state.mkdir(parents=True, exist_ok=True)
        run_bed(build_bed(resumed, paths_for(resumed),
                          {**settings, "title": "restart_round_trip_rewritten"},
                          state_dir, args.state_year, 1, 1,
                          state_day=split_day,
                          save_dir=whole_state), args.ranks, tables)
        failures = compare_states(state_dir, whole_state,
                                  "written", "rewritten")
        headline = (f"restart round trip at year {args.state_year} day "
                    f"{split_day}: written against rewritten, no simulated day "
                    "in between")
        verdict = ("A difference here is state the write and read of a state "
                   "file does not carry. It cannot be a member Soil::serialize "
                   "omits entirely -- such a member is absent from both files "
                   "-- so it is one the deserialize path reads and then "
                   "overwrites.")

    else:
        # Both arms write the state for the END of the same simulated day. The
        # whole arm reaches it uninterrupted; the split arm reaches it one day
        # after resuming from the state written the day before.
        whole_state.mkdir(parents=True, exist_ok=True)
        run_bed(build_bed(whole, paths_for(whole),
                          {**settings, "title": "restart_one_day_whole"},
                          whole_state, args.state_year, 1, 0,
                          state_day=split_day + 1), args.ranks, tables)
        run_bed(build_bed(resumed, paths_for(resumed),
                          {**settings, "title": "restart_one_day_split"},
                          state_dir, args.state_year, 1, 0,
                          state_day=split_day), args.ranks, tables)
        continued = bed_root / "continued"
        continued_state = bed_root / "state_continued"
        continued_state.mkdir(parents=True, exist_ok=True)
        run_bed(build_bed(continued, paths_for(continued),
                          {**settings, "title": "restart_one_day_continued"},
                          state_dir, args.state_year, 1, 1,
                          state_day=split_day,
                          save_year=args.state_year,
                          save_day=split_day + 1,
                          save_dir=continued_state), args.ranks, tables)
        failures = compare_states(whole_state, continued_state,
                                  "whole", "split")
        headline = (f"one simulated day either side of a restart at year "
                    f"{args.state_year} day {split_day}: state at the end of "
                    f"day {split_day + 1}, whole against split")
        verdict = ("The two arms are in different simulated states one day "
                   "after the split. One day of arithmetic separates the branch "
                   "from this comparison, so what differs is what that day did "
                   "differently -- state the serializer did not carry reaching "
                   "state it did.")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "biosphere/scripts/verify_lpj_restart_continuity.py",
        "mode": mode,
        "nyear": args.nyear, "state_year": args.state_year,
        "state_day": split_day if mode != "annual" else None,
        "ranks": args.ranks, "npatch": args.npatch,
        "wetlands_active": active,
        "tables_compared": list(tables) if mode == "annual" else [],
        "continuous": not failures,
        "failures": failures,
    }
    out = GENERATED / "lpj_restart_continuity.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")

    print(headline)
    for failure in failures:
        print(f"  {failure}")
    print(f"\nwrote {rel(out)}")
    if failures:
        print(f"\n{verdict}")
        sys.exit(1)


if __name__ == "__main__":
    main()
