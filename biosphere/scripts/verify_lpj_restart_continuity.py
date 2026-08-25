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

WHY IT IS A YEAR BOUNDARY AND NOT MID-YEAR. LPJ-GUESS serializes exactly once,
at the end of simulation year `state_year - 1`, and a restarted run sets
`date.year = state_year` and begins at day 0. `framework/framework.cpp` has no
other save point, and a restarted run can never reach `state_year - 1` again, so
neither a mid-year stop nor ExoPlaSim's zero-step round trip is expressible
here: there is no second state file to compare the first against. The year
boundary is what the model has, and it is enough, because `state_year` is the
first simulated year every unserialized quantity can differ in. Giving this
model a mid-year save point, and with it the zero-step round trip that names a
lost quantity directly instead of reporting a year of divergence, is
WORLD-FUJ4.

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

from _paths import (GENERATED, GUESS_BINARY, GUESS_SOURCE,  # noqa: E402
                    PROJECT_ROOT, RUNS)

import run_lpj_guess  # noqa: E402
import wetland_gate  # noqa: E402

# The serialization block plib takes as the later declaration. `state_path` is
# absolute for the reason every other path in a generated instruction file is:
# each rank chdirs into its own directory before reading anything.
SERIALIZATION = """
! Restart continuity fixture. save_state writes at the END of year
! state_year - 1; restart jumps date.year to state_year and starts at day 0.
state_year {state_year}
save_state {save_state}
restart {restart}
state_path "{state_path}"
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
              state_year: int, save_state: int, restart: int) -> Path:
    """One runnable directory, and the instruction file it runs."""
    bed.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(paths["pfts_source"], bed / "vesper_pfts.ins")
    for extra in (GUESS_SOURCE / "data" / "ins").glob("*.ins"):
        if not (bed / extra.name).exists():
            shutil.copyfile(extra, bed / extra.name)
    instruction = bed / "run.ins"
    instruction.write_text(
        run_lpj_guess.build_instruction(paths, settings)
        + SERIALIZATION.format(state_year=state_year, save_state=save_state,
                               restart=restart,
                               state_path=state_dir.resolve()))
    return instruction


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
    for path in (whole, resumed):
        if path.exists():
            raise SystemExit(f"{path} exists; remove it before re-running")
    state_dir.mkdir(parents=True, exist_ok=True)

    def paths_for(bed: Path) -> dict:
        return {"driver": Path(args.driver).resolve(),
                "soilmap": Path(args.soilmap).resolve(),
                "pfts": (bed / "vesper_pfts.ins").resolve(),
                "pfts_source": Path(args.pfts).resolve()}

    settings = {"nyear": args.nyear, "npatch": args.npatch,
                "nfix_a": 0.234, "nfix_b": -0.172, "ifbvoc": 0,
                "outputs": tables, **wetland_gate.switches(active)}

    run_bed(build_bed(whole, paths_for(whole),
                      {**settings, "title": "restart_continuity_whole"},
                      state_dir, args.state_year, 1, 0), args.ranks, tables)
    run_bed(build_bed(resumed, paths_for(resumed),
                      {**settings, "title": "restart_continuity_resumed"},
                      state_dir, args.state_year, 0, 1), args.ranks, tables)

    failures = compare(whole, resumed, args.state_year, tables)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "biosphere/scripts/verify_lpj_restart_continuity.py",
        "nyear": args.nyear, "state_year": args.state_year,
        "ranks": args.ranks, "npatch": args.npatch,
        "wetlands_active": active,
        "tables_compared": list(tables),
        "continuous": not failures,
        "failures": failures,
    }
    out = GENERATED / "lpj_restart_continuity.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")

    print(f"{args.nyear} simulated years whole against a resume at year "
          f"{args.state_year}; {len(tables)} tables compared")
    for failure in failures:
        print(f"  {failure}")
    print(f"\nwrote {rel(out)}")
    if failures:
        print("\nThe resumed run is not the run it continues. Every difference "
              "above is simulated state the serializer does not carry.")
        sys.exit(1)


if __name__ == "__main__":
    main()
