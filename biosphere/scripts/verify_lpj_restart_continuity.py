#!/usr/bin/env python3
"""Is a continued LPJ-GUESS run the same experiment as the run it continues?

    python biosphere/scripts/verify_lpj_restart_continuity.py --nyear 12 --state-year 8
    python biosphere/scripts/verify_lpj_restart_continuity.py --self-test

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

`--self-test` is a fourth thing and needs neither model nor forcing: it checks
the two integers `run_lpj_guess.py --save-state` and `--continue-from` decide
against the arithmetic `framework/framework.cpp` performs on them, so a
continuation that would silently skip or repeat a simulated year is refused
before anything is bought. The three modes above test the MODEL's serializer;
that tests the RUNNER's wiring, and neither covers the other.

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
import re
import shutil
import struct
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


# The VESPDRV8 driver layout, as `build_lpj_driver.py` writes it and
# `vendor/lpj-guess/modules/vesperinput.cpp` reads it. Restated here for one
# purpose only -- taking a subset of the cells out of a driver file -- and the
# two consumers of the layout are checked against each other by
# `slice_driver`, which refuses a file whose length the layout does not
# explain.
DRIVER_MAGIC = b"VESPDRV8"
DRIVER_HEADER = "<6i2d"          # cells, intervals, year length, years, subdaily, pad, CO2, Ndep
DRIVER_PROVENANCE_BYTES = 64
DRIVER_INTERVAL_BYTES = 4 * 8    # start, end, duration, and the fourth double


def slice_driver(driver: Path, wanted: list[tuple[float, float]],
                 out: Path) -> list[tuple[float, float]]:
    """Write a driver file holding only the named cells, and say which it took.

    WHY A SUBSET IS A LEGITIMATE BED HERE, and it is worth stating because a
    subset of a coupled model would not be. LPJ-GUESS simulates each gridcell
    independently: no cell reads another cell's state, and this project's
    stochastic streams are keyed by COORDINATE rather than by traversal order
    or rank (`derive_stochastic_seed`, DEMO-5). So a cell in a two-cell driver
    integrates the same trajectory it integrates in a sixteen-hundred-cell one,
    and a restart defect that reaches that cell reaches it either way.

    WHAT IT BUYS. The whole grid at the spin-up floor this fixture needs costs
    about thirteen minutes per arm on the shared host; four cells cost seconds.
    That is the difference between a measurement taken once and a member-by-
    member search, which is what naming lost state actually takes.

    WHAT IT COSTS, and why the verdict is written elsewhere. A subset is a
    weaker statement than the grid: a defect confined to cells it does not hold
    is invisible to it. `main` therefore writes a subsetted run's report to a
    different path, and the production gate in `run_lpj_guess.py` keeps reading
    only the whole-grid one.
    """
    raw = driver.read_bytes()
    if raw[:8] != DRIVER_MAGIC:
        raise SystemExit(
            f"{rel(driver)} does not carry the {DRIVER_MAGIC.decode()} magic, "
            "so this is not a driver file this fixture can subset.")
    head = struct.calcsize(DRIVER_HEADER)
    (ncells, nintervals, year_length, nyears, subdaily,
     _pad, co2, ndep) = struct.unpack(DRIVER_HEADER, raw[8:8 + head])
    offset = 8 + head + DRIVER_PROVENANCE_BYTES
    provenance = raw[8 + head:offset]
    offset += nyears * nintervals * DRIVER_INTERVAL_BYTES
    span = nintervals * nyears
    body = len(raw) - offset
    if ncells < 1 or body % ncells:
        raise SystemExit(
            f"{rel(driver)} has {body} bytes of cell records for {ncells} "
            "cells, which is not a whole number of records. The layout this "
            "script restates does not describe this file; re-read "
            "build_lpj_driver.py.")
    cell_bytes = body // ncells
    # lon, lat, soilcode+pad, four retention numbers, one usable share per soil
    # layer, then the four forcing arrays.
    layers = (cell_bytes - 2 * 8 - 2 * 4 - 4 * 8 - 4 * span * 8) // 8
    if layers < 1:
        raise SystemExit(
            f"{rel(driver)} has {cell_bytes}-byte cell records, too short for "
            f"{span} forcing samples. The layout does not describe this file.")

    coords = []
    for index in range(ncells):
        at = offset + index * cell_bytes
        lon, lat = struct.unpack("<dd", raw[at:at + 16])
        coords.append((lon, lat))

    chosen: list[int] = []
    taken: list[tuple[float, float]] = []
    for lon, lat in wanted:
        best = min(range(ncells),
                   key=lambda i: (coords[i][0] - lon) ** 2 + (coords[i][1] - lat) ** 2)
        if (abs(coords[best][0] - lon) > 0.5 or abs(coords[best][1] - lat) > 0.5):
            raise SystemExit(
                f"no cell in {rel(driver)} is within half a degree of "
                f"{lon},{lat}; the nearest is {coords[best][0]},{coords[best][1]}")
        if best in chosen:
            raise SystemExit(
                f"{lon},{lat} names the same cell as one already asked for, "
                f"{coords[best][0]},{coords[best][1]}")
        chosen.append(best)
        taken.append(coords[best])

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as handle:
        handle.write(DRIVER_MAGIC)
        handle.write(struct.pack(DRIVER_HEADER, len(chosen), nintervals,
                                 year_length, nyears, subdaily, 0, co2, ndep))
        handle.write(provenance)
        handle.write(raw[8 + head + DRIVER_PROVENANCE_BYTES:offset])
        for index in chosen:
            at = offset + index * cell_bytes
            handle.write(raw[at:at + cell_bytes])
    # THE BUILD IDENTITY TRAVELS WITH THE SLICE. A subset of a build's driver is
    # that build's forcing over fewer cells, and `provenance.require_build`
    # refuses an unstamped driver outright where the pairing matters. Copying
    # the sidecar is what makes a sliced bed usable by run_lpj_guess.py at all.
    sidecar = driver.with_name(driver.stem + "_provenance.json")
    if sidecar.is_file():
        shutil.copyfile(sidecar, out.with_name(out.stem + "_provenance.json"))
    return taken


def parse_cells(text: str) -> list[tuple[float, float]]:
    """`lon,lat` pairs separated by semicolons, as the output tables print them."""
    cells = []
    for piece in text.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        parts = piece.split(",")
        if len(parts) != 2:
            raise SystemExit(
                f"--cells takes lon,lat pairs separated by ';', not {piece!r}")
        try:
            cells.append((float(parts[0]), float(parts[1])))
        except ValueError:
            raise SystemExit(f"--cells: {piece!r} is not a lon,lat pair")
    if not cells:
        raise SystemExit("--cells was given no cells")
    return cells


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
    state = {
        "restart": bool(restart), "save_state": bool(save_state),
        "state_year": state_year, "state_day": state_day,
        "save_year": state_year if save_year is None else save_year,
        "save_day": state_day if save_day is None else save_day,
        "state_path": str(state_dir.resolve()),
        "save_path": str((state_dir if save_dir is None
                          else save_dir).resolve()),
    }
    # The serialization block is `run_lpj_guess.serialization_block`, reached
    # through `build_instruction`, so this fixture and the production runner
    # cannot drift on what `state_year` means while both keep parsing.
    instruction.write_text(run_lpj_guess.build_instruction(
        paths, {**settings, "state": state}))
    # What this arm ASKED FOR, beside the instruction that asks it. check_instants
    # holds the model's own log against this rather than against a second copy
    # of the arithmetic.
    (bed / "asked_state.json").write_text(json.dumps(state, indent=2) + "\n")
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


def parsed_instants(bed: Path, ranks: int) -> dict:
    """What the model says it parsed, read back off its own log.

    THE DEFECT THIS EXISTS FOR. `state_day -1` and `save_day -1` are the year
    boundary, and they are the only negative integers any instruction file in
    this model declares. `libraries/plib` stored an integer as
    `(int)(num + 0.5)`, which truncates toward zero, so both arrived as 0: every
    restart taken here saved at the end of day 0 of `state_year` and resumed at
    day 1 of it, an arbitrary-day restart, while the instruction file, this
    fixture and the runner all said year boundary. Nothing refused it, because
    0 is in range and every number downstream was computed from 0 consistently.

    A fixture that only restates the rule cannot catch that: it would compute
    -1 from its own arithmetic and agree with itself. So `framework.cpp` prints
    what it PARSED and what it DERIVED, and this reads that back and compares it
    with what the instruction file asked for. That is a check with a wrong
    answer rather than one that can only differ.
    """
    instants: dict = {}
    for rank in range(1, ranks + 1):
        log = bed / f"run{rank}" / "guess.log"
        if not log.is_file():
            continue
        text = log.read_text(errors="replace")
        found: dict = {}
        for line in text.splitlines():
            head = line.strip()
            if head.startswith("Restart instants:"):
                parts = head.replace("Restart instants:", "").split()
                found["state_day"] = int(parts[1])
                found["save_day"] = int(parts[3])
            elif head.startswith("resume: state covers year"):
                parts = head.split()
                found["resume"] = (int(parts[4]), int(parts[6]))
            elif head.startswith("save: state written covers year"):
                parts = head.split()
                found["save"] = (int(parts[5]), int(parts[7]))
        if found:
            return found
    return instants


def check_instants(bed: Path, ranks: int, label: str,
                   year_length: int) -> list[str]:
    """Every way the model's restart instants differ from the ones asked for."""
    asked = json.loads((bed / "asked_state.json").read_text())
    if not (asked["restart"] or asked["save_state"]):
        return []
    seen = parsed_instants(bed, ranks)
    if not seen:
        return [f"{label}: the model logged no restart instants, so what it "
                "parsed from the instruction file is not known. Rebuild: "
                "framework.cpp prints them whenever a run restarts or saves."]
    failures = []
    for name in ("state_day", "save_day"):
        if seen.get(name) != asked[name]:
            failures.append(
                f"{label}: the instruction file declares {name} "
                f"{asked[name]} and the model parsed {seen.get(name)}")
    wanted = resume_and_save_instants(asked, year_length)
    if "resume" in seen and seen["resume"] != wanted["state_covers"]:
        failures.append(
            f"{label}: the model resumes from a state covering {seen['resume']} "
            f"and the block asks for {wanted['state_covers']}")
    if "save" in seen and seen["save"] != wanted["state_written_covers"]:
        failures.append(
            f"{label}: the model writes a state covering {seen['save']} and the "
            f"block asks for {wanted['state_written_covers']}")
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


RUN_ID = re.compile(r"^(lpj_[0-9a-f]{32})$")


def reduced_pfts(source: Path, target: Path, spinup: int) -> None:
    """A copy of the PFT file with the derived spin-up replaced by a short one.

    The runner reads `nyear_spinup` out of the PFT file it imports, and the
    generated one declares the DERIVED floor: thousands of simulated years, and
    the whole grid behind them. That is the right number for a run that buys a
    record and the wrong one for a check on the runner's PLUMBING, which is what
    this mode tests -- state directory creation, the continuation block in the
    manifest, the refusals, and whether a continuation's rows are the rows the
    uninterrupted run has. None of that is a function of the spin-up length.

    The floor stays where it is; this writes a copy into the bed. The copy is a
    bed, not a run: nothing it produces is a result about this world.
    """
    text = source.read_text(encoding="utf-8")
    replaced, count = re.subn(r"(?m)^nyear_spinup\s+\d+",
                             f"nyear_spinup {spinup}", text)
    if count != 1:
        raise SystemExit(
            f"{rel(source)} declares nyear_spinup {count} times and this needs "
            "exactly one to replace")
    target.write_text(replaced, encoding="utf-8")


def invoke_runner(bed: Path, driver: Path, pfts: Path, ranks: int,
                  npatch: int, nyear: int, extra: list[str],
                  expect_failure: bool = False) -> tuple[str | None, str]:
    """One `run_lpj_guess.py` invocation, and the run id it made."""
    command = [sys.executable,
               str(PROJECT_ROOT / "biosphere" / "scripts" / "run_lpj_guess.py"),
               "--nyear", str(nyear), "--npatch", str(npatch),
               "--ranks", str(ranks), "--driver", str(driver),
               "--pfts", str(pfts)] + extra
    result = subprocess.run(command, cwd=PROJECT_ROOT, capture_output=True,
                            text=True)
    output = result.stdout + result.stderr
    (bed / "runner.log").write_text(
        (bed / "runner.log").read_text() if (bed / "runner.log").is_file() else "")
    with (bed / "runner.log").open("a") as handle:
        handle.write(f"$ {' '.join(command)}\n{output}\n")
    if expect_failure:
        return (None if result.returncode else "UNEXPECTED SUCCESS"), output
    if result.returncode:
        return None, output
    for line in output.splitlines():
        match = RUN_ID.match(line.strip())
        if match:
            return match.group(1), output
    return None, output


def runner_check(bed_root: Path, driver: Path, pfts_source: Path, ranks: int,
                 nyear: int, spinup: int, tables: tuple[str, ...]) -> list[str]:
    """Does the PRODUCTION runner continue a run, against a running model?

    THE PART OF WORLD-GQYP THAT WAS NEVER MEASURED. `--save-state` and
    `--continue-from` are wired through `run_lpj_guess.py`, the instruction file
    it derives, the run manifest and five refusals, and `--self-test` holds the
    ARITHMETIC against `framework/framework.cpp`. None of the PLUMBING had ever
    been exercised against a running model: whether the state directory is
    created where the manifest says it is, whether a continuation's manifest
    names its parent, whether the recorded state hashes match the files on disk
    at the moment of the continuation, and whether the continued record is the
    record the uninterrupted run has.

    Three runs, all on the same bed: a PARENT that saves, a CONTINUATION from
    it, and an uninterrupted CONTROL twice as long. The continuation's years
    have to be the control's years, exactly, row for row. Then a fourth run that
    has to be REFUSED: a continuation at a patch count the parent never ran is
    a different experiment resumed from someone else's state.
    """
    bed_root.mkdir(parents=True, exist_ok=True)
    pfts = bed_root / "runner_pfts.ins"
    reduced_pfts(pfts_source, pfts, spinup)

    failures: list[str] = []
    parent, output = invoke_runner(bed_root, driver, pfts, ranks, 1, nyear,
                                   ["--save-state", "--label",
                                    "runner_check_parent"])
    if parent is None:
        return [f"the parent run failed:\n{output[-2000:]}"]
    parent_dir = RUNS / parent
    manifest = json.loads((parent_dir / "run_manifest.json").read_text())
    saved = manifest.get("saved_state")
    if not saved:
        failures.append(
            f"{parent} was asked for --save-state and its manifest records no "
            "saved_state block")
    elif not (Path(saved["dir"]) / "meta.bin").is_file():
        failures.append(
            f"{parent}'s manifest names {saved['dir']} as its state directory "
            "and there is no state file in it")
    elif saved["covers_year"] != spinup + nyear - 1:
        failures.append(
            f"{parent}'s manifest says its state covers year "
            f"{saved['covers_year']}; a run of {nyear} retained years behind a "
            f"{spinup}-year spin-up ends at {spinup + nyear - 1}")

    child, output = invoke_runner(bed_root, driver, pfts, ranks, 1, nyear,
                                  ["--continue-from", parent, "--label",
                                   "runner_check_child"])
    if child is None:
        return failures + [f"the continuation failed:\n{output[-2000:]}"]
    child_manifest = json.loads((RUNS / child / "run_manifest.json").read_text())
    continuation = child_manifest.get("continuation")
    if not continuation:
        failures.append(f"{child} continues {parent} and its manifest records "
                        "no continuation block, so the chain has no provenance")
    elif continuation.get("parent") != parent:
        failures.append(
            f"{child}'s manifest names {continuation.get('parent')} as its "
            f"parent and it continued {parent}")

    control, output = invoke_runner(bed_root, driver, pfts, ranks, 1,
                                    2 * nyear, ["--label",
                                               "runner_check_control"])
    if control is None:
        return failures + [f"the uninterrupted control failed:\n{output[-2000:]}"]

    # The years the continuation is responsible for, in the control's own
    # numbering: it resumes where the parent stopped and adds its own record.
    first_year = spinup + nyear
    failures += compare(RUNS / control, RUNS / child, first_year, tables)

    # THE REFUSAL. A continuation at a patch count the parent never ran resumes
    # a different experiment from someone else's state, and the guard exists to
    # stop it. A guard that has never refused anything is not known to refuse.
    refused, output = invoke_runner(bed_root, driver, pfts, ranks, 3, nyear,
                                    ["--continue-from", parent, "--label",
                                     "runner_check_refusal"],
                                    expect_failure=True)
    if refused is not None:
        failures.append(
            "a continuation at npatch 3 from an npatch 1 parent was NOT "
            "refused")
    elif "npatch" not in output:
        failures.append(
            "the continuation at npatch 3 was refused and the refusal does not "
            f"name npatch, so it may have been refused for another reason:\n"
            f"{output[-1000:]}")
    return failures


def resume_and_save_instants(state: dict, year_length: int) -> dict:
    """The two instants `framework/framework.cpp` computes from a state block.

    Mirrored from lines 214-218 rather than inferred from the parameter names,
    and `self_test` refuses if those lines have moved out from under this. The
    resume instant is the FIRST day the resuming run simulates; the save instant
    is the LAST day the written state covers.
    """
    resume_year = state["state_year"] - 1 if state["state_day"] < 0 else state["state_year"]
    resume_day = year_length - 1 if state["state_day"] < 0 else state["state_day"]
    save_year = state["save_year"] - 1 if state["save_day"] < 0 else state["save_year"]
    save_day = year_length - 1 if state["save_day"] < 0 else state["save_day"]
    first_year, first_day = ((resume_year + 1, 0) if resume_day == year_length - 1
                             else (resume_year, resume_day + 1))
    return {"state_covers": (resume_year, resume_day),
            "first_simulated": (first_year, first_day),
            "state_written_covers": (save_year, save_day)}


def self_test() -> None:
    """Does the production runner's state block resume where its parent stopped?

    NO MODEL, NO FORCING. `run_lpj_guess.py` decides two integers -- which
    simulated year a saving run's state covers, and which one a continuing run
    resumes at -- and an off-by-one in either is the failure this whole row
    exists to avoid: a run that silently skips or repeats a simulated year while
    reporting the spin-up it was asked for. The identity has a right answer. The
    parent's state covers its own last simulated year, and the continuation's
    first simulated day is the day after that, with no day repeated and none
    skipped.

    It is checked through `run_lpj_guess.build_instruction`, so what is tested
    is the text the model will actually parse and not a second copy of the rule.

    THE SPIN-UP IS THE PART THAT BITES, and this test exists because it did.
    The model counts `date.year` from zero THROUGH `nyear_spinup`, so a run of
    `--nyear N` ends at simulated year `nyear_spinup + N - 1`. The first version
    of `--save-state` computed its save point from `nyear` alone and named a
    year in the middle of the run -- with the derived floor in
    `vesper_pfts.ins` that is thousands of simulated years early, and the state
    would have been written, accepted and continued from without anything
    objecting. So the spin-up here is READ FROM THE PFT FILE the runs actually
    import, never from a constant written into this test: a test that carries
    its own copy of the number cannot catch the number being wrong.
    """
    year_length = model_year_days(yaml.safe_load(CONFIG.read_text()))
    pfts = GENERATED / "vesper_pfts.ins"
    if not pfts.is_file():
        raise SystemExit(
            f"{rel(pfts)} is absent, so the spin-up this test has to reckon "
            "with is not known. Build it with build_vesper_pfts.py.")
    spinup = run_lpj_guess.spinup_years(pfts)

    framework = (GUESS_SOURCE / "framework" / "framework.cpp").read_text(
        encoding="utf-8")
    for expression in (
            "const int resume_year = state_day < 0 ? state_year - 1 : state_year;",
            "const int save_point_year = save_day < 0 ? save_year - 1 : save_year;"):
        if expression not in framework:
            raise SystemExit(
                "framework/framework.cpp no longer computes\n  "
                f"{expression}\nso resume_and_save_instants mirrors arithmetic "
                "the model does not do. Re-read framework.cpp and update it.")

    failures: list[str] = []
    parent_nyear, child_nyear = 1253, 1253
    parent_total = spinup + parent_nyear
    child_total = parent_total + child_nyear

    parent = {"restart": False, "save_state": True,
              "state_year": parent_total, "state_day": -1,
              "save_year": parent_total, "save_day": -1,
              "state_path": "/parent/state", "save_path": "/parent/state"}
    child = {"restart": True, "save_state": True,
             "state_year": parent_total, "state_day": -1,
             "save_year": child_total, "save_day": -1,
             "state_path": "/parent/state", "save_path": "/child/state"}

    parent_instants = resume_and_save_instants(parent, year_length)
    child_instants = resume_and_save_instants(child, year_length)

    # A run simulates years 0 .. nyear_spinup + nyear - 1, so the state a saving
    # run writes has to cover the last day of the last of those.
    if parent_instants["state_written_covers"] != (parent_total - 1, year_length - 1):
        failures.append(
            f"a run of {parent_nyear} retained years behind a {spinup}-year "
            f"spin-up writes a state covering "
            f"{parent_instants['state_written_covers']}, not the last day of "
            f"year {parent_total - 1}")
    if child_instants["first_simulated"] != (parent_total, 0):
        failures.append(
            f"the continuation's first simulated day is "
            f"{child_instants['first_simulated']}, not day 0 of year "
            f"{parent_total}: it would "
            + ("repeat" if child_instants["first_simulated"][0] < parent_total
               else "skip") + " simulated years")
    if child_instants["state_written_covers"] != (child_total - 1, year_length - 1):
        failures.append(
            f"the continuation writes a state covering "
            f"{child_instants['state_written_covers']}, not the last day of "
            f"year {child_total - 1}, so a third run could not continue it")

    # THE RUNNER'S OWN ARITHMETIC, not a restatement of it. The dicts above say
    # what a correct state block contains; these are the ones run_lpj_guess.py
    # actually builds from a parent manifest, and the two have to agree.
    built_parent = run_lpj_guess.state_block(
        spinup, parent_nyear, Path("/parent"), None)
    parent_manifest = {
        "run_id": "lpj_parent", "physical": {"nyear": parent_nyear,
                                             "nyear_spinup": spinup},
        "saved_state": {"covers_year": parent_total - 1},
    }
    built_child = run_lpj_guess.state_block(
        spinup, child_nyear, Path("/child"), parent_manifest,
        parent_state=Path("/parent/state"))
    for label, built, wanted in (("--save-state", built_parent, parent),
                                 ("--continue-from", built_child, child)):
        for key in ("state_year", "state_day", "save_year", "save_day"):
            if built[key] != wanted[key]:
                failures.append(
                    f"{label} builds {key} {built[key]}, and the instant a "
                    f"correct block names is {wanted[key]}")

    refusals = run_lpj_guess.continuation_refusals(
        parent_manifest, {}, {"nyear_spinup": spinup}, child_nyear)
    if any("save point" in r for r in refusals):
        failures.append(
            "continuation_refusals rejects a parent whose recorded save point "
            f"is the end of its own run: {refusals}")

    # The text the model parses, not a third copy of the rule.
    settings = {"title": "self_test", "nyear": child_nyear, "npatch": 1,
                "nyear_spinup": parent_total,
                "root_seed": 1, "nfix_a": 0.0, "nfix_b": 0.0, "ifbvoc": 0,
                "outputs": ("cmass.out",), "state": child,
                **wetland_gate.switches(False)}
    text = run_lpj_guess.build_instruction(
        {"driver": Path("/d"), "soilmap": Path("/s"), "pfts": Path("/p")},
        settings)
    for line in (f"state_year {parent_total}", "state_day -1",
                 f"save_year {child_total}", "restart 1", "save_state 1",
                 f"nyear_spinup {parent_total}", f"nyear {child_nyear}",
                 'state_path "/parent/state"', 'save_path "/child/state"'):
        if line not in text:
            failures.append(f"the generated instruction file omits {line!r}")
    if run_lpj_guess.build_instruction(
            {"driver": Path("/d"), "soilmap": Path("/s"), "pfts": Path("/p")},
            {**settings, "state": None}).count("state_path") != 0:
        failures.append("a run that neither saves nor restarts declares a "
                        "state_path, so the model would look for a state file")

    for failure in failures:
        print(f"  {failure}")
    if failures:
        raise SystemExit(
            f"\n{len(failures)} failures. A continuation that resumes at the "
            "wrong simulated year reports the spin-up it was asked for and "
            "integrates a different one.")
    print(f"restart instants agree over a {year_length}-day simulation year "
          f"and the {spinup}-year spin-up {rel(pfts)} declares: a run of "
          f"{parent_nyear} retained years writes a state covering year "
          f"{parent_total - 1}, a continuation resumes at day 0 of year "
          f"{parent_total} and writes its own covering year {child_total - 1}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true",
                        help="the state-block arithmetic against the model's "
                             "own, with no model and no forcing, and exit")
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
    parser.add_argument("--nyear-spinup", type=int, default=1,
                        help="the simulated years in front of --nyear. ONE by "
                             "default, which is the model's own minimum and is "
                             "what makes this fixture runnable at all: the PFT "
                             "file declares the DERIVED spin-up floor, which is "
                             "thousands of simulated years, and inheriting it "
                             "turns a twelve-year bed into a "
                             "twelve-thousand-year one. One also puts all but "
                             "the first simulated year in the output tables, "
                             "which is what the annual mode compares. It leaves "
                             "the CENTURY accelerator window empty, so this "
                             "fixture does not exercise the accelerator's own "
                             "restart state.")
    parser.add_argument("--runner", action="store_true",
                        help="instead: exercise the PRODUCTION runner's state "
                             "plumbing end to end -- a parent that saves, a "
                             "continuation from it, an uninterrupted control "
                             "twice as long, and a continuation at a patch "
                             "count the parent never ran, which has to be "
                             "refused. Needs --cells: it runs run_lpj_guess.py, "
                             "which runs the whole driver it is given.")
    parser.add_argument("--cells", type=str, default=None,
                        help="run only these cells, as lon,lat pairs separated "
                             "by ';' -- the coordinates the output tables "
                             "print. The driver is sliced into the bed and both "
                             "arms read the slice. Cells are independent in "
                             "LPJ-GUESS and this project's stochastic streams "
                             "are keyed by coordinate, so a cell integrates the "
                             "same trajectory either way; what a subset loses is "
                             "reach, so its verdict is written to a separate "
                             "report and does not gate a continuation.")
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

    if args.self_test:
        self_test()
        return

    # Every year here is a `date.year`, which counts from zero THROUGH the
    # spin-up: the run covers 0 .. nyear_spinup + nyear - 1 and
    # `commonoutput.cpp:785` writes no annual row before nyear_spinup. So the
    # split has to leave something on both sides AND land where the annual mode
    # has output to compare.
    last_year = args.nyear_spinup + args.nyear - 1
    if not args.nyear_spinup < args.state_year <= last_year:
        raise SystemExit(
            f"--state-year {args.state_year} has to fall inside the run and "
            f"after the spin-up: this run covers simulated years 0 to "
            f"{last_year} and writes output from {args.nyear_spinup}. A "
            "restart at the first output year continues nothing and one after "
            "the last compares nothing.")

    if sum((args.one_day, args.round_trip, args.runner)) > 1:
        raise SystemExit(
            "--one-day, --round-trip and --runner are different questions: what "
            "a simulated day did differently either side of a restart, what the "
            "write and read of a state file does not carry, and whether the "
            "production runner continues a run at all. Ask one.")

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
            "\n\nBuild the model and the forcing and run this again. This is an "
            "absent measurement, not a pass.")

    # THE MODEL'S OWN PRECONDITION ON THE SPIN-UP, asked here rather than
    # discovered from an abort one second into a run that has already taken the
    # host lock. `parameters.cpp:1402` refuses nyear_spinup <= freenyears, and
    # freenyears is 200 in the generated PFT file rather than the 100 the
    # shipped instruction files carry, so a bed sized from the wrong one parses
    # and then dies.
    freenyears = run_lpj_guess.declared_int(
        Path(args.pfts), "freenyears",
        why="the smallest spin-up this fixture may ask for is not known")
    if args.nyear_spinup <= freenyears:
        raise SystemExit(
            f"--nyear-spinup {args.nyear_spinup} is not above the "
            f"{freenyears} simulated years {rel(Path(args.pfts))} declares as "
            "freenyears, and the model refuses that outright: nitrogen "
            "limitation switches on after freenyears and the CENTURY "
            "accelerator window is derived from the difference. Ask for more "
            f"than {freenyears}.")

    # THE BINARY MUST CONTAIN THE MODEL IN THIS TREE. This fixture spawns
    # LPJ-GUESS and its verdict is about the serialization code that ran, so a
    # binary the tree has moved under makes the verdict describe source that
    # was never executed -- and continuity is precisely what a changed
    # serialization block breaks. Same refusal `run_lpj_guess.py` makes, from
    # the same reader. world-w62x.
    from build_lpj_guess import verify as verify_binary
    stale = verify_binary(GUESS_BINARY)
    if stale:
        raise SystemExit(
            "\n".join(stale) + "\n\nThe restart continuity fixture has NOT "
            "been run: its verdict would be about source this binary does not "
            "contain. Rebuild with biosphere/scripts/build_lpj_guess.py.")

    # The retained tables. The wetland set is included only when the gate has
    # granted activation, on exactly the terms run_lpj_guess.py includes it:
    # a table the run does not write cannot be compared, and asking for one
    # would make this fixture fail for a reason that is not a restart defect.
    declaration = wetland_gate.read_declaration()
    active = wetland_gate.granted(declaration)[0]
    tables = run_lpj_guess.OUTPUTS + (
        tuple(declaration["acceptance"]["retained_outputs"]) if active else ())

    bed_root = args.bed or (RUNS / "restart_continuity")
    bed_root.mkdir(parents=True, exist_ok=True)
    whole, resumed = bed_root / "whole", bed_root / "resumed"
    state_dir = bed_root / "state"
    whole_state = bed_root / "state_whole"
    for path in (whole, resumed, bed_root / "continued"):
        if path.exists():
            raise SystemExit(f"{path} exists; remove it before re-running")
    state_dir.mkdir(parents=True, exist_ok=True)

    cells = parse_cells(args.cells) if args.cells else None
    if cells is not None:
        if args.ranks > len(cells):
            raise SystemExit(
                f"--ranks {args.ranks} against {len(cells)} cells: a rank with "
                "no cells writes no output and the merge would be short. Ask "
                f"for at most {len(cells)} ranks.")
        sliced = bed_root / "driver_subset.bin"
        taken = slice_driver(Path(args.driver), cells, sliced)
        print(f"driver sliced to {len(taken)} cells: "
              + "; ".join(f"{lon},{lat}" for lon, lat in taken))
        args.driver = sliced

    def paths_for(bed: Path) -> dict:
        return {"driver": Path(args.driver).resolve(),
                "soilmap": Path(args.soilmap).resolve(),
                "pfts": (bed / "vesper_pfts.ins").resolve(),
                "pfts_source": Path(args.pfts).resolve()}

    settings = {"nyear": args.nyear, "nyear_spinup": args.nyear_spinup,
                "npatch": args.npatch,
                "root_seed": 20260828,
                "nfix_a": 0.234, "nfix_b": -0.172, "ifbvoc": 0,
                "outputs": tables, **wetland_gate.switches(active)}

    # The default split day. The middle of the simulation year is chosen before
    # any result is seen and for a stated reason: day 0 is where every annual
    # accumulator resets and the last day is where they flush, so both are days
    # on which a member that is lost the rest of the year looks carried.
    year_length = model_year_days(yaml.safe_load(CONFIG.read_text()))
    split_day = (args.state_day if args.state_day is not None
                 else year_length // 2)
    mode = ("runner" if args.runner
            else "round-trip" if args.round_trip
            else "one-day" if args.one_day else "annual")

    if mode == "runner":
        if cells is None:
            raise SystemExit(
                "--runner needs --cells. It invokes run_lpj_guess.py, which "
                "simulates every cell in the driver it is handed, and the "
                "runner's plumbing is not a function of how many cells ran.")
        # args.driver is the SLICE by now; args.pfts is still the source the
        # bed's reduced copy is made from.
        failures = runner_check(bed_root, Path(args.driver).resolve(),
                                Path(args.pfts).resolve(),
                                args.ranks, args.nyear, args.nyear_spinup,
                                run_lpj_guess.OUTPUTS)
        report = {
            "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "generator": "biosphere/scripts/verify_lpj_restart_continuity.py",
            "mode": mode, "nyear": args.nyear,
            "nyear_spinup": args.nyear_spinup,
            "state_year": None, "state_day": None,
            "ranks": args.ranks, "npatch": 1,
            "wetlands_active": active,
            "tables_compared": list(run_lpj_guess.OUTPUTS),
            "binary_sha256": run_lpj_guess.sha256(GUESS_BINARY),
            "cells": [[lon, lat] for lon, lat in taken],
            "continuous": not failures,
            "failures": failures,
        }
        out = GENERATED / "lpj_runner_continuation.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2) + "\n")
        print(f"production runner end to end: a parent of {args.nyear} retained "
              f"years behind a {args.nyear_spinup}-year spin-up, a continuation "
              f"of {args.nyear} more, and a {2 * args.nyear}-year control")
        for failure in failures:
            print(f"  {failure}")
        print(f"\nwrote {rel(out)}")
        if failures:
            print("\nThe production runner does not continue a run. "
                  "--continue-from must not be used to buy a record.")
            sys.exit(1)
        return


    if mode == "annual":
        run_bed(build_bed(whole, paths_for(whole),
                          {**settings, "title": "restart_continuity_whole"},
                          state_dir, args.state_year, 1, 0), args.ranks, tables)
        run_bed(build_bed(resumed, paths_for(resumed),
                          {**settings, "title": "restart_continuity_resumed"},
                          state_dir, args.state_year, 0, 1), args.ranks, tables)
        failures = (check_instants(whole, args.ranks, "the uninterrupted arm",
                                   year_length)
                    + check_instants(resumed, args.ranks, "the resumed arm",
                                     year_length)
                    + compare(whole, resumed, args.state_year, tables))
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
        failures = (check_instants(whole, args.ranks, "the writing arm",
                                   year_length)
                    + check_instants(resumed, args.ranks, "the rewriting arm",
                                     year_length)
                    + compare_states(state_dir, whole_state,
                                     "written", "rewritten"))
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
        failures = (check_instants(whole, args.ranks, "the uninterrupted arm",
                                   year_length)
                    + check_instants(resumed, args.ranks, "the splitting arm",
                                     year_length)
                    + check_instants(continued, args.ranks, "the continued arm",
                                     year_length)
                    + compare_states(whole_state, continued_state,
                                     "whole", "split"))
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
        "nyear": args.nyear, "nyear_spinup": args.nyear_spinup,
        "state_year": args.state_year,
        "state_day": split_day if mode != "annual" else None,
        "ranks": args.ranks, "npatch": args.npatch,
        "wetlands_active": active,
        "tables_compared": list(tables) if mode == "annual" else [],
        # WHICH EXECUTABLE THIS VERDICT IS ABOUT. A continuity verdict is a
        # statement about serialization code, so it travels with the binary
        # that ran it; `run_lpj_guess.py` refuses --continue-from unless a
        # verdict here says `continuous` for the binary it is about to run.
        "binary_sha256": run_lpj_guess.sha256(GUESS_BINARY),
        # WHICH CELLS THIS VERDICT IS ABOUT. Null is the whole driver, which is
        # the only thing the production gate accepts; a list is a subset bed
        # built to localise a defect and its report is written elsewhere, so a
        # cheap pass over four cells can never be read as a pass over the grid.
        "cells": ([[lon, lat] for lon, lat in taken]
                  if cells is not None else None),
        "continuous": not failures,
        "failures": failures,
    }
    out = (GENERATED / ("lpj_restart_continuity_subset.json" if cells is not None
                        else "lpj_restart_continuity.json"))
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
