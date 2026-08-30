#!/usr/bin/env python3
"""Verify that a high-cadence segment's raw stream survives a failed conversion.

WORLDBUILDING CONTEXT: Vesper is a fictional planet and this is engineering
work on the simulation of it. Everything below is a property of
`continue_exoplasim.py`'s file handling, checked on fixtures. It runs no model,
compiles nothing and spawns no process.

WHY IT EXISTS. A high-cadence orbit is 2 GB of raw at T21 and 15 GB at T42, it
costs a model run to produce, `.gitignore` keeps model output out of the
repository on purpose, and nothing else in the tree recovers it. So the one
property that matters is that the raw stream is still findable, and reported,
after the postprocessing that reads it has failed.

WHAT WENT WRONG, and each case below is one half of it. On run_67323a923013's
orbit 127, pyburn raised on the raw, ExoPlaSim's `integritycheck` quietly
re-ran it under `example.nl` and left a twelve-bin average of the REGULAR
variable set at `highcadence/MOST_HC.00127.nc`, and the segment reported
success. The runner then printed `high-cadence raw kept: []` while the 2.0 GB
raw sat at the run root, because it globbed `highcadence/` for a file the model
never puts there. The extractor later read 1463 samples out of that raw, which
is the number `plasim.f90` says it wrote and twelve is not.

Run it:

    python exoplasim/scripts/verify_high_cadence_rescue.py

Exits 0 when every case holds and 1 on the first that does not.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))

import continue_exoplasim as cx                        # noqa: E402


YEAR = 42
STAMP = {"run_id": "run_fixture", "source_build": "fixture-build",
         "segment_start_year_index": YEAR, "orbits": 1, "interval_steps": 4,
         "runsteps_per_orbit": 5850, "expected_samples_per_orbit": 1463,
         "executable_sha256": None, "purpose": "diagnostic",
         "started_utc": "2026-01-01T00:00:00+00:00"}
PAYLOAD = b"not a real Fortran stream, but it is the bytes that must survive\n" * 64


def fixture_run(tmp: Path, raw_in: str) -> Path:
    """A run directory with the raw stream placed where `raw_in` says.

    `root` is where the model leaves it (`__init__.py:891`), which is also
    where it stays on the SUCCESS path because the move into `highcadence/`
    matches only the netCDF. `crashed` is where `_crash()` puts the whole
    working directory when the postprocessor's exception propagates.
    """
    run_dir = tmp / "run_fixture"
    (run_dir / "highcadence").mkdir(parents=True)
    name = cx.high_cadence_raw_name(YEAR)
    if raw_in == "root":
        (run_dir / name).write_bytes(PAYLOAD)
    elif raw_in == "crashed":
        crashed = tmp / "run_fixture_crashed"
        crashed.mkdir()
        (crashed / name).write_bytes(PAYLOAD)
    else:                                              # pragma: no cover
        raise ValueError(raw_in)
    return run_dir


def write_netcdf(path: Path, times: int) -> None:
    with Dataset(path, "w") as ds:
        ds.createDimension("time", times)
        ds.createDimension("lat", 2)
        ds.createDimension("lon", 2)
        var = ds.createVariable("spd", "f4", ("time", "lat", "lon"))
        # netCDF4 reshapes the array it is handed, which numpy 2.5 deprecates.
        # Scoped to this line because it is the library's assignment and not
        # this file's, and a blanket filter in a check hides the next one.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            var[:] = np.zeros((times, 2, 2), dtype="f4")


def case_sample_count_is_the_models_own() -> list[str]:
    """The expected sample count reproduces plasim.f90's write condition exactly.

    A check whose right answer is wrong is worse than no check, so this does not
    trust the closed form: it counts the timesteps the model's own condition
    would fire on, over a grid of run lengths, orbit counts and intervals.

    The one anchor from a real run is pinned separately. `runsteps_per_orbit`
    5850 at interval 4 must give 1463, which is what
    `aeolian/scripts/extract_high_cadence_wind.py` actually pulled out of
    run_67323a923013's MOST_HC.00127.
    """
    problems = []

    def brute(runsteps: int, orbits: int, interval: int) -> int:
        # plasim.f90:776 nhcstp = 1 at the top of each call, :896 increments it,
        # :841-843 write when nhcstp >= hcstartstep .and. nhcstp < hcendstep
        # .and. mod(nhcstp - hcstartstep, hcinterval) == 0.
        start, end = cx.HIGH_CADENCE_START_STEP, runsteps * orbits
        return sum(1 for n in range(1, runsteps + 1)
                   if start <= n < end and (n - start) % interval == 0)

    for runsteps in (12, 100, 5850, 11700):
        for orbits in (1, 2, 5):
            for interval in (1, 3, 4, 7):
                want = brute(runsteps, orbits, interval)
                got = cx.expected_high_cadence_samples(runsteps, orbits, interval)
                if got != want:
                    problems.append(
                        f"expected_high_cadence_samples({runsteps}, {orbits}, "
                        f"{interval}) is {got}; counting the model's own write "
                        f"condition gives {want}")
    if cx.expected_high_cadence_samples(5850, 1, 4) != 1463:
        problems.append(
            "the T21 orbit this project has actually run wrote 1463 samples at "
            "interval 4 and the formula no longer says so, so the check has "
            "come loose from the one measurement that anchors it")
    return problems


def case_the_conversion_substitution_is_refused() -> list[str]:
    """A twelve-bin file under the high-cadence name is refused; a full one is not.

    Twelve is not an arbitrary wrong number: it is the regular output's bin
    count, and it is what ExoPlaSim's integritycheck fallback leaves behind when
    pyburn fails on the raw and it re-runs under `example.nl`. This is the case
    that turns that silent substitution into a refusal.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        substituted = tmp / "MOST_HC.00042.nc"
        write_netcdf(substituted, 12)
        try:
            cx.validate_high_cadence(substituted, 1463)
        except RuntimeError as exc:
            if "12 time samples" not in str(exc):
                problems.append(f"the refusal does not name what it found: {exc}")
        else:
            problems.append(
                "a twelve-bin file passed as the high-cadence product. That is "
                "the integritycheck fallback's output wearing the right name, "
                "and it is the failure this check exists for")

        good = tmp / "MOST_HC.00043.nc"
        write_netcdf(good, 1463)
        try:
            cx.validate_high_cadence(good, 1463)
        except RuntimeError as exc:
            problems.append(f"a file with the model's own sample count was "
                            f"refused: {exc}")

        missing = tmp / "MOST_HC.00044.nc"
        try:
            cx.validate_high_cadence(missing, 1463)
        except RuntimeError:
            pass
        else:
            problems.append("an absent conversion passed the check")
    return problems


def _check_rescue(run_dir: Path, where: str) -> list[str]:
    problems = []
    found_years = cx.discover_high_cadence_years(run_dir)
    if found_years != [YEAR]:
        problems.append(f"with the raw {where}, discovery returned "
                        f"{found_years} and the fixture holds orbit {YEAR}")
    entries = cx.preserve_high_cadence_raw(run_dir, range(YEAR, YEAR + 1), STAMP)
    if len(entries) != 1:
        return [f"with the raw {where}, the rescue found {len(entries)} "
                f"streams and the fixture has exactly one. An empty list is "
                f"the defect this exists to remove."]
    entry = entries[0]
    landed = run_dir / "highcadence" / cx.high_cadence_raw_name(YEAR)
    if not landed.is_file():
        problems.append(f"with the raw {where}, the rescue reported a stream "
                        f"that is not at {landed}")
        return problems
    if landed.read_bytes() != PAYLOAD:
        problems.append(f"with the raw {where}, the preserved file is not the "
                        f"bytes the model wrote")
    if entry["bytes"] != len(PAYLOAD):
        problems.append(f"with the raw {where}, the record says "
                        f"{entry['bytes']} bytes and the file is {len(PAYLOAD)}")
    if Path(entry["path"]).name != landed.name:
        problems.append(f"with the raw {where}, the reported path "
                        f"{entry['path']} does not name the file that exists")
    for key in ("sha256", "run_id", "source_build", "expected_samples_per_orbit"):
        if not entry.get(key):
            problems.append(f"with the raw {where}, the record carries no "
                            f"{key}; an unstamped path is a file nobody can "
                            f"attribute to a run")
    sidecar = run_dir / "highcadence" / "high_cadence_raw.json"
    if not sidecar.is_file():
        problems.append(f"with the raw {where}, no high_cadence_raw.json was "
                        f"written beside the stream. The run manifest is the "
                        f"other record and `_crash()` moves it away, so the "
                        f"sidecar is what survives that")
    elif json.loads(sidecar.read_text())[0]["file"] != landed.name:
        problems.append(f"with the raw {where}, the sidecar does not name the "
                        f"preserved stream")
    # The rescue runs on both paths and may run twice; a second call must not
    # lose the file it already moved.
    again = cx.preserve_high_cadence_raw(run_dir, range(YEAR, YEAR + 1), STAMP)
    if len(again) != 1 or not landed.is_file():
        problems.append(f"with the raw {where}, rescuing a second time lost "
                        f"the stream. The failure path calls it after the "
                        f"success path already has")
    return problems


def case_raw_at_the_run_root_is_found() -> list[str]:
    """The success path: the model leaves the raw at the run root and always has.

    The control is the glob this replaced. `highcadence/` holds the netCDF and
    never the raw, because ExoPlaSim's move matches `MOST_HC.NNNNN*.nc`, so a
    search there for an extensionless file returns nothing on a segment that
    worked perfectly. That is what printed an empty list.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir = fixture_run(Path(tmpdir), "root")
        write_netcdf(run_dir / "highcadence" / f"MOST_HC.{YEAR:05d}.nc", 1463)
        superseded = [p for p in (run_dir / "highcadence").glob("MOST_HC.[0-9]*")
                      if p.suffix == ""]
        if superseded:
            problems.append(
                "the superseded glob found the raw in highcadence/, so this "
                "case no longer reproduces the layout the model produces and "
                "proves nothing about the report that was empty")
        problems += _check_rescue(run_dir, "at the run root")
    return problems


def case_raw_in_the_crashed_sibling_is_found() -> list[str]:
    """The failure path: `_crash()` has moved the whole run directory away.

    `exoplasim/__init__.py:1509` is `mv <run>/* <run>_crashed/`, so after a
    postprocessor exception the raw is neither where the model left it nor
    where anything looked, the run directory is empty, and the only record of
    the stream's location is a sibling directory name. This puts it back under
    the run id, which is the handle everything else uses.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir = fixture_run(Path(tmpdir), "crashed")
        if (run_dir / cx.high_cadence_raw_name(YEAR)).exists():
            problems.append("the fixture left a raw in the run directory; "
                            "`_crash()` empties it and this case is about "
                            "what survives that")
        problems += _check_rescue(run_dir, "in the crashed sibling")
        sidecar = run_dir / "highcadence" / "high_cadence_raw.json"
        if not sidecar.is_file():
            return problems
        entry = json.loads(sidecar.read_text())[0]
        if "crashed" not in entry.get("recovered_from", ""):
            problems.append(
                "the record does not say the stream was recovered from the "
                "crashed sibling. Where a file came from is half of what makes "
                "the rescue auditable")
    return problems


def case_nothing_found_is_said_out_loud() -> list[str]:
    """A segment that asked for a high-cadence stream and has none says so.

    The empty list was indistinguishable from "there is nothing to report", and
    that is exactly what made a 2 GB file look absent for two days.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir = Path(tmpdir) / "run_fixture"
        (run_dir / "highcadence").mkdir(parents=True)
        entries = cx.preserve_high_cadence_raw(
            run_dir, range(YEAR, YEAR + 1), STAMP)
        if entries:
            problems.append("the rescue invented a stream where there is none")
        if (run_dir / "highcadence" / "high_cadence_raw.json").exists():
            problems.append("a sidecar was written for no stream, so a later "
                            "reader would find a record of nothing")
    return problems


def main() -> None:
    cases = [
        ("the expected sample count is the model's own write condition",
         case_sample_count_is_the_models_own),
        ("a substituted twelve-bin conversion is refused",
         case_the_conversion_substitution_is_refused),
        ("the raw at the run root is found and stamped",
         case_raw_at_the_run_root_is_found),
        ("the raw in the crashed sibling is found and stamped",
         case_raw_in_the_crashed_sibling_is_found),
        ("no stream is reported as no stream, not as an empty list",
         case_nothing_found_is_said_out_loud),
    ]
    failed = 0
    for name, run in cases:
        problems = run()
        if problems:
            failed += 1
            print(f"[ FAIL ] {name}", flush=True)
            for p in problems:
                print(f"         {p}", flush=True)
        else:
            print(f"[  ok  ] {name}", flush=True)
    print(f"\n{len(cases)} cases, {failed} failed")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
