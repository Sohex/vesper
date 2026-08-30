#!/usr/bin/env python3
"""Verify that a high-cadence segment's raw stream survives a failed conversion.

WORLDBUILDING CONTEXT: Vesper is a fictional planet and this is engineering
work on the simulation of it. Everything below is a property of
`continue_exoplasim.py`'s file handling or of `pyburn`'s request dispatch,
checked on fixtures. It runs no model, compiles nothing and spawns no process.

WHY IT EXISTS. A high-cadence orbit is 2 GB of raw at T21 and 15 GB at T42, it
costs a model run to produce, `.gitignore` keeps model output out of the
repository on purpose, and nothing else in the tree recovers it. So the one
property that matters is that the raw stream is still findable, and reported,
after the postprocessing that reads it has failed.

THE SECOND SUBJECT IS THE POSTPROCESSOR THE RESCUE CALLS. A rescued raw is
worth nothing if the conversion that reads it answers the wrong question, and
the conversion is one `pyburn.dataset` call over a list of variable codes. The
dispatch cases below hold that list to meaning one thing: a code, that code as
a string and the variable's name are one request and must come back identical,
every derivation the file carries can be asked for by the code it is written
for, and a request the raw cannot support is reported rather than raised or
silently dropped.

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

import hashlib
import json
import os
import re
import struct
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


# A T21 raw stream, which is the grid this project runs and the one whose
# record layout is transcribed below. The main header carries (nlon, nlat,
# nlev, ntru) in words 5 to 8; a spectral record is one per level and NSP long,
# and a gridpoint record is nlat*nlon long. `readallvariables` builds its
# ENTIRE time axis by counting code-139 records, so one per sample is what
# makes the file readable at all.
FIXTURE_NLAT, FIXTURE_NLON, FIXTURE_NLEV, FIXTURE_NTRU = 32, 64, 10, 21
FIXTURE_NSP = (FIXTURE_NTRU + 1) * (FIXTURE_NTRU + 2)


def _fortran_record(header: tuple[int, ...], payload) -> bytes:
    """One Fortran sequential record: a marker-wrapped header, then the data."""
    head = struct.pack("<8i", *header)
    body = np.asarray(payload, dtype="<f4").tobytes()
    return (struct.pack("<i", len(head)) + head + struct.pack("<i", len(head))
            + struct.pack("<i", len(body)) + body + struct.pack("<i", len(body)))


def synthetic_spectral_raw(path: Path, nsamples: int = 3,
                           seed: int = 20260830,
                           omit: tuple[int, ...] = ()) -> Path:
    """Write a raw stream pyburn can read, carrying the spectral wind source.

    WHY THIS IS SYNTHESISED RATHER THAN SLICED. The property under test is that
    key order does not change a derived field, and it holds for any spectral
    input at all -- the transform is linear and the derivation reads only the
    divergence and vorticity records. A slice of a real orbit would tie this
    file to `exoplasim/runs/`, which `.gitignore` keeps out of the repository,
    so the check would pass or fail on whether a 2 GB untracked artifact
    happened to be on the machine. Everything here is fixed by `seed`.

    The values are not physical and are not meant to be: nothing downstream
    reads them, and a check that needed them to be physical would be checking
    the model rather than pyburn's dispatch.

    `omit` drops the named codes from the stream. It exists for one case --
    a raw with no code-152 record, which is a file pyburn cannot read at all
    rather than a request it cannot answer -- and writing that as a subtraction
    from the readable fixture is what makes the two comparable.
    """
    rng = np.random.default_rng(seed)
    sigmah = np.linspace(0.05, 1.0, FIXTURE_NLEV)
    zsig = np.zeros(FIXTURE_NLAT * FIXTURE_NLON, dtype="f4")
    zsig[:FIXTURE_NLEV] = sigmah
    parts = [_fortran_record(
        (333, 0, 20260830, 0, FIXTURE_NLON, FIXTURE_NLAT,
         FIXTURE_NLEV, FIXTURE_NTRU), zsig)]
    for step in range(nsamples):
        spectral = lambda: rng.normal(0.0, 1.0e-5, FIXTURE_NSP)
        # Surface geopotential and log surface pressure, one spectral record
        # each; `dataset` reads the second before it reaches the request.
        if 129 not in omit:
            parts.append(_fortran_record(
                (129, 0, 20260830, step, FIXTURE_NSP, 1, 0, step), spectral()))
        lnps = np.zeros(FIXTURE_NSP)
        # Mode zero alone, so the surface is uniform at 1000 hPa. The
        # coefficient reaches the grid unscaled and pyburn wants log(Pa).
        lnps[0] = np.log(1.0e5)
        if 152 not in omit:
            parts.append(_fortran_record(
                (152, 0, 20260830, step, FIXTURE_NSP, 1, 0, step), lnps))
        for code in (130, 155, 138):  # air temperature, divergence, vorticity
            if code in omit:
                continue
            for _level in range(FIXTURE_NLEV):
                parts.append(_fortran_record(
                    (code, 1, 20260830, step, FIXTURE_NSP, 1, 0, step),
                    spectral()))
        parts.append(_fortran_record(
            (139, 0, 20260830, step, FIXTURE_NLON, FIXTURE_NLAT, 0, step),
            250.0 + rng.normal(0.0, 5.0, FIXTURE_NLAT * FIXTURE_NLON)))
    path.write_bytes(b"".join(parts))
    return path


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


class FakePyburn:
    """Enough of pyburn to exercise the re-conversion's decisions, and no more.

    The real one needs a Fortran stream and a spectral transform. What is under
    test here is not the transform -- `case_a_derived_first_request_converts`
    runs the real one on a synthetic raw for that -- it is WHICH CODES ARE
    ASKED FOR, what comes back, and what happens to a file already sitting
    under the name the conversion writes to. So this records the request and
    writes a netCDF holding what the request implies.
    """

    ilibrary = {"131": ["ua"], "132": ["va"], "259": ["spd"]}

    def __init__(self, samples: int, drop: str | None = None):
        self.samples = samples
        # A code the fake accepts and then does not produce, which is what the
        # real `dataset` does with a requested variable it can neither find in
        # the raw nor derive: it says so in the log and returns what it has.
        self.drop = drop
        self.requested: list[str] | None = None
        self.written: list[str] | None = None

    def dataset(self, filename, variablecodes, **kwargs):
        self.requested = list(variablecodes)
        data = {"time": [np.arange(self.samples), None]}
        for code in variablecodes:
            if str(code) == self.drop:
                continue
            data[self.ilibrary[str(code)][0]] = [
                np.zeros((self.samples, 2, 2)), None]
        return data

    def netcdf(self, rdataset, filename, logfile=None):
        # RETURNS THE OPEN DATASET, like the real one: pyburn.netcdf hands back
        # a Dataset it has not closed, and `postprocess` closes it on the next
        # line. A fake that closed it for the caller would let a caller that
        # forgets pass here and hash an unflushed file in production.
        self.written = sorted(k for k in rdataset if k != "time")
        nc = Dataset(filename, "w")
        nc.createDimension("time", self.samples)
        nc.createDimension("lat", 2)
        nc.createDimension("lon", 2)
        for name in self.written:
            nc.createVariable(name, "f4", ("time", "lat", "lon"))
        return nc


def _with_fake_pyburn(fake):
    """Put `fake` where `from exoplasim import pyburn` will find it."""
    import exoplasim as exo_module
    previous = getattr(exo_module, "pyburn", None)
    exo_module.pyburn = fake
    return exo_module, previous


def reconvert_fixture(tmp: Path, existing_samples: int | None,
                      samples: int = 1463):
    """A run whose raw is preserved and whose conversion may be a substitute."""
    run_dir = tmp / "run_fixture"
    (run_dir / "highcadence").mkdir(parents=True)
    (run_dir / "highcadence" / cx.high_cadence_raw_name(YEAR)).write_bytes(PAYLOAD)
    if existing_samples is not None:
        write_netcdf(run_dir / "highcadence" / f"MOST_HC.{YEAR:05d}.nc",
                     existing_samples)
    fake = FakePyburn(samples)
    exo_module, previous = _with_fake_pyburn(fake)
    try:
        record = cx.reconvert_high_cadence(
            run_dir, YEAR, samples, {"radius": 1.2, "gravity": 12.81,
                                     "gascon": 287.0})
    finally:
        exo_module.pyburn = previous
    return run_dir, fake, record


def case_the_substituted_conversion_is_kept_as_evidence() -> list[str]:
    """A twelve-bin substitute is moved aside, not overwritten.

    It is the only record of what the segment reported, and a consumer that
    read it read a twelve-bin average of the regular variable set. Deleting it
    would leave the corrected file and no way to say what had been there.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir, _fake, record = reconvert_fixture(Path(tmpdir), 12)
        aside = (run_dir / "highcadence" /
                 f"MOST_HC.{YEAR:05d}.{cx.SUBSTITUTED_SUFFIX}.nc")
        if not aside.is_file():
            problems.append(
                f"the substituted conversion is gone; it should be at "
                f"{aside.name}, because it is the evidence of what was reported")
        else:
            with Dataset(aside) as nc:
                if len(nc.dimensions["time"]) != 12:
                    problems.append(
                        "the file moved aside is not the twelve-bin one")
        if record.get("moved_aside") is None:
            problems.append("the record does not name the file moved aside, so "
                            "the substitution leaves no trace in the manifest")
        if not record.get("reconverted"):
            problems.append("the record does not say a conversion was run")
        with Dataset(run_dir / "highcadence" / f"MOST_HC.{YEAR:05d}.nc") as nc:
            if len(nc.dimensions["time"]) != 1463:
                problems.append("the corrected conversion did not take the "
                                "high-cadence file's name")
    return problems


def case_a_valid_conversion_is_not_redone() -> list[str]:
    """A file already holding the model's sample count is left alone.

    The raw is 2 GB and the transform is minutes; re-running it over a correct
    product would also replace a file whose hash something may already cite.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir, fake, record = reconvert_fixture(Path(tmpdir), 1463)
        if fake.requested is not None:
            problems.append("the conversion ran again over a file that already "
                            "held the samples the model wrote")
        if record.get("reconverted"):
            problems.append("the record claims a conversion that did not happen")
        if any((run_dir / "highcadence").glob(f"*{cx.SUBSTITUTED_SUFFIX}*")):
            problems.append("a valid conversion was moved aside as a substitute")
    return problems


def case_the_recorded_hash_describes_the_file() -> list[str]:
    """The sha256 in the record is the sha256 of the file on disk.

    `pyburn.netcdf` RETURNS the Dataset it opened and does not close it --
    `postprocess` closes it on the following line. A caller that hashes the
    path before that close records the digest of a file still missing its
    final metadata, and the recorded provenance then names a file that no
    longer exists in that form. That happened, on the real orbit 127.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir, _fake, record = reconvert_fixture(Path(tmpdir), None)
        out = run_dir / "highcadence" / f"MOST_HC.{YEAR:05d}.nc"
        digest = hashlib.sha256(out.read_bytes()).hexdigest()
        if record.get("sha256") != digest:
            problems.append(
                f"the record says {record.get('sha256')} and the file hashes "
                f"to {digest}. The conversion was hashed before it was closed")
        if record.get("bytes") != out.stat().st_size:
            problems.append("the recorded size is not the file's size either")
    return problems


def case_the_substitution_survives_a_second_pass() -> list[str]:
    """Re-running the rescue does not erase that a conversion was substituted.

    The second pass finds a valid conversion and has no way to derive that one
    was ever replaced, so a record rebuilt from scratch reports no
    substitution. That is the same trap `recovered_from` is guarded against,
    and it would quietly delete the only durable evidence every time anyone
    re-ran the command.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir, _fake, first = reconvert_fixture(Path(tmpdir), 12)
        if first.get("moved_aside") is None:
            problems.append("the first pass recorded no substitution to keep")
        fake = FakePyburn(1463)
        exo_module, previous = _with_fake_pyburn(fake)
        try:
            second = cx.reconvert_high_cadence(
                run_dir, YEAR, 1463, {"radius": 1.2, "gravity": 12.81,
                                      "gascon": 287.0})
        finally:
            exo_module.pyburn = previous
        if second.get("moved_aside") != first.get("moved_aside"):
            problems.append(
                f"the second pass reports moved_aside "
                f"{second.get('moved_aside')!r} where the first reported "
                f"{first.get('moved_aside')!r}. Re-running the rescue must not "
                f"erase that the conversion had been substituted")
    return problems


def case_the_request_is_the_declared_set() -> list[str]:
    """The conversion asks for HIGH_CADENCE_CODES, unpadded and unreordered.

    It used to prepend `pyburn.tscode` and drop that variable again, because
    `pyburn.dataset` closed its per-key loop logging a name that only the arm
    finding a code ALREADY IN THE RAW binds; `ua`, `va` and `spd` are each
    derived from the spectral divergence and vorticity, so a request opening
    with one of them raised UnboundLocalError before anything was written.
    That is world-2jj3, and `case_a_derived_first_request_converts` below holds
    the repair against the real pyburn. What is held here is that the caller no
    longer compensates for it: a workaround left standing beside a fix is a
    trap for whoever reads the two together.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        _run_dir, fake, _record = reconvert_fixture(Path(tmpdir), None)
        declared = [str(code) for code in cx.HIGH_CADENCE_CODES]
        if fake.requested != declared:
            problems.append(
                f"the request was {fake.requested} and the declared set is "
                f"{declared}. A code added to the head of the request, or an "
                f"order other than the declared one, is a workaround for a "
                f"defect that is fixed")
        if fake.written != sorted(FakePyburn.ilibrary[c][0] for c in declared):
            problems.append(
                f"the product carries {fake.written} and the declared set is "
                f"{sorted(FakePyburn.ilibrary[c][0] for c in declared)}. A "
                f"conversion that quietly adds a variable is the same class of "
                f"defect as one that quietly substitutes them")
    return problems


def case_a_missing_declared_field_is_refused() -> list[str]:
    """A product short of a declared field is refused rather than written.

    `pyburn.dataset` REPORTS a requested variable it cannot produce and returns
    what it has, which is right for the namelist a run postprocesses under --
    that lists every code the model can write, not the codes this
    configuration wrote. It is not right for a three-code request whose whole
    purpose is the gust distribution, so the completeness check belongs to the
    caller that knows what it asked for.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir = Path(tmpdir) / "run_fixture"
        (run_dir / "highcadence").mkdir(parents=True)
        (run_dir / "highcadence" /
         cx.high_cadence_raw_name(YEAR)).write_bytes(PAYLOAD)
        fake = FakePyburn(1463, drop="259")
        exo_module, previous = _with_fake_pyburn(fake)
        try:
            cx.reconvert_high_cadence(run_dir, YEAR, 1463,
                                      {"radius": 1.2, "gravity": 12.81,
                                       "gascon": 287.0})
        except RuntimeError as exc:
            if "spd" not in str(exc):
                problems.append(f"the refusal does not name the missing "
                                f"variable: {exc}")
        else:
            problems.append(
                "a conversion missing spd was accepted and written. That is "
                "the same class of substitution as the twelve-bin average: a "
                "file under the high-cadence name that is not the declared "
                "product")
        finally:
            exo_module.pyburn = previous
    return problems


def case_a_derived_first_request_converts() -> list[str]:
    """The real pyburn, on a request whose FIRST code is derived.

    THE ONE CASE HERE THAT RUNS THE ACTUAL SPECTRAL TRANSFORM, on a synthetic
    raw this file writes, because the property under test is pyburn's and a
    fake cannot hold it. The right answer is available without a model: key
    order is not an input to any derivation, so a request whose first code is
    derived must return exactly what the same request returns reordered. That
    is an identity, and it is what `world-2jj3` broke -- `dataset` closed its
    per-key loop reading a name that the `ua`, `va`, `spd` and `wa` arms never
    bind, so the derived-first request raised UnboundLocalError before writing
    anything and the reordered one did not.

    THE LOG IS CHECKED SEPARATELY, AND AGAINST A DIFFERENT REPAIR. Initialising
    the name ahead of the loop stops the traceback and leaves the line
    reporting whichever variable was collected before, which on the declared
    request alone reads as correct: every variable in it is dimensioned by time
    on axis 0, so a stale shape gives the right count under the right name. So
    the second request below adds a code pyburn can neither find in the raw nor
    derive. Nothing is stored for it, and under that repair the log announces
    it as collected at the previous variable's shape -- a line naming a variable
    the dataset does not contain. Requiring the log to name what the dataset
    carries, and only that, is what separates the two repairs.
    """
    from exoplasim import pyburn

    problems = []
    declared = [str(code) for code in cx.HIGH_CADENCE_CODES]
    names = [pyburn.ilibrary[code][0] for code in declared]
    marker = str(pyburn.tscode)
    with tempfile.TemporaryDirectory() as tmpdir:
        raw = synthetic_spectral_raw(Path(tmpdir) / "MOST_HC.00042")
        log = Path(tmpdir) / "pyburn.log"
        constants = {"radius": 1.2, "gravity": 12.81, "gascon": 287.0}
        try:
            derived_first = pyburn.dataset(
                str(raw), declared, mode="grid", zonal=False,
                substellarlon=180.0, physfilter=False, logfile=str(log),
                **constants)
        except Exception as exc:
            return [f"the declared request {declared} raised "
                    f"{type(exc).__name__}: {exc}. Its first code is derived "
                    f"from the spectral divergence and vorticity, and key "
                    f"order is not an input to any derivation"]
        reordered = pyburn.dataset(
            str(raw), [marker] + declared, mode="grid", zonal=False,
            substellarlon=180.0, physfilter=False, logfile=None, **constants)
        for name in names:
            if name not in derived_first:
                problems.append(f"the declared request produced no {name}")
            elif not np.array_equal(derived_first[name][0], reordered[name][0]):
                problems.append(
                    f"{name} differs between the declared request and the "
                    f"same request led by code {marker}. Key order changed a "
                    f"derived field")
        # A code that IS in pyburn's library, is absent from any high-cadence
        # raw, and has no derivation arm. Requesting it is legal and produces
        # nothing, which is the ordinary state of a dozen entries in every
        # namelist this project postprocesses under.
        underivable = "145"
        absent = pyburn.ilibrary[underivable][0]
        log.write_text("", encoding="utf-8")
        with_absent = pyburn.dataset(
            str(raw), declared + [underivable], mode="grid", zonal=False,
            substellarlon=180.0, physfilter=False, logfile=str(log),
            **constants)
        if absent in with_absent:
            return problems + [
                f"code {underivable} produced a {absent} out of a raw that "
                f"carries no such record. The fixture no longer holds the "
                f"property this case is built on"]
        collected = [line.split()[2] for line in
                     log.read_text(encoding="utf-8").split("\n")
                     if line.startswith("Collected variable:")]
        if collected != names:
            problems.append(
                f"the log reports {collected} collected and the dataset "
                f"carries {names}. A line naming a variable the dataset does "
                f"not hold is reporting another variable's shape under that "
                f"variable's name")
    return problems


PYBURN_SOURCE = (Path(__file__).resolve().parents[2] / "vendor" / "exoplasim"
                 / "exoplasim" / "pyburn.py")
MODEL_SOURCE = (Path(__file__).resolve().parents[2] / "vendor" / "exoplasim"
                / "exoplasim" / "plasim" / "src")


def _named(dataset: dict) -> dict:
    """The variables in a pyburn dataset, without the coordinate axes."""
    return {name: entry for name, entry in dataset.items()
            if name not in ("lat", "lon", "lev", "levp", "time")}


def _arm_inputs(function: str) -> tuple[dict[str, tuple[str, ...]], list[str]]:
    """The raw codes each derivation arm of `function` indexes, from the source.

    Read out of the text rather than out of the running module because the
    thing being checked is the DECLARATION against the CODE, and a table
    derived from the code it is supposed to constrain constrains nothing.
    Both spellings an arm uses are resolved: a literal `rawdata["142"]` and a
    `rawdata[str(divcode)]` through the module's own constant.

    THE PARSER'S OWN ASSUMPTION IS RETURNED RATHER THAN ASSUMED. Every arm
    dispatches on `key==str(<name>code)` against a module-level integer, and a
    reader of this function has to take that on trust; an arm written any other
    way would simply not be seen, and a check that cannot see a defect passes.
    So dispatch lines this cannot read come back as the second return value,
    and the cases that call it report them as failures of their own.
    """
    from exoplasim import pyburn

    lines = PYBURN_SOURCE.read_text(encoding="utf-8").split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith(f"def {function}("))
    end = next(i for i, l in enumerate(lines[start + 1:], start + 1)
               if l.startswith("def "))
    inputs: dict[str, set[str]] = {}
    unreadable: list[str] = []
    current: tuple[str, ...] = ()
    for line in lines[start:end]:
        opened = re.match(r"            (?:el)?if key==(.*?):", line)
        if opened:
            current = tuple(str(getattr(pyburn, name)) for name in
                            re.findall(r"str\((\w+code)\)", opened.group(1)))
            if not current:
                unreadable.append(f"{function}: {line.strip()}")
            for code in current:
                inputs.setdefault(code, set())
        if not current:
            # Everything before the first arm: the coordinate axes and the
            # already-present branch, which are not a derivation's inputs.
            continue
        for indexed in re.findall(r"rawdata\[([^\]]+)\]", line):
            through = re.fullmatch(r"str\((\w+code)\)", indexed)
            literal = re.fullmatch(r'"(\d+)"', indexed)
            if through:
                read = str(getattr(pyburn, through.group(1)))
            elif literal:
                read = literal.group(1)
            else:
                unreadable.append(f"{function}: rawdata[{indexed}]")
                continue
            for code in current:
                inputs[code].add(read)
    return ({code: tuple(sorted(reads, key=int))
             for code, reads in inputs.items()}, unreadable)


def _codes_the_model_writes() -> set[str]:
    """Every literal code the model hands to a record writer, from the source."""
    codes = set()
    for source in sorted(MODEL_SOURCE.glob("*.f90")):
        for call in re.finditer(r"call write(?:gp|sp|scalar)\((.*)\)\s*$",
                                source.read_text(encoding="utf-8"), re.M):
            depth, field, args = 0, "", []
            for char in call.group(1):
                if char == "(":
                    depth += 1
                elif char == ")":
                    depth -= 1
                if char == "," and depth == 0:
                    args.append(field.strip())
                    field = ""
                else:
                    field += char
            args.append(field.strip())
            if len(args) >= 3 and args[2].isdigit():
                codes.add(args[2])
    return codes


def case_the_three_spellings_of_a_request_agree() -> list[str]:
    """A code, that code as a string and the variable's name are one request.

    `dataset`'s docstring promises all three and its per-key loop used to
    resolve them in a branch on `type(key)==int` that normalised only the
    string side. An integer reached the derivation dispatch still an integer,
    compared unequal to every `key==str(...code)` arm, and produced NOTHING:
    `dataset(raw,[131,132,259])` returned an empty dataset where the same
    request spelled with strings returned `ua`, `va` and `spd`. world-hs5z.

    The identity needs no model and no expected values: the spelling of a
    request is not an input to any derivation, so all three must return the
    same variables with bit-identical arrays, or all three must refuse. It is
    asked of pyburn's WHOLE table rather than of the three codes the row
    named, because the defect was in the one line every request goes through.
    """
    from exoplasim import pyburn

    problems = []
    codes = list(pyburn.ilibrary.keys())
    with tempfile.TemporaryDirectory() as tmpdir:
        raw = str(synthetic_spectral_raw(Path(tmpdir) / "MOST_HC.00042"))
        constants = {"radius": 1.2, "gravity": 12.81, "gascon": 287.0,
                     "mode": "grid", "zonal": False, "substellarlon": 180.0,
                     "physfilter": False, "logfile": None}
        asked = {
            "the integer codes": [int(code) for code in codes],
            "the codes as strings": codes,
            "the variable names": [pyburn.ilibrary[code][0] for code in codes],
        }
        got = {}
        for how, request in asked.items():
            try:
                got[how] = _named(pyburn.dataset(raw, request, **constants))
            except Exception as exc:
                problems.append(f"the whole library requested as {how} raised "
                                f"{type(exc).__name__}: {exc}")
        if len(got) < len(asked):
            return problems
        reference = got["the codes as strings"]
        if not reference:
            return problems + ["the fixture produced no variables at all, so "
                               "this case can no longer tell the spellings apart"]
        for how, produced in got.items():
            if sorted(produced) != sorted(reference):
                problems.append(
                    f"requested as {how} the dataset carries "
                    f"{sorted(produced)} and requested as the codes as strings "
                    f"it carries {sorted(reference)}. The spelling of a "
                    f"request is not an input to any derivation")
                continue
            for name in reference:
                if not np.array_equal(produced[name][0], reference[name][0]):
                    problems.append(
                        f"{name} differs between the request spelled as {how} "
                        f"and the same request spelled as the codes as strings")
    return problems


def case_the_advanced_request_keeps_the_callers_key() -> list[str]:
    """`advancedDataset` reads each variable's options off the key it was given.

    Its per-variable options live in the CALLER's dict, keyed however the
    caller spelled the request. The resolution used to rebind the key to the
    numeric code and then look the options up under that, so the name-keyed
    form this function's own docstring shows raised `KeyError '139'`, and the
    integer form got past the lookup and fell through the derivation dispatch
    to return nothing. world-hs5z.

    Checked on an option that CHANGES THE ANSWER rather than on the call
    returning: `zonal` collapses the longitude axis, so a request that reached
    its options has one fewer dimension than a request that did not, and a
    stub dict would pass the first test and fail this one.
    """
    from exoplasim import pyburn

    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        raw = str(synthetic_spectral_raw(Path(tmpdir) / "MOST_HC.00042"))
        constants = {"radius": 1.2, "gravity": 12.81, "gascon": 287.0,
                     "mode": "grid", "substellarlon": 180.0,
                     "physfilter": False, "logfile": None}
        # One code present in the raw and one derived from it, so both branches
        # of the loop are asked. 139 is written by the model; 259 is not.
        asked = {
            "the variable names": {"ts": {"zonal": True}, "spd": {"zonal": True}},
            "the codes as strings": {"139": {"zonal": True}, "259": {"zonal": True}},
            "the integer codes": {139: {"zonal": True}, 259: {"zonal": True}},
        }
        plain = None
        for how, request in asked.items():
            try:
                got = _named(pyburn.advancedDataset(raw, request, **constants))
            except Exception as exc:
                problems.append(f"the request spelled with {how} raised "
                                f"{type(exc).__name__}: {exc}")
                continue
            if sorted(got) != ["spd", "ts"]:
                problems.append(
                    f"the request spelled with {how} produced {sorted(got)}, "
                    f"not the ts and spd it asked for")
                continue
            if plain is None:
                plain = _named(pyburn.advancedDataset(
                    raw, {"139": {}, "259": {}}, **constants))
            for name in ("ts", "spd"):
                if got[name][0].ndim >= plain[name][0].ndim:
                    problems.append(
                        f"{name} came back with {got[name][0].ndim} dimensions "
                        f"under zonal=True and {plain[name][0].ndim} without "
                        f"it, spelled with {how}. The per-variable options were "
                        f"not read off the caller's own key")
    return problems


def case_a_library_entry_added_at_run_time_is_addressable() -> list[str]:
    """A row put into `ilibrary` after import answers to its name as well.

    `run_exoplasim.register_energy_diagnostic_codes` makes codes 360-387 and
    460-487 requestable by inserting rows into `pyburn.ilibrary` at call time,
    which is the project's own way of teaching the postprocessor a code the
    vendored table does not carry. `slibrary` is built from `ilibrary` once at
    import and does not follow, so without a fallback the code answers and the
    name does not -- the same request meaning two things, which is the defect
    the resolution exists to remove.

    The row here is invented and removed again; it is a lookup fixture, not a
    claim that any such variable exists.
    """
    from exoplasim import pyburn

    problems = []
    code = "9001"
    if code in pyburn.ilibrary:                        # pragma: no cover
        return [f"code {code} is a real library entry now; this case needs "
                f"one that is not"]
    pyburn.ilibrary[code] = ["fixturevar", "a_lookup_fixture", "1"]
    try:
        by_code, meta_code = pyburn._resolvekey(code)
        by_name, meta_name = pyburn._resolvekey("fixturevar")
        if (by_code, meta_code) != (by_name, meta_name):
            problems.append(
                f"code {code} resolves to {(by_code, meta_code)} and its name "
                f"to {(by_name, meta_name)}. A row added after import is "
                f"addressable by one spelling and not the other")
    except Exception as exc:
        problems.append(f"a library row added after import does not resolve: "
                        f"{type(exc).__name__}: {exc}")
    finally:
        del pyburn.ilibrary[code]
    return problems


def case_every_derivation_is_addressable() -> list[str]:
    """Every derivation arm can be reached by the code it is written for.

    An arm dispatched on `key==str(<name>code)` whose code is not an
    `ilibrary` key is unreachable: the head of the per-key loop raises
    "Unknown variable code requested" before the dispatch, and the arm itself
    opens with `ilibrary[key]` and would raise again if it were reached.
    Codes 270, 271 and 275 were three such arms. world-e9yz.

    Both directions of the repair are held here. Deleting an arm satisfies
    this; so does adding the library row. What it refuses is the state that
    was there -- code carrying a derivation nothing can ask for.
    """
    from exoplasim import pyburn

    problems = []
    numbers = dict(re.findall(r"^(\w+code)\s*=\s*(\d+)",
                              PYBURN_SOURCE.read_text(encoding="utf-8"), re.M))
    for function in ("dataset", "advancedDataset"):
        found, unreadable = _arm_inputs(function)
        problems += [f"this check cannot read the dispatch at {line}, so an "
                     f"arm written that way is invisible to it"
                     for line in unreadable]
        for code in found:
            if code not in pyburn.ilibrary:
                named = sorted(n for n, v in numbers.items() if v == code)
                problems.append(
                    f"{function} has a derivation arm for code {code} "
                    f"({', '.join(named) or 'unnamed'}) and it is not an "
                    f"ilibrary key, so no request can reach it")
    names = [entry[0] for entry in pyburn.ilibrary.values()]
    for name in sorted({n for n in names if names.count(n) > 1}):
        problems.append(
            f"two library codes are both named {name}; a dataset is keyed on "
            f"the name, so one of them overwrites the other")
    return problems


def case_no_derivation_shadows_a_code_the_model_writes() -> list[str]:
    """No arm derives a code the model itself puts in the raw.

    When the model writes a code, the loop's already-present branch takes it
    and the arm below can never run, so an arm for such a code is dead. Worse,
    it is dead for a reason a reader cannot see: the arm and the library row
    disagree about what the code MEANS, and the library row is the one the
    output is labelled with.

    That was 268 and 269. `outmod.f90` writes glacier elevation as 268 and
    ground-plus-glacier elevation as 269 on every stream it writes, and
    `ilibrary` names them `icez` and `netz` to match; burn7 numbers shortwave
    and longwave net atmospheric radiation there, and the arms carried burn7's
    meaning. On a raw missing those records the arms fired and stored a
    radiative flux under the name `icez` in units of m2 s-2. world-e9yz.

    The raw is the authority for what a code means in it, which is what makes
    this checkable at all: the model source is the other half of the pair.
    """
    from exoplasim import pyburn

    written = _codes_the_model_writes()
    if not written:
        return ["no code was parsed as written by the model; the record "
                "writers moved or changed shape"]
    problems = []
    for function in ("dataset", "advancedDataset"):
        found, unreadable = _arm_inputs(function)
        problems += [f"this check cannot read the dispatch at {line}, so an "
                     f"arm written that way is invisible to it"
                     for line in unreadable]
        for code in sorted(found, key=int):
            if code in written:
                problems.append(
                    f"{function} derives code {code}, which the model writes "
                    f"into the raw itself. The arm is unreachable and it and "
                    f"ilibrary's {pyburn.ilibrary[code][0]} disagree about "
                    f"what code {code} is")
    return problems


def case_the_declared_derivation_inputs_are_the_arms_own() -> list[str]:
    """`_DERIVATION_INPUTS` says exactly what the arms read out of the raw.

    The table is what lets a request the raw cannot support be reported
    instead of raised, and it is only worth that if it tracks the arms. Both
    directions: a table entry naming a record the arm does not read suppresses
    a derivation that would have worked, and an arm reading a record the table
    does not name is the bare `KeyError` back again.
    """
    from exoplasim import pyburn

    problems = []
    for function in ("dataset", "advancedDataset"):
        found, unreadable = _arm_inputs(function)
        problems += [f"this check cannot read {line}, so what that arm reads "
                     f"out of the raw is not held to the table"
                     for line in unreadable]
        declared = {code: tuple(sorted(codes, key=int))
                    for code, codes in pyburn._DERIVATION_INPUTS.items()}
        for code in sorted(set(found) | set(declared), key=int):
            if code not in declared:
                problems.append(
                    f"{function} has an arm for code {code} and "
                    f"_DERIVATION_INPUTS does not, so its inputs are never "
                    f"checked and a raw without them raises KeyError")
            elif code not in found:
                problems.append(
                    f"_DERIVATION_INPUTS declares inputs for code {code} and "
                    f"{function} has no arm for it")
            elif found[code] != declared[code]:
                problems.append(
                    f"code {code}: {function} reads {found[code]} out of the "
                    f"raw and _DERIVATION_INPUTS declares {declared[code]}")
    return problems


def case_an_underivable_request_is_reported_not_raised() -> list[str]:
    """A derivation the raw cannot support says so and the conversion goes on.

    Every arm indexes `rawdata` directly, so a code the raw does not carry came
    back as a bare `KeyError('142')` from the middle of the per-key loop --
    naming neither the variable requested nor the reason, and taking every
    variable that would have followed with it. The namelist a run
    postprocesses under lists every code the model was COMPILED to be able to
    write, and the default request when no namelist is given is every key in
    `ilibrary`, so a request for a derivation this configuration cannot
    support is ordinary rather than exceptional. Same disposition as a code
    with no arm at all, for the same reason.

    Held on the pairing: the unsupported code is reported and SKIPPED, the
    codes around it still come back, and the report names the record that was
    missing. A repair that simply swallowed the error would pass the first
    two and fail the third.
    """
    from exoplasim import pyburn

    problems = []
    # 260 is precipitation, derived from codes 142 and 143. The fixture is a
    # spectral wind stream and carries neither.
    with tempfile.TemporaryDirectory() as tmpdir:
        raw = str(synthetic_spectral_raw(Path(tmpdir) / "MOST_HC.00042"))
        log = Path(tmpdir) / "pyburn.log"
        constants = {"radius": 1.2, "gravity": 12.81, "gascon": 287.0,
                     "mode": "grid", "zonal": False, "substellarlon": 180.0,
                     "physfilter": False}
        try:
            got = _named(pyburn.dataset(raw, ["139", "260", "259"],
                                        logfile=str(log), **constants))
        except Exception as exc:
            return [f"a request carrying code 260, whose inputs the raw does "
                    f"not carry, raised {type(exc).__name__}: {exc}. The "
                    f"request is ordinary and the codes beside it are "
                    f"derivable"]
        if "pr" in got:
            return ["code 260 produced a pr out of a raw that carries no "
                    "precipitation records. The fixture no longer holds the "
                    "property this case is built on"]
        for name in ("ts", "spd"):
            if name not in got:
                problems.append(
                    f"the request produced no {name}. A code the raw cannot "
                    f"support took the codes beside it down with it")
        report = [line for line in log.read_text(encoding="utf-8").split("\n")
                  if line.startswith("NOT COLLECTED") and " 260 " in line]
        if not report:
            problems.append(
                "nothing in the log names code 260. A derivation dropped "
                "without a word is the failure this replaced")
        elif not all(code in report[0] for code in ("142", "143")):
            problems.append(
                f"the report for code 260 is {report[0].strip()!r} and does "
                f"not name codes 142 and 143, which are the records it wanted")
    return problems


def case_the_ecological_codes_are_refused_out_loud() -> list[str]:
    """`refuse_eco_codes` fires, names the codes, and says where they are read.

    A refusal nobody has watched refuse is a comment. This is the pairing:
    a postprocessor code list carrying an ecological code is rejected and a
    list without one passes, so the check fails both if the refusal goes and
    if it starts refusing everything.

    `scripts/smoke_test.py:check_every_written_code_is_named_or_refused` holds
    the other half statically -- that codes 600 to 628 are DECLARED as refused
    rather than merely missing from pyburn's table. This holds that the
    declaration is executable: the refusal fires, names the code, and says
    where the stream IS read, so a reader who hits it is told what to do
    instead.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_exoplasim                                # noqa: E402

    problems = []
    ordinary = [139, 131, 132, 259]
    try:
        run_exoplasim.refuse_eco_codes(ordinary, "FIXTURE")
    except Exception as exc:
        problems.append(f"an ordinary code list was refused: "
                        f"{type(exc).__name__}: {exc}")
    carrying = ordinary + [610]
    try:
        run_exoplasim.refuse_eco_codes(carrying, "FIXTURE")
        problems.append(
            "a code list carrying 610 was accepted. That stream is written to "
            "its own unit with its own interval bounds and is read directly; "
            "asking pyburn for it stops the postprocessor naming neither the "
            "code nor the reason")
    except ValueError as exc:
        said = str(exc)
        if "610" not in said:
            problems.append(f"the refusal does not name the code it refused: "
                            f"{said}")
        if "compare_eco_streams" not in said:
            problems.append(
                f"the refusal does not say where the stream IS read, so a "
                f"reader is told no and not told what to do: {said}")
    return problems


def case_one_code_carries_one_record_shape() -> list[str]:
    """A code whose records disagree in length is refused, not reshaped.

    `readallvariables` keeps the FIRST record's header per code and
    `refactorvariable` reshapes the whole joined array by it, so two different
    fields sharing a code do not fail -- they come back as one variable with
    the wrong time axis and the wrong values, under the name the library
    supplies. That is the worst of the outcomes available, and it is what a
    diagnostic-block collision produces.

    The pairing is the point: two EXTRA records of the same length under one
    code must still read, because that is what every multi-level field and
    every extra output step is.
    """
    from exoplasim import pyburn

    problems = []
    scalar = _fortran_record((51, 0, 20260830, 0, 1, 1, 0, 0), np.zeros(1))
    spectral = _fortran_record((51, 0, 20260830, 0, FIXTURE_NSP, 1, 0, 0),
                               np.zeros(FIXTURE_NSP))
    with tempfile.TemporaryDirectory() as tmpdir:
        base = synthetic_spectral_raw(Path(tmpdir) / "MOST_HC.00042").read_bytes()

        same = Path(tmpdir) / "same.raw"
        same.write_bytes(base + scalar + scalar)
        try:
            pyburn.readfile(str(same))
        except Exception as exc:
            problems.append(
                f"two records of the SAME length under code 51 were refused: "
                f"{type(exc).__name__}: {exc}. That is what every multi-level "
                f"field looks like")

        mixed = Path(tmpdir) / "mixed.raw"
        mixed.write_bytes(base + scalar + spectral)
        try:
            pyburn.readfile(str(mixed))
            problems.append(
                "a code carrying a 1-long record and a 484-long record was "
                "read. The two are joined and reshaped by the first one's "
                "dimensions, so the result is a variable with the wrong time "
                "axis and the wrong values under a name the library supplies")
        except Exception as exc:
            said = str(exc)
            if "51" not in said:
                problems.append(f"the refusal does not name the code that "
                                f"collided: {said}")
    return problems


def case_the_readers_own_prerequisite_says_so() -> list[str]:
    """A raw with no log surface pressure is refused as unreadable, by name.

    Code 152 is read once, before any requested variable is looked at, and the
    half- and full-level pressures, the surface pressure, the horizontal
    gradients and every derivation that uses a pressure are built from it. So
    its absence is not a variable that could not be produced -- the disposition
    `_logcollected` and `_DERIVATION_INPUTS` settle for that -- it is a file
    that cannot be read at all, whatever was requested. It came back as a bare
    `KeyError('152')` raised from before the loop, which reads as a request
    having gone wrong.

    Held on the message rather than on the raising, because raising was never
    the problem: the refusal has to be distinguishable from the ordinary case
    by someone reading a log.
    """
    from exoplasim import pyburn

    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        raw = synthetic_spectral_raw(Path(tmpdir) / "MOST_HC.00042", omit=(152,))
        constants = {"radius": 1.2, "gravity": 12.81, "gascon": 287.0,
                     "mode": "grid", "zonal": False, "substellarlon": 180.0,
                     "physfilter": False, "logfile": None}
        try:
            pyburn.dataset(str(raw), ["139"], **constants)
            problems.append(
                "a raw carrying no code-152 record was read. Every pressure "
                "the postprocessor reports is built from it")
        except KeyError as exc:
            problems.append(
                f"a raw with no code-152 record raised a bare KeyError({exc}), "
                f"which names a dictionary key and not the file's defect")
        except Exception as exc:
            said = str(exc)
            if "152" not in said:
                problems.append(f"the refusal does not name code 152: {said}")
            if "prerequisite" not in said:
                problems.append(
                    f"the refusal does not separate a file that cannot be read "
                    f"from a variable that could not be produced: {said}")
    return problems


def case_the_planet_is_read_from_the_run() -> list[str]:
    """Radius comes back in Earth radii, which is what pyburn wants.

    `Model.configure` writes PLARAD in METRES and pyburn's `radius` is in Earth
    radii. Passing metres scales every wind by 6.371e6; that has happened, in
    `extract_high_cadence_wind.py`, and was caught only because the answer was
    compared against a climatology.
    """
    problems = []
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir = Path(tmpdir) / "run_fixture"
        run_dir.mkdir()
        (run_dir / "planet_namelist").write_text(
            " &planet_nl\n PLARAD = 7645464.0 \n GA = 12.81 \n"
            " GASCON = 287.0171275001432 \n /END\n", encoding="utf-8")
        got = cx.run_planet_constants(run_dir)
        if abs(got["radius"] - 1.2) > 1e-9:
            problems.append(
                f"radius came back as {got['radius']}, not 1.2 Earth radii. "
                f"PLARAD is metres and pyburn's radius is Earth radii")
        if got["gravity"] != 12.81 or got["gascon"] != 287.0171275001432:
            problems.append(f"gravity or gascon is not the run's: {got}")
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
        ("the substituted conversion is kept as evidence",
         case_the_substituted_conversion_is_kept_as_evidence),
        ("a conversion that already holds the model's samples is not redone",
         case_a_valid_conversion_is_not_redone),
        ("the conversion requests the declared set and nothing else",
         case_the_request_is_the_declared_set),
        ("a conversion missing a declared field is refused",
         case_a_missing_declared_field_is_refused),
        ("pyburn converts a request whose first code is derived",
         case_a_derived_first_request_converts),
        ("the recorded hash describes the file that was written",
         case_the_recorded_hash_describes_the_file),
        ("a re-run does not erase that the conversion was substituted",
         case_the_substitution_survives_a_second_pass),
        ("the conversion reads the planet off the run's own namelist",
         case_the_planet_is_read_from_the_run),
        ("a code, its string and its name are one request",
         case_the_three_spellings_of_a_request_agree),
        ("advancedDataset reads each variable's options off the caller's key",
         case_the_advanced_request_keeps_the_callers_key),
        ("a library row added at run time answers to its name too",
         case_a_library_entry_added_at_run_time_is_addressable),
        ("every derivation is reachable by the code it is written for",
         case_every_derivation_is_addressable),
        ("no derivation shadows a code the model writes itself",
         case_no_derivation_shadows_a_code_the_model_writes),
        ("the declared derivation inputs are the ones the arms read",
         case_the_declared_derivation_inputs_are_the_arms_own),
        ("a derivation the raw cannot support is reported, not raised",
         case_an_underivable_request_is_reported_not_raised),
        ("the ecological stream's codes are refused out loud",
         case_the_ecological_codes_are_refused_out_loud),
        ("one code carries one record shape, or the read is refused",
         case_one_code_carries_one_record_shape),
        ("the reader's own missing prerequisite says it is one",
         case_the_readers_own_prerequisite_says_so),
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
