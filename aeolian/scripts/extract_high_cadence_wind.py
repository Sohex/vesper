#!/usr/bin/env python3
"""Pull the near-surface wind out of a raw high-cadence file, in chunks.

    python aeolian/scripts/extract_high_cadence_wind.py \\
        exoplasim/runs/<run>/MOST_HC.00066

DUST-5 wants a gust distribution: the shape of the near-surface wind speed
histogram, sampled fast enough to resolve a 30-hour day. `continue_exoplasim.py
--high-cadence` produces the raw file, one sample every fourth timestep. This
turns it into a small netCDF holding the bottom model level alone.

Chunked, checkpointed and reporting progress throughout, because this is a
15 GB input and `notes/large-data.md` says that is not optional.

## Why this exists: pyburn cannot read a high-cadence file at all

Not "not efficiently". At all, and for a reason unrelated to size.

`pyburn.readallvariables` builds its time axis by counting records of code 139,
and reshapes every variable to `(ntimes, nlev, dim1)`. The high-cadence stream
writes its GRIDPOINT fields twice per timestep and its SPECTRAL fields once, so
`ntimes` comes out at twice the true sample count and every spectral variable
fails to reshape. That would have happened on the original file too, after
however many hours, so the run that was killed at 73 minutes was never going to
produce anything.

The winds are the spectral half: `outmod.f90:hcadencesp` writes divergence (155)
and vorticity (138), and `ua`, `va` and `spd` are all derived from them by a
transform. So this keeps the spectral records, plus exactly ONE code 139 per
timestep to serve as the time marker, and drops the rest. That fixes the axis
and shrinks a timestep from 350 records to 43 on the way.

## And it could not have been done in one pass either

`pyburn.readfile` opens with `fbuffer = fb.read()`. Measured 2026-08-17: 15.3 GB
of raw, 18.9 MB resident per sample, 27.6 GB after 73 minutes. That per-sample
figure is a full decoded T42 record at ~280 field-levels, so **pyburn holds the
entire decoded record set whatever you ask it to emit** -- restricting the output
codes bought nothing. Hence chunks, even after the filter above cuts a timestep
to 43 records.

Reimplementing the spectral transform to avoid all this would trade a resource
problem for a correctness risk, which is the wrong trade. This is the same
pyburn doing the same transform on a bounded working set.

## The raw format, as far as this needs it

Fortran unformatted sequential: every record is a 4-byte length, the payload,
then the same length again. A 32-byte record is a header of eight int32,

    (code, level, date, time, dim1, dim2, dim3, year_length)

and the record after it is that field's payload. The file opens with a single
grid-descriptor record, code 333, carrying nlon, nlat, nlev and the spectral
truncation; every chunk needs a copy of it or pyburn cannot size anything.
Timesteps are contiguous and are delimited by `(date, time)` changing.

## What comes out, and what happens if it dies

`<run>/highcadence/MOST_HC.NNNNN_wind.nc`, holding `spd`, `ua` and `va` at the
bottom model level only, on the full time axis. Roughly 290 MB against 15.3 GB
of raw, which is the point: the levels above the surface are not what a
saltation threshold reads.

Each chunk is written into that file as it completes and recorded in a sidecar
state file, so **re-running the same command resumes rather than restarting.**
Nothing is cleaned up until the last chunk lands.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import CONFIG, PROJECT_ROOT  # noqa: E402

# 139 is not wanted for itself. pyburn builds its ENTIRE time axis by counting
# occurrences of code 139 -- `readallvariables` appends to `variables["time"]`
# on every one -- and its assembly loop leaves `variable` unbound for any code
# in the file that is not also in this list. So 139 has to be both present in
# the chunk and requested here.
WIND_CODES = ["139", "131", "132", "259"]
FIELDS = ("spd", "ua", "va")
GRID_DESCRIPTOR_CODE = 333
HEADER_BYTES = 32

# Chosen against the measured 18.9 MB per decoded sample, so a chunk's working
# set lands near 1.5 GB. This is the ceiling the script exists to enforce, not a
# tuning knob: raise it and you are back to the failure it was written for.
DEFAULT_CHUNK_SAMPLES = 80


def rss_gb() -> float:
    try:
        with open(f"/proc/{os.getpid()}/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1e6
    except OSError:
        pass
    return float("nan")


def scan_records(path: Path):
    """Byte extents of the grid descriptor and of every sample.

    Returns `(preamble, samples)` where samples is a list of
    `[date, time, start, end]` bracketing every record belonging to one sample.
    """
    preamble = None
    main_header = None
    samples: list[list[int]] = []
    with open(path, "rb") as fh:
        while True:
            pos = fh.tell()
            raw = fh.read(4)
            if len(raw) < 4:
                break
            n = int(np.frombuffer(raw, "<i4")[0])
            body = fh.read(n)
            fh.read(4)
            if n != HEADER_BYTES:
                raise SystemExit(
                    f"{path}: expected a {HEADER_BYTES}-byte header at byte "
                    f"{pos}, found a {n}-byte record. Not a PlaSim raw file, "
                    "or truncated.")
            code, _lev, date, tim = (int(x) for x in np.frombuffer(body, "<i4")[:4])
            raw = fh.read(4)
            if len(raw) < 4:
                raise SystemExit(f"{path}: header at {pos} has no payload; truncated")
            m = int(np.frombuffer(raw, "<i4")[0])
            fh.seek(m + 4, 1)
            end = fh.tell()
            if code == GRID_DESCRIPTOR_CODE:
                preamble = (pos, end)
                main_header = np.frombuffer(body, "<i4")
                continue
            if samples and samples[-1][0] == date and samples[-1][1] == tim:
                samples[-1][3] = end
            else:
                samples.append([date, tim, pos, end])
    if preamble is None:
        raise SystemExit(f"{path}: no code {GRID_DESCRIPTOR_CODE} grid descriptor")
    return preamble, samples, main_header


TIME_MARKER_CODE = 139        # surface temperature, and pyburn's clock


def spectral_dim1(preamble_header) -> int:
    """A spectral record's first dimension, NESP, from the grid descriptor.

    Derived rather than declared: it is `(ntru + 1) * (ntru + 2)`, so it moves
    with the truncation and a hardcoded 1892 would silently be T42-only.
    """
    ntru = int(preamble_header[7])
    return (ntru + 1) * (ntru + 2)


def write_chunk(src: Path, dst: Path, preamble, group,
                spectral_dim1_value: int) -> None:
    """The grid descriptor, then the records that make a readable time axis.

    Keeps every spectral record, because that is where the winds are, and
    exactly one code 139 per timestep, because that is what pyburn counts to get
    `ntimes`. Dropping the second copy of each gridpoint field is the whole fix:
    with both present the axis comes out doubled and nothing reshapes.
    """
    with open(src, "rb") as fh, open(dst, "wb") as out:
        fh.seek(preamble[0])
        out.write(fh.read(preamble[1] - preamble[0]))
        for _date, _time, start, end in group:
            fh.seek(start)
            seen_marker = 0
            while fh.tell() < end:
                head_len = int(np.frombuffer(fh.read(4), "<i4")[0])
                head = fh.read(head_len)
                fh.read(4)
                body_len = int(np.frombuffer(fh.read(4), "<i4")[0])
                body = fh.read(body_len)
                fh.read(4)
                code, dim1 = (int(x) for x in np.frombuffer(head, "<i4")[[0, 4]])
                keep = dim1 == spectral_dim1_value
                if code == TIME_MARKER_CODE:
                    seen_marker += 1
                    keep = seen_marker == 1
                if keep:
                    n32, m32 = np.int32(head_len).tobytes(), np.int32(body_len).tobytes()
                    out.write(n32 + head + n32)
                    out.write(m32 + body + m32)


def open_output(path: Path, n_time: int, lat, lon) -> None:
    """Create the output file sized for the whole run, so chunks can fill it."""
    with Dataset(path, "w") as ds:
        ds.createDimension("time", n_time)
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        for name in FIELDS:
            v = ds.createVariable(name, "f4", ("time", "lat", "lon"), zlib=True)
            v.units = "m s-1"
            v.level = "bottom model level"
        ds.note = ("Bottom model level only, extracted in chunks because pyburn "
                   "holds the whole decoded record set in memory. See "
                   "aeolian/scripts/extract_high_cadence_wind.py.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", type=Path, help="a MOST_HC.NNNNN raw file")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--chunk", type=int, default=DEFAULT_CHUNK_SAMPLES)
    ap.add_argument("--restart", action="store_true",
                    help="discard the checkpoint and start over")
    args = ap.parse_args()

    args.raw = args.raw.resolve()
    if not args.raw.is_file():
        raise SystemExit(f"no raw file at {args.raw}")
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    planet = config["planet"]

    out = args.output or (args.raw.parent / "highcadence" / f"{args.raw.name}_wind.nc")
    out.parent.mkdir(parents=True, exist_ok=True)
    work = args.raw.parent / f".{args.raw.name}.hcwind"
    state_path = work / "state.json"
    if args.restart:
        shutil.rmtree(work, ignore_errors=True)
        out.unlink(missing_ok=True)
    work.mkdir(exist_ok=True)

    print(f"scanning {args.raw.name}  ({args.raw.stat().st_size / 1e9:.1f} GB)", flush=True)
    t0 = time.time()
    preamble, samples, main_header = scan_records(args.raw)
    nesp = spectral_dim1(main_header)
    groups = [samples[i:i + args.chunk] for i in range(0, len(samples), args.chunk)]
    offsets, at = [], 0
    for g in groups:
        offsets.append(at)
        at += len(g)
    print(f"  {len(samples)} samples in {len(groups)} chunks of up to {args.chunk}"
          f"  ({time.time() - t0:.0f} s to scan)", flush=True)

    state = {"done": []}
    if state_path.is_file() and out.is_file():
        state = json.loads(state_path.read_text())
        if state.get("samples") != len(samples):
            print("  checkpoint is for a different file; starting over", flush=True)
            state = {"done": []}
            out.unlink(missing_ok=True)
        elif state["done"]:
            print(f"  RESUMING: {len(state['done'])}/{len(groups)} chunks already done",
                  flush=True)

    import exoplasim.pyburn as pyburn
    # EARTH RADII, not metres. `Model.configure` passes `radius=radius_earth`
    # straight through to pyburn, whose own default is 1.0. Passing metres scales
    # every wind by 6.371e6 and yields a global mean of 4.7e7 m/s, which is
    # exactly how this was caught: the answer was checked against the snapshot
    # climatology's 7.356 m/s before being used for anything.
    radius_earth = float(planet["radius_earth"])
    started = time.time()
    todo = [k for k in range(len(groups)) if k not in state["done"]]

    for i, k in enumerate(todo):
        chunk_raw = work / f"chunk{k:04d}"
        chunk_nc = work / f"chunk{k:04d}.nc"
        write_chunk(args.raw, chunk_raw, preamble, groups[k], nesp)
        pyburn.postprocess(
            str(chunk_raw), str(chunk_nc), logfile=str(work / "pyburn.log"),
            variables=WIND_CODES, mode="grid",
            timeaverage=False, times=None, interpolatetimes=False,
            radius=radius_earth, gravity=float(planet["gravity_m_s2"]),
            gascon=287.0)
        with Dataset(chunk_nc) as ds:
            slab = {v: np.asarray(ds[v][:, -1, :, :], dtype="f4") for v in FIELDS}
            lat, lon = np.asarray(ds["lat"][:]), np.asarray(ds["lon"][:])
        if not out.is_file():
            open_output(out, len(samples), lat, lon)
        a = offsets[k]
        with Dataset(out, "a") as ds:
            for name in FIELDS:
                ds[name][a:a + slab[name].shape[0]] = slab[name]
        chunk_raw.unlink(missing_ok=True)
        chunk_nc.unlink(missing_ok=True)

        state["done"] = sorted(state["done"] + [k])
        state["samples"] = len(samples)
        state_path.write_text(json.dumps(state) + "\n", encoding="utf-8")

        elapsed = time.time() - started
        rate = (i + 1) / elapsed
        left = (len(todo) - i - 1) / rate if rate > 0 else float("nan")
        print(f"  chunk {len(state['done'])}/{len(groups)}  "
              f"{len(groups[k])} samples  "
              f"{elapsed / 60:.1f} min elapsed  ~{left / 60:.1f} min left  "
              f"rss {rss_gb():.1f} GB", flush=True)

    shutil.rmtree(work, ignore_errors=True)
    with Dataset(out, "a") as ds:
        try:
            ds.source_raw = str(args.raw.relative_to(PROJECT_ROOT))
        except ValueError:
            ds.source_raw = str(args.raw)
        ds.samples = len(samples)
    print(f"\nwrote {out}  ({out.stat().st_size / 1e6:.0f} MB, {len(samples)} samples)",
          flush=True)
    print("the raw file is NOT deleted; remove it once the fit is done", flush=True)


if __name__ == "__main__":
    main()
