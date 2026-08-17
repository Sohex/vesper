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

## Why this exists rather than calling pyburn once

Calling `pyburn.postprocess` on the whole file does not work at this cadence.
Measured 2026-08-17 on a single T42 orbit: 15.3 GB of raw, 1462 samples, and
pyburn had reached 27.6 GB resident after 73 minutes with nothing written before
it was killed. The cost is 18.9 MB per sample, which is very close to a full
decoded T42 record at ~280 field-levels, so **pyburn holds the entire decoded
record set in memory whatever you ask it to emit.** Restricting the output codes
to three winds bought nothing at all.

The fix is not to reimplement pyburn. The raw file carries SPECTRAL divergence
and vorticity, not gridpoint winds -- `outmod.f90:hcadencesp` writes codes 155
and 138 -- so `ua`, `va` and `spd` are all derived through a spectral transform,
and hand-rolling that to save memory would trade a resource problem for a
correctness risk. Instead this splits the raw file on timestep boundaries and
hands pyburn one chunk at a time: the same code doing the same transform with a
bounded working set.

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

WIND_CODES = ["131", "132", "259"]
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
                continue
            if samples and samples[-1][0] == date and samples[-1][1] == tim:
                samples[-1][3] = end
            else:
                samples.append([date, tim, pos, end])
    if preamble is None:
        raise SystemExit(f"{path}: no code {GRID_DESCRIPTOR_CODE} grid descriptor")
    return preamble, samples


def write_chunk(src: Path, dst: Path, preamble, group) -> None:
    """The grid descriptor, then one contiguous run of samples, copied blockwise."""
    with open(src, "rb") as fh, open(dst, "wb") as out:
        fh.seek(preamble[0])
        out.write(fh.read(preamble[1] - preamble[0]))
        start, end = group[0][2], group[-1][3]
        fh.seek(start)
        remaining = end - start
        while remaining > 0:
            block = fh.read(min(1 << 24, remaining))
            if not block:
                break
            out.write(block)
            remaining -= len(block)


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
    preamble, samples = scan_records(args.raw)
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
    radius_m = float(planet["radius_earth"]) * 6371000.0
    started = time.time()
    todo = [k for k in range(len(groups)) if k not in state["done"]]

    for i, k in enumerate(todo):
        chunk_raw = work / f"chunk{k:04d}"
        chunk_nc = work / f"chunk{k:04d}.nc"
        write_chunk(args.raw, chunk_raw, preamble, groups[k])
        pyburn.postprocess(
            str(chunk_raw), str(chunk_nc), logfile=str(work / "pyburn.log"),
            variables=WIND_CODES, mode="grid",
            timeaverage=False, times=None, interpolatetimes=False,
            radius=radius_m, gravity=float(planet["gravity_m_s2"]), gascon=287.0)
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
        ds.source_raw = str(args.raw.relative_to(PROJECT_ROOT))
        ds.samples = len(samples)
    print(f"\nwrote {out}  ({out.stat().st_size / 1e6:.0f} MB, {len(samples)} samples)",
          flush=True)
    print("the raw file is NOT deleted; remove it once the fit is done", flush=True)


if __name__ == "__main__":
    main()
