"""Benchmark and verify pyburn's raw reader against a synthetic PlaSim output file.

WHY A SYNTHETIC FILE. The raw model output is deleted once pyburn has turned it
into netCDF, so by the time anyone wants to re-measure the reader there is
nothing left to measure it on, and producing a real one costs a model run. What
the reader's cost actually depends on is record count and record size, and both
are reproducible exactly: `--make` writes a file with the geometry of one T42
orbit at `NLOWIO = 0` -- 77 single-level codes and 36 ten-level codes over 185
output steps, 80,845 records of 128x64 float32, 2.65 GB.

The file is a valid Fortran sequential unformatted stream and reads correctly
through unmodified pyburn, which is what makes it usable as a bench rather than
a mock.

WHAT IT CHECKS. `--verify` reads the same file with the current reader and with
a reference copy of an older one, and requires every variable to agree in value,
shape AND dtype, plus the whole of `readfile()`'s refactored output. The right
answer is zero differing fields, so it is a test rather than a comparison of two
formulations: a reader change that alters any number fails it.

    python exoplasim/scripts/bench_pyburn_read.py --make /tmp/raw.bin
    python exoplasim/scripts/bench_pyburn_read.py --bench /tmp/raw.bin
    python exoplasim/scripts/bench_pyburn_read.py --verify /tmp/raw.bin --against OLD_PYBURN.py

`notes/audits/pyburn-postprocessing-cost.md` carries the measurements this
produced and what they changed.
"""

from __future__ import annotations

import argparse
import importlib.util
import struct
import sys
import time
from pathlib import Path

import numpy as np

import _paths  # noqa: F401  -- anchors paths and puts lib/ on sys.path

sys.path.insert(0, str(_paths.MODEL_SRC.parent))

# One T42 orbit at NLOWIO = 0. These are the dimensions the audit measured, not
# tuning knobs: change them and the numbers stop being comparable to the ones
# recorded there.
NLON, NLAT, NLEV, NTRU = 128, 64, 10, 42
NSTEPS = 185
CODES_2D = [139] + list(range(200, 276))   # 139 is what pyburn builds its time axis from
CODES_3D = list(range(300, 336))


def write_synthetic(path: Path, nsteps: int = NSTEPS) -> int:
    """Write a raw output file with the record structure pyburn expects.

    Each record is a Fortran sequential unformatted record with 4-byte markers:
    a marker giving the header length, eight int32 of header, the marker
    restated, then the same around the data payload.
    """
    ngrid = NLON * NLAT
    nrec = 0
    with open(path, "wb", buffering=1 << 22) as out:

        def record(header, databytes):
            hb = struct.pack("<8i", *header)
            marker = struct.pack("<i", 32)
            out.write(marker + hb + marker)
            dmarker = struct.pack("<i", len(databytes))
            out.write(dmarker + databytes + dmarker)

        # The main header record. Its dim1*dim2 has to equal its own word count,
        # because pyburn infers word length from the ratio of the two.
        zsig = np.zeros(ngrid, dtype="<f4")
        zsig[:NLEV] = np.linspace(0.1, 1.0, NLEV)
        record([0, 0, 0, 0, NLON, NLAT, NLEV, NTRU], zsig.tobytes())

        field = np.random.default_rng(0).random(ngrid).astype("<f4").tobytes()
        for t in range(nsteps):
            nstep = t * 32
            for code in CODES_2D:
                record([code, 0, 20260819, 0, NLON, NLAT, nstep, 0], field)
                nrec += 1
            for code in CODES_3D:
                # header[1] is the level index; pyburn reads a first record of 1
                # as the marker of a multi-level variable.
                for lev in range(1, NLEV + 1):
                    record([code, lev, 20260819, 0, NLON, NLAT, nstep, 0], field)
                    nrec += 1
    return nrec


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def bench(path: Path) -> None:
    from exoplasim import pyburn

    t = time.perf_counter()
    buf = open(path, "rb").read()
    print(f"file into memory     {time.perf_counter() - t:6.2f} s   ({len(buf) / 1e9:.2f} GB)")

    t = time.perf_counter()
    headers, variables = pyburn.readallvariables(buf)
    print(f"readallvariables     {time.perf_counter() - t:6.2f} s")
    print(f"  {len(variables) - 3} codes, {len(variables['time'])} times, "
          f"dtype {np.asarray(variables[str(CODES_3D[0])]).dtype}")

    t = time.perf_counter()
    pyburn.readfile(str(path))
    print(f"readfile             {time.perf_counter() - t:6.2f} s")


def verify(path: Path, against: Path) -> int:
    """Compare the current reader with a reference copy. Zero differences required."""
    from exoplasim import pyburn

    reference = _load("pyburn_reference", against)
    buf = open(path, "rb").read()

    t = time.perf_counter()
    ref_headers, ref_vars = reference.readallvariables(buf)
    t_ref = time.perf_counter() - t
    t = time.perf_counter()
    headers, variables = pyburn.readallvariables(buf)
    t_new = time.perf_counter() - t
    print(f"readallvariables   reference {t_ref:6.2f} s   current {t_new:6.2f} s"
          f"   {t_ref / t_new:.1f}x")

    failures = []
    if set(ref_vars) != set(variables):
        failures.append(f"variable set differs: {set(ref_vars) ^ set(variables)}")
    if ref_headers != headers:
        failures.append("headers differ")
    for key in set(ref_vars) & set(variables):
        a, b = np.asarray(ref_vars[key]), np.asarray(variables[key])
        if a.shape != b.shape:
            failures.append(f"{key}: shape {a.shape} vs {b.shape}")
        elif a.dtype != b.dtype:
            failures.append(f"{key}: dtype {a.dtype} vs {b.dtype}")
        elif not np.array_equal(a, b):
            failures.append(f"{key}: values differ")
    print(f"variables compared   {len(ref_vars)}")

    ref_data, data = reference.readfile(str(path)), pyburn.readfile(str(path))
    for key in set(ref_data) & set(data):
        a, b = np.asarray(ref_data[key]), np.asarray(data[key])
        if a.shape != b.shape or a.dtype != b.dtype or not np.array_equal(a, b):
            failures.append(f"readfile {key}: differs")
    print(f"readfile keys        {len(ref_data)}")

    if failures:
        print(f"\nFAIL: {len(failures)} differences")
        for line in failures[:10]:
            print("   " + line)
        return 1
    print("\nidentical in value, shape and dtype")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--make", metavar="PATH", type=Path,
                    help="write a synthetic raw file with one orbit's geometry")
    ap.add_argument("--steps", type=int, default=NSTEPS,
                    help=f"output steps to write (default {NSTEPS}, one orbit)")
    ap.add_argument("--bench", metavar="PATH", type=Path,
                    help="time the current reader on a raw file")
    ap.add_argument("--verify", metavar="PATH", type=Path,
                    help="compare the current reader against --against on a raw file")
    ap.add_argument("--against", metavar="PYBURN_PY", type=Path,
                    help="reference copy of pyburn.py to compare against")
    args = ap.parse_args()

    if args.make:
        nrec = write_synthetic(args.make, args.steps)
        size = args.make.stat().st_size
        print(f"{args.make}: {nrec} data records, {args.steps} steps, {size / 1e9:.2f} GB")
    if args.bench:
        bench(args.bench)
    if args.verify:
        if not args.against:
            ap.error("--verify needs --against pointing at a reference pyburn.py")
        return verify(args.verify, args.against)
    if not (args.make or args.bench or args.verify):
        ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
