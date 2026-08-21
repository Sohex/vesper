#!/usr/bin/env python
"""Compare two model restarts record by record, and say what kind of difference it is.

A change to the MPI decomposition cannot be checked by a checksum: each process
sums the Legendre transform over the latitudes it holds and the partial sums are
completed by a reduce-scatter, so changing WHICH latitudes a process holds
regroups a floating point sum and the last bits move. What the check has to
distinguish is a regrouped sum from a different computation, and those differ by
many orders of magnitude, not by a factor.

So every record is reported in one of three states:

  identical   byte for byte. Fields that are only scattered in and gathered
              back out -- the land-sea mask, the orography -- MUST land here,
              because no arithmetic touches them and a permutation applied in
              one direction and not the other shows up nowhere else.
  rounding    differs, but by no more than `--tol` relative to the record's own
              RMS. That is what regrouping a sum of NLAT terms produces.
  DIFFERENT   differs by more than that.

`plasim_status` is sequential unformatted Fortran: a 16-character name record
then a data record, repeatedly. Reals are 8 byte because the model is compiled
with -fdefault-real-8; a data record whose length is not a multiple of 8 is
read as 4-byte integers.
"""
import argparse
import struct
import sys
from pathlib import Path


def records(path: Path):
    """(name, payload) pairs, in file order."""
    raw = path.read_bytes()
    out, pos, name = [], 0, None
    while pos < len(raw):
        (n,) = struct.unpack_from("<i", raw, pos)
        body = raw[pos + 4:pos + 4 + n]
        (m,) = struct.unpack_from("<i", raw, pos + 4 + n)
        if m != n:
            raise SystemExit(f"{path}: record markers disagree at byte {pos} "
                             f"({n} vs {m}); this is not a sequential "
                             f"unformatted file written by this model")
        pos += 8 + n
        if n == 16 and all(32 <= c < 127 for c in body):
            name = body.decode().strip()
            continue
        out.append((name if name is not None else "<unnamed>", body))
        name = None
    return out


def as_numbers(body: bytes):
    if len(body) % 8 == 0:
        return struct.unpack(f"<{len(body)//8}d", body), "real8"
    return struct.unpack(f"<{len(body)//4}i", body), "int4"


def relative(a, b) -> float:
    """Largest elementwise difference, scaled by the record's own RMS.

    Scaled by the RMS rather than elementwise, because a field with values
    passing through zero makes an elementwise relative difference meaningless
    and reports 1.0 for a perfectly good rounding difference.
    """
    n = len(a)
    rms = (sum(x * x for x in a) / n) ** 0.5
    worst = max(abs(x - y) for x, y in zip(a, b))
    if rms == 0.0:
        return 0.0 if worst == 0.0 else float("inf")
    return worst / rms


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("left", type=Path)
    ap.add_argument("right", type=Path)
    ap.add_argument("--tol", type=float, required=True,
                    help="largest relative difference, against the record's own RMS, "
                         "that still counts as a regrouped sum rather than a "
                         "different computation. Declare it before running.")
    ap.add_argument("--exact", action="append", default=[],
                    help="a record name that must be BIT IDENTICAL, repeatable. "
                         "Use it for fields that are only scattered and gathered.")
    ap.add_argument("--quiet", action="store_true",
                    help="report only the records that are not identical")
    args = ap.parse_args()

    la, lb = records(args.left), records(args.right)
    if [n for n, _ in la] != [n for n, _ in lb]:
        print("FAIL: the two restarts do not hold the same records in the same order")
        return 1

    worst_name, worst_val = None, 0.0
    n_same = n_round = n_diff = 0
    exact_failed = []
    for (name, a), (_, b) in zip(la, lb):
        if a == b:
            n_same += 1
            if not args.quiet:
                print(f"[ identical ] {name}")
            continue
        if len(a) != len(b):
            print(f"[ DIFFERENT ] {name}: record lengths differ")
            n_diff += 1
            continue
        va, kind = as_numbers(a)
        vb, _ = as_numbers(b)
        if kind == "int4":
            print(f"[ DIFFERENT ] {name}: integer record differs")
            n_diff += 1
            continue
        r = relative(va, vb)
        if r > worst_val:
            worst_name, worst_val = name, r
        if name in args.exact:
            exact_failed.append(name)
        if r <= args.tol:
            n_round += 1
            print(f"[ rounding  ] {name}: {r:.3e}")
        else:
            n_diff += 1
            print(f"[ DIFFERENT ] {name}: {r:.3e}")

    print()
    print(f"{len(la)} records: {n_same} identical, {n_round} at rounding scale, "
          f"{n_diff} beyond it")
    if worst_name is not None:
        print(f"worst: {worst_name} at {worst_val:.3e} relative, tolerance {args.tol:.0e}")

    for name in args.exact:
        if name not in [n for n, _ in la]:
            print(f"FAIL: --exact named {name}, which is not a record in these files")
            return 1
    if exact_failed:
        print(f"FAIL: these must be bit identical and are not: {', '.join(exact_failed)}")
        print("      They carry no arithmetic between the scatter and the gather,")
        print("      so a difference means the permutation is not its own inverse.")
        return 1
    if n_diff:
        print("FAIL: differences beyond a regrouped sum.")
        return 1
    print("PASS: every difference is at the scale of a regrouped sum.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
