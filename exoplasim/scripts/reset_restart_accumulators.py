#!/usr/bin/env python3
"""Zero the ACCUMULATOR records in a PlaSim restart, leaving state untouched.

    python exoplasim/scripts/reset_restart_accumulators.py IN.rest OUT.rest

WHY. `run_exoplasim.py --restart-from` hands a donor run's restart to
`model.configure(restartfile=...)` and PlaSim reads it WHOLE -- including the
partial accumulation the donor happened to be carrying when it stopped, and the
counters that go with it. The seeded run then opens mid-window: it adds its own
steps to somebody else's sum and divides by its own count.

The damage is one output record and it is measurable on a field that cannot
vary. Under `NLOWIO = 1` the land mask reads 0.95341 in orbit 0 bin 0 and
exactly 1.0 in every other bin of every orbit; a cold-started low-I/O run reads
1.0 throughout, and a clean-I/O run accumulates nothing so shows nothing. Both
conditions are needed, which is what identifies the cause. TASKS.md CLIM-31,
`exoplasim/notes/first-output-bin.md`.

WHY ZEROING IS SAFE, which is the whole argument for this being a fix rather
than a hack. Every record named below is one the model itself sets to zero at an
interval boundary -- `outmod.f90:2435`, `seamod.f90:212`, and the equivalents in
icemod and oceanmod. So a restart with them zeroed describes a run sitting
exactly at the start of an accumulation window, which is a state PlaSim occupies
several times an orbit. Nothing here is prognostic: temperatures, ice, snow and
soil water are left byte-identical, and the file's length and record order do
not change.

THE FORMAT is sequential Fortran unformatted: each record is a 4-byte length, a
payload, and a 4-byte trailing length. Records alternate between a 16-character
NAME and the data it labels, which is what `restart_ini` walks to build its
table. Zeroing the payload bytes gives 0 for an integer and +0.0 for an IEEE
float alike, so the same operation serves both without knowing the type.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import struct

# Reset to zero by the model at every interval boundary. Counters first, then
# the arrays they normalise. `darea` is NOT here and must never be: it looks
# like an accumulator by its name and is grid area weights, which are state.
ACCUMULATOR_RECORDS = {
    # counters
    "naccu", "naccua", "naccuice", "naccuo", "naccuoce", "naccuout",
    # outmod
    "atsa", "atsama",
    # outmod, the accumulated hurricane indices. Zeroed at the same interval
    # boundary as everything above (outmod.f90:2490) and written to the restart
    # unconditionally (plasim.f90:958), so they are here whether or not
    # NSTORMDIAG is on. notes/audits/dormant-exoplasim-modules.md finding 4.
    "agpi", "aventi", "alaav", "ampoti", "avrmpi", "acapen", "alnb", "achim",
    # seamod
    "cheata", "clhdta", "clhfla", "clwfla", "cpmea", "cprsa", "croffa",
    "cshdta", "cshfla", "cswfla", "ctauxa", "ctauya", "cust3a",
    # icemod
    "xcfluxa", "xcfluxna", "xcfluxra", "xcpmea", "xcroffa", "xfluxca",
    "xflxicea", "xheata", "ximelta", "xofluxa", "xqmelta", "xscflxa",
    "xsmelta", "xstoia", "xtsfluxa",
    # oceanmod
    "ydssta", "yfldoa", "yfssta", "yheata", "yifluxa", "yqhda",
}


def read_records(data: bytes):
    """Every sequential-unformatted record as (start, end, payload)."""
    out, pos, n = [], 0, len(data)
    while pos < n:
        if pos + 4 > n:
            raise ValueError(f"truncated record header at byte {pos}")
        (length,) = struct.unpack("<i", data[pos:pos + 4])
        body = pos + 4
        tail = body + length
        if length < 0 or tail + 4 > n:
            raise ValueError(f"record at byte {pos} claims {length} bytes, "
                             "which does not fit; this is not a little-endian "
                             "sequential unformatted file")
        (check,) = struct.unpack("<i", data[tail:tail + 4])
        if check != length:
            raise ValueError(f"record at byte {pos} has header {length} and "
                             f"trailer {check}; they must match")
        out.append((pos, tail + 4, data[body:tail]))
        pos = tail + 4
    return out


def reset(src: Path, dst: Path) -> dict:
    raw = src.read_bytes()
    records = read_records(raw)
    out = bytearray(raw)
    zeroed, seen, name = [], [], None
    for start, end, payload in records:
        if len(payload) == 16:
            candidate = payload.decode("ascii", "replace").strip()
            # A 16-byte record is a name only if it looks like one. A data
            # record of four 4-byte values is also 16 bytes, so the check
            # matters: names are printable and unpadded on the right.
            if candidate and all(32 <= b < 127 for b in payload):
                name = candidate
                seen.append(name)
                continue
        if name is not None and name in ACCUMULATOR_RECORDS:
            body = start + 4
            out[body:body + len(payload)] = b"\x00" * len(payload)
            zeroed.append((name, len(payload)))
        name = None
    dst.write_bytes(bytes(out))
    return {"records": len(records), "names": seen, "zeroed": zeroed}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("destination", type=Path)
    ap.add_argument("--list", action="store_true",
                    help="print every record name found and exit without writing")
    args = ap.parse_args()
    if args.list:
        for i, n in enumerate(
                x for x in reset(args.source, Path("/dev/null"))["names"]):
            print(f"{i:4d}  {n}")
        return
    info = reset(args.source, args.destination)
    if not info["zeroed"]:
        raise SystemExit(
            f"{args.source} carries none of the {len(ACCUMULATOR_RECORDS)} "
            "accumulator records this knows about. Either the restart predates "
            "them or the names have moved; refusing to claim a reset that did "
            "not happen.")
    print(f"{info['records']} records, {len(info['names'])} named; "
          f"zeroed {len(info['zeroed'])}:")
    for n, size in info["zeroed"]:
        print(f"  {n:12} {size:9d} bytes")
    if args.source.stat().st_size != args.destination.stat().st_size:
        raise SystemExit("the copy changed size, which a payload zeroing "
                         "cannot do; something is wrong")


if __name__ == "__main__":
    main()
