"""Is a run taken in pieces the same experiment as the run taken whole?

Worldbuilding frame: a correctness check on the Vesper climate model's
ecological output stream. Nothing here is about the simulated planet.

    python exoplasim/scripts/compare_eco_streams.py WHOLE SEG1 SEG2 [SEG3 ...]

Each argument is one stream file. The model writes `plasim_eco` in place, and a
run driven through `run_exoplasim.py` or `continue_exoplasim.py` has it moved
aside per orbit as `MOST_ECO.NNNNN`; a bed that invokes the executable directly
leaves it under its own name.

WHAT IT IS FOR. The ecological stream (EFOR-2, `outmod.f90:ecogp`) reduces the
model's own timestep state to the intervals the biosphere integrates over. A
stream that is not restart-exact silently makes a continuation a different
experiment: the interval the segment boundary falls inside covers a shorter
span than it declares, and the seam lands in the sequence at exactly the place
a replay protocol later tests for one.

THE IDENTITY, and it has a right answer. Concatenate the blocks of the segments
and they must reproduce the blocks of the same run taken whole: the same number
of blocks, in the same order, over the same declared intervals, with the same
payloads. Everything this reports is a comparison against that identity rather
than a plausibility judgement, and it separates the three ways it can fail:

  STRUCTURE   a different number of blocks, or blocks over different intervals.
              This is the stream's own defect: the interval counter did not
              survive the restart, so a partial interval was thrown away or
              double counted. It is what `naccueco` is serialized for.
  BOUNDS      the same blocks over intervals that disagree about their own
              start, end or duration. Also the stream's own.
  PAYLOAD     the same blocks over the same intervals, carrying different
              numbers. This is INHERITED: it says the model's state did not
              reproduce across the restart, and the stream is reporting that
              faithfully. It cannot be fixed in the stream.

The third is not hypothetical on this model. `notes/audits/ecological-stream-restart-continuity.md`
carries the measurement and the run recipe.

RECORD LAYOUT. Fortran unformatted sequential, as `writegp` and `writescalar`
leave it: a 4-byte length, the payload, the same length again. The file opens
with two records that are not a block, the code-333 grid descriptor and the
sigma table. Every record after that is a header of eight int32 followed by its
payload, and a block is the run of records from a code-600 record up to the
next one.
"""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

import numpy as np

from _paths import COMPONENT_ROOT  # noqa: F401

# The stream's own codes, from outmod.f90:ecogp. The three scalars are the
# interval and are what makes a block self-describing.
INTERVAL_START, INTERVAL_END, INTERVAL_DURATION = 600, 601, 602
FIELD_NAMES = {
    610: "tas", 611: "ts", 612: "ps", 613: "hus", 614: "wind",
    615: "swdown", 616: "swup", 617: "swnet", 618: "lwnet", 619: "lwup",
    620: "czen", 621: "pr", 622: "prsn", 623: "prc", 624: "evap",
    625: "tasmax", 626: "tasmin", 627: "tsmax", 628: "tsmin",
}
# Position 7 of the header is `nstep - nstep1`, the step index WITHIN the run.
# It differs between a whole run and its segments BY CONSTRUCTION and carries no
# information the interval scalars do not carry better, so it is excluded rather
# than reported as a difference.
RUN_RELATIVE_STEP = 6


def raw_records(path: Path) -> list[bytes]:
    data = path.read_bytes()
    out, i = [], 0
    while i < len(data):
        if i + 4 > len(data):
            raise SystemExit(f"{path}: truncated record marker at byte {i}")
        n = struct.unpack("<i", data[i:i + 4])[0]
        out.append(data[i + 4:i + 4 + n])
        i += 4 + n + 4
    return out


def blocks(path: Path) -> list[dict]:
    """Every interval in one stream file, as header/payload pairs."""
    records = raw_records(path)[2:]
    if len(records) % 2:
        raise SystemExit(f"{path}: odd number of records after the file header, "
                         f"so it is not header/payload pairs")
    out: list[dict] = []
    for j in range(0, len(records), 2):
        head = struct.unpack("<8i", records[j])
        if head[0] == INTERVAL_START:
            out.append({"records": []})
        if not out:
            raise SystemExit(f"{path}: a record of code {head[0]} appears "
                             f"before any interval start (code "
                             f"{INTERVAL_START})")
        out[-1]["records"].append((head, records[j + 1]))
    for b in out:
        b["by_code"] = {h[0]: p for h, p in b["records"]}
    return out


def scalar(block: dict, code: int) -> float | None:
    payload = block["by_code"].get(code)
    if payload is None:
        return None
    return struct.unpack("<d" if len(payload) == 8 else "<f", payload)[0]


def interval(block: dict) -> tuple:
    return (scalar(block, INTERVAL_START), scalar(block, INTERVAL_END),
            scalar(block, INTERVAL_DURATION))


def compare(whole: list[dict], pieces: list[dict]) -> tuple[list[str], int]:
    """Report lines, and the number of checks that failed."""
    lines, failed = [], 0

    if len(whole) != len(pieces):
        lines.append(f"STRUCTURE  FAIL  the whole run wrote {len(whole)} "
                     f"intervals and the segments wrote {len(pieces)} between "
                     f"them. The interval counter did not survive the restart.")
        return lines, 1
    lines.append(f"STRUCTURE  ok    {len(whole)} intervals either way")

    bad = [k for k, (a, b) in enumerate(zip(whole, pieces))
           if interval(a) != interval(b)]
    if bad:
        failed += 1
        lines.append(f"BOUNDS     FAIL  {len(bad)} of {len(whole)} intervals "
                     f"declare different bounds")
        for k in bad[:5]:
            lines.append(f"                 interval {k}: whole "
                         f"{interval(whole[k])} segments {interval(pieces[k])}")
    else:
        lines.append(f"BOUNDS     ok    every interval declares the same "
                     f"start, end and duration in absolute seconds")

    codes = sorted(set(FIELD_NAMES) & set(whole[0]["by_code"]))
    first_diff = None
    worst: dict[int, float] = {}
    for k, (a, b) in enumerate(zip(whole, pieces)):
        for code in codes:
            pa, pb = a["by_code"][code], b["by_code"][code]
            if pa == pb:
                continue
            va = np.frombuffer(pa, dtype="<f4")
            vb = np.frombuffer(pb, dtype="<f4")
            d = float(np.max(np.abs(va - vb)))
            worst[code] = max(worst.get(code, 0.0), d)
            if first_diff is None:
                first_diff = k
    if first_diff is None:
        lines.append("PAYLOAD    ok    every interval carries identical values")
    else:
        failed += 1
        lines.append(f"PAYLOAD    FAIL  intervals 0 to {first_diff - 1} are "
                     f"identical and interval {first_diff} is the first to "
                     f"differ")
        lines.append("                 INHERITED, not the stream's: identical "
                     "structure and bounds over differing values means the "
                     "model's state did not reproduce across the restart.")
        for code in sorted(worst, key=lambda c: -worst[c])[:8]:
            lines.append(f"                 {FIELD_NAMES[code]:>7}  worst "
                         f"absolute difference {worst[code]:.6e}")

    # Headers, excluding the run-relative step index.
    head_bad = 0
    for a, b in zip(whole, pieces):
        for (ha, _), (hb, _) in zip(a["records"], b["records"]):
            if [x for k, x in enumerate(ha) if k != RUN_RELATIVE_STEP] != \
               [x for k, x in enumerate(hb) if k != RUN_RELATIVE_STEP]:
                head_bad += 1
    if head_bad:
        failed += 1
        lines.append(f"HEADERS    FAIL  {head_bad} records disagree on code, "
                     f"level, date, time or grid")
    else:
        lines.append("HEADERS    ok    code, level, date, time and grid agree "
                     "on every record")
    return lines, failed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("whole", type=Path,
                        help="the stream from the run taken in one piece")
    parser.add_argument("segments", type=Path, nargs="+",
                        help="the streams from the same run taken in pieces, "
                             "in order")
    args = parser.parse_args()

    whole = blocks(args.whole)
    pieces: list[dict] = []
    for path in args.segments:
        pieces.extend(blocks(path))

    print(f"whole    {args.whole}")
    for p in args.segments:
        print(f"segment  {p}")
    print()
    lines, failed = compare(whole, pieces)
    for line in lines:
        print(line)
    print()
    print("the segments reproduce the whole run" if not failed
          else f"{failed} of the four checks failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
