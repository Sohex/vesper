#!/usr/bin/env python3
"""Compare two ExoPlaSim restarts record by record.

    python exoplasim/scripts/diff_restarts.py A/plasim_status B/plasim_status
    python exoplasim/scripts/diff_restarts.py A B --top 20 --json out.json

Worldbuilding. Vesper is an invented planet and this script is about the
simulation of it: two checkpoint files from a toy climate model. Nothing here
is a claim about the world.

**A sha256 answers "are these the same" and nothing else.** When the answer is
no -- which it is for any A/B this project runs, since CLIM-44 established that
even a rank-count change moves the answer -- the next question is always the
same: WHICH records, and by how much. Round-off in 83 of 199 records is a
different finding from one record wrong by a factor, and the hash cannot tell
them apart.

This was written three times inline before it was written once here: for
CLIM-44's rank comparison, for the `zsolars` check inside it, and again for
CLIM-45's reduction identity. Each time the interesting output was the same
table.

## The format, which is why this is short

A restart is a flat sequence of Fortran sequential-access records. Each named
quantity is TWO records: a 16-byte name, then its data. So the file is walked
by reading a leading length word, taking that many bytes, checking the trailing
length word matches, and treating any 16-byte record as the name of whatever
follows. There is no index and no header.

`restartmod.f90` writes it and `reseek` reads it back by scanning for names, so
a record's size cannot dislodge any other record -- which is what makes it safe
to compare two files whose records differ in length, as they do across the
`zsolars` fix.

## What it reports

Records present in one file and not the other, then the differing ones ranked
by `max|a - b| / max|a|`, which is the measure that separates round-off from a
real change. Integer and seed records are compared but not ranked, since a
relative difference on a counter means nothing.

**Reading the number.** This project's own scale, measured: 1e-16 to 1e-14 is
arithmetic reordering, and CLIM-45's reduction identity sits there at one
timestep; 1e-10 after tens of steps is that same seed amplified by the model's
own chaos, which is where CLIM-44 found an 8-rank against 16-rank comparison;
anything at 1e-3 or above is a physics change and should be treated as one
until shown otherwise.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import PROJECT_ROOT  # noqa: E402
import restart_format  # noqa: E402

def read_records(path: Path) -> dict:
    """{name: raw bytes}, in file order. `restart_format` owns the framing."""
    return restart_format.payloads(path)


def compare(a: bytes, b: bytes) -> dict | None:
    """Relative and absolute difference, or None if the record is not float64."""
    if len(a) != len(b) or len(a) % 8 or len(a) < 8:
        return None
    x = np.frombuffer(a, dtype="<f8")
    y = np.frombuffer(b, dtype="<f8")
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        return None
    scale = max(float(np.abs(x).max()), float(np.abs(y).max()), 1e-300)
    return {"elements": int(x.size),
            "differing": int((x != y).sum()),
            "max_abs": float(np.abs(x - y).max()),
            "max_rel": float(np.abs(x - y).max() / scale),
            "field_max": float(np.abs(x).max())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("a", type=Path, help="a restart, or a run directory holding one")
    ap.add_argument("b", type=Path)
    ap.add_argument("--top", type=int, default=12, help="how many rows to print")
    ap.add_argument("--json", type=Path, help="write the full comparison here")
    args = ap.parse_args()

    paths = []
    for p in (args.a, args.b):
        p = p / "plasim_status" if p.is_dir() else p
        if not p.is_file():
            raise SystemExit(f"{p} is not a file")
        paths.append(p)

    ra, rb = (read_records(p) for p in paths)
    only_a = sorted(set(ra) - set(rb))
    only_b = sorted(set(rb) - set(ra))
    common = [k for k in ra if k in rb]
    identical = [k for k in common if ra[k] == rb[k]]

    rows, unranked = [], []
    for name in common:
        if ra[name] == rb[name]:
            continue
        stats = compare(ra[name], rb[name])
        if stats is None:
            unranked.append({"name": name,
                             "bytes_a": len(ra[name]), "bytes_b": len(rb[name])})
        else:
            rows.append({"name": name, **stats})
    rows.sort(key=lambda r: -r["max_rel"])

    print(f"{paths[0]}\n{paths[1]}\n")
    print(f"  records: {len(ra)} and {len(rb)}, {len(common)} in common")
    print(f"  byte-identical: {len(identical)}   differing: {len(common) - len(identical)}")
    if only_a:
        print(f"  only in the first : {', '.join(only_a)}")
    if only_b:
        print(f"  only in the second: {', '.join(only_b)}")
    if unranked:
        print(f"  not float64, compared but not ranked: "
              f"{', '.join(r['name'] for r in unranked)}")
    if rows:
        print(f"\n  {'record':14} {'max rel':>11} {'max abs':>12} {'elements differing':>20}")
        for r in rows[:args.top]:
            print(f"  {r['name']:14} {r['max_rel']:11.3e} {r['max_abs']:12.4e} "
                  f"{r['differing']:>10}/{r['elements']}")
        if len(rows) > args.top:
            print(f"  ... {len(rows) - args.top} more, --top to see them")
        worst = rows[0]["max_rel"]
        reading = ("arithmetic reordering" if worst < 1e-13 else
                   "round-off amplified by the model's own chaos" if worst < 1e-6 else
                   "A PHYSICS CHANGE until shown otherwise")
        print(f"\n  worst {worst:.3e} -- {reading}")
    else:
        print("\n  every common record is byte-identical")

    if args.json:
        args.json.write_text(json.dumps(
            {"a": str(paths[0]), "b": str(paths[1]),
             "records": [len(ra), len(rb)], "identical": len(identical),
             "only_in_a": only_a, "only_in_b": only_b,
             "unranked": unranked, "differing": rows}, indent=2) + "\n",
            encoding="utf-8")
        try:
            shown = args.json.relative_to(PROJECT_ROOT)
        except ValueError:
            shown = args.json
        print(f"wrote {shown}")


if __name__ == "__main__":
    main()
