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
conditions are needed, which is what identifies the cause. CLIM-31,
`exoplasim/notes/first-output-bin.md`.

WHY ZEROING IS SAFE, which is the whole argument for this being a fix rather
than a hack. Every record it touches is one the model itself puts back to a
clean value at an interval boundary -- `outmod.f90:outreset`, `seamod.f90:212`,
and the equivalents in icemod and oceanmod. So a restart it has written
describes a run sitting exactly at the start of an accumulation window, which
is a state PlaSim occupies several times an orbit. Nothing prognostic is
touched: temperatures, ice, snow and soil water are left byte-identical, and
the file's length and record order do not change.

WHICH RECORDS, AND WHY NOT A LIST HERE. `restart_schema.py` says which records
are accumulators and what the model resets each one TO, and that column is
checked against the model source rather than asserted. A list maintained here
covered 49 of the 114 accumulators a restart carries and had drifted without
anything going red -- including `aprl`, `aprc`, `assol`, `atsol` and the six
spectral accumulators, which is to say most of what CLIM-31 was about.

CLEAN IS NOT ALWAYS ZERO, which is why the value comes from the schema and not
from this file's idea of an accumulator. `tempmin` resets to 1.0e3 and `atsami`
to 1.0e10, both running MINIMA whose zeroed form would report 0 K for the rest
of the run. `asndch` and `aanrho` are reset by the model nowhere at all -- the
first on purpose, the second not -- so they are left exactly as they were and
named in the report.

THE FORMAT is sequential Fortran unformatted; `restart_format.py` owns it.
"""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import restart_format  # noqa: E402
import restart_schema  # noqa: E402
from _paths import MODEL_SRC  # noqa: E402

MODEL_SOURCE = MODEL_SRC / "plasim" / "src"


def clean_payload(record, value: float, real_bytes: int) -> bytes:
    """The bytes of a record holding `value` in every element."""
    count = record.nbytes // real_bytes
    if count * real_bytes != record.nbytes:
        raise restart_schema.ConversionError(
            f"'{record.name}' is {record.nbytes} bytes, which is not a whole "
            f"number of {real_bytes}-byte reals")
    code = "<f4" if real_bytes == 4 else "<f8"
    return struct.pack(f"<{count}{'f' if real_bytes == 4 else 'd'}",
                       *([value] * count)) if count else b""


def reset(src: Path, dst: Path) -> dict:
    """Every accumulator put back to the value the model resets it to.

    Returns what it did, per record, for the caller's manifest.
    """
    records = restart_format.read(src)
    geometry, real_bytes = restart_schema.describe(records)
    inventory = restart_schema.inventory_from_source(MODEL_SOURCE)

    unknown = sorted({r.name for r in records} - set(restart_schema.POLICY))
    if unknown:
        raise restart_schema.ConversionError(
            f"{src} holds records the schema does not name: "
            f"{', '.join(unknown)}. Add them to restart_schema.POLICY rather "
            "than letting a record nobody has classified through a tool whose "
            "whole job is knowing which records are accumulators.")

    out, zeroed, sentinels, left = [], [], [], []
    for rec in records:
        pol = restart_schema.POLICY[rec.name]
        if pol.semantic != restart_schema.ACCUMULATOR:
            out.append(rec)
            continue
        if pol.model_reset == "zero":
            out.append(restart_format.Record(name=rec.name,
                                             payload=b"\x00" * rec.nbytes,
                                             offset=rec.offset))
            zeroed.append((rec.name, rec.nbytes))
        elif pol.model_reset == "sentinel":
            out.append(restart_format.Record(
                name=rec.name,
                payload=clean_payload(rec, pol.reset_value, real_bytes),
                offset=rec.offset))
            sentinels.append((rec.name, pol.reset_value))
        else:
            out.append(rec)
            left.append(rec.name)
    restart_format.write(dst, out, overwrite=True)
    return {"records": len(records), "geometry": geometry.label,
            "real_bytes": real_bytes,
            "names": [r.name for r in records], "zeroed": zeroed,
            "sentinels": sentinels, "never_reset_by_the_model": left}


def self_test() -> int:
    """Prove the reset against the model's own reset routines.

    A right answer rather than a comparison: `outreset` and its per-module
    equivalents say what every accumulator's clean value is, so a restart this
    tool has written must hold exactly those values and nothing else may have
    moved. The controls are the four records whose clean value is not zero --
    if the tool ever goes back to zeroing by name, two running minima land on
    0 K and this fails.
    """
    import tempfile

    runs = sorted(Path(__file__).resolve().parents[1].glob("runs/run_*/plasim_restart"))
    if not runs:
        print("self-test: no run directory holds a plasim_restart")
        return 1
    donor = runs[0]
    inventory = restart_schema.inventory_from_source(MODEL_SOURCE)
    resets = restart_schema.model_resets_from_source(MODEL_SOURCE)

    failures = []
    with tempfile.TemporaryDirectory(prefix="reset_selftest_") as d:
        out = Path(d) / "reset"
        info = reset(donor, out)
        before = {r.name: r for r in restart_format.read(donor)}
        after = {r.name: r for r in restart_format.read(out)}

        if [*before] != [*after]:
            failures.append("the record set or its order changed")
        if donor.stat().st_size != out.stat().st_size:
            failures.append("the file changed size, which a payload rewrite "
                            "cannot do")

        touched = 0
        for name, rec in after.items():
            pol = restart_schema.POLICY[name]
            if pol.semantic != restart_schema.ACCUMULATOR:
                if rec.payload != before[name].payload:
                    failures.append(f"{name} is not an accumulator and moved")
                continue
            kind, value = restart_schema.derived_model_reset(
                name, inventory, resets)
            if kind == "zero":
                if rec.payload != b"\x00" * rec.nbytes:
                    failures.append(f"{name} is reset to zero by the model "
                                    "and was not zeroed")
                touched += 1
            elif kind == "sentinel":
                want = clean_payload(rec, value, info["real_bytes"])
                if rec.payload != want:
                    failures.append(
                        f"{name} is reset to {value} by the model and holds "
                        "something else; a zeroed running minimum reports 0 K "
                        "for the rest of the run")
                touched += 1
            else:
                if rec.payload != before[name].payload:
                    failures.append(f"{name} is reset by the model nowhere "
                                    "and was changed anyway")

    print(f"self-test on {donor.parent.name}/{donor.name}, "
          f"{info['geometry']} at {info['real_bytes']}-byte reals")
    print(f"  {len(info['zeroed'])} accumulators zeroed, "
          f"{len(info['sentinels'])} put back to a nonzero clean value "
          f"({', '.join(f'{n}={v:g}' for n, v in info['sentinels'])}), "
          f"{len(info['never_reset_by_the_model'])} left alone "
          f"({', '.join(info['never_reset_by_the_model'])})")
    print(f"  {touched} records checked against the model's own reset routines")
    for f in failures:
        print(f"  FAIL {f}")
    if failures:
        return 1
    print("  every accumulator holds the value the model resets it to, and "
          "nothing else moved")
    return 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path, nargs="?")
    ap.add_argument("destination", type=Path, nargs="?")
    ap.add_argument("--list", action="store_true",
                    help="print every record name found and exit without writing")
    ap.add_argument("--self-test", action="store_true",
                    help="check a written restart against the model's own "
                         "reset routines, including the four records whose "
                         "clean value is not zero")
    args = ap.parse_args()
    if args.self_test:
        raise SystemExit(self_test())
    if args.list:
        for i, n in enumerate(
                x for x in reset(args.source, Path("/dev/null"))["names"]):
            print(f"{i:4d}  {n}")
        return
    info = reset(args.source, args.destination)
    if not info["zeroed"]:
        raise SystemExit(
            f"{args.source} carries no accumulator record the schema names. "
            "Either the restart predates them or the names have moved; "
            "refusing to claim a reset that did not happen.")
    print(f"{info['records']} records at {info['geometry']}; "
          f"zeroed {len(info['zeroed'])}, "
          f"{len(info['sentinels'])} put back to a nonzero clean value, "
          f"{len(info['never_reset_by_the_model'])} left as they were")
    for n, v in info["sentinels"]:
        print(f"  {n:12} <- {v:g}")
    for n in info["never_reset_by_the_model"]:
        print(f"  {n:12} left: the model resets it nowhere")
    if args.source.stat().st_size != args.destination.stat().st_size:
        raise SystemExit("the copy changed size, which a payload zeroing "
                         "cannot do; something is wrong")


if __name__ == "__main__":
    main()
