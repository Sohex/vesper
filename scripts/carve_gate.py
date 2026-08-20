#!/usr/bin/env python3
"""Which open tasks touch a step upstream of the carve.

    python scripts/carve_gate.py

`orogen` is gated on nothing outstanding still moving an artifact the carve
verdict is computed from. That is not a ceremony about irreversibility -- a
wrong verdict is recoverable, because a build is replaced wholesale rather than
edited -- it is the ordinary condition that the inputs have stopped moving.

## Why this is its own script

Answering it needs two things: the graph, to know what is upstream of
`carve_list`, and `TASKS.md`, to know what is open. Those belong to different
tools, and the direction of the dependency is the whole point.

`pipeline.py` is the GRAPH. It answers what exists and what must run to reach a
target, and its own docstring states the rule it lives by: consistency answers
"do these describe the same world", the graph answers "what exists and what must
run", and asking either to do the other's job would give two answers to one
question. It used to read `TASKS.md` as well, which made a prose tracker an
input to the pipeline planner: editing a status cell changed what the planner
reported, and a regex over human-written markdown sat in the middle of it. The
gate is a TRACKER question that needs graph help, so it reads the graph -- never
the other way round.

The practical cost of the old arrangement was not tidiness. Because the gate was
printed by `--status`, every row filed against a step read as a blocker on that
step, including rows that were only saying "this work happens here". The gate
could not distinguish "a finding still moves this artifact" from "this step has
not run yet", and the pipeline not having run is not a blocker on running it.

## The margin, which is deliberate

What this CANNOT do is decide about a task that names no step, or one whose
relevance is a matter of degree. Those are reported as a RESIDUAL to be
arbitrated rather than silently counted or silently ignored, because a gate that
hides its own uncertainty is a gate nobody should trust. The graph narrows the
judgement; it does not replace it.

A task naming a step the graph does not have is worse than either, because it is
INVISIBLE: neither blocking nor residual. That is how a step rename would disarm
the gate, so it is reported loudly and separately.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import pipeline  # noqa: E402  the graph, and only in this direction

TASKS = ROOT / "TASKS.md"


def open_tasks() -> tuple[dict, list]:
    """Open tasks that name a step, and the residual that names none.

    OPENNESS IS READ FROM THE STATUS COLUMN, not from where a row sits. This
    used to slice between the `## Open` and `## Closed` headings, which stopped
    working when TASKS.md went to one table per prefix with closed ids left in
    place as stubs -- every stub would have counted as a residual, and the gate
    would have reported 142 tasks to arbitrate. Reading the status is also
    simply more robust: it does not care about headings, ordering, or which
    table a row is in.
    """
    text = TASKS.read_text(encoding="utf-8")
    attached: dict[str, list[str]] = {}
    residual: list[str] = []
    for line in text.splitlines():
        if not line.startswith("| ") or line.startswith(("| id", "| ---")):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 4 or cells[3].lower().startswith(("done", "wontfix")):
            continue
        tid = cells[0]
        m = re.search(r"\[step:\s*([a-z_0-9, ]+)\]", line)
        if m:
            for s in (x.strip() for x in m.group(1).split(",")):
                attached.setdefault(s, []).append(tid)
        else:
            residual.append(tid)
    return attached, residual


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()

    graph = pipeline.load()
    by_id = pipeline.steps_by_id(graph)
    attached, residual = open_tasks()

    gate = pipeline.upstream("carve_list", by_id) | {"carve_list"}
    blocking = {s: t for s, t in attached.items() if s in gate}
    unknown = {s: t for s, t in attached.items() if s not in by_id}

    if unknown:
        print("NAMES A STEP THIS GRAPH DOES NOT HAVE, so the gate cannot see it:")
        for s, ids in sorted(unknown.items()):
            print(f"  {s:24s} {', '.join(ids)}")
        print("  Repoint the [step: ...] marker in TASKS.md, or add the row to "
              "config/pipeline.yaml.\n")

    print("THE CARVE GATE: open tasks touching a step upstream of carve_list")
    if blocking:
        for s, t in sorted(blocking.items()):
            print(f"  {s:24s} {', '.join(t)}")
    else:
        print("  none attached")

    print(f"\nRESIDUAL, naming no step and therefore not counted either way: "
          f"{len(residual)}")
    if residual:
        print("  " + " ".join(residual))
        print("  These are for a person to arbitrate. The gate is narrowed by "
              "the graph, not decided by it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
