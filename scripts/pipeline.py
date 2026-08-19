#!/usr/bin/env python3
"""Plan the pipeline: what exists, what must run, and what blocks a target.

    python scripts/pipeline.py --status            # what is present, what is not
    python scripts/pipeline.py --plan carve_list   # ordered steps to reach a target
    python scripts/pipeline.py --register          # the artifact register

The graph is `config/pipeline.yaml` and that file is its only source. This script
reads it; `WORKFLOW.md` explains it and references its step ids. Nothing restates
it, because two representations of one graph is the duplication both documents
exist to prevent.

## What this does NOT do

It does not run anything. Planning and executing are separated deliberately:
five steps in the graph cost hours, and a script that could start one by accident
is worse than no script. `--plan` prints; you run.

It also does not check whether artifacts AGREE with each other. That is
`check_consistency.py`, which compares terrain hashes, provenance stamps and
conventions, and it is the better tool for it. The division is: consistency
answers "do these describe the same world", this answers "what exists and what
must run to reach X". Asking either to do the other's job would give two answers
to one question.

## The carve gate, and the margin in it

`orogen` is gated on nothing outstanding touching a step upstream of
`carve_list`. That is computable from the graph plus the step each open task
names, and `--status` computes it. What it CANNOT do is decide about a task that
names no step, or one whose relevance is a matter of degree. Those are reported
as a residual to be arbitrated rather than silently counted or silently ignored,
because a gate that hides its own uncertainty is a gate nobody should trust.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT / "config" / "pipeline.yaml"
TASKS = ROOT / "TASKS.md"


def load() -> dict:
    return yaml.safe_load(GRAPH.read_text(encoding="utf-8"))


def active_build() -> str:
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    return cfg.get("source_build") or cfg.get("model", {}).get("source_build", "")


def steps_by_id(graph: dict) -> dict:
    return {s["id"]: s for s in graph["steps"]}


def resolve(path: str, build: str) -> Path:
    return ROOT / path.replace("{build}", build).replace("{name}", "k25v")


def present(step: dict, build: str) -> tuple[bool | None, list[str]]:
    """Does every artifact this step declares exist? Presence only, by design.

    Returns None where presence is NOT DECIDABLE from the filesystem. A climate
    run is named by a UUID and has no fixed path, so the honest answer is "ask
    INDEX.json", not a guess. The first version of this used INDEX.json as
    evidence for each run step, and since four steps write that one file, three
    runs reported as done because a fourth had happened.
    """
    if step.get("presence") == "by_index":
        return None, []
    missing = [w for w in step.get("writes", []) if not resolve(w, build).exists()]
    return (not missing), missing


def upstream(target: str, by_id: dict, seen: set | None = None) -> set:
    """Every step that must precede `target`. Cycle-safe: the graph is cyclic."""
    seen = set() if seen is None else seen
    for need in by_id[target].get("needs", []):
        if need not in seen:
            seen.add(need)
            upstream(need, by_id, seen)
    return seen


def order(target: str, by_id: dict) -> list[str]:
    """Depth-first postorder, which is a valid run order and tolerates cycles."""
    out, mark = [], set()

    def visit(node: str, stack: frozenset) -> None:
        if node in mark or node in stack:
            return                      # already placed, or a cycle edge
        for need in by_id[node].get("needs", []):
            visit(need, stack | {node})
        mark.add(node)
        out.append(node)

    visit(target, frozenset())
    return out


def task_steps() -> tuple[dict, list]:
    """Open tasks that name a step, and the residual that names none."""
    text = TASKS.read_text(encoding="utf-8")
    body = text[text.index("## Open"):text.index("## Closed")]
    attached, residual = {}, []
    for line in body.splitlines():
        if not line.startswith("| ") or line.startswith(("| id", "| ---")):
            continue
        tid = line.strip("|").split("|")[0].strip()
        m = re.search(r"\[step:\s*([a-z_0-9, ]+)\]", line)
        if m:
            for s in (x.strip() for x in m.group(1).split(",")):
                attached.setdefault(s, []).append(tid)
        else:
            residual.append(tid)
    return attached, residual


def cmd_register(graph: dict) -> int:
    by_id = steps_by_id(graph)
    consumers: dict[str, list[str]] = {}
    for s in graph["steps"]:
        for need in s.get("needs", []):
            consumers.setdefault(need, []).append(s["id"])
    print("| artifact | step | read by |")
    print("| --- | --- | --- |")
    for s in graph["steps"]:
        reads = ", ".join(consumers.get(s["id"], [])) or "terminal"
        for w in s.get("writes", []) or ["(no artifact)"]:
            print(f"| `{w}` | `{s['id']}` ({s['script']}) | {reads} |")
    return 0


def cmd_status(graph: dict) -> int:
    build = active_build()
    by_id = steps_by_id(graph)
    print(f"active build: {build}\n")
    miss_total = 0
    for s in graph["steps"]:
        ok, missing = present(s, build)
        if not s.get("writes"):
            continue
        if ok is None:
            print(f"[ ask idx ] {s['id']}  (UUID-named; INDEX.json is the record)")
        elif ok:
            print(f"[ present ] {s['id']}")
        else:
            miss_total += 1
            print(f"[ MISSING ] {s['id']}  ({len(missing)} of "
                  f"{len(s['writes'])}: {missing[0]})")
    attached, residual = task_steps()
    unknown = {s: t for s, t in attached.items() if s not in by_id}
    gate = upstream("carve_list", by_id) | {"carve_list"}
    blocking = {s: t for s, t in attached.items() if s in gate}
    print(f"\n{len(graph['steps'])} steps, {miss_total} with a missing artifact")
    if unknown:
        # A task naming a step the graph does not have is INVISIBLE to the gate:
        # it is neither blocking nor residual, so it is silently uncounted. That
        # is how a step rename disarms the gate, so say it loudly rather than
        # dropping the row.
        print("\nNAMES A STEP THIS GRAPH DOES NOT HAVE, so the gate cannot see it:")
        for s, ids in sorted(unknown.items()):
            print(f"  {s:24s} {', '.join(ids)}")
        print("  Repoint the [step: ...] marker in TASKS.md, or add the row here.")
    print("\nTHE CARVE GATE: open tasks touching a step upstream of carve_list")
    if blocking:
        for s, t in sorted(blocking.items()):
            print(f"  {s:24s} {', '.join(t)}")
    else:
        print("  none attached")
    print(f"\nRESIDUAL, naming no step and therefore not counted either way: "
          f"{len(residual)}")
    if residual:
        print("  " + " ".join(residual))
        print("  These are for a person or the assistant to arbitrate. The gate is\n"
              "  narrowed by the graph, not decided by it.")
    return 0


def cmd_plan(graph: dict, target: str, force: bool) -> int:
    by_id = steps_by_id(graph)
    if target not in by_id:
        raise SystemExit(f"no step '{target}'. Known: {', '.join(sorted(by_id))}")
    build = active_build()
    seq = order(target, by_id)
    print(f"to make `{target}` current, in order:\n")
    hours, undecidable = 0, 0
    for sid in seq:
        s = by_id[sid]
        ok, _ = present(s, build)
        if ok is True and not force and sid != target:
            print(f"  [ skip    ] {sid:26s} artifact present")
            continue
        # Everything below here is a step the plan may have to run, and an
        # undecidable one is NOT an exception to that: it is the case where the
        # filesystem cannot say, which is a reason to cost it and print its gate
        # rather than a reason to pass over it. Every hours-costing step in the
        # graph is `by_index`, so treating undecidable as skippable made the
        # hour count structurally zero and silenced every gate attached to a run.
        mark = "HOURS" if s["cost"] == "hours" else s["cost"]
        hours += s["cost"] == "hours"
        if ok is None:
            undecidable += 1
            print(f"  [ check   ] {sid:26s} {mark:8s} UUID-named: has INDEX.json "
                  f"one on this build?")
        else:
            print(f"  [ run     ] {sid:26s} {mark:8s} {s['script']}")
        if s.get("gate"):
            print(f"              GATE: {' '.join(s['gate'].split())}")
    print(f"\n{len(seq)} steps in the closure, {hours} of them cost hours.")
    if undecidable:
        print(f"{undecidable} of the steps listed are UUID-named, so whether they "
              f"still need to run is\nINDEX.json's answer and not the "
              f"filesystem's. The hour count assumes they do.")
    print("This prints a plan. It does not run anything.")
    return 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--plan", metavar="STEP")
    ap.add_argument("--register", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="with --plan, list steps whose artifacts already exist")
    args = ap.parse_args()
    graph = load()
    if args.register:
        raise SystemExit(cmd_register(graph))
    if args.plan:
        raise SystemExit(cmd_plan(graph, args.plan, args.force))
    raise SystemExit(cmd_status(graph))


if __name__ == "__main__":
    main()
