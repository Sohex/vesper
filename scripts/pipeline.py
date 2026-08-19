#!/usr/bin/env python3
"""Plan the pipeline: what exists, what must run, and what blocks a target.

    python scripts/pipeline.py --status            # what is present, what is not
    python scripts/pipeline.py --plan carve_list   # ordered steps to reach a target
    python scripts/pipeline.py --register          # the artifact register
    python scripts/pipeline.py --purge orogen      # what a change to X makes worthless
    python scripts/pipeline.py --purge orogen --execute      # ... and delete it

The graph is `config/pipeline.yaml` and that file is its only source. This script
reads it; `WORKFLOW.md` explains it and references its step ids. Nothing restates
it, because two representations of one graph is the duplication both documents
exist to prevent.

## What this does NOT do

It does not run anything. Planning and executing are separated deliberately:
five steps in the graph cost hours, and a script that could start one by accident
is worse than no script. `--plan` prints; you run.

`--purge` is the one thing here that touches the filesystem, and it only ever
deletes. It is the same separation seen from the other side: running a step is
expensive and is yours to start, but working out WHAT a change makes worthless
is a graph question, and doing it by hand is how a keep-list goes wrong. It is a
dry run until `--execute`, exactly like `scripts/archive_runs.py`.

## Why purge is a graph question

CLAUDE.md rule 7: an upstream change makes everything below it WORTHLESS rather
than stale, so the question after any change is "what is now worthless". That is
answerable only from a complete graph, which is what `config/pipeline.yaml` is
for -- and it is the reason `writes` has to be the complete artifact set rather
than a marker file. A row naming only its report purges its report and leaves
the .sra beside it, then says it succeeded.

Two edges are not `needs` and purge uses both:

`reads_export` is the step consuming `source/{build}/`. Four steps do it and all
four declared `needs: []`, which read as "depends on nothing" when it meant
"depends on the terrain and nothing else in this pass". `--purge orogen` seeds
from those four rather than from `orogen`'s own `needs`.

And the traversal STOPS AT `orogen` rather than passing through it. Loop A is a
cycle -- carve_list feeds orogen, orogen's export feeds hydrography -- so
unrestricted reachability from any climate step reaches the generation step and
then everything, including `source/{build}/` itself. That answer is wrong, not
merely alarming: a new climatology does not invalidate the terrain it was
computed on. The cycle closes across an ITERATION boundary and "what is worthless
now" is a within-pass question, so the generation step is where it is cut.

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
import shutil
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "lib"))
from provenance import INERT_CONFIG_KEYS, artifact_drift  # noqa: E402
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


def expand(path: str, build: str) -> list[Path]:
    """Every file a `writes` entry currently stands for.

    ONE resolver for presence and for purge, deliberately. Two would drift, and
    the direction they would drift in is the dangerous one: a purge that expands
    a pattern the presence check does not would delete something `--status` had
    just called absent.

    Three forms, matching what `config/pipeline.yaml` documents. A glob returns
    its matches, so an unmatched glob is legitimately empty. A directory (written
    with a trailing slash) returns the files under it, because a step that owns a
    directory owns what accumulates in it -- per-run convergence products are
    named by a UUID and cannot be listed any other way. Anything else is a single
    path, returned whether or not it exists; the caller decides what absence
    means.
    """
    if "*" in path:
        pattern = path.replace("{build}", build).replace("{name}", "k25v")
        return sorted(ROOT.glob(pattern))
    target = resolve(path, build)
    if path.endswith("/"):
        return sorted(q for q in target.rglob("*") if q.is_file())
    return [target]


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
    # Through expand(), so a glob is present when it matches anything and a
    # directory when it holds anything. Testing resolve().exists() called an
    # empty directory present, which is what a purged step looks like.
    missing = [w for w in step.get("writes", [])
               if not any(q.exists() for q in expand(w, build))]
    return (not missing), missing


def stale(step: dict, build: str, config: dict) -> list[str]:
    """Config drift under a step's own outputs. Empty when nothing moved.

    PRESENCE IS NOT CURRENCY, and `present` above says so in its own name. A
    staged surface field is regenerated from the seed, the config and the code,
    so a config change makes it worthless rather than merely old -- but it stays
    on disk, and a planner that reads the filesystem calls that done. That is
    how the star-reweighted canopy of SPEC-5 sat unstaged through two commits
    while every document described it: `CLAUDE.md` rule 7 promises the ordering
    regenerates derived artifacts as a side effect, and the promise held only
    for steps whose output is judged by something other than existing.

    Checked only for steps with a TRACED inert set in
    `lib/provenance.py:INERT_CONFIG_KEYS`. A step absent from there is not
    checked, because guessing an inert set would flag every step on any
    parameter change, which teaches re-running generators to silence a message.
    """
    inert = INERT_CONFIG_KEYS.get(step["id"])
    if inert is None:
        return []
    out = []
    for w in step.get("writes", []):
        drift = artifact_drift(resolve(w, build), config, inert)
        if drift:
            out += [f"{Path(w).name}: {d}" for d in drift]
    return out


def upstream(target: str, by_id: dict, seen: set | None = None) -> set:
    """Every step that must precede `target`. Cycle-safe: the graph is cyclic."""
    seen = set() if seen is None else seen
    for need in by_id[target].get("needs", []):
        if need not in seen:
            seen.add(need)
            upstream(need, by_id, seen)
    return seen


def downstream(seed: str, by_id: dict) -> set:
    """Every step whose output a change to `seed` makes worthless.

    Reverse reachability over `needs`, plus the export edge: `orogen` reaches
    every step marked `reads_export`, because those consume `source/{build}/`
    without declaring a `needs` on the step that writes it, and must not declare
    one (`config/pipeline.yaml` says why against the key).

    CUT AT `orogen`, which is what makes the answer finite and correct rather
    than merely finite. Loop A is a cycle, so from any climate step the walk
    reaches carve_list, then orogen, then -- through the export -- everything,
    and it would offer to delete the terrain because a climatology moved. The
    generation step is the iteration boundary: things computed ON a build are
    downstream of it, and the build is not downstream of them.

    The seed is never in its own result, so `--purge X` deletes what X invalidates
    and never X itself. For `orogen` that is the difference between clearing the
    tree and deleting the export.
    """
    out: set[str] = set()
    frontier = ({s for s, v in by_id.items() if v.get("reads_export")}
                if seed == "orogen" else {seed})
    if seed == "orogen":
        out |= frontier
    while frontier:
        nxt = set()
        for sid, step in by_id.items():
            if sid in out or sid == seed or sid == "orogen":
                continue          # the cut, and the seed-excludes-itself rule
            if frontier & set(step.get("needs", [])):
                nxt.add(sid)
        out |= nxt
        frontier = nxt
    return out


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
        # The export edge is not a `needs` edge and would otherwise be invisible
        # here, which showed `source/{build}/` as read by nobody -- the exact
        # misreading `reads_export` was added to correct.
        if s.get("reads_export"):
            consumers.setdefault("orogen", []).append(s["id"])
    print("| artifact | step | read by |")
    print("| --- | --- | --- |")
    for s in graph["steps"]:
        reads = ", ".join(consumers.get(s["id"], [])) or "terminal"
        for w in s.get("writes", []) or ["(no artifact)"]:
            print(f"| `{w}` | `{s['id']}` ({s['script']}) | {reads} |")
    return 0


def cmd_status(graph: dict) -> int:
    build = active_build()
    config = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    stale_total = 0
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
            drift = stale(s, build, config)
            if drift:
                stale_total += 1
                print(f"[  STALE  ] {s['id']}  (config moved under it: "
                      f"{drift[0]}{'; ...' if len(drift) > 1 else ''})")
            else:
                print(f"[ present ] {s['id']}")
        else:
            miss_total += 1
            print(f"[ MISSING ] {s['id']}  ({len(missing)} of "
                  f"{len(s['writes'])}: {missing[0]})")
    attached, residual = task_steps()
    unknown = {s: t for s, t in attached.items() if s not in by_id}
    gate = upstream("carve_list", by_id) | {"carve_list"}
    blocking = {s: t for s, t in attached.items() if s in gate}
    print(f"\n{len(graph['steps'])} steps, {miss_total} with a missing "
          f"artifact, {stale_total} present but built from an older config")
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
    config = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
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
        if ok and stale(s, build, config):
            print(f"  [ run     ] {sid:26s} {s.get('cost', ''):8s} "
                  f"{s.get('script', '')}   <- present but the config moved "
                  f"under it")
            continue
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


def cmd_purge(graph: dict, target: str, execute: bool) -> int:
    by_id = steps_by_id(graph)
    if target not in by_id:
        raise SystemExit(f"no step '{target}'. Known: {', '.join(sorted(by_id))}")
    build = active_build()
    doomed = downstream(target, by_id)

    print(f"active build: {build}")
    print(f"worthless if `{target}` changes: {len(doomed)} of "
          f"{len(graph['steps'])} steps\n")

    # Runs are named by a UUID and have no path in the graph, so purge cannot
    # find them and must not pretend the tree is clear without them. They also
    # have their own extract-then-delete path, which verifies the archive before
    # removing anything -- reimplementing that here would be a second copy of the
    # one thing in this repository that must not have two.
    runs = sorted(s for s in doomed if by_id[s].get("presence") == "by_index")

    total, plan = 0, []
    for step in graph["steps"]:            # graph order, so the tree reads top-down
        sid = step["id"]
        if sid not in doomed:
            continue
        files = []
        for w in step.get("writes", []):
            files += [q for q in expand(w, build) if q.is_file()]
        if not files:
            continue
        size = sum(q.stat().st_size for q in files)
        total += size
        plan.append((sid, files))
        print(f"  {sid:26s} {len(files):>4} files  {size/1e6:>9.1f} MB")

    print(f"\n{'deleting' if execute else 'would delete'} "
          f"{sum(len(f) for _, f in plan)} files, {total/1e9:.2f} GB")

    def note_runs() -> None:
        if runs:
            print(f"\nNOT COVERED, and not a gap this can close: {', '.join(runs)}")
            print("  Runs are UUID-named, so the graph has no path for them, and\n"
                  "  they are deleted by extracting their identity first:\n"
                  "      python scripts/archive_runs.py --include-live --execute")

    if not execute:
        note_runs()
        print("\ndry run; pass --execute to delete")
        return 0

    for sid, files in plan:
        for q in files:
            # missing_ok: expand() listed the tree a moment ago, and two steps
            # can legitimately name the same file. Racing a concurrent session
            # is not worth aborting a half-done purge over.
            q.unlink(missing_ok=True)
        # A directory a step OWNS goes whole, subdirectories included. Unlinking
        # its files and then rmdir()ing left the empty per-run subdirectories
        # behind, and an empty directory reads as a present artifact to anything
        # checking existence rather than contents. A directory a step merely
        # writes INTO is shared and is never named this way.
        for d in step_writes_dirs(by_id[sid], build):
            if d.is_dir():
                shutil.rmtree(d)
        print(f"  purged {sid}")

    print(f"\nfreed {total/1e9:.2f} GB")
    print("`world_state.json` and `exoplasim/runs/INDEX.json` are generated from\n"
          "what is on disk and now describe a tree that has moved. Regenerate:\n"
          "    python exoplasim/scripts/index_runs.py\n"
          "    python scripts/world_state.py")
    note_runs()
    return 0


def step_writes_dirs(step: dict, build: str) -> list[Path]:
    """The directories a step declares outright, which it therefore owns."""
    return [resolve(w, build) for w in step.get("writes", []) if w.endswith("/")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--plan", metavar="STEP")
    ap.add_argument("--register", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="with --plan, list steps whose artifacts already exist")
    ap.add_argument("--purge", metavar="STEP",
                    help="delete every artifact a change to STEP makes worthless. "
                         "Excludes STEP's own output, and never crosses `orogen`")
    ap.add_argument("--execute", action="store_true",
                    help="with --purge, actually delete; default is a dry run")
    args = ap.parse_args()
    graph = load()
    if args.register:
        raise SystemExit(cmd_register(graph))
    if args.purge:
        raise SystemExit(cmd_purge(graph, args.purge, args.execute))
    if args.plan:
        raise SystemExit(cmd_plan(graph, args.plan, args.force))
    raise SystemExit(cmd_status(graph))


if __name__ == "__main__":
    main()
