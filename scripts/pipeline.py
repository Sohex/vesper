#!/usr/bin/env python3
"""Plan the pipeline: what exists, what must run, and what blocks a target.

    python scripts/pipeline.py --status            # what is present, what is not
    python scripts/pipeline.py --plan carve_list   # ordered steps to reach a target
    python scripts/pipeline.py --register          # the artifact register
    python scripts/pipeline.py --purge orogen      # what a change to X makes worthless
    python scripts/pipeline.py --purge orogen --execute      # ... and delete it
    python scripts/pipeline.py --plan orogen --no-maps       # ... without the map frames

The graph is `config/pipeline.yaml` and that file is its only source. This script
reads it; the pipeline chapters under `docs/src/pipeline/` explain it and reference its step ids. Nothing restates
it, because two representations of one graph is the duplication both documents
exist to prevent.

## What this does NOT do

It does not run anything. Planning and executing are separated deliberately:
several steps in the graph cost hours, and a script that could start one by accident
is worse than no script. `--plan` prints; you run.

`--purge` is the one thing here that touches the filesystem, and it only ever
deletes. It is the same separation seen from the other side: running a step is
expensive and is yours to start, but working out WHAT a change makes worthless
is a graph question, and doing it by hand is how a keep-list goes wrong. It is a
dry run until `--execute`, exactly like `scripts/archive_runs.py`.

## The map after every step that can move it

A plan carries a `MAP:` line after each step from which the map is reachable --
`map_affecting` below derives that set from the same graph, so it needs no
second list. The point is a series: one frame per state of the world, so the
terrain, the lakes and the biomes can be watched changing through a pass rather
than only at the end. `maps/snapshot.py` keys a frame on what it was drawn from,
so a step that moves none of the map's four inputs costs a manifest row and no
render. `--no-maps` drops the lines.

## Why purge is a graph question

CLAUDE.md rule 7: an upstream change makes everything below it WORTHLESS rather
than stale, so the question after any change is "what is now worthless". That is
answerable only from a complete graph, which is what `config/pipeline.yaml` is
for -- and it is the reason `writes` has to be the complete artifact set rather
than a marker file. A row naming only its report purges its report and leaves
the .sra beside it, then says it succeeded.

A purge starts at the seed's OWN output and works down. A change to a step
invalidates what that step wrote before anything else, so a walk that only went
downstream answered "nothing" at every step nothing else `needs`, with a success
exit; it also made the completeness rule on `writes` untestable at every leaf,
since a leaf could name its whole output set and still purge none of it.

Two edges are not `needs` and purge uses both:

`reads_export` is the step consuming `source/{build}/`. Four steps do it,
each declaring only the `needs` of its own pass -- which read as "depends on
nothing" when it meant "depends on the terrain". `--purge orogen` seeds
from those four rather than from `orogen`'s own `needs`.

And the traversal STOPS AT `orogen` rather than passing through it. Loop A is a
cycle -- carve_list feeds orogen, orogen's export feeds hydrography -- so
unrestricted reachability from any climate step reaches the generation step and
then everything, including `source/{build}/` itself. That answer is wrong, not
merely alarming: a new climatology does not invalidate the terrain it was
computed on. The cycle closes across an ITERATION boundary and "what is worthless
now" is a within-pass question, so the generation step is where it is cut. That
cut holds for `orogen` seeded from itself too: `source/` is read-only, so
`--purge orogen` clears what the generation invalidates and is the one seed that
does not take its own output with it.

It also does not check whether artifacts AGREE with each other. That is
`check_consistency.py`, which compares terrain hashes, provenance stamps and
conventions, and it is the better tool for it. The division is: consistency
answers "do these describe the same world", this answers "what exists and what
must run to reach X". Asking either to do the other's job would give two answers
to one question.

## What this does NOT read: the issue tracker

This module answers graph questions and nothing else. It does not read the
tracker, so nothing recorded in `bd` can change what it reports.

It used to, to print a carve gate -- "open tasks touching a step upstream of
`carve_list`" -- and that was wrong twice over. A regex over human-written
markdown sat inside the pipeline planner, so editing a status cell changed what
the planner said about the pipeline. And the gate answered a different question
from the one it named: whether a row still moves the verdict depends on what
closing it would CHANGE, which is the row's prose, while a `[step: <id>]` marker
records where the work is filed. No traversal turns the second into the first,
so the output looked decided while deciding nothing. The gate is a judgement
made by READING the open issues before carving; `config/pipeline.yaml` states
it against the `orogen` step.
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


def load() -> dict:
    return yaml.safe_load(GRAPH.read_text(encoding="utf-8"))


def active_build() -> str:
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    return cfg.get("source_build") or cfg.get("model", {}).get("source_build", "")


def steps_by_id(graph: dict) -> dict:
    return {s["id"]: s for s in graph["steps"]}


def spectrum_name() -> str:
    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    return str(cfg.get("radiation", {}).get("stellar_spectrum", "k25v"))


def model_resolution() -> str:
    """The configured ladder rung, checked against its own grid dimensions.

    `{res}` and `{res_lower}` in a `writes` entry stand for this. Every
    ExoPlaSim input family is named and rooted by resolution -- `inputs/t21/`
    holding `orogen_T21_surf_0172.sra` -- and the graph carried `t42` and
    `T42` as literals, so it declared artifacts that do not exist and declared
    none of the ones that do. SPAT-2.
    """
    import sys as _sys
    import yaml
    if str(ROOT / "lib") not in _sys.path:
        _sys.path.insert(0, str(ROOT / "lib"))
    import rungs
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    return rungs.model_grid(cfg)[0]


def _substitute(path: str, build: str) -> str:
    res = model_resolution()
    return (path.replace("{build}", build)
                .replace("{name}", spectrum_name())
                .replace("{res_lower}", res.lower())
                .replace("{res}", res))


def resolve(path: str, build: str) -> Path:
    return ROOT / _substitute(path, build)


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
        return sorted(ROOT.glob(_substitute(path, build)))
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

    STRICTLY downstream: the seed is not in its own result. `purge_set` below is
    what `--purge` deletes and it adds the seed back, because a step's own output
    is the first thing a change to that step invalidates. The two are kept apart
    so the cut at `orogen` is stated once, against reachability, where the cycle
    is.
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


def purge_set(seed: str, by_id: dict) -> set:
    """Every step `--purge seed` deletes: what `seed` invalidates, plus `seed`.

    A step's own output is the first thing a change to that step makes worthless,
    so it is the one case the tool must always get right. Reverse reachability
    alone gets it wrong at every leaf: a step nothing else `needs` has an empty
    downstream set, and `--purge` reported success having deleted none of the
    step's own artifacts. `writes` is required to name the COMPLETE output set
    precisely so a purge deletes exactly that, and a leaf could satisfy that rule
    perfectly while purging nothing. world-ysd1.

    `orogen` is the one exemption, BY NAME and for rule 7: its output is
    `source/{build}/`, which is read-only. A build is retired by adding its
    successor and stubbing it under `archive/builds/`, never by deleting the
    export in place, and `--purge orogen` means "clear what this generation
    invalidates" rather than "delete the terrain". The cut in `downstream` stops
    every OTHER seed from reaching the export; this stops the export's own seed.
    """
    out = downstream(seed, by_id)
    if seed != "orogen":
        out.add(seed)
    return out


def purge_writes(seed: str, graph: dict, by_id: dict) -> list[tuple[str, list[str]]]:
    """The `writes` entries `--purge seed` would expand, per step, in graph order.

    Declared entries, not files: what is on disk is the tree's business and this
    is the graph's. `cmd_purge` expands these and `smoke_test.py` asserts over
    them, so the plan a check reads is the plan the delete walks.
    """
    doomed = purge_set(seed, by_id)
    return [(s["id"], list(s.get("writes", [])))
            for s in graph["steps"] if s["id"] in doomed]


MAP_STEPS = ("basemap", "projections")


def map_affecting(by_id: dict) -> set:
    """Every step whose output can change the map.

    Ancestors over `needs`, plus the export edge: the map reads
    `source/{build}/` and declares `reads_export`, so `orogen` can change it and
    so can everything upstream of `orogen`. Loop A is a cycle, so that is most
    of loop A -- correctly. A carve verdict changes the next terrain, and the
    next terrain is the picture.

    NOT cut at `orogen` the way `downstream` is. That cut exists because "what
    is now worthless" is a within-pass question; this one is "what can move the
    picture", and it crosses the iteration boundary by design.

    A step in this set is one the map is rendered after. That does not mean the
    picture changes there: `maps/snapshot.py` keys a frame on the identity of
    what it was drawn from, so a step that moves none of those four inputs is
    recorded against the frame that already exists rather than rendering a
    second copy of it.
    """
    out: set[str] = set()
    frontier = set(MAP_STEPS)
    while frontier:
        nxt: set[str] = set()
        for sid in frontier:
            step = by_id.get(sid, {})
            nxt |= set(step.get("needs", []))
            if step.get("reads_export"):
                nxt.add("orogen")
        nxt -= out | set(MAP_STEPS)
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
    print(f"\n{len(graph['steps'])} steps, {miss_total} with a missing "
          f"artifact, {stale_total} present but built from an older config")
    return 0


def cmd_plan(graph: dict, target: str, force: bool, maps: bool = True) -> int:
    config = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    by_id = steps_by_id(graph)
    if target not in by_id:
        raise SystemExit(f"no step '{target}'. Known: {', '.join(sorted(by_id))}")
    build = active_build()
    seq = order(target, by_id)
    affecting = map_affecting(by_id) if maps else set()
    print(f"to make `{target}` current, in order:\n")
    hours, undecidable, snapshots = 0, 0, 0
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
        # rather than a reason to pass over it. Most hours-costing steps in the
        # graph are `by_index`, so treating undecidable as skippable silenced
        # most of the hour count and every gate attached to a run.
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
        if sid in affecting:
            snapshots += 1
            print(f"              MAP:  python maps/snapshot.py --step {sid}")
    print(f"\n{len(seq)} steps in the closure, {hours} of them cost hours.")
    if undecidable:
        print(f"{undecidable} of the steps listed are UUID-named, so whether they "
              f"still need to run is\nINDEX.json's answer and not the "
              f"filesystem's. The hour count assumes they do.")
    if snapshots:
        print(f"{snapshots} of them can move the map, so a frame is rendered after "
              f"each.\nA step that moves none of the four inputs the map is drawn "
              f"from is recorded\nagainst the frame that already exists rather "
              f"than redrawing it; --no-maps drops\nthe lines entirely.")
    print("This prints a plan. It does not run anything.")
    return 0


def cmd_purge(graph: dict, target: str, execute: bool) -> int:
    by_id = steps_by_id(graph)
    if target not in by_id:
        raise SystemExit(f"no step '{target}'. Known: {', '.join(sorted(by_id))}")
    build = active_build()
    doomed = purge_set(target, by_id)

    print(f"active build: {build}")
    own = (f"`{target}` included"
           if target in doomed
           else f"`{target}`'s own output is the read-only export and stays")
    print(f"worthless if `{target}` changes: {len(doomed)} of "
          f"{len(graph['steps'])} steps, {own}\n")

    # Runs are named by a UUID and have no path in the graph, so purge cannot
    # find them and must not pretend the tree is clear without them. They also
    # have their own extract-then-delete path, which verifies the archive before
    # removing anything -- reimplementing that here would be a second copy of the
    # one thing in this repository that must not have two.
    runs = sorted(s for s in doomed if by_id[s].get("presence") == "by_index")

    total, plan = 0, []
    for sid, writes in purge_writes(target, graph, by_id):   # graph order, top-down
        files = []
        for w in writes:
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
                  "      python scripts/archive_runs.py --include-live\n"
                  "  (that flag treats EVERY run on disk as dead, not only these;\n"
                  "  it is a dry run until --execute -- review its plan and spare\n"
                  "  survivors with --keep first)")

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
            if not d.is_dir():
                continue
            owners = [s["id"] for s in graph["steps"]
                      for w in s.get("writes", []) if w.endswith("/")
                      and resolve(w, build) == d]
            outside = [o for o in owners if o not in doomed]
            if outside:
                print(f"  kept {d} (also owned by {', '.join(outside)}, "
                      "which this purge does not reach)")
                continue
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
                    help="delete every artifact a change to STEP makes worthless, "
                         "STEP's own output included. Never crosses `orogen`, and "
                         "never deletes the export")
    ap.add_argument("--execute", action="store_true",
                    help="with --purge, actually delete; default is a dry run")
    ap.add_argument("--no-maps", action="store_true",
                    help="with --plan, leave out the map snapshot after each "
                         "step that can move the picture")
    args = ap.parse_args()
    graph = load()
    if args.register:
        raise SystemExit(cmd_register(graph))
    if args.purge:
        raise SystemExit(cmd_purge(graph, args.purge, args.execute))
    if args.plan:
        raise SystemExit(cmd_plan(graph, args.plan, args.force,
                                  maps=not args.no_maps))
    raise SystemExit(cmd_status(graph))


if __name__ == "__main__":
    main()
