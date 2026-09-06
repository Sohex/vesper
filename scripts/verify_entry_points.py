#!/usr/bin/env python3
"""Every script with an entry point starts: it imports, and its parser builds.

    python scripts/verify_entry_points.py
    python scripts/verify_entry_points.py --jobs 1          # serial
    python scripts/verify_entry_points.py exoplasim/scripts # a subtree only
    python scripts/verify_entry_points.py --self-test       # the classifier only

WHAT IT ASSERTS. Every `.py` under the script directories that has a `__main__`
block STARTS, in an interpreter of its own: its module body executes, so its
import graph resolves on THIS interpreter, and where it builds an argparse
parser the parser is constructed without raising. The second is not hypothetical
-- a strict per-build default raises while the parser is being built when the
build directory is missing, and that is the intended failure.

STARTING IS NOT RUNNING, AND THE PROBE IS CHOSEN PER SCRIPT SO THAT IT IS NOT.
`--help` short-circuits only where there is a parser to short-circuit it: on a
script that never builds one, `--help` is an argument nothing reads and the
script does its whole job. Twenty-seven entry points here are in that class, and
starting all of them once rewrote NINE TRACKED ARTIFACTS -- gate reports and
measurement products whose generated timestamp moved while every other byte
stayed put. A gate that mutates tracked state is a gate an agent cannot run
freely, and the churn lands in the next `git add`. So the probe is picked from
the source: `help` where the source constructs an `ArgumentParser`, `import`
where it does not. `import` executes the module body under a `__name__` that is
not `__main__`, which is the whole of what "does this start" asks and none of
what running it costs. `scripts/_entry_probe.py` is the child that does it.

A REFUSAL IS NOT A BROKEN IMPORT, AND THE EXIT STATUS CANNOT TELL THEM APART.
Several scripts here refuse on purpose: one will not mix a second build's
landform composition into a roughness, another reports an open seam. They start
perfectly well and exit non-zero, and reporting that as an entry point that does
not start teaches a reader to ignore this gate -- which costs exactly the
failure rule 8 promoted it to a gate for. The child therefore reports a CLASS
rather than a status, because inside the child the difference is an exception
type: `ImportError` for a graph that no longer resolves, `SystemExit` with a
non-zero code for a script that chose its own exit, anything else for a
statement or a parser that raised. Two of those fail this gate and the third is
reported beside them. A per-script allowlist would decide the same question by
naming names, and would go stale silently the moment a script changed its mind.
`--self-test` builds one script of each class in a temporary directory and
checks the classifier against them, including that the two which would write
when RUN do not write when STARTED. It also checks the static half -- that a
script the `help` probe is chosen for builds its parser before it does anything,
which is what makes that probe safe on it -- and every entry point is held to
that in a full pass. Both halves run in a full pass, so neither can rot
unexercised.

WHY IT IS NOT IN `smoke_test.py` ANY MORE. `smoke_test.py` is the gate every
session runs before every commit, so its cost is multiplied by every commit in
the project; six agents paying it at once is a measurable share of this host's
load. This pass is a fresh interpreter per entry point, each importing numpy,
netCDF4 and whatever else the module reaches for at import time, and its own
opt-out flag used to say it "dominates runtime". A smoke test answers "is this tree
coherent" from static reads in seconds. This is a real check with real value and
it is a PRE-BUILD and PRE-PUSH check, not a per-commit one, so it has its own
entry point and its own place in `CLAUDE.md` rule 8.

WHY NOT A STATIC PASS INSTEAD. Because a static pass answers a DIFFERENT
question, and saying otherwise would be the cheap gate pretending to be the
expensive one. `smoke_test.py` parses every module, runs pyflakes F821 over it,
and asks CPython's own symbol table which loaded names nothing binds; all three
are about NAME BINDING within a source text. None of them can say whether the
import graph resolves against the packages installed in this venv, whether a
module-level statement raises, or whether an argparse default blows up while
being evaluated. Importing is a dynamic act and only an import answers it.

WHY A FRESH INTERPRETER PER SCRIPT, when one process importing all of them would
be far cheaper. Isolation, and it is load-bearing rather than tidy. Many scripts
here begin by inserting `lib/` or their own directory on `sys.path`; run them in
one process and the second script's `import rungs` succeeds because the first
one edited the path, not because the second one is correct. Module basenames
also repeat across the script directories, so `sys.modules` would hand the
second one the first one's module. A single process would report a pass this
tree has not earned, which is worse than the cost it saves. It stays one child
per entry point: choosing the probe rather than adding one does not raise the
count.

PARALLELISM IS DELIBERATELY MODEST. This host runs several agents at once, so
this takes a quarter of the logical cores and no more; `--jobs` overrides it.
The work is subprocesses, so threads carry it and the GIL is not in the way.
Anything heavier than this runs under `scripts/lock_and_run` -- see CLAUDE.md.
This does not, because a quarter of the cores for a few seconds is not the kind
of load that lock exists to serialise. It does, however, verify that wrapper:
`lock_and_run` carries no `.py` extension, so it falls outside the sweep below
and this is the only gate that starts it.

IT IS NOT THE COMPILE GATE'S SIBLING BY ACCIDENT.
`exoplasim/scripts/verify_model_compiles.py` is the same shape for the model
Fortran -- a process per translation unit, too expensive per commit, run before
paying for a build. Rule 8 names them together.
"""

from __future__ import annotations

import argparse
import ast
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap

# ONE CORE PER PROBE. `--help` is not a numerical workload, but most of these
# scripts import numpy on the way to their parser, and numpy's bundled OpenBLAS
# opens a spinning thread pool sized to the logical core count. Multiplied by
# the worker count that is the whole machine, for nothing. Sized when the
# library loads, so it has to be in the environment the child starts with.
PROBE_ENV = dict(os.environ)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    PROBE_ENV.setdefault(_v, "1")
PROBE_ENV.setdefault("OMP_WAIT_POLICY", "passive")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from smoke_test import SCRIPT_DIRS  # noqa: E402  -- one declaration of the tree

ROOT = Path(__file__).resolve().parents[1]
PROBER = Path(__file__).resolve().parent / "_entry_probe.py"

# Entry points that must not be started at all: `rebuild_binaries.py` recompiles
# the model, which neither probe below makes cheap or side-effect free.
SKIP = {"rebuild_binaries.py"}

# The classes `_entry_probe.py` reports back. Two of them fail this gate; the
# third is a script that started and said no, which is not this gate's business.
FAILING = ("import", "start")
REPORTED = ("refusal",)


def probe_mode(source: str) -> str:
    """`help` where there is a parser to short-circuit on, `import` otherwise.

    Read out of the source on every pass, so a script that gains or loses its
    parser moves class the same day rather than waiting for a list to be edited.
    """
    return "help" if "ArgumentParser" in source else "import"


def _builds_parser_first(body: list[ast.stmt]) -> bool:
    """Is the first thing this block does the construction of its parser?

    The statements a parser is allowed to follow are the ones that cannot do
    work: imports, a docstring, a `global` declaration, and the
    `sys.path.insert` a script uses to reach its own directory.
    """
    for stmt in body:
        text = ast.unparse(stmt)
        if "ArgumentParser" in text:
            return True
        if isinstance(stmt, (ast.Import, ast.ImportFrom, ast.Global)):
            continue
        if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Constant):
            continue
        if isinstance(stmt, ast.Expr) and text.startswith("sys.path.insert("):
            continue
        return False
    return False


def help_short_circuits(path: Path, source: str) -> str:
    """`--help` on this script reaches its parser before it does any work.

    The `help` probe is safe BECAUSE argparse exits at `--help` before the rest
    of the start path runs, and that is a property of the SCRIPT rather than of
    argparse: a parser built after a read, a solve or a write short-circuits
    nothing that matters, and starting such a script would run all of it. This
    gate chooses that probe, so this gate is where the choice is checked. A
    static read, so it fails on the source rather than after the artifact has
    already been rewritten. "" is a pass.

    Building the parser behind a helper is not a defect in the script; it is a
    shape this check cannot see through, and a probe that cannot be shown to be
    safe is not one to run. Move the construction into the start path.
    """
    blocks = []
    for n in ast.parse(source).body:
        if isinstance(n, ast.FunctionDef) and n.name == "main":
            blocks.append(("main()", n.body))
        if isinstance(n, ast.If) and "__main__" in ast.unparse(n.test):
            blocks.append(("its __main__ block", n.body))
    if any(_builds_parser_first(body) for _, body in blocks):
        return ""
    where = blocks[0][0] if blocks else "this script"
    return (f"{path.name}: {where} does not build its ArgumentParser before "
            f"anything else, so --help cannot be shown to short-circuit "
            f"before the work")


def entry_points(roots: list[Path]) -> list[tuple[Path, str]]:
    """Every script under `roots` worth starting, with the probe it takes.

    `_`-prefixed files are not entry points: they are the path helpers imported
    by their neighbours, and `_entry_probe.py`, which is this gate's own child.
    A file with no `__main__` has no entry point to check.
    """
    out = {}
    for d in roots:
        if not d.is_dir():
            continue
        for f in d.glob("*.py"):
            if f.name in SKIP or f.name.startswith("_"):
                continue
            src = f.read_text(encoding="utf-8", errors="replace")
            if "__main__" not in src:
                continue
            out[f] = probe_mode(src)
    return sorted(out.items())


def classify(path: Path, mode: str) -> dict:
    """Start one script in an interpreter of its own and name what happened.

    Returns the child's verdict: `kind` in ok / import / start / refusal. A
    child that dies without one -- a timeout, a signal, a crash in an extension
    -- did not start, and says so under the class its silence belongs to.
    """
    argv = [sys.executable, str(PROBER), "--mode", mode, str(path)]
    try:
        r = subprocess.run(argv, capture_output=True, text=True, timeout=120,
                           cwd=ROOT, env=PROBE_ENV)
    except subprocess.TimeoutExpired:
        return {"kind": "start", "detail": "no verdict within 120 s", "mode": mode}
    lines = [l for l in r.stdout.splitlines() if l.strip()]
    try:
        return json.loads(lines[-1])
    except (IndexError, ValueError):
        tail = (r.stderr.strip().splitlines() or ["(no output)"])[-1]
        return {"kind": "start",
                "detail": f"the probe returned {r.returncode} without a "
                          f"verdict: {tail}",
                "mode": mode}


def lock_wrapper(verbose: bool) -> list[str]:
    """`scripts/lock_and_run --self-test`, which is the only thing that runs it.

    The wrapper has no `.py` extension, deliberately, because the invocation is
    a command prefix a caller types. That puts it outside `entry_points()`, so
    without this it would be the one script in the tree nothing ever starts --
    and it is the script every expensive run goes through. A lock that has
    quietly stopped excluding is invisible until two integrations land on the
    host at once.

    Its own arms are behavioural rather than static: streams and exit status
    through, a second command excluded until the first finishes, a nested
    wrapper running instead of deadlocking, a dead wrapper's claim cleared and a
    hand-rolled one left alone. It runs against a private lock path and never
    touches `/tmp/world.lock`.
    """
    tool = ROOT / "scripts" / "lock_and_run"
    if not tool.is_file():
        return ["scripts/lock_and_run is missing, and every expensive run in "
                "this tree is documented to go through it"]
    if not os.access(tool, os.X_OK):
        return ["scripts/lock_and_run is not executable, so the documented "
                "invocation does not run"]
    r = subprocess.run([sys.executable, str(tool), "--self-test"],
                       capture_output=True, text=True, timeout=180, cwd=ROOT,
                       env=PROBE_ENV)
    if verbose:
        for line in r.stdout.splitlines():
            print(f"  {line}")
    if r.returncode == 0:
        return []
    failed = [l.strip() for l in r.stdout.splitlines() if l.strip().startswith("FAIL")]
    return [f"scripts/lock_and_run --self-test: {f}" for f in failed] or [
        "scripts/lock_and_run --self-test failed without naming an arm: "
        + (r.stderr.strip().splitlines() or ["(no output)"])[-1]]


def worktree_links_self_test(verbose: bool) -> list[str]:
    """`scripts/check_worktree_links.py --self-test`.

    Starting that script proves it imports; it does not prove the ledger
    comparison still discriminates. The self-test builds a throwaway repository
    of the shape the hazard needs -- a generated ignored file beside tracked
    source, a wholly-ignored directory, files either side of the hash budget --
    adds A REAL GIT WORKTREE to it, and answers nine questions, including the
    two that decide whether the guard is usable at all: an identical-bytes
    rewrite must NOT report, and re-linking must re-baseline so a legitimate
    regeneration is accepted rather than refused.

    Here rather than in `smoke_test.py` because it creates a worktree and
    spawns git, which is the line rule 8 draws between the two tiers, and
    because this file already runs the one other self-test of that shape.
    """
    r = subprocess.run([sys.executable, "scripts/check_worktree_links.py", "--self-test"],
                       capture_output=True, text=True, timeout=180, cwd=ROOT,
                       env=PROBE_ENV)
    if verbose:
        for line in r.stdout.splitlines():
            print(f"  {line}")
    if r.returncode == 0:
        return []
    failed = [l.strip() for l in r.stdout.splitlines() if l.strip().startswith("FAIL")]
    return [f"scripts/check_worktree_links.py --self-test: {f}" for f in failed] or [
        "scripts/check_worktree_links.py --self-test failed without naming a case: "
        + (r.stderr.strip().splitlines() or ["(no output)"])[-1]]


def export_refuses_an_existing_build(verbose: bool) -> list[str]:
    """`vendor/orogen/tools/export-planet.mjs` refuses to overwrite an export.

    Rule 7's one enforcement point. `source/` is what the world is
    reconstructed from, git tracks none of its payload, and the writer is the
    vendored Node exporter -- so the refusal lives in the exporter, where it
    fires in the main checkout and in a worktree alike, before the generation
    that is almost all of what an export costs. A guard that fired only through
    a Python wrapper would be bypassed by the invocation `source/README.md`
    documents, which is the exporter itself.

    HERE RATHER THAN IN `smoke_test.py` because it spawns Node and, in the
    negative arm, generates a planet: that is rule 8's line between the tiers.
    It is also the tier that matters, since the thing it protects is paid for
    just before a build.

    BOTH ARMS, and each decides by an observable that answers the question
    asked rather than by an exit status. The positive arm requires the refusal
    to name the directory it refused; the negative arm requires an ordinary
    export into an empty directory to COMPLETE and leave a manifest, because a
    guard that refused every path would pass the positive arm alone and would
    stop every export in the project. 2,000 regions and one T21 grid put the
    pair at about a second.
    """
    import shutil

    exporter = ROOT / "vendor" / "orogen" / "tools" / "export-planet.mjs"
    if not exporter.is_file():
        return [f"{exporter.relative_to(ROOT)} is missing, and it is the only "
                "thing that refuses an overwrite of an existing source/ build"]
    node = shutil.which("node")
    if node is None:
        return ["node is not on this host, so the one refusal that enforces "
                "rule 7 at the writer cannot be exercised. See "
                "docs/src/reference/environment.md"]
    # Named separately, because without it the negative arm below fails on a
    # missing import and reads as a broken refusal. `link_worktree.py` links it
    # -- it is an install and not a build of tracked source -- so its absence
    # means an unlinked worktree rather than a defect in the exporter.
    if not (exporter.parent.parent / "node_modules").exists():
        return ["vendor/orogen/node_modules is absent, so the exporter cannot "
                "run and its refusal cannot be exercised. Run "
                "`python scripts/link_worktree.py` in a worktree, or "
                "`npm i delaunator` under vendor/orogen"]

    cwd = exporter.parent.parent
    bad: list[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        held = Path(tmp) / "already-a-build"
        held.mkdir()
        (held / "manifest.json").write_text("{}", encoding="utf-8")
        r = subprocess.run(
            [node, "tools/export-planet.mjs", "--seed", "1", "--regions", "2000",
             "--grid", "T21", "--no-raw", "--quiet", "--out", str(held)],
            cwd=cwd, capture_output=True, text=True, timeout=300, env=PROBE_ENV)
        if verbose:
            print(f"  refusal arm: exit {r.returncode}")
        if str(held) not in r.stderr or "already holds an export" not in r.stderr:
            bad.append("export-planet.mjs did not refuse an --out directory "
                       "that already holds a manifest, so nothing stops an "
                       f"existing source/ build being overwritten: {r.stderr.strip()[-200:]}")
        elif (held / "README.txt").exists() or (held / "grid").exists():
            bad.append("export-planet.mjs refused and wrote into the directory "
                       "anyway")

        fresh = Path(tmp) / "new-build"
        r = subprocess.run(
            [node, "tools/export-planet.mjs", "--seed", "1", "--regions", "2000",
             "--grid", "T21", "--no-raw", "--quiet", "--out", str(fresh)],
            cwd=cwd, capture_output=True, text=True, timeout=300, env=PROBE_ENV)
        if verbose:
            print(f"  ordinary arm: exit {r.returncode}")
        if not (fresh / "manifest.json").is_file():
            bad.append("an ordinary export into an empty directory did not "
                       "write a manifest, so the refusal above is refusing "
                       f"everything: {(r.stderr.strip() or '(no output)')[-200:]}")
    return [f"export-planet.mjs: {b}" for b in bad]


def verify(roots: list[Path], jobs: int, verbose: bool) -> tuple[list[str], list[str]]:
    """Start every entry point under `roots`. Returns (failures, refusals)."""
    files = entry_points(roots)
    helps = sum(1 for _, m in files if m == "help")
    print(f"{len(files)} entry points, {jobs} at a time "
          f"({helps} started with --help, {len(files) - helps} imported only)")
    unsafe = [help_short_circuits(f, f.read_text(encoding="utf-8",
                                                 errors="replace"))
              for f, m in files if m == "help"]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        verdicts = list(pool.map(lambda fm: classify(*fm), files))
    failures = [f"[probe] {u}" for u in unsafe if u]
    refusals: list[str] = []
    for (f, _), v in zip(files, verdicts):
        line = f"{f.relative_to(ROOT)}: {v['detail']}"
        if v["kind"] in FAILING:
            failures.append(f"[{v['kind']}] {line}")
        elif v["kind"] in REPORTED:
            refusals.append(line)
    if verbose:
        for (f, mode), v in zip(files, verdicts):
            mark = {"ok": " ok ", "refusal": "said no"}.get(v["kind"], "FAIL")
            print(f"  {mark:>7}  {mode:6s}  {f.relative_to(ROOT)}")
    return failures, refusals


# One script per class the classifier claims to separate, plus the two whose
# whole point is that starting them must not RUN them. Each row is (filename,
# source, expected class, expected side effect). The side effect column is the
# property world-ujcj is about: a script that writes when started leaves
# `<name>.written` beside itself, and the two rows that would write are the two
# the gate must reach without writing.
SELF_TEST_CASES = [
    ("broken_import.py", """
        import a_package_this_tree_has_never_had

        def main():
            print(a_package_this_tree_has_never_had)

        if __name__ == "__main__":
            main()
     """, "import", False),
    ("refuses_at_module_level.py", """
        import sys
        from pathlib import Path

        raise SystemExit("refusing: these band shares belong to another build")

        if __name__ == "__main__":
            pass
     """, "refusal", False),
    ("raises_at_module_level.py", """
        TABLE = {"a": 1}
        MISSING = TABLE["b"]

        if __name__ == "__main__":
            print(MISSING)
     """, "start", False),
    ("gate_refuses_in_main.py", """
        from pathlib import Path

        def main():
            Path(__file__).with_suffix(".written").write_text("ran\\n")
            raise SystemExit(1)

        if __name__ == "__main__":
            main()
     """, "ok", False),
    ("healthy_parser.py", """
        import argparse
        from pathlib import Path

        def main():
            ap = argparse.ArgumentParser(description="a well-formed parser")
            ap.add_argument("--n", type=int, default=1)
            a = ap.parse_args()
            Path(__file__).with_suffix(".written").write_text(f"{a.n}\\n")

        if __name__ == "__main__":
            main()
     """, "ok", False),
    ("parser_default_raises.py", """
        import argparse
        from pathlib import Path

        def strict_default():
            raise FileNotFoundError("no build directory to default to")

        def main():
            ap = argparse.ArgumentParser()
            ap.add_argument("--build", default=strict_default())
            ap.parse_args()

        if __name__ == "__main__":
            main()
     """, "start", False),
    ("refuses_before_parser.py", """
        import argparse

        def main():
            raise SystemExit("refusing: the SPAT-5 seam is open")
            argparse.ArgumentParser().parse_args()

        if __name__ == "__main__":
            main()
     """, "refusal", False),
]

# The static half: does `--help` on this source reach the parser before the
# work? Each row is (name, source, expected complaint substring; "" for a pass).
SHORT_CIRCUIT_CASES = [
    ("parser_first.py", """
        import argparse
        from pathlib import Path

        def main():
            ap = argparse.ArgumentParser()
            Path("out").write_text("ran")

        if __name__ == "__main__":
            main()
     """, ""),
    ("parser_after_a_read.py", """
        import argparse
        from pathlib import Path

        def main():
            mesh = Path("mesh").read_text()
            ap = argparse.ArgumentParser()
            ap.parse_args()

        if __name__ == "__main__":
            main()
     """, "does not build its ArgumentParser before anything else"),
    ("parser_in_the_guard.py", """
        def work():
            return 1

        if __name__ == "__main__":
            import argparse
            import sys
            sys.path.insert(0, ".")
            argparse.ArgumentParser().parse_args()
     """, ""),
    ("parser_out_of_reach.py", """
        import argparse

        def build():
            return argparse.ArgumentParser()

        def main():
            build().parse_args()

        if __name__ == "__main__":
            main()
     """, "does not build its ArgumentParser before anything else"),
]


def self_test(verbose: bool) -> list[str]:
    """The classifier against one script of each class, in a temporary tree.

    Every arm can fail and says how: a broken import that came back anything but
    `import` would mean this gate no longer catches the failure rule 8 promoted
    it to a gate for, and a refusal that came back `import` would mean it is
    back to calling a gate's correct no a broken module. The two `ok` rows carry
    a second assertion -- that the script did not WRITE -- because "started" and
    "ran" are the distinction the probe exists to hold apart, and a row that
    started passing by running would still report `ok`.
    """
    problems = []
    with tempfile.TemporaryDirectory(prefix="entry-probe-") as tmp:
        d = Path(tmp)
        for name, body, expect, writes in SELF_TEST_CASES:
            f = d / name
            src = textwrap.dedent(body).lstrip()
            f.write_text(src, encoding="utf-8")
            mode = probe_mode(src)
            v = classify(f, mode)
            wrote = f.with_suffix(".written").exists()
            ok = v["kind"] == expect and wrote == writes
            if verbose:
                print(f"  {' ok ' if ok else 'FAIL'}  {mode:6s}  {name} -> "
                      f"{v['kind']}{' and wrote' if wrote else ''}")
            if v["kind"] != expect:
                problems.append(f"--self-test {name}: expected class "
                                f"{expect}, got {v['kind']} ({v['detail']})")
            if wrote != writes:
                problems.append(f"--self-test {name}: starting it "
                                f"{'wrote' if wrote else 'did not write'} "
                                f"beside itself, which is not what this probe "
                                f"is for")
    for name, body, expect in SHORT_CIRCUIT_CASES:
        src = textwrap.dedent(body).lstrip()
        got = help_short_circuits(Path(name), src)
        ok = (expect in got) if expect else not got
        if verbose:
            print(f"  {' ok ' if ok else 'FAIL'}  static  {name} -> "
                  f"{'short-circuits' if not got else got.split(': ', 1)[-1]}")
        if not ok:
            problems.append(f"--self-test {name}: expected "
                            f"{expect or 'no complaint'}, got "
                            f"{got or 'no complaint'}")
    return problems


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("roots", nargs="*", type=Path,
                    help="directories to check (default: every script "
                         "directory smoke_test.py declares)")
    ap.add_argument("--jobs", type=int, default=None,
                    help="workers (default: a quarter of the logical cores, "
                         "so a shared host keeps most of itself)")
    ap.add_argument("--verbose", action="store_true",
                    help="print each entry point with its verdict")
    ap.add_argument("--self-test", action="store_true",
                    help="check the classifier against one script of each "
                         "class, and stop")
    a = ap.parse_args()
    if a.self_test:
        problems = self_test(verbose=True)
        print(f"{len(problems)} classifier arms failed" if problems
              else "the classifier separates every class")
        for p in problems:
            print(p, file=sys.stderr)
        raise SystemExit(1 if problems else 0)

    roots = [ROOT / r for r in a.roots] if a.roots else SCRIPT_DIRS
    jobs = a.jobs if a.jobs else max(1, (os.cpu_count() or 4) // 4)
    problems, refusals = verify(roots, jobs, a.verbose)
    # This gate's own two self-checks belong to the repository's `scripts/`, not
    # to any directory that happens to be named one: `<component>/scripts` is
    # spelled the same and must not drag them in.
    if any(r.resolve() == ROOT / "scripts" for r in roots):
        problems += self_test(a.verbose)
        problems += lock_wrapper(a.verbose)
        problems += worktree_links_self_test(a.verbose)
        problems += export_refuses_an_existing_build(a.verbose)
    if refusals:
        print(f"\n{len(refusals)} entry points started and refused, which is "
              f"not this gate's business:")
        for r in refusals:
            print(f"  {r}")
    for p in problems:
        print(p, file=sys.stderr)
    print(f"\n{len(problems)} entry points do not start" if problems
          else "\nevery entry point starts")
    raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()
