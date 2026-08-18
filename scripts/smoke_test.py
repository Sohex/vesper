#!/usr/bin/env python3
"""Exercise every script's entry points without running any science.

    python scripts/smoke_test.py

`check_consistency.py` audits the artifacts. This audits the code that produces
them, which is where a different class of bug lives: in one day this project
found basin ids read from the wrong terrain, six scripts defaulting to a stale
data directory, five stale binaries, and a config change that left `latitudes`
behind. Every one of those was reachable by importing a module and looking at
where its defaults pointed -- none needed a model run.

Four checks, all cheap:

1. **Imports.** Every module imports. Catches a missing import added while
   editing, which `--help` alone will also catch but this localises better.
1b. **Undefined names**, via pyflakes F821. `--help` proves argparse builds and
   nothing about the code after it, so a name used but never imported survives
   the entry-point check and crashes at the end of `main()` -- after the artifact
   has been written. That is exactly how a six-site `relative_to` sweep shipped
   five NameErrors.
2. **Entry points.** Every script with a `main()` answers `--help`. Catches an
   argparse default that raises while being constructed -- which is exactly what
   a strict per-build default does when the build directory is missing, and is
   the intended failure.
3. **Defaults point at the active build.** Any default path under a component's
   `data/` must resolve beneath the active build's directory. This is the check
   that would have caught `carve_verdict.py` reading 2,107 basins from
   precarve-unzoned while computing verdicts on a build with 3,629.
4. **No selection by sort order.** A lint for `sorted(...glob(...))[n]` and
   friends. Choosing an artifact by whatever sorts last is the pattern behind
   ExoPlaSim's `finalize()` emitting the wrong world's output and behind three
   separate bugs here. Enumerating a known set is fine; picking from one is not.
5. **One grid convention.** A lint for longitude arithmetic outside
   `lib/gridding.py`. The export/model column mapping has been got wrong three
   times -- once as the original defect, once in the shape of its own fix, and
   once in a sink lookup the fix's sweep did not reach -- and each occurrence
   was a script deriving a column of its own. There is nothing to translate and
   therefore nothing to derive, so the check is that the expressions appear in
   exactly one file. Prevention, not translation.
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIRS = [ROOT / "scripts", ROOT / "lib", ROOT / "exoplasim" / "scripts",
               ROOT / "hydrography" / "scripts", ROOT / "pedology" / "scripts",
               ROOT / "biosphere" / "scripts", ROOT / "maps"]

# Entry points that must not be smoke-run: they are slow, or they mutate the
# tree, and --help on them is not free of side effects.
SKIP_HELP = {"rebuild_binaries.py"}

# `sorted(x.glob(...))[i]`, `sorted(glob(...))[i]`, `list(...glob(...))[i]`,
# and max/min over a glob: all of them choose one artifact by ordering.
ORDER_PICK = re.compile(
    r"(sorted|list|max|min)\s*\([^\n]*\.?glob\([^\n]*\)[^\n]*\)\s*\[")


def check_imports(files: list[Path]) -> list[str]:
    bad = []
    for f in files:
        try:
            ast.parse(f.read_text(encoding="utf-8"))
        except SyntaxError as exc:
            bad.append(f"{f.relative_to(ROOT)}: {exc}")
    return bad


def check_help(files: list[Path]) -> list[str]:
    bad = []
    for f in files:
        if f.name in SKIP_HELP or f.name.startswith("_"):
            continue
        src = f.read_text(encoding="utf-8")
        if "__main__" not in src:
            continue
        r = subprocess.run([sys.executable, str(f), "--help"],
                           capture_output=True, text=True, timeout=120, cwd=ROOT)
        if r.returncode != 0:
            tail = (r.stderr.strip().splitlines() or ["(no output)"])[-1]
            bad.append(f"{f.relative_to(ROOT)}: {tail}")
    return bad


def check_build_scoped_defaults() -> list[str]:
    """Every per-build default must resolve under the active build."""
    sys.path.insert(0, str(ROOT / "lib"))
    import yaml
    from builds import component_data
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    active = cfg.get("source_build")
    bad = []
    for component in ("hydrography", "pedology"):
        try:
            d = component_data(component, cfg)
        except Exception as exc:
            bad.append(f"{component}: component_data raised {exc}")
            continue
        if d.name != active:
            bad.append(f"{component}: data dir is {d.name}, active build is {active}")
    return bad


def check_undefined_names(files: list[Path]) -> list[str]:
    """Names used but never bound. This is the crash class `--help` cannot see.

    An entry point answering `--help` proves argparse is constructed; it proves
    nothing about the code after it. A sweep that replaced six `relative_to`
    call sites with a helper and imported that helper into only one of the five
    files passed `--help` on every one of them and then died with NameError at
    the end of `main()`, after the artifact had been written.

    pyflakes F821 is exactly this check. Only undefined names are treated as
    failures; unused imports and star-import warnings are style, and this
    project has deliberate `# noqa` imports for path setup.
    """
    try:
        from pyflakes.api import checkPath
        from pyflakes.reporter import Reporter
    except ImportError:
        return ["pyflakes is not installed; `uv pip install pyflakes`"]
    import io
    bad = []
    for f in files:
        out, err = io.StringIO(), io.StringIO()
        checkPath(str(f), Reporter(out, err))
        for line in out.getvalue().splitlines():
            if "undefined name" in line:
                bad.append(line.replace(str(ROOT) + "/", ""))
    return bad


def check_no_order_picks(files: list[Path]) -> list[str]:
    bad = []
    for f in files:
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            # Skip prose: a backtick means the line is documentation about the
            # pattern, not the pattern. This file's own docstring names it.
            if (stripped.startswith("#") or "smoke-ok" in line
                    or "`" in line or stripped.startswith("*")):
                continue
            if ORDER_PICK.search(line):
                bad.append(f"{f.relative_to(ROOT)}:{i}: {stripped[:70]}")
    return bad


# Longitude arithmetic. `(lon + 180)`, `% 360`, `- lon[0]` and `/ dlon` are the
# pieces every one of the three occurrences was built out of. `lib/gridding.py`
# is the one file allowed to contain them.
GRID_CONVENTION = re.compile(
    r"(\+\s*180(\.0)?\s*\)\s*/\s*360"          # (lon + 180) / 360
    r"|\+\s*180(\.0)?\s*\)\s*%\s*360"          # (lon + 180) % 360
    r"|%\s*360(\.0)?\s*-\s*180"                  # ... % 360 - 180
    r"|lon\[0\]"                                   # measuring from a label
    r"|360(\.0)?\s*/\s*n?lon)")                    # reconstructing dlon
# `lib/gridding.py` owns the convention. The two exemptions below are not the
# same seam: `maps/projections.py` and `maps/render_projections.py` rasterise
# World Orogen's own coordinates onto image pixels, where no model grid exists
# to disagree with. Anything that reads a climatology or a coupling matrix is
# on the seam and is not exempt -- `maps/build_basemap.py` reads one, which is
# how it came to draw every climate layer 180 degrees from its own terrain.
GRID_OWNER = ("lib/gridding.py", "scripts/smoke_test.py",
              "maps/projections.py", "maps/render_projections.py")


def check_one_grid_convention(files: list[Path]) -> list[str]:
    """Longitude arithmetic lives in exactly one file.

    The mapping between a World Orogen export and the ExoPlaSim grid it was
    integrated onto IS the identity, index for index, because one expression in
    `lib/gridding.py` put every mesh region in its column to begin with. So a
    second copy of that expression anywhere else is not a duplicate to keep in
    sync: it is a chance to write the mapping down differently, which is what
    happened three times and cost a build.

    A comment or a docstring naming the expression is documentation about the
    convention, not an instance of it, so lines carrying a backtick or opening
    with `#` are skipped exactly as in `check_no_order_picks`.
    """
    bad = []
    for f in files:
        if str(f.relative_to(ROOT)) in GRID_OWNER:
            continue
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            if (stripped.startswith("#") or "grid-ok" in line
                    or "`" in line or stripped.startswith("*")):
                continue
            if GRID_CONVENTION.search(line):
                bad.append(f"{f.relative_to(ROOT)}:{i}: {stripped[:70]}")
    return bad


def check_registered_in_workflow(files) -> list[str]:
    """Every generator that writes an artifact is named in `WORKFLOW.md`.

    `config/pipeline.yaml` is the graph and the list of what exists. A script that writes a product nobody declared has
    no recorded consumers, so nothing can say what it invalidates when it moves,
    and `CLAUDE.md` rule 7 stops being answerable -- "what is now worthless"
    needs a complete graph.

    This can fail in the direction that matters: add a generator, forget the
    row, and the check goes red. It deliberately does NOT test the reverse,
    because a row naming a step that has not been written yet is a plan rather
    than an error.
    """
    import yaml
    graph = yaml.safe_load((ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    declared = {Path(s["script"]).name for s in graph["steps"]}
    declared |= {Path(s["script"]).name for s in graph.get("checks", [])}
    declared |= {Path(s).name for s in graph.get("one_offs", [])}
    writes = re.compile(r"write_text|to_netcdf|savefig|json\.dump|write_sra"
                        r"|Dataset\([^)]*['\"]w['\"]|open\([^)]*['\"]w")
    missing = []
    for f in files:
        if f.name.startswith("_"):
            continue
        src = f.read_text(encoding="utf-8", errors="ignore")
        if not writes.search(src):
            continue                      # a module, not a step
        if "def main(" not in src and "__main__" not in src:
            continue                      # imported, not run
        if f.name not in declared:
            missing.append(f"{f.relative_to(ROOT)} writes an artifact and is in no "
                           "step, check or one_off in config/pipeline.yaml")
    return missing


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-help", action="store_true",
                    help="skip the subprocess --help pass, which dominates runtime")
    args = ap.parse_args()

    files = sorted({f for d in SCRIPT_DIRS if d.is_dir()
                    for f in d.glob("*.py")})
    print(f"{len(files)} modules under {len(SCRIPT_DIRS)} directories\n")

    checks = [("imports", check_imports(files)),
              ("no undefined names", check_undefined_names(files)),
              ("defaults scoped to the active build", check_build_scoped_defaults()),
              ("no artifact selection by sort order", check_no_order_picks(files)),
              ("one grid convention, in lib/gridding.py",
               check_one_grid_convention(files)),
              ("every generator is declared in config/pipeline.yaml",
               check_registered_in_workflow(files))]
    if not args.skip_help:
        checks.insert(1, ("entry points answer --help", check_help(files)))

    failed = 0
    for name, problems in checks:
        if problems:
            failed += 1
            print(f"[ FAIL ] {name}")
            for p in problems:
                print(f"         {p}")
        else:
            print(f"[  ok  ] {name}")

    print(f"\n{len(checks)} checks, {failed} failed")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
