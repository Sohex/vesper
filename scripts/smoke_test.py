#!/usr/bin/env python3
"""Is this tree coherent? Every check is a static read, and none of them runs science.

    python scripts/smoke_test.py

`check_consistency.py` audits the artifacts. This audits the code that produces
them, which is where a different class of bug lives: in one day this project
found basin ids read from the wrong terrain, six scripts defaulting to a stale
data directory, five stale binaries, and a config change that left `latitudes`
behind. Every one of those was reachable by reading a module and looking at
where its defaults pointed -- none needed a model run.

THE COST BOUND IS PART OF WHAT THIS FILE IS. It is the gate every session runs
before every commit, so its cost is multiplied by every commit in the project,
and this host runs several sessions at once. So a check belongs here when it is
a parse, a text read, a config comparison, or a small synthetic fixture; a check
that spawns a process per unit of the tree does not, however true it is. Two
such passes used to run here, both behind their own opt-out flag -- and the
`--help` flag's own help text said it "dominates runtime", which was the file
admitting the problem rather than fixing it. They are now their own gates,
neither weakened and neither deleted:

* **`scripts/verify_entry_points.py`** -- every script with a `__main__` answers
  `--help` in an interpreter of its own, about 136 of them. It catches an
  argparse default that raises while being constructed, which is what a strict
  per-build default does when the build directory is missing, and an import
  graph that no longer resolves. Nothing static answers either question, and the
  checks below do NOT stand in for it: they are all about name binding inside a
  source text.
* **`exoplasim/scripts/verify_model_compiles.py`** -- `gfortran -fsyntax-only`
  over every translation unit `plasim/CMakeLists.txt` names. Check 12f below
  catches the one class that has actually cost a build, statically and in the
  commit that makes it; the compile gate is the authority on whether the source
  is legal at all.

Both run before a build or a push, where a per-unit process is the cheap half of
what is about to be paid. `CLAUDE.md` rule 8 names them.

The checks, all cheap (plus registered-script existence and the two purge
properties -- never reaching the terrain, and always reaching the seed's own
writes -- run from `main()` with the rest):

1. **Every module parses.** An `ast.parse` of every source, which localises a
   syntax error to the file that has it. It does not IMPORT: nothing here does,
   and `verify_entry_points.py` is what answers that question.
1b. **Undefined names**, via pyflakes F821. A name used but never imported is
   bound by nothing and crashes at the end of `main()` -- after the artifact has
   been written. That is exactly how a six-site `relative_to` sweep shipped five
   NameErrors, and no entry-point check can see it: `--help` proves argparse was
   constructed and nothing about the code after it.
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
6. **Every generator is declared in `config/pipeline.yaml`.** An artifact no step
   writes has no derivable consumers, so nothing can say what it invalidates.
7. **Every step is named in its component README.** A step that works and is
   invisible gets reimplemented beside itself.
8. **`local_slope_deg` reproduces an analytic gradient**, on a synthetic
   icosphere. The tolerance is set between two plausible implementations rather
   than picked: see the docstring on the check itself.
9. **The convergence window follows the declared purposes.** Synthetic segment
   manifests against `segments.py:production_window`, including the live case: a
   diagnostic tail that the old trailing-window rule would have averaged into a
   verdict.
10. **No imported module name is rebound.** `import climatology` at the top
   and `climatology = <a Path>` in `main()` leaves every `climatology.foo()`
   after it calling a Path method. The name is defined and the import is live,
   so neither the import check nor pyflakes sees it.
11. **A resume refuses a rewritten spectrum file.** The config names the
   spectrum and the model reads the file, so a config comparison cannot see
   `k25v.dat` regenerated in place. CONS-3.
12b. **SHTns is asked for the one grid mode that runs no timing race.**
   `SHT_GAUSS` and `SHT_GAUSS_FLY` both benchmark algorithm variants at startup
   and keep the winner, which sets the transform functions and so the restart
   hash. Archive CLIM-44 is that, paid for once. CLIM-74 measured the
   alternatives at about 0.1% of runtime, so the tempting edit is all cost.
12c. **pyburn's compiled extensions load on this interpreter.** `pyfft` and
   `pyfft991` are untracked f2py build artifacts tagged with a CPython ABI, so
   rebuilding the venv on a new interpreter invalidates them and every run then
   fails at postprocessing -- reported by ExoPlaSim as the MODEL crashing.
12d. **The transform gates run the configured spectral filter.** The SHTns
   drivers set `filterkappa` and `nfilterexp` as literals while claiming to run
   the model's configuration; `filter_power` moved and they did not.
12e. **The rung-to-dimension table is `lib/rungs.py` and nowhere else.** The
   path lint tests for a rung inside something path-shaped, so a `case "$res"
   in T21) nlat=32` and a `{"T21": (32, 64)}` are out of its scope by
   construction -- and those are what the surviving SPAT-2 copies were. Both
   rung lints read shell scripts as well as Python, which `d.glob("*.py")`
   made impossible.
12f. **No continuation marker has been dropped from the model Fortran.** A
   text check for the class that cost the most: `landmod.f90`'s namelist gained
   a line and the line above it kept its terminator, so the statement ended
   early and `,snowcovz` began a new one. See the check for why a grep is the
   right shape for this and a compile is not a substitute.
12h. **No name is loaded that nothing binds**, in this project's Python AND in
   the Python the model ships. Asked of CPython's own symbol table, so it is
   the language's resolver that answers rather than a pattern. `world-ro6` swept
   the model Fortran for unreferenced procedures and applied the same sweep to
   `pyburn.py`, where a Fortran call graph cannot see a Python caller; the run
   after it integrated ten orbits and died in the postprocessor. Narrower than
   check 1b on purpose -- see the check.
12i. **The spectral tail slope is fitted above the roundoff floor.** A
   synthetic spectrum with a tail of known slope sitting on a floating-point
   floor. `spectral_tail.py` used to cut the spectrum at `m = NTRU` and fit
   through the dead top, so its reported tail slope was a property of the floor;
   the check's right answer is the slope put in.
13. **The tools `environment.md` names are actually on this host.** That
   document sends a reader to `ncdump`, NCO, `h5diff` and `yq` rather than a
   Python session, and nothing else checks the claim is true. Both
   directions, so the document and the check cannot drift apart.
14. **Every unclosed issue carries exactly one `batch:<n>`.** The batch label is
   the tracker's ordering -- what a change to a row makes worthless -- and it is
   maintained by hand, so a newly filed row silently sits outside the order and
   `bd ready --label batch:2` quietly under-reports. A text read of the tracked
   export, which is the copy every other checkout sees.
15. **Every pedology bracket contains the value it brackets.** The declarations
   under `pedology/config/` bracket nearly every number they carry, and nothing
   read the brackets: a value edited outside its own, or a bracket whose ends
   were swapped, reached a soil build with nothing objecting. Four dispositions,
   because a bracket on a LEVEL no key states and a bracket the carbonate
   equilibrium DERIVES are enforced elsewhere, and what is checked for those two
   is that they still are. `pedology/scripts/outgassing_gate.py` runs the same
   job for its own declaration at the artifact tier; every input to this one is
   a config read.
"""

from __future__ import annotations

import os

# ONE CORE, BEFORE ANYTHING IMPORTS NUMPY. Nothing in this file is a numerical
# workload: the two numeric checks fit a plane on a 2,562-vertex icosphere and
# run an AR(1) recursion, and neither is BLAS-bound. numpy's bundled OpenBLAS
# nonetheless opens a pool sized to the logical core count -- 32 here -- and
# OpenBLAS threads SPIN after a parallel region rather than sleeping, so the
# gate was observed at 1129% CPU and kept burning cores through its own
# teardown. That is what made a static gate the thing several sessions were
# waiting on at once, and there is no reading of it in which a coherence check
# should take eleven cores. Set before the import because a pool is sized when
# the library loads and cannot be resized afterwards.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("OMP_WAIT_POLICY", "passive")

import argparse
import ast
import builtins
import importlib.util
import inspect
import math
from pathlib import Path
import re
import shutil
import subprocess
import symtable
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIRS = [ROOT / "scripts", ROOT / "lib", ROOT / "exoplasim" / "scripts",
               ROOT / "hydrography" / "scripts", ROOT / "pedology" / "scripts",
               ROOT / "biosphere" / "scripts", ROOT / "maps"]

# `sorted(x.glob(...))[i]`, `sorted(glob(...))[i]`, `list(...glob(...))[i]`,
# and max/min over a glob: all of them choose one artifact by ordering.
ORDER_PICK = re.compile(
    r"(sorted|list|max|min)\s*\([^\n]*\.?glob\([^\n]*\)[^\n]*\)\s*\[")


def check_modules_parse(files: list[Path]) -> list[str]:
    bad = []
    for f in files:
        try:
            ast.parse(f.read_text(encoding="utf-8"))
        except SyntaxError as exc:
            bad.append(f"{f.relative_to(ROOT)}: {exc}")
    return bad


def check_purge_never_reaches_the_terrain() -> list[str]:
    """`--purge` must never offer to delete the export, from any seed.

    This has a right answer rather than merely a different one, which is what
    makes it a test: loop A is a cycle, so reachability from a climate step comes
    back round through `carve_list` to `orogen` and then, via the export, to
    everything. The cut at the generation step is the only thing stopping a purge
    seeded anywhere in the loop from proposing `source/{build}/`. Assert it from
    EVERY seed rather than the one that motivated it, because the next step added
    inside loop A gets the guarantee for free or not at all.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    import pipeline
    graph = pipeline.load()
    by_id = pipeline.steps_by_id(graph)
    bad = []
    for seed in by_id:
        if "orogen" in pipeline.downstream(seed, by_id):
            bad.append(f"--purge {seed} reaches orogen, so it would offer to "
                       f"delete the export")
        # Asserted over the DELETE set rather than reachability, because those
        # are two sets now: `purge_set` adds the seed's own output back, so the
        # generation step is the one seed that can put `source/{build}/` into a
        # purge plan without any traversal reaching it. world-ysd1.
        for sid, writes in pipeline.purge_writes(seed, graph, by_id):
            for w in writes:
                if w.startswith("source/"):
                    bad.append(f"--purge {seed} would delete {w} (written by "
                               f"{sid}); source/ is read-only")
    # The export edge is what seeds `--purge orogen`; a typo in it silently
    # empties that seed rather than erroring, and the purge would report success
    # having deleted nothing.
    readers = [s["id"] for s in graph["steps"] if s.get("reads_export")]
    if not readers:
        bad.append("no step declares reads_export, so --purge orogen is a no-op")
    return bad


def check_purge_covers_the_seeds_own_writes() -> list[str]:
    """`--purge X` must delete every path X's own `writes` names.

    An identity over the whole graph, and the right answer is known before the
    walk runs: `config/pipeline.yaml` requires `writes` to be the COMPLETE set of
    a step's outputs, and a change to a step invalidates that set first. So the
    plan `--purge X` builds contains X's rows, entry for entry, for every X.

    It fails on the walk this replaced. Reverse reachability alone returned an
    empty set at every step nothing else `needs`, so three registered leaves each
    purged none of their own output and exited zero -- and the completeness rule
    on `writes` was unenforceable exactly there, since naming every artifact
    changed nothing about what got deleted. world-ysd1.

    `orogen` is the one exemption and it is asserted in the opposite direction,
    in the terrain check above: its `writes` is `source/{build}/`, and rule 7
    makes that read-only.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    import pipeline
    graph = pipeline.load()
    by_id = pipeline.steps_by_id(graph)
    bad = []
    for seed, step in by_id.items():
        if seed == "orogen":
            continue
        planned = dict(pipeline.purge_writes(seed, graph, by_id))
        covered = set(planned.get(seed, []))
        for w in step.get("writes", []):
            if w not in covered:
                bad.append(f"--purge {seed} does not delete {w}, which {seed} "
                           f"itself writes")
    return bad


# The one `<component>/data/` path in the graph that is not per-build, exempt BY
# NAME and with a reason: the dust optics table is a property of the mineralogy
# and the stellar spectrum, and no terrain enters it. A pattern would exempt the
# next flat write too, which is the bug this check exists to catch.
BUILD_SCOPED_DATA_EXEMPT = {
    "exoplasim/data/dust": "dust optics are per-mineralogy, not per-terrain",
}


def check_build_scoped_defaults() -> list[str]:
    """Every `<component>/data/` path the graph writes is namespaced by the build.

    `lib/builds.py:component_data` says what a flat `<component>/data/` costs:
    it holds whichever build was active when it was last written, so a script
    reading it pairs one terrain's rows with another's columns, and four scripts
    did. `world_state.py` reported 2,107 basins against a build that had 2,540.

    THE COMPONENT LIST COMES FROM `config/pipeline.yaml`, not from a literal
    here. It was a literal -- `("hydrography", "pedology")` -- and it asserted
    that `component_data(c, cfg).name` equals `cfg["source_build"]`, which is the
    string `component_data` builds the path out of. That comparison cannot fail:
    it was the resolver agreeing with itself, and it left minerals and maps, both
    of which write per-build data, unexamined. world-60x0.

    What can fail now is a real disagreement between two independent
    declarations. The graph says where a step writes; the resolver says where a
    reader looks. A step that writes a flat `<component>/data/x.nc` fails, and so
    does a resolver that lands somewhere other than the directory the graph
    namespaces. Adding a component to the graph puts it under this check with no
    edit here, which is the half the literal list got wrong.

    `strict=True` is deliberately NOT passed. It refuses a path whose directory
    does not exist, and a build whose per-build data has not been generated yet
    is the resting state of the tree rather than a defect; this asks where a
    reader would look, not whether the artifact is there.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    import yaml
    from builds import component_data
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    graph = yaml.safe_load((ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    active = str(cfg.get("source_build", ""))
    if not active:
        return ["config/planet.yaml declares no source_build"]

    bad, components = [], set()
    for step in graph.get("steps", []):
        for w in step.get("writes", []):
            m = re.match(r"^([A-Za-z0-9_]+)/data/(.+)$", str(w))
            if not m:
                continue
            component, rest = m.group(1), m.group(2)
            head = rest.split("/")[0]
            if any(str(w).startswith(k + "/") for k in BUILD_SCOPED_DATA_EXEMPT):
                continue
            if head != "{build}":
                bad.append(f"{step['id']} writes {w}: {component}/data/ is not "
                           f"namespaced by the build, so it holds whichever "
                           f"build was active when it was last written")
                continue
            components.add(component)

    if not components:
        bad.append("no step in config/pipeline.yaml writes a per-build "
                   "<component>/data/{build}/ path, so this check reads nothing")
    for component in sorted(components):
        try:
            d = component_data(component, cfg)
        except Exception as exc:
            bad.append(f"{component}: component_data raised {exc}")
            continue
        want = ROOT / component / "data" / active
        if d != want:
            bad.append(f"{component}: component_data resolves to {d}, but the "
                       f"graph writes that component's data under {want}")
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
    """Every generator that writes an artifact is registered in `config/pipeline.yaml`.

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


def check_registered_paths_exist() -> list[str]:
    """Every script `config/pipeline.yaml` names actually exists.

    The register is only usable if its paths resolve: CLAUDE.md leans on it to
    answer "what is now worthless", and a row pointing at nothing answers that
    question wrongly and silently. This found `exoplasim/scripts/
    extract_high_cadence_wind.py`, registered under `one_offs` while the script
    lived in `aeolian/scripts/` -- a path left behind when the script moved, with
    the real one registered three lines below it.

    All three registers are held to it -- `steps`, `checks` and `one_offs` --
    because a gate row pointing at nothing is the worse case of the two: a step
    that will not run is noticed the next time someone runs it, while a gate
    that will not run reports nothing and reads as a pass.
    """
    import yaml
    graph = yaml.safe_load(
        (ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    named = ([(s["script"], f"step {s['id']}") for s in graph["steps"]]
             + [(c["script"], f"check {c['id']}") for c in graph.get("checks", [])]
             + [(o, "one_offs") for o in graph.get("one_offs", [])])
    return [f"{path} is named by {where} in config/pipeline.yaml and does not exist"
            for path, where in named if not (ROOT / path).exists()]


def check_documented_in_component(files) -> list[str]:
    """Every pipeline step is named in its component's README.

    "An undocumented component is not complete" is a rule in CLAUDE.md and this
    is what makes it hold. A fresh session reads the component README to learn
    what a directory does; a step that works and is invisible there gets
    reimplemented beside itself, which is how this project came to have four
    copies of a path resolver.

    It checks the NAME is present, not that the description is any good, because
    the second is not mechanisable and the first catches the actual failure:
    a script that nobody wrote down at all.

    `checks` rows are entries here for the same reason `one_offs` are. A gate
    is run by a person who went looking for it, so a gate absent from its
    component README is a gate nobody will think to run.
    """
    import yaml
    graph = yaml.safe_load(
        (ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    missing = []
    # `one_offs` are checked too, and that omission is why this was extended.
    # `run_albedo_bracket.sh` is registered there rather than as a step, so the
    # steps-only version passed while the script that re-derives the flux on new
    # terrain was named in no README at all. A tool being run occasionally is a
    # reason to write it down, not a reason not to.
    entries = ([(Path(s["script"]), "a pipeline step") for s in graph["steps"]]
               + [(Path(c["script"]), "registered under checks")
                  for c in graph.get("checks", [])]
               + [(Path(o), "registered under one_offs")
                  for o in graph.get("one_offs", [])])
    for script, what in entries:
        readme = ROOT / script.parts[0] / "README.md"
        if not readme.is_file():
            continue                     # no component README to be absent from
        if script.name not in readme.read_text(encoding="utf-8"):
            missing.append(f"{script} is {what} and is not named in "
                           f"{script.parts[0]}/README.md")
    return missing





def check_uncovered_trailing_orbits_refuse() -> list[str]:
    """Orbits past the manifest's last segment refuse rather than defaulting.

    THE LIVE CASE. A continuation killed partway leaves its orbits on disk and
    no segment saying what they were, because the record is appended when the
    block finishes. `low_io_orbits` then reads them as low-I/O -- correct when
    "unlabelled" means a run older than segments, wrong when it means
    UNRECORDED -- so a block of clean orbits becomes invisible and the
    I/O-regime guard has nothing to refuse. The guard walked around rather than
    failing.

    Class 17, and the pairing is the point: a LEADING gap must still pass,
    because orbit 0 is uncovered on older runs and every orbit is uncovered on
    runs that predate segments entirely.
    """
    import json as _json
    import tempfile
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    from segments import production_window, refuse_orbits_no_segment_covers

    bad = []

    def run_dir(tmp, segments):
        d = Path(tmp)
        if segments is not None:
            (d / "run_manifest.json").write_text(
                _json.dumps({"segments": segments}), encoding="utf-8")
        return d

    def seg(a, b, low_io=True):
        return {"start_year_index": a, "end_year_index": b,
                "purpose": "spinup", "low_io": low_io}

    with tempfile.TemporaryDirectory() as tmp:
        # The killed continuation: 110 orbits, segments stop at 69.
        d = run_dir(tmp, [seg(0, 69)])
        try:
            refuse_orbits_no_segment_covers(d, 110)
            bad.append("40 uncovered trailing orbits did not refuse")
        except RuntimeError:
            pass
        # And it refuses through the door callers actually use.
        try:
            production_window(d, 110, 42)
            bad.append("production_window accepted a run with a trailing gap")
        except RuntimeError as exc:
            if "recorded nowhere" not in str(exc):
                bad.append(f"production_window refused for the wrong reason: {exc}")
        # The paired positive: the same manifest with the orbits it claims.
        try:
            refuse_orbits_no_segment_covers(d, 70)
        except RuntimeError as exc:
            bad.append(f"a fully covered run was refused: {exc}")

    with tempfile.TemporaryDirectory() as tmp:
        # A LEADING gap is legitimate and must pass: orbit 0 uncovered.
        d = run_dir(tmp, [seg(1, 79)])
        try:
            refuse_orbits_no_segment_covers(d, 80)
        except RuntimeError as exc:
            bad.append(f"a leading gap was refused: {exc}")

    with tempfile.TemporaryDirectory() as tmp:
        # No manifest at all: every run made before segments existed.
        d = run_dir(tmp, None)
        try:
            refuse_orbits_no_segment_covers(d, 80)
        except RuntimeError as exc:
            bad.append(f"a run with no manifest was refused: {exc}")
    return bad


def check_io_regime_window() -> list[str]:
    """A verdict window is refused when it spans a change of I/O regime.

    THE LIVE FAILURE. This project's bootstrap ran 70 orbits of low-I/O spinup
    and then a clean-I/O tail. Fitted across the join, the integrated
    autocorrelation time of the per-orbit mean surface temperature came back at
    9.08 orbits; fitted on each side it is 2.00 and 1.00. The whole difference
    was a +0.1616 K step at the boundary, read as memory the planet does not
    have -- and tau is what sizes the window and prices the production span.

    Class 17: every case has a right answer, and every refusal is paired with
    the positive that proves the setup was real rather than the manifest merely
    existing. Synthetic manifests in a temp directory; it touches no run.
    """
    import json as _json
    import tempfile
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    from segments import io_regime_changes, production_window

    bad = []

    def run_dir(tmp: str, segments) -> Path:
        d = Path(tmp)
        (d / "run_manifest.json").write_text(
            _json.dumps({"segments": segments}), encoding="utf-8")
        return d

    def seg(a, b, low_io):
        return {"start_year_index": a, "end_year_index": b,
                "purpose": "spinup", "low_io": low_io}

    def case(name, got, want):
        if got != want:
            bad.append(f"{name}: got {got!r}, expected {want!r}")

    def raises(name, fn) -> None:
        try:
            fn()
        except RuntimeError:
            return
        bad.append(f"{name}: expected a refusal and got a window")

    with tempfile.TemporaryDirectory() as tmp:
        d = run_dir(tmp, [seg(0, 79, False)])
        # The positive that proves the setup: one clean regime throughout and
        # the window is the plain tail, so the refusals below are caused by the
        # join and not by `low_io` being declared at all.
        case("one regime throughout takes the tail",
             production_window(d, 80, 10), (70, 79))
        case("one regime has no change", io_regime_changes(d, range(80)), [])

    with tempfile.TemporaryDirectory() as tmp:
        # The bootstrap's own shape: a low-I/O approach, then a clean tail.
        d = run_dir(tmp, [seg(0, 69, True), seg(70, 81, False)])
        case("the join is found", io_regime_changes(d, range(82)), [70])
        raises("a window spanning the join is refused",
               lambda: production_window(d, 82, 21))
        # Paired positive: the same run judged inside the clean block is fine,
        # so the refusal is about the join and not about the run's length.
        case("a window inside the clean block is accepted",
             production_window(d, 82, 12), (70, 81))
        # And the boundary is exactly where it should be: one orbit more
        # reaches back into the low-I/O block.
        raises("one orbit past the join is refused",
               lambda: production_window(d, 82, 13))

    with tempfile.TemporaryDirectory() as tmp:
        # A run that switched BACK. Both joins are reported, and the window is
        # judged on the first one it meets rather than on the segment list.
        d = run_dir(tmp, [seg(0, 9, True), seg(10, 19, False),
                          seg(20, 29, True)])
        case("both joins are found",
             io_regime_changes(d, range(30)), [10, 20])
        case("a subrange reports only its own joins",
             io_regime_changes(d, range(20, 30)), [])

    with tempfile.TemporaryDirectory() as tmp:
        # A segment with no `low_io` key is low-I/O, which is what every run
        # made before 2026-08-17 was. So an undeclared block beside a declared
        # clean one IS a join, and the unsafe default is what makes it one.
        d = run_dir(tmp, [{"start_year_index": 0, "end_year_index": 9,
                           "purpose": "spinup"}, seg(10, 19, False)])
        case("an undeclared block joins a clean one",
             io_regime_changes(d, range(20)), [10])
    return bad


def check_io_step_measurement() -> list[str]:
    """`assess_convergence.io_step_at_join` measures the step it names.

    The relaxation fit spans the approach by construction, so it may cross a
    join the verdict window may not. What makes that admissible is measuring
    the step against the criterion the asymptote decides, so this holds the
    measurement to a series whose step is known by construction.
    """
    import json as _json
    import tempfile
    import numpy as _np
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    from assess_convergence import io_step_at_join

    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        (d / "run_manifest.json").write_text(_json.dumps({"segments": [
            {"start_year_index": 0, "end_year_index": 19, "purpose": "spinup",
             "low_io": True},
            {"start_year_index": 20, "end_year_index": 29, "purpose": "spinup",
             "low_io": False}]}), encoding="utf-8")
        # Flat at 288.0, then flat at 288.5: the step is 0.5 K exactly and no
        # trend can be confused for it.
        series = _np.array([288.0] * 20 + [288.5] * 10)
        join, step = io_step_at_join(series, d, 0, 29)
        if join != 20:
            bad.append(f"join: got {join}, expected 20")
        if abs(step - 0.5) > 1e-12:
            bad.append(f"step: got {step}, expected 0.5")
        # A range that stops short of the join has no join and no step, which
        # is the case that must not report a spurious one.
        join, step = io_step_at_join(series, d, 0, 19)
        if (join, step) != (None, 0.0):
            bad.append(f"no join: got {(join, step)}, expected (None, 0.0)")
    return bad





def check_resume_adopts_the_manifest_thread_count() -> list[str]:
    """A continuation takes `ncpus` from the run, not from the config file.

    `run_exoplasim.py --ncpus N` writes N into the loaded config before the
    manifest is stamped, so the manifest carries N and the config FILE still
    declares whatever it declares. Comparing them refused a run that should
    have been extended, with no flag able to say otherwise, and cost a re-run
    of a paired soil-thermal experiment. world-q4gh.

    ADOPTION AND NOT AN EXEMPTION: a continuation integrates on the executable
    the run directory already holds, and ExoPlaSim compiles one per
    (resolution, layers, threads), so a config declaring 16 against a manifest
    of 8 must resolve the p8 binary. Dropping the key from the comparison
    instead would resolve a p16 binary into a run whose earlier orbits are p8.

    The check reads the source rather than driving `main`, because driving it
    needs a run directory, a binary and a restart. Class 17 is satisfied by the
    third case: the ADOPTION direction is what is asserted, not merely that the
    key is mentioned.
    """
    src = (ROOT / "exoplasim" / "scripts" / "continue_exoplasim.py").read_text(
        encoding="utf-8")
    bad = []
    if "manifest_ncpus" not in src:
        bad.append("continue_exoplasim.py does not read the manifest's ncpus, "
                   "so a run prepared with --ncpus cannot be resumed at all")
        return bad
    adopt = src.index("manifest_ncpus")
    drift = src.index("drift = config_drift(")
    if adopt > drift:
        bad.append("the thread count is adopted AFTER the drift comparison, "
                   "so the comparison still refuses the run it should extend")
    # The direction: the manifest's value goes INTO the config, never the
    # other way, or the binary resolver reads a count the run does not have.
    if 'config["model"]["ncpus"] = manifest_ncpus' not in src:
        bad.append("the manifest's thread count is not written into the "
                   "config, so the binary resolver can still read the file's")
    if "INERT_CONFIG_KEYS" in src and '"ncpus"' in src.split(
            "INERT_CONFIG_KEYS")[1][:400]:
        bad.append("ncpus looks like it was added to INERT_CONFIG_KEYS; that "
                   "waives the check instead of adopting the value, and lets a "
                   "p16 binary resolve into a run whose orbits are p8")
    return bad


def check_production_span_is_self_limiting() -> list[str]:
    """The span comes from the run's own criteria, not from twenty tau.

    THE MEASUREMENT THIS REPLACES A RULE WITH. Tau on this model is a property
    of the WINDOW rather than of the process -- 1.00 at twenty orbits, 6.30 at
    a hundred and forty on one clean block, every reading supported -- so a
    span of twenty tau never closes: each span bought raises the tau that
    prices the next. The standard error of the mean it exists to correct has
    already converged, 0.0201 K at twenty orbits against 0.0255 at a hundred
    and forty. `exoplasim/notes/memory-time-and-the-production-span.md`.

    Class 17: each case has a right answer and each is driven by a synthetic
    report, so it tests the rule and not today's runs.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    import run_lengths as rl

    bad = []

    def report(storage, offset, floor):
        return {"resolving_power": {
            "window_orbits_for_storage_criterion": storage,
            "window_orbits_for_offset_criterion": offset,
            # The same numbers with the memory taken out, which must NOT be
            # read: a rule that swept the block by prefix would take them and
            # under-buy by exactly the memory correction.
            "window_orbits_for_storage_criterion_if_independent": 1.0,
            "window_orbits_for_offset_criterion_if_independent": 1.0,
            "window_orbits_for_offset_criterion_prices": "prose, not a number",
            "required_window_is_a_lower_bound": floor}}

    got = rl.production_span_from_report(report(10.0, 30.0, False))
    if got != (30.0, False, "window_orbits_for_offset_criterion"):
        bad.append(f"the binding criterion: got {got!r}")

    # THE OFFSET ROW'S QUALIFIER. Its own `..._prices` field says it prices the
    # verdict only where `offset_statistic_source` is the drift fallback; on the
    # exponential-fit path the statistic is the FIT's half width and this number
    # is computed from the expected relaxation instead. Taking it anyway read
    # 42.8 orbits on a run whose statistic was six times its target.
    fit = report(10.0, 30.0, False)
    fit["resolving_power"]["offset_statistic_source"] = "exponential_fit"
    got = rl.production_span_from_report(fit)
    if got[2] != "window_orbits_for_storage_criterion":
        bad.append(f"the offset row was used on the fit path: got {got!r}")
    if got[1] is not True:
        bad.append("dropping a criterion did not make the answer a floor, so a "
                   "caller would take the remaining criteria as sufficient")
    # And the drift fallback keeps it, or the qualifier would just disable the
    # row rather than reading it.
    drift = report(10.0, 30.0, False)
    drift["resolving_power"]["offset_statistic_source"] = "drift_fallback"
    if rl.production_span_from_report(drift)[2] != "window_orbits_for_offset_criterion":
        bad.append("the offset row was dropped on the drift-fallback path too")
    got = rl.production_span_from_report(report(40.0, 30.0, False))
    if got[0] != 40.0 or got[2] != "window_orbits_for_storage_criterion":
        bad.append(f"the other criterion binding: got {got!r}")
    if rl.production_span_from_report(report(10.0, 30.0, True))[1] is not True:
        bad.append("a lower-bound span was not reported as a floor")

    # A report predating the block cannot answer, and must say so rather than
    # returning a number: buying too few orbits is the one direction this must
    # not fail silently in.
    try:
        rl.production_span_from_report({"resolving_power": {}})
        bad.append("a report with no required-window keys returned a span")
    except RuntimeError:
        pass

    total, floor, _ = rl.commissioning_orbits_from_report(70.0, report(10.0, 30.0, False))
    if total != 100.0 or floor is not False:
        bad.append(f"the approach is not added: got {(total, floor)!r}")

    # The a-priori rule still exists and still means twenty tau, so the two can
    # be compared and the note's claim that it overbuys stays checkable.
    if rl.production_span_orbits(6.3) != 20.0 * 6.3:
        bad.append("production_span_orbits is no longer twenty tau")
    return bad


def check_no_write_through_a_symlink() -> list[str]:
    """`sra.py:write_sra` refuses to write a staged field through a link.

    THE LIVE HAZARD. `link_worktree.py` links `exoplasim/inputs/<rung>/` entry
    by entry so a worktree can READ the staged fields to prepare a run, and the
    four builders that write that directory resolve their default output to it.
    So a generator run in a worktree writes STRAIGHT THROUGH into the main
    checkout and replaces the fields a commissioning run there is reading. One
    agent queued exactly that rebuild and cancelled it before it took the lock.

    Class 17: every case has a right answer, and the refusals are paired with
    the positive proving an ordinary write still works. Temp directories only;
    it touches no staged field.
    """
    import tempfile
    import numpy as _np
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    from sra import write_sra

    bad = []
    field = _np.zeros((4, 8))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        donor = tmp / "main"
        donor.mkdir()
        (donor / "f.sra").write_text("donor", encoding="ascii")

        try:
            write_sra(tmp / "plain.sra", 174, field)
        except SystemExit as exc:
            bad.append(f"an ordinary path was refused: {exc}")

        link = tmp / "linked.sra"
        link.symlink_to(donor / "f.sra")
        try:
            write_sra(link, 174, field)
            bad.append("a symlinked file was written through")
        except SystemExit:
            pass

        linked_dir = tmp / "inputs"
        linked_dir.symlink_to(donor, target_is_directory=True)
        try:
            write_sra(linked_dir / "g.sra", 174, field)
            bad.append("a path under a symlinked directory was written through")
        except SystemExit:
            pass

        # The point of the guard: the donor is untouched by either refusal.
        if (donor / "f.sra").read_text(encoding="ascii") != "donor":
            bad.append("the donor file was modified despite the refusal")
        if (donor / "g.sra").exists():
            bad.append("a file was created in the donor directory")
    return bad


def check_production_window() -> list[str]:
    """`segments.py:production_window` picks the window a declaration implies.

    Class 17: every case has a right answer that is not a matter of taste, and
    every negative is paired with the positive proving the setup was real. The
    live failure it reproduces is the third one -- a three-orbit diagnostic tail
    on a differently patched binary, which the trailing-window rule would have
    averaged into a convergence verdict.

    Synthetic manifests in a temp directory; it touches no run.
    """
    import json as _json
    import tempfile
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    from segments import production_window

    bad = []

    def run_dir(tmp: str, segments) -> Path:
        d = Path(tmp)
        if segments is not None:
            (d / "run_manifest.json").write_text(
                _json.dumps({"segments": segments}), encoding="utf-8")
        return d

    def seg(a, b, purpose):
        return {"start_year_index": a, "end_year_index": b, "purpose": purpose}

    def case(name, got, want):
        if got != want:
            bad.append(f"{name}: got {got!r}, expected {want!r}")

    def raises(name, fn) -> None:
        try:
            fn()
        except RuntimeError:
            return
        bad.append(f"{name}: expected a refusal and got a window")

    with tempfile.TemporaryDirectory() as tmp:
        # No manifest at all: every orbit is production, so the window is the
        # plain tail. This is the pre-existing behaviour and every run made
        # before purposes were declared depends on it.
        case("no manifest takes the tail",
             production_window(run_dir(tmp, None), 80, 10), (70, 79))

    with tempfile.TemporaryDirectory() as tmp:
        # 1 to 79 and not 1 to 76: the segment has to COVER the orbits, or
        # `refuse_orbits_no_segment_covers` refuses the run before the window
        # is computed. The three-orbit trailing gap this used to carry was
        # incidental to what the case asserts and the case is unchanged
        # without it -- but a trailing gap is now a defect in its own right,
        # because it is what a killed continuation leaves and it hides an
        # I/O regime nobody recorded.
        d = run_dir(tmp, [seg(1, 79, "spinup")])
        # The positive that proves the setup: with the tail declared spinup the
        # answer is still the tail, so the next case's shift is caused by the
        # declaration and not by the manifest merely existing.
        case("a spinup tail takes the tail",
             production_window(d, 80, 10), (70, 79))

    with tempfile.TemporaryDirectory() as tmp:
        d = run_dir(tmp, [seg(1, 76, "spinup"), seg(77, 79, "diagnostic")])
        case("a diagnostic tail is dropped",
             production_window(d, 80, 10), (67, 76))

    with tempfile.TemporaryDirectory() as tmp:
        d = run_dir(tmp, [seg(1, 71, "spinup"), seg(72, 74, "diagnostic"),
                          seg(75, 79, "spinup")])
        raises("a diagnostic hole inside the window is refused",
               lambda: production_window(d, 80, 10))
        # Paired positive: the same run assessed over a window that clears the
        # hole is fine, so the refusal is about the hole and not the manifest.
        case("a window clearing the hole is accepted",
             production_window(d, 80, 5), (75, 79))

    with tempfile.TemporaryDirectory() as tmp:
        d = run_dir(tmp, [seg(0, 9, "diagnostic")])
        raises("a run with no production orbits is refused",
               lambda: production_window(d, 10, 5))

    with tempfile.TemporaryDirectory() as tmp:
        d = run_dir(tmp, [seg(8, 9, "diagnostic")])
        raises("too few production orbits for the window is refused",
               lambda: production_window(d, 10, 10))
    return bad


def check_spectrum_guard() -> list[str]:
    """A resume refuses a spectrum file rewritten under the same name.

    CONS-3. The config names the spectrum and the model reads the FILE, so a
    config comparison cannot see `k25v.dat` regenerated in place. Every negative
    here is paired with the positive proving the manifest was otherwise
    acceptable, or the refusal would prove nothing.
    """
    import copy
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    import yaml
    from run_exoplasim import stellar_spectrum_digest, require_stellar_spectrum

    bad = []
    config = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    digest = stellar_spectrum_digest(config)
    if not digest or not digest.get("sha256") or not digest.get("hr_sha256"):
        return ["stellar_spectrum_digest returned no digest for the configured "
                "spectrum, so nothing below tests anything"]

    manifest = {"stellar_spectrum_digest": copy.deepcopy(digest)}
    if require_stellar_spectrum(manifest, config):
        bad.append("an unchanged spectrum was reported as a backfill")

    for key in ("sha256", "hr_sha256"):
        tampered = {"stellar_spectrum_digest": copy.deepcopy(digest)}
        tampered["stellar_spectrum_digest"][key] = "0" * 64
        try:
            require_stellar_spectrum(tampered, config)
            bad.append(f"a changed {key} was resumed across")
        except RuntimeError:
            pass

    # The backfill path prints its warning; swallow it so the check's own
    # output is the only thing on stdout.
    import contextlib
    import io
    unstamped: dict = {}
    with contextlib.redirect_stdout(io.StringIO()):
        backfilled = require_stellar_spectrum(unstamped, config)
    if not backfilled:
        bad.append("an unstamped run did not report that it was backfilled")
    if unstamped.get("stellar_spectrum_digest") != digest:
        bad.append("the backfill did not land in the manifest")
    return bad


def check_input_stamp_guard() -> list[str]:
    """A re-derived input file is detected under an unchanged name.

    world-nvs2, and the same shape as `check_spectrum_guard` above: the config
    stamp compares parsed configuration values and is blind to a DERIVED FILE a
    generator read, because that file's name does not move when its own
    generator re-runs. Every negative here is paired with the positive proving
    the record was otherwise acceptable, or the refusal would prove nothing.
    """
    import json
    import tempfile
    from provenance import artifact_input_drift, input_drift, input_stamp

    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        read = root / "derived.json"
        read.write_text('{"shape": 1}\n', encoding="utf-8")
        absent = root / "never_written.json"

        # The stamp keys a path relative to the PROJECT root, so a path outside
        # it -- which this temporary directory is -- is keyed absolute. Read the
        # keys back rather than assume either shape; `input_drift` joins them
        # onto its `root` and an absolute key survives that join unchanged.
        stamp = input_stamp([read, absent])["source_inputs"]
        if stamp.get(str(read)) is None:
            bad.append("a file that exists was stamped as absent")
        if str(absent) not in stamp:
            bad.append("an absent input was dropped from the stamp rather "
                       "than recorded as absent")

        if input_drift(stamp, root):
            bad.append("an unchanged input was reported as moved")

        read.write_text('{"shape": 2}\n', encoding="utf-8")
        if not input_drift(stamp, root):
            bad.append("a rewritten input was not reported")

        absent.write_text("{}\n", encoding="utf-8")
        if len(input_drift(stamp, root)) != 2:
            bad.append("an input that appeared after the build was not reported")

        read.unlink()
        if not any("missing" in line for line in input_drift(stamp, root)):
            bad.append("a deleted input was not reported as missing")

        # And the artifact-facing wrapper's three-valued answer: None for a
        # record with no stamp at all, which is UNOBSERVABLE and not "current".
        unstamped = root / "unstamped_report.json"
        unstamped.write_text(json.dumps({"generator": "x"}) + "\n",
                             encoding="utf-8")
        if artifact_input_drift(unstamped, root) is not None:
            bad.append("an artifact with no input stamp was reported as checked")
        stamped = root / "stamped_report.json"
        stamped.write_text(json.dumps({"source_inputs": stamp}) + "\n",
                           encoding="utf-8")
        if not artifact_input_drift(stamped, root):
            bad.append("artifact_input_drift did not see the drift its own "
                       "input_drift reports")
    return bad


def check_donor_surface_guard() -> list[str]:
    """`--restart-from` refuses a frozen surface and admits a reread one.

    THE TWO ERRORS THIS SITS BETWEEN. Too loose and a run seeded from a donor
    silently discards a restaged `dwmax` or `dalbcl` and reproduces its parent,
    because `landini` takes those from the restart when `nrestart > 0`. Too
    strict and it refuses over a code the model REREADS from the `.sra` at every
    start, which a restart cannot freeze in either direction -- and that costs a
    cold start for nothing. PHYS-15's paired arm is the case that found the
    second: `nwetsoil = 0` drops codes 1742, 1750 and 1760 from the staged set,
    so the arm refused on "codes differ" against a donor that staged them, over
    three files it never opens and the donor never carried.

    Driven against the REAL config and the real staged files, with synthesised
    donor manifests, so the code sets are the ones a run would actually stage
    rather than a fixture's idea of them. Every negative is paired with the
    positive proving the manifest was otherwise acceptable.
    """
    import copy
    import yaml
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    import run_exoplasim as rx
    import restart_surface

    config = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    codes = rx.intended_surface_codes(config)
    reread = set(rx.REREAD_SURFACE_CODES)
    frozen = sorted(codes - reread)
    staged_reread = sorted(codes & reread)
    if not frozen:
        return ["this config stages no restart-frozen surface code, so the "
                "refusing half of the guard has no sample"]

    def manifest(from_file, hashes=None):
        present = [c for c in from_file if rx.surface_sra(config, c).is_file()]
        h = {str(c): rx.file_sha256(rx.surface_sra(config, c)) for c in present}
        h.update(hashes or {})
        return {"surface_fields": {"from_file": sorted(from_file)},
                "surface_field_sha256": h}

    bad = []

    # POSITIVE. A donor staging exactly what this run stages is acceptable.
    why = rx.donor_surface_reason(manifest(codes), config)
    if why is not None:
        bad.append(f"a matching donor was refused: {why}")

    # NEGATIVE. A donor missing a restart-FROZEN code is refused, because the
    # run would inherit the donor's field and discard the staged one.
    if rx.donor_surface_reason(manifest(set(codes) - {frozen[0]}), config) is None:
        bad.append(f"a donor missing frozen code {frozen[0]} was accepted")

    # NEGATIVE. A frozen code whose CONTENT moved is refused.
    if rx.donor_surface_reason(
            manifest(codes, {str(frozen[0]): "0" * 64}), config) is None:
        bad.append(f"a donor whose frozen code {frozen[0]} changed content "
                   "was accepted")

    # NEGATIVE. No hashes at all: content changes cannot be ruled out.
    if rx.donor_surface_reason(
            {"surface_fields": {"from_file": sorted(codes)}}, config) is None:
        bad.append("a donor recording no surface hashes was accepted")

    # THE REPAIR, BOTH DIRECTIONS. A reread code is transparent: the model looks
    # for the file at every start rather than for a restart record, so neither
    # its presence nor its content can be frozen into one.
    if staged_reread:
        code = staged_reread[0]
        why = rx.donor_surface_reason(manifest(set(codes) - {code}), config)
        if why is not None:
            bad.append(f"a donor that did not stage reread code {code} was "
                       f"refused: {why}")
        why = rx.donor_surface_reason(manifest(codes, {str(code): "0" * 64}),
                                      config)
        if why is not None:
            bad.append(f"a donor whose reread code {code} changed content was "
                       f"refused: {why}")
        # And the arm that found it: the same donor against a config that turns
        # the moisture term off, which drops all three codes from the set.
        off = copy.deepcopy(config)
        off["surface"]["soil_albedo_moisture"]["enabled"] = False
        dropped = codes - rx.intended_surface_codes(off)
        if dropped != set(rx.WET_ALBEDO_SURFACE_CODES):
            bad.append(f"turning the moisture term off dropped {sorted(dropped)}, "
                       "not the saturated albedo pair")
        why = rx.donor_surface_reason(manifest(codes), off)
        if why is not None:
            bad.append("the nwetsoil = 0 arm was refused against a donor that "
                       f"staged the pair: {why}")

    # THE CLAIM THE EXEMPTION RESTS ON, checked where it is made rather than
    # trusted: every reread code must be one `restart_surface` registers as
    # never reaching a restart.
    stray = sorted(reread - set(restart_surface.REREAD_EVERY_START))
    if stray:
        bad.append(f"codes {stray} are exempted from the donor guard without "
                   "restart_surface.REREAD_EVERY_START saying the model rereads "
                   "them")
    return bad


def check_staged_surface_build_guard() -> list[str]:
    """A staged surface field from another build is refused at the read.

    world-xgtj. `exoplasim/inputs/<rung>/orogen_<RUNG>_surf_<code>.sra` is keyed
    by the RUNG alone while `surface_albedo` rewrites it per BUILD, so the path
    cannot say which build's field is in it. Every negative below is paired with
    the positive proving the fixture was otherwise acceptable, and the fixture
    is built from the REGISTRY's own hashes so it cannot pass by describing a
    build that does not exist.
    """
    import hashlib
    import json
    import tempfile
    import yaml
    from orogen import _KNOWN_TERRAIN_HASHES
    from provenance import staged_surface_field

    config = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    rung = str(config["model"]["resolution"]).upper()
    import builds as _builds
    active = _builds.terrain_hash(config)
    other = next((h for h in _KNOWN_TERRAIN_HASHES if h != active), None)
    if other is None:
        return ["the registry holds only one terrain hash, so a cross-build "
                "read cannot be constructed and nothing below tests anything"]

    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        rung_dir = root / "exoplasim" / "inputs" / rung.lower()
        rung_dir.mkdir(parents=True)
        sra = rung_dir / f"orogen_{rung}_surf_0174.sra"
        sra.write_text("self-test background albedo\n", encoding="ascii")
        report = rung_dir / "albedo_report.json"

        def stamp(terrain, codes=(174, 175, 176, 212)):
            report.write_text(json.dumps(
                {"terrain_hash": terrain, "codes": list(codes)}) + "\n",
                encoding="utf-8")

        stamp(active)
        try:
            rec = staged_surface_field(174, config, root=root)
            if rec["terrain_hash"] != active:
                bad.append("the record does not carry the build it read")
            if rec["sha256"] != hashlib.sha256(sra.read_bytes()).hexdigest():
                bad.append("the record's sha256 is not the file's")
        except SystemExit as exc:
            bad.append(f"a field staged from the active build was refused: {exc}")

        stamp(other)
        try:
            staged_surface_field(174, config, root=root)
            bad.append("a field staged from another build was read without a "
                       "declaration")
        except SystemExit:
            pass

        # A declaration NAMES the build, so it admits exactly that one.
        other_name = _KNOWN_TERRAIN_HASHES[other]["name"]
        try:
            rec = staged_surface_field(174, config, root=root,
                                       for_build=other_name)
            if rec["declared_cross_build"] != other_name:
                bad.append("a declared cross-build read did not record the "
                           "declaration")
        except SystemExit as exc:
            bad.append(f"a correctly declared cross-build read was refused: {exc}")
        stamp(active)
        try:
            staged_surface_field(174, config, root=root, for_build=other_name)
            bad.append("a declaration for one build admitted another build's field")
        except SystemExit:
            pass
        try:
            staged_surface_field(174, config, root=root, for_build="no-such-build")
            bad.append("an unregistered build name was accepted as a declaration")
        except SystemExit:
            pass

        # An unstamped field is UNOBSERVABLE, not current.
        report.unlink()
        try:
            staged_surface_field(174, config, root=root)
            bad.append("a staged field with no provenance beside it was read")
        except SystemExit:
            pass
        # And a code no record names is unstamped even when the directory has one.
        stamp(active, codes=(173,))
        try:
            staged_surface_field(174, config, root=root)
            bad.append("a code no record names was read off another code's stamp")
        except SystemExit:
            pass
    return bad


def check_slope_fit() -> list[str]:
    """`lib/orogen.py:local_slope_deg` reproduces a gradient it can get wrong.

    Two cases on a synthetic icosphere, because the real mesh costs 1.7 GB and
    ten seconds and this has to be cheap enough to run every time:

    * a CONSTANT surface must give exactly zero, which catches a sign or an
      index error that a smooth field would hide;
    * a linear ramp `f = a*z` must reproduce its analytic surface gradient
      `a*sqrt(1-z^2)/R` to 2%, which is what catches the estimator being biased.

    THE SECOND CAN FAIL, and it is worth saying what would make it. The first
    implementation of this field took the steepest drop to a neighbour, and that
    estimator scores 2.7% at the median and 25.8% at worst on this same sphere,
    because a gradient rarely points exactly at a neighbour. The plane fit scores
    0.08% and 1.7%. So the tolerance sits between two implementations that both
    look reasonable, which is the only kind of threshold worth writing down.

    An icosphere rather than a Fibonacci sphere, and that was not free either: on
    a Fibonacci sphere with nearest-neighbour adjacency, 5.6% of vertices exceed
    2% and all of them sit near the spiral's poles, where the point pattern and
    therefore the neighbour geometry degrade. That is a defect of the test mesh
    and it would have been read as a defect of the fit.

    The expectation is stated in the mesh's OWN cartesian frame on purpose.
    Writing it in latitude is how this was first got wrong -- the export is
    y-up, so `z` is not `sin(lat)`, and a correct estimator was blamed for a 30%
    error that lived entirely in the frame. See `source/README.md`.
    """
    import numpy as np
    sys.path.insert(0, str(ROOT / "lib"))
    from orogen import Export

    golden = (1 + 5 ** 0.5) / 2
    verts = [np.array(v, dtype=float) for v in
             [(-1, golden, 0), (1, golden, 0), (-1, -golden, 0), (1, -golden, 0),
              (0, -1, golden), (0, 1, golden), (0, -1, -golden), (0, 1, -golden),
              (golden, 0, -1), (golden, 0, 1), (-golden, 0, -1), (-golden, 0, 1)]]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
             (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
             (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(4):
        cache: dict = {}
        split = []

        def midpoint(a: int, b: int) -> int:
            key = (min(a, b), max(a, b))
            if key not in cache:
                verts.append((verts[a] + verts[b]) / 2)
                cache[key] = len(verts) - 1
            return cache[key]

        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            split += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = split

    pos = np.array(verts)
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)
    adjacent = [set() for _ in pos]
    for a, b, c in faces:
        for i, j in ((a, b), (b, c), (c, a)):
            adjacent[i].add(j)
            adjacent[j].add(i)
    off = np.concatenate([[0], np.cumsum([len(a) for a in adjacent])]).astype(np.int64)
    lst = np.concatenate([sorted(a) for a in adjacent]).astype(np.int64)

    n, radius, amplitude = len(pos), 7000.0, 100.0

    class Stub:
        n_regions = n
        adjacency = (off, lst)
        manifest = {"planet": {"radiusKm": radius}}

        def __init__(self, elevation):
            self._elevation = elevation

        def field(self, name):
            return {"x": pos[:, 0], "y": pos[:, 1], "z": pos[:, 2],
                    "elevation_km": self._elevation}[name]

    fit = Export.local_slope_deg.func
    problems = []

    flat = float(np.abs(fit(Stub(np.full(n, 0.7)))).max())
    if flat != 0.0:
        problems.append(f"a constant surface has slope {flat:.3e} deg, "
                        f"must be exactly 0")

    ramp = fit(Stub(amplitude * pos[:, 2]))
    expect = np.degrees(np.arctan(
        amplitude * np.sqrt(np.clip(1 - pos[:, 2] ** 2, 0, 1)) / radius))
    live = expect > 0.1                       # away from the ramp's own poles
    worst = float(np.abs(ramp[live] / expect[live] - 1).max())
    if worst > 0.02:
        problems.append(f"a linear ramp is reproduced to only {worst:.1%}; a "
                        f"steepest-drop estimator scores 25.8% here and the "
                        f"plane fit 1.7%, so this is the difference between them")
    return problems


def check_shtns_init_is_deterministic() -> list[str]:
    """The model asks SHTns for the ONE grid mode that runs no timing race.

    CLIM-74 measured what the alternatives cost and the answer was determinism.
    `SHT_GAUSS` benchmarks its algorithm variants at startup and keeps the
    winner, and `SHT_GAUSS_FLY` is not the deterministic floor it reads as:
    `shtns_set_grid_auto` sets `quick_init` only for `quick_init`, `reg_fast`
    and `reg_poles`, so `gauss_fly` falls through to `choose_best_sht` exactly
    as `gauss` does. Probed at T170, both return a different pick vector on
    every run; `SHT_QUICK_INIT` runs zero races.

    A different transform function is a different rounding order, so the model's
    restart hash follows the winner of that race. Archive CLIM-44 is this model
    failing to be reproducible and it was paid for once. The choice is therefore
    not a tuning preference, and this check is here because the tempting edit --
    swapping in `SHT_GAUSS` or `SHT_AUTO` to let the library pick something
    faster -- reintroduces it silently and costs about 0.1% of runtime at best.
    `exoplasim/notes/shtns-algorithm-selection.md` has the measurements.
    """
    src = ROOT / "vendor/exoplasim/exoplasim/plasim/src/shtnsmod.f90"
    problems = []
    if not src.is_file():
        return problems
    racing = ("SHT_GAUSS_FLY", "SHT_AUTO", "SHT_REG_DCT", "SHT_REG_FAST",
              "SHT_REG_POLES", "SHT_GAUSS")
    seen = False
    for n, line in enumerate(src.read_text(errors="replace").splitlines(), 1):
        code = line.split("!")[0]
        if "klay" not in code or "=" not in code:
            continue
        seen = True
        for name in racing:
            # SHT_GAUSS is a prefix of SHT_GAUSS_FLY, so match the token.
            if re.search(rf"\b{name}\b", code):
                problems.append(
                    f"shtnsmod.f90:{n} asks SHTns for {name}, which runs a "
                    f"startup timing race; its winner sets the transform "
                    f"functions and so the restart hash. Archive CLIM-44.")
        if not re.search(r"\bSHT_QUICK_INIT\b", code):
            problems.append(
                f"shtnsmod.f90:{n} sets the SHTns grid mode without "
                f"SHT_QUICK_INIT, which is the only mode measured to run no "
                f"timing race at all.")
    if not seen:
        problems.append("shtnsmod.f90 has no klay assignment; the SHTns grid "
                        "mode is no longer where this check looks for it")
    return problems


def check_compiled_extensions() -> list[str]:
    """The f2py extensions pyburn needs are importable by THIS interpreter.

    `pyburn.readfile` imports `exoplasim.pyfft` (pyburn.py:662) to build the
    Gaussian grid, so a missing or mis-tagged extension makes every raw model
    file unreadable. The `.so` files are untracked build artifacts compiled for
    one CPython ABI; git tracks only `pyfft.f90` and `pyfft991.f90`. Rebuilding
    the venv on a new interpreter therefore invalidates them silently, and
    nothing notices until a run finishes and its postprocessing fails.

    It cost a pair of diagnostic runs to find, and the failure does not look
    like itself: `__init__.py:1150` turns any postprocessing exception into
    `_crash()`, which moves the working directory to `<run>_crashed/` and
    reports "ExoPlaSim has crashed or begun producing garbage". The model had
    integrated perfectly.

    Imports rather than inspecting filenames, because the question is whether
    THIS interpreter can load them, not whether a file of about the right name
    is on disk.
    """
    problems = []
    if not (ROOT / "vendor/exoplasim/exoplasim").is_dir():
        return problems
    import importlib
    for name in ("exoplasim.pyfft", "exoplasim.pyfft991"):
        try:
            importlib.import_module(name)
        except Exception as exc:
            problems.append(
                f"{name} will not import ({type(exc).__name__}), so pyburn "
                f"cannot read raw model output and every run will fail at "
                f"postprocessing. Rebuild it from its .f90 with numpy.f2py "
                f"against this interpreter.")
    return problems


def check_no_control_patch() -> list[str]:
    """No deliberate corruption is sitting in the committed model source.

    Several checks work by patching a defect INTO `plasim/src`, building, and
    requiring the result to fail -- that is what gives them teeth. They restore
    the tree with `git checkout` when they finish.

    That is safe until something commits while one is running. It happened:
    `verify_threaded_numerics.sh` had `lo = mypid * NHOR + 1` replaced by an
    offset that makes every thread's grid band overlap its neighbour's, a
    `git add -A` swept it into a commit, and the script's own trap then restored
    the tree to the corrupted commit. The model then crashed in the radiation on
    every bed and every seed, and looked for all the world like a regression in
    the change being tested -- it cost a long bisect, because the corruption was
    in the baseline too.

    So every control patch leaves a marker, and this refuses to let one be
    committed. It cannot catch a patch that only DELETES code, which is why the
    rule is that a control adds its marker first.
    """
    src = ROOT / "vendor/exoplasim/exoplasim/plasim/src"
    problems = []
    if not src.is_dir():
        return problems
    for path in sorted(src.glob("*.f90")):
        for n, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if "CONTROL PATCH IN PROGRESS" in line or "! CONTROL:" in line:
                problems.append(f"{path.name}:{n} carries a control patch marker; "
                                f"a verification script's deliberate corruption "
                                f"has been left in the source")
    return problems


def check_omp_directive_length() -> list[str]:
    """No OpenMP directive line in the model runs past the fixed-form limit.

    A `!$omp threadprivate(...)` list that overruns the line length is
    TRUNCATED, and the names past the cut are simply not threadprivate. Nothing
    fails to compile and nothing warns at the volume this build prints; what
    happens is that every thread shares one copy of a variable that was meant to
    be private, the last writer wins, and the model stops being reproducible run
    to run.

    That is exactly what happened when ten grid pointers were added to the list
    in `plasimmod.f90`: the line reached 142 characters, the tail was dropped,
    and two runs of one binary gave different restarts. It took a bisect against
    a known-good restart sha to find, because the symptom -- nondeterminism in a
    threaded build -- points at the new code rather than at a directive that
    reads correctly in the editor.

    The bar is 132 because these are `.f90` files, gfortran compiles them as
    FREE form, and 132 columns is gfortran's free-form default. It is a
    property of the language rather than of any flag, which is why
    `-ffixed-line-length-132` was measured to be a no-op and left out of
    `config/planet.yaml`'s `f90_opts` (`notes/audits/model-build-flags.md`).
    The check is over DIRECTIVES only: ordinary over-length lines elsewhere in
    the vendored sources are upstream's and are not silent in the same way.
    """
    src = ROOT / "vendor/exoplasim/exoplasim/plasim/src"
    problems = []
    if not src.is_dir():
        return problems
    for path in sorted(src.glob("*.f90")):
        for n, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            stripped = line.lstrip()
            if stripped.startswith("!$omp") and len(line) > 132:
                problems.append(f"{path.name}:{n} OpenMP directive is "
                                f"{len(line)} characters, over the 132 limit")
def _fortran_sources() -> list[Path]:
    """Every Fortran source this project holds: the model, and the gates' probes.

    The probe and gate sources under `exoplasim/scripts/` are compiled too, by
    the `verify_*.sh` arms, and a dropped continuation is the same defect there.
    """
    dirs = [ROOT / "vendor/exoplasim/exoplasim/plasim/src",
            ROOT / "exoplasim" / "scripts"]
    return sorted(f for d in dirs if d.is_dir() for f in d.glob("*.f90"))


def _without_trailing_comment(text: str) -> str:
    """`text` up to the first `!` that is not inside a character literal."""
    quote = ""
    for i, c in enumerate(text):
        if quote:
            if c == quote:
                quote = ""
        elif c in "'\"":
            quote = c
        elif c == "!":
            return text[:i]
    return text


def check_no_dropped_continuation() -> list[str]:
    """A line that continues a statement follows a line that says it continues.

    THE INSTANCE. `landmod.f90`'s `namelist/landmod_nl/` gained `,snowcovz` on a
    new line, and the line above it -- until then the statement's last -- kept no
    `&`. The statement therefore ended one name early and the added line began a
    NEW statement with `&`, which gfortran rejects as an invalid character in a
    name. `landmod.o` was the only object that failed, that stopped the link, and
    no binary could be built from the tip at all. It survived a merge and five
    parallel batches with every check green. world-adts.

    WHY A TEXT CHECK WHEN THERE IS ALSO A COMPILE GATE. Because the two fail
    differently and only one of them is free. This runs over every Fortran source
    the project holds, including the ones no current build compiles, and it names
    the line rather than the symbol the parser choked on. The compile gate is the
    authority on whether the source is legal; this is the one that catches the
    class in the commit that makes it, and it needs no compiler, no module
    directory and no SHTns.

    WHAT IT ASSERTS. These are `.f90` files and gfortran compiles them as FREE
    form, where a continued statement is marked by a trailing `&` and the
    continuation may repeat one at its start. So a source line whose first
    non-blank character is `&` or `,` must follow a source line that ends in `&`,
    and a chain that ends at the last line of a file is unterminated. Comment and
    blank lines between continuation lines are legal and are skipped rather than
    breaking the chain. `!$omp` directives carry their own chain and are tracked
    separately: a dropped `&` in a `threadprivate` list drops names off it
    silently, which is the failure `check_omp_directive_length` exists for
    reached by the other route.
    """
    problems = []
    for path in _fortran_sources():
        # One chain for statements, one for `!$omp` directives.
        continued = {False: False, True: False}
        opened_at = {False: 0, True: 0}
        for n, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
            stripped = raw.strip()
            if not stripped:
                continue
            directive = stripped[:5].lower() == "!$omp"
            if stripped.startswith("!") and not directive:
                continue
            body = stripped[5:] if directive else _without_trailing_comment(stripped)
            body = body.strip()
            if not body:
                continue
            if body[0] in "&," and not continued[directive]:
                what = "!$omp directive" if directive else "statement"
                problems.append(
                    f"{path.name}:{n} continues a {what} and the line above it "
                    f"has no '&', so this begins a new one: {stripped[:60]!r}")
            continued[directive] = body.endswith("&")
            opened_at[directive] = n
        for directive, still_open in continued.items():
            if still_open:
                what = "!$omp directive" if directive else "statement"
                problems.append(
                    f"{path.name}:{opened_at[directive]} ends the file with a "
                    f"{what} continued onto nothing")
    return problems


def check_diag_writes_are_answered() -> list[str]:
    """No write to the shared diagnostics unit is reachable by every thread.

    `nud` is unit 6 and `opendiag` opens it once, under a root guard, so it is
    the diagnostics file for the WHOLE thread team; `plasim.f90` then wraps
    `mpstart` through `mpstop` in one `!$omp parallel`, so every statement in
    the model runs on NPRO threads unless something on the path tests
    `mypid == NROOT`. An unguarded `write(nud,...)` is therefore NPRO
    interleaved copies in `plasim_diag`, in an order that is not the same twice.

    THE GAP THIS CLOSES. world-0ihs read all 928 sites by hand, found thirteen
    that every thread reached and fixed them, and left the convention resting on
    the next author reading the audit note. The model gains write sites; a
    convention nothing enforces is a convention until someone is in a hurry.

    `exoplasim/scripts/lint_diag_writes.py` holds the passes and the argument
    for each of the four answers it accepts. It is a text parse of `plasim/src`,
    no build and no run, and it answers a set of reduced fixtures on every
    invocation before it reports on the tree.
    """
    script = ROOT / "exoplasim" / "scripts" / "lint_diag_writes.py"
    if not script.is_file():
        return [f"{script.relative_to(ROOT)} is gone, and it is what holds the "
                "guard convention on the shared diagnostics unit"]
    r = subprocess.run([sys.executable, str(script)],
                       capture_output=True, text=True, cwd=ROOT)
    if r.returncode == 0:
        return []
    return [(r.stdout.strip() or r.stderr.strip() or "(no output)")]


def check_masked_only_locals_are_answered() -> list[str]:
    """No procedure local has all its definitions inside a `where`.

    A `where` selects which elements an ASSIGNMENT stores; it does not restrict
    which elements of a right-hand side the compiler EVALUATES, and a
    vectorising compilation evaluates all of them. So an automatic local with no
    initialiser, written only inside a mask and then read, is read on every lane
    the mask discarded as whatever the stack held -- and the declared
    `-ffpe-trap=invalid,zero,overflow` turns that into SIGFPE. Where the read
    sits under a different mask or none, the indeterminate lane is not discarded
    at all: it is stored.

    THE GAP THIS CLOSES. world-d016 found the class in `tands`, `mktsoil` and
    `mkdca` while re-deriving the masked-division population, fixed those three,
    and left the class itself unenumerated. Nothing in the tree could say
    whether there were others, and the answer was seventeen more in `rainmod`
    alone.

    `exoplasim/scripts/lint_masked_locals.py` holds the pass, the five
    directions it over-reports in and the two it cannot see, and the argument
    for why it is not a rule inside `lint_implicit_save.py`. It is a text parse
    of `plasim/src`, no build and no run, and it answers twenty reduced fixtures
    on every invocation before it reports on the tree.
    """
    script = ROOT / "exoplasim" / "scripts" / "lint_masked_locals.py"
    if not script.is_file():
        return [f"{script.relative_to(ROOT)} is gone, and it is what says "
                "whether a local's only definition is inside a mask"]
    r = subprocess.run([sys.executable, str(script)],
                       capture_output=True, text=True, cwd=ROOT)
    if r.returncode == 0:
        return []
    return [(r.stdout.strip() or r.stderr.strip() or "(no output)")]


# Command-line tools `docs/src/reference/environment.md` sends a reader to, mapped
# to the Arch package shipping each. The package belongs in the failure message
# because the tool name is usually not the package name: looking for a binary
# called `graphviz`, or supposing `h5diff` needs a package of its own, are both
# mistakes made while writing that section.
DOCUMENTED_TOOLS = {
    "ncdump": "netcdf",
    "ncks": "nco", "ncra": "nco", "ncdiff": "nco", "ncwa": "nco",
    "ncatted": "nco",
    "h5diff": "hdf5",
    "yq": "go-yq",
    "gdalinfo": "gdal", "ogrinfo": "gdal",
    "dot": "graphviz",
    "valgrind": "valgrind",
    "ncdu": "ncdu",
}

# Packages the same document names, which live in the venv rather than on the host.
DOCUMENTED_MODULES = ("dask", "flox", "bottleneck")


def check_documented_tools() -> list[str]:
    """Every tool `environment.md` names resolves, and every tool named here is in it.

    `environment.md` directs a reader to `ncdump`, NCO, `h5diff` and `yq`
    instead of a Python session, and carries the traps that come with them: that
    `ncwa` cannot area-weight this grid, and that the time axis is not a
    calendar. A document naming a tool the host lacks sends that reader to write
    the Python session anyway, having first wasted the lookup.

    Checked in BOTH directions. Forward catches a host that has changed under
    the document. Backward catches the document dropping a tool this list still
    asserts, which would leave the check guarding a claim nobody makes. Neither
    direction catches a tool ADDED to the document and not to this list; that is
    the seam, and the list is here rather than parsed out of the prose because
    deciding by regex over English what counts as a named tool is a worse
    failure than the one it would prevent.

    These read artifacts rather than produce them, so no pipeline step imports
    one and a failure here does not mean a run would be wrong. It means the
    documentation is.
    """
    doc = ROOT / "docs" / "src" / "reference" / "environment.md"
    if not doc.is_file():
        return [f"{doc.relative_to(ROOT)} is missing"]
    text = doc.read_text(encoding="utf-8")
    # First token of every inline code span: `ncdump -h` names ncdump.
    named = {span.split()[0] for span in re.findall(r"`([^`]+)`", text)
             if span.split()}

    problems = []
    for tool, package in sorted(DOCUMENTED_TOOLS.items()):
        if shutil.which(tool) is None:
            problems.append(
                f"environment.md names `{tool}` but it is not on PATH; "
                f"install the {package} package, or drop it from the document")
        if tool not in named:
            problems.append(
                f"`{tool}` is asserted here but environment.md no longer names it")
    for module in DOCUMENTED_MODULES:
        if importlib.util.find_spec(module) is None:
            problems.append(
                f"environment.md names {module} but it does not import; "
                f"it is a venv package, so check requirements.txt too")
        if module not in named:
            problems.append(
                f"{module} is asserted here but environment.md no longer names it")
    return problems


# An enumeration of the ladder -- "T21/T42/T85/T127/T170" -- is prose about the
# ladder, not a path, and the slashes in it are what made it look like one.
# Stripped before the path test rather than exempted after it.
RUNG_LADDER = re.compile(
    r"\b[tT](?:21|31|42|63|85|106|127|170)"
    r"(?:/[tT](?:21|31|42|63|85|106|127|170))+\b")
# A ladder rung appearing in something PATH-SHAPED: next to a separator, or in
# a filename with an extension. `T42` in prose, as a dict key in a per-rung
# table, or as a command-line default is not this; `inputs/t42/` and
# `orogen_T42_surf_0129.sra` are.
RUNG_IN_PATH = re.compile(
    r"""(?:^|[-/_'"(\[])[tT](?:21|31|42|63|85|106|127|170)(?:[-/_.'"),\]]|$)""")
# The registry itself, the build matrix that declares which rungs have
# binaries, and the reproducibility recipe that names one executable. A mapping
# rather than a tuple, so every exemption has to say why it is one.
RUNG_OWNER = {
    "lib/rungs.py": "the registry itself",
    "lib/gridding.py": "the one grid convention",
    "scripts/smoke_test.py": "this lint",
    "exoplasim/scripts/rebuild_binaries.py":
        "the build matrix, which declares which rungs have binaries",
    "exoplasim/scripts/shtns_variant_sweep.py": "names the rungs it sweeps",
    "exoplasim/scripts/reproducibility_matrix.py":
        "the recipe, which names one executable",
}


def check_no_rung_literal_in_a_path(files: list[Path]) -> list[str]:
    """No artifact path carries a resolution literal.

    SPAT-2. The ladder is T21/T42/T85/T127/T170 and the pipeline named T42 and
    T85 in the paths it wrote to, so a run at any other rung wrote its fields
    into a directory named for a rung it was not. It was not hypothetical: the
    graph declared `exoplasim/inputs/t42/orogen_T42_surf_0129.sra` while the
    staged file on disk was `inputs/t21/orogen_T21_surf_0129.sra`, so every
    ExoPlaSim input row in `config/pipeline.yaml` pointed at a file that did
    not exist and declared none of the ones that did.

    The rung belongs in the path -- an input family IS per-resolution -- but it
    has to arrive as a parameter. `{res}` and `{res_lower}` in the graph,
    `rungs.model_grid(config)` in a script, `rung_of_latitudes` where the
    artifact itself says how big it is.

    Prose is skipped exactly as in the other line lints, and the registry, the
    build matrix and the reproducibility recipe are exempt: those NAME rungs on
    purpose rather than building a path out of one.
    """
    bad = []
    for f in list(files) + [ROOT / "config" / "pipeline.yaml"]:
        rel = str(f.relative_to(ROOT))
        if rel in RUNG_OWNER:
            continue
        for i, line in enumerate(f.read_text(encoding="utf-8",
                                             errors="ignore").splitlines(), 1):
            stripped = line.strip()
            if (stripped.startswith("#") or "rung-ok" in line
                    or "`" in line or stripped.startswith("*")):
                continue
            probe = RUNG_LADDER.sub("<ladder>", line)
            if "/" not in probe and not re.search(r"\.(sra|nc|json|rest|x)\b", probe):
                continue                      # not path-shaped
            if RUNG_IN_PATH.search(probe):
                bad.append(f"{rel}:{i}: {stripped[:78]}")
    return bad



# The ladder spelled as a MAPPING rather than as a path: `T21) nlat=32`, or
# `{"T21": 32}`, or a bare `(32 64 128 192 256)`. A dict from rung to dimensions
# is not path-shaped, so `check_no_rung_literal_in_a_path` was never going to
# see it, and four of these survived the sweep that created lib/rungs.py.
#
# Exempt by NAME and with a reason, never by pattern: an exemption that matches
# a shape exempts the next copy too.
RUNG_TABLE_EXEMPT = {
    "lib/rungs.py": "the registry itself",
    "scripts/smoke_test.py": "this lint",
    "exoplasim/scripts/restart_convert_selftest.py":
        "enumerates SOURCE-TARGET grid pairs for the converter, which is a "
        "property of the conversion rather than a rung-to-dimension mapping; "
        "it maps no rung name to anything",
}


def check_no_rung_table_outside_rungs(files: list[Path]) -> list[str]:
    """The rung-to-dimension mapping is `lib/rungs.py` and nowhere else.

    SPAT-2 again, and the half its first lint could not reach. That lint tests
    for a rung inside something PATH-SHAPED, so a `case "$res" in T21) nlat=32`
    or a `{"T21": (32, 64)}` is out of scope by construction -- and those are
    what the surviving copies are. `lib/rungs.py`'s own docstring says it exists
    because "four scripts each carried their own copy of the mapping, and two of
    the four were missing rungs the others had"; four more were still carrying
    one when this was written, and one of them cited SPAT-2 in the comment
    directly above its own copy.

    Two shapes, because the copies come in two:

    - A rung literal on the same line as one of THAT RUNG's own dimensions, its
      truncation, latitudes or longitudes. Two such lines in a file is a table;
      one is a single case and might be a filename.
    - Three or more of the ladder's latitude counts as standalone integers on
      one line, which is the ladder written in dimensions with the rungs left
      off. `verify_gauss_weights.sh` had exactly that.

    Shell scripts are included, which the path lint's `d.glob("*.py")` never
    was.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    import rungs

    dims = {}
    for rung in rungs.RUNGS:
        nlat, nlon, ntru = rungs.geometry(rung)
        dims[rung.upper()] = {str(nlat), str(nlon), str(ntru)}
    nlats = {str(rungs.geometry(r)[0]) for r in rungs.RUNGS}
    pair = re.compile(r"\b[tT](21|31|42|63|85|106|127|170)\b")
    number = re.compile(r"(?<![\w.])\d{2,3}(?![\w.])")

    bad = []
    for f in files:
        rel = str(f.relative_to(ROOT))
        if rel in RUNG_TABLE_EXEMPT:
            continue
        hits = []
        for i, line in enumerate(f.read_text(encoding="utf-8",
                                             errors="ignore").splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("#") or "rung-table-ok" in line:
                continue
            numbers = set(number.findall(line))
            paired = [f"T{m.group(1)}" for m in pair.finditer(line)
                      if dims["T" + m.group(1)] & numbers]
            if paired:
                hits.append(f"{rel}:{i}: {', '.join(sorted(set(paired)))} beside "
                            f"its own dimensions -- {stripped[:60]}")
            elif len(numbers & nlats) >= 3:
                hits.append(f"{rel}:{i}: {len(numbers & nlats)} of the ladder's "
                            f"latitude counts -- {stripped[:60]}")
        # One line is a single case; two is a table.
        if len(hits) >= 2 or any("latitude counts" in h for h in hits):
            bad.extend(hits)
    return bad


def check_configured_grid() -> list[str]:
    """`resolution`, `latitudes` and `longitudes` are one fact, not three.

    `config/planet.yaml` carries all three because `read_sra` validates staged
    surface files against the latter two, so a resolution changed without them
    refuses its own inputs. That makes a stale pairing possible and silent
    until something reads a file. The rung table derives the dimensions and
    `rungs.model_grid` refuses a disagreement; this is where that runs without
    anyone having to call it.
    """
    import yaml
    sys.path.insert(0, str(ROOT / "lib"))
    import rungs
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    try:
        rungs.model_grid(cfg)
    except RuntimeError as exc:
        return [str(exc)]
    return []


def check_ladder_restatements() -> list[str]:
    """Every restatement of the ladder, against `lib/rungs.py`.

    The registry has carried this check since SPAT-2 and NOTHING RAN IT, which
    is the same shape of defect as the copies it exists to catch: a check with
    no caller certifies nothing, and `plasimmod.f90` names it in a comment
    directly above its own copy of the ladder. WORLD-J37 wired it in here.

    Three restatements, none removable: CMake has no Python, the vendored
    package must stay importable knowing nothing about this repository, and the
    Fortran comment is what someone editing `NLAT_ATM` by hand reads.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    import rungs
    return rungs.check_restatements(ROOT)


def check_snow_conductivity_restatements() -> list[str]:
    """Every restatement of the snow conductivity relation, against `lib/snow.py`.

    The climate column and the ecology column modelled the same snow with two
    different relations -- Fourteau et al. (2021) Eq. (18) in `landmod` and
    Sturm et al. (1997) in `soil.cpp` -- and at the declared snow density they
    differed by close to a factor of two, so one snowfall insulated one model's
    soil about twice as well as the other's. WORLD-GJOV made it one relation.

    They cannot be held together by matching NUMBERS: the vegetation model's
    snow density is prognostic across a compaction range and the climate
    model's is a single namelist key, so equal values at one density is a
    coincidence at one point of two curves. It is the RELATION that is shared,
    `lib/snow.py` is the one declaration of it, and the two restatements are
    literals only because a Fortran model and a C++ model cannot import a
    Python module at runtime. This is what makes a restatement honest;
    `biosphere/scripts/snow_thermal_gate.py` runs it too, on the register side.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    import snow
    return snow.check_restatements(ROOT)


def check_ladder_timestep_declarations() -> list[str]:
    """The three timestep quantities, each against what it restates.

    WORLD-F997: "the timestep at rung X" named three different things in this
    tree. They are three QUANTITIES rather than a disagreement, and the fix is
    that each is defined once in `lib/rungs.py` and checked against its source
    here -- the measured ceilings against the probe grid they are read out of,
    the escalation route against the chapter that decides it.

    A route table that agrees with no chapter and a ceiling table that agrees
    with no measurement are exactly how the T170 row came to declare a step the
    grid marks as starting and dying. WORLD-6QRR.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    import rungs
    return (rungs.check_stability_ceilings(ROOT)
            + rungs.check_timestep_restatements(ROOT))


def check_configured_timestep() -> list[str]:
    """`model.timestep_minutes` is a step the route runs this rung at.

    The defect this refuses: the step was one scalar with no dependence on
    `model.resolution`, so changing the rung left the step where it was and the
    first T42 arm ran at T21's 45 minutes and died in its forty-seventh orbit.
    WORLD-TD3, WORLD-J37.

    The DECLARED configuration is judged on both quantities, so an off-route
    step fails here even when it is below the ceiling. A diagnostic arm is
    launched from a generated config rather than from this file and reports
    instead of refusing; `lib/rungs.py:timestep_problems` carries the split.

    THE CONTROLS are below and are what make this a test rather than a
    description. There are two because there are two ways to be off the ladder
    and each has its own right answer: a step above the rung's measured ceiling,
    and a step the route does not run this rung at even though the rung could
    take it. The second is the defect above, and it is the one a single
    ceiling-shaped control would miss.
    """
    import yaml
    sys.path.insert(0, str(ROOT / "lib"))
    import rungs
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    try:
        active, rung = rungs.configured_timestep(cfg)
    except RuntimeError as exc:
        return [str(exc)]
    problems = rungs.timestep_problems(rung, active)

    # NEGATIVE CONTROL ONE: above the ceiling. T85 is measured to refuse at 60,
    # so a configuration declaring it has to fail here.
    if not rungs.timestep_problems("T85", 60.0):
        problems.append(
            "the timestep check does not refuse T85 at dt 60, which is above "
            "the step T85 is measured to refuse at: a check that cannot fail "
            "is not a check")
    # NEGATIVE CONTROL TWO: below every ceiling and off the route. This is the
    # defect's own shape -- the step is perfectly safe and simply is not the one
    # the ladder runs this rung at -- and a ceiling-only check passes it.
    if not rungs.timestep_problems("T21", 30.0):
        problems.append(
            "the timestep check does not refuse T21 at dt 30, which is well "
            "under T21's ceiling and is not on the route: a check that only "
            "sees ceilings cannot catch a rung run at another rung's step")
    # And the other direction: the declared route must be accepted, or the
    # control above is passing on a check that refuses everything.
    for rung_ok, dt_ok in rungs.ESCALATION_ROUTE:
        if rungs.timestep_problems(rung_ok, dt_ok):
            problems.append(
                f"the timestep check refuses {rung_ok} at dt {dt_ok}, which is "
                "on the escalation route it is meant to accept")
    return problems


def check_gate_filter_matches_config() -> list[str]:
    """A transform gate runs the spectral filter `config/planet.yaml` declares.

    `verify_shtns_equivalence.f90` asserts "THE CONFIGURATION UNDER TEST IS THE
    ONE THE MODEL RUNS" and then sets `filterkappa` and `nfilterexp` as
    literals. `filter_power` moved from 8 to 16 in config and the three drivers
    did not, so the gate certified a filter the model had stopped using and
    said nothing about it.

    The strength is not cosmetic to the verdict. Both arms of the comparison
    carry the same `skspgp(n)`, so it divides out of a per-mode ratio, but the
    gate's error is a field norm: the filter reweights the residual's spectrum
    against a denominator the low modes own.
    `exp(-kappa*x**16)/exp(-kappa*x**8)` peaks at `exp(kappa/4)`, which is 7.39
    at kappa 8 and falls at `n/NTRU = (1/2)**(1/8) = 0.917` -- the mid-to-high
    band where the SHTns residual is largest. A gate left at gamma 8 sees that
    band at a seventh of the amplitude the shipped configuration gives it.

    Both directions, so config and the drivers cannot drift apart in either.
    """
    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    model = cfg.get("model", {})
    want = {"filterkappa": ("filter_kappa", float(model["filter_kappa"])),
            "nfilterexp": ("filter_power", float(model["filter_power"]))}
    assign = re.compile(r"^\s*(filterkappa|nfilterexp)\s*=\s*([0-9.eE+-]+)\s*(!.*)?$")
    bad, seen = [], set()
    for f in sorted((ROOT / "exoplasim" / "scripts").glob("*.f90")):
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            m = assign.match(line)
            if not m:
                continue
            key, got = m.group(1), float(m.group(2))
            cfg_key, cfg_val = want[key]
            seen.add(key)
            if got != cfg_val:
                bad.append(f"{f.relative_to(ROOT)}:{i}: {key} = {m.group(2)}, "
                           f"but config/planet.yaml model.{cfg_key} is {cfg_val:g}")
    for key, (cfg_key, _) in want.items():
        if key not in seen:
            bad.append(f"no driver under exoplasim/scripts sets {key}; the "
                       f"config.{cfg_key} comparison has nothing to check")
# Settings that configure() re-copies over on every continuation, so a run that
# does not REAPPLY them integrates its later segments on the compiled defaults.
# The file records three that were lost this way; the list is here so a fourth
# cannot be added to run_exoplasim alone.
REDECLARED_ON_CONTINUE = (
    "declare_hyperdiffusion", "declare_energy_fixer", "declare_robert_filter",
    "declare_dynamics_only", "declare_conversion_time_level",
    "declare_dealias_conversion",
)

# The same rule for declarations the CALLER parameterises rather than reading
# from the config, so their call sites cannot be matched on `(model, config)`.
# They are on exactly the same footing -- configure() overwrites the namelist
# either way -- and only the shape of the check differs.
REDECLARED_ON_CONTINUE_ANY_ARGS = (
    "set_low_io", "declare_ecological_stream",
)


def check_continuation_redeclares_everything() -> list[str]:
    """Whatever a prepare declares, a continuation declares again.

    `model.configure()` re-copies the shipped namelists over the configured
    ones every time it runs, so a setting written at prepare time and not
    rewritten on continuation silently reverts partway through a run. That has
    now happened three times in this file's history -- the stellar spectrum,
    the shortwave gas weights, and hyperdiffusion.

    Hyperdiffusion was the expensive one, because its fallback is not "off" but
    a DIFFERENT OPERATOR: `plasim.f90:1443` gives T21 ndel 2 where the config
    derives 4, grad^4 instead of grad^8. The production baseline integrated 84
    of its 85 orbits that way. world-1nz.
    """
    cont = (ROOT / "exoplasim" / "scripts" / "continue_exoplasim.py").read_text(
        encoding="utf-8")
    prep = (ROOT / "exoplasim" / "scripts" / "run_exoplasim.py").read_text(
        encoding="utf-8")
    bad = []
    for name in REDECLARED_ON_CONTINUE:
        if f"def {name}(" not in prep:
            bad.append(f"{name} is asserted here and run_exoplasim.py no "
                       "longer defines it")
        elif f"{name}(model, config)" not in cont:
            bad.append(f"run_exoplasim declares {name} and "
                       "continue_exoplasim.py never reapplies it, so a "
                       "continued segment reverts to the compiled default")
    for name in REDECLARED_ON_CONTINUE_ANY_ARGS:
        if f"def {name}(" not in prep:
            bad.append(f"{name} is asserted here and run_exoplasim.py no "
                       "longer defines it")
        elif f"{name}(" not in cont.replace(f"def {name}(", ""):
            bad.append(f"run_exoplasim declares {name} and "
                       "continue_exoplasim.py never reapplies it, so a "
                       "continued segment reverts to the compiled default")
    return bad


def check_per_level_namelist_keys_cover_every_level() -> list[str]:
    """A namelist key backed by an NLEV array is written for every level.

    `ndel`, `tdissd`, `tdissz`, `tdisst` and `tdissq` are declared `(NLEV)` in
    plasimmod, and a Fortran namelist scalar assigns ELEMENT ONE and leaves the
    rest. Written as scalars they reached the model top and nine levels in ten
    kept whatever `readnl` had preset. The model echoes what it read, and it
    read `NDEL=4, 9*2` at T21: world-1nz's grad^8 on one level and grad^4 on
    the other nine, with humidity damped 7.4 times too hard there.

    At T42 the same scalar landed in the wrong UNIT. `readnl`'s `NTRU==42`
    branch fills the arrays in seconds, `dayseccheck` decides days-against-
    seconds from MAXVAL over the whole array, and the preset's 65664 suppressed
    the conversion element one needed -- leaving the top level damped 86400
    times too hard and nothing else touched. world-td3.

    This CALLS the writer against the real config and reads what it would put
    in the namelist, rather than matching its source text: the first version
    matched text and passed a scalar written through a helper.
    """
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    import yaml
    import run_exoplasim

    per_level = ("NDEL", "TDISSD", "TDISSZ", "TDISST", "TDISSQ")
    config = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    nlev = int(config["model"]["layers"])

    class _Recorder:
        def __init__(self):
            self.written = {}

        def _edit_namelist(self, namelist, key, value):
            self.written[key] = value

    def _count(value: str) -> int:
        """How many array elements this namelist value actually supplies."""
        total = 0
        for item in str(value).split(","):
            item = item.strip()
            if not item:
                continue
            total += int(item.split("*", 1)[0]) if "*" in item else 1
        return total

    bad = []
    recorder = _Recorder()
    run_exoplasim.declare_hyperdiffusion(recorder, config)
    for key in per_level:
        if key not in recorder.written:
            bad.append(f"declare_hyperdiffusion no longer writes {key}, which "
                       "this check asserts is declared for every level")
            continue
        supplied = _count(recorder.written[key])
        if supplied != nlev:
            bad.append(
                f"declare_hyperdiffusion writes {key} = "
                f"{recorder.written[key]!r}, which supplies {supplied} of "
                f"{nlev} levels -- the rest keep the model's own preset")
    return bad


def check_autocorrelation_estimator() -> list[str]:
    """The one standard error over a series with memory recovers a known answer.

    `lib/autocorrelation.py` decides whether two runs are at the same
    equilibrium and whether one has converged, and before it existed both of
    those took their error bar from the raw sample count. A test of it has to
    be able to FAIL, so it is run on synthetic AR(1) series whose integrated
    autocorrelation time is known in closed form, tau = (1+r)/(1-r), and
    against the empirical spread of many independent window means -- the
    quantity the estimator claims to predict.

    EVERY TOLERANCE HERE WAS FIXED BEFORE THE ESTIMATOR WAS RUN ON A MODEL
    SERIES, and the seed is fixed so the check is a check and not a lottery.
    The controls are the specific wrong answers: an independent series must not
    come back with memory, a trended series must be refused as non-stationary,
    and the naive count must be shown to understate the error it was being used
    for.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import numpy as np
        import autocorrelation as ac
    except ImportError as exc:
        return [f"lib/autocorrelation.py does not import: {exc}"]

    TAU_RTOL = 0.10          # tau, from a long series, to a tenth
    SE_RTOL = 0.15           # predicted window-mean SE against the empirical spread
    rng = np.random.default_rng(20260825)

    BURN = 500                # discarded so every series starts stationary

    def ar1(n, r, count=1):
        """`count` independent AR(1) series of length n, unit marginal variance.

        Generated across the ensemble rather than one at a time: the recursion
        is over time and vectorised over the count, which is what keeps this a
        check somebody will actually leave in the gate.
        """
        e = rng.standard_normal((count, n + BURN))
        x = np.zeros((count, n + BURN))
        s = np.sqrt(1.0 - r * r)
        for i in range(1, n + BURN):
            x[:, i] = r * x[:, i - 1] + s * e[:, i]
        return x[:, BURN:]

    bad = []
    for r in (0.0, 0.3, 0.6, 0.8):
        want = (1.0 + r) / (1.0 - r)
        got = ac.integrated_time(ar1(20000, r)[0])["tau"]
        if abs(got - want) > TAU_RTOL * want:
            bad.append(f"AR(1) at r={r} has tau {want:.3f} and the estimator "
                       f"returned {got:.3f}")

    # THE CLAIM ITSELF: sigma * sqrt(tau / n) is the standard deviation of a
    # window mean. Checked against the spread of 1500 independent windows.
    for r, n in ((0.6, 20), (0.8, 20), (0.6, 60)):
        tau = (1.0 + r) / (1.0 - r)
        empirical = float(np.std(ar1(n, r, 1500).mean(axis=1), ddof=1))
        predicted = float(np.sqrt(tau / n))
        naive = float(1.0 / np.sqrt(n))
        if abs(predicted - empirical) > SE_RTOL * empirical:
            bad.append(f"at r={r}, n={n} a window mean's spread is "
                       f"{empirical:.4f} and sigma*sqrt(tau/n) predicts "
                       f"{predicted:.4f}")
        # The control: the count-based form must be visibly wrong here, or
        # this check is not testing anything the project did not already have.
        if naive > 0.75 * empirical:
            bad.append(f"at r={r}, n={n} the count-based standard error "
                       f"{naive:.4f} is not detectably below the true "
                       f"{empirical:.4f}, so this case proves nothing")

    # A trend is not variability, and a tau taken across one describes the
    # approach. The guard must refuse the second series and accept the first.
    flat = ar1(200, 0.6)[0]
    if not ac.stationary_enough(flat)["stationary"]:
        bad.append("the stationarity guard refused a stationary AR(1) series")
    trended = flat + 0.05 * np.arange(flat.size)
    if ac.stationary_enough(trended)["stationary"]:
        bad.append("the stationarity guard accepted a series with a trend "
                   "larger than its own scatter")

    # The window that a target standard error needs, against the relation it
    # inverts. An identity, so it is exact rather than tolerant.
    n = ac.samples_for_standard_error(sigma=0.07, tau=4.2, target=0.02)
    if abs(0.07 * np.sqrt(4.2 / n) - 0.02) > 1e-12:
        bad.append("samples_for_standard_error does not invert its own relation")
    return bad


def check_restart_schema_covers_the_model() -> list[str]:
    """Every restart record the model writes has a policy, with the right reset.

    THIS IS THE ANSWER TO "what do I have to update if I add a restart record".
    `exoplasim/scripts/restart_schema.py` is the one place, and this makes the
    tree go red rather than leaving it to be found by whoever next converts a
    restart or seeds a run.

    It fails in the direction that matters and in both directions at once. A
    `put_restart_*` call the policy does not name is an unknown record, and
    `convert_restart.py` treats one as a hard error rather than passing it
    through -- because a record nobody has classified has no defensible
    behaviour across a resolution change. A policy entry no call site writes is
    a record that has been removed from the model and left behind here. And the
    reset column is checked against `outreset` and its per-module equivalents,
    which is what stops "an accumulator's clean value is zero" from silently
    becoming false: it already is for four of them.

    `reset_restart_accumulators.py` reads the same policy, so a gap here is a
    seeded run opening mid-window on the donor's partial accumulation.

    IT ALSO CHECKS THE RESOLVER, not only the names. A shape may be written in
    a symbol no `Geometry` knows -- `dwatcl(NHOR,NLSOILWX)` was, from the day
    the land column landed -- and the name check cannot see it, because the
    name is covered and it is the LENGTH that cannot be predicted. That gap
    reached the tree through `convert_restart.py --self-test`, which nothing
    runs; `check_every_shape_resolves` needs no restart file and runs here.

    It also holds every accumulator's DECLARED initial value against its reset.
    A cold run's first output window accumulates from the declaration with no
    reset before it, so the two disagreeing makes the first record of every
    cold start a different quantity from the rest of the run. That is what
    `atsami` was: declared at 0.0, reset to 1.0e10, and a running minimum, so
    the first record reported a minimum of 0 K.
    """
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    try:
        import restart_schema
    except ImportError as exc:
        return [f"exoplasim/scripts/restart_schema.py does not import: {exc}"]
    src = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src"
    if not src.is_dir():
        return [f"{src} is missing; the vendored model source moved"]
    return (restart_schema.check_policy_covers_source(src)
            + restart_schema.check_every_shape_resolves(src)
            + restart_schema.check_first_window_matches_the_rest(src))


def _fortran_call_args(text: str) -> list[str]:
    """Split one Fortran argument list on top-level commas."""
    parts, depth, current = [], 0, []
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(char)
    parts.append("".join(current).strip())
    return parts


def _codes_written_to_unit(src: Path, unit: str) -> set[int]:
    """Every service-format code the model writes to one output unit.

    Parsed from the writers themselves. `writegp(unit, field, code, level)` and
    `writesp(unit, field, code, level, scale, offset)` put the code third;
    `writescalar(unit, value, code)` puts it last. A call whose unit is a
    variable rather than a literal belongs to whichever stream the caller passes
    and is not counted for either.
    """
    codes: set[int] = set()
    for path in sorted(src.glob("*.f90")):
        for match in re.finditer(r"call\s+write(gp|sp|scalar)\s*\((.*)\)\s*$",
                                 path.read_text(), flags=re.M):
            args = _fortran_call_args(match.group(2))
            if len(args) < 3 or args[0] != unit:
                continue
            raw = args[2] if match.group(1) in ("gp", "sp") else args[-1]
            if re.fullmatch(r"\d+", raw):
                codes.add(int(raw))
    return codes


def _codes_pyburn_derives(pyburn_source: str) -> set[int]:
    """Every code pyburn has a derivation branch for.

    The branches test `key==str(<name>code)` against module-level integers, so
    the two together give the set. A code in neither the model's writers nor
    here is one pyburn's dispatch marks derived, matches no branch, and drops
    without a word.
    """
    numbers = dict(re.findall(r"^(\w+code)\s*=\s*(\d+)", pyburn_source,
                              flags=re.M))
    return {int(numbers[name])
            for name in re.findall(r"key==str\((\w+code)\)", pyburn_source)
            if name in numbers}


def check_requested_codes_are_produced() -> list[str]:
    """Every postprocessor code asked for is one something can answer.

    A code in `REGULAR_CODES` or `SNAPSHOT_CODES` is satisfied two ways: the
    model writes it to that stream's unit, or `pyburn` derives it from codes
    that are written. A code that is neither is a request that LOOKS satisfied.
    `pyburn.dataset` finds it absent from the raw data, sets `derived=True`,
    falls past every branch and drops it silently, so the run finishes clean and
    the product is missing a field a consumer read off the list and expected.

    This is `world-dy7a`, and 168 was not alone: 163, 171 and 238 were in both
    lists on the same terms, and none of the four appears in any climatology.
    The list is a contract with whoever writes against it, which is why the
    failure has to be at the point the list is edited rather than at the point
    someone looks for the field.

    Both directions are NOT checked. A code the model writes that no list
    requests is an ordinary and deliberate state -- the model writes far more
    than any product carries -- and `world-j0az` weighs those one at a time.
    """
    src = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src"
    pyburn = ROOT / "vendor" / "exoplasim" / "exoplasim" / "pyburn.py"
    if not src.is_dir():
        return [f"{src} is missing; the vendored model source moved"]
    if not pyburn.is_file():
        return [f"{pyburn} is missing; the postprocessor moved"]
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    try:
        import run_exoplasim
    except ImportError as exc:
        return [f"exoplasim/scripts/run_exoplasim.py does not import: {exc}"]

    derived = _codes_pyburn_derives(pyburn.read_text())
    problems = []
    for name, unit, requested in (
            ("REGULAR_CODES", "40", run_exoplasim.REGULAR_CODES),
            ("SNAPSHOT_CODES", "140", run_exoplasim.SNAPSHOT_CODES)):
        written = _codes_written_to_unit(src, unit)
        if not written:
            problems.append(f"no code was parsed as written to unit {unit}; "
                            f"the writers moved or changed shape")
            continue
        for code in requested:
            if code not in written and code not in derived:
                problems.append(
                    f"{name} asks for code {code} and nothing writes it to "
                    f"unit {unit} or derives it in pyburn, so the product "
                    f"drops it without a word")
    return problems


# The two steps whose script deliberately resolves NEITHER climatology, and what
# names the file instead. Both take it from the caller, which is a legitimate
# third answer and not a gap: the drift this guards against is a script that
# silently resolves the OTHER one. A new entry here is a decision that has to be
# argued, which is the point of writing it down rather than inferring it.
CLIMATOLOGY_FROM_THE_CALLER = {
    "analyze_climatology": "--label is required and both product names derive "
                           "from it; a fixed default here reported 280.9 K for "
                           "a 293.8 K world",
    "ice_mask": "--climatology is required=True",
    "ocean_forcing": "the forcing contract takes an explicit pass artifact; "
                       "resolving config's current baseline could select the "
                       "newer side of a chronological pass pair",
    "ocean_loop_convergence": "the assessor takes both consecutive pass "
                               "artifacts from the caller; there is no single "
                               "current climatology in a pairwise test",
}


RESOLVERS = {"bootstrap_climatology_path", "climatology_path",
             "best_available_climatology", "climatology_path_for_state"}

# Which resolvers each declared need permits. The need is what must EXIST; the
# resolver is what the step READS, and the two are different statements, which
# is why one need admits two resolvers and the other admits one.
PERMITTED_RESOLVERS = {
    "bootstrap_climatology": {"bootstrap_climatology_path",
                              "best_available_climatology",
                              "climatology_path_for_state"},
    "baseline_climatology": {"climatology_path"},
}


def check_climatology_needs_match_call_sites() -> list[str]:
    """Each step resolves a climatology `config/pipeline.yaml` permits it.

    THERE ARE TWO CLIMATOLOGIES AND THE GRAPH HAS ALWAYS SAID WHICH IS WHICH.
    The bootstrap is the run on terrain-only surface fields and exists to
    produce the climate the derived fields are built FROM; the baseline is the
    run on those fields once they exist. `config/planet.yaml` carried ONE key,
    so `lib/paths.py:climatology_path` returned the same file to all eleven
    steps and the distinction reached nothing -- a step wanting the bootstrap
    read an artifact that cannot exist on a first pass.

    THE EDGE AND THE READ ARE DIFFERENT STATEMENTS, and that is why a
    `bootstrap_climatology` step has two permitted resolvers. The edge says
    which artifact must EXIST before the step can run, and for these steps it is
    the bootstrap: they run on a first pass, when no baseline exists at all.
    What the step should READ is decided by whether its answer depends on the
    climate STATE. Where it does, `best_available_climatology` returns the
    baseline once one is named and the bootstrap before that, and says which --
    a step pinned to the bootstrap forever is right on the first pass and
    holding the loop back on every pass after. Where the answer depends only on
    the model calendar or the grid, `bootstrap_climatology_path` is right and a
    best-available call would make the step newly need a baseline for nothing.

    A `baseline_climatology` step still gets exactly one resolver:
    `climatology_path`. Its edge already says a baseline must exist, so there is
    no earlier stage for it to prefer and nothing to choose between.

    Calling more than one resolver in one script is a failure whatever the
    declaration, because then the step reads two stages of the world at once.

    THE STAMP IS PART OF THE PROPERTY. A script that calls
    `best_available_climatology` must also write `climatology_stage` into what
    it produces: the choice between stages is defensible because it is recorded,
    so an unstamped product would leave a reader unable to tell a first-pass
    artifact from a later one. `lib/paths.py` carries that argument in full.

    WHAT IT DOES NOT COVER, said here so nobody reads it as more: a step that
    declares neither key and reads a climatology anyway is invisible to this,
    because the graph does not state which one it wants. Adding a key to those
    steps is what would make them checkable.
    """
    import yaml
    graph = yaml.safe_load(
        (ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    problems = []
    for step in graph.get("steps") or []:
        needs = [n for n in (step.get("needs") or [])
                 if n in PERMITTED_RESOLVERS]
        if not needs:
            continue
        step_id = step["id"]
        if len(needs) > 1:
            problems.append(
                f"step {step_id} needs both climatologies; they are two "
                f"artifacts and a step reading both reads two worlds")
            continue
        script = ROOT / str(step["script"])
        if not script.is_file():
            continue                    # check_registered_paths_exist owns this
        source = script.read_text(encoding="utf-8")
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue                    # check_modules_parse owns this
        called = set()
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            f = node.func
            name = f.id if isinstance(f, ast.Name) else (
                f.attr if isinstance(f, ast.Attribute) else None)
            if name in RESOLVERS:
                called.add(name)
        permitted = PERMITTED_RESOLVERS[needs[0]]
        wrong = called - permitted
        if "climatology_path_for_state" in called and step_id != "soil":
            problems.append(
                f"step {step_id} uses climatology_path_for_state(), which is "
                "reserved for soil's required --state contract")
        if wrong:
            problems.append(
                f"step {step_id} declares `needs: {needs[0]}` and "
                f"{step['script']} calls {sorted(wrong)[0]}(), which is not "
                f"one of {sorted(permitted)}. The bootstrap and the baseline "
                f"are different artifacts, not two versions of one")
        elif len(called) > 1:
            problems.append(
                f"step {step_id}'s {step['script']} calls {sorted(called)}. "
                f"One script resolves ONE climatology, or the step reads two "
                f"stages of the world at once")
        elif not called and step_id not in CLIMATOLOGY_FROM_THE_CALLER:
            problems.append(
                f"step {step_id} declares `needs: {needs[0]}` and "
                f"{step['script']} resolves no climatology at all. Call one of "
                f"{sorted(permitted)} from lib/paths.py, or take the path from "
                f"the caller and say so in CLIMATOLOGY_FROM_THE_CALLER with "
                f"the argument")
        if "best_available_climatology" in called:
            problems.extend(stamp_problems(step_id, step["script"], tree))
    return problems


def stamp_problems(step_id: str, script: str, tree: ast.AST) -> list[str]:
    """The stage the resolver returned is the stage the product records.

    Stronger than "the key is present", and that is the point: a script that
    stamps the string `bootstrap` or a variable from somewhere else satisfies
    presence and tells the reader something false. The property is a
    pass-through -- the second element of what `best_available_climatology`
    returned has to be what a `climatology_stage` entry is set to somewhere in
    the script. One such entry is enough; several products may carry the value
    on from the first.
    """
    stage_names = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        call = node.value
        if not isinstance(call, ast.Call):
            continue
        f = call.func
        name = f.id if isinstance(f, ast.Name) else (
            f.attr if isinstance(f, ast.Attribute) else None)
        if name != "best_available_climatology":
            continue
        for target in node.targets:
            if isinstance(target, ast.Tuple) and len(target.elts) == 2 \
                    and isinstance(target.elts[1], ast.Name):
                stage_names.add(target.elts[1].id)
    if not stage_names:
        return [f"step {step_id}'s {script} calls "
                f"best_available_climatology() without unpacking the stage it "
                f"returns. It returns `(path, stage)` and the stage is not "
                f"optional: it is what the product records"]
    stamped = False
    for node in ast.walk(tree):
        value = None
        if isinstance(node, ast.Dict):
            for key, val in zip(node.keys, node.values):
                if isinstance(key, ast.Constant) \
                        and key.value == "climatology_stage":
                    value = val
                    if any(isinstance(n, ast.Name) and n.id in stage_names
                           for n in ast.walk(value)):
                        stamped = True
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Attribute) \
                        and target.attr == "climatology_stage":
                    if any(isinstance(n, ast.Name) and n.id in stage_names
                           for n in ast.walk(node.value)):
                        stamped = True
    if stamped:
        return []
    return [f"step {step_id}'s {script} reads the best available climatology "
            f"and never records {sorted(stage_names)[0]} as "
            f"`climatology_stage`. The choice between the two stages is "
            f"defensible because it is stamped; a product carrying no stage, "
            f"or a stage that is not the one resolved, cannot be told from a "
            f"first-pass one"]


def check_best_available_climatology_resolves_by_stage() -> list[str]:
    """The best-available resolver returns the right stage in all three cases.

    THREE RIGHT ANSWERS, not three things that merely differ. With a baseline
    named it must return the baseline and say `baseline`; with the baseline
    null it must return the bootstrap and say `bootstrap`; with NEITHER named
    it must raise rather than guess. The third is the one that makes this not a
    fallback: a fallback has somewhere to go when it runs out of inputs, and
    this has nowhere.

    Driven against a synthetic config root rather than the project's own, so
    the check states the property instead of restating today's
    `config/planet.yaml`. The paths need not exist: the resolver chooses on
    what config DECLARES, which is itself part of the property -- choosing on
    what happens to be on disk is exactly the silent degradation the no-fallback
    rule is about.
    """
    import tempfile
    sys.path.insert(0, str(ROOT / "lib"))
    import paths as paths_lib

    cases = [
        ("baseline named", "baseline_climatology: a/base.nc\n"
                           "bootstrap_climatology: a/boot.nc\n",
         "a/base.nc", "baseline"),
        ("baseline null", "baseline_climatology: null\n"
                          "bootstrap_climatology: a/boot.nc\n",
         "a/boot.nc", "bootstrap"),
    ]
    problems = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "config").mkdir()
        cfg = root / "config" / "planet.yaml"
        for label, text, want_path, want_stage in cases:
            cfg.write_text(text, encoding="utf-8")
            try:
                got = paths_lib.best_available_climatology(root=root)
            except SystemExit as exc:
                problems.append(
                    f"best_available_climatology raised with {label}: {exc}")
                continue
            if got.path != root / want_path:
                problems.append(
                    f"with {label} the resolver returned {got.path}, not "
                    f"{root / want_path}")
            if got.stage != want_stage:
                problems.append(
                    f"with {label} the resolver reported stage "
                    f"{got.stage!r}, not {want_stage!r}")
        # NEITHER NAMED. Raising is the answer; returning anything at all here
        # would be the resolver inventing a world.
        cfg.write_text("baseline_climatology: null\n"
                       "bootstrap_climatology: null\n", encoding="utf-8")
        try:
            got = paths_lib.best_available_climatology(root=root)
        except SystemExit:
            pass
        else:
            problems.append(
                f"with neither climatology named the resolver returned "
                f"{got.path} at stage {got.stage!r} instead of raising. That "
                f"is the fallback the no-fallback rule forbids")
        # An override is labelled by comparison, not assumed to be either.
        cfg.write_text("baseline_climatology: a/base.nc\n"
                       "bootstrap_climatology: a/boot.nc\n", encoding="utf-8")
        for given, want_stage in ((root / "a" / "base.nc", "baseline"),
                                  (root / "a" / "boot.nc", "bootstrap"),
                                  (root / "a" / "elsewhere.nc", "unnamed")):
            got = paths_lib.best_available_climatology(given, root=root)
            if got.stage != want_stage:
                problems.append(
                    f"--climatology {given.name} was labelled {got.stage!r}, "
                    f"not {want_stage!r}")
    return problems


def check_no_shadowed_imports(files: list[Path]) -> list[str]:
    """A name bound by `import X` is never rebound to something else.

    `build_soil.py` did `import climatology` at module level and then
    `climatology = args.climatology or climatology_path()` inside main(), so
    every `climatology.annual_mean_of(...)` after it called a method on a
    PosixPath. The module was still imported and the name was still defined,
    which is why the import check and pyflakes F821 both pass on it -- the
    failure is a type, not a name, and it surfaces only when the function is
    actually run. Two scripts carried it, and two others had already dodged it
    with an alias (`import climatology as clim`, `as climatology_lib`), which
    is the fix and is why the rule is stated as "never rebound" rather than
    "beware": the alias makes the collision impossible instead of remembered.

    The rebinding is function-LOCAL, so it does not corrupt other functions --
    only uses inside the shadowing scope break. The check ignores scope anyway
    and reports any rebinding at all, because a name that means a module in one
    function and a path in the next is worth an alias regardless.
    """
    problems = []
    for path in files:
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.asname is None and "." not in a.name:
                        imported.add(a.name)
        if not imported:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store) \
                    and node.id in imported:
                problems.append(
                    f"{path.relative_to(ROOT)}:{node.lineno} rebinds `{node.id}`, "
                    f"which is an imported module; import it under an alias")
            # Parameters shadow too, and that is the worse case: the whole
            # function body sees the wrong object, not just the lines after an
            # assignment. `thermostat_efficiency.read_weathering(climatology:
            # Path)` called `climatology.annual_mean_of` on its own argument.
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                a = node.args
                for arg in (a.posonlyargs + a.args + a.kwonlyargs
                            + ([a.vararg] if a.vararg else [])
                            + ([a.kwarg] if a.kwarg else [])):
                    if arg.arg in imported:
                        problems.append(
                            f"{path.relative_to(ROOT)}:{arg.lineno} parameter "
                            f"`{arg.arg}` of {node.name}() shadows an imported "
                            f"module; import it under an alias")
    return problems

# The Python the model itself ships, which this project imports and runs. The
# name check below reads it; the other lints do not, because they are about how
# THIS project keeps its artifacts and the vendored package does not keep any.
VENDOR_PY_DIRS = [ROOT / "vendor" / "exoplasim" / "exoplasim"]

# Bound by the interpreter in every module namespace, so they are never
# "unbound" however the source reads.
MODULE_DUNDERS = frozenset({"__file__", "__name__", "__doc__", "__package__",
                            "__spec__", "__loader__", "__builtins__",
                            "__debug__", "__path__", "__all__", "__dict__"})

# symtable makes a scope for a deferred annotation (PEP 649) and for the type
# parameters of a generic. Neither has an ast node to report against and both
# resolve lazily, so the sweep does not descend into them.
_UNREPORTED_SCOPES = frozenset({"annotation", "type alias", "type parameters",
                                "type parameter", "TypeVar bound"})


def _module_level_bindings(top: symtable.SymbolTable) -> set[str]:
    """Names bound at module scope, asked of the resolver rather than the text.

    A `def` inside a module-level `for` binds a module global just as a `def` in
    column zero does, and reading `tree.body` alone misses it -- which this check
    reported as an unbound name on its own first run.
    """
    return {s.get_name() for s in top.get_symbols()
            if s.is_assigned() or s.is_imported()}


def _implicit_bindings(path: Path, tree: ast.AST) -> tuple[set[str], list[str]]:
    """Module-scope names bound by something a scope walk cannot see.

    Two of them, both decidable:

    * `from X import *` binds whatever X exports. The check resolves X against
      the file's own package directory and reads that module's module-level
      bindings. A star import it cannot resolve is REPORTED rather than waved
      through, because silently widening the allowed set is how a check of this
      kind stops being able to fail.
    * `import pkg.sub` inside `pkg/__init__.py` binds `sub` as well as `pkg`.
      Importing a submodule sets it as an attribute of the parent package, and
      inside the package's own `__init__` the parent's attributes ARE the module
      globals. `exoplasim/__init__.py` reaches `gcmt` and `pyburn` that way.
    """
    names, problems = set(), []
    package = path.parent.name if path.name == "__init__.py" else None
    for node in tree.body:
        if isinstance(node, ast.Import) and package is not None:
            for a in node.names:
                parts = a.name.split(".")
                if a.asname is None and len(parts) > 1 and parts[0] == package:
                    names.add(parts[1])
        elif isinstance(node, ast.ImportFrom):
            if not any(a.name == "*" for a in node.names):
                continue
            target = (node.module or "").split(".")[-1]
            source = path.parent / f"{target}.py"
            if not source.is_file():
                problems.append(
                    f"{path.relative_to(ROOT)}:{node.lineno} `from {node.module} "
                    f"import *` cannot be resolved, so the names it binds are "
                    f"unknown and this file cannot be checked")
                continue
            try:
                names |= _module_level_bindings(symtable.symtable(
                    source.read_text(encoding="utf-8"), str(source), "exec"))
            except SyntaxError as exc:
                problems.append(f"{source.relative_to(ROOT)}: {exc}")
    return names, problems


def _scope_nodes(tree: ast.AST) -> dict[tuple[str, int], ast.AST]:
    """(symtable scope name, line) -> the ast node that made the scope."""
    out: dict[tuple[str, int], ast.AST] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            out[(node.name, node.lineno)] = node
        elif isinstance(node, ast.Lambda):
            out[("lambda", node.lineno)] = node
        elif isinstance(node, (ast.GeneratorExp, ast.ListComp,
                               ast.SetComp, ast.DictComp)):
            for kind in ("genexpr", "listcomp", "setcomp", "dictcomp"):
                out.setdefault((kind, node.lineno), node)
    return out


def check_no_unbound_names(files: list[Path]) -> list[str]:
    """No name is loaded that nothing in an enclosing scope binds.

    A dead-code sweep has to be run with the language's own resolver. `world-ro6`
    swept the model Fortran for unreferenced procedures and then applied the same
    sweep to `pyburn.py`, where a Fortran call graph cannot see a Python caller:
    it removed `readallvariables`, which `readfile` calls, and the next run
    integrated ten orbits and died in the postprocessor with NameError. Four of
    the five names it removed really were dead, which is why it looked right.

    So this check asks CPython's own symbol table, not a pattern: for every
    scope, a name that is referenced and resolves to a module global that nothing
    binds is a name that can only ever raise NameError. That is decidable and it
    is what a sweep must not be able to create. It found five in `pyburn.py`
    besides the `readallvariables` break, three of them misspellings
    (`vairable`, `gpvuar`, `rottlgridvar`) of names in the same scope -- which
    means those branches had never executed here or upstream.

    It is deliberately NARROWER than pyflakes F821, which also reports a name
    read before its binding in the same scope. That is a real smell but it is
    not always a bug: `pyburn.dataset()` carries `theta` from one iteration of
    its variable loop into the next, and it works. This check has one answer and
    can only fail when the name is unreachable, so it can be pointed at vendored
    source that this project does not otherwise hold to its own style.
    """
    problems: list[str] = []
    for path in files:
        try:
            src = path.read_text(encoding="utf-8")
            tree = ast.parse(src, str(path))
            top = symtable.symtable(src, str(path), "exec")
        except SyntaxError as exc:
            problems.append(f"{path.relative_to(ROOT)}: {exc}")
            continue
        bound = _module_level_bindings(top) | MODULE_DUNDERS
        implicit, unresolved = _implicit_bindings(path, tree)
        problems.extend(unresolved)
        if unresolved:
            continue
        bound |= implicit
        # `global x` in any scope, with an assignment there, binds x at module
        # scope even though no module-level statement does.
        stack = [top]
        while stack:
            scope = stack.pop()
            stack.extend(scope.get_children())
            bound |= {s.get_name() for s in scope.get_symbols()
                      if s.is_declared_global()}
        nodes = _scope_nodes(tree)
        found: set[tuple[int, str, str]] = set()
        stack = [top]
        while stack:
            scope = stack.pop()
            if scope.get_type() in _UNREPORTED_SCOPES:
                continue
            stack.extend(scope.get_children())
            unbound = set()
            for sym in scope.get_symbols():
                name = sym.get_name()
                if not sym.is_referenced() or hasattr(builtins, name):
                    continue
                if name in bound:
                    continue
                if scope.get_type() == "module":
                    if not (sym.is_assigned() or sym.is_imported()):
                        unbound.add(name)
                elif sym.is_global():
                    unbound.add(name)
            if not unbound:
                continue
            node = tree if scope.get_type() == "module" else \
                nodes.get((scope.get_name(), scope.get_lineno()))
            if node is None:
                found.add((scope.get_lineno(), ", ".join(sorted(unbound)),
                           scope.get_name()))
                continue
            for sub in ast.walk(node):
                if isinstance(sub, ast.Name) and isinstance(sub.ctx, ast.Load) \
                        and sub.id in unbound:
                    found.add((sub.lineno, sub.id, scope.get_name()))
        for lineno, name, scope_name in sorted(found):
            problems.append(f"{path.relative_to(ROOT)}:{lineno} loads `{name}`, "
                            f"which nothing binds; in {scope_name}()")
    return problems
def check_tail_fit_stops_above_roundoff() -> list[str]:
    """`spectral_tail.py`'s tail slope is the slope it was given, not the floor's.

    A synthetic T42 spectrum: an inertial range, a tail of KNOWN slope below the
    fit band, and every wavenumber sitting on a floating-point roundoff floor
    with the dead band above the truncation carrying only that floor. The right
    answer is the slope put in, which is what makes this a test rather than a
    comparison -- and it is one the script's earlier shape cannot pass, because
    it cut the spectrum at `m = NTRU` and fitted straight through the floor.

    The tolerance is 0.1 in the slope against an effect of up to 17: fitting to
    the truncation returns about -23 whatever it is given below that, since the
    floor anchors the bottom of the fit and the true tail cannot be seen through
    it. So the check has margin of two orders over its own resolution, and a
    reimplementation that quietly extends the band again fails it by a wide one.
    """
    import numpy as np
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    try:
        from spectral_tail import FLOOR_MARGIN, slope
    except ImportError as exc:
        return [f"exoplasim/scripts/spectral_tail.py does not import: {exc}"]
    n, nlon, floor_value, tol = 42, 128, 2.6e-15, 0.1
    m = np.arange(nlon // 2 + 1)
    hi = max(3, n // 3)
    rng = np.random.default_rng(0)
    bad = []
    for true_tail in (-25.0, -30.0, -35.0, -40.0):
        mm = np.maximum(m, 1).astype(float)
        inertial = 1e-1 * mm ** -2.04
        tail_law = inertial[hi] * (mm / hi) ** true_tail
        full = np.maximum(np.where(m <= hi, inertial, tail_law),
                          floor_value * rng.uniform(0.4, 1.6, len(m)))
        full[0] = 0.0
        spec = full[:n + 1]
        floor = float(np.median(full[n + 1:]))
        live = [int(k) for k in np.arange(len(spec))[1:]
                if spec[k] > FLOOR_MARGIN * floor]
        got = slope(spec, hi, min(n, max(live)))
        if not abs(got - true_tail) <= tol:
            bad.append(f"a tail of {true_tail:+g} fitted above the roundoff "
                       f"floor came back as {got:+.2f}")
    return bad

def check_slab_capacity_follows_the_run() -> list[str]:
    """The convergence slab takes its sea water from the RUN, not from the model.

    `CRHOS` and `CPS` are `icemod_nl` keys and every run this project writes
    declares them, so a run integrated before a change to `icemod.f90` must be
    assessed at the capacity it used. Read once at import, the pair silently
    re-verdicts every run older than the change: the capacity sets
    `tau_expected`, which sets the fallback remaining offset, which is what
    `OFFSET_TOLERANCE_K` fails a run on. The model's own CPS has already moved
    by 4.54 per cent once.

    Two synthetic run directories, one overriding the pair and one declaring
    nothing. The right answers are known in closed form -- the override's own
    product, and the compiled declaration -- and a module-level constant cannot
    give both, which is what makes this a test rather than a comparison.
    """
    import tempfile
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import assess_convergence
        import sea_water
    except ImportError as exc:
        return [f"assess_convergence.py or lib/sea_water.py does not import: {exc}"]
    bad = []
    compiled = sea_water.constants()
    mld = assess_convergence._MLD
    with tempfile.TemporaryDirectory() as tmp:
        # A run that declares its own sea water, as run_exoplasim.py writes it.
        override = Path(tmp) / "with_namelist"
        override.mkdir()
        crhos, cps = 1011.0, 3777.0
        (override / "icemod_namelist").write_text(
            f" &icemod_nl\n CPS = {cps}\n CRHOS = {crhos}\n /END\n", encoding="utf-8")
        got, water = assess_convergence.slab_heat_capacity(override)
        want = mld * crhos * cps
        if abs(got - want) > 1e-6 * want:
            bad.append(f"a run declaring CRHOS={crhos} CPS={cps} was assessed at "
                       f"{got:.6e} J/m2/K, not its own {want:.6e}")
        if water.get("CPS") != cps:
            bad.append(f"the assessment recorded CPS={water.get('CPS')}, not the "
                       f"run's {cps}")
        # A run that declares nothing falls back to the compiled model.
        bare = Path(tmp) / "no_namelist"
        bare.mkdir()
        got_bare, _ = assess_convergence.slab_heat_capacity(bare)
        want_bare = mld * compiled["CRHOS"] * compiled["CPS"]
        if abs(got_bare - want_bare) > 1e-6 * want_bare:
            bad.append(f"a run declaring no sea water was assessed at {got_bare:.6e} "
                       f"J/m2/K, not the compiled {want_bare:.6e}")
        # The two must differ, or the check above proves nothing.
        if abs(got - got_bare) <= 1e-6 * want_bare:
            bad.append("the override and the compiled default gave the same "
                       "capacity, so this check cannot see the difference it exists "
                       "to catch")
    return bad


def check_melting_point_follows_the_run() -> list[str]:
    """The state closure takes the melting point from the RUN, not from a literal.

    `tmelt` is pumamod's and a `planet_nl` key, so a run can set it and
    `icemod` takes it through `iceini` rather than holding a compile-time copy.
    It weights the sea-ice part of the rebuilt mixed-layer heat content, so a
    literal that no longer tracks the model rescales the storage term the
    convergence criterion passes on, and rescales it without failing.

    Two synthetic run directories with known answers in closed form -- the
    override's own value and the planet module's `planet_ini` assignment -- plus
    the absence of a module-level copy, which is the form the defect took.
    """
    import tempfile
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import close_state_energy
        import sea_water
    except ImportError as exc:
        return [f"close_state_energy.py or lib/sea_water.py does not import: {exc}"]
    bad = []
    if hasattr(close_state_energy, "TMELT"):
        bad.append("close_state_energy.py carries a module-level TMELT again; the "
                   "melting point is sea_water.melting_point(run_dir)'s to state")
    compiled = sea_water.melting_point()["TMELT"]
    with tempfile.TemporaryDirectory() as tmp:
        override = Path(tmp) / "with_namelist"
        override.mkdir()
        tmelt = 260.5
        (override / "planet_namelist").write_text(
            f" &planet_nl\n GA = 12.81\n TMELT = {tmelt}\n /END\n", encoding="utf-8")
        got = sea_water.melting_point(override)["TMELT"]
        if got != tmelt:
            bad.append(f"a run declaring TMELT={tmelt} was closed at {got}")
        bare = Path(tmp) / "no_namelist"
        bare.mkdir()
        got_bare = sea_water.melting_point(bare)["TMELT"]
        if got_bare != compiled:
            bad.append(f"a run declaring no melting point was closed at {got_bare}, "
                       f"not the planet module's {compiled}")
        if got == got_bare:
            bad.append("the override and the planet module gave the same melting "
                       "point, so this check cannot see the difference it exists "
                       "to catch")
    # It must be the melting point and not the sea-water freezing point: the
    # stale comment this check replaces named the wrong one of the two.
    if abs(compiled - sea_water.constants()["TFREEZE"]) < 1.0:
        bad.append("sea_water.melting_point() is returning something within a "
                   "kelvin of TFREEZE; those are two different quantities")
    return bad


def check_cloud_tables_match_the_papers() -> list[str]:
    """The shortwave cloud tables in `radmod.f90` are the ones the papers print.

    `swr`'s computed-cloud branch reads Stephens et al. (1984) Tables 1(a),
    1(b) and 1(c) by bilinear interpolation, 324 entries of Fortran source that
    nothing else in the tree constrains. A slipped digit in any of them is a
    silently different cloud: no compile fails, no run refuses, and the model
    reflects a little more or a little less near infrared for the rest of its
    life. `exoplasim/scripts/stephens_tables_vs_fits.py` holds an independent
    transcription of both papers, already checked against the papers' own
    claims about which of their tables the 1984 revision altered, and reads the
    model source back against it.

    It also checks the conservative floor `swr` applies to the co-albedo, which
    is the one place the table can be exactly zero and the two-stream solution
    divides by it. Both are pure reads: no build, no run, no artifact.
    """
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import numpy as np
        import stephens_tables_vs_fits as tables
    except ImportError as exc:
        return [f"exoplasim/scripts/stephens_tables_vs_fits.py does not "
                f"import: {exc}"]
    bad = []
    for what, run in (
            ("the transcription of the two papers", tables.check_transcription),
            ("radmod's port of the 1984 tables", tables.check_model_source),
            ("the co-albedo floor against Stephens Eq. (1)",
             lambda: tables.check_conservative_floor(
                 math.sqrt(np.finfo(np.float64).eps)))):
        try:
            run()
        except SystemExit as exc:
            bad.append(f"{what}: {exc}")
    return bad


def check_run_length_derivation() -> list[str]:
    """A declared run length is derived from a timescale, and from the right one.

    `lib/run_lengths.py` turns two timescales into two lengths, and the whole
    value of it is that the two are kept apart: a commissioning span is bought
    in the memory time and a settling block in the relaxation time. This checks
    the arithmetic against answers known in closed form, checks that the two
    tolerances that must be one number ARE one number, and checks that the
    derived lengths still reproduce the operational experience they were built
    to agree with -- because a derivation that no longer lands where experience
    does is either wrong or has found something, and either way it must not pass
    silently.
    """
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import assess_convergence
        import run_lengths
    except ImportError as exc:
        return [f"lib/run_lengths.py or assess_convergence.py does not import: {exc}"]
    bad = []
    if run_lengths.SETTLING_RESIDUAL_K != assess_convergence._OFFSET_TOLERANCE_K:
        bad.append(
            f"a settling block decays to {run_lengths.SETTLING_RESIDUAL_K} K while "
            f"the offset criterion allows {assess_convergence._OFFSET_TOLERANCE_K} K; "
            "they are one number and a settling residual the next verdict can see "
            "is not settled")
    # Closed form: n = tau * ln(A0 / residual). At A0 = e * residual it is
    # exactly tau, which is an identity rather than a comparison.
    tau = 7.0
    want = tau
    got = run_lengths.settling_orbits(math.e * run_lengths.SETTLING_RESIDUAL_K, tau)
    if abs(got - want) > 1e-9:
        bad.append(f"a perturbation e times the residual should settle in one tau, "
                   f"{want}; got {got}")
    # A perturbation already below the residual has nothing to decay.
    if run_lengths.settling_orbits(0.5 * run_lengths.SETTLING_RESIDUAL_K, tau) != 1.0:
        bad.append("a perturbation smaller than the residual was given more than "
                   "one orbit to decay")
    # Closed form: the span is the multiple, and the commissioning length adds
    # the approach to it rather than replacing it.
    if run_lengths.production_span_orbits(3.0) != 3.0 * run_lengths.PRODUCTION_SPAN_TAU_MULTIPLE:
        bad.append("the production span is not the declared multiple of tau")
    if run_lengths.commissioning_orbits(11.0, 3.0) != 11.0 + run_lengths.production_span_orbits(3.0):
        bad.append("a commissioning length is not its approach plus its span")
    # The relaxation bracket is READ rather than declared, so its absence is a
    # refusal from the module itself. A gate reports that as a failed check and
    # not as a traceback: a check that dies looks like a broken gate rather than
    # a tree that is missing an artifact it needs.
    try:
        relaxation = run_lengths.tau_relaxation_orbits_bracket(ROOT)
    except RuntimeError as exc:
        bad.append(str(exc))
        relaxation = None
    for name, bracket in (("memory", run_lengths.TAU_MEMORY_ORBITS_BRACKET),
                          ("relaxation", relaxation)):
        if bracket is None:
            continue
        low, high = bracket
        if not 0 < low <= high:
            bad.append(f"the {name}-time bracket {bracket} is not an ordered "
                       "pair of positive times")
    # THE DERIVED LENGTHS MUST STILL LAND ON THE EXPERIENCE. A reconvergence has
    # repeatedly taken ten to twenty orbits after a step change worth about half
    # a kelvin, and the settling bracket for that perturbation has to cover it.
    if relaxation is not None:
        low, high = run_lengths.settling_bracket(0.5, root=ROOT)
        if not (low <= 20.0 and high >= 10.0):
            bad.append(f"a 0.5 K step settles in {low:.1f} to {high:.1f} orbits "
                       "by derivation, which no longer overlaps the ten to "
                       "twenty this project has repeatedly seen")
    return bad


def _convergence_report(scatter: float, tau: float, relaxation: float,
                        purpose: str = "production",
                        tau_resolved: bool = True) -> dict:
    """The parts of a convergence report the two bound checks read."""
    return {
        "assessed_purpose": purpose,
        "resolving_power": {
            "temperature_residual_scatter_k": scatter,
            "temperature_residual_tau_orbits": tau,
            "tau_span_supports_the_estimate": tau_resolved,
        },
        "metrics": {"relaxation_orbits_expected": relaxation},
    }


def _convergence_tree(tmp: str, reports: dict) -> Path:
    """A repository root carrying nothing but convergence reports."""
    import json as _json
    root = Path(tmp)
    directory = root / "exoplasim" / "analysis" / "convergence"
    directory.mkdir(parents=True, exist_ok=True)
    for name, report in reports.items():
        (directory / name).write_text(_json.dumps(report), encoding="utf-8")
    return root


def check_convergence_bounds_are_re_read() -> list[str]:
    """The window's nominal inputs are held to the reports they must bound.

    `notes/audits/frozen-derived-quantities.md` tier 1: the scatter and the
    memory time size every commissioning run on the ladder, and nothing carried
    a measurement back to them. They are UPPER BOUNDS on a series that still
    drifts, so a report reading BELOW them is the expected state and is not the
    defect -- the defect was that nothing could tell that state from a bound
    nobody had looked at since the series moved.

    Class 17: every case has a right answer, and every refusal is paired with
    the positive that proves the setup was real. The anchor cases are the ones
    the audit asks for, because they are what separates "still holds" from
    "never re-examined"; the bound cases are what a real falsification looks
    like. Synthetic reports in a temp directory; the live tree is checked last.
    """
    import tempfile
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import assess_convergence
        import run_lengths
    except ImportError as exc:
        return [f"assess_convergence.py or lib/run_lengths.py does not import: {exc}"]

    bad = []
    scatter_anchor_name = Path(
        assess_convergence.CONVERGENCE_BOUND_ANCHOR["artifact"]).name
    memory_anchor_name = Path(
        run_lengths.MEMORY_BRACKET_ANCHOR["artifact"]).name
    anchor_scatter = assess_convergence.CONVERGENCE_BOUND_ANCHOR["orbit_scatter_k"]
    anchor_tau = run_lengths.MEMORY_BRACKET_ANCHOR["observation"]
    tau_top = max(run_lengths.TAU_MEMORY_ORBITS_BRACKET)
    scatter_top = assess_convergence.NOMINAL_ORBIT_SCATTER_K
    relaxation_top = assess_convergence.NOMINAL_RELAXATION_ORBITS

    def anchor(scatter=anchor_scatter, tau=anchor_tau, **overrides):
        report = _convergence_report(scatter, tau, relaxation_top - 1.0)
        report["resolving_power"].update(overrides)
        return report

    def anchors(scatter_overrides=None, memory_overrides=None):
        """Both independent evidence files, even when their runs differ."""
        reports = {
            scatter_anchor_name: anchor(
                tau=tau_top * 0.5, **(scatter_overrides or {})),
            memory_anchor_name: anchor(
                scatter=scatter_top * 0.5, **(memory_overrides or {})),
        }
        return reports

    def case(name, root, want_scatter_problem, want_tau_problem):
        got = assess_convergence.check_convergence_bounds(root)
        if bool(got) != bool(want_scatter_problem):
            bad.append(f"{name}: check_convergence_bounds returned {got!r}, "
                       f"expected {'a refusal' if want_scatter_problem else 'no problem'}")
        got = run_lengths.check_memory_bracket(root)
        if bool(got) != bool(want_tau_problem):
            bad.append(f"{name}: check_memory_bracket returned {got!r}, "
                       f"expected {'a refusal' if want_tau_problem else 'no problem'}")

    # THE POSITIVE THAT PROVES THE SETUP. The anchor reads exactly what was
    # recorded and every bound has headroom, so the refusals below are caused by
    # what each case changes and not by the fixture merely existing.
    with tempfile.TemporaryDirectory() as tmp:
        case("an anchor reading what was recorded passes",
             _convergence_tree(tmp, anchors()), False, False)

    # THE ANCHOR IS GONE: no report at all. Both bounds are unexamined, and an
    # unexamined bound is not a held one.
    with tempfile.TemporaryDirectory() as tmp:
        case("a missing anchor is a refusal",
             _convergence_tree(tmp, {}), True, True)

    # THE ANCHOR HAS MOVED: the run was extended or re-assessed. This is the arm
    # that closes the circle the audit named -- a longer series changes the
    # reading, and the reading refuses the declaration derived from the old one.
    with tempfile.TemporaryDirectory() as tmp:
        case("a moved anchor scatter is a refusal",
             _convergence_tree(tmp, anchors(scatter_overrides={
                 "temperature_residual_scatter_k": anchor_scatter * 0.5})),
             True, False)
    with tempfile.TemporaryDirectory() as tmp:
        case("a moved anchor memory time is a refusal",
             _convergence_tree(tmp, anchors(memory_overrides={
                 "temperature_residual_tau_orbits": anchor_tau * 0.5})),
             False, True)

    # A READING BELOW THE BOUND IS NOT A DEFECT, which is the case the audit is
    # explicit about: both declared values are upper bounds on a drifting
    # series, so a smaller measurement is what they predict.
    with tempfile.TemporaryDirectory() as tmp:
        case("a second report reading well below both bounds passes",
             _convergence_tree(tmp, anchors() | {
                 "run_0000_convergence.json": _convergence_report(
                     scatter_top * 0.1, tau_top * 0.1, 1.0)}),
             False, False)

    # THE BOUND IS VIOLATED. A settled production report above either bound is
    # the bound being wrong rather than old.
    with tempfile.TemporaryDirectory() as tmp:
        case("a settled report above the scatter bound is a refusal",
             _convergence_tree(tmp, anchors() | {
                 "run_0000_convergence.json": _convergence_report(
                     scatter_top * 1.1, tau_top * 0.5, 1.0)}),
             True, False)
    with tempfile.TemporaryDirectory() as tmp:
        case("a settled report above the memory bound is a refusal",
             _convergence_tree(tmp, anchors() | {
                 "run_0000_convergence.json": _convergence_report(
                     scatter_top * 0.5, tau_top * 1.1, 1.0)}),
             False, True)
    with tempfile.TemporaryDirectory() as tmp:
        case("a report deriving more relaxation than the window is sized on "
             "is a refusal",
             _convergence_tree(tmp, anchors() | {
                 "run_0000_convergence.json": _convergence_report(
                     scatter_top * 0.5, tau_top * 0.5, relaxation_top * 1.1)}),
             True, False)

    # THE FILTERS ARE REAL, and each is proved by the same numbers passing on one
    # side of it and refusing on the other. A diagnostic report prices an A/B
    # arm's perturbed orbits; a report whose window cannot fix its own tau says
    # so. Neither is a reading of the settled trajectory.
    for label, kwargs in (("a diagnostic report", {"purpose": "diagnostic"}),
                          ("a report whose window cannot fix its own tau",
                           {"tau_resolved": False})):
        with tempfile.TemporaryDirectory() as tmp:
            case(f"{label} above both bounds is not a refusal",
                 _convergence_tree(tmp, anchors() | {
                     "run_0000_convergence.json": _convergence_report(
                         scatter_top * 2.0, tau_top * 2.0, 1.0, **kwargs)}),
                 False, False)
        with tempfile.TemporaryDirectory() as tmp:
            case(f"the same numbers in a settled production report ARE a "
                 f"refusal, against {label}",
                 _convergence_tree(tmp, anchors() | {
                     "run_0000_convergence.json": _convergence_report(
                         scatter_top * 2.0, tau_top * 2.0, 1.0)}),
                 True, True)

    # A REPORT WITH NO PURPOSE FIELD predates the field. It is production unless
    # its filename says otherwise, and both readings are proved here.
    with tempfile.TemporaryDirectory() as tmp:
        old = _convergence_report(scatter_top * 2.0, tau_top * 2.0, 1.0)
        del old["assessed_purpose"]
        case("an old report with no purpose field is read as production",
             _convergence_tree(tmp, anchors() | {
                 "run_0000_convergence.json": old}),
             True, True)
    with tempfile.TemporaryDirectory() as tmp:
        old = _convergence_report(scatter_top * 2.0, tau_top * 2.0, 1.0)
        del old["assessed_purpose"]
        case("an old report named diagnostic is read as one",
             _convergence_tree(tmp, anchors() | {
                 "run_0000_convergence_diagnostic.json": old}),
             False, False)

    # AND THE LIVE TREE. Everything above proves the instrument; this is what it
    # says about the declarations actually in force.
    bad.extend(assess_convergence.check_convergence_bounds(ROOT))
    bad.extend(run_lengths.check_memory_bracket(ROOT))
    return bad


def check_relaxation_bracket_is_read() -> list[str]:
    """The relaxation-time bracket is READ from its producer, never declared.

    `notes/audits/frozen-derived-quantities.md` tier 1 calls this the cleanest
    instance in the audit: the producer existed, the consumer existed, and the
    wire between them was a human. So the test is that there is no number to
    freeze -- `lib/run_lengths.py` refuses rather than falling back when the
    artifact is absent or malformed, and takes whatever a present one says.
    """
    import json as _json
    import tempfile
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import run_lengths
    except ImportError as exc:
        return [f"lib/run_lengths.py does not import: {exc}"]

    bad = []
    relative = run_lengths.RELAXATION_CEILING_ARTIFACT

    def tree(tmp: str, payload) -> Path:
        root = Path(tmp)
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if payload is not None:
            path.write_text(_json.dumps(payload), encoding="utf-8")
        return root

    def refuses(name, root) -> None:
        try:
            run_lengths.tau_relaxation_orbits_bracket(root)
        except RuntimeError:
            return
        bad.append(f"{name}: expected a refusal and got a bracket")

    # The positive: an artifact carrying a bracket IS the bracket, to the digit.
    with tempfile.TemporaryDirectory() as tmp:
        root = tree(tmp, {"fitted_bracket_orbits": [3.25, 17.5]})
        got = run_lengths.tau_relaxation_orbits_bracket(root)
        if got != (3.25, 17.5):
            bad.append(f"the emitted bracket [3.25, 17.5] was read as {got!r}")
        settling = run_lengths.settling_bracket(math.e * run_lengths.SETTLING_RESIDUAL_K,
                                                root=root)
        # Closed form: at a perturbation e times the residual a settling block
        # is exactly one tau, so the settling bracket IS the tau bracket.
        if settling != (3.25, 17.5):
            bad.append("a settling block for a perturbation of e times the "
                       f"residual should be the tau bracket itself; got {settling!r}")

    with tempfile.TemporaryDirectory() as tmp:
        refuses("no artifact is a refusal, not a default", tree(tmp, None))
    with tempfile.TemporaryDirectory() as tmp:
        refuses("an artifact with no bracket is a refusal",
                tree(tmp, {"verdict": "untested", "fitted_bracket_orbits": None}))
    with tempfile.TemporaryDirectory() as tmp:
        refuses("a disordered bracket is a refusal",
                tree(tmp, {"fitted_bracket_orbits": [17.5, 3.25]}))
    with tempfile.TemporaryDirectory() as tmp:
        refuses("a non-positive bracket is a refusal",
                tree(tmp, {"fitted_bracket_orbits": [0.0, 3.25]}))

    # AND THE PRODUCER STILL EMITS ONE THIS TREE CAN READ. A bracket nobody can
    # read is the frozen state with the number deleted, so the absence of the
    # artifact is reported as a failed check rather than raised through the gate.
    try:
        low, high = run_lengths.tau_relaxation_orbits_bracket(ROOT)
    except RuntimeError as exc:
        bad.append(str(exc))
    else:
        if not 0 < low <= high:
            bad.append(f"this tree emits a relaxation bracket of ({low}, {high})")
    return bad


def check_fit_tail_fraction_is_stated_once() -> list[str]:
    """The exponential fit's tail fraction is one statement, not three.

    THE TEST THAT CAN FAIL: move the constant, and every consumer must move with
    it. A copy is exactly what would not -- that is what a copy is -- so the
    span the assessment records, the span `check_relaxation_ceiling.py`
    reconstructs for an artifact too old to carry one, and the fit's own mask
    are all driven from the one name and read back.
    """
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    try:
        import assess_convergence
        import check_relaxation_ceiling
    except ImportError as exc:
        return [f"assess_convergence.py or check_relaxation_ceiling.py does "
                f"not import: {exc}"]

    bad = []
    if check_relaxation_ceiling.FIT_TAIL_FRACTION is not assess_convergence.FIT_TAIL_FRACTION:
        bad.append("check_relaxation_ceiling.py carries its own tail fraction "
                   "rather than the assessment's")
    # A default is bound once, at definition, so a constant spelled as one is a
    # snapshot that cannot follow the name it was copied from.
    default = inspect.signature(
        assess_convergence.approach_to_equilibrium).parameters["fraction"].default
    if default is not None:
        bad.append(f"approach_to_equilibrium's fraction defaults to {default!r}, "
                   "which snapshots the constant at definition")

    original = assess_convergence.FIT_TAIL_FRACTION
    try:
        # Closed form at both settings: orbits are 0..n-1, the mask is
        # `orbit >= (n-1) * fraction`, so a hundred orbits keep 65 at 0.35 and
        # 50 at 0.5. Both are counted rather than one, so a consumer frozen at
        # the original value fails the second and not the first.
        for fraction, want in ((0.35, 65), (0.5, 50)):
            assess_convergence.FIT_TAIL_FRACTION = fraction
            got = assess_convergence.fit_span_orbits(100)
            if got != want:
                bad.append(f"at a tail fraction of {fraction} the fit sees "
                           f"{got} of 100 orbits, not {want}")
            got = check_relaxation_ceiling.recorded_or_reconstructed_span(
                {"completed_orbits": 100})
            if got != want:
                bad.append(f"at a tail fraction of {fraction} "
                           f"check_relaxation_ceiling.py reconstructs a span of "
                           f"{got} orbits, not {want}; it is not reading the "
                           "assessment's constant")
    finally:
        assess_convergence.FIT_TAIL_FRACTION = original
    # A RECORDED SPAN IS TAKEN AS RECORDED and not reconstructed, which is what
    # lets an artifact written at another fraction still be read correctly.
    got = check_relaxation_ceiling.recorded_or_reconstructed_span(
        {"completed_orbits": 100, "metrics": {"relaxation_fit_span_orbits": 7}})
    if got != 7:
        bad.append(f"a recorded fit span of 7 orbits was reconstructed as {got}")
    return bad


def check_commissioning_evidence_is_re_read() -> list[str]:
    """Every commissioning row is re-read from the run record it was copied from.

    `notes/audits/frozen-derived-quantities.md` tier 1: `_check_ceilings` and
    `_check_route` validate the rungs and the ordering, and nothing re-read an
    orbit count or a verdict. The counts permit a route step, so a phantom row
    licenses an advance on evidence nothing can check.

    Class 17: the record is the right answer and the table is the copy. Every
    refusal below is paired with the positive that proves the fixture was real,
    and the table is substituted rather than the check reimplemented, so what is
    tested is the function the gate calls.
    """
    import json as _json
    import tempfile
    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import rungs
    except ImportError as exc:
        return [f"lib/rungs.py does not import: {exc}"]

    bad = []

    def tree(tmp: str, entries: list) -> Path:
        root = Path(tmp)
        path = root / "exoplasim" / "runs" / "INDEX.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_json.dumps({"runs": entries}), encoding="utf-8")
        return root

    def record(run, orbits, resolution="T21", status="quasi_equilibrated"):
        return {"run_id": run, "orbits_on_disk": orbits, "status": status,
                "physical": {"resolution": resolution}}

    def case(name, table, entries, want_problem) -> None:
        original = rungs.COMMISSIONING_EVIDENCE
        with tempfile.TemporaryDirectory() as tmp:
            rungs.COMMISSIONING_EVIDENCE = table
            try:
                got = rungs.check_commissioning_evidence(tree(tmp, entries))
            finally:
                rungs.COMMISSIONING_EVIDENCE = original
        if bool(got) != bool(want_problem):
            bad.append(f"{name}: got {got!r}, expected "
                       f"{'a refusal' if want_problem else 'no problem'}")

    row = {"verdict": "endured", "orbits": 50, "run": "run_aaaaaaaaaaaa"}
    table = {("T21", 45.0): row}

    case("a row matching its record passes",
         table, [record("run_aaaaaaaaaaaa", 50)], False)
    case("an orbit count the record contradicts is a refusal",
         table, [record("run_aaaaaaaaaaaa", 49)], True)
    case("a row filed under a rung the record contradicts is a refusal",
         table, [record("run_aaaaaaaaaaaa", 50, resolution="T42")], True)
    case("endurance on a run the record calls failed is a refusal",
         table, [record("run_aaaaaaaaaaaa", 50, status="failed")], True)
    case("a row whose run is in no record and does not say so is a refusal",
         table, [record("run_zzzzzzzzzzzz", 50)], True)
    case("the same row saying its record is gone passes",
         {("T21", 45.0): dict(row, record_gone="searched every index; absent")},
         [record("run_zzzzzzzzzzzz", 50)], False)
    case("a row claiming its record is gone while it is there is a refusal",
         {("T21", 45.0): dict(row, record_gone="searched every index; absent")},
         [record("run_aaaaaaaaaaaa", 50)], True)
    case("a row that names no run is a refusal",
         {("T21", 45.0): {"verdict": "endured", "orbits": 50}},
         [record("run_aaaaaaaaaaaa", 50)], True)

    blew = {"verdict": "blew_up", "orbits": 46, "run": "run_aaaaaaaaaaaa"}
    case("a blow-up on a run the record calls failed passes",
         {("T21", 45.0): blew},
         [record("run_aaaaaaaaaaaa", 46, status="failed")], False)
    case("a blow-up on a run the record calls equilibrated is a refusal",
         {("T21", 45.0): blew},
         [record("run_aaaaaaaaaaaa", 46, status="quasi_equilibrated")], True)
    case("a verdict that is neither of the two is a refusal",
         {("T21", 45.0): dict(row, verdict="settled")},
         [record("run_aaaaaaaaaaaa", 50)], True)

    # A TREE WITH NO RUN RECORDS AT ALL is not a disagreement, the same rule
    # `check_stability_ceilings` applies to an absent probe grid.
    original = rungs.COMMISSIONING_EVIDENCE
    with tempfile.TemporaryDirectory() as tmp:
        rungs.COMMISSIONING_EVIDENCE = table
        try:
            got = rungs.check_commissioning_evidence(Path(tmp))
        finally:
            rungs.COMMISSIONING_EVIDENCE = original
    if got:
        bad.append(f"a tree with no run records was read as a disagreement: {got!r}")

    bad.extend(rungs.check_commissioning_evidence(ROOT))
    return bad


# Every bracket pedology declares, and what each one is a bracket ON. The table
# is written out rather than discovered by key name, for the reason
# `pedology/scripts/outgassing_gate.py` gives about its own: a check that
# guesses which keys are brackets stops checking the moment someone names one
# differently, and stops silently. `_unregistered_pedology_brackets` below is
# the other half -- it walks the files and reports a bracket this table does not
# name -- so the two cannot drift apart.
#
# Four dispositions, and they are not interchangeable:
#
#   value    the file states the bracket and the number inside it. The number
#            has to be inside the bracket.
#   pair     the same, with the two ends carried as `<key>_minimum` and
#            `<key>_maximum` beside the value rather than as a two-element list.
#   level    the bracket is on a quantity NO key in the file states, because it
#            is a land mean over a build's own soil map. Nothing static can
#            evaluate it, so what is checked here is that the named generator
#            still READS the key: the run-time refusal is the enforcement, and a
#            refusal someone deletes is what this catches.
#   derived  the bracket is the sentinel `derived`, and the named module fills
#            it from `config/planet.yaml`. The sentinel has to be intact, the
#            module has to still name the key, and the value it brackets is
#            checked against the bracket that module computes.
#
# (kind, label, config, value path, bracket path, enforcer)
PEDOLOGY_CONFIGS = {
    "pedogenesis": Path("pedology/config/pedogenesis.yaml"),
    "surface_classes": Path("pedology/config/surface_classes.yaml"),
}
PEDOLOGY_BRACKETS = (
    ("value", "texture.clay_conversion", "pedogenesis",
     ("texture", "clay_conversion"), ("texture", "clay_conversion_bracket"), None),
    ("value", "texture.sand_to_silt_loss_ratio", "pedogenesis",
     ("texture", "sand_to_silt_loss_ratio"),
     ("texture", "sand_to_silt_loss_ratio_bracket"), None),
    ("value", "regolith.maximum_depth_m", "pedogenesis",
     ("regolith", "maximum_depth_m"), ("regolith", "maximum_depth_bracket_m"), None),
    ("value", "regolith.dry_erosion_baseline", "pedogenesis",
     ("regolith", "dry_erosion_baseline"),
     ("regolith", "dry_erosion_baseline_bracket"), None),
    ("value", "ph.parent_by_category.evaporite", "pedogenesis",
     ("ph", "parent_by_category", "evaporite"),
     ("ph", "parent_bracket_evaporite"), None),
    ("value", "ph.leaching_slope", "pedogenesis",
     ("ph", "leaching_slope"), ("ph", "leaching_slope_bracket"), None),
    ("value", "ph.gibbsite_buffer_ph", "pedogenesis",
     ("ph", "gibbsite_buffer_ph"), ("ph", "gibbsite_buffer_ph_bracket"), None),
    ("value", "water.volumetric_capacity_by_texture.sand", "pedogenesis",
     ("water", "volumetric_capacity_by_texture", "sand"),
     ("water", "volumetric_capacity_by_texture_bracket", "sand"), None),
    ("value", "water.volumetric_capacity_by_texture.silt", "pedogenesis",
     ("water", "volumetric_capacity_by_texture", "silt"),
     ("water", "volumetric_capacity_by_texture_bracket", "silt"), None),
    ("value", "water.volumetric_capacity_by_texture.clay", "pedogenesis",
     ("water", "volumetric_capacity_by_texture", "clay"),
     ("water", "volumetric_capacity_by_texture_bracket", "clay"), None),
    ("value", "water.volumetric_capacity_organic", "pedogenesis",
     ("water", "volumetric_capacity_organic"),
     ("water", "volumetric_capacity_organic_bracket"), None),
    ("value", "catena.frost_production_bonus", "pedogenesis",
     ("catena", "frost_production_bonus"),
     ("catena", "frost_production_bonus_bracket"), None),
    ("value", "catena.slope_transport", "pedogenesis",
     ("catena", "slope_transport"), ("catena", "slope_transport_bracket"), None),
    ("value", "catena.slope_fines_loss", "pedogenesis",
     ("catena", "slope_fines_loss"), ("catena", "slope_fines_loss_bracket"), None),
    ("value", "catena.maximum_fines_loss", "pedogenesis",
     ("catena", "maximum_fines_loss"), ("catena", "maximum_fines_loss_bracket"), None),
    ("pair", "andisol.volumetric_capacity_allophane", "pedogenesis",
     ("andisol", "volumetric_capacity_allophane"),
     (("andisol", "volumetric_capacity_allophane_minimum"),
      ("andisol", "volumetric_capacity_allophane_maximum")), None),
    ("value", "surface_cover.loess.deposition_g_m2_yr", "surface_classes",
     ("surface_cover", "loess", "deposition_g_m2_yr"),
     ("surface_cover", "loess", "deposition_bracket"), None),
    # The one the pair `maximum_depth_m` and `erosion_coefficient_per_relief_m`
    # is jointly constrained by. Neither key's own bracket expresses it, so a
    # sweep that moves them independently can satisfy both and violate this.
    ("level", "regolith.regolith_depth_bracket_m", "pedogenesis",
     None, ("regolith", "regolith_depth_bracket_m"),
     "pedology/scripts/build_soil.py"),
    ("derived", "ph.parent_bracket_silicate", "pedogenesis",
     None, ("ph", "parent_bracket_silicate"), "pedology/scripts/carbonate_ph.py"),
    ("derived", "ph.parent_bracket_carbonate", "pedogenesis",
     None, ("ph", "parent_bracket_carbonate"), "pedology/scripts/carbonate_ph.py"),
    ("derived", "ph.endorheic_alkalinity_bonus", "pedogenesis",
     ("ph", "endorheic_alkalinity_bonus"),
     ("ph", "endorheic_alkalinity_bonus_bracket"),
     "pedology/scripts/carbonate_ph.py"),
)


_MISSING = object()


def _dig(doc, path):
    """Follow a key path, returning a sentinel rather than raising."""
    node = doc
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return _MISSING
        node = node[key]
    return node


def _unregistered_pedology_brackets(docs: dict, table=None) -> list[str]:
    """Brackets in the files that `PEDOLOGY_BRACKETS` does not name.

    The half of the pairing that keeps the table honest. Two shapes count as a
    bracket: a leaf whose key says `bracket`, and a `<key>_minimum` /
    `<key>_maximum` pair sitting beside a `<key>` that they bound. A clip is
    neither -- `ph.minimum` and `ph.maximum` bound the OUTPUT and have no value
    of their own beside them -- so the sibling requirement is what tells the two
    apart rather than a list of exceptions.
    """
    registered = set()
    for kind, _label, config, _value, bracket, _enforcer in (
            PEDOLOGY_BRACKETS if table is None else table):
        if kind == "pair":
            registered.update((config,) + p for p in bracket)
        else:
            registered.add((config,) + bracket)

    found = []
    for config, doc in docs.items():
        def walk(node, path=()):
            if isinstance(node, dict):
                for key, child in node.items():
                    if (key.endswith("_minimum") and key[:-8] + "_maximum" in node
                            and key[:-8] in node):
                        found.append((config,) + path + (key,))
                        found.append((config,) + path + (key[:-8] + "_maximum",))
                    walk(child, path + (key,))
            elif path and "bracket" in path[-1]:
                found.append((config,) + path)
        walk(doc)
    return [f"{'.'.join(p[1:])} in {PEDOLOGY_CONFIGS[p[0]]} is a bracket that "
            f"scripts/smoke_test.py:PEDOLOGY_BRACKETS does not name, so nothing "
            f"checks it" for p in sorted(set(found)) if p not in registered]


def _pedology_bracket_problems(docs: dict, derived: dict | None,
                               table=None) -> list[str]:
    """Every declared bracket, checked against the value it is a bracket on.

    `derived` maps a derived bracket's key path to the two numbers
    `carbonate_ph.py` computes for it, or is None where that module refused --
    which is itself reported, because a bracket nothing can evaluate is the
    frozen state with the number deleted.
    """
    table = PEDOLOGY_BRACKETS if table is None else table
    bad = list(_unregistered_pedology_brackets(docs, table))

    for kind, label, config, value_path, bracket_path, enforcer in table:
        doc = docs[config]
        where = PEDOLOGY_CONFIGS[config]

        if kind == "derived":
            declared = _dig(doc, bracket_path)
            if declared != "derived":
                bad.append(f"{label}'s bracket in {where} is {declared!r} and not "
                           f"the string `derived`. It is a function of "
                           f"config/planet.yaml's pCO2_bar through {enforcer}, "
                           f"and a number here cannot learn that pCO2 moved")
                continue
            source = ROOT / enforcer
            if not source.is_file() or bracket_path[-1] not in source.read_text(
                    encoding="utf-8"):
                bad.append(f"{label}'s bracket is declared `derived` and "
                           f"{enforcer} does not name {bracket_path[-1]}, so "
                           f"nothing fills it")
                continue
            bracket = None if derived is None else derived.get(bracket_path[-1])
            if bracket is None:
                bad.append(f"{label}'s bracket is declared `derived` and this "
                           f"tree cannot evaluate it; see the refusal above")
                continue
        elif kind == "pair":
            ends = [_dig(doc, p) for p in bracket_path]
            if not all(isinstance(e, (int, float)) for e in ends):
                bad.append(f"{label} in {where} is bracketed by "
                           f"{ends!r}, which is not two numbers")
                continue
            bracket = [float(e) for e in ends]
        else:
            declared = _dig(doc, bracket_path)
            if (not isinstance(declared, (list, tuple)) or len(declared) != 2
                    or not all(isinstance(x, (int, float)) for x in declared)):
                bad.append(f"{label} in {where} is bracketed by {declared!r}, "
                           f"which is not two numbers")
                continue
            bracket = [float(x) for x in declared]

        low, high = bracket
        if not low <= high:
            bad.append(f"{label} in {where} is bracketed [{low}, {high}], whose "
                       f"low end is above its high end")
            continue

        if kind == "level":
            # Nothing here states the value, so what is checkable is that the
            # generator that CAN state it still reads the key.
            source = ROOT / enforcer
            if not source.is_file() or bracket_path[-1] not in source.read_text(
                    encoding="utf-8"):
                bad.append(f"{label} brackets a level no key in {where} states, "
                           f"and {enforcer} no longer names "
                           f"{bracket_path[-1]}, so nothing evaluates it on any "
                           f"run")
            continue

        if value_path is None:
            # A derived bracket that brackets no key of its own: what it bounds
            # is `parent_by_category`, and `carbonate_ph.resolve` refuses a
            # parent outside it on every run.
            continue

        value = _dig(doc, value_path)
        if not isinstance(value, (int, float)):
            bad.append(f"{label} in {where} is {value!r}, which is not a number "
                       f"its bracket [{low}, {high}] can contain")
            continue
        if not low <= float(value) <= high:
            bad.append(f"{label} in {where} is {value} and its own bracket is "
                       f"[{low}, {high}]. A value its bracket does not contain "
                       f"is one of the two, not both")
    return bad


def check_steering_weight_is_layer_mass() -> list[str]:
    """The aerosol steering wind's vertical weight, against answers it can miss.

    All three aerosols advect on the wind of the layers at or below a declared
    sigma, and all three used to contract that level axis with a plain mean
    while their comments said the contraction was mass-weighted. The model's
    layer thicknesses differ by a factor of four, so the two are different
    quantities and a reader had no way to know which one the number was.

    `aeolian/scripts/build_dust.py:check_steering_weights` is the check and this
    runs it: the weights sum to one, a constant column comes back as that
    constant, `levp` is refused, and -- the control -- the layer-mass weights
    DIFFER from a plain mean on levels of unequal thickness, which is what says
    the first two are looking at the weight rather than past it.

    Every component keeps its own `_paths`, and an earlier check in this file
    has already imported one of them, so a plain import here would hand
    `build_dust` the wrong component's module. The cached entry is set aside
    for the length of this check and put back, which keeps this a static read
    rather than the subprocess the alternative would cost.
    """
    sys.path.insert(0, str(ROOT / "lib"))
    sys.path.insert(0, str(ROOT / "aeolian" / "scripts"))
    shadowed = sys.modules.pop("_paths", None)
    try:
        import build_dust
        return build_dust.check_steering_weights()
    finally:
        sys.modules.pop("_paths", None)
        if shadowed is not None:
            sys.modules["_paths"] = shadowed


def check_pedology_values_are_inside_their_brackets() -> list[str]:
    """Every bracket pedology declares contains the value it is a bracket on.

    `world-9ctm` gave the regolith, texture and pH blocks of
    `pedology/config/pedogenesis.yaml` the brackets they lacked, and nothing read
    them: a value edited outside its own bracket, or a bracket whose ends were
    swapped, reached a soil build with nothing objecting. That is the frozen
    state CLAUDE.md names -- a declared quantity no loop re-derives and no check
    contradicts -- and the repair is the one `pedology/scripts/outgassing_gate.py`
    already runs for its own declaration, in the per-commit tier because every
    input to it is a config read.

    Three of the four dispositions have a right answer that can fail on a value
    nobody looked at: a bracket that is not two numbers, ends the wrong way
    round, and a value outside its own bracket. The fourth, `level`, brackets
    the land-mean regolith depth that `maximum_depth_m` and
    `erosion_coefficient_per_relief_m` produce jointly; no key states it, so
    what this asserts is that `build_soil.py` still reads the key and refuses on
    it at run time.
    """
    import copy
    import yaml
    docs = {}
    for name, path in PEDOLOGY_CONFIGS.items():
        full = ROOT / path
        if not full.is_file():
            return [f"{path} does not exist"]
        docs[name] = yaml.safe_load(full.read_text(encoding="utf-8"))

    # The derived brackets, from the module that owns them rather than from a
    # number written down anywhere.
    sys.path.insert(0, str(ROOT / "pedology" / "scripts"))
    sys.path.insert(0, str(ROOT / "lib"))
    derived = None
    try:
        import carbonate_ph
        planet = yaml.safe_load(
            (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
        block = carbonate_ph.derived_ph_block(
            float(planet["atmosphere"]["pCO2_bar"]))
        derived = {k: v for k, v in block.items() if isinstance(v, (list, tuple))}
    except Exception as exc:                                   # noqa: BLE001
        return [f"the derived pH brackets cannot be evaluated: {exc!r}"]

    bad = _pedology_bracket_problems(docs, derived)

    # AND THE CHECK ITSELF CAN FAIL. Each fixture is one declaration broken in
    # one way, and each asserts the sentence it has to produce. A gate whose
    # failure path is never exercised reads as a pass for the same reason an
    # unregistered one does.
    def mutated(path, value, config="pedogenesis"):
        d = {k: copy.deepcopy(v) for k, v in docs.items()}
        node = d[config]
        for key in path[:-1]:
            node = node[key]
        node[path[-1]] = value
        return d

    fixtures = (
        ("a value outside its own bracket",
         mutated(("ph", "leaching_slope"), 3.0), "does not contain"),
        ("a bracket whose ends are the wrong way round",
         mutated(("catena", "slope_transport_bracket"), [8.0, 2.0]),
         "low end is above its high end"),
        ("a bracket that is not two numbers",
         mutated(("regolith", "dry_erosion_baseline_bracket"), "wide"),
         "not two numbers"),
        ("a min/max pair its value sits outside",
         mutated(("andisol", "volumetric_capacity_allophane"), 0.9),
         "does not contain"),
        ("a derived bracket restated as a number",
         mutated(("ph", "parent_bracket_silicate"), [5.5, 8.2]),
         "not the string `derived`"),
        ("a value outside the bracket the equilibrium derives",
         mutated(("ph", "endorheic_alkalinity_bonus"), 9.0), "does not contain"),
        ("a bracket in another file its value sits outside",
         mutated(("surface_cover", "loess", "deposition_g_m2_yr"), 500.0,
                 "surface_classes"), "does not contain"),
    )
    for label, broken, expected in fixtures:
        said = _pedology_bracket_problems(broken, derived)
        if not any(expected in s for s in said):
            bad.append(f"the fixture {label!r} was not caught: expected a "
                       f"problem saying {expected!r}, got {said!r}")

    # The two dispositions whose enforcement is a file rather than a number.
    # Neither can be broken by editing the declaration, so the fixture moves the
    # enforcer instead: a generator that no longer names the key is the failure
    # each of them exists to catch.
    elsewhere = "pedology/scripts/_paths.py"   # a real file that names no bracket
    for kind, expected in (("level", "no longer names"),
                           ("derived", "nothing fills it")):
        table = tuple((k, lab, cfg, val, br, elsewhere if k == kind else enf)
                      for k, lab, cfg, val, br, enf in PEDOLOGY_BRACKETS)
        said = _pedology_bracket_problems(docs, derived, table)
        if not any(expected in s for s in said):
            bad.append(f"the fixture 'a {kind} bracket whose enforcer no longer "
                       f"reads it' was not caught: expected a problem saying "
                       f"{expected!r}, got {said!r}")
    return bad


def check_every_open_issue_is_batched() -> list[str]:
    """Every unclosed issue carries exactly one `batch:<n>` label.

    `docs/src/practice/conventions.md` records what the label means: where
    `step:<id>` says WHERE work lands, `batch:<n>` says WHEN, ordered by what a
    change to the row makes worthless under rule 7. Unlike the step marker it is
    not derivable from anything -- it is a judgement taken by reading the row,
    the same judgement the carve gate is -- so nothing regenerates it and a row
    filed without one simply falls out of the order. That is a silent failure of
    exactly the shape this file exists for: `bd ready --label batch:2` returns a
    shorter list and says nothing about what it left out.

    Read from `.beads/issues.jsonl` rather than by shelling out to `bd`. Two
    reasons and either alone decides it. This gate's cost is multiplied by every
    commit in the project, and a subprocess per invocation is what moved two
    other passes out of this file. And the export is the copy every OTHER
    checkout sees, so it is the right thing to hold to the convention: a label
    that exists only in one machine's Dolt working set is not yet a record.

    The cost of reading the export rather than the database is that a `bd label`
    write is invisible here until the export is regenerated, which the pre-commit
    hook does. So this check can lag by the length of one session and cannot
    report a false FAILURE from that lag -- only a false pass, which the next
    commit closes.
    """
    import json

    export = ROOT / ".beads" / "issues.jsonl"
    if not export.exists():
        return [f"{export.relative_to(ROOT)} is missing: the tracker's export is "
                f"the only record of the issue set outside one machine's database"]

    known = {f"batch:{n}" for n in range(1, 10)}
    unbatched: list[str] = []
    multiple: list[str] = []
    unknown: list[str] = []
    undeclared_decision: list[str] = []
    for line in export.read_text().splitlines():
        if not line.strip():
            continue
        try:
            issue = json.loads(line)
        except json.JSONDecodeError as exc:
            return [f"{export.relative_to(ROOT)} is not valid JSONL: {exc}"]
        if issue.get("_type") != "issue" or issue.get("status") == "closed":
            continue
        labels = set(issue.get("labels") or [])
        batches = {label for label in labels if label.startswith("batch:")}
        if not batches:
            unbatched.append(issue["id"])
        elif len(batches) > 1:
            multiple.append(f"{issue['id']} carries {', '.join(sorted(batches))}")
        if batches - known:
            unknown.append(f"{issue['id']} carries {', '.join(sorted(batches - known))}")
        # The narrower claim implies the wider one: a decision that gates the
        # current pass is a decision. Carrying only the narrow label hides the
        # row from the query that is meant to list everything the author owes.
        if "decision-blocks-loop-a" in labels and "needs-decision" not in labels:
            undeclared_decision.append(issue["id"])

    def _some(ids: list[str], limit: int = 12) -> str:
        # A wholesale loss would otherwise print hundreds of lines and bury the
        # rest of the gate's verdicts, which is what a gate exists not to do.
        head = ", ".join(ids[:limit])
        return head if len(ids) <= limit else f"{head} and {len(ids) - limit} more"

    problems = []
    if unbatched:
        problems.append(f"{len(unbatched)} unclosed issues carry no batch:<n> "
                        f"label: {_some(sorted(unbatched))}")
    if multiple:
        problems.append(f"{len(multiple)} carry more than one, so the order is "
                        f"ambiguous for them: {_some(sorted(multiple))}")
    if unknown:
        problems.append(f"batch labels outside batch:1 through batch:9, which "
                        f"conventions.md does not define: {_some(sorted(unknown))}")
    if undeclared_decision:
        problems.append(f"{len(undeclared_decision)} carry decision-blocks-loop-a "
                        f"without needs-decision, so they are absent from the "
                        f"list of what the author owes: {_some(sorted(undeclared_decision))}")
    return problems


# Where the bin centres sit in each signature. `annual_mean_of` takes a Dataset
# and reads `time` off it itself, so it is not here and cannot carry the defect.
CENTRES_ARGUMENT = {"bin_weights": 0, "infer_ntimes": 0,
                    "annual_mean": 1, "masked_mean": 1}
# Everything that MANUFACTURES an axis rather than reading one off a file.
SYNTHESISED_AXIS = {"arange", "linspace", "range", "indices", "ogrid", "mgrid"}


def _climatology_call_sites() -> list[tuple]:
    """Every call into `lib/climatology.py`, as (path, lineno, function, args).

    Its own file list, and not `SCRIPT_DIRS`, because `SCRIPT_DIRS` does not
    reach `aeolian/`, `ocean/`, `minerals/` or `analysis/` -- and the two call
    sites that carried the defect below were in one of them, invisible to every
    file-based check in this file.

    Names are resolved rather than matched: `annual_mean` is a climatology
    function in the modules that import it and a LOCAL function taking weights
    in `aeolian/scripts/build_sea_salt.py`, and a text match cannot tell those
    apart.
    """
    # PRUNED rather than filtered: `.claude/worktrees` holds whole checkouts of
    # this repo and `.venv` holds the installed world, so a walk that descends
    # and then discards costs minutes and reports every other session's copy of
    # a call site as if it were this tree's.
    prune = {".git", ".venv", ".claude", "vendor", "archive", "__pycache__",
             "node_modules", "build", ".mypy_cache", ".pytest_cache"}
    definitions = (ROOT / "lib" / "climatology.py").resolve()
    paths = []
    for parent, dirs, names in os.walk(ROOT):
        dirs[:] = sorted(d for d in dirs if d not in prune)
        paths.extend(Path(parent) / n for n in sorted(names)
                     if n.endswith(".py"))
    # This file is excluded because its own fixture below calls the weighting
    # on an index ON PURPOSE: that call is the control that proves the two
    # spellings differ, and a scan that flagged it would be flagging the test.
    gate = Path(__file__).resolve()
    sites = []
    for path in paths:
        if path.resolve() in (definitions, gate):
            continue
        try:
            source = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        if "climatology" not in source:
            continue                    # cheaper than parsing the reference trees
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue                    # check_modules_parse owns this
        modules, direct = set(), {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.name.split(".")[-1] == "climatology":
                        modules.add(a.asname or a.name)
            elif isinstance(node, ast.ImportFrom):
                if (node.module or "").split(".")[-1] == "climatology":
                    for a in node.names:
                        direct[a.asname or a.name] = a.name
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            f, name = node.func, None
            if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) \
                    and f.value.id in modules:
                name = f.attr
            elif isinstance(f, ast.Name) and f.id in direct:
                name = direct[f.id]
            if name in CENTRES_ARGUMENT:
                sites.append((path, node.lineno, name, node.args))
    return sites


def check_bin_weights_take_bin_centres() -> list[str]:
    """The bin weighting is handed a time coordinate, never a bin index, and it
    refuses the axes that cannot answer.

    `lib/climatology.py` recovers how many raw records each bin of a pyburn
    climatology holds from the SPACING of the bin centres, and weights the
    annual mean by those counts. Hand it `np.arange(nbin)` and the spacing is
    uniform, which asserts that the bins are even instead of measuring them; a
    uniform weight vector is an ordinary array, and the annual mean it produces
    is the one the project had before the module existed.

    THIS FAILS SILENTLY IN BOTH DIRECTIONS, which is why it is a gate and not a
    review note. On a climatology whose bins happen to be nearly even the wrong
    answer is nearly right, so no result looks odd; and the call is right by
    coincidence rather than by construction, so it stops being right when the
    cadence or the window changes and still nothing objects.

    Three parts, and the last two are what make the first a test rather than a
    lint:

    1. **No call site synthesises the axis.** Resolved through the imports,
       so the local `annual_mean(field, weights)` in `build_sea_salt.py` is not
       confused with the climatology function of the same name.
    2. **The weighting inverts real centres.** A fixture reconstructs the bin
       centres pyburn's own binning implies for a record count, over a sweep of
       counts, and requires that the weighting inverts them back to that record
       count and returns the counts themselves. The CONTROL is uniformity:
       wherever the bin count does not divide the record count the answer has
       to be measurably apart from equal weights, or passing an index would be
       undetectable. Without that control a weighting that had collapsed to
       uniform would pass half of this and prove nothing. A declaration of the
       record count is checked against the centres too, so one from another run
       is refused rather than reweighting the file to fit.
    3. **The axes that carry no answer are refused.** An index is one, and an
       evenly spaced axis over an even bin count is the other: [1, 2, 1, 2, ...]
       over twelve bins stamps exactly the centres [3, 3, 3, ...] stamps at
       twice the write interval, so the file cannot say which it is and the two
       weight neighbouring bins a factor of two apart. The sweep holds 18 and 30
       records in twelve bins inside it for that reason, and asserts both the
       refusal and the answer a declared record count buys.
    """
    problems = []
    for path, lineno, name, args in _climatology_call_sites():
        pos = CENTRES_ARGUMENT[name]
        if len(args) <= pos:
            continue                    # passed by keyword, or a bad call
        arg = args[pos]
        made = sorted({
            n.func.id if isinstance(n.func, ast.Name) else n.func.attr
            for n in ast.walk(arg)
            if isinstance(n, ast.Call)
            and (n.func.id if isinstance(n.func, ast.Name)
                 else getattr(n.func, "attr", None)) in SYNTHESISED_AXIS})
        if made:
            problems.append(
                f"{path.relative_to(ROOT)}:{lineno} calls {name}() on "
                f"{made[0]}(), which is a bin INDEX. The weighting reads the "
                f"spacing of the bin centres, so evenly spaced integers assert "
                f"that the bins are even instead of measuring them, and it "
                f"returns equal weights without raising. Pass the file's own "
                f"`time` variable")

    sys.path.insert(0, str(ROOT / "lib"))
    try:
        import numpy as np
        import climatology as clim
    except ImportError as exc:
        return problems + [f"lib/climatology.py does not import: {exc}"]

    DT = 32.0            # a write interval; the recovery is scale-free in it
    T0 = 31.0            # and offset-free, so neither is a project number
    checked = 0
    alternating = 0
    for nbin in (4, 5, 12):
        uniform = np.full(nbin, 1.0 / nbin)
        # THE INDEX CONTROL. A bin index is refused at the boundary, so what is
        # asserted is the refusal. Both spellings are here because an index
        # scaled by a constant is still an index and still evenly spaced.
        for index in (np.arange(nbin), np.arange(nbin) * 5.0):
            try:
                clim.bin_weights(index)
            except ValueError:
                pass
            else:
                problems.append(
                    f"lib/climatology.py takes {index.tolist()} as the bin "
                    f"centres of a {nbin}-bin climatology. That is a bin "
                    f"INDEX: the weighting reads the spacing of real centres, "
                    f"so evenly spaced integers from the time origin assert "
                    f"that the bins are even instead of measuring them")
        for ntimes in range(nbin, 250):
            counts = clim.counts_for(ntimes, nbin)
            if counts.min() <= 0:
                continue
            # pyburn stamps each bin with the MEAN of the raw timestamps in it,
            # so with a regular stream the gap between two centres is the mean
            # of their two record counts.
            starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
            centres = T0 + DT * (starts + (counts - 1) / 2.0)
            gaps = np.diff(centres)
            recoverable = not np.allclose(gaps, gaps[0], rtol=1e-9, atol=1e-9)
            even = bool(np.all(counts == counts[0]))
            checked += 1
            if not recoverable:
                # THE EVENLY SPACED AXES. Counts that alternate about a constant
                # pair sum, [1, 2, 1, 2, ...] and the families above it, stamp
                # the same centres as an EQUAL binning of a different record
                # count: 18 records in 12 bins written every 320 timesteps and
                # 36 written every 160 give the identical twelve numbers. So the
                # file carries neither the record count nor the weights, and the
                # two answers are a factor of two apart on neighbouring bins.
                # The weighting refuses that rather than picking one, and
                # answers only when the caller declares what the run wrote.
                declared = clim.bin_weights(centres, ntimes)
                if not np.allclose(declared * ntimes, counts, rtol=0, atol=1e-9):
                    problems.append(
                        f"a declared {ntimes} records in {nbin} bins weights "
                        f"them {(declared * ntimes).tolist()}, not by their "
                        f"record counts {counts.tolist()}")
                if not even:
                    alternating += 1
                    # The refusal has to protect something: an alternating
                    # binning sits exactly half a record per bin from uniform,
                    # so returning uniform here is a real error and not a
                    # rounding one.
                    spread = float(np.max(np.abs(declared - uniform)))
                    if spread + 1e-12 < 0.5 / ntimes:
                        problems.append(
                            f"{nbin} bins holding {counts.tolist()} are within "
                            f"{spread:.3e} of uniform, so this case proves "
                            f"nothing about the refusal")
                if nbin % 2 == 0:
                    # An even bin count always admits both shapes, so no answer
                    # is derivable from the centres alone.
                    try:
                        got = clim.bin_weights(centres)
                    except ValueError:
                        pass
                    else:
                        problems.append(
                            f"{nbin} bins holding {counts.tolist()} have evenly "
                            f"spaced centres, which both a multiple of {nbin} "
                            f"records and an odd multiple of {nbin // 2} "
                            f"records produce, and the weighting answered "
                            f"{got.tolist()} instead of refusing")
                else:
                    # An odd bin count admits only the equal shape, so the
                    # weights are determined even though the record count is
                    # not, and refusing would be over-strict.
                    w = clim.bin_weights(centres)
                    if not np.allclose(w, uniform, rtol=0, atol=1e-12):
                        problems.append(
                            f"{nbin} bins holding {counts.tolist()} admit only "
                            f"equal weights, and the weighting returned "
                            f"{w.tolist()}")
                continue
            got = clim.infer_ntimes(centres)
            if got != ntimes:
                problems.append(
                    f"bin centres built from {ntimes} records in {nbin} bins "
                    f"invert to {got} records")
                continue
            w = clim.bin_weights(centres)
            if abs(float(w.sum()) - 1.0) > 1e-12:
                problems.append(f"weights for {ntimes} records in {nbin} bins "
                                f"sum to {float(w.sum()):.12f}, not 1")
            # The weights ARE the record counts, so scaling them back by the
            # record count has to return integers.
            if not np.allclose(w * ntimes, counts, rtol=0, atol=1e-9):
                problems.append(
                    f"weights for {ntimes} records in {nbin} bins are not the "
                    f"record counts {counts.tolist()}: {(w * ntimes).tolist()}")
            # A declaration is checked against the centres rather than trusted,
            # so the count from another run or another I/O regime is refused.
            for wrong in (ntimes + 1, ntimes * 2):
                if np.allclose(clim.counts_for(wrong, nbin) / wrong,
                               counts / ntimes, rtol=0, atol=1e-12):
                    continue            # same weights, so nothing to refuse
                try:
                    clim.bin_weights(centres, wrong)
                except ValueError:
                    pass
                else:
                    problems.append(
                        f"centres from {ntimes} records in {nbin} bins accept "
                        f"a declared {wrong} records, so a declaration is "
                        f"taken on trust rather than checked against the file")
            # THE CONTROL, and the whole point. An uneven split moves at least
            # one whole record between bins, so the answer from the centres and
            # the answer an index would have given are at least half a record
            # apart. If they were not, passing an index would be undetectable.
            spread = float(np.max(np.abs(w - uniform)))
            if spread < 0.5 / ntimes:
                problems.append(
                    f"{nbin} bins hold {counts.tolist()} records and the "
                    f"weighting is within {spread:.3e} of uniform, so passing "
                    f"an index here would be undetectable")
    if checked < 300:
        problems.append(f"the sweep only reached {checked} record counts; it is "
                        f"meant to cover a few hundred and something has "
                        f"narrowed it")
    if alternating < 2:
        problems.append(f"the sweep reached {alternating} alternating bin "
                        f"patterns; 18 and 30 records in 12 bins are the cases "
                        f"the weighting has to refuse and they have gone out "
                        f"of the fixture")
    return problems


def main() -> None:
    argparse.ArgumentParser(
        description="The fast static gate: every check here is a read, a parse "
                    "or a small synthetic fixture. The two passes that spawn a "
                    "process per unit are scripts/verify_entry_points.py and "
                    "exoplasim/scripts/verify_model_compiles.py, which run "
                    "before a build or a push. Takes no arguments.").parse_args()

    files = sorted({f for d in SCRIPT_DIRS if d.is_dir()
                    for f in d.glob("*.py")})
    # The rung lints read shell too: SPAT-2's copies did not all land in Python,
    # and `d.glob("*.py")` made every shell offender invisible to both of them.
    shell_files = sorted({f for d in SCRIPT_DIRS if d.is_dir()
                          for f in d.glob("*.sh")})
    vendor_files = sorted({f for d in VENDOR_PY_DIRS if d.is_dir()
                           for f in d.glob("*.py")})
    print(f"{len(files)} modules under {len(SCRIPT_DIRS)} directories\n")

    checks = [("every module parses", lambda: check_modules_parse(files)),
              ("no undefined names", lambda: check_undefined_names(files)),
              ("every per-build data path is namespaced by the build",
               lambda: check_build_scoped_defaults()),
              ("no artifact selection by sort order", lambda: check_no_order_picks(files)),
              ("one grid convention, in lib/gridding.py",
               lambda: check_one_grid_convention(files)),
              ("every generator is declared in config/pipeline.yaml",
               lambda: check_registered_in_workflow(files)),
              ("every registered script exists",
               lambda: check_registered_paths_exist()),
              ("every step is named in its component README",
               lambda: check_documented_in_component(files)),
              ("purge never reaches the terrain, from any seed",
               lambda: check_purge_never_reaches_the_terrain()),
              ("purge deletes every path the seed's own writes names",
               lambda: check_purge_covers_the_seeds_own_writes()),
              ("local_slope_deg reproduces an analytic gradient",
               lambda: check_slope_fit()),
              ("a re-derived input file is caught under an unchanged name",
               lambda: check_input_stamp_guard()),
              ("a staged surface field from another build is refused",
               lambda: check_staged_surface_build_guard()),
              ("a donor restart is refused for a frozen field, not a reread one",
               lambda: check_donor_surface_guard()),
              ("the convergence window follows the declared purposes",
               lambda: check_production_window()),
        ("a resume adopts the run's own thread count",
               lambda: check_resume_adopts_the_manifest_thread_count()),
        ("the production span comes from the run's own criteria",
               lambda: check_production_span_is_self_limiting()),
        ("a staged field is never written through a symlink",
               lambda: check_no_write_through_a_symlink()),
        ("orbits no segment covers refuse rather than defaulting",
               lambda: check_uncovered_trailing_orbits_refuse()),
        ("a verdict window is refused across an I/O-regime change",
               lambda: check_io_regime_window()),
        ("the I/O step at a join is measured against the offset criterion",
               lambda: check_io_step_measurement()),
              ("no control patch is left in the model source",
               lambda: check_no_control_patch()),
              ("SHTns is asked for the one mode that runs no timing race",
               lambda: check_shtns_init_is_deterministic()),
              ("pyburn's compiled extensions load on this interpreter",
               lambda: check_compiled_extensions()),
              ("no OpenMP directive line is truncated",
               lambda: check_omp_directive_length()),
              ("no continuation marker is dropped in the model Fortran",
               lambda: check_no_dropped_continuation()),
              ("no diagnostics write is reachable by every thread",
               lambda: check_diag_writes_are_answered()),
              ("no model local has all its definitions inside a where",
               lambda: check_masked_only_locals_are_answered()),
              ("each step resolves a climatology its needs permit",
               lambda: check_climatology_needs_match_call_sites()),
              ("the best-available resolver names the stage it returned",
               lambda: check_best_available_climatology_resolves_by_stage()),
              ("no imported module name is rebound",
               lambda: check_no_shadowed_imports(files)),
              ("no name is loaded that nothing binds, model Python included",
               lambda: check_no_unbound_names(files + vendor_files)),
              ("the spectral tail slope is fitted above the roundoff floor",
               lambda: check_tail_fit_stops_above_roundoff()),
              ("the restart schema covers every record the model writes",
               lambda: check_restart_schema_covers_the_model()),
              ("every postprocessor code requested is one something produces",
               lambda: check_requested_codes_are_produced()),
              ("the autocorrelation estimator recovers a known answer",
               lambda: check_autocorrelation_estimator()),
              ("the bin weighting is handed bin centres, never a bin index",
               lambda: check_bin_weights_take_bin_centres()),
              ("the transform gates run the configured spectral filter",
               lambda: check_gate_filter_matches_config()),
              ("a continuation redeclares what a prepare declared",
               lambda: check_continuation_redeclares_everything()),
              ("every per-level namelist key is written for every level",
               lambda: check_per_level_namelist_keys_cover_every_level()),
              ("no artifact path carries a resolution literal",
               lambda: check_no_rung_literal_in_a_path(files + shell_files)),
              ("the rung-to-dimension table is lib/rungs.py and nowhere else",
               lambda: check_no_rung_table_outside_rungs(files + shell_files)),
              ("the configured resolution matches its own grid dimensions",
               lambda: check_configured_grid()),
              ("every restatement of the ladder agrees with lib/rungs.py",
               lambda: check_ladder_restatements()),
              ("the ceiling and the route agree with what they restate",
               lambda: check_ladder_timestep_declarations()),
              ("both columns model snow from lib/snow.py's one relation",
               lambda: check_snow_conductivity_restatements()),
              ("the configured timestep is one the route runs this rung at",
               lambda: check_configured_timestep()),
              ("a resume refuses a rewritten spectrum file",
               lambda: check_spectrum_guard()),
              ("the convergence slab's sea water follows the run",
               lambda: check_slab_capacity_follows_the_run()),
              ("the state closure's melting point follows the run",
               lambda: check_melting_point_follows_the_run()),
              ("a declared run length is derived from the right timescale",
               lambda: check_run_length_derivation()),
              ("the window's nominal inputs still bound the reports",
               lambda: check_convergence_bounds_are_re_read()),
              ("the relaxation bracket is read from its producer",
               lambda: check_relaxation_bracket_is_read()),
              ("the fit's tail fraction is stated once",
               lambda: check_fit_tail_fraction_is_stated_once()),
              ("every commissioning row is re-read from its run record",
               lambda: check_commissioning_evidence_is_re_read()),
              ("the model's shortwave cloud tables are the papers' tables",
               lambda: check_cloud_tables_match_the_papers()),
              ("the tools environment.md names are on this host",
               lambda: check_documented_tools()),
              ("every pedology bracket contains the value it brackets",
               lambda: check_pedology_values_are_inside_their_brackets()),
              ("the aerosol steering wind is weighted by layer mass",
               lambda: check_steering_weight_is_layer_mass()),
              ("every unclosed issue carries exactly one batch:<n>",
               lambda: check_every_open_issue_is_batched())]
    # Run and REPORT one at a time, rather than evaluating the list and then
    # printing it. A gate that prints nothing until it is finished cannot be
    # told apart from a gate that has hung, and stdout is block-buffered into a
    # pipe, so the flush is the half that makes it true.
    failed = 0
    for name, run in checks:
        problems = run()
        if problems:
            failed += 1
            print(f"[ FAIL ] {name}", flush=True)
            for p in problems:
                print(f"         {p}", flush=True)
        else:
            print(f"[  ok  ] {name}", flush=True)

    print(f"\n{len(checks)} checks, {failed} failed")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
