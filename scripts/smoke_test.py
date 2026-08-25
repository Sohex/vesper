#!/usr/bin/env python3
"""Exercise every script's entry points without running any science.

    python scripts/smoke_test.py

`check_consistency.py` audits the artifacts. This audits the code that produces
them, which is where a different class of bug lives: in one day this project
found basin ids read from the wrong terrain, six scripts defaulting to a stale
data directory, five stale binaries, and a config change that left `latitudes`
behind. Every one of those was reachable by importing a module and looking at
where its defaults pointed -- none needed a model run.

The checks, all cheap (plus registered-script existence and the
purge-never-reaches-the-terrain property, run from `main()` with the rest):

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
12g. **The model source compiles**, front end only, under the declared flags.
   Every other check of the Fortran here is a grep or a parse, so a tree that no
   binary could be built from passed all of them. `--skip-compile` opts out.
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
"""

from __future__ import annotations

import argparse
import ast
import builtins
import importlib.util
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
        doomed = pipeline.downstream(seed, by_id)
        if "orogen" in doomed:
            bad.append(f"--purge {seed} reaches orogen, so it would offer to "
                       f"delete the export")
        if seed in doomed:
            bad.append(f"--purge {seed} includes {seed} itself")
    # The export edge is what seeds `--purge orogen`; a typo in it silently
    # empties that seed rather than erroring, and the purge would report success
    # having deleted nothing.
    readers = [s["id"] for s in graph["steps"] if s.get("reads_export")]
    if not readers:
        bad.append("no step declares reads_export, so --purge orogen is a no-op")
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
        d = run_dir(tmp, [seg(1, 76, "spinup")])
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


def check_model_source_compiles() -> list[str]:
    """The model source passes gfortran's front end, under the declared flags.

    THE GAP THIS CLOSES. Every other check of the Fortran in this file is a grep
    or a text parse, so the property that makes the source usable at all was
    checked by nothing until somebody ran `rebuild_binaries.py` by hand -- and
    rule 4 means that is the most expensive moment to find out. A tree no binary
    could be built from reported 24 of 24 green here. world-adts.

    WHAT IT COSTS, in the units that survive a contended host: 36 translation
    units, one `gfortran -fsyntax-only` each, against a full build's same set
    compiled AND optimised at `-O2 -march=znver4 -funroll-loops` and then linked
    against SHTns and FFTW. The front end reads the source, resolves the `use`
    graph and writes module files; there is no code generation and no link. That
    is the same shape of cost as the `--help` pass this file already runs over
    every script, and `--skip-compile` opts out of it the same way.

    IT DOES NOT BUILD A BINARY, deliberately. Whether an executable is CURRENT is
    rule 4's question and `check_consistency.py` answers it; whether the source
    compiles is this one. `exoplasim/scripts/verify_model_compiles.py` holds the
    gate itself, including why the flag line has to be the declared one.
    """
    script = ROOT / "exoplasim" / "scripts" / "verify_model_compiles.py"
    if not script.is_file():
        return [f"{script.relative_to(ROOT)} is gone, and it is what checks "
                "that the model source compiles"]
    r = subprocess.run([sys.executable, str(script)],
                       capture_output=True, text=True, cwd=ROOT)
    if r.returncode == 0:
        return []
    return [(r.stderr.strip() or r.stdout.strip() or "(no output)")]


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

    THE CONTROL is below and is what makes this a test rather than a
    description: the same check on T42 at T21's step has to fail.
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

    # THE NEGATIVE CONTROL. world-td3's configuration exactly: the rung moved to
    # T42 and the step left at T21's 45. If this passes, the check above cannot
    # fail and is worth nothing.
    if not rungs.timestep_problems("T42", 45.0):
        problems.append(
            "the timestep check does not refuse T42 at dt 45, which is the "
            "configuration world-td3 ran and blew up: a check that cannot fail "
            "is not a check")
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
            + restart_schema.check_first_window_matches_the_rest(src))


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

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-help", action="store_true",
                    help="skip the subprocess --help pass, which dominates runtime")
    ap.add_argument("--skip-compile", action="store_true",
                    help="skip the gfortran -fsyntax-only pass over the model source")
    args = ap.parse_args()

    files = sorted({f for d in SCRIPT_DIRS if d.is_dir()
                    for f in d.glob("*.py")})
    # The rung lints read shell too: SPAT-2's copies did not all land in Python,
    # and `d.glob("*.py")` made every shell offender invisible to both of them.
    shell_files = sorted({f for d in SCRIPT_DIRS if d.is_dir()
                          for f in d.glob("*.sh")})
    vendor_files = sorted({f for d in VENDOR_PY_DIRS if d.is_dir()
                           for f in d.glob("*.py")})
    print(f"{len(files)} modules under {len(SCRIPT_DIRS)} directories\n")

    checks = [("imports", check_imports(files)),
              ("no undefined names", check_undefined_names(files)),
              ("every per-build data path is namespaced by the build",
               check_build_scoped_defaults()),
              ("no artifact selection by sort order", check_no_order_picks(files)),
              ("one grid convention, in lib/gridding.py",
               check_one_grid_convention(files)),
              ("every generator is declared in config/pipeline.yaml",
               check_registered_in_workflow(files)),
              ("every registered script exists",
               check_registered_paths_exist()),
              ("every step is named in its component README",
               check_documented_in_component(files)),
              ("purge never reaches the terrain, from any seed",
               check_purge_never_reaches_the_terrain()),
              ("local_slope_deg reproduces an analytic gradient",
               check_slope_fit()),
              ("a re-derived input file is caught under an unchanged name",
               check_input_stamp_guard()),
              ("a staged surface field from another build is refused",
               check_staged_surface_build_guard()),
              ("the convergence window follows the declared purposes",
               check_production_window()),
              ("no control patch is left in the model source",
               check_no_control_patch()),
              ("SHTns is asked for the one mode that runs no timing race",
               check_shtns_init_is_deterministic()),
              ("pyburn's compiled extensions load on this interpreter",
               check_compiled_extensions()),
              ("no OpenMP directive line is truncated",
               check_omp_directive_length()),
              ("no continuation marker is dropped in the model Fortran",
               check_no_dropped_continuation()),
              ("no imported module name is rebound",
               check_no_shadowed_imports(files)),
              ("no name is loaded that nothing binds, model Python included",
               check_no_unbound_names(files + vendor_files)),
              ("the spectral tail slope is fitted above the roundoff floor",
               check_tail_fit_stops_above_roundoff()),
              ("the restart schema covers every record the model writes",
               check_restart_schema_covers_the_model()),
              ("the transform gates run the configured spectral filter",
               check_gate_filter_matches_config()),
              ("a continuation redeclares what a prepare declared",
               check_continuation_redeclares_everything()),
              ("every per-level namelist key is written for every level",
               check_per_level_namelist_keys_cover_every_level()),
              ("no artifact path carries a resolution literal",
               check_no_rung_literal_in_a_path(files + shell_files)),
              ("the rung-to-dimension table is lib/rungs.py and nowhere else",
               check_no_rung_table_outside_rungs(files + shell_files)),
              ("the configured resolution matches its own grid dimensions",
               check_configured_grid()),
              ("every restatement of the ladder agrees with lib/rungs.py",
               check_ladder_restatements()),
              ("the ceiling and the route agree with what they restate",
               check_ladder_timestep_declarations()),
              ("the configured timestep is one the route runs this rung at",
               check_configured_timestep()),
              ("a resume refuses a rewritten spectrum file",
               check_spectrum_guard()),
              ("the tools environment.md names are on this host",
               check_documented_tools())]
    if not args.skip_help:
        checks.insert(1, ("entry points answer --help", check_help(files)))
    if not args.skip_compile:
        checks.append(("the model source compiles under the declared flags",
                       check_model_source_compiles()))

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
