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
11. **`TASKS.md` header counts match their tables.** The count is what the
   one-table-per-prefix layout exists to make readable in one place, and a
   hand-maintained number goes stale the way every other one in this project
   has.
12. **A resume refuses a rewritten spectrum file.** The config names the
   spectrum and the model reads the file, so a config comparison cannot see
   `k25v.dat` regenerated in place. CONS-3.
13. **The tools `environment.md` names are actually on this host.** That
   document sends a reader to `ncdump`, NCO, `h5diff` and `yq` rather than a
   Python session, and nothing else checks the claim is true. Both
   directions, so the document and the check cannot drift apart.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
from pathlib import Path
import re
import shutil
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

    Cheap, and it fails the moment a script is renamed without the register
    following, which is the whole point of having one.
    """
    import yaml
    graph = yaml.safe_load(
        (ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    named = ([(s["script"], f"step {s['id']}") for s in graph["steps"]]
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


def check_task_counts() -> list[str]:
    """Each TASKS.md prefix's header count matches the table under it.

    The header line is the answer to "what is left here" and "what may I number
    next", and its convention says the table is what it is computed FROM -- so
    the two disagreeing means one of the two questions is being answered wrong.
    It has drifted once already: CLIM read `4 open of 30 issued` against a table
    of 31 rows with 3 open, both halves stale from the same commit, which is
    what a hand-maintained count does. Openness is read from the status
    column, not from where a row sits.
    """
    text = (ROOT / "TASKS.md").read_text(encoding="utf-8")
    problems, prefix, header = [], None, None
    issued = open_now = 0

    def settle() -> None:
        if not prefix:
            return
        if header is None:
            problems.append(f"{prefix} has no `N open of M issued.` line")
        elif header != (open_now, issued):
            problems.append(
                f"{prefix} says {header[0]} open of {header[1]} issued; "
                f"its table has {open_now} open of {issued}")

    for line in text.splitlines():
        m = re.match(r"^## ([A-Z]+) --", line)
        if m:
            settle()
            prefix, header, issued, open_now = m.group(1), None, 0, 0
            continue
        m = re.match(r"^(\d+) open of (\d+) issued\.", line)
        if m and prefix:
            header = (int(m.group(1)), int(m.group(2)))
            continue
        if not line.startswith("| ") or line.startswith(("| id", "| ---")):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 4 or not re.fullmatch(r"[A-Z]+-[0-9a-z]+", cells[0]):
            continue
        issued += 1
        if not cells[3].lower().startswith(("done", "wontfix")):
            open_now += 1
    settle()
    if prefix is None:
        problems.append("no prefix sections found; has TASKS.md changed shape?")
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

    The bar is 132, which is what `most_compiler_omp` passes as
    `-ffixed-line-length-132`, and the check is over DIRECTIVES only: ordinary
    over-length lines elsewhere in the vendored sources are upstream's and are
    not silent in the same way.
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
               check_registered_in_workflow(files)),
              ("every registered script exists",
               check_registered_paths_exist()),
              ("every step is named in its component README",
               check_documented_in_component(files)),
              ("purge never reaches the terrain, from any seed",
               check_purge_never_reaches_the_terrain()),
              ("local_slope_deg reproduces an analytic gradient",
               check_slope_fit()),
              ("the convergence window follows the declared purposes",
               check_production_window()),
              ("TASKS.md counts match their tables", check_task_counts()),
              ("no control patch is left in the model source",
               check_no_control_patch()),
              ("no OpenMP directive line is truncated",
               check_omp_directive_length()),
              ("no imported module name is rebound",
               check_no_shadowed_imports(files)),
              ("a resume refuses a rewritten spectrum file",
               check_spectrum_guard()),
              ("the tools environment.md names are on this host",
               check_documented_tools())]
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
