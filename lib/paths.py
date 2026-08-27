"""Path helpers shared across components, and the guards on the climatology.

`rel` exists because of a specific recurring bug. There are TWO climatologies:
`bootstrap_climatology_path` resolves the run on terrain-only surface fields,
which the derived fields are built FROM, and `climatology_path` the baseline run
on those fields once they exist. `best_available_climatology` resolves the more
DETERMINED of the two and says which it returned, for the steps whose answer
depends on the climate state rather than on the calendar or the grid.
`config/pipeline.yaml` says which climatology each step needs to EXIST, and
`scripts/smoke_test.py:check_climatology_needs_match_call_sites` holds the call
sites to it. `require_clean_io` and `require_configured_grid` are what all three
must survive.

`Path.relative_to` RAISES when the path is not under the given root. Every script
here prints "wrote <path>" relative to the project root, and that print happens
*after* the artifact has been written. So passing `--output` to somewhere outside
the tree -- /tmp, a scratch directory, an absolute path typed by hand -- writes
the file correctly and then dies with a ValueError, leaving an artifact on disk
and a traceback on the terminal that looks like the run failed.

`pedology/scripts/build_soil.py` hit this and left a soil map with no provenance
record. It was then fixed at one of its two call sites and not the other, which
is the more instructive half of the story: the fix has to be applied everywhere
the pattern appears, not just where it was first observed. See
`docs/src/practice/failure-modes.md`.

Provenance records have the same problem for a different reason: a path recorded
as absolute is not portable, and a path that raises is worse than either.
"""

from __future__ import annotations

from netCDF4 import Dataset
from pathlib import Path
from typing import NamedTuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]

BASELINE = "baseline"
BOOTSTRAP = "bootstrap"
UNNAMED = "unnamed"


class Climatology(NamedTuple):
    """A climatology WITH the stage of determination it came from.

    `stage` is `baseline`, `bootstrap`, or `unnamed` for a file handed in on the
    command line that is neither of the two `config/planet.yaml` names. It exists
    so the choice reaches the product's provenance record: a reader has to be
    able to tell a first-pass artifact from a later one without re-deriving it.
    """

    path: Path
    stage: str


def rel(path: Path | str, root: Path | None = None) -> str:
    """Path relative to the project root, or absolute if it lies outside.

    Never raises. Use for anything printed to a human or written into a
    provenance record.
    """
    p = Path(path)
    base = Path(root) if root is not None else PROJECT_ROOT
    try:
        return str(p.relative_to(base))
    except ValueError:
        return str(p)


def bootstrap_climatology_path(root: Path | None = None) -> Path:
    """The climatology taken on TERRAIN-ONLY surface fields.

    A SECOND KEY BECAUSE THERE ARE TWO CLIMATOLOGIES AND THE GRAPH HAS ALWAYS
    SAID WHICH IS WHICH. `config/pipeline.yaml` is the register of which steps
    declare `needs: bootstrap_climatology` and which declare
    `baseline_climatology`; it is not restated here, because a list in two
    places is a list that goes out of step. Config carried one key for both, so
    every one of those steps resolved to the same file and the graph's
    distinction reached nothing.
    `scripts/smoke_test.py:check_climatology_needs_match_call_sites` is what
    holds each step's script to its own declaration.

    THE ORDER IS WHY IT MATTERS. The bootstrap exists to produce the climatology
    the derived surface fields are built FROM, and the baseline is the run on
    those fields. A step that needs the bootstrap and reads the baseline is
    asking for an artifact that does not exist at all on a first pass.

    WHAT THIS IS FOR, now that a third resolver exists. Call this where the
    answer does NOT depend on the climate state -- a phase of the model
    calendar, a grid, a land mask -- because the two climatologies carry the
    same one and pinning to the bootstrap keeps the step runnable on a first
    pass without asking for anything it cannot use. Where the answer DOES
    depend on the state, `best_available_climatology` is the resolver: on a
    later pass the bootstrap is no longer the best available, only the
    earliest, and its docstring carries the argument.

    Same no-fallback rule as `climatology_path`, and for the same reason: a
    fallback returns a plausible number from a different world instead of an
    error. Same grid guard, applied here so callers inherit it.
    """
    import yaml
    project = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    config = yaml.safe_load(
        (project / "config" / "planet.yaml").read_text(encoding="utf-8"))
    declared = config.get("bootstrap_climatology")
    if not declared:
        raise SystemExit(
            "config/planet.yaml has no `bootstrap_climatology`. Name one there "
            "or pass --climatology; there is deliberately no fallback. This is "
            "the TERRAIN-ONLY climatology the derived surface fields are built "
            "from, and it is not `baseline_climatology`.")
    path = project / declared
    if path.is_file():
        require_configured_grid(path, config)
    return path


def climatology_path(name: str | None = None, root: Path | None = None) -> Path:
    """The climatology every downstream component is driven from.

    Read from `config/planet.yaml`'s `baseline_climatology`, with no fallback,
    because a fallback here is the most expensive kind of bug this project has:
    it returns a plausible number computed from a different world instead of an
    error. The default it replaced pointed at `climatology_s096`, which is
    pre-carve terrain under the superseded `k2` spectrum and the surface the
    antipodal carve verdict was taken from.

    It lives in `lib/` because copies of it did not stay in step. pedology and
    biosphere were migrated to the config key; `surface_water.py` kept a private
    module constant, and `build_surface_albedo.py` and
    `build_surface_soil_water.py` kept an argparse default, all three still
    naming the superseded directory. Two of those only surfaced when a
    baseline re-run ran them for the first time in months.

    Returns the FILE. Callers must not rebuild the name from a directory: the
    label is chosen per product, so a bootstrap climatology is
    `bootstrap_regular_climatology.nc` and reconstructing
    `baseline_regular_climatology.nc` beside it finds nothing.

    `name` still accepts a directory under `exoplasim/analysis/` for the old
    layout, so existing callers that pass one keep working.

    THE GRID IS CHECKED HERE. A climatology's name says nothing about its rung,
    so changing `model.resolution` used to neither rename nor invalidate nor
    refuse the declared file, and seven consumers then took the grid FROM the
    file and computed on whatever rung it happened to carry.
    `require_configured_grid` is what refuses that, and it is applied here so
    every caller inherits it instead of two builders having their own copy.
    """
    import yaml
    project = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    config = yaml.safe_load(
        (project / "config" / "planet.yaml").read_text(encoding="utf-8"))
    if name is not None:
        path = (project / "exoplasim" / "analysis" / name
                / "baseline_regular_climatology.nc")
    else:
        declared = config.get("baseline_climatology")
        if not declared:
            raise SystemExit(
                "config/planet.yaml has no `baseline_climatology`. Name one "
                "there or pass --climatology; there is deliberately no "
                "fallback.")
        path = project / declared
    # Only when it is there: a file that does not exist is the caller's error to
    # report, and several of them say something more useful about it than this
    # could.
    if path.is_file():
        require_configured_grid(path, config)
    return path


def climatology_stage(path: Path, root: Path | None = None) -> str:
    """Which of the two configured climatologies a given file IS.

    `baseline`, `bootstrap`, or `unnamed` when it is neither -- which is what a
    hand-typed `--climatology` is, and saying so is better than picking the
    nearer of the two. Compared by resolved path, because the config names are
    repo-relative and a caller's argument may be anything.
    """
    import yaml
    project = Path(root) if root is not None else PROJECT_ROOT
    config = yaml.safe_load(
        (project / "config" / "planet.yaml").read_text(encoding="utf-8"))
    resolved = Path(path).resolve()
    for stage, key in ((BASELINE, "baseline_climatology"),
                       (BOOTSTRAP, "bootstrap_climatology")):
        declared = config.get(key)
        if declared and (project / declared).resolve() == resolved:
            return stage
    return UNNAMED


def best_available_climatology(override: Path | None = None,
                               root: Path | None = None) -> Climatology:
    """The MOST DETERMINED climatology that exists, and which one that is.

    The baseline when `config/planet.yaml` names one, the bootstrap when it does
    not, and a refusal when neither key is set. `override` is the caller's
    `--climatology`, passed straight through and labelled by comparison so one
    call site covers both routes.

    WHAT IT IS FOR. Bootstrapping and looping move a world from its least
    determined state to its most determined one, so a step reads the best input
    available AT THE POINT IT RUNS and not the one its first pass happened to
    have. A step pinned to the bootstrap forever is right on the first pass,
    where the bootstrap is genuinely the best that exists, and wrong on every
    pass after, where it is merely the earliest. `docs/src/pipeline/loops.md`
    argues it; `pedology/scripts/build_soil.py` is the case that established it,
    taking its runoff from a terrain-only climatology that carries no lakes on
    any iteration while taking its vegetation from a run that had already
    reached the baseline.

    THIS IS NOT THE FALLBACK THE NO-FALLBACK RULE FORBIDS, and the difference is
    not a matter of degree. That rule exists because a fallback returns a
    PLAUSIBLE NUMBER COMPUTED FROM A DIFFERENT WORLD instead of an error: the
    default `climatology_path` used to carry named pre-carve terrain under a
    superseded spectrum, so a caller that had lost its climatology got numbers
    from a planet that no longer existed and nothing said so. The bootstrap and
    the baseline are the SAME world at two stages of determination -- same
    build, same terrain hash, same config, checked by `require_build` at every
    call site -- so neither answer is from somewhere else. Three further things
    separate them, and all three have to hold or this would be a fallback:

      1. The choice is made on what config DECLARES, never on what happens to be
         on disk. A named baseline that is missing raises through the caller's
         own existence check rather than quietly degrading to the bootstrap,
         which is exactly the failure the rule is about.
      2. The choice is STAMPED. Every product built through this records its
         stage, so a first-pass artifact is distinguishable from a later one by
         reading it rather than by re-deriving it.
      3. There is no third option. With neither key set this raises, and it
         raises naming both keys rather than reaching for a default.

    WHERE IT MUST NOT BE USED. A quantity that depends on the model calendar or
    on the grid rather than on the climate STATE is equally determined at either
    stage, and calling this for one of those makes a step newly dependent on a
    baseline being named while telling it nothing.
    `biosphere/scripts/build_vesper_header.py` fits a solstice offset against
    solar declination and `aeolian/scripts/build_dust_source_fields.py` takes
    only the grid, so both stay on `bootstrap_climatology_path`.

    THE GRAPH EDGE IS UNCHANGED, and it is a different statement. A step here
    still declares `needs: bootstrap_climatology` in `config/pipeline.yaml`,
    because that edge says which artifact must EXIST before the step can run and
    the bootstrap is the one that must: the step runs on a first pass, when no
    baseline exists at all. The edge is the ordering constraint; this resolver is
    what the step should READ. `scripts/smoke_test.py` holds both.
    """
    import yaml
    project = Path(root) if root is not None else PROJECT_ROOT
    config = yaml.safe_load(
        (project / "config" / "planet.yaml").read_text(encoding="utf-8"))
    if override is not None:
        path = Path(override).resolve()
        stage = climatology_stage(path, project)
    else:
        declared = config.get("baseline_climatology")
        stage = BASELINE
        if not declared:
            declared = config.get("bootstrap_climatology")
            stage = BOOTSTRAP
        if not declared:
            raise SystemExit(
                "config/planet.yaml names neither `baseline_climatology` nor "
                "`bootstrap_climatology`. This step reads the best available "
                "of the two and there is deliberately no third option: name "
                "one there or pass --climatology.")
        path = project / declared
    # Only when it is there, for the reason `climatology_path` gives: a file
    # that does not exist is the caller's error to report, and several callers
    # say something more useful about it than this could.
    if path.is_file():
        require_configured_grid(path, config)
    return Climatology(path, stage)


def require_configured_grid(climatology: Path, cfg: dict | None = None,
                            root: Path | None = None) -> None:
    """Refuse a climatology whose grid is not the configured rung's.

    The rung is a property of the run that produced the file and appears
    nowhere in its name, so nothing stopped a T21 climatology from driving a
    T42 configuration. What that produces is not an error but a plausible
    number computed on the wrong world, which is this project's most expensive
    failure shape: `check_consistency.py` built its land-sea mask this way and
    could therefore pass a T21 coupling for a T42 run.

    `lib/rungs.py:model_grid` supplies the configured dimensions, and refuses a
    `resolution` paired with another grid's `latitudes` before this compares
    anything, so a stale config cannot be what the climatology is judged
    against.
    """
    import yaml
    if cfg is None:
        project = Path(root) if root is not None else Path(__file__).resolve().parents[1]
        cfg = yaml.safe_load(
            (project / "config" / "planet.yaml").read_text(encoding="utf-8"))
    import sys
    lib = str(Path(__file__).resolve().parent)
    if lib not in sys.path:
        sys.path.insert(0, lib)
    import rungs
    rung, nlat, nlon = rungs.model_grid(cfg)
    with Dataset(climatology) as ds:
        got = (len(ds.dimensions["lat"]), len(ds.dimensions["lon"]))
    if got != (nlat, nlon):
        raise SystemExit(
            f"{rel(climatology)} is {got[0]}x{got[1]} and config/planet.yaml "
            f"is {rung}, {nlat}x{nlon}. A climatology carries no rung in its "
            "name, so this is the only thing between a run at one resolution "
            "and a product computed on another. Name a climatology built at "
            f"{rung}, or move model.resolution to the rung this one is.")


def require_clean_io(climatology: Path) -> None:
    """Refuse a climatology built from PlaSim's low-I/O output path.

    Such a product is wrong in two ways at once, and both are artifacts of that
    path rather than properties of the binning. The first output record of every
    orbit carries no boundary layer, so bottom-level wind reads about 7.5x the
    other bins and humidity 27% low. And the binned `spd` is additionally
    vector-cancelled, because the model accumulates the components over each
    output interval before anything averages them.

    Measured on this run: under `NLOWIO = 1`, orbit 65 gives a first-bin wind
    ratio of 9.0 and a binned `spd` 1.14x the mean of instantaneous speeds. Under
    `NLOWIO = 0`, orbits 70 and 76 give first-bin ratios of 1.03 and 0.99 and
    binned/snapshot ratios of 1.002 and 0.999. **Both defects vanish together.**

    So this refuses rather than corrects. Correcting was the earlier approach and
    the corrections are themselves wrong once the defect is gone: reading wind
    from the snapshot product costs a factor of 32 in samples for nothing, and
    scaling the binned wind up by 1.55 would be a 55% error on clean output.

    A product with no `low_io` attribute predates the stamp and is assumed
    tainted, which is the honest default -- every climatology built before
    2026-08-17 was.
    """
    with Dataset(climatology) as ds:
        low_io = getattr(ds, "low_io", None)
    if low_io is None or int(low_io) != 0:
        raise SystemExit(
            f"{climatology} was built from orbits run with PlaSim's low-I/O "
            "accumulation (or predates the stamp), so its wind and humidity are "
            "wrong. Rebuild it from NLOWIO = 0 orbits; see "
            "exoplasim/notes/first-output-bin.md.")


def snapshot_beside(regular: Path) -> Path:
    """The instantaneous-sample product beside a given regular climatology.

    Still produced and still useful -- it carries orbital phase, which the binned
    product does not -- but no longer a workaround for anything. See
    `require_clean_io`.

    IT TAKES THE REGULAR PRODUCT RATHER THAN RESOLVING ONE. This resolved the
    baseline itself, which made it wrong twice over. `build_dust.py` and
    `build_sea_salt.py` are its only callers and both are `bootstrap_climatology`
    steps, so the snapshot it handed them came from the other run; and a caller
    passing `--climatology` got its regular field overridden and its wind tail
    still read from the configured one, which is half an escape hatch. Both
    callers now pass the regular product they actually resolved, so the pair
    cannot come from two different runs and one flag moves both.
    """
    if "_regular_climatology.nc" not in regular.name:
        # Without this the substitution is a no-op and the "snapshot" returned
        # is the binned product itself, which is exactly the field the caller
        # asked for a snapshot INSTEAD of: a 12-bin mean has averaged away the
        # wind tail the fit is measuring.
        raise SystemExit(
            f"{rel(regular)} is not named `<label>_regular_climatology.nc`, so "
            "the snapshot beside it cannot be derived from its name. Name the "
            "sample file explicitly.")
    snapshot = regular.with_name(
        regular.name.replace("_regular_climatology.nc", "_snapshot_climatology.nc"))
    if not snapshot.is_file():
        raise FileNotFoundError(f"no snapshot climatology at {snapshot}")
    return snapshot
