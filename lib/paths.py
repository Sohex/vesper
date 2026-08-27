"""Path helpers shared across components, and the guards on the climatology.

`rel` exists because of a specific recurring bug. There are TWO climatology
resolvers because there are two climatologies: `bootstrap_climatology_path` for
the run on terrain-only surface fields, which the derived fields are built FROM,
and `climatology_path` for the baseline run on those fields once they exist.
`config/pipeline.yaml` says which of the two each step needs, and
`scripts/smoke_test.py:check_climatology_needs_match_call_sites` holds the call
sites to it. `require_clean_io` and `require_configured_grid` are what both must
survive.

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

PROJECT_ROOT = Path(__file__).resolve().parents[1]


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

    A SECOND KEY BECAUSE THERE ARE TWO CLIMATOLOGIES AND SEVEN STEPS WANT THE
    OTHER ONE. `config/pipeline.yaml` already distinguishes them -- `surface_water`,
    `groundwater`, `soil`, `dust`, `sea_salt`, `volcanic_sulfate` and
    `vesper_header` all declare `needs: bootstrap_climatology` -- while
    `carve_verdict`, `ice_mask`, `lpj_driver` and `error_budget` declare
    `baseline_climatology`. Config carried one key for both, so every one of
    those eleven resolved to the same file and the graph's distinction reached
    nothing.

    THE ORDER IS WHY IT MATTERS. The bootstrap exists to produce the climatology
    the derived surface fields are built FROM, and the baseline is the run on
    those fields. A step that needs the bootstrap and reads the baseline is
    asking for an artifact that does not exist yet on a first pass, and on a
    later pass it silently reads a climate produced by the fields it is
    supposed to be producing.

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
    snapshot = regular.with_name(
        regular.name.replace("_regular_climatology.nc", "_snapshot_climatology.nc"))
    if not snapshot.is_file():
        raise FileNotFoundError(f"no snapshot climatology at {snapshot}")
    return snapshot
