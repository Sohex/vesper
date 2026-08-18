"""Path display helpers shared across components.

One function, and it exists because of a specific recurring bug.

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
`notes/failure-modes.md`.

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
    re-baseline ran them for the first time in months.

    Returns the FILE. Callers must not rebuild the name from a directory: the
    label is chosen per product, so a bootstrap climatology is
    `bootstrap_regular_climatology.nc` and reconstructing
    `baseline_regular_climatology.nc` beside it finds nothing.

    `name` still accepts a directory under `exoplasim/analysis/` for the old
    layout, so existing callers that pass one keep working.
    """
    import yaml
    project = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    if name is not None:
        return (project / "exoplasim" / "analysis" / name
                / "baseline_regular_climatology.nc")
    config = yaml.safe_load(
        (project / "config" / "planet.yaml").read_text(encoding="utf-8"))
    declared = config.get("baseline_climatology")
    if not declared:
        raise SystemExit(
            "config/planet.yaml has no `baseline_climatology`. Name one there "
            "or pass --climatology; there is deliberately no fallback.")
    return project / declared


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


def snapshot_climatology_path(name: str | None = None,
                              root: Path | None = None) -> Path:
    """The instantaneous-sample product beside the configured climatology.

    Still produced and still useful -- it carries orbital phase, which the binned
    product does not -- but no longer a workaround for anything. See
    `require_clean_io`.
    """
    regular = climatology_path(name, root=root)
    snapshot = regular.with_name(
        regular.name.replace("_regular_climatology.nc", "_snapshot_climatology.nc"))
    if not snapshot.is_file():
        raise FileNotFoundError(f"no snapshot climatology at {snapshot}")
    return snapshot
