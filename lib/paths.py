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


def snapshot_sibling(regular: Path) -> Path:
    """The instantaneous-sample product beside a binned climatology.

    Two fields have to be read from here rather than from the binned product:
    wind speed and specific humidity.

    Wind wants it on the physics. A turbulent flux is driven by the mean of the
    SPEED, and the binned `spd` is much closer to the speed of a time-mean
    vector, which cancels where the wind reverses: 5.14 m/s per orbit against
    7.50 from instantaneous samples over the same orbits.

    Humidity wants it because of a defect. PlaSim's low-I/O path wrote a corrupt
    first output record per model call -- in these same two fields -- and a
    binned mean carries it. Runs from 2026-08-17 set `NLOWIO = 0` and no longer
    produce it, but every climatology built before that does; see
    `exoplasim/notes/first-output-bin.md`. The two errors happened to run
    opposite ways, so the uncorrected binned wind was within 3.1% of this and
    dropping the bad record alone would have been 33% low.

    Taken from the given path rather than from config, so a `--climatology`
    override stays self-consistent.
    """
    regular = Path(regular)
    stem = regular.name
    if "_regular_climatology.nc" not in stem:
        raise ValueError(
            f"{regular} is not a *_regular_climatology.nc, so its snapshot "
            "sibling cannot be named; pass the binned product")
    snapshot = regular.with_name(
        stem.replace("_regular_climatology.nc", "_snapshot_climatology.nc"))
    if not snapshot.is_file():
        raise FileNotFoundError(
            f"no snapshot climatology at {snapshot}. Wind and humidity are read "
            "from instantaneous samples, not from the binned mean; rebuild the "
            "climatology so both products exist.")
    return snapshot


def snapshot_climatology_path(name: str | None = None,
                              root: Path | None = None) -> Path:
    """`snapshot_sibling` of the configured climatology."""
    return snapshot_sibling(climatology_path(name, root=root))
