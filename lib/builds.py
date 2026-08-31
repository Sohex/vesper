"""Resolve which World Orogen build a component should read.

`source/` holds one namespaced directory per build, because builds now multiply
faster than they can be swapped in place: a carve verdict, a lithology split and
any future iteration each produce a new terrain, and each one has downstream
artifacts computed against it. Overwriting `source/` would make a result's
provenance depend on when it was computed rather than on what it was computed
from.

Every build carries `manifest.hashes.finalElevation`, so the build a result came
from is recoverable from the result itself. The name is for humans; the hash is
the identity.
"""

from __future__ import annotations

from pathlib import Path

import yaml

# `maps/frames.py` reaches this module as `from lib import builds`, where lib is
# not on sys.path as a directory, so a bare sibling import fails there. Anchor
# it the way every script anchors its own imports.
import sys as _sys
if str(Path(__file__).resolve().parent) not in _sys.path:
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
import rungs  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT_ROOT / "source"
CONFIG = PROJECT_ROOT / "config" / "planet.yaml"


def _config(config: dict | None = None) -> dict:
    return config if config is not None else yaml.safe_load(
        CONFIG.read_text(encoding="utf-8"))


def is_build(path: Path) -> bool:
    """A build directory is one holding at least one grid export.

    Tested by shape rather than by naming a rung: `exoplasim-T42` was the
    literal here, so a build that had every other grid and not that one was not
    a build. SPAT-2.
    """
    return Path(path).is_dir() and any(Path(path).glob("exoplasim-T*"))


def available() -> list[str]:
    return sorted(p.name for p in SOURCE.iterdir() if is_build(p))


# What one finished grid export holds, from the recipe in `source/README.md`:
# the field catalogue, the CF netCDF, the binary field directory and the
# per-export README. `raw/` is deliberately not here -- exactly one export per
# build carries it and `mesh_export_of` is the check for that.
_EXPORT_CONTENTS = ("manifest.json", "planet.nc", "README.txt")


def incomplete_exports(root: Path) -> list[str]:
    """Why a build directory is not a finished set of exports; empty if it is.

    DERIVED RATHER THAN LISTED, and that is the whole point of it.
    `scripts/archive_builds.py` held the completeness test as five literal
    directory names -- T21, T42, T63, T85 and `grid-512x256`. The recipe in
    `source/README.md` has since retired T63 and added T127 and T170, so the
    ACTIVE build matched none of that list and could never have been archived
    however superseded it became, while the older build passed on a T63 nothing
    emits any more.

    A build cannot be asked which grids the recipe of its own day emitted,
    because nothing records that. What it CAN be asked is whether each export
    it carries is finished, and that is the property "still being generated"
    actually names: an export is written file by file, so a directory missing
    its manifest or its netCDF is one being written now. A Gaussian export must
    also name a rung in `lib/rungs.py`, so a misspelled or half-named directory
    is refused rather than counted.

    AN EXPORT DIRECTORY CARRYING NO DATA IS ABSENT, NOT UNFINISHED. The
    exporter writes `README.txt` from a template that needs nothing computed,
    so a directory holding that and nothing else is a grid whose export never
    produced anything -- indistinguishable from a grid the recipe of that day
    never emitted, and with no payload to lose either way.
    `source/precarve-craton/exoplasim-T127` is one, and calling it unfinished
    would block a build that is otherwise written through. What blocks is a
    directory with DATA in it and a piece missing, which is what an export
    being written now looks like.

    This is not the only refusal in front of a deletion, and it was never meant
    to be. `lib/orogen.py`'s registry rejects any terrain hash a human has not
    registered, `mesh_export_of` rejects a build with no `raw/` carrier or with
    two, and the caller skips the configured build, so a build being generated
    right now is refused several times over before this is consulted.
    """
    root = Path(root)
    reasons = []
    exports = sorted(d for d in root.iterdir()
                     if d.is_dir() and (d.name.startswith("exoplasim-")
                                        or d.name.startswith("grid-")))
    if not exports:
        return [f"{root.name} carries no grid export"]
    for d in exports:
        if not any(f.name != "README.txt" for f in d.iterdir()):
            continue  # no data written: absent, not unfinished
        if d.name.startswith("exoplasim-"):
            try:
                resolution_of(d)
            except RuntimeError as exc:
                reasons.append(str(exc))
                continue
        missing = [f for f in _EXPORT_CONTENTS if not (d / f).is_file()]
        if not (d / "grid").is_dir() or not any((d / "grid").iterdir()):
            missing.append("grid/")
        if missing:
            reasons.append(f"{d.name} is missing {', '.join(missing)}")
    return reasons


def build_root(config: dict | None = None) -> Path:
    """Directory of the configured build."""
    name = _config(config).get("source_build")
    if not name:
        raise RuntimeError(
            f"config/planet.yaml has no source_build. Available: {available()}")
    path = SOURCE / name
    if not is_build(path):
        raise RuntimeError(
            f"source_build {name!r} is not a build directory. Available: {available()}")
    return path


def mesh_export(config: dict | None = None) -> Path:
    """The export directory carrying `raw/`: the build's MESH STORAGE CARRIER.

    The native mesh is written once per build rather than once per grid,
    because it is the same mesh for every grid. It happens to be stored under
    `exoplasim-T42`, and that name has been read as "the T42 export" by
    consumers that meant "the export the raw mesh is in" -- which is how a T42
    literal ended up in code paths that have nothing to do with T42.

    So the carrier is IDENTIFIED rather than named: whichever export holds
    `raw/`. That is also a check, because a build with none and a build with
    two are both wrong and were both previously invisible.
    """
    return mesh_export_of(build_root(config))


def mesh_export_of(root: Path) -> Path:
    """The mesh storage carrier under one build directory. See `mesh_export`."""
    root = Path(root)
    carriers = sorted(d for d in root.glob("exoplasim-T*") if (d / "raw").is_dir())
    if not carriers:
        raise RuntimeError(
            f"no export under {root} carries raw/. The native mesh is stored "
            "once per build and every mesh consumer reads it from there; "
            "without it nothing can integrate from the mesh.")
    if len(carriers) > 1:
        raise RuntimeError(
            f"{len(carriers)} exports under {root} carry raw/ "
            f"({', '.join(d.name for d in carriers)}). The mesh is one thing "
            "per build and two copies are two chances to read the older one.")
    return carriers[0]


def mesh_export_at(path: Path) -> Path:
    """The mesh carrier for a directory that is either the carrier or the build.

    A caller naming an export on a command line has one of two directories:
    the one World Orogen's `--out` wrote, which holds `raw/` directly, or a
    build directory under `source/`, whose carrier is a subdirectory. Both are
    real and a measurement script is handed both, so resolving them here is
    what lets a usage example name the BUILD. That is the whole point: naming
    the carrier means spelling `exoplasim-T42` in a path that has nothing to do
    with T42, which is the literal `mesh_export` above exists to keep out.
    """
    path = Path(path)
    return path if (path / "raw").is_dir() else mesh_export_of(path)


def grid_export(config: dict | None = None, resolution: str | None = None,
                *, build: str | None = None) -> Path:
    """Export whose grid matches the configured model resolution.

    `build` names another build's directory, for a consumer processing a build
    that is not the active one. The resolution is checked against the ladder
    registry, so a typo asks for a grid that cannot exist rather than a
    directory that silently is not there.
    """
    cfg = _config(config)
    res = (resolution or str(cfg["model"]["resolution"])).upper()
    rungs.geometry(res)                      # refuses anything off the ladder
    root = (SOURCE / build) if build else build_root(cfg)
    return root / f"exoplasim-{res}"


def component_data(component: str, config: dict | None = None,
                   *, strict: bool = False) -> Path:
    """Per-build data directory for a component, e.g. `hydrography/data/<build>`.

    Drainage, basins and coupling matrices are properties of a terrain, so they
    are namespaced by build exactly as `source/` is. Four scripts were written
    against the flat `<component>/data/` before that was true, and each one would
    happily pair one terrain's rows with another's columns: `world_state.py`
    reported 2,107 basins against a build that had 2,540, and `carve_verdict.py`
    raised an IndexError only because the counts happened to differ. Had they
    matched, it would have computed a wrong answer in silence.

    **Never falls back to the flat directory.** It used to, "so older layouts
    keep working", and that fallback is the original bug wearing the shape of its
    fix: it turns a missing input into a silent read of whichever build happened
    to be active when the flat directory was last written. Returning a path that
    does not exist is strictly better -- the caller fails on a missing file,
    naming it, instead of succeeding against the wrong terrain.

    `strict=True` additionally refuses to return the path at all, raising here
    rather than at the caller's first read. Use it where the failure should
    surface at argument-parse time, which is most places.
    """
    cfg = _config(config)
    name = str(cfg.get("source_build", ""))
    named = PROJECT_ROOT / component / "data" / name
    if strict and not named.is_dir():
        raise RuntimeError(
            f"no per-build data at {named}; build it for {name!r}. There is no "
            "fallback: the flat directory holds whichever build was active when "
            "it was last written, and reading it would pair one terrain's rows "
            "with another's columns.")
    return named


def terrain_hash(config: dict | None = None) -> str:
    """finalElevation of the configured build. The identity a result belongs to."""
    import json
    manifest = mesh_export(config) / "manifest.json"
    return json.loads(manifest.read_text(encoding="utf-8"))["hashes"]["finalElevation"]


def soilmap(config: dict | None = None) -> Path:
    """Pedology's soil map for the configured build AND RUNG.

    Soil texture derives from lithology, so this is per-build like everything
    else downstream of a terrain. It was a flat `pedology/data/soilmap.txt` with
    no build in the name, which is the same trap that caught four hydrography
    scripts with one fewer script in it -- and this one feeds LPJ-GUESS, which
    would have grown a biosphere on another planet's soil without complaint.

    The rung is in the name for the same reason the build is. The map is one
    row per LAND CELL of a climate grid -- 1,019 rows at T21 -- so it is a
    property of the support as much as of the terrain, and a name without the
    rung in it means a run at another rung overwrites it and every consumer
    that has not been re-run reads another grid's soil. SPAT-2.
    """
    res = rungs.model_grid(_config(config))[0]
    return component_data("pedology", config) / f"soilmap_{res}.txt"


def land_column_states(config: dict | None = None) -> Path:
    """The land column property contract's per-cell states, for the build AND RUNG.

    What `pedology/scripts/land_column_properties.py --emit` writes and what
    every consumer of the simulated soil's hydraulic states reads. Per build and
    per rung for exactly the reasons `soilmap` is: the states are derived from
    the soil map's texture and regolith on a climate grid's land cells, so a
    name without both in it lets one rung's states reach another rung's model.
    """
    res = rungs.model_grid(_config(config))[0]
    return component_data("pedology", config) / f"land_column_states_{res}.txt"


def resolution_of(grid_dir: Path) -> str:
    """Resolution implied by a grid export directory, e.g. `exoplasim-T42` -> T42.

    The grid IS the resolution, so anything writing a surface field should take
    its resolution from the grid it integrated onto rather than from the config.
    The two are the same whenever `--grid` is left at its default and differ
    exactly when someone builds for another resolution -- which is the case that
    matters, because the output path and the SRA filename both encode it. Reading
    the resolution from config while reading the data from `--grid` writes T21
    fields into files named T42, in the T42 directory, with no error.
    """
    name = Path(grid_dir).name
    if "-" in name:
        tag = name.rsplit("-", 1)[1].upper()
        if tag.startswith("T") and tag[1:].isdigit():
            # AGAINST THE LADDER, not against the spelling. T-and-digits is a
            # shape, not a rung: `exoplasim-T99` parsed cleanly and three
            # surface-field generators wrote a whole family under it, so
            # `inputs/t99/orogen_T99_surf_0129.sra` existed for a grid that
            # cannot be built. `grid_export` above already refuses this way and
            # this is the same rule. world-yop.
            rungs.geometry(tag)
            return tag
    raise RuntimeError(
        f"cannot read a resolution from grid directory {name!r}; expected "
        "something like 'exoplasim-T42'")
