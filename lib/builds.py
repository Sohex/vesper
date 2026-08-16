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

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT_ROOT / "source"
CONFIG = PROJECT_ROOT / "config" / "planet.yaml"


def _config(config: dict | None = None) -> dict:
    return config if config is not None else yaml.safe_load(
        CONFIG.read_text(encoding="utf-8"))


def available() -> list[str]:
    return sorted(p.name for p in SOURCE.iterdir()
                  if p.is_dir() and (p / "exoplasim-T42").is_dir())


def build_root(config: dict | None = None) -> Path:
    """Directory of the configured build."""
    name = _config(config).get("source_build")
    if not name:
        raise RuntimeError(
            f"config/planet.yaml has no source_build. Available: {available()}")
    path = SOURCE / name
    if not (path / "exoplasim-T42").is_dir():
        raise RuntimeError(
            f"source_build {name!r} is not a build directory. Available: {available()}")
    return path


def mesh_export(config: dict | None = None) -> Path:
    """Export carrying `raw/`. The mesh is identical across grids in a build."""
    return build_root(config) / "exoplasim-T42"


def grid_export(config: dict | None = None, resolution: str | None = None) -> Path:
    """Export whose grid matches the configured model resolution."""
    cfg = _config(config)
    res = (resolution or str(cfg["model"]["resolution"])).upper()
    return build_root(cfg) / f"exoplasim-{res}"


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

    Falls back to the flat directory when no per-build one exists, so older
    layouts keep working; pass `strict=True` to refuse that fallback.
    """
    cfg = _config(config)
    name = str(cfg.get("source_build", ""))
    root = PROJECT_ROOT / component / "data"
    named = root / name
    if named.is_dir():
        return named
    if strict:
        raise RuntimeError(
            f"no per-build data at {named}; build it for {name!r} rather than "
            f"falling back to {root}, which holds whichever build was active "
            "when it was last written")
    return root


def terrain_hash(config: dict | None = None) -> str:
    """finalElevation of the configured build. The identity a result belongs to."""
    import json
    manifest = mesh_export(config) / "manifest.json"
    return json.loads(manifest.read_text(encoding="utf-8"))["hashes"]["finalElevation"]


def soilmap(config: dict | None = None) -> Path:
    """Pedology's soil map for the configured build.

    Soil texture derives from lithology, so this is per-build like everything
    else downstream of a terrain. It was a flat `pedology/data/soilmap.txt` with
    no build in the name, which is the same trap that caught four hydrography
    scripts with one fewer script in it -- and this one feeds LPJ-GUESS, which
    would have grown a biosphere on another planet's soil without complaint.
    """
    return component_data("pedology", config) / "soilmap.txt"


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
            return tag
    raise RuntimeError(
        f"cannot read a resolution from grid directory {name!r}; expected "
        "something like 'exoplasim-T42'")
