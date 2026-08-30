"""Read BIO-11's one rootable-surface artifact without spatial fallback."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

import builds
import rungs
from gridding import require_same_rows


def rootable_path(config: dict) -> Path:
    rung = rungs.model_grid(config)[0]
    return builds.component_data("biosphere", config) / f"rootable_fraction_{rung}.nc"


def read_rootable_partition(config: dict, lat, lon, path: Path | None = None,
                            land=None) -> tuple[dict[str, np.ndarray], dict]:
    """Return rootable, solved-water and dry-barren land-relative fractions."""
    from netCDF4 import Dataset

    target = Path(path) if path is not None else rootable_path(config)
    if not target.is_file():
        raise SystemExit(
            f"{target} does not exist. Build BIO-11's rootable_fraction for "
            "this build and rung; barren land and solved water are not zero.")
    with Dataset(target) as ds:
        build = str(getattr(ds, "vesper_source_build", ""))
        if build != str(config.get("source_build", "")):
            raise SystemExit(
                f"{target} belongs to source build {build!r}, not "
                f"{config.get('source_build')!r}")
        file_lat = np.asarray(ds["lat"][:], dtype=float)
        file_lon = np.asarray(ds["lon"][:], dtype=float)
        fields = {
            "rootable": np.asarray(ds["f_rootable"][:], dtype=float),
            "water": np.asarray(ds["f_nonrootable_water"][:], dtype=float),
            "barren": np.asarray(ds["f_nonrootable_barren"][:], dtype=float),
        }
        identity = str(getattr(ds, "vesper_spatial_contract_identity", ""))
    require_same_rows(file_lat, lat, what=f"{target}'s rows and its consumer")
    shape = (len(lat), len(lon))
    if file_lon.size != len(lon) or any(field.shape != shape
                                        for field in fields.values()):
        raise SystemExit(
            f"{target}'s partition does not have consumer shape {shape}. "
            "Equal-looking coordinates do not repair "
            "a support mismatch.")
    finite = np.logical_and.reduce([np.isfinite(field) for field in fields.values()])
    for name, field in fields.items():
        if np.any(finite & ((field < -1e-7) | (field > 1.0 + 1e-7))):
            raise SystemExit(f"{target} carries {name} outside [0,1]")
    if np.any(finite & (np.abs(sum(fields.values()) - 1.0) > 2e-6)):
        raise SystemExit(f"{target}'s rootable/water/barren partition does not sum to one")
    if land is not None:
        land = np.asarray(land, dtype=bool)
        if land.shape != shape:
            raise SystemExit("the consumer's land mask does not match the rootable grid")
        if np.any(land & ~finite):
            raise SystemExit(
                f"{target} has {int(np.sum(land & ~finite))} model-land cells "
                "with no native-mesh support. SPAT-5 must name their ownership; "
                "BIO-11 does not fill them from a nearest region.")
    return {name: np.clip(field, 0.0, 1.0)
            for name, field in fields.items()}, {
        "path": str(target.relative_to(builds.PROJECT_ROOT)),
        "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "spatial_contract_identity": identity,
    }


def read_rootable(config: dict, lat, lon, path: Path | None = None,
                  land=None) -> tuple[np.ndarray, dict]:
    """Return the land-relative rootable fraction and its provenance digest."""
    partition, provenance = read_rootable_partition(
        config, lat, lon, path, land=land)
    return partition["rootable"], provenance
