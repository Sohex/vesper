"""Reader for World Orogen data exports.

The export is a manifest plus flat little-endian binaries. This wraps that in
something that reads a field by name, keeps the mesh adjacency in CSR form, and
refuses to silently hand back a field from a build other than the one asked for.

Nothing here is specific to hydrography; other components should use it too.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
import json
from pathlib import Path

import numpy as np

from _paths import MESH_EXPORT

# Elevation conversion changed meaning in this build: below-sea-level land used
# to take the bathymetric branch and read ten times too deep. Refuse to run
# against an export that predates the fix rather than producing quiet nonsense.
_KNOWN_TERRAIN_HASHES = {
    "821aa71b37a7beda0b59398c7f005b91531050000ca46660d0f724cdb3f401a3":
        "2026-08 build: over-erosion fixed, sub-sea-level land at 1.0 km/unit",
}

OCEAN, LAND, INLAND_WATER = 0, 1, 2


@dataclass(frozen=True)
class Basin:
    """One preserved closed basin, as the catalogue describes it.

    Elevations carry two conventions. Unsuffixed keys are the generator's model
    parameter; `*Km` keys are physical kilometres converted through the land
    branch. Only the physical ones are used here.

    `hypsometry` in the catalogue is measured on the NATURAL (pre-conditioning)
    terrain, so it overstates what the finished terrain can hold. Use
    `build_hydrography.py`, which recomputes it on the final terrain.
    """

    index: int
    id: str
    sink: int                 # region index of the sink, on the final terrain
    sink_elevation_km: float
    spill_elevation_km: float
    spills_into: str
    natural_area_km2: float
    natural_volume_km3: float
    final_flooded_area_km2: float
    final_volume_km3: float
    catchment_area_km2: float

    @property
    def natural_hypsometry_is_usable(self) -> bool:
        """The natural curve is only safe where erosion barely touched the rim."""
        return self.final_volume_km3 >= 0.95 * self.natural_volume_km3


class Export:
    """A World Orogen export directory."""

    def __init__(self, root: Path | None = None, *, require_known_build: bool = True):
        self.root = Path(root or MESH_EXPORT)
        manifest_path = self.root / "manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"No manifest at {manifest_path}")
        self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.terrain_hash = self.manifest["hashes"]["finalElevation"]
        if require_known_build and self.terrain_hash not in _KNOWN_TERRAIN_HASHES:
            raise RuntimeError(
                f"Export terrain hash {self.terrain_hash[:16]} is not a build this "
                "code has been checked against. Elevation and basin-catalogue "
                "conventions have changed between builds; verify before overriding "
                "with require_known_build=False."
            )

        if self.manifest.get("raw") is None:
            raise RuntimeError(
                f"{self.root} has no raw/ mesh; hydrography needs the native mesh"
            )
        self._raw_fields = {f["name"]: f for f in self.manifest["raw"]["fields"]}
        self._cache: dict[str, np.ndarray] = {}

    # -- scalars ---------------------------------------------------------

    @property
    def n_regions(self) -> int:
        return int(self.manifest["numRegions"])

    @property
    def radius_km(self) -> float:
        return float(self.manifest["planet"]["radiusKm"])

    @property
    def surface_area_km2(self) -> float:
        return float(self.manifest["planet"]["surfaceAreaKm2"])

    # -- fields ----------------------------------------------------------

    def field(self, name: str) -> np.ndarray:
        """Read a raw-mesh field by name, memoised."""
        if name in self._cache:
            return self._cache[name]
        spec = self._raw_fields.get(name)
        if spec is None:
            raise KeyError(
                f"{name!r} is not in this export. Available: "
                f"{', '.join(sorted(self._raw_fields)[:12])}..."
            )
        a = np.fromfile(self.root / spec["path"], dtype=spec["dtype"])
        expected = int(np.prod(spec["shape"]))
        if a.size != expected:
            raise RuntimeError(f"{name}: read {a.size} values, manifest says {expected}")
        self._cache[name] = a
        return a

    def __getattr__(self, name: str) -> np.ndarray:
        # Convenience: export.elevation_km rather than export.field("elevation_km").
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self.field(name)
        except KeyError as exc:
            raise AttributeError(name) from exc

    # -- mesh ------------------------------------------------------------

    @cached_property
    def adjacency(self) -> tuple[np.ndarray, np.ndarray]:
        """CSR neighbour structure as (offsets, neighbours)."""
        mesh = self.manifest["raw"]["mesh"]
        off = np.fromfile(self.root / mesh["adjOffset"]["path"], dtype=mesh["adjOffset"]["dtype"])
        lst = np.fromfile(self.root / mesh["adjList"]["path"], dtype=mesh["adjList"]["dtype"])
        if off.size != self.n_regions + 1:
            raise RuntimeError(f"adjOffset has {off.size} entries, expected {self.n_regions + 1}")
        return off, lst

    def neighbours(self, region: int) -> np.ndarray:
        off, lst = self.adjacency
        return lst[off[region]:off[region + 1]]

    # -- basins ----------------------------------------------------------

    @cached_property
    def basins(self) -> list[Basin]:
        out = []
        for i, b in enumerate(self.manifest["basins"]["preserved"]):
            fp = b["finalPreserved"]
            fc = b["finalCatchment"]
            out.append(Basin(
                index=i,
                id=b["id"],
                sink=int(fp["sink"]),
                sink_elevation_km=float(fp["sinkElevationKm"]),
                spill_elevation_km=float(fp["spillElevationKm"]),
                spills_into=b["spillsInto"],
                natural_area_km2=float(b["areaKm2"]),
                natural_volume_km3=float(b["volumeKm3"]),
                final_flooded_area_km2=float(fp["floodedAreaKm2"]),
                final_volume_km3=float(fp["volumeKm3"]),
                catchment_area_km2=float(fc["areaKm2"]),
            ))
        return out

    # -- provenance ------------------------------------------------------

    def provenance(self) -> dict:
        return {
            "export_root": str(self.root),
            "terrain_hash": self.terrain_hash,
            "terrain_build": _KNOWN_TERRAIN_HASHES.get(self.terrain_hash, "unrecognised"),
            "seed": self.manifest["seed"],
            "num_regions": self.n_regions,
            "planet_radius_km": self.radius_km,
            "surface_area_km2": self.surface_area_km2,
        }
