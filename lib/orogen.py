"""Reader for World Orogen data exports.

The export is a manifest plus flat little-endian binaries. This wraps that in
something that reads a field by name, keeps the mesh adjacency in CSR form, and
refuses to silently hand back a field from a build other than the one asked for.

Project-level: every component reads the export through this, so conventions
like the terrain-hash allowlist and the two land definitions are enforced once.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
import json
from pathlib import Path

import numpy as np

from builds import mesh_export as _configured_mesh_export

# Elevation conversion changed meaning in this build: below-sea-level land used
# to take the bathymetric branch and read ten times too deep. Refuse to run
# against an export that predates the fix rather than producing quiet nonsense.
_KNOWN_TERRAIN_HASHES = {
    "821aa71b37a7beda0b59398c7f005b91531050000ca46660d0f724cdb3f401a3":
        "2026-08 build: over-erosion fixed, sub-sea-level land at 1.0 km/unit",
    # WITHDRAWN. The terrain move in this build was not an improvement: a basin
    # protection floor was captured pre-erosion and stored absolute, letting the
    # carve take back divides that erosion had raised, with the same array used
    # as the assertion baseline so it never tripped. Reverted upstream; the
    # regenerated export is bit-identical to 821aa71b again. Kept here only so an
    # export from that window is recognised rather than silently accepted.
    "27b7479aa486f5dacebccb0c638ff839a60a98e617c2437229600ef0bacf32ec":
        "WITHDRAWN 2026-08 build: basin-protection floor drift, superseded by "
        "821aa71b. Do not use.",
    # Evaporite split into salt crust (0.50) and playa clastics (0.30). The two
    # have different erodibility, 3.50 against 2.80, so stream power sees a
    # different surface and the terrain moved on both planets. The basin
    # catalogue hash did not move, because detection still runs on the
    # pre-conditioning surface, so per-basin work computed against 821aa71b
    # still resolves.
    "26fc76914da14289ff26f15a130192bd84d59031098569adb66186ffdabb28b7":
        "2026-08 precarve-zoned: threshold selection, crust/fill lithology split",
    # SUPERSEDED. The verdict this build applied was decided on climate read
    # 180 degrees out: the coupling matrix numbers its columns on the Orogen
    # grid's -180..180 and an ExoPlaSim climatology on 0..360, so every basin
    # integrated its antipode's rainfall. Of its 1,522 carves, 850 are not
    # justified by the climate meant to justify them and cannot be un-cut.
    # Kept registered so results computed from it stay readable and traceable.
    "3899a0c57d1eee2f47ba9054c218a171a7aa4532e2437c9104070c2c3dfaece6":
        "2026-08 carved-zoned: iteration-1 carve verdict applied, crust/fill "
        "split. SUPERSEDED by 010f2143 -- its verdict used antipodal climate.",
    # SUPERSEDED. The same first pass recomputed after the longitude fix, so its
    # verdict is the right one: 1,089 carve, 170 marginal, 2,370 preserved. The
    # lithology under it was not. Correct as a build of the model as it then
    # stood, but that
    # model decided which deposit sat on top by the order the branches were
    # typed in, so closed-basin fill lost to three separate rules: fold-belt
    # exhumation inside orogens, and basins on oceanic crust or flood basalt
    # never reaching the endorheic branch at all. 203 preserved basins had no
    # fill cell anywhere and were reaching ExoPlaSim as vegetated land.
    "010f214397338008ae28de1d8ecce93486d7ce8232e1cfb554006cfebbc14b6f":
        "2026-08 carved-zoned-v2: corrected iteration-1 carve verdict, "
        "crust/fill split. SUPERSEDED by 5bed5549 -- basin fill was being "
        "overwritten by cover-chain branch order.",
    # Same verdict, same catalogue, cover chain now a declared table walked in
    # order with closed-basin fill first. 99.77% of preserved-basin area carries
    # fill against 74.6%, and no basin is left without any. The first build on
    # which terrain, lithology and the drainage verdict are mutually consistent.
    "5bed5549315da14b22275fea51a0b6f5b34d79cdf2237c9380e8471e0b431c78":
        "2026-08 carved-zoned-v4: corrected verdict, cover-chain fix, "
        "crust/fill split",
}

# Basin ids are computed on the pre-conditioning surface, so they survive a
# carve iteration. This hash is the check: if it matches, a carve verdict
# computed against an earlier export still refers to the same basins.
_KNOWN_CATALOGUE_HASHES = {
    "2d1f8e57c26b60c608deaa62bd3c44b7da04a4ee5bd518249447630095d96a98":
        "2026-08 catalogue, 3629 preserved from 81904 detected",
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
        self.root = Path(root) if root is not None else _configured_mesh_export()
        manifest_path = self.root / "manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"No manifest at {manifest_path}")
        self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.terrain_hash = self.manifest["hashes"]["finalElevation"]
        self.catalogue_hash = self.manifest["hashes"].get("basinCatalogue")
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
            "catalogue_hash": self.catalogue_hash,
            "catalogue_known": self.catalogue_hash in _KNOWN_CATALOGUE_HASHES,
            "drainage_hypothesis": self.manifest["basins"].get("drainageHypothesis"),
            "selection_source": self.manifest["basins"].get("selectionSource"),
            "seed": self.manifest["seed"],
            "num_regions": self.n_regions,
            "planet_radius_km": self.radius_km,
            "surface_area_km2": self.surface_area_km2,
        }
