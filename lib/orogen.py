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
        {"name": "precarve-unzoned", "note":
         "2026-08 build: over-erosion fixed, sub-sea-level land at 1.0 km/unit"},
    # WITHDRAWN. The terrain move in this build was not an improvement: a basin
    # protection floor was captured pre-erosion and stored absolute, letting the
    # carve take back divides that erosion had raised, with the same array used
    # as the assertion baseline so it never tripped. Reverted upstream; the
    # regenerated export is bit-identical to 821aa71b again. Kept here only so an
    # export from that window is recognised rather than silently accepted.
    "27b7479aa486f5dacebccb0c638ff839a60a98e617c2437229600ef0bacf32ec":
        {"name": "withdrawn-basin-floor-drift", "note":
         "WITHDRAWN 2026-08 build: basin-protection floor drift, superseded by "
        "821aa71b. Do not use."},
    # Evaporite split into salt crust (0.50) and playa clastics (0.30). The two
    # have different erodibility, 3.50 against 2.80, so stream power sees a
    # different surface and the terrain moved on both planets. The basin
    # catalogue hash did not move, because detection still runs on the
    # pre-conditioning surface, so per-basin work computed against 821aa71b
    # still resolves.
    "26fc76914da14289ff26f15a130192bd84d59031098569adb66186ffdabb28b7":
        {"name": "precarve-zoned", "note":
         "2026-08 precarve-zoned: threshold selection, crust/fill lithology split"},
    # SUPERSEDED. The verdict this build applied was decided on climate read
    # 180 degrees out: the coupling matrix numbers its columns on the Orogen
    # grid's -180..180 and an ExoPlaSim climatology on 0..360, so every basin
    # integrated its antipode's rainfall. Of its 1,522 carves, 850 are not
    # justified by the climate meant to justify them and cannot be un-cut.
    # Kept registered so results computed from it stay readable and traceable.
    "3899a0c57d1eee2f47ba9054c218a171a7aa4532e2437c9104070c2c3dfaece6":
        {"name": "carved-zoned", "note":
         "2026-08 carved-zoned: iteration-1 carve verdict applied, crust/fill "
        "split. SUPERSEDED by 010f2143 -- its verdict used antipodal climate."},
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
        {"name": "carved-zoned-v2", "note":
         "2026-08 corrected iteration-1 carve verdict, crust/fill split. "
         "SUPERSEDED by 5bed5549 -- basin fill was being overwritten by "
         "cover-chain branch order."},
    # SUPERSEDED. Cover chain became a declared table walked in order with
    # closed-basin fill first, which fixed the lithology. Its verdict integrated
    # ExoPlaSim's mrro as catchment runoff,
    # which is river-routed net divergence rather than local generation, so it
    # under-carved: 749 more basins overflow under precipitation minus
    # evaporation. Recoverable rather than wrong, since every basin it cut should
    # have been cut, and the carve list carried those forward.
    "5bed5549315da14b22275fea51a0b6f5b34d79cdf2237c9380e8471e0b431c78":
        {"name": "carved-zoned-v4", "note":
         "2026-08 corrected verdict, cover-chain fix, crust/fill split. "
         "SUPERSEDED by a4d204f6 -- its verdict read mrro as catchment runoff."},
    # Pass 2: 1,838 carved of 3,629, of which 1,089 carried forward from v4 and
    # 749 decided against P - E. Endorheic land falls 60.10% to 43.06% and
    # closed-basin fill 16.5% to 12.35%, because a net-divergence field
    # under-reports exactly the large catchments the criterion turns on.
    "a4d204f6e3e706597ff1470d064b457768243483ff884e9b574f055774de3ce9":
        {"name": "carved-zoned-v5", "note":
         "2026-08 iteration-2 verdict on P - E runoff, crust/fill split"},
    # Pre-carve base at the corrected gravity, 12.81 m/s2, and the first build
    # with both km bugs fixed: orog_* is real kilometres and the basin
    # catalogue's ...Km fields carry the relief scaling.
    #
    # THIS HASH IS THE SAME AT 10.1989 AND 12.81 -- verified upstream by
    # generating both. Gravity is not in any hash, so the allowlist alone cannot
    # tell two gravities apart and `Export` checks `planet.gravityMS2` against
    # the config separately. Do not treat a hash match as sufficient here.
    "974ceb78fccfd54c3a842cf4ae8d7ae79029816abd9b4fbf90409165090b0ec1":
        {"name": "precarve-zoned-g1281", "note":
         "2026-08 pre-carve base at g = 12.81, crust/fill split, orography and "
         "basin-km unit fixes"},
    # Same seed, same sliders, same gravity. What changed is lithology: two
    # erodibility values corrected against the literature (carbonate 1.30 ->
    # 0.45, schist 1.10 -> 0.45) and the class spread compressed to about 4x
    # with --lithology-strength 0.682, since Orogen ran 14x against Moosdorf's
    # global index of 3.2x and Zondervan's fluvially expressed contrast of ~4x.
    #
    # The land/sea mask is BIT-IDENTICAL to 974ceb78 -- zero cells of 2,500,001
    # differ -- and so are basinCatalogue, basinMembership and basinsPreserved,
    # so a carve verdict computed against the old export still refers to these
    # basins. Mean land elevation moves -3.0 m. What does move is finished basin
    # geometry: finalPreserved volumeKm3 by 8.1% in the median, which is what the
    # carve criterion tests, so hydrography and the verdict must be recomputed.
    "e2b510660d1cbe999d237734969e9ac9bf8636b6d85a04f513a4bc0d1f92d321":
        {"name": "precarve-zoned-g1281-erod4x", "note":
         "2026-08 pre-carve base at g = 12.81 with corrected carbonate/schist "
         "erodibility and the class spread compressed to ~4x "
         "(lithologyStrength 0.682)"},
    # Same seed, sliders and gravity again. What changed is that the arc and
    # forearc lithology rules FIRE AT ALL. Three of them were unreachable: they
    # tested `subductFactor`, which is a property of a cell's own crust and is
    # high on the DOWNGOING slab, while the rules used it as a proximity measure
    # for the overriding plate. Melange had the identical bug. They are now keyed
    # on backArcDist with a dip-modulated arc-trench gap, as a band about the
    # volcanic front rather than a disc from the trench.
    #
    # Melange goes 0.00% -> 6.17% of land and arc andesite appears for the first
    # time, so this build exposes rock classes no previous one ever did. Their
    # albedo and erodibility had therefore never been checked; both were audited
    # here, albedo held and erodibility did not (arc_andesite 0.90 -> 0.50,
    # rift_bimodal 0.95 -> 0.50, arc_basalt 0.85 -> 0.65, and the acid-vs-basic
    # ordering had been inverted against Moosdorf).
    #
    # PRE-CARVE. No carve verdict has been determined for this terrain.
    # `rift_bimodal` is 0% here and that is geography, not a bug: continental
    # rifts need two adjacent continental SUPER plates, each continent is its own
    # super plate, and none of the ten are adjacent. The same fact means there is
    # no continent-continent collision anywhere on this world.
    # Regenerated from the recipe in source/README.md after the previous export
    # was deleted unconsumed. Same hash to the bit, which is the check that the
    # recipe reproduces: the `--code` path and the explicit-slider path give
    # identical terrain.
    #
    # NOTE the hash does not cover the SCHEMA. This export renames
    # `surface_rock` to `substrate_class`, and nothing in any hash moves,
    # exactly as nothing moves when gravity changes. Two exports can share a
    # terrain hash and disagree on field names, so read the manifest field list
    # rather than assuming one from the hash.
    # Craton weight no longer multiplied by (1 - basin). That factor said a
    # craton under sedimentary cover is not a craton, which is backwards --
    # Earth's cratons are mostly platform -- and cover is modelled separately,
    # so it double-counted. Cratonic BASEMENT goes 1.41% of land to 3.46%.
    # Exposed gneiss barely moves, because these cratons carry thick basin
    # cover; that is the shield-versus-platform distinction, not a failure.
    #
    # The terrain hash moves even though this is a lithology change, because
    # craton weight also feeds CRATON_AMP_SUPPRESS and therefore relief.
    "e931b0d9947a232210a3cc11ba82003b5da15765d86e87956bc44c220932760e":
        {"name": "precarve-craton", "note":
         "2026-08 pre-carve base at g = 12.81; craton basin double-count "
         "removed, arc/forearc rules reachable, substrate_class naming"},
    "2e06d17682075096ab2e0093eb0f5618417cdee38383f1ca153646995bffad5a":
        {"name": "precarve-substrate", "note":
         "2026-08 pre-carve base at g = 12.81; arc/forearc rules reachable, arc "
         "erodibility and albedo grounded, surface_rock renamed substrate_class"},
    # NOT ACTIVATABLE AS IT STANDS, and registered anyway so the resolution audit
    # that measured it resolves. Same planet code, radius, gravity and lithology
    # strength as precarve-craton, differing ONLY in --regions: 10,000,005 against
    # 2,500,001, a mean edge of 7.59 km against 15.19. It carries all five
    # Gaussian grids on the T21/T42/T85/T127/T170 ladder, every one verified to
    # the same terrain hash, so the reason not to point source_build at it is no
    # longer the grid set. It was generated at glacialErosion 0.8, whose ice
    # mask PHYS-13 finds is placed by an
    # Earth-calibrated latitude threshold blind to this planet's obliquity,
    # spectrum and rotation, so a generation meant to be commissioned should
    # settle that first. `notes/audits/orogen-resolution.md` is what it was for.
    "ab0d679bd81360cd30fb67a3ce13e726b4aa1bb9ba7e4afb601c10aa13f9a526":
        {"name": "precarve-craton-10m", "note":
         "2026-08 resolution reference at 4x the region count: measurement "
         "artifact for the resolution audit, all five grids, ice mask not yet fixed"},
}

# Basin ids are computed on the pre-conditioning surface, so they survive a
# carve iteration. This hash is the check: if it matches, a carve verdict
# computed against an earlier export still refers to the same basins.
_KNOWN_CATALOGUE_HASHES = {
    "2d1f8e57c26b60c608deaa62bd3c44b7da04a4ee5bd518249447630095d96a98":
        "2026-08 catalogue, 3629 preserved from 81904 detected",
    # Same 3,629 basins, same ids, same order, same per-region membership --
    # `basinsPreserved` and `basinMembership` are byte-identical to the previous
    # catalogue and the id sets were diffed upstream. The hash moved only because
    # the catalogue's ...Km fields now carry the relief scaling they always
    # should have. So carve verdicts computed against 2d1f8e57 still apply.
    "bc84109168789f519b53fc8197ddc976be93ecd8326b61a8086de85e638343cc":
        "2026-08 catalogue, same 3629 basins, ...Km fields now relief-scaled",
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

        # Gravity is NOT in the terrain hash, and this is not an oversight in the
        # hash: Orogen runs its whole pipeline in model units and applies the 1/g
        # relief scaling only at the model-unit-to-km conversion on export. So two
        # builds at different gravities are bit-identical in `finalElevation`,
        # `basinCatalogue` and `params` alike, and differ only in
        # `manifest.planet.gravityMS2` and the `elevation_km` it scales.
        #
        # For us that means the terrain hash is not sufficient identity. A build
        # generated at another gravity would pass the allowlist silently while
        # every vertical quantity derived from it was off by the ratio -- 31% at
        # the gravities this project has used. Gravity is therefore checked
        # separately, against the config that declares it.
        self.gravity_m_s2 = float(self.manifest["planet"]["gravityMS2"])
        self.relief_scale = self.manifest["planet"].get("reliefScale")

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

    # Fields the fork has renamed. A bare AttributeError on an old name sends the
    # reader looking for a typo; naming the replacement costs one dict.
    RENAMED = {
        "surface_rock": "substrate_class",
        "surface_rock_pre_erosion": "substrate_class_pre_erosion",
    }

    def __getattr__(self, name: str) -> np.ndarray:
        # Convenience: export.elevation_km rather than export.field("elevation_km").
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self.RENAMED:
            raise AttributeError(
                f"{name!r} was renamed to {self.RENAMED[name]!r}. It is the top "
                f"of the cover/basement stack -- consolidated lithology plus "
                f"basin fill -- and was never what is at the surface, which is "
                f"derived downstream.")
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

    @cached_property
    def local_slope_deg(self) -> np.ndarray:
        """Dip of the least-squares plane through each region and its
        neighbours, in degrees.

        This is a REGIONAL dip, not a hillslope, and the distinction is the
        whole of its correct use. Neighbour spacing on this mesh is about 20 km,
        so what it measures is the tilt of the landscape across adjacent cells;
        a talus at its 34 degree repose angle is two orders of magnitude below
        anything the mesh can carry, and `vendor/orogen/tools/README.md` says so
        under scarps. Read it for "is this a low-relief surface or is it
        mountainous", never for "how steep is this slope".

        A PLANE FIT rather than the steepest finite difference to a neighbour.
        On a smooth analytic ramp the two barely differ -- the max-drop estimator
        is 1.2% low at the median, 11% at its 5th percentile, because the
        gradient rarely points exactly at a neighbour -- so this is not chosen
        for accuracy on smooth ground. It is chosen because on ROUGH ground the
        two measure different things: a max drop reports local roughness and a
        fit reports the tilt the roughness sits on, and it is the tilt that
        Freyssinet et al. (2005) p. 685 write when they give the "overall dips of
        1 to 5 degrees" of a planation surface. On this world's real terrain the
        two disagree by a factor of about two in the upper percentiles, so the
        choice is not free and it is made by what the source means.

        Built from `elevation_km`, so the planet's 1/g relief scaling is already
        inside it. Displacements are tangent-plane projections of the export's
        own unit-sphere coordinates times the manifest radius.
        """
        off, lst = self.adjacency
        degree = np.diff(off)
        if not (degree >= 3).all():
            raise RuntimeError(
                "a region has fewer than three neighbours, so the tangent-plane "
                "fit below is underdetermined for it")

        elev = self.field("elevation_km").astype(np.float64)
        pos = np.stack([self.field("x"), self.field("y"), self.field("z")],
                       axis=1).astype(np.float64)
        radius = float(self.manifest["planet"]["radiusKm"])

        src = np.repeat(np.arange(self.n_regions, dtype=np.int64), degree)
        dst = lst.astype(np.int64)

        # An ARBITRARY orthonormal tangent basis per region, not east/north.
        # Only the magnitude of the fitted gradient is returned and that is
        # invariant to the choice, so the basis is picked for robustness: cross
        # the normal with whichever axis it is least aligned with, which cannot
        # degenerate. A true east/north basis does, on the two regions this mesh
        # puts on the rotation axis.
        normal = pos
        ref = np.zeros_like(normal)
        ref[np.arange(len(normal)), np.argmin(np.abs(normal), axis=1)] = 1.0
        east = np.cross(normal, ref)
        east /= np.linalg.norm(east, axis=1)[:, None]
        north = np.cross(normal, east)

        # Chord to the neighbour, projected into the tangent plane at src.
        # Accumulated one cartesian component at a time: the stacked form needs
        # a 15M x 3 float64 temporary and roughly doubles peak memory for no
        # gain, and this reader is opened inside scripts that are already
        # holding a climatology.
        u = np.zeros(src.size)
        v = np.zeros(src.size)
        for k in range(3):
            step = (pos[dst, k] - pos[src, k]) * radius
            u += step * east[src, k]
            v += step * north[src, k]
            del step
        rise = elev[dst] - elev[src]

        seg = off[:-1].astype(np.int64)
        suu = np.add.reduceat(u * u, seg)
        suv = np.add.reduceat(u * v, seg)
        svv = np.add.reduceat(v * v, seg)
        sur = np.add.reduceat(u * rise, seg)
        svr = np.add.reduceat(v * rise, seg)
        del u, v, rise

        det = suu * svv - suv * suv
        scale = suu * svv
        bad = det <= 1e-12 * np.maximum(scale, 1e-30)
        if bad.any():
            raise RuntimeError(
                f"{int(bad.sum())} regions have collinear neighbours, so the "
                f"plane fit is singular there")
        gu = (svv * sur - suv * svr) / det
        gv = (suu * svr - suv * sur) / det
        return np.degrees(np.arctan(np.hypot(gu, gv))).astype(np.float32)

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
            # The NAME, which is a directory-safe build id, not the note.
            # These were prose strings and `build_hydrography.py` files its
            # products under this value, so a build whose note began with
            # something other than its own name produced a directory called
            # "2026-08 precarve-zoned-g1281-erod4x: pre-carve base at g = ...".
            # Registry entries are {name, note} now so the two cannot be
            # confused.
            "terrain_build": (_KNOWN_TERRAIN_HASHES.get(self.terrain_hash) or {})
                             .get("name", "unrecognised"),
            "terrain_build_note": (_KNOWN_TERRAIN_HASHES.get(self.terrain_hash) or {})
                                  .get("note"),
            "catalogue_hash": self.catalogue_hash,
            "catalogue_known": self.catalogue_hash in _KNOWN_CATALOGUE_HASHES,
            "drainage_hypothesis": self.manifest["basins"].get("drainageHypothesis"),
            "selection_source": self.manifest["basins"].get("selectionSource"),
            "seed": self.manifest["seed"],
            "num_regions": self.n_regions,
            "planet_radius_km": self.radius_km,
            "surface_area_km2": self.surface_area_km2,
        }
