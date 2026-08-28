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
#
# An entry is {name, note} plus an OPTIONAL `refusal`. REGISTERED and ACTIVATABLE
# are two different questions and this key is what separates them. Registration
# says an export with this hash is recognised, which is what keeps a result
# already computed from it readable; `refusal` says the tree may not POINT at it,
# and carries the reason as a sentence. Its absence is the only statement that a
# build may be named as `config/planet.yaml`'s `source_build`, and
# `scripts/check_consistency.py` is what enforces that.
#
# The key exists because the refusals used to be prose: "NOT ACTIVATABLE AS IT
# STANDS" in a comment above one entry, "Do not use." at the end of another's
# note, "SUPERSEDED by <hash>" inside a third's. No gate can read any of that, so
# config named a build the registry refused for as long as nobody noticed, and
# batch P0C measured the basin catalogue on two builds because it could not tell
# which was active. Add a `refusal` whenever a note gains the word WITHDRAWN,
# SUPERSEDED or NOT ACTIVATABLE; those three are the whole trigger, and a build
# that acquires one and does not get the key is invisible to the guard again.
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
        "821aa71b. Do not use.",
         "refusal":
         "withdrawn: the basin-protection floor was captured pre-erosion and "
         "stored absolute, so the carve took back divides erosion had raised. "
         "Reverted upstream, and the regenerated export is bit-identical to "
         "821aa71b, which is the build to name instead."},
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
        "split. SUPERSEDED by 010f2143 -- its verdict used antipodal climate.",
         "refusal":
         "superseded by 010f2143: 850 of its 1,522 carves were decided on the "
         "antipode's rainfall and cannot be un-cut."},
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
         "cover-chain branch order.",
         "refusal":
         "superseded by 5bed5549: closed-basin fill lost to cover-chain branch "
         "order, and 203 preserved basins reach ExoPlaSim as vegetated land."},
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
         "SUPERSEDED by a4d204f6 -- its verdict read mrro as catchment runoff.",
         "refusal":
         "superseded by a4d204f6: its verdict read mrro as catchment runoff, "
         "which is net divergence rather than local generation, so 749 basins "
         "that overflow under P - E are still standing on this terrain."},
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
    # Refused for activation -- the `refusal` below is the statement of that, and
    # the rest of this comment is its argument -- and registered anyway so the
    # resolution audit
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
         "artifact for the resolution audit, all five grids, ice mask not yet fixed",
         "refusal":
         "generated at glacialErosion 0.8, whose ice mask PHYS-13 finds is placed "
         "by an Earth-calibrated latitude threshold blind to this planet's "
         "obliquity, spectrum and rotation. A generation meant to be commissioned "
         "settles that first; regenerating gives a new terrain hash and an entry "
         "with no refusal. Registered so notes/audits/orogen-resolution.md, which "
         "measured this terrain, still resolves."},
    # NO `refusal`, and that is the whole point of this entry. It is the first
    # generation taken with a view to being commissioned rather than measured:
    # `--glacial 0`, so the ice is an honest null awaiting the mask PHYS-13's
    # route supplies, rather than the Earth-calibrated latitude ramp that
    # refuses precarve-craton-10m above. Same planet code, seed, radius, gravity
    # and lithology strength as the two builds above; against
    # precarve-craton-10m it differs in glacialErosion, in the basin floor, and
    # in the generator changes since -- the depth floor compared in physical
    # kilometres (WORLD-YRIL) and the relief curve's branch above the shape
    # function's domain.
    #
    # `--basin-min-area 850` is passed EXPLICITLY and has to be. The planet
    # code's basin slider maps onto a ladder with no 850 rung and it outranks
    # BASIN_MIN_AREA_KM2 unless a basin flag is given, so two generations were
    # taken at 1000 before the exporter was made to refuse the disagreement.
    # `source/README.md` carries the recipe and the argument.
    #
    # Its ANCESTRY is not settled terrain: LITH-26 replaces the shelf substrate
    # classification and moves terrain when it lands, and it lands on a later
    # iteration because it needs a climatology this build has to produce first.
    # `notes/audits/orogen-first-pass-gate.md` is the reading that says nothing
    # outstanding changes the FIRST export, which is the question the orogen
    # step's gate asks.
    "2004672995403ead35e5ccbee400ae10cf3549ca662363c629c798707766c0c9":
        {"name": "canonical-10m-base", "note":
         "2026-08 first pass at 10M without ice: --glacial 0, basin floor 850, "
         "all five ladder grids plus grid-512x256, pre-carve"},
    # The SECOND pass, and the first this project has taken with both of the
    # generation's inputs supplied. `canonical-10m-base` was commissioned to a
    # baseline; that climatology produced the carve list and the ice mask; both
    # are consumed AT GENERATION, which is why this is a new build and not an
    # alteration of the one above.
    #
    # WHAT MOVED AND WHAT DID NOT, from the two manifests' hash blocks:
    # `preConditioningElevation` and `preErosionElevation` are BYTE-IDENTICAL to
    # the build above, which is what the seed and region count being unchanged
    # requires. `finalElevation`, `basinsPreserved`, `basinMembership`,
    # `surfaceRock` and `erodibility` all moved, because the carve and the ice
    # are erosion-loop inputs and the loop exposes different rock.
    #
    # `--glacial 0.3` and NOT the 0.8 the pre-carve builds carried. With the
    # mask supplying placement, the slider does only the erosion-rate half of
    # the job it used to do, and 0.8 was chosen for the half it no longer does.
    # It is undeclared either way and world-sr77 carries that, deferred.
    # The THIRD pass, and the first whose carve list came from a real
    # intersection: 195 basins both bounding climates agreed to cut, out of a
    # union of 485 where they disagreed about 290. carve1's list was taken the
    # same way; what is new is that the disagreement is now most of the union,
    # which is what a second pass on the marginal population looks like.
    #
    # Its ice mask is written from carve1's own baseline climatology rather
    # than carve1's donor: 307,118 regions against 313,790, because the carved
    # world is warmer. Same seed, same region count, same 850 km2 basin floor,
    # same --glacial 0.3, and the same basinCatalogue hash as both builds
    # before it, so a verdict computed against any of the three still refers
    # to the same basins.
    "f496ae9fd749da592cd2f0d588f44819fcee488922e249b09b91f9310ddfc3ec":
        {"name": "canonical-10m-carve2", "note":
         "2026-08 third pass at 10M: 195 basins cut and 40 notched out of the "
         "4657 carve1 left standing, ice from carve1's baseline freezing "
         "height, --glacial 0.3, basin floor 850"},
    "4884dc8a6120bf5661f0e6fefbbf8d6bc30c78d7d2f1de6ce5a6c93073ff5afe":
        {"name": "canonical-10m-carve1", "note":
         "2026-08 second pass at 10M: the first carve, 4115 basins cut and 60 "
         "notched at their saddles out of a catalogue of 8772, with ice placed "
         "from the baseline climatology's freezing height rather than the "
         "latitude ramp -- --ice-mask, --glacial 0.3, basin floor 850"},
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
    # canonical-10m-base. A DIFFERENT basin population, not the same one
    # re-measured: no carve verdict computed against a catalogue above carries
    # over, and none exists to carry. Selected at the 850 floor on ice-free
    # terrain, where the binding floor is minCells -- 12 cells over a mean
    # 73.449 km2 cell is 881.39 km2, above the declared 850, so the area floor
    # does not bind at this region count. `source/README.md` put that crossover
    # at 10.37M regions before the build was taken.
    # SHARED BY canonical-10m-base AND canonical-10m-carve1, byte-identical
    # across the carve, which is this key's whole purpose demonstrated rather
    # than asserted. The ids are computed on the pre-conditioning surface, so a
    # carve iteration cannot move them, and the two builds' manifests carry the
    # same `preConditioningElevation` hash as well. A verdict computed against
    # this catalogue therefore still refers to these basins on either build.
    #
    # `preserved` is NOT part of this hash and differs between them: 8772 on
    # the base, where no preserve list was given, and 4657 on carve1, where the
    # list carved 4115 and notched 60 more. The catalogue is what was selected;
    # the preserved set is what survived a water balance.
    "35c922c3b481af4bf7b814a6719e05969db0223b9a01895beef01cf200a73b28":
        "2026-08 catalogue at the 850 floor with minCells binding, 8772 "
        "selected, shared by canonical-10m-base, canonical-10m-carve1 and "
        "canonical-10m-carve2 -- three builds and two carves apart, and "
        "byte-identical across all of them",
}

# The registry is keyed by hash because the hash is the identity, but config and
# `source/` both address a build by NAME, so the guard has to be able to go that
# way too. Two entries sharing a name would make that lookup return whichever
# came first, and it would also put two terrains in one `<component>/data/<name>`
# directory, so the collision is worth more than a wrong answer here.
_BUILDS_BY_NAME: dict[str, dict] = {}
for _hash, _entry in _KNOWN_TERRAIN_HASHES.items():
    if _entry["name"] in _BUILDS_BY_NAME:
        raise RuntimeError(
            f"two registry entries are both named {_entry['name']!r}; a build "
            "name addresses a source/ directory and a per-build data directory, "
            "so it has to be unique")
    _BUILDS_BY_NAME[_entry["name"]] = _entry
del _hash, _entry


def registry_entry(terrain_hash: str | None = None, *,
                   name: str | None = None) -> dict | None:
    """The registry's record for a build, by hash or by name; None if unknown.

    Prefer the hash: a name is a label someone chose and a hash is what the
    export is. The name route exists for the case where there is nothing on disk
    to hash, which is a legitimate state -- a build is disposable until a climate
    run has consumed it, so between a generator change and the next generation
    `source/` is deliberately empty while `config/planet.yaml` still names one.
    """
    if terrain_hash is not None:
        entry = _KNOWN_TERRAIN_HASHES.get(terrain_hash)
        if entry is not None:
            return entry
    if name is not None:
        return _BUILDS_BY_NAME.get(name)
    return None


def activation_refusal(terrain_hash: str | None = None, *,
                       name: str | None = None) -> str | None:
    """Why this build may not be named as `source_build`, or None if it may be.

    This is about ACTIVATION only. Reading a refused build's export stays legal
    and has to: `notes/audits/orogen-resolution.md` measured one, and a result is
    only readable while the build it came from is still recognised here.
    """
    entry = registry_entry(terrain_hash, name=name)
    return (entry or {}).get("refusal")


OCEAN, LAND, INLAND_WATER = 0, 1, 2


def plane_fit_slope_deg(off, lst, pos, elev, radius) -> np.ndarray:
    """Dip of the least-squares plane through each generator and its neighbours.

    The mesh arithmetic of `Export.local_slope_deg`, lifted out of the class so
    that a mesh which is not an Orogen export can be given the SAME estimator
    rather than a second one that looks like it. The Earth calibration harness
    builds its own Fibonacci mesh and its own Voronoi adjacency, and a slope
    computed there by a different rule would make a comparison against this
    world's index a comparison of two estimators.

    `off`/`lst` are a symmetric CSR neighbour structure, `pos` the generators on
    the unit sphere as `(n, 3)`, `elev` their elevation and `radius` the sphere
    radius IN THE SAME LENGTH UNIT as `elev`; the returned dip does not depend
    on which unit that is, only on the two agreeing. Read `local_slope_deg`'s
    docstring for what the number means and what it may not be used for.
    """
    off = np.asarray(off, dtype=np.int64)
    lst = np.asarray(lst, dtype=np.int64)
    pos = np.asarray(pos, dtype=np.float64)
    elev = np.asarray(elev, dtype=np.float64)
    n = off.size - 1
    degree = np.diff(off)
    if not (degree >= 3).all():
        raise RuntimeError(
            "a region has fewer than three neighbours, so the tangent-plane "
            "fit below is underdetermined for it")

    # An ARBITRARY orthonormal tangent basis per region, not east/north. Only
    # the magnitude of the fitted gradient is returned and that is invariant to
    # the choice, so the basis is picked for robustness: cross the normal with
    # whichever axis it is least aligned with, which cannot degenerate.
    normal = pos
    ref = np.zeros_like(normal)
    ref[np.arange(n), np.argmin(np.abs(normal), axis=1)] = 1.0
    east = np.cross(normal, ref)
    east /= np.linalg.norm(east, axis=1)[:, None]
    north = np.cross(normal, east)

    src = np.repeat(np.arange(n, dtype=np.int64), degree)
    dst = lst
    u = np.zeros(src.size)
    v = np.zeros(src.size)
    for k in range(3):
        step = (pos[dst, k] - pos[src, k]) * float(radius)
        u += step * east[src, k]
        v += step * north[src, k]
        del step
    rise = elev[dst] - elev[src]

    seg = off[:-1]
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
    # Depression relief on the two surfaces, in the generator's dimensionless
    # elevation parameter. `natural_spill_depth` is what detectBasins measured
    # on the pre-conditioning surface and is the basis Orogen's `retain`
    # allowance is a fraction of; `final_spill_depth` is what survived. The two
    # differ on almost every basin, which is why both are carried.
    natural_spill_depth: float
    final_spill_depth: float
    final_spill_depth_km: float
    retained_fraction: float | None

    @property
    def natural_hypsometry_is_usable(self) -> bool:
        """The natural curve is only safe where erosion barely touched the rim."""
        return self.final_volume_km3 >= 0.95 * self.natural_volume_km3

    @property
    def still_closed(self) -> bool:
        """Is there still a depression on the finished terrain?

        PRESERVED IS NOT THE SAME AS STILL CLOSED, and this is the field that
        separates them. Preservation is a promise about ONE agent: the
        drainage-enforcement carve will not lower a protected divide below its
        floor, and `assertDividesNotLowered` proves it pass by pass. It is not
        a promise about the depression, because hydraulic, thermal and glacial
        erosion, ridge sharpening and soil creep all run over protected cells
        by design, and the sink is deliberately left unprotected so the basin
        floor can erode. On both registered builds the median preserved basin
        keeps about four fifths of its relief, a sixth of them gain relief, and
        a handful keep none at all.

        A basin that keeps none impounds nothing: it is published as preserved,
        with `retain: 1` and a hypsometry curve measured on a surface the
        finished terrain no longer has. Consumers that need an impoundment --
        a lake balance, a capacity, a carve verdict -- must test this rather
        than membership of the preserved set.

        notes/audits/basin-catalogue-floor.md carries the measurement.
        """
        return self.final_spill_depth > 0.0


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
        return plane_fit_slope_deg(
            off, lst,
            np.stack([self.field("x"), self.field("y"), self.field("z")], axis=1),
            self.field("elevation_km").astype(np.float64),
            float(self.manifest["planet"]["radiusKm"]))

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
                natural_spill_depth=float(b["natural"]["spillDepth"]),
                final_spill_depth=float(fp["spillDepth"] or 0.0),
                final_spill_depth_km=float(fp["spillDepthKm"] or 0.0),
                retained_fraction=(None if fp.get("retainedFraction") is None
                                   else float(fp["retainedFraction"])),
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
