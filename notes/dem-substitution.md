# Substituting a standard DEM for the Orogen export

*Written 2026-08-21. Numbers measured on `precarve-craton` (terrain
`e931b0d9`) and, where the fine support matters, on `precarve-craton-10m`
(`ab0d679b`).*

Vesper's geography comes from a fork of World Orogen, and every other component
is a consumer of that export. This note asks what would have to arrive alongside
a standard digital elevation model -- an Earth product, or any other externally
sourced raster -- for the rest of the pipeline to keep running. It doubles as an
inventory of what the tree actually reads from `source/`, which is the part
worth keeping whether or not anyone ever makes the substitution.

The short answer is that the elevation field is the smaller half of the export.
A DEM answers `elevation_km` and nothing else. Three bodies of information have
no elevation representation at all, a fourth is computable from a DEM but is
nowhere computed by the consumers, and one loop in the pipeline loses its sink.

## What has no elevation representation

### A land/sea class, not a sea-level test

`surface_class` is authoritative (CLAUDE.md rule 1) and its entire content is
the distinction elevation cannot carry: dry closed-basin floor below the datum
is land. The sign test is `land_mask`, and on `precarve-craton` the two disagree
over 1.91% of the planet, 48,119 regions, of which 4,560 sit outside every
preserved basin and so are invisible to the `land_mask | is_endorheic`
reconstruction as well. On the 10M fine support the same disagreement is 3.31%
and 332,768 regions, 16,406 of them outside the catalogue, because a finer mesh
resolves more small enclosed depressions.

So a DEM needs an ocean-connectivity mask beside it -- ocean means connected to
the world ocean, not below the datum -- plus an explicit decision about which
sub-datum depressions are dry. On Earth that decision is a known coastline
product and the cells at issue are the Qattara and Dead Sea floors. Without it,
`exoplasim/scripts/build_boundary_conditions.py`, all of hydrography, both
pedology surface products and both prospectivity scripts have no input.

### Lithology, as a map and as a class table

`substrate_class`, `basement_rock`, `cover_rock` and `cover_thickness`, plus
`manifest.lithology.rockClasses`, which carries per-class code, albedo and
erodibility (20 classes in this build). This is the largest dependency after
elevation. It feeds soil texture and phosphorus in `pedology/scripts/build_soil.py`,
the surface albedo fields, the surface term of aerodynamic roughness, dust
emission, the LPJ-GUESS driver, brine paths, solute routing, groundwater
hydraulic properties, and both minerals products. `erodibility`, `rock_albedo`
and `scarp_potential` are derived from the class table rather than measured, so
the table is the thing that has to arrive.

The Earth-side equivalent is a lithological map such as GLiM, resampled to the
same support, with per-class albedo and erodibility values checked the way this
project checks them rather than inherited from whatever the source shipped.

### Tectonic setting, including exhumation

`plate_index`, `plate_is_ocean`, `r_boundaryType`, `r_stress`,
`r_subductFactor`, `r_t_craton`, `r_t_foldBelt`, `backArcDist` and
`erosionDelta`. `minerals/scripts/build_prospectivity.py` is built directly on
craton weight, fold-belt weight, propagated stress and cumulative exhumation.

`erosionDelta` is the one with no substitute. It is how much rock has been
removed above the present surface, and a present-day elevation field does not
record it. Orogen has no time axis (`docs/src/reference/no-time-axis.md`) but it
does carry the cumulative amount. A DEM substitution either imports an
exhumation estimate from elsewhere or prospectivity loses its depth axis and
becomes a surface-lithology field only.

## What is computable but is nowhere computed

These are not missing information. They are missing artifacts: derivable from a
raster by standard methods, but no consumer derives them, because every consumer
expects them present.

**Support geometry.** `cell_area` summing to 4 pi R^2, the `x`/`y`/`z` unit
vectors in the generator's y-up frame, `lat`/`lon`, and mesh adjacency in CSR
under `manifest.raw.mesh`. `hydrography/scripts/drainage.py` and
`build_hydrography.py` walk that CSR directly; `lib/gridding.py:land_weighted`
needs the area and a region-to-cell map.

**The resolution ladder.** Five Gaussian grids from T21 to T170 plus the uniform
`grid-512x256` for cartography, each with `grid_cell_area` and
`gauss_weights.bin`, resampled by the export's own rules: mean for continuous
fields, nearest for categorical. Underneath them the pipeline wants a support
finer than the model grid, because `build_surface_roughness.py` takes its
orographic roughness from the distribution of elevation inside a T42 cell.
That is about 305 regions per cell globally on `precarve-craton` and about 1,221
on the 10M reference. A DEM at or coarser than the model grid cannot supply it;
a 30 m DEM supplies it with room to spare.

**Drainage and the basin catalogue.** `drain_to`, `drainage_terminal`,
`flow_accumulation`, `ocean_basin_index`, `basin_index`, `is_endorheic`, and
`manifest.basins` with the `preserved[]` hypsometry curves, spill elevations,
volumes, `resolution.avgEdgeKm`, `drainageConsistency` and `counts` (3,621
preserved basins here, 9,419 on the fine support). Priority-flood depression
analysis is exactly what a DEM toolchain does, and hydrography already refloods
independently, using the export's `drainage_terminal` as a cross-check rather
than as input. But `lib/orogen.py:Basin` and `hydrography/scripts/lake_balance.py`
read the catalogue out of the manifest, so it has to exist there.

**Bathymetry, mostly not needed.** ExoPlaSim runs a slab ocean, and
`build_boundary_conditions.py` writes only a land mask and land geopotential,
zeroing elevation wherever the cell is not land. A land-only DEM therefore costs
the climate component nothing. Ocean depth enters only the ocean-basin and sill
bookkeeping in `manifest.hydrology`.

## The planet block, and the relief decision

`manifest.planet` supplies `radiusKm`, `gravityMS2`, `surfaceAreaKm2` and
`reliefScale`; `lib/orogen.py` reads gravity from it and checks it against
`config/planet.yaml`.

The sharp point is that `elevation_km` already carries the 1/g scaling for this
planet's gravity, while an Earth DEM's relief was built at Earth gravity on an
Earth radius. Substituting one for the other forces an explicit choice: rescale
the amplitude and reproject onto Vesper's radius, or declare the DEM
authoritative and change `config/planet.yaml` to match. Rule 2 does not allow
leaving the two in disagreement, and `derive()` raises when gravity and mass do
not agree. Whichever way it goes, the choice is a decision to record, not a
detail of the import script.

## Packaging and identity

`Export.__init__` refuses a build whose `hashes.finalElevation` is not in
`lib/orogen.py:_KNOWN_TERRAIN_HASHES`, and the basin catalogue has its own
allowlist beside it. A DEM therefore has to be packaged as a build before
anything will read it: a `manifest.json` carrying the field catalogue and the
`hashes`, `planet`, `lithology`, `basins`, `hydrology` and `landSeaMask` blocks,
a registry entry saying what the build is, a row in `source/README.md`, and
resolution through `lib/builds.py`. `scripts/check_consistency.py` then enforces
one terrain hash across every artifact exactly as it does now.

Rule 7 changes shape under the substitution. The durable set currently stops at
the planet code, the seed and the recipe, because the export is regenerable from
them. An imported DEM is not regenerable, so the raster itself joins the durable
set, along with the projection and rescaling that produced the export from it,
and the archive rules have to say so.

## What stops existing

Loop A has no sink. `carve_verdict` to `carve_list` to `orogen_generate` exists
because Orogen can regenerate terrain with basins cut; a fixed raster cannot be
regenerated. Either the carve becomes a conditioning pass that edits the DEM in
place, which is ordinary hydrological conditioning and keeps the loop in a
different form, or the world is a single iteration and those steps plus the
endmember bracket leave `config/pipeline.yaml`.

Gone with no substitute: `elevation_pre_erosion`, which
`analysis/orogen_resolution.py` reads to measure the terrain-information floor,
and `scarp_potential`, which `analysis/orogen_resolution_controls.py` reads as
its control. `elevation_pre_conditioning` goes too, and so does the
`climateComputed` block; nothing currently reads either.

## Verdict

The substitution is bounded. Beside the raster it needs an ocean-connectivity
mask that keeps dry sub-datum basins as land, a lithology map with a class table
carrying albedo and erodibility, a tectonic map with craton, fold-belt, boundary
type and exhumation, the planet parameters with an explicit relief-rescaling
decision, and a packaging layer that emits all of it on a fine native support
and the five Gaussian grids with a manifest and a registered terrain hash.
Drainage and the basin catalogue are computable, but they still have to be
written into the manifest. What does not survive is the carve loop and the
resolution diagnostics that read Orogen's intermediate surfaces.
