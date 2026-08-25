# The endorheic basin catalogue's floor: where it is, and what sets it

Vesper is a generated world. This note is about World Orogen's selection of
which modelled depressions enter the endorheic basin catalogue of an export of
its terrain, and every number here is a property of that generated terrain.

The catalogue matters out of proportion to its size. The carve list is drawn
from the preserved set, loop A closes by Orogen consuming that list, and the
generator regenerates rather than edits, so whatever decides which depressions
are catalogued decides the terrain of every later pass.

Measured on 2026-08-24 against the two registered builds that carry a payload,
`precarve-craton` at 2,500,001 regions and `precarve-craton-10m` at 10,000,005.
BOTH ARE PRE-CARVE, so every basin number below is what is physically sitting in
`source/` before any carve verdict has been applied: a limit, not a state. The
measurement is `hydrography/scripts/catalogue_floor.py`, registered under
`one_offs` in `config/pipeline.yaml`, which writes
`hydrography/analysis/catalogue_floor.json`.

## The instrument, and why it is one

`selectBasins` applies three floors and then a nesting rule that drops any
depression lying inside an already-selected one. Counting how many depressions
fail each floor overstates every floor, because the nesting rule then removes
much of what a relaxed floor would have admitted. So the measurement relaxes one
floor and re-runs the WHOLE selection, and reports the preserved set that comes
out.

That is only worth reading if the re-run is the selection Orogen performed.
Three controls, each able to fail, and the script exits nonzero on the first two:

- REPRODUCTION. Reapplying the published criteria and the nesting rule to the
  published catalogue returns the published preserved set, id for id, on both
  builds: 3,621 and 9,419, with no id on either side alone.
- SENSITIVITY. Halving every floor must admit strictly more and doubling every
  floor strictly fewer. Reproduction alone is not evidence, because a routine
  that ignored the criteria and echoed the published list would pass it. A
  single floor is deliberately not asserted on: a floor another floor dominates
  is inert on that build in the relaxing direction, and that inertness is the
  measurement rather than a fault. On `precarve-craton` the area floor is inert
  in exactly that direction, which is the finding below.
- RESIDUE. The below-sea-level counts derived here from `surface_class`,
  `land_mask` and `is_endorheic` equal the generator's own `landSeaMask` counts
  in the manifest, 48,119 / 43,559 and 332,768 / 316,362. Two independent paths
  to one number.

## The mesh floor stops binding between the two builds

Relaxing `minCells` and re-selecting:

| build | regions | mean cell km2 | preserved | with minCells relaxed | added |
| --- | --- | --- | --- | --- | --- |
| precarve-craton | 2,500,001 | 293.80 | 3,621 | 8,084 | +123.3% |
| precarve-craton-10m | 10,000,005 | 73.45 | 9,419 | 9,427 | +0.1% |

At the smaller build the mesh is keeping out more basins than the catalogue
contains. At the larger one it keeps out eight. That is the crossover
`notes/audits/orogen-resolution.md` predicted near 8.8M regions, confirmed here
against the preserved set rather than against a cell-area arithmetic, and it
settles the resolution half of the question: at 10M regions the catalogue is not
resolution-limited in area.

Two things the raw rejection counts get wrong, and both are the same mistake in
different clothes. Before nesting, 18 depressions at 10M clear the depth and
area floors and fail on cells; after nesting only 8 of them would ever have been
preserved, because the other 10 sit inside a basin already selected. And the
manifest's `resolution.bindingFloor` is computed from the MEAN cell area over
the sphere, so it is a statement about the typical depression rather than a
guarantee about any particular one: local cell area varies with the mesh jitter,
and those 18 are depressions above `effectiveMinAreaKm2` that the cell floor
rejected anyway, in the region of the mesh where cells run large. The field now
says which it is.

## The declared area floor is what binds at 10M, and it is physical

Relaxing `minAreaKm2` at `precarve-craton-10m` adds 1,105 basins and about
9.9e5 km2; at `precarve-craton` it adds nothing at all, because the cell floor
there already excludes everything the area floor would. So the two builds are on
opposite sides of the crossover, and the larger one is governed by a fixed
physical size in real km2. That is the sense in which the catalogue can be made
physical, and at 10M regions in area it already is.

## The depth floor is neither physical nor resolution-set

`BASIN_MIN_DEPTH_KM` is compared against a depression's depth in the model's
DIMENSIONLESS elevation parameter, not in kilometres. The comparison is
deliberate: that is the currency the terrain noise the floor exists to reject
lives in. What is not deliberate is that the constant, the manifest key
`resolution.minDepthKm` and the published `selectionCriteria.minDepthKm` all
name a physical depth, so a consumer reasons about a 50 m floor that was never
enforced.

The model-unit-to-km curve is quartic above sea level, linear below it, and
SATURATES at model elevation 1. One threshold is therefore many physical depths.
Measured over the depressions sitting within 2% of the threshold:

| build | physical depth the threshold demanded, km |
| --- | --- |
| precarve-craton | 0.000 to 0.485, median 0.062 |
| precarve-craton-10m | 0.000 to 0.492, median 0.085 |

so the same declared number admits depressions a few metres deep near sea level
and demands hundreds of metres on high ground. The consequence in the preserved
set, which is what the carve list is drawn from:

| build | preserved | shallower than the declared floor | share |
| --- | --- | --- | --- |
| precarve-craton | 3,621 | 773 | 21.3% |
| precarve-craton-10m | 9,419 | 2,077 | 22.1% |

The saturation is the sharp end of it. A depression whose sink and spill both
sit above model elevation 1 has a published physical depth of exactly zero and
still clears the floor. `hydrography/data/*/basins.nc` carries the downstream
consequence: 5 preserved basins on `precarve-craton` and 16 on
`precarve-craton-10m` have `capacity_km3` of exactly zero, which is a preserved
closed basin that can impound no water at all. The clamp is not confined to
basins: 0.26% and 0.25% of land by area sits at or above model elevation 1 on
the two builds, and every cell of it is published at the same `elevation_km`.
That last figure belongs to whoever owns the relief curve, not to the catalogue.

Applying the same declared number to `depthKm` instead is a counterfactual, not
a proposal, because it changes the preserved set and therefore the terrain:

| build | as built | with the floor in km | drops | admits |
| --- | --- | --- | --- | --- |
| precarve-craton | 3,621 | 2,879 | 773 | 31 |
| precarve-craton-10m | 9,419 | 7,471 | 2,077 | 129 |

Both directions are populated, so this is not a floor that is merely too loose.
It is loose near sea level and tight on high ground, and the sign of the error
is set by where the depression sits rather than by how deep it is.

Relaxing the depth floor entirely at 10M adds 2,521 basins against the area
floor's 1,105 and the cell floor's 8. So at the region count this project has
settled on, the floor that decides the catalogue is the one whose units are
wrong, and no further refinement of the mesh moves it.

## What the catalogue misses, in rule 1's terms

CLAUDE.md rule 1 says land comes from `surface_class`, never from `land_mask`,
and that reconstructing it as `land_mask | is_endorheic` misses the depressions
too small to enter the basin catalogue. This is how much, and the answer is an
AREA, because a share of land is an area share. Cell areas are extensive and are
summed; the count share is a different statistic and is given beside it.

| build | below sea level and dry | outside any preserved basin | of land | by count |
| --- | --- | --- | --- | --- |
| precarve-craton | 1.404e7 km2 | 1.286e6 km2 | 0.405% | 9.48% |
| precarve-craton-10m | 2.436e7 km2 | 1.157e6 km2 | 0.364% | 4.93% |

The missed AREA is converged and the missed COUNT is not. Four times the region
count multiplies the number of missed pieces by 3.6 while leaving their total
area within 10% of where it was, because every newly resolved depression falls
under the declared 1000 km2 area floor and none is kept. So the floor rule 1
warns about sits at roughly 0.4% of land area, it is a property of the declared
area criterion rather than of the mesh, and refining the generation does not
close it. Lowering `minAreaKm2` would; nothing else will.

The deepest ground outside every preserved basin is -0.513 km on
`precarve-craton-10m`, against -0.591 km inside one, so what the catalogue
misses is not shallow ground. It is ordinary dry closed-basin floor in pieces
below the declared area floor.

## Verdict

- Resolution-limited: TRUE on `precarve-craton`, where the mesh keeps out more
  basins than the catalogue holds. FALSE on `precarve-craton-10m`, where it
  keeps out eight.
- Physics-limited: TRUE at 10M in area, where a fixed km2 floor governs.
- Neither, in depth: the depth floor is compared in model units on both builds
  and governs the catalogue at 10M. Making it physical is a change to the
  selection criteria, which changes the terrain, so it is a loop A decision and
  belongs with the generation that carries it rather than with this note.
