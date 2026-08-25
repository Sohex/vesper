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

## The depth floor was neither physical nor resolution-set

Both builds in `source/` were selected with `BASIN_MIN_DEPTH_KM` compared
against a depression's depth in the model's DIMENSIONLESS elevation parameter,
not in kilometres, while the constant, the manifest key
`resolution.minDepthKm` and the published `selectionCriteria.minDepthKm` all
name a physical depth. A consumer of either build therefore reasons about a
50 m floor that was never enforced, and the numbers below are what that cost.

`selectBasins` now compares the floor against `b.depthKm`, and
`basinResolutionContext` publishes `minDepthComparedIn` so the currency is
read from the manifest rather than inferred. Neither build inherits that: the
comparison is a selection criterion, so it reaches the catalogue only through
a generation. The counterfactual table below is what the next one carries,
measured on this terrain.

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
still clears the floor: 2 preserved basins on `precarve-craton` are in that
state and none on `precarve-craton-10m`. The clamp is not confined to basins:
0.26% and 0.25% of land by area sits at or above model elevation 1 on the two
builds, and every cell of it is published at the same `elevation_km`. That last
figure belongs to whoever owns the relief curve, not to the catalogue.

A preserved basin that can impound no water is a SEPARATE fault with a separate
cause, and the depth floor is not it. `hydrography/data/*/basins.nc` gives
`capacity_km3` of exactly zero for 5 preserved basins on `precarve-craton` and
16 on `precarve-craton-10m`, and on the larger build not one of those 16 is
shallow: their published `depthKm` runs from 0.046 to 1.644 km and 15 of the 16
clear the declared floor in kilometres. What they share is
`finalPreserved.retainedFraction` of exactly zero -- the depression the
catalogue measured on the pre-conditioning surface does not exist on the final
terrain, because the drainage conditioning breached its rim. It is not confined
to those: the median preserved basin on both builds retains about 0.80 of its
relief, 7,702 of 9,419 retain less than all of it, and 18 retain none. Making
the depth floor physical does not remove them, and neither build carries a
carve list, so this is the conditioning and not a verdict.

MEASURED 2026-08-25 from `basins.preserved[].finalPreserved` in both manifests
and `capacity_km3` in `hydrography/data/<build>/basins.nc`.

Applying the same declared number to `depthKm` instead does not reach either
build, because it changes the preserved set and therefore the terrain. This is
what the next generation carries, measured on this one:

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
  and governs the catalogue at 10M. `selectBasins` now compares it in
  kilometres, which makes it physical, but a selection criterion reaches a
  catalogue only through a generation -- so on everything currently in
  `source/` the verdict above still stands, and the change lands with the
  generation `world-q5ig` carries.

# Preserved is not the same as still closed

Measured on 2026-08-25 against the same two builds, from
`basins.preserved[].finalPreserved` in each manifest and from `capacity_km3` in
`hydrography/data/<build>/basins.nc`. Both builds are pre-carve and neither
carries a preserve or carve list -- `drainageHypothesis` is empty on both -- so
nothing here is a declared verdict being honoured or broken.

The catalogue's preserved set is the input to the carve list, and the carve list
is what closes loop A. So what "preserved" entitles a basin to is not a naming
question.

## What the conditioning does to a preserved basin

| build | preserved | retained < 1 | retained > 1 | retained = 0 | median retained |
| --- | ---: | ---: | ---: | ---: | ---: |
| precarve-craton | 3,621 | 2,988 | 633 | 3 | 0.796 |
| precarve-craton-10m | 9,419 | 7,702 | 1,717 | 18 | 0.804 |

The median preserved basin comes out of a generation with about four fifths of
the relief the catalogue measured on it. About a sixth come out DEEPER. A
handful come out with none: 18 at 10M, of which 16 carry `capacity_km3` of
exactly zero in `basins.nc`, which is an independent instrument -- that file
recomputes hypsometry from its own priority flood over the finished terrain
rather than reading the catalogue's curve.

It is not the depth floor and not the height-curve saturation. Of the 18, 16
have a published `depthKm` clearing the declared 0.05 km floor in kilometres,
running up to 1.644 km, so `world-yril` making the floor physical leaves every
one of them in the catalogue.

## Protection is not leaking, and the rim comes down anyway

Preservation is a promise about ONE agent. `buildBasinProtection` builds a
`noLower` divide band and a per-cell `carveAllowance`, `priorityFloodCarve`
clamps its carve kernel against them, and `assertDividesNotLowered` proves pass
by pass that no protected divide went below its floor. Both call sites pass
protection, and the invariant is a throw rather than a warning. That contract is
kept.

It is a narrow contract, and `buildBasinProtection` says so: "Ordinary
hydraulic, thermal and glacial erosion still run over protected cells. Only the
drainage-enforcement carving is held off". Three further things follow from the
pipeline's own ordering, and together they are the whole of the effect:

- **The sink is deliberately unprotected.** `buildBasinProtection` clears
  `noLower` at every sink so the basin floor can erode, which is what lets a
  sixth of the population deepen.
- **Erosion runs between the carve passes.** Each `priorityFloodCarve` call
  re-snapshots its floor from the elevation it finds, so erosion lowering a rim
  between two passes is invisible to the invariant by construction.
- **Two passes after the carve are unprotected entirely.** `sharpenRidges` and
  `applySoilCreep` take no `protection` argument and run after
  `erodeComposite`, and `remeasureBasins` measures after both.

Attributing the loss confirms it is the rim rather than the floor that moves.
On `precarve-craton-10m` the final spill sits below the natural one on 8,993 of
9,419 preserved basins and above it on 426; the sink rises on 1,660 and falls on
3,755; and among the basins that lost relief the rim accounts for the whole of
the loss on the median basin. A breached rim and a raised floor are both
present, the rim dominates, and neither is the carve.

So the answer to "is a median retained fraction of 0.80 the intended behaviour"
is yes, in the sense that every mechanism producing it is one the pipeline
deliberately leaves enabled. `vendor/orogen/tools/test-basins.mjs` asserts the
carve contract and passes because the carve contract holds.

## What is not honest is the publication

The basins that keep nothing are published with `selectedBy: 'threshold'`,
`retain: 1`, and a `hypsometry` curve computed on the pre-conditioning surface.
Every one of those three describes a depression the finished terrain does not
have. Nothing on the artifact distinguishes a preserved basin that still
impounds from one that does not, and until 2026-08-25 nothing on the Python side
could: `lib/orogen.py`'s `Basin` carried no retention field at all, and
`basins.nc` carried none either, so a consumer asking "is this a closed basin"
got the preserved set and no way to narrow it.

That is now carried rather than argued: `Basin.still_closed`,
`basins.nc:has_impoundment`, `natural_spill_depth` and `retained_fraction`, and
`BasinSet.has_impoundment` in `lake_balance.py`.

## The consequence for a carve list, and it is not the zeros

The zeros are 16 basins in 9,419 and they carve, which is the right instruction
reached by the wrong route: `depth_at_spill_m` floors at 1 m, so any nonzero
incision drives their retain to 0. `export_carve_list.py` now counts them and
says so rather than letting the floor pass for a verdict.

The larger consequence is on the marginal class, and it is a basis mismatch
rather than a defect in either component.

**Hydrography computes retain against the FINISHED depression, in metres.**
`retain = clip(1 - cut/depth_at_spill, 0, 1)`, and `depth_at_spill_m` comes from
`basins.nc`, whose spill and sink are both measured on the conditioned terrain,
in kilometres through the height curve.

**Orogen spends it against the NATURAL depression, in model units.**
`buildBasinProtection` sets `carveAllowance = (1 - retain) * selected[i].depth`,
where `depth` is what `detectBasins` measured on the pre-conditioning surface,
in the dimensionless elevation parameter. It has to be: protection is built
before erosion, so the finished depth does not exist yet, and the carve operates
on model elevations.

Both statements are locally correct and neither artifact says which depth it
means. The two bases differ by exactly `retainedFraction`, and the height curve
is quartic on land, so the gap is neither small nor a constant:

| build | natural/final depth, model units | the same in km |
| --- | ---: | ---: |
| precarve-craton | 1.254 | 1.754 |
| precarve-craton-10m | 1.242 | 1.686 |

Expressed as the physical metres of rim Orogen is permitted to cut against the
metres hydrography intended, over every preserved basin on `precarve-craton-10m`:

| retain sent | actual / intended, 5th pct | median | 95th pct |
| ---: | ---: | ---: | ---: |
| 0.90 | 0.35 | 2.69 | 8.68 |
| 0.50 | 0.25 | 2.12 | 8.39 |
| 0.10 | 0.75 | 1.69 | 6.88 |

The median instruction buys about twice the incision it asks for and the
population spans a factor of thirty, so it is not a scale factor to divide out.
Carried through to the impoundment that would survive, using this build's own
geometry as the estimate of the next generation's:

| intended retain | delivered, 5th pct | median | 95th pct |
| ---: | ---: | ---: | ---: |
| 0.75 | 0.51 | 1.00 | 1.00 |
| 0.50 | 0.16 | 0.83 | 1.00 |
| 0.25 | 0.01 | 0.47 | 1.00 |

The medians at 1.00 are the other half of it: where erosion has already taken
the rim below the allowance floor, the allowance does nothing at all and the
basin keeps whatever erosion left. So a marginal instruction is not
systematically overspent or underspent; it is not transmitted.

This bites the marginal class alone. Retain 1 gives an allowance of zero and is
bit-identical to no allowance, and retain 0 clears `noLower` outright; both
cross exactly. The marginal class is the one the retain machinery exists for --
`hydrography/notes/retain-fraction.md` calls the notched through-flowing valley
"one of the few landforms this project decides rather than inherits" -- and it
is the only class the crossing does not carry.

The estimate above holds this build's geometry still and asks what an allowance
would have done to it. The next generation's terrain is not this one's, so the
numbers are the magnitude of the mismatch rather than a prediction of a
particular basin. What is a property of the conditioning rather than of the
carve list, and so does transfer, is the natural-to-final ratio in the table
above.

## The basis the interface now declares, and what the conversion closes

The carve list's retain is declared to be a fraction of the depression's
NATURAL relief in the generator's dimensionless elevation parameter. That is
the basis `buildBasinProtection` spends `(1 - retain)` against and the only
depth that exists when protection is built, and it is now stated on the format
in `vendor/orogen/tools/README.md`, on the two functions that spend it in
`vendor/orogen/js/basins.js`, and in the carve list's own header.
`export_carve_list.py:to_natural_relief_basis` puts the verdict into it before
writing, multiplying the incision fraction by `retainedFraction`.

Measured on 2026-08-25 from `basins.preserved[].finalPreserved.retainedFraction`
and `basins.preserved[].natural.spillDepth` in each manifest, as the ratio of
the incision Orogen is permitted in MODEL UNITS to the incision intended:

| build | ratio before, p05 / median / p95 | after, p05 / median / p95 |
| --- | ---: | ---: |
| precarve-craton | 0.75 / 1.26 / 2.71 | 1.0000 / 1.0000 / 1.0000 |
| precarve-craton-10m | 0.73 / 1.24 / 2.65 | 1.0000 / 1.0000 / 1.0000 |

The before column is `1 / retainedFraction` and does not depend on the retain
sent; the after column is exact on every basin whose instruction the mechanism
can express, which is the identity `export_carve_list.py --selftest` asserts.
The remaining spread is the two ends the mechanism cannot express, and both are
countable rather than argued.

**Saturation, and it is confined to the deep end of the marginal class.** An
intended cut deeper than the whole natural relief cannot be an allowance, and
the closest instruction Orogen has is carve. It requires
`retainedFraction > 1 / (1 - retain)`, so it reaches only the basins the
conditioning left DEEPER than the catalogue measured them, and only where the
verdict is already cutting most of the depression. On `precarve-craton-10m`, of
9,401 basins with a measurable ratio: none saturate at retain 0.90, 94 at 0.50
and 1,103 at 0.10. On `precarve-craton`, of 3,618: none, 35 and 381. That is a
verdict class moving, so it is reported per pass in the sidecar's
`carve_list_basis.verdict_class_moved` rather than absorbed.

**The residual is a bracket because the ratio is an estimate.** The conversion
uses this generation's per-basin `retainedFraction` as the estimate of the
next generation's, which is defensible because the ratio is a property of the
conditioning rather than of the carve list. Its population spread is what
bounds that: `retainedFraction` runs 0.369 / 0.797 / 1.340 at the 5th, 50th and
95th percentile on `precarve-craton` and 0.378 / 0.804 / 1.372 on
`precarve-craton-10m`, so a basin whose next-generation ratio lands at the
population's tails instead of its own median gets between 0.46 and 1.68 times
the intended incision. That is the honest residual, and it is a sixth of the
factor of thirty the mismatch carried.
