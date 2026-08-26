# The relief curve's domain: what the clamp at model elevation 1 was

Vesper is a generated world. This note is about the curve World Orogen uses to
turn its dimensionless elevation parameter into kilometres, and every number
here is a property of that generator and of the two registered builds of its
terrain.

The question `world-34gu` asks is whether `elevToHeightKm`'s `Math.min(elev, 1)`
is a deliberate physical ceiling on relief or an artefact of a curve fitted on
[0, 1]. It is the second, and the curve says so without needing an argument
about geology.

## What the curve is

The land branch is a scale times a shape function:

    h(t) = RELIEF_KM_PER_UNIT * s(t),   s(t) = t^4 (5 - 4t),   s'(t) = 20 t^3 (1 - t)

with `RELIEF_KM_PER_UNIT` = 6. s is a Hermite interpolant on [0, 1]. Its four
defining conditions are all at the two endpoints: s(0) = s'(0) = s''(0) = 0,
which is what gives the extensive flatlands, and s(1) = 1 with s'(1) = 0, which
is the peak compression the source comment names. A polynomial pinned that way
is a statement about the unit interval, not about the terrain's range.

Three things follow, and each is checkable rather than arguable:

- s'(t) is NEGATIVE for t > 1, so the polynomial turns over. It falls back
  through 6 km immediately above 1, reaches zero at t = 1.25, and is below sea
  level beyond. `Math.min(elev, 1)` was the guard against that inversion.
- The generator already names [0, 1] as the shape's domain elsewhere.
  `applyDetailNoise` in `js/terrain-post.js` refuses to invert the same
  polynomial at `elev >= 0.99`, with the reason recorded on the line: above
  elevation 1 the formula corrupts downstream erosion through NaN propagation.
- The physical ceiling on land relief exists in this codebase and is a different
  object. It is the 1/g factor in `scaledHeightKm`, out of sigma/(rho g), which
  is applied AFTER this curve and is argued from crustal strength. Reading the
  domain edge as a second ceiling double-counts one and misattributes the other.

Nothing bounds the model's elevation parameter at 1 either. `applyFinalShaping`
normalises land by RANK between the land field's own minimum and maximum, peak
compression is a power rather than a clamp, and detail noise, domain warping,
erosion, ridge sharpening and soil creep all run afterwards.

## What the clamp did, measured

MEASURED 2026-08-26 against the two registered builds that carry a payload,
`precarve-craton` at 2,500,001 regions and `precarve-craton-10m` at 10,000,005,
from the T42 export's raw mesh arrays in region order. Land is `surface_class`,
and area shares are area-weighted by `cell_area`. BOTH BUILDS ARE PRE-CARVE.

| build | max elevation, model units | land at or above 1, by area | by count |
| --- | ---: | ---: | ---: |
| precarve-craton | 1.3687 | 0.2594% | 2,953 cells |
| precarve-craton-10m | 1.4280 | 0.2486% | 11,251 cells |

The parameter runs 37% and 43% past the shape's domain, so this is a substantial
excursion rather than an edge case, and the share is converged across a fourfold
change in region count. Both builds published a maximum land `elevation_km` of
4.5933 km, identical to four decimals, because both were the clamp times the
relief scale.

What the polynomial would have done unclamped, on those cells: 2,949 of 2,953
and 11,243 of 11,251 would publish BELOW the join value, and 49 and 134 would
publish below sea level. The clamp was load-bearing.

## The branch, and what it changes

Above the domain the shape is replaced by its argument, s(t) := t, so the branch
is `RELIEF_KM_PER_UNIT * elev`. It introduces no constant the curve did not
already have: it is continuous at t = 1 where s(1) = 1 already, its gradient is
the curve's own mean gradient across its domain, and it is strictly increasing
for every t. It is applied before `reliefScale`, so the terrain stays
gravity-invariant.

What it costs is a kink. The gradient at the join steps from 0 on the left,
which is exactly the peak-compression condition s'(1) = 0, to 6 km per model
unit on the right. The same function already carries a kink of that shape at sea
level, where the land branch arrives with slope 0 and the ocean branch leaves
with 10 km per unit, and the kink replaces a gradient of exactly zero over all
the terrain it affects.

Recomputed from the same `elevation` arrays:

| build | max land elevation_km, clamped | with the branch | cells whose height moves | distinct heights among them, before / after |
| --- | ---: | ---: | ---: | ---: |
| precarve-craton | 4.5933 km | 6.2867 km | 2,953 | 1 / 2,952 |
| precarve-craton-10m | 4.5933 km | 6.5594 km | 11,251 | 1 / 11,188 |

The median cell in that band gains 239 m and 206 m. Both maxima land under
Earth's own maximum relief divided by this planet's gravity ratio, which is a
check on the branch rather than a criterion it was chosen to meet: the branch
was derived from the curve before any of this was computed.

Every number above is reproducible from `lib/orogen.py`'s
`Export.field('elevation')`, `Export.field('surface_class')` and
`Export.field('cell_area')` on the T42 export of each build, evaluating the two
curves on the land cells and area-weighting the shares.

## It moves the terrain, through selection

`selectBasins` compares `BASIN_MIN_DEPTH_KM` against `selectionDepthKm`, which
is this curve applied to a depression's spill and its sink at reference gravity.
The preserved set it returns reaches `r_elevation` through `inciseOutlets` and
`buildBasinProtection`, both inside `runPostProcessing` and both before erosion.
So the curve is not only a publication rule.

Over the whole catalogue, counting depressions with a sink or a spill above
model elevation 1:

| build | catalogue | touched by the branch | cross the 0.05 km depth floor upward | downward |
| --- | ---: | ---: | ---: | ---: |
| precarve-craton | 81,871 | 199 | 108 | 0 |
| precarve-craton-10m | 286,534 | 667 | 390 | 0 |

The move is one-directional by construction: the branch only raises heights
above the join, and a spill is at or above its sink, so a depression's depth can
only grow. Those counts are candidates on the depth floor alone; the area floor,
the cell floor and the nesting rule still apply, so the preserved set gains at
most that many.

Basin ids are unaffected, being built from the sink region index and the member
set on the pre-conditioning surface, so a carve verdict keyed by id still
resolves. `hashes.basinCatalogue` does move, because `depthKm`,
`sinkElevationKm` and `spillElevationKm` are published fields inside it.

`computeScarpPotential` also reads the curve, and `scarp_potential` moves with
it, but that call runs after `runPostProcessing` and feeds no erosion.

## What this is NOT the cause of

`hydrography/data/<build>/basins.nc` gives `capacity_km3` of exactly zero for 5
preserved basins on `precarve-craton` and 16 on `precarve-craton-10m`. The clamp
accounts for 2 of the 5 and none of the 16.

MEASURED 2026-08-26 from `basins.preserved[]` in each manifest and from
`sink_elevation_km`, `spill_km` and `capacity_km3` in `basins.nc`. On
`precarve-craton` exactly 2 preserved basins have both sink and spill above
model elevation 1, and those 2 are the export's only preserved basins with
`depthKm` of exactly zero; on `precarve-craton-10m` there are none. The
remaining 3 and 16 have their sink and their spill at ordinary heights, from
0.003 to 2.671 km, equal to each other to the precision the file carries, over
areas of 47 to 528 km2 -- one or two mesh cells on whichever build each belongs
to.

`notes/audits/basin-catalogue-floor.md` reached the same split independently and
attributes those basins to `finalPreserved.retainedFraction` of zero: the
drainage conditioning breached the rim, so the depression the catalogue measured
does not exist on the finished terrain. That is a separate fault with a separate
cause, and the branch here does not touch it.

## Verdict

The clamp was an artefact of the shape function's domain, not a physical ceiling
on relief. `js/color-map.js` now carries a branch above 1, the derivation is
stated where the curve is declared, and
`vendor/orogen/tools/test-basins.mjs` asserts strict monotonicity from 0.001 to
2.0 -- an assertion the clamped implementation fails at 1.001.

Neither build in `source/` inherits it. The curve reaches a build only through a
generation, and this project's terrain changes are batched into one.
