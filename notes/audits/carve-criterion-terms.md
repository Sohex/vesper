# Audit: the four terms of the carve criterion, and the clamp that was deciding it

*Measured 2026-08-17 on `precarve-craton` under `baseline_regular_climatology.nc`,
after HYD-13 and PHYS-5. Unlike the two audits it follows, this one FIXED what it
found; every number below is measured on the artifacts, and the state before each
fix is given beside the state after.*

Four open tasks all moved the same quantity, so they were settled together:
HYD-15 the retain mapping, HYD-14 the seasonal rectification, HYD-11 the dry
column, HYD-12 the tolerance. A fifth thing turned up in the middle of it, was
larger than any of them, and is finding 1.

Findings are tagged **[numeric]** where computed here and **[inspection]** where
read out of code.

---

## 1. A clamp, not the climate, was deciding 73% of the overflowing basins

**[numeric]** `carve_verdict.py` floored the open-water estimate at the model's
own land evaporation:

    penman = np.maximum(penman, evap)

Its stated reason was that a saturated surface cannot evaporate less than the
moisture-limited ground beside it, so a Penman below `evap` had to be Penman
failing. It bound on **2,472 of 4,105 land cells, 60.2%**, and where a whole
catchment was floored the aridity index came out at exactly -1.0 and the basin
carved regardless of its climate: **1,252 of the 2,543 basins with catchment
runoff**, 73% of the overflowing set. The floor was added to fix 43 basins.

**It was not catching a failure. It was catching the wet regime.**

| | where the floor bound | where it did not |
| --- | ---: | ---: |
| mean soil wetness | 0.597 | 0.121 |
| mean model evaporation, mm/day | 1.731 | 1.310 |
| mean Penman, mm/day | 1.064 | 3.212 |

It fired on 1,324 of the 1,347 land cells at full wetness. The mechanism is
roughness, and it is in the boundary conditions this project supplies: land
carries a roughness field with a median `z0` of 0.521 m against open water's
1.5e-4 m, which at this reference height is a transfer coefficient **6.4 times
larger**. A smooth lake in a rough, wet, vegetated landscape genuinely
evaporates less than the land around it. The nesting the floor defended -- that
`wet` is always a lower bound on lake evaporation -- holds only where the ground
is moisture-limited, and reverses where it is not.

**The floor also contaminated the check that was supposed to catch it.** Over
ocean, `max(penman, evap)` is bounded below by the model's own answer, so the
ratio could not come out under 1 whatever Penman did. Unfloored it reads 1.0845;
floored it read 1.0909, and 1.0073 once the level fix below was in. A validation
that cannot fall below 1 is not measuring the thing it names.

The floor is removed. `wet` is now described as a sensitivity whose side depends
on the regime, not as a bound, and the three-way partition in
`export_carve_list.py` that assumed nesting is replaced by an agreement count --
it had an assertion that the three classes summed to the catalogue, and that
assertion failed the moment the floor went, which is the check working.

Found by the error-budget session, which reconstructed the criterion from the
carve list and noticed `aridity_index_penman` was exactly -1.0 on 1,252 basins.

## 2. Penman read its two humidity terms at two different heights

**[numeric]** `es_a` was saturation at `tas`, the 2 m temperature, while `e_air`
came from `hus` at the lowest model level, which `SIGMA_LOWEST` puts of order
300 m up. A Penman deficit is `e_s(T_a) - e_a` at ONE reference height, and the
transfer coefficient this scheme builds is derived over `z_ref`, so that height
is the lowest model level.

The 2 m air is **1.85 K warmer** than the lowest level over land, worth **13.7%**
on `e_s`. Taking the saturation term low and warm and the actual term high and
dry inflates the deficit, one-signed.

Evaluating everything at the lowest model level -- temperature, humidity, wind
and pressure -- is the fix, and the ocean is where it can be tested, because an
ocean cell already is the surface Penman is written for and the model gives it
water's own roughness:

| Penman against the model's own ocean evaporation | ratio | land mean, mm/day |
| --- | ---: | ---: |
| levels mixed, as it stood | 1.0845 | 2.974 |
| one reference level | **0.9672** | 2.866 |

The error falls from 8.45% to 3.28% and changes sign. It is not tuned to 1.000
and must not be: the residual is the method's, and part of it is known -- the
diurnal range is `maxt - mint` at 2 m, applied at 300 m where the real cycle is
damped, which biases Penman UP and so cannot explain a low ratio.

There is a second reason the model level is the right height rather than merely a
consistent one. A sub-grid lake does not sit under the 2 m air the model reports,
which is a diagnostic of the dry ground beside it.

## 3. The sub-grid dry column is real, is much smaller than it was priced at, and stays a bracket

**[numeric]** `missed-couplings.md` finding 2 priced "Penman over a dry land
column" at 10-18% of land-mean Penman and 144 to 235 basins. Both halves of that
were inflated, for two separate reasons.

Relative humidity at the lowest model level is already 68.3% over land, so a
floor at 70% barely bites. Measured on the corrected scheme:

| assumed column relative humidity | land-mean Penman, mm/day | change | basins carved |
| --- | ---: | ---: | ---: |
| as computed | 2.866 | -- | 1790 |
| at least 70% | 2.863 | -0.1% | 1791 |
| at least 80% | 2.745 | -4.2% | 1799 |
| at least 90% | 2.586 | -9.8% | 1825 |

So the item is worth **1 to 35 basins**, not 144 to 235. The gap has two causes,
and both are worth stating because the same arithmetic will be attempted again:
the level fix in finding 2 has already removed most of what the dry-column
bracket was measuring, and the 144-per-10% sensitivity it was multiplied by was
itself measured with the floor of finding 1 in place, which inflates any
evaporation sensitivity by about a factor of four.

**It stays a bracket and is not applied.** A lake tens of kilometres across grows
an internal boundary layer deeper than this reference level, so the air over a
real lake is wetter than the cell mean -- but how much is a fetch problem, and a
humidity floor chosen to move the answer would be a knob. `penman_open_water`
takes `column_relative_humidity` for exactly this bracket and no product passes
it.

## 4. Clamping runoff per season would count the soil store twice

**[numeric]** HYD-14 asked for a bracket on the seasonal rectification in
catchment runoff, on the argument that the annual mean of `P - E` charges a
catchment's dry-season deficit against its wet-season supply. `grid-convention-and-runoff.md`
finding 3 measured the three treatments and put the truth inside [161, 266] mm/yr.

**The upper end is not a bound. It is a double count, and the model's own soil
water says so cell by cell.**

Over one annual cycle at steady state a land cell's storage returns to where it
started, so `annual(P - E)` equals the runoff that cell generated, exactly. It is
an identity, not an estimate. The negative bins are the store being drawn down;
the water that refills it in the wet season is water that did not run off. Adding
the negative bins back counts the wet season twice.

Land means per Vesper year, on the current baseline:

| | mm |
| --- | ---: |
| annual mean of `P - E`, the conserved answer | 63.3 |
| sum of the positive bins, which per-bin clamping would use | 114.6 |
| the dry-bin deficit it discards | 51.3 |
| seasonal range of the model's own soil water plus snow | 62.8 |

Cell by cell, the ratio of that storage range to that deficit has a **median of
0.988**, with quartiles 0.95 to 1.07. The store IS the deficit. That is the
check: a bucket cannot supply water it never stored, and it cannot generate
negative runoff.

Reconstructing the model's own local generation as `P - E - dW/dt` confirms it
from the other side. Its annual mean is `annual(P - E)` by construction, and
clamping THAT per bin gives 76.0 mm against 114.6 for clamping `P - E` -- the
residual being the 12-bin centred difference smoothing sub-bin cycling, which
also makes the reconstruction go negative on 39.6% of bins, an impossibility for
a bucket and therefore a measure of the smoothing rather than of the water.

**The clamp stays after aggregation** rather than moving per cell. A negative
annual `P - E` over pure land is not physical here -- the mask is binary, so
there are no part-ocean cells, and this model's routed runoff never re-enters the
evaporating bucket, so no cell can re-evaporate imported water. It appears on
1,597 of 4,105 land cells with a mean of -3.1 mm/yr and a worst case of -39,
which is the residual non-periodicity of a five-orbit window. Clamping a
zero-mean noise term at zero biases the estimate up by half of it; doing so at
the catchment level, where the noise has already averaged down, biases it least.
Per-cell clamping would move the catchment mean by 0.6% and reclassify 407
basins, for an operation whose own bias is three times the change it makes.

**The seasonal rectification that does exist is in the lake, not the catchment,
and it is small here.** A lake near its spill can be pushed over it by a wet
season even when its annual balance sits below -- and carving is irreversible, so
that ratchets, by exactly the argument `docs/src/pipeline/sequencing.md` A2 makes for the stellar
cycle. Measured: of the 1,752 basins below their spill at annual equilibrium, the
seasonal storage swing exceeds the volume gap to the spill on **8**, 0.5%. Median
swing-to-gap ratio 0.004, ninetieth percentile 0.046.

HYD-14 is therefore closed without a change to the denominator. The criterion
already uses the conserved quantity.

## 5. The sill's slope is the basin's depth, and at n = 1 the depth cancels

**[numeric]** HYD-15. The retain mapping folded `S^n` into the coefficient, which
asserts that a sill's gradient is independent of the depth behind it.

`retain-fraction.md` recorded that the slope "is not in the export". The saddle's
own slope is not, and could not be: a saddle is a saddle, and the cross-divide
drop between `spill_region` and `spill_exit_region` on this build has a median of
0.1 m with half the basins negative. But stream power is not written about the
divide, it is written about the channel the overflow runs DOWN, and that channel
IS in the export. Following the drainage ten receiver steps from
`spill_exit_region` covers a median of 135 km, a dozen mesh regions.

Measured that way, over the 2,937 basins whose outflow path is at least 50 km:

| | |
| --- | ---: |
| `r(log S, log depth_at_spill)` | **0.735** |
| `d log S / d log depth`, ordinary least squares | 1.04 |
| the same by reduced major axis | 1.41 |
| `r(log S, log Q)`, the confound check | 0.072 |
| 5-step against 20-step estimate of S, same basins | r = 0.903, median ratio 1.00 |

The outlet gradient and the basin's depth are the same relief counted twice. With
`S ~ depth^p`, the effective depth exponent in `retain` is `n*p - 1`, which at
Orogen's own n = 1 and p = 1.04 is **+0.02**: the depth cancels. That is a
different mapping, not a different coefficient, which is why it had to be settled
rather than bracketed.

It is also new information rather than a proxy for the water: against discharge
the gradient measures r = 0.072.

`S` is now measured per basin and normalised by the population's geometric mean,
so the coefficient keeps the meaning its Earth calibration gives it.

## 6. The incision coefficient was stale, and being a literal is why

**[numeric]** The coefficient is calibrated by matching the density of Earth's
standing through-flowing impounded basins. **How many basins overflow is a
property of the climate**, so the target moves whenever the verdict does -- and a
literal cannot.

The published 161 was fitted against a verdict taken before HYD-13. Against the
verdict that replaced it, 161 leaves **55** standing basins within 35 degrees
where the Earth density asks for 33.6. Recalibrating with no other change gives
235.9, outside the published 128-218 bracket.

It is now solved on every run by bisection, from `EARTH_STANDING_BASINS = 15` over
`EARTH_BAND_LAND_MKM2 = 78.9`, against this world's own land area in the same
band, which the export puts at 176.9 Mkm2 for a target of 33.6.

Decomposed on the pre-floor-removal verdict, so that the mapping change is
separated from the evaporation changes:

| | coefficient | carve | marginal | preserve | standing < 35 deg |
| --- | ---: | ---: | ---: | ---: | ---: |
| as published, no slope term | 161.0 | 1598 | 123 | 1900 | 55 |
| recalibrated, still no slope term | 235.9 | 1648 | 73 | 1900 | 33 |
| slope term, n = 0.5 | 94.8 | 1628 | 93 | 1900 | 33 |
| **slope term, n = 1, adopted** | **102.9** | **1608** | **113** | **1900** | **34** |
| slope term, n = 2 | 330.2 | 1601 | 120 | 1900 | 33 |

n = 1 is Orogen's own and is not a free choice here: GRAV-4 records that it is
baked into the Braun-Willett closed form the generator solves, and the post-hoc
1/g relief scaling is only correct at that value.

**The counts barely move and the composition changes completely**, which is the
result, and it is the same shape the previous revision of this mapping produced.
The marginal class had a median depth at spill of 1,360 m against a population
median of 127; it now has 126, which is the population. The depth selection is
gone, exactly as the algebra says it should be. What selects a marginal basin now
is a flat outlet: median outlet gradient 0.00015 against the population's 0.00068,
four and a half times flatter. 180 basins change side.

The marginal class does not collapse to the 30 of the discharge-only mapping,
because the outlet gradient has more spread than depth does. A basin is now
marginal because its overflow leaves down a gentle channel, which is the right
reason for a notched valley to survive.

## 7. The tolerance stands at 0.25, and it was the wrong thing to be asked about

**[inspection]** HYD-12 asked whether `TOLERANCE = 0.25` is still right now that
Penman's ocean error is measured rather than hardcoded. It is, and the reasoning
is unchanged in form: the tolerance is set by the uncertainty that dominates, and
Penman's method error is the smallest term in the list, not the largest.

| term | size on lake evaporation |
| --- | --- |
| the assumed biosphere | 3.7 to 7.1 K |
| dust | 10 to 20% |
| the sub-grid dry column, finding 3 | 0 to 9.8% |
| Penman's own method error, finding 2 | 3.28% |

Tying the tolerance to the method error would let the smallest term set the scale
for the largest.

The framing is what was wrong. `retain_margin` has never decided a single basin
and cannot: `retain` is the larger of it and the incision value, the margin is
positive exactly when the basin does not overflow, and that is exactly when the
incision value is already 1. Re-verified under the new mapping, 0 of 3,621. The
term that WAS deciding basins by fiat was the floor in finding 1, on 73% of the
overflowing set, and it was one line away in the same file.

## 8. Two more instances of the antipodal grid defect

**[numeric]** GRID-1 consolidated the export/model column convention into
`lib/gridding.py`. Applying the same invariant that caught HYD-13 to every
remaining site found two more, both of which survived that sweep:

`surface_water.py:per_basin_forcing` sampled each basin's lake precipitation and
open-water evaporation at its sink with `col = mod(round((sink_lon - lon[0]) /
dlon), nlon)`, measuring an Orogen longitude from the CLIMATOLOGY's first label.
It is 65 columns out of 128. **52.47% of basin sinks landed on cells the model
calls ocean, against 3.04% index for index**, and the median relative change in
the sampled field is 38 to 46%. Every lake on this planet was taking its
precipitation and its evaporation from its antipode. It was named in
`grid-convention-and-runoff.md` finding 1 and was not fixed with the rest.

`maps/build_basemap.py:upsample` resampled climate fields onto the render grid by
interpolating against the climatology's own 0..360 longitude labels while the
render grid runs -180..180. It is 64 columns out, so every climate-derived layer
on the basemap -- the biome tint, the ice, the sea ice -- was drawn half a world
from the terrain beneath it. Nothing computed reads the basemap, so no result was
wrong; the map was.

`scripts/smoke_test.py` now lints for longitude arithmetic outside
`lib/gridding.py`, with `maps/projections.py` and `maps/render_projections.py`
exempt because they rasterise Orogen's own coordinates onto image pixels and no
model grid is involved. Anything that reads a climatology or a coupling matrix is
on the seam and is not exempt, which is how the basemap was caught.

---

## The combined effect, which is not the sum of the parts

The four evaporation and runoff changes move the same denominator, and three move
it the same way. Measured on the ratio form of the criterion, which decides the
2,559 basins with catchment runoff:

| | carve | ocean ratio | land Penman, mm/day |
| --- | ---: | ---: | ---: |
| as it stood: levels mixed, floor in | 1721 | 1.0909 | 3.256 |
| one reference level, floor still in | 1770 | 1.0073 | 3.205 |
| levels mixed, floor removed | 1739 | 1.0845 | 2.974 |
| **one level, floor removed: adopted** | **1790** | **0.9672** | **2.866** |

The two deltas measured separately are +49 and +18; measured together they are
+69, so they do not sum, and 69 basins differ between the first row and the last.

**The floor gets worse under the level fix, not better.** With the level fix in
and the floor kept, the number of basins whose index sits at exactly -1.0 rises
from 1,252 to 1,499: a correction that lowers Penman drives more of the catchment
onto the clamp, so the clamp absorbs the correction. Either had to be settled
before the other could be measured at all.

Through the full carve list, which also decides the 1,062 basins with no
catchment runoff, and with the mapping and coefficient of findings 5 and 6:

| | carve | marginal | preserve |
| --- | ---: | ---: | ---: |
| the list currently on disk | 1598 | 123 | 1900 |
| after all of the above | **1938** | **95** | **1588** |

243 of the new carves are basins with no catchment runoff at all, whose own lake
surface gains more precipitation than it evaporates. There were none before,
because the floor put lake evaporation at the land rate. They sit at a median
latitude of 70.3 degrees, where a lake evaporates 185 mm/yr against 269 mm/yr of
precipitation: cold basins fill. That is the case the discharge form of the test
was rewritten to handle, and this is the first climate in which it fires.

## The no-runoff population

HYD-17 (archive/tasks.md) refuted the window hypothesis: the count is noise
about a constant across one- to ten-orbit windows. What is not periodic is the
soil store; the measured land-mean drift is 8.6 mm per bin.

---

## Tasks

Tracked in `TASKS.md`, not restated here.

| finding | id |
| --- | --- |
| 1. the land-rate floor | `HYD-11` |
| 2. the level mixing | `HYD-11` |
| 3. the sub-grid dry column | `HYD-11` |
| 4. seasonal rectification is a double count | `HYD-14` |
| 5. the sill slope is the basin's depth | `HYD-15` |
| 6. the incision coefficient was stale | `HYD-15` |
| 7. the tolerance | `HYD-12` |
| 8. two more antipodal grid sites | `GRID-1` |
