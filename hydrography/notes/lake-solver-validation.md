# Validating the lake solver against Earth

HYD-4. Registered 2026-08-17, **before any data was fetched**, in the same way
`pedology/scripts/validate_against_earth.py` and the albedo bracket's criteria
were. What follows above the line is the test; what follows below is what it
found.

## Why this test and not another

The lake solver is the least-checked product in this pipeline, and it is load
bearing twice over. It sets how much of the land is under water, which is worth
-0.0139 on land-mean albedo and therefore reaches the climate; and it is the only
thing that could break the degeneracy between salt crust extent and salt crust
albedo, which `exoplasim/notes/parameter-decisions.md` records as unresolvable
from the climate side alone.

Penman was validated the same way and that is the standard being matched here:
applied to ocean cells, which already are open water, it reproduced the model's
own evaporation to 1.7%. The lake solver has had nothing equivalent.

## What is under test, and what is not

The solver has three parts and only one of them is testable on Earth.

**Testable: the equilibrium relation.** A terminal lake grows until net
evaporation over its surface consumes what its catchment delivers:

    R * (A_catch - A_lake) = (E_lake - P_lake) * A_lake

which rearranges to a prediction that needs no hypsometry at all:

    A_lake / A_catch = R / (R + E_lake - P_lake)

**Not testable here: the hypsometric inversion**, which converts a solved area to
a level and a volume through curves measured on this world's terrain. Earth
basins have their own hypsometry and testing it would test the data, not the
code.

**Not testable here: the overflow cascade**, for the same reason. What can be
said about it is already known -- it reports whether it reached a fixed point,
and a non-converged solution is refused rather than used.

## The prediction, and the thresholds

Stated before fetching anything.

The relation above is scale free, so the test is of the *ratio* rather than of an
area. Earth's terminal lakes span four orders of magnitude in size, so a relation
that is right in form should hold across all of them or fail visibly.

- **PASS** if the median absolute log-ratio of predicted to observed lake area is
  below 0.30, i.e. predictions typically within a factor of 2.
- **MARGINAL** between 0.30 and 0.60, a factor of 2 to 4.
- **FAIL** above 0.60, or if the error correlates with basin size, aridity or
  lake area, since any of those would mean the relation is missing a term rather
  than being noisy.

A factor of 2 sounds loose and is not. The quantity that reaches the climate is
lake *area* as a fraction of land, and the albedo consequence of getting it wrong
by a factor of 2 is smaller than the playa albedo uncertainty already carried at
1.99 W/m2. What would matter is an order of magnitude, or a bias.

## What will make this noisy, stated in advance

**Earth's terminal lakes are mostly not in equilibrium**, and the ones people
have measured best are the ones humans have most altered. The Aral Sea and Lake
Chad are the two most cited endorheic basins on Earth and both are dominated by
twentieth-century irrigation withdrawal; their present areas measure water
management, not climate. Where a pre-development balance is published, that is
the number to use, and the choice has to be made per basin and stated.

**Lake area is not a smooth function of climate.** A basin's hypsometry can be
flat, so a shallow terminal lake swings enormously in area for a small change in
volume. Lake Eyre is the extreme case, dry in most years and 9,500 km2 in a wet
one. A single-year area is the wrong observable; a long-term mean is the right
one and is not always available.

**Runoff over the catchment is the least certain input**, because for most of
these basins it is itself modelled rather than gauged.

None of that is a reason to skip the test. It is a reason to report the scatter
and to say which basins carry it.

## Method

1. Choose basins on whether a defensible water balance exists, before looking at
   whether the prediction works for them.
2. Take catchment area, lake area, catchment runoff, lake evaporation and lake
   precipitation from the literature, one source per basin where possible.
3. Compute the predicted ratio from the relation above and compare.
4. Report every basin, including the ones that miss.

---

## Results

**FAIL, on the criterion registered this morning.** Median absolute log-ratio of
predicted to observed lake area is **0.622** against a pass at 0.30 and a fail
above 0.60, over 145 scored lakes. Only 30% land within a factor of two. Run by
`hydrography/scripts/validate_lake_solver.py`; the per-lake table is in
`hydrography/analysis/lake_solver_validation.json`.

The registered criteria also failed it a second way, which matters more than the
median: **the error correlates with runoff at +0.760** and with catchment area at
+0.467. The threshold said a correlation means a missing term rather than noise,
and it does.

### The failure is two-sided, and each side names its own missing term

By runoff quartile, which is a post-hoc split and labelled as one:

| catchment runoff | n | median abs log10 | bias |
| --- | ---: | ---: | ---: |
| 0.0 - 5.7 mm/yr | 37 | 0.888 | **-0.888** |
| 5.7 - 21.6 | 37 | 0.374 | +0.052 |
| 21.6 - 67.5 | 37 | 0.443 | +0.428 |
| 67.5 - 467.6 | 37 | 0.734 | **+0.734** |

**The wet tail over-predicts by about fivefold, and the term it is missing is the
spill cap.** The relation lets a lake grow until evaporation consumes its supply.
A real basin stops at its sill and passes the rest downstream, so its area is set
by hypsometry rather than by evaporation. `lake_balance.solve` HAS that cap; the
relation tested here does not, because the whole point of the relation was to be
testable without hypsometry.

**The dry tail under-predicts by about eightfold, and the terms it is missing are
groundwater and history.** The relation says a nearly dry catchment supports a
nearly absent lake. Observed lakes in those basins are far larger, and two
reasons are already in this project's own reading: Issyk-Kul takes roughly a
third of its input as groundwater recharge, which no surface balance sees, and
Earth's arid terminal lakes are substantially relicts of wetter climates rather
than bodies in equilibrium with today's. The design section above warned about
exactly that, in advance, and it turns out to dominate a quarter of the sample.

### What this does and does not say about the solver

**It does not validate it.** That was the point of running it and the answer is
no.

**It does not cleanly invalidate it either**, and saying so is not special
pleading: the wet tail fails on a term the solver implements and the test could
not, and the dry tail fails on terms neither has and Earth's own lakes violate.
What the test bounds is the equilibrium RELATION, not the code.

**Where the relation works is the middle**, at a bias of +0.05 across the second
quartile, and even there the scatter of 0.374 misses the registered pass
threshold. So the honest statement is that this relation is good to about a
factor of two in the regime where a lake is genuinely evaporation-limited, and
not usable outside it.

**The regime that matters on Vesper is the one this test cannot probe.** Of 2,465
water-holding basins there, 1,737 sit at their spill: geometry-limited, where the
cap does the work and the equilibrium relation is not what sets the area. Only
728 are evaporation-limited, which is the population this test speaks to.

So HYD-4 closes with the solver still uncertified, and with something better than
the unqualified "never validated against anything" it started with: a measured
bound on the relation it iterates, a measured direction of failure at each
extreme, and the knowledge that most of this world's lakes are set by a cap
rather than by that relation.

### What would test the cap

Hypsometry. A basin's area-at-spill against its solved area is the quantity, and
Earth's endorheic basins do have measured hypsometry in places -- the Great Basin
pluvial lakes and the Tibetan closed basins are both mapped. That is a different
and larger piece of work than this one, and it is the only route to testing the
half of the solver this world actually exercises.

---

## HYD-7: the spill cap, and what a 15 km representation costs

Measured 2026-08-17, on Copernicus DEM 90 m from the AWS Registry of Open Data.
HYD-4 could not reach the cap because it needed hypsometry; this does, by asking
the question the cap actually turns on.

### The question, narrowed

The cap engages when a basin's storage is exhausted, so what has to be right is
`area_at_spill` and `capacity`. Those come from flooding the terrain at the mesh
scale of about 15 km. The test is therefore not "is the algorithm correct" but
**what does representing a real basin at 15 km do to its measured storage**, and
in which direction.

Two real endorheic basins, deliberately in different settings: the Qaidam in
Tibet, a broad pan, and the Great Salt Lake basin in the Basin and Range, a
graben. Each DEM is used at its native 90 m and then block-averaged to the mesh
scale, and the flooded area and volume compared level by level.

| depth above floor | Qaidam area | Qaidam volume | Great Salt area | Great Salt volume |
| ---: | ---: | ---: | ---: | ---: |
| 25 m | 0.00 | 0.00 | 0.63 | 0.55 |
| 50 m | 0.00 | 0.00 | 0.76 | 0.64 |
| 100 m | 0.65 | 0.51 | 0.76 | 0.71 |
| 200 m | 0.72 | 0.66 | 0.79 | 0.77 |
| 400 m | 1.03 | 0.87 | 1.00 | 0.87 |

Ratios are coarse over native, so below 1 means the mesh scale sees less.

### The result, and it is one-signed

**A 15 km representation understates storage, by 13% at 400 m of depth and by 30
to 50% in the shallow range, and on a broad pan it can miss shallow flooding
entirely.** The Qaidam shows literally no flooded area at 25 or 50 m above its
floor, because block-averaging fills the depression's own floor with the ridges
around it.

That direction matters more than the magnitude. Understated storage means a given
water supply fills a basin sooner, and understated `area_at_spill` raises the
critical aridity index, which is the threshold a basin overflows against. **Both
channels push the same way: toward spilling, and therefore toward over-carving.**

**The shallow end is where this world lives.** Vesper's basins have a median mean
depth at spill of 46.5 m, which sits in the 25-to-100 m band where the deficit is
worst rather than in the 400 m band where it converges.

### What it does not say

Vesper's terrain is generated at mesh scale, so it never had sub-grid depressions
to lose. This is not a bug in `build_hydrography.py` and it is not a correction
to apply. It is a statement about how much storage a real world of this relief
would have that the represented one does not, and it belongs in the error budget
as a one-signed bias on the carve verdict.

### The shape check, which passes

Separately, and cheaply: the area-volume scaling of the solved lakes matches
Earth's form. Fitting `V = c A^k` over lakes above 10 km2,

| population | n | exponent k | median mean depth |
| --- | ---: | ---: | ---: |
| Vesper solved lakes, all | 2032 | 1.259 | 34.4 m |
| Vesper, evaporation-limited only | 296 | 1.628 | 8.7 m |
| Earth endorheic terminal lakes | 146 | 1.168 | 3.0 m |
| Earth, all natural lakes above 10 km2 | 14591 | 1.208 | 6.7 m |

The exponent is the shape: how fast area grows as a basin fills, which is what
decides when the cap engages. Ours is 1.26 against Earth's 1.17 to 1.21, so the
basins fill like real ones do. The level differs -- our lakes are deeper for
their area -- and that is confounded between a real property of this terrain and
Earth's terminal lakes being shrunken relicts, which HYD-4 measured independently
as an eightfold under-prediction in the dry tail, in the same direction.

Do not compare Vesper's basins-at-spill against Earth's modern lakes: at spill
they are basins filled to the brim, at 0.799 and 46.5 m, and Earth's endorheic
lakes are mostly nothing like full. That comparison was run first and is a
category error.
