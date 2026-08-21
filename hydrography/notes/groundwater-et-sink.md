# GW-15: the internal sink, its design, and a prediction made before the test

Worldbuilding. Vesper is an invented planet; this designs a term for the
groundwater solver written for it, and tests that term against Earth.

Written 2026-08-20 after GW-3 missed both its thresholds. **The prediction in
section 4 is recorded before the term was implemented**, so the test can fail.

## 1. What GW-3 established

A steady-state water table whose only exits are the coast and surface seepage
cannot hold itself below the land surface across a dry continent. It pinned at
the surface on 95.5% of Australian land and 63.0% of Vesper's, and scored no
skill at all against 70,119 bores: residual standard deviation 40.06 m against
an observed 40.05, because a constant explains nothing.

The reason needs no model. For a confined aquifer draining recharge `R` over a
distance `L`, head above the outlet goes as `R L^2 / 2T`. At Australia's 8 mm/yr
over 1,000 km, holding the table 250 m down needs about 600 km of aquifer.
Recharge in a continental interior cannot reach the sea, so in steady state it
has nowhere to go but up.

**Every exit the model has is at the edge of the domain. Real arid groundwater
leaves in the middle.**

## 2. The term

Groundwater evapotranspiration: where the water table is shallow enough for
roots and the capillary fringe to reach it, water leaves upward. This is
channel 2 of the four in `groundwater-scoping.md` section 2, and it was filed as
GW-5, a one-signed omission. GW-3 reclassifies it as a precondition.

Shah, Nachabe and Ross (2007) is the source and it settles the FORM, which is
the part that matters here. Their finding, from Richards-equation simulations
across soil textures and land covers: **"The decline of ET with DTWT is better
simulated by an exponential decay function than the commonly used linear
decay."** MODFLOW's EVT package uses the linear form; this does not.

Their extinction depths run from 0.18 m for bare sand to 1.86 m for clay under
forest, and those are shallow. That number is what makes the prediction below
falsifiable rather than obvious.

    ET(d) = ET_max * exp(-d / lambda)

for water table depth `d` below the surface. The sink is head-dependent, so the
matrix stops being fixed and the problem is nonlinear again -- but MONOTONE, in
the direction that helps: the deeper the table, the weaker the sink, so a cell
cannot oscillate the way the exponential transmissivity did under GW-9.

## 3. Why this changes the shape of the answer, not just its offset

Setting the sink equal to the local recharge gives the equilibrium depth
directly:

    R = ET_max * exp(-d / lambda)      =>      d = lambda * ln(ET_max / R)

**The depth becomes a function of recharge.** That is the mechanism Fan et al.
(2013) name as the regional control -- "regions of deep WTD correspond to
regions of low recharge, the great deserts of the world stand out" -- and it is
exactly what the current model cannot produce, because it has no sink whose
strength depends on where the table sits.

The logarithm is the catch. Australian recharge spans 1 to 869 mm/yr, a factor
of 900, and a logarithm compresses that to a factor of 13 in depth.

## 4. The prediction, recorded before implementing

With a spatially constant `ET_max` of order 1,500 mm/yr, so that the ONLY source
of spatial variation is recharge:

1. **The surface pinning largely goes.** The at-surface fraction falls from
   95.5% to something under 20%.
2. **The model acquires skill it currently has none of.** Residual standard
   deviation falls below the observed 40.05 m, which the present constant output
   cannot do at any offset.
3. **And it still misses the declared thresholds**, because the deep tail
   survives: observed depth reaches 49 m at the 95th percentile, and
   `lambda ln(ET_max/R)` at Shah's lambda of order a metre cannot reach that
   without a recharge below about 1e-19 mm/yr.

**So the expected verdict is: necessary, insufficient, and diagnostic.** If (1)
and (2) fail, the term is not the missing physics and this note is wrong. If (3)
succeeds -- if the deep tail IS reproduced -- then the tail was never evidence
for the relict-water argument and that argument should be withdrawn.

## 5. What it would still not reach

Australian groundwater carries residence times of 1e4 to 1e6 years and is not in
equilibrium with modern recharge. A steady-state model cannot represent a water
table still draining from a wetter past, and `docs/src/reference/no-time-axis.md`
forbids the alternative. That is the same objection HYD-4 raised against the lake
solver's dry tail, arriving in the same place by a different route, and no sink
fixes it.

`ET_max` should ultimately be the Penman estimate `carve_verdict.py` already
computes and already shares with the lake solver and the carve criterion, so
that one evaporation rule governs lakes, carving and the water table rather than
three. Using a constant here is a deliberate simplification for the test in
section 4, whose whole purpose is to isolate recharge as the source of variation.

---

# The result: one prediction held, one failed, and the failure is the finding

Run 2026-08-20 on the GW-3 Earth harness, everything else held.

**The inert path is bit-identical.** With `et_max_m_s` unset, `depth_m` is
`np.array_equal` to the pre-change result and closure's ET door reads exactly
0.0. The term cannot affect a run that does not ask for it.

`ET_max` is a spatially constant 1,500 mm/yr, deliberately, so that the ONLY
source of spatial variation is recharge. Closure holds throughout, 1e-15 to
1e-13, with the third door carrying real water.

| lambda | at surface | median depth | 95th pct | residual mean | residual sd |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 m | 0.0% | 2.59 m | 3.5 m | -12.61 m | 40.05 m |
| 1.0 m | 0.0% | 5.18 m | 6.9 m | -10.49 m | 40.05 m |
| 2.0 m | 0.0% | 10.36 m | 13.5 m | **-6.23 m** | 40.06 m |

against observed median 5.89 m, 95th percentile 50.9 m, standard deviation
40.05 m, and thresholds of |mean| <= 8.92 and sd <= 24.56.

## Against the prediction

**1. Pinning goes: HELD, and by more than predicted.** 95.5% to 0.0%, against a
predicted "under 20%". The sink does exactly the job it was added for.

**2. The model acquires skill: FAILED.** The residual standard deviation is
40.05 m, which is the observed standard deviation, unchanged from the 40.06 m
the surface-pinned model scored. Pearson correlation between model and
observation is **+0.031**, so the model explains **0.1% of the variance**.

**3. The thresholds: the mean one PASSES at lambda = 2 m**, at -6.23 m against
8.92. The standard deviation misses at every lambda, and it is the one that
matters.

**So the headline is not the one predicted.** The prediction expected the miss to
come from a missing deep tail. It does not. The whole distribution is compressed:

| | 5th pct | median | 95th pct | ratio |
| --- | ---: | ---: | ---: | ---: |
| model | 2.91 m | 4.31 m | 5.56 m | **1.9** |
| observed | 0.90 m | 5.89 m | 50.88 m | **56.5** |

## Why, and it is not a bug

The mechanism works exactly as designed, and the correlations prove it:

- model against recharge: Spearman **-0.977**
- observed against recharge: Spearman **-0.156**

`d = lambda ln(ET_max / R)` makes the model a near-deterministic function of
recharge, which is what the equation says it should be. **The real water table
is not a function of recharge.** Australian recharge spans 1 to 869 mm/yr, but a
logarithm turns a factor of 900 into a factor of 13, and the recharge actually
present is concentrated enough that the model spans less than a factor of 2.

Observed depth tracks cell elevation slightly better than recharge, at Spearman
+0.222, and Fan et al. (2013)'s own finding is that terrain dominates water
table depth at local scales -- "the well-articulated gradient is the topography
from valley to ridge, spanning decameters to kilometers". **At 15.19 km a cell
holds no valley and no ridge, so the model has no mechanism for the variation
the observations are made of.** That is GW-6, and this is the second and sharper
measurement of it: the first said the model was constant at zero, this says it
is still nearly constant when the physics that removes the pinning is correct.

## What stands

**The sink stays.** It is right physics from a real source, it removes a state
the Earth data says is wrong, it costs nothing when off, and closure now has
three doors and still balances. Calling it necessary was correct.

**Calling it sufficient would have been wrong, and the prediction said so** --
though for the wrong reason, and that is worth more than a prediction that
happened to land. The failure is not a missing tail. It is that a 15 km cell
cannot carry the terrain signal that sets a water table, so no sink, no
thickness and no conductivity recovers the variance.

`ET_max` as a constant is the obvious next refinement and it will not change
this. Using the Penman field would give ET_max the spatial variation of
evaporative demand, which on this continent is smooth and broad, while the
observed variation is at the scale of valleys. It would move the mean and leave
the standard deviation where it is.

---

# GW-16: what the sink does to the carve result

Run 2026-08-20 once the solver converged. The sink-off case reproduces
1798 / 1881 / 83 / 83 / 0 exactly, so everything below is the sink and nothing
else.

| lambda | median shift | carving | flipped | to carve | to hold |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 m | 0.00105 | 1,916 | 118 | 118 | 0 |
| 1.0 m | 0.00117 | 1,925 | 127 | 127 | 0 |
| 2.0 m | 0.00140 | 1,930 | 132 | 132 | 0 |
| none, GW-4 | 0.00073 | 1,881 | 83 | 83 | 0 |

against a surface-only baseline of 1,798.

**The count goes UP, and I expected it to go down.** The reasoning was that
seepage is what handed a zero-runoff basin its supply, so a sink that removes
water before it reaches the surface should hand over less. That is true of
seepage and false of the answer, because it ignores what pinning was doing to
the flow.

A water table pinned at the surface has its head fixed at the terrain. It cannot
develop a gradient of its own, so there is almost nothing for lateral exchange to
respond to and seepage does all the work locally. Freeing it is what lets
groundwater actually flow: the free set goes from 399,882 cells to 624,228, and
the median absolute exchange rises from 0.00073 of a basin's recharge to 0.0014.
**The sink increases exchange by letting the water table detach and move**, and
that outweighs the water it removes.

So GW-4's 83 was an UNDERESTIMATE rather than the upper bound GW-16 was opened
to test. The direction is unchanged and now holds across the whole bracket:
every flip is toward carving and none toward holding, at every lambda.

## What the earlier 69 flips to hold were

An intermediate run reported 162 to 164 flips with 69 of them toward holding.
Those were entirely an artefact of `measure_carve_effect` comparing recharge
against SEEPAGE, which with a sink differ by the evaporation. Attributing that
difference to exchange made basins look like they were losing water underground
when they were losing it to the air. Once the comparison isolates the exchange,
`give + qg` against `give`, the holds vanish and the count settles.

That is worth keeping because closure never saw it: the global books balanced at
9.2e-13 throughout. The water was accounted for and attributed to the wrong door,
and only a derived quantity being physically absurd -- a redistribution moving
90% of every basin's water -- gave it away.

## Still bounded by the same thing

GW-3 says this water table has no skill against real observations, and nothing
here changes that. What the bracket shows is that the carve direction survives
both the sink and the lambda range, not that the count is right. It remains a
count of basins crossing a threshold in a model whose depth field is known not
to reproduce the one quantity it can be tested against.

---

# Why the depth range collapses, and what each downstream consumer needs

## The model has no local baselevel

`solve()` has exactly one fixed-head boundary: `is_ocean`, at sea level. That is
the whole of it, and it is the range problem.

With the sink dominant the equilibrium depth is `d = lambda ln(ET_max / R)`. A
logarithm turns Australia's factor of 900 in recharge into a factor of 13 in
depth, so at lambda of order a metre the entire field lives between 2.7 and
6.9 m. Reaching 100 m from this mechanism would need `ET_max/R = e^100`. **The
sink structurally cannot produce a deep water table.**

Real water tables are deep for a different reason: the land stands high above a
NEARBY drainage baselevel and the water drains away sideways. That needs a
baselevel within reach, and ours is the sea, up to 1,000 km off, which no
plausible transmissivity reaches. So the lateral term never competes and the
local balance wins everywhere.

It is also why mesh distance to the sea is the best single predictor of observed
depth at +0.268: it is a crude proxy for height above baselevel, and the only
part of that mechanism this geometry still carries.

**The baselevels exist already.** `surface_water.nc` carries `lake` and
`discharge_m3_s` per region and `regions.nc` carries the `receiver` network.
A river or lake IS a fixed head at the water surface. Using them collapses the
drainage length from 1,000 km to something like 15 to 50 km, where the arithmetic
that condemned the coastal boundary becomes survivable: holding 20 m of head over
15 km at 8 mm/yr needs `T` near 1.4e-3, which at Gleeson's conductivity is about
2 km of aquifer -- thick, but a real sedimentary basin rather than the 600 km the
coastal version demanded.

## MIN-6 and SURF-7 fail differently, and only one is closed

Measured on the depth-consistent Australian bores.

**MIN-6 fails on RANGE, and no amount of skill fixes it.** Supergene enrichment
keys on the water table being deep -- Reich and Vasconcelos put leached caps at
several hundred metres "particularly when the water table was deep enough". The
model's 95th percentile over land is 6.86 m, and it puts 1.1% of land past 30 m
where the observations put 8.5%. A perfectly skilful model confined to 2-7 m
still could not drive that rule.

**SURF-7 fails on CLASSIFICATION, and is not closed.** It wants a mask, "at or
near a water table, or at zones of groundwater outflow", so the test is
discrimination rather than a value:

| threshold | observed base rate | AUC |
| --- | ---: | ---: |
| <= 1 m | 4.3% | 0.573 |
| <= 2 m | 17.2% | 0.546 |
| <= 5 m | 45.9% | 0.520 |
| <= 10 m | 68.0% | 0.540 |

Barely above the 0.5 of no discrimination, and too weak to place a surface class
on -- but strongest at the tightest threshold, which is where SURF-7 operates.

**And the route for it is a saturated-area FRACTION, not a cell flag.** Given a
cell-mean water table depth and the distribution of elevations WITHIN the cell,
the fraction lying below the table is the discharge area. That is the standard
topographic-index construction, and the sub-grid hypsometry it needs is the same
quantity GRAV-6 says the downscaling machinery has to persist anyway. It would
also replace a genuinely unstable approximation: "at surface" is currently
all-or-nothing per 294 km2 cell, which is why that fraction swings between 63%
and 14.5% on a formulation change.

The fraction inherits the cell-mean depth's error, so it is necessary and not
sufficient, the same shape as the sink. It turns skill in the mean into a usable
mask; it does not manufacture skill.

---

# GW-17 and GW-18: the range is fixed, the pattern is not

Run 2026-08-20. River cells are derived rather than assumed: a priority flood on
the mesh, steepest descent on the filled surface, and flow accumulation, which is
the same construction `build_hydrography.py` uses for Vesper. A cell above the
accumulation threshold takes a fixed head at its own surface.

| configuration | 95th pct depth | r^2 | share of the 0.857 ceiling |
| --- | ---: | ---: | ---: |
| no rivers, D = 100 m | 6.9 m | 0.0114 | 1.3% |
| rivers >= 5e3 km2, D = 100 m | 6.8 m | 0.0190 | 2.2% |
| rivers >= 5e3 km2, D = 2000 m | **52.4 m** | 0.0173 | 2.0% |
| rivers >= 1e4 km2, D = 2000 m | **59.6 m** | 0.0166 | 1.9% |

against an observed 95th percentile of 42 m on depth-consistent bores.

**The range problem is solved.** With local baselevels and a basin-scale
thickness the model reaches 52 to 60 m at the 95th percentile where it reached
6.9 m before. That was MIN-6's blocker, and it was a blocker about what values
the field can take rather than about skill.

**The skill problem is untouched**, at about 2% of the ceiling. The model now
produces deep water tables and puts them in the wrong places.

## And the mechanism I blamed is not the one the observations follow

The subdued-replica argument says depth should track height above the nearest
drainage baselevel. That is now computable, by following the drainage path to
the first river cell. It does not:

| predictor of observed cell-mean depth | Spearman |
| --- | ---: |
| height above the nearest river | **+0.064** |
| cell elevation | +0.218 |
| log flow accumulation | -0.095 |
| mesh distance to the sea | +0.268 |

Height above the nearest river is the WEAKEST of them, on a mesh where that
height has a median of 56.9 m and a 95th percentile of 320 m, so it is not for
want of range in the predictor.

**So the diagnosis in the section above was half right and the half it got wrong
is the important half.** Missing local baselevels did explain the range
collapse, and adding them does fix it. They do not explain the missing skill,
because at 15.19 km the valley-to-divide structure that mechanism lives on is at
or below the cell size: real drainage spacing is 5 to 20 km, so an interfluve is
one or two cells wide and the model has nowhere to put a subdued replica.

What predicts observed depth at this scale is weak and diffuse -- no single field
exceeds 0.27, and a flexible fit over all of them reaches R^2 = 0.1435 on
depth-consistent bores. There is resolvable signal, it is a quarter of the
variance rather than the 0.86 the cell decomposition allows, and it is not
carried by one mechanism that a physical model could be given.
