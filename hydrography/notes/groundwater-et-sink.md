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
