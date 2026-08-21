# GW-3: the criterion, declared before the score

Worldbuilding. Vesper is an invented planet; this note is about validating the
groundwater solver written for it against measurements made on Earth.

**Written 2026-08-20, BEFORE any model output was compared to any observation.**
That ordering is the whole point of the document: `docs/src/practice/conventions.md`
requires a threshold fixed before results are seen, and a criterion chosen after
the run it judges is not a criterion.

## What is being tested

The same solver, unchanged, on Earth. It reads seven fields from an export --
radius, x, y, z, region count, cell area, surface class -- and `Geometry` derives
its own Delaunay from the generators, so nothing is reimplemented for this run.

The mesh is Orogen's own `generateFibonacciSphere`, ported faithfully with its
RNG and the project's own `jitter = 0.75`, sized to Vesper's exact cell area of
293.80 km2, giving 1,736,112 regions at Earth's radius. It reproduces Vesper's
operator error: relative RMS against the analytic Legendre eigenvalue of 0.0006,
0.1062, 0.1186 and 0.1235 at l = 1 to 4, against Vesper's 0.0217, 0.1082, 0.1202
and 0.1243. **The discretisation is therefore not a free parameter in what
follows**, and a mismatch cannot be blamed on the mesh.

Australia is modelled as an island with an ocean boundary, which is what it is.
Land comes from ETOPO 2022 at 60 arc-seconds and totals 7,748,094 km2 against
Australia's 7,688,000, the 0.8% excess being what a 15 km discretisation of a
coastline gives.

## The benchmark is Fan's own performance, not a number I chose

Fan et al. (2013) Table S1 reports her residual statistics against the same
compilation this is scored on. For **Australia and Asia**, over 51,126 cells
carrying observations:

| forcing | mean residual | standard deviation |
| --- | ---: | ---: |
| CLM recharge | -5.08 m | 20.85 m |
| Doll-Fiedler recharge | -8.92 m | 24.56 m |

Her global totals are -1.62 m and 17.91 m.

## The criterion

Residual is defined as `model_depth - observed_depth` in metres, positive where
the model puts the water table too deep, evaluated **at observation locations
only** and never model-everywhere against observations-everywhere.

**PASS requires both:**

1. standard deviation of the residual **<= 24.56 m**
2. absolute mean residual **<= 8.92 m**

Those are Fan's own Australia-and-Asia figures under the worse of her two
recharge forcings. The criterion is therefore "no worse than the published state
of the art in the same region", which is falsifiable and is not a number chosen
to be reachable.

**The comparison is generous to Fan and strict here, deliberately.** Her figures
are at 30 arc-seconds, about 1 km; this runs at 15.19 km, which is 15 times
coarser in the direction that matters, because Fan's own finding is that terrain
signal dominates water table depth at local scales. Matching her at this
resolution would be surprising. Failing to is the expected outcome and would
still be informative.

## What a failure would and would not mean

`T = K D` is NOT Fan's formulation. She uses an exponential decay of permeability
with depth, which this project abandoned under GW-9 because a metre-scale
e-folding length has no cell-mean value at 15.19 km. **So a miss confounds
formulation with implementation, and saying which is part of the result rather
than an excuse for it.** Three things separate them:

- The reduction, closure, uniqueness and trace identities already pass, so the
  implementation solves the equations it states.
- The operator matches Vesper's to within 2% at l = 2 to 4, so it is not the
  discretisation.
- A residual that correlates with terrain slope or with recharge points at the
  formulation; one that does not points at the inputs.

Fan's own residuals correlate with slope at -0.49 globally, so she has this
problem too, and that number is the thing to compare against rather than zero.

## Recorded before the run

Recharge is Berghuijs et al. (2022), 30 arc-second, in mm/yr. Permeability is
GLHYMPS `Permeability_no_permafrost`, whose values over Australian land come out
at the 5th, 50th and 95th percentiles as -15.20, -13.00 and -11.80 in log10 k --
exactly Gleeson (2011)'s siliciclastic, unconsolidated and carbonate class means,
which is the same table `hydrography/config/groundwater.yaml` uses for Vesper. The
two worlds therefore share their permeability parameterisation exactly.

Aquifer thickness is the same `D` used on Vesper. Conductivity is `K = k rho g /
mu` at Earth's gravity, which is 0.766 times Vesper's for the same rock, and
getting that wrong would make this a different model rather than the same one.

Observations are the 75,321 Australian sites in
`hydrography/data/earth_validation/aus_wtd_sites.csv`, assembled from the
Australian Groundwater Explorer with the datum reconciliation recorded in
`docs/src/reference/external-data.md`.

---

# The result: MISS on both thresholds

Scored 2026-08-20 at 70,119 of the 75,321 Australian sites, the remainder
falling on cells the solver excludes as non-conductive or ocean.

| | observed | model |
| --- | ---: | ---: |
| median | 5.89 m | 0.00 m |
| mean | 14.86 m | 0.12 m |
| standard deviation | 40.05 m | 2.28 m |

| residual, model minus observed | value | threshold | |
| --- | ---: | ---: | --- |
| mean | **-14.74 m** | <= 8.92 | MISS |
| standard deviation | **40.06 m** | <= 24.56 | MISS |

**The model has no skill, not merely a bias.** Its residual standard deviation
of 40.06 m is the observed standard deviation of 40.05 m, because the model
output is very nearly a constant: it puts the water table at the surface at
99.2% of scored sites, where 5.9% of observations are shallower than one metre.
A prediction that explains none of the variance is what a constant does.

## It is not the discretisation, the permeability or the implementation

Each was excluded before the run rather than after the answer:

- The Earth operator reproduces Vesper's relative RMS error against the analytic
  eigenvalue to within 2% at l = 2 to 4.
- GLHYMPS over Australian land returns log10 k of -15.20, -13.00 and -11.80 at
  the 5th, 50th and 95th percentiles, which ARE Gleeson (2011)'s class means and
  are the same numbers `groundwater.yaml` gives Vesper.
- The reduction, closure, uniqueness and trace identities all pass, and this run
  converged in five passes to a residual of 2.2e-16 with zero leak.

## It is the formulation, and no aquifer thickness fixes it

Sweeping the one free parameter, with everything else held:

| thickness | median T | at surface | median depth | 95th pct depth |
| ---: | ---: | ---: | ---: | ---: |
| 100 m (declared) | 7.5e-05 | 95.5% | 0.00 m | 0.00 m |
| 1 km | 7.5e-04 | 78.7% | 0.00 m | 39 m |
| 10 km | 7.5e-03 | 38.8% | 7.96 m | 143 m |
| 100 km | 7.5e-02 | 12.0% | 69.10 m | 322 m |
| 1,000 km | 7.5e-01 | 2.0% | 178.45 m | 517 m |

against an observed median of 5.89 m and 95th percentile near 49 m.

The observed MEDIAN is reproduced at a thickness of about 10 km, which is a
hundred times the declared value and not a defensible aquifer. **And even there
the DISTRIBUTION is wrong**: the 95th percentile comes out at 143 m against 49 m
observed. No thickness matches both ends, so this is not a parameter that was set
badly. It is the wrong shape of answer.

The reason is a one-dimensional estimate that needs no model. For a confined
aquifer of transmissivity `T` draining recharge `R` over a distance `L`, the head
above the outlet goes as `R L^2 / 2T`. Australian recharge is 8 mm/yr, and the
interior sits about 1,000 km from any coast, so holding the water table 250 m
below a 300 m land surface needs `T` near 0.5 m2/s, which at Gleeson's
conductivity is an aquifer some 600 km thick. **Recharge in the interior of a dry
continent cannot reach the sea.** In the model it therefore has nowhere to go but
up, and the water table rises until it seeps.

## The missing physics is an internal sink, and it is already on the books

`groundwater-scoping.md` section 2 lists four channels through which a water
table changes an answer here. The second is "evaporation from a shallow water
table in discharge zones, which is a sink the single bucket cannot represent",
and it is GW-5, filed as a one-signed omission. **This result reclassifies it
from a refinement to a precondition**: without a sink between the recharge cell
and the coast, steady state has only one place to put the water.

Real arid groundwater loses recharge to evapotranspiration from shallow tables
and to internal discharge over tens of kilometres, not the thousand it would
take to reach the sea. Australian groundwater also carries residence times of
1e4 to 1e6 years and is not in equilibrium with modern recharge at all, which is
the same objection HYD-4 raised against the lake solver's dry tail and the same
one `no-time-axis` prevents this project from answering.

## What this says about the Vesper result, which is the point of running it

`groundwater_report.json` records `at_surface_fraction_land` of **0.630** and a
median land water table depth of **0.00 m**. Vesper shows the same pathology,
milder because its endorheic basins put internal sinks far closer than a coast,
but present: most of its land also has the water table at the surface.

**So GW-4's 83 carve flips rest on seepage from a surface-pinned water table,
and this test says a surface-pinned water table is wrong.** What survives is the
direction, which never depended on the depth: a basin with no surface runoff can
acquire water from outside its surface catchment, and the mechanism is a
zero-crossing rather than a magnitude. What does not survive is the count. It
should be read as an upper bound until an internal sink exists, because
overstating seepage overstates the supply that lifts a basin over the threshold.

---

# Amendment: the observation filter, declared before the rescore

Written 2026-08-20, before the model was scored against a filtered set.

## The set was measuring two different things

`aus_wtd_sites.csv` now carries `bore_depth_m`, and it should have from the
start. A bore screened well below the water table measures a POTENTIOMETRIC HEAD
in a confined aquifer, which is a different surface. This set runs to 129 m at
the 95th percentile, and Australia contains the Great Artesian Basin.

Decomposing the observed variance within and between 15.19 km cells:

| kept | within-cell | ceiling on R^2 at this cell size |
| --- | ---: | ---: |
| all bores | 32.1% | 0.679 |
| depth <= 30 m | 21.7% | 0.783 |
| depth <= 10 m | **14.3%** | **0.857** |

And what resolvable geography can predict rises with it. Held-out R^2 on
cell-mean depth from elevation, local relief, mesh distance to the sea, recharge
and permeability:

| kept | linear | 12x12 binned on distance x elevation |
| --- | ---: | ---: |
| all bores | +0.071 | +0.122 |
| depth <= 30 m | +0.153 | +0.130 |
| depth <= 10 m | **+0.280** | +0.190 |

**So the original GW-3 verdict understated the ceiling and partly scored the
model against the wrong quantity.** It does not overturn it -- the model scored
r^2 = 0.001 against a ceiling of 0.679, and 0.001 against 0.857 is no better --
but the target was wrong and the ceiling was wrong, and both matter for what
comes next.

## The filter

**Primary: `bore_depth_m <= wtd_m + 20`.** A bore whose bottom is within 20 m of
the water level it reports cannot be open to a deeper aquifer, so it is
measuring the water table. This is a statement about construction against
measurement, not about the model, and it is chosen on that ground rather than on
which cut flatters the score.

Reported alongside, as sensitivities: absolute cuts at 30 m and 10 m, and the
unfiltered set.

## The thresholds do not move, and the comparison changes meaning

PASS stays at residual standard deviation <= 24.56 m and |mean| <= 8.92 m, from
Fan's own Australia-and-Asia figures.

**But Fan did not filter, and could not: her four columns are latitude,
longitude, elevation and depth, with no construction data at all.** Her scatter
therefore contains the same mixed-surface noise this filter removes. So a
filtered score is NOT like-for-like against her numbers -- some of any
improvement is removed observation noise rather than a better model.

Both are reported for that reason. The unfiltered score is the one comparable to
Fan; the filtered score is the better test of the model. Quoting the filtered
number against her threshold and calling it a pass would be scoring an easier
exam.

---

# The rescore, and why its PASS is not a pass

Run 2026-08-20 against the declared filters. **The verdict does not change: the
model has no skill.** What changed is that the thresholds stopped being able to
say so, and the reason is worth more than the numbers.

## A data defect the first score ran straight past

**3,917 sites, 5.2%, record a water level below the bottom of their own bore.**
A median of 23.9 m of water in a 2.6 m bore; at worst 1,036 m in a bore 90 m
deep. Those are keying and unit errors, not hydrology, and they dominated the
spread: excluding them takes the observed standard deviation from 38.8 m to
18.0 m and the 95th percentile from 49.4 m to 40.5 m.

The declared filter did not catch them, and the sensitivity cuts selected FOR
them. `bore_depth <= wtd + 20` is satisfied whenever the bore is far SHALLOWER
than the water level, so the impossible records passed it; and "bore depth under
10 m" is exactly the set where a spurious deep water level is most likely, which
is why the observed standard deviation ROSE to 60 m under the strictest cut. A
filter bounded on one side only is not a filter.

`aus_wtd_sites.csv` now carries `depth_consistent`. 57,093 of 75,321 sites pass
it.

## The scores

| set | n | observed sd | residual mean | residual sd | r^2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| unfiltered, no sink | 70,119 | 40.05 m | -14.74 m | 40.06 m | 0.0007 |
| unfiltered, with sink | 70,119 | 40.05 m | -10.49 m | 40.05 m | 0.0010 |
| consistent, with sink | 53,410 | 18.36 m | -7.14 m | 18.26 m | 0.0114 |
| consistent + shallow, with sink | 37,501 | 11.20 m | -3.79 m | **11.07 m** | 0.0315 |

The last two rows meet both of Fan's thresholds. **They are not passes.**

**The residual standard deviation equals the observed standard deviation at
every level.** 40.06 against 40.05, 18.26 against 18.36, 11.07 against 11.20.
That is the signature of a constant, and it is invariant to how the set is
filtered. A model that predicts the mean everywhere scores 11.20 m on the last
row and clears a 24.56 m bar comfortably.

**Fan's thresholds only mean something against a set as noisy as hers.** She did
not filter and could not, her four columns carrying no construction data at all,
so 24.56 m is the scatter of a real model against observations that include this
same 5% of impossible records. Transplanting that number onto a set whose own
spread is 11 m turns it from a test into a formality. The declaration written
before this run said quoting a filtered number against her threshold would be
scoring an easier exam; that is exactly what the last two rows are.

**So the criterion for a filtered set has to be a SKILL criterion, not an
absolute scatter.** Against the cleaned set the ceiling is R^2 near 0.86 and the
model reaches 0.03.

## What did improve, and it is real

Cleaning the observations and adding the sink together take r^2 from 0.0007 to
0.0315, a factor of 45. That is a real gain and it is still three percent of the
variance. It says the earlier "no skill whatsoever" was partly the observations'
fault and mostly the model's, which is a more useful statement than either alone.

---

# GW-19: the criterion is replaced, and this one cannot be passed by a constant

The scatter bars above are withdrawn as a PASS test. They stay recorded as what
Fan achieved, which is still the useful comparison, but they cannot decide
anything here for a reason the rescore demonstrated: **residual standard
deviation equalled observed standard deviation at every filter level** -- 40.06
against 40.05, 18.26 against 18.36, 11.07 against 11.20 -- so a model predicting
the mean everywhere scores the observations' own spread and clears a 24.56 m bar
on any set quiet enough. The bar moves with the data; the failure does not.

## The replacement

Scored on depth-consistent bores, `0 < wtd_m <= bore_depth_m`, at observation
locations.

    R^2 = 1 - sum (model - observed)^2 / sum (observed - mean observed)^2

**PASS requires R^2 > 0.07.**

That is not a round number chosen for comfort. It is what the single best
resolvable predictor achieves on its own: mesh distance to the sea correlates
with observed cell-mean depth at Spearman +0.268, and a naive monotone
regression on it explains about 7% of the variance. **A physical model that
cannot beat "how far is it to the coast" has not earned its machinery**, and the
bar is therefore set by what the geography gives away for free rather than by
what seems achievable.

Two reference points recorded beside it, neither of them the bar:

- **0.857** is the ceiling. Between-cell variance is 85.7% of the total on
  depth-consistent bores, so no model at 15.19 km can exceed it however good.
- **0.1435** is what a flexible statistical fit over every resolvable field --
  elevation, local relief, distance to the sea, recharge, permeability --
  reaches on held-out cells. That is the practical target: a physical model
  should approach what a regression on its own inputs can do.

## Where the model stands against it

| configuration | R^2 | verdict |
| --- | ---: | --- |
| no sink, D = 100 m, no rivers | 0.0042 | MISS |
| sink, D = 100 m, no rivers | 0.0114 | MISS |
| sink, rivers, D = 2000 m | 0.0173 | MISS |

The best configuration reaches a quarter of the bar and 6% of the practical
target. **The verdict is unchanged and is now stated in a measure that filtering
cannot flatter.**

---

## The resolution confound, separated: it is the formulation

This note declared, before any scoring, that "a miss confounds formulation with
implementation, and saying which is part of the result rather than an excuse for
it", and that running at 15.19 km against Fan's roughly 1 km was "15 times
coarser in the direction that matters". That was a promissory note. It is now
paid, by running the same case at three cell sizes with
`hydrography/scripts/earth_calibration.py`. Measured 2026-08-20, on
depth-consistent bores, with GW-15's sink and GW-17's baselevels ON.

| mean edge | regions | bores | cells | ceiling | model sd | Pearson | R2 on cell means |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 15.19 km | 1,736,112 | 53,410 | 3,276 | 0.6335 | 1.91 m | +0.0703 | -0.365 |
| 10.74 km | 3,473,013 | 54,921 | 4,682 | 0.6894 | 3.35 m | +0.0665 | -0.357 |
| 7.60 km | 6,935,660 | 55,009 | 6,548 | 0.7408 | 4.33 m | +0.0711 | -0.344 |

Halving the cell edge does two things and fails to do the third. The CEILING
rises from 0.6335 to 0.7408, so more of the observed variance genuinely becomes
resolvable. The model's own spread more than doubles, 1.91 m to 4.33 m, so it is
producing more structure. And the correlation does not move at all: +0.0703,
+0.0665, +0.0711, which is a Pearson r2 near 0.005 at every resolution. The
model gains variance without gaining agreement. The extra structure is in the
wrong places.

The dependence structure says why, and it is the clearest result here:

| mean edge | rho(model, recharge) | rho(obs, recharge) | rho(model, elevation) | rho(obs, elevation) |
| --- | --- | --- | --- | --- |
| 15.19 km | -0.902 | -0.272 | -0.047 | +0.125 |
| 10.74 km | -0.931 | -0.279 | -0.020 | +0.161 |
| 7.60 km | -0.923 | -0.268 | -0.003 | +0.199 |

The model is very nearly a monotone function of recharge at every cell size,
Spearman -0.90 to -0.93, and the real water table is not, at -0.27. Against
elevation the two move in opposite directions as the mesh refines: the
observations' terrain signal STRENGTHENS, +0.125 to +0.161 to +0.199, exactly as
Fan's finding that terrain dominates locally predicts once cells are small
enough to show it, while the model's decays to nothing, -0.047 to -0.003.

So the resolvable terrain signal is real, it is increasingly available, and the
model is increasingly not using it. Refining the mesh under this formulation
produces a finer recharge map.

**The verdict is formulation, not resolution.** With `T = K D` at constant
thickness and a strong evapotranspiration sink, the steady state is set locally,
by the balance between recharge and what the sink removes at that depth, and the
lateral term is too weak for terrain to organise the water table across a
catchment. Fan's exponential decay is not a refinement of that, it is the
mechanism that couples the two: transmissivity falls with depth, so a deeper
table conducts less, which is what makes the surface follow topography with a
damped amplitude rather than following recharge.

GW-9 abandoned the exponential decay on the argument that a metre-scale
e-folding length has no cell-mean value at 15.19 km. That argument was never
measured, and this result puts it back on the table: whatever the e-folding
length means at a cell mean, the formulation replacing it does not reproduce the
dependence structure of the thing it models, at any of the three cell sizes
tried. The next experiment is the cheap one, the same mesh with a depth-decaying
transmissivity, and it needs no new data.

Two things this does NOT show. It does not show the implementation is wrong: the
reduction, closure, uniqueness and trace identities pass, the operator matches
Vesper's to 2%, and the solve converges with nothing pinned. And it does not
show the depth range is wrong: the model's median depth is 5.12, 5.16 and 5.18 m
across the three, stable and physically ordinary. What is wrong is WHERE the
water is, not how much of it there is or how the equations are solved.

### A stale result, caught by rebuilding

The harness that produced the original figures was written inline and lost
(GW-20). Recovering it from the session transcript showed it predates GW-15's
sink and GW-17's baselevels, and run as recovered it pins 93.8% of cells at the
surface with a median depth of 0.00 m, which is the pathology
`config/pipeline.yaml`'s gate names in as many words. The archived score was
therefore measured on a configuration this component no longer uses. Restoring
the current physics moves the baseline from R2 -0.579 to -0.365 and Pearson
+0.026 to +0.070, which changes no verdict but is the honest number. A result
whose harness is not kept is a result that quietly ages out of agreement with
the model it judged.

## What GW-21 found, including two corrections to the section above

The section above blamed the formulation and named Fan's depth-decaying
transmissivity as the fix. Checking the record first was the right move and it
withdrew that: `groundwater.py:transmissivity` already carries a quantitative
argument against it, `exp(h/f)` being convex so a cell mean carries
`exp(sigma^2/2f^2)`, which at 100 m of sub-grid relief against Fan's 0.95 m `f`
is `exp(5000)`, and `notes/water-table-convergence.md` records that Picard on it
LIMIT-CYCLES because `T` moves by a factor of e per e-folding length. Two
solvers died there. The experiment proposed was the one that had already failed.

What was run instead, on the cached 15.19 km Earth mesh:

**The sink does set the depth locally, exactly as its closed form says.** With
`E(h) = et_max exp(-(z-h)/lambda)` the balance is
`d = lambda ln(et_max A / supply)`, which names no neighbour. Median depth is
linear in `lambda`: 5.12, 10.23, 25.51, 50.86, 126.35 m at 1, 2, 5, 10 and 25 m.
`rho(model, recharge)` is pinned at -0.90 for every one of them, and `et_max`
barely registers because it enters logarithmically.

**Transmissivity magnitude, not its depth dependence, controls whether terrain
appears.** Raising `D` from 100 m to 20 km moves `rho(model, elevation)` from
-0.047 to +0.466 and `rho(model, recharge)` from -0.902 to -0.245, which brackets
the observed -0.272 and +0.125. The dependence structure is therefore reachable
with the formulation already in hand. Without any sink the model is not
recharge-dominated at all, -0.106, and carries +0.082 of terrain signal; it is
simply pinned at the surface on 91.4% of cells, because at `T = K * 100 m` the
aquifer cannot carry its recharge to a baselevel.

**And none of it improves agreement.** Across `D` over a 200-fold range, `lambda`
over 25-fold, and the sink on or off, Pearson on cell means stays between +0.058
and +0.101. Matching the dependence structure does not produce a match.

### The correction that matters: rank against magnitude

This note has been quoting Pearson. On rank the model is not skill-free:
Spearman is +0.276 at the current setting and +0.291 at its best. Reporting
"no skill" was too strong and is withdrawn.

What that rank skill IS, though, settles the question:

| field or model | Pearson | Spearman |
| --- | --- | --- |
| recharge alone | -0.118 | -0.272 |
| elevation alone | +0.072 | +0.125 |
| log10 K alone | +0.042 | +0.052 |
| the model, current settings | +0.070 | +0.276 |
| the model, best of everything tried | +0.101 | +0.291 |

The model's rank agreement is its recharge field's rank agreement. Solving a
groundwater equation on top of that input recovers what the input already had
and adds no more than 0.02. Meanwhile R^2 on cell means is negative everywhere
and grows worse as `D` rises, -0.365 to -0.749 to -2.318 to -7.025, because the
model's spread overshoots: 45 m and 75 m of standard deviation against the
observations' 28.5.

So the model orders cells about as well as its recharge forcing does, and places
the magnitudes badly, and the extremes worst of all, which is why Pearson sits
so far below Spearman.

### What this means, and what it does not

It is not a resolution gap: the sweep above shows the correlation flat from
15.19 to 7.60 km. It is not the exponential transmissivity: that is ruled out
twice over and was never run. And it is not obviously a defect in the flow
solver, whose identities pass and whose dependence structure can be tuned onto
the observed values without agreement following.

The reading this leaves is that the INPUTS do not determine the answer. Recharge,
GLHYMPS permeability at 15 km and a cell-mean elevation carry, between them,
about the rank skill the model achieves, and the independent statistical result
recorded above -- 0.1435 held-out R^2 from a flexible fit over every resolvable
field -- says the room above that is modest against a within-cell ceiling of
0.63. Those two numbers are different metrics and should not be equated, but
they point the same way.

That is a harder finding than a formulation bug, because no reformulation of a
flow equation adds information its forcing does not contain. It also means the
honest use of this field is unchanged and is what `hydrography/README.md`
already says: the per-basin exchange direction survives, and a per-cell depth
does not mean anything.

## Which surface is the depth measured from, and why both answers are wrong

A model depth and an observed depth are only comparable if they hang from the
same surface, and they do not. The model's is below the DEM's CELL MEAN. The
bore's is below the ground at the bore, and bores sit in valleys -- the bias Fan
states herself and `docs/src/reference/external-data.md` records, that wells are
where people are, which is valleys and oases. Measured on the 48,552 Australian
bores that fall on conductive cells, the siting is real: the DEM at a bore sits
a mean of 8.38 m and a median of 1.91 m BELOW its cell mean.

Correcting it is not as simple as using the bore's own elevation, and the first
attempt is worth recording because the failure is instructive. Scoring
`bore_elevation - head` gives R2 -18.0 and a model spread of 107 m against the
observations' 28. The reason is a datum and resolution mismatch, not physics:
the surveyed bore elevation and the ETOPO value AT THE SAME POINT differ with a
standard deviation of 51.97 m, median -0.28 m, so there is no systematic offset
to remove and only scatter, which lands whole on the residual. **Do not mix an
observation's own elevation with a head anchored to a different DEM.**

Taking both from the DEM removes that, and the result is the interesting one:

| depth measured below | Pearson | Spearman | R2 | model sd |
| --- | --- | --- | --- | --- |
| the cell mean | +0.0689 | +0.2734 | -0.354 | 3.96 m |
| the DEM at the bore | +0.1278 | +0.2004 | -3.122 | 44.36 m |
| the observations | | | | 28.19 m |

Pearson nearly doubles and everything else gets worse, because the model's
spread goes to 44 m against an observed 28. The two rows are two assumptions
about what happens INSIDE a cell and neither is true. Scoring below the cell
mean assumes the depth is constant across the cell, so the water table copies
the terrain exactly. Scoring below the DEM at the bore adds the full within-cell
surface variation, standard deviation 36.24 m, to a single flat head, which
assumes the water table is level across the cell and the depth carries all the
relief. A real water table is a subdued replica of the topography and sits
between the two.

That is the sub-grid drainage term of GW-23 seen from the scoring side rather
than the physics side, and it is now BRACKETED rather than argued: the correct
within-cell behaviour lies between a model spread of 3.96 m and one of 44.36 m,
and the observations say 28.19 m. Any parameterisation of sub-grid venting has
those three numbers to hit, which is a sharper target than "add a sink".

### GW-23 tested at its bounding limit, and a defect found on the way

The sub-grid drain has a free end: a Robin condition `Q = C max(0, h - z_valley)`
needs a solver term, but `C -> infinity` is just a fixed head at the cell's own
valley floor, which the existing `fixed_head_m` supplies. That BOUNDS what the
whole family can do, and it can be run on the cached mesh in a second.

22,284 of 26,372 Australian land cells have 10 m or more of within-cell relief
to vent into, at a median of 30 m. Draining them:

| configuration | Pearson | R2 | model sd |
| --- | --- | --- | --- |
| cell-mean surface, no drain | +0.0703 | -0.365 | 1.91 m |
| cell-mean surface, drain to DEM minimum | -0.0575 | -5.069 | 45.85 m |
| cell-mean surface, drain to mean - 2 sd | -0.0885 | -8.975 | 65.03 m |
| DEM-at-bore surface, no drain | **+0.1341** | -3.053 | 36.18 m |
| DEM-at-bore surface, drain to DEM minimum | +0.0492 | -3.025 | 43.32 m |
| the observations | | | 28.19 m |

Two things fall out and they point in opposite directions.

**The drain does not help, at either end of the scoring convention.** Below the
cell mean it inverts the correlation, +0.0703 to -0.0575, because it makes the
model deepest exactly where relief is greatest and the observations do not do
that. Below the DEM at the bore it still costs, +0.1341 to +0.0492. Median depth
goes to 25.5 m against an observed 12.9. `C -> infinity` over-drains, and since
the other end of the family is the undrained baseline, the family is bracketed
by two configurations without skill. A finite `C` would have to beat both ends,
which nothing here suggests.

**The scoring surface was a defect, and fixing it is the largest single
improvement in this whole investigation.** Simply hanging the depth from the DEM
at the bore rather than from the cell mean takes Pearson from +0.0703 to
+0.1341, nearly double, with no change to the model at all. That is a
measurement correction, not a modelling one, and it had been suppressing the
score from the start.

It is not skill. R2 stays at -3.05 because the spread goes to 36.18 m against
28.19 observed, which is the over-correction the section above brackets: hanging
every bore's depth from its own ground while the head stays flat across the cell
gives the depth all of the relief. The best configuration measured anywhere in
this note is Pearson +0.1341, against a statistical fit's 0.1435 and a ceiling of
0.6335.

One bug worth recording because it is a class rather than an instance. The first
version of the surface option changed the per-bore array and left the per-cell
aggregation reading `sol["depth"]`, so `--surface` ran, reported, and changed
nothing: identical R2 and Pearson to four decimal places while the model spread
moved from 3.96 m to 36.18. An option that appears to work and silently does
nothing is worse than one that fails. The aggregation now averages the model
over the same bores as the observations, which is what it should have done.

### GW-25: the surface correction is real, and it should not be the default

The two scorings are the ends of one family. Within a cell a water table is a
subdued replica of the terrain, `h(x) = h_cell + alpha (z(x) - z_cell)`, so the
predicted depth at a bore is `cell_depth + (1 - alpha) delta` with
`delta = z_bore - z_cell` taken from the DEM at both ends. `alpha = 1` is the
cell-mean scoring, depth constant across the cell; `alpha = 0` is the
bore-surface scoring, head flat across the cell. Sweeping it, on 3,276 cells:

| alpha | model sd | Pearson | Spearman | R2 |
| --- | --- | --- | --- | --- |
| 1.00, cell mean | 3.87 m | +0.0703 | +0.2760 | -0.365 |
| 0.90 | 5.99 m | +0.1399 | +0.2627 | -0.404 |
| 0.80 | 9.77 m | +0.1439 | +0.2441 | -0.494 |
| 0.36 | 28.56 m | +0.1365 | +0.2127 | -1.497 |
| 0.00, bore surface | 44.28 m | +0.1341 | +0.2036 | -3.053 |
| the observations | 28.49 m | | | |

The three measures disagree about where to stand, which is the first sign that
`alpha` is not being chosen by the physics. Pearson peaks near 0.8, Spearman and
R2 are both best at 1.0, and the one rule declared BEFORE looking at any
correlation -- match the observed spread of 28.49 m -- picks 0.36, which is
worse than either end on two measures out of three.

The attribution settles it. Pearson jumps at the first step off `alpha = 1` and
then flatlines, which is what happens when a term is being added that carries
its own signal rather than revealing the model's:

| predictor of observed cell-mean depth | Pearson | Spearman |
| --- | --- | --- |
| the model alone | +0.0703 | +0.2760 |
| `delta` alone, where the bore sits within its cell | +0.1290 | +0.1849 |
| model plus `delta` | +0.1341 | +0.2036 |

**`delta` out-predicts the model.** Where a bore sits inside its own cell
explains more of the observed depth, on Pearson, than the groundwater solve
does. So the near-doubling from the surface correction is mostly `delta`'s
predictive power folded into the model's score, and reporting it as the model's
number would be crediting the physics with a geometric covariate. A joint least
squares on both reaches R2 +0.0208 -- positive, unlike the -0.365 of a direct
comparison, because fitting an intercept and a scale is a different question --
against a flexible statistical fit's 0.1435 and a ceiling of 0.6335.

So the DEFAULT STAYS `cell-mean`, and `--surface dem-at-bore` stays available
and labelled. The defect was real and worth finding: the model's depth and the
bore's hang from surfaces a mean 8.4 m apart. But correcting it does not recover
model skill, it imports a covariate, and the honest report of the model is the
one that does not.

What the correction DOES establish is a sharper statement of GW-21's finding.
That row said the model adds nothing over its recharge forcing. This adds that
it is also out-predicted by pure geometry -- where the bore sits relative to its
cell mean, which involves no groundwater physics at all.

## The 0.28 was never measured, and it had been the target

This note quoted 0.28 four times as what "a flexible statistical fit over every
resolvable field reaches on held-out cells", `hydrography/README.md` quoted it,
`earth_calibration.py` PRINTED it beside every score, and GW-21 reasoned from
it. Reproducing the computation it came from, recovered from the session
transcript and re-run on the current artifacts, it does not exist. The original
run reported a linear held-out R2 of +0.0703 over six predictors and +0.1215 for
a nonparametric 12 by 12 binning on hops-to-sea against elevation. Nothing in it
produced 0.28.

Re-measured now on 2,939 multi-bore cells, target `log(1 + cell-mean depth)`,
the same 50/50 split and seed:

| fit | held-out R2 |
| --- | --- |
| linear, six predictors | +0.0711 |
| nonparametric 12 by 12 bins, hops by elevation | +0.0955 |
| gradient boosting, all seven predictors | +0.1411 |
| **gradient boosting, the physical model EXCLUDED** | **+0.1435** |
| gradient boosting, the physical model ALONE | -0.0467 |

So the honest number for what the resolvable fields support is **0.1435**, and
every citation has been corrected to it. A figure that no computation produced,
carried in four documents and printed beside every score, is worse than no
figure: it reads as measured, it was used to judge the model, and it set the
target twice as high as the evidence allows.

Two results in that table matter more than the correction.

**The model alone has a NEGATIVE held-out R2**, -0.0467, on the log target. It
is worse than predicting the mean depth everywhere.

**Adding the model to a flexible fit makes the fit slightly WORSE**, +0.1435 to
+0.1411. Given every resolvable field, a learner does better ignoring the
groundwater solve than using it. That is a stronger statement than GW-21's, and
it is the one to quote: the model does not merely fail to add to its recharge
forcing, it is information-negative against its own inputs.

### What this does to GW-21

GW-21 concluded the INPUTS are the limit rather than the formulation. That
survives and is sharpened at both ends. The inputs really are weak: 0.1435
against a between-cell ceiling of 0.6335, so four fifths of what the cell means
could in principle carry is not in elevation, relief, distance to the sea,
recharge or permeability at this resolution. And the model does not reach even
that weak bar, sitting at -0.0467 alone. Both halves are true at once, and the
earlier phrasing implied the second could not be, because it measured the model
against a target that was never real.

## GW-22: the Australian verdict was REGIONAL, and the model has a regime

Australia was chosen first because it is the arid analogue for a world of
endorheic basins. That was the right call for relevance and the wrong one for
diagnosis: it is the regime where this model degenerates. Scored on the United
States, at the same 15.19 km mesh, the same solver and the same settings:

| set | bores | cells | ceiling | Pearson on cell means | model sd | observed sd |
| --- | --- | --- | --- | --- | --- | --- |
| Australia, all | 53,410 | 3,276 | 0.6335 | +0.0703 | 1.91 m | 18.36 m |
| United States, all | 679,924 | 17,763 | 0.5310 | +0.1614 | 34.02 m | 34.26 m |
| United States, CONFIRMED UNCONFINED | 71,265 | 6,574 | 0.8079 | **+0.2596** | 21.03 m | 29.22 m |
| United States, confirmed confined | 47,613 | 4,675 | 0.8071 | +0.3176 | 13.43 m | 32.54 m |

On bores the USGS labels unconfined -- a water table, which is what the solver
computes, rather than the potentiometric head a confined bore reads -- the
correlation is nearly four times Australia's. Australia can only guess at
confinement from bore depth; `aquifer_type_code` states it.

**The model has a regime, and that is the finding.** In Australia its spread is
1.91 m against an observed 18.36: the evapotranspiration sink wins over the
lateral term everywhere, the steady state collapses to the local balance
`d = lambda ln(et_max A / supply)`, and the flow solution contributes almost no
variance. In the United States the spread is 34.02 m against an observed 34.26,
because there the recharge is large enough and the aquifers transmissive enough
that lateral flow sets the head and the solve is doing real work. Sweeping
`lambda` across its whole declared bracket moves the United States spread from
20.94 m to 21.28 m while moving the median from 1.92 m to 7.66 m, which is the
same statement from the other side: there the variance is not the sink's.

### Does it meet the bar?

The declared bar is R2 = 0.07. Direct R2 is negative because the model runs
shallow, a median 5.96 m against 9.50, so the honest test is whether the PATTERN
carries once that bias is removed. Fitting a two-parameter affine correction on
half the cells and scoring the other half, at `lambda` = 2.0, the top of its
declared physical bracket, over 20 random splits:

| | |
| --- | --- |
| held-out R2 | +0.0683, standard deviation 0.0138 |
| range over splits | +0.0443 to +0.0954 |
| splits clearing 0.07 | 10 of 20 |
| Pearson over all cells | +0.2709, so rho^2 = 0.0734 |

**It straddles the bar rather than clearing it.** A single favourable split gave
+0.0748 and reporting that alone would have been a pass on a coin toss; twenty
say the answer is 0.068 give or take 0.014. Against Australia's rho^2 of 0.005
this is a real and reproducible signal, and it is not skill in any strong sense:
the fitted slope is +0.155 on an intercept of 18.96 m, so the calibration shrinks
the model's own variation nearly sixfold and most of the prediction is a
constant.

### What this means for Vesper, PROVISIONALLY, and why it cannot be more

`docs/src/reference/external-data.md` chose Australia because this world is
dominated by arid endorheic basins, and the temptation is to conclude that
Vesper therefore sits in the regime where this model degenerates. That
conclusion is not available yet, and the reason is loop A.

Measured on the CURRENT state: the sink takes a median 97.0% of a land cell's
recharge, and more than 90% of it over 54.1% of land area. On its face that is
the unfavourable regime almost everywhere. But every input to that number comes
from an uncommissioned state.

- The forcing is the BOOTSTRAP climatology. `config/planet.yaml` carries
  `baseline_climatology: null`, and this project's own vocabulary says a
  bootstrap run's numbers are not the baseline. Recharge is what sets the sink
  fraction, and recharge is exactly what a baseline would move.
- The terrain is PRE-CARVE. `lib/orogen.py` records what a carve verdict did on
  the build where one was applied: endorheic land fell from 60.10% to 43.06%.
  So "dominated by arid endorheic basins" is itself a property of the open loop,
  and closing it cut that dominance by more than a quarter on the one occasion
  it has been closed.

So the honest statement is conditional. On the current bootstrap forcing and
pre-carve terrain, most of Vesper's land is sink-dominated and its depth field
is a recharge map there. Whether that survives commissioning is unknown, and
`sink_fraction` is the instrument for asking again rather than the answer. It
should be re-read after every baseline, which costs nothing because the solve
already computes it.

That also cuts at the Australia-first decision, which rested on the same
uncommissioned premise. The reasoning was sound when written and it is still the
right first choice for relevance; it is simply not yet established that the
regime it was chosen to represent is the regime this planet ends in.

### What `sink_fraction` does and does not separate

It is the share of a cell's recharge that groundwater evapotranspiration removes.
Near 1 the depth is the local balance and the flow solve contributed nothing.
Below 1 something else took the water, and there are TWO somethings, which the
Vesper run makes obvious: 37.1% of land has a sink fraction under 0.5 and a
median depth of 0.00 m. Those cells are pinned at the surface and shedding their
recharge as SEEPAGE, not carrying it laterally. A low sink fraction therefore
means "the sink did not set this depth", which is weaker than "the flow solution
did". Reading it as a flow-dominated fraction overstates what it says, and the
pinned flag is what separates the two.

The usable rule it does give is a three-way one, and it needs both fields.
`sink_fraction` near 1 says the depth is the local balance, a recharge map in a
water table's units. Low with `at_surface` set says the cell is pinned and
seeping, where the depth is 0 by construction and carries no information either.
Low with `at_surface` clear is the only case where the lateral term did the work,
and it is the case GW-22 showed correlates against real bores. Anything reading
`depth_m` should ask which of the three it is holding, per cell, which is what
shipping the two fields together is for.
