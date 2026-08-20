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
