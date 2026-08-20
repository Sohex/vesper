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
