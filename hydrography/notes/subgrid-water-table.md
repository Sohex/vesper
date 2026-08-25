# What is sub-grid about the water table, and what a cell-scale parameter may claim

This is worldbuilding. Vesper is an invented super-Earth; everything here is
about the steady-state groundwater model written for it, the mesh it is solved
on, and the Earth measurements it has been scored against.

Everything measured here was measured on 2026-08-24, on `precarve-craton-10m`
where a build is named and on `precarve-craton` where the two are compared. The
Earth figures quoted are carried from `earth-calibration-criterion.md` and
`groundwater-et-sink.md` and are not re-derived.

GW-6 is not a defect and is not a task waiting for effort. It is a STANDING
CONSTRAINT on what may be written into this component, and this note states it
as a rule with a decision procedure, applies that procedure to the three
parameters that were queued behind it, and records the arithmetic each verdict
rests on. `groundwater-scoping.md` section 5 raised it, `groundwater-et-sink.md`
and `earth-calibration-criterion.md` measured it five times, and none of them
says what a future parameter is allowed to claim. That is what is missing and
what this supplies.

## 1. The finding, in one paragraph, because it is settled

Fan et al. (2013) put the well-articulated water table gradient at the valley to
ridge transect, decameters to kilometres. A region on this project's mesh is
7.60 km on the active build and 15.19 km on the one before it, so that transect
is inside a cell. Five separate measurements say it stays there. The solved
depth field is very nearly a monotone function of its recharge forcing, Spearman
-0.90 to -0.93, where the real one is not, -0.27. Halving the cell edge raises
the resolvable ceiling and more than doubles the model's own spread while moving
the correlation not at all, +0.0703, +0.0665, +0.0711. A generation at four
times the region count adds at most 12% to relief at any separation from 5 to
200 km and the pre-erosion substrate is converged outright, so a finer generator
does not hold the missing texture either. Draining every cell to its own DEM
minimum, which is the bounding limit of the whole sub-grid-drain family,
inverts the correlation rather than improving it. And where a bore sits inside
its own cell out-predicts the groundwater solve on Pearson, +0.1290 against
+0.0703, so the strongest single term available at this scale is geometry with
no physics in it at all.

The one place the model does work is where lateral flow rather than the local
sink sets the head: on United States bores labelled unconfined it reaches
Pearson +0.2599 with a spread of 20.72 m against an observed 28.94. So the gap
is regional, not universal, and `sink_fraction` in `water_table.nc` is what
tells the two apart per cell.

## 2. The rule this licenses

**A cell-scale parameter may not be justified by a mechanism whose length scale
is below the cell, and the test for that is whether the parameterisation has a
CELL MEAN.**

The test is arithmetic and not a judgement. Expand any parameterisation `P(h)`
about the cell-mean head with within-cell spread `sigma`:

    E[P(h)] = P(hbar) + (1/2) P''(hbar) sigma^2 + ...

- `P` LINEAR in the head: the correction is exactly zero at every `sigma`. The
  parameterisation may be evaluated at a cell mean at any resolution, and the
  cell size never enters.
- `P` convex or concave: there is a correction, its size must be COMPUTED, and
  it must be compared against the size of the effect the term is being added
  for. A term whose cell-mean correction is larger than the effect it carries is
  not a coarse version of the mechanism; it is a different number.

Two worked cases, and they are three orders of magnitude apart:

**Fan's exponential decay, `T = A exp(h/f)`.** Convex, so the cell mean carries
`exp(sigma^2 / 2 f^2)`. Fan's curve reaches `f = 0.95 m` on steep bedrock. The
elevation spread among a region and its neighbours on this world's land is a
median 20.4 m on the active build and 33.8 m on the one before it, with the
upper quartile at 117 and 177 m -- and that is relief BETWEEN cells, which GW-9
is explicit is a different quantity from relief WITHIN one and is used here only
as an order of magnitude. At `sigma / f` of order 20 to 180 the factor is
`exp(200)` to `exp(5000)`. INADMISSIBLE, and no accuracy of implementation
recovers it.

**The unconfined saturated column, `T = K max(h - z_bottom, b_min)`.** Linear
above the floor, so the Jensen term is exactly zero there; the floor is the only
place the form bends at all. For a Gaussian head with mean column `m` and spread
`sigma`, the relative excess of `E[max]` over `max(E)` is

| `m / sigma` | relative excess |
| ---: | ---: |
| 0.5 | 39.6% |
| 1 | 8.33% |
| 2 | 0.424% |
| 3 | 0.0127% |
| 5 | 1.1e-8 |

So the unconfined form is admissible wherever the saturated column is about
three times the sub-grid relief and no longer. The column is the thickness less
the depth, so the worst case is a thin aquifer with a deep table, and the
thickness floor is its upper bound: at the 100 m floor against a median 20 to
34 m of relief the ratio is at best 3 to 5 and the correction at best under
0.02%, while against the roughest quartile's 117 to 177 m the ratio falls below
one and the correction reaches tens of percent. **That is the constraint GW-6
puts on GW-18's floor**, and it is a reason to keep the floor at Gleeson's 100 m
rather than lower it, stated in advance of any run. It is also the reason the
run must report `at_transmissivity_floor`: a cell whose column has reached
`min_saturated_thickness_m` is on a numerical bound, and its depth is a lower
bound rather than a value.

**The second half of the rule.** Sub-grid information may reach a cell-scale
parameter ONLY as a statistic of the distribution the cell contains -- a
fraction, a rank, a quantile -- and never as a resolved gradient or as a
position within the cell. GW-25 is why the second clause is there: folding a
bore's position within its own cell into the model's score nearly doubled its
Pearson, and the term doing the work was the geometry, not the physics. A model
that consumes a within-cell position is scoring a covariate.

## 3. GW-24: depth-dependent transmissivity, and what it is worth

The form is `T = K (h - z_bottom)`, the textbook unconfined Dupuit-Forchheimer
case, and the current `T = K D` is its confined approximation. It is real
physics, the config has always said so in as many words, and under
`docs/src/practice/failure-modes.md` class 16 the fact that it does not improve
an agreement is not a reason to leave it out. It is implemented, it passes the
checks below, and it is OFF by default. The reason is size, and the size is
measurable rather than arguable.

### What it is worth in metres of head

For one-dimensional drainage of recharge `R` over a length `L` to a fixed head,
the confined mound above the outlet is `R L^2 / 2 K D`, and the unconfined one
is `D (sqrt(1 + 2 H_c / D) - 1)`. For `H_c` small against `D` that is
`H_c (1 - H_c / 2D)`, so

    the unconfined correction, as a fraction of the head, is H / 2D

which is the fraction of the saturated column the water table has drained. That
is not an estimate that needs a solve: it is `depth / D` per cell, both of which
the CONFINED run already writes. `build_groundwater.py` now reports it as
`saturated_column_drained`, and the decision to enable `--unconfined` rests on
that number from the run in hand rather than on this paragraph.

Anchored on the Australian case the Earth comparator is scored against -- 8 mm/yr
of recharge, 15 km to a local baselevel, Gleeson's median conductivity of
7.5e-7 m/s:

| aquifer thickness | confined mound over 15 km | unconfined correction |
| ---: | ---: | ---: |
| 100 m | 381 m, so the table pins | not a perturbation at all |
| 365 m, this world's sourced median | 104 m | 14%, about 15 m |
| 1,000 m | 38 m | 1.9%, about 0.7 m |
| 2,000 m, GW-17's range fix | 19 m | 0.48%, about 0.09 m |

**Against the instrument.** The solver resolves head to about a micrometre: its
uniqueness identity is declared at 1e-9 relative and measured at 0.000e+00, and
the water balance closes at 1.7e-16 against a declared 1e-10. So the tolerance
is not what limits this. What limits it is the model's own between-cell spread,
1.91 m in Australia and 20.72 m in the United States, against observed spreads
of 18.36 and 28.94 m. At the 2 km thickness GW-17 needed to reproduce the
observed depth range the correction is 0.09 m, which is 0.4% of the model's own
spread and cannot move a correlation. At 365 m it is 15 m, which is comparable
to the spread -- but there the confined mound is already 104 m over 15 km, far
above the median relief the terrain offers over that distance, so the table pins
and the evapotranspiration sink rather than the flow solve sets the depth.

**And the MESH is not the instrument either.** GW-8 records the discrete
operator at about 12% relative RMS against its analytic eigenvalue, and that has
been read as a noise floor on any field solved with it, which would put this
whole table out of reach before it was measured. It is not one. That 12% is a
TRUNCATION error and does not converge at all: `notes/mesh-geometry.md` fits an
order of 0.004 in the cell size over a sixty-four-fold refinement of Orogen's
own generator. The solved head's error is a different quantity, converges at
second order, and is 2.8e-06 relative on `precarve-craton`. Every entry in the
table above is orders above that, so the discretisation is not what decides
this.

**Nor does the ledger refuse it.** `config/land_water_ledger.yaml` declares
`transient_saturated_storage` unrepresentable and owned by PLHY-4, and that
absence is seasonal head, aquifer storage and a transient capillary or baseflow
response. `T = K (h - z_bottom)` evaluated at the steady head is still an
equilibrium and carries no storage term, so the groundwater store's declared
limit is untouched. The same holds for GW-18, whose thickness is a static field.

**So the term is large only in the configuration the model is already known to
be wrong in, and negligible in the one that reproduces the observed range.**
That is the finding, and it is why the default does not move.

### What it costs, which is not nothing

`T` depends on the head, so:

- The matrix is no longer fixed. The problem stops being a LINEAR
  complementarity problem, and the argument that licensed `--uniqueness-check`
  as an IDENTITY -- one symmetric positive definite matrix, therefore exactly
  one solution -- lapses. Under the unconfined form that check is a MEASUREMENT
  against `UNCONFINED_HEAD_RELATIVE`, declared at 1e-6 before any unconfined run.
- The static dry set survives, but only because of the floor. A face whose
  saturated column reached zero would stop conducting, so which cells connect to
  recharge or to the sea would depend on the head, and the connectivity query
  that made the uniqueness identity pass bit-identically would stop being
  static. `min_saturated_thickness_m` is strictly positive for that reason and
  not for a numerical one.
- The Picard iteration's rate is governed by the same number as the effect. From
  `groundwater.py --dupuit-test`, on a one-dimensional aquifer:

  | saturated column contrast across the domain | passes to the bar |
  | ---: | ---: |
  | 1.12 | 8 |
  | 1.60 | 12 |
  | 5.07 | 18 |
  | 24.9 | not converged in 60 |

  Where the unconfined form barely changes the answer it converges in eight
  passes; where it changes the answer a great deal it does not converge in the
  budget this component runs. One number sets both.

### The checks, declared before they were run

`groundwater.py --dupuit-test`, on a synthetic one-dimensional aquifer with an
analytic answer and nothing about this world in it.

**The CONFINED path first, against its own analytic answer.** The two forms now
share every line of `solve` except the reassembly, and a shared line is a shared
failure, so the guard against the unconfined arrival quietly changing the
confined model is a check rather than a claim. One pass and a relative error of
2e-15 to 2e-14 at 100 m, 2 km and 20 km of thickness.

**The Kirchhoff identity, declared exact at 1e-12.** Over a FLAT aquifer base
the arithmetic-mean face thickness makes the discrete unconfined flux exactly
`g K (u_j - u_i)` for `u = b^2/2`, so the nonlinear scheme and the linear one in
`u` are the same matrix. Measured at 4.6e-15 and 1.1e-14 at 100 and 200 cells.
PASS.

**The analytic Dupuit parabola, declared at 1e-3 relative.** Measured at 3.6e-12
and 1.4e-11. PASS, and by nine orders, because the three-point scheme reproduces
a quadratic exactly at cell centres and the Kirchhoff identity makes the problem
quadratic.

**A SLOPING aquifer base is a different problem, and this is the number that
matters for GW-18.** With the base tilted 1 in 1,000 -- 100 m of base relief over
100 km, which is less than this world's terrain does -- the Kirchhoff identity
departs by 17.6% and the solution departs from the flat-base parabola by 15.4%.
Reported and not judged: it is measuring the size of the term the transform
drops, not an error. A real aquifer base follows the terrain, so no Kirchhoff
transform exists for the case this model would actually run, and the exactness
above is the instrument rather than the solver.

**A round-off floor on the residual bar, which the confined model never met.**
The declared bar is 1e-12 of total recharge on the free-cell water balance.
Under the confined form the direct solve makes that balance exact whatever the
transmissivity. Under the unconfined form the balance is re-evaluated with the
conductance the NEW head implies, so it is a difference of face fluxes formed
from heads whose magnitude is the aquifer thickness while their difference is a
fraction of a metre, and the cancellation puts a floor at about

    eps * sum_faces Trans (|h_src| + |h_dst|) / total recharge

`solve` computes that floor, uses it as the bar where it exceeds 1e-12, and
REFUSES to report a convergence at all where it exceeds the closure bar of
1e-10, because such a run could not be certified. On the one-dimensional case
the observed plateau sat at 0.3 and 0.7 of the un-inflated bound at 400 m and
2 km of saturated column, so the bound is the right shape. On a real mesh the
denominator is a planet's recharge and the floor is many orders below the bar;
the guard exists so that a configuration which crosses it says so.

### What must be run before this can be turned on

Nothing here has touched the real mesh, and the risk in this change has always
been convergence rather than algebra. Two solvers died on the previous attempt
at a depth-dependent transmissivity.

    python hydrography/scripts/build_groundwater.py --unconfined --uniqueness-check
    python hydrography/scripts/build_groundwater.py --unconfined --reduction-test
    python hydrography/scripts/earth_calibration.py   # with the unconfined arm

and the questions they answer, in order: does it converge at all on 10 million
regions within the pass budget; does the reduction identity still hold bitwise
at zero permeability; do the two active-set trajectories land within 1e-6
relative; and does the dependence structure that GW-21 showed is reachable by
transmissivity MAGNITUDE alone become reachable in a way that also moves the
agreement, which is the thing raising `D` uniformly never did.

## 4. GW-18: a sourced thickness, and what it does to GW-17's range fix

The constant 100 m is Gleeson's own "on the order of 100 m" for the depth his
lithology maps describe. GW-17 then found that a uniform 2 km with local
baselevels takes the 95th percentile depth from 6.9 m to 52 m against an
observed 42 m, which fixed the RANGE. That 2 km is not sourced from anything.

**The only sediment thickness this project holds is the export's own
`cover_thickness`**, the surviving thickness of the sedimentary and volcanic
cover over basement. It is a STATE and not a history, which is why
`docs/src/reference/no-time-axis.md` does not refuse it: it is the thickness
that is there now, not an accumulation over a duration.

Measured on the active build, floored at Gleeson's 100 m: 35.9% of land sits on
the floor because it carries no cover at all, the land median is 365 m, the 95th
percentile 998 m, the maximum 2,950 m, and the area-weighted mean 380 m. **Only
0.4% of land reaches the 2 km that GW-17 applied everywhere.**

So the honest statement is not that GW-18 sources the range. It is that a
sourced thickness CONTRADICTS the assumed one over 99.6% of the land, and puts
2 km only in the deepest basins. The prediction, recorded in
`config/groundwater.yaml` before the run: the depth range shrinks, and cratonic
cells move from lateral-flow control back towards the local
recharge-and-evapotranspiration balance, which `sink_fraction` reports per cell.
GW-18's own row already expected no gain in PATTERN -- no resolvable field
predicts observed depth above Spearman 0.27, and height above the nearest river,
the mechanism thickness would strengthen, is the weakest at 0.064 -- so a loss of
range with no gain in pattern is a real possible outcome and the run is what
decides.

**The floor is not adjustable downward.** Section 2's table is why: at 100 m
against this world's sub-grid relief the unconfined form's cell mean is good to
0.02% on median ground and to tens of percent on the roughest quartile. A
thinner floor puts the parameterisation into the range where it has no cell mean
at all, which is the same disqualification Fan's `f` earned.

**And the operator carries a thickness field at the same relative accuracy it
carries anything else.** A spatially varying `D` is a multiplier on the
transmissivity spanning a factor of thirty across this world's land, which is an
order-one perturbation to the operator and not a small one.
`notes/mesh-geometry.md` scores the operator on a transmissivity perturbation
against an analytic answer and finds the perturbation carried with a relative
error of 0.12 OF THE PERTURBATION, the same at an amplitude of 0.01 as at 0.2,
because the confined assembly is linear in `T`. So the discretisation does not stand between
GW-18 and a measurable answer. A run does.

**It costs nothing structurally.** A spatially varying `D` under the confined
form leaves `T` independent of the head, so the matrix stays fixed, the problem
stays a linear complementarity problem, and every identity this component is
certified by is untouched. That is the whole difference between GW-18 and GW-24,
and it is why they are not the same change even though they read the same field:
GW-18 hands the solver a multiplier, GW-24 hands it a datum, and only the second
changes the class of the problem.

## 5. GW-26: the index, what transports, and what does not

A saturated FRACTION is the quantity GW-6 leaves reachable where a per-cell
DEPTH is not, and `build_topographic_index.py` computes it. Three things about
it are decisions rather than implementation.

**A fraction is a fraction of a population, and there is only one population.**
`f_sat_max` is a rank statistic over the terrain inside a cell. The only such
population this project holds is the mesh regions inside a climate-grid cell, so
`f_sat_max` lives at the GRID and there is no per-region saturated fraction.
Manufacturing one would need a sub-region hypsometry, which is exactly what
section 1 says does not exist. That is a constraint on WET-2, SURF-7 and LSHY-6,
which are the three rows waiting on this: what is available to them is a
climate-grid-cell fraction, not a native-mesh one.

**And the population is an AREA, not a region count.** `f_sat_max` multiplies
into a saturated fraction OF A CELL, so what it ranks over is that cell's land
area. CLIMBER-X's tabulated CDF is over equal-area DEM pixels, where the two
coincide; the land regions of this mesh are not equal-area and span a factor of
5.2 from their 5th to their 95th percentile, at a coefficient of variation of
0.43. A region count weights a sliver and a full cell alike, and it does so
twice: in the share, and in the cell mean the share is taken about, which is the
class boundary. Measured on the active build, moving both to area weights shifts
the cell-mean index by a median of 0.19 -- the same order as the 0.23 between
the two slope arms, which this note already ships as a bracket -- and moves
`f_sat_max` by a median of 0.004 at T21 and by up to 0.30 in the worst cell.
`lib/gridding.py` owns the operators: `cell_moments` for the area-weighted mean
and the within-cell spread, `cell_fraction` for the categorical area share.

**The sub-grid population is real and is worth having.** On the active build,
aggregated onto the configured grid, the index's spread WITHIN a cell has a
median of 2.22 against a spread of cell means BETWEEN cells of 1.82. More of the
index's variance is inside the cells than between them, so what a cell-mean
depth throws away is the larger half. That is the case for computing it at all.

**What the three waiting rows get, which is world-d9u4's decision.** They take
the grid-cell fraction and change quantity, which is the third of the three
routes that row named. Section 2's rule decides it rather than a preference
does: sub-grid information reaches a cell-scale parameter only as a statistic of
the cell's own distribution. Parameterising a sub-region hypsometry from the
finer build's statistics puts a distribution BELOW the cell, which section 1
says is not there; letting each consumer downscale by its own stated rule is the
same invention moved to where three rows would each make it differently. Only
the third route stays inside what the terrain supports, and it is the route
`groundwater_access` already took for the same quantity class under PLHY-5,
where a per-cell root-access depth was replaced by a per-cell AREA.

Per row, and two of the three turn out not to need anything built:

- **LSHY-6 is already served.** It asks for the native-mesh to climate-grid
  crossing to preserve the mosaic and explicitly rejects a cell-mean water table
  depth as an area proxy. A climate-grid saturated AREA is what it asked for;
  the fraction arrives as one of its response-tile shares, alongside PLHY-5's
  accessible area. Its dependence on a native-mesh fraction is deleted, not met.
- **SURF-7 never needed this quantity.** A groundwater discharge mask is a
  statement about where the solved water table meets the surface, and that is
  RESOLVED: `water_table.nc` carries `at_surface` and `seepage_m3_s` per region
  from the solve. `f_sat` is a statistic about the terrain a cell does NOT
  resolve, so it is not the finer version of that mask and cannot be substituted
  for it. SURF-7 takes the per-region discharge fields it already has, and its
  dependence on this row is deleted.
- **WET-2 is the only row that genuinely wanted a below-mesh fraction**, and it
  is the one that changes quantity. Its resolved classes -- open lake, river and
  playa water, and dry mineral soil -- stay per region, from `surface_water.nc`
  and the depth field. Its saturated non-inundated mineral class is unresolved
  and crosses as a climate-grid AREA share. Mutual exclusivity, which that row
  requires, is then a constraint stated at two supports rather than one: the
  resolved classes partition a cell's land area between them, and the saturated
  share is taken out of what they leave. That is a real cost of the decision and
  it belongs to WET-2 to carry.

**The absolute scale does not transport, and the rank statistic does.**
CLIMBER-X keys three wetland schemes on absolute index values: a CDF tabulated
on integer bins 1 to 15, `cti_mean_crit = 5.5` below which no wetland forms, and
a cut at 14 above which the maximum wetland fraction is set to zero. Those are
calibrated against an index computed on Earth at about a kilometre. The index
computed here has a land mean of 15.9 and a cell mean above 14 on 83.3% of the
cells that hold enough regions to estimate one, so applied literally that rule
would return no wetland almost everywhere -- a statement about two indices'
scales, not about this world. The scale is not fixable either: `a` carries a
length, so the whole distribution shifts with the mesh, and the measured shift
between this project's own two builds, 7.60 against 15.19 km, is 0.86 on the
plane-fit slope arm and 0.80 on the receiver-drop arm, both close to `ln 2`.
`f_sat_max` is a share of a cell's own population above that cell's own mean, so
it survives any additive shift; an absolute threshold does not. The config
refuses the absolute thresholds by name and the script writes no product that
keys on one.

**The slope is the same wrong slope GW-9 found**, one level removed: the index
takes `tan(beta)` and the mesh's slope estimators are regional dips rather than
hillslope gradients, which understates the gradient and therefore overstates the
index, one-signed. Both estimators are computed and shipped, and they differ by
0.23 in the land mean -- 15.91 for the plane fit against 16.14 for the receiver
drop -- so the slope choice is NOT what puts this index off CLIMBER-X's scale.
The cell size is. What the two arms DO differ on is how much land they push onto
the slope floor, 9.0% against 16.6%, and the receiver drop is higher because a
depression-filled surface is exactly flat inside a filled pit.

**Nothing may consume `f_sat` until it is scored.** The bar is declared in
`hydrography/config/topographic_index.yaml`, before any fraction was computed:
the area under the ROC curve for an observed water table within a metre of the
surface must exceed 0.573, which is what the model's cell-mean depth reaches
alone. The attribution is exact rather than argued, because `f_sat` is strictly
monotone in the depth for a fixed `f_sat_max`: ranking cells by `f_sat` with a
constant `f_sat_max` reproduces ranking them by depth exactly, so any gain is
`f_sat_max`'s and nothing else's. A miss means the index adds nothing to the
depth it multiplies, and the fraction should be withdrawn rather than shipped
with a caveat.

The estimator question GW-26 flagged is unchanged and is not settled here.
PALADYN and ClimaLand both estimate a grid-cell mean water table from column
water content rather than from a lateral solve, so the column estimator is the
field's normal practice and this project's Dupuit-Forchheimer solution is the
departure needing the argument. Nothing in this note decides that; it is a
decision about which estimator the world's numbers rest on.
