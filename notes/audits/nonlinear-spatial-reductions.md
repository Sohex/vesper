# Which mesh-to-grid reductions the order of operations changes, and by how much

**Measured:** the roughness arms, sections 3 and 5, on 2026-08-26 on
`canonical-10m-base`, terrain hash `20046729`; the pedology and saturation arms,
sections 1, 2, 4 and 6, on 2026-08-24 and 2026-08-25 on `precarve-craton-10m`,
terrain hash `ab0d679b`. Both builds carry 10,000,005 regions and every arm is
measured at every rung of the T21/T42/T85/T127/T170 ladder.
`analysis/spatial_reduction_gap.py` is the measurement and
`analysis/spatial_reduction_gap.json` its output; the artifact holds one build
at a time and currently holds the later one, so the sections measured on the
earlier build are quoted here and not reproducible from it until they are
re-run.

This is worldbuilding. Vesper is an invented super-Earth; every quantity below
is a modelled field of that planet or a parameterisation this pipeline applies
to it. Nothing was simulated: no World Orogen generation, no ExoPlaSim run, no
soil solve, no biosphere.

`biosphere/notes/spatial-support-ecological-aggregation-audit.md` finding 4
lists five paths where a cell mean is handed to a nonlinear law and does not
size any of them, so it cannot say which deserve an operator. This sizes them.
The procedure is `hydrography/notes/subgrid-water-table.md` section 2, unchanged:
expand about the cell mean, COMPUTE the correction, and compare it against the
size of the effect the term carries. Its second clause binds every case here --
sub-grid information reaches a cell-scale parameter only as a statistic of the
cell's own distribution, never as a resolved gradient or a within-cell position.

## The bars, and why there is not one of them

Each case is judged against the instrument its own consuming step already
reports itself against, declared before the measurement ran. A single percentage
across six different quantities would have been a preference rather than a
criterion.

| case | instrument |
| --- | --- |
| roughness through the exchange coefficient | the roughness builder's own `z_ref` bracket over the liquid-water span, which is the ignorance it already declares |
| parent texture through weathering | a factor 1.35 in the weathering intensity, Dunne's `S_y.x` of 0.13 log units carried in `pedology/config/pedogenesis.yaml` |
| erodibility through the regolith depth law | `regolith.minimum_depth_m`, the thinnest profile the pedogenesis model distinguishes |
| subgrid elevation through Clausius-Clapeyron | 1% of the saturation vapour pressure, declared rather than sourced. It is the one bar here that is not an instrument, and it has not been moved: there is nothing sourced to move it to, and a bar placed after the result it judges is not a criterion |
| the staged albedo pair through the model's wetting | the accepted baseline's 0.12 W m-2 state-storage tolerance, which `config/partial_surface.yaml` selected the tile operator against. An albedo error reaches it as `delta_alpha * S_down` |
| the aeolian roughness mixture | each class's own declared `bracket` in `aeolian/config/dust.yaml`, the level uncertainty its measurement set supports, converted into a factor on the power of `u*` the arm reports |

Two further cases are not Jensen gaps at all and are reported as what they are:
a calibration that moves with the support, and a statistic that carries a length
from the mesh.

## Result, in one table

| reduction | consuming step | verdict | gap, in the consumer's units |
| --- | --- | --- | --- |
| erodibility mixed before the regolith depth law | `soil` | **MATERIAL, by an order of magnitude** | 0.19 to 0.28 m of regolith at the land mean, against a 0.02 m bar, on 63 to 87% of land area |
| subgrid elevation through saturation vapour pressure | `boundary_conditions`, then evaporation | **NO VERDICT, and the whole of it turns on one number** | 0.32% to 12.0% of `e_sat` against a 1% bar; the bar is crossed once the land-mean lapse rate passes 3.95 to 7.07 K per km, against this planet's dry adiabat of 12.75 |
| surface roughness averaged as a length | `surface_roughness` | **MATERIAL, in the land mean and over most of land area** | -10.3% to -8.1% at the land mean against a 5.1% instrument, past the instrument on 81% of land area at T21 falling to 53% at T170, and -24% at the land mean in the top decile of barren share |
| parent texture mixed before weathering | `soil` | **NOT material, and the premise is wrong** | the clay total is exactly 0.0 at every intensity and every rung; sand and silt trade at most 0.0047 of the land mean, inside Dunne's own scatter everywhere |
| the orographic term across the ladder | `surface_roughness` | DERIVED, and flat where it used to be a calibration | the land mean spreads by 0.9% from T21 to T170, where the solved coefficient it replaced moved by a factor of 1.98 |
| `subgrid_slope` as an elevation spread over the mesh spacing | `soil` | a length carried from the MESH, not a gap | doubled between this project's two builds while the spread it is built from moved 1.6%; the run is now declared, and 11.8% of regolith depth at the land mean rests on the change |
| the staged dry and saturated albedo pair through the model's mixing | `surface_albedo`, then `landmod:wetalb` | **MATERIAL over the middle of the wetting range, and no two-field repair exists** | 0.147 W m-2 of global absorbed shortwave against a 0.12 W m-2 bar, past it from 0.10 to 0.47 of the saturation range and peaking on the model's own evaporation knee |
| the aeolian roughness of a mixed erodible patchwork | `build_dust`, then the in-model emission | NOT material, and the operation was right for the wrong reason | exact for MB95's drag partition, 4.9% on `u*` cubed against the 62.8% the class's own declared bracket carries |

## 1. Erodibility mixed before the regolith depth law. The large one

The law is

    depth = maximum_depth * P / (P + E),   E proportional to erodibility

which is convex in the erodibility, so the mean of the depths lies above the
depth at the mean. `pedology/scripts/build_soil.py` area-mixed every rock class's
erodibility inside a cell and applied the law once to the mixture; this
measurement is what moved it onto `lib/gridding.py`'s `cell_expectation`, which
evaluates the law on each region's own rock and averages the depths, and the
soil report's `regolith_aggregation` block now carries the difference on the map
it built. Erodibility is a property of the rock a region is made of and
spans 0.25 for quartzite to 3.5 for evaporite in the export's own lithology
table, a factor of fourteen, and the area-weighted land mean is 1.70. Both
orders of the reduction are therefore well defined and the difference between
them is the operator's and nothing else's.

Everything the law needs except the erodibility collapses into one ratio,
`rho = E / P` at the cell mean, so no climatology is needed to
size this: `rho` is swept over the decades the law can reach.

| `rho` | depth at the mixed erodibility | land-mean gap | land area above the 0.02 m bar |
| ---: | ---: | ---: | ---: |
| 0.01 | 4.95 m | 0.0002 m | 0% |
| 0.1 | 4.55 m | 0.016 m | 43% |
| 0.5 | 3.33 m | 0.152 m | 78% |
| 1.0 | 2.50 m | 0.262 m | 80% |
| 3.0 | 1.25 m | 0.334 m | 81% |
| 10 | 0.45 m | 0.211 m | 80% |
| 100 | 0.05 m | 0.031 m | 70% |

At T42. The gap is ONE-SIGNED across the whole sweep and every rung: mixing
erodibility first always makes the soil thinner than the cell's own mixture of
rocks would be. That is a pre-registrable invariant for SPAT-8 -- this operator
cannot move the other way -- and it is the one place in this audit where a
universal sign is available, because the curvature of a single saturating
function does not change sign the way a coupled climate response does.

Refinement reduces it and does not remove it. The land-mean gap at `rho = 1`
runs 0.278, 0.262, 0.234, 0.209, 0.187 m across T21 to T170, so T170 still
carries 9.4 times the bar, and the land area above the bar only falls from 87%
to 63%. A finer climate grid is not a route out of this one.

## 2. Subgrid elevation through the model's saturation vapour pressure

`build_boundary_conditions.py` averages elevation over a cell's land and hands
the model one height. Temperature is linear in height through a lapse rate, so
that reduction is EXACT for temperature and there is no Jensen term to find.
Saturation vapour pressure is convex in temperature over the whole liquid-water
span, so the cell's own spread of heights raises the mean saturation deficit
above the one computed at the mean height, and evaporation is what consumes
that.

### Two things the first measurement got wrong, corrected 2026-08-25

Both were corrected before the arm was re-run and the correction was committed
first, because a criterion fixed after the result it judges is not a criterion.
The bar was NOT touched.

**The lapse ceiling was Earth's.** The bracket ran to 9.8 K per km, which is
Earth's `g/cp`. This planet's dry adiabat is 12.75 K per km, and `lib/lapse.py`
already treats it as the ceiling a measured environmental rate must fall below,
raising when it does not. That is the ceiling the bracket now takes. The
correction WIDENS the bracket and moves the arm further from a verdict, which
is the test that it is a correction and not a preference.

**The function was not the model's.** The arm evaluated an idealised
Clausius-Clapeyron with a constant latent heat of 2.5e6 J/kg. The model
evaluates Magnus-Teten through `plasimmod.f90:ra1s/ra2s/ra4s` with two
coefficient sets and a phase switch at `tmelt`, and the quantity being measured
is a curvature, which the two functions do not share. At the cold end of the
reference bracket a cell's high ground falls below `tmelt`, where the model's
coefficients are steeper and the gap is consequently smaller than the liquid
branch alone reports. The arm now reads the coefficients out of `p_earth.f90`
in the manner `lib/sea_water.py` established rather than copying them.

**The floor has nothing under it.** 4.0 K per km is declared, not bounded.
`lib/lapse.py` states that the measured environmental rate sits BELOW the moist
adiabatic rate evaluated at the window-mean state, which is why it deliberately
declines to floor there; nothing at this step keeps the land-mean rate away
from zero, and the gap goes to zero with it. So the low end of the bracket is a
choice and the arm reports it as one.

### The corrected bracket

| rung | 4.0 K/km, 313 K | 4.0 K/km, 273 K | 12.75 K/km, 313 K | 12.75 K/km, 273 K |
| --- | ---: | ---: | ---: | ---: |
| T21 | 0.75% | 1.03% | 6.83% | 12.03% |
| T42 | 0.61% | 0.81% | 5.79% | 10.43% |
| T85 | 0.47% | 0.60% | 4.65% | 8.48% |
| T127 | 0.38% | 0.47% | 3.86% | 7.04% |
| T170 | 0.32% | 0.37% | 3.26% | 5.91% |

Land-mean relative gap in `e_sat`. The bar is crossed on 10.6% to 56.8% of land
area depending on the corner, and the 95th-percentile cell reaches 47.6% at the
steepest corner, so the land mean understates what individual cells carry.

### The verdict, and the one number it turns on

NO VERDICT. The gap crosses the bar inside the bracket, so the measurement has
neither passed nor failed and reporting it as either would be a preference.

What the correction bought is that the undecidedness is now known to lie on ONE
axis. The reference-temperature bracket is worth a factor of about 1.8 across
the whole ladder and never reaches the bar on its own at the declared floor;
the lapse-rate bracket is worth a factor of 12 to 16 and carries the crossing.
So the arm reports the lapse rate at which the land mean reaches the bar, which
converts the deferral into a comparison a climatology settles in one step.

| rung | critical lapse, 273 K | critical lapse, 313 K | fraction of the dry adiabat |
| --- | ---: | ---: | ---: |
| T21 | 3.95 | 4.64 | 0.31 to 0.36 |
| T42 | 4.35 | 5.16 | 0.34 to 0.40 |
| T85 | 4.90 | 5.85 | 0.38 to 0.46 |
| T127 | 5.40 | 6.48 | 0.42 to 0.51 |
| T170 | 5.87 | 7.07 | 0.46 to 0.55 |

K per km, against a dry adiabat of 12.75 K per km. The decision procedure is
fixed here, before the climatology that answers it exists: measure the
land-area-weighted environmental lapse rate with
`lapse.environmental_lapse_k_per_km` on the first accepted baseline
climatology, and compare it against this table at the rung in use. Above it the
reduction needs a sub-grid orographic term in the surface evaporation; at or
below it the reduction is admissible and no operator is needed.

Refinement is not a route out. The critical rate rises by only a factor 1.49
from T21 to T170 while the bracket it is being compared against spans a factor
of 3.2, so a finer rung moves the answer by less than the ignorance does. That
is where section 5 lands on the roughness for a different reason, and it is
why this belongs to SPAT-5's partial-surface decision rather than to a rung
choice.

## 3. Surface roughness averaged as a length. Material, and it grew

`build_surface_roughness.py` gave every mesh region a roughness from its land
cover, area-averaged the LENGTHS over a cell's land, and handed the model one
`z0`. The model then takes `ce = k^2 / ln(z_ref/z0)^2` and the turbulent flux is
linear in `ce`, so what the cell owes the atmosphere is the area mean of `ce`
over its own surfaces, not `ce` of the mean length.

The gap is in the land mean and not only in a tail. Measured on
`canonical-10m-base` on 2026-08-26, the land-mean gap runs -10.3%, -9.8%,
-9.1%, -8.6%, -8.1% from T21 to T170, against an instrument -- the step's own
`z_ref` bracket -- of 5.1%. Twice the instrument at the coarsest rung and still
past it at the finest.

**Its size depends on how rough the field is, and that is why this measurement
moved.** The first pass of it was taken while the field's land mean was solved
onto Earth's, which put an orographic term of several metres under every land
cell and swamped the two-decade contrast between barren ground and canopy; the
gap then read a few tenths of a per cent and the verdict was "not material in
the land mean". Against a derived land mean, where the orographic term is a
small addition to the cover roughness rather than the whole of it, the same
reduction decides a first-order quantity. The operator was already the right
one. What changed is what it is worth.

The land area whose gap exceeds the instrument runs 81%, 74%, 64%, 58%, 53%
across the ladder. **That share falls with refinement**, the opposite of what
the first pass reported, and for the same reason: a coarser cell mixes more
surface classes, so more cells carry the contrast at all, and a fine enough cell
is one class and has no gap to have.

Stratifying by the cell's barren share finds where the field's purpose lives. In
the top decile of barren share at T21 the land-area-weighted gap is -24% and the
1st-percentile cell reaches -34%; those cells are the flat closed-basin floors
the carve verdict integrates evaporation over. Averaging the lengths puts a
canopy minority in charge of a playa cell's exchange coefficient and overstates
its evaporation by a third.

**Corrected.** The builder reduces in `ce` and writes the length that reproduces
the cell's area mean of it. The height that average is taken at is not free and
is not the model's lowest level: Mason (1988) shows the average that reproduces
the correct area-mean surface stress is taken at the BLENDING HEIGHT, and his
Eq (14) puts it an order of magnitude lower on this mesh. `notes/audits/
tuned-values.md` section 9 carries that and the derivation of the orographic
term beside it.

## 4. Parent texture mixed before weathering. The premise does not hold

`build_soil.py:weather_texture` mixes every rock class's parent texture and then
converts weatherable minerals to clay, and finding 4 names it as an
aggregate-then-process defect. It is not one, and the arithmetic says why. The
clay total is

    clay + clip(1 - quartz - clay, 0, 1) * clay_yield * (1 - exp(-clay_conversion * W))

which is AFFINE in the mixed quartz and clay shares wherever the clip does not
bind, and on this world's lithology it never binds. So mixing first and
weathering first give the same clay to the last bit: measured at 0.0 land-mean
gap and 0.0 maximum cell gap, at every intensity in the sweep and every rung.

The split of that conversion between sand and silt is not affine -- it goes
through a ratio and two `minimum` clamps -- and it does move, by up to 0.0047 of
the land mean and 0.0127 in the worst cell at the highest intensity, sand
gaining exactly what silt loses. Against the instrument, what Dunne's own
scatter on the runoff exponent does to the same cell's sand, that is inside the
law's published error at every intensity: 0.0047 against 0.0055 at `W = 6`, and
0.0021 against 0.0222 at `W = 1`.

**The nonlinearity in pedogenesis is real and it is entirely in `W(q, T)`**, a
power of runoff times an exponential of temperature, both of which are cell-mean
CLIMATE quantities. There is no sub-grid population of either -- the climate is
solved at the grid -- so GW-6's constraint applies unchanged: what a cell-scale
parameter may consume is a statistic of a distribution the cell contains, and
for runoff and temperature this project does not hold one. Sizing that term
needs a climatology and a within-cell distribution of it, and neither exists.
The texture reduction itself is admissible at a cell mean and needs no operator.

## 5. The orographic term across the ladder. Was a calibration, is now derived

Not a Jensen gap, and no longer a calibration either. `build_surface_roughness.py`
used to solve a constant relating subgrid relief to a roughness so the land mean
landed on ExoPlaSim's own `dz0land`, and it solved it on whichever grid it was
building for: the constant rose by a factor of 1.98 from T21 to T170 while the
relief it multiplied fell by 5.96. Two rungs built with the defaults therefore
differed by their terrain and by their calibration at once, and finding 8's
convergence question could not be answered from them.

The relation the constant sat in was the defect. Turbulent orographic form drag
is quadratic in the subgrid SLOPE, not linear in relief amplitude, so a
coefficient multiplying an elevation spread has to carry the missing horizontal
scale inside itself -- and the horizontal scale of a within-cell spread IS the
cell. That is why it moved with the rung. The orographic term is now derived from
the mesh's own plane-fit slope through Wood and Mason (1993) Eq (33) and
Beljaars et al. (2004) Eq (6), with no constant to solve;
`notes/audits/tuned-values.md` section 9 carries the derivation.

Measured on `canonical-10m-base` on 2026-08-26, the derived land-mean roughness
across the whole ladder:

| rung | land-mean `z0` | ratio to T42 | median slope variance | median subgrid relief |
| --- | ---: | ---: | ---: | ---: |
| T21 | 0.4800 m | 0.9974 | 3.46e-4 | 153.57 m |
| T42 | 0.4813 m | 1.0000 | 1.44e-4 | 74.77 m |
| T85 | 0.4826 m | 1.0027 | 5.31e-5 | 42.29 m |
| T127 | 0.4836 m | 1.0049 | 2.72e-5 | 31.09 m |
| T170 | 0.4845 m | 1.0067 | 1.66e-5 | 25.58 m |

The spread is 0.9% end to end. These are the arm's own aggregate-then-process
land means, so they are not the field's own -- the builder reduces in `ce` at the
blending height and reports that per rung -- but the comparison across rungs is
like for like and it is the one this section asks about.

The slope variance falls by 21 across the ladder and the land mean does not
follow it, because the drag is a property of the land and the mesh rather than of
the cell: each cell averages the same per-region drag over a different number of
regions. **SPAT-8's constraint is discharged.** A ladder comparison no longer has
to pass one coefficient to every rung, because there is no coefficient; a derived
land mean that moved with the rung would say the scheme was still carrying the
support inside it, and this table is that check.

## 6. `subgrid_slope` measured the mesh, not the gradient

Not a Jensen gap either. `build_soil.py:subgrid_slope` divided the within-cell
elevation spread by the MESH spacing and handed the result to the catena term as
a gradient. The spread is a legitimate cell statistic; the divisor was a length
that belongs to the mesh, so the statistic carried the mesh's resolution into a
quantity the pedogenesis model reads as terrain.

Measured on the same T42 cells, on this project's two builds of the same planet
at the same seed:

| build | regions | mesh spacing | median spread | median `tan beta` | median catena divisor |
| --- | ---: | ---: | ---: | ---: | ---: |
| `precarve-craton` | 2,500,001 | 15.19 km | 73.06 m | 0.004809 | 0.9811 |
| `precarve-craton-10m` | 10,000,005 | 7.60 km | 71.88 m | 0.009464 | 0.9635 |

The elevation spread the statistic is built from moves by 1.6%, which is the
terrain converging as `notes/audits/orogen-resolution.md` says it does. The
inferred gradient DOUBLES, 1.968, because the divisor halved.

### Which mechanism, measured 2026-08-25

A statistic that shifts with the region count has one of three mechanisms and
they take different repairs, so the shift was decomposed rather than assumed:
`analysis/subgrid_slope_support.py`, on the same two builds, at T42, over the
4,629 cells both builds give at least 30 land regions, against the 1.15x
transport bar `orogen-resolution.md` fixed for this class.

| test | what it separates | p50 | p75 | p90 | p95 | p99 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| spread over the mesh spacing, fine/coarse | the shift itself | 2.128 | 2.126 | 2.138 | 2.082 | 2.009 |
| residual once the spacing ratio is divided out | a LENGTH in the definition | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| spread on a seeded 1-in-4 subsample of the fine build | a population bias | 0.968 | 0.994 | 1.001 | 0.999 | 1.007 |
| spread, fine/coarse | the terrain's own support dependence | 1.064 | 1.063 | 1.069 | 1.041 | 1.005 |

The residual is exactly one at every quantile: the explicit divisor is the whole
of the shift. That is the compound topographic index's mechanism, whose `a`
carries a length and which moves by about `ln 2` between the same two builds. It
is NOT the population bias that disqualified a minimum-over-a-ball relief form,
and it is NOT the self-affine mechanism `computeScarpPotential`'s one-edge
gradient has, where the shift drifts from 1.408 at p50 to 1.695 at p99 and
cannot be divided out. The distinction is what makes the repair a normalisation:
over a declared run the same quantiles agree to within 1.069x, inside the bar,
because the numerator is a within-cell statistic of a fixed cell and the terrain
has converged.

So `catena.gradient_baseline_km` declares the run, at 30 km, above Orogen's
measured ~20 km terrain-information floor. It is a unit for the spread rather
than a length the statistic samples at, which is why it can sit below the 90 km
the scarp relief term needs: what has to converge here is the within-cell
spread, and it does. The magnitude is not the terrain's hillslope gradient and never was:
a T42 cell is hundreds of kilometres across, real catenas run at 100 m, and
`slope_transport` is declared against real hillslope gradients rather than
fitted to this distribution. What the field carries is the pattern.

The consequence is in the catena divisor and is not negligible. At the land mean
over those cells, `1 + slope_transport * tan(beta)` was 1.0777 on the coarse
build against 1.1642 on the fine one; over the declared run it is 1.0393 and
1.0416, agreeing to 0.2%. Regolith depth at the land mean on the current build
therefore rises by 11.8%, which is far above the 0.02 m the pedogenesis model
distinguishes. The soil map staged before this change is worthless rather than
stale.

## 7. The wet-soil albedo mixing. Pre-registered before the measurement

### The ocean's zenith fit is not an aggregate-then-process composition here

`radmod.f90`'s ocean albedo is strongly convex in the cosine of the zenith
angle -- the ECHAM-3 branch is `min(0.05/(mu0 + 0.15), 0.15)` -- and the model
does hold a branch that hands it a LONGITUDINAL MEAN of that cosine, taken over
a latitude row. That branch is `ndcycle == 0` and this world does not run it:
`ndcycle`'s compiled default is 1, nothing in `config/planet.yaml` or in
`run_exoplasim.py` assigns it, and `radini` is the only place it could be
assigned. So the zenith cosine reaching the fit is instantaneous and there is no
aggregation in front of it to have an order. The composition is real in the
source and absent from this configuration, which is a different verdict from
"small", and it becomes live the moment a run turns the diurnal cycle off.

### The staged fields

`build_surface_albedo.py` is affine at the crossing everywhere except here, and
that is worth saying first because it is the larger part of the finding. Every
nonlinear step in the builder -- the band split through `band_shapes`, the
per-class wetting ratios, the derived evaporite threshold, the lake paint -- is
evaluated PER MESH REGION and reduced afterwards, so the six staged albedo
fields are area-weighted means of quantities the law has already been applied
to. The builder's own recombination identity, `z1*175 + z2*176 == 174`, is the
check that says so and it is asserted before the write. Snow albedo, sea-ice
albedo and the ocean's zenith fit are all evaluated inside the compiled model on
grid-scale prognostics with no sub-grid population at all, so GW-6's constraint
applies to them unchanged and no operator is warranted.

One composition is not affine. The builder stages a DRY albedo and a SATURATED
albedo per cell, each an area-weighted mean over the cell's own lithology, and
`landmod.f90:wetalb` mixes them per cell through `wet_soil_albedo`, which is
Sadeghi, Jones and Philpot (2015) Eqs (4) and (13). That mixing is linear in the
Kubelka-Munk transform `r = (1-R)^2 / (2R)` and therefore concave in the albedo,
so the cell's own mixture of rocks does not give the mixed albedo of the cell's
mean rock. The sub-grid population is real and is the mesh's `substrate_class`,
which spans two decades of albedo between a dark basalt and a playa.

**The bar.** The consuming step is the surface energy balance, and its declared
instrument is the accepted baseline's 0.12 W m-2 state-storage tolerance --
`config/partial_surface.yaml` selected the tile operator against it and
`ocean/config/transport_loop.yaml` carries it as the transport loop's own
controlling scalar. An albedo error reaches that balance as
`delta_F = delta_alpha * S_down`, so the bar is the albedo difference worth
0.12 W m-2 at the accepted baseline's own land-mean downward shortwave, read
from the climatology rather than assumed.

**The invariants, pre-registered.**

- The gap is EXACTLY ZERO at both ends of the saturation axis. `sadeghi_mix`
  returns the staged dry field at `s = 0` and the staged wet field at `s = 1` on
  their own branches, and both are area-weighted means, so the two orders agree
  to the last bit there. A nonzero gap at either end indicts the arm.
- The gap is exactly zero in a cell of one substrate class, at every saturation.
  This is `cell_mean`'s constant-preservation invariant reached through the
  mixing, and it is falsifiable rather than asserted.
- The gap is ONE-SIGNED in saturation for a fixed cell, because the mixing is
  concave in the albedo pair and the curvature of the Kubelka-Munk inverse does
  not change sign on the physical albedo range. The sign is therefore
  pre-registrable HERE, unlike in the ocean arms, for the same reason the
  regolith depth law's was: a single saturating function has one curvature.
- Refinement must shrink it, because in the limit of one region per cell the
  two orders are the same computation.

### The result

Measured 2026-08-30 on `canonical-10m-carve2`, terrain hash `f496ae9f`, at every
rung, with the watts taken on the accepted baseline climatology, which is T21.

**The one-class control passes.** A cell holding a single substrate class gives
the same answer in both orders to 1.8e-15 in albedo against a 1e-12 bar, at
every rung, and the gap is exactly zero at both ends of the saturation axis
because `sadeghi_mix` returns the staged fields there on their own branch. That
is the arm that could have failed.

| rung | land cells | of one class | band 1 land-mean gap | band 2 | worst cell |
| --- | ---: | ---: | ---: | ---: | ---: |
| T21 | 1,639 | 147 | -1.974e-03 | -1.766e-03 | 1.55e-02 |
| T42 | 5,754 | 848 | -1.821e-03 | -1.641e-03 | 1.65e-02 |
| T85 | 20,719 | 5,283 | -1.587e-03 | -1.450e-03 | 1.68e-02 |
| T127 | 44,392 | 15,032 | -1.398e-03 | -1.297e-03 | 1.68e-02 |
| T170 | 76,168 | 30,657 | -1.246e-03 | -1.174e-03 | 1.68e-02 |

All at the sweep's worst saturation, 0.2853 of the 0 to 0.7609 range the model's
declared endpoints allow. The gap is NEGATIVE, so the modelled wet ground is
brighter than the cell's own mixture of rocks would be and the surface absorbs
less than it should. That is concavity, as pre-registered.

**The pre-registered sign holds on the population and not on every cell.** The
invariant said the gap does not change sign along the saturation axis within a
cell, because the curvature of a single saturating function does not change. It
is violated on 7 of the 1,492 mixed cells at T21 and on 62 of 45,511 at T170 --
0.47 per cent falling to 0.14 -- and that is reported rather than rounded away.
The mechanism is that the mixing is concave in the PAIR and the pair does not
move together: a cell whose two classes have very different wetting ratios has
a saturated spread of a different shape from its dry one, and the curvature of
the composition can change sign between them. The claim was too strong as
stated; what survives it is the land mean, which is one-signed at every
saturation and every rung.

**Refinement shrinks the land mean and does not touch the worst cell.** The land
mean falls by a factor 1.6 from T21 to T170 while the worst cell rises slightly,
because a finer grid has more cells of one class -- 9.0 per cent at T21 against
40.2 at T170 -- and the mixed cells that remain are no less mixed.

**MATERIAL over the middle of the wetting range.** The gap reaches 0.147 W m-2
of global-mean absorbed shortwave against the accepted baseline's 0.12 W m-2
storage tolerance, and stays past it from about 0.10 to about 0.47 of the
saturation range:

| saturation | 0.000 | 0.095 | 0.190 | 0.285 | 0.380 | 0.476 | 0.571 | 0.666 | 0.761 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W m-2 | 0.000 | -0.105 | -0.142 | **-0.147** | **-0.137** | -0.120 | -0.099 | -0.076 | -0.054 |

Taken on the accepted baseline's own downward shortwave per band, land-area
weighted at 79.3 and 92.6 W m-2, and weighted by each cell's land area over the
whole planet so the number is directly comparable with a global storage
tolerance.

**The peak sits on the model's own evaporation knee, which is what settles the
materiality rather than a guess about how wet this world is.** `landmod.f90:181`
puts `drhsfull` at 0.4 of capacity, the fraction above which the wetness factor
reaches one, and `config/planet.yaml` maps a full surface layer to a saturation
of 0.7609. So the knee is at a saturation of 0.304 and the gap peaks at 0.285.
The modelled soil crosses that point every time it wets or dries through the
limiter's ceiling, which is the transition the limiter exists to describe.

### The pre-registered repair does not work, and the reason is a boundary condition

The pre-registration said the mixing is linear in the Kubelka-Munk transform, so
staging the albedo whose transform is the cell's area mean of the per-region
transforms would make the composition exact at every saturation with no
residual. That is wrong at first order and the test is what says so.

| saturation | truth | staged in albedo | staged in the transform |
| ---: | ---: | ---: | ---: |
| 0.000 | 0.275000 | +0.0 | -1.14e-01 |
| 0.285 | 0.239999 | +1.0e-03 | -1.07e-01 |
| 0.761 | 0.199963 | +5.0e-04 | -9.7e-02 |

On a cell of two classes at equal area with albedos 0.10 and 0.45. Staging in
the transform is a hundred times worse than the defect it was meant to remove,
and it fails hardest at zero saturation.

**The two staged fields are not free parameters.** Code 174 IS the dry albedo:
`sadeghi_mix` returns it unchanged at zero saturation and the radiation reads
the same field on dry ground, so it is pinned to the cell's area-mean dry
albedo by a boundary condition that has nothing to do with wetting. The
saturated pair is pinned at the other end by the same argument. There is
therefore no two-field staging that is exact across the saturation axis: what
the model would need is the cell's own curvature, and two numbers do not carry
it.

The mixing is linear in the transform of ITS ARGUMENT and the inverse transform
is applied afterwards, which is the step the pre-registration missed. Reducing
in the transform makes the argument exact and leaves the inverse's own Jensen
term, which is convex where the measured one is concave -- so the two available
stagings bracket the truth from opposite sides rather than one of them reaching
it.

**The disposition is therefore a model change or a declared bracket, and not a
builder change.** Either `wetalb` reads a third staged field carrying the cell's
curvature, or the measured gap is carried as declared model-form error with the
table above as its bracket. Nothing in `build_surface_albedo.py` can fix it,
which is why this finding does not move that builder.

### The model change is priced, and it clears the bar

The choice between those two was not takeable on the numbers above, because
nobody had measured whether the third staged field works. The repair this audit
pre-registered for the same defect was wrong by two orders, so a proposed form
is priced before it is adopted. `analysis/spatial_reduction_gap.py`'s arm now
carries three forms against the same truth, and the criterion -- inside the
0.12 W m-2 storage tolerance at every saturation the modelled column reaches --
was fixed before they were run.

| saturation | 0.000 | 0.095 | 0.190 | 0.285 | 0.380 | 0.476 | 0.571 | 0.666 | 0.761 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| as built | 0.000 | -0.105 | -0.142 | **-0.147** | -0.137 | -0.120 | -0.099 | -0.076 | -0.054 |
| staged in the transform | 0.857 | 0.872 | **0.876** | 0.873 | 0.866 | 0.858 | 0.849 | 0.840 | 0.831 |
| three points | 0.000 | -0.029 | -0.023 | -0.005 | -0.025 | -0.039 | **-0.042** | -0.038 | -0.030 |

Global-mean absorbed shortwave in W m-2, at T21 on the accepted baseline's own
downward shortwave. The three-point form peaks at 0.042 and is inside the bar
at every saturation; the defect peaks at 0.147 and is outside it over three of
the nine points. The albedo side falls with it at every rung: the band 1
land-mean gap from 1.97e-03 to 5.58e-04 at T21 and from 1.25e-03 to 3.42e-04 at
T170, and the worst cell from 1.68e-02 to about 4.6e-03 across the ladder.

**The refuted repair is worse on the mesh than on the two-class cell that
refuted it.** Staging in the transform is 0.83 to 0.88 W m-2 across the whole
range -- six times the defect, seven times the bar, and 0.857 at the zero
saturation the defect gets exactly right. Its worst cell grows with refinement,
from 8.6e-02 at T21 to 2.6e-01 at T170, where every other form's shrinks.

**Why three staged points and not a fitted shape parameter.** Sadeghi's `sigma`
is `s_dry / s_sat` and a per-cell value chosen to force agreement at one
saturation would be the residual of a fit with no derivation to carry to another
planet. The third albedo is not: it is the exact area mean of the per-region
mixed albedo, the same reduction the other two staged fields are, so all three
are exact area means of per-region quantities and the composition is exact at
three saturations by construction.

**The third point's saturation is derived.** `landmod.f90`'s `drhsfull` is the
fill fraction above which the evaporation limiter's wetness factor reaches one,
and `wetalb` maps a fill fraction onto saturation through `skinsrad` and
`skinsrfc`. The knee is the image of one under the other, read from the model
source and from `config/planet.yaml` by the measurement rather than written into
it. That the gap also peaks near there is a result and not the reason.

## 8. The aeolian roughness mixture. Right operation, wrong reason, and immaterial

`aeolian/config/dust.yaml` collapses a patchwork of erodible surfaces to one
roughness per lithology class and then one per cell, as a geometric mean. Two
laws read that single number and they are not nonlinear in the same way, so one
reduction cannot be right for both by assumption. Measured 2026-08-30 by
`aeolian/scripts/roughness_mixing_order.py`; the artifact is
`aeolian/analysis/roughness_mixing_order.json`.

**The drag partition is AFFINE in `ln z0` and the geometric mean is therefore
EXACT for it.** Marticorena and Bergametti (1995) give
`feff = 1 - ln(z0/z0s) / ln(a (X/z0s)^b)`, whose denominator carries no `z0`, so
the area mean of `ln z0` reproduces the area-mean drag efficiency to the last
bit. That is an identity rather than a comparison and it is the arm that can
fail: the two orders agree to 1.7e-16 on the widest mixture in the table and to
exactly zero on the rest, against a 1e-12 bar, and the clip binds on no
endmember of either mixture. Its premise is the clip, and a table whose classes
clipped would break the affinity without breaking the arithmetic, which is why
the arm reports the clip separately from the residual.

**The reason the config gave is not that reason, and the difference is not
pedantic.** It read "the drag goes as 1/ln(z/z0) and it is ln(z0) that
averages". The same sentence carries to the exchange coefficient, where
`ce = k^2 / ln(z_ref/z0)^2` and section 3 above measures 8 to 10 per cent of
land-mean error from averaging in the wrong quantity. A correct operation
resting on a reason that generalises wrongly is one copy-paste from a defect,
and the config and `build_dust.py:mosaic_scalar_z0` now carry the affine
argument instead.

**The friction velocity is not affine in `ln z0`, and that residual is what the
arm sizes.** The same value is written as the namelist `DUSTZ0` the in-model
friction velocity divides by, and `u* = k U / ln(z_ref/z0)` goes as `1/L` with
`L = ln(z_ref/z0)`. `1/L^p` is convex in `L` for every positive `p`, so the area
mean of `u*^p` is at or above `u*^p` at the mean `L` and the geometric mean
UNDERSTATES it. The sign is pre-registrable and one-signed here, for the reason
the regolith depth law's is: a single convex function of one variable does not
change curvature.

| mixture | spread in `ln z0` | `u*` | `u*^2` | `u*^3` |
| --- | ---: | ---: | ---: | ---: |
| `playa_clastic`'s three landform bands | 3.005 | 1.0078 | 1.0238 | **1.0488** |
| `evaporite` against `playa_clastic`, worst share | 1.025 | 1.0015 | 1.0047 | **1.0093** |

Process-then-aggregate over aggregate-then-process, at the shorter end of the
reference-height bracket, which is where the gap is largest. The reference
height is `lapse.reference_height_m` over the same 273.15 to 313.15 K liquid
span the roughness builder brackets over, and it runs 106.2 to 121.7 m on this
planet.

**NOT MATERIAL, by an order of magnitude, against the instrument the step
already declares.** Each class carries a `bracket` on its own `z0` in
`aeolian_z0_by_class_m`, which is the level uncertainty its measurement set
supports. On `playa_clastic` that bracket is worth a factor 1.63 on `u*^3`,
against the 1.049 the order of the mixture is worth. In logs the declared
ignorance is 10.2 times the gap. The between-class mixture is narrower still and
is bounded by the band case rather than argued to be: the two erodible classes
span 1.025 in `ln z0` where the three bands span 3.005, and the arm reports both
spreads so the bound is checkable.

**The band shares come from a build this project no longer holds, and the
verdict does not turn on them.** `playa_roughness_mix.py` refuses to run on the
configured build because the shares are an area-weighted measurement over the
raw mesh and cannot be carried across a generation, so the artifact carries
`precarve-craton` and `precarve-craton-10m` and neither is `canonical-10m-carve2`.
The two it does carry give 1.042 and 1.049 on `u*` cubed, a spread of 0.7 per
cent where the instrument is 62.8, and the endmember roughnesses the mixture is
built from are measured surfaces rather than terrain and do not move at all. A
share vector on the current build cannot bring 1.05 within reach of 1.63.

The class mixture is measured over a declared share sweep rather than the
export's own per-cell shares, and that is a deliberate limit rather than an
omission. Per-cell shares would move the gap between the pure endpoints, where
it is exactly one, and the equal-share case, where it is largest; the arm
already evaluates the largest, so a per-cell measurement can only lower the
answer. It would sharpen a number that is already an order of magnitude inside
its instrument.

## What follows

- The regolith depth law is the one reduction here that needs an expectation
  operator, and it is in `pedology/`. `lib/gridding.py:cell_expectation` is the
  operator; nothing in this batch applied it, because `pedology/` was owned
  elsewhere while this was measured.
- The saturation-deficit term cannot be settled without a climatology, and the
  bracket is reported rather than collapsed. What IS settled is the decision
  procedure: the critical lapse rate per rung is fixed above, so the first
  accepted baseline climatology answers the arm with one comparison rather than
  a re-measurement. It is not evidence for or against a rung.
- The texture reduction needs no operator and finding 4's bullet about it should
  be read as withdrawn: the aggregation is affine and the nonlinearity is in a
  climate variable with no sub-grid population.
- A convergence comparison across the ladder must fix the roughness
  coefficient. `subgrid_slope` transports over its declared run, but its
  magnitude is a spread over a declared length rather than a hillslope
  gradient, so a threshold anchored on measured hillslopes does not belong on
  it at any region count.
- The roughness reduction is corrected in place, and every staged
  `orogen_*_surf_0173.sra` predates the correction.
