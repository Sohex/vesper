# Which mesh-to-grid reductions the order of operations changes, and by how much

**Measured:** 2026-08-24, on `precarve-craton-10m`, terrain hash
`ab0d679b`, 10,000,005 regions, at every rung of the T21/T42/T85/T127/T170
ladder. `analysis/spatial_reduction_gap.py` is the measurement and
`analysis/spatial_reduction_gap.json` its output.

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
| subgrid elevation through Clausius-Clapeyron | 1% of the saturation vapour pressure, declared rather than sourced |

Two further cases are not Jensen gaps at all and are reported as what they are:
a calibration that moves with the support, and a statistic that carries a length
from the mesh.

## Result, in one table

| reduction | consuming step | verdict | gap, in the consumer's units |
| --- | --- | --- | --- |
| erodibility mixed before the regolith depth law | `soil` | **MATERIAL, by an order of magnitude** | 0.19 to 0.28 m of regolith at the land mean, against a 0.02 m bar, on 63 to 87% of land area |
| subgrid elevation through saturation vapour pressure | `boundary_conditions`, then evaporation | **MATERIAL over the upper half of its bracket** | 0.35% to 7.7% of `e_sat`, against a 1% bar; crosses the bar inside the bracket |
| surface roughness averaged as a length | `surface_roughness` | **NOT material in the land mean; material on a growing tail** | 0.29% to -0.01% at the land mean against a 7.0% instrument, but past the instrument on 1.9% of land area at T21 rising to 4.2% at T170, and reaching 30% on the closed-basin floors |
| parent texture mixed before weathering | `soil` | **NOT material, and the premise is wrong** | the clay total is exactly 0.0 at every intensity and every rung; sand and silt trade at most 0.0047 of the land mean, inside Dunne's own scatter everywhere |
| the orographic coefficient re-solved per rung | `surface_roughness` | a support-dependent CALIBRATION, not a gap | moves by a factor 1.98 across the ladder while the relief it multiplies falls by 5.96 |
| `subgrid_slope` as an elevation spread over the mesh spacing | `soil` | a length carried from the MESH, not a gap | doubled between this project's two builds while the spread it is built from moved 1.6%; the run is now declared, and 11.8% of regolith depth at the land mean rests on the change |

## 1. Erodibility mixed before the regolith depth law. The large one

`pedology/scripts/build_soil.py` area-mixes every rock class's erodibility
inside a cell and then applies

    depth = maximum_depth * P / (P + erosion_weight * E),   E proportional to erodibility

which is convex in the erodibility, so the mean of the depths lies above the
depth at the mean. Erodibility is a property of the rock a region is made of and
spans 0.25 for quartzite to 3.5 for evaporite in the export's own lithology
table, a factor of fourteen, and the area-weighted land mean is 1.70. Both
orders of the reduction are therefore well defined and the difference between
them is the operator's and nothing else's.

Everything the law needs except the erodibility collapses into one ratio,
`rho = erosion_weight * E / P` at the cell mean, so no climatology is needed to
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

## 2. Subgrid elevation through Clausius-Clapeyron

`build_boundary_conditions.py` averages elevation over a cell's land and hands
the model one height. Temperature is linear in height through a lapse rate, so
that reduction is EXACT for temperature and there is no Jensen term to find.
Saturation vapour pressure is exponential in temperature, so the cell's own
spread of heights raises the mean saturation deficit above the one computed at
the mean height, and evaporation is what consumes that.

Bracketed twice, because neither factor is known before a climatology exists:
the lapse rate over 4.0 to 9.8 K per km, and the reference air temperature over
the liquid-water span.

| rung | 4.0 K/km, 313 K | 4.0 K/km, 273 K | 9.8 K/km, 313 K | 9.8 K/km, 273 K |
| --- | ---: | ---: | ---: | ---: |
| T21 | 0.82% | 1.41% | 4.54% | 7.71% |
| T42 | 0.67% | 1.15% | 3.81% | 6.60% |
| T85 | 0.52% | 0.90% | 3.04% | 5.34% |
| T127 | 0.42% | 0.74% | 2.51% | 4.44% |
| T170 | 0.35% | 0.62% | 2.11% | 3.76% |

Land-mean relative gap in `e_sat`. The bar is crossed on 11% to 56% of land area
depending on the corner. So this one CROSSES ITS BAR INSIDE ITS BRACKET and is
not resolved here: at the dry-adiabatic, cold corner it is material at every
rung, and at the moist, warm corner it is below the bar at T85 and above. What
settles it is a climatology, not a finer grid, and the honest statement is the
bracket rather than a number. The 95th-percentile cell reaches 26% at the wet
corner, so the land mean understates what individual cells carry.

## 3. Surface roughness averaged as a length. Not material, except where it is

`build_surface_roughness.py` gave every mesh region a roughness from its land
cover, area-averaged the LENGTHS over a cell's land, and handed the model one
`z0`. The model then takes `ce = k^2 / ln(z_ref/z0)^2` and the turbulent flux is
linear in `ce`, so what the cell owes the atmosphere is the area mean of `ce`
over its own surfaces, not `ce` of the mean length.

In the land mean the two agree: 0.29% at T21 falling through 0.18%, 0.08%,
0.03% to -0.01% at T170, against an instrument -- the step's own `z_ref` bracket
-- of 7.0%. Twenty-four times inside it at the worst rung. By the declared bar
that is noise.

It is not noise everywhere. The 1st-percentile cell carries -13% at T21 and
-22% at T170, and the land area whose gap exceeds the instrument runs 1.9%,
2.4%, 3.2%, 3.9%, 4.2% across the ladder. **That share grows with refinement**,
which is the opposite of what a reader would assume, and the reason is that a
finer cell is more often DOMINATED by one surface with a minority of the other,
which is exactly where a logarithm parts company with its argument.

Stratifying by the cell's barren share finds the tail where the field's purpose
lives. In the top decile of barren share the 1st-percentile cell reaches -31%,
and those cells are the flat closed-basin floors the carve verdict integrates
evaporation over. Averaging the lengths puts a canopy minority in charge of a
playa cell's exchange coefficient and overstates its evaporation by up to a
third.

**Corrected.** The builder now reduces in `ce` and writes the length that
reproduces the cell's area mean of it. The inversion needs `z_ref`, which is
unknown before a climatology, so it is taken at the midpoint of the same
liquid-water bracket the report already declares -- and the whole bracket is
worth 1.5e-4 of the land mean, measured, so the anchor choice is not a free
parameter in disguise. At T21 the land mean stays at the anchored 2.0 m by
construction, the roughest land cell falls from 7.572 to 7.526 m and the
smoothest from 0.003781 to 0.003772 m.

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

## 5. The orographic coefficient is re-solved on every rung

Not a Jensen gap. `build_surface_roughness.py` solves the constant relating
subgrid relief to a roughness so that the land mean lands on ExoPlaSim's own
`dz0land`, and it solves it on whichever grid it is building for.

| rung | coefficient | ratio to T42 | median subgrid relief |
| --- | ---: | ---: | ---: |
| T21 | 0.004088 | 0.820 | 141.0 m |
| T42 | 0.004985 | 1.000 | 66.6 m |
| T85 | 0.006149 | 1.234 | 38.3 m |
| T127 | 0.007147 | 1.434 | 28.5 m |
| T170 | 0.008093 | 1.623 | 23.7 m |

The relief falls by 5.96 across the ladder, as it must -- a smaller cell holds
less of the terrain's variance -- and the coefficient rises by 1.98 to hold the
anchor. Their product still falls, from 0.576 to 0.192 m, so the orographic
contribution to `z0` is genuinely smaller at fine support AND the constant that
sets it has moved. Two rungs built with the defaults therefore differ by their
terrain and by their calibration at once, and finding 8's convergence question
cannot be answered from them.

The solve stays, because it is the anchoring argument the field rests on and
removing it would move the global roughness ExoPlaSim was tuned against.
`--orographic-coefficient` supplies one value to a whole ladder comparison
instead, the default still solves, and a single-rung build is unchanged. SPAT-8
must pass it.

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

## What follows

- The regolith depth law is the one reduction here that needs an expectation
  operator, and it is in `pedology/`. `lib/gridding.py:cell_expectation` is the
  operator; nothing in this batch applied it, because `pedology/` was owned
  elsewhere while this was measured.
- The saturation-deficit term cannot be settled without a climatology. It is
  bracketed and the bracket is reported; it is not evidence for or against a
  rung.
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
