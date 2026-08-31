# The ocean's support crossing: which reductions the order of operations changes

This is worldbuilding. Vesper is an invented super-Earth; the bathymetry,
coastline and ocean below are modelled fields of that planet, and every bar
quoted is a tolerance this project has already declared for one of them.
Nothing is simulated here: no World Orogen generation, no ExoPlaSim run, no
ocean run.

`notes/audits/ocean-and-marine-biosphere.md` section 8b says the ocean has
another support crossing and that "narrow connections and partial coastal cells
can change circulation even when cell-mean depth is unchanged". Section 8d says
"reached the ocean" is a boundary condition whose shelf and coastal geometry the
delivery scales with. Neither sizes anything. This sizes them, on terrain alone,
because bathymetry and coastline are terrain and not climate.

`notes/audits/nonlinear-spatial-reductions.md` is the same audit on the land
fields and its procedure is unchanged here: expand about the cell mean, COMPUTE
the correction, and compare it against the size of the effect the term carries.
Its second clause binds every case below -- sub-grid information reaches a
cell-scale parameter only as a statistic of the cell's own distribution, never
as a resolved gradient or a within-cell position.

`analysis/ocean_support_reduction.py` is the measurement and
`analysis/ocean_support_reduction.json` its output.

## Result, in one table

Measured 2026-08-30 on `canonical-10m-carve2`, terrain hash `f496ae9f`,
10,000,005 regions, at every rung of the T21/T42/T85/T127/T170 ladder and on the
36 x 36 GOLDSTEIN grid OCN-3's shipped configurations use. Ocean is 56.75 per
cent of the planet's area by `surface_class`.

| arm | verdict | in the bar's units |
| --- | --- | --- |
| ocean volume through the cell-mean depth | CONTROL, PASSES | relative residual 0 to 1.7e-16 against 1e-12, at every support |
| the horizon-area gap on a uniform ocean | CONTROL, PASSES | exactly zero, at every support and every horizon |
| the gap is confined to cells the horizon cuts | CONTROL, PASSES | every horizon at every support, 150 cases |
| ocean area on either side of a depth horizon | **MATERIAL, by one to two orders** | +12 to +49 times the bar at its positive peak and -8 to -42 at its negative one, on every support including T170 |
| deep-water connectivity | **MATERIAL, and it does not converge** | the detached area disagrees with the mesh's by up to 83 Mkm2, and the disagreement neither shrinks nor keeps its sign along the ladder |
| the binary 0.5 wet mask's own area error | a second mechanism, reported separately | 3.3 to 5.0 times the bar on the Gaussian ladder and 12.3 on the GOLDSTEIN grid, at zero horizon |

## 1. The three controls

**Volume.** The cell-mean depth reproduces the mesh's ocean volume to between 0
and 1.7e-16 relative at every support, against a 1e-12 bar. That is the affine
identity and it is the arm that could have failed: a binning error, a wrong area
weight or a mismatched population all break it and none of them breaks the arms
that follow in a way a reader would notice.

**Constant preservation.** On a uniform ocean the horizon-area gap is exactly
zero -- not small, zero -- at every support and every horizon in the sweep.

**Containment.** The measured gap never exceeds the ocean area held in cells the
horizon actually cuts, at any horizon on any support. That bound is computed
from `lib/gridding.py:cell_quantiles` at probabilities 0 and 1, which are the
cell's shallowest and deepest region exactly, so the check is against the
distribution rather than against a second estimate of the mean.

## 2. The area on either side of a depth horizon

`depth_only` holds the cell's true ocean area fixed and moves only the test onto
the cell mean, so it isolates the threshold's own term; `binarised` adds the 0.5
wet mask the atmosphere applies, so the pair separates the two mechanisms
instead of reporting their sum.

| support | positive peak | at | negative peak | at | sign changes between |
| --- | ---: | ---: | ---: | ---: | --- |
| T21 | +45.4 | 1.00 km | -34.3 | 2.50 km | 1.75 and 2.00 km |
| T42 | +26.3 | 1.00 km | -19.9 | 2.25 km | 1.75 and 2.00 km |
| T85 | +18.0 | 1.25 km | -13.6 | 2.25 km | 1.75 and 2.00 km |
| T127 | +14.0 | 1.25 km | -10.6 | 2.50 km | 1.75 and 2.00 km |
| T170 | +12.0 | 1.25 km | -8.4 | 2.50 km | 1.75 and 2.00 km |
| GOLDSTEIN 36 x 36 | +49.2 | 1.00 km | -41.6 | 2.50 km | 1.75 and 2.00 km |

Multiples of the 0.001-of-planetary-area bar, on the `depth_only` arm. The
positive peak at T21 is 33.4 million square kilometres of ocean the cell-mean
support says is deeper than a kilometre and the terrain says is not.

**The sign is an outcome and it changes inside the sweep, exactly where the
pre-registration said it could.** Above 1.9 km the cell mean is shallower than
the horizon over most straddling cells and the coarse support loses deep water;
below it the mean is deeper and the coarse support invents it. No universal sign
was assumed and none exists: the crossing sits at the land-mean ocean depth,
which is where a mean has to change which side of a threshold it falls on.

**Refinement shrinks it and does not remove it.** The positive peak runs 45.4,
26.3, 18.0, 14.0, 12.0 times the bar from T21 to T170, a factor of 3.8 across
the whole ladder against a bar it still exceeds twelvefold at the finest rung.
The GOLDSTEIN candidate is worse than T21. A finer ocean grid is not a route out
of this one.

**The mask is a second mechanism and it does not shrink the same way.** At zero
horizon there is no threshold on depth at all, so the whole `binarised` gap is
the 0.5 wet cut: 3.3 to 5.0 times the bar across the ladder and 12.3 on the
GOLDSTEIN grid. `analysis/coastline_threshold_cost.json` already prices the LAND
side of that same cut; this is its ocean-side counterpart and nothing carried it
before.

## 3. Connectivity, and the one arm that does not converge

Deep water at a horizon is reported as the area NOT attached to the largest
connected body, because that is what a circulation cares about: a basin the
terrain isolates and the grid attaches is a basin the model ventilates when the
planet does not.

| support | 1.00 | 1.25 | 1.50 | 1.75 | 2.00 | 2.25 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T21 | 2.0 | 12.5 | 158.9 | 138.8 | 92.0 | 52.7 |
| T42 | 14.3 | 16.2 | 135.4 | 143.4 | 103.0 | 63.4 |
| T85 | 15.5 | 16.1 | 50.8 | 144.0 | 106.2 | 68.5 |
| T127 | 16.7 | 17.1 | 52.9 | 139.3 | 99.5 | 71.0 |
| T170 | 8.1 | 18.6 | 149.1 | 139.5 | 101.0 | 72.6 |
| GOLDSTEIN 36 x 36 | 2.8 | 6.8 | 166.7 | 138.9 | 81.1 | 43.1 |
| **the mesh** | **20.8** | **26.0** | **83.7** | **138.9** | **105.3** | **78.9** |

Detached area in millions of square kilometres, by horizon in km.

**The disagreement does not shrink with refinement and does not keep its sign.**
At 1.00 km every support understates the detached area, T170 by a factor of 2.6.
At 1.50 km every support except T85 and T127 overstates it, T170 by 65 million
square kilometres and the GOLDSTEIN grid by 83. T85 and T127 land near the
mesh's answer at 1.50 km and T170, the finest rung on the ladder, does not.

**That is not a support artifact and it is not a Jensen term, and the
pre-registration's dichotomy did not anticipate it.** P5 said a gap that fails
to shrink is carrying a length or a population bias. The mechanism here is a
third one: the quantity is discontinuous. Between 1.50 and 1.75 km the mesh's
own largest connected body falls from 161 to 54 million square kilometres and
the count of separately connected bodies above the area bar rises from 12 to 20.
The world ocean's deep water percolates apart in that interval, and an area
measured across a percolation transition moves by the size of a basin when the
support moves a passage by one cell. Refinement cannot make that continuous; it
only moves where the transition sits.

The consequence for the support is stronger than a convergence statement, not
weaker. A support cannot be chosen by taking the finest rung available, because
the finest rung is not the closest here. What a support has to reproduce is the
horizon at which the terrain's deep water separates and which bodies it
separates into, and that is a property the cell-mean bathymetry does not carry
at any resolution measured.

## 4. What the capacity bar says, and the one comparison it is waiting for

The loop's controlling scalar is `ocean/config/transport_loop.yaml`'s
`heat_transport_change`, 0.12 W m-2 over ocean area. On this world's ocean area
that is 0.050 PW of heat transport: a passage whose capacity the reduction moves
by a fraction `f` moves the heat it carries by `f` times its own transport, so
the material fractional error is 0.050 PW divided by that transport.

| the passage carries | a capacity error is material above |
| ---: | ---: |
| 0.1 PW | 50.0% |
| 0.25 PW | 20.0% |
| 0.5 PW | 10.0% |
| 1.0 PW | 5.0% |
| 2.0 PW | 2.5% |

The decision procedure is fixed here, before the ocean run that answers it
exists: take each passage's own heat transport from the first accepted ocean
run and compare the measured fractional capacity error against this table. That
converts the deferral into one comparison rather than a re-measurement, which is
the disposition the land audit's saturation-deficit arm took with the lapse
rate and for the same reason -- a bracket over an unknown circulation would be
wider than every effect it is meant to size.

## The pre-registration

Everything in this section was written and committed before the measurement ran.
The bars are the ones the consuming decisions already declare; none was chosen
after a result, and none is a percentage invented for this note.

### The bars, and what each one is an instrument for

| arm | bar | where the bar comes from |
| --- | --- | --- |
| ocean volume through the cell-mean depth | 1e-12 relative | `lib/remap.py:CLOSURE_TOLERANCE`, the float64 round-off floor this project closes every conservative crossing against. This arm is a CONTROL: the reduction is exact and the bar exists to catch an implementation error, not a Jensen term |
| ocean area on either side of a depth horizon | 0.001 of planetary area | `ocean/config/transport_loop.yaml`, criterion `sea_ice_area`. It is the project's one declared bar on a planetary AREA FRACTION, set for the area whose climate role -- an albedo and insulation contrast switched by a threshold on a continuous field -- is the closest analogue a shelf/deep partition has |
| deep-water connectivity | 0.001 of planetary area, categorical in form | the same bar, applied to the area of water whose connectivity CLASS differs between the two orders. The form is `transport_loop.yaml`'s `carve_drivers` criterion: a topology change is counted, not averaged |
| passage transport capacity | 0.12 W m-2 of ocean-area-mean heat transport | `ocean/config/transport_loop.yaml`, criterion `heat_transport_change`, which is the loop's own controlling scalar and the same 0.12 W m-2 storage tolerance `config/partial_surface.yaml` selected the tile operator against |

The last bar is stated in the arm's own units at the point of reporting. A
passage carrying an ocean heat transport `H` and whose capacity the reduction
moves by a fraction `f` moves the heat by `f * H`, and the bar is reached at

    f * H = 0.12 W m-2 * A_ocean

so the material fractional capacity error is `0.12 * A_ocean / H`. `A_ocean` is
measured here; `H` is not knowable without an ocean run, so the arm publishes
the critical fraction against `H` and the first accepted ocean run answers it in
one comparison rather than a re-measurement. That is the same disposition the
land audit's saturation-deficit arm took with the lapse rate, and for the same
reason: a bracket over an unknown circulation would be wider than the effect it
is meant to size.

### The invariants, with no universal sign for refinement

**P1. Ocean volume is affine in depth, so the cell-mean reduction is EXACT for
it at every support.** A nonzero volume gap indicts the binning or the areas,
never the physics. FALSIFIABLE, and it is the arm that can fail.

**P2. Constant preservation.** A uniform depth field reduces to itself in every
covered cell, so the area-on-either-side-of-a-horizon gap is exactly zero on a
uniform ocean at every support and every horizon. This is `cell_mean`'s own
invariant exercised through this arm rather than asserted about it.

**P3. The horizon-area gap is confined to STRADDLING cells.** A cell whose
whole fine depth distribution lies on one side of the horizon gives the same
answer in both orders, so the gap is bounded above by the ocean area held in
cells the horizon cuts. A gap exceeding that bound is a binning error.

**P4. No sign is assumed, in either arm, and the two arms have opposite
mechanisms.** For the horizon area the gap is one sign where the cell mean is
deeper than the horizon and the other where it is shallower, so the global sign
is an outcome and not a prediction. For connectivity, coarsening acts in two
directions at once: a cell mean is at least the cell minimum, which DEEPENS
bottlenecks and opens passages, while a passage narrower than a cell can fall
below the 0.5 ocean-fraction the atmosphere binarises at, which CLOSES it. Which
dominates, and at what depth, is the measurement.

**P5. Refinement.** In the limit of one mesh region per cell the two orders are
the same computation, so every gap here must approach zero along the
T21/T42/T85/T127/T170 ladder. A gap that does not shrink is carrying a support
artifact -- a length or a population bias, `subgrid_slope`'s mechanism in the
land audit's section 6 -- and not a Jensen term. This is the check that
separates the two, and it is pre-registered because a gap that grows with
refinement would otherwise be read as evidence about resolution.

### What the arms do NOT cover, declared before the measurement

The delivery law itself. Section 8d's terminal export has no marine side: ANUT-6
routes to the coast and OCN-13 begins there, and neither exists. What is
measured here is the geometry that any such law would read -- shelf area,
coastal ocean area, connectivity -- and not the delivery. An arm sizing the
nonlinearity of a law that has not been written would be sizing a preference.

## What the ocean support has to carry, and why a better mean is not the answer

The three quantities this audit asks of one cell's bathymetry behave differently
and the difference is structural rather than a matter of accuracy.

**Volume is affine in depth.** The area-weighted mean reproduces it exactly, at
every support, and the control arm above is that identity rather than a
comparison.

**The area on either side of a depth horizon is a threshold on the
distribution.** A mean is one number and a threshold reads a distribution, so no
single depth per cell can carry it. That is not a statement about the mean being
a poor estimator: it is that the two quantities are not functions of each other.
Choosing the median, or the depth that puts the cell on the right model level,
or the deepest point, moves which horizon the cell is right at and leaves it
wrong at every other one.

**Connectivity is a property of the graph the cells form** and is not a
per-cell quantity at all. It depends on which cells are wet, which is the
binarisation, and on which are deep enough at a horizon, which is the threshold
above. Both of its inputs are the two the mean does not carry.

So the support's disposition is not a choice between reductions. It is a choice
between ONE NUMBER and a DISTRIBUTION, and the operator for the second already
exists: `lib/gridding.py:cell_quantiles` is the area-weighted quantile table a
sub-grid hypsometry is, and `area_fraction_above` reads the share above a
horizon back out of it. `hydrography/scripts/build_spatial_support.py` already
carries a per-cell hypsometry on the LAND support for the same reason. The ocean
support does not exist yet, so nothing has defaulted; what this audit says is
what it must not default to.

**The pipeline's gate as written cannot catch this.** `config/pipeline.yaml`'s
`ocean_support` step requires that bathymetry "preserve extensive quantities on
their declared wet support". Volume is the extensive quantity, and the cell mean
preserves it exactly, so a support built entirely on cell-mean depths passes
that gate and fails every threshold a consumer reads from it. The gate needs the
distribution stated as a requirement rather than implied by conservation.


## What follows

- The ocean support must carry a per-cell depth DISTRIBUTION and not a cell
  depth. `lib/gridding.py:cell_quantiles` is the operator and OCN-11 is the
  step; nothing has defaulted yet, because the support does not exist.
- `config/pipeline.yaml`'s `ocean_support` gate asks for extensive preservation,
  which the cell mean already gives exactly. It now asks for the distribution as
  well, because the arm that could have failed did not and the arms that matter
  are the ones conservation cannot see.
- A support cannot be chosen by taking the finest rung. Section 3 measures the
  finest rung on the ladder landing further from the mesh than two coarser ones,
  because the quantity crosses a percolation transition. What a candidate has to
  reproduce is the horizon at which this world's deep water separates and the
  bodies it separates into.
- The 0.5 wet mask's own ocean-side area error is 3.3 to 12.3 times the same
  bar. It is a mask question rather than a bathymetry one, so it is carried
  beside SPAT-5's land side in `analysis/coastline_threshold_cost.json`, where
  every candidate rule is now priced on both sides of the coastline and
  `config/partial_surface.yaml`'s retain verdict names both halves of its cost.
  The Gaussian ladder's figures there reproduce this audit's exactly, computed
  from the land ledger's two one-signed halves rather than by binning the mesh,
  which is a second route to the same number. The 36 x 36 GOLDSTEIN figure is
  that grid's own binarisation and stays here: it is OCN-11's mask and not the
  one ExoPlaSim applies.
- The two sides do not select the same rule. The area-conserving threshold is
  the only one of the five under the planetary-area bar on the ocean side, at
  0.003 to 0.09 times it against the retained rule's 3.3 to 5.0, and it is worse
  on the below-datum land-volume closure the retain verdict rests on. The tile
  model is selected, so the binary mask's remaining job is topology, and a
  topology criterion is what section 3 measures the support against.
- The passage-capacity arm is not deferred, it is fixed: the critical fractional
  error against a passage's own heat transport is tabulated above, and the first
  accepted ocean run answers it in one comparison.
