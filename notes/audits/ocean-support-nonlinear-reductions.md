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
