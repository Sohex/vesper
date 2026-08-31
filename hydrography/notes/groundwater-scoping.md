# Groundwater: what a water table would change, and what it would not

This is worldbuilding. Vesper is a fictional super-Earth and everything below is
about a simulation of it: a terrain generator, a toy climate model, and the
water balance this project computes between them. No sentence here is about the
real world except where it names a source measured on Earth.

Scoping for a steady-state water table on the region mesh. It answers the
question that comes before building one: what does the pipeline get wrong today
that a water table would get right, and what would it still get wrong.

## 1. What the pipeline assumes now

Every water balance in this project is a surface balance. `surface_water.py`
takes runoff as P-E over land, accumulates it down the drainage network, and
solves lake area against open-water evaporation. `carve_verdict.py` integrates
the same fields over surface catchments. ExoPlaSim's land surface underneath
both is a single bucket: `landmod.f90:1098` makes runoff literally the overflow
above `dwmax`, with no water table and no drainage flux.

Fan (2019) writes the catchment budget as

    dS/dt = P - ET - Qr - Qg

for stream outflow `Qr` and net groundwater outflow `Qg`, and says of the last
term that "a common assumption is that the Qg term is negligible". That is the
assumption this pipeline makes, and it makes it everywhere, without stating it.

## 2. The term is not negligible in the regime this world occupies

This is the part that decides whether the work is worth doing, and it does not
rest on any Vesper measurement.

Fan (2019) section 3.1 puts the conditions for a leaky catchment as small
catchment size, position at either end of a steep regional gradient, deep
permeable substrate extending beyond the catchment, and **drier climate**. Her
Figure 2 is the mechanism: under a humid case at 381 mm/yr of recharge the water
table stays high enough to intercept local stream beds, local flow systems close,
and groundwater divides form near the topographic divides, so catchments are
self-contained. Under an arid case at 25 mm/yr the water table falls below the
stream beds in the headwaters, recharge joins the regional flow system and
resurfaces in the lower basin, and **the groundwater divides no longer exist**.
Higher substrate permeability widens the same gap.

So the surface catchment stops being the water catchment precisely in the dry,
high-relief, closed-basin setting that this world is mostly made of, and
precisely for the basins the carve verdict decides. The assumption is weakest
where it is load-bearing.

Fan's Hypothesis 1 bounds it from the other side and is the reason not to
overclaim: aggregate the catchments along a regional gradient into one larger
catchment and the export from the upper cancels the import into the lower, so
the large catchment is self-contained again. At steady state groundwater
redistributes water between basins; it does not create or destroy it. A basin
whose groundwater divide follows its surface divide gains nothing from this
model at all.

That leaves four channels through which a water table changes a Vesper number,
and they are the whole of it:

1. Groundwater divides that do not follow surface divides, and regional flow
   that crosses them. Fan's export/import pair.
2. Evaporation from a shallow water table in discharge zones, which is a sink
   the single bucket cannot represent because `drhs` throttles on
   `dwatc/dwmax` and the cell dries out.
3. Leakage from lakes through permeable substrate.
4. Baseflow, which redistributes discharge in time and changes no annual total.

Channels 1 and 2 have opposite signs at a given basin. Which one wins there is
not knowable before the solve, and that is the honest statement of what this
model predicts.

## 3. The inputs all already exist

A steady-state water table needs terrain, recharge and permeability.

**Terrain** is the region mesh, with the depression-filled surface and the
`receiver` network `build_hydrography.py` already resolves.

**Recharge** is the field `surface_water.py` already builds from the baseline
climatology, read through the coupling matrix. It crosses the export/ExoPlaSim
grid boundary, so it is subject to CLAUDE.md rule 3 and to the half-cell
longitude offset in hydrography's known approximations; `basin_means` already
handles both and is the route to reuse.

**Permeability** is the one that would normally have to be invented, and does
not. Gleeson et al. (2011) Table 1 gives geometric-mean log permeability for
five combined hydrolithologies keyed on the Duerr et al. (2005) lithologic
classes:

| hydrolithology | log k (m2) | sigma | n |
| --- | ---: | ---: | ---: |
| carbonate | -11.8 | 1.5 | 47 |
| volcanic | -12.5 | 1.8 | 33 |
| unconsolidated | -13.0 | 2.0 | 113 |
| crystalline | -14.1 | 1.5 | 17 |
| siliciclastic sedimentary | -15.2 | 2.5 | 20 |

Those are the same classes `pedology/config/pedogenesis.yaml` already maps
Orogen's lithology onto for Hartmann et al. (2014) phosphorus: PA, PB, VA, VB,
MT, SS, SC, SU, EP, with Moosdorf and Hartmann as shared authors on both
compilations. The permeability table therefore drops onto an existing mapping
rather than requiring a new one, on the precedent that mapping already set.

Two things about that table are load-bearing:

**It is scale-independent over the range this mesh uses.** Gleeson section 3
paragraph 8 tests permeability against hydrolithologic unit length and finds no
dependence for crystalline or coarse-grained unconsolidated units, concluding
the geometric mean "can represent large areas (5-100 km in length)". A region
here is 230 km2, about 15 km across, which is inside that range. The single
exception Gleeson names is carbonate, where permeability does increase with
scale through karst, so carbonate terrain is a one-signed underestimate here.

**Evaporite has no assigned permeability.** Duerr's EV class is in Gleeson's
"not assigned" row along with water bodies and ice. This world has a great deal
of playa and salt crust, so that is a real hole and it is over the surfaces the
brine and duricrust rules care about most. It is a missing number, not missing
machinery, and it must be carried as unassigned rather than filled by the
nearest class.

The class spread is 3.4 orders of magnitude from carbonate to siliciclastic,
and the within-class sigma is 1.5 to 2.5 orders. So a water table depth from
this table is a bracket, not a value, and the sigma column supplies the bracket
directly.

## 4. Gravity carries into hydraulic conductivity, and permeability does not

Permeability `k` is a property of the pore geometry. Hydraulic conductivity is
`K = k rho g / mu`, so Gleeson's `k` transfers to Vesper unchanged while `K`
scales with gravity: 12.81 against 9.81 m/s2 makes conductivity 1.31 times
Earth's for the same rock. Higher conductivity drains the same recharge under a
smaller head gradient, so the modelled water table is flatter, deeper under
ridges and more concentrated at discharge zones than Earth's would be under the
same climate and rock.

This is the pattern `notes/audits/orogen-gravity.md` and the DUST-6 finding
exist to catch, it is free to get right at the start, and it is expensive to
retrofit once a calibration has absorbed it.

## 5. The resolution limit, stated before anything is built

Fan et al. (2013) separate three scales of control, and they do not survive the
mesh equally.

**Regional recharge control is resolvable.** "Regions of deep WTD correspond to
regions of low recharge, the great deserts of the world stand out", against
shallow tables under the high-recharge tropical swamps. Vesper's aridity
gradients are far wider than a region.

**Basin-scale convergence is resolvable, and is the channel that matters here.**
Fan's own example is "arid basins where groundwater convergence from surrounding
mountains maintains valley ecosystems (oases) otherwise absent", which she calls
the larger examples of lateral convergence. That is a whole closed basin, which
is the object this component already catalogues.

**Valley-scale texture is not.** Fan puts the well-articulated gradient "from
valley to ridge spanning decameters to kilometres", and says that at local
scales terrain signals dominate and override climate boundaries. At 15.19 km a
region cannot hold a valley-to-ridge transect, so the texture Fan measures as
locally dominant is sub-grid here and stays sub-grid. This is the same argument
GRAV-6 makes for glacial valleys, and it has the same conclusion: the process is
real, it is not resolvable at this mesh, and claiming it would be precision
theatre.

Sea level is Fan's dominant global driver, and Vesper has one, so the coastal
band of shallow water table comes free and correct.

## 6. What it unblocks downstream

**Pedology has a field declared empty and waiting for exactly this.**
`pedology/config/surface_classes.yaml` refuses to place Ullyott and Nash's
groundwater silcrete, with the reason written in the config: it sits at or near
a water table and this project models none. Fenske et al. (2025) put calcrete
formation with it, "by the precipitation of dissolved groundwater calcite under
dry conditions". Both become placeable from a mean depth field and a discharge
mask.

The rest of Fenske does not become reachable, and the reason is a project rule
rather than an effort budget. Their model hardens a layer at a depth set by the
RANGE of water table fluctuation over a characteristic timescale, and they put
duricrust formation time at 10^5 years or longer, citing Tardy for crusts
"mostly monogenic, at least millions if not tens of millions of years old".
`docs/src/reference/no-time-axis.md` forbids asking any component here for that
duration. A steady-state solve gives a mean depth and no fluctuation range at
all. The nearest reachable substitute is the range over the stellar cycle, which
is the same machinery SURF-4 already asks for on lake occupancy, and it is a
57-year band rather than Fenske's timescale, so it would be a different quantity
wearing the same name. Note also that Fenske puts gypcrete outside this
mechanism entirely, as an evaporitic blanket formed by surface and air processes
in hyper-arid settings, so the existing placement of that one is untouched.

**Minerals had ONE rule keyed on a water table through a rainfall proxy, and it
was deleted rather than replaced.** Supergene copper carried
`maximum_precipitation_mm_yr: 500.0` with the reasoning that copper is flushed
once percolation runs year round and the admission, in the same comment, that a
thick leached cap needs a deep water table rather than a dry climate as such.
MIN-4 removed it on 2026-08-18 because Sillitoe (2005) p. 736 states there is no
wet bound on rainfall at all, so what stood there was contradicted by the
canonical source and not merely unsourced.
`minerals/config/downstream_prospectivity.yaml` now carries no rainfall term
standing in for a depth. Bauxite's and nickel laterite's precipitation
thresholds are their own sources' rainfall criteria, applied by Price et al.
(1997) and Butt and Cluzel (2013) to gridded climate fields as rainfall, so they
are not proxies either.

**What a depth field would still have to clear there is RANGE, not skill.**
Reich and Vasconcelos (2015) p. 306 put leached caps at several hundred metres
"particularly when the water table was deep enough", which is the depth the
mechanism is keyed on; `groundwater-et-sink.md` measures the modelled 95th
percentile against it. A depth field also does not supply the descent RATE,
which is a duration and is refused by the same rule as above.

**Hydrography gets the `Qg` term** for the carve criterion, per basin. This is
the largest consequence and the one to be slowest about, because the carve
verdict feeds loop A.

## 7. What it does not reach

**The relict half of HYD-4.** `hydrography/notes/lake-solver-validation.md`
found the dry runoff quartile under-predicting observed lake area about eightfold
and named two missing terms: groundwater, and Earth's arid terminal lakes being
substantially relicts of wetter climates rather than equilibria. A water table
reaches the first and cannot reach the second, which is a disequilibrium history
and therefore a duration.

**And it must not be tuned against that number.** The case for the term is that
Fan measures it on Earth and derives it from theory, not that adding it would
close a gap in this project's own comparison. A correct term that worsens an
agreement is information; `docs/src/practice/failure-modes.md` class 16.
Channels 1 and 2 of section 2 pull in opposite directions, so the sign of the
effect on lake area is not decided in advance, and any version of this model
whose parameters were chosen after seeing that comparison is not a model.

**The climate feedback.** Groundwater-supported evapotranspiration in drylands
is real and first-order, and reaching it means either a water table lower
boundary in `landmod.f90`, which drags CLAUDE.md rule 4 and a full binary
rebuild behind it, or a fork of the LPJ-GUESS soil column. Neither belongs in a
first pass. Raising `dwmax` in convergence zones through the existing
`build_surface_soil_water.py` step looks like a cheap substitute and is not one:
a deeper bucket changes storage and timing, not the availability floor a water
table sets, so it would move the answer without representing the mechanism.

## 7b. GW-9 settled: the e-folding length does not survive a 15 km cell

Measured and argued 2026-08-20, after both solver schemes failed to converge.

The transmissivity is `T = A(x) exp(h/f)` for an e-folding length `f` set by
Fan's slope curve, `f = a / (1 + b * slope)`, with the fork's constants giving a
bedrock floor of `20 / (1 + 125 * 0.16) = 0.95 m` at the slope cap. **Two
separate things are wrong with that here, and they are not the same defect.**

**The slope is the wrong slope.** Fan fits that curve against a hillslope
gradient resolved at about 1 km. What this project has is `local_slope_deg`, a
plane fit through each region and its neighbours at 15.19 km, and
`minerals/README.md` documents it as an upper bound on REGIONAL DIP, which is the
quantity it was built for. Terrain gradient falls with the scale it is measured
over, so a regional dip understates the hillslope gradient at the same place, and
feeding it to Fan's curve returns an `f` that is too large. One-signed, toward a
deeper and flatter water table.

**And the floor is incoherent at any slope.** This is the one that broke the
solver, and it is a modelling failure rather than a numerical one. Take the cell
mean of the transmissivity when the water table varies within the cell with
spread `sigma`:

    E[exp(h/f)] = exp(hbar/f) * exp(sigma^2 / 2 f^2)

With sub-grid relief of order 100 m and `f` of 1 m, that factor is `exp(5000)`.
It is not a correction to be applied; it is the statement that **the
parameterisation cannot be evaluated at a cell mean head at all** when `f` is far
below the sub-grid relief. The solver's measured `|z|/f` reaching 4588, against
an `exp` argument limit near 709 in double precision, is the same fact arriving
as an overflow.

Setting `f` of order `sigma` makes that factor `exp(0.5)`, which is an ordinary
order-one term. So the requirement is not a tuned value but a scale:

**The effective e-folding length for a cell cannot be smaller than the relief
that cell contains.** `f` is a soil and weathering profile depth, a genuine
1 km-scale quantity; a 15.19 km control volume averaging over hundreds of metres
of relief has no business carrying it. Fan's own numbers are not wrong, they are
being asked a question at three orders of magnitude from where they were fitted.

### What this means for GW-1

The blocker is not the iteration scheme. Picard limit-cycled and Kirchhoff
overflowed, and both are symptoms of a stiffness that a coherent `f` would not
produce. Newton on the full residual might well converge -- **and it would be
converging a formulation that applies a metre-scale parameter at 15 km, which is
worse than not converging, because the answer would look like a result.**

The sub-grid relief per cell is exactly what GRAV-6 already parks glacial erosion
behind: its recommendation is that the process "belongs with the downscaling
machinery, which has to persist sub-grid hypsometry anyway". The same quantity
settles this. That is not a coincidence, and it is not a reason to defer this one
alongside it: a water table reaches the carve verdict, the surface classes and
the mineral rules, all of which sit in loops that run before any downscaling.

## 8. Where it belongs

In `hydrography/`, as a peer of `surface_water.py`, not a sibling component. It
solves on the same mesh, over the same `receiver` network and the same basin
catalogue, reads recharge through the same coupling matrices, and would need a
cross-component provenance stamp under rule 5 for data that is one solve away if
it sat anywhere else. `surface_water.py` is already the climate-dependent member
of a directory whose README calls the rest of itself climate-independent, so a
climate-dependent groundwater solve is not an anomaly there.

It is therefore a loop A step, downstream of a climatology and upstream of
`carve_verdict`, and it needs its row in `config/pipeline.yaml` in the commit
that creates it.

## 9. The tests, declared before the code

**The reduction identity, which is the strong one.** Drive permeability to zero
and the solver must reproduce the current surface-only lake solution exactly,
not approximately. Every cell's recharge stays where it fell, `Qg` is zero for
every basin, and `surface_water.nc` comes back bit-comparable.

**Closure.** At convergence, total recharge equals total discharge across
seepage, baseflow to the network and lake leakage, to machine precision. This is
a conservation law and it either holds or the discretisation is wrong.

**The divide test.** Under spatially uniform permeability and a terrain-following
water table, groundwater catchments must reproduce surface catchments. If they
do not, the solver is wrong and not the world.

**Earth.** The same code on Earth topography, Earth recharge and GLHYMPS
permeability, scored against Fan et al. (2013)'s compilation of 1,603,781 well
sites. The threshold is declared before the data is scored, per the standing
convention, and the scoring must account for a sampling bias Fan states
explicitly: observations favour valleys and oases, wells are not sited at
random, and her own model reads deeper than the observations in arid regions for
that reason. A comparison that ignores it will fail the model for being right.
