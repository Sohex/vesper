# Hydrography

Turns the World Orogen geography into drainage: where water goes, which basins
it collects in, and how much each can hold. Everything here except the lake
solver's forcing is climate-independent and can be built before ExoPlaSim runs.

```bash
python hydrography/scripts/build_hydrography.py   # ~13 s, climate-independent
python hydrography/scripts/lake_balance.py        # solver smoke test and sweep
python hydrography/scripts/surface_water.py       # ~3 s, needs a climatology
python hydrography/scripts/groundwater.py         # the discretisation checks
python hydrography/scripts/build_groundwater.py   # the water table, needs a climatology
python hydrography/scripts/land_water_ledger.py   # the store and flux ownership contract
python hydrography/scripts/build_groundwater_access.py   # that water table on the climate grid
```

## Why this component exists

Upstream World Orogen forced every land cell to drain to the sea. The fork stops
doing that, so closed basins survive into the export, and it deliberately does
not decide what happens next: which basins hold lakes is a precipitation-versus-
evaporation balance, and that belongs downstream.

What arrives is therefore unrouted. `drain_to` is raw steepest descent, and
`drainage_terminal` is -2 for anything whose chain ends in neither the ocean nor
a preserved basin sink. On this planet that is **most of the land**, draining
into unpreserved pits, most of them a single mesh cell of noise -- the exact
share is per build and is in `world_state.json`. Integrating
precipitation over catchments without resolving that would silently discard most
of the land's water.

## What it does

**Resolves drainage** with a priority flood over the 2.5M-region mesh, filling
the noise pits while keeping the surviving basins as genuine terminals. How many
survive is a property of the build and is in `world_state.json`; it named a
count and a build here, and both were several iterations out of date.
Every land region ends up assigned to the world ocean or to exactly one basin.
The fill depth has a median of a few metres, which is the scale that confirms
these were noise rather than landforms; the counts are in the hydrography
report for the build they were computed on.

**Recomputes hypsometry on the finished terrain.** The catalogue's
`hypsometry` is measured on the natural, pre-conditioning surface and overstates
what the basins can hold: the finished terrain holds substantially less of the
catalogue's
figure. Level, area and volume curves are rebuilt from the flood
and are exact at their sample points, verified against brute force to 4e-15
relative on volume and exactly on area.

**Finds spill levels and targets** on the finished terrain. Most basins overflow
into another basin rather than to the ocean, so filling has to be solved as a
cascade rather than basin by basin.

**Couples to the climate grid.** `coupling_*.nc` gives each basin's catchment
area per grid cell as a sparse matrix, which is what you integrate P-E over.
Built for T42 and T85.

## Products

Everything under `data/` is **per build**, written to `data/<source_build>/`.
The paths below are shown relative to that directory. There is no flat
fallback and `component_data` is strict about it: a flat `data/basins.nc`
would hold whichever build was active when it was last written, and pairing
one terrain's rows with another's columns is a failure this project has
already had.

| File | Contents |
| --- | --- |
| `data/<build>/regions.nc` | per mesh region: resolved terminal, depression-filled surface |
| `data/<build>/basins.nc` | per basin: hypsometric curves, spill level and target, catchment area, capacity |
| `data/<build>/coupling_<grid>.nc` | sparse basin-by-grid-cell catchment areas |
| `data/<build>/hydrography_report.json` | diagnostics, river mouths, marginal seas, provenance |
| `data/<build>/surface_water.nc` | per region: lake, lake depth, river discharge; per basin: solved area, level, volume, overflow |
| `data/<build>/topographic_index_<grid>.nc` | per region: the compound topographic index on both slope arms; per grid cell: its mean, its within-cell spread, and `f_sat_max` |
| `analysis/lake_balance_sweep.json` | solver sensitivity under placeholder forcing |
| `analysis/surface_water_report.json` | the solved water balance and its forcing |
| `analysis/topographic_index_report.json` | the index's distribution, the scale measurement, and the score declared before it |
| `analysis/land_water_ledger_report.json` | the land water ledger's graph, its closure fixtures, its graph mutations, and every store that still has a shadow copy |

## The land water ledger

`land_water_ledger.py` is a contract, not a solver. It reads
`config/land_water_ledger.yaml` and answers one question about any proposed
coupling: which store is being debited, by which component, over which
interval. It refuses the answer "all of them".

The pipeline already has a water balance and it is not the same thing. `P - E`
closes at equilibrium, the groundwater solver conserves its own recharge
against its own seepage, and the basin code adds only signed exchange to
catchment supply. Each is right on its own terms, and together they let one
interval's positive `P - E` be spent three times: as ExoPlaSim's bucket
overflow, as this component's catchment runoff, and as the groundwater
builder's recharge.

The ledger books water mass in kilograms, absolute, with one source and one
destination per term. Its invariant is that water held in reservoirs plus water
at terminal nodes equals water held at the start plus water across boundary
nodes, and a second invariant beside it says each evaporative component is
debited by exactly one term from exactly one store -- which the mass identity
does not imply, because two models withdrawing the same millimetre from the
same store conserve mass and are still wrong.

Both arms are exercised on every invocation: eleven closure fixtures, ten built
to be wrong in a named way and required to be rejected, and twelve mutations of
the declaration itself that the graph check is required to catch. A fixture or
a mutation that does not get its verdict is a defect in the checker and exits
non-zero.

Every term carries the `undeclared` sentinel today, so the ledger is DEFINED
and does not CLOSE. `--strict` is the arm that refuses, and it refuses while any
term is undeclared, while any store still has a shadow copy in a second
component, or while any term's producer supplies only an annual mean to a store
that changes inside the year. `notes/land-water-ledger.md` is the contract and
the argument.

## The lake solver

`lake_balance.py` is machinery, not a result. At equilibrium

    r (C - A) + p A = e A      so      A = r C / (e - p + r)

for catchment area C, lake area A, runoff depth r over dry land, and
precipitation p and evaporation e over open water. Overflow cascades downstream
and the system iterates to a fixed point.

Running it directly sweeps uniform placeholder forcing, which exercises the
solver and shows the sensitivity. Under 10 to 200 mm/yr of runoff against 400 to
1600 mm/yr of lake evaporation, lake area lands between 0.21% and 8.7% of the
planet. That range is a property of the terrain, not a prediction.

## Surface water under a real climate

`surface_water.py` supplies the forcing the solver was waiting for, from the
`climatology_s096` baseline, and accumulates the same water down the drainage
network to get rivers. This is the first thing in the project to decide
`surface_class == 2`, which World Orogen deliberately leaves empty.

Current figures -- how many basins hold water, lake area, largest river -- are in
`world_state.json` and in the run's own report, not here. They change with every
iteration and every baseline re-run, and this document is about how the
solver works.

Against Earth, which is the only calibration available, the ratios are the part
worth carrying: this world has roughly twice the land at somewhat lower runoff
depth, so total river discharge is of the same order, and the largest river is
comparable to the Amazon. Lakes and endorheic drainage are both several times
Earth's share, and by a similar factor -- which is the consistency check, since
they are two views of the same aridity.
Those last two are computed by different routes and agreeing is a check rather
than a coincidence.

Two choices in the forcing are worth stating, because neither is the obvious
field:

**Runoff is P-E, not `mrro`.** At steady state they are the same quantity:
whatever falls on land and does not evaporate has to leave. The baseline run's
global water budget closes to a small fraction of a percent and the land's
surplus matches the sea's deficit, so the budget is trustworthy. Its `mrro`
diagnostic accounts for a fraction of that surplus, because it is river-routed
net divergence rather than local generation, so the diagnostic is the unreliable
one and P-E is what a water
balance can be built on. Both are recorded in the report.

**Lake evaporation is Penman, shared with the carve verdict.** Over land the
model's `evap` is actual evapotranspiration, throttled by soil moisture, which
is the wrong quantity for open water: using it would make every lake as dry as
the desert it sits in. `carve_verdict.penman_open_water` answers the right
question instead, and using the same function here means the lakes and the
terrain they sit in are judged by one rule.

It is also the estimate that can be checked. Applied to ocean cells, which
already are open water and which the model gives water's own roughness, Penman
reproduces the model's own evaporation to within a few percent, against 7.8% for
a Priestley-Taylor estimate at alpha 1.26. The run's own ratio is computed by
`validate_over_ocean` and written into `carve_verdict.json`; it is not repeated
here, because every prose copy of it in this repository was a copy of a hardcoded
constant.


**Penman is NOT floored at the model's land rate.** It was, on the reasoning that
a saturated surface cannot evaporate less than the moisture-limited ground beside
it. That is false here: land carries a roughness field with a median `z0` of
0.521 m against open water's 1.5e-4, which is a transfer coefficient 6.4 times
larger, so a smooth lake in a rough wet landscape genuinely evaporates less than
the land around it. The floor bound on 60.2% of land cells and decided 73% of the
overflowing basins by clamp; `notes/audits/carve-criterion-terms.md` finding 1.

One caveat stands about the validation itself: the ocean is where air is
near-saturated and wind is well resolved, which is the opposite of an inland arid
basin, so the measured ratio is an upper bound on the estimate's accuracy in the
places that decide the verdict.

### Basin coupling and longitude conventions

The coupling matrix inherits the Orogen grid, columns -180 to 180; an
ExoPlaSim climatology numbers its own from 0 to 360, so indexing one with the
other's numbering reads each basin's climate from near its antipode while
shapes, latitudes and plausibility all survive. `coupling_*.nc` carries the
`cell_lat` and `cell_lon` it was built on, `basin_means` takes the longitude
axis of the fields it is given and remaps columns before indexing, and it
refuses a coupling file too old to state its own convention. The check: ask
the coupling to average a field of known longitude and compare the answer to
where the basins are. It is cheap; rerun it whenever either grid changes.

Lakes are painted by filling each basin to its solved *area*, in ascending order
of the flooded surface, rather than by thresholding on the solved level. A mesh
region is 230 km2, so thresholding gives every basin with a smaller lake a 15 km
blob, and a wet landscape of many small basins comes out speckled with lakes
that are mostly rounding error. Filling by area reproduces the solved total to
0.4% and simply does not draw a lake too small to reach a whole region.

Rivers come from the priority flood's own discovery pointer, now kept as
`receiver` in `regions.nc`. Steepest descent on the filled surface would not do:
a filled depression is flat and has no downhill neighbour, so whole tributaries
would be dropped. Accumulation peels leaves off that tree rather than sorting on
elevation, for the same reason.

**Overflow is routed onto the mesh**, which needs the saddle and not just its
height. `spill_levels` was taking the minimum over a basin's exit edges and
keeping only the value; it now keeps the argument too, so `basins.nc` carries
`spill_region` and `spill_exit_region` either side of the saddle. Without them a
basin pinned at its spill had a known outflow and nowhere to start it, and the
largest flows on the planet were missing from the network: routing the
overflowing basins multiplied the biggest river's discharge severalfold.

There is no double counting, because the solver has already resolved the
cascade: a basin's overflow is its final equilibrium value with everything
upstream included, and the paths are disjoint segments, one basin's saddle to
the next one's sink.


## The water table

`build_groundwater.py` solves a steady-state, vertically integrated, unconfined
water table on the region mesh and reports the term every other balance in this
component drops. `groundwater.py` is the machinery: the Voronoi discretisation,
the transmissivity, and the complementarity solve.

Fan (2019) writes the catchment budget as `dS/dt = P - ET - Qr - Qg` and notes
that dropping the net groundwater term `Qg` is the common assumption. It is the
one made everywhere here, and Fan's own Figure 2 puts its failure in this
world's regime: at high recharge the water table intercepts the stream beds,
local flow systems close and the groundwater divides sit near the topographic
ones; at low recharge the table drops below the beds and the groundwater divides
stop existing. Dry, high-relief and closed-basin is where the assumption is
weakest and where the carve verdict is decided.

`hydrography/notes/groundwater-scoping.md` is the argument: what a water table
changes, what it does not reach, and the four channels through which the term is
non-zero at all. Two of them pull in opposite directions at a given basin, so
this model does not predict the sign of its own effect on lake area, and it must
not be tuned against one.

**The steady state is not a simplification chosen for cost.** An equilibrium
water table is the form `docs/src/reference/no-time-axis.md` leaves reachable.
An aquifer with a residence time is not, and neither is Fenske et al. (2025)'s
duricrust model, which hardens a layer over the RANGE of water table
fluctuation across 10^5 years or more.

**What it reads.** Terrain and lithology from the export, recharge as P-E
through the same coupling convention `surface_water.py` uses, and
`config/groundwater.yaml` for the subsurface: Gleeson et al. (2011) permeability
per hydrolithology, on the same Duerr class mapping pedology already uses for
Hartmann phosphorus, over a constant saturated thickness.

**Transmissivity is `T = K D`, and Fan's exponential decay with depth is NOT
used.** It was, and it could not be evaluated on this mesh: the cell mean of
`T = A exp(h/f)` carries a factor `exp(sigma^2/2f^2)`, which at 100 m of
sub-grid relief against a metre-scale e-folding length is `exp(5000)`. That is
not a large correction, it is a statement that the parameterisation has no
cell-mean value at 15.19 km. `D` of order 100 m is the depth Gleeson's own maps
represent, so conductivity and the thickness it multiplies are one paper's
statement about one rock at one scale. GW-9.

**Two terms exist and are OFF unless asked for**, each bit-identical to a run
without it:

- `--et-lambda` gives the water table a sink, groundwater evapotranspiration
  exponential in depth on Shah et al. (2007), whose finding is that the decline
  with depth is exponential rather than the linear form MODFLOW's EVT uses.
  `ET_max` is the Penman field the lakes and the carve verdict already take.
  Without it the table pins at the surface over most of the land, because the
  only exits are the coast and seepage and recharge in a continental interior
  reaches neither.
- `fixed_head_m` imposes a head on river and lake cells, so a water surface IS
  the water table there. Without local baselevels the sea is the only fixed head
  and no defensible transmissivity reaches it.

**Gravity enters once, and correctly.** Permeability is pore geometry and
carries over from Earth unchanged; hydraulic conductivity is `k rho g / mu` and
is this world's. The config tabulates permeability and never conductivity, and
the solver takes gravity from `config/planet.yaml`.

## What this field may be used for, and what it may not

**It is not a wetland mask.** Wetland extent depends on native-mesh relief,
seasonal persistence, routed inundation and groundwater discharge; thresholding
an ExoPlaSim/LPJ cell-mean depth would erase all four and can double-count the
solved lakes, rivers and playas. WET-2 therefore derives mutually exclusive
native-mesh persistent-peat, saturated-soil, seasonal-inundation, open-water and
dry-soil fractions before aggregation. WET-3 makes those states consume the
same mass-conserving groundwater/surface-water exchange as PLHY-4 rather than
LPJ-GUESS's fixed wetland runon or saturation source.

**The Earth result is regime-dependent, not a blanket failure.** Australia was
the first test because it resembles the proposed arid, endorheic Vesper state.
There the current solver reaches Pearson +0.0703 on cell means and produces only
1.91 m of spread against 18.36 m observed: the evapotranspiration sink takes the
recharge locally and the flow solve contributes almost no spatial variance.

GW-22 then tested the same solver and 15.19 km mesh against **73,451 USGS bores
explicitly labelled unconfined**. It reaches Pearson +0.2599, with model spread
20.72 m against 28.94 m observed. An affine correction fitted on half the cells
and scored on the other half gives held-out R2 +0.0646 with standard deviation
0.0216 over twenty splits, clearing the declared 0.07 bar on nine. The signal is
real and reproducible; its uncertainty straddles the declared bar rather than
supporting a universal pass/fail label.

**What changed between continents is measurable.** Where `sink_fraction` is
near one, depth is chiefly the local recharge--ET balance. Where it is low and
`at_surface` is set, the cell is pinned and seeping. Low `sink_fraction` with
`at_surface` clear is the regime where lateral flow set the head, and it is the
regime that showed skill against the US unconfined wells. `water_table.nc`
ships both fields so every depth consumer can distinguish these cases.

Vesper's final regime is not yet known. On the current bootstrap climatology
and pre-carve terrain the sink takes a median 97.0% of recharge, but both
recharge and basin structure change when the climate/terrain loop closes.
Re-read `sink_fraction` after each accepted baseline rather than transporting
either the Australian failure or the US result wholesale. Uses that need a
per-cell depth should preserve this regime flag and the permeability bracket;
wetland area additionally needs the native-mesh seasonal treatment above.

**The depth field is a bracket, not a value.** Gleeson's within-class spread is
1.5 to 2.5 orders of magnitude. `--sigma -1` and `--sigma +1` shift every class
by its own standard deviation and are the arms; the central run alone
overstates what is known.

### Crossing to the climate grid: `build_groundwater_access.py`

`water_table.nc` is per REGION and the biosphere is driven per climate-grid
CELL, so something has to reduce ten million regions to a few thousand cells.
Averaging the depth is the reduction that suggests itself and it destroys the
thing groundwater is in the model for: root access and discharge vegetation
depend on the FRACTION of a cell with a shallow table, and a mean erases a
narrow riparian strip in a dry cell while manufacturing a broad one from a
single wet corner. PLHY-5.

So the access criterion is evaluated on the native regions, where the depth
lives, and **what crosses to the grid is an AREA**. No hypsometry is used and
none is needed -- the within-cell distribution is in hand, and reconstructing it
would be an inference where there is a measurement. GRID-2's hypsometry is for
criteria that are terrain-conditioned and have no native evaluation.

`sink_fraction`, `at_surface`, the lithologies with no assigned permeability and
the cells where the sink fraction is undefined each cross as their own share of
the cell's land, so the regime flags this README says every depth consumer must
preserve are preserved by construction rather than by discipline. The
permeability bracket crosses too: `--water-table` once per member writes every
fraction along a member axis labelled by each file's own `sigma`.

Access depths are SWEPT and never chosen. 1.5 m is the vendored LPJ-GUESS soil
column, `SOILDEPTH_UPPER + SOILDEPTH_LOWER`, and is therefore the depth below
which the biosphere as it stands has no roots at all; 5 and 20 m are declared
model-form brackets for a phreatophytic access that model cannot represent. The
artifact carries all three and the consumer argues for one.

`--selftest` runs the reduction against six invariants on synthetic input, with
no build and no solve. The one that matters is the negative control: a cell that
is deep everywhere except a narrow shallow strip must return the strip's area
share, which a depth-averaging reduction would fail while every conservation
check still passed.

**Evaporite has no permeability and is not given one.** Gleeson puts it in the
"not assigned" row with water and ice. Those regions are removed from the
conductive network and their recharge leaves as local seepage, which is what the
surface-only balance already does with it. The alternative policies are in the
config with what each would be claiming.

### Transmissivity is `K D`, and the depth decay it replaced

`T = K D`: Gleeson conductivity over a constant aquifer thickness. `K` and `D`
come from one source at one scale -- Gleeson et al. (2011) put the
permeabilities at 5-100 km and the lithology maps carrying them at "the shallow
subsurface (on the order of 100 m)", against a region about 15 km across.

Transmissivity therefore does not depend on the head, the matrix is fixed, and
the problem is a box-constrained linear complementarity problem that converges
in tens of passes with no relaxation and no damping.

**What this costs, and it is real.** Aquifer thickness no longer depends on
terrain slope, so regolith does not thin on steep ground. Transmissivity no
longer falls as the water table drops, which is the confined approximation and
overstates flow in dry ground. This is the coarser model, taken because it is
the one the sources support at this resolution.

**What it replaced.** Fan et al. (2007)'s exponential decay at a slope-dependent
e-folding length has no cell-mean value here: `exp(h/f)` is convex, so the mean
over a cell whose water table varies by `sigma` carries `exp(sigma^2/2f^2)`, and
100 m of sub-grid relief against the 0.95 m `f` Fan's curve reaches on steep
bedrock makes that `exp(5000)`. Fan's `f` is a hillslope closure fitted over
1.25 km cells; mixing it into a 15 km regional permeability was the original
error. Two solvers were built on it and neither converged.
`notes/water-table-convergence.md` carries both, and they are kept because they
are what a future attempt would otherwise repeat.

### Two options on the same field, and only one changes the problem

`config/groundwater.yaml` carries an aquifer thickness that can be a FIELD and a
transmissivity that can be UNCONFINED. They read the same number and they are not
the same change. Both are off by default and
`notes/subgrid-water-table.md` carries the arithmetic and the verdicts.

**`aquifer.thickness_source: cover_thickness`** takes the saturated thickness
from the export's own surviving cover over basement, floored at Gleeson's 100 m,
so the range is sourced rather than assumed. It is a THICKNESS and not a
history, which is why `docs/src/reference/no-time-axis.md` does not refuse it.
`T` stays independent of the head, the matrix stays fixed, and every identity
this component is certified by is untouched. What it does NOT do is source
GW-17's uniform 2 km: only a small fraction of this world's land carries cover
that deep, so a sourced thickness contradicts that assumption nearly everywhere
and should shrink the depth range rather than explain it. The prediction is in
the config, written before the run.

**`aquifer.transmissivity: unconfined`** is `T = K (h - z_bottom)` with the base
at surface minus that same thickness. It is real physics -- the confined form
overstates flow in dry ground and always said so -- and it is LINEAR in the head,
so unlike Fan's exponential it has an exact cell mean at any resolution. But the
correction it makes is `depth / 2D`, the fraction of the saturated column the
table has drained, and `build_groundwater.py` reports that fraction from the
CONFINED run so the decision rests on a measured number. And it costs: `T`
depends on the head, so the problem stops being a LINEAR complementarity
problem, `--uniqueness-check` stops being an identity and becomes a measurement
against a bar declared in `groundwater.py`, and the outer iteration's
convergence rate is governed by the same number as the size of the effect.
`groundwater.py --dupuit-test` checks the scheme against an analytic
one-dimensional aquifer and prices that trade; nothing may enable this without
the convergence run the note names.

### What is checked, and what missed

`groundwater.py` run directly applies the discrete operator to Legendre
polynomials, which are eigenfunctions of the Laplace-Beltrami operator on the
sphere with an eigenvalue fixed by geometry. It is the only check here that a
wrong face width cannot survive: closure holds for any symmetric weights at all,
and the first face-width estimate was 178x off while passing everything else.
The criterion was declared at 10% relative RMS; `l = 1` passes and `l = 2` to
`4` MISS at about 12%. The error is distributed truncation rather than a few bad
faces, so it is the first-order scheme's own error on an irregular mesh and it
is the water table's mesh-scale noise floor. It is reported rather than tuned
away, and the numbers are in the report.

`--reduction-test` drives permeability to zero, which must reproduce the
surface-only balance exactly rather than closely: every cell returns its own
recharge as seepage and every basin's `Qg` is zero. It passes bitwise.

`--divide-test` does two things. It reports how far the groundwater catchments
agree with the surface ones under a table following the FILLED surface, which is
a measurement at **79.3%** and not a bar: `terminal` comes from a priority
flood's discovery pointer and a flux trace is a steepest-descent rule, and
hydrography's own README says why the flood cannot use descent. The two
disagree on identical terrain by construction. A 95% bar was declared for it
before that was understood and it MISSED; it is left standing as a miss.

It then runs the check that CAN be exact. Handed a flux field whose only
outgoing flux at each cell is the face to that cell's own surface `receiver`,
the trace machinery must reproduce `terminal` on every land region, with no
groundwater involved at all. **It does, at 100.0000% of land area.** That check
earned its place immediately: it caught `groundwater_receiver` reading the
module's sign convention backwards and tracing every cell to the neighbour it
receives most water from, which is to say following the water uphill.

`--uniqueness-check` re-solves from the opposite initial active set. The matrix
is symmetric positive definite, so the complementarity problem has exactly one
solution and any two trajectories must reach it: an identity, not a comparison.
**It passes at 0.000e+00** -- bit-identical head fields from either direction.
Getting there required the dry set to become a static property of the graph
rather than something discovered mid-iteration; before that it missed by 12.76 m
because the two trajectories were marking different cells dry and so solving
slightly different problems.

`--operator-noise` perturbs every face coefficient to ask what GW-8's own
truncation error does to the answer. It is a sensitivity harness, not a model
parameter.

### Two cell areas, and they are not the same

The mesh is a spherical Voronoi tessellation, so its dual is the convex hull of
the region centroids and every face width follows from the circumcentres. The
exact spherical areas that construction implies sum to `4 pi R^2` to eight
figures. The export's own `cell_area` is a different quantity: it sums to 0.068%
more than the sphere, and per region the two disagree by more than 1% over most
of the mesh.

So the solver uses both, deliberately. Fluxes divide by the Voronoi area their
own faces bound, because a finite-volume divergence taken over any other area is
not consistent and the operator check fails by a factor. Water VOLUMES use the
export's `cell_area`, because every other component computes them that way and
the reduction identity has to be exact. Whether `cell_area` should be what it is
belongs upstream in the exporter and is not decided here.

## The topographic index, and the saturated fraction

`build_topographic_index.py` computes `ln(a / tan beta)` per mesh region from
the export's own upslope contributing area and a local slope, and per GRID CELL
the share of that cell's regions whose index exceeds the cell mean. It is
terrain only: no climate, no water table, no solve. GW-26.

**It answers a different question from the water table, deliberately.**
`notes/subgrid-water-table.md` states the constraint the depth field runs into
and the rule it licenses: the valley-to-ridge gradient is below this mesh and
stays there, and sub-grid information may reach a cell-scale parameter only as a
statistic of the distribution the cell contains. TOPMODEL carries the terrain
control statistically and returns a saturated FRACTION rather than a depth,
which is what a wetness classification, a discharge mask and a hydrologic mosaic
consume anyway. The form is
`f_sat = min(f_sat_max * exp(-f_grad * z_wt), 1)`, SIMTOP after Niu et al.
(2005), implemented independently in PALADYN, ClimaLand and CLIMBER-X.

**A fraction is a fraction of a population, and there is only one.** The
population is the mesh regions inside a climate-grid cell, so `f_sat_max` lives
at the grid and there is NO per-region saturated fraction here. A region has no
sub-population and inventing one is exactly what the note refuses. The artifact
is therefore one per (build, grid), on the coupling matrices' precedent.

**`f_grad` is declared and bracketed, never fitted, and the bracket is a
convention.** CLIMBER-X writes `exp(-f_wtab * w_table)`; ClimaLand writes
`exp(-f_over/2 * z_wt)`, where the same numeral means half as much per metre.
Mixing them is a silent factor of two, and the bracket is exactly that factor.

**Absolute thresholds are refused and no product keys on one.** The published
ones -- a CDF on integer bins 1 to 15, a critical cell mean of 5.5, a cut at 14 --
are calibrated against an index computed on Earth at about a kilometre. This
index carries a length in `a`, so the whole distribution shifts with the mesh,
measurably so between this project's own two builds, and it sits far above that
range. `f_sat_max` is a share of a cell's population above that cell's own mean,
so it survives the offset; a threshold does not.

**Nothing may consume `f_sat` until it is scored.** The bar is in
`config/topographic_index.yaml`, declared before any fraction was computed, and
the report's `consumers_licensed` says whether it has been met. The closure is
strictly monotone in the depth for a fixed `f_sat_max`, so any discrimination it
gains over the depth alone is attributable to the terrain half and to nothing
else -- which is what makes the score a test rather than a comparison.

## The Earth comparator

Earth's endorheic fraction is **about one-fifth of its land surface**: Wang et
al. (2018), "Recent global decline in endorheic basin water storages", Nature
Geoscience 11, 926-932, measures 31.8 million km2 across 48,813 landlocked
watersheds from 15-arcsecond HydroSHEDS. The units match ours: Wang's figure
is catchment area that drains internally, the same measure as this component's
endorheic share of land -- not the basin-floor area, a quantity several times
smaller and easy to confuse. This world is more endorheic by roughly a factor
of two after carving, and closer to four before it; that follows from the fork
preserving closed basins, and it is the fact this component exists to handle.
The share falls with each carve, and it is the number most likely to be quoted
from a stale per-build product: these files are per-build, and a figure
computed against one terrain says nothing about another.


**That figure is the hyper-arid limit, not a property of the world.** It assumes
no basin ever overflows. A basin that overflows year on year incises its outlet,
and over the 1e4 to 1e6 years a landscape needs to relax, that drains the lake
and the depression stops existing. So a basin pinned at its spill is a transient,
not a landscape state, and the endorheic share depends on how many basins the
climate keeps overflowing:

The dependence is steep and monotone: as the aridity threshold rises from humid
to hyper-arid, the number of basins that carve falls by more than two orders of
magnitude while the endorheic share of land rises from under a tenth to roughly
three quarters. The counts are not reproduced here, because they belong to
whichever terrain produced them and the SHAPE is the point being made; no
registered step emits the sweep today, so re-derive it from `carve_verdict.py`
at varied aridity when a build needs it.

Regenerating it means running
`carve_verdict.py` against the carved terrain, and that is an iteration-2 carve
decision rather than a documentation chore, so it is left for whoever takes that
decision.

The decision variable is climate-free per basin and is stored as
`critical_aridity_index` in `basins.nc`. A basin overflows exactly when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

so the index is pure geometry and the climate supplies one number per basin.
Median here is 3.03 on the carved terrain, meaning the typical basin overflows
unless evaporative demand over open water exceeds about three times the runoff
depth. Earth's
surviving endorheic basins sit well above that, which is why Earth keeps only
about one-fifth. `carve_verdict()` in `lake_balance.py` applies it.

## What this component cannot fix

**We cannot carve the terrain here.** Fixing a rim means changing elevation,
which belongs upstream in World Orogen: incision depends on the `erodibility`
field, and a second copy of the terrain in this component would desynchronise
ExoPlaSim's orography from the maps and from everything else reading `source/`.

What this component can do is decide *which* basins should carve, which is what
`carve_verdict.py` and `export_carve_list.py` are for. The generator takes it:
`--preserve-basins FILE` accepts a retain fraction per basin, where 1 keeps the
rim, 0 carves it and anything between cuts a notch and tapers it over the divide
band. `--carve-basins FILE` is subtractive. Both are documented in the fork's
`tools/README.md`.

**`carve_list.json` records the criterion's inputs, not only its verdict.**
Alongside each basin's retain fraction it writes the three catchment means the
overflow test is built from -- `precipitation_km_per_year`,
`lake_evaporation_km_per_year` and `land_evaporation_km_per_year` -- in the same
km/year depth as `runoff_km_per_year`. They are there so a consumer can
re-evaluate the test under a perturbation without restating the formula.
`scripts/error_budget.py` is that consumer, and it used to recover P and E by
inverting the recorded aridity index and evaporation margin instead: exact
algebra, a second copy of the criterion, and wrong the first time the criterion
moved. Under `runoff_source p_minus_e` the recorded runoff is
`max(P - E_land, 0)` by construction, which is the identity to check before
trusting them.

The count of retain-0 entries in `carve_list.json` matches
`carvedByRetainZero` in the resulting build's manifest, which is the check
that the list was consumed as written. Iteration 2 is a run, not a generator
change.

The terrain still is not in equilibrium: many basins fill to their spill under
this climate, which is the water balance saying their outlets should have been
cut. They hold a large share of both the land and the lake area on the map, so
what is drawn is largely a landscape that has not relaxed yet. That is the case
for another carve iteration, and the counts are in `world_state.json` and
`surface_water.nc` rather than here, because they move with every pass.

## Known approximations

- **The surface balances take net groundwater flow as zero.** Runoff is P-E
  over land, accumulated down the drainage network, and a basin's catchment is
  its surface catchment. Fan (2019) writes the budget as
  `dS/dt = P - ET - Qr - Qg` and notes that dropping the last term is the
  common assumption; it is the one `surface_water.py` and `carve_verdict.py`
  make, and it is weakest in the dry closed-basin regime this world is mostly
  made of. `build_groundwater.py` measures the term rather than removing it,
  and nothing downstream consumes it yet.
- **A few overflow paths disagree with the cascade.** Most overflowing basins
  have a saddle opening into exactly the basin the solver routes them to. A
  handful are cycle-collapsed: their target was reassigned to the cycle's
  primary while their saddle still leads where it physically leads, so the
  routed river and the solver's bookkeeping name different destinations.
- **The two grids are offset half a cell in longitude.** Orogen's exported grid
  centres its columns at -180 + (k + 0.5) dlon; an ExoPlaSim climatology centres
  its own at k dlon. Remapping by nearest column leaves up to 1.4 degrees of
  slip, about 185 km at the equator. Whether the boundary conditions handed to
  ExoPlaSim carry the same offset is an `exoplasim/` question and has not been
  checked here.
- **Catchment forcing is binned, not integrated.** A basin's runoff is a
  catchment-area-weighted mean over whole grid cells, so a catchment that
  straddles a cell gets that cell's value entire. Comparing a basin's solved
  overflow against the discharge the mesh delivers to its sink, the ratio has a
  median of 0.99, exactly as it should for a lake passing on everything but its
  own evaporation, and a ninetieth percentile of 1.13. That spread is this
  binning.
- **Merged basins.** 205 basins were collapsed into a neighbour to break spill
  cycles, which arise where basins share a saddle and each names the other as
  its outlet. Physically they merge into one lake at that level. The survivor's
  hypsometry still describes only its own depression, so a merged group filled
  past the shared saddle holds more water than the curve says.
- Catchments are assigned to grid cells by region centre, matching the
  exporter's categorical resampling rule, so a cell straddling a divide is
  attributed whole.
- The solver finds equilibrium, not a seasonal cycle. Basins with large storage
  relative to annual throughput will lag; nothing here models that.

## One-off tools

- `scripts/build_earth_wtd_sites.py` -- one-off: assembles Australian bore water table depths into one site table for GW-3, the first leg of the solver's Earth comparator. Reads the Australian Groundwater Explorer's per-state download, whose `level_<state>.csv` sits BESIDE the geodatabase rather than inside it. Reconciles three datums against an identity that can fail, and reports the quality flag rather than filtering on it, because the codes are per-agency and filtering would select which state survives. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/earth_calibration.py` -- one-off: the GW-3 Earth comparator end to end, mesh, ETOPO, GLHYMPS permeability, recharge, drainage, solve and score, with `--edge-km` or `--regions` setting the mesh so the same case can be run at more than one cell size. The original was written inline and lost, which is GW-20 and is why this exists as a script: the one number that judges this component was recorded and not reproducible. Stages cache under `data/earth_validation_cache/edge<E>km/` and are skipped when present; the cache is gitignored and regenerable, the result is `analysis/earth_calibration.json`. Runs with the sink and the river baselevels ON, as the pipeline gate requires. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/build_us_wtd_sites.py` -- one-off: the United States half of the Earth comparator's observations, GW-22, companion to `build_earth_wtd_sites.py`. `external-data.md` says take Australia FIRST and the US second; this is the second leg, and it exists because GW-21 concluded the limit on the Australian score is the input fields rather than the formulation or the mesh, which is a claim resting on one region. Fetches USGS `monitoring-locations` and `field-measurements` at parameter code 72019, depth to water in feet below land surface, paging blind on the `next` link because `numberMatched` is absent from every response. The US set can be BETTER than the Australian one rather than merely bigger: `aquifer_type_code` separates confined from unconfined directly, where Australia can only infer it from bore depth, and a `Static` qualifier marks readings not taken while pumping, which `external-data.md` notes Fan's four columns cannot support. Writes `data/earth_validation/us_wtd_sites.csv`. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/fetch_recharge.py` -- one-off: pulls a window of the global recharge grid by BYTE RANGE rather than downloading it. `RechargeTotal.nc` is a classic netCDF3 whose one variable is a contiguous big-endian float64 array, so whole rows can be ranged out of a multi-gigabyte file; the script carries the data offset and the exact windows, which is the part `docs/src/reference/external-data.md` cannot hold compactly. Windows are deliberately WIDER than the region they serve wherever a bbox edge would cut land, because a cell with no coverage is not land, so it is ocean, so it is a fixed head at sea level. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/validate_lake_solver.py` -- one-off: checks the lake solver against real endorheic basins on Copernicus DEM tiles, which is where HYD-7's mesh-scale storage deficit was measured. Registered under `one_offs` in `config/pipeline.yaml`; it generates nothing the pipeline reads.
