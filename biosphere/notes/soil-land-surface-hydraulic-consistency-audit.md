# Soil and land-surface hydraulic consistency audit

**Recorded:** 2026-08-21
**Scope:** the water and heat path from ExoPlaSim's land surface through
pedology, hydrography and groundwater into the vendored LPJ-GUESS-CNP tree.
This was a source, artifact and pipeline audit only. No ExoPlaSim or
LPJ-GUESS simulation was run.

## Result

The current system does not have one terrestrial water column. It has two
independent fast water balances and a third subsurface balance:

- ExoPlaSim carries one liquid land bucket, a separate snowpack and five soil
  temperature layers. Its bucket controls latent heat and produces routed
  runoff when it overflows.
- LPJ-GUESS reconstructs rain, snow, melt, interception, soil liquid and ice,
  evapotranspiration and runoff from the climate driver. Its 15 water layers
  are parameterized independently and do not debit ExoPlaSim's bucket.
- the groundwater solver treats positive annual `P-E` as recharge, removes
  groundwater ET and returns seepage/baseflow on the hydrography mesh. Its
  state and fluxes do not presently bound either fast water balance.

This is more than differing model detail. The same atmospheric surplus is
currently allowed to be ExoPlaSim surface runoff, catchment runoff and
groundwater recharge without an infiltration/drainage partition. At the same
time, the two vegetation-capable models can evaporate separate stores. A
closed annual `P-E` budget therefore does not prove that the coupled system
conserves water or represents its timing.

The preferred central architecture is one restart-exact, climate-side fast
land water and energy column. Pedology supplies a single vertically explicit
hydraulic and thermal property artifact. ExoPlaSim owns precipitation
interception, snow and soil phase, infiltration, storage, latent heat, surface
runoff and drainage at atmospheric timesteps. LPJ-GUESS consumes the same
state and flux history and returns vegetation boundary controls such as LAI,
canopy capacity, stomatal conductance and root demand. The groundwater model
receives drainage and returns one signed lower-boundary exchange. Every
withdrawal is debited once.

That architecture is a recommendation, not a commitment to a maximal Richards
solver. The current bucket, a parsimonious multilayer column and a
gradient-driven column are testable model hypotheses. The conservation core,
state ownership and property contract should be fixed before choosing among
them. This follows Clark et al.'s separation of conservation structure from
flux parameterization and spatial representation.

## Sources read

The source trace covered ExoPlaSim `landmod.f90`, `fluxmod.f90` and restart
calls in `plasim.f90`; `build_surface_soil_water.py`; the pedology soil
generator; surface-water and groundwater builders; the ecological driver;
LPJ-GUESS `soilinput.cpp`, `soilwater.cpp`, `soil.cpp` and
`vesperinput.cpp`; and the relevant pipeline and README contracts.

Four additional primary sources were fetched into `references/`, read, and
recorded in `references/INDEX.md`:

- Cosby et al. (1984), the texture regressions used by LPJ-GUESS;
- Lawrence and Slater (2008), the organic-soil mixing used by the fork;
- Best et al. (2011), a coupled land-surface precedent that makes surface
  partitioning, vertical flow, roots, freezing, heat and lower boundaries one
  conservation problem; and
- Clark et al. (2015), the SUMMA framework for testing alternative process and
  scaling hypotheses inside a common conservation core.

The earlier plant-hydraulics audit supplies the direct LPJ-GUESS-RE and
LPJ-GUESS--ParFlow precedents. Earth regressions, profile shapes and fitted
parameters in all of these sources are provenance and model-form evidence,
not measurements of Vesper.

## Findings

### 1. One pedology artifact becomes two incompatible vadose-zone soils

Pedology writes texture, organic carbon, bulk density, regolith depth, an
available-water capacity, weathered-bedrock fraction, andic fraction and
phosphate fixation. ExoPlaSim's project-added surface builder reads only
`awc`, converts millimetres to metres and installs it as `dwmax`.

LPJ-GUESS does not consume that capacity. `soilinput.cpp` independently derives
saturation, matric potential, wilting point and field capacity from the Cosby
texture regressions, then derives an empirical percolation coefficient.
`iforganicsoilproperties 1` makes it redistribute total soil carbon down a
fixed Earth profile and mix mineral and one ideal organic soil following
Lawrence and Slater. `vesperinput.cpp` subsequently scales those independently
derived capacities by regolith depth and weathered-bedrock fraction.

The paths disagree by construction:

- pedology uses declared sand/silt/clay capacity endmembers plus additive
  organic and allophane terms; LPJ uses Cosby retention parameters and ignores
  pedology's AWC, bulk density and andic hydraulic effect;
- pedology and LPJ calculate organic volume and its vertical distribution by
  different rules; and
- ExoPlaSim sees regolith depth through total AWC but no weathered-bedrock
  storage, while LPJ adds capped sub-regolith storage without changing its
  physical layer or root geometry.

A no-simulation calculation on the checked-in 4,105-cell pre-carve soil map
quantifies the seam. With its zero initial soil carbon, applying LPJ's active
Cosby equations and the current Vesper regolith/bedrock scaling gives a median
capacity of 145.3 mm against pedology's 84.6 mm. The median LPJ/pedology ratio
is 1.61, its 10th--90th percentile range is 0.90--2.55, and the median absolute
difference is 48.6 mm. The fields correlate (Pearson 0.793) because both start
from related texture and depth, but they are not the same capacity. These are
dated artifact diagnostics, not current-world results.

Cosby et al. support texture as a useful predictor of mean hydraulic behavior,
but emphasize that large residual variance remains within texture classes.
Lawrence and Slater show that organic matter changes hydraulic and thermal
properties together, and their results depend on vertical organic state. A
canonical artifact therefore needs retention and conductivity parameters,
uncertainty and layer geometry; one clipped AWC number is not sufficient.
This does not mean copying a vadose-zone conductivity into the groundwater
solver. Aquifer permeability/transmissivity lives at a different material and
spatial scale and remains a separate groundwater property, joined to the land
column through an explicit contact and lower-boundary contract.

### 2. AWC is being asked to stand in for an entire land column

ExoPlaSim's active hydrology is one prognostic liquid store `dwatc`. Each
timestep adds the net water flux, calls any excess over `dwmax` runoff, and
clips the store between zero and capacity. Evaporation is throttled by the
store divided by 40% of capacity. The model has no vertical liquid layers,
canopy interception, root extraction, drainage flux or lower hydraulic
boundary in this path.

Its five soil-temperature layers use fixed global heat capacity and thermal
conductivity. They do not respond to texture, organic/andic material, liquid
water, ice or regolith/bedrock geometry. Snow is an energy-bearing state, but
soil-water phase does not constrain infiltration or runoff.

Changing `dwmax` was a useful project integration step: it proved that
pedology can alter climate storage and runoff. It cannot make a scalar AWC
encode infiltration capacity, unsaturated conductivity, matric potential,
freeze/thaw, capillary rise, vertical redistribution or drainage timing. Nor
does a deeper bucket represent groundwater access; it changes storage and
memory without imposing a water-table boundary.

Best et al. provide a relevant precedent rather than a prescription. JULES
puts interception and throughfall, surface runoff/infiltration, root extraction,
Richards flow, phase-aware soil heat and drainage inside one coupled surface
balance. It also shows that apparently technical choices about routing excess
water upward or downward change snowmelt runoff, subsequent soil moisture and
atmospheric temperature. Those choices must be exposed as hypotheses at
Vesper's grid and gravity, not hidden in a capacity map.

### 3. Positive `P-E` currently has mutually inconsistent owners

The hydrography surface-water path rejects ExoPlaSim `mrro` as an incomplete
routed diagnostic and computes annual runoff as `max(P-E, 0)`. Pedology uses
the same quantity for weathering and regolith formation. The groundwater
builder then uses that same positive `P-E` as recharge. Meanwhile ExoPlaSim
has already routed overflow from its land bucket as runoff.

At equilibrium, `P-E` is a sound total terrestrial export constraint. It does
not say how export divides between infiltration-excess or saturation-excess
surface runoff, soil drainage, lateral interflow, groundwater recharge,
baseflow and changing storage. Treating the total as each component at once is
double ownership, even if the groundwater solver later conserves its own
recharge and the basin code carefully adds only signed groundwater
redistribution to catchment supply.

The pipeline needs one flux ledger with one owner for every store and crossing:
precipitation, canopy and snow storage, liquid and ice by soil layer,
evaporation, transpiration, surface runoff, drainage/recharge, signed
groundwater exchange, baseflow and open-water storage. Pedology may consume an
accepted runoff or drainage climatology for weathering, but it should not
invent the partition from the atmospheric residual.

### 4. ExoPlaSim and LPJ-GUESS evaporate different water histories

LPJ-GUESS reconstructs rainfall/snowfall, snowmelt, interception, infiltration,
layer water, plant uptake, evapotranspiration and runoff from the ecological
weather input. None of those withdrawals debits `dwatc`, so the latent heat in
the climate and the water stress in vegetation can be based on different
stores. Conversely, vegetation currently changes climate albedo and will
change roughness, but LAI, roots, interception and stomatal state do not control
ExoPlaSim evaporation.

This is the open hydrological biosphere--climate loop. Passing LPJ soil water
to ExoPlaSim after the fact would not close it: atmospheric latent heat and
precipitation partitioning are fast processes and must use the same state when
the flux is evaluated. The central owner should therefore be the climate-side
land column, or a shared kernel called by it. LPJ returns vegetation boundary
conditions and demand; the owner computes the accepted flux and returns the
fulfilled uptake by source. This prevents both models from withdrawing the
same millimetre.

The current pipeline is sequential rather than an online coupler. It therefore
cannot use vegetation calculated from the future climate segment to control
that same segment. LSHY-4 must choose either a shared/online operator or an
explicit fixed-point outer iteration that exchanges chronology-bound controls
and demonstrates convergence of both vegetation and hydrologic fluxes. Replaying
LPJ controls into a changed climate without convergence would only move the
inconsistency by one iteration.

Jia et al.'s LPJ-GUESS--ParFlow coupling is a useful ownership precedent: the
hydrologic model supplies state and LPJ returns a boundary flux. It is not a
reason to import ParFlow or its Earth configuration. Here, latent heat makes
the atmospheric land scheme the natural fast owner, with the independent
groundwater solver supplying its lower boundary.

### 5. Snow, frozen soil and thermal state are reconstructed twice

ExoPlaSim produces snow and evolves soil temperature, while LPJ-GUESS
independently partitions precipitation around an Earth-derived temperature
threshold, melts its own snow, evolves soil temperature and ice, and reduces
available liquid water under freezing. The present forcing does not carry the
producer's snow water equivalent, melt, soil temperature, liquid/ice state or
runoff/drainage chronology.

The two columns can consequently disagree on when water reaches soil, whether
it is liquid, whether frozen pores impede movement and whether melt becomes
runoff. Because freeze/thaw also changes heat capacity and latent heat, this is
an energy-consistency problem as well as a biological input problem.

The ecological forcing contract should carry authoritative interval-bounded
land state and fluxes from the accepted climate column. LPJ may retain an
independent snow/soil reconstruction only as a named structural control. The
canonical soil artifact must also provide moisture-, phase- and
composition-dependent thermal properties or a declared reduced approximation.

### 6. Both active lower-boundary treatments are limited

The standard 15-layer LPJ path moves deep water only on days for which rain
plus snowmelt is at least 0.1 mm. It has no dry-day drainage, upward capillary
flow or physical groundwater boundary. LPJ-GUESS-RE demonstrates that this can
preserve an artificial deep dry-season reservoir and that lower-boundary
choice materially changes vegetation. PLHY-2 already owns that finding.

The consistency work changes how it should be resolved. Porting a more complex
LPJ-only column while leaving ExoPlaSim's bucket active would improve one of
two balances but make their structural mismatch larger. The corrected LPJ
bucket and LPJ-GUESS-RE remain useful labelled controls. The central coupled
case should consume the shared land-column state and the single groundwater
exchange defined under PLHY-4.

The groundwater side is also an equilibrium model. It supplies a steady-state
water table and annual-mean balance, not seasonal head, aquifer storage or a
transient capillary/baseflow response. A first coupled bracket may use a fixed
head or orbit-mean signed exchange and iterate the mean recharge and ET to
closure. It must not manufacture a seasonal groundwater response from that
artifact. PLHY-4 and LSHY-3 need to decide whether that equilibrium boundary is
adequate or a transient saturated-storage state is required by the responses
being claimed.

### 7. Cell-mean capacity erases the hydrologic mosaic

`build_surface_soil_water.py` area-averages land and lake capacities into one
T42 `dwmax`. That makes the lake fraction a wetter local bucket, but routed
river water, seepage and baseflow cannot enter an open-water store and return
to the atmosphere. A cell containing a dry rootable upland and a small
groundwater-fed playa therefore behaves as one medium-wet surface.

This is not fixed by choosing a better mean. Runoff generation, groundwater
access, lake evaporation and root access are nonlinear and occur on different
fractions. The hydrologic representation needs conservative response tiles or
an equivalent subgrid distribution built from BIO-11/BIO-17 surface fractions,
PLHY-5 native-mesh groundwater access and GRID-2 hypsometry. Mean water-table
depth or mean AWC must not substitute for an accessible-area fraction.

### 8. The current interfaces are project scaffolding, not immutable constraints

Several awkward features are ours: `soilmap.txt`'s extra depth/AWC columns,
the scalar `dwmax` builder, the regolith capacity rescaling, the duplicated
`P-E` use and the 12-bin ecological driver. They established end-to-end wiring
but can be replaced. The vendored ExoPlaSim bucket is also editable source,
not an external service contract. The audit should not optimize a new design
around preserving those choices.

The checked-in pedology report also contains a stale note saying LPJ-GUESS
does not consume regolith depth. The generator has been corrected to describe
capacity scaling without physical-layer or root-geometry changes. The derived
report was deliberately not hand-edited or regenerated during this audit.

The constraints worth preserving are scientific and operational: exact mass
and energy ownership, restart continuity, explicit units and time bounds,
Vesper gravity in potential/head conversions, conservative scale crossings,
and the ability to reduce a richer implementation to named simpler cases.

## Recommended implementation order

1. Define the canonical soil hydraulic/thermal property artifact and compare
   both current consumers against it without changing a model.
2. Define the state/flux ownership ledger and its surface--vadose--groundwater
   boundaries.
3. Implement the restart-exact climate-side land column behind a selectable
   bucket reduction, then carry its state and flux chronology through EFOR.
4. Close vegetation demand/fulfilled uptake and groundwater lower-boundary
   exchange without a second withdrawal.
5. Add phase/thermal consistency and subgrid hydrologic response tiles.
6. Only then run the registered structural cases and propagate their spread.

The first two stages and deterministic fixtures need no climate or
LPJ-GUESS execution. Timing measurements, climate response cases and
ecological comparisons still require explicit permission.

## Acceptance boundary

Before a coupled result is interpreted, the implementation must demonstrate:

- identical declared hydraulic/thermal properties at each consumer boundary,
  including units, gravity convention, layer geometry and uncertainty case;
- exact closure of precipitation to storage change, ET components, runoff and
  drainage, and of drainage plus signed groundwater exchange to groundwater
  storage/ET/baseflow;
- one and only one debit for canopy evaporation, soil evaporation,
  transpiration and groundwater uptake;
- rain/snow and liquid/ice partition closure, with water and energy conserved
  across freeze/thaw;
- restart equivalence through active precipitation, melt, drainage and
  vegetation-demand intervals;
- conservative native-mesh to the accepted climate support area and flux
  aggregation under SPAT-1 through SPAT-4; and
- analytical reductions in which the multilayer column reproduces the
  registered single-bucket case, impermeable and free-drainage limits behave
  correctly, and zero groundwater exchange gives the no-groundwater case.

Annual `P-E` closure, a plausible LAI map or agreement of two cell-mean water
capacities is not sufficient: each can coexist with compensating flux and
state errors exposed above.
