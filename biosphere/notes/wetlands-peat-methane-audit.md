# Wetlands, peat and methane audit

Worldbuilding. Vesper is invented; every wetland area, peat stock and methane
quantity below is a model input or prediction, not an observation. This audit
was written 2026-08-21 from the vendored LPJ-GUESS-CNP source, the groundwater,
hydrography and atmospheric interfaces, and the primary papers listed at the
end. No LPJ-GUESS build or simulation and no climate run was performed.

## Executive finding

LPJ-GUESS contains a useful northern-peatland process model, but it is not a
planetary wetland or methane model ready to enable. Four distinct problems are
hidden behind `run_peatland` and `ifmethane`:

1. where peatlands, seasonally inundated soils and open waters exist;
2. how groundwater, routing, precipitation and evapotranspiration conserve
   water while producing a daily water table;
3. how peat carbon, vegetation and redox produce, oxidize and transport CH4;
4. how all land and aquatic sources plus soil and atmospheric sinks determine
   an atmospheric abundance and radiative feedback.

The current Vesper configuration enables none of them. That is the correct
fail-closed state. The dormant path requires an externally prescribed
`PEATLAND` land-cover fraction; classifies every such stand north of 40 N as a
detailed peatland and every stand south of it as a simplified inundated mineral
wetland; fixes a 1.5 m Earth peat column and Earth gas boundary conditions; and
has no connection to Vesper groundwater or routed surface water.

The direct conservation blocker this audit found has since been repaired, under
WORLD-C4J8, and finding 3 below records both what it was and what stands in its
place. The low-latitude wetland infiltration path filled every soil layer to
capacity regardless of available rain; `ifsaturatewetlands` controlled only
whether the created water was recorded for a later attempt to subtract it from
runoff, not whether the water was added; and even with the switch on an unmet
remainder could remain as an external source rather than routed runon. The path
is now bounded by the rain that arrived, and the one daily exchange that would
deliver routed and groundwater water to a wetland stand is PLHY-4's and does
not exist, so the wetland has a NAMED absence where it had a source term.

What has to be declared and repaired before any of those switches move is in
`biosphere/notes/wetland-activation-contract.md`, and
`biosphere/scripts/wetland_gate.py` holds the refusal.

The practical route is to derive a seasonally resolved sub-grid wetness
classification from Vesper's hydrography, native-mesh groundwater state and
topographic distribution; make one mass-conserving groundwater/surface-water
exchange own the daily water table; retain the existing detailed methane model
as a northern-Earth structural benchmark rather than a geographic branch; add
the omitted dry-soil sink and inland-water source brackets; and close surface
fluxes through an offline oxidant/lifetime calculation before replacing the
currently prescribed 1.6 ppmv CH4. Dynamic peat history is a separate,
millennial-timescale problem and cannot be invented from the repository's
no-time-axis equilibrium premise.

## What exists now

### Dormant activation

- `vendor/lpj-guess/data/ins/landcover.ins` sets `run_peatland 0`.
- `data/ins/global.ins` sets `ifmethane 0`, `ifsaturatewetlands 0` and
  `wetland_runon 0`; it does not import `wetlandpfts.ins`, and methane outputs
  are not retained by the Vesper harness.
- `modules/soilmethane.cpp` runs only for a `PEATLAND` stand, with emissions
  beginning only in the final 100 simulation orbits of spin-up.
- `framework/externalinput.cpp` reads a `PEATLAND` fraction from an external
  land-cover file. It does not predict wetland or peatland extent. With fixed
  land cover and multiple active classes it can instead assign equal fractions,
  which is not a scientific fallback for this world.

### Two unrelated models selected by Earth latitude

`Stand::is_highlatitude_peatland_stand()` is exactly `PEATLAND && lat >= 40.0`;
`is_true_wetland_stand()` is `PEATLAND && lat < 40.0`. This is not absolute
latitude: a Vesper cell at 60 S follows the low-latitude inundated-soil path.
The branch controls hydrology, vegetation stress, decomposition and methane,
not merely a default parameter.

North of the line, the code uses the Wania et al. peat hydrology and detailed
CH4 diffusion, plant transport and ebullition model. South of it, a saturated
mineral-wetland path converts a fixed fraction of heterotrophic respiration
directly to atmospheric CH4. The original papers imposed northern geographical
domains because their detailed formulation was developed and evaluated there;
that limitation is evidence for a model-form bound, not a planetary classifier.
This discovery prompted the biosphere-wide inventory in
`implicit-earth-assumptions.md` finding 8 and BIO-29. The same rule applies here
as elsewhere: latitude may locate the cell for geometry and diagnostics, but it
may not determine wetland type or any other ecological process regime.

### Peat and gas state

The detailed path fixes three 0.1 m acrotelm layers over twelve 0.1 m catotelm
layers. The catotelm is always saturated. Fixed porosity, water-table response,
root distribution and a 7.5 kg C m-2 acrotelm threshold stand in for peat
structure; the bulk CENTURY C-N-P pools are not vertically indexed. Slow and
passive pools are simply assigned the anaerobic modifier on the assumption that
they occupy the catotelm. Accumulated peat does not change column depth,
hydraulic properties or the acrotelm/catotelm boundary.

Methane production takes daily bulk heterotrophic respiration, distributes it
over the fixed root profile and multiplies it by anoxia and a fixed CH4:CO2
ratio. The detailed ratio is 0.085. The simplified ratio is 0.027 and its source
comment says it was changed "to match global emissions." Carbon emitted as CH4
is subtracted from the soil CO2 flux, and the detailed routine contains useful
carbon checks; those conservation features should be preserved.

The detailed transport path includes diffusion, live-aerenchymatous-plant
transport and ebullition. It fixes gravity at 9.81 m s-2, pressure at 101325 Pa,
atmospheric CH4 at 1.7 micro-atm, O2 at 209000 micro-atm and 10 m wind at zero.
It uses 24 h to convert gas-transfer velocity and 0.01 Earth day gas substeps.
Ebullition is forced directly to the atmosphere even when the computed water
table would place bubbles below an unsaturated layer. Snow at or above the
threshold stops diffusion, and live leaf biomass/phenology removes plant
transport too early to represent dead tiller venting and winter/thaw pulses.

Restart state is incomplete. `Soil::serialize` retains `wtp`, `awtp` and
several yesterday gas stores, and does not retain `Wtot`, `wtd`, `stand_water`,
`mwtp`, `Frac_ice` or `rootfrac`. `Frac_ice` is the current ice fraction and is
distinct from the `Frac_ice_yesterday` that is retained; `rootfrac` is
initialized on day zero only. An arbitrary-day restart therefore needs an
exact-continuity fixture before either peat hydrology or methane output can be
accepted. The annual water-table update itself tests
`date.day == Date::MAX_YEAR_LENGTH`. `Date::next()` resets `day` to 0 on the
last day of the last month, so `day` never exceeds `MAX_YEAR_LENGTH - 1` and
that block is unreachable on the original and the Vesper calendar alike. `awtp`
therefore holds the 0.0 it is initialized to for the whole run, and across
restarts because it is serialized. The consequence reaches the carbon:
`update_acrotelm_co2` interpolates the acrotelm CO2 concentration linearly in
`awtp` between the atmospheric value at -300 mm and the pore-water value at 0,
so `awtp == 0` pins the simulated acrotelm to the pore-water concentration
everywhere and always.

### Vesper interfaces

The groundwater solver already supplies a much better physical seam than
`wetland_runon`: recharge, local baselevels, gravity-correct hydraulic
conductivity, native roughly 15 km regions and a `sink_fraction` diagnostic.
Its per-cell depth is not itself an adequate wetland map, and the LPJ/ExoPlaSim
cell is coarser still. Wetland extent depends on sub-grid relief, routing,
seasonality and persistence, not mean water-table depth alone.

Hydrography already resolves basin membership, routing, solved lakes, rivers,
playas and surface classes. These can distinguish persistent open water from
land that is saturated or seasonally inundated, provided the native-mesh
fractions are aggregated rather than thresholded after aggregation. PLHY-4 and
PLHY-5 already require one mass-conserving plant/soil/groundwater exchange and
native-mesh preservation. Wetlands must consume that same flux; they cannot
make a second withdrawal or runon term.

`config/planet.yaml` prescribes 1.6 ppmv CH4 from a fixed-modern-Earth-biogenic-
flux photochemical calculation. ExoPlaSim can radiatively consume that
abundance, but neither predicts atmospheric CH4 nor transports and oxidizes a
biosphere flux. BVOC-6 already owns the reduced Vesper oxidant and methane-
lifetime seam.

## Findings

### 1. Activation needs a source, state, output and acceptance contract

The existing off switches correctly prevent an ungrounded result. Enabling
`run_peatland` without an input fraction can assign equal land-cover shares;
enabling `ifmethane` then turns only the `PEATLAND` share into a source. The
standard Vesper PFT file does not define wetland plants, and the run harness
does not retain the pathway, water-table, stock or closure diagnostics needed
to interpret it.

The activation gate must require accepted wetland-class fractions, daily
hydrology, planet gas/pressure/gravity inputs, C-N-P-water closure and exact
restart continuity. Missing inputs must fail rather than fall back to equal
fractions, a fixed runon or the 40 N branch. BIO-28 remains the small planetary-
constant guard; the WET tasks below make the full pathway fail closed.

### 2. Wetland extent is a first-order prediction, not one land-cover column

WETCHIMP found that internally derived Earth wetland areas differed by nearly a
factor four and that area disagreement propagated directly into emission
uncertainty. Kleinen et al. likewise found peat carbon strongly dependent on
minimum, mean or maximum inferred extent. The problem is larger on Vesper,
where there is no observed wetland mask to prescribe.

One `PEATLAND` fraction collapses states with different carbon and gas physics:
persistent peat-forming land, saturated but non-inundated mineral soil,
seasonally flooded land, open lakes/rivers/playas and dry mineral soil. Vesper
needs mutually exclusive native-mesh fractions and transition rules. The
fractions must close to the rootable and open-water surface areas already used
by BIO-11 and hydrography, with explicit prevention of double counting.
Sub-grid topographic index or hypsometry may provide a bounded extent model,
but no latitude threshold or coarse-cell mean depth may select the physics.

### 3. Low-latitude wetland hydrology created water; the path is now bounded

As audited, `initial_infiltration()` had a low-latitude `PEATLAND` stand compute
every layer's saturation deficit and add that full deficit. It did so even when
`rain_melt` was smaller: `rain_melt` was clamped to zero and the whole deficit
added regardless. `ifsaturatewetlands` only recorded the amount, and recorded the
whole `total_potential` rather than the created part, so the record over-stated
even where the rain could have supplied the water. `hydrology()` then tried to
subtract the record from runoff; where runoff was too small the remainder
survived as `awetland_water_added`. With the switch off the same water was added
and not recorded at all. Three further defects sat on the same path: a
`soiltype.runon = wetland_runon` scalar in place of a routed quantity, an annual
water-table average guarded on `date.day == Date::MAX_YEAR_LENGTH` which
`Date::next()` never produces, and a `Soil::serialize` that omitted `Wtot`,
`wtd`, `stand_water`, `mwtp`, `Frac_ice` and `rootfrac`.

All four are repaired under WORLD-C4J8 and
`biosphere/notes/wetland-activation-contract.md` section 2 describes what stands
in their place: infiltration bounded by `min(soil.rain_melt, total_potential)`
and distributed in proportion to each layer's deficit, no runon term at all
until PLHY-4's exchange declares one, the annual average on
`Date::MAX_YEAR_LENGTH - 1`, and the six members serialized with
`biosphere/scripts/verify_lpj_restart_continuity.py` as the behavioural check.

What is NOT repaired is the model this path needs. The detailed peat path is
still not closed: its own source comment notes that its evapotranspiration
treatment makes the hydrologic cycle impossible to close. Wania et al. identify
missing groundwater as a primary reason fen water tables and vegetation can be
wrong and state that bog/fen separation requires it. There must be one daily
water ledger linking precipitation, snow, surface routing, PLHY uptake/ET,
groundwater recharge/discharge, storage and runoff; the wetland receives actual
routed or groundwater water and returns all losses, and saturation is a state
outcome, never a source term. `hydrography/config/land_water_ledger.yaml` is the
graph that ledger closes on and PLHY-4 owns the exchange, which is why
`hydrology.water_ledger` is the one field of the hydrology block still refusing.

### 4. Peat stock and vertical structure cannot be inferred from one equilibrium run

Peat formation is a competition between litter inputs, vertical decomposition,
water table, temperature, nutrients, fire/drainage and millennia of history.
Wania et al. show strong dependence on spin-up duration and peatland age;
Kleinen et al. use a 50-year persistence criterion for peat-forming area and
Holocene history for stock accumulation. The fixed 1.5 m active column neither
predicts total peat depth nor couples accumulated carbon back to hydraulics.

This repository intentionally has no fabricated history. It therefore cannot
claim a unique present peat inventory. It can model active-layer process fluxes
and carry a declared equilibrium/age/initial-stock bracket, or accept a future
history scenario as a separate artifact. A vertical peat extension must share
SDEC-3's C-N-P/redox coordinate, conserve all three elements as area and depth
change, and distinguish the methane-active column from deeper stored peat.

### 5. Wetland vegetation is an Earth taxon bundle, not a Vesper trait set

`wetlandpfts.ins` defines Earth boreal/tundra woody PFTs, a moss, and wet grasses
described by Carex, Eriophorum, Juncus and Typha. Fixed water-table thresholds,
snow/GDD limits, root distributions, aerenchyma, tiller geometry and peat
formation traits determine both vegetation and gas transport.

These belong in PCAR-5's covarying trait registry. The registry must distinguish
at least bog/fen nutrient strategies, peat moss, emergent aerenchymatous plants,
wet mineral-soil vegetation and non-vegetated inundation; couple roots and
nutrient acquisition to PLHY/SDEC; and retain the current Earth PFTs only as a
named transplanted-Earth bracket. Vesper traits must not be created by copying
Earth taxon values or tuning methane totals.

### 6. Methane production is a globally tuned respiration split

The code contains no methanogen substrate or electron-acceptor state. It uses
bulk CENTURY respiration as the carbon source, the fixed root profile as a
vertical allocator and constant CH4:CO2 ratios as production. Wania et al. made
the detailed ratio adjustable because anaerobic observations span a wide
range; their parameter set was fitted at seven northern sites. Spahni et al.'s
simplified ratio intentionally combines production, oxidation and transport
into an emission factor for broad Earth regions. The current fork's 0.027 value
is explicitly another global-emission fit.

The existing detailed and simple paths are defensible Earth structural bounds,
not Vesper central values and not latitude regimes. The central interface should
consume vertically resolved substrate, temperature, saturation/redox and
electron-acceptor state from SDEC-3/SDEC-8, expose production separately from
oxidation, and preserve the current carbon ledger. A reduced ratio model can
remain as a broad, named model-form bracket.

### 7. Gas transport imports Earth atmosphere, gravity, clock and phenology

Pressure and gravity enter bubble solubility and the hydrostatic ebullition
threshold. Atmospheric O2 and CH4 set the diffusive boundary. Wind controls
air-water exchange. Rotation and absolute seconds set velocities and substeps.
Living aerenchymatous biomass controls plant transport. Every one is presently
an Earth constant or an unresolved upstream state.

BIO-22/BIO-23/BIO-28 must supply absolute time, local pressure, gravity and the
declared atmospheric gas state; ExoPlaSim forcing can supply wind. PCAR-5 must
supply plant transport traits and phenology. The transport bracket must include
winter and thaw pulses, dead-tiller persistence, pressure/wind-driven
ebullition, and oxidation of bubbles crossing unsaturated layers. Kallingal et
al.'s LPJ-GUESS v4.1 calibration found strong equifinality among production,
oxidation and transport parameters and identifies missing wind, pressure,
winter emissions and peaks as material limitations. No single fitted pathway
partition should be privileged.

### 8. Wetlands alone do not close the methane surface budget

The current module emits only from `PEATLAND`. Spahni et al.'s global model had
to represent northern peatlands, inundated wetlands, wet mineral soils and dry
soil uptake separately. Curry's diffusion-reaction model shows that aerobic
near-surface soils are a material, moisture- and porosity-dependent CH4 sink.
Rosentreter et al. synthesize large and highly skewed lake, river, reservoir and
other aquatic sources and warn explicitly about ecosystem-area uncertainty and
double counting. Saunois et al. show that bottom-up natural-source sums already
exceed atmospheric top-down constraints in some budgets.

Vesper needs a complete labelled surface ledger: persistent peatland;
inundated/wet mineral soil; dry mineral-soil oxidation; lakes, rivers and
playas/open water; fire; and any omitted geological or biological source kept
as an external uncertainty rather than silently zeroed. Aquatic sources may
begin as broad area/temperature/productivity/transport brackets, but they must
use the same mutually exclusive surface fractions as finding 2.

### 9. Wetland and peat transitions are not generic land-cover change

LPJ can read changing `PEATLAND` fractions, but its transition machinery was
built for land use. Donor peatland vegetation can pass through harvest logic,
and the CNP fork's transfer structure moves C and N litter/SOM while omitting P
stocks and litter P in the receiving-pool code. It does not form peat, rewet a
column, preserve age/depth structure or release the correct C-N-P-water stocks
when a wetland contracts.

Expansion, contraction, drainage, rewetting and permafrost/thaw therefore need
explicit stock-conserving transition semantics. Peat fire stays under FIRE-6
and FIRE-7: `Patch::has_fires()` currently disables all fire on `PEATLAND`, a
safe dormant baseline but not a claim that Vesper peat cannot burn.

### 10. A surface flux does not determine atmospheric CH4

Atmospheric abundance depends on transport, OH and other sinks, and feedbacks
between CH4 and oxidizing capacity. Saunois et al. place atmospheric OH as the
dominant but uncertain sink and diagnose source-sink closure with atmospheric
growth and inversions, not by dividing a wetland flux by an Earth constant.
Vesper's K-star UV, water vapour, ozone, BVOCs, fire/soil NOx and circulation all
affect the lifetime.

The accepted seasonal surface ledger should feed BVOC-6's reduced
chemistry/transport bracket and return a global and, if resolved, latitudinal
CH4 state plus lifetime and imbalance. Only a converged result may replace the
fixed 1.6 ppmv radiative input. The climate--hydrology--wetland--CH4--oxidant--
climate loop needs hashes and an explicit convergence test; it cannot silently
mix a Vesper source with the fixed-modern-Earth-flux photochemical abundance.

### 11. Acceptance requires area, water, carbon and pathway diagnostics together

A plausible global CH4 total can arise from compensating wetland area,
production ratio, oxidation and transport errors. WETCHIMP and the LPJ-GUESS
MCMC study both demonstrate that aggregate agreement does not identify those
components. An accepted artifact must therefore retain native and aggregated
surface fractions; daily water table and routed/groundwater exchange; peat
depth/age bracket and C-N-P stocks; substrate production, oxidation, diffusion,
plant and ebullition fluxes; dry-soil uptake and aquatic sources; and water,
C-N-P and atmospheric methane residuals.

Pre-register matched no-methane, source-only, wetland-extent, production/redox,
transport and atmospheric-feedback brackets. Source and no-simulation
conservation/restart fixtures can proceed now. Any LPJ-GUESS, chemistry or
climate experiment with non-negligible compute still requires explicit
permission.

## Recommended order

1. Repair the free-water, annual-water-table and restart correctness defects
   while keeping peat and methane disabled.
2. Define mutually exclusive native-mesh surface fractions and one
   mass-conserving groundwater/surface-water exchange.
3. Register wetland vegetation and active-peat C-N-P traits, then establish the
   no-history peat-stock bracket.
4. Separate production, oxidation and transport; add dry-soil uptake and
   inland-water source brackets; emit one closed surface ledger.
5. Close that ledger through the reduced oxidant/lifetime model before changing
   atmospheric CH4 or climate.
6. Only then execute preregistered matched cases, with explicit permission.

## Sources read

- Wania, Ross and Prentice (2009a), *Integrating peatlands and permafrost into a
  dynamic global vegetation model: 1. Evaluation and sensitivity of physical
  land surface processes*, `10.1029/2008GB003412`.
- Wania, Ross and Prentice (2009b), *Integrating peatlands and permafrost into a
  dynamic global vegetation model: 2. Evaluation and sensitivity of vegetation
  and carbon cycle processes*, `10.1029/2008GB003413`.
- Wania et al. (2010), *Implementation and evaluation of a new methane model
  within a dynamic global vegetation model: LPJ-WHyMe v1.3.1*,
  `10.5194/gmd-3-565-2010`.
- Spahni et al. (2011), *Constraining global methane emissions and uptake by
  ecosystems*, `10.5194/bg-8-1643-2011`.
- Kleinen, Brovkin and Schuldt (2012), *A dynamic model of wetland extent and
  peat accumulation: results for the Holocene*, `10.5194/bg-9-235-2012`.
- Melton et al. (2013), *Present state of global wetland extent and wetland
  methane modelling: conclusions from a model inter-comparison project
  (WETCHIMP)*, `10.5194/bg-10-753-2013`.
- Kallingal et al. (2024), *Optimising CH4 simulations from the LPJ-GUESS model
  v4.1 using an adaptive Markov chain Monte Carlo algorithm*,
  `10.5194/gmd-17-2299-2024`.
- Curry (2007), *Modeling the soil consumption of atmospheric methane at the
  global scale*, `10.1029/2006GB002818`.
- Rosentreter et al. (2021), *Half of global methane emissions come from highly
  variable aquatic ecosystem sources*, `10.1038/s41561-021-00715-2`.
- Saunois et al. (2020), *The Global Methane Budget 2000-2017*,
  `10.5194/essd-12-1561-2020`.
