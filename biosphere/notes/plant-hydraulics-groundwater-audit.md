# Plant hydraulics, drought and groundwater audit

This note records the second follow-up biosphere investigation: the active
LPJ-GUESS-CNP water-stress path from atmospheric demand through soil and roots
to vegetation, and its presently absent connection to Vesper's groundwater
model. Fire and soil biogeochemistry are outside scope except where water
stress changes mortality or nutrient uptake. No LPJ-GUESS build or simulation
was run.

## Scope and sources read

The source audit covered `modules/canexch.cpp`, `modules/soil.cpp`,
`modules/soilwater.cpp`, `modules/vegdynam.cpp`, the PFT/root state in
`framework/guess.h`, `modules/vesperinput.cpp`, `data/ins/global.ins`, the
Vesper driver, the pipeline graph, and the groundwater solver and artifact
contract. The primary sources below were fetched into `references/`, read, and
recorded in `references/INDEX.md`:

- Jackson et al. (1996), the Earth-biome root profile used by the active PFTs;
- Shah, Nachabe and Ross (2007), the source of the groundwater-ET extinction
  curve already used by Vesper's groundwater solver;
- Anderegg et al. (2016), cross-species hydraulic traits and drought mortality;
- Papastefanou et al. (2024), hydraulic-failure mortality in LPJ-GUESS-HYD;
- Meyer et al. (2025), the current LPJ-GUESS-HYD formulation, sensitivity and
  evaluation;
- Verbruggen et al. (2025), LPJ-GUESS-RE soil hydrology and aquifer boundary;
  and
- Jia et al. (2026), a two-way ParFlow--LPJ-GUESS coupling.

These studies establish model forms and failure modes. Their Earth species,
soil and aquifer parameters are not a Vesper calibration.

## Active Vesper configuration

`global.ins` selects the 15-layer soil, `wateruptake "rootdist"`, and Jackson
root distributions. Each layer is 100 mm, making a fixed 1.5 m active column.
The PFT-specific `root_beta` values are Earth-global calibrations; any root
fraction below 1.5 m is placed in the fifteenth layer. Drought-limited
establishment is off, and its common fallback `drought_tolerance` is deliberately
near zero. Wetland saturation is also off.

The Vesper adapter scales layer water capacities at the modelled
regolith--bedrock contact, including a capped bedrock-water fraction, but does
not change layer geometry or root distribution. The LPJ driver does not read
`water_table.nc`, and the pipeline contains no dependency between groundwater
and either `lpj_driver` or `lpj_run`.

## Findings

### 1. The active drought model is a supply--demand scalar, not plant hydraulics

For the selected `rootdist` uptake mode, each layer can supply the lesser of
its plant-available water volume and `emax * rootdist`. Until the available
volume becomes limiting, uptake does not decline continuously with water
content or soil water potential. Total supply is compared with atmospheric
demand; a shortage reduces canopy conductance and assimilation. The daily
supply/demand ratio becomes `wscal`, whose annual mean also controls subsequent
root-to-leaf allocation and raingreen phenology.

There are no soil, root, xylem or leaf water potentials; hydraulic resistance,
plant water storage, embolism/cavitation, hydraulic redistribution, adaptive
rooting or hydraulic-failure mortality. Drought can reduce growth and thereby
raise the existing five-orbit growth-efficiency mortality, but that delayed,
source-limited route is not equivalent to hydraulic mortality. Papastefanou et
al. demonstrate that adding hydraulic failure can decouple carbon losses from
carbon gains, while also warning that the result is sensitive to the mortality
mapping and omits groundwater uptake and redistribution. Anderegg et al. find
hydraulic safety traits predictive across species but explicitly do not offer a
first-principles mortality predictor.

This leaves a defensible empirical standard model, provided it is labelled and
bracketed, but not a mechanistic drought-survival claim. It is PLHY-1 and the
hydraulic model-form comparison is PLHY-6; drought mortality must be assessed
jointly with DEMO-4 rather than added as an independent certainty.

### 2. Atmospheric drought is not yet represented by the forcing contract

The current adapter supplies temperature, precipitation, shortwave radiation
and diurnal temperature range, but discards the climatology's pressure,
humidity, wind and net longwave radiation. BIO-23 already records the resulting
errors in equilibrium evapotranspiration, air pressure and the psychrometric
term. Meyer et al. show the practical distinction: standard LPJ-GUESS gives
species nearly identical evapotranspiration responses to vapor-pressure deficit,
where LPJ-GUESS-HYD uses VPD directly and resolves an isohydric--anisohydric
continuum. BIO-13's rainfall chronology is equally load-bearing because monthly
rain spread smoothly over days cannot produce the same dry spells.

The active canopy code also labels its natural-vegetation water-stress
conductance equation as faulty and says the crop-only fix has no natural-
vegetation equivalent. That upstream TODO is on the exact path used here and
must be resolved or bounded before interpreting drought response. Separately,
`ifdroughtlimitedestab 0` means the run admits no direct drought filter on
establishment; simply turning it on would not fix this because the generic PFTs
carry a near-zero placeholder rather than supported Vesper values.

Completing BIO-13 and BIO-23 is therefore a prerequisite, not an optional
refinement to hydraulics. PLHY-1 records the remaining standard-model decision,
regression and diagnostic work.

### 3. The soil bucket can preserve an artificial deep dry-season reservoir

The standard 15-layer scheme sets `soil.percolate` only on days with at least
0.1 mm of rain plus snowmelt. Consequently its percolation and base drainage
paths do not run through dry days. Verbruggen et al. identify this same default
condition as the cause of abrupt inter-layer gradients and water retained in
deep layers through a dry season, capable of supporting evergreen trees for
the wrong reason. The code has no upward capillary flow and no lower hydraulic
boundary.

LPJ-GUESS-RE replaces the bucket movement with water-potential-gradient flow,
allows upward movement, variable soil depth, and bedrock, free-drainage or
prescribed-aquifer bottom boundaries. It improved upper-layer soil moisture in
the tested drylands but did not solve the deeper profile, remained sensitive to
pedotransfer functions, and assumed vertically homogeneous hydraulic
properties. It also cost about five times the default in the authors' experience;
the aquifer boundary added another factor of ten and worsened mass-balance
error. Thus a wholesale RE merge is a useful bracket, not a cheap prerequisite
to every run. PLHY-2 first isolates the dry-day artifact and defines a lower-
boundary interface, then treats RE as a pre-registered model-form comparison.

### 4. Root depth, regolith depth and groundwater access are three different quantities

The Jackson equation describes cumulative root biomass in observed Earth
biomes. The implementation truncates it at 1.5 m and adds all residual mass to
the bottom layer, so that layer represents both 1.4--1.5 m roots and every root
the empirical curve placed deeper. Its `root_beta` values were recalibrated to
match an older two-layer parameterization, not Vesper plants. Verbruggen et al.
and Jia et al. independently identify fixed shallow root depth and static
profiles as limiting groundwater-dependent vegetation and phreatophytes.

Vesper has already improved one axis: actual regolith depth scales water
capacity. But root allocation remains unchanged through thin regolith, while
sub-bedrock capacity is forced into soil-texture layers and capped at the
overlying soil's capacity because higher values break the model's porosity
assumption. A numerical capacity floor also keeps nominal water in otherwise
non-storage layers. None of this expresses a taproot, weathered-rock porosity,
root plasticity, or a water table below 1.5 m.

PLHY-3 separates soil-profile depth, rootable depth, rooting strategy and
weathered-bedrock storage. It must remove the bottom-layer residual
conflation and register a deep-root/adaptive-root bracket rather than infer
alien root architecture directly from the water table.

### 5. Groundwater supply cannot be added to LPJ without removing the existing sink

Vesper's groundwater solve already removes
`ETmax * exp(-depth / lambda)` from recharge. Shah et al. support an exponential
rather than linear extinction curve, but show that its scale depends jointly
on soil and land cover and that very shallow groundwater connects saturated and
vadose domains. The current Vesper term deliberately collapses capillary soil
evaporation and plant uptake into one bulk sink. Giving LPJ an additional
groundwater supply while leaving that sink in place would spend the same water
twice.

`water_table.nc` also warns that `depth_m` is not uniformly a lateral-flow
prediction. Where `sink_fraction` approaches one, depth is principally the
local recharge--ET balance; where `at_surface` is set, the head is pinned.
Existing Earth tests find strong regional variation in cellwise skill. These
flags are validity diagnostics, not plant-availability parameters.

The first safe interface is therefore a flux contract: expose groundwater ET
at native hydrography resolution, partition the one conserved withdrawal among
bare-soil evaporation and PFT uptake, return unused demand to the groundwater
balance, and replace rather than supplement the solver's bulk sink. A later
soil-column bottom boundary can use depth where the artifact is certified, but
raw mean depth is not an acceptable shortcut. PLHY-4 defines this coupling;
GW-5 remains the separate climate-land-surface feedback of the same conserved
flux.

### 6. The spatial and iteration seams are as important as the local equation

Groundwater is solved on approximately 15 km mesh regions; LPJ is driven on a
roughly 300 km T42 cell. Root access and discharge vegetation depend on the
fraction of a cell with shallow groundwater, not its mean depth. Averaging the
depth can erase riparian and oasis habitat or manufacture broad access from a
narrow discharge zone. The existing coupling matrices provide the region-to-
cell mapping; GRID-2's mesh-under-cell hypsometry can supply terrain-conditioned
fractions, while GW-6 records the unresolved within-region terrain. PLHY-5
therefore requires native-mesh integration of accessible area and flux,
retaining `sink_fraction`, `at_surface` and uncertainty regimes through
aggregation.

The temporal graph is also open. Groundwater currently consumes bootstrap
climate in loop A, while soil and LPJ run later without groundwater; its
recharge does not come from LPJ drainage, and vegetation ET never returns to
the groundwater solve. Jia et al. provide the relevant architectural precedent:
they replace LPJ soil moisture with hydrologic-model state daily and return
LPJ precipitation minus ET as a boundary flux. Their result is explicitly
two-way and still suffers fixed roots and coarse-grid oversaturation. Vesper
does not need ParFlow to adopt the invariant: exactly one component owns each
water store and every exchanged flux is removed from its source. PLHY-4 must
close recharge, uptake/evaporation, seepage and storage over a declared
iteration rather than connect two steady-state artifacts one way.

### 7. LPJ-GUESS-HYD is the right structural bracket, not a drop-in truth

LPJ-GUESS-HYD adds plant water potentials, root/stem/leaf resistances,
vulnerability curves, VPD-coupled stomatal regulation and cavitation mortality.
Meyer et al. find improved drought ET but no meaningful GPP improvement over
standard LPJ-GUESS in their evaluation. Of seven new hydraulic parameters,
`psi50` and the maximum soil-to-leaf potential difference dominate much of the
response; the values are constrained for Earth species, not the generic Vesper
PFTs. The model still omits plant capacitance, turgor-limited growth, drought
leaf shedding, root-interface failure and downstream biotic mortality.

The scheme also explicitly includes the gravitational potential `rho * g * h`.
Unlike the active standard model, a port must obtain gravity from Vesper rather
than retain Earth's constant: at Vesper's higher gravity the lift cost of a tall
canopy is correspondingly larger. Plant height, vulnerability and conductance
therefore form a coupled trait prior; importing each independently would create
hydraulic strategies no organism exhibits.

PLHY-6 registers standard versus HYD as a no-retuning structural sensitivity,
with Vesper gravity and a coherent trait registry. PLHY-7 retains the daily
forcing, layer, uptake, groundwater-source, water-potential, cavitation and
mortality diagnostics needed to distinguish a real mechanism from compensating
errors. The comparison is not executable until its forcing, water balance,
soil/root and demographic prerequisites are closed, and it still requires
explicit permission to run LPJ-GUESS.

## Interpretation boundary and order

Until PLHY-4 closes, the active LPJ result is rainfall-fed equilibrium potential
natural vegetation; it does not predict groundwater-dependent vegetation,
oases, riparian belts or drought mortality. The safe order is to finish the
daily forcing and acceptance contracts, isolate the standard bucket and root
defects, define the conserved native-mesh groundwater flux exchange, and only
then compare standard versus RE and HYD model forms on pre-registered cases.
No raw water-table-depth feed and no additive groundwater ET route is physically
closed.

PLHY-1 through PLHY-7 record that work in `TASKS.md`. Source/fixture and
diagnostic plumbing can proceed without a model run; every measured sensitivity
or production execution remains subject to the explicit compute permission.
