# 3. The pipeline, step by step

### 3.1 Geography

World Orogen generates the terrain from a planet code, which encodes the seed
and every slider, plus a carve list. The code is in
`source/worldorogen_seed.txt` and is not a build: it survives across builds,
while a build is what one pass of the generator produced from it.

Before a spatial producer is meaningful, `spatial_support_gate` validates the
versioned vocabulary in `config/spatial_support.yaml`. A support is not a shape:
its identity includes geometry and coordinates, native area or volume, any
effective fraction, time support, field semantics, aggregation operator and
nonlinear order, and provenance. Native terrain regions, Gaussian atmosphere
cells, GOLDSTEIN cells, hydrologic units, soil units, separate land and marine
ecological response units, LPJ grid cells and demographic patches, and
comparison supports are distinct kinds. In particular,
an LPJ patch explicitly says it is a stochastic demographic sample and not a
spatial response unit. `lib/spatial_support.py` is the validator; SPAT-10 owns
applying it to every conversion and production artifact.

The aggregation vocabulary names one operator that does not reduce a cell to a
number. A field declared `distribution_quantiles` under
`area_weighted_distribution` carries the area-weighted quantile table
`lib/gridding.py:cell_quantiles` computes, and must state the probability
vector it is resolved at, because a table without its probabilities is a block
of numbers in the field's units and nothing more. The pair is enforced in both
directions: a distribution cannot declare a mean, and the operator cannot carry
a single-valued field. That is what lets a support state that it preserved a
within-cell distribution rather than having already thrown it away, which is
what the ocean's bathymetry needs and what conservation alone cannot see -- the
gate compares the vocabulary against the reductions `lib/gridding.py` exposes,
so a reduction that module gains without a term to name it fails
`spatial_support_gate`.

`spatial_conversion_gate` is that application envelope. It binds a conversion
to the source and destination contract digests, shapes and coordinates; carries
separate area, volume, water, salt, energy and C-N-P closure ledgers; and
requires constant, identity, reduction, round-trip, operator-order and
common-support tests. A quantity that does not apply still has a named reason.
Nearest fallbacks, changed surface ownership or connectivity, discarded modes,
and unmapped extensive stores are lists even when empty, so a finite field can
never hide what the operator dropped.

Verify a build by `manifest.hashes.finalElevation`. A seed alone does not
identify a planet, because fixes to the generator change the terrain under a
fixed seed. Which build is current is in `world_state.json`; `lib/orogen.py`
holds the registry, refuses any terrain hash it has not been checked against,
and notes what was wrong with each superseded build.

The fork exports Gaussian grids directly off the mesh at several truncations,
plus a uniform 512x256 for mapping, plus PNG maps; `source/README.md` has the
table. Only the T42 export carries `raw/`, the native mesh, and that mesh is
the same for every export.

### 3.2 Boundary conditions

Both the land mask and the topography are integrated from the native mesh, not
sampled from a gridded export. The export resamples categorical fields by the
region containing the cell centre, which at T42 discards the coastline, and
averages continuous fields over all regions in a cell, which drags a coastal
cell's albedo toward open water's 0.06.

Land comes from `surface_class`, never from `land_mask`. The two disagree over
dry closed-basin floor below sea level -- `manifest.landSeaMask` has the size
of the disagreement -- and `land_mask` would flood it (CLAUDE.md rule 1).

The atmosphere boundary receives both the binary ownership mask on code 172
and the native-mesh subaerial area share on code 1720. ExoPlaSim still converts
its legacy `yls` mask to zero or one at 0.5, so the physical exchange remains
binary until the tile kernel is complete. The immutable fraction is separately
bound as `dlf`, `ylf` and `xlf` in the atmosphere/land, ocean and sea-ice state
owners, and all three aliases are re-read on cold and restart starts. The
compiled `spat5_tile_combine` primitive bypasses absent-tile arithmetic at the
pure endpoints and preserves identical tile inputs bit-for-bit. Fourteen
restart-required records retain land/ocean temperature, humidity, wetness,
roughness and broadband/two-band albedo publications. Land/glacier and
ocean/ice now consume separate turbulent and radiative bundles, after which
the shared diagnostics are restored through that combine. Both bundles are
still seeded from the one binary evaluation; independent tile evaluation is
the remaining physics seam. `coastline_threshold_cost` measures
that rounding on all five exported rungs rather than reading its partly
cancelling net. On the active carve2 export, land dropped into ocean falls from
11.58% of mesh land at T21 to 3.67% at T170, while water promoted to land falls
from 10.56% to 2.93%.
`partial_surface_gate` checks the mutually exclusive class partition and the
declared measurement-before-model-change decision. `coastline_flux_bracket`
refuses until an accepted baseline exists on this build, then bounds the
missing tile rate from same-latitude observed cells. It is explicitly a
fixed-climate model-form bracket, not a coupled tile result. The downstream
`partial_surface_decision_gate` compares its surface-energy bound with the
accepted run's already-declared state-storage tolerance, verifies that code
1720 is generated, staged and re-read by every surface-state owner on every
start, and refuses until the selected exchange and separate tile states exist
in the model source.

Eight surface fields are always supplied: topography (129), land mask (172),
subaerial area fraction (1720), roughness (173), broadband and two-band albedo
(174, 175, 176) and forest fraction (212). Three more are conditional on model
keys: soil water capacity
(229, `model.soil_water_source`), the dust radiative source (1811,
`model.dust_source`), and the dust emission fields (1801,
`model.dust_emission`).

Roughness is measured from land cover and derived from this world's own
subgrid slope rather than asserted uniform, so its land mean is an output and
not an anchor; every other surface field falls back to a
declared uniform, because Earth's maps are tied to Earth's continents.
`exoplasim/README.md` has the construction.

Lakes enter through those same fields rather than through the mask: the solved
lake extent is composited per mesh region into the albedo, since almost every
lake on this planet is far below the grid, and it is one of the larger single
terms in the land-mean albedo because the cells carrying water are the bright
playa and salt crust -- the measured value is `lakes.delta` in
`exoplasim/inputs/t42/albedo_report.json`, which the error budget reads rather
than transcribes. Lakes deliberately get neither water's roughness (no water
column's heat capacity to go with it) nor a soil-water representation (routed
river water never re-enters the evaporating bucket).
`exoplasim/notes/lake-representation.md` has both arguments.

### 3.3 Climate

ExoPlaSim 3.4.2, on the resolution, rank count, timestep, layer count and
slab-ocean depth that `config/planet.yaml`'s `model:` block declares, with
interactive sea ice, glaciers enabled, and a measured stellar spectrum
rather than a blackbody.

`ocean_column_response` records why the dormant multilayer slab is not the
route beyond that surface tier. Its vertical solver is cheap and viable, but
the source overwrites the final declared layer with `mldepth`, applies every
surface watt to level 1, and at the default diffusivity raises the unit seasonal
SST response by 2.79 times for a converged ten-layer subdivision of the same
50 m. That is a material, profile-dependent model change without penetrative
shortwave. The active resolved tier remains the separately contracted
cGENIE/GOLDSTEIN loop.

`ocean_shortwave_penetration` prices the companion optical change without
activating it. Hale--Querry extinction integrated against k25v puts 69.53% of
the clear-water shortwave in the top metre, against 57.73% for a 5772 K Planck
spectrum on the same grid, and leaves 3.19% below a hypothetical 50 m column.
The operator closes energy and exactly reduces to existing surface deposition,
but OCN-1 left it no accepted ExoPlaSim layer destination. Biological
attenuation stays with OCN-14 and reflected ocean colour with OCN-7.

`ocean_carbon_gate` keeps the marine carbon feedback deliberately absent while
`config/planet.yaml` prescribes atmospheric CO2. It checks every ExoPlaSim run
entry point, refuses a carbon-bearing ocean return, and records the ownership,
conservation, permission and convergence requirements that must replace this
boundary before air-sea CO2 exchange can enter a new pipeline loop. It runs no
model and does not treat a diagnostic ocean carbon inventory as an atmospheric
reservoir.

That spectrum is `k25v`, built from BT-Settl, replacing a shipped file that
was the star K2-18 rather than the spectral type K2; snow, ice and glacier
albedos were 0.10 to 0.17 too dark under it, which weakened the very feedback
the glacier and stellar-cycle machinery exists to resolve. Everything before
the re-baseline carries that bias, direction known; state it wherever those
numbers are quoted. See `exoplasim/notes/stellar-spectrum-audit.md`.

A run is spun up cheap and finished clean: spin-up segments carry PlaSim's
low-I/O accumulation and the orbits a climatology is built from do not,
because an accumulation cannot be undone while instantaneous records can be
averaged and still answer questions about variance and extremes. The clean
regime is the default, so a forgotten flag costs time rather than data, and a
climatology from low-I/O orbits is refused rather than warned about.
`exoplasim/README.md` has the invocation and
`exoplasim/notes/first-output-bin.md` the evidence.

A run directory names everything that changes the answer -- resolution, flux,
CO2, rotation, obliquity, eccentricity, the physics switches, the spectrum,
and a digest of every surface input. Anything physical that is not in the name
is a collision waiting for the run that changes it; the spectrum and the flux
precision each nearly caused one before being added.

`model.energy_diagnostics` adds PlaSim's 28-term energy decomposition on codes
360-387, registered at run time because the postprocessor does not ship them
and patching the vendored tree is what a reinstall silently undoes. It needs a
settled run: a segment taken one orbit off a restart sits several W/m2 out of
balance, because `configure()` resets surface fields a restart does not carry.

Convergence is a fixed six-part test, not a judgement; the criteria and the
derivation of each threshold are in `assess_convergence.py`, and runs that
miss are labelled, not rounded. The test window is the last `--window`
PRODUCTION orbits: segments declare what they were for, and orbits run to
measure the model rather than the planet are dropped from the tail and refused
inside the window.

**A window that cannot resolve a threshold is not testing it.** The model has
variability at the window's own timescale, so consecutive orbits are not
independent samples and a criterion whose statistic has a standard error
comparable to its own threshold flips rather than discriminates. Every
assessment reports what its window can resolve and the window each criterion
would need at that run's own scatter and autocorrelation, so the number is
produced rather than swept for; `exoplasim/notes/convergence-lengths.md`
carries the arithmetic, the bracket it lands in, and the operational lengths of
the approach itself. A separate five-orbit window with 32 snapshots per orbit,
run after the pass/fail decision, forms the climatology.

### 3.4 Offline ocean circulation

The adopted circulation is cGENIE/GOLDSTEIN, wrapped by `ocean/`; ExoPlaSim is
still the only atmosphere and remains authoritative for sea ice. The first
ocean pass begins from an explicitly zero transport return. Later passes build
one chronological atmosphere/land/ice forcing bundle, integrate the ocean,
assess its heat and freshwater/salt closure, and conservatively return ocean
heat convergence on surface code 903 plus surface velocity for sea-ice
advection.

Two contracts stand between a climatology and a run. OCN-10 owns time bounds,
units, signs, masks, state-versus-flux semantics, and one owner for every heat,
water and salt term. OCN-11 owns the wet area and volume, partial coasts,
bathymetry and connectivity, river mouths, and mappings between the independent
ocean grid and every atmosphere rung. Their registered entry points refuse
until their Beads dependencies are accepted; the pipeline therefore cannot
mistake a provisional file for a valid coupling artifact.

The forcing procedure, provenance chain and across-pass exit are declared in
`ocean/config/transport_loop.yaml` and checked by
`ocean/scripts/transport_loop_gate.py`. The executed return channel has a
separate proof in `exoplasim/scripts/verify_ocean_flux_channel.py`; the static
gate does not substitute for that model run.

### 3.5 Hydrography

The export deliberately does not route water -- routing is a hydrology
decision and the exporter leaves it downstream. `build_hydrography.py`
resolves it with a priority flood that keeps the preserved basins as
terminals, rebuilds each basin's hypsometry on the *finished* terrain (the
catalogue's curves are measured on the pre-conditioning surface and overstate
capacity substantially), and writes the sparse basin-by-grid-cell coupling
matrix that climate is integrated over. `hydrography/README.md` has the
construction.

### 3.6 The carve verdict

A basin that overflows year on year incises its outlet, and over 1e4 to 1e6
years that drains the lake and the depression stops existing. So a basin
pinned at its spill is a transient, not a landscape state.

The test is geometry on one side and climate on the other:

    (E - P) / runoff  <=  catchment / area_at_spill - 1

The right-hand side ships in `basins.nc` as `critical_aridity_index`. The
left-hand side is the climatology integrated over each catchment.

`E` is evaporation from open water, which the model does not have there. It
is estimated with the Penman combination equation using water's albedo and
roughness, and validated by applying the same calculation to ocean cells,
which *are* open water; that validation ratio is computed on every run and
written into `carve_verdict.json`. **It is not floored at the model's land
evaporation, and must not be**: this world's land is aerodynamically rougher
than open water, so a smooth lake in a rough wet landscape evaporates less
than the ground around it, and a floor asserting otherwise once decided most
of the overflowing basins by clamp rather than by climate.

### 3.7 Pedology

`pedology/` weathers the lithology into soil under the climate, which is the
step neither the climate model nor the vegetation model does. Parent material
sets which minerals are available, climate sets how far they have been
converted, relief and erosion set how much regolith survives, and the
biosphere sets the organic fraction.

The weathering law, the texture model and their calibrations are
`pedology/README.md`'s to describe. What crosses component boundaries is the
convention: **runoff means `P - E`, not the model's `mrro`.** `mrro` is
river-routed net water flux, negative in places and non-zero over ocean;
using it understates land runoff by 6.6x. See
`exoplasim/notes/water-and-energy-closure.md`.

**Volcanism enters here as a process, not only as a rock class.** Andic soil
properties need ongoing ejecta, because volcanic glass is metastable and a
surface that stops receiving ash weathers past it. The arc classes are the
only place in this pipeline where volcanism is known to be *continuing*; flood
basalt and ocean-island provinces are reported undetermined for want of an
eruption-age field rather than counted or assumed absent. The second gate is
leaching, at the literature's own boundary: below about 1500 mm of
precipitation the same parent material weathers to halloysite, an ordinary
clay soil. The consequence on this world is that most volcanic terrain is
*not* andisol, and the reason matters -- andic material FIXES phosphorus
rather than supplying it, so getting this wrong would have inverted a nutrient
result rather than merely scaling one.

**The same weathering also produces solute fluxes, and they leave the
component.** `weathering_fluxes.py` computes CO2 consumption and dissolved
silica per lithology from concentrations times runoff. Silica is the supply
side of silcrete and diatomite, and the endorheic share is what matters there:
silica reaching the ocean is diluted into an enormous reservoir, while silica
reaching a closed basin concentrates until it saturates. CO2 is the carbon
cycle, and [section 4](loops.md) says why that loop is left open.

**It also derives what the surface has BECOME, which is a different question
from what the profile is.** `build_surface_classes.py` emits two independent
axes -- `surface_cover`, what wind and light see, and `duricrust`, what is
cementing at or below the surface -- because the two do not compete for the
same physical position, and resolving them in one precedence chain once cost a
build its basin fill. Two of its inputs are deliberately external: the dust
deposition field decides loess and pavement (the aeolian roughness bracket is
worth a factor of 40 there, so the answer is a bracket), and the per-basin
chemical divide decides which duricrust can form at all.
`pedology/notes/derived-surface-classes.md` has the design.

Every Earth calibration lives in `pedology/config/` with its source or an
explicit statement that it is declared; nothing in the scripts hardcodes a
Vesper number, so porting the component to another world is a config change.

### 3.8 The biosphere

`biosphere/` runs the vendored LPJ-GUESS CNP v1.0 fork, with phosphorus
limitation deliberately disabled until its Vesper inputs are ready. The in-tree
port adapts the model to this world's calendar and astronomy; it carries no
planetary numbers itself but reads a generated `vesper.h` derived from
`config/planet.yaml`. The calendar decisions -- 24-hour model days, rescaled
degree-day limits, and why -- are `biosphere/README.md`'s.
What the pipeline has to know is the LOCK: the year is compiled into the
model, and the year is a function of the semimajor axis, so the axis has to be
locked before the biosphere is built against it and later flux changes go
through luminosity instead (loop C).

LPJ's spatial population is the BIO-11 rootable surface, not every binary land
cell. `build_rootable_fraction.py` partitions native-mesh land into solved
water, dry barren substrate, and rootable ground on each atmosphere rung. Fully
non-rootable cells are omitted from the driver; partial cells retain one
environment but every extensive NPP/carbon total, soil-carbon feedback and
prediction denominator uses the same effective rootable area. The artifact is
support- and build-stamped, and a model-land cell with no native support is a
SPAT-5 refusal rather than a nearest-cell fill.

Forcing arrives as one binary driver file carrying however many years of
climate it was built with, and the model cycles through them. One year is a
fixed climate; several are how a variable star reaches the biosphere, and a
single repeating year cannot represent one at all.

The stochastic model has one declared root in
`biosphere/config/stochastic_seeds.yaml`. `run_lpj_guess.py` writes that root
into both the instruction file and run manifest. LPJ-GUESS then derives a
separate restart-serialized Park-Miller state from cell coordinates, stand,
replicate patch and process. MPI rank and traversal order are absent from the
key, and establishment, fire occurrence, fire mortality, background mortality
and disturbance cannot advance one another's sequences. The
`stochastic_seed_gate` checks fixed vectors and rank/order fixtures without
launching the model.

Downstream consumers use one `lib/lpj_output.py` equilibrium statistic rather
than choosing the greatest year number. The declared window is ten complete
forcing cycles—ten years under a repeated one-year climate. It is accepted only
when every cell has every annual row and the sequence of cycle means is no
longer trending. Reports separate annual temporal spread, fixed-patch root-seed
spread, and across-patch-count spread; an absent ensemble is labelled
`not_measured`, not silently treated as zero.

For albedo feedback those equilibrium FPCs are conditional on rootable ground.
The modelled compositor replaces only BIO-11's rootable substrate contribution;
the solved-water and dry-barren contributions retain their native-mesh weights,
and forest code 212 receives the whole-cell rootable tree share. This is why a
partly flooded cell is not equivalent to multiplying its already mixed albedo
by one minus canopy cover.

### 3.9 Economic minerals

`minerals/` emits ore prospectivity per deposit type, in two artifacts split
by genesis and therefore by lifetime: the tectonic and magmatic half is a pure
function of the export and survives a re-run of the baseline; the weathering,
drainage and brine half reads a climatology, a lake solution and the per-basin
chemical divide, so it carries a climate in its identity as well as a terrain.
Design in `docs/src/reference/economic-minerals.md`; mechanics in `minerals/README.md`.
Two structural rules matter to the rest of the pipeline, and both are
expensive to retrofit.

**Minerals are an overlay, never a lithology.** `substrate_class` sets
erodibility, bare-rock albedo, soil texture and solute chemistry, and a
deposit is orders of magnitude smaller than a mesh cell, so admitting ore as a
rock class would move the planet's energy balance on the strength of a mine.
The layer is prospectivity per cell -- it does NOT place deposits, which
belong with the downscaling pass -- read only by what comes after climate, and
consequently rebuildable without invalidating a climate run. Even the
tectonic half is COMPUTED downstream, from the export rather than inside the
generator, so a value that does not exist during generation cannot feed back
into erodibility or albedo: the overlay rule is structural instead of
remembered.

**The split is by genesis** because whatever concentrates a deposit has to be
MODELLED wherever the deposit is placed: arcs, fold belts, LIPs and cratons
are Orogen's; climate, runoff and the chemical divide are not. Supergene
enrichment needs both, and is the case that shows why this is a split rather
than a handover.
