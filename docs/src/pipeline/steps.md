# 3. The pipeline, step by step

### 3.1 Geography

World Orogen generates the terrain from a planet code, which encodes the seed
and every slider, plus a carve list. The code is in
`source/worldorogen_seed.txt` and is not a build: it survives across builds,
while a build is what one pass of the generator produced from it.

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

Seven surface fields are always supplied: topography (129), land mask (172),
roughness (173), broadband and two-band albedo (174, 175, 176) and forest
fraction (212). Three more are conditional on model keys: soil water capacity
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

### 3.4 Hydrography

The export deliberately does not route water -- routing is a hydrology
decision and the exporter leaves it downstream. `build_hydrography.py`
resolves it with a priority flood that keeps the preserved basins as
terminals, rebuilds each basin's hypsometry on the *finished* terrain (the
catalogue's curves are measured on the pre-conditioning surface and overstate
capacity substantially), and writes the sparse basin-by-grid-cell coupling
matrix that climate is integrated over. `hydrography/README.md` has the
construction.

### 3.5 The carve verdict

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

### 3.6 Pedology

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

### 3.7 The biosphere

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

Forcing arrives as one binary driver file carrying however many years of
climate it was built with, and the model cycles through them. One year is a
fixed climate; several are how a variable star reaches the biosphere, and a
single repeating year cannot represent one at all.

### 3.8 Economic minerals

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
