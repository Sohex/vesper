# Vesper: the worldbuilding pipeline

How geography, climate, water and life are computed for this planet, in what
order, and why the order is not a straight line.

**This document is CANONICAL for the pipeline's REASONING**: what each stage is
for, why the order is what it is, why the loops close, and the traps. The GRAPH
itself -- every step, what it writes, what must precede it -- is
`config/pipeline.yaml`, and the two do not overlap. This document references step
ids and does not restate them, because one graph written in two places is the
duplication both files exist to prevent.

The rule that follows, and it is a rule rather than an aspiration: **an artifact
that no step in `config/pipeline.yaml` generates does not exist.** If a script
writes something that is in no step, either the step is missing or the product is
one nobody should be reading. Both are defects, and the second is worse, because
an unregistered artifact has no derivable consumers and so no way to know what it
invalidates when it changes. `smoke_test.py` fails when a generator is absent
from the graph. `CLAUDE.md` rule 7 depends on all of this: "what is now
worthless" is only answerable from a graph that is complete.

Vesper is a super-Earth orbiting a K2.5V dwarf: larger than Earth, higher
gravity, a longer day, and more obliquity. Every one of those is declared in
`config/planet.yaml`, which is the only place they are written down. They are
stable rather than immutable -- the gravity has been corrected once, and the
prose copies of it that were scattered through this repo all went stale at that
moment while continuing to read as current.

The flux, the year, the active build and the mean surface temperature do move,
every iteration. They are in `world_state.json`, which is generated; see section
5.

---

## 0. Vocabulary, because two of these collide

Written 2026-08-17 after "re-baseline" was used to mean one thing and read as
another, and the two differ by roughly a day of work and a whole regeneration of
the terrain. These are the project's words; use them and not synonyms.

**build** -- one World Orogen export, namespaced under `source/`. Identified by
its terrain hash, never by its name. Superseded means wrong, not merely old.

**iteration** -- one turn of loop A in section 6: carve list to Orogen, terrain
regenerated, hydrography and boundary conditions rebuilt on it, climate re-run,
verdict taken again. An iteration produces a NEW BUILD. This is the expensive
thing, and it is what "starting again from Orogen" means.

**bootstrap run** -- the first climate run on a new terrain, made with only the
surface fields that are pure functions of terrain. It exists to produce the
climatology that the remaining fields need. **Its numbers are not the baseline.**

**baseline run** -- the second climate run of an iteration, on the full surface
fields: lakes, the lake compositing in the albedo, and soil water capacity, all
built from the bootstrap's climatology. The climatology of a baseline run is what
downstream components read.

**re-run the baseline** -- a NEW baseline run on the SAME terrain, because
something changed that is not the terrain: a model patch, a boundary field, a
configuration value. This is not an iteration and Orogen is not involved.

  Do not say "re-baseline". It reads as an iteration and means this, and that
  gap has already cost one misunderstanding.

**segment** -- a contiguous block of orbits added to an existing run by
`continue_exoplasim.py`. Runs are made of segments; each records its I/O regime
and its purpose, and the purpose is DECLARED by the caller with `--purpose`
rather than inferred from the flags. `spinup` and `post_equilibrium_climatology`
are the run's own trajectory; `diagnostic` is orbits run to measure the model
rather than the planet, and those are kept out of convergence windows and
climatologies. See `exoplasim/scripts/segments.py`.

**carve verdict** -- the finding: which basins overflow, per basin, with its
evidence. `hydrography/analysis/carve_verdict.json`.

**carve list** -- the artifact Orogen consumes, `carve_list.txt`, one retain
fraction per basin. The verdict is a conclusion; the list is an instruction.

### What a re-run of the baseline actually costs

The distinction is worth a table, because most of the pipeline does not move when
only the climate does:

| | redone? | why |
| --- | --- | --- |
| Orogen, `source/` | no | the terrain has not changed |
| drainage, basins, coupling | no | pure functions of terrain |
| land mask, topography, roughness | no | pure functions of terrain |
| model binaries | only if it was a patch | see CLAUDE.md rule 4 |
| the climate run | yes | that is the point |
| climatology | yes | follows the run |
| lakes, albedo 174-176, soil, code 229 | yes | all are built from a climatology |
| a second run on those rebuilt fields | usually | this is where the loop closes |
| carve verdict, dust, everything downstream | yes | below the climatology |

At the measured 88 s of model time per orbit, a settling run off an existing
near-equilibrium plus ten clean orbits is a couple of hours of wall clock, and
the loop usually wants two of them. An ITERATION is far more, because it
regenerates the terrain and everything that is a function of it.

## 1. The components

```
config/planet.yaml   Canonical planet, star, orbit, atmosphere. Every component reads it.
source/              World Orogen exports. Canonical, read-only.
lib/                 Shared readers: orogen.py (the export), gridding.py (mesh to grid).
hydrography/         Drainage, catchments, basin capacity, lake balance, carve verdict.
exoplasim/           Boundary conditions, climate integrations, climatology.
pedology/            Weathers lithology into soil, and into solute fluxes: CO2, silica, phosphorus.
biosphere/           LPJ-GUESS: vegetation, leaf area, carbon, PFT composition.
minerals/            Ore prospectivity per deposit type. Reads the export; never feeds climate.
aeolian/             Offline dust: emission, transport, deposition, optical depth.
                     Reads a climatology and the lake solution; feeds pedology's
                     loess and the radiation's dust forcing. NOT a leaf: it is
                     downstream of one climate and upstream of the next.
analysis/            Project-level products that belong to no single component:
                     the error budget, the dust optics and forcing.
maps/                Rendering. Terminal: nothing reads its output. Its generators
                     live in the component root rather than a scripts/ directory,
                     which is why an audit that globbed */scripts missed them.
```

Each component owns its own `scripts/`, `data/` or `analysis/`, and reads
`config/planet.yaml` and `source/`. Paths resolve from the file location, not the
working directory, so anything runs from anywhere.

## 2. The data flow

```
World Orogen (fork)
   |  seed + planet code -> terrain, lithology, closed basins, hydrology
   v
source/<build>/exoplasim-T21|T42|T85 + grid-512x256 + maps
   |  raw/ mesh lives in the T42 export and is identical for all of them
   v
lib/gridding.py            integrate mesh fields onto any model grid
   |                        |
   |                        v
   |                  exoplasim/build_boundary_conditions.py -> land mask, topography
   |                  exoplasim/build_surface_albedo.py      -> albedo, forest fraction
   |                  exoplasim/build_surface_soil_water.py  -> dwmax  (off by default)
   |                        |
   |                        v
   |                  ExoPlaSim spin-up
   |                        |
   |                  build_climatology.py -> averaged climatology
   |                                       -> per-orbit climatologies (--per-year)
   |                                       -> climate series, per bin per orbit
   |                        |
   |          +-------------+-------------+
   |          |                           |
   v          v                           v
hydrography/            pedology/build_soil.py        hydrography/carve_verdict.py
  drainage,               weathers lithology            which basins overflow
  catchments,             under the climate         <-- integrates climate over
  hypsometry,               |                           catchments
  coupling matrix           v
   |                  biosphere/build_lpj_driver.py -> one binary, N years
   |                        |
   |                        v
   |                  biosphere/run_lpj_guess.py  (LPJ-GUESS, MPI)
   |                        |
   |                        +--> cpool.out ---> back to build_soil.py
   |                        |                   (soil and biosphere iterate)
   |                        v
   |                  fpc.out -> build_surface_albedo.py --mode modelled
   |                        |     albedo and forest fraction from what grew
   |                        v
   |                  back to ExoPlaSim
   v
carve list -> back to World Orogen -> new terrain
```

Three loops close in that diagram and a fourth is cut across iterations;
section 4 says why each has to be what it is.

**Two branches hang off it that the diagram does not draw**, because they would
turn one picture into four. Both are in the register below.

```
climatology + surface_water.nc -> aeolian/build_dust.py -> dust_baseline.{nc,json}
   |                                        |
   |   deposition -> pedology (loess, phosphorus)
   |   optical depth -> exoplasim/dust_optics.py -> dust_aerofile.py -> the model
   |                 -> exoplasim/dust_forcing.py -> carve_verdict.py --dust
   v
source/ + soil + drainage -> minerals/build_prospectivity.py            (terrain only)
                          -> minerals/build_downstream_prospectivity.py (carries a climate)
```

## 2b. The artifact register

**The register is `config/pipeline.yaml`, and it is not reproduced here.** That
file declares every step, what it writes, and what must precede it; this document
explains why the order is what it is and what each stage is for. One graph in one
place, referenced from the other, because a table restating it would go stale the
way every hand-maintained list of current values in this project has -- and this
section used to be exactly that table.

    python scripts/pipeline.py --register     # artifact, step, and what reads it
    python scripts/pipeline.py --status       # what is present, and the carve gate
    python scripts/pipeline.py --plan <step>  # the ordered steps to reach a target

Consumers are derived rather than declared: what reads an artifact is the set of
steps that `need` the step which writes it. Declaring both would be the same
duplication one level down.

`pipeline.py` plans and never runs. Five steps in the graph cost hours, and a
script that could start one by accident is worse than no script. It also does not
check whether artifacts AGREE -- that is `check_consistency.py`, and asking
either to do the other's job would give two answers to one question.

## 3. The pipeline, step by step

### 3.1 Geography

World Orogen generates the terrain from a planet code, which encodes the seed and
every slider, plus a carve list. The code is in `source/worldorogen_seed.txt` and
is not a build: it survives across builds, while a build is what one pass of the
generator produced from it. The exports carry the seed and the mesh region count
in each `manifest.json`.

Verify a build by `manifest.hashes.finalElevation`. A seed alone does not
identify a planet, because fixes to the generator change the terrain under a
fixed seed. Which build is current is in `world_state.json`; `lib/orogen.py`
holds the registry, and it refuses any terrain hash it has not been checked
against, with a note on each superseded build saying what was wrong with it.
Several builds are superseded rather than merely older, and they stay registered
so results computed from them remain readable and datable.

The fork exports Gaussian grids directly off the mesh at T21, T42 and T85, plus a
uniform 512x256 for mapping, plus PNG maps. Only the T42 export carries `raw/`,
the native mesh, and that mesh is the same for every export.

### 3.2 Boundary conditions

Both the land mask and the topography are integrated from the native mesh, not
sampled from a gridded export. The export resamples categorical fields by the
region containing the cell centre, which at T42 discards the coastline, and
averages continuous fields over all regions in a cell, which drags a coastal
cell's albedo toward open water's 0.06.

Land comes from `surface_class`, never from `land_mask`. The two disagree over
dry closed-basin floor below sea level -- `manifest.landSeaMask` has the size of
the disagreement -- and `land_mask`
would flood it.

Seven surface fields are supplied. Topography (129), land mask (172), roughness
(173), broadband and two-band albedo (174, 175, 176) and forest fraction (212),
with soil water capacity (229) when `model.soil_water_source` is set.

Roughness replaces a uniform `dz0land = 2.0 m` that asserted forest-scale
roughness over salt crust and playa, which are closed-basin floors and flat by
construction. Its land mean is anchored to that same 2.0 m, so the global value
the model was tuned against does not move and only the distribution does. The
orographic half is measured as the standard deviation of elevation among the
~610 mesh regions inside each cell, which is what a 2.5M-region mesh is for and
what a gridded elevation field cannot give. Everything else falls back
to a uniform namelist default, which is declared rather than accidental: Earth's
roughness and vegetation maps are tied to Earth's continents and would be
meaningless here.

Lakes enter through those same fields rather than through the mask. The solved
lake extent is applied per mesh region and then integrated, so each cell receives
an area-weighted composite -- almost every lake on this planet is far below the
grid. It is worth -0.0139 on land-mean albedo, since the cells carrying water are
the bright playa and salt crust. Roughness (173) *is* supplied, but from land
cover and subgrid relief rather than from lake extent: giving a lake cell water's
roughness without a water column's heat capacity would cut turbulent exchange 10x
and leave it decoupled and hot. And no setting of 229 can
sustain a lake, because routed river water never re-enters the evaporating
bucket. See `exoplasim/notes/lake-representation.md`.

### 3.3 Climate

ExoPlaSim 3.4.2, T42, 16 MPI ranks, 45-minute timestep, 10 layers, 50 m slab
ocean, interactive sea ice, glaciers enabled, and a measured stellar spectrum
rather than a blackbody.

That spectrum is `k25v`, built from BT-Settl and interpolated to 4965 K, and it
replaces a wrong one that every run of the first three eras inherited.
ExoPlaSim's `k2.dat` is the star K2-18, an M2.5V at about 3450 K, not the
spectral type K2; the package ships no K dwarf spectrum at all. Snow, ice and
glacier albedos were consequently 0.10 to 0.17 too dark everywhere, which
weakened the very feedback the glacier and stellar-cycle machinery exists to
resolve. The model confirms the fix from its own log: the energy fraction below
0.75 microns is more than three times larger under `k25v` than the 0.11588 that
`k2` gives, and `lib/stellar.py` is where the current share comes from.

Everything before the re-baseline carries that bias. It is not invalidated in
kind and the direction is known, but it should be stated wherever it is quoted.
See `exoplasim/notes/stellar-spectrum-audit.md`.

A run is spun up cheap and finished clean. Spin-up segments carry PlaSim's
low-I/O accumulation and the orbits a climatology is built from do not, because
the two write different things: an accumulation cannot be undone, while
instantaneous records can be averaged into the same twelve bins and still answer
questions about variance, extremes and single records. Nothing downstream changes
shape either way. The clean regime is the default and the cheap one is opt-in, so
a forgotten flag costs time rather than data, and asking for a climatology from
low-I/O orbits is refused rather than warned about. `exoplasim/README.md` has the
invocation and `exoplasim/notes/first-output-bin.md` the evidence.

A run directory names everything that changes the answer: resolution, flux, CO2,
rotation, obliquity, eccentricity, the physics switches, the spectrum, and a
digest of every surface input. Two of those were added after they had already
nearly caused a collision. The spectrum marker was missing when the k25v
re-baseline would have landed in the completed k2 run's directory, and the flux
was rounded to hundredths, so 0.945 and 0.94 both resolved to `s094`. Only the
geography digest separated them, which was luck. Anything physical that is not in
the directory name is a collision waiting for the run that changes it.

Enabling `model.energy_diagnostics` adds PlaSim's 28-term energy decomposition on
codes 360-387, which was the instrument for the residual between the top of the
atmosphere and the surface. It closed the decomposition and thereby ruled itself
out: the residual is an offset in the reported top-of-atmosphere net radiation,
which the 28 terms do not see. The postprocessor does
not ship those codes, so `run_exoplasim.py` registers them at run time rather
than patching the vendored tree, which a reinstall would silently undo. The
residual needs a settled run: a segment taken one orbit off a restart sits several
W/m2 out of balance, because `configure()` resets the surface fields a restart
does not carry.

Convergence is a fixed six-part test, not a judgement: temperature drift below
0.05 K per orbit, top-of-atmosphere and surface balance trends below 0.05 W/m2
per orbit, sea-ice drift below 0.001 of planetary area per orbit, energy storage
in the prognostic state below 0.12 W/m2, and extrapolated remaining approach
below 0.15 K. The exact criteria and the derivation of each threshold are in
`assess_convergence.py`. Runs that miss are labelled, not rounded.

The test window is the last `--window` PRODUCTION orbits. Segments declare what
they were for, and orbits run to measure the model rather than the planet are
dropped from the tail and refused inside the window.

A separate five-orbit window with 32 snapshots per orbit, run after the pass/fail
decision, forms the climatology.

### 3.4 Hydrography

The export deliberately does not route water: `drain_to` is raw steepest descent,
and on this planet most of the land drains into unpreserved single-cell pits --
the share is per build and is in `world_state.json`. Routing
is a hydrology decision and the exporter leaves it downstream.

`build_hydrography.py` resolves it with a priority flood over the 2.5M-region
mesh, which fills the noise pits while keeping the preserved basins as
terminals. It then rebuilds each basin's hypsometry on the *finished* terrain,
because the catalogue's curves are measured on the pre-conditioning surface and
overstate capacity substantially; the measured ratio is in the hydrography report
for the build it was computed on. Finally it writes a sparse basin-by-grid-cell
coupling matrix, which is the interface climate is integrated over.

### 3.5 The carve verdict

A basin that overflows year on year incises its outlet, and over 1e4 to 1e6 years
that drains the lake and the depression stops existing. So a basin pinned at its
spill is a transient, not a landscape state.

The test is geometry on one side and climate on the other:

    (E - P) / runoff  <=  catchment / area_at_spill - 1

The right-hand side ships in `basins.nc` as `critical_aridity_index`. The
left-hand side is the climatology integrated over each catchment.

`E` is evaporation from open water, which the model does not have there. It is
estimated with the Penman combination equation using water's albedo and roughness,
evaluated wholly at the lowest model level, and validated by applying the same
calculation to ocean cells, which *are* open water and which the model also gives
water's roughness. That ratio is computed on every run and written into
`carve_verdict.json`; it is not quoted here, because the two times it was quoted
in prose it was quoting a hardcoded constant.

**It is not floored at the model's land evaporation, and must not be.** This
world's land is several times aerodynamically rougher than open water, so a
smooth lake in a rough wet landscape evaporates less than the ground around it.
A floor asserting otherwise stood for a day and decided 73% of the overflowing
basins by clamp rather than by climate.

### 3.6 Pedology

`pedology/` weathers the lithology into soil under the climate, which is the step
neither the climate model nor the vegetation model does. Parent material sets
which minerals are available, climate sets how far they have been converted,
relief and erosion set how much regolith survives, and the biosphere sets the
organic fraction.

Weathering intensity is the Walker-Hays-Kasting form, a power of runoff times an
exponential of temperature, normalised so Earth's land mean is 1. Texture comes
from mixing parent materials and converting weatherable primary minerals to clay,
with quartz tracked separately because it never becomes clay: that single fact is
why granite and basalt diverge under identical climate.

Runoff means `P - E`, not the model's `mrro`. `mrro` is river-routed net water
flux, so it is negative in places and non-zero over ocean; using it understates
land runoff by 6.6x. See `exoplasim/notes/water-and-energy-closure.md`.

**Volcanism enters here as a process, not only as a rock class.** Andic soil
properties need ongoing ejecta, because volcanic glass is metastable and a
surface that stops receiving ash weathers past it. The arc classes are the only
place in this pipeline where volcanism is known to be *continuing*, since they
are placed as a band about the volcanic front; flood basalt and ocean-island
provinces are reported undetermined for want of an eruption-age field rather
than counted or assumed absent. The second gate is leaching, at the literature's
own boundary: below about 1500 mm of precipitation the same parent material
weathers to halloysite instead, which is an ordinary clay soil. The consequence
on this world is that most volcanic terrain is *not* andisol, and the reason
matters -- andic material FIXES phosphorus rather than supplying it, so getting
this wrong would have inverted a nutrient result rather than merely scaling one.

**The same weathering also produces solute fluxes, and they leave the
component.** `weathering_fluxes.py` computes CO2 consumption and dissolved silica
per lithology from concentrations times runoff. Silica is the supply side of
silcrete and diatomite, and the endorheic share is what matters there: silica
reaching the ocean is diluted into an enormous reservoir, while silica reaching a
closed basin concentrates until it saturates. CO2 is the carbon cycle, and
section 4 says why that loop is left open.

**It also derives what the surface has BECOME, which is a different question from
what the profile is.** `build_surface_classes.py` emits two independent fields on
the mesh: `surface_cover`, what wind and light see, and `duricrust`, what is
cementing at or below the surface. Two axes rather than one chain, because a
duricrust forms inside a profile and loess sits on top of one, so they do not
compete for the same physical position and resolving them together would answer a
question that has none by branch order -- which is how `carved-zoned-v2` lost
basin fill from 203 basins. Precedence exists only within an axis and every
region records which rules it also satisfied.

Two of its inputs are outside pedology and both are deliberate. The dust
deposition field decides loess and pavement, and the aeolian roughness bracket is
worth a factor of 40 there, so the answer is a bracket rather than a value. And
the per-basin chemical divide decides which duricrust can form at all, because
gypcrete needs Ca surviving to saturation and the alkaline path removes it early.

Every Earth calibration lives in `pedology/config/pedogenesis.yaml` with its
source, and every surface-class threshold in `pedology/config/surface_classes.yaml`
with its source or an explicit statement that it is declared. Nothing in the
scripts hardcodes a Vesper number, so porting the component to another world is a
config change.

### 3.7 The biosphere

`biosphere/` runs LPJ-GUESS 4.1.1, patched for this world's calendar and
astronomy. The patch carries no planetary numbers itself: it points the model at
a generated `vesper.h`, which `build_vesper_header.py` derives from
`config/planet.yaml`, because the year length is a function of the stellar flux.

The model steps in 24-hour days with a year rounded from the orbital period. That
is why the semimajor axis has to be locked before the biosphere is built against
it, and why later flux changes go through luminosity instead: the period is a
function of the flux, so moving the flux the other way silently invalidates a
compiled model. That keeps every per-day
rate constant calibrated against the absolute time it was calibrated against and
confines the error to daylength, where its sign is known. Stepping in real 30-hour
Vesper days would put 25.9% into respiration, decomposition and phenology alike.

Annual degree-day limits are rescaled by the same orbit-derived factor, because a
Vesper year reaches about half the annual GDD of a 365-day one and Earth thresholds
would otherwise exclude every tree for reasons unrelated to the climate.

Forcing arrives as one binary driver file carrying however many years of climate
it was built with, and the model cycles through them. One year is a fixed
climate. Several are how a variable star reaches the biosphere; a single
repeating year cannot represent one at all.

## 3.8 Economic minerals

`minerals/` emits ore prospectivity per deposit type, in two artifacts split by
genesis and therefore by lifetime. The tectonic and magmatic half is a pure
function of the export and survives a re-run of the baseline; the weathering,
drainage and brine half reads a climatology, a lake solution and the per-basin
chemical divide, so it carries a climate in its identity as well as a terrain.
Design in `notes/economic-minerals.md`. Two rules govern where it attaches, and
both are expensive to retrofit.

**Minerals are an overlay, never a lithology.** `substrate_class` sets
erodibility and therefore terrain, bare-rock albedo and therefore climate, soil
texture and therefore the biosphere, and solute chemistry and therefore the
weathering fluxes. A mesh cell is 15 km across and a porphyry is one or two, so
admitting ore as a rock class would move the planet's energy balance on the
strength of a mine. The layer is prospectivity per cell, read only by what comes
after climate, and consequently rebuildable without invalidating a climate run.

It emits a continuous 0-1 prospectivity field and does NOT place individual
deposits. A porphyry is one to two kilometres against a 15.19 km cell and a vein
is under a hundredth of one, so a discrete deposit is invisible at every
resolution this pipeline runs at. Deposits belong with the downscaling pass,
which is where glacial overdeepening goes for the same reason -- if a thing is
smaller than a cell, it is not this pipeline's to place.

The split is by genesis: whatever concentrates a deposit has to be MODELLED
wherever the deposit is placed. Tectonic and magmatic types depend on arcs, fold
belts, LIPs and cratons, which Orogen models; weathering, drainage and brine
types depend on climate, runoff and the chemical divide, which it does not.
Supergene enrichment needs both, and is the case that shows why this is a split
rather than a handover.

Note that the tectonic half is nonetheless COMPUTED downstream, from the export
rather than inside the generator. Every input but one is already exported, and
keeping the arithmetic outside means a value that does not exist during
generation cannot feed back into erodibility or albedo -- the overlay rule above
becomes structural instead of remembered -- while a rule change costs no rebuild.

**Exhumation is what earns the layer**, rather than it being a recolouring of the
rock map. A porphyry forms one to five kilometres down and is destroyed by deep
erosion; an orogenic gold system forms five to fifteen down and is revealed by
it. The same erosion that removes one exposes the other, and Orogen tracks it
through `cover_thickness` and `erosionDelta`. So porphyry prospectivity rejects
the 59.5% of the arc belt stripped to its granodiorite root, where the porphyry
level went with the section above it.

## 4. Why this is not a straight line

Three quantities each depend on the other two.

**Drainage depends on climate.** Which basins survive is a water balance.

**Climate depends on drainage.** Closed-basin fill is a large minority of this
planet's land and the brightest thing on it -- salt crust and playa clastics
against a much darker land mean. Carve the basins and the world gets darker. This
is also the channel through which a lithology bug reached the climate, twice, so
measure it on the surface the model actually sees rather than on the rock table.

**Climate depends on the biosphere, and the biosphere on climate.** Bare rock and
a vegetated surface differ in land albedo by enough to be worth several kelvin,
and the two reach any given design mean at *non-overlapping* stellar fluxes. That is the load-bearing fact: no single flux is robust to the vegetation
question, so the orbit and the biosphere are one choice, not two.

The two interact rather than adding. Vegetation paints everything that can carry
a canopy at a single value, so it masks bare-rock variation but not the barren
classes; a lithology change confined to closed-basin fill therefore moves the
*vegetated* albedo roughly twice as far as it moves the bare one. Quote the
vegetated figure when the question is what the climate will do.

**That coupling is radiative and aerodynamic, and it is not hydrological.** Three
fields carry the vegetation state into ExoPlaSim: albedo, roughness and forest
fraction. Albedo is the priced channel above. Roughness is a real second one,
since it sets the turbulent exchange coefficient and so the evaporation rate.
Forest fraction is narrower than it sounds -- with SIMBA off it only mixes the
snow albedo, so it too is radiative. **What is absent is transpiration.**
ExoPlaSim's land surface is a single bucket with no stomatal control, no LAI
dependence and no rooting depth, so a vegetated cell and a bare cell holding the
same soil water evaporate identically apart from roughness. Soil water capacity
does come back from `pedology/`, which carries the biosphere's soil carbon, but
that is a channel through the soil rather than through the plant, and it has been
measured near-inert.

That absence is a limitation rather than a footnote, because of what the channel
would have been worth. `pedology/README.md` records, from Lapides et al. (2024),
that giving LPJ-GUESS a bedrock vadose zone raises annual transpiration by 100 to
150 mm, which is the same order as this world's entire land runoff. Runoff is the
denominator of the carve criterion and is a small residual of two much larger
fluxes, so an error in land evaporation arrives in it several times magnified.
The sign is not obvious either way: stomatal closure under stress cuts
evaporation below what a wet bucket gives, while deep roots reaching water the
bucket cannot hold raise it above. So the biosphere-climate coupling is priced in
kelvin and unpriced in millimetres, and the carve is decided in millimetres. The
evidence, and the reading of the model source that establishes the absence, are
in `notes/audits/missed-couplings.md`, finding 4.

**Soil depends on the biosphere, and the biosphere on soil.** Texture, pH and
regolith depth are weathering products of lithology under a climate, but the
organic fraction is what the vegetation leaves behind, and it changes the bulk
density and water-holding capacity the vegetation then grows in. `pedology/`
therefore iterates against `biosphere/` rather than running once before it.

The loop is therefore: assume, compute, feed back, repeat.

**There is a fourth loop, and it is currently cut rather than closed: dust.**
Emission reads the wind and the surface, the burden sets an optical depth, the
optical depth changes the radiation, and the radiation changes the wind. The
offline chain in `aeolian/` runs that path exactly once per iteration and stops,
which is what makes it PRESCRIBED rather than interactive: the dust the model
sees is the dust the previous climate produced. That is a defensible
approximation for one iteration deep and it is not one for a converged answer,
because the reopening test in `notes/dust.md` fired -- a prescribed field cannot
respond to the winds the dust itself changes. DUST-3 closes it by putting
emission in the model. Until then the cut is the reason `build_dust.py` sits
BELOW the climatology in the register and surface code 1811 sits ABOVE the next
run, which reads as a contradiction in the diagram and is really one loop drawn
across two iterations.

**One loop is deliberately left open: the carbon cycle.** `config/planet.yaml`
fixes CO2 at 450 ppm, and nothing in this project solves the carbonate-silicate
balance that would set it. That is a defensible choice for a snapshot climate --
CO2 was chosen alongside the flux to land the target temperature -- but it should
be read as an assumption rather than a result, because the pipeline now measures
what the assumption costs. Silicate weathering converts CO2 into buriable
alkalinity, so a weathering rate is also a statement about the outgassing the
world needs to hold its atmosphere steady.

Closing this loop is not a matter of adding a step. It needs a weathering law in
pCO2 and a climate response to it, which is a coupled calculation across pedology
and exoplasim rather than an analysis in either, and the weathering
concentrations currently in use carry no CO2 dependence at all. What exists is a
bound on the size of the assumption, not its resolution.

One caution about reading any of it. The budget is a sum over lithologies, and
a class that is a small share of the land can still decide the total if its
solute chemistry is extreme -- which is how an unexamined mapping for a rock
that was 0.00% of land came to supply nearly half the planet's CO2 drawdown once
a tectonic bug was fixed and the class appeared. `weathering_fluxes.py` now
reports per-class contributions for that reason. Check which class is on top
before quoting the total.

Two things about it are worth holding on to. The requirement is driven by **land
area, not by weathering intensity**: this world weathers less per unit area than
Earth, being drier, and needs more outgassing anyway because it has roughly twice
the land. A big-land planet is a high-outgassing planet or it is a cold one. And
only the exorheic share joins the marine carbonate feedback that stabilises CO2 --
endorheic alkalinity still buries carbon, on its own basin floor, but it is
decoupled from the loop that regulates. `pedology/scripts/thermostat_efficiency.py`
measures that share.

**The verdict map is antitone, not monotone, and that changes what the loop
does.** Carving removes closed-basin fill, the brightest lithology, so the land
darkens, the world warms, open-water evaporation rises, and basins that were
marginal would now stay closed. A larger carve set produces a *smaller* next
verdict. An antitone map does not approach a fixed point from one side; it
oscillates, and the successive verdicts bracket the answer rather than converging
onto it.

That is a better procedure than a one-sided approach, because a bracket is
measurable. Take the verdict at both bounding climates -- the cold, bright,
bare-rock end and the warm, dark, vegetated end -- and carve only the
intersection. Everything between the two is the marginal set *by construction*
rather than by a tolerance chosen after the fact, and the width of the bracket is
the honest uncertainty on the carve.

Over-carving is a budget item, not a lost landscape. A build is regenerated from
the planet code plus a verdict rather than edited, so an over-carve costs a
terrain, hydrography and boundary-condition rebuild and nothing else. That is how
`carved-zoned` was abandoned once its verdict turned out to have been computed on
antipodal climate: its carves could not be un-cut *within that build*, and the
build was replaced wholesale.

On iteration 2, the already-carved set should be re-evaluated against the new
climate and the number that would no longer have carved reported. That number is
the overshoot, and it is the honest measure of how much the first pass cost.

## 5. Where the pipeline currently stands

**`world_state.json` in the project root, generated by `scripts/world_state.py`.**
It is regenerated from the artifacts after anything that changes a build, a run
or a verdict, and it is the only place current values are written down.

They are not repeated here on purpose. This document describes how the pipeline
works and why it is shaped the way it is; those change slowly. Which build is
active, how many basins survive, what the world's mean temperature is and how
much of the land is closed-basin fill change every iteration, and prose restating
them is wrong within the day. A number appears in this document only when it is a
decision, a threshold, an identity, or when the magnitude carries an argument
that fails without it.

What does not change between iterations, and is worth knowing before reading the
generated state:

- **The land fraction is fixed.** Carving changes where water leaves a basin, not
  where the coast is, so the land/sea split is identical across every build.
  Anything that moves it is a bug.
- **This world is far more endorheic than Earth**, by a large factor rather than
  a little, in every iteration so far. That is the fact the hydrography component
  exists to handle, and it is why lakes and evaporite matter to the climate here
  when they would be a detail on Earth.
- **Carving reduces closed-basin fill, and fill is the brightest thing on the
  land.** So each carve iteration darkens the world and warms it, which is the
  feedback that makes this a loop rather than a sequence.
- **Vegetation masks bare-rock variation but not the barren classes.** A
  lithology change confined to closed-basin fill therefore moves the *vegetated*
  land albedo roughly twice as far as it moves the bare one. Bare-rock figures
  systematically understate what the climate will see.
- **The basin catalogue is stable across builds.** Basin ids are computed on the
  pre-conditioning surface, so a verdict computed against one build still refers
  to the same basins in the next. `manifest.hashes.basinCatalogue` is the check.

## 5b. Where the planet sits, and why

The flux was set by a derivation rather than inherited, on 2026-08-16. Recorded
here because the previous target had propagated through six files with no
statement of where it came from, and that is exactly the failure this document
exists to prevent.

**The old target was 290-293 K and had no derivation.** It traced to a constant
in `compare_albedo_bracket.py` whose docstring justified it as "the design range
this project has been aiming at since the first sweep, unchanged". That is a
number justified by its own persistence. It had reached `CLAUDE.md`, this file,
`config/planet.yaml` three times, `world_state.json`, and the convergence
criterion in `assess_convergence.py`, which chose its 0.15 K tolerance *because*
the band was 3 K wide.

**The mean is now chosen for habitability by latitude band.** Summer and winter
temperature per band were measured on three converged runs and projected across
candidate means. The trade is tropics against poles, and the tropics win on area:
cooling moves half the land out of sustained heat stress at the cost of a tenth
of it going from harsh to extreme. Flux 0.945 puts the tropics near +33 C in
their warmest month rather than +36, keeps most land in Earth-like conditions,
and leaves the polar margins severe but small.

**Glaciers are decoupled from the mean, and that is the load-bearing finding.**
Cooling is close to useless for making them. Summer amplification is nearly flat
with latitude, 0.67 to 1.11 K per K of global mean, while winter amplification
runs 0.92 to 3.30. So cooling buys brutal winters and barely touches the summers
that control ablation. Freezing a polar summer needs about 19 K of global
cooling, which would put polar winter near -94 C.

Glaciers come from relief instead. Peaks reach 4.59 km and a summit needs only
about 1.37 km above its grid cell to sit at freezing in summer, so every band
poleward of 30 degrees has peaks below freezing year round, with snowfall around
275 mm/yr water equivalent at 50-60 degrees. The equilibrium line moves roughly
700 m of elevation over the stellar cycle, so they advance and retreat visibly on
a generational rhythm.

**The poles are the worst place for glaciers on this world**, which inverts the
terrestrial intuition and does so from first principles. At 32 degrees obliquity
the pole receives 0.530 of the stellar constant as daily mean insolation at
summer solstice against the equator's 0.270 -- nearly double -- and both polar
caps are land, so there is no ocean buffer and no inherited ice to reflect it
away. Polar summers reach +39 C, hotter than the tropics, while polar winters
reach -31 C: a seasonal range of 70 K, beyond anything terrestrial. Earth's poles
are cold in summer largely *because* ice is already there. Vesper never
established that feedback.

The cold-summer band is therefore 50-60 degrees, which is where the glaciers are.
That band being at 50-60 rather than 45 or 65 is this world's particular
continents rather than a general rule; the inversion itself is mechanism.

**The cycle pushes the mean up, very slightly.** This was predicted the other
way round and measured wrong. The expectation was that the ice-albedo feedback
would amplify cooling, making T(f) concave and pulling the cycle-mean below the
static mean. Three quasi-equilibrated T42 points say the opposite: the slope is
196 K per unit flux across 0.9125-0.945 and 209 across 0.945-0.968, so **T(f) is
convex** over the range the cycle actually spans. Those three sit on a superseded
terrain and the curvature is what transfers, not the intercept; the slope itself
is canonical in `lib/sensitivity.py`, measured on the active build, and nothing
should quote 196 or 209 as a conversion factor.

The reason is that there is not enough ice for the ice-albedo term to lead. Sea
ice runs 6.2% to 0.2% across the interval, while the Planck response alone goes
as 1/(4 sigma T^3) and therefore grows with temperature, and the water-vapour
feedback grows with it too. Both push convex and both outweigh ice here.

The magnitude is negligible either way -- about +0.02 to +0.05 K before damping --
so nothing downstream changes. The mechanism is corrected because it was stated
as a reason, and a reason that is backwards will be reused.

One thing here is NOT settled and must not be carried as though it were. The
glacier result is unmodelled: the GCM reports `glac = 0` everywhere because it
cannot see a 4.6 km peak inside a 300 km cell, so the question is entirely
sub-grid.

`notes/glacier-rough-pass.md` now does that sub-grid integration properly rather
than with band-mean elevations, and **confirms the mechanism**: the area-weighted
mean glacier latitude is 47 to 57 degrees in every case and falls as the world
cools, because relief rather than latitude sets where ice survives. It remains a
rough pass -- an assumed lapse rate, a temperature criterion with no mass
balance, and a proxy offset for the 0.945 climatology -- so its areas are an
upper bound and are not properties of the world. The detailed treatment belongs
with the downscaling and sub-grid sampling machinery, not in the climate loop.

## 6. What happens next

Three nested loops and then a resolution change. The order matters in places
where it is not obvious, so those places are called out rather than left to be
rediscovered.

**A. The terrain loop.** Carve list to World Orogen, regenerate the terrain,
rebuild hydrography and boundary conditions on it, re-run the T42 baseline, take
the verdict again. Hydrography has to be rebuilt *after* the carve and before the
climate, because the carve changes the drainage the climate is integrated over.

Concretely, from a build in `source/` with nothing derived from it yet. Each
step reads the one before, and `check_consistency.py` reports which are still
missing at any point:

```bash
python exoplasim/scripts/rebuild_binaries.py --verify # CLAUDE.md rule 4. FIRST, always
python scripts/check_consistency.py                   # rule 8
python scripts/smoke_test.py

python hydrography/scripts/build_hydrography.py       # drainage, basins, coupling
python exoplasim/scripts/build_boundary_conditions.py # land mask, topography
python exoplasim/scripts/build_surface_albedo.py      # lithology albedo, no lakes yet
python exoplasim/scripts/build_surface_roughness.py   # z0
                                                      # soil_water_source: uniform
python exoplasim/scripts/run_exoplasim.py             # BOOTSTRAP run; hours
python exoplasim/scripts/assess_convergence.py <run>  # NOT optional: it writes the
                                                      # run's status, and nothing may
                                                      # treat a run as settled without it
python exoplasim/scripts/build_climatology.py <run>   # then set baseline_climatology
python exoplasim/scripts/index_runs.py                # INDEX.json is the only record

python hydrography/scripts/surface_water.py           # lakes, now a climate exists
python exoplasim/scripts/build_surface_albedo.py --lakes <surface_water.nc>
python pedology/scripts/build_soil.py                 # soil, with no biosphere yet
                                                      # soil_water_source: pedology
python exoplasim/scripts/build_surface_soil_water.py  # 229, from that soil
python aeolian/scripts/build_dust.py                  # needs the climatology AND the lakes
python exoplasim/scripts/dust_optics.py               # only if the spectrum moved
python exoplasim/scripts/dust_aerofile.py             # the model's aerofile
python exoplasim/scripts/build_surface_dust.py        # 1811, only if ndustrad = 1

python exoplasim/scripts/run_exoplasim.py             # the BASELINE; hours
python exoplasim/scripts/assess_convergence.py <run>
python exoplasim/scripts/build_climatology.py <run>   # repoint baseline_climatology
python exoplasim/scripts/analyze_climatology.py <...> # Koppen, biomes, the report
python exoplasim/scripts/run_stellar_cycle.py         # A2 below; hours. NOT optional

python exoplasim/scripts/dust_forcing.py              # the surface forcing the verdict reads
python hydrography/scripts/carve_verdict.py           # the verdict, weighted per A2
python hydrography/scripts/export_carve_list.py       # the list Orogen consumes

python scripts/error_budget.py                        # re-rank against the new state
python scripts/world_state.py                         # LAST: it reads everything above
```

`python scripts/pipeline.py --plan carve_list` prints that list against the
current state, with what is already present marked as skippable, so it does not
have to be followed from the top every time. `--status` answers the other half:
**which open tasks touch a step upstream of the carve.** That is the carve gate
-- not a ceremony about irreversibility, but the ordinary condition that nothing
outstanding still moves an artifact the verdict is computed from. A task names
its step as `[step: <id>]`; one that names none is reported as a residual and
arbitrated rather than counted or ignored, because the graph narrows the
judgement and does not replace it.

The three at the top and the two at the bottom are the ones most often skipped
and the ones that cost most when they are. `rebuild_binaries.py --verify` is
first because a `.venv` reinstall silently discards every patch and nothing else
in the list would notice. `world_state.py` is last because it reads every
artifact above it, and running it earlier records a state that no longer holds.

Everything below the climatology in the register -- the derived surface classes,
brine paths, the phosphorus budget, weathering fluxes, both prospectivity fields
-- is regenerated after the baseline and before anything quotes it. They are not
in the list above because they gate nothing in the loop; they are consumers, and
`CLAUDE.md` rule 7 governs them: when the climatology moves they are worthless,
and regenerating them is a step rather than a task.

**The cycle run is in that list because leaving it to prose lost it once.** A2
said it belongs before the verdict, nothing tracked it, and a full verdict was
taken on the mean climate before anyone noticed. It is now CYC-1 as well.

**The first climate run on a new terrain is a bootstrap, and its numbers are not
the baseline.** Three of the surface fields a run consumes cannot be built
without a climatology: the lakes, the lake compositing inside the albedo (174 to
176), and soil water capacity (229), which comes from a soil that has to be
weathered under a climate. Everything else -- topography, land mask, roughness,
forest fraction, the lithology half of the albedo -- is a pure function of the
terrain. So the loop is entered by running the model on the fields that do not
need it, and the run exists to produce the climatology the rest need. This is
written here because discovering it one field at a time costs a run each time.

**`model.soil_water_source` has to say `uniform` for the bootstrap and `pedology`
for the baseline, and the flip goes between them.** With it set to `pedology`,
`run_exoplasim.py` requires surface code 229 and refuses to start without it,
which is the one field the bootstrap exists to make possible. Flipping the key
changes a value `continue_exoplasim.py` compares against the run's manifest, and
it will not resume across that, so the bootstrap has to be finished -- converged,
seasonal window run, climatology built -- before the key moves. Every other
surface field the config asks for is a pure function of the terrain and is built
before the bootstrap starts.

**Lakes and soil water have to be in place before the run whose climatology the
verdict uses, not merely before the verdict.** Both change evaporation, which is
the numerator of the carve criterion, so a verdict taken on a climate that lacked
them is a verdict on the wrong evaporation. Their weights are not equal and the
ordering is worth its cost for one of them: lakes are worth -0.0139 on land-mean
albedo, because the cells carrying water are the bright playa and salt crust,
while the soil bucket has been measured at 1.06x on runoff for a 3.75-fold change
in capacity. Build both, and expect the lakes to be what moves the answer.

Set `baseline_climatology` in `config/planet.yaml` once a climatology exists --
until then every consumer raises rather than guessing -- and repoint it at the
baseline once that run finishes, so nothing downstream reads the bootstrap.

Seasonal snapshots are written by default, and `analyze_climatology.py` needs the
orbital phase they carry. A segment run with `--no-seasonal-output` has to be
extended before it can produce a climatology, which has happened once.

**The flux is re-derived on a new terrain, not carried.** Terrain moves lithology
exposure, relief and therefore land albedo, so the flux that lands the design mean
moves with it. Take two converged points spanning the target and interpolate;
never extrapolate from one, which this project has paid for twice. A point counts
only if it clears the remaining-offset test in `assess_convergence.py` rather than
the drift criteria alone: a T21 point at 0.90 once passed every drift test with
+0.348 K of approach still to run, and a point entering a bracket a third of a
kelvin low bends the slope. Which knob moves the flux is in loop C below, and it
depends on whether the calendar is locked.

**Change one thing at a time between runs, and prefer that to hitting a
temperature.** Re-deriving the flux is right when the flux is the answer you
want; it is wrong when it destroys a comparison you need more. The bootstrap and
the baseline differ in their surface fields -- lakes, the lake compositing in the
albedo, soil water -- and if they also differ in flux then nothing measures what
those fields are worth, and the bracket's slope stays measured in a regime the
model has left. Hold the flux, measure the surface, and move the flux afterwards
on a slope measured with the surface in place.

The temperature target is a weaker constraint than it looks, which is what makes
this trade easy. The design mean came out of the habitability-by-latitude
derivation in section 5b rather than being specified, so it is a property of the
terrain it was derived on. And it cannot be held to better than the terms already
outstanding: run-to-run spread on a converged pair is 0.23 K, and dust is priced
at 0.5 to 2 K of cooling that nothing has computed yet.

That is three converged runs before the first verdict: two for the flux bracket
and one for the baseline. The bootstrap is one of the bracket runs rather than a
fourth, since any converged climatology will do for fields that are themselves
about to be rebuilt.

**Know which base you are on, because a pre-carve build restarts the count.** The
first verdict taken on one is iteration 1, whatever a previous line of builds had
already carved, and the overshoot measurement in section 4 applies only where the
terrain being re-verdicted is itself the carved one.

**A3. A forcing change costs a PREDICTION and an A/B, not an iteration.**
Rewritten 2026-08-18. What stood here told you to land one forcing change per
iteration so that each could be attributed. That was wrong, and it was wrong in a
way worth recording, because the reasoning looked sound.

**Attribution cannot gate anything, so it must not be priced as though it could.**
A correct term goes in because it exists; section 7's physics-is-not-a-knob rule
already settles inclusion, and no attribution result can reverse it. The old rule
therefore spent runs buying information that could never change a decision.
`CLAUDE.md` rule 7 removes the usual fallback too: a wrong verdict is recoverable,
because a build is replaced wholesale rather than edited, so "we could not
attribute it and the carve is permanent" is not an argument either.

**Attribution still pays in exactly one place: diagnosing a surprise.** The
physics rule depends on it -- if a correct term makes an agreement worse, that is
information about the implementation, the comparison, or a second error
cancelling the first -- and telling those apart means knowing which term moved.
But that cost is CONTINGENT. Serializing iterations pays it in full every time,
in converged runs of hours, against a risk that materialises occasionally.

**So buy it the cheap way. Testing a change is not the same as needing an
iteration.** A converged run measures a RESPONSE in kelvin and is expensive. A
short A/B measures a FORCING in W/m2 and is not: branch two arms off one common
restart, run a few orbits, and difference the flux diagnostics. Convert to kelvin
afterwards with the canonical slope in `lib/sensitivity.py` rather than by waiting
for the model to equilibrate. Every term can be attributed this way while the
whole bundle still costs ONE converged run.

Four things make that A/B trustworthy, and none of them is optional:

- **Every forcing change lands with a quantitative prediction of its own effect**,
  stated before it is run, with what result would mean "wrong". That is section
  7's rule that a check needs a right answer, applied to physics rather than to
  code. PHYS-6 is the worked example: 2.5 W/m2 here against 1.75 solar-weighted,
  checked against Earth's measured 1.5-2.5.
- **Both arms use the SAME BINARY and differ only by a namelist key.** This is
  what the no-op-until-enabled convention is for: `h2osww` defaults to 1.0 and
  `ndustrad` to 0 precisely so that one executable can run both arms. Two
  binaries would confound the term with the rebuild, and the low-I/O patch
  changes the restart layout, so arms built either side of it cannot share a
  restart at all.
- **Both arms branch from ONE restart.** Run-to-run spread on a converged pair is
  0.23 K, which is larger than several of the terms being tested; a shared
  initial condition removes most of it and turns the comparison into a paired one.
- **The segments are labelled as diagnostics**, so a short A/B tail never enters a
  convergence window or a climatology. That is CLIM-9, and until it lands this
  practice is unsafe rather than merely untidy.

**Bundle the terms, then check the SUM against the sum of predictions.** If they
agree, no attribution was needed and none was bought. If they disagree, bisect,
and only into the subset whose predictions were soft. What bundling genuinely
risks is `notes/failure-modes.md` class 15, two errors that nearly cancel; the
mitigation is the per-term prediction above, not serialization, because two
cancelling errors hide just as well in a serial sequence that nobody predicted
the size of.

**What survives from the old rule is a budgeting fact, not an ordering one.**
`build_surface_albedo.py` has a `modelled` mode that takes tree cover from an
LPJ-GUESS `fpc.out` instead of asserting a uniform vegetated endmember, and it is
the mode this project should end up in. The two land-albedo endmembers are 15 to
19 W/m2 apart in absorbed flux, against 21 W/m2 for the entire 0.85-to-0.95
stellar sweep that produced a 33 K range. So adopting it is a forcing change
comparable to the largest this project varies on purpose, and **it requires the
flux to be re-derived**. That is a cost to budget for, and it is the same cost
whether or not a carve shares the iteration.

The corollary is the answer to "when does the biosphere run". Not before the
carve, because its driver is built from a climatology and a soil the carve
replaces. After the re-baseline, on final terrain, where its output can be
adopted deliberately in the iteration after that.

**A4. Dust is not carve-neutral, and the carve waits on one climate run.** Its
global-mean forcing is negligible in kelvin; the local effect is not. The
question splits in two and only one half is settled.

The LAKE half is done. Dust dims the lake, Penman evaporation falls, and more
basins overflow: a dust-free verdict UNDER-carves by about 23 basins, which is
small and in the recoverable direction. That was settled offline by putting the
dust SURFACE forcing through `carve_verdict.py --dust-forcing`, and the surface
balance is the operative one -- the top-of-atmosphere sign over bright basin fill
is energy retained in the atmosphere, not delivered to the ground.

The CATCHMENT half is larger, runs the other way, and cannot be reached offline,
because it is a precipitation response and not a surface energy balance.
Suppressed precipitation cuts the runoff the criterion divides by, which
OVER-carves, and over-carving is the irreversible direction. It needs one
prescribed-dust climate run, specified with its gates and a prediction in
`aeolian/notes/prescribed-dust-run.md`. The ORDER here is a dependency and not a
serialization: the verdict divides by a runoff that dust moves, so it has to be
taken on a climatology that already has dust in it. Run it, take the verdict on
its climatology, then carve. A3 no longer asks for a separate iteration and never
should have; what it asks for is the prediction that run already carries.

**Do not record a carve list as dust-independent.** `notes/dust.md` has the
numbers and `aeolian/analysis/dust_runoff_sensitivity.json` has the conversion
from a precipitation change into basins.

**B. The soil and biosphere loop, at T42.** For a given climate:

1. `pedology/scripts/build_soil.py`, with no biosphere on the first pass.
2. `biosphere/scripts/build_lpj_driver.py`, then `run_lpj_guess.py`.
3. `build_soil.py --soil-carbon <run>/cpool.out`, then LPJ-GUESS again.
4. Repeat 3 until the criteria in `pedogenesis.yaml` are met.

Check whether step 3 moves anything before assuming it needs iterating:
LPJ-GUESS computes its own soil carbon internally, so the pedology organic
feedback may be second-order.

**A2. The stellar cycle, and where it belongs in the order.** The cycle run is
NOT a final flourish. It has to come after a baseline that is FINAL, not merely
converged, because a cycle is variance about a mean and centring it on the wrong
one describes the wrong world -- but it belongs BEFORE the carve verdict, for a
reason that is easy to miss.

  "Final" is the operative word and it is stricter than "converged". Anything
  still outstanding that moves the mean -- a radiative correction, a decision to
  carry dust -- has to land first, or the cycle describes a world that will not
  exist. A converged run at the wrong mean is exactly as useless as an
  unconverged one. Run the cycle last of the things that change the climate.

Carving is irreversible. Outlet incision does not undo when the warm phase
returns, so every cold, wet excursion that pushes a basin to overflow carves it
permanently. **The terrain therefore ratchets toward the state implied by the
cycle's wet extreme, not its mean.** Taking the verdict on the mean climate
systematically under-carves, because overflow is a threshold process and the wet
phase contributes disproportionately. Bedrock incision is a 10^3 to 10^5 year
process against cycles of 11 and 57 Earth years, so it is the integral over many
cycles that matters rather than any single one -- which means the verdict wants a
climate somewhere between the mean and the wet extreme, weighted by time spent
overflowing, not the mean alone.

**The two components ratchet differently, and that is the reason there are two.**
The medium one is too fast for ice to follow and too fast for an outlet to
incise, so it contributes only through the tail of its distribution. The long one
is slow enough that glaciers equilibrate and long enough that a wet excursion is
sustained, so it is the one that actually cuts. And because 57/11 is
non-commensurate the deepest minima differ in depth rather than repeating, so
successive advances reach different distances and the landscape ends up recording
which past minima were severe. A cycle run therefore has to be long enough to
sample that: its length is set in periods of the LONG component, not the medium
one.

The cycle run also measures the damping factor, which is currently known only as
a range of 0.24 to 0.6 -- a factor of 2.5 on every temperature excursion derived
from a flux amplitude. That single number converts any future amplitude choice
into a climate without another run, so measuring it once is worth more than the
run that measures it.

Sequence, then: baseline at the chosen mean, cycle run centred on it, verdict on
a cycle-informed climate, and only then the terrain loop.

**C. The vegetation-climate loop.** `build_surface_albedo.py --mode modelled`
turns the run's foliar cover into surface albedo and forest fraction, then the
climate runs again on it. Two checks belong here and neither is optional.

*The flux.* If the modelled land albedo differs much from the assumed value the
world may leave the design mean derived in section 5b. **Move the flux by
changing the star's luminosity, not the orbit.** The year length depends on the
semimajor axis, and moving the orbit changes it; the year is compiled into
LPJ-GUESS, so that would force a rebuild and a driver regeneration. With the
semimajor axis locked, F = L/a^2, so luminosity is the free parameter and the
calendar does not move.

The knobs swap at the lock, and the lock is what makes this a rule rather than a
preference. Before the calendar is fixed, the semimajor axis is the better knob:
it keeps the star itself untouched, so `k25v` stays valid, where covering the
same flux range with luminosity moves the effective temperature far enough to
leave the window the spectrum was interpolated in. Once LPJ-GUESS is compiled
against a year, that knob costs a rebuild and a driver regeneration, and
luminosity is the only cheap one left. A bracket taken on the axis varies two
things rather than one, because the year moves with it, so the measured slope
carries a small seasonality effect through the calendar. That is second order for
an annual mean, and it is stated rather than implied.

Mind the scaling, and mind the spectrum. L scales **1:1** with F at fixed orbit,
not 2:1 as an earlier revision of this document said, and at fixed radius
L ~ T^4, so the effective temperature moves as F^(1/4):

| flux change | luminosity | effective temperature | band-1 fraction |
| --- | --- | --- | --- |
| 1.5% | 1.5% | 18 K | ~0.9% |
| 4% | 4% | 49 K | ~2.3% |

`k25v` is interpolated to 4965 K between the 4900 and 5000 K BT-Settl models, and
that 100 K grid step is worth 4.7% in the fraction of flux below 0.75 um. So
small adjustments are free, but **past about 2 to 3% in luminosity, re-run
`build_stellar_spectrum.py` at the new effective temperature**, and past roughly
+35 K the target leaves the 4900-5000 K bracket entirely and the pinned SVO grid
points have to change with it. The spectrum exists to get snow and ice albedo
right, so letting it drift silently would undo the reason it was built.

*The carve verdict.* It was taken on assumed vegetation. Real vegetation changes
the climate, which changes evaporation over catchments, which can change the
verdict. Re-run it; if basins flip, back to loop A.

**D. The resolution change.** Only after A, B and C have settled.

1. Climate at T85, with the T42 vegetation regridded as its boundary condition.
2. Soil and biosphere at T85 on that climatology, loop B again.
3. Climate at T85 again with T85 vegetation, unless step 2's vegetation turns out
   close to the regridded T42 field, which is a cheap comparison worth making
   first.

**The biosphere never decides the resolution.** LPJ-GUESS gridcells are
independent columns, so its cost is linear in land-cell count and trivial at
either resolution -- tens of minutes against a couple of hours on 16 ranks,
where the climate model's cost is what actually decides.
Running it at T85 on T42 forcing would resolve detail that is not in its input.
Expect T85 to lower total NPP, and treat that as a resolution bias rather than a
result: productivity saturates with water, so averaging the forcing before the
model sees it inflates the answer.

**E. The cycle reaching the biosphere, last, on the settled world.** This is not
a second decision about when to run the cycle: A2 says when, and why it has to
precede the verdict. What is left for the end is letting LPJ-GUESS see it, on a
world whose terrain and climate have stopped moving. The amplitudes and periods
are the two components in `config/planet.yaml`, and the run has to be long enough
to sample the long one, which A2 also covers. Build the climatology with
`--per-year`, pass the sequence to
`build_lpj_driver.py --climatology y0.nc y1.nc ...`, and the biosphere sees the
cycle rather than its average. That matters because productivity responds
annually and saturates, so a run on the cycle mean over-predicts it: 4.8% at the
one cell measured. It does *not* reach survival thresholds, which LPJ-GUESS gates
on a twenty-year mean of coldest-month means and which therefore smooth the cycle
away almost exactly. Then the regional products, biomes, Koppen and lake maps.

Three things about the calendar and the resolution that are easy to get wrong,
each worked through in `biosphere/README.md`:

- **Spin-up is counted in simulation years**, and a simulation year is about
  half an Earth one, so the shipped 500 buys roughly half the absolute time
  Earth practice assumes. Year counts are scaled UP and annual sums scaled DOWN
  by the reciprocal factor; the two directions are opposite, both lists are
  named in `build_vesper_pfts.py`, and the factor is derived there from the
  configured orbit rather than written anywhere.
- **T85 redistributes precipitation, it does not merely resolve it.** Steeper
  relief means stronger orographic ascent and deeper rain shadows, so lee basins
  dry and windward coasts wet. That reaches the carve verdict, which should be
  re-taken at T85 rather than carried over from T42.
- **The 30-hour day widens the real diurnal range and LPJ-GUESS cannot see it.**
  `dtr` reaches only the biogenic VOC scheme, and every cold limit runs through a
  twenty-year mean, so there is no daily-minimum mortality to trigger.

## 7. Conventions this project holds to

Provenance travels with every artifact: config, input hashes, software versions
and the terrain hash go into each run manifest and each analysis report.

Thresholds are fixed before results are seen. The albedo bracket's convergence
and marginality criteria were written into the script's docstring before any run
finished, and were then applied against a result that cleared one of them by
0.05 K.

Estimates that cannot be verified are bracketed rather than guessed, and the
bracket is reported. This is how albedo, evaporation and the carve verdict are
all handled.

Claims are checked against the artifact rather than the documentation. Several
findings in this project's history came from comparing two products that were
supposed to agree and finding they did not.
