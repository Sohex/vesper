# Vesper: the worldbuilding pipeline

How geography, climate, water and life are computed for this planet, in what
order, and why the order is not a straight line.

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

## 1. The components

```
config/planet.yaml   Canonical planet, star, orbit, atmosphere. Every component reads it.
source/              World Orogen exports. Canonical, read-only.
lib/                 Shared readers: orogen.py (the export), gridding.py (mesh to grid).
hydrography/         Drainage, catchments, basin capacity, lake balance, carve verdict.
exoplasim/           Boundary conditions, climate integrations, climatology.
pedology/            Weathers lithology into soil, and into solute fluxes: CO2, silica, phosphorus.
biosphere/           LPJ-GUESS: vegetation, leaf area, carbon, PFT composition.
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

Three loops close in that diagram, and section 4 says why each has to.

## 3. The pipeline, step by step

### 3.1 Geography

World Orogen generates the terrain from a planet code, which encodes seed and
every slider. The current build is `01eshm059lt0b9mpgro2y83t`: seed 16236323,
2,500,001 mesh regions, closed basins preserved, lithology-modulated erosion.

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

Land comes from `surface_class`, never from `land_mask`. The two disagree by 1.9%
of the planet, all of it dry closed-basin floor below sea level, and `land_mask`
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
bucket. See `notes/lake-representation.md`.

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
resolve. The model confirms the fix from its own log: energy fraction below
0.75 microns is 0.38438 under `k25v` against 0.11588 under `k2`.

Everything before the re-baseline carries that bias. It is not invalidated in
kind and the direction is known, but it should be stated wherever it is quoted.
See `exoplasim/notes/stellar-spectrum-audit.md`.

A run directory names everything that changes the answer: resolution, flux, CO2,
rotation, obliquity, eccentricity, the physics switches, the spectrum, and a
digest of every surface input. Two of those were added after they had already
nearly caused a collision. The spectrum marker was missing when the k25v
re-baseline would have landed in the completed k2 run's directory, and the flux
was rounded to hundredths, so 0.945 and 0.94 both resolved to `s094`. Only the
geography digest separated them, which was luck. Anything physical that is not in
the directory name is a collision waiting for the run that changes it.

Enabling `model.energy_diagnostics` adds PlaSim's 28-term energy decomposition on
codes 360-387, which is the instrument for the constant -0.455 W/m2 that does not
close between the top of the atmosphere and the surface. The postprocessor does
not ship those codes, so `run_exoplasim.py` registers them at run time rather
than patching the vendored tree, which a reinstall would silently undo. The
residual needs a settled run: a segment taken one orbit off a restart sits several
W/m2 out of balance, because `configure()` resets the surface fields a restart
does not carry.

Convergence is a fixed six-part test, not a judgement: temperature drift below
0.05 K per orbit, top-of-atmosphere and surface balance trends below 0.05 W/m2
per orbit, sea-ice drift below 0.001 of planetary area per orbit, and mean
absolute imbalances below 0.5 W/m2. Runs that miss are labelled, not rounded.

A separate five-orbit window with 32 snapshots per orbit, excluded from the
pass/fail decision, forms the climatology.

### 3.4 Hydrography

The export deliberately does not route water: `drain_to` is raw steepest descent,
and on this planet most of the land drains into unpreserved single-cell pits --
the share is per build and is in `world_state.json`. Routing
is a hydrology decision and the exporter leaves it downstream.

`build_hydrography.py` resolves it with a priority flood over the 2.5M-region
mesh, which fills the noise pits while keeping the 3,629 preserved basins as
terminals. It then rebuilds each basin's hypsometry on the *finished* terrain,
because the catalogue's curves are measured on the pre-conditioning surface and
overstate capacity by about 1.5x. Finally it writes a sparse basin-by-grid-cell
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
and validated by applying the same calculation to ocean cells, which *are* open
water: 3.736 mm/day against the model's own 3.672, a ratio of 1.017.

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

Every Earth calibration lives in `pedology/config/pedogenesis.yaml` with its
source. Nothing in the scripts hardcodes a Vesper number, so porting the
component to another world is a config change.

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
181-day year reaches half the annual GDD of a 365-day one and Earth thresholds
would otherwise exclude every tree for reasons unrelated to the climate.

Forcing arrives as one binary driver file carrying however many years of climate
it was built with, and the model cycles through them. One year is a fixed
climate. Several are how a variable star reaches the biosphere; a single
repeating year cannot represent one at all.

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
and the two reach the 290 to 293 K design target at *non-overlapping* stellar
fluxes. That is the load-bearing fact: no single flux is robust to the vegetation
question, so the orbit and the biosphere are one choice, not two.

The two interact rather than adding. Vegetation paints everything that can carry
a canopy at a single value, so it masks bare-rock variation but not the barren
classes; a lithology change confined to closed-basin fill therefore moves the
*vegetated* albedo roughly twice as far as it moves the bare one. Quote the
vegetated figure when the question is what the climate will do.

**Soil depends on the biosphere, and the biosphere on soil.** Texture, pH and
regolith depth are weathering products of lithology under a climate, but the
organic fraction is what the vegetation leaves behind, and it changes the bulk
density and water-holding capacity the vegetation then grows in. `pedology/`
therefore iterates against `biosphere/` rather than running once before it.

The loop is therefore: assume, compute, feed back, repeat.

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
static mean. Three converged points say the opposite: the slope is 196 K per unit
flux across 0.9125-0.945 and 209 across 0.945-0.968, so **T(f) is convex** over
the range the cycle actually spans.

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

**B. The soil and biosphere loop, at T42.** For a given climate:

1. `pedology/scripts/build_soil.py`, with no biosphere on the first pass.
2. `biosphere/scripts/build_lpj_driver.py`, then `run_lpj_guess.py`.
3. `build_soil.py --soil-carbon <run>/cpool.out`, then LPJ-GUESS again.
4. Repeat 3 until the criteria in `pedogenesis.yaml` are met.

Check whether step 3 moves anything before assuming it needs iterating:
LPJ-GUESS computes its own soil carbon internally, so the pedology organic
feedback may be second-order.

**A2. The stellar cycle, and where it belongs in the order.** The cycle run is
NOT a final flourish. It has to come after a converged baseline at the chosen
mean, because a cycle is variance about a mean and centring it on the wrong one
describes the wrong world -- but it belongs BEFORE the carve verdict, for a
reason that is easy to miss.

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
world may leave the 290-293 K design band. **Move the flux by changing the star's
luminosity, not the orbit.** The year length depends on the semimajor axis, and
moving the orbit changes it; the year is compiled into LPJ-GUESS, so that would
force a rebuild and a driver regeneration. With the semimajor axis locked,
F = L/a^2, so luminosity is the free parameter and the calendar does not move.

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
independent columns, so its cost is linear in cell count and trivial either way:
4,106 land cells at T42 and 16,489 at T85, 23 against 94 minutes on 16 ranks.
Running it at T85 on T42 forcing would resolve detail that is not in its input.
Expect T85 to lower total NPP, and treat that as a resolution bias rather than a
result: productivity saturates with water, so averaging the forcing before the
model sees it inflates the answer.

**E. The stellar cycle, last, on the settled world.** 0.91 to 1.01 S-Earth over 8
Earth years. Build the climatology with `--per-year`, pass the sequence to
`build_lpj_driver.py --climatology y0.nc y1.nc ...`, and the biosphere sees the
cycle rather than its average. That matters because productivity responds
annually and saturates, so a run on the cycle mean over-predicts it: 4.8% at the
one cell measured. It does *not* reach survival thresholds, which LPJ-GUESS gates
on a twenty-year mean of coldest-month means and which therefore smooth the cycle
away almost exactly. Then the regional products, biomes, Koppen and lake maps.

Three things about the calendar and the resolution that are easy to get wrong,
each worked through in `biosphere/README.md`:

- **Spin-up is counted in simulation years**, so the shipped 500 is 247 Earth
  years. Year counts are scaled up by 2.022 while annual sums are scaled down by
  0.4946; the two directions are opposite and both lists are named in
  `build_vesper_pfts.py`.
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
