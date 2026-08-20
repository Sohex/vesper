# Vesper: the worldbuilding pipeline

How geography, climate, water and life are computed for this planet, in what
order, and why the order is not a straight line.

**This document is CANONICAL for the pipeline's REASONING**: what each stage is
for, why the order is what it is, why the loops close, and the traps. The GRAPH
itself -- every step, what it writes, what must precede it -- is
`config/pipeline.yaml`, and the two do not overlap. This document references
step ids and does not restate them, because one graph written in two places is
the duplication both files exist to prevent.

The rule that follows, and it is a rule rather than an aspiration: **an
artifact that no step in `config/pipeline.yaml` generates does not exist.** If
a script writes something that is in no step, either the step is missing or the
product is one nobody should be reading. Both are defects, and the second is
worse, because an unregistered artifact has no derivable consumers and so no
way to know what it invalidates when it changes. `smoke_test.py` fails when a
generator is absent from the graph. `CLAUDE.md` rule 7 depends on all of this:
"what is now worthless" is only answerable from a graph that is complete.

Vesper is a super-Earth orbiting a K2.5V dwarf: larger than Earth, higher
gravity, a longer day, and more obliquity. Every one of those is declared in
`config/planet.yaml`, which is the only place they are written down. The flux,
the year, the active build and the mean surface temperature move every
iteration; they are in `world_state.json`, which is generated. See section 5.

---

## 0. What a re-commissioning costs

**The words are in `CLAUDE.md` under "Vocabulary" and are not repeated here.**
What belongs here is the cost, because the cost is what makes the distinctions
worth drawing. The rows below are a RE-COMMISSIONING: the same build,
everything below the climatology redone. A bare re-run of the baseline is
cheaper still, and is only honest when nothing in the "yes" rows is reachable
from whatever changed.

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

**Price an orbit in WALL CLOCK, not in model time.** What a run waits on is the
model plus postprocessing plus the ocean and ice stream writes, and quoting the
model half alone understates a run by about a third. Per-segment model-only
figures are `native_runtime_seconds` in each run's manifest; the wall figure is
the gap between consecutive `MOST_REST.NNNNN`. At T42 on 16 ranks an orbit
costs about 140 s end to end -- a ceiling, measured before the second pyburn
fix landed, which should bring it nearer 110 s once a run made after the change
measures it. So a settling run off a near-equilibrium restart plus ten clean
orbits is a couple of hours, and a commissioning usually wants two of them. An
ITERATION is far more, because the generation in front of it regenerates the
terrain and every field that is a pure function of it.

**That affordability assumes the pyburn fixes are resident.** Unpatched,
postprocessing is roughly three quarters of the orbit rather than an overhead
on it, and a commissioning goes from about four hours to about twelve.
`notes/audits/pyburn-postprocessing-cost.md` has the measurements.

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
turn one picture into four. Both are in the register.

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

**The register is `config/pipeline.yaml`, and it is not reproduced here** --
one graph in one place, referenced from the other.

    python scripts/pipeline.py --register     # artifact, step, and what reads it
    python scripts/pipeline.py --status       # what is present, and what is not
    python scripts/pipeline.py --plan <step>  # the ordered steps to reach a target
    python scripts/pipeline.py --purge <step> # what a change to that step makes worthless

Consumers are derived rather than declared: what reads an artifact is the set
of steps that `need` the step which writes it, plus the steps marked
`reads_export` for the export itself. Declaring both directions would be the
same duplication one level down.

`pipeline.py` plans and never RUNS A STEP. Several steps in the graph cost
hours, and a script that could start one by accident is worse than no script.
It also does not check whether artifacts AGREE -- that is
`check_consistency.py`, and asking either to do the other's job would give two
answers to one question.

`--purge` is the exception and it only ever deletes; it is a dry run until
`--execute`. Rule 7 makes "what is now worthless" the question after any
change, and that is a graph question -- doing it by hand produces a keep-list
rather than an answer. Two properties to know before using it. It never
crosses `orogen`: loop A is a cycle, so unrestricted reachability from any
climate step would come back round and offer to delete the terrain, which is
wrong rather than merely alarming -- a new climatology does not invalidate the
build it was computed on. And it does not touch climate RUNS, which are
UUID-named and have no path in the graph; it names them and points at
`scripts/archive_runs.py`, which extracts a run's identity and verifies the
archive before deleting anything.

## 3. The pipeline, step by step

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
of the disagreement -- and `land_mask` would flood it.

Seven surface fields are supplied. Topography (129), land mask (172),
roughness (173), broadband and two-band albedo (174, 175, 176) and forest
fraction (212), with soil water capacity (229) when `model.soil_water_source`
is set.

Roughness is measured from land cover and subgrid relief rather than
asserted uniform, with its land mean anchored to the model's tuned default so
only the distribution moves; every other surface field falls back to a
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

ExoPlaSim 3.4.2, T42, 16 MPI ranks, 45-minute timestep, 10 layers, 50 m slab
ocean, interactive sea ice, glaciers enabled, and a measured stellar spectrum
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
inside the window. A separate five-orbit window with 32 snapshots per orbit,
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
cycle, and section 4 says why that loop is left open.

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

`biosphere/` runs LPJ-GUESS 4.1.1, patched for this world's calendar and
astronomy; the patch carries no planetary numbers itself but reads a generated
`vesper.h` derived from `config/planet.yaml`. The calendar decisions -- 24-hour
model days, rescaled degree-day limits, and why -- are `biosphere/README.md`'s.
What the pipeline has to know is the LOCK: the year is compiled into the
model, and the year is a function of the semimajor axis, so the axis has to be
locked before the biosphere is built against it and later flux changes go
through luminosity instead (loop C).

Forcing arrives as one binary driver file carrying however many years of
climate it was built with, and the model cycles through them. One year is a
fixed climate; several are how a variable star reaches the biosphere, and a
single repeating year cannot represent one at all.

## 3.8 Economic minerals

`minerals/` emits ore prospectivity per deposit type, in two artifacts split
by genesis and therefore by lifetime: the tectonic and magmatic half is a pure
function of the export and survives a re-run of the baseline; the weathering,
drainage and brine half reads a climatology, a lake solution and the per-basin
chemical divide, so it carries a climate in its identity as well as a terrain.
Design in `notes/economic-minerals.md`; mechanics in `minerals/README.md`.
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

## 4. Why this is not a straight line

Three quantities each depend on the other two.

**Drainage depends on climate.** Which basins survive is a water balance.

**Climate depends on drainage.** Closed-basin fill is a large minority of this
planet's land and the brightest thing on it -- salt crust and playa clastics
against a much darker land mean. Carve the basins and the world gets darker.
This is also the channel through which a lithology bug reached the climate,
twice, so measure it on the surface the model actually sees rather than on the
rock table.

**Climate depends on the biosphere, and the biosphere on climate.** Bare rock
and a vegetated surface differ in land albedo by enough to be worth several
kelvin, and the two reach any given design mean at *non-overlapping* stellar
fluxes. That is the load-bearing fact: no single flux is robust to the
vegetation question, so the orbit and the biosphere are one choice, not two.

The two interact rather than adding. Vegetation paints everything that can
carry a canopy at a single value, so it masks bare-rock variation but not the
barren classes; a lithology change confined to closed-basin fill therefore
moves the *vegetated* albedo roughly twice as far as it moves the bare one.
Quote the vegetated figure when the question is what the climate will do.

**That coupling is radiative and aerodynamic, and it is not hydrological.**
Albedo, roughness and forest fraction carry the vegetation state into the
model; **what is absent is transpiration**, because the land surface is a
single bucket with no stomatal control, LAI dependence or rooting depth. The
absence is a limitation rather than a footnote: the channel would have been
worth the same order as this world's entire land runoff, runoff is the
denominator of the carve criterion, and the sign is not obvious either way.
So the biosphere-climate coupling is priced in kelvin and unpriced in
millimetres, and the carve is decided in millimetres. Evidence and the
model-source reading that establishes the absence:
`notes/audits/missed-couplings.md`, finding 4.

**Soil depends on the biosphere, and the biosphere on soil.** Texture, pH and
regolith depth are weathering products of lithology under a climate, but the
organic fraction is what the vegetation leaves behind, and it changes the bulk
density and water-holding capacity the vegetation then grows in. `pedology/`
therefore iterates against `biosphere/` rather than running once before it.

The loop is therefore: assume, compute, feed back, repeat.

**There is a fourth loop, and it is currently cut rather than closed: dust.**
Emission reads the wind and the surface, the burden sets an optical depth, the
optical depth changes the radiation, and the radiation changes the wind. The
offline chain in `aeolian/` runs that path exactly once per iteration and
stops, which is what makes it PRESCRIBED rather than interactive: the dust the
model sees is the dust the previous climate produced. That is defensible one
iteration deep and not for a converged answer -- the reopening test in
`notes/dust.md` fired. DUST-3 closes it by putting emission in the model.
Until then, `build_dust.py` sits BELOW the climatology in the register and
surface code 1811 sits ABOVE the next run, which reads as a contradiction in
the diagram and is really one loop drawn across two iterations.

**One loop is deliberately left open: the carbon cycle.** `config/planet.yaml`
fixes CO2, and nothing in this project solves the carbonate-silicate balance
that would set it -- a defensible choice for a snapshot climate, to be read as
an assumption rather than a result. The pipeline measures what the assumption
costs (a weathering rate is also a statement about the outgassing the world
needs to hold its atmosphere steady), and what exists is a bound on the size
of the assumption, not its resolution; closing the loop needs a weathering law
in pCO2 and a climate response to it, a coupled calculation across pedology
and exoplasim rather than an analysis in either.

Two things to hold on to, and one caution. The requirement is driven by **land
area, not weathering intensity** -- a big-land planet is a high-outgassing
planet or it is a cold one. Only the exorheic share joins the marine carbonate
feedback that stabilises CO2; `pedology/scripts/thermostat_efficiency.py`
measures the share. And the budget is a sum over lithologies in which a tiny
class with extreme solute chemistry can decide the total, so
`weathering_fluxes.py` reports per-class contributions: check which class is
on top before quoting the total.

**The verdict map is antitone, not monotone, and that changes what the loop
does.** Carving removes closed-basin fill, the brightest lithology, so the
land darkens, the world warms, open-water evaporation rises, and basins that
were marginal would now stay closed. A larger carve set produces a *smaller*
next verdict. An antitone map does not approach a fixed point from one side;
it oscillates, and successive verdicts bracket the answer rather than
converging onto it.

That is a better procedure than a one-sided approach, because a bracket is
measurable. Take the verdict at both bounding climates -- the cold, bright,
bare-rock end and the warm, dark, vegetated end -- and carve only the
intersection. Everything between the two is the BRACKETED set *by
construction* rather than by a tolerance chosen after the fact, and the width
of the bracket is the honest uncertainty on the carve.

**`bracketed`, not `marginal`.** `marginal` names a per-basin landform in the
carve list; `bracketed` names our uncertainty. The two were briefly one word
and it leaked into an exporter; the definitions and the third neighbour,
`disputed`, are in CLAUDE.md's vocabulary.

**Both arms run at ONE flux, and it is the design flux.** Decided 2026-08-19.
The alternative -- each endmember at whatever flux keeps IT in the design
range -- collapses the bracket: those two worlds sit at the same global mean
by construction, so the intersection lands close to either verdict alone and
the reported width understates the real uncertainty, which something
downstream will quote. It is also the reading that matches what is actually
unknown: the orbit is chosen, not uncertain, while the vegetation state is
genuinely unknown when the verdict is taken, because the carve happens before
LPJ-GUESS has ever run and the verdict inherits an assumed biosphere.
Bracket the thing you do not know.

**The flux bracket itself runs on the VEGETATED branch only.** That is the
functional world, the one this project has chosen and the one the design flux
is defined against. Running it on both branches would be deriving two orbits
for a world that has one.

**The bare-rock arm is a BOUND, not a world.** Nobody claims Vesper sits at
the design flux without a biosphere; that run exists to produce the cold-end
verdict and nothing else. Say so wherever its numbers appear, because a reader
meeting its climatology in this repository will otherwise take it for a
description.

**Expect the bracket to come out wide, and do not read width as failure.** The
endmember spread widens as the world cools and sea ice grows back to amplify
it -- the T21 work measured it nearly doubling across a tenth of a unit of
flux -- and the bare-rock arm sits at the cold end of exactly that behaviour.
A large bracketed set is the correct answer to a genuinely uncertain question,
and if the bare arm grows glaciers where the vegetated one does not, that is a
finding about the bound rather than a reason to move the arm.

Over-carving is a budget item, not a lost landscape. A build is regenerated
from the planet code plus a verdict rather than edited, so an over-carve costs
a terrain, hydrography and boundary-condition rebuild and nothing else. That
is how `carved-zoned` was abandoned once its verdict turned out to have been
computed on antipodal climate: its carves could not be un-cut *within that
build*, and the build was replaced wholesale.

On iteration 2, the already-carved set should be re-evaluated against the new
climate and the number that would no longer have carved reported. That number
is the overshoot, and it is the honest measure of how much the first pass
cost.

## 5. Where the pipeline currently stands

**`world_state.json` in the project root, generated by `scripts/world_state.py`.**
It is regenerated from the artifacts after anything that changes a build, a
run or a verdict, and it is the only place current values are written down.
They are not repeated here on purpose: which build is active, how many basins
survive, and what the mean temperature is change every iteration, and prose
restating them is wrong within the day.

What does not change between iterations, and is worth knowing before reading
the generated state:

- **The land fraction is fixed.** Carving changes where water leaves a basin,
  not where the coast is, so the land/sea split is identical across every
  build. Anything that moves it is a bug.
- **This world is far more endorheic than Earth**, by a large factor rather
  than a little, in every iteration so far. That is the fact the hydrography
  component exists to handle, and it is why lakes and evaporite matter to the
  climate here when they would be a detail on Earth.
- **Each carve iteration darkens the world and warms it**, and **vegetated
  albedo moves roughly twice as far as bare-rock albedo** under a lithology
  change confined to closed-basin fill -- both are section 4's couplings, and
  the second means bare-rock figures systematically understate what the
  climate will see.
- **The basin catalogue is stable across builds.** Basin ids are computed on
  the pre-conditioning surface, so a verdict computed against one build still
  refers to the same basins in the next. `manifest.hashes.basinCatalogue` is
  the check.

## 5b. Where the planet sits, and why

The flux was set by a derivation rather than inherited, on 2026-08-16.
Recorded because the previous target had propagated through six files with no
statement of where it came from -- a number justified by its own persistence
-- and had reached the convergence criterion before anyone asked.

**The mean is chosen for habitability by latitude band.** Summer and winter
temperature per band were measured on converged runs and projected across
candidate means; the trade is tropics against poles, and the tropics win on
area -- a lower flux moves half the land's warm-season monthly means back
inside the design comfort band at the cost of a tenth of it going from harsh
to extreme, keeping most land in Earth-like conditions and leaving the polar
margins severe but small. The derivation is CODIFIED: `derive_design_flux.py`,
step `design_flux`, writes `exoplasim/analysis/design_flux.json` with its
thresholds declared in the script before it runs. Re-run it on a new terrain
rather than reinventing it, and know its first finding: the
comfort-maximizing rule alone prefers a much lower flux, and the recorded
choice is optimal only under a cap on the cold-extreme land fraction -- the
quantity "severe but small" gestured at and never fixed. The cap must be
DECLARED in advance by the next re-derivation; the purged artifact inferred it
from the answer, which is what CLIM-30 exists to undo.

**Glaciers are decoupled from the mean, and that is the load-bearing
finding.** Cooling is close to useless for making them: summer amplification
is nearly flat with latitude while winter amplification is several times
larger, so cooling buys brutal winters and barely touches the summers that
control ablation. Freezing a polar summer needs about 19 K of global cooling,
which would put polar winter near uninhabitable.

Glaciers come from relief instead. Summits a kilometre or so above their grid
cell sit at freezing in summer, so every band poleward of 30 degrees has peaks
below freezing year round, and the equilibrium line moves roughly 700 m of
elevation over the stellar cycle, so they advance and retreat visibly on a
generational rhythm.

**The poles are the worst place for glaciers on this world**, which inverts
the terrestrial intuition and does so from first principles. At 32 degrees
obliquity the summer-solstice pole receives nearly double the equator's daily
mean insolation, both polar caps are land, so there is no ocean buffer and no
inherited ice to reflect it away, and polar summers come out hotter than the
tropics while polar winters are severe: a seasonal range beyond anything
terrestrial. Earth's poles are cold in summer largely *because* ice is already
there; Vesper never established that feedback. The cold-summer band, and
therefore the glaciers, sit in the mid-latitudes -- which band is this world's
particular continents rather than a general rule; the inversion itself is
mechanism.

**The cycle pushes the mean up, very slightly.** Predicted the other way
round and measured wrong: ice-albedo feedback was expected to make T(f)
concave, and converged points say it is convex over the range the cycle spans
-- there is not enough ice for the ice term to lead, while the Planck response
and the water-vapour feedback both grow with temperature and both push convex.
The magnitude is negligible; the mechanism is corrected because a reason that
is backwards will be reused. The flux-to-kelvin slope is canonical in
`lib/sensitivity.py`, and the superseded per-segment values that first
measured the curvature must not be quoted as conversion factors.

One thing here is NOT settled and must not be carried as though it were. The
glacier result is unmodelled: the GCM reports `glac = 0` everywhere because it
cannot see a mountain peak inside a 300 km cell, so the question is entirely
sub-grid. `notes/glacier-rough-pass.md` does the sub-grid integration properly
and **confirms the mechanism** -- the area-weighted mean glacier latitude sits
in the mid-latitudes in every case and falls as the world cools, because
relief rather than latitude sets where ice survives -- but it remains a rough
pass with an assumed lapse rate and no mass balance, so its areas are an upper
bound and are not properties of the world. The detailed treatment belongs with
the downscaling and sub-grid sampling machinery, not in the climate loop.

## 6. What happens next

Three nested loops and then a resolution change. The order matters in places
where it is not obvious, so those places are called out rather than left to be
rediscovered.

**A. The terrain loop, one turn of which is an ITERATION.** A generation on
the carve list, then the commissioning of the build it produces: hydrography
and boundary conditions rebuilt on the new terrain, bootstrap, derived fields,
baseline, verdict. Hydrography has to be rebuilt *after* the carve and before
the climate, because the carve changes the drainage the climate is integrated
over.

The block below is the COMMISSIONING, and it is the same block whether the
build is new or is being re-commissioned; a re-commissioning skips the
generation and starts from a build already in `source/`. Each step reads the
one before, and `check_consistency.py` reports which are still missing at any
point:

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
current state, with what is already present marked skippable.

The other half is **the carve gate: nothing outstanding may still move an
artifact the verdict is computed from.** It is not a ceremony about
irreversibility -- a wrong verdict is recoverable, because a build is replaced
wholesale -- and it is not a computation. **Read `TASKS.md` and decide.** The
open set is around a dozen rows and each says what it is waiting on.

It was briefly computed, as "open tasks touching a step upstream of
`carve_list`", and the arithmetic is why that was abandoned: 58% of the graph
is upstream of the carve, so the filter selected 10 of 11 open rows. The single
row it dropped was already marked blocked on the downscaling pass. A filter
that reproduces its input is not a gate, and dressing the judgement as a script
made it look decided when nothing had been decided. A task's `[step: <id>]`
marker stays, as the annotation it always was: where the work lands, not a
verdict on whether it blocks.

The three at the top and the two at the bottom are the ones most often skipped
and the ones that cost most when they are: `rebuild_binaries.py --verify`
because nothing else in the list notices a stale binary, `world_state.py` last
because it reads every artifact above it, and running it earlier records a
state that no longer holds.

Everything below the climatology in the register -- the derived surface
classes, brine paths, the phosphorus budget, weathering fluxes, both
prospectivity fields -- is regenerated after the baseline and before anything
quotes it. They are not in the list above because they gate nothing in the
loop; they are consumers, and rule 7 governs them: when the climatology moves
they are worthless, and regenerating them is a step rather than a task.

**The cycle run is in that list because leaving it to prose lost it once.** A2
said it belongs before the verdict, nothing tracked it, and a full verdict was
taken on the mean climate before anyone noticed. It is now CYC-1 as well.

**The first climate run on a new terrain is a bootstrap, and its numbers are
not the baseline.** Three of the surface fields a run consumes cannot be built
without a climatology: the lakes, the lake compositing inside the albedo (174
to 176), and soil water capacity (229), which comes from a soil weathered
under a climate. Everything else is a pure function of the terrain. So the
loop is entered by running the model on the fields that do not need it, and
the run exists to produce the climatology the rest need. Written here because
discovering it one field at a time costs a run each time.

**`model.soil_water_source` has to say `uniform` for the bootstrap and
`pedology` for the baseline, and the flip goes between them.** With it set to
`pedology`, `run_exoplasim.py` requires surface code 229 and refuses to start
without it -- the one field the bootstrap exists to make possible. Flipping
the key changes a value `continue_exoplasim.py` compares against the run's
manifest, and it will not resume across that, so the bootstrap has to be
finished -- converged, seasonal window run, climatology built -- before the
key moves.

**Lakes and soil water have to be in place before the run whose climatology
the verdict uses, not merely before the verdict.** Both change evaporation,
the numerator of the carve criterion, so a verdict taken on a climate that
lacked them is a verdict on the wrong evaporation. Their weights are not
equal -- `analysis/error_budget.json` ranks them -- so build both, and expect
the lakes to be what moves the answer.

**Leave `baseline_climatology` null until the BASELINE exists, and pass the
bootstrap's climatology explicitly with `--climatology` to each field builder
that needs it.** Every consumer raises rather than guessing while it is null,
which is the designed state and not a gap. Naming the bootstrap there instead
would hand it to every defaulting consumer -- `carve_verdict.py` among them --
and a verdict defaulted onto a climate with no lakes and no soil water is a
verdict on the wrong evaporation. `config/planet.yaml` says the same at the
key; repoint it at the baseline's REGULAR climatology once that run finishes.

Seasonal snapshots are written by default, and `analyze_climatology.py` needs
the orbital phase they carry. A segment run with `--no-seasonal-output` has to
be extended before it can produce a climatology, which has happened once.

**The flux is re-derived on a new terrain, not carried.** Terrain moves land
albedo, so the flux that lands the design mean moves with it. Take two
converged points spanning the target and interpolate; never extrapolate from
one, which this project has paid for twice, and count a point only if it
clears the remaining-offset test in `assess_convergence.py` rather than the
drift criteria alone. Which knob moves the flux is in loop C, and depends on
whether the calendar is locked.

**Change one thing at a time between runs, and prefer that to hitting a
temperature.** The bootstrap and the baseline differ in their surface fields,
and if they also differ in flux then nothing measures what those fields are
worth: hold the flux, measure the surface, and move the flux afterwards on a
slope measured with the surface in place. The temperature target is a weaker
constraint than it looks -- it is itself a property of the terrain it was
derived on, and it cannot be held to better than the terms already
outstanding.

**The run count to a carve list.** Corrected 2026-08-19; it used to say three,
written before section 4's intersection carve was costed. The bracket needs
BOTH bounding climates, and the cold one is a commissioning of its own because
its lakes, albedo compositing and soil water have to be built from a bare-rock
climatology; it cannot inherit the vegetated ones.

| | converged runs |
| --- | ---: |
| flux bracket, two points spanning more than 5 K, vegetated branch | 2 |
| vegetated baseline at the design flux | 1 |
| bare-rock arm: its own bootstrap, then its baseline | 2 |
| stellar cycle, per A2 | 1 |

Six. The bootstrap is one of the bracket points rather than a seventh, since
any converged climatology will do for fields that are themselves about to be
rebuilt. Three of the six pair naturally and one cannot: the two bracket
points are independent, the vegetated baseline and the bare-rock bootstrap are
independent, and so is each arm's baseline given its own derived fields. The
cycle run is the exception and A2 says why: it has to follow a baseline that
is FINAL, not merely converged. Everything between the runs is functionally
free -- nearly every other step in the graph costs minutes or seconds -- so
plan in runs and ignore the rest.

**Seed the arms rather than cold-starting them.** `run_exoplasim.py
--restart-from` changes only the spin-up path and not the equilibrium, so the
bare-rock arm starts from the vegetated baseline's restart and relaxes across
the endmember gap instead of from nothing. Expect that relaxation to be slow:
a step change in forcing relaxes on a longer timescale than a settled run
drifts, measured at several times the drift fit on the one run that has done
it.

**Know which base you are on, because a pre-carve build restarts the count.**
The first verdict taken on one is iteration 1, whatever a previous line of
builds had already carved, and the overshoot measurement in section 4 applies
only where the terrain being re-verdicted is itself the carved one.

**A3. A forcing change costs a PREDICTION and an A/B, not an iteration.**
Rewritten 2026-08-18; what stood here said to land one forcing change per
iteration so each could be attributed. But attribution cannot gate anything --
a correct term goes in because it exists, per section 7's physics-is-not-a-knob
rule, and a wrong verdict is recoverable because a build is replaced wholesale
-- so the old rule spent converged runs buying information that could never
change a decision. Attribution pays in exactly one place, diagnosing a
surprise, and that cost is CONTINGENT; serializing iterations pays it in full
every time against a risk that materialises occasionally.

So buy it the cheap way. A converged run measures a RESPONSE in kelvin and is
expensive; a short A/B measures a FORCING in W/m2 and is not: branch two arms
off one common restart, run a few orbits, difference the flux diagnostics, and
convert to kelvin with the canonical slope in `lib/sensitivity.py`. Every term
can be attributed this way while the whole bundle still costs ONE converged
run.

Four things make that A/B trustworthy, and none is optional:

- **Every forcing change lands with a quantitative prediction of its own
  effect**, stated before it is run, with what result would mean "wrong". That
  is section 7's rule that a check needs a right answer, applied to physics.
- **Both arms use the SAME BINARY and differ only by a namelist key.** This is
  what the no-op-until-enabled convention is for. Two binaries would confound
  the term with the rebuild, and the low-I/O patch changes the restart layout,
  so arms built either side of it cannot share a restart at all.
- **Both arms branch from ONE restart.** Run-to-run spread on a converged pair
  is larger than several of the terms being tested; a shared initial condition
  turns the comparison into a paired one.
- **The segments are labelled as diagnostics**, so a short A/B tail never
  enters a convergence window or a climatology.

**Bundle the terms, then check the SUM against the sum of predictions.** If
they agree, no attribution was needed and none was bought. If they disagree,
bisect, and only into the subset whose predictions were soft. What bundling
genuinely risks is `notes/failure-modes.md` class 15, two errors that nearly
cancel; the mitigation is the per-term prediction above, not serialization,
because two cancelling errors hide just as well in a serial sequence nobody
predicted the size of.

**What survives from the old rule is a budgeting fact.** `build_surface_albedo.py`
has a `modelled` mode that takes tree cover from an LPJ-GUESS `fpc.out`
instead of asserting a uniform vegetated endmember, and it is the mode this
project should end up in. The two land-albedo endmembers are as far apart in
absorbed flux as the largest sweep this project varies on purpose, so adopting
it is a forcing change of the first rank and **it requires the flux to be
re-derived**. That is a cost to budget for, and it is the same cost whether or
not a carve shares the iteration. The corollary is the answer to "when does
the biosphere run": not before the carve, because its driver is built from a
climatology and a soil the carve replaces; after the re-baseline, on final
terrain, where its output can be adopted deliberately in the iteration after
that.

**A4. Dust is not carve-neutral, and the carve waits on one climate run.** Its
global-mean forcing is negligible in kelvin; the local effect is not. The
question splits in two and only one half is settled.

The LAKE half is done. Dust dims the lake, Penman evaporation falls, and more
basins overflow: a dust-free verdict UNDER-carves, modestly and in the
recoverable direction. That was settled offline by putting the dust SURFACE
forcing through `carve_verdict.py --dust-forcing`; the surface balance is the
operative one, because the top-of-atmosphere sign over bright basin fill is
energy retained in the atmosphere, not delivered to the ground.

The CATCHMENT half is larger, runs the other way, and cannot be reached
offline, because it is a precipitation response and not a surface energy
balance: reduced rainfall cuts the runoff the criterion divides by, which
OVER-carves, and over-carving is the irreversible direction. It needs one
prescribed-dust climate run, specified with its gates and a prediction in
`aeolian/notes/prescribed-dust-run.md`. The ORDER here is a dependency and not
a serialization: the verdict divides by a runoff that dust moves, so it has to
be taken on a climatology that already has dust in it. Run it, take the
verdict on its climatology, then carve. **Do not record a carve list as
dust-independent.** `notes/dust.md` has the numbers and
`aeolian/analysis/dust_runoff_sensitivity.json` the conversion from a
precipitation change into basins.

**B. The soil and biosphere loop, at T42.** For a given climate:

1. `pedology/scripts/build_soil.py`, with no biosphere on the first pass.
2. `biosphere/scripts/build_lpj_driver.py`, then `run_lpj_guess.py`.
3. `build_soil.py --soil-carbon <run>/cpool.out`, then LPJ-GUESS again.
4. Repeat 3 until the criteria in `pedogenesis.yaml` are met.

Check whether step 3 moves anything before assuming it needs iterating:
LPJ-GUESS computes its own soil carbon internally, so the pedology organic
feedback may be second-order.

**A2. The stellar cycle, and where it belongs in the order.** The cycle run is
NOT a final flourish. It has to come after a baseline that is FINAL, not
merely converged -- a cycle is variance about a mean, and centring it on the
wrong one describes a world that will not exist, so anything still outstanding
that moves the mean has to land first. But it belongs BEFORE the carve
verdict, for a reason that is easy to miss.

Carving is irreversible: outlet incision does not undo when the warm phase
returns, so every cold, wet excursion that pushes a basin to overflow carves
it permanently. **The terrain therefore ratchets toward the state implied by
the cycle's wet extreme, not its mean**, and a verdict taken on the mean
climate systematically under-carves, because overflow is a threshold process
and the wet phase contributes disproportionately. Bedrock incision is slow
against both cycle periods, so it is the integral over many cycles that
matters: the verdict wants a climate between the mean and the wet extreme,
weighted by time spent overflowing.

**The two components ratchet differently, and that is the reason there are
two.** The medium one is too fast for ice to follow and too fast for an outlet
to incise, so it contributes only through the tail of its distribution. The
long one is slow enough that glaciers equilibrate and long enough that a wet
excursion is sustained, so it is the one that actually cuts. And because the
two periods are non-commensurate, the deepest minima differ in depth rather
than repeating, so the landscape ends up recording which past minima were
severe. A cycle run therefore has to be long enough to sample that: its length
is set in periods of the LONG component, not the medium one.

The cycle run also measures the damping factor, currently known only as a
factor-of-2.5 range, which converts any future amplitude choice into a climate
without another run -- measuring it once is worth more than the run that
measures it.

Sequence, then: baseline at the chosen mean, cycle run centred on it, verdict
on a cycle-informed climate, and only then the terrain loop.

**C. The vegetation-climate loop.** `build_surface_albedo.py --mode modelled`
turns the run's foliar cover into surface albedo and forest fraction, then the
climate runs again on it. Two checks belong here and neither is optional.

*The flux.* If the modelled land albedo differs much from the assumed value
the world may leave the design mean derived in section 5b. **Move the flux by
changing the star's luminosity, not the orbit.** The year length depends on
the semimajor axis, and the year is compiled into LPJ-GUESS, so moving the
orbit forces a rebuild and a driver regeneration. With the axis locked,
F = L/a^2, so luminosity is the free parameter and the calendar does not move.

The knobs swap at the lock, and the lock is what makes this a rule rather than
a preference. Before the calendar is fixed, the semimajor axis is the better
knob: it keeps the star itself untouched, so `k25v` stays valid, where
covering the same flux range with luminosity moves the effective temperature
far enough to leave the window the spectrum was interpolated in. Once
LPJ-GUESS is compiled against a year, that knob costs a rebuild and luminosity
is the only cheap one left. A bracket taken on the axis varies two things
rather than one, because the year moves with it, so the measured slope carries
a small seasonality effect through the calendar; that is second order for an
annual mean, and it is stated rather than implied.

Mind the scaling, and mind the spectrum. L scales **1:1** with F at fixed
orbit, not 2:1 as an earlier revision of this document said, and at fixed
radius L ~ T^4, so the effective temperature moves as F^(1/4):

| flux change | luminosity | effective temperature | band-1 fraction |
| --- | --- | --- | --- |
| 1.5% | 1.5% | 18 K | ~0.9% |
| 4% | 4% | 49 K | ~2.3% |

`k25v` is interpolated between BT-Settl models 100 K apart, and that grid step
is worth 4.7% in the fraction of flux below 0.75 um. So small adjustments are
free, but **past about 2 to 3% in luminosity, re-run
`build_stellar_spectrum.py` at the new effective temperature**, and past
roughly +35 K the target leaves the interpolation bracket entirely and the
pinned SVO grid points have to change with it. The spectrum exists to get snow
and ice albedo right, so letting it drift silently would undo the reason it
was built.

*The carve verdict.* It was taken on assumed vegetation. Real vegetation
changes the climate, which changes evaporation over catchments, which can
change the verdict. Re-run it; if basins flip, back to loop A.

**D. The resolution change.** Only after A, B and C have settled.

1. Climate at T85, with the T42 vegetation regridded as its boundary
   condition.
2. Soil and biosphere at T85 on that climatology, loop B again.
3. Climate at T85 again with T85 vegetation, unless step 2's vegetation turns
   out close to the regridded T42 field, which is a cheap comparison worth
   making first.

**The biosphere never decides the resolution.** LPJ-GUESS gridcells are
independent columns, so its cost is linear in land-cell count and trivial at
either resolution; the climate model's cost is what decides. Running it at
T85 on T42 forcing would resolve detail that is not in its input. Expect T85
to lower total NPP, and treat that as a resolution bias rather than a result:
productivity saturates with water, so averaging the forcing before the model
sees it inflates the answer.

**E. The cycle reaching the biosphere, last, on the settled world.** A2 says
when the cycle RUNS; what is left for the end is letting LPJ-GUESS see it.
Build the climatology with `--per-year` and pass the sequence to
`build_lpj_driver.py --climatology y0.nc y1.nc ...`. That matters because
productivity responds annually and saturates, so a run on the cycle mean
over-predicts it; it does *not* reach survival thresholds, which LPJ-GUESS
gates on a twenty-year mean and which therefore smooth the cycle away almost
exactly. Then the regional products, biomes, Koppen and lake maps.

Three things about the calendar and the resolution that are easy to get wrong,
each worked through in `biosphere/README.md`:

- **Spin-up is counted in simulation years**, and a simulation year is about
  half an Earth one, so the shipped default buys roughly half the absolute
  time Earth practice assumes. Year counts are scaled UP and annual sums
  scaled DOWN by the reciprocal factor; the two directions are opposite, both
  lists are named in `build_vesper_pfts.py`, and the factor is derived there
  from the configured orbit rather than written anywhere.
- **T85 redistributes precipitation, it does not merely resolve it.** Steeper
  relief means stronger orographic ascent and deeper rain shadows, so lee
  basins dry and windward coasts wet. That reaches the carve verdict, which
  should be re-taken at T85 rather than carried over from T42.
- **The 30-hour day widens the real diurnal range and LPJ-GUESS cannot see
  it.** `dtr` reaches only the biogenic VOC scheme, and every cold limit runs
  through a twenty-year mean, so there is no daily-minimum plant mortality to
  trigger.

## 7. Conventions this project holds to

Provenance travels with every artifact: config, input hashes, software
versions and the terrain hash go into each run manifest and each analysis
report.

Thresholds are fixed before results are seen. The albedo bracket's convergence
and marginality criteria were written into the script's docstring before any
run finished, and were then applied against a result that cleared one of them
by 0.05 K.

Estimates that cannot be verified are bracketed rather than guessed, and the
bracket is reported. This is how albedo, evaporation and the carve verdict are
all handled.

Claims are checked against the artifact rather than the documentation. Several
findings in this project's history came from comparing two products that were
supposed to agree and finding they did not.
