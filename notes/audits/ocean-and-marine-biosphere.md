# The ocean: what the model already has, what a marine biosphere is worth, and what carving does to the delivery number

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from several real disciplines: **Vesper is a fictional planet and this is
engineering work on the simulation of it.** Words like ocean, phytoplankton,
chlorophyll, salinity and productivity name modelled quantities of an invented
world throughout. Nothing here is Earth science and nothing here is biology.

Audited 2026-08-21, commissioned to answer whether this project should go beyond
ExoPlaSim's slab ocean and what marine ecosystem model would fit if it did. The
question turned out to have three separable halves that the candidate model
names conflate: the COUPLING CLASS, the CIRCULATION HOST, and the ECOSYSTEM
TIER. Each is answered below.

Four findings are about what the vendored model already contains and nobody has
read; one is a pre-carve number read as a state, in a config block that is
load-bearing; and one is a pricing that settles the coupling class. Section 9 is
a second reading, which added the model class the first one skipped. The tasks
they became are named at the end.

Numbers are labelled with what they were measured on. The climate figures come
from the BOOTSTRAP climatology, which is the only one on this build, so read
them as scale rather than as the baseline.

---

## 1. PlaSim's slab is a degenerate case of a code path that is already general

`vendor/exoplasim/exoplasim/plasim/src/oceanmod.f90:15`:

    parameter(NLEV_OCE = 1)           ! Number of Layers

Every ocean array in the module is dimensioned on it -- `dlayer`, `vdiffkl`,
`hdiffk`, `ysst`, `ymld` -- and line 255 branches on `NLEV_OCE > 1` to run
vertical diffusion between layers, with `vdiffkl` the per-layer coefficient
already in `oceanmod_nl`. Line 229 sets `dlayer(NLEV_OCE) = mldepth`, so the
single-layer case is the general case collapsed.

**A multi-layer column ocean is a compile-time constant and two layer choices,
not new physics code.** It is the cheapest thing in this document that acts on a
term `config/planet.yaml` already declares: `mixed_layer_depth_m: 50.0` is
ExoPlaSim's default arriving unnamed, it sets the seasonal amplitude, and this
world's orbit is 182.8 Earth days, so one slab damps the simulated seasonality
about twice as hard as Earth's ocean does.

It changes the binary, so CLAUDE.md rule 4 applies in full: `NLEV_OCE` is a
source constant, and every one of the compiled executables goes stale together.

## 2. There is already a channel for a prescribed ocean heat transport, and it is wired at both ends

`oceanmod_nl` carries `nfluko`, 0 for none, 1 for heat-budget, 2 for newtonian.
At `nfluko == 1`, `oceanmod.f90:323` reads a field:

    call mpsurfgp('yfsst',yfsst,NHOR,14)

`surfmod.f90:227` maps that name to **surface code 903**, and `code_surf_file`
resolves it to `N064_surf_0903.sra` at T42. That is the same `.sra` machinery
this project already writes codes 129, 172, 173, 174, 175, 176, 212, 229 and
1811 through. The field is monthly, `NHOR x 14`, in W/m2. `getflxco` interpolates
it to the timestep and `addfc` applies it, under its own header comment:

    add oceanic flux correction (prescribed advection)

`addfc` distinguishes the sea ice cases internally, and `mkiflx` handles the
residual into the ice model.

**The other end is written too.** `exoplasim/scripts/close_ocean_energy.py`
postprocesses the `ocean_output` stream and reads code 903 among its terms. Both
ends existed and both were null: no code 903 field had ever been written, and
the closure carried the term in the list of things it required to be identically
zero. Section 11 is the end-to-end verification and what it changed.

What this means for cost. An ocean heat transport reaches the model as ONE new
pipeline step writing one `.sra`, plus one namelist key on the
`KEY@namelist` mechanism `run_exoplasim.py` already has. What it does NOT supply
is the field itself, which is finding 6's problem.

`analysis/error_budget.json` records no-q-flux as STRUCTURAL with "no cheap
version exists". The cheap version of the BOUND is `nhdiff`, which
`notes/audits/absent-and-inherited-physics.md` finding 1 established. The cheap
version of the TERM is this, and it was not known when that finding was written.

## 3. LSG was vendored, and it was not the switch it looked like

`vendor/exoplasim/exoplasim/lsg/src/lsgmod.f90` was 9,313 lines of a real
three-dimensional geostrophic ocean. `oceanmod_nl` carries `nlsg` and `naomod`,
the atmosphere-to-ocean step ratio, and `most.c` knew how to build it. It looked
like the answer. Four things in the source said otherwise.

- `most.c:2415`, comment and code together: "LSG works currently only with T21
  PlaSim", and it FORCED `Resolution = RES_T21`. Production here is T42 with T85
  in loop D.
- `plasim/src/cpl.f90` hardcoded the coupler grids as `nxa=64, nya=32` on the
  atmosphere side and `nxo=72, nyot=76` on the ocean side, and the planet as
  `parameter(radea=6.371E6)`. This world's radius is 1.20 Earth.
- The LSG grid is a fixed 72 x 76 x 22 present-day Earth configuration, and
  `lsg/dat/` shipped Earth bathymetry, an Earth runoff map and Earth restart
  fields to go with it.
- `oceanmod.f90:437` is a hard error: LSG coupling requires
  `n_days_per_year = 360`. That is a count of SIDEREAL DAYS PER ORBIT and is this
  planet's 146, not a setting, so the requirement can never be met here.

**Adopting LSG meant forking it to a depth comparable to the LPJ-GUESS port and
then running the atmosphere at T21.**

**The four constraints were accepted and the source is deleted.** world-6ak and
world-cmz removed the whole `lsg/` tree, `plasim/src/cpl.f90` and `most.c` as
trees nothing builds, reads or ships. `plasim/CMakeLists.txt:183-184` still links
a four-line `src/lsgmod.f90` stub and `cpl_stub.f90`, so the `nlsg > 0` branch
compiles and refuses; `oceanmod.f90:418-434` carries the note at that branch.
The bodies are recoverable from git if the verdict is ever revisited, and this
section is the list of what a revisit would have to re-derive rather than
re-discover.

The candidate this project actually carries for the resolved tier is
`vendor/cgenie`, under OCN-3. It is not an adopted component: nothing reads it
and it does not build where it stands.

## 4. The exorheic fraction is a pre-carve LIMIT, and the drainage geometry gives no ceiling at all

CLAUDE.md rule 9 names this trap and it is easy to walk into: a pre-carve build
is a limit, not a state.

`hydrography/data/precarve-craton/hydrography_report.json` gives
`ocean_draining_fraction_of_land = 0.23777`, measured on `precarve-craton`. That
is the number with ZERO carving, which is what is physically sitting in
`source/` at the start of every cycle. Carving only ever moves it up on a given
terrain, so it is a FLOOR for this build; a new generation can create new
basins, so it is not a monotone trend across iterations.

Walking the spill graph in `basins.nc`, measured on `precarve-craton`: all 3,621
basins have an exit -- none has `spill_region < 0` -- and every `spill_target`
chain terminates at -1, the world ocean. **The all-carve limit is therefore
exactly 1.0000, and the geometry constrains nothing above the floor.**

Carving where `critical_aridity_index >= t` and following the chain through
carved basins only, so that a carved basin spilling into an uncarved one stays
endorheic. The last column indexes exorheic land per unit ocean area against
Earth at 87 percent exorheic over a 0.29 land fraction, which is 0.3554:

| case | exorheic land | per unit ocean area | x Earth |
| --- | ---: | ---: | ---: |
| no carve, what is in `source/` now | 0.2373 | 0.180 | 0.51 |
| t = 5.0 | 0.3286 | 0.250 | 0.70 |
| t = 2.75, the median `critical_aridity_index` | 0.4277 | 0.325 | 0.91 |
| t = 2.0 | 0.5254 | 0.399 | 1.12 |
| t = 1.0 | 0.7536 | 0.573 | 1.61 |
| all-carve limit | 1.0000 | 0.760 | 2.14 |

The no-carve figure here is 0.2373 against the report's published 0.23777,
because it is computed as `land - sum(catchment_km2)` over the report's own land
denominator; the difference is area bookkeeping and moves nothing.

**This is the registered prediction for any future marine productivity work**,
in place of a single number: riverine delivery to this world's ocean, per unit
ocean area, lies between about half Earth's and about twice Earth's, and where
it lands inside that span is a carve-verdict question end to end. A marine
result quoted against the pre-carve figure alone is quoted against the floor.

**The same number is read as a state in a load-bearing config block.**
`config/planet.yaml`'s `ocean` block argues its salinity direction from "the
exorheic share of the planet's surface is under half Earth's per unit ocean
area", which is true only at the no-carve end of the table above. The block is
DECLARED and carries Earth's 34.7 psu, so nothing numeric moves today, but that
sentence is the justification for bracketing the freezing point, and the
direction it derives -- fresher ocean, warmer freezing point, more simulated sea
ice -- does not survive the middle of the table.

Two things make the solute question harder than the water question, and neither
is settled here. Carving turns a basin from a solute SINK into a solute SOURCE,
because what accumulated on its floor becomes available to a through-flowing
river. And `brine_paths.py`, the evaporite lithology classes and the phosphorus
routing in `pedology/scripts/phosphorus_budget.py` are all constructed on the
pre-carve terrain, so they describe the floor case as well.

## 5. What a marine biosphere is worth as a climate input, priced

Priced against the BOOTSTRAP climatology, measured 2026-08-21: planetary albedo
0.2631 and 0.806 K per W/m2 of top-of-atmosphere forcing, through
`lib/sensitivity.py` at its 202 K per flux ratio slope, which `verify()` passes.

Four pathways by which a modelled ocean ecosystem could reach this world's
climate. One is priceable, two are structurally absent from the model, and one
is gated on a switch that is off.

**5a. Chlorophyll darkening the simulated sea surface.** `seamod.f90:13` sets
`albsea = 0.069`, and `seamod.f90:155-158` applies `doceanalb(1)` and
`doceanalb(2)` as two-band namelist SCALARS everywhere `dls < 0.5`. There is no
spatial ocean albedo field and no surface code for one, so this pathway costs a
new surfcode as well as a field.

The simulated pigment signal is confined to the water-leaving part of the
reflectance; the Fresnel term does not see it. Bracketing the Earth-Sun
water-leaving albedo change at 0.0005 to 0.005, DECLARED rather than sourced,
and the surface-to-planetary attenuation at 0.25 to 0.5, which is the error
budget's 0.5 with the factor-of-two `lib/sensitivity.py` refuses to absorb:

**0.014 to 0.28 K**, over an ocean that is 0.5683 of this world's surface.

The top of that span is below the smallest item currently resident in the error
budget, lakes composited into albedo at 0.61 K and playa albedo at -0.52 K.

It is already discounted for this host, and by a computable amount. Pigment
absorption sits in the blue and the red, and the water-leaving signal dies past
about 0.7 um, so the whole effect scales with the band-1 flux share.
`world_state.json` gives 0.3824 for the k25v spectrum against 0.5374 for a
5772 K Planck at the same 0.75 um split: **this star gives a blue-and-red
pigment 0.712 of its Earth-Sun leverage before anything else is considered.**

**5b. Chlorophyll heating the simulated mixed layer.** This is the pathway that
carries the larger effect on Earth, and it is structurally absent.
`oceanmod.f90:mksst` adds the whole atmospheric heat flux to level 1:

    zsst(:,1)=zsst(:,1)+(yheat(:)+yfldo(:))*zcpsdt/ymld(:,1)

Lower levels are reached only by `vdiffkl`. So nothing resolves where in the
water column the shortwave is absorbed, **and finding 1 does not fix it**:
adding layers gives the column a vertical structure but still deposits all the
shortwave in the top one. A penetration profile is separate code. Unpriceable
until it is written.

**5c. A biogenic sulfur pathway to cloud droplets.** `radmod.f90` has no droplet
number, no effective radius and no aerosol-cloud term of any kind, and
`aeromod.f90` is the FFSL transport scheme only; the dust work added a direct
radiative effect and nothing indirect. There is no term to hook a biogenic
aerosol source onto. This costs a module rather than a field.

**5d. Ocean carbon.** Zero leverage as configured, because `pCO2_bar` is
prescribed and `config/planet.yaml` already says to read it as an assumption.
Conditional on a prognostic carbon cycle, which would also need the
carbonate-silicate side that `pedology/scripts/weathering_fluxes.py` half
covers, air-sea repartitioning is worth, at the same 0.806 K per W/m2 and with
the broadband scheme's measured 0.4 of a line-by-line CO2 answer as the low arm:

| repartitioning | at 0.4x scheme | at 1.0x |
| --- | ---: | ---: |
| 50 ppmv | 0.18 K | 0.45 K |
| 100 ppmv | 0.35 K | 0.87 K |
| 200 ppmv | 0.63 K | 1.59 K |

That is budget-scale. **It is the only pathway on which a modelled marine
ecosystem is a first-order control of this world's climate, and it is gated on a
switch this project has not thrown.**

**The verdict.** A marine biosphere is a LEAF of the pipeline, in the sense
`docs/src/pipeline/components.md` uses for the minerals overlay, with one named
condition: it becomes a climate input the day pCO2 becomes prognostic, and not
before.

The consequence worth more than the pricing is a schedule split. The ocean
PHYSICS is inside loop A, because a heat transport field moves the global mean,
which moves the design flux, which moves the carve verdict; it lands before a
carve or it invalidates one. The ocean BIOLOGY is outside loop A, so it can be
built, re-run and explored without costing a re-commissioning.

## 6. The model recommendation

**Coupling class: offline full-flux, on the aeolian pattern.** The ocean model
reads a climatology and writes fields the next climate run is forced with. It is
not synchronously coupled to the atmosphere.

The reason is the resolution ladder, and it is sharper than "a deep ocean is
expensive". A coupled intermediate-complexity atmosphere-ocean model spends 98
percent of its time in the atmosphere; the frictional-geostrophic ocean is
nearly free. So what a dynamic ocean costs is not the ocean. It is being obliged
to integrate the ATMOSPHERE for the thousands of model years the ocean needs,
and that price is paid at whatever rung the atmosphere is running. Offline
coupling is the only class in which the ocean's cost does not climb the
T21/T42/T85/T127/T170 ladder with the atmosphere, and this project stops at the
first rung that passes SPAT-8 rather than at a fixed resolution
(`docs/src/pipeline/costs.md`).

**Two other coupling architectures were considered and eliminated**, recorded
here so neither is re-opened. An offline path supplying only wind stress, winds
and albedo, letting an energy-moisture balance atmosphere compute the heat
fluxes with its own parameterised transport, keeps the partition tunable and
needs no loop -- but the tuning is of an Earth-calibrated second atmosphere on a
30-hour, 32-degree world, which is class 16. A synchronous path with gearing
solves the partition inside the model, and does not reach the top of the ladder:
its shipped form is single-task and compile-time fixed to T21 with ten levels,
its matched-grid convention puts a serial r^4 ocean at 512 x 256, and the only
lever that would make it affordable at a high rung is the one that turns it back
into offline coupling. The arithmetic is in `notes/external-model-survey.md`
section 59.

**Circulation host: cGENIE, and the resolved tier is off the board.** Two
eliminations and one candidate.

The NEMO family goes on one argument: PISCES, MEDUSA and PlankTOM are NEMO
configurations, so adopting any of them means adopting NEMO, whose grid and
configuration machinery are harder to move off Earth than the alternatives'.

**MITgcm goes on the coupling-class argument, applied to itself.** This is
recorded rather than left implicit, so the next reader does not re-open it. The
objection is NOT that MITgcm is cluster software; a coarse global configuration
runs on a workstation, and the runtime parameters for sphere radius, gravity and
rotation rate are exactly the axis this world moves along. The objection is that
the property this project needs from an ocean model is an EQUILIBRIUM, and a
primitive-equation model reaches one only after thousands of model years at a
timestep set by resolved dynamics. That is the same argument that chose offline
coupling at the top of this section, and it does not stop applying because the
model in question is the one being considered. The connectivity MITgcm would buy
is real, but it is only available at a resolution whose spin-up is unaffordable,
so it is not a fallback: it is a benefit priced out of reach. A coarse MITgcm
resolves no more than cGENIE does and still cannot be equilibrated.

That leaves **cGENIE**, which reaches biogeochemical steady state in a day on one
core, already has a published and archived ExoPlaSim coupling, carries the closed
carbon cycle finding 5d needs, and hardcodes planetary radius in four places
while its transport and carbon parameters are fitted to Earth observations.

It is a candidate, not an adoption. Section 9 reads its source and cost the way
finding 3 read LSG's, and OCN-3 establishes viability and measured cost on this
machine before anything is committed to it.

**The consequence of dropping the resolved tier is that connectivity now rides
entirely on resolution.** Section 8b's straits, sills and partial coasts were the
one thing the other candidate was for. With it gone they become a constraint on
how fine a grid cGENIE must run, not a choice between hosts, and if that grid
turns out to be unaffordable the honest outcome is a declared
unresolved-connectivity bracket rather than a quiet acceptance of whatever
36 x 36 represents. That is what makes OCN-18 load-bearing rather than an
optimisation.

**Ecosystem tier: trait-based, and the host decides which one.** MARBL is the
better-engineered library by some distance, being deliberately driver-agnostic
with a documented coupling interface. Its ecology is a small set of fixed Earth
plant functional types, two of which are specific Earth evolutionary inventions.
That is the same problem `biosphere/README.md` already records for LPJ-GUESS:
"the shipped plant functional types are Earth's ... should be declared that way
rather than presented as a prediction". Darwin and ECOGEM are both trait-based,
so in both the modelled community composition is an output of the trait space
rather than an input, which on a world whose nutrient supply is bracketed across
a factor of four (finding 4) is the difference between a prediction and an
assumption. That property does not separate them, and the host
decides which is reachable: Darwin is an MITgcm package and leaves the board with
it, so the tier is ECOGEM, with driver-agnostic MARBL as the fallback wherever
the trait space turns out not to be portable.

An upper-trophic model is forced by a lower-trophic field and cannot precede
one. It is the right tool when this project wants modelled nekton biomass, and
it is not the first step.

## 7. The forcing trap, and the loop it creates

An ocean model forced by a climatology that a zero-transport slab produced, and
restored toward that climatology's own sea surface temperature, returns a heat
transport of zero by construction. The offline model would be reproducing the
answer it was meant to supply.

The forcing procedure therefore has to be declared before the first run rather
than discovered from it: wind stress and bulk formulae, salinity restored
weakly, and temperature NOT restored. It also creates a loop -- transport field,
new baseline, new forcing, new transport field -- and every loop in this project
carries a declared exit predicate in `config/pipeline.yaml` rather than a
stopping habit. That is OCN-5, and it is a decision rather than a measurement.

The hand-over itself has a rule, taken from the eliminated synchronous
architecture, where it is a working arrangement rather than an inference. Terms
that depend on ocean temperature are recomputed live against it -- saturation
specific humidity and hence latent heat, net longwave, and sensible heat --
while only the transfer coefficients and the downward radiative and moisture
fields are replayed. Naive replay of net heat flux severs the negative feedback
that holds sea surface temperature. Evaporation is deliberately not recomputed,
for moisture conservation, with rescaling precipitation and runoff considered
and rejected. `notes/external-model-survey.md` section 59f.

## 7b. Sea ice: ExoPlaSim stays authoritative, and the ocean returns a velocity rather than a state

Adopting an offline circulation host lifts OCN-21's block as a side effect, since
the block is that a slab ocean has no surface velocity to advect ice with and a
circulation model produces exactly that. It also creates a question the heat
transport decision did not have, because both components own sea ice: ExoPlaSim's
`icemod` must run, being the atmosphere's lower boundary, and the circulation
candidate ships `genie-goldsteinseaice` beside it.

**ExoPlaSim's `icemod` is authoritative, and the return path carries the ocean's
surface velocity, not its ice.** Ice stays prognostic in the atmosphere and gains
an advection term driven by a field the ocean is already producing. That is the
same shape as the heat-transport decision -- the ocean returns a transport agent
and the atmosphere keeps the state.

The alternative was available and is rejected. `icemod` already has the channel
for it: `read_ice_surface` reads climatological sea surface temperature, ice
cover and ice thickness as 14-month surface fields, `nice = 0` prescribes ice
outright, and `nfluko = 1` diagnoses an ice flux correction from the mismatch
against the prescribed thickness and clamps cover to it. Four arguments against
using it, the first of which is decisive on its own:

- **It severs the feedback this world's climate question is about.** Ice
  prescribed from a previous pass cannot respond to the atmosphere within a run,
  so a cooling excursion grows no ice. Albedo-driven bistability at the cold end
  is what the design flux, the extreme-cold land fraction and the stellar cycle's
  grand minima all turn on.
- **`nfluko = 1` is a correction toward a target**, which is a knob in the sense
  `docs/src/practice/failure-modes.md` class 16 names. It can bound the term; it
  cannot be used to make the ice come out right.
- **The lag is structural rather than incidental.** Heat-flux convergence
  tolerates arriving from a previous pass because it is a smooth, large-scale
  integral quantity. Ice fraction is a threshold field at a moving margin, and a
  stale margin puts the wrong albedo at exactly the latitudes where the feedback
  lives.
- **It spends what the offline decision bought.** Ice returned from the ocean
  arrives on the ocean grid and must be remapped onto every rung of the ladder,
  and a fraction crossing a threshold is the worst case for conservative
  remapping across a factor-of-twelve scale gap.

The published synchronous coupling does the opposite -- ice temperature, height,
fraction and albedo all come FROM the ocean -- and that is correct for a
synchronous arrangement, where there is no lag. The argument does not transfer.

**What this costs**: `icemod` has no advection at all, so the term has to be
written. The reference is small -- a second-order explicit transport on two
prognostic fields using upper-ocean velocity, with a diffusion term beside it and
no ice momentum equation -- but it is a change under `vendor/exoplasim`, so
CLAUDE.md rule 4 applies. Until it exists the term stands declared in
`scripts/error_budget.py` rather than corrected.

**One consequence for OCN-10's contract**: sea ice is not in the return set. The
returned fields are heat-flux convergence and ocean surface velocity. Ice
fraction, thickness and albedo remain ExoPlaSim's, and the ocean receives them as
forcing rather than supplying them.

## 8. Integration with the terrestrial and biosphere audits

The first version of this audit was written beside, rather than through, the
land-biosphere audits.  Reading the two together exposes seven additional
interfaces.  None changes the choice to scope an offline circulation model and
a trait-based ecosystem.  They change what those words must mean before either
can produce an interpreted result.

### 8a. A climatology is not yet a conservative ocean forcing artifact

EFOR-1 through EFOR-8 establish that a forcing file needs explicit interval,
clock, field and conservation semantics.  The ocean needs the same discipline,
but not the same artifact: it consumes atmosphere-ocean momentum, heat and
freshwater exchanges, routed runoff and sea-ice freeze/melt rather than the
land-column state LPJ consumes.  Wind stress, net and penetrative shortwave,
longwave, sensible and latent heat, precipitation minus evaporation, river
discharge, salt and ice freshwater all need one owner, time basis, support and
sign convention.  A weak salinity restoring term can remain a named numerical
control; it cannot silently replace the freshwater and salt budgets being
investigated.  OCN-10 owns this boundary.

### 8b. The ocean has another support crossing, not a privileged grid

T42 is an operating convention rather than a constraint, while the 10M Orogen
export is the reference support that fully resolves the generated terrain.
The ocean grid therefore cannot be chosen once and called the planet.  It needs
a versioned map from the 10M bathymetry and coastline into wet fractions,
depth/volume, shelves, straits, sills and routed river mouths, plus conservative
maps to every accepted T21/T42/T85/T127/T170 atmosphere.  Narrow connections
and partial coastal cells can change circulation even when cell-mean depth is
unchanged.  OCN-11 extends SPAT-1 through SPAT-5 and SPAT-8/10 to this crossing;
OCN-3's coarse configuration is a candidate whose adequacy must be measured,
not a fixed resolution.

### 8c. An external ocean model does not escape the implicit-Earth audit

Runtime radius, gravity and rotation parameters are necessary but not
sufficient.  Calendar and rate constants, reference pressure, equation of
state, geothermal and tidal assumptions, vertical and lateral closures,
salinity and initial hydrography can still import Earth.  In an ecosystem,
fixed Redfield relations, nutrient inventories, temperature functions,
sinking speeds, irradiance/PAR conversions, latitude classifiers and supplied
community traits do the same.  Latitude remains valid for geometry and
Coriolis; it is not a water-mass, productivity or ecological regime.  OCN-12
applies BIO-22/BIO-29's clock and latitude rules and the repository's
gravity/pressure audit shape before any external ocean or ecosystem model is
accepted.  Section 9d is that audit's sharpest case, because an EMIC's skill is
substantially in a calibration fitted to Earth.

**That audit is `notes/audits/ocean-tier-implicit-earth.md`**, which applies
`model-earth-centrism.md`'s method and vocabulary to both candidate tiers and
answers OCN-12 as scoping: the inventory, and the preconditions a tier would
have to meet before it could be chosen. Its sharpest single item is not a
constant sweep would reach either. `gem_carbchem.f90:99` turns depth into
pressure by dividing by ten, which is seawater density times EARTH gravity with
neither on the line, and it reaches every carbonate equilibrium constant, so the
carbonate compensation depth sits substantially shallower than the model would
place it. The ecosystem half is `notes/audits/ecosystem-tier-ecogem-marbl.md`,
whose section 13 closes the silicon axis.

### 8d. "Reached the ocean" is a boundary condition, not marine nutrition

ANUT-6 routes terrestrial dissolved and particulate export to the coast;
ANUT-3/4 provide particulate and fixed-N deposition, and ANUT-5 defines the
accepted species/time interface. None specifies
what happens after arrival.  A marine calculation needs an initial and forcing
ledger for at least C, N, P, S, Fe and Si, with salt, alkalinity, DIC and
oxygen as coupled chemical state: dissolved/particulate and bioavailable/refractory
forms, estuarine/coastal retention, dust dissolution, hydrothermal and benthic
sources, uptake, remineralisation, sinking, burial, fixation, denitrification
and redox/anoxia. Any resulting marine CH4 or N2O is a named source to the
shared atmospheric trace-gas ledger, not permission for the ocean component to
set an atmospheric abundance independently. Carving changes both river delivery
and release of stored basin solutes, so the no-carve and all-carve values are
boundary brackets rather than marine productivity endmembers. OCN-13 begins at
ANUT's terminal-export boundary and closes the internal ocean ledger without
making ANUT own marine biogeochemistry.

### 8e. Trait-based is not yet Vesper physiology

A trait-based ecosystem makes community composition emergent only inside the
trait space it is given.  OCN-14 must declare that space with covarying strategies and retain a
transplanted-Earth case as a labelled bracket rather than tuning traits to a
desired productivity map.  It also imports PCAR-1/2 and BIO-25's lessons:
rates need an absolute-time versus orbital-phase registry; light limitation
needs K-star spectrum/window-weighted photon supply and water-column spectral
attenuation; and sinking, buoyancy and mixing responses must use Vesper gravity
and the accepted circulation state.  Nutrient limitation consumes OCN-13's
species and cannot infer general sufficiency from a C-N-P subset.

### 8f. Marine output needs acceptance, and shared feedbacks need one owner

OCN-15 gives the ocean and its biosphere an accepted artifact: source hashes,
support and interval identity; volume/area, heat, freshwater, salt and elemental
closure; circulation and overturning diagnostics; light, limitation,
productivity, biomass, export and redox state; restart continuity; common-
support convergence; and declared null/reduced cases.  A plausible global
productivity total cannot compensate for a leaking or resolution-dependent
ocean.

Two feedbacks cross existing ownership.  OCN-8 should contribute marine sulfur
or organic precursor and particle properties to BVOC-9's shared aerosol-number
and cloud bracket rather than designing a second cloud module.  Ocean carbon
must remain fail-closed while atmospheric CO2 is prescribed.  If that premise
is reopened, OCN-16 joins air-sea exchange and carbonate chemistry to ANUT-6,
OCN-13 and VOLC-8's weathering/outgassing bracket; only then does ocean biology
leave the pipeline's leaf tier and enter a climate convergence loop.

### 8g. The existing marine substrate already has a static latitude classifier

This is not only a future external-model problem. Orogen's `shelfClass()`
currently selects `carbonate` below 30 degrees absolute latitude and
`shelf_clastic` above it. Its own comments call latitude a proxy for sea-surface
temperature and promise downstream reclassification, but no downstream step or
task existed. Temperature alone would still be incomplete: carbonate state,
biological production and terrigenous supply are the mechanisms the class is
standing in for. LITH-26 owns a climate/ocean-state-derived reclassification on
the 10M support and makes Orogen's latitude result a labelled bootstrap only.
OCN-11 supplies support, ANUT-6 supplies sediment delivery and OCN-13/14 supply
chemistry and biology; none should independently rewrite the rock map.

## 9. cGENIE read at source: cost, constants, coupling and ceiling

Section 6 selects the EMIC tier on the coupling-class argument. This section is
the reading behind that selection: what the model costs, what it hardcodes,
what its published ExoPlaSim coupling supplies, and where its resolution stops.
It is to cGENIE what finding 3 is to LSG, and it is the evidence OCN-3 through
OCN-4 and OCN-17 through OCN-20 are drawn from.

### 9a. The cost argument selects this tier

Section 6 chose offline coupling because equilibrating a deep ocean takes
thousands of model years against orbits priced in wall-clock hours, and then
selected a host whose expense is the reason equilibration is unaffordable. The
EMIC tier dissolves that tension rather than routing around it. Three costs,
each measured by the model's own authors on Earth configurations:

| source | configuration | cost |
| --- | --- | --- |
| Ridgwell et al. (2007) | GENIE-1, 12 active ocean biogeochemical tracers, CO2-climate feedback on | better than 1000 model years per 2.4 GHz CPU hour |
| Edwards and Marsh (2005) | frictional-geostrophic ocean, EMBM atmosphere, dynamic-thermodynamic sea ice | 20 kyr of integration in about a day on a PC |
| Capirala and Olson (2026) | 36 x 36 x 16, to physical and biogeochemical steady state at 20,000 model years | 24 to 48 hours on a single CPU core |

What that buys is not only a cheaper run. Finding 4 brackets nutrient delivery
across a factor of four, and at this price that bracket is an ensemble rather
than a choice. It also makes OCN-5's loop something that can be run to its exit
predicate and measured, instead of argued about in advance of a host nobody can
afford to iterate.

### 9b. The ExoPlaSim-to-cGENIE coupling is built, published and archived

This is the part that makes the tier more than an alternative worth naming.
Liu et al. (2024) and Capirala and Olson (2026) both drive cGENIE with
ExoPlaSim, one way: ExoPlaSim v3.0.6 at T21 with 10 layers, run 100 model years
to radiative balance at about 10 hours on 8 cores, supplying annually averaged
wind fields and zonal planetary albedo, regridded to cGENIE's 36 x 36 grid.

Both halves are archived. The model is `cgenie.muffin` v0.9.50 under an MIT
licence at `10.5281/zenodo.10798347`; the regridding path is
`acapirala/exoplasim_genie_regrid`, CC-BY-4.0, concept DOI
`10.5281/zenodo.10802839` resolving to v1.2, which carries
`make_exoplasimtogenie_input` plus topography converters in both directions.

**It arrives with a knob, and the knob is the reason to read it before
adopting it.** The regrid README instructs the user to set the wind stress
scaling `ea_11` and `go_13` to 2.0 for ExoPlaSim v3.0.6 and 2.6 for v3.3.0, and
to move `bg_par_gastransfer_a` from 0.715 to 1.1 alongside. That is a
multiplier on a physical forcing whose value is selected by which version of
the atmosphere produced the field. Under this project's conventions it is a
knob and not a coupling: it must be derived, bracketed or declared, never
inherited. `vendor/exoplasim` is a fork of neither release, so neither pair of
numbers transfers, and the quantity the scaling stands in for has to be
identified before a field crosses this boundary at all.

### 9c. Finding 3's treatment, applied to cgenie.muffin v0.9.50

Finding 3 rejected LSG by reading four constraints out of its source. The same
reading here returns a different answer, and the difference is the whole reason
this tier is worth a price rather than a rejection.

Already parameterised:

- **Rotation.** `genie-goldstein/src/fortran/initialise_goldstein.F:381-385`
  branches: where the solar and sidereal day lengths differ, `fsc = 4*pi/sidaylen`,
  and otherwise it falls back to the Earth literal `2*7.2921e-5`. `sidaylen` and
  `sodaylen` are namelist reals in `ini_gold_nml` (line 207) and in
  `ini_embm_nml` (`initialise_embm.F:217`). The source comment dates the change
  and names its author: "CL (01/15/24) : used sidereal day length for Coriolis
  effect scaling factor". Liu et al.'s modification is in the tree, not in a
  branch.
- **Calendar.** `yearlen` is a namelist real and `syr = yearlen * sodaylen`
  (`initialise_goldstein.F:292`), so the year in seconds follows from
  days-per-year and the length of a day rather than from 365.25. Vesper's orbit
  and its 30-hour day reach that pair directly. Two separate Earth constants are
  NOT covered by it and are compile-time: `global_daysperyear = 365.25` at
  `genie-main/genie_control.f90:190`, and `conv_yr_d = 365.25` at
  `gem_cmn.f90:585`.
- **Ocean depth and grid.** `dsc = par_dsc`, so maximum depth is configuration;
  `imax`, `jmax` and `kmax` are build-time, and muffingen writes `GENIENX`,
  `GENIENY` and `GOLDSTEINNLEVS`, so 36 x 36 x 16 is a configuration rather than
  a constant.

Not lifted:

- **Planetary radius, in four independent places.** Three bare literal
  assignments, `initialise_goldstein.F:376`, `initialise_embm.F:459` and
  `genie-goldsteinseaice/src/fortran/initialise_seaice.F:246`, each
  `rsc = 6.37e6`; and one compile-time
  `REAL,PARAMETER::const_rEarth = 6.37E+06` at `gem_cmn.f90:802`. The parameter
  is the serious one. It sets every grid-cell area and volume in BIOGEM,
  ATCHEM, SEDGEM, ROKGEM, ECOGEM and GEMLITE, and it scales the overturning
  streamfunction. At 1.20 Earth radii every area is out by 1.44, and the
  biogeochemical inventory built on those areas with them.
- **Gravity**, `gsc = 9.81` at `initialise_goldstein.F:386`, entering the
  density scale `rhosc = rh0sc*fsc*usc*rsc/gsc/dsc`.

So the count is four Earth constants across five lines, two of which upstream
has already turned into namelist inputs, against LSG's four constraints of
which none is liftable short of a fork the size of the LPJ-GUESS port. That is
a different order of problem. It is NOT a finding that changing the remaining
two is correct: `rsc` and `gsc` enter derived scale factors, and whether the
frictional-geostrophic closure and its fitted transport parameters still mean
anything at 1.20 radii is exactly what OCN-12 has to decide rather than assume.

### 9d. The Earth content that no constant sweep reaches

Both descriptive papers state plainly that the model's parameters were fitted to
Earth observations, and this is the exposure that matters more than the
constants.

Edwards and Marsh (2005) present the transport parameters as "a first attempt at
tuning a 3-D climate model by a strictly defined procedure": a 1,000-member
ensemble scored against observed meridional overturning and Atlantic heat
transport, with the warning that "single-parameter sensitivity studies can
therefore be misleading". Ridgwell et al. (2007) calibrated the ocean carbon
cycle by assimilating three-dimensional observed phosphate and alkalinity with
an ensemble Kalman filter, and offer a global export production of 8.9 PgC/yr
and CaCO3 export of 1.2 PgC/yr as evidence of the fit.

An EMIC's skill is substantially IN its calibration, which is the opposite of a
primitive-equation model whose planetary parameters are runtime values and whose
closures are argued rather than fitted. Section 6 dropped that alternative on
cost, so this exposure is not a comparison any more: it is simply the cost of the
remaining candidate, worse in kind and not merely in count, and none of it appears
in a grep for constants. OCN-12 owns it, and this is the case it has to be
sharpest about.

One measured point cuts the other way, and it belongs here because it is the
obvious objection to the tier. Edwards and Marsh find that model errors "are
reduced only moderately by a doubling of resolution". Section 8b's concern is
still real, because straits, sills and partial coasts are not representable at
10 degrees of longitude whatever the error statistics do. But the presumption
that a coarse frictional-geostrophic ocean is wrong in proportion to its
coarseness is not supported by its own authors' measurement, and OCN-11 should
test that rather than assume it.

### 9e. What this tier buys that the resolved tier does not

**Finding 5d, and OCN-16.** This audit's own pricing makes ocean carbon the only
pathway on which a modelled marine ecosystem is a first-order control of this
world's climate, and parks it behind prescribed pCO2. BIOGEM, SEDGEM and ROKGEM
(Colbourn et al., 2013) are a closed carbon cycle with carbonate chemistry,
sediment burial and terrestrial weathering. A circulation host with an ecosystem
package bolted on supplies the ocean half and neither the sediment nor the
weathering side, which is one more reason the resolved tier was not the answer to
this. If OCN-16 is ever opened, this is the tier that answers it, and
`pedology/scripts/weathering_fluxes.py` already covers part of what ROKGEM does.

**Trait-based ecology at EMIC cost.** ECOGEM (Ward et al., 2018) resolves an
arbitrary number of plankton populations with traits assigned by size and
functional group at runtime; its reference configuration is 16 populations
across eight size classes in two functional types, and 1.1 adds a diatom group
(Naidoo-Bagwell et al., 2024). The property section 6 used to prefer a
trait-based tier over MARBL, that community composition emerges from a declared
trait space rather than arriving as fixed types, is a property of ECOGEM, so
dropping the resolved host costs nothing on this axis. MARBL remains the fallback
precisely because it is driver-agnostic.

### 9f. What it costs that the resolved tier does not

**A second atmosphere.** cGENIE ships EMBM, and the published coupling
prescribes wind stress, wind speed and planetary albedo into it rather than
replacing it. The modelled land and the modelled ocean would then be reading
different atmospheres. Section 8a's one-owner rule does not forbid that, but it
does forbid leaving it unstated: OCN-10 must name which heat, water and
momentum terms EMBM owns and which arrive prescribed.

**The forcing trap is untouched.** Both published couplings are one way, and
this project wants the transport back as surface code 903. Section 7 and OCN-5
apply unchanged, and if anything the trap is sharper here, because a model this
cheap makes it easy to run the loop before its exit predicate has been declared.

**MATLAB.** muffingen and the ExoPlaSim regrid script are both MATLAB. Whether
either runs under Octave is unverified, and it is a host question before it is a
modelling one.

### 9g. 36 x 36 x 16 is a convention, and the connector is not what caps it

The published coupling regrids onto 36 x 36 with 16 levels, which reads like a
design point of the connector. It is not. Across the 547 configurations shipped
with v0.9.50:

| resolution | configs |
| --- | ---: |
| 36 x 36 x 16 | 376 |
| 18 x 18 x 16 | 65 |
| 36 x 36 x 17 | 58 |
| 36 x 36 x 8 | 22 |
| 18 x 18 x 8 | 15 |
| 12 x 12 x 8 | 4 |
| 48 x 40 x 16 | 2 |
| 36 x 36 x 32 | 1 |

A non-square higher horizontal grid and a doubled vertical grid both exist as
working configurations. `maxi`, `maxj` and `maxk` are preprocessor defines
(`ocean.cmn:32-33`) that muffingen writes into the build options
(`muffingen.m:1580-1584`), and `par_dsc` sets maximum depth from the namelist.
The grid is genuinely parametric and 36 x 36 x 16 is what 376 configurations
happen to use.

**What caps it is the ocean component, not the regridder.**

- **It is serial.** Every OpenMP directive in the tree is commented out:
  `tstepo.F:74` reads `c!$omp parallel`, `biogem.f90:718` and `:749` are
  commented the same way, and `genie.F:484` is inside a `c$$$` block. The
  single-core costs in 9a are therefore not a choice about how the runs were
  configured, and there are no cores to spend on a finer grid.
- **The barotropic solve runs every ocean timestep and scales badly.**
  `invert` factorises the streamfunction matrix once at
  `initialise_goldstein.F:1761`, but `ubarsolv` is called from `goldstein.F:429`
  each step, and it is a forward elimination and a back substitution over the
  banded factors: each is an outer loop of `imax*(jmax+1)` with an inner width
  of at most `imax+1`. Under uniform horizontal refinement by a factor r that is
  r^3 per call, and the band storage `gap(mpxi*mpxj, 2*mpxi+3)` grows as r^3 with
  it. The one-off factorisation is worse, roughly r^4.
- **The timestep is a namelist integer, not a stability result.**
  `initialise_goldstein.F:515` sets `tv = sodaylen*yearlen/(nyear*tsc)` and
  assigns it uniformly to every level; `nyear` defaults to 100. There is no
  Courant number, no stability test and no diagnostic anywhere in the component,
  and the source's own comment at `tstipo.F:24` says the variable-timestep option
  "prevents convergence and obscures instabilities". Refining the grid means
  raising `nyear` by hand, which multiplies the per-year cost by a further factor
  of about r, and forgetting to raise it produces an unstable run rather than a
  complaint.

Taken together the barotropic cost per model year goes roughly as r^4 on one
core. Doubling to 72 x 72 is therefore of order 16x, which turns 9a's 24 to 48
hours into something between two and four weeks for the same 20,000-year
spin-up. **That estimate is read off the loop structure and is NOT a
measurement**, which is exactly why OCN-3 asks for wall clock on this machine:
the tier's whole attraction is its price, and its price is a property of
36 x 36 rather than of the model.

Two further ceilings are worth knowing before anyone tries.
`GOLDSTEINMAXISLES` defaults to 10 (`ocean.cmn:24-25`) and caps the island
count at compile time, while a finer coastline resolves more islands, so
muffingen's `.paths` and `.psiles` generation is the part of the connector that
genuinely does need checking at higher resolution. And `dt` is derived from
`sodaylen` and `yearlen`, so adopting this world's calendar changes the absolute
ocean timestep at fixed `nyear` before any grid change is considered, with
nothing in the component checking the result.

### 9h. There is no disabled parallelism to re-enable, and threads do not touch the solver

The commented OpenMP in 9g reads like working parallelism somebody switched off.
It is not. Both blocks are abandoned plumbing experiments.

In `tstepo.F:74-90` the directives are `!$OMP SECTIONS` around two calls to
`tstepo_flux_t` operating on duplicate copies of the tracer state, beside the one
real `call tstepo_flux()`. In `biogem.f90:718-735` the sections call
`sub_wasteCPUcycles1` and `sub_wasteCPUcycles2`, which are real routines in
`biogem_lib.f90:2261-2299`. That is somebody testing whether OpenMP could be
plumbed through the build at all, not a parallel decomposition of the physics.
**Making the component threaded is new work, not an un-commenting**, and the
absence of any partitioning of the tracer arrays is the reason.

One free result falls out of reading it. `tstepo.F:32-35` declares `ts_t1`,
`ts1_t1`, `rho_t1`, `ts_t2`, `ts1_t2` and `rho_t2`, and lines 62-68 fill all six
on every ocean timestep. Nothing reads them: the only consumers are the
commented calls. That is six full array copies per timestep serving dead code,
and it grows with resolution exactly as the rest does. Deleting it is a
correctness-neutral saving available before any parallelisation question is
settled.

**Threads and the solver are complementary, and neither substitutes for the
other.**

- `ubarsolv` will not thread. Its forward elimination
  `gb(j) = gb(j) - ratm(j,j-i)*gb(i)` and its back substitution both carry a
  loop-carried dependence; a banded triangular solve is sequential along the
  band by construction. Threading the ocean does not make the barotropic step
  faster.
- What threads well is the other half. BIOGEM's `do n=1,n_vocn` loop over the
  vectorised ocean columns is close to embarrassingly parallel, and the tracer
  loops in `tstepo_flux` are 3D sweeps. That is where the tracer count buys
  work worth dividing.
- What the barotropic solve wants is an algorithmic change. The band structure
  is a consequence of the lexicographic ordering `k = i + j*n` at
  `invert.f:32`, which fixes the bandwidth at `imax`. A fill-reducing ordering,
  or an iterative treatment of what is after all an elliptic streamfunction
  equation, attacks the r^3 directly rather than dividing it. The complication
  to price first is the island machinery: `ratm`, `psisl`, `erisl` and
  `matinv_gold` all assume the current factorisation, and `GOLDSTEINMAXISLES`
  bounds them.

**Which of the two pays more is unmeasured, and the two arguments are
different.** The barotropic solve dominates GROWTH with resolution; it may still
be a minority of RUNTIME at 36 x 36, where twelve or more BIOGEM tracers are
doing the bulk of the work. Those are not the same claim and this audit
establishes neither. A profile at the shipped grid decides the order, and it is
cheap next to either piece of work.

The acceptance criterion already exists and should be used rather than invented.
`genie-knowngood/` ships per-component NetCDF output for four configurations,
with a changelog, which is the right answer an optimisation has to reproduce.
Under this project's conventions that is what makes the work testable rather
than merely benchmarkable.

Patching performance into cGENIE makes it a maintained fork on the ExoPlaSim and
LPJ-GUESS pattern. That is the expected end state for an external model here
rather than an argument against one, and
`docs/src/reference/design-intent.md` records the reason: these models were not
built for this project's subject or its runtime constraints, so the question
about a candidate is whether its physics is worth the fork rather than whether it
can be adopted unmodified. What a fork does NOT do is answer the physics
question, which is why OCN-3 still precedes OCN-19 and OCN-20: performance work
is spent on a host that has been chosen, not on one that is being evaluated.

## 10. What this audit did NOT establish

Recorded so the next reader knows the edges.

- Nothing has been built on this machine. cGENIE's source and licence have been
  read at v0.9.50 and its costs are quoted from its authors' papers on Earth
  configurations, but nothing has been compiled or timed HERE, which is the only
  number OCN-3 accepts.
- MITgcm was dropped on the coupling-class argument rather than on a measurement
  of it. Nobody here has timed a coarse global MITgcm configuration or its
  spin-up, and the claim that its connectivity is priced out of reach is an
  inference from the equilibrium requirement, not a benchmark. It is a cheap
  thing to falsify if a future reader wants to: the number that would reopen it
  is model years per wall-clock day at a resolution that resolves the straits
  section 8b names.
- Section 9 does not establish that cGENIE's remaining Earth constants can be
  changed correctly. It establishes where they are and how many there are.
  `notes/audits/ocean-tier-implicit-earth.md` takes that further and stops in a
  specific place: it establishes that `rsc` and `const_rEarth` are two
  independent statements of the radius that the geochemistry and the circulation
  read separately, so a port that changes one silently disagrees with the other,
  and it extends OCN-17's declared bracket from the wind stress scaling alone to
  the seven-parameter tuned block that ships beside it. Whether the
  frictional-geostrophic closure survives 1.20 radii remains unanswered, and
  that audit adds why it is hard to answer: the model's own error function
  scores against Earth observational fields and is inert here.
- Section 9g's r^4 cost growth is an estimate from the shape of `ubarsolv`'s
  loops, not a timing. No cGENIE configuration has been compiled or run here at
  any resolution, and the two to four week figure for a 72 x 72 spin-up should be
  treated as an order of magnitude that OCN-18 exists to replace.
- Nothing here checks what muffingen produces above 36 x 36. That it CAN write
  larger dimensions is read from `muffingen.m:1580-1584`; whether the island and
  path generation stays correct there is unknown, and so is whether the tool runs
  outside MATLAB.
- Section 9h does not establish which half of the ocean component dominates
  runtime at 36 x 36. It establishes that the growth argument and the runtime
  argument are different and that nobody here has measured either. OCN-19's
  profile is the first thing that should happen if this host is selected.
- Whether a fill-reducing reordering or an iterative solve is the better answer
  is not established either, and the island machinery that constrains both has
  been located but not priced.
- The wind-stress scaling in the published ExoPlaSim coupling, 2.0 or 2.6
  depending on the atmosphere's version, is recorded as a knob that has to be
  explained. What physical quantity it stands in for has NOT been determined,
  and neither paper is read closely enough here to say whether it is documented
  anywhere but the regrid README.
- The water-leaving albedo span in 5a is DECLARED, not sourced. It carries the
  pricing's whole uncertainty and no paper in `references/INDEX.md` backs it.
- Every climate figure is on the bootstrap climatology. The 5a and 5d numbers
  move when a baseline exists, though not by enough to change 5a's verdict.
- The carve thresholds in finding 4's table are a sweep over the basins' own
  `critical_aridity_index` and are NOT a carve verdict. A verdict integrates a
  climate over catchments and is `hydrography/scripts/carve_verdict.py`'s to
  produce.
- Nothing here re-examines `nhdiff`, which is the cheapest bound on the missing
  transport and was settled separately.

---

## 11. The prescribed heat-transport channel, verified end to end

Finding 2 established that the channel exists at both ends. This section is the
verification that it carries a field correctly, run BEFORE anything computes a
real one, because a first real field arriving through an unverified channel
cannot be separated from a wrong channel carrying a right field. Section 7 makes
this the channel the adopted offline architecture depends on, so it is the one
piece of ocean work that had to come first.

The instrument is `exoplasim/scripts/verify_ocean_flux_channel.py`, registered
under `one_offs`. It stages a field and its prediction into a run directory and
then answers each criterion against the stream that run wrote. Every criterion
and every tolerance is a module constant declared in that file's source, and
they were written before any run existed.

### 11a. The field, and why its answer is known in advance

    f(phi, lambda) = A * (P2(mu) + 0.3 mu) * (1 + 0.4 cos(lambda + pi/3))

with `mu = sin(phi)` and `A` 20 W/m2. Four properties, each of which is a
prediction the channel could have failed:

- **Constant in time.** One month is written to the `.sra` and
  `surfmod.f90:get_surf_array` expands a single record to all fourteen by copy.
  `getflxco`'s month weights sum to one whatever the calendar says, so `yfsst2`
  equals the prescribed field at every ocean step and the accumulation `yfssta`
  equals it exactly. This deliberately verifies the CHANNEL and not the monthly
  interpolation: the calendar cannot influence the answer, so a failure cannot
  be attributed to it.
- **Orthogonal to the constant.** `P2 + 0.3 mu` is degree two in `mu` and the
  grid's latitudes and weights are Gauss-Legendre, which integrates degree two
  exactly, and the zonal factor has longitudinal mean one. So the area-weighted
  global integral is zero analytically. A prescribed advection that does not
  integrate to zero is a heat source.
- **Asymmetric in both coordinates.** The `0.3 mu` term breaks the equatorial
  symmetry `P2` alone would have and the phase breaks the zonal one, so a
  north-south flip, a longitudinal roll, a sign error and a transpose are all
  DIFFERENT fields. The comparison names the failure instead of reporting a
  difference.
- **Nonzero over land.** `getflxco` interpolates over the whole grid and only
  `addfc` restricts to `yls < 1`, so the diagnostic must carry the field over
  land as well.

### 11b. What was run

Measured 2026-08-24. A diagnostic segment at T42 on the current model source,
built at eight threads with `build_model.py --no-publish`, in a scratch
directory outside the run registry: cold start, `nfluko = 1`, 512 model steps at
a 45-minute step, `nout` at its default of 32, giving 16 ocean-stream records of
86,400 s each. `nhdiff` 0, `nlsg` 0, `NLEV_OCE` 1, so codes 904, 905 and 906
were required to be identically zero and were. The companion sea surface
temperature climatology at code 169, which `oceanini` refuses `nfluko` without,
was uniform at 290 K so that the ocean is ice-free everywhere and the delivery
identity is testable across the whole of it.

The segment is not in `exoplasim/runs/INDEX.json` and its output is not a
climatology input. It is instrumentation.

### 11c. What the criteria returned

| criterion | result | bar |
| --- | --- | --- |
| reported field equals written field, every record, every cell | 1.90e-6 W/m2 | 1.0e-4 |
| area-weighted global integral, written field | 1.98e-8 W/m2 | 1.0e-6 |
| area-weighted global integral, worst record | 2.05e-8 W/m2 | 1.0e-4 |
| reported field over land equals written field | 1.90e-6 W/m2, on 4116 land cells | 1.0e-4 |
| delivery: scatter about one fitted proportionality | 6.49e-5 of the largest change | 1.0e-3 |
| delivery: fitted record interval against the declared one | 3.17e-7 | 1.0e-3 |

The reproduction figure is the float32 resolution of the service stream on a
20 W/m2 field and nothing else. The field's integral before it is written is
-2.75e-14 W/m2, which is the quadrature; 1.98e-8 is what the `.sra`'s
five-decimal write leaves of it, and the bar for that is derived from the
quantisation rather than chosen.

**The orientation check discriminates by five to seven orders of magnitude.**
The nearest wrong field differs from the reported one by 0.51 W/m2 (a
one-cell eastward roll), then 14.68 (a quarter turn), 16.79 (a north-south
flip), 20.76 (a half turn) and 72.67 (a sign reversal), against 1.90e-6 for the
field as written.

**The delivery test is the sharper half.** On an ice-free ocean cell `mksst` and
`addfc` between them give
`CRHOS * CPS * mld * (SST_k - SST_{k-1}) = (yheata_k + yfssta_k) * dt_record`,
so across every such cell and every record the temperature change must be ONE
constant times the total flux. Fitting that constant over 4076 cells and 15
record pairs, 61,140 samples, leaves a maximum residual of 3.05e-5 K against a
largest change of 0.469 K, and returns a record interval of 86,400.027 s against
a declared 86,400. That recovers the model's own output interval to 0.03 s in a
day from the sea surface temperature response alone, which is the strongest
statement available that the field is not merely reported but delivered.

### 11d. What the verification caught

**The fitted record interval missed by 4.75 per cent on the first pass, and the
model was right.** `close_ocean_energy.py` and `close_state_energy.py` both
carried `CRHOS = 1030.0` and `CPS = 4180.0` as literals attributed to
`oceanmod.f90`. `oceanmod` does not own them: `icemod.f90` declares both as
`icemod_nl` keys and passes them to `oceanini` so the two modules cannot hold
different sea water, and its `CPS` is sea water's specific heat at S = 34.7 and
its freezing point, 3990.34 J/kg/K, not fresh water's 4180 at about 25 C. Every
slab heat capacity built from the copied pair was too large by their ratio.

Nothing failed, because there was nothing for it to fail against.
`assess_convergence.py` imports those two names to build `SLAB_HEAT_CAPACITY`,
which sets the expected relaxation time, which sets the remaining offset its
convergence verdict is judged on. It had already replaced an uncited 3990 with
the copied 4180 on the argument that the model's own value is the right one, so
the correction moved the number AWAY from the model. That is the whole failure
mode: a residual quietly rescaled, in a quantity that decides a pass.

`lib/sea_water.py` now reads all four of the salinity quartet from
`icemod.f90`'s declaration and then from the run's own `icemod_namelist`, and
the three scripts read it from there.

**Two operational facts, recorded because a hand-run segment is how this
verification is repeated.** The model needs a large stack: at T42 on eight
threads it terminates without a Fortran backtrace under the default limit, and
runs under `OMP_STACKSIZE` and an unlimited stack rlimit, which is what
`run_exoplasim.py` sets and a direct invocation does not. And a run directory
from an earlier segment is not a namelist source: `NGUIDBG` in `plasim_nl` and
`ZETA` in `carbonmod_nl` no longer exist in the model, and each is a hard
runtime error on read. Neither is written by anything current.

### 11e. What this changes, and what it does not

`close_ocean_energy.py`'s D4 and D6 now carry the flux correction:
D4 is `CRHOS*CPS*mld*d(SST)/dt = yheat + yfsst` and D6 is
`hfns + yfsst = CRHOS*CPS*mld*d(SST)/dt`. On every run made so far `yfsst` is
identically zero, so both evaluate to what they did; the term was removed from
the list of things the closure requires to be zero, because refusing a run that
carries a flux correction would make the closure unusable on exactly the
configuration this work exists to produce.

**A globally zero-integral field does not have a zero integral over the ocean,
and the ocean is what `addfc` applies it to.** On this build's land mask the
field above integrates to -0.919 W/m2 over ocean cells while integrating to zero
over the globe. So the forcing contract OCN-10 owns has to state which support
the zero-integral condition holds on, and the answer is the ocean: a field
constructed to have zero global integral and then applied only to water is a net
heat source or sink of the ocean-area-weighted mean of its land part. This is a
property of the geography and it changes with every build and every carve.

**What this section does NOT establish.**

- The monthly interpolation is deliberately not exercised. The field is constant
  in time so that the calendar cannot affect the answer. A time-varying field is
  a separate verification and the identity available for it is weaker: whatever
  the calendar, `getflxco` returns a convex combination of two months, so the
  interpolated value must lie inside the monthly envelope.
- The ice branches of `addfc` are not exercised. The segment is ice-free by
  construction, so branches (b) and (c) and the residual path into `mkiflx` are
  untested, and the delivery identity was evaluated only on cells with no ice in
  any record.
- Nothing here says a heat transport field is right. It says the channel carries
  the field it is given, without loss, without a mask, without an orientation
  error, and delivers it to the slab at the declared heat capacity and interval.
  The field itself is finding 6's problem and OCN-5's loop.
- The segment is at one rung. The instrument takes `--rung` and the argument is
  rung-independent, but only T42 has been run.

---

## Tasks

OCN-1 through OCN-20. Findings 1, 2, 5a, 5b and 5c are the
exploration rows; finding 3 is recorded here and became no row, being a
constraint rather than work; finding 4's config half is OCN-9 and its prediction
half is this document. Section 8 supplies the cross-component contracts and
acceptance rows OCN-10 through OCN-16, plus LITH-26 for the existing static
latitude shelf classifier. Section 9 supplies OCN-3's single-host viability and cost row and
OCN-4's ecosystem tier, its coupling half is OCN-17, its resolution ceiling is
OCN-18, and its profile-then-optimise split is OCN-19 and OCN-20.
