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
load-bearing; and one is a pricing that settles the coupling class. The tasks
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
already postprocesses the `ocean_output` stream and lists code 903 among the
terms it reads, recording that it is "identically zero in this configuration,
and checked to be". So the input channel and the instrument that verifies it
both exist, and both are currently null.

What this means for cost. An ocean heat transport reaches the model as ONE new
pipeline step writing one `.sra`, plus one namelist key on the
`KEY@namelist` mechanism `run_exoplasim.py` already has. What it does NOT supply
is the field itself, which is finding 6's problem.

`analysis/error_budget.json` records no-q-flux as STRUCTURAL with "no cheap
version exists". The cheap version of the BOUND is `nhdiff`, which
`notes/audits/absent-and-inherited-physics.md` finding 1 established. The cheap
version of the TERM is this, and it was not known when that finding was written.

## 3. LSG is vendored, and it is not the switch it looks like

`vendor/exoplasim/exoplasim/lsg/src/lsgmod.f90` is 9,313 lines of a real
three-dimensional geostrophic ocean. `oceanmod_nl` already carries `nlsg` and
`naomod`, the atmosphere-to-ocean step ratio, and `most.c` knows how to build it.
It looks like the answer. Four things in the source say otherwise.

- `most.c:2415`, comment and code together: "LSG works currently only with T21
  PlaSim", and it FORCES `Resolution = RES_T21`. Production here is T42 with T85
  in loop D.
- `plasim/src/cpl.f90` hardcodes the coupler grids as `nxa=64, nya=32` on the
  atmosphere side and `nxo=72, nyot=76` on the ocean side, and the planet as
  `parameter(radea=6.371E6)`. This world's radius is 1.20 Earth.
- The LSG grid is a fixed 72 x 76 x 22 present-day Earth configuration, and
  `lsg/dat/` ships Earth bathymetry, an Earth runoff map and Earth restart
  fields to go with it.
- `oceanmod.f90:334` is a hard error: LSG coupling requires
  `n_days_per_year = 360`. This world's orbit is 182.8 Earth days at the
  current flux, and its day is 30 hours.

**Adopting LSG means forking it to a depth comparable to the LPJ-GUESS port and
then running the atmosphere at T21.** Recorded here so that the next reader who
finds `lsgmod.f90` does not re-discover the four constraints.

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

**Coupling class: offline, on the aeolian pattern.** The ocean model reads a
climatology and writes fields the next climate run is forced with. It is not
synchronously coupled to the atmosphere. This is the only class that is
affordable, because equilibrating a deep ocean takes thousands of model years
and an ExoPlaSim orbit is priced in wall-clock hours
(`docs/src/pipeline/costs.md`).

**Circulation host: MITgcm.** PISCES, MEDUSA and PlankTOM are NEMO
configurations, so adopting any of them means adopting NEMO, whose grid and
configuration machinery are harder to move off Earth than MITgcm's. MITgcm takes
the sphere radius, gravity and the rotation rate as runtime parameters, which is
exactly the axis this world moves along.

**Ecosystem tier: Darwin.** MARBL is the better-engineered library by some
distance, being deliberately driver-agnostic with a documented coupling
interface, and it is the right fallback if the host question turns out
differently. Its ecology is a small set of fixed Earth plant functional types,
two of which are specific Earth evolutionary inventions. That is the same
problem `biosphere/README.md` already records for LPJ-GUESS: "the shipped plant
functional types are Earth's ... should be declared that way rather than
presented as a prediction". Darwin is trait-based, so the modelled community
composition is an output of the trait space rather than an input, which on a
world whose nutrient supply is bracketed across a factor of four (finding 4) is
the difference between a prediction and an assumption.

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
An MITgcm grid therefore cannot be chosen once and called the planet.  It needs
a versioned map from the 10M bathymetry and coastline into wet fractions,
depth/volume, shelves, straits, sills and routed river mouths, plus conservative
maps to every accepted T21/T42/T85/T127/T170 atmosphere.  Narrow connections
and partial coastal cells can change circulation even when cell-mean depth is
unchanged.  OCN-11 extends SPAT-1 through SPAT-5 and SPAT-8/10 to this crossing;
OCN-3's coarse configuration is a candidate whose adequacy must be measured,
not a fixed resolution.

### 8c. MITgcm and Darwin do not escape the implicit-Earth audit

Runtime radius, gravity and rotation parameters are necessary but not
sufficient.  Calendar and rate constants, reference pressure, equation of
state, geothermal and tidal assumptions, vertical and lateral closures,
salinity and initial hydrography can still import Earth.  In an ecosystem,
fixed Redfield relations, nutrient inventories, temperature functions,
sinking speeds, irradiance/PAR conversions, latitude classifiers and supplied
community traits do the same.  Latitude remains valid for geometry and
Coriolis; it is not a water-mass, productivity or ecological regime.  OCN-12
applies BIO-22/BIO-29's clock and latitude rules and the repository's
gravity/pressure audit shape before either external model is accepted.

### 8d. “Reached the ocean” is a boundary condition, not marine nutrition

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

Darwin makes community composition emergent only inside the trait space it is
given.  OCN-14 must declare that space with covarying strategies and retain a
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

## 9. What this audit did NOT establish

Recorded so the next reader knows the edges.

- MITgcm has not been built on this host, and neither its build nor Darwin's
  licence, tracer cost or trait-space requirements have been checked against
  anything. Section 6 is a recommendation from the coupling and Earth-content
  arguments, not from a measurement.
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

## Tasks

OCN-1 through OCN-16 in `TASKS.md`. Findings 1, 2, 5a, 5b and 5c are the
exploration rows; finding 3 is recorded here and became no row, being a
constraint rather than work; finding 4's config half is OCN-9 and its prediction
half is this document. Section 8 supplies the cross-component contracts and
acceptance rows OCN-10 through OCN-16, plus LITH-26 for the existing static
latitude shelf classifier.
