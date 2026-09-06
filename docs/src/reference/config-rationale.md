# Why config/planet.yaml says what it says

The arguments behind `config/planet.yaml`, keyed by the setting each belongs
to; the file says what each value is and points here. Where a block disagrees
with `exoplasim/notes/parameter-decisions.md` or an audit, this file is the
current rationale and governs; the note is the dated derivation and stays as
written.

## The statuses

Every setting carries one, and it records what STANDS BEHIND the value rather
than how old it is. Four atomic statuses are defined at the top of
`config/planet.yaml`, and the two lists have to agree. Compound labels build
on them where a block is not one thing: MIXED and PARTLY DETERMINED mark
blocks whose sub-keys carry different atomic statuses, TRANSITIVE marks
operational values that are facts about the run rather than the world, and
DECLARED ABSENT is DECLARED applied to a deliberate null.

DETERMINED, measured or derived by an artifact that exists, or decided in a way
that is now load-bearing; downstream work may rely on it. DECLARED, a stated
position on something this project cannot compute or has not computed
reproducibly, written down so it is not read back as a result. PROVISIONAL, a
placeholder or an interim value that WILL move. DERIVED, computed by a named
script from other values in the file, so regenerated rather than edited.

The distinction is not cosmetic and it decays in one direction only:
`stellar_cycle` has carried an amplitude and a centre since before there was a
baseline flux to centre on, and both were read back later as though they had
been chosen. They had not. A value that has never been decided should not be
indistinguishable from one that has -- and a value whose evidence has since been
deleted should not be either, which is why a status is demoted in the same pass
that removes what it rested on.

## `source_build`

DETERMINED. Which World Orogen build every component reads. source/ is namespaced by build
because they now multiply faster than they can be swapped in place, and a
result's provenance should depend on what it was computed from rather than on
when it was computed.

lib/orogen.py is the registry and holds the hash, the note and what superseded
what. No list of builds belongs here: it goes stale the moment a build lands,
and a stale list beside the live setting reads as though it were checked.

## `baseline_climatology`

Which climatology every downstream component reads by default. Named here
rather than hardcoded in a _paths helper, because it was hardcoded once -- to a
climatology computed on pre-carve terrain under the wrong stellar spectrum --
and six scripts across pedology and biosphere took that superseded default
silently.

The key is either a path or null, and null is a real state rather than an
omission: a carve moves `source_build`, every climatology on the terrain it
replaces becomes a property of a world this project no longer holds, and the
key goes back to null until a commissioning produces one. `climatology_path()`
raises on null rather than falling back, which is the point -- a stale default
returns a plausible number from the wrong world, and that is worse than an
error. Read the current state from `config/planet.yaml`; the comment above the
key there is where the reason for the current one lives.

It takes the REGULAR climatology of the BASELINE run and never the bootstrap's.
The bootstrap is the first run on a terrain and exists so that lakes, the lake
compositing in the albedo and the pedology soil water field can be built at
all; the run whose climatology the carve verdict uses is the one made
afterwards, with those fields in place. A step that will take either reads
`best_available_climatology`, which names the stage it chose.

A named climatology is not automatically the ACTIVE build's climate.
check_consistency.py compares terrain hashes and reports the mismatch; that is
where to find out, not a comment here that has to be kept true by hand.

## `planet`

```
planet:
```


DETERMINED. Radius and gravity are Orogen's and canonical; mass follows from
gravity, not the other way round. Rotation, obliquity and eccentricity chosen.

## `gravity_m_s2`

```
gravity_m_s2: 12.81
```

Surface gravity is canonical, not derived. This is the value World Orogen
was run with, and it already scaled the terrain's maximum relief as 1/g, so
the geography in source/ cannot be separated from it. Mass follows from it:
mass_earth = (g/g_earth) * radius_earth^2, and `derive()` raises if the two
keys below disagree, so they cannot drift apart. The formula is written here
and the value is not.

## `rotation_direction`

```
rotation_direction: prograde
```

The IAU cartographic convention makes rotation direction and longitude sense
one declaration rather than two: `W`, the angle from the ascending node to the
prime meridian, increases with time for a prograde rotator and decreases for a
retrograde one, and planetographic longitude is positive west for a prograde
rotator and positive east for a retrograde one. So a world that has not said
which way it turns cannot say which way its longitudes run either, and the
three keys here are written together for that reason.

Prograde is read off the code rather than chosen, and it is structural in three
independent places.

`radmod.f90:solang` builds the solar hour angle as a term rising with the time
of day plus a term rising with the grid column, so the column at local noon
falls as the day advances: the subsolar point tracks toward decreasing
longitude and the planet turns toward increasing longitude. That is eastward
rotation on an east-positive axis.

`plasimmod.f90` initialises `plavor`, the planetary vorticity coefficient, to a
positive constant, so the Coriolis parameter is positive in the northern
hemisphere and the simulated flow is Earth-handed: cyclonic circulation turns
counterclockwise north of the equator.

ExoPlaSim reaches the model as `ROTSPD = 1/rotationperiod`, and
`run_exoplasim.py` derives the rotation period from `rotation_hours`, which is
signless and positive. A retrograde world is therefore not expressible on this
path at all: it would need a negative rotation period, which the day-count
arithmetic beside the namelist edit would then mangle. Every run this project
has integrated, archived included, was prograde, and none of them could have
been anything else.

## `longitude_positive`

```
longitude_positive: east
```

`maps/projections.py:unit` places a point at
`(cos lat sin lon, sin lat, cos lat cos lon)` and the World Orogen export
recovers longitude as `atan2(x, z)` about the same y-up pole, so increasing
longitude is a positive right-handed rotation about the north pole. East is
therefore both the direction of increasing longitude and the direction the
planet turns, and the two agree because `rotation_direction` is prograde.

Positive east on a prograde rotator is the IAU planetocentric system, which is
the one this project wants: planetocentric latitude is what a sphere has, and
the planetographic pairing for a prograde rotator would be positive west. The
mixed pairing, planetocentric latitude with planetographic longitude sense, is
explicitly not IAU-approved.

This declaration changes no join and is not a datum for one. Longitude labels
differ between the export and ExoPlaSim's own output by construction, and the
remedy is CLAUDE.md rule 3: map by index, share one coordinate source. Declaring
a longitude sense makes an axis legible; it does not make a longitude match
safe, and nothing here weakens the rule that forbids the match.

## `prime_meridian`

```
prime_meridian: orogen_export_zero
```

This world carries no fixed observable surface feature that a meridian could be
anchored to, which is the case the IAU report provides for directly: where there
is none, the expression for `W` defines the prime meridian rather than being
corrected to fit a landmark. What is available is a phase, and the phase is
already fixed.

The World Orogen mesh frame supplies it. Every mesh region carries unit-sphere
cartesian coordinates in the generator's own y-up frame, and the export's
longitude is `atan2(x, z)`, so the meridian through `x = 0, z > 0` is the frame's
zero. `maps/render_projections.py:draw_graticule` already draws that meridian and
the antimeridian heavier, in a pass over multiples of 180 degrees, and it draws
them on the export's own longitudes rather than on any model output's labels.

Naming the export frame is better than writing a `W` expression, and better than
picking a landform. It needs no celestial reference frame, which this project
does not have; it is reproducible from the planet code and the seed, both of
which are in the durable set of CLAUDE.md rule 7, so the meridian survives a
regeneration that the terrain itself does not; and it is what every artifact in
the tree was already drawn and integrated against, so the declaration costs
nothing downstream.

The value names a datum rather than an offset, deliberately. An offset in degrees
would read as a knob, and nothing in the pipeline applies one; moving the prime
meridian would mean adopting a different datum and reworking everything that
reads the export frame. The Working Group's standing emphasis is the reason to
write it down now: any change to this value changes the longitude of every point
on the planet, so a longitude system is chosen once and not revisited.

## `star`

```
star:
```


DETERMINED. Spectral type, activity level, and the surface ultraviolet that
follows from them are all decided; see the notes inside the block.

## `spectral_type`

```
spectral_type: K2.5V
```

Provisional empirical midpoint between K2V and K3V in the
Pecaut--Mamajek dwarf sequence (see exoplasim/notes/parameter-decisions.md).

## `activity`

```
activity: active
```

DETERMINED: this star is chromospherically ACTIVE, of the kind epsilon
Eridani represents -- young, spotted, and with ultraviolet well above a
quiescent K dwarf's. Chosen to match the intense starspot cycle this world is
built around, and it settles a decision three other quantities were waiting
on: the ozone column, the surface ultraviolet environment and the Lacis-
Hansen ultraviolet weight are all functions of activity, so all three now
take their active-star values from Segura et al. (2003), whose K2V IS
epsilon Eridani observed by IUE.

## `surface_uv_relative_to_earth`

```
surface_uv_relative_to_earth: 0.4
```

DETERMINED, and now measured rather than chosen. Segura et al. (2003),
"Ozone Concentrations and Ultraviolet Fluxes on Earth-Like Planets Around
Other Stars", Astrobiology 3, 689-708, models a K2V host explicitly -- very
nearly this star -- with an Earth-like O2 atmosphere:

ozone column at 1 PAL O2   6.64e18 cm-2 against the Sun's 8.36e18, so 0.794x
surface UV-B               "about 0.4 times Earth's flux"

So this world is substantially BETTER protected from ultraviolet than Earth,
not comparably. The star's lower ultraviolet output wins against the thinner
ozone column it produces; the two do not cancel.



This closes a dimension instead of opening one. Surface UV is the product of
how much the star emits and how much the ozone column absorbs, and neither is
known here: BT-Settl is a photospheric model and carries no chromosphere,
which is where a K dwarf's 200-320 nm flux mostly comes from, and the ozone
column is prescribed as Earth's rather than computed. Pinning the product
means neither factor has to be determined separately.

It is defensible rather than merely convenient because ozone is itself
UV-produced, so the column tracks the incident ultraviolet and buffers the
surface against the star's activity. That is what lets this world have a
dramatic bolometric activity cycle and a terrestrial surface UV environment
at the same time -- the two would otherwise pull against each other, since
the cycle wants an active star and the biosphere wants a quiet one.

What this does NOT settle is ozone's radiative effect on the climate, which
is a different quantity from its shielding of the surface; pricing it is CLIM-32.

## `metallicity`

```
metallicity: 0.0
```

DECLARED 2026-08-18. Solar, which is unremarkable for a K2.5V star; it earns
an entry as the star's most powerful otherwise-undeclared parameter.

The leverage is what earns it a line. Moving half a dex along the BT-Settl grid
moves the shortwave band-1 share by +0.0114 at [M/H] -0.5 or -0.0084 at +0.5,
which is five times what the resampler bias fixed alongside it was worth and
twenty times what rounding the surface gravity to the grid costs. That share
weights every two-band snow, sea ice, glacier and ground albedo, so it reaches
the surface energy balance directly. Measured in `notes/audits/stellar-spectrum-oracle.md`.

Surface gravity is deliberately NOT declared here, and the contrast is the
point. It is not a free parameter: it follows from `mass_solar`,
`luminosity_solar` and `effective_temperature_k` through R/Rsun = sqrt(L)
(Tsun/Teff)^2, and `build_stellar_spectrum.py:derive_log_g` computes it, snaps
it to the 0.5 dex grid and refuses to build if the snap lands anywhere but the
grid point its pinned SVO record ids serve. Declaring it would create a fourth
copy of a quantity the other three already fix, free to drift away from them
silently. Metallicity has no such derivation, which is exactly why it has to be
stated.

## `orbit`

```
orbit:
```


MIXED. `earth_solar_constant_w_m2` is a physical constant and
`sweep_flux_earth` is a sweep specification; both are DETERMINED.
`baseline_flux_earth` is CHOSEN, and `longitude_vernal_equinox_degrees` is
DECLARED. Each is argued below.

## `baseline_flux_earth`

```
baseline_flux_earth: 0.945
```

CHOSEN. The habitability-by-latitude derivation of
`docs/src/pipeline/state.md` section 5b is codified as `derive_design_flux.py`
(step `design_flux`), which declares its thresholds in advance. It has been RUN
on this terrain and it REFUSES: four land bands have a negative warmest-bin
response between the two measured flux points, so the per-band projection has
nothing to scale by. That refusal is not one another run lifts, because the
sign is a property of how this world's seasons respond rather than of where the
two points sit, so a design decision is what settles the number.

The record is `exoplasim/analysis/design_flux.json` with `basis: chosen`,
written by the same script and only after it caught its own refusal. It carries
the refusal with the offending bands, the finding the choice rests on in
`notes/audits/design-flux-two-point-response.md`, and what would reopen the
derivation: a projection method that can represent a band whose seasonal range
contracts with warming. `scripts/check_consistency.py` refuses a chosen record
short of any of those, and `docs/src/pipeline/sequencing.md` loop C re-derives
the flux on every new terrain, so the choice is re-affirmed per terrain rather
than inherited.

It is an INPUT and not a result: it fixes the semi-major axis below, and that
orbit is compiled into LPJ-GUESS. What the flux and the biosphere share is a
FINDING rather than an artifact, and that is the part the decision rests on:

The flux and the biosphere are one choice, not two. The vegetated and
bare-rock windows for the design band DO NOT OVERLAP, so this flux is
habitable *because* the world is vegetated rather than independently of it.
That non-overlap is the load-bearing fact and is why the number is here.

This fixes the orbit. a = sqrt(L/F) = 0.585173 AU and the year is 182.801
days, and that year is compiled into LPJ-GUESS. Subsequent flux changes
should therefore move the star's luminosity and leave the orbit alone.
L scales 1:1 with F at fixed a, and Teff as F^(1/4), so past 2 to 3% in
luminosity the spectrum needs rebuilding at the new temperature; see
`docs/src/pipeline/sequencing.md` loop C.

## `longitude_vernal_equinox_degrees`

```
longitude_vernal_equinox_degrees: 102.7
```

DECLARED, 2026-08-19, after `notes/audits/inherited-earth-constants.md`
finding 4: the number is ExoPlaSim's Earth default, and nothing about this
world determines it.

What it sets is the phase of perihelion against the equinoxes and solstices,
in ExoPlaSim's convention the true longitude of the vernal equinox measured
from perihelion. Nothing about this world determines that phase: Orogen has no
time axis, no precession is modelled, and the phase is a snapshot choice of
exactly the kind pCO2 is. So it cannot be DETERMINED, only declared.

What the choice is worth. At e = 0.02 the perihelion-to-aphelion flux ratio is
((1+e)/(1-e))^2 = 1.083, 8.3% peak to peak over every orbit -- larger than the
stellar cycle's total declared amplitude of 6.0%, and unlike the cycle it is
hemispherically antisymmetric: whichever hemisphere's summer falls near
perihelion gets brighter summers and darker winters, the other the reverse. At
102.7 degrees the geometry is Earth-like, perihelion falling near the southern
summer solstice, so the south is the amplified hemisphere. That lands on the
season `docs/src/pipeline/state.md` section 5b establishes as the operative one for the ice line, which is
why the value deserves a name rather than a default.

Why it is KEPT at 102.7 rather than moved: the phase is already compiled into
LPJ-GUESS through the solstice offset in `vesper.h`, and every existing run and
climatology carries it, so any alternative costs a biosphere rebuild, a driver
regeneration and a re-run for a quantity nothing constrains. A symmetric
alternative exists if one is ever wanted -- perihelion at an equinox (0 or 180)
makes the two hemispheres' seasons equal in amplitude -- and choosing it would
be a worldbuilding decision about the simulated seasons, not a correction.

## `atmosphere`

```
atmosphere:
```


DETERMINED for the composition, DECLARED for the CO2. The partial pressures sum
to exactly 1 bar by construction; 450 ppm is chosen -- and was CHECKED,
2026-08-16, though the check no longer has an artifact. See the last paragraph.

450 ppm is prescribed. Nothing in this project solves the carbonate-silicate
balance that would set it, and closing that loop is deliberately not on the
roadmap: it would replace a prescribed CO2 with a prescribed OUTGASSING rate,
which is no better constrained, at the cost of a weathering law in pCO2 and
several coupled climate integrations. One free parameter would move, not go.

What was done instead is a check, and 450 ppm passes it with room.

At steady state a planet's outgassing equals its own silicate weathering, so
the ratio of this world's silicate CO2 consumption to Earth's IS the ratio of
outgassing it needs. That comparison avoids absolute outgassing estimates,
which span a factor of several depending on whether metamorphic and
diagenetic sources are counted. `pedology/scripts/weathering_fluxes.py`
computes it.

Holding 450 ppm at this climate asks for about 1.07x Earth's outgassing.
Radiogenic heat production scales with mass at fixed composition and this
planet is 1.881 Earth masses, so the first-order expectation is near 1.9x.
The requirement sits below the supply with a margin around 1.7x.

Read it as a FLOOR, not a target: supplying more breaks nothing, it settles
the thermostat at a higher CO2 and a warmer state, which is a different world
rather than an inconsistent one. What the check rules out is the opposite
case, where 450 ppm would have needed an implausible outgassing rate and the
world would quietly draw down to something colder than every run assumes.

Two limits on the check itself. It is a plausibility bound rather than a
derivation -- melt production depends on spreading rate and mantle
temperature, higher gravity compresses the melting column, and 43% land means
less ocean basin than Earth, so the sign of the net correction is not claimed.
And it is computed on a climatology and terrain that both predate the current
build, so the requirement moves on the baseline re-run. The margin is wide, not
unlimited. Re-run it after.

As of 2026-08-19 that check has no artifact at all: `weathering_fluxes.json` and
the climatology it was taken on were deleted with everything else downstream of
the export. The argument above stands as an argument, and the number it produced
does not. `pedology/scripts/weathering_fluxes.py` is the step that restores it.

## `pN2_bar`

```
pN2_bar: 0.78025
```

Partial pressures sum to exactly 1 bar. H2O remains interactive; it is not
included as a fixed surface partial pressure.

## `radiation`

```
radiation:
```


DETERMINED. The k25v spectrum replaced ExoPlaSim's k2.dat, which is the star
K2-18, an M2.5V, not a K dwarf. Two-band surface albedo is required under a red
spectrum and is on.

## `two_band_albedo`

```
two_band_albedo: true
```

Separate visible/near-IR surface albedos, which matter under a red spectrum:
snow albedo splits 0.959/0.614 at maximum and 0.496/0.289 at minimum.

## `stellar_spectrum`

```
stellar_spectrum: k25v
```

`k25v`, in exoplasim/inputs/stellarspectra/, built by
exoplasim/scripts/build_stellar_spectrum.py from the BT-Settl (CIFIST2011)
grid, interpolated log-linearly to the declared 4965 K between the 4900 and
5000 K models at log g 4.5. It puts 0.382 of shortwave below 0.75 um.

This was `k2` through every run of the first three eras, chosen believing
ExoPlaSim ships a measured K2 *dwarf* spectrum matching this star's K2.5V. It
does not: `k2.dat` is the star K2-18, an M2.5V at about 3450 K. The filename
names the object, not the spectral type, and the package ships no K dwarf at
all. The model's own log settles it, reporting 0.1159 of shortwave below
0.75 um where this star should put 0.382.

It matters because nstarfile takes precedence over nstartemp and the spectrum
sets the weighting for every snow, ice and glacier albedo. Those ran 0.10 to
0.17 too dark throughout, weakening the very feedback the glacier and
stellar-cycle machinery exists to resolve. See exoplasim/notes/stellar-spectrum-audit.md
for what that biases and by how much; results predating this remain valid in
kind, with the direction of the error known.

## `ocean`

```
ocean:
```

MIXED, and the block the model reaches through more keys than it looks like.
`horizontal_diffusion` and `horizontal_diffusivity_m2_s` are the cheapest bound
on the missing ocean heat transport and exist to be BRACKETED rather than
tuned; the coefficient is realised at the value it declares, which it was not
before world-mll, so an arm measured earlier carries the wrong label on its
diffusivity.

`salinity_psu` is DECLARED, and it reaches the model through FOUR compiled
constants and not one: the freezing point, sea water's density, sea water's
specific heat and the heat of fusion of sea ice. `run_exoplasim.py` derives the
first three from it and `sea_ice_fusion_j_kg` declares the fourth, because that
one depends on the ice's own salinity and temperature and `icemod.f90` carries
neither as a variable. A salinity bracket therefore moves the mixed-layer heat
capacity and the snow-ice flooding threshold as well as the freezing point.

`cold_start` is DECLARED and is what a run with no restart begins from. This
world has no sea surface temperature or sea-ice climatology and cannot have
one -- an SST field is what the model PRODUCES -- so the cold start is a stated
latitude profile, hemispherically symmetric by construction, and `icemod.f90`
refuses a cold start that declares nothing rather than constructing a field
from a sentinel. The full argument is in the block's own comment; the reason it
is not merely cosmetic is that a hysteresis probe is only meaningful if both
sides of it can be reached, and the constructed field this replaced put every
arm on the iced side.

## `surface`

```
surface:
```


PARTLY DETERMINED. Sea ice and glaciers are decided. `mixed_layer_depth_m` is
NOT: 50 m is a default, it sets seasonal amplitude, and this world's year is
half Earth's so it damps seasonality about twice as hard as Earth's ocean does.
`scripts/error_budget.py` books it as a structural item; pricing it is CLIM-33.

`sea_ice_max_thickness_m` and `sea_ice_lead_closing_m` are DECLARED and are the
two Earth Arctic lengths inside the sea-ice model. The first is negative, which
means no maximum: it is not only a clamp, since the model zeroes the conductive
heat flux once ice reaches it, and the process it stood in for is ice export,
which this model does not have and which the error budget already carries as a
declared structural term. The second is the whole of the lead
parameterisation and sets how much ice must grow before a cell goes white for
albedo and roughness; it is exposed at ExoPlaSim's value rather than
recalibrated, because nothing here can recalibrate it.

## `glaciers`

```
glaciers:
```

Persistent snow converts to glacier, and glacier thickness is added to the
surface geopotential. That orographic term is the point: it is the feedback
that makes an ice sheet grow into its own cold, and without it a cold branch
gets snow albedo but no elevation response. On a world testing albedo-driven
bistability, leaving it off would answer a different question.

initial_height_m -1 places no initial ice, so glaciers only appear where snow
genuinely persists year-round; the module cannot manufacture them. It does
not model ice flow, so continental ice-sheet extent is underestimated.

## `model`

```
model:
```


TRANSITIVE, except where noted. Resolution, layers, ranks and timestep are
operational choices that change with what is being run, not properties of the
world -- they moved from T21/8 to T42/16 between the flux bracket and the
bootstrap. The ozone scale and band weights inside the block ARE determined and
say so.

## `resolution`

```
resolution: T42
```

T42 for the iteration-1 baseline. The bracket ran at T21, which was right
for a several-kelvin question. This pass feeds the carve verdict and the
first LPJ-GUESS input, both of which will change and be redone, so the
more expensive grid is not warranted yet. This does not make T42 a production
constraint. The supported ladder is T21/T42/T85/T127/T170; after the current
support settles, CLIM-52 can map its restart into the next rung and Loops A--C
must settle again there. SPAT-8 chooses the first rung at which the coupled
decision quantities converge, with the ~10M-region, 7.60 km Orogen mesh as the
fine aggregation reference. A converted restart is a spin-up accelerator, not
an equilibrium carried unchanged across resolution.

## `ncpus`

```
ncpus: 16
```

16 physical cores on the Ryzen 9 7950X3D. Not 32: those are SMT threads and
MPI ranks on sibling threads contend for the same FPU and cache, which does
not help a compute-bound spectral model. NLAT must divide by ncpus, so at
T42 the valid choices are 1, 2, 4, 8, 16, 32, 64; 16 leaves 4 latitudes per
rank. Spectral transforms need global transposes, so scaling is sublinear
and may turn over before 16 -- measured against the 8-core baseline below.
RESOLUTION-DEPENDENT, and not simply "all the physical cores". Measured on
this machine, one orbit:

T42, 16 ranks   111 s      T21, 16 ranks    75 s
T21,  8 ranks    53 s

At T21 each of 16 ranks holds about 128 gridpoints, and the global transposes
a spectral model needs cost more than the arithmetic they carry, so half the
ranks run 1.4x faster and leave 8 cores free for a second run. At T42 there
is enough work per rank for 16 to pay. NLAT must divide by ncpus: 32 at T21
allows 1, 2, 4, 8, 16, 32.

## `uniform_land_surface`

```
uniform_land_surface: true
```

Only topography (0129) and the land mask (0172) are supplied. ExoPlaSim's
configure() clears every other surface .sra when a landmap is given, so
roughness, albedo, field capacity, forest fraction and the glacier mask all
take uniform namelist defaults. That is deliberate: Earth's surface maps are
tied to Earth's continents and would be meaningless on this geography.
Declared here so the fallback is a decision, not an accident.

## `energy_diagnostics`

```
energy_diagnostics: false
```

nenergy in plasim_nl, which the Python API does not expose, so
run_exoplasim.py edits the namelist directly. Adds the 28 energy-budget
terms on codes 360-387 to the regular output, to locate which one carries
the residual between the top of the atmosphere and the surface. None does:
the decomposition closes to 0.02 W/m2 and the residual is an offset in the
reported top-of-atmosphere net radiation, and the source is the adiabatic
spectral step.

ON, because `energy_fixer` is driven by `denergy26` and `denergy27` and cannot
run without them. CLIM-1 also requires `close_term_energy.py` to re-measure on
the first low-I/O-off block of any T85 run, since the quantity is resolution-
and timestep-dependent rather than carried from T42.

`2` rather than `true` asks additionally for the CONVERSION DECOMPOSITION, a
control for world-0ov. It prints the adiabatic conversion's reference half both
as the semi-implicit scheme applies it and as it stands at time t, in the same
arithmetic and units as `denergy02`, so the displacement between the two can be
read directly. It costs four extra spectral transforms a timestep and is for a
diagnostic arm, not for production. `true` and `1` are the same setting.

## `conversion_time_level`

```
conversion_time_level: false
```

`nconvtime` in plasim_nl. The adiabatic reference conversion's advective half is
explicit in `calcgp` at time t and its divergence half is the `tkp*c` part of
`tau`, applied on `sdt`. This takes the divergence half back to the state at t so
the two halves meet, and it is world-0ov's first repair route.

IT CHANGES WHAT THE MODEL INTEGRATES. It leaves that half out of the
semi-implicit treatment in the temperature equation while the divergence solve
still treats the temperature implicitly, so the timestep it is stable at is its
own question.

AND THE ANSWER IS NOT THE EXPLICIT GRAVITY-WAVE TIMESTEP. That limit,
`dt < a / (c sqrt(N(N+1)))` with `c = sqrt(R T0 / (1 - kappa))`, is a necessary
condition; the speed in it is within 2 percent of the largest eigenvalue of the
model's own semi-implicit vertical structure matrix, so it is not mis-derived.
The mode this term destabilises is not a gravity wave: `sdt - sd` is the second
time difference, O(dt^2) for a smooth mode and exactly `-2 sd` for the leapfrog
computational mode, so the term feeds that mode, the Robert-Asselin filter is
the only thing damping it, and the boundary is a function of `PNU` as well as of
the rung. At `PNU = 0` there is no stable timestep at all.

The model therefore measures rather than estimates: `plasim.f90`'s
`conversion_time_amplification` iterates its own linearised adiabatic step at
the configured rung, timestep, vertical grid, reference temperature, Robert
coefficient and damping, and refuses a configuration whose fastest mode grows.
`exoplasim/scripts/conversion_time_stability.py` computes the same map without
building or running the model, which is where a caller finds the boundary before
buying a run.

OFF, because the price is a much shorter timestep at every rung and the
operation the sink actually comes from is not yet named (world-pkf).

## `robert_filter`

```
robert_filter: 0.1
```

`PNU` in planet_nl, the leapfrog time filter's coefficient, and the key is in the
PLANET namelist rather than the model one because that is where the model
declares it. Absent leaves what `p_earth.f90`'s `planet_ini` sets, which is 0.1
and not the 0.0 `plasimmod.f90` declares.

Present because the filter had been eliminated as a candidate for the adiabatic
energy sink on the strength of that declaration. Varying it from 0.02 to 0.25
moves the sink by 3.5 percent, so the elimination survives -- but by measurement
now rather than by a misread default.

## `ozone_scale`

```
ozone_scale: 0.794
```

Scales ExoPlaSim's prescribed Earth ozone column, which no part of the model
derives from the host star. 0.794 is Segura et al.
(2003) Table 1:  6.64e18 cm-2 for a K2V host at 1 PAL O2 against the
Sun's 8.36e18. Set through radmod_nl directly, like NENERGY and STARFILE,
since the Python API does not expose it.

This corrects only the column. The absorption coefficients remain fitted to a
solar-shaped band 1, and a 4965 K star puts proportionally less of its band-1
flux in the 200-350 nm Hartley and Huggins bands, so absorption per unit
ozone is still overestimated. See exoplasim/notes/ozone.md.

## `ozone_uv_weight`

```
ozone_uv_weight: 0.335
```

Spectral re-weighting of the Lacis & Hansen (1974) ozone absorptances, which
radmod.f90 uses and which give absorptance as a fraction of TOTAL INCIDENT
SOLAR flux. Each of their three terms therefore carries the Sun's share of
flux in the band it represents, and applying them unchanged to a K dwarf puts
solar band weights on a non-solar spectrum.

The ultraviolet weight is MEASURED, from Segura et al. (2003) Table 2, which
reports incoming ultraviolet at the planet for Sun and K2V hosts placed at
equal total insolation. Summing UV-C, UV-B and the 315-350 nm part of UV-A
gives 17.92 W/m2 against the Sun's 53.44, so 0.335.

That table is epsilon Eridani, observed by IUE: a young, chromospherically
ACTIVE K2V, which is the star this world has been declared to have. So this
is the right weight rather than a bound.

A 4965 K blackbody gives 0.469, HIGHER than the observed 0.335: ultraviolet
line blanketing in a real stellar atmosphere removes more flux than the
chromosphere puts back, so a blackbody overestimates a cool star's
ultraviolet rather than bounding it from below.

The visible weight stays computed, 0.3358 / 0.3673 = 0.914, measured from the
k25v spectrum itself, where line blanketing is already represented and the
chromosphere contributes nothing.

Both default to 1.0 upstream, which reproduces Lacis & Hansen exactly, so the
patch leaves any solar-host run bit-identical.

The same argument applies to water vapour and to CO2, and those are
`h2o_sw_weight` and `co2_sw_weight` below.

## `h2o_sw_weight`

The same spectral re-weighting in the term that is ten times larger than ozone's.
Lacis & Hansen's water vapour absorptance is their Eq. 21, a fit to Yamamoto
(1962), and Yamamoto defines it as a fraction of the SOLAR constant; `radmod.f90`
divides it by `zsolar2` and multiplies it back by band-2 flux, so the two cancel
and a K dwarf gets the Sun's absorbed fraction of total flux. Sets `H2OSWW`,
which defaults to 1.0 and reproduces Lacis & Hansen exactly.

Derived in `exoplasim/notes/shortwave-water-vapour.md`, which also carries the
bracket, the checks that partly failed and what survives them, and the reason
the value belongs to the k25v spectrum rather than to the blackbody `solarini`
builds when `NSTARFILE` is 0.

## `h2o_sw_level`

A SECOND correction to the same water vapour term and a separate decision from
the weight above, which is why it is a separate key. `h2o_sw_weight` is a
star-over-Sun ratio, so an error in the absolute level of Lacis and Hansen
Eq. 21 divides straight out of it and nothing else in the scheme touches it.
Sets `H2OSWL`, which defaults to 1.0 and is Eq. 21 unmodified.

Eq. 21 is low, and four independent constructions of the quantity say so: this
project's reconstruction from Howard's bands, HITRAN2020 correlated-k on the same
path, and Ramaswamy and Freidenreich (1992), who scored Lacis and Hansen against
a line-by-line reference and found them under-absorbing. None of those carries a
water vapour continuum, and neither does Eq. 21, so all of them are floors: the
continuum absorbs in the WINDOWS between the bands, which a band sum has no
absorption in at all.

**The bracket is what the continuum is not known to.** Its ends are arms to run
rather than an error bar on a settled number, and closing it needs a line-by-line
calculation with a continuum on this path, which is a correlated-k bundle this
project does not have. `exoplasim/notes/corrk-cross-check.md` derives the value
and the bracket and states what each source contributed.

There is no shortwave continuum term to set instead. `radmod.f90` has a water
vapour continuum coefficient, `th2oc`, and it is in `lwr` only; `swr` has ozone
in band 1 and water vapour in band 2 and nothing else, so the continuum reaches
the shortwave through this key or not at all.

## `co2_sw_weight`

```
co2_sw_weight: 1.510
```

The third of the same family, and the one that is a NEW ABSORBER rather than a
re-weighting: `swr` has no shortwave CO2 at all, because Lacis & Hansen did not
parameterise it and the port is faithful. Sets `CO2SWW`.

**Its default is 0.0 and not 1.0, and that is not an inconsistency.** The other
weights scale an absorptance the scheme already has, so 1.0 is the solar value
and also the no-op. This one scales an absorptance the scheme does not have, so
0.0 is the no-op and 1.0 is the solar-weighted term. Anything that treats the
four keys as interchangeable will invert this one.

Derived in `exoplasim/notes/shortwave-co2.md`, which also carries the 2.7 um
overlap decision, the Earth-column check the derivation had to pass, and the
statement that this correction and `h2o_sw_weight`'s have the same sign.

## `land_longwave_emissivity` and `sea_longwave_emissivity`

```
land_longwave_emissivity: <the lithology-weighted land mean>
sea_longwave_emissivity: 0.98
```

`lwr` builds a per-cell surface emissivity from these two, one for the land
fraction of a cell and one for the rest, and it wrote them as one literal until
they were named: the modelled land emitted as a perfect blackbody by
construction and everything the land-sea mask did not call land took 0.98, with
no comment, no unit and no source anywhere in the tree.

THE LAND VALUE IS DERIVED. `analysis/rock_emissivity.py` integrates one minus
the directional-hemispherical reflectance of the ECOSTRESS spectra against a
Planck function, per Orogen rock class, and the key is that table area-weighted
over the active build's land. Only hemispherical measurements are read, because
Kirchhoff needs the whole hemisphere and a bidirectional reflectance gives an
upper bound on emissivity rather than an estimate of it. Every class carries a
solid-to-particulate preparation bracket, and the two arms are the sensitivity
pair to run rather than a guessed one.

ONE SCALAR IS ENOUGH FOR THE LAND, AND THAT IS MEASURED.
`analysis/emissivity_contrast.py` asks what a per-cell field out of the same
lithology map would buy OVER a scalar set at the field's own land mean, in the
units the surface energy balance reports. Surface net longwave is
`-eps*(sigma*Ts^4 - LWdown)` and is exactly linear in the emissivity, so the
answer needs no run. The criterion is 1.4 W/m2, the top of the model's own dry
adiabatic energy sink, and the field misses it by a factor of three at T21 and
by a factor of eight at the bound where the spectra's out-of-band behaviour is
a blackbody. The blackbody land surface it replaces was worth several times the
whole contrast a field would have bought, one-signed on every land cell. So the
number was the defect and the shape was not;
`notes/audits/surface-longwave.md` carries the measurement and `world-vhhs`
re-runs it at a finer rung, where cell averaging removes less of the contrast.

The sea value is still declared rather than derived. Sea water is 0.985 to
0.99 in the thermal window, so 0.98 is slightly low, and the bracket to run an
arm over is 0.98 to 0.99.

Both are written into every run's `radmod_namelist` unconditionally and checked
there by `verify_staged_namelists`, because a continuation that drops them
silently returns the modelled surface to the compiled literal -- which for the
land is now a blackbody rather than the value this file carries.

## The ozone profile's two lengths, `BO3` and `CO3`

This file states neither. `run_exoplasim.py:ozone_profile_lengths_m` reads
upstream's own compiled `bo3` and `co3` out of `radmod.f90` and scales both by
`lib/lapse.py:pressure_length_ratio_to_earth`, so the two halves of the
derivation are read where they live and the product is written down nowhere.

`mko3` places its synthetic ozone profile with a logistic in geometric height,
centred at `BO3` with width `CO3`, both in metres. It builds the height
coordinate it compares them against hypsometrically from the model's own
`gascon` and `ga`, correctly; upstream's two constants are Earth's and do not
follow, so on a higher-gravity world the whole profile is pushed toward the
model top.

Height is the wrong invariant. The modelled photochemical maximum is set by
pressure-like conditions, ultraviolet optical depth and three-body
recombination density, so pressure is what to hold fixed, and holding it means
scaling both by `(gascon/ga)` over Earth's own -- the same factor the heights
already carry. That ratio moves with the declared composition and gravity, so a
number in this file would be a copy that could not learn either had changed. Doing so reproduces Earth's own layer-by-layer share of the
ozone column on this world's sigma levels to machine precision, which is the
form of the claim that could have failed;
`exoplasim/notes/ozone.md` carries the recursion and the numbers.

This is a REDISTRIBUTION and not a column change. `ozone_scale` is the column.
The budget effect at the surface is second-order and the top model layer's
shortwave heating rate is where it shows, which is the layer setting the
modelled tropopause temperature and static stability.

Written into `radmod_namelist` by `configure_otherargs`, so they are reapplied
on every continuation, and checked there by `verify_staged_namelists`. They are
NOT routed through `configure(ozone=...)`: that path also rewrites `A0O3`,
`A1O3`, `ACO3` and `TOFFO3` from the same dict, and those four are a separate
question about profile SHAPE that `exoplasim/notes/ozone.md` records and this
does not touch.

## `energy_diagnostics_3d`

```
energy_diagnostics_3d: false
```

The same 28 terms per level, codes 460-487. The column totals established
that no single term carries the residual and that the decomposition closes,
which is what ruled the gridpoint physics out. Per-level is what showed the
large-scale condensation lead to be a phase-booking difference of a flat
13.3% at every sub-freezing level rather than a leak.

OFF as of 2026-08-19 with the column set, and this is the half that cost
something: carrying a level axis, it was about two thirds of an orbit's
output bytes on its own. Both findings above are recorded, which is what
makes the field disposable rather than the diagnostic being unwanted. Back
on for T85 on the same condition as `energy_diagnostics`.

## `roughness_source`

```
roughness_source: lithology
```

Aerodynamic roughness from land cover and subgrid relief, code 173, by
build_surface_roughness.py. Replaces the uniform dz0land = 2.0 m, which
asserts forest-scale roughness over the salt-crust and playa share of land
(`world_state.json` has the current figure). Those are closed-basin floors, flat by construction, and they are the
cells the carve verdict integrates evaporation over: the default gives them
about 8x the turbulent exchange a real playa surface has. The reference height
that ratio is taken over is `lib/lapse.py:reference_height_m`, at this planet's
gravity; the builder's report brackets it, because the height is linear in the
lowest-level air temperature and the field is built before any run measures one.

NVEG is 0 in these runs, so landmod.f90:409 takes dz0 = dz0clim directly and
173 is the field that matters; SIMBA's separate vegetation and orographic
terms never run.

The land mean is DERIVED and not anchored. The orographic term is the
turbulent form drag this world's own subgrid slope exerts, through Wood and
Mason (1993) Eq (33) and Beljaars et al. (2004) Eq (6), with no free
coefficient in it; the cell average is of drag coefficients at Mason (1988)'s
blending height. The land mean of PlaSim's own Earth boundary dataset is
reported beside it as a comparison, and the builder's report carries the
distance, both ends of the declared bracket and the land means the field
spans. `notes/audits/tuned-values.md` section 9 has the derivation and what
each paper settled.

## `soil_water_source`

```
soil_water_source: pedology
```

`uniform` for a BOOTSTRAP, `pedology` for the baseline, and the flip goes
between them. Soil water capacity comes from the pedology soil map as code 229, and that
soil has to be weathered under a climatology, so on a terrain that has none
yet the field cannot exist and run_exoplasim refuses to start without it.
The flip cannot happen mid-run either, since continue_exoplasim compares the
config against the run's manifest key by key and will not resume across a
changed value, so the bootstrap has to be
finished before it moves. `docs/src/pipeline/sequencing.md` loop A has the order.

It has to be `pedology` BEFORE the run whose climatology the verdict uses,
not merely before the verdict: the land-mean derived bucket depth is several
times shallower than the uniform 0.5 m one it replaces, which changes
evaporation, which changes P - E, which is the numerator of the carve criterion.
The depth is a property of the soil map, so it moves with the build and the
rung; `exoplasim/scripts/build_surface_soil_water.py` writes it as
`land_mean_capacity_m` onto the provenance sidecar beside the staged field.

## `land_albedo_source`

```
land_albedo_source: vegetated
```

Background land albedo, from lithology rather than ExoPlaSim's uniform 0.22.
Rock classes run 0.10 for basalt to 0.50 for evaporite, so a uniform value is
a planetary-albedo error of a few hundredths.

The world is taken to be vegetated. What a canopy reflects is `vegetation_albedo`
and no longer a constant here; the classes that cannot carry one are listed in
`barren_rock_classes` and keep their own albedo. The bracket showed this is not
a detail: the bare-rock and vegetated endmembers differ by 3.7-7.1 K and reach a
given design mean at non-overlapping fluxes, so the biosphere and the orbit are
one choice, not two.

What the land mean comes out at is per build and belongs in the albedo report
beside the .sra files, not here. It moves with the lithology, which is why
the flux is re-derived on every new terrain rather than carried.

## `lithology_albedo_overrides`

```
lithology_albedo_overrides:
  playa_clastic:
    albedo: 0.23
    replaces: 0.19
```

DETERMINED, and the key exists because one class is worth overriding. Rock
class albedo is the generator's: `vendor/orogen/js/lithology.js` carries a
value per class, salt crust among them, and this project reads them rather than
restating them.

The exception is `playa_clastic`. It is the largest single albedo lever on this
terrain, because it is the biggest lithological share of the land -- an earlier
override on `evaporite` was retired when the generator zoned that class into
crust and clastics separately, and the lever moved with the area. The
generator's value comes from Post et al. (2000), 52 pyranometer measurements
over 0.3-2.8 um: the right KIND of measurement on the wrong population, because
it is a field mean over all soils rather than over saline desert soils.

The override is that field level corrected by how much brighter saline desert
soils are than the soil population in ECOSTRESS, measured within one library so
the laboratory-over-field offset cancels, then re-weighted to this star.
`analysis/playa_albedo.py` is the derivation and
`notes/audits/orogen-lithology.md` is the argument.

`replaces` is the guard rather than documentation: it names the value the
override expects to find, and refuses to fire if the generator has redefined the
class underneath it. That is what retired the previous override rather than
leaving it to be noticed. If this class moves again, the flux moves with it.

## `barren_rock_classes`

```
barren_rock_classes: [evaporite, playa_clastic]
```

DETERMINED. Rock classes that cannot carry vegetation, so they keep their own albedo in
`vegetated` mode rather than being handed the canopy value. Nothing roots in
salt crust and nothing much roots in playa mud, and between them they are
a large share of this planet's land -- `world_state.json` says how large on the
active build -- so treating either as vegetated would be a first-order error in
the energy balance.

This started as a single hardcoded reference to evaporite and became wrong
the moment Orogen split that class into crust and clastics. Driven from
config now so the next addition is a config change rather than a silent one.

## `geography_land_threshold`

```
geography_land_threshold: 0.5
```


DETERMINED. The land fraction at or above which a model gridcell counts as
land. The export is a mesh and the model is a grid, so every cell arrives with a
fractional land area and something has to make it binary;
`build_boundary_conditions.py` and `build_surface_albedo.py` both threshold at
this value, and they read the same key so that the land mask and the albedo
field cannot disagree about which cells are land.

0.5 is majority-rule and is a threshold rather than a measurement, which is why
it is a decision and lives here. What it costs is a coastline effect: cells are
hundreds of kilometres across, so it rounds partial coasts either into the sea
or onto the land, and the two scripts sharing one key is what keeps that
rounding consistent rather than correct.

It is also the one value that makes the rounding SINGLE. `oceanmod.f90`
hard-binarises `yls` at 0.5 whatever this builder writes, so a builder
threshold anywhere else means the mask written and the rounding the model would
apply to a fractional field disagree, and the cell is rounded twice by two
rules. That is not a preference between values; it is the difference between
one rounding and two.

The size of the effect and the price of every alternative are in
`notes/audits/coastline-threshold-cost.md`, and
`build_boundary_conditions.py:coastline_ledger` reproduces both in every
boundary-condition report rather than leaving them in prose. Three rules were
weighed and the note carries the arithmetic: moving the threshold so land AREA
is conserved, exempting cells that hold below-datum land from the threshold,
and leaving it here. The exemption is the only form that keeps
`source/README.md`'s first rule intact, and it is refused because at the rungs
this world is run at it misassigns several times more surface than it rescues.
The route out of the residue is the RUNG, which is the only term in the ledger
that falls fast with support.

This is separate from the land mask rule that matters most in this project.
Land for ANALYSIS comes from `surface_class` in the export and never from
`land_mask`, because the two disagree over dry closed-basin floor below sea
level; see `source/README.md`. This key is about gridding a fraction, not about
which of those two definitions is right.

## `stellar_cycle`

```
stellar_cycle:
```


PARTLY DETERMINED. The activity level is chosen; see `star.activity`.

PROVISIONAL, 2026-08-16. Periods and amplitudes are CHOSEN -- two components
at 11 and 57 Earth years, 2.5% and 3.5% peak-to-peak -- and centred on the
baseline flux. The block stays PROVISIONAL rather than DETERMINED because
nothing has been run with it yet: what the first cycle run has to return
before it settles is the real damping factor, which is analytic and
Planck-only here; whether the cycle-mean temperature sits where the convexity
argument says it does; and whether the aligned excursion is tolerable.
Amplitudes may want revising afterwards, and revising them is expected rather
than a failure.

The amplitude is a stipulated property of the star -- a worldbuilding choice,
defensible for a young active K dwarf, not an estimate -- and the period is
stated in absolute time rather than local years, since a stellar dynamo does
not know about this planet's orbit.

The mechanism -- a secular luminosity change against a spot and facular cycle
-- is decided for modelling purposes: the cycle is applied as a grey
multiplier, because a spot-driven swing of this size moves the band-1 fraction
by 0.007 and is worth under 0.02 K at this world's ice cover. That
measurement, and the ice fraction at which it would stop holding, are in
exoplasim/notes/parameter-decisions.md. The physics that made it a real
question: a secular change moves effective temperature as F^(1/4) and needs
the spectrum rebuilt past a couple of percent, while a spot cycle is strongly
non-grey and concentrated in the blue -- the band snow and ice albedo respond
to, and the whole reason a measured spectrum is used at all -- so the two give
different cryosphere responses at the same flux.

Convergence for a cycle run is undefined in drift terms: drift-based criteria
assume an approach to a steady state, and a forced cycle has none. Validate
against periodicity instead, comparing a year against the same phase one cycle
later.

## `components`

```
components:
```

Two superposed sinusoidal components of bolometric variation about the
baseline flux, applied at every radiation timestep by
patches/exoplasim-3.4.2-star-cycle.patch. Two rather than one, because
epsilon Eridani genuinely has both a short and a long cycle, and because the
two do different jobs here: the medium one is climatic and the long one is
geomorphic.

The periods are NOT epsilon Eri's values. That star is a reference for what a
K2.5V of this activity can do, not a template to copy. What was taken from it
is the plausible amplitude envelope -- 4.3% median, 2.4-7.2% at 95%,
bolometric. The 6% total here is inside that and above its median.

The ratio is deliberately non-integer. At 57/11 = 5.18 the phase relationship
returns only after about 314 Earth years, so grand minima differ in depth
rather than repeating identically. That is the point: glaciers advancing
further in the deeper minima leave moraines at different distances, so the
landscape records which past minima were severe. An exact ratio erases that.

Damping below is ANALYTIC and Planck-only, from a slab thermal timescale of
1.17 Earth years. Measuring the real value is what the first cycle run is for.

Per-component kelvin consequences are not recorded here. The canonical local
sensitivity is `lib/sensitivity.py`; corrected magnitudes are in
`exoplasim/notes/parameter-decisions.md`. THE AMPLITUDES BELOW ARE THE
DECISION.
