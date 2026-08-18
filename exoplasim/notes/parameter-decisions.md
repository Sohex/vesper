# Baseline parameter decisions

## Source geography

Everything in this section describes the equirectangular PNG maps that the
completed runs were built from. Those maps came from a 510k-region World Orogen
build and are the files now in `../source/maps/`. The canonical geography has
since moved to full data exports under `../source/exoplasim-*/`, from a
2.5M-region build with preserved endorheic basins and lithology-modulated
erosion, which supply Gaussian grids directly and make the remapping below
unnecessary. See the project `CLAUDE.md`.

The files in `source/` are canonical and remain untouched. World Orogen's
exporter documented and implemented these encodings:

- 16384 x 8192 equirectangular maps, north at the top and the antimeridian at
  the left/right edges (0 degrees longitude at map center).
- Full heightmap: linear 16-bit grayscale from -5 km to +6 km.
- Land-only heightmap: linear 16-bit grayscale from 0 to +6 km, with ocean
  black.
- Land mask: opaque RGBA, exactly white for land and black for ocean.

The mask is flat-shaded by Orogen region while heightmaps use smooth vertex
interpolation. Narrow coastal disagreements are therefore expected. The
conversion uses the explicit mask for land/sea and the land-only heightmap for
physical land elevation. Bathymetry is not needed by ExoPlaSim's slab ocean.

The conversion is conservative in spherical area. It integrates the source
mask and land elevation over T42 cells defined by Gaussian quadrature latitude
weights, then uses a 50% land threshold. Topography is the land-area-weighted
mean elevation in retained land cells. Physical height is converted to the
geopotential expected by SRA code 129 using this planet's gravity.

## Planet and star

Surface gravity is declared, not calculated from mass and radius. World Orogen
ran with `10.1989 m s-2` and scaled the terrain's maximum relief as 1/g, so that
value is a property of the geography rather than a free parameter, and mass
follows from it: `1.4976 = (10.1989 / 9.80665) * 1.20^2` Earth masses.

The earlier revision of this file declared 1.50 Earth masses and calculated
`10.2152604 m s-2`, which is 0.16% away from the gravity the terrain was
actually built under. Every run predating this change used the calculated value.
`derive()` now reads `planet.gravity_m_s2` and refuses to run if
`planet.mass_earth` disagrees with it.

The provisional star is an empirical K2V/K3V midpoint based on the Pecaut--
Mamajek dwarf sequence: 0.80 solar masses, log luminosity -0.49 (0.32359 solar
luminosities), 4965 K, and an implied 0.769 solar radii. This is internally
consistent with the requested mid-K dwarf. The earlier 0.82 solar masses,
0.45 solar luminosities, and 4900--5000 K combination implies an oversized
star and is closer to an evolved or earlier-K object.

Orbit distance and period are derived separately for every flux case. The
model uses 1361 W m-2 as modern Earth flux.

| Flux | Stellar flux | Semimajor axis | Orbital period |
| --- | ---: | ---: | ---: |
| 0.85 S-Earth | 1156.85 W m-2 | 0.617008 AU | 197.9198 Earth days |
| 0.90 S-Earth | 1224.90 W m-2 | 0.599624 AU | 189.6145 Earth days |
| 0.95 S-Earth | 1292.95 W m-2 | 0.583631 AU | 182.0794 Earth days |

## Atmosphere and model representation

Partial pressures sum to 1.0 bar: N2 0.78025 bar, O2 0.21 bar, Ar 0.0093 bar,
and CO2 0.00045 bar (450 ppmv). ExoPlaSim uses the full mixture to derive mean
molecular weight, gas constant, and total pressure. Only CO2 and interactive
H2O are explicit radiative absorbers; Rayleigh scattering remains Earthlike.
O2 and Ar are not chemistry or radiation species. Ozone is a prescribed
Earthlike climatology, not photochemistry.

The stellar spectrum is represented by ExoPlaSim's 4965 K blackbody option,
not by a detailed K-dwarf atmosphere spectrum. Separate visible and near-IR
surface albedos are enabled; this matters for the weaker snow/ice reflectivity
under a redder K-dwarf spectrum.

The baseline uses a 50 m thermodynamic slab ocean, interactive sea ice, ten
atmospheric levels, double precision, and T42. A two-way exponential physics
filter is enabled because the custom topography has sharp small-scale features.

## Non-Earth rotation/calendar workaround

ExoPlaSim 3.4.2's wrapper sets the physical `SIDEREAL_YEAR` correctly, but if
`rotationperiod != 1` it then resets `N_DAYS_PER_YEAR` from a hard-coded
360-Earth-day run. It also defaults each Python-level "year" to 360 Earth days.
Run generation must override:

- `N_RUN_STEPS` to one physical orbital period, and
- `N_DAYS_PER_YEAR` to the nearest number of 30-hour sidereal rotations in
  that orbit.

The physical orbit and run-step count control seasonal forcing. Later seasonal
analysis will use orbital phase from snapshot output rather than assume Earth
calendar month lengths.

ExoPlaSim 3.4.2 also accepts integer postprocessor codes according to its API
documentation but silently omits derived variables (including total
precipitation, vector winds, and net-energy diagnostics) when integers are
passed. Run generation passes the same codes as strings and validates the
required NetCDF schema after every run.

Although the postprocessor library names codes 165/166 as 10 m vector winds,
ExoPlaSim 3.4.2 has neither a native output producer nor a derivation branch
for them. Prevailing-flow maps therefore use the lowest model level
(sigma about 0.983), labeled as such, together with native surface stress.

## Equilibration

ExoPlaSim provides `runtobalance()`, whose default test is less than 0.05 W m-2
per modeled year of drift in both top-of-atmosphere and surface energy balance
over a long baseline. For this workflow, the first smoke run is not interpreted.
Production runs retain annual diagnostics for global surface temperature,
TOA and surface net energy, ocean/slab drift, and sea ice. Equilibrium requires
small multi-decadal trends rather than a single quiet year, followed by a
separate multi-orbit climatology window.

For the T42 worldbuilding baseline, the explicit ten-orbit convergence test
requires temperature drift below 0.05 K/orbit; TOA and surface-balance trends
below 0.05 W m-2/orbit; sea-ice drift below 0.001 of planetary area/orbit; and
mean absolute TOA and surface imbalances below 0.5 W m-2. These criteria are
stricter than merely observing one balanced year, but are not presented as a
publication-grade coupled-ocean equilibrium test.

The 0.90-S-Earth baseline passed after 46 modeled orbits (indices 0--45). Over
the final ten-orbit test window the temperature trend was +0.0246 K/orbit, TOA
balance trend -0.0154 W m-2/orbit, surface-balance trend -0.0203 W m-2/orbit,
and sea-ice trend -0.000309 of planetary area/orbit. Mean TOA and surface
imbalances were -0.262 and +0.0216 W m-2. Five new orbits (46--50), not used in
the pass/fail decision, were then retained with 32 snapshots per orbit and
averaged into the analysis climatology.

## Climate interpretation conventions

Astronomical seasons are assigned by simulated solar longitude: 0--90 degrees
is northern spring, 90--180 northern summer, 180--270 northern autumn, and
270--360 northern winter. The model's orbital phase variables are used rather
than Earth calendar dates.

Köppen--Geiger precipitation thresholds were developed for Earth-year totals.
Because this planet's orbit is roughly half an Earth year, applying them to
precipitation per local orbit would make identical physical rain rates appear
about twice as dry. The diagnostic therefore annualizes modeled precipitation
rates to 365.2425 days while retaining the simulated seasonal temperature and
rainfall cycle. The resulting classes and biome groupings are heuristic ecological
interpretations and do not model soils, nutrient limitation, fire, vegetation
competition, or ocean heat transport.

PlaSim output code 182 uses negative values for upward evaporation. Analysis
plots positive evaporation as `-evap` and computes `P-E` as `pr + evap` in the
native sign convention. Very small negative runoff values created by spectral
ringing are clipped to zero for display but remain untouched in NetCDF source
products.

## Completed flux sweep

The 0.95-S-Earth endpoint passed the same six fixed convergence criteria after
46 spin-up orbits (indices 0--45). Its five independent climatology orbits
(46--50) average 292.079 K, 2.904 mm/day precipitation, -0.372 W m-2 TOA net
radiation, +0.034 W m-2 surface heat flux, and 0.00264 planetary sea-ice
fraction. This case naturally satisfies the provisional 290--293 K design
target at 450 ppm CO2, so no greenhouse tuning is warranted before regional
interpretation.

The 0.85-S-Earth endpoint was integrated much longer because its ice-rich
branch relaxed slowly. At index 75, its final ten-orbit trends were
-0.000040 K/orbit temperature, +0.001998 W m-2/orbit TOA balance,
+0.01018 W m-2/orbit surface balance, and -0.000168 planetary sea-ice
fraction/orbit. Five criteria passed, but the mean TOA balance was
-0.504154 W m-2, 0.004154 W m-2 outside the fixed absolute 0.5 threshold.
It is therefore recorded as quasi-equilibrated rather than a strict pass.
Climatology indices 76--80 average 258.797 K, 1.113 mm/day precipitation,
-0.573 W m-2 TOA net radiation, -0.113 W m-2 surface heat flux, and 0.3254
planetary sea-ice fraction.

The sweep's 22 K jump from 0.85 to 0.90 and 11 K jump from 0.90 to 0.95,
together with the collapse in sea ice from 32.5% to 6.3% to 0.26%, indicate a
strong nonlinear ice-albedo response in this model/configuration. Three points
are sufficient to identify the broad regimes but not to locate transition
thresholds precisely or justify linear interpolation.

## Time-varying stellar-cycle experiments

The first stellar-cycle comparison holds the established 0.90-S-Earth orbit
fixed (0.599624 AU and 189.6145 Earth days) and varies the star's
disk-integrated bolometric irradiance sinusoidally. Both experiments use an
8.000-Julian-year period, equal to 15.410 local orbital periods. The extreme
case spans 0.85--0.95 S-Earth and the control spans 0.88--0.92 S-Earth. Each
begins from the final equilibrated 0.90 restart at mean flux with irradiance
increasing. Four full cycles are required before interpretation, with cycles
three and four compared at equal phase to assess periodic convergence.

ExoPlaSim does not natively represent starspots. A pinned, project-local patch
to ExoPlaSim 3.4.2 applies the sinusoid at every radiation timestep using the
absolute timestep stored in restart files. This avoids discontinuous annual
forcing changes and phase resets. The stellar spectral partition remains the
4965 K blackbody configuration throughout: these runs isolate bolometric
forcing and do not represent spot-induced color changes, UV variability,
flares, or rotational spot modulation.

The 0.85--0.95 range is treated as an intentionally extreme active-star
worldbuilding experiment rather than typical behavior for a quiet 6--8 Gyr K
dwarf. Its period is observationally ordinary for K dwarfs, but its 11.1%
peak-to-peak bolometric variation would imply unusually extensive and variable
spot coverage. The 0.88--0.92 case tests whether any large response depends on
that extreme amplitude.

Both integrations completed 62 local orbits, or 4.024 stellar cycles. Over the
final two cycles, the extreme case averages 280.73 K and its local-orbit mean
temperature spans 276.68--284.69 K; the control averages 280.96 K and spans
279.33--282.64 K. Their seasonally adjusted temperature responses are 7.96 K
and 3.24 K peak-to-peak, respectively, with temperature maxima about 1.44
Earth years after flux maximum. The slab ocean therefore prevents either run
from approaching the static 0.85 and 0.95 equilibrium endpoints during each
half-cycle.

Equal-stellar-phase comparison of cycles three and four gives temperature
phase-curve RMSEs of 0.069 K (extreme) and 0.153 K (control). Final-two-cycle
mean TOA imbalances are -0.463 and -0.543 W m-2. The response is settled enough
for worldbuilding interpretation but is described as quasiperiodic rather
than a strict periodic radiative equilibrium. The orbital year is not an
integer divisor of the activity cycle, so seasonal state and activity-cycle
phase do not exactly repeat together.

## Land surface fields are uniform, deliberately

ExoPlaSim's `configure()` runs `rm workdir/*.sra` whenever a landmap or topomap
is supplied (`__init__.py:2938`) and then writes back only those two. Every other
surface field is therefore absent from the run directory, and `surfmod.f90:125`
reads each one with `inquire(exist=)` and skips it silently when it is missing.

The completed runs consequently used uniform land-surface properties: 2.0 m
roughness, 0.22 background albedo, Earth field capacity, 0.5 forest fraction and
no glacier mask. `landmod.f90:304-308` presets those before attempting the read,
so the fallback is a set of namelist values rather than zeros, and nothing about
the output is corrupt.

This is now a declared choice rather than an accident. `model.uniform_land_
surface` in `config/planet.yaml` must be true, and `surface_field_report()`
refuses to run otherwise and records the fallback values in every run manifest.

Albedo is the exception, and is now supplied. See below.

Keeping the rest uniform is the defensible option. Earth's roughness, albedo and
vegetation maps are tied to Earth's continents; imprinting them on this geography
would place forests and deserts by coincidence of coordinates. A geographic land
surface should come from this project's own biome work, not from the shipped
Earth climatology.

One consequence for resolution: the wipe happens regardless of truncation, so the
fact that ExoPlaSim ships surface data only up to T42 does not constrain us. A
T63 or T85 run would default to exactly the same uniform values a T42 run already
does. Resolution can be chosen on cost and on what the geography justifies.

## Run directories are keyed to the geography

`finalize()` selects output as `sorted(glob("MOST*"))[-1]` (`__init__.py:1501`),
so a shorter run in a directory that still holds a longer previous run copies the
older world's final year out under the new world's name. `crashtolerant` has the
same exposure through `cp MOST_REST.%05d plasim_restart`.

We do not call `finalize()`, and the prepare step refuses a directory containing
any of `MOST.*`, `MOST_REST.*`, `MOST*DIAG*`, `snapshots/` or `highcadence/`
regardless of `--force-prepare`. Beyond that, `run_id` now ends in a digest of the
two boundary-condition SRA files, so a change of geography lands in a different
directory by construction rather than by vigilance. Runs made before that change
keep their old names and are not resumable, which is correct: they predate both
the current geography and the current gravity.

## Background land albedo comes from lithology

Uniform albedo was the one part of the uniform land surface that could not be
justified. ExoPlaSim's `albland` default is 0.22; area-weighted bare-rock albedo
over this planet's land is **0.315**, ranging from 0.10 for basalt to 0.50 for
evaporite. Written as a planetary figure that is an albedo error of +0.041, or
about -12.5 W/m2 of absorbed flux at 0.90 S-Earth before snow, cloud and
vegetation feedbacks. For comparison the entire 0.85-to-0.95 flux sweep spans
roughly 21 W/m2 absorbed, so a uniform 0.22 is a first-order error rather than a
refinement.

The spread is not incidental to this world. Evaporite covered 20.8% of the land
on the pre-carve terrain this was computed on; closed-basin fill is 16.5% on the
current `carved-zoned-v4`, split 1.9% salt crust and 14.5% playa clastics,
second behind schist at 22.4% and ahead of intracratonic clastics at 18.9%, and
evaporite is the brightest class in the table at 0.50. (The manifest's
`compositionLand` says 18.6%, because it measures against `land_mask` and so
drops the dry sub-sea-level basin floors where playa fill concentrates. The
figures here are computed from `substrate_class`, then named `surface_rock`, against `surface_class`.) That is a direct consequence of the
endorheic drainage: closed basins accumulate playa and salt-pan fill. A planet
whose land mostly does not drain to the sea is a brighter planet, and the two
facts come from the same place.

`build_surface_albedo.py` writes codes 174, 175 and 176 from the export's
`rock_albedo` field, whose own description names it "a surface boundary condition
for the climate stage". With `NSIMPLEALBEDO=0` the radiation uses the two-band
pair 175/176; 174 is written as well so the broadband diagnostic agrees rather
than silently keeping 0.22. Both bands get the same value, because we have one
albedo per rock class and no spectral split; that reproduces single-band
behaviour while still varying geographically.

Staging order matters. `configure()` clears `workdir/*.sra` whenever a landmap is
given, so these files are copied in *after* that call, by `stage_surface_extras`.
Anything staged earlier is silently deleted.

### The caveat, which is not small

This is substrate albedo, not land-surface albedo. A vegetated Earth-like surface
sits nearer 0.12 to 0.20. Writing bare rock produces a genuinely bare-rock
planet, which is honest when no vegetation model has run, but it will be cold
relative to a vegetated world, and that cold bias propagates into whatever
vegetation model consumes the output.

### Bracket it before iterating

The gap between the two physical endmembers is larger than it first looks. Bare
rock gives a land mean of 0.315. A fully vegetated surface, with evaporite left
bare because nothing grows on a salt pan, gives 0.223. That is 12 to 16 W/m2 of
absorbed flux depending on the vegetation albedo assumed, against **21 W/m2 for
the entire 0.85-to-0.95 stellar sweep, which produced a 33 K range**.

So land albedo is not a bootstrap detail to be tuned away. It is comparable to
the largest forcing this project has deliberately varied, and it acts through a
positive feedback: darker means warmer means more vegetation means darker. This
world has already shown a strong ice-albedo nonlinearity over that same forcing
range, with sea ice going 32.5% to 6.3% to 0.26% across the sweep. A forcing of
similar size through a second positive feedback is the standard setup for more
than one stable state.

If that is the case, iterating from a single starting albedo does not find the
answer. It finds the basin of attraction the starting point was in. The method
therefore is to run both endmembers as static climates first, with no vegetation
model in the loop, and compare:

- If the two climates land within a few K, the vegetation feedback is a
  correction. Pick the branch, iterate twice, done.
- If they diverge sharply, or one glaciates and the other does not, the world is
  bistable. That is a first-class worldbuilding finding rather than a problem,
  and the state is then chosen deliberately and on the record.

Only then run the vegetation loop, warm-starting within a branch, which is
correct because it follows that branch, and never across branches, which would
destroy the test.

One coupling to keep in view: 20.2% of the land is evaporite, and evaporite forms
in closed basins. If the carve verdict removes many basins, there is less
evaporite and the world is darker. So the drainage hypothesis, the albedo and the
climate are a three-way loop, not two separate two-way loops, and all three are
gated on the same first climate pass.

`model.land_albedo_source` chooses which way to be wrong:

- `lithology` writes bare rock as exported. Physically what the surface is before
  anything grows on it, and cold. This is the default, because the first pass
  exists to feed a vegetation model and a bare planet is what that model should
  be handed.
- `vegetated` is the opposite endmember: vegetated ground everywhere except
  evaporite. Land mean 0.223. Physical, warm, and the other half of the bracket.
- `scaled` keeps the lithology pattern but rescales the land mean to
  `--target-mean`. Retained for experiments and **not recommended in
  production**: it pins the mean to a guess at the answer, which is circular, and
  it leaves every cell somewhere physically wrong, making basalt darker than any
  real rock and evaporite far too dark for a salt pan.
- `uniform` accepts ExoPlaSim's 0.22 and writes nothing.

Part of the lithology signal survives vegetation regardless: nothing grows on a
salt pan, so the bright evaporite basins stay bright. That argues for carrying
the pattern through the loop rather than flattening it.

## Boundary conditions are built from the mesh, not the maps

`build_boundary_conditions.py` replaces `convert_orogen.py`. Two things changed.

Land comes from `surface_class`, not from the PNG mask's `elevation > 0` test.
The two disagree by 1.9% of the planet, all of it dry basin floor below sea
level, and the old convention put it under water.

Both fields are integrated from the native 2.5M-region mesh rather than sampled
from the gridded export. The export resamples categorical fields, `surface_class`
among them, by taking the region containing the cell centre; at T42 a cell holds
roughly 300 regions, so a point sample discards the coastline. Land fraction per
cell is now the area-weighted fraction of its regions that are land, thresholded
at 0.5, and topography is the land-area-weighted mean elevation over land regions
only, so ocean depths never drag a coastal cell down.

Dry basin floors keep negative elevation. 119 T42 land cells sit below sea level,
the deepest at -58 m; the mesh reaches -562 m but a 90,000 km2 cell averages the
depth away. Flattening them to zero would undo the preservation the pipeline is
built around.

Against the superseded files, which came from the 510k-region map PNGs: land
fraction 0.4284 against 0.4137, topography maximum 5,101 m against 6,001 m, and
119 below-sea-level land cells where there were none.

## Precision stays at 8 bytes

ExoPlaSim's `precision` argument sets `-fdefault-real-8` or leaves the compiler
at 4-byte reals. Its own documentation says 4-byte "may run slightly faster, but
possibly at the cost of reduced stability". We keep 8, for reasons specific to
what this project measures rather than a general preference.

The convergence criteria are the binding constraint. They require mean absolute
TOA and surface imbalances below 0.5 W m-2 and trends below 0.05 W m-2 per orbit.
Annual means are accumulated by naive summation into plain `real` arrays in
`plasimmod.f90`, so those accumulators take the compiled precision. Summing of
order 10^3 to 10^4 samples of a ~300 W m-2 flux at single precision carries a
worst-case accumulation error around 0.03 to 0.2 W m-2, which is the same order
as the trend threshold and a large fraction of the absolute one. The 0.85-S-Earth
endpoint already missed its threshold by 0.004 W m-2; that is the margin we work
in, and single precision does not resolve it reliably.

The immediate experiment makes it worse. Bracketing the albedo endmembers is a
comparison between two separate integrations, and numerical noise does not cancel
in a difference taken across independent runs.

Single precision would be a reasonable tool for exploratory scanning, where the
question is a several-kelvin difference rather than a fraction of a watt. It is
the wrong tool for the equilibration diagnostics that gate decisions here.

### The trap, which applies whichever precision is chosen

ExoPlaSim names the binary `most_plasim_t<res>_l<layers>_p<ncpus>.x`. The `p` is
the CPU count, not the precision, and nothing in the path records what it was
compiled with. `Model.__init__` reuses whatever exists at that path unless
`recompile=True`. So changing `model.precision_bytes` alone changes nothing, and
a binary built once at 4 bytes silently serves every later run that believes it
is at 8.

Two guards. `run_id` now carries `r8` or `r4`, so runs at different precision
cannot share a directory. And every run manifest records the executable's SHA-256,
so a swapped binary is visible in provenance even though the filename is not.
Neither prevents the mismatch; they make it auditable. Changing precision means
deleting the binary or passing `recompile=True`, and restart files are not
portable across the change because the Fortran record layout differs.

## Cores, and the measured cost of an orbit

16 MPI ranks, one per physical core on a Ryzen 9 7950X3D. Not 32: those are SMT
siblings and MPI ranks on paired threads contend for the same FPU and cache,
which does not help a compute-bound spectral model. `NLAT` must divide by the
rank count, so at T42 the choices are 1, 2, 4, 8, 16, 32 and 64, and 16 leaves
four latitudes per rank.

Measured, not assumed. At 8 ranks the completed sweep averaged 1.86 to 2.70
minutes per orbit across the three flux cases. The first 16-rank run took 1 minute
44 seconds wall clock for the whole prepare, compile check, staging, one orbit and
postprocessing, at 1411% CPU. So a 50-orbit spin-up is about two hours, not the
days this note previously implied.

Changing the rank count changes the executable name, so ExoPlaSim rebuilds
automatically; the T42 build takes about 15 seconds.

## First orbit on the new world

Prepared and run end to end against the regenerated geography, the mesh-derived
boundary conditions and lithology albedo. Against the old world's first orbit,
which is the only fair comparison since neither is equilibrated:

| | new | old |
| --- | ---: | ---: |
| surface temperature | 263.77 K | 265.29 K |
| 2 m air temperature | 263.45 K | 264.98 K |
| precipitation | 0.970 mm/day | 1.003 mm/day |
| TOA net radiation | 9.43 W/m2 | 16.24 W/m2 |
| surface heat flux | 23.14 W/m2 | 28.80 W/m2 |
| planetary sea ice | 0.1035 | 0.0933 |
| planetary albedo | 0.2551 | not recorded |

Cooler, brighter and closer to balance, which is the expected direction for a
substrate albedo of 0.315 replacing a uniform 0.22. Internally consistent too:
evaporation of 0.987 mm/day implies 28.6 W/m2 of latent heat against a reported
`hfls` of -29.5 W/m2, and precipitation nearly matches evaporation as it should
over a full orbit.

Note that `pr` and `evap` are in m s-1, not mm/day. The conversion is 86400 x 1000.

## Glaciers are enabled

`NGLACIER` was 0, ExoPlaSim's default, through the first prepare. It is now 1,
with `GLACELIM` 2.0 m water equivalent and `ICESHEETH` -1.

The reason is not the albedo, which snow already supplies: `dalb` is blended
toward the snow value continuously with depth whether or not the glacier module
runs. It is the orography. `glaciermod` keeps a lithographic and a glacier
orography and sets the surface geopotential to their sum, so an ice sheet raises
the ground it sits on and grows into its own cold. That is the feedback that
decides whether a cold branch runs away, and the experiment in progress is
specifically about whether this world has more than one stable state. Running it
with the ice-sheet feedback switched off would answer a different question.

Snow does not accumulate without limit in the absence of the module. `newsnow.f90`
caps `dsnowz` at 3000 m water equivalent, which is ice-sheet scale, so the
difference is not runaway mass but the missing elevation response.

The module is conservative. `ICESHEETH` -1 places no initial ice, so a glacier
appears only where snow survives a full model year, and orbit 0 is bit-comparable
with the module on or off: 263.79 K against 263.77 K, glacier fraction exactly
zero, mean land snow 9 mm. Whether it ever fires is an outcome, not an
assumption. After one orbit the deepest land snow is 0.730 m, on a tropical
summit at 9.8 S rather than at a pole, with the land mean growing 16 mm per
orbit, so the 2 m threshold is reachable within a 50-orbit spin-up.

Two limitations to carry. The module does not move ice, so continental ice-sheet
extent is underestimated, as its own documentation says. And `groundoro == 0` is
treated as a sentinel in the initialisation, but only on the `ICESHEETH >= 0`
branch, so the 119 land cells we place below sea level are untouched.

Glacier albedo is spectrally split, 0.745 below 0.75 um and 0.431 above, which is
a further reason the two-band configuration matters under a 4965 K spectrum.

`run_id` gains a `_glac` marker when the module is on, so a glaciers-on and
glaciers-off pair cannot share a run directory.

## Full audit of ExoPlaSim configuration

`configure()` takes 102 parameters. We set 26; the rest sit at defaults. Audited
against the source rather than the documentation, because several defaults are
wrong for this world in ways nothing announces.

### Changed as a result

**Stellar spectrum.** ExoPlaSim ships `stellarspectra/k2.dat` and this star is
K2.5V, so the 4965 K blackbody was replaced by the measured spectrum.
`radmod.f90:801` gives `nstarfile` precedence over `nstartemp`, and the spectrum
sets the weighting behind every snow, ice, glacier and surface albedo. It is not
a refinement. Overall snow albedo falls from 0.562 to 0.395, because a real K
dwarf puts far more flux above 0.75 um than a blackbody of the same effective
temperature does, and snow is dark there. The below-0.75 um band barely moves,
0.745 to 0.754; the above-0.75 um band falls 0.431 to 0.347. At orbit 0 the
planet is 2.6 K warmer, absorbs 29.1 against 9.2 W m-2, and carries 0.075 against
0.104 sea-ice fraction. The blackbody was biasing this world cold and icy in
exactly the quantity the current experiment is designed to measure.

**Forest fraction.** `dforest` was at its uniform 0.5 default while the
background albedo asserted bare rock. It is not decorative: `landmod.f90:383-392`
blends snow albedo between forested and unforested endpoints by `dforest`, so the
default was darkening snow as though half the land were canopy. It is now written
as code 212 by `build_surface_albedo.py`, tracking the albedo mode: 0 for
`lithology`, 0.5 off-evaporite for `vegetated`.

### Corrected reasoning on glaciers

Enabling the glacier module is still right, but not for the reason recorded when
it was switched on. `landmod.f90:839-845` hard-clips snow at `dsmax`, whose
default is **5.0 m water equivalent**, and line 424 pins a glacier cell at exactly
that. The glacier orography contribution is therefore capped at about five
metres and is inert. PlaSim cannot grow an ice sheet, and raising `maxsnow` would
not fix it, because without ice flow the accumulation zone would thicken into a
tower rather than spread into a sheet.

What the module does contribute is albedo persistence. A glacier cell uses
`albgmin`/`albgmax` instead of the snow curve, so it holds a minimum albedo of
0.745 below 0.75 um where snow decays to 0.501 as it approaches melting. That
hysteresis is a genuine bistability mechanism, and it is what the setting buys.

State the limitation plainly when interpreting: this experiment tests albedo
bistability including glacier albedo hysteresis and excluding ice-sheet growth.

### Left at defaults, deliberately

- `wetsoil` False. It adjusts land albedo by soil moisture and is tuned to Earth
  observations, so it is the wrong tool under a K dwarf, and it would fight the
  lithology albedo we supply.
- `vegetation` False, so `NVEG=0` and SimBA is off. Vegetation is LPJ-GUESS's job
  downstream. Note `NCVEG` in `vegmod_namelist` is the growth accelerator, not the
  on/off switch; reading it as the switch is an easy mistake.
- `snowicealbedo` None, which is what lets the spectrum compute the albedos above.
  Setting it would override them with a single number.
- `co2weathering`, `evolveco2` False. CO2 is fixed at 450 ppm by design; a
  silicate-weathering feedback is a separate experiment.
- `aquaplanet`, `desertplanet`, `drycore`, `aerosol`, `synchronous` all False.
- `modeltop` None, giving `PTOP` 5000 Pa, appropriate for a 1 bar atmosphere.

### Left at defaults, but worth knowing

- `soilwatercap` None, so `WSMAX` is Earth's 0.5 m uniformly. This is the bucket
  capacity the lake work will want to set per cell from basin hypsometry.
- `cpsoil` None and `soildepth` 1.0, so land heat capacity is uniform at
  2.4e6 J m-3 K-1. The export carries per-class rock densities, so this could be
  made lithology-dependent later; it is second-order next to albedo.
- `stormclim` False. Turning it on yields tropical-cyclone diagnostics, which is
  a worldbuilding output rather than a physics correction.
- Diffusion and physics-filter coefficients are all at resolution-dependent
  defaults. We enable the filter but tune nothing.

### Two traps found in the process

`starfile` and `starfilehr` are `character(len=80)` in `radmod.f90`, and
`configure()` absolutises whatever path it is handed. A venv path is easily
longer: ours was 90 characters. Fortran truncates at 80, emits a namelist
*warning* rather than an error, then dies with an end-of-file inside `readdat`
naming a path that does not exist. `stage_stellar_spectrum` copies the two files
into the run directory and rewrites the namelist to bare names, which the model
resolves against its own working directory.

T63 additionally requires `pyfft991` to be compiled by hand, per the ExoPlaSim
tutorial. T42 and T85 do not.

### Model biases to carry into interpretation

From the ExoPlaSim tutorial, and worth repeating wherever results are read: no
deep-ocean circulation, so high latitudes run somewhat too cold, which matters
directly for the glaciation question; monsoons come out weak; there is a bias
toward Mediterranean patterns; and small islands and peninsulas come out drier
than they should.

## Albedo-endmember bracket, T21, six cases

Two land-surface endmembers at three fluxes, 50 orbits each, glaciers enabled,
K2 spectrum, mesh-derived boundary conditions.

| flux | mode | ts K | pr mm/day | TOA W/m2 | sea ice | albedo |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 0.90 | bare rock | 273.55 | 1.815 | -0.63 | 11.10% | 0.278 |
| 0.90 | vegetated | 280.68 | 2.207 | -0.44 | 4.75% | 0.189 |
| 0.95 | bare rock | 285.67 | 2.487 | -0.33 | 2.09% | 0.212 |
| 0.95 | vegetated | 290.05 | 2.803 | -0.33 | 0.39% | 0.155 |
| 1.00 | bare rock | 294.41 | 3.157 | -0.31 | 0.04% | 0.188 |
| 1.00 | vegetated | 298.12 | 3.455 | -0.12 | 0.00% | 0.143 |

### The result that matters

The land-surface assumption is worth 7.14 K at 0.90, 4.38 K at 0.95 and 3.71 K at
1.00. It narrows as the sea ice that amplifies it disappears; by 1.00 both
branches are ice-free and the residual 3.71 K is land albedo alone.

Interpolating to equilibrium, the flux that puts this world in the 290 to 293 K
design range is **0.9516 to 0.9697 S-Earth if vegetated, and 0.9765 to 0.9936 if
bare rock**. Those windows do not overlap, by 0.0068 in flux.

The reason is structural rather than a numerical accident: the endmember spread
near the target, 3.7 to 4.4 K, is wider than the 3 K target range itself.
Whenever that holds, no single flux can place both branches inside the range.

So the orbit cannot be chosen independently of the biosphere. Either the flux is
picked and the world's habitability is contingent on its vegetation state, which
is a defensible and rather appealing worldbuilding fact, or the vegetation state
is resolved first and the flux follows from it. What cannot be done is to pick a
flux that is robust to the question.

### What the bracket does not show

It does not show bistability. Albedo is a fixed input here because vegetation is
not interactive, so the six runs are separately forced problems with separately
correct answers. A wide spread is sensitivity. Bistability would require one set
of boundary conditions admitting two stable states, and establishing it needs the
vegetation feedback actually closed, or at minimum a run with `vegetation=2`.

### Caveats

Every case still carries a negative TOA balance, from -0.63 to -0.12 W/m2, so
each reported temperature is an upper bound on its own equilibrium. The
equilibrium estimates above apply the measured 0.94 K per W/m2 to that residual.
Only three of the six pass the full convergence test unaided; the rest miss on
mean TOA alone. They should be extended before any of this feeds a carve verdict.

Glaciers formed in none of the six. Mean snow never exceeded 7 mm, so the module
stayed inert throughout and this bracket measures snow and sea-ice albedo
feedback only. Enabling it cost nothing and bought nothing here.

Flux sensitivity falls as ice runs out: +12.12 K across 0.90 to 0.95 on the bare
rock branch, then +8.74 K across 0.95 to 1.00, with sea ice going 11.1% to 2.1%
to 0.04%. The old sweep's +11.2 K per 0.05 step sits between the two, so the
model's flux response is unchanged; only the offset moved.

## Baseline chosen: vegetated at 0.96 S-Earth

From the bracket, the vegetated branch reaches the 290 to 293 K target between
0.9516 and 0.9697 S-Earth. 0.96 sits mid-range and gives about **291.4 K**.

The choice is a package, not two independent settings. The bare-rock window,
0.9765 to 0.9936, does not overlap the vegetated one, so 0.96 is habitable
*because* the world is vegetated. A world at this orbit that lost its biosphere
would fall to roughly 286 K, and one that gained a biosphere at 0.98 would rise
past 297 K. That coupling is a property of the world worth keeping in view rather
than an artefact to be tidied away.

### The stellar cycle at 0.91 to 1.01

Centred on the baseline, 10.4% peak-to-peak, against the earlier experiment's
11.1%. Still an intentionally active star: this amplitude implies unusually
extensive and variable spot coverage for a quiet 6 to 8 Gyr K dwarf, and should
be presented as a deliberate choice rather than as typical behaviour.

Static endpoints on the vegetated branch are about 282.2 K at 0.91 and 299.7 K at
1.01, a 17.5 K span. The earlier cycle work measured the 50 m slab damping a
33.3 K static span to a 7.96 K seasonally adjusted response, a ratio of 0.239. If
that carries over, the prediction here is about **4.2 K peak-to-peak, roughly
289.3 to 293.5 K**, which keeps the world inside the habitable band across the
whole cycle. That is the reason to centre at 0.96 rather than at an endpoint.

The 1.01 endpoint is extrapolated one step beyond the 1.00 run, so it is the
least supported number in the set. The prediction should be treated as a
hypothesis the cycle runs will test, not as a result.

Cycle amplitudes live in `config/planet.yaml` rather than as literals in
`run_stellar_cycle.py`. They were hardcoded while the baseline was 0.90; leaving
them there would have applied the old 0.85-0.95 range about the new centre
without complaint.

**Superseded 2026-08-16** by the two-component design below. The single-sinusoid
`stellar_cycle.cases` representation carried its own minimum and maximum, so it
also carried its own mean, which could disagree with `orbit.baseline_flux_earth`;
the script had to check that they matched. The component representation states
amplitudes about the baseline and has no second mean to disagree with.

### Two components, 11 and 57 Earth years

*Decided 2026-08-16, against a 0.945 baseline.*

One sinusoid was a placeholder. Epsilon Eridani, the reference for what a K2.5V
of this activity level can do, genuinely carries both a short and a long cycle at
once, and the two do different jobs on this world:

| component | period | peak-to-peak | damping | response |
| --- | --- | --- | --- | --- |
| medium | 11 Earth yr | 2.5% | 0.83 | ~3.1 K, climatic |
| long | 57 Earth yr | 3.5% | 0.99 | ~5.2 K, geomorphic |

The periods are **not** epsilon Eri's measured values, and copying them would
have been a category error: that star bounds the plausible *amplitude* envelope
(4.3% median, 2.4 to 7.2% at 95%, bolometric), not the periods this world's
narrative wants. The 6% total is inside that envelope and above its median.

The damping is analytic and Planck-only, from a slab thermal timescale of
tau = C/(4 sigma T^3) = 1.17 Earth years, giving 1/sqrt(1 + (2 pi tau / P)^2).
Real damping includes the ice-albedo feedback and will be larger in the cold
phase. Measuring it is what the first cycle run is for.

**The ratio is deliberately non-commensurate.** At 57/11 = 5.18 the phase
relationship returns only after about 314 Earth years, so successive grand minima
differ in depth rather than repeating identically. That is the point rather than
a detail: glaciers advance further in the deeper minima and leave moraines at
different distances, so the landscape records which past minima were severe. An
integer ratio erases that record.

Note that the aligned envelope is approached rather than attained. The
cycle-mean temperature sits very slightly *above* the static mean at the same
flux, because T(f) is measured convex over this range; see the curvature section
below. The offset is under 0.05 K.

### The cold-regime slope, measured 2026-08-16

Two 60-orbit T42 runs on `precarve-zoned-g1281` bracketing the 0.945 baseline.
Both runs have since been archived and their output deleted; the ids below
resolve under `archive/runs/`, which keeps each run's manifest, convergence
assessment and climate series. The measurement stands -- a flux slope turns on
sea ice and the Planck response, not on which terrain produced it.

| flux | run | mean T | fitted asymptote | sea ice |
| --- | --- | --- | --- | --- |
| 0.9125 | `run_5aed450f3972` | 282.489 K | 282.491 K | 6.174% |
| 0.968 | `run_1dbb75d05aca` | 293.741 K | 293.663 K | 0.206% |

**201 K per unit flux ratio**, or 2.01 K per 0.01 — from the fitted asymptotes.
Run means give 203 and drift-implied endpoints 206, so call it 201 to 206.

**Neither run formally converged, and both are one criterion short.** The cold
run fails `abs_mean_toa_lt_0.5_w_m2` at -0.653 W/m2 and is still cooling at
-0.028 K/orbit, with a drift-implied remaining offset of -0.269 K. The warm run
fails `extrapolated_offset_lt_0.15_k`, though its remaining offset is only
-0.078 K; its relaxation fit is degenerate, which is what actually trips the
criterion. Both are *quasi-equilibrated* in the sense this project already uses
for the 0.85-flux case, and the slope is quoted with that caveat rather than
without it.

**This supersedes the 290.3 K figure for flux 0.945**, and it supersedes it by
exactly the mechanism this file already warns about. 290.3 K came from
extrapolating down from the 0.968 run using a sensitivity near 149 K per unit
flux, measured where the world is nearly ice-free. The 0.9125-to-0.968 interval
crosses the ice transition -- sea ice grows thirtyfold across it -- so the true
interval sensitivity is 35% larger. The chord now gives **289.0 K at 0.945**.

The chord was expected to *underestimate* the interior point, on the reasoning
that the ice-albedo feedback makes dT/df larger at low flux and T(f) therefore
concave. **Measured, that is wrong** -- see the curvature section below. The
chord's 289.0 K turned out to be accurate to 0.03 K, not a lower bound.

The consequence for the cycle is that its predicted response was understated:

| component | peak-to-peak | static span | damped (analytic) |
| --- | --- | --- | --- |
| medium, 11 yr | 2.5% | 5.03 K | 4.18 K |
| long, 57 yr | 3.5% | 7.05 K | 6.98 K |
| aligned envelope | 6.0% | 12.08 K | **11.15 K** |

against the 8.3 K written down when the components were chosen. The amplitudes
are still the decision; only their consequence moved. Note also that the aligned
envelope is an extreme reached only near the 314-year recurrence: the root-sum-
square of the two damped semi-amplitudes is about 2.9 K, so a typical excursion
is roughly +/-2.9 K rather than +/-5.6 K.

### The cycle is grey, and that is measured rather than assumed

The patched `radmod.f90` scales `gsol0` and leaves the Lacis-Hansen band split
alone, so the cycle is a bolometric multiplier with a fixed 4965 K spectral
partition. A real spot-driven cycle is not grey: the star reddens at minimum.

Sized before building anything. For a photosphere at 4965 K and spots at 4000 K,
(T_spot/T_phot)^4 = 0.4213, so a 6% peak-to-peak bolometric swing needs the spot
covering fraction to move by 10.4 percentage points. Band-1 fractions (below the
0.75 um split) are 0.4127 for the photosphere and 0.2690 for the spots, so about
f = 0.10 `zsolar1` runs from 0.4097 at flux maximum to 0.4026 at minimum. The
cold phase is redder by **0.0071**.

Both endmembers are blackbodies, which is the right model for a spot but not for
a photosphere: the BT-Settl spectrum this world actually declares puts 0.384383
below the split rather than 0.4127. The argument here is a difference between
two phases and survives the offset unchanged; the absolute level does not, and
`lib/stellar.py` is where it comes from.

Carried through to what it does, via the effective ice albedo shift of +0.0025:

| sea ice | planetary albedo | temperature |
| --- | --- | --- |
| 0.21% | +0.000005 | 0.0011 K |
| 3.60% | +0.000089 | 0.0183 K |
| 10.0% | +0.000248 | 0.0509 K |

Under 0.02 K at this world's ice cover. **Runtime spectral interpolation is not
warranted**, and neither is regenerating the spectrum per phase.

This is the k2-to-k25v correction repeating exactly: fixing a genuinely wrong
stellar spectrum moved absorbed shortwave by 0.04 W/m2 here, because the spectral
partition acts on snow and ice and there is almost none of either. The general
form is in `notes/failure-modes.md` -- a correction's size depends on how much of
the surface it acts on.

**The condition under which this stops holding is stated so it can be checked**:
it is negligible because sea ice is a few percent. On a cold branch at 10% ice it
is 0.05 K, still small; at the sort of ice cover a much colder world carries it
would not be. Anything that moves this world onto a substantially icier state
should re-derive the table above rather than inheriting the conclusion.

### Resolution still to settle

The bracket ran at T21, which was the right call for a question about a several
kelvin spread. The production baseline feeds hydrography, the carve verdict and
eventually biomes, and those want regional detail: the completed sweep was T42,
the coupling matrices in `hydrography/data/` are built for T42 and T85, and none
exists for T21. Recommend T42 for the baseline unless there is a reason not to.

### A stale albedo label, corrected

The vegetated endmember was described as 0.197 land-mean albedo in several
places. The correct figure is **0.223**. The 0.197 came from an early grid-based
calculation in which an "evaporite" grid cell carried an area-averaged albedo
rather than the class value of 0.50, so the bright fraction was understated.

This is a labelling error, not a modelling one. The mesh-based path was already
in place before the T21 bracket ran, so every measured temperature in that
bracket used 0.223 and all six results stand unchanged. Only the number quoted
alongside them was wrong, which also means the endmember separation is 0.091
rather than 0.117 in land albedo, and the forcing gap is nearer 12 to 16 W/m2
than 15 to 19.

## The 0.96 baseline climate

First described climate for the current world. Converged on all six criteria;
five-orbit climatology with 32 snapshots per orbit.

| | |
| --- | ---: |
| surface temperature | 292.97 K |
| 2 m air temperature | 292.23 K |
| precipitation | 2.949 mm/day |
| TOA net radiation | -0.351 W/m2 |
| surface heat flux | +0.090 W/m2 |
| planetary sea ice | 0.133% |
| land area fraction | 42.84% |

Koppen classes by land area, largest first: BWh 19.4%, Dfb 16.2%, BSh 11.5%,
Dfa 6.8%, Af 6.7%, Cfa 5.7%, Aw 5.5%, Dfc 5.1%, Cfb 4.8%, Csa 4.4%.

Broad biome interpretation: continental temperate and mixed forest 28.3%, desert
19.4%, steppe and semidesert 12.0%, temperate and subtropical forest 11.5%,
seasonal tropical woodland and savanna 9.0%, everwet tropical forest 6.7%, boreal
5.9%, Mediterranean woodland 4.7%, monsoon tropical forest 2.2%.

**31.4% of land is desert or steppe.** That is a dry world, and it is the same
story the hydrography tells independently: 76% of land drained to a closed basin
before carving, and evaporite covered 20.8% of it. Three separate calculations,
Koppen from the climate, drainage from the terrain, and lithology from tectonic
history, agreeing that this planet's interiors do not reach the sea.

### What this climatology is not

It is the climate of a superseded surface. It was computed on the pre-carve
terrain, so it still carries the 3,629 preserved basins including the 1,522 the
verdict has since carved. And it used the exported evaporite albedo of 0.50,
which is now a declared override at 0.40 pending the zoned split upstream. The
two corrections move in opposite directions, roughly +5.0 and -2.2 W/m2 on the
pre-carve fractions, so the next baseline is not obviously warmer or cooler and
the estimate is not worth trusting far enough to skip measuring it.

Quote these numbers as the 0.96 pre-carve baseline, never as Vesper's climate.

## Where the albedo uncertainty actually lives

The evaporite constant was the largest single lever in this planet's energy
balance. Orogen's crust/fill split did not remove it, it moved it, and it is
worth knowing exactly where it went.

Per 0.10 of albedo error, on the exports as delivered:

| class | share of land | pre-carve | carved |
| --- | ---: | ---: | ---: |
| `playa_clastic` | 18.80% | 2.48 W/m2 (2.34 K) | 1.50 W/m2 (1.41 K) |
| `evaporite` (crust) | 2.11% | 0.28 W/m2 (0.26 K) | 0.18 W/m2 (0.17 K) |

Nine to one in favour of the fill. Almost the entire remaining sensitivity sits
on the 0.30 assigned to playa clastics, which was chosen the same way the old
0.50 was.

That number is also a mixture: "Playa mud / alluvial fan fill" covers light clay
playa at roughly 0.30-0.35 and desert-varnished fan gravel and pavement at
roughly 0.15-0.25. A plausible mix spans 0.25 to 0.33, worth 1.99 W/m2 or 1.87 K
pre-carve, which is 40% of what the whole split was worth.

It is recorded and not overridden. 0.30 is defensible as a central value, and
picking a different one to make a temperature land would be exactly the error the
split was meant to fix, one level down.

### A degeneracy neither side can break

Crust area and crust albedo are degenerate in the effective value: from the
climate side, 10.6% of basin floor at 0.50 and 20% at 0.40 are the same planet.
Nothing in the terrain geometry separates them. Only an independent measurement
of crust extent would, and that is flooding frequency, which is a water balance
and therefore ours, and our lake solver is the least-validated thing in this
pipeline. So neither side can close it, and both have declined to invent a
constraint. Left open and labelled.

### Where this stops being worth doing

Each split replaces one guessed number with two better-founded ones, but it does
not recurse usefully forever. At some depth the dominant uncertainty stops being
"which facies is this" and becomes "what does bare rock of this type reflect
under a 4965 K spectrum", and neither the terrain model nor the climate model has
anything to say about that. One more level is probably worth it. The level after
that is not.

## The carved-zoned baseline, and a calibration error worth more than the result

Ran at 0.96 S-Earth on `carved-zoned`, vegetated, glaciers on, K2 spectrum.
Converged on all six criteria after 70 orbits.

| | baseline (`precarve-unzoned`) | carved-zoned |
| --- | ---: | ---: |
| land-mean albedo prescribed | 0.2240 | 0.1722 |
| surface temperature | 292.93 K | **295.18 K** |
| precipitation | 2.954 | 3.128 mm/day |
| planetary albedo | 0.1519 | 0.1276 |
| sea ice | 0.133% | 0.018% |
| mean TOA | -0.496 | -0.464 W/m2 |

**Predicted 299.4 K. Measured 295.18 K.** The albedo arithmetic was right to
0.0004 in land-mean albedo; the climate response was wrong by a factor of three.

### The sensitivity used for every extrapolation was measured in the wrong regime

0.94 K per W/m2 came from the T21 bracket's 0.90 to 0.95 step, where sea ice fell
from 11.1% to 2.1%. That step was amplified by the ice-albedo feedback. This run
sits at 0.018% sea ice with nothing left to lose, so the same forcing buys far
less warming:

    prescribed forcing        +6.80 W/m2
    measured warming          +2.25 K
    implied sensitivity        0.331 K per W/m2

    Planck response alone      0.312 K per W/m2

The realised response is within 6% of the bare Planck value, which says the net
non-Planck feedback in this state is close to zero: no ice left, and whatever
water vapour adds is offset by lapse rate and cloud.

So there are two sensitivities on this planet and they differ threefold. Near the
ice transition it is about 0.94 K per W/m2; in the ice-free state above roughly
293 K it is about 0.33. Every extrapolation in this project that used the first
number outside the ice-transition regime overestimated by about three times,
including the +7.06 K forecast for this run and the 0.925 flux estimate derived
from it.

Record which regime a sensitivity was measured in whenever one is quoted. The
project has repeatedly warned against extrapolating across the ice-albedo
transition and then did exactly that by carrying a coefficient measured inside it
into a state without it.

### What it means for the flux

At 0.33 K per W/m2 and an ice-free planetary albedo of 0.128, a 0.05 change in
stellar flux is worth about 14.8 W/m2 absorbed and therefore about 4.9 K, or 98 K
per unit flux, against the 242 K per unit flux the bracket suggested. Returning
from 295.18 K to the middle of the 290-293 K band needs roughly -0.037 in flux,
so about 0.92.

That number is itself an extrapolation, now with a coefficient measured in the
right regime but still only one point. It should be measured, not trusted.

---

## Superseded numbers in this file, and what replaced them

This file records decisions with the evidence available when they were taken, so
the numbers below are left where they are rather than rewritten. What follows is
the index of which ones have since moved and why.

**Land albedo 0.223 (vegetated) and 0.314 (bare).** Both were measured on the
pre-carve terrain. On `carved-zoned-v4` they are **0.179 and 0.276**. Two changes
moved them: carving removed bright closed-basin fill, and the v4 lithology fix
restored fill that the cover chain had been overwriting. Note the two move by
different amounts -- the v4 fix was +0.0034 bare against +0.0077 vegetated --
because vegetation paints everything that can carry a canopy at one value, so it
masks bare-rock variation but not the barren classes. When a lithology change is
confined to basin fill, expect the vegetated figure to move about twice as far.

**Evaporite at 20.8% of land.** That was one class on the pre-carve terrain. It
is now split, and closed-basin fill is 16.5% of land on `carved-zoned-v4`: 1.9%
salt crust, 14.5% playa clastics.

**Playa clastics at 0.30 albedo, and the 0.25-0.33 plausible mix around it.**
The generator moved that class to **0.19** and cites Post et al. (2000) for it:
52 pyranometer measurements over 0.3-2.8 um with a mean of 0.189, which is the
same quantity measured in the field rather than a mixture argued from facies.
See `vendor/orogen/js/lithology.js`.

The sensitivity section above is therefore an underestimate now, in two ways at
once. The class is a larger share of land on `precarve-craton` than it was when
that was written, and it sits below the bottom of the 0.25-0.33 range the
mixture argument produced, so the lever is both longer and pulled further. It
remains the largest single lever on this planet's energy balance.

Note the two groundings disagree and both are in `references/INDEX.md`:
Henderson-Sellers and Wilson (1983) put salt playas and light sand deserts at
0.28-0.44, against Post's all-soils mean of 0.189. The generator took the field
pyranometer measurement, which is the right *kind* of number; whether US
agricultural soils are the right *population* for playa mud and desert-varnished
fan gravel is not settled here.

**The albedo bracket's flux windows, 0.952-0.970 vegetated and 0.977-0.994
bare.** Measured at T21 on pre-carve terrain. The vegetated window has since been
measured directly on the current terrain from three converged T42 points, giving
192.2 K per unit flux ratio and a baseline of **0.945**. The bracket's *implied*
sensitivity of about 167 K per unit flux was right; a later 0.331 K/W/m2 figure,
measured across an albedo step in a nearly ice-free state, was not, and using it
predicted 291.9 K for a run that converged at 287.47 K.

**The measured stellar spectrum.** Every run before the k25v baseline re-run used
`k2.dat`, which is the star K2-18, an M2.5V at about 3450 K, not a K dwarf. The
correction changed absorbed shortwave by 0.04 W/m2 on this world, because it acts
on snow and ice and there is almost none at 291 K. It is not null on the cold
branch. See `stellar-spectrum-audit.md`.

**Any carve verdict before the longitude fix.** The coupling matrix numbers its
columns -180 to 180 and a climatology numbers its own 0 to 360, so every basin
integrated its antipode's climate. See `../../hydrography/README.md`.
