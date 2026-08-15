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
Because this planet's orbit is 189.6 Earth days, applying them to precipitation
per local orbit would make identical physical rain rates appear about twice as
dry. The diagnostic therefore annualizes modeled precipitation rates to
365.2425 days while retaining the simulated seasonal temperature and rainfall
cycle. The resulting classes and biome groupings are heuristic ecological
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

The spread is not incidental to this world. Evaporite covers 20.8% of the land,
second behind schist at 22.4% and ahead of intracratonic clastics at 18.9%, and
evaporite is the brightest class in the table at 0.50. (The manifest's
`compositionLand` says 18.6%, because it measures against `land_mask` and so
drops the dry sub-sea-level basin floors where playa fill concentrates. The
figures here are computed from `surface_rock` against `surface_class`.) That is a direct consequence of the
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
bare because nothing grows on a salt pan, gives 0.197. That is 15 to 19 W/m2 of
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
  evaporite. Land mean 0.197. Physical, warm, and the other half of the bracket.
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
