# Ecological climate forcing audit

**Recorded:** 2026-08-21
**Scope:** the path from ExoPlaSim state and diagnostics through climatology
products and `build_lpj_driver.py` into the LPJ-GUESS-CNP process consumers.
This was a source and pipeline audit only. No ExoPlaSim or LPJ-GUESS simulation
was run.

## Result

The current 12-bin input is project scaffolding, not a scientific or upstream
model constraint. The `VESPDRV` driver format, `DRIVER_BINS = 12`, the four climate arrays and
their interpolation were all added by this project. ExoPlaSim integrates at a
45-minute step, and LPJ-GUESS already has a daily input path plus a subdaily
temperature/radiation path. There is no reason to preserve the binary or its
bin count in a replacement design.

The scientifically preferred path is to retain one coherent chronological
weather sequence from ExoPlaSim, reduce it at the producer with explicit time
bounds and process-specific operators, and let LPJ consume those intervals
directly. The existing averaged climatology remains useful for maps, climate
diagnostics and equilibrium summaries. It should not be the source of weather
for a nonlinear ecological model.

This audit therefore does not propose a larger monthly file or a more elaborate
monthly weather generator. It proposes a separate ecological forcing artifact.
Its transport format should be selected after the field, time and conservation
contract is fixed. A self-describing, chunked NetCDF-like product is an obvious
candidate, but that is an engineering choice rather than a premise of the
science.

## 1. The twelve bins and binary are entirely artificial

`run_exoplasim.py` asks pyburn to turn each orbit into 12 regular time-bin
means. `build_lpj_driver.py` then refuses any other bin count, converts
precipitation rates to bin totals, and writes temperature, total precipitation,
net surface shortwave and diurnal range to the driver file. `vesperinput.cpp` repeats
the 12-bin check and calls LPJ-GUESS's monthly interpolation functions. Those
functions generate smooth quasi-daily temperature and radiation and spread
each bin's precipitation total smoothly across its days.

None of these constraints comes from the CNP fork. The stock input interface
requires daily temperature, precipitation and radiation. Its CF input module
also has paths for daily relative humidity, wind and temperature extrema, and
the framework supports subdaily temperature and radiation vectors. The present
adapter proved that a Vesper grid, calendar and soil can run end to end. It did
not establish 12 seasonal bins as an adequate forcing representation.

The replacement should not version the next `VESPDRV` magic around a new fixed
number of bins. It should carry explicit interval bounds in seconds, orbital
position and local solar phase. A consumer can then aggregate or subdivide according to the
process it is actually executing.

## 2. Chronology is discarded before LPJ-GUESS sees it

The default `baseline_climatology` averages corresponding bins across all
selected post-equilibrium orbits. The pipeline gives that averaged product to
`lpj_driver`. Interannual weather variability, the ordering of wet and dry
spells, co-occurring heat and drought, and the relationship among wind,
humidity, rain and radiation are gone at that point.

`build_climatology.py --per-year` preserves orbit-to-orbit differences, but
each output is still the same 12-bin reduction. Passing several such files to
the current driver proves that the adapter can cycle multiple forcing years. It
does not preserve daily weather within any of them. The README's earlier
description of several files as "interannual variability" was consequently too
strong.

This matters even for a stationary baseline. Sitch et al. explicitly required
interannually varying spin-up climate because fire and other water-limited
responses occur preferentially in dry years. Gerten et al. added stochastic
wet-day precipitation because interception, soil water, transpiration and
runoff do not respond to a monthly total as though its timing were irrelevant.
The active code has daily phenology, GDD/chilling, interception, snow,
infiltration, soil water, stomatal stress, respiration and fire operators. A
mean seasonal cycle is not a neutral input to those nonlinear processes.

The pipeline needs a distinct `ecological_forcing` product made from consecutive
accepted model orbits before climatological averaging. The averaged
climatology and the chronological forcing can be derived from the same run, but
neither should masquerade as the other.

## 3. Existing clean and snapshot output are not daily forcing

The clean-I/O regular stream contains about 182 records per orbit, but these are
instantaneous atmospheric samples written once per 24 Earth hours by default.
They are not contiguous 24-hour means or precipitation totals. Moreover,
Vesper's rotation is 30 hours, exactly 40 model steps against the 32-step output
interval. The samples therefore visit only five local solar phases. CLIM-11
measured the resulting fixed alias in an annual energy diagnostic.

The 32 seasonal snapshots are also instantaneous samples. Averaging snapshots
across orbits is useful for a phase composite, not for reconstructing the order
of events. Increasing pyburn's output bin count cannot recover the intervals
between either set of samples.

ExoPlaSim already accumulates precipitation, fluxes and extrema internally at
every timestep. A compact ecological stream can use the same live fields and
write only the required surface diagnostics at exact boundaries. This avoids
both bad alternatives: accepting 12 bins because they are small, or writing the
model's entire three-dimensional raw payload every timestep. The latter is
known to be about 15 GB per T42 orbit even when the postprocessed field list is
only three winds, because the current raw path has no field selection.

The ecological stream must be restart-exact. Partial interval accumulators,
extrema and event state have to be serialized, or a model-call boundary will
create a false weather boundary of the same class as the accumulator defects
already documented for the general output path.

## 4. Cadence follows processes, not a universal bin count

The first implementation boundary can remain LPJ-GUESS's existing 24-hour
ecological step. That is an explicit model-form boundary, not a claim that a
Vesper solar day is 24 hours. Each interval must be recorded in absolute
seconds. Since 24 hours is 0.8 rotations, its local solar phase advances each
step and must remain visible.

Some processes need structure within that interval. LPJ-GUESS already loops
over subdaily temperature and radiation for canopy exchange. PCAR-1 requires
that path to integrate light and dark periods on the 30-hour rotation. BVOC
leaf temperature, fire weather and a convection-derived lightning diagnostic
also require synchronous subdaily state. Those fields can share timestamps but
need not force every daily hydrology and biogeochemistry operator onto the
45-minute GCM timestep.

The required producer-side reductions are:

| cadence | quantities | reduction |
| --- | --- | --- |
| every forcing interval | start/end time, duration, orbital position, local solar phase, grid and source identity | exact coordinates, never inferred from record number |
| 24-hour ecological interval | 2 m temperature, pressure, near-surface humidity, wind vector and speed, shortwave and longwave components | duration-weighted means plus named extrema where consumers use them |
| 24-hour ecological interval | large-scale rain, convective rain and snowfall | time integrals, retaining phase and convective partition |
| selected subdaily intervals | temperature, downwelling shortwave or spectral photon supply, humidity and wind | synchronous duration-weighted state with no independent interpolation |
| convection timestep to daily | upstream lightning potential and convective precipitation | synchronous diagnostic integral, event count and relevant maximum as FIRE-1 defines |
| static or slowly varying | cell area, land/rootable fractions, coordinates and surface identity | provenance-bound fields, not repeated as weather |

This is a minimum contract, not a mandate to store every listed field twice.
For example, exact subdaily samples can generate a daily mean and extrema while
the artifact builder verifies their closure.

## 5. Field meaning has to be corrected at the same seam

The current four-field payload hides several semantic choices.

First, `rss` is net downward surface shortwave. It is appropriate for a surface
energy balance, but not automatically the incident photon supply to a canopy:
it already includes reflection by the surface state ExoPlaSim was run with.
ExoPlaSim also retains upward surface shortwave, so downwelling shortwave can be
carried explicitly. The forcing contract should expose incident, upward and net
terms or a conservative equivalent. PCAR-1, PCAR-2 and BIO-25 then own how a
canopy turns incident energy into absorbed photons; the adapter should not hide
that decision by naming net surface energy "insolation".

Second, ExoPlaSim separately diagnoses large-scale precipitation, convective
precipitation and snowfall. The adapter adds the first two and discards the
third, after which LPJ-GUESS repartitions all precipitation from daily mean air
temperature at 0 C. That can turn a mixed or transient precipitation interval
into all rain or all snow. The ecological artifact should retain the modelled
phase, and the LPJ input seam should accept rain and snow explicitly while
closing their sum.

Third, the regular products already contain or can derive pressure, humidity,
wind, net longwave, temperature extrema and solar geometry, but the driver
discards them. BIO-23 owns their use in evapotranspiration and gas partial
pressures; FIRE-2 owns the coherent fire-weather sequence. They should receive
one shared source sequence rather than separate interpolated copies.

Finally, all variables must remain synchronous. Richardson's weather generator
conditioned temperature and solar radiation on wet versus dry state and retained
serial and cross-correlations because independently plausible marginal series
do not make plausible weather. Hempel et al. demonstrate the same failure mode
for independently bias-corrected variables. Vesper does not need Earth bias
correction, and these papers are not a license to apply one. They establish why
the native ExoPlaSim joint chronology should be preserved.

## 6. A forcing replay protocol is still required

LPJ-GUESS spin-up is much longer than the post-equilibrium climate segment that
is economical to archive. Some forcing therefore has to repeat. The transparent
central case is cyclic replay of one consecutive, accepted ExoPlaSim block,
with the full daily and subdaily chronology intact. It must begin and end at
declared orbital boundaries, and the seam between the last and first interval
must be tested rather than assumed harmless.

One repeated block understates uncertainty in rare droughts, wet spells and
compound extremes. The accepted biosphere result should therefore carry a
pre-registered forcing ensemble: different consecutive climate blocks when
available, plus block order or phase sensitivities that do not break seasons.
For a stellar-cycle run, the physical cycle order is not shuffleable and must
be replayed whole.

A Vesper-trained stochastic weather generator is a fallback if the archived
chronology is too short for the needed sensitivity. It would need held-out
validation of wet/dry spell lengths, precipitation amounts, temperature and
radiation extremes, autocorrelation and cross-variable dependence. Richardson's
first-order wet/dry model is evidence that these properties matter, not a
ready-made parameterization for Vesper. Earth GWGEN observations and ISIMIP
bias adjustment are not admissible substitutes for missing alien weather.

## 7. Acceptance is conservation plus information retention

The ecological forcing artifact should fail before LPJ-GUESS if it has missing
or overlapping intervals, a restart discontinuity, inconsistent grids or
source builds, impossible phase partitions, non-finite values, or a duration
that does not tile the declared orbit block.

For extensive quantities, sums from model timesteps through daily forcing and
through the replay block must close. For intensive quantities, duration-weighted
means must reproduce the source. Temperature extrema must bound the means;
rain plus snow and convective plus large-scale partitions must close by their
declared definitions. The acceptance report should also compare distributions,
spell lengths and cross-variable relationships against the source stream so an
artifact can conserve annual means while still destroying the weather.

These checks can be developed against small synthetic fixtures and retained raw
records without launching either model. Production timing and response tests
remain subject to explicit permission.

## Task boundary

EFOR-1 through EFOR-8 record the new producer-to-consumer path. They do not
duplicate the process work already assigned elsewhere:

- BIO-13 owns rainfall-event effects in LPJ hydrology;
- BIO-23 owns pressure, radiation, humidity and wind in EET and gas physics;
- FIRE-1 and FIRE-2 own lightning and fire-weather semantics;
- PCAR-1 owns subdaily photosynthesis on the 30-hour rotation;
- LSHY-4 and LSHY-5 own the shared land column, vegetation exchange and phase
  semantics carried by EFOR-8;
- BIO-22 owns absolute-time conversion of ecological rates; and
- BIO-20 owns the cheap adapter/consumer regression boundary.

EFOR supplies those tasks with one chronological, conservative source instead
of letting each reconstruct weather independently.
