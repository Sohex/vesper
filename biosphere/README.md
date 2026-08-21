# Biosphere

> Figures in this document are illustrative of method, and were measured on
> builds and climates that have since moved. Current values live in
> `world_state.json`, generated from the artifacts. See the convention in
> `CLAUDE.md`.



Replaces the assumed vegetation with a modelled one. LPJ-GUESS takes the
ExoPlaSim climatology and returns leaf area, carbon and plant functional type
composition per gridcell, which becomes the surface albedo and forest fraction
that the next climate run is forced with.

The model is mechanically integrated and has run end to end on real Vesper
cells at smoke scale. A follow-up source audit exposed ecological-calendar and
time-base correctness work below that interface; the smoke result is therefore
an integration result, not yet a meaningful biosphere result.
An interpretable full run still needs the surface-area, forcing, aggregation and
acceptance contracts tracked as BIO-11 through BIO-15 and the newly exposed
calendar, time-base and forcing corrections in BIO-21 through BIO-25.  The audit
behind the original porting choices is in `notes/lpj-guess-porting-audit.md`, the
remaining modelling gaps are evidenced in `notes/modelling-gap-audit.md`, the
implicit Earth assumptions below the port are in
`notes/implicit-earth-assumptions.md`, the soil decomposition, C-N-P and
pedology/groundwater seams are audited in
`notes/soil-decomposition-biogeochemistry-audit.md`, plant physiology and carbon
allocation are audited in `notes/plant-physiology-carbon-allocation-audit.md`,
BVOC emissions, secondary organic aerosol and atmospheric coupling are audited
in `notes/bvoc-soa-atmospheric-coupling-audit.md`, wetlands, peat and methane
are audited in `notes/wetlands-peat-methane-audit.md`, and
`notes/productivity-prediction.md` registers what the answer should be before the
model can contradict it.

## Why this component exists

The albedo bracket left the biosphere as the largest unresolved lever on this
world. Bare rock gives a land-mean albedo of 0.276 and a vegetated surface 0.179,
worth several kelvin, and the two endmembers reach the design band then in force
at non-overlapping stellar fluxes. The bracket's windows (0.977-0.994 bare against
0.952-0.970 vegetated) were measured on pre-carve terrain at T21; the flux has
since been measured directly on the current terrain and the vegetated baseline
is 0.945. The project picked vegetated, which means the current baseline is
habitable *because* the world is assumed to be vegetated. That assumption is
currently a constant in `config/planet.yaml` (`land_albedo_source: vegetated`,
0.15 on anything that can carry a canopy). This component is what turns it into
a result.

## State

| item | state |
| --- | --- |
| model obtained | `mateusdp/LPJ-GUESS-NTD`, tag `LPJ-GUESS-CNP_v1.0`, commit `b368b893`; MPL-2.0 notices restored from verified 4.1.1 |
| build | vendored CNP tree plus Vesper port; build verification awaits an available compute window |
| smoke test | stock 4.1.1 port: bundled 3-cell demo, 550 years, 73 s, expected PFTs |
| Earth-assumption audit | follow-up complete; BIO-21 through BIO-28 track cross-cutting non-fire assumptions, PCAR-1 through PCAR-10 track plant physiology/allocation, SDEC-1 through SDEC-9 track soil decomposition/biogeochemistry, BVOC-1 through BVOC-10 track volatile carbon through chemistry/aerosol/climate, and WET-1 through WET-11 track wetlands, peat and the methane source-to-atmosphere loop |
| productivity prediction | registered, unscored |
| calendar and astronomy port | mechanical calendar and orbital geometry applied; natural phenology still has unreachable Earth dates under BIO-21 |
| PFT degree-day rescale | generated from the orbit, 500 -> 247 gdd5min_est; the remaining annual-rate semantics are BIO-22 |
| input module | `vesperinput`, runs end to end, splits across MPI ranks |
| soil | from `pedology/`, loop closing at smoke scale |
| run harness | written; records inputs, binary and model identity in its manifest |
| albedo and forest feedback | modelled mode exists; rootable/lake and spectral corrections are BIO-17 and BIO-18 |
| aerodynamic feedback | modelled roughness is open as BIO-16 |
| BVOC/SOA feedback | LPJ source is present but off; carbon closure, PFT traits, reduced atmospheric chemistry/transport, direct optics and cloud effects are separated under BVOC-1 through BVOC-10 |
| wetlands/peat/methane | LPJ's northern-Earth peat/CH4 source is present but off; current low-latitude saturation creates water and the full extent, groundwater, peat-stock, source/sink and atmospheric closure is separated under WET-1 through WET-11 |
| full run | harness exists, but BIO-11 through BIO-15 and BIO-21 through BIO-25 must close before its output is interpreted; prior estimate ~25 min on 16 ranks |

The model is vendored at `vendor/lpj-guess/` as a git subtree from the CNP fork.
The Vesper calendar, astronomy and input-module changes live directly in that
tree, so the source being reviewed is the source being compiled and run.
Phosphorus limitation remains off until Vesper's soil-P fields and a replacement
productivity prediction land together; this is a source normalization, not a
silent change from the established C-N experiment.

```bash
cmake -S vendor/lpj-guess -B vendor/lpj-guess/build \
  -DCMAKE_BUILD_TYPE=Release -DUNIT_TESTS=OFF
cmake --build vendor/lpj-guess/build --parallel 16
```

`framework/vesper.h` is still generated and ignored; it contains the planetary
numbers. The committed port only includes it and wires the in-tree
`modules/vesperinput.*` into the build. See "Running it" below for the
generators.

## The three things that decide whether this is credible

**The calendar.** Vesper's year is a few hundred Earth days and its day is 30
hours, so nothing about LPJ-GUESS's 365 x 24 h grid survives contact. Settled and
patched: 24-hour steps, a year length derived from `lib/orbit.py` and rounded to
whole days, which keeps every per-day rate constant
calibrated against the absolute time it was calibrated against and confines the
error to daylength alone. Verified against the unpatched model on identical
forcing, where annual evapotranspiration falls to 0.492 of its former value
against an expected 0.496.

**The PAR fraction.** LPJ-GUESS assumes half of shortwave is photosynthetically
active, which multiplies straight into productivity; the wrong spectrum is a
factor of four on it (0.078 from the mislabelled `k2.dat`, the M2.5V star
K2-18, against 0.309). `FRADPAR` is derived from the BT-Settl K2.5V spectrum
at `exoplasim/inputs/stellarspectra/k25v` by `build_vesper_header.py`, which
writes the value and its derivation into `generated/vesper.h` and
`generated/vesper_provenance.json`. See
`exoplasim/notes/stellar-spectrum-audit.md`.

**The PFTs.** The shipped plant functional types are Earth's, and their
bioclimatic limits are Earth calibrations. Keeping them is defensible as an
Earth-analogue biosphere and should be declared that way rather than presented as
a prediction. Their degree-day thresholds must be rescaled to the length of a
Vesper year, and this is not optional: running the patched model on Earth's own
demo data collapses boreal needleleaf and temperate broadleaf to grass, because
Earth `gdd5min` cannot be met in a Vesper year. `build_vesper_pfts.py` does it,
deriving the factor from the configured orbit. It deliberately leaves
`phengdd5ramp` alone, which is a within-season accumulation already in absolute
time, and scaling that would be a real error.

There is a second PFT question the literature answers more sharply than expected.
Earth's 400-700 nm photosynthetic window is an accident of our star, and Lehmer
et al. 2021 predict peak pigment absorbance around a K2V at 675, 711 and 746 nm.
A 400-750 nm window on this spectrum gives 0.99 of Earth's photon flux, so a
50 nm redward shift would almost exactly cancel the dimmer, redder star. The
21% oxygen atmosphere caps how far that can go, since anything past about 800 nm
needs anoxygenic photosynthesis. See the prediction note.

## What the answer is expected to be

Registered in full, with tolerances and known biases, in
`notes/productivity-prediction.md`. The short version: a dimmer, redder star
under a clearer sky delivers 0.84 to 0.89 of Earth's PAR photons per square
metre of land, and a warmer, wetter climate more than repays it, so **per unit
area this world should be productive within about 30% of Earth**. It has 2.11
times Earth's land area, so **total NPP should land near 2.1 times Earth's**.

The line that actually decides the pipeline is not either of those. It is tree
cover, predicted at 55% of land: the climatology that produced these numbers
assumed a vegetated surface, and if LPJ-GUESS returns substantially less canopy
than that, the climate is not one that biosphere would sustain and the loop turns
again.

## Soil comes from `pedology/`

`vesperinput` takes `file_soilmap`, the map `pedology/scripts/build_soil.py`
emits, and uses LPJ-GUESS's own richer `SoilInput` path: sand, clay, silt,
organic carbon, pH, bulk density and C:N per cell rather than one of ten texture
codes. The driver file still carries a soil code per cell, derived from
lithology alone through `SOIL_CODE_BY_ROCK`, and that is the fallback when no
soil map is given. Which one was used is printed at the top of the run, because
the difference is invisible in the output otherwise.

The fallback is parent material, not soil, and the distinction matters: the same
granite gives coarse grus in a cold arid place and deep kaolinitic clay in a wet
tropical one. Prefer the soil map.

That component closes a loop with this one: soil organic matter is a product of
the biosphere, and it changes the bulk density and water-holding capacity the
biosphere then grows in. Iteration 0 runs on mineral soil; every iteration after
feeds `cpool.out` back into `build_soil.py`. See `pedology/README.md` for the
convergence criteria, which are fixed in advance, and for the finding that
the runoff-versus-precipitation choice is a 3.41x sensitivity on weathering
intensity rather than a bracket.

## Interannual forcing, and the stellar cycle

The driver file carries however many years of climate it was built with and
`vesperinput` cycles through them. One year is a fixed climate. Several are how a
variable star reaches the biosphere:

```bash
python biosphere/scripts/build_lpj_driver.py \
    --climatology year0.nc year1.nc ... yearN.nc
```

Each file contributes one year, in the order given, and the wrap means spin-up
sees the whole sequence rather than one arbitrary phase of it.

This began as a single repeating year, which was not a missing feature but a
wrong assumption: a fixed climatology cannot represent a variable star at all,
and the model would have shown a flat line no matter how the star behaved.

Verified with three years at -3, 0 and +3 K on one cell. Annual NPP locks to the
forcing period exactly and responds strongly:

| phase | offset | NPP kgC/m2 |
| --- | --- | --- |
| 0 | -3 K | 0.508 |
| 1 | 0 K | 0.377 |
| 2 | +3 K | 0.191 |

Two things worth reading off that. The response is steep, nearly halving per 3 K,
because this cell is water-limited and warming raises evaporative demand. And the
mean over the cycle, 0.359, is **4.8% below** the value at the mean climate,
0.377. That is the concavity argument made concrete: running a variable star on
its average climate over-predicts productivity, and the error is one-sided.

### Correction: the cycle does not reach survival thresholds

An earlier version of this file argued that the cycle would cross PFT
cold-survival limits on the 13.8% of land lying within 2.1 K of one, and that
thresholds do not average. **That is wrong for this model**, and the reason is
worth knowing.

`tcmin_surv` is not compared against a cold year. `vegdynam.cpp:153` tests it
against `climate.mtemp_min20`, which `driver.cpp:617-627` builds as the **mean of
the last twenty years' coldest monthly means**. That is a deliberate choice in
LPJ-GUESS: one hard winter does not extirpate a simulated species, so the model
smooths before it removes a cohort. Twenty simulation years is 9.9 Earth years, which spans about
1.2 periods of an 8-Earth-year stellar cycle, so the window averages the cycle
almost exactly.

So the cycle reaches the vegetation through **productivity**, which responds
annually and is where the 4.8% concavity effect measured above lives, and not
through survival. Biome boundaries move because growth changes, not because
plants freeze.

The 13.8% figure is still the right measure of how much land is climatically
marginal. It is not a measure of what this model will do with it.

## Resolution: T42 first, T85 only behind the climate

LPJ-GUESS gridcells are independent columns. There is no lateral flow between
them and no communication after the initial split, so cost is exactly linear in
cell count and resolution buys no dynamics, only a finer view of the forcing.

| | land cells | CPU | wall on 16 ranks |
| --- | --- | --- | --- |
| T42 | 4,106 | 6.2 h | 23 min |
| T85 | 16,489 | 25.0 h | 94 min |

**The biosphere is never the reason to choose a resolution.** Ninety-four
minutes is nothing beside a T85 ExoPlaSim equilibrium, and the pipeline already
follows `config.model.resolution`, with T85 exports present for every build. So
the biosphere should simply match whatever the climate ran at.

Running LPJ-GUESS at T85 on T42 forcing would be worse than pointless: the extra
cells would carry interpolated climate, so the model would resolve a detail that
is not in its input, and the output would look sharper than the information
behind it.

**Expect T85 to lower total NPP, and treat that as a resolution bias rather than
a result.** Productivity saturates with water, so it is concave, and averaging
the forcing before the model sees it inflates the answer. The same effect is
already measured on the Miami side of `productivity-prediction.md`: applied at
land means it gives 603 gC/m2/yr and applied per gridcell 443, a factor of 0.73
purely from resolving heterogeneity. T42 to T85 is a much smaller step than that,
but it points the same way, and it means a T42 and a T85 answer are not directly
comparable without saying so.

**The concavity is not the only thing that moves, and the other one is not a
bias.** T85 resolves higher, steeper relief: `CLAUDE.md` records it recovering
the full 5,769 m of mesh relief against T42's 5,101. Sharper orography means
stronger forced ascent on windward slopes and deeper rain shadows behind them, so
precipitation redistributes rather than merely smoothing differently. Expect
windward coasts wetter and lee basins drier, and expect that to run past the
biosphere into the carve verdict, since a drier rain-shadow basin is less likely
to overflow and more likely to survive as a basin. A T85 pass is therefore not
only a finer picture of the same world; parts of it are a different water
balance, and the carve verdict should be re-taken rather than assumed to carry
over.

## Spin-up is in simulation years, and that halves it

`nyear_spinup 500` reads like an absolute statement and is not. A simulation year
is about half an Earth year, so the shipped 500 gives this world roughly half the
absolute vegetation and soil development that Earth practice assumes.

`build_vesper_pfts.py` therefore scales year *counts* UP by the reciprocal of the
factor it scales annual *sums* DOWN by, deriving both from the configured orbit.
Confusing those two directions would be worse than doing neither, so both lists
are named in that script and the rescaled values are written into
`generated/vesper_pfts.ins` and its provenance:

| parameter | shipped | scaled? | why |
| --- | --- | --- | --- |
| `nyear_spinup` | 500 | up | time to reach steady state |
| `distinterval` | 100 | up | disturbance return time |
| `freenyears` | 100 | up | time to build an N pool before N limits |
| `estinterval` | 5 | no | counted in growing seasons, not absolute time |

Slow soil carbon does not need integrating for all of that: `ifcentury 1` solves
the equilibrium pool sizes analytically, accumulating running means between 70%
and 80% of the spin-up and solving at the end of that window
(`guess.h:3030-3034`). What the longer spin-up buys is enough absolute time for
the *vegetation* to reach steady state before that solve happens, which slow
forest succession needs.

The cost is a doubled run, 23 minutes to about 46 on 16 ranks. Not a
consideration.

## The 30-hour day widens the diurnal range, and this model cannot see it

A 30-hour rotation gives longer daytime heating and longer nighttime cooling than
Earth's, so the real diurnal temperature range is wider, and on marginal ground a
night could dip below a freezing threshold that the daily mean never approaches.

The forcing carries this correctly: `dtr` in the driver comes from ExoPlaSim's
own `maxt` and `mint`, which are timestep extrema and so include the full 30-hour
trough.

**LPJ-GUESS has nowhere to put it.** `climate.dtr` is read in exactly one place,
`bvoc.cpp:276`, for leaf temperature in the biogenic VOC scheme, which is off.
There is no daily-minimum plant mortality anywhere in the model; every cold limit runs
through `mtemp_min20`, a twenty-year mean of monthly means. So a wider diurnal
range is physically real here and radiatively present in the climate, and the
vegetation model is structurally blind to it.

Recorded rather than worked around. Representing it would mean adding a frost
plant-mortality mechanism that LPJ-GUESS does not have, which is a much larger change
than this project needs, and the data is already in the driver if it is ever
wanted.

## Cost

5.46 s of CPU per gridcell for 530 years at `npatch 5`, measured
(`notes/lpj-guess-porting-audit.md`), so 4,106 land cells at T42 are about
6.2 CPU-hours and 23 minutes of wall on 16 ranks -- the resolution table
above. The shipped demo configuration runs nearer 24 s per cell; the
difference is the patch count. The vegetation loop is not the expensive half
of this iteration; ExoPlaSim is.

## Running it

Everything flux-dependent is generated, never written down, because the year
length is a function of `orbit.baseline_flux_earth` and moves whenever the flux
does. Regenerate all three after any orbit change, and rebuild: the year length
sizes arrays at compile time, so a stale binary is silently wrong rather than
failing. The driver file records the year length it was built for and
`vesperinput` refuses a mismatch, which is the backstop for exactly that.

```bash
python biosphere/scripts/build_vesper_header.py   # vesper.h, installed into the tree
python biosphere/scripts/build_vesper_pfts.py     # degree-day limits rescaled
python biosphere/scripts/build_lpj_driver.py      # climate + soil codes + gridlist
cmake --build vendor/lpj-guess/build --parallel 16
```

### Fire is GLOBFIRM, with flux and occurrence diagnostics

`run_lpj_guess.py` writes `firemodel "GLOBFIRM"` after its `import` of
`vesper_pfts.ins`, which still carries the shipped `firemodel "BLAZE"`. The later
declaration wins, so the generated instruction file is what decides this and the
PFT file is not worth editing. BLAZE wants a SimFIRE input built from Earth
observations, and it forces `weathergenerator "GWGEN"`, which wants sub-daily
statistics this world does not have; the shipped demo makes the same two
substitutions for the same reason. The mismatch is not silent: LPJ-GUESS aborts
on BLAZE with a non-GWGEN generator, so a run that starts is a run on GLOBFIRM.

GLOBFIRM writes its carbon loss in the `Fire` column of `cflux.out` and its
inferred return time plus burned fraction in `firert.out`. The harness retains
both. BLAZE has the richer daily/monthly burned-area and SIMFIRE analysis
outputs, but those are separate from GLOBFIRM's annual occurrence diagnostic.

One gate worth knowing when a short run reports no fire at all. With
`iftwolayersoil 0`, which `iforganicsoilproperties` requires,
`vegdynam.cpp:1386` accumulates fire-season days only after model year 100, from
a hard-coded constant that is not rescaled alongside `nyear_spinup` and
`distinterval`. The spin-up is far longer than that, so it never reaches a
reported year, but it does mean fire is off for the first tenth of it.

## One-off tools

- `scripts/score_prediction.py` -- one-off: scores a productivity prediction against an LPJ-GUESS run, the machinery behind BIO-2's nitrogen bracket. Registered under `one_offs` in `config/pipeline.yaml`; it generates nothing the pipeline reads.
