# Biosphere

> Figures in this document are illustrative of method, and were measured on
> builds and climates that have since moved. Current values live in
> `world_state.json`, generated from the artifacts. See the convention in
> `CLAUDE.md`.



Replaces the assumed vegetation with a modelled one. LPJ-GUESS takes the
ExoPlaSim climatology and returns leaf area, carbon and plant functional type
composition per gridcell, which becomes the surface albedo and forest fraction
that the next climate run is forced with.

The port is complete and runs end to end on real Vesper cells at smoke scale.
What remains before a full run is a current climatology. The audit behind the
porting choices is in `notes/lpj-guess-porting-audit.md`, and
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
| model obtained | LPJ-GUESS 4.1.1, MPL-2.0, Zenodo 8065737, md5 verified |
| build | clean on GCC 16.1.1, MPI on, netCDF off (not installed, not needed) |
| smoke test | bundled 3-cell demo, 550 years, 73 s, expected PFTs |
| Earth-assumption audit | complete, see the note |
| productivity prediction | registered, unscored |
| calendar and astronomy patch | written, applied, verified |
| PFT degree-day rescale | generated from the orbit, 500 -> 247 gdd5min_est |
| input module | `vesperinput`, runs end to end, splits across MPI ranks |
| soil | from `pedology/`, loop closing at smoke scale |
| run harness | not written; runs so far are hand-assembled in scratch |
| albedo and forest feedback | not written, and it is the component's purpose |
| full run | blocked on a current climatology; ~25 min on 16 ranks |

The model is not in this repository. It lives at
`/home/cfutro/git/lpj-guess/guess_4.1` beside `ExoPlaSim` and
`planet_heightmap_generation`, for the same reason those do.

```bash
cd /home/cfutro/git/lpj-guess/guess_4.1
patch --forward --strip=1 --directory=. < <world>/biosphere/patches/lpj-guess-4.1.1-vesper.patch
cd ../build && cmake ../guess_4.1 -DCMAKE_BUILD_TYPE=Release && make -j16
```

`--reverse` restores the pristine 4.1.1 tree. The patch itself contains no
planetary numbers: it points `framework/guessmath.h` at a generated `vesper.h`
and wires `vesperinput` into the build. Copy `src/vesperinput.*` into `modules/`
alongside it. See "Running it" below for the generators.

## The three things that decide whether this is credible

**The calendar.** Vesper's year is 180.655 Earth days and its day is 30 hours,
so nothing about LPJ-GUESS's 365 x 24 h grid survives contact. Settled and
patched: 24-hour steps, 181-day year, which keeps every per-day rate constant
calibrated against the absolute time it was calibrated against and confines the
error to daylength alone. Verified against the unpatched model on identical
forcing, where annual evapotranspiration falls to 0.492 of its former value
against an expected 0.496.

**The PAR fraction.** LPJ-GUESS assumes half of shortwave is photosynthetically
active, which multiplies straight into productivity. Deriving the right value for
this star turned up a bug in the climate runs: ExoPlaSim's `k2.dat` is the star
K2-18, an M2.5V at about 3450 K, not the spectral type K2. Using it would have
put the PAR fraction at 0.078 instead of 0.309, a factor of four on
productivity. Resolved: a correct BT-Settl K2.5V spectrum now exists at
`exoplasim/inputs/stellarspectra/k25v`, and `FRADPAR` is 0.40. See
`exoplasim/notes/stellar-spectrum-audit.md` for what it means upstream.

**The PFTs.** The shipped plant functional types are Earth's, and their
bioclimatic limits are Earth calibrations. Keeping them is defensible as an
Earth-analogue biosphere and should be declared that way rather than presented as
a prediction. Their degree-day thresholds must be rescaled by 0.4946, and this is
not optional: running the patched model on Earth's own demo data collapses boreal
needleleaf and temperate broadleaf to grass, because Earth `gdd5min` cannot be
met in 181 days. `build_vesper_pfts.py` does it, deriving the factor from the
configured orbit; `gdd5min_est` 500 becomes 247. It deliberately leaves
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
weathering intensity is currently bracketed by a factor of 14.6 on the climate
model's runoff.

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
LPJ-GUESS: one hard winter does not extirpate a species, so the model smooths
before it kills. Twenty simulation years is 9.9 Earth years, which spans about
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
is 0.4946 Earth years, so the shipped 500 gives this world **247 Earth years** of
vegetation and soil development where Earth practice assumes 500.

`build_vesper_pfts.py` therefore scales year *counts* up by 2.022, the reciprocal
of the factor it scales annual *sums* down by. Confusing those two directions
would be worse than doing neither, so both lists are named in that script:

| parameter | shipped | rescaled | why |
| --- | --- | --- | --- |
| `nyear_spinup` | 500 | 1011 | time to reach steady state |
| `distinterval` | 100 | 202 | disturbance return time |
| `freenyears` | 100 | 202 | time to build an N pool before N limits |
| `estinterval` | 5 | 5 | counted in growing seasons, not absolute time |

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
There is no daily-minimum mortality anywhere in the model; every cold limit runs
through `mtemp_min20`, a twenty-year mean of monthly means. So a wider diurnal
range is physically real here and radiatively present in the climate, and the
vegetation model is structurally blind to it.

Recorded rather than worked around. Representing it would mean adding a frost
mortality mechanism that LPJ-GUESS does not have, which is a much larger change
than this project needs, and the data is already in the driver if it is ever
wanted.

## Cost

4106 land gridcells at T42, about 24 s each for 550 years, so roughly 1.7 hours
across 16 cores. The vegetation loop is not the expensive half of this
iteration; ExoPlaSim is.

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
cd /home/cfutro/git/lpj-guess/build && make -j16
```

### Fire is GLOBFIRM, and `cflux.out` is the only place it shows

`run_lpj_guess.py` writes `firemodel "GLOBFIRM"` after its `import` of
`vesper_pfts.ins`, which still carries the shipped `firemodel "BLAZE"`. The later
declaration wins, so the generated instruction file is what decides this and the
PFT file is not worth editing. BLAZE wants a SimFIRE input built from Earth
observations, and it forces `weathergenerator "GWGEN"`, which wants sub-daily
statistics this world does not have; the shipped demo makes the same two
substitutions for the same reason. The mismatch is not silent: LPJ-GUESS aborts
on BLAZE with a non-GWGEN generator, so a run that starts is a run on GLOBFIRM.

GLOBFIRM writes no `firert.out` and no burned area. Those are BLAZE-only, so the
`Fire` column of `cflux.out` is the whole diagnostic, and without it in the
harness's output list the burning is visible only as a mortality in `cmass` and
`dens` that no output explains.

One gate worth knowing when a short run reports no fire at all. With
`iftwolayersoil 0`, which `iforganicsoilproperties` requires,
`vegdynam.cpp:1386` accumulates fire-season days only after model year 100, from
a hard-coded constant that is not rescaled alongside `nyear_spinup` and
`distinterval`. The spin-up is far longer than that, so it never reaches a
reported year, but it does mean fire is off for the first tenth of it.
