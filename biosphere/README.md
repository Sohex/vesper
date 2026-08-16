# Biosphere

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
world. Bare rock gives a land-mean albedo of 0.314 and a vegetated surface 0.223,
worth 3.7 to 7.1 K, and the two endmembers reach the 290-293 K design band at
non-overlapping stellar fluxes: 0.977-0.994 bare against 0.952-0.970 vegetated.
The project picked 0.96 and vegetated, which means the current baseline is
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

For the declared 0.91-1.01 cycle, the predicted swing is about 4.2 K
peak-to-peak after slab damping, so plus or minus 2.1 K on a monthly mean.
Against the PFT cold-survival thresholds, **13.8% of land sits within that
distance of one**, where the cycle can kill a PFT outright and a mean-climate run
would never see it. Survival limits are thresholds, so that is a bias rather than
noise, and it falls on biome boundaries rather than on global totals.

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
