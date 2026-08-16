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
| input module | `vesperinput`, built and run end to end on real cells |
| soil texture from lithology | mapped, but it is parent material not soil; see below |
| full run | blocked on a current climatology; 3-cell shakedown passes |

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

## Soil is currently lithology, and that is a known shortcut

`build_lpj_driver.py` maps World Orogen rock classes onto LPJ soil texture codes
through a table, `SOIL_CODE_BY_ROCK`. That table is parent material, not soil.
Soil texture is what a rock weathers *into under a climate*, and the same granite
gives coarse grus in a cold arid place and deep kaolinitic clay in a wet tropical
one. The mapping ignores that, and so assigns the same texture to both.

Everything a real pedogenesis stage needs already exists in this pipeline, which
is what makes the shortcut conspicuous. Jenny's five soil-forming factors are
climate, organisms, relief, parent material and time: ExoPlaSim has the first,
LPJ-GUESS itself produces the second, Orogen has relief, erodibility and
lithology for the third and fourth, and the fifth is a modelling choice.

Three things the current table gets demonstrably wrong:

- **No climate dependence at all**, as above. This is the first-order error.
- **Salinity and sodicity are absent.** 12.65% of carved-zoned land is evaporite
  or playa, and the only handling is the blunt rule that nothing roots there.
  Salt-affected soil is a gradient, not a binary.
- **Regolith depth is not represented.** Orogen knows erosion rate and exhumation;
  a thin soil over bedrock holds far less water than a deep one, and LPJ-GUESS's
  soil codes cannot say so.

**Not built yet, deliberately.** The next run should measure how much this
matters before anything is built on top of it: run once with the lithology
mapping and once with uniform medium soil, and report the spread in NPP. If it is
small the mapping is adequate and a pedology stage is a refinement; if it is
large the stage is required and the current numbers carry its uncertainty. That
is the same bracket-before-building discipline the albedo and evaporation
questions got.

If it is required, it belongs as a sibling component reading lithology, climate,
relief and drainage, and feeding LPJ-GUESS's richer `SoilInput` interface, which
takes sand, clay, organic carbon, pH, C:N and bulk density per cell rather than
one of ten codes. Note that it closes another loop: soil organic matter is a
product of the biosphere that the biosphere then grows in.

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
