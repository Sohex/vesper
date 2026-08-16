# Biosphere

Replaces the assumed vegetation with a modelled one. LPJ-GUESS takes the
ExoPlaSim climatology and returns leaf area, carbon and plant functional type
composition per gridcell, which becomes the surface albedo and forest fraction
that the next climate run is forced with.

Nothing in here runs yet. This directory holds the port audit and the built
model; see `notes/lpj-guess-porting-audit.md` for what was found and what order
the remaining work goes in, and `notes/productivity-prediction.md` for what the
answer is expected to be, registered before the model can contradict it.

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
| calendar patch | not written |
| input module | not written |
| soil texture from lithology | not designed |
| first run | blocked on a `carved-zoned` climatology |

The model is not in this repository. It lives at
`/home/cfutro/git/lpj-guess/guess_4.1` beside `ExoPlaSim` and
`planet_heightmap_generation`, for the same reason those do.

```bash
cd /home/cfutro/git/lpj-guess/build && cmake ../guess_4.1 -DCMAKE_BUILD_TYPE=Release && make -j16
```

## The three things that decide whether this is credible

**The calendar.** Vesper's year is 180.655 Earth days and its day is 30 hours,
so nothing about LPJ-GUESS's 365 x 24 h grid survives contact. The audit
recommends stepping in 24-hour days with a 181-day year, which keeps every
per-day rate constant calibrated against the absolute time it was calibrated
against and confines the error to daylength alone.

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
a prediction. Their degree-day thresholds do need rescaling, because a 181-day
year accumulates about half the annual degree-days a 365-day year does for the
same temperatures.

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

## Cost

4106 land gridcells at T42, about 24 s each for 550 years, so roughly 1.7 hours
across 16 cores. The vegetation loop is not the expensive half of this
iteration; ExoPlaSim is.
