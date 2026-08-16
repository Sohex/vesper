# Vesper: the worldbuilding pipeline

How geography, climate, water and life are computed for this planet, in what
order, and why the order is not a straight line.

Vesper is a super-Earth: 1.2 Earth radii, surface gravity 10.1989 m/s2, a 30-hour
day, 32 degrees of obliquity, e = 0.02, orbiting a K2.5V dwarf at 0.96 S-Earth
with a 180.7-day year. Mean surface temperature 292.9 K.

---

## 1. The components

```
config/planet.yaml   Canonical planet, star, orbit, atmosphere. Every component reads it.
source/              World Orogen exports. Canonical, read-only.
lib/                 Shared readers: orogen.py (the export), gridding.py (mesh to grid).
hydrography/         Drainage, catchments, basin capacity, lake balance, carve verdict.
exoplasim/           Boundary conditions, climate integrations, climatology.
pedology/            Weathers lithology into soil texture, pH and organic content.
biosphere/           LPJ-GUESS: vegetation, leaf area, carbon, PFT composition.
```

Each component owns its own `scripts/`, `data/` or `analysis/`, and reads
`config/planet.yaml` and `source/`. Paths resolve from the file location, not the
working directory, so anything runs from anywhere.

## 2. The data flow

```
World Orogen (fork)
   |  seed + planet code -> terrain, lithology, closed basins, hydrology
   v
source/<build>/exoplasim-T21|T42|T85 + grid-512x256 + maps
   |  raw/ mesh lives in the T42 export and is identical for all of them
   v
lib/gridding.py            integrate mesh fields onto any model grid
   |                        |
   |                        v
   |                  exoplasim/build_boundary_conditions.py -> land mask, topography
   |                  exoplasim/build_surface_albedo.py      -> albedo, forest fraction
   |                  exoplasim/build_surface_soil_water.py  -> dwmax  (off by default)
   |                        |
   |                        v
   |                  ExoPlaSim spin-up
   |                        |
   |                  build_climatology.py -> averaged climatology
   |                                       -> per-orbit climatologies (--per-year)
   |                                       -> climate series, per bin per orbit
   |                        |
   |          +-------------+-------------+
   |          |                           |
   v          v                           v
hydrography/            pedology/build_soil.py        hydrography/carve_verdict.py
  drainage,               weathers lithology            which basins overflow
  catchments,             under the climate         <-- integrates climate over
  hypsometry,               |                           catchments
  coupling matrix           v
   |                  biosphere/build_lpj_driver.py -> one binary, N years
   |                        |
   |                        v
   |                  biosphere/run_lpj_guess.py  (LPJ-GUESS, MPI)
   |                        |
   |                        +--> cpool.out ---> back to build_soil.py
   |                        |                   (soil and biosphere iterate)
   |                        v
   |                  fpc.out -> build_surface_albedo.py --mode modelled
   |                        |     albedo and forest fraction from what grew
   |                        v
   |                  back to ExoPlaSim
   v
carve list -> back to World Orogen -> new terrain
```

Three loops close in that diagram, and section 4 says why each has to.

## 3. The pipeline, step by step

### 3.1 Geography

World Orogen generates the terrain from a planet code, which encodes seed and
every slider. The current build is `01eshm059lt0b9mpgro2y83t`: seed 16236323,
2,500,001 mesh regions, closed basins preserved, lithology-modulated erosion.

Verify a build by `manifest.hashes.finalElevation`. A seed alone does not
identify a planet, because fixes to the generator change the terrain under a
fixed seed. The current terrain is `821aa71b37a7...`.

The fork exports Gaussian grids directly off the mesh at T21, T42 and T85, plus a
uniform 512x256 for mapping, plus PNG maps. Only the T42 export carries `raw/`,
the native mesh, and that mesh is the same for every export.

### 3.2 Boundary conditions

Both the land mask and the topography are integrated from the native mesh, not
sampled from a gridded export. The export resamples categorical fields by the
region containing the cell centre, which at T42 discards the coastline, and
averages continuous fields over all regions in a cell, which drags a coastal
cell's albedo toward open water's 0.06.

Land comes from `surface_class`, never from `land_mask`. The two disagree by 1.9%
of the planet, all of it dry closed-basin floor below sea level, and `land_mask`
would flood it.

Six surface fields are supplied. Topography (129), land mask (172), broadband and
two-band albedo (174, 175, 176) and forest fraction (212). Everything else falls
back to a uniform namelist default, which is declared rather than accidental:
Earth's roughness and vegetation maps are tied to Earth's continents and would be
meaningless here.

### 3.3 Climate

ExoPlaSim 3.4.2, T42, 16 MPI ranks, 45-minute timestep, 10 layers, 50 m slab
ocean, interactive sea ice, glaciers enabled, and a measured stellar spectrum
rather than a blackbody.

That spectrum is `k25v`, built from BT-Settl and interpolated to 4965 K, and it
replaces a wrong one that every run of the first three eras inherited.
ExoPlaSim's `k2.dat` is the star K2-18, an M2.5V at about 3450 K, not the
spectral type K2; the package ships no K dwarf spectrum at all. Snow, ice and
glacier albedos were consequently 0.10 to 0.17 too dark everywhere, which
weakened the very feedback the glacier and stellar-cycle machinery exists to
resolve. The model confirms the fix from its own log: energy fraction below
0.75 microns is 0.38438 under `k25v` against 0.11588 under `k2`.

Everything before the re-baseline carries that bias. It is not invalidated in
kind and the direction is known, but it should be stated wherever it is quoted.
See `exoplasim/notes/stellar-spectrum-audit.md`.

A run directory names its spectrum, as it already named its geography and its
physics switches. Two runs differing only in spectrum are different climates, and
without the marker a re-baseline would have landed in the completed run's
directory. Directories predating the marker carry none and are all `k2`.

Convergence is a fixed six-part test, not a judgement: temperature drift below
0.05 K per orbit, top-of-atmosphere and surface balance trends below 0.05 W/m2
per orbit, sea-ice drift below 0.001 of planetary area per orbit, and mean
absolute imbalances below 0.5 W/m2. Runs that miss are labelled, not rounded.

A separate five-orbit window with 32 snapshots per orbit, excluded from the
pass/fail decision, forms the climatology.

### 3.4 Hydrography

The export deliberately does not route water: `drain_to` is raw steepest descent,
and on this planet 63% of the land drains into 220,649 unpreserved pits. Routing
is a hydrology decision and the exporter leaves it downstream.

`build_hydrography.py` resolves it with a priority flood over the 2.5M-region
mesh, which fills the noise pits while keeping the 3,629 preserved basins as
terminals. It then rebuilds each basin's hypsometry on the *finished* terrain,
because the catalogue's curves are measured on the pre-conditioning surface and
overstate capacity by about 1.5x. Finally it writes a sparse basin-by-grid-cell
coupling matrix, which is the interface climate is integrated over.

### 3.5 The carve verdict

A basin that overflows year on year incises its outlet, and over 1e4 to 1e6 years
that drains the lake and the depression stops existing. So a basin pinned at its
spill is a transient, not a landscape state.

The test is geometry on one side and climate on the other:

    (E - P) / runoff  <=  catchment / area_at_spill - 1

The right-hand side ships in `basins.nc` as `critical_aridity_index`. The
left-hand side is the climatology integrated over each catchment.

`E` is evaporation from open water, which the model does not have there. It is
estimated with the Penman combination equation using water's albedo and roughness,
and validated by applying the same calculation to ocean cells, which *are* open
water: 3.736 mm/day against the model's own 3.672, a ratio of 1.017.

### 3.6 Pedology

`pedology/` weathers the lithology into soil under the climate, which is the step
neither the climate model nor the vegetation model does. Parent material sets
which minerals are available, climate sets how far they have been converted,
relief and erosion set how much regolith survives, and the biosphere sets the
organic fraction.

Weathering intensity is the Walker-Hays-Kasting form, a power of runoff times an
exponential of temperature, normalised so Earth's land mean is 1. Texture comes
from mixing parent materials and converting weatherable primary minerals to clay,
with quartz tracked separately because it never becomes clay: that single fact is
why granite and basalt diverge under identical climate.

Runoff means `P - E`, not the model's `mrro`. `mrro` is river-routed net water
flux, so it is negative in places and non-zero over ocean; using it understates
land runoff by 6.6x. See `exoplasim/notes/water-and-energy-closure.md`.

Every Earth calibration lives in `pedology/config/pedogenesis.yaml` with its
source. Nothing in the scripts hardcodes a Vesper number, so porting the
component to another world is a config change.

### 3.7 The biosphere

`biosphere/` runs LPJ-GUESS 4.1.1, patched for this world's calendar and
astronomy. The patch carries no planetary numbers itself: it points the model at
a generated `vesper.h`, which `build_vesper_header.py` derives from
`config/planet.yaml`, because the year length is a function of the stellar flux.

The model steps in 24-hour days with a 181-day year. That keeps every per-day
rate constant calibrated against the absolute time it was calibrated against and
confines the error to daylength, where its sign is known. Stepping in real 30-hour
Vesper days would put 25.9% into respiration, decomposition and phenology alike.

Annual degree-day limits are rescaled by the same orbit-derived factor, because a
181-day year reaches half the annual GDD of a 365-day one and Earth thresholds
would otherwise exclude every tree for reasons unrelated to the climate.

Forcing arrives as one binary driver file carrying however many years of climate
it was built with, and the model cycles through them. One year is a fixed
climate. Several are how a variable star reaches the biosphere; a single
repeating year cannot represent one at all.

## 4. Why this is not a straight line

Three quantities each depend on the other two.

**Drainage depends on climate.** Which basins survive is a water balance.

**Climate depends on drainage.** Evaporite forms in closed basins and is 20.8% of
this planet's land at albedo 0.50, the brightest class there is. Carve the basins
and the world gets darker.

**Climate depends on the biosphere, and the biosphere on climate.** Bare rock
gives a land-mean albedo of 0.314, a vegetated surface 0.223. That difference is
worth 3.7 to 7.1 K, and the two reach the 290 to 293 K design target at
*non-overlapping* stellar fluxes: 0.977 to 0.994 bare, 0.952 to 0.970 vegetated.
No single flux is robust to the question, so the orbit and the biosphere are one
choice, not two.

**Soil depends on the biosphere, and the biosphere on soil.** Texture, pH and
regolith depth are weathering products of lithology under a climate, but the
organic fraction is what the vegetation leaves behind, and it changes the bulk
density and water-holding capacity the vegetation then grows in. `pedology/`
therefore iterates against `biosphere/` rather than running once before it.

The loop is therefore: assume, compute, feed back, repeat.

Monotone carving guarantees the loop *terminates*, since basins are only ever
removed. It does not guarantee it lands in the right place, and it has a
direction: carving removes evaporite, which is the brightest lithology, so the
land darkens, the world warms, open-water evaporation rises, and basins that were
marginal would have stayed closed. Iteration 1 therefore carves at the coolest,
brightest state available and cannot take any of it back, so the pipeline
systematically over-carves. Orogen's applied export bears this out directly:
evaporite fell from 20.8% of land to 12.65%.

On iteration 2, the already-carved set should be re-evaluated against the new
climate and the number that would no longer have carved reported. That number is
the overshoot, and it is the honest measure of how much the first pass cost.

## 5. Where the pipeline currently stands

| stage | state |
| --- | --- |
| terrain | `821aa71b`, 3,629 preserved basins, no carve list applied |
| drainage | resolved, 99,085 pits filled at a median 8.5 m |
| hypsometry | rebuilt on the finished terrain |
| coupling | T42 and T85 matrices built |
| climate | T42 at 0.96 S-Earth, vegetated, converged on all six criteria, 292.88 K |
| carve verdict | first pass complete: 1,522 of 3,629 carve, 235 partial, 1,872 preserved |
| pedology | built, closing the loop through soil carbon; weathering 0.50 on `P - E`, bracketed to 2.81 if driven by precipitation |
| biosphere | LPJ-GUESS ported, driven, parallel, harnessed and scored; awaiting a current climatology for the first full run |
| stellar cycle | deferred to last |

Selected results:

- Land is 43.17% of the surface by `surface_class`.
- Before carving, 76% of land drains to a closed basin, against roughly 13% on
  Earth. After the first carve verdict that falls to 50.3% under the Penman
  estimate, or 33.6% under the land-evaporation sensitivity. Orogen's applied
  export measures 55.15%.
- 1,248 basins have zero catchment runoff and survive as dry salt pans rather
  than lakes. That agrees independently with evaporite being 20.8% of the land.
- Surviving lakes total about 1% of the planet's surface.
- Precipitation 2.95 mm/day, sea ice 0.12%, planetary albedo 0.152.

## 6. What happens next

Three nested loops and then a resolution change. The order matters in places
where it is not obvious, so those places are called out rather than left to be
rediscovered.

**A. The terrain loop.** Carve list to World Orogen, regenerate the terrain,
rebuild hydrography and boundary conditions on it, re-run the T42 baseline, take
the verdict again. Hydrography has to be rebuilt *after* the carve and before the
climate, because the carve changes the drainage the climate is integrated over.

**B. The soil and biosphere loop, at T42.** For a given climate:

1. `pedology/scripts/build_soil.py`, with no biosphere on the first pass.
2. `biosphere/scripts/build_lpj_driver.py`, then `run_lpj_guess.py`.
3. `build_soil.py --soil-carbon <run>/cpool.out`, then LPJ-GUESS again.
4. Repeat 3 until the criteria in `pedogenesis.yaml` are met.

Check whether step 3 moves anything before assuming it needs iterating:
LPJ-GUESS computes its own soil carbon internally, so the pedology organic
feedback may be second-order.

**C. The vegetation-climate loop.** `build_surface_albedo.py --mode modelled`
turns the run's foliar cover into surface albedo and forest fraction, then the
climate runs again on it. Two checks belong here and neither is optional.

*The flux.* If the modelled land albedo differs much from the assumed 0.223 the
world may leave the 290-293 K design band. **Move the flux by changing the star's
luminosity, not the orbit.** The year length depends on the semimajor axis, which
is derived as sqrt(L/F), so changing L and F together leaves the orbit and the
calendar untouched: a 4% flux change needs 2% of luminosity, which is 0.5% of
effective temperature, inside the spectral-type uncertainty already declared.
Moving the orbit instead changes the year length, which is compiled into
LPJ-GUESS and would force a rebuild and a driver regeneration.

*The carve verdict.* It was taken on assumed vegetation. Real vegetation changes
the climate, which changes evaporation over catchments, which can change the
verdict. Re-run it; if basins flip, back to loop A.

**D. The resolution change.** Only after A, B and C have settled.

1. Climate at T85, with the T42 vegetation regridded as its boundary condition.
2. Soil and biosphere at T85 on that climatology, loop B again.
3. Climate at T85 again with T85 vegetation, unless step 2's vegetation turns out
   close to the regridded T42 field, which is a cheap comparison worth making
   first.

**The biosphere never decides the resolution.** LPJ-GUESS gridcells are
independent columns, so its cost is linear in cell count and trivial either way:
4,106 land cells at T42 and 16,489 at T85, 23 against 94 minutes on 16 ranks.
Running it at T85 on T42 forcing would resolve detail that is not in its input.
Expect T85 to lower total NPP, and treat that as a resolution bias rather than a
result: productivity saturates with water, so averaging the forcing before the
model sees it inflates the answer.

**E. The stellar cycle, last, on the settled world.** 0.91 to 1.01 S-Earth over 8
Earth years. Build the climatology with `--per-year`, pass the sequence to
`build_lpj_driver.py --climatology y0.nc y1.nc ...`, and the biosphere sees the
cycle rather than its average. This matters beyond totals: about 14% of land sits
within one cycle's swing of a PFT cold-survival threshold, and thresholds do not
average. Then the regional products, biomes, Koppen and lake maps.

## 7. Conventions this project holds to

Provenance travels with every artifact: config, input hashes, software versions
and the terrain hash go into each run manifest and each analysis report.

Thresholds are fixed before results are seen. The albedo bracket's convergence
and marginality criteria were written into the script's docstring before any run
finished, and were then applied against a result that cleared one of them by
0.05 K.

Estimates that cannot be verified are bracketed rather than guessed, and the
bracket is reported. This is how albedo, evaporation and the carve verdict are
all handled.

Claims are checked against the artifact rather than the documentation. Several
findings in this project's history came from comparing two products that were
supposed to agree and finding they did not.
