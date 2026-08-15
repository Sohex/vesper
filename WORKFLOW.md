# Vesper: the worldbuilding pipeline

How geography, climate, water and life are computed for this planet, in what
order, and why the order is not a straight line.

Vesper is a super-Earth: 1.2 Earth radii, surface gravity 10.1989 m/s2, a 30-hour
day, 32 degrees of obliquity, e = 0.02, orbiting a K2.5V dwarf at 0.96 S-Earth
with a 189.6-day year. Mean surface temperature 292.9 K.

---

## 1. The components

```
config/planet.yaml   Canonical planet, star, orbit, atmosphere. Every component reads it.
source/              World Orogen exports. Canonical, read-only.
lib/                 Shared readers: orogen.py (the export), gridding.py (mesh to grid).
hydrography/         Drainage, catchments, basin capacity, lake balance, carve verdict.
exoplasim/           Boundary conditions, climate integrations, climatology.
```

Each component owns its own `scripts/`, `data/` or `analysis/`, and reads
`config/planet.yaml` and `source/`. Paths resolve from the file location, not the
working directory, so anything runs from anywhere.

## 2. The data flow

```
World Orogen (fork)
   |  seed + planet code -> terrain, lithology, closed basins, hydrology
   v
source/exoplasim-T21|T42|T85 + grid-512x256 + maps
   |  raw/ mesh lives in the T42 export and is identical for all of them
   v
lib/gridding.py            integrate mesh fields onto any model grid
   |                        |
   |                        v
   |                  exoplasim/build_boundary_conditions.py -> land mask, topography
   |                  exoplasim/build_surface_albedo.py      -> albedo, forest fraction
   |                        |
   |                        v
   |                  ExoPlaSim spin-up -> climatology
   |                        |
   v                        v
hydrography/build_hydrography.py     hydrography/carve_verdict.py
   drainage, catchments,             which basins overflow
   hypsometry, coupling matrix   <-- integrates climate over catchments
   |
   v
carve list -> back to World Orogen -> new terrain
```

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
ocean, interactive sea ice, glaciers enabled, measured K2 stellar spectrum rather
than a blackbody.

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
| biosphere | assumed, not modelled |
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

1. Send the carve list to World Orogen; regenerate the terrain.
2. Rebuild hydrography and boundary conditions on the new terrain.
3. Re-run the T42 baseline; check whether the verdict holds.
4. Run LPJ-GUESS on the climatology to replace the assumed biosphere, then feed
   the real vegetation back as albedo and forest fraction.
5. Once terrain and biosphere have settled, one T85 equilibrium for the regional
   products: biomes, Koppen, lake maps.
6. The stellar cycle last, 0.91 to 1.01 S-Earth over 8 Earth years, on the
   settled world.

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
