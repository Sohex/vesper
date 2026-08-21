# Reading a World Orogen export

The format reference is `vendor/orogen/tools/README.md`, which is authoritative
and should be read before writing anything that consumes this directory. This
file is the consuming side: the conventions, the traps, and the fields whose
names do not mean what they look like.

**Nothing here is a current value.** Counts, fractions and elevations belong to
whichever build is in hand; read them from `manifest.json`, from
`world_state.json`, or by measuring. What is written down here is which field to
trust and why.

`lib/orogen.py` wraps most of this and refuses a build it has not been checked
against. Prefer it to reading the files directly.


Five exports of the same planet per build, all from seed 16236323, and all
carrying the same `manifest.hashes.finalElevation` within a build. The REGION
COUNT is a property of the build rather than of the project and is in
`manifest.numRegions`: `precarve-craton` is 2,500,001 at a 15.19 km mean edge,
and 10,000,004 is the default for anything generated from now on, for the reason
under the recipe below. Check that hash against `lib/orogen.py` before trusting any number quoted
about the terrain: the seed and parameters alone do not identify a build, because
fixes to the generator change the terrain under a fixed seed, and one such fix
(over-erosion) moved mean land elevation by a factor of four. Anything derived
from a superseded terrain is not comparable, and the registry note on each build
says what moved.

| Directory | Grid | Notes |
| --- | --- | --- |
| `exoplasim-T21/` | 64×32 Gaussian | grid only |
| `exoplasim-T42/` | 128×64 Gaussian | the only one with `raw/` (native 2.5M-region mesh) |
| `exoplasim-T63/` | 192×96 Gaussian | grid only |
| `exoplasim-T85/` | 256×128 Gaussian | grid only |
| `grid-512x256/` | 512×256 uniform | grid only; for mapping/visualisation |

Each has `manifest.json` (the field catalogue -- path, dtype, shape, units,
description, plus semantic blocks for lithology, basins, hydrology, plates),
`README.txt`, `grid/<field>.bin`, and `planet.nc`. The manifest is ~66 MB, mostly
the basin hypsometry catalogue; index into it rather than dumping it.

`planet.nc` is the easy path (CF-1.8; `manifest.json` is the field catalogue
and carries the count). Use the `.bin` files when you
need `raw/` or `gauss_weights.bin`.

## Choosing a land mask -- use `surface_class`, and only `surface_class`

Two land definitions exist and they disagree:

| Source | Definition |
| --- | --- |
| `surface_class == 1` | **authoritative**; ocean means connected to the world ocean |
| `land_mask` | elevation-sign test, `elevation_km > 0` |

The disagreeing regions are dry closed-basin floor lying below sea level, which
`land_mask` would flood. Preserving that terrain is the whole point of the fork.
Both fields carry descriptions in the manifest saying so and pointing at each
other, and **`manifest.landSeaMask` states the disagreement in numbers** -- read
it from there, so it is checkable rather than folklore and cannot go stale here.

**Do not reconstruct it as `land_mask | is_endorheic`.** Most of the disagreeing
regions sit inside a preserved basin, but not all: the rest are smaller enclosed
depressions that `fixupTopology` kept as genuinely not-sea but that never cleared
the basin selection thresholds, so they are absent from the basin catalogue and
carry `basin_index == -1`. The naive union misses exactly those and silently
floods them. `surface_class` is the only correct source.

`surface_class == 2` (`inland_water`) is **empty**, by design rather than
oversight -- and filling it is now done downstream: `hydrography/scripts/surface_water.py`
solves lake extent and `build_surface_albedo.py --lakes` carries it into the
climate as a composite albedo. Orogen measures basin geometry but never decides water levels -- that
is a precipitation-versus-evaporation balance and belongs downstream. A large fraction of
this planet is flagged `is_endorheic`, and it falls with each carve iteration. Filling those basins is our job,
using the hypsometry curves in `manifest.basins.preserved[]`, and it feeds back
into climate through albedo and evaporation.

## Masks: take them from `planet.nc`, not from an image

`surface_class` in `planet.nc` is the mask. `build_boundary_conditions.py`
integrates it from the native mesh, which is the current path and the correct
one.

An image path exists and nothing current uses it. If you are working from the
PNGs anyway,
`orogen-surfacemask-*.png` is the right one -- white land, grey inland water,
black ocean, so thresholding at `> 0` gives land/sea, and grey appears only once
water levels are assigned downstream. `orogen-landmask-*.png` keeps its old
`elevation > 0` meaning rather than being silently redefined, so it is the
*wrong* mask and floods the dry closed basins.

## Elevation conventions below sea level

`elevToHeightKm` takes its land flag from `surface_class`, so below-sea-level
land converts at 1.0 km per unit against the ocean's 10, and dry basin floors
agree with the catalogue. Every build registered in `lib/orogen.py` carries
this convention; an export the registry does not know may predate it, and its
dry floors read ten times too deep.

**The basin catalogue renamed keys in the same pass, and it is a breaking
change.** Unsuffixed keys (`sinkElevation`, `depth`, `spillElevation`) are the
generator's model parameter; `Km` and `Km3` keys are genuinely physical,
converted through the land branch. Volumes and hypsometry are integrated in
physical height. Read the suffixed keys unless you specifically want model units,
and do not assume an unsuffixed key is kilometres.

Note also that the top-level catalogue entry describes the **natural**
pre-conditioning basin, while `finalPreserved` describes the finished terrain.
They differ by a large factor -- the natural surface floods far more area at
spill than the finished one does. `hypsometry` is on the natural terrain, which
is why `hydrography/` recomputes it on the finished one. Compare the two in the
manifest rather than trusting a figure quoted anywhere.

## Other gotchas

- **`x`, `y`, `z` are y-up, so the polar axis is `y` and NOT `z`.** The mesh
  carries unit-sphere cartesian coordinates in the generator's own three.js
  frame, and they are exact: `lat = degrees(arcsin(y))` and
  `lon = degrees(atan2(x, z))`, both reproducing the `lat`/`lon` fields to
  within float32. Assuming the usual z-up convention silently rotates the
  planet 90 degrees and puts the pole on the equator. It is safe to use `x, y,
  z` for anything rotation-invariant -- great-circle distances between regions,
  dot products, tangent-plane fits -- which is what
  `hydrography/scripts/export_carve_list.py` and `lib/orogen.py:local_slope_deg`
  both do. Mixing them with `lat`/`lon` in the same expression is where this
  bites. Found 2026-08-18 while validating a slope field against an analytic
  ramp: the ramp was built on `z`, the expectation was written in latitude, and
  the estimator was blamed for a 30% error that was entirely in the frame.
- **`elevation` is not kilometres.** It is the generator's internal shaping
  parameter (nonlinear hypsometric curve; 0.5 ≈ 1.1 km, 1.0 = 6 km; ocean linear
  at 10 km/unit). Use `elevation_km` for physical orography.
- **Weighting.** `grid_cell_area` in the gridded output is the true cell area and
  sums exactly to 4πR², so it is a correct area weight on every grid.
  `grid/gauss_weights.bin` is equivalent on the Gaussian grids and `cos(lat)` on
  the uniform one; all three agree to ~2e-4. `raw/cell_area.bin` is the mesh
  region area. (Earlier exports called the gridded field `cell_area` and it was
  a mesh diagnostic, not a cell area -- that trap is fixed, but any code written
  against an older export needs checking.)
- **Distance fields are in cell hops, not km.** Convert with
  `avgEdgeKm = π × R / √numRegions`, and `manifest.basins.resolution.avgEdgeKm`
  has it computed for the build in hand. Read it from there.
- **Never match by longitude between the export and ExoPlaSim output.**
  `source/<build>/exoplasim-*/planet.nc` labels longitudes from −178.5938;
  ExoPlaSim's own output labels them from 0. Same grid, different labels, and the
  correct mapping is **by index**: for these files that is a roll of zero, and
  the two axes are also half a cell apart, so do not "fix" the offset by rolling
  until the labels line up. Only the coordinate axes disagree; the cells are the
  same cells.

  This is *not* a claim that the export's gridded `surface_class` and the mask we
  integrate for ExoPlaSim agree cell for cell. They do not, in both directions,
  because `build_boundary_conditions.py` integrates from the native mesh while
  the export emits by the region containing the cell centre.
  That difference is the reason the mesh integration exists and is expected. Anything that keys on lon/lat across that
  boundary silently matches zero cells, which has now happened three times on
  three different scripts. Share one coordinate source -- in practice the
  climatology, since the LPJ-GUESS driver and the pedology soil map are both
  built from it -- and never reconstruct one.
- **`manifest.planet` rotation/obliquity/eccentricity are Earth defaults**
  (23.93 h, 23.44°, 0.0167), not this world's. Orogen does not consume them, so
  they were never overridden. **`config/planet.yaml` is authoritative** for
  everything except radius and gravity.
- `manifest.planetRadiusKm` at the top level says 6371; `manifest.planet.radiusKm`
  says the real 7645.2. Trust the latter.
- **`lithology.compositionLand` is measured against `land_mask`**, numerator and
  denominator both, so it omits the dry sub-sea-level basin floors entirely. The
  distortion is not uniform: evaporite gains 1.17x on area against 1.046x for
  land overall, so its share is understated. That is
  the expected direction, because playa fill accumulates in exactly the closed
  basins `land_mask` excludes. Compute composition from `substrate_class` and
  `surface_class` rather than quoting the table.

## Generating a build

**A build is disposable until a climate run has consumed it.** Nothing derived
exists until hydrography and ExoPlaSim have been run on it, so before that point
the right response to any change in the generator is to regenerate rather than
to migrate. After that point it is not, because runs and verdicts start
depending on it and superseding it costs the whole downstream chain.

So what has to survive is the recipe, not the export. This is it:

```bash
# The planet code carries the seed and every terrain slider. It does NOT carry
# radius, gravity, or lithology strength, so those three are passed explicitly
# and are the ones to check against config/planet.yaml.
cd vendor/orogen
COMMON="--code 01eshm059lt0b9mpgro2y83t \
        --radius 7645.2 --gravity 12.81 \
        --lithology-strength 0.682 --regions 10000004 \
        --netcdf --quiet"

node --max-old-space-size=32000 tools/export-planet.mjs $COMMON \
     --grid T42 --out ../../source/<build>/exoplasim-T42          # only this one keeps raw/
for G in T21 T85 T127 T170; do
  node --max-old-space-size=32000 tools/export-planet.mjs $COMMON \
       --grid $G --no-raw --out ../../source/<build>/exoplasim-$G
done
node --max-old-space-size=32000 tools/export-planet.mjs $COMMON \
     --grid 512x256 --no-raw --out ../../source/<build>/grid-512x256
```

**The ladder is T21/T42/T85/T127/T170**, and it is the same list as
`exoplasim/scripts/rebuild_binaries.py`'s `MATRIX` on purpose: a truncation with
an export and no binary is terrain nothing can run, and a binary with no export
has no geography to run on. T63 is no longer emitted; `precarve-craton` carries
one because it predates the ladder.

`T170` and `grid-512x256` are both 512 by 256 and are NOT interchangeable: the
truncation emits Gauss-Legendre latitudes, which is what a spectral model wants,
and the other is uniform in latitude. Keep both; the shape they share is a
coincidence of the truncation table.

**`--regions 10000004` is the default, and it is a measured choice rather than a
preference.** `notes/audits/orogen-resolution.md` has the argument; the three
numbers that decide it are these. Orogen's terrain noise sits at fixed physical
wavelengths and the finest it DESIGNS is about 20 km, so `elevation_pre_erosion`
gains only 8% of semivariance per lag doubling below 10 km against 26% at 50 km;
a 7.60 km mesh oversamples that floor by 2.6x and a finer one manufactures
erosion texture at the mesh scale rather than resolving terrain. The endorheic
basin catalogue CONVERGES there: preserved basins go 3,621, 6,345, 9,419, 9,649
across 2.5M, 5M, 10M and 25M regions, so the step past 10M adds 2.4% where the
one before it added 48%, because `minCells` has stopped binding and the declared
1000 km2 floor governs. And the cost is superlinear: 1,093 s and 9.0 GB at 10M
against 3,875 s and 21.1 GB at 25M.

At 10,000,004 the T42 export with `raw/` takes about eighteen minutes and 9 GB,
and the four `--no-raw` grids are quicker. The heap has to be raised from the
14336 the 2.5M build used. Then register the terrain hash in `lib/orogen.py` and
point `source_build` at it.

Add `--preserve-basins FILE` to apply a carve list; without it the build is
pre-carve, which is what a first pass on new geography wants.

The code decodes to seed 16236323, 2,500,001 regions, 100 plates, 10 continents
and the erosion and shaping sliders. An explicit flag beats the decoded value,
which is why `--regions` above overrides the count and the rest of the code still
applies. Verify against `manifest.params` after
generating: `lithologyStrength` in particular must read 0.682 and not the 1 the
code decodes to, because it is a fork addition the code predates.

## `source/maps/`

The equirectangular PNGs are gitignored payload, regenerated by the command
below rather than kept.

**A rendering is only valid for the gravity it was made at.** The heightmap
carries the `1/g` relief scaling, so a render made at one gravity reads high or
low by that ratio against a mesh at another, and nothing in the image says which
gravity produced it. The land MASK is unaffected: a mask does not scale. Check
the render against `elevation_km` from the mesh before trusting a height read
off an image.

That is why the command takes gravity and radius from `config/planet.yaml`
instead of having them typed. A planet code encodes sliders only -- never radius
or gravity -- so the code alone does not pin a rendering.

To render them:

```bash
# Gravity and radius are read from the config rather than typed, because a typed
# copy of either renders the maps for a planet this project no longer has. Run
# these two from the repo root, before cd-ing.
GRAV=$(python -c "import yaml;print(yaml.safe_load(open('config/planet.yaml'))['planet']['gravity_m_s2'])")
RAD=$(python -c "import yaml;print(yaml.safe_load(open('config/planet.yaml'))['planet']['radius_earth']*6371)")

cd vendor/orogen
node --max-old-space-size=12288 tools/export-maps.mjs \
    --code 01eshm059lt0b9mpgro2y83t \
    --radius $RAD --gravity $GRAV --width 16384
```

## Gravity is applied at export, and only at export

Orogen runs its whole pipeline in model units; the 1/g relief scaling is applied
at the model-unit-to-km conversion on export. So two builds differing only in
gravity are **bit-identical in every hash** -- `finalElevation`, `basinCatalogue`
and `params` alike -- and differ only in `manifest.planet.gravityMS2` and the
`elevation_km` it scales.

- The terrain hash is not sufficient identity for us. `lib/orogen.py` reads
  `gravityMS2` separately and `scripts/check_consistency.py` checks it against
  the config, because a build from another gravity would otherwise pass the
  allowlist while every vertical quantity was off by the ratio.
- **Basin ids and carve verdicts survive a gravity change.** The catalogue is
  bit-identical, so an existing verdict replays. Whether it *should* is a
  separate question -- the climate driving the water balance moves -- but nothing
  forces the loop to restart from zero.
- Erosion also ran in model units, so a higher-gravity planet does not get
  steeper-slope collapse. The landscape's shape is identical; only the vertical
  scale changes.
- The catalogue's own `depthKm`, `volumeKm3` and hypsometry do **not** carry the
  scaling, while `elevation_km` does. Anything mixing the two is comparing
  verticals that differ by `reliefScale`. Our hypsometry is rebuilt from
  `elevation_km` so the carve criterion is safe; the one exception is a reported
  diagnostic, which says so.
- `orog_mean/std/min/max` are declared `units: 'km'` and are neither scaled nor
  converted through the hypsometric curve -- they are raw model units. We do not
  consume them. Anything that starts to must convert them first.

