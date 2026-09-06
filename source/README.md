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


Five Gaussian exports of the same planet per build, all from seed 16236323, and
all carrying the same `manifest.hashes.finalElevation` within a build. The
REGION COUNT is a property of the build rather than of the project and is in
`manifest.numRegions`: `precarve-craton` is the older 2.5M-region build at a
15.19 km mean edge, while `precarve-craton-10m` is the current fine-support
reference at about 10M regions and 7.60 km. The latter fully samples Orogen's
measured ~20 km terrain-information floor and is the support to use for new
aggregation/resolution work. Check the terrain hash against `lib/orogen.py`
before trusting any number quoted
about the terrain: the seed and parameters alone do not identify a build, because
fixes to the generator change the terrain under a fixed seed, and one such fix
(over-erosion) moved mean land elevation by a factor of four. Anything derived
from a superseded terrain is not comparable, and the registry note on each build
says what moved.

| Directory | Grid | Notes |
| --- | --- | --- |
| `exoplasim-T21/` | 64×32 Gaussian | grid only |
| `exoplasim-T42/` | 128×64 Gaussian | the only one with `raw/` (native mesh storage carrier) |
| `exoplasim-T85/` | 256×128 Gaussian | grid only |
| `exoplasim-T127/` | 384×192 Gaussian | grid only |
| `exoplasim-T170/` | 512×256 Gaussian | grid only |
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

**TWO thresholds decide this and they interact.** `selectBasins` keeps a
depression only if it clears BOTH `minAreaKm2`, a declared landform floor, and
`minCells`, a resolution floor of 12 mesh cells; the manifest publishes
`effectiveMinAreaKm2 = max(minAreaKm2, minCells * cell)` and `bindingFloor`
naming which one bound. So a statement about which floor governs is a statement
about a BUILD, and the two builds in `source/` sit on opposite sides of the
crossover: 293.80 km2 a cell at 2.5M regions puts the cell floor at 3,526 km2,
and 73.45 km2 at 10M puts it at 881.4.

**`minAreaKm2` is 850 km2** (WORLD-BVL9, 2026-08-26), which is inert at 10M
regions and binds at 10,369,311 and above -- 3.7 per cent finer. The catalogue
is therefore resolution-set at today's build and physics-set at any refinement
of it. Read `bindingFloor` from the manifest rather than assuming either.

**Builds generated before that change carry the previous 1000 km2 floor**, which
bound by 13 per cent at 10M, and a selection criterion reaches a catalogue only
through a generation. What the 850 buys at the next one, measured beforehand in
`notes/audits/basin-catalogue-floor.md`: about 1,105 more basins and 9.9e5 km2,
86 per cent of the 1.157e6 km2 that sat outside every preserved basin under the
old floor -- ordinary dry basin floor in small pieces, its deepest ground at
-0.513 km against -0.591 km inside a preserved basin.

**Lowering it further does nothing on its own.** Below the cell floor `minCells`
takes over, and `minCells` is what gives a depression enough cells for
`basinHypsometry` to build the capacity curve the carve criterion compares fill
against, so moving it is a trade rather than a tightening.

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
- **`elevation` is not kilometres, and it is not bounded by 1.** It is the
  generator's internal shaping parameter (nonlinear hypsometric curve; 0.5 ≈ 1.1
  km, 1.0 = 6 km; ocean linear at 10 km/unit). The land curve is the Hermite
  shape `t^4(5-4t)` on [0,1] and is LINEAR at 6 km per unit above it, because
  the polynomial turns over past its domain and would fall below sea level by
  1.25 -- WORLD-34GU replaced the clamp that used to guard against that with a
  branch, so heights above 1 are now distinct rather than published at one
  value. The parameter reaches 1.369 and 1.428 on the two registered builds.
  `notes/audits/relief-curve-domain.md` has the derivation. Use `elevation_km`
  for physical orography.
- **Weighting: three area quantities, three questions, and only one of them is
  a grid weight.** `raw/cell_area.bin` is the MESH region's dual area and is the
  right weight for reducing mesh regions onto a grid cell, which is what
  `lib/gridding.py`'s `cell_*` operators do with it; averaged onto a grid it
  carries no information about the cell and is not a grid weight.
  `grid/grid_cell_area.bin` is the area of the grid CELL and is the weight for
  area-averaging a gridded field. `grid/gauss_weights.bin` is that same
  partition on a Gaussian grid expressed as weights summing to 2, and it is
  what the spectral model's own global budget is taken over, so on those grids
  the two are one fact: `grid_cell_area = R² × Δlon × gauss_weight`.
  `cos(lat)` is a fourth thing and is none of them; it is not a partition of the
  sphere at all.
  A Gaussian row is a quadrature abscissa and not a cell centre, so a partition
  built on midpoints between rows is a DIFFERENT partition of the sphere. It
  closes to 4πR² just as exactly and it reports a global mean the model does not
  take: it is wider than the quadrature interval by 22 per cent in the polar row
  at every truncation, and that does not shrink with resolution. Exports built
  before 2026-08-25 carry that partition in `grid/grid_cell_area.bin` and in
  `planet.nc`'s `grid_cell_area`; on those, take `grid/gauss_weights.bin`, or
  `lib/gridding.py:gaussian_grid(nlat, nlon).cell_area(radius_m)`, which
  constructs it and never reads the field.
  (Earlier exports called the gridded field `cell_area` and it was a mesh
  diagnostic, not a cell area -- that trap is fixed, but any code written
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
# ONE invocation, one generation, six grids. The flags before the first --grid
# describe the PLANET and are global; each --grid opens an output target and the
# export options after it belong to that target alone, which is how T42 keeps
# raw/ while the rest do not.
#
# The planet code carries the seed and every terrain slider. It does NOT carry
# radius, gravity, or lithology strength, so those three are passed explicitly
# and are the ones to check against config/planet.yaml.
#
# --basin-min-area is passed for a DIFFERENT reason: the code DOES carry a basin
# slider, and it outranks terrain-config.js unless a basin flag is given. Its
# ladder has no 850 rung, so the fork's floor can only be stated here. The
# exporter refuses rather than silently taking the code's 1000; see below.
cd vendor/orogen
node --max-old-space-size=32000 tools/export-planet.mjs \
    --code 01eshm059lt0b9mpgro2y83t \
    --radius 7645.2 --gravity 12.81 \
    --lithology-strength 0.682 --regions 10000004 \
    --basin-min-area 850 \
    --netcdf --quiet \
    --grid T42     --out ../../source/<build>/exoplasim-T42 \
    --grid T21     --out ../../source/<build>/exoplasim-T21  --no-raw \
    --grid T85     --out ../../source/<build>/exoplasim-T85  --no-raw \
    --grid T127    --out ../../source/<build>/exoplasim-T127 --no-raw \
    --grid T170    --out ../../source/<build>/exoplasim-T170 --no-raw \
    --grid 512x256 --out ../../source/<build>/grid-512x256   --no-raw
```

**An output directory that already holds an export is refused.** The exporter
checks every `--out` for a `manifest.json` or a `planet.zip` before the
generation starts and stops with exit 2, naming the directory. `source/` is
read-only by rule 7 and git tracks none of this payload, so a directory
rewritten in place has nothing to fall back on -- and the rewrite need not
announce itself, since re-exporting the same code onto a different grid leaves a
directory whose terrain hash still matches every reader's record while the bytes
have moved. There is no override flag. Export to a NEW build name; if a
directory holds a failed export, remove it yourself, having looked at what is in
it.

**One invocation, because the grid does not decide the terrain.** Six exports of
one build used to run six generations of the same planet, and generation is
almost all of what an export costs: over the six-grid run of 2026-08-26 the grid
that also wrote a 4.6 GB `raw/` payload was the FASTEST of the six, so writing
the export is the small remainder. Every target of one invocation therefore
carries the same `manifest.hashes.finalElevation` by construction rather than by
repetition, and `--ice-mask` and `--preserve-basins`, which are consumed at
generation, shape the single terrain all six are written from.

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
area floor governs. And the cost is superlinear: 1,093 s and 9.0 GB at 10M
against 3,875 s and 21.1 GB at 25M.

**Those four counts were measured under the 1000 km2 floor and the floor is now
850**, which moves the crossover from 8.81M regions to 10.37M. So at 10M the
cell floor is once again the binding one and the convergence argument above is
evidence about the OLD floor: the counts under 850 are higher -- about 1,105
higher at 10M -- and whether they still converge across the same four region
counts is not measured. Re-take it from the next generation's catalogues rather
than carrying these numbers across the change. What does NOT move is the
mesh argument the paragraph rests on, which is about semivariance and erosion
texture and has no floor in it.

At 10,000,004 the six-grid build ran in 2h23m as six generations, each export
between 22 and 26 minutes and about 9 GB. As one generation it is bounded below
by 33 minutes -- six exports at total T with a shared generation G cost T - 5G,
and G cannot exceed the fastest single export -- and the measured split at
250,001 regions, where generation is 91 per cent of an export, puts it near 45.
The heap stays at 32000: the bundle for one target is released before the next
is built, so the peak is the mesh plus the largest single target, which measured
20 per cent above a per-grid invocation's peak at 250,001 regions rather than
six times it. Then register the terrain hash in `lib/orogen.py` and point
`source_build` at it.

Add `--preserve-basins FILE` to apply a carve list; without it the build is
pre-carve, which is what a first pass on new geography wants.

### Ice is supplied, not inferred

Glacial erosion carves where `vendor/orogen/js/glacial-ice.js` puts the ice, and
its built-in placement is a latitude-and-elevation ramp that consults no
temperature at all -- not the spectrum, not the obliquity, not the rotation. At
the code's `glacialErosion` the ice line is put there by the erosion slider.

`--ice-mask FILE` replaces that with a climatology's answer. The mask is one
float32 per mesh region in region order plus a `FILE.json` sidecar, produced by

```bash
python analysis/ice_mask_freezing_height.py \
    --climatology <the commissioned climatology> \
    --grid-export source/<build>/exoplasim-<rung> \
    --write-mask <path>/ice_mask.f32
```

`--build` names the export the MESH comes from and defaults to the configured
one; only the T42 export keeps `raw/`, so a climatology at another truncation
needs `--grid-export` to name the export whose GRID it was run on. The criterion
is the warmest-bin surface temperature lapse-corrected to each region's own
elevation, because the model forms no permanent land ice at any truncation this
project will run -- `glac` is a fact about the grid, not about the planet.

The mask is matched BY REGION INDEX. A mesh is a pure function of the seed and
the region count, so the sidecar carries both and the generator refuses a
mismatch before it builds anything. Never match this by coordinate.

So a build on new geography takes two passes, and the first is deliberately
without ice: generate with `--glacial 0`, an honest null rather than an
Earth-calibrated guess, commission that build to a baseline, write the mask from
its climatology, and regenerate with `--ice-mask`. It cannot be done in one pass
and it cannot be bolted on afterwards -- glacial, hydraulic and thermal erosion
share one iteration loop and a mid-loop priority flood cuts outlets through the
depressions glaciation makes, so a later glacial pass would leave them undrained
and move both the drainage network and the basin catalogue.

The code decodes to seed 16236323, 2,500,001 regions, 100 plates, 10 continents
and the erosion and shaping sliders. An explicit flag beats the decoded value,
which is why `--regions` above overrides the count and the rest of the code still
applies. Verify against `manifest.params` after
generating: `lithologyStrength` in particular must read 0.682 and not the 1 the
code decodes to, because it is a fork addition the code predates. Verify
`basins.selectionCriteria.minAreaKm2` too, for the opposite reason -- there the
code has an answer and it is the WRONG one.

**THE BASIN FLOOR IS THE ONE THE CODE ARGUES WITH.** `basinAreaFromSlider` maps
the code's slider onto `BASIN_AREA_LADDER_KM2`, which has rungs at 300, 1000 and
3000 and none at 850, and that value outranks `BASIN_MIN_AREA_KM2` unless one of
the basin flags is given. Taking the rung splits a generation in half:
`js/elevation.js` reads `BASIN_MIN_AREA_KM2` DIRECTLY for its resolvable test
while `js/basins.js` selects the catalogue at whatever the exporter passes, so
the mesh decides what it resolves at 850 while the catalogue is chosen at 1000.
Because basin preservation clamps the carve kernel, that reaches the finished
terrain and not only the catalogue. `export-planet.mjs` refuses the
disagreement rather than resolving it silently, and names both numbers; passing
one basin flag skips the code's whole basin block, which is safe here only
because `preserveBasins` defaults to true. That last point is the trap the
`BASIN_FLAGS` comment in the exporter is about, and it is why the list names
every basin flag rather than one.

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

An ice mask does not change that. It is a climate input, not a gravity one, and
the glacial altitude gate stays dimensionless on purpose: the relief ceiling
goes as 1/g because a crustal root fails at sigma/(rho g), and a dry-adiabatic
freezing height goes as 1/g because the adiabat is g/cp, so the two cancel and
the dimensionless form is the gravity-invariant one. Two builds differing only in
gravity stay bit-identical with or without a mask.

- The terrain hash is not sufficient identity for us. `lib/orogen.py` reads
  `gravityMS2` separately and `scripts/check_consistency.py` checks it against
  the config, because a build from another gravity would otherwise pass the
  allowlist while every vertical quantity was off by the ratio.
- **Basin ids and carve verdicts survive a gravity change.** The catalogue is
  bit-identical, and so is the PRESERVED set and the order it is enumerated in,
  which is what `basin_index` depends on. The two keys that decide those --
  `selectionDepthKm` for the depth floor and the sort volume -- are measured at
  reference gravity for exactly this reason: `selectBasins` runs during
  generation, and what it returns is notched and protected into the terrain. So
  an existing verdict replays. Whether it *should* is a separate question -- the
  climate driving the water balance moves -- but nothing forces the loop to
  restart from zero.
- Erosion also ran in model units, so a higher-gravity planet does not get
  steeper-slope collapse. The landscape's shape is identical; only the vertical
  scale changes.
- The catalogue's own `depthKm`, `volumeKm3`, spill and sink heights and
  hypsometry carry the scaling, exactly as `elevation_km` does, so a water
  balance may take its levels from one and its surface from the other. The one
  depth that does not is the generator's internal `selectionDepthKm`, which is
  the criterion the depth floor is compared against rather than a measurement,
  and which is kept out of the catalogue so that adding a field does not move
  `hashes.basinCatalogue`. `basins.resolution.minDepthComparedIn` names the
  currency, and undoing the scaling on the published spill and sink heights
  recovers the number, which is what
  `hydrography/scripts/catalogue_floor.py` does.
- `orog_mean/std/min/max` ARE physical kilometres, on the same vertical scale as
  `elevation_km`, with the hypsometric curve and the 1/g relief scaling both
  applied. This entry used to say the opposite -- raw model units, do not consume
  -- which was true of builds generated before 2026-08-16 and is true of neither
  build in `source/`. `vendor/orogen/tools/README.md` is authoritative on the
  format and records the fix; gridded `orog_mean` is bit-identical to
  `elevation_km`, which is how the correction was checked rather than read.
  WHAT THEY ARE STILL NOT is a statement about LAND: the four are taken over
  EVERY region in a cell, so a coastal cell's mean sits below its own coast and
  a peak excess measured against it is inflated by seabed rather than by terrain.
  For a land-population statistic take
  `hydrography/data/{build}/support_exoplasim-{res}.nc`, which names the
  population of every quantity it carries and holds the sub-grid hypsometry as
  well (SPAT-3).
