# Data export

Headless planet generation and full-fidelity data export, for feeding
World Orogen output into ExoPlaSim, LPJ-GUESS and anything else downstream.

Everything the simulation computes is exported — not just elevation and a land
mask, but the tectonic internals: plate Euler poles and velocities, boundary
classification, stress magnitude and orientation, subduction factor, every
distance field, terrain archetype weights, and the elevation decomposition that
built the final heightmap.

## Setup

```bash
npm install delaunator
```

That's the only runtime dependency. There is no build step.

## Usage

```bash
# Defaults: 100k regions, 24 plates, 1024x512 grid, raw + grid output
node tools/export-planet.mjs --out out/world-01

# Reproducible, higher resolution, with NetCDF
node tools/export-planet.mjs --seed 12345 --regions 250000 --grid 1024x512 --netcdf

# Just what ExoPlaSim needs, on a T42-ish grid, no climate simulation
node tools/export-planet.mjs --seed 12345 --grid 128x64 --no-climate \
    --only geometry,elevation --netcdf

# See what fields exist
node tools/export-planet.mjs --list-fields

# Full option list
node tools/export-planet.mjs --help
```

The same export is available in the browser under **Export → Export Raw Data**,
which downloads a `.zip` with identical contents (no NetCDF — that's CLI only).

## Output layout

```
out/world-01/
  manifest.json          field catalogue: path, dtype, shape, units, description
  README.txt             short reading guide, generated per export
  raw/<field>.bin        one value per mesh region
  raw/mesh_*.bin         mesh topology (triangles, halfedges, CSR adjacency)
  grid/<field>.bin       equirectangular raster, row-major
  grid/lat.bin, lon.bin  grid cell-centre coordinates (float64)
  planet.nc              optional, --netcdf: gridded fields as NetCDF classic
```

All `.bin` files are flat little-endian arrays with **no header**. The manifest
is the only thing you need to read them.

```python
import json, numpy as np

m = json.load(open('manifest.json'))
f = next(x for x in m['grid']['fields'] if x['name'] == 'elevation')
elev = np.fromfile(f['path'], dtype=f['dtype']).reshape(f['shape'])
# elev[0] is the northernmost row; elev[:, 0] is longitude -180.
```

With `--netcdf` you can skip all that:

```python
import xarray as xr
ds = xr.open_dataset('planet.nc')
ds.elevation.plot()
```

## Conventions

**Coordinates.** `lat = asin(y)`, `lon = atan2(x, z)` from the region unit
vector. The grid is cell-centre registered: row 0 is the northernmost row, and
column 0 is centred at longitude −180 + half a cell.

**Elevation** is in kilometres relative to sea level.

**There are two land definitions and they are not interchangeable.**

| definition | field | what it is |
|---|---|---|
| `elevation > 0` | `land_mask` | the erosion domain — the cells the stream-power solve touches |
| not ocean-connected, not under a lake | `surface_class == 1` | the **subaerial surface**, which includes dry closed-basin floors below sea level |

They differ by the dry sub-sea-level basin floors: 1.91% of Vesper, and those
cells are 100% endorheic and disproportionately evaporite. Dividing by the wrong
one is a real bias, not a rounding difference — the endorheic share of land is
76.2% against `surface_class` and 79.7% against `land_mask`.

**Every published fraction names its denominator**: `basins.drainageConsistency`
carries `landDenominator` / `landAreaKm2` / `fractionOfLand`, and
`lithology.compositionDenominator` does the same. `manifest.landSeaMask`
quantifies the gap. If you compute your own fraction, use `surface_class` unless
you specifically mean the erosion domain.

**Resampling.** Continuous fields are area-weighted means of every region whose
centre falls in the grid cell, falling back to the containing region for cells
no centre lands in. This matters when downsampling: a 250k-region mesh onto a
128×64 grid puts ~30 regions in each cell, and nearest-sampling would discard
all but one of them. Categorical fields (plate index, boundary type, Köppen
class, masks) are always nearest-region sampled — averaging a class index is
meaningless. `--grid-method nearest` forces nearest everywhere.

**Raw is lossless.** `raw/` is the simulation's native irregular Delaunay sphere
mesh. Region positions are in `raw/lat.bin` / `raw/lon.bin` / `raw/x|y|z.bin`,
per-cell areas in `raw/cell_area.bin` (km², summing to the sphere's area to
within ~0.1%), and the topology files let you rebuild the dual mesh and walk
neighbours. Use this if you want to do your own conservative regridding.

**Integer types in NetCDF** are all widened to `NC_INT`. NetCDF-3's `NC_BYTE` is
signed with no portable unsigned marker, and byte/short variables need trailing
padding for odd element counts that some readers mishandle. A few extra bytes on
mask fields is worth not having to think about it.

## Field groups

`--only` takes a comma-separated list of these:

| Group | Contents |
| --- | --- |
| `geometry` | lat, lon, x/y/z, cell area |
| `elevation` | final elevation, pre-erosion elevation, land mask |
| `plates` | plate index/id, ocean flag, super-plate index, surface velocity, plate-physics diagnostics |
| `tectonics` | stress + direction, subduction factor, boundary type, all distance fields, terrain archetype weights |
| `terrain` | the elevation decomposition — each stage's contribution to the final heightmap |
| `lithology` | surface / basement / cover rock class, cover thickness, erodibility, scarp potential, bare-rock albedo |
| `hydrology` | flow accumulation, drainage receiver, marginal-sea index |
| `orography` | sub-grid elevation statistics (grid output only) |
| `basins` | closed-basin membership, endorheic mask, surface class, inland water, pre-conditioning elevation, drainage terminals |
| `climate` | wind, ocean currents, precipitation, temperature, Köppen |

Climate fields require climate to have run (omit `--no-climate`). Super-plate
fields require `--plates >= 8`. Anything the pipeline produces that isn't in the
catalogue is still exported, under group `other`, with a generic descriptor.

## Endorheic basins

Upstream Orogen forces every land cell to drain to the sea. Three separate
passes did the forcing — `fixupTopology` raised every below-sea-level closed
depression to 5 m, and two carve-biased priority floods cut outlets through the
rest, taking hundreds of metres out of a rim. This fork preserves closed basins
instead.

**Orogen does not decide which basins are endorheic.** That is a water balance
— catchment precipitation against evaporation over the lake — and climate is
ExoPlaSim's job. What Orogen guarantees is that the *geometry* survives, and
that you get everything needed to run the balance yourself. Three things are
kept deliberately distinct:

| | what it is | who decides |
| --- | --- | --- |
| Potential basin geometry | every closed depression the terrain has | measured here, all of it exported |
| Preserved endorheic terrain | the subset whose rims the drainage conditioning may not breach | a signal-vs-noise threshold, here |
| Actual lake water level | how much water is in a basin | **you**, or a downstream model — never inferred |

The only threshold Orogen applies is *is this depression a landform the mesh
resolves, or is it noise?* — `--basin-min-depth` (50 m default, above the
detail-noise amplitude), `--basin-min-area` (1,000 km² default, an absolute
physical number so the same planet keeps the same basins at any region count),
and `--basin-min-cells` (12, which only binds at coarse resolutions).

```bash
# what depressions does this planet have?
node tools/export-planet.mjs --seed 12345 --list-basins

# keep a basin that falls under the size floors
node tools/export-planet.mjs --seed 12345 --preserve-basin orogen-basin-rgys-6fp8pp

# fill two of them to specific surface elevations, in km
node tools/export-planet.mjs --seed 12345 \
    --water-level orogen-basin-rgys-6fp8pp=0.10 \
    --water-level orogen-basin-rog9-j2e2kp=0.05

# vanilla behaviour: everything drains to the ocean
node tools/export-planet.mjs --seed 12345 --no-basins
```

`manifest.json` carries a `basins` section with, per preserved basin: the
`hypsometry` curve (level → flooded area → volume, 24 samples from sink to
spill), `finalCatchment.areaKm2` measured on the *finished* terrain rather than
the pre-conditioning one and routed through filled depressions so every land
cell reaches a terminal, the spill elevation and what it spills into, and
`finalPreserved.retainedFraction` — how much of the basin's natural depth
survived erosion. On a typical planet preserved basins retain ~75–80% of their
natural spill depth against ~3% with `--no-basins`.

Every `...Km` / `...Km3` key in that section is on the same vertical scale as
`elevation_km`, relief scaling included, so a water balance can take its levels
from the basin record and its surface from the grid without mixing scales.
**Builds before 2026-08-16 got this wrong** on non-Earth-gravity planets: the
basin conversions skipped the 1/g scaling while `elevation_km` applied it, so
the two were apart by `1/reliefScale` — 4% at 1.04 g, 31% at 1.31 g. Unsuffixed
keys remain the dimensionless model parameter and are unaffected, as is
`retainedFraction`, which is a ratio of model-unit depths.

### Which sink to group by

Three exported fields describe one routing tree: `drain_to` (immediate
receiver), `drainage_terminal` (where the path ends) and `flow_accumulation`
(area passing through). All three are steepest descent on the depression-filled
surface, so they cross filled flats and every land cell resolves to a terminal.

Group `drainage_terminal` by **`finalCatchment.rootSink`**, which equals
`finalPreserved.sink`. That reproduces `finalCatchment.areaKm2` exactly, and
`manifest.basins.drainageConsistency` publishes the check (`disagreeing: 0`,
`worstRelativeError`, and both totals) so a consumer can verify it rather than
take it on trust.

Do **not** group by the top-level `sink`. That is the catalogue sink, measured on
the pre-conditioning surface, and erosion moves a basin's low point: on a
2.5M-region planet 878 of 3629 preserved basins ended with a deeper cell
elsewhere in the member set. Both sinks are published — `rootSink` for routing,
`catalogueSink`/`sink` for identity — because they answer different questions.

Basin ids look like `orogen-basin-r<sink>-<hash>` and are stable across runs of
the same seed and parameters. They are **not** stable across region counts,
because a region index is a property of the mesh — match basins across
resolutions by centroid and area instead.

Ordinary hydraulic, thermal and glacial erosion still act on preserved basins;
only the drainage-enforcement carving is held off. Sediment reaching a preserved
sink leaves the elevation model the same way sediment reaching the sea does —
real endorheic basins do aggrade, but the deposition rule has no sediment budget
or timescale, so modelling infill is left to whoever has one, using the exported
hypsometry.

## Map images

`tools/export-maps.mjs` renders the same equirectangular PNGs the browser's
Export Map button produces, headlessly. It rasterizes the same dual-mesh
triangles with the same colour functions and the same antimeridian handling —
barycentric interpolation reproduces the GPU's Gouraud shading on the smooth
heightmap types.

```bash
# All four standard maps at full size, from a planet code
node --max-old-space-size=12288 tools/export-maps.mjs \
    --code 09oa7lj8kek17v5639gvhe --width 16384

# Just a couple, from a seed
node tools/export-maps.mjs --seed 12345 --regions 250000 --width 4096 \
    --types color,landmask
```

`--code` supplies the seed and every slider, so a code round-trips to the same
planet the browser would build. Types: `color`, `biome`, `satellite`-style
`koppen`, `heightmap`, `landheightmap`, `landmask` (`biome` and `koppen` need
climate, so omit `--no-climate`-style flags — they turn climate on
automatically). Filenames match the browser's: `orogen-<type>-<name>.png`.

Output formats match too: RGBA-8 for the colour types and the land mask,
16-bit grayscale for the two heightmaps. A 16384×8192 set takes roughly 90
seconds end to end and needs `--max-old-space-size` raised, since the raw RGBA
raster alone is 537 MB.

**A map of a carved planet needs the carve list.** `--preserve-basins FILE` and
`--carve-basins FILE` take the same drainage hypothesis `export-planet.mjs` does
and mean the same thing. A planet code cannot encode one, so rendering a code
without the list its data export was built from silently draws a different
terrain — same filename, same code, different planet. The map manifest closes
that hole: `terrain.finalElevation` is the same hash the data export publishes,
so a map can be *checked* against the export it belongs with rather than assumed
to match it, and `terrain.basins` records the hypothesis (`selectionSource`,
counts, and a `listHash` identifying the exact list). Filenames still come from
the planet code, which the list is not in — pass `--name` and the tool says so.

Rendering is CPU-only and faithful to the GPU path: against a browser-exported
reference at 16384×8192, `--no-basins` output matches to a mean absolute
difference of 5 parts in 65535, the residual being rasterizer fill-rule and
interpolation differences at cell edges. Add `--no-basins` when you want to
reproduce upstream drainage; without it the terrain reflects preserved
endorheic basins and will legitimately differ.

## Lithology

Rock type is a product of tectonic history, so it is derived here rather than
downstream — every input already exists in the tectonics bundle, and a separate
tool would have to reimplement Orogen's own vocabulary to get at it. More
importantly, erosion can then respond to it.

**Two layers.** A crystalline or oceanic `basement` from deep tectonics, and a
sedimentary or volcanic `cover` of some thickness on top. Surface rock is the
cover where it survives and the basement where erosion has stripped it, so
orogen cores exhume and cratons strip to shield without needing a stratigraphic
column — which this model could not support anyway, having no time axis.

**Closed basins are zoned, not uniform.** A closed basin traps both a dissolved
and a clastic load, and they do not settle in the same place: the clastic load
drops at the margin as playa mud and alluvial fans, while the dissolved load
travels to the lowest ground and precipitates as it evaporates. So basin fill is
two classes — `evaporite` (salt crust) within
`LITHO_SALT_CRUST_DEPTH_FRAC` of the basin's own relief above its sink, and
`playa_clastic` over the rest. On Vesper that is ~10% salt to ~90% clastic by area.

This matters far beyond the rock label. `evaporite` carries albedo 0.50, near the
brightest natural land surface there is, and it used to be applied to the whole
basin under the name "Evaporite / playa fill". On a planet with a fifth of its
land in closed basins that single constant was worth **~2.8 W/m² per 0.10 of
error** in the downstream energy balance — larger than most terms anyone was
modelling deliberately. Zoning replaces it with two albedos over a measured area
split: effective 0.32 rather than 0.50.

The zoning knob is a fraction of each basin's own relief, so it means the same
thing at any region count and on any planet, and it is a geometric proxy for
flooding frequency — *not* a water balance. How often a basin actually floods is
climate. Where the low ground is, is terrain.

**A closed basin inside an orogen keeps its fill.** Cover thins toward zero as
fold-belt intensity rises, because an active orogen is being exhumed rather than
buried — but that rule does not apply to the closed-basin branch, which is where
the exhumed material goes. Altiplano, Qaidam and Tarim hold kilometres of fill
inside the most actively deforming belts on Earth. Orogen applied the thinning to
basin fill for one release, and the failure mode is worth knowing because it was
silent: the cells kept `cover_rock` = `playa_clastic`/`evaporite` while
`cover_thickness` went to 0, so `substrate_class` reported orogenic basement and any
check on the cover class still passed. It cost 15% of preserved-basin area on
Vesper and 0.0027 of mean land rock albedo, which ExoPlaSim integrates directly.

**And a closed basin on oceanic crust keeps its fill too.** Same reasoning, other
branch: the fill is the youngest deposit, so it is the cover, and the marine
parent material stays in `basement` as `morb`. Nothing is lost, because this is a
two-layer model. Before this, 326 preserved basins on oceanic crust surfaced as
pelagic ooze and carbonate platform and 163 of them — wholly oceanic — held no
basin fill at all, which is not a thing an endorheic basin can be. It is only
0.94% of Vesper's land and +0.0004 of mean land rock albedo. The reason it was
worth fixing anyway is that downstream selects the barren surface by **class
name** (`evaporite` / `playa_clastic`), so those basins were being handed to the
climate model as vegetated land — a far larger step than the albedo numbers
suggest, and one decided by which branch happened to fire.

**And on a flood-basalt province, for the third time.** 33 basins, 629,000 km².
By then the model was answering the same question two ways depending on crust
type — a basin on *oceanic* LIP got fill, on *continental* LIP got basalt — which
is what settled it.

**So the order is now a table, not a chain.** `COVER_SEQUENCE` in
`js/lithology.js` lists cover assignments youngest-first, and `classifyLithology`
walks it: the table *is* the chain, so there is no second copy to drift. All three
bugs were the same one — "which deposit is on top" is a chronology, this model has
no time axis, and the answer was being taken from the order the branches were
typed in. Closed-basin fill is first, because a basin with no outlet keeps filling
after volcanism and orogeny have stopped. Adding a cover class means adding an
entry and declaring what it `produces`; a test fails if a class reaches a surface
without one. Together the three cost 15% of preserved-basin area and 0.0034 of
mean land rock albedo.

If you consume basin extent, read `basin_index` / `is_endorheic`; if you consume
the substrate, read `substrate_class`. They answer different questions, and
after these fixes they agree over closed basins: 99.8% of preserved-basin area on
Vesper is fill, and no preserved basin lacks it entirely. The residual is cells
whose veneer erosion has stripped to basement, which is the model working.

**Erosion responds to it.** The `erodibility` field is a relative stream-power
multiplier, **mean-normalised to 1 over land**. That normalisation is the whole
contract: lithology decides *where* erosion goes, never *how much* of it there
is, so the Hydraulic Erosion slider keeps the calibration it always had. Granite
and quartzite hold up ridges; shale, molasse and evaporite strip away. Cover is
stripped as material is removed, and a cell swaps to its basement erodibility
once the cover is gone.

**The normalisation is a property of the field when it is computed, not of the
field as shipped, so do not rescale the delivered field against an assumed mean
of 1.** `buildErodibility` normalises over land at the point it runs, which is
before cover stripping; every cell that later strips takes its basement value
and the land mean moves off 1. It moves by a few percent, in either direction,
and by an amount that depends on the planet, so it cannot be corrected for
either. It is also sensitive to which definition of land you use, since
`surface_class` and `elevation > 0` disagree over dry closed-basin floor below
sea level. The field ships on the relative scale the contract describes and is
meant to be consumed as it stands; a consumer that renormalises it introduces a
few percent of error for no gain.

```bash
node tools/export-planet.mjs --list-rocks            # the rock table
node tools/export-planet.mjs --seed 1 --lithology-strength 0.5   # half contrast
node tools/export-planet.mjs --seed 1 --no-lithology  # uniform erodibility
```

`--lithology-strength 0` is exactly uniform erodibility and reproduces upstream
erosion; 1 is the full contrast implied by the rock table.

### Escarpments

`scarp_potential` (0–1) marks where an escarpment belongs. It is **not
geometry** — a cliff is a sub-kilometre feature and cells here are tens of km
across, so the landform lives entirely below the grid. The field is for a
higher-resolution downstream pass or a renderer placing scarps.

Three factors multiplied: erodibility contrast across the cover/basement
interface, proximity to where that interface daylights (cover surviving on one
side, stripped on the other), and enough local gradient to hang a cliff on.
Relief is measured against land neighbours only, so the continental slope does
not light up every coastline.

Because cover is always the weaker layer in this model (sedimentary and volcanic
cover runs 0.65–3.5 on the erodibility scale; crystalline basement 0.35–0.45),
these are **stripped-edge scarps** — plateau margins and shield/cover boundaries
of the Great Escarpment type — rather than resistant-caprock cuestas, which
would need a third layer. Read it as "an escarpment belongs here", not as a
claim about which way it faces.

Typically ~13% of land is nonzero and ~1.5% exceeds 0.25.

**Why not just give talus angle a rock dependence?** Because the talus threshold
in the thermal erosion step is not an angle of repose and cannot be one. Its
units make the slider's 0.8–1.2 correspond to a real angle of 0.007°–0.011°,
which 80–86% of land neighbour pairs already exceed. Raising it would not help
either: on a 150k-region planet the median neighbour pair sits at 0.037°, the
99th percentile at 0.37°, and the steepest pair anywhere at 16.9°, while a
genuine 34° repose angle would need a model-slope of ~4300. Talus slopes are two
orders of magnitude below what a 50 km grid can represent. `scarp_potential`
works with that limit instead of against it. The full measurement is recorded in
a comment on the thermal step in `js/terrain-post.js`.

**This is parent material, not soil.** Soil is CLORPT — climate, organisms,
relief, parent material, time — and two of those arrive only after ExoPlaSim. A
soil model here would have to invent a climate. Combine this rock map with real
ExoPlaSim output in a later stage to get the soil texture LPJ-GUESS wants.

One classification has a climate flavour and is flagged as such in the manifest:
`carbonate` vs `shelf_clastic` is assigned by absolute latitude as a proxy for
sea-surface temperature. Latitude is pure geometry so it stays inside Orogen's
remit, but reclassify those cells downstream if you have real SST.

## Non-Earth planets

Orogen used to hardcode Earth. It no longer does, which matters because
ExoPlaSim exists to simulate arbitrary planets — otherwise the pipeline would be
modelling the climate of a super-Earth whose topography was built for Earth's
gravity.

```bash
node tools/export-planet.mjs --planet Vesper --radius 7400 --mass-radius 1.6,1.16 \
    --rotation 31.2 --obliquity 18.5 --grid T42 --netcdf
```

Two parameters actually change the terrain:

- **`--radius`** sets every physical length scale — cell area, mean edge length,
  basin area floors, scarp gradients, the km-per-hop conversion for the distance
  fields.
- **`--gravity`** (or `--mass-radius M,R` in Earth units) scales maximum relief
  as **1/g**. Crustal strength caps the load a mountain root can carry before it
  spreads under its own weight, so a high-gravity world has subdued topography
  and a low-gravity one has Olympus Mons.

Gravity is applied where model units become kilometres, not inside the terrain
pipeline, and that has a consequence worth knowing before you build an
allowlist: **the reproducibility hashes over the pipeline's arrays do not move
with gravity.** Two builds differing only in `--gravity` share
`finalElevation`, `preConditioningElevation`, `preErosionElevation` and
`params` (which carries the sliders — radius and gravity are in
`manifest.planet`, not there). Terrain shape, erosion and basin selection do
not respond to g at all; only the vertical scale of the km fields does.

So a downstream reader keying on `finalElevation` alone **cannot tell two
gravities apart**, and would accept a build whose every vertical quantity is on
a different scale. Check `manifest.planet.gravityMS2` alongside the hash.

Two useful corollaries: basin ids and carve verdicts survive a gravity change
untouched, so a water-balance loop does not restart; and `basinCatalogue` *does*
move, because its `...Km` fields are correctly scaled.

Everything else — rotation period, obliquity, eccentricity, solar constant — is
carried in the manifest for ExoPlaSim and **not consumed by Orogen**, because it
drives climate. The manifest's `planet.usedByOrogen` and
`planet.passedThroughForDownstream` say which is which.

Earth values are the defaults and are exact, so a run that names no planet is
bit-identical to one from before this existed.

## Spectral (Gaussian) grids

Spectral models do not use equally spaced latitudes — their rows sit at
Gauss–Legendre latitudes, which differ from equal-angle by more than a degree at
mid-latitudes. Pass a truncation instead of a size and the grid is emitted
directly onto those rows:

```bash
node tools/export-planet.mjs --grid T42 --netcdf     # 128x64 Gaussian
```

`T21 T31 T42 T63 T85 T106 T127 T170` are recognised. Quadrature weights are
written to `grid/gauss_weights.bin`. Emitting the Gaussian grid straight off the
mesh removes the lossy equirect intermediate you would otherwise have to regrid
through.

## Sub-grid orography

Resampling gives a cell its mean elevation and nothing else, but a GCM wants to
know what the terrain does *below* the grid. The `orography` group adds, per
grid cell: `orog_std` (the usual roughness input for surface drag and orographic
gravity-wave schemes), `orog_min`, `orog_max`, `orog_count`, plus
`orog_anisotropy` and `orog_angle` from the eigenvalues of the slope covariance
tensor — 0 for a perfectly ridged cell, 1 for isotropic, with the principal axis
in degrees anticlockwise from east.

Cells no region centre fell into report `orog_count` 0 and NaN statistics; that
only happens when the grid is finer than the mesh. `--no-subgrid` skips the lot.

`orog_mean`, `orog_std`, `orog_min` and `orog_max` are physical kilometres, on
the same vertical scale as `elevation_km` — hypsometric curve and the planet's
1/g relief scaling both applied. **Builds before 2026-08-16 got this wrong**:
these four were computed straight from the dimensionless model elevation while
declaring km, so `orog_mean` was bit-identical to the `elevation` field and a
4.56 km peak was reported as 1.14 "km". The error was not a constant factor, so
it cannot be corrected downstream — regenerate instead. Check
`manifest.hashes.finalElevation` against a build you trust, and note that hash
does *not* move with gravity (see below).

## Rivers and marginal seas

`flow_accumulation` is upstream contributing area in **km²**, not a cell count,
so it compares across resolutions and planets. Flow terminates at open ocean and
at preserved basin **sinks** — so a terminal lake's sink carries the whole inflow
its water balance needs, not just its own cell. `manifest.hydrology` lists the
largest river mouths with their discharge areas — the freshwater flux weighting a
coastal ocean model wants, and the answer to "where are the great rivers".

`coastalOutletCount` is **every** land cell discharging to the sea, which is a
coastline length in cells, not a river count — on Vesper that is 70,000 outlets
with a median discharge of 339 km², i.e. one cell's own area. Use
`dischargeDistribution` to pick your own threshold: 1,155 outlets exceed
10,000 km² and 35 exceed 100,000 km².

`ocean_basin_index` marks submerged cells not connected to the world ocean, and
`manifest.hydrology.oceanBasins` gives each one its **sill depth**. That number
sets the basin's character: a deep sill exchanges freely, a shallow one restricts
it into a Mediterranean-type (evaporitic) or Black Sea-type (stratified, anoxic
at depth) sea. It is not recoverable from bathymetry without this computation.

## Driving the carve decision from a water balance

Orogen measures depression geometry; whether a basin is endorheic is a water
balance and belongs downstream. Two flags let that verdict come back in:

```bash
# 1. baseline export -> compute the verdict downstream
node tools/export-planet.mjs --code <code> --grid T42 --out run1

# 2. feed the verdict back. the preserved set becomes exactly this list
node tools/export-planet.mjs --code <code> --grid T42 \
    --preserve-basins verdict.txt --out run2
```

`--preserve-basins FILE` replaces threshold selection outright — the preserved
set is exactly what the file names. `--carve-basins FILE` is subtractive:
threshold selection minus the listed basins. The explicit form is usually
clearer, since it states the hypothesis rather than a floors-minus-subtraction.

File format, one basin per line:

```
# drainage hypothesis, iteration 3
orogen-basin-r1a2-xyz              # closed: rim kept intact
orogen-basin-r9zz-abc      0.45    # overflows: rim cut by 55% of basin relief
```

A trailing number is a **retain fraction** in [0,1]. 1 keeps the rim; 0 is the
same as not listing the basin at all; between the two cuts a notch at the saddle
so a basin that only just overflows becomes a through-flowing valley with a
residual lake rather than a fully trenched one. A JSON array of ids, or of
`{id, retain}` objects, is accepted too.

### Ids are stable across the loop

Basin ids are computed from the **pre-conditioning surface**, which is produced
before any preserve or carve decision is taken. They therefore do not change
when the carve set changes — verified by test across preserved sets of 387, 67
and 14 basins on the same seed, all giving an identical catalogue hash. A
verdict list computed against one export stays valid against the next, which is
what lets the loop converge.

Check `hashes.basinCatalogue` to be sure: if it matches the export your verdicts
were computed against, the ids still mean the same thing. It changes if the
seed, generation parameters or region count change — all of which do alter the
terrain and therefore the catalogue.

`manifest.basins` echoes the hypothesis in `selectionSource`,
`drainageHypothesis` and `incision`, so an export states on its face which
drainage assumption produced it. Two exports of one seed now legitimately
differ, and the elevation hash alone will not say why.

## Reproducibility

`manifest.json` carries a `hashes` block: SHA-256 over the final elevation, the
pre-conditioning elevation, the pre-erosion elevation, the basin catalogue, the
preserved-basin id list, per-region basin membership, and the parameter set.
Same seed and parameters give the same hashes, in the browser and headless
alike, so a downstream run can assert it is looking at the planet it thinks it
is.

## Tests

```bash
node --test tools/test-basins.mjs
node --test tools/test-lithology.mjs
node --test tools/test-integration.mjs
```

26 basin tests covering non-mutating detection, basin measurement on a synthetic
crater, deterministic ids, the selection floors, the divide-protection
invariant, hypsometry monotonicity, water-level validation and clamping,
final-terrain catchments, surface classification, and end-to-end determinism.
23 lithology tests covering the rock table, classification validity, the
mean-normalisation invariant, exhumation, scarp potential (bounds, and that it
requires contrast, a daylighting interface and relief together), and the promise
that lithology redistributes erosion without rescaling it.
22 integration tests covering planetary parameters, Gauss–Legendre grids,
sub-grid statistics, hydrology, and planet-code round-tripping — including the
invariant that Earth defaults leave output bit-identical.

## Interpreting the tectonic fields

`r_boundaryType` is `0` interior, `1` convergent, `2` divergent, `3` transform
(also listed in `manifest.boundaryTypes`).

Distance fields (`dist_*`, `riftDist`, `ridgeDist`, `backArcDist`, …) are in
**cell hops**, not kilometres, because that's how the simulation uses them.
Convert with the mean edge length: `avgEdgeKm = π × 6371 / √numRegions`. The
resolution-dependent band widths the simulation used are in
`manifest.tectonicBands`, so you can recover what each band meant physically.

`manifest.plates[]` gives each plate its Euler pole, angular velocity (before
and after the physics pass, as `omega` and `omegaAfter`), centroid, area
fraction, crustal density, and the drag / slab-pull / ridge-push / mantle-flow
diagnostics. Surface velocity at a point **p** is `omega · (pole × p)`; the
per-region decomposition into east/north components is already exported as
`plate_velocity_east` / `plate_velocity_north`.

`omega` is in the generator's arbitrary angular units, not rad/Myr. Direction
and relative magnitude are meaningful; absolute rate is not.

## Files

- `export-planet.mjs` — the CLI
- `lib/netcdf-write.mjs` — minimal NetCDF classic (CDF-1 / CDF-2) writer
- `../js/pipeline.js` — the shared headless generation pipeline
- `../js/data-export.js` — the shared field registry, resampler and serializer
