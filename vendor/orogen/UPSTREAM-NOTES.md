# Upstream notes (was CLAUDE.md)

**Renamed on vendoring.** This file is the mainline Orogen project's own agent
direction. Its goals diverge substantially from what this world uses the
generator for, so it is deliberately NOT named CLAUDE.md: it must not be picked
up as instructions inside this repository. The governing instructions are
`/CLAUDE.md` at the repo root.

Kept because it documents the fork's intent, and because renaming it rather than
deleting it keeps `git subtree pull` from re-adding it silently as a new file.


## What this fork is

This is a **personal fork of World Orogen**, repurposed as the first stage of a
worldbuilding pipeline:

```
World Orogen  →  ExoPlaSim  →  LPJ-GUESS  →  (downstream worldbuilding)
```

World Orogen's job here is to produce a **planet's solid surface and its
tectonic history**, and to hand every number it computed to the next stage.
It is no longer a web toy that happens to have an export button; it is a data
generator that happens to have a viewer attached.

**This fork optimises for one user (the repo owner) and one workflow.** Upstream
World Orogen's priorities do not apply here. Do not weigh changes against
"would other users like this" — the answer is irrelevant.

## Priorities (in order)

1. **Physical accuracy and internal consistency.** Terrain, tectonics and
   geology should be as defensible as we can make them. Where upstream chose a
   cheap approximation because it "looks right", we are willing to replace it
   with something slower and more correct. Where a field has real physical
   units, expose them and document them.
2. **Completeness of data output.** Nothing the simulation computes should stay
   internal. If a stage derives a field, that field should be reachable by the
   exporter with units and a description. New simulation code must extend the
   export path in the same change.
3. **Reproducibility.** Same seed + same parameters ⇒ same planet, in the
   browser and headless, forever. Never introduce unseeded randomness, wall-clock
   dependence, or platform-dependent iteration order into the pipeline.
4. **Artistic plausibility.** Still matters — the output should look like a real
   planet. But it no longer outranks accuracy; when they conflict, accuracy wins
   and we accept a less pretty map.

## Explicit non-goals

These were upstream priorities. They are **not** constraints here, and code may
be changed in ways that degrade or remove them:

- **Instant feedback / generation speed.** A generate taking minutes is fine.
  Do not reject an approach because it is slow. Do not add caching, deferral or
  approximation *purely* for responsiveness. (Don't be gratuitously wasteful
  either — just stop treating latency as a design constraint.)
- **Climate simulation.** ExoPlaSim owns climate. The existing wind / ocean
  current / precipitation / temperature / Köppen code is kept only because it is
  already written and its outputs are harmless to export. **Do not invest in it,
  tune it, or extend it.** Do not let climate constrain a terrain or tectonics
  decision. The Earth-match tuning suite in `tuning/climate/` is frozen — no
  need to re-run `evaluate.mjs` after changes, and a regression there is not a
  blocker.
- **Mobile support.** Desktop only. Ignore the responsive bottom-sheet layout,
  touch targets, `state.isTouchDevice` branches, and pinch-to-zoom. Do not add
  mobile equivalents for new UI. Don't gratuitously break what exists, but never
  spend effort on it.
- **Onboarding, marketing and discoverability.** The tutorial modal, the What's
  New modal, the survey prompt, SEO/AISEO meta tags, JSON-LD, `llms.txt`,
  `sitemap.xml`, the hidden `<main>` crawler block — all inert. Do not update
  them when features change.
- **Ease of use for newcomers.** The UI can assume the operator knows the model.
  Dense controls, raw field names and debug layers are fine.
- **Backwards compatibility with upstream.** Planet codes, saved URLs and the
  upstream file layout may be broken when there's a reason to.

## Architecture notes specific to this fork

**One pipeline, two hosts.** `js/pipeline.js` holds the generation pipeline and
has no DOM, no Worker globals and no CDN imports. Both the browser worker
(`js/planet-worker.js`) and the headless CLI (`tools/export-planet.mjs`) call
`runGeneratePipeline()`. **Never fork the pipeline** — if a stage needs to change,
change it there so both hosts stay identical. The caller supplies Delaunator via
`setDelaunator()` (CDN in the browser, `node_modules` in Node).

**Internals are outputs.** `assignElevation()` returns a `tectonics` bundle
(`buildTectonicsBundle` in `js/elevation.js`) carrying every per-region field
stages 1–3 produce: stress and its orientation, subduction factor, boundary
classification, every distance field, terrain archetype weights, and the
resolution-dependent band widths needed to interpret the hop-count fields. The
worker forwards it, plus a flattened per-plate table (Euler poles, angular
velocity, centroid, area, density, physics diagnostics) and the super-plate
table. When you add a per-region field to any stage, add it to that bundle and
give it an entry in `FIELD_META` in `js/data-export.js`.

**The exporter is the contract.** `js/data-export.js` is DOM-free and shared by
the browser download button and the CLI. It produces:
- `raw/` — one value per mesh region, plus mesh topology (CSR adjacency,
  triangles, halfedges) so the dual mesh can be rebuilt downstream.
- `grid/` — the same fields on a regular equirectangular lat/lon grid.
  Continuous fields are **area-weighted means** of the regions in each cell
  (correct when downsampling a fine mesh to a coarse ExoPlaSim grid);
  categorical fields are nearest-region sampled.
- `manifest.json` — dtype, shape, units and description for every field.

**Endorheic basins are preserved, not decided.** `js/basins.js` detects every
closed depression on the pre-conditioning surface, and the drainage conditioning
is forbidden to breach the rims of the ones selected. Four things used to destroy
them and all four are now gated on a protection mask: `fixupTopology` (stage 13
of `assignElevation`), the two carve-biased `priorityFloodCarve` passes, and the
stream-power solve's "no cell may sit below its receiver" clamp, which raised
every endorheic sink to its uphill pit-handling receiver.

Orogen deliberately does **not** decide which basins are endorheic — that is a
water balance, and climate is ExoPlaSim's. The only threshold applied here is
signal-vs-noise ("is this a landform the mesh resolves?"), which is why the
floors are absolute physical numbers (`BASIN_MIN_*` in `terrain-config.js`)
rather than anything resolution- or count-relative. Keep these three separate,
they are separate on purpose: *potential basin geometry*, *preserved endorheic
terrain*, and *actual lake water level* (supplied by the caller, never inferred).

**The carve decision comes from downstream.** `--preserve-basins FILE` /
`--carve-basins FILE` accept a per-basin verdict with an optional retain
fraction. Basin ids are computed from the pre-conditioning surface, so they are
invariant under the carve decision — that is what makes the water-balance loop
converge, and there is a test asserting it. Do not move detection after
selection; it would break the loop.

Partial retain requires an ACTIVE cut (`inciseOutlets`), not merely a lowered
protection floor. With a terminal root the flood has no deficit at that basin, so
permitting a deeper carve changes nothing — the first implementation was inert
and measured identical to full preservation. Same failure mode as the divide
ring. If a knob does not move a measured number, it is not a feature.

Catchments route on the FILLED surface via the flood's own visitation tree.
Raw steepest descent stranded 57% of land in local pits, making
`finalCatchment.areaKm2` under half the real contributing area while being
documented as the area to integrate precipitation over. The flood tree *is*
steepest descent on the filled surface (the heap pops by fill value, so the
neighbour that first reaches a cell is its lowest-fill neighbour) and it
additionally crosses flats, where the raw gradient is zero.

`drain_to`, `drainage_terminal` and `flow_accumulation` are three views of that
one tree and come from one call. Building them separately produced two live
contradictions: `computeRivers` was handed the basin MEMBER mask as its terminal
set, so flow stopped at rims and no sink anywhere had lake inflow; and
catchments were rooted at the catalogue sink while the record published the
remeasured one, so 878 of 3629 basins accumulated nothing at their published
sink and the two products disagreed by 18% of land. Root at the sink of the
surface being measured — `opts.sinks` — and never at the catalogue sink for
anything downstream of erosion. `auditCatchmentConsistency` publishes the check
in the manifest; `summariseBasins` throws if the two records name different cells.

Basin divide protection is a per-cell **allowance below the running elevation**
(`carveAllowance`), never an absolute floor. Protection is built pre-erosion and
the carve runs during it, so an absolute floor is a pre-erosion height: where
erosion had raised a divide the carve could take it back down, invisibly to the
invariant guard, which uses the same array as its baseline. It also moved terrain
on every planet with no retain fraction set anywhere. A knob must not move the
baseline — there is a test asserting retain=1 is bit-identical to no allowance.

`priorityFloodCarve` asserts its own invariant via `assertDividesNotLowered` —
if a protected divide is ever lowered it throws rather than silently shipping a
breached rim, which is exactly how an earlier attempt at this failed. Run
`node --test tools/test-basins.mjs` after touching drainage, erosion, or basins.

**Lithology is parent material, and erosion responds to it.** `js/lithology.js`
classifies rock from the tectonics bundle as two layers — `basement` from deep
tectonics plus a `cover` veneer of some thickness — and builds an `erodibility`
field that scales stream-power K per cell. That field is **mean-normalised to 1
over land**: lithology decides where erosion goes, never how much of it there is,
so the erosion sliders keep their calibration. Erosion strips cover and swaps a
cell to basement erodibility once it is gone (`updateExhumation`). Rock classes
and their relative erodibilities live in `ROCK_CLASSES` in that module, following
the `koppen.js` precedent; the thresholds are `LITHO_*` in `terrain-config.js`.

**The cover order is a stratigraphic claim, so it is a table and not an
if/else chain.** `COVER_SEQUENCE` in `js/lithology.js` lists cover assignments
youngest-first and `classifyLithology` walks it; the table *is* the chain, not a
comment describing one. This exists because the same bug shipped three times in
that file: an orogen exhuming the closed basin inside it, a closed basin on
oceanic crust surfacing as pelagic ooze, and one on flood basalt surfacing as
basalt. Each was the question "which deposit is on top?" — a chronology, in a
model with no time axis — being answered by the order someone typed the branches
in. Together they cost 15% of preserved-basin area and 0.0034 of mean land rock
albedo, which ExoPlaSim integrates.

**Closed-basin fill is first in that sequence and that placement is the load-
bearing one.** A closed basin is defined by the condition that keeps operating
after every episodic process stops — no outlet, so whatever arrives stays.
Volcanism resurfaces a landscape once; a basin fills continuously afterwards.
Adding a cover class means adding an entry with its `produces` list, and a test
fails if a class reaches a surface undeclared. Don't answer a "what is on top"
question by moving a branch; answer it in the table, where the order is visible.

A time axis would be the general fix and is **not** warranted: these were solved
by a rule about which processes persist, not by dating anything, and real
chronology means forward-integrating tectonics and sediment budgets — a
landscape evolution model, not this one. Revisit only if a downstream stage
needs sediment thickness, subsidence history or paleogeography; nothing in
ExoPlaSim or LPJ-GUESS asks for age. The seam if it ever happens is
`updateExhumation`, which already integrates one process through erosion.

**Talus angle is not lithology-dependent, on purpose.** The `talusSlope` in the
thermal erosion step is not an angle of repose: its units make the slider's
0.8–1.2 a real angle of 0.007°–0.011°, already exceeded by 80–86% of land
neighbour pairs, so the step is linear diffusion with a vestigial floor. A real
34° repose angle needs a model-slope of ~4300 while the steepest neighbour pair
on a 150k-region planet is 16.9° — talus slopes are two orders of magnitude
below what the grid can represent. Don't "fix" this by raising the threshold or
varying it by rock; the measurement and reasoning are in a comment on the
thermal step. Where escarpments belong is reported by `scarp_potential` instead,
which is a sub-grid marker, not geometry.

**Rock albedo is consumed downstream as physics, so its constants are not free.**
`ROCK_CLASSES.albedo` feeds the climate stage's energy balance directly. The
`evaporite` value of 0.50 was picked for plausibility and applied to every cell of
every closed basin under the name "Evaporite / playa fill" — a clean-halite
reflectance over a surface that is mostly not halite. On a planet with a fifth of
its land in closed basins that was worth ~2.8 W/m² per 0.10 of error, more than
the terms downstream was modelling on purpose, and it made a converged surface
temperature conditional on a number nobody had chosen for radiative accuracy.
Basin fill is now zoned by `saltCrustMask` into `evaporite` (crust, in the sump)
and `playa_clastic` (margins), effective albedo 0.32. Before changing any albedo,
assume something downstream is integrating it, and say so rather than retuning
quietly — a silent shift moves every climate result without announcing itself.

This is **not soil** — soil is CLORPT and needs climate, so it belongs to a stage
after ExoPlaSim. Do not add one here. The single climate-flavoured call is
carbonate-vs-clastic by latitude, flagged in the manifest as overridable
downstream. Run `node --test tools/test-lithology.mjs` after touching erosion or
rock classification.

**Planets are not assumed to be Earth.** `js/planet-params.js` owns radius,
gravity and the orbital metadata. Only two of them touch terrain: `radiusKm`
sets every physical length scale (threaded as an optional argument through
`geometry.js` and its consumers — never as module state, because tests run
several pipelines in one process), and `gravityMS2` scales maximum relief as
1/g. Rotation, obliquity, eccentricity and insolation are
carried for ExoPlaSim and must NOT be wired into terrain — they drive climate.
Earth defaults are exact and the regression test asserts that an unnamed planet
is bit-identical to one from before this existed; keep it that way.

**Gravity is a unit conversion, not a terrain input, and that is a sharp edge.**
`applyFinalShaping` takes one argument, the elevation array; `elevation.js`
never sees a planet. The 1/g scaling is applied where model units become
kilometres — `elevToHeightKm(...) * reliefScale`. So the pipeline, and every
reproducibility hash over its arrays, is **gravity-independent**: two builds
differing only in gravity share `finalElevation`, `preConditioningElevation`,
`basinCatalogue` and `params` (which carries sliders, not radius or gravity).
Terrain shape, erosion and basin selection do not respond to g at all; only the
vertical scale of the km fields does. Consequences worth stating because they
have already bitten:

- A downstream terrain-hash allowlist **cannot** tell two gravities apart. Check
  `manifest.planet.gravityMS2` alongside the hash.
- Carve verdicts and basin ids survive a gravity change untouched.
- Every km-labelled quantity must apply `reliefScale`, or it silently belongs to
  a different planet than the one beside it in the manifest. Both known misses
  are fixed — the basin catalogue's `...Km` fields (threaded through
  `detectBasins`/`attachHypsometry`/`remeasureBasins`) and the `orog_*` sub-grid
  statistics (now computed from `elevation_km`, not the model array). Positive
  heights scale and negative ones do not; `basins.js heightKm` and
  `data-export.js elevKm` must stay identical.

Making g act on the terrain itself (steeper slopes collapsing under higher
gravity) would be the physically richer model and is **not** what any of this
does. It would change every hash and invalidate every existing carve verdict, so
it is a deliberate future decision, not a bug to quietly fix.

**Grids can be Gaussian.** `--grid T42` emits Gauss–Legendre latitudes directly
off the mesh rather than equirectangular, because spectral models run on those
rows and regridding through equirect is a lossy extra step. `makeGrid()` in
data-export.js owns the distinction; row binning goes through `latitudeEdges` /
`rowForLatitude` so it is correct for both.

**The planet code is data-driven now.** `encodePlanetCode` derives its packing
from the same field table `decodeFormat` reads. Hand-maintaining a parallel list
of RADICES indices is how the two silently drifted apart and dropped a slider;
don't reintroduce that. Adding a slider means extending `SLIDERS`, `RADICES` and
the `DECODE_FORMATS[BASE_LEN]` field list, bumping `BASE_LEN`, and giving the
previous length its own entry with `defaults` that reproduce the behaviour those
older planets actually had.

**Coordinate convention.** The whole app derives geography as
`lat = asin(y)`, `lon = atan2(x, z)`. This is *not* the frame
`generateFibonacciSphere` built internally (that one is z-up). Everything
user-visible and everything exported uses the y-up convention. Don't "fix" one
without fixing all of them.

**Elevation units** are km relative to sea level.

**Two land definitions, and published fractions must name which.** `land_mask`
is `elevation > 0` — the erosion domain, the cells the stream-power solve
touches. `surface_class == 1` is the subaerial surface, which additionally
includes the dry closed-basin floors below sea level (1.9% of a planet with
preserved endorheic drainage). Those cells are entirely endorheic and
disproportionately evaporite, so an elevation-sign denominator biases both the
endorheic share (76.2% → 79.7%) and `lithology.compositionLand` (evaporite
20.8% → 18.6%) — non-uniformly, and against the classes most characteristic of
the drainage regime. Use `surface_class` for anything published; use
`land_mask` only where the erosion domain is genuinely what is meant, and say so
(`buildErodibility` is the legitimate case). Every fraction in the manifest
carries its `landDenominator`.

## Key rules

**Scale invariance still applies.** The result must look and measure the same
regardless of the region count (2K to 2.5M). Never use raw cell-hop counts or
neighbour-displacement magnitudes without scaling by resolution:

- Smoothing passes target a physical distance:
  `Math.max(minPasses, Math.round(targetKm / avgEdgeKm))` where
  `avgEdgeKm = (π × 6371) / √numRegions`. Never a bare `smooth(mesh, field, 5)`.
- Multipliers on neighbour-displacement quantities must normalise by
  `avgEdgeRad = π / √numRegions`.
- BFS hop thresholds are `Math.round(targetKm / avgEdgeKm)`, not fixed integers.
- Thresholds already in physical units (degrees latitude, km, °C, mm) need no
  scaling.
- When in doubt: "if I double numRegions, does this value change meaning?"

**Every new field gets exported.** Adding a per-region array to the pipeline
without adding it to the tectonics bundle / `FIELD_META` is an incomplete change.
Undocumented fields still get exported, but under a generic descriptor — that's
a safety net, not the goal.

**Determinism.** All randomness goes through `makeRng(seed)` / `SimplexNoise(seed)`.
No `Math.random()` in the pipeline. No iteration over a `Set`/object whose order
depends on insertion history in a way that affects output — where it already
happens (plate seed sets), keep it stable.

**Never define a set as the complement of one someone else owns.** `rock !=
evaporite` meaning "barren" was correct when written, passes any schema or hash
check, and broke silently the moment the evaporite class was split — playa fill
inherited a vegetation canopy because it joined a complement by default. The
defence is not validation, it is stating the set: `LAND_HEIGHT_BRANCH_CLASSES`
lists which surface classes take the land branch rather than testing `!== ocean`,
and there is a test that fails if a surface class is added without someone
deciding its branch. Same idea as publishing `carvedByRetainZero` instead of
letting an absence imply it — say the thing rather than let it be inferred.

**Terrain tunables** live in `js/terrain-config.js`. Prefer adding a named
constant there over hardcoding. (`js/climate-config.js` is the climate
equivalent, and is frozen — see non-goals.)

**Update `README.md`** when controls, the pipeline, the export format, the CLI
flags, or the file layout change. `tools/README.md` documents the export
format and CLI; keep it accurate, it's the interface to the rest of the
pipeline. Skip the tutorial modal, the What's New modal and the SEO files.

**Planet code encoding** (`js/planet-code.js`) packs slider values into a base36
string. If a slider's range, step or count changes, update `SLIDERS`, `RADICES`,
`encodePlanetCode`/`decodePlanetCode`, and the slider wiring in `js/main.js`.
Breaking old codes is acceptable; silently mis-decoding them is not.

## Running things

```bash
# Headless generate + full data export (needs: npm i delaunator)
node tools/export-planet.mjs --seed 12345 --regions 250000 --grid 1024x512 --netcdf

# List the exportable field catalogue / this planet's closed basins
node tools/export-planet.mjs --list-fields
node tools/export-planet.mjs --seed 12345 --list-basins

# Regression suites — run after any drainage/erosion/basin/lithology change
node --test tools/test-basins.mjs
node --test tools/test-lithology.mjs
node --test tools/test-integration.mjs

# Browser (no build step)
python3 -m http.server 8000     # then open http://localhost:8000
```
