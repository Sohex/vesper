# World Orogen

A browser-based procedural planet generator that creates realistic terrestrial planets with tectonic plate simulation, elevation modeling, and interactive editing. Uses native ES modules with no build step required.

[![Live Site](https://img.shields.io/badge/Try_it-orogen.studio-brightgreen)](https://orogen.studio/) ![Three.js](https://img.shields.io/badge/Three.js-0.160.0-blue) ![No Build](https://img.shields.io/badge/build-none-green)

> **Note — this is a personal fork.** Upstream World Orogen is concept art for
> planets; this fork is the first stage of a worldbuilding data pipeline
> (World Orogen → ExoPlaSim → LPJ-GUESS). It prioritises physical accuracy,
> completeness of data output, and reproducibility over generation speed,
> onboarding and mobile support, and it treats the built-in climate simulation
> as frozen legacy — climate is ExoPlaSim's job. See [CLAUDE.md](CLAUDE.md) for
> the full priority list, and [tools/README.md](tools/README.md) for the data
> export format.

## Philosophy

The generator's job is to produce a planet's solid surface and its tectonic
history, and to hand every number it computed to the next stage. Terrain,
tectonics and geology should be as defensible as we can make them; where a
cheap approximation and a correct-but-slower model disagree, the correct one
wins. Nothing the simulation computes stays internal.

The viewer is still there, and the interactive plate editing still works — but
the product is the data, not the picture.

## Guiding Principles

1. **Physical accuracy and internal consistency** — defensible tectonics and geology, with real units where they exist.
2. **Completeness of data output** — every computed field is exportable, documented, and unit-labelled.
3. **Reproducibility** — same seed and parameters produce the same planet, in the browser and headless.
4. **Artistic plausibility** — the output should still look like a real planet, but it no longer outranks accuracy.

Explicit non-goals: generation speed, climate simulation quality, mobile
support, onboarding, and backwards compatibility with upstream.

## Features

- **Fibonacci sphere meshing** with Voronoi cell tessellation via Delaunay triangulation
- **Tectonic plate simulation** — farthest-point seed placement with top-3 jitter, round-robin flood fill with directional growth bias, growth-rate governor, compactness penalty to prevent spindly shapes, multi-pass boundary smoothing, and fragment reconnection
- **Ocean/land assignment** — farthest-point continent seeding, round-robin growth with separation guarantees, trapped sea absorption, targeting ~30% land coverage
- **Collision detection** — convergent, divergent, and transform boundary classification with density-based subduction modeling; dual-layer super plate system groups same-type plates into ~20 tectonic units for broad orogenic belts blended 50/50 with fine-grained individual plate orogeny
- **Elevation generation** — three distance fields (mountain/ocean/coastline) combined via harmonic-mean formula, stress-driven uplift, asymmetric mountain profiles, continental shelf/slope/abyss profiles, foreland basins, plateau formation, and rift valleys with graben profiles
- **Ocean floor features** — mid-ocean ridges at divergent boundaries, deep trenches at subduction zones, fracture zones at transform boundaries, back-arc basins behind subduction zones
- **Island arcs** — volcanic island chains at ocean-ocean convergent boundaries with ridged noise shaping
- **Hotspot volcanism** — dual-component mantle plume model (broad thermal swell + volcanic peak) with drift-trail island chains, domain-warped shape distortion, drift-direction elongation, summit calderas on active domes, radial rift-zone ridges, age-dependent volcanic texture, and per-hotspot variation in strength/decay/spacing
- **Terrain post-processing** — noise-based domain warping (FBM simplex noise with greedy mesh walk) to deform the elevation field for organic coastlines and mountain ridges, independently controllable bilateral smoothing to blend harsh BFS distance-field boundaries, glacial erosion that carves fjords, U-shaped valleys, and lake basins at high latitudes and altitudes via latitude-driven ice flow with drainage accumulation, priority-flood pit resolution with canyon carving (Barnes et al. algorithm that ensures every land cell drains to the ocean, carving dramatic canyons through mountain saddle points rather than filling basins), iterative implicit stream power hydraulic erosion (Braun-Willett style) that carves self-reinforcing river valleys with automatic sediment deposition in flat receivers, thermal erosion that softens ridges via talus-angle material transport, ridge sharpening that accentuates mountain ridgelines, and always-on soil creep (Laplacian diffusion) that rounds off hillslopes
- **Coastal roughening** — fractal noise with active/passive margin differentiation, domain warping for bays/headlands, and offshore island scattering
- **3D globe rendering** with atmosphere rim shader, translucent water sphere, terrain displacement, and starfield
- **Equirectangular map projection** with antimeridian wrapping
- **Interactive editing** — Ctrl-click plates to mark them for reshaping (multi-select with visual tinting), then click Rebuild to apply all changes at once. Ctrl-click again to undo a pending selection. Press Escape to cancel all pending edits
- **Seasonal wind simulation** — pressure-driven wind patterns with a longitude-varying ITCZ that tracks the thermal equator (~5° over ocean, up to 15-20° over continents), Gaussian pressure bands (subtropical highs, subpolar lows, polar highs), land/sea thermal contrast for monsoon-like pressure reversals, elevation barometric effects, and Coriolis-deflected geostrophic wind with natural cross-equatorial flow reversal. Computed for both summer and winter seasons.
- **Ocean surface currents** — rule-based geographic gyre simulation driven by wind belts (trade winds, westerlies, polar easterlies) with a longitude-varying ITCZ equatorial countercurrent. Continental shelves are classified as western or eastern boundaries via coast-normal BFS, producing subtropical gyres (CW in NH, CCW in SH) with western boundary intensification (Gulf Stream, Kuroshio effect) and weaker eastern boundary return flow. Detects circumpolar channels for unobstructed eastward currents (Antarctic Circumpolar Current). Currents are colored by heat transport: red = warm poleward flow, blue = cold equatorward flow, black = zonal (neutral). Computed for both summer and winter seasons.
- **Precipitation** — blended dual-model approach: a complex moisture advection simulation is blended with a fast heuristic zonal model (blend weight and all climate constants are tuned against real-Earth Köppen data; see `js/climate-config.js`). The advection model simulates wind-driven moisture transport from coasts with six mechanisms: ITCZ convective uplift, frontal convergence, orographic rain/shadow, lee cyclogenesis, polar-front precipitation, and subtropical high suppression. The heuristic model provides smooth latitude-based patterns (ITCZ wet belt, subtropical dry belt, mid-latitude recovery, polar dryness) modulated by continentality and orographic effects. Blending the two reduces splotchiness while preserving terrain-informed detail and strengthening subtropical desert formation (~20–35°). Visualized on a brown (dry) → green (moderate) → blue (wet) color ramp. Computed for both summer and winter seasons.
- **Map type switcher** — first-class Terrain / Satellite / Climate / Heightmap tabs with color legends for each view
- **On-demand climate** — optional deferred climate computation; skip climate during generation for faster terrain iteration, compute it on demand when needed
- **Detailed visualization** — twenty-six selectable inspection layers organized by category (Geology, Atmosphere, Ocean, Climate, Elevation) for viewing each component in isolation. Wind/pressure layers show directional wind arrows, ocean current layers show current arrows colored by heat transport, on both globe and map views. Precipitation layers use a brown→green→blue ramp showing dry to wet regions.
- **Heightmap import** — bring your own equirectangular B&W heightmap (Earth, Mars, hand-drawn maps) onto a 3D globe. Black pixels become ocean, brighter pixels become higher land. The import page (`/import`) runs full climate simulation (wind, precipitation, temperature, K&ouml;ppen) on your imported terrain, with optional terrain sculpting (smoothing, erosion, ridge sharpening). Supported formats: PNG, JPEG, WebP.
- **Non-Earth planets** — radius and surface gravity are parameters, not assumptions. Radius sets every physical length scale; gravity scales maximum relief as 1/g, so a high-gravity world has subdued topography. Rotation, obliquity, eccentricity and insolation are carried through for the climate stage.
- **Spectral grid output** — `--grid T42` emits Gauss–Legendre latitudes directly, matching what a spectral climate model runs on instead of forcing a lossy regrid from equirectangular.
- **Sub-grid orography** — per-cell elevation standard deviation, extremes, and slope anisotropy/orientation, the standard inputs for surface roughness and orographic gravity-wave drag.
- **Rivers and marginal seas** — flow accumulation as real upstream area in km², river mouths with discharge areas for freshwater flux, and sill depths for every marginal sea (which is what separates a Mediterranean from a Black Sea).
- **Lithology and differential erosion** — surface rock type derived from tectonic history as a two-layer basement/cover model, with erosion responding to rock strength. Granite and quartzite hold up ridges while shale and evaporite strip away; the erodibility field is mean-normalised so lithology redistributes erosion rather than rescaling it. Cover is stripped as terrain erodes, exhuming orogen cores and cratonic shields. Which deposit sits on top is a declared table (`COVER_SEQUENCE`), walked in youngest-first order rather than left to the order of an if/else chain — closed-basin fill outranks volcanism and orogeny, because a basin with no outlet keeps filling after the episodic processes have stopped. Also exports `scarp_potential`, marking where escarpments belong — a cliff is sub-grid at planetary cell sizes, so this is a hint for a higher-resolution pass rather than geometry. Parent material only — soil needs climate and belongs downstream.
- **Endorheic basins** — closed depressions keep their rims instead of being carved open to the sea. Depression geometry, spill hierarchy, hypsometry curves and final-terrain catchments are catalogued and exported; whether a basin actually holds water is left to the downstream water balance, and lake levels are supplied rather than inferred. See [tools/README.md](tools/README.md).
- **Full data export** — every field the simulation computes, not just the heightmap: plate Euler poles and surface velocities, boundary classification (convergent/divergent/transform), stress magnitude and orientation, subduction factor, all distance fields, terrain archetype weights, and the per-stage elevation decomposition. Exported on the native mesh (lossless, with topology) and/or resampled to an equirectangular grid, as flat binaries with a JSON manifest — plus NetCDF from the headless CLI. See [tools/README.md](tools/README.md).
- **Headless generation** — `tools/export-planet.mjs` runs the identical pipeline in Node with no browser, for scripted and reproducible runs.
- **Map export** — download high-resolution equirectangular PNGs (color terrain, satellite biome, climate/Köppen, B&W heightmap, land-only heightmap, or B&W land mask) at configurable widths up to 65536px with tiled rendering. **Export All** downloads Satellite, Climate, Heightmap, and Land Mask in one click, auto-computing climate if needed.

## Quick Start

Serve the project with any local HTTP server (required for ES modules):

```bash
# Python
python3 -m http.server 8000

# Or Node.js
npx serve .
```

Then open **http://localhost:8000** in your browser. No dependencies to install, no build step.

Click **Build New World** to create a new random planet. The button changes color and label based on what you've adjusted:
- **Build New World** (blue) — generates a fresh planet with a new random seed
- **Rebuild** (amber) — re-renders the current planet at a new detail/roughness level without changing continent shapes
- **Regenerate** (red) — creates new tectonic plates when the Plates or Continents slider has changed

### Navigation

A top navigation bar connects the two pages:
- **Generate** (`/`) — procedural planet generation with tectonic plates, erosion, and climate
- **Import** (`/import`) — import your own equirectangular B&W heightmap, view it on a 3D globe, and run climate simulation. Black (0) = ocean, brighter = higher elevation. Supports PNG, JPEG, and WebP.

### Sharing Planets

Every generated planet produces a **planet code** (shown below the Build button) that encodes the random seed, all slider values, and any plate edits. An unedited planet is 21 characters; plate edits (applied via Rebuild) extend the code to include the toggled plates. Older codes (13–18 characters) from previous versions are still supported — missing sliders default to their current default values. To share a planet:

- **Copy** the code with the copy button and send it to someone
- **Load** a code by pasting it into the planet code field and clicking Load (or pressing Enter). The Load button turns blue when a new code is ready to apply.
- **URL sharing** — the code is also stored in the URL hash (e.g. `#a7f3kq9xp2b`), so you can share the full URL directly. Opening a URL with a valid hash auto-loads that planet, including any plate edits.

## Controls

### Shape Your World

Core world parameters that control the planet's structure (changing these requires a full rebuild):

| Control | Range | Default | Description |
|---------|-------|---------|-------------|
| Detail | 5,000 – 2,560,000 | 204,000 | Number of Voronoi cells on the sphere. Only affects rendering resolution — continent shapes are stable across detail levels (generated on a fixed ~20K reference grid) |
| Irregularity | 0 – 1 | 0.75 | Randomization of Fibonacci point positions |
| Plates | 4 – 120 | 80 | Number of tectonic plates |
| Continents | 1 – 10 | 4 | Target number of separate landmasses |
| Roughness | 0 – 0.5 | 0.40 | Fractal noise magnitude for terrain roughness |
| Continent Size Variety | 0 – 1 | 0.35 | How much continent sizes vary — 0 keeps continents similar in size, 1 allows a mix of large and small landmasses |
| Land Coverage | 0 – 1 | 0.3 | Percentage of the planet covered by land. Low values create ocean worlds, high values create desert worlds. Above 40% coverage, precipitation is progressively dampened to simulate reduced oceanic moisture |

### Terrain Sculpting

Post-processing passes that refine the terrain (collapsed by default — the defaults produce good results). These do not require a full rebuild; adjusting any slider lights up the **Reapply** button at the bottom of this section — click it to reapply only the sculpting passes on the current planet.

| Control | Range | Default | Description |
|---------|-------|---------|-------------|
| Terrain Warp | 0 – 1 | 0.75 | Domain warping — deforms the elevation field using noise to produce organic, squiggly coastlines and mountain ridges |
| Smoothing | 0 – 1 | 0.10 | Blends harsh terrain boundaries from tectonic generation |
| Glacial Erosion | 0 – 1 | 0.50 | Ice-age sculpting — carves fjords, U-shaped valleys, and lake basins at high latitudes and altitudes via latitude-driven ice flow |
| Hydraulic Erosion | 0 – 1 | 0.50 | Iterative stream-power erosion — resolves endorheic basins via priority-flood canyon carving, then carves river valleys and dendritic drainage networks, with automatic sediment deposition in flat receivers |
| Thermal Erosion | 0 – 1 | 0.10 | Slope-driven material transport — softens ridges and creates natural talus slopes |
| Ridge Sharpening | 0 – 1 | 0.50 | Accentuates mountain ridgelines — pushes peaks further above their surroundings for more dramatic terrain |

### Climate

Global climate offsets that adjust temperature and precipitation without a full rebuild. Changing these triggers a fast climate-only recompute.

| Control | Range | Default | Description |
|---------|-------|---------|-------------|
| Temperature | -15 – 15 | 0 | Global temperature offset in °C — positive makes the planet warmer, negative colder. Climate zones shift accordingly |
| Precipitation | -1 – 1 | 0 | Global precipitation scale — positive makes the planet wetter, negative drier. Affects desert and rainforest distribution |

### Auto Climate

Climate simulation (wind, ocean currents, precipitation, temperature, Köppen classification) runs automatically during generation when detail is ≤ 300K regions. Above 300K, climate is skipped for faster terrain iteration and computed on demand when switching to a climate-dependent view.

### Visual Options

- **Map Type** — segmented Terrain / Satellite / Climate / Heightmap tabs for quick switching between the four most common visualizations. Each tab shows a color legend:
  - **Terrain** — elevation color ramp from deep ocean through sea level to mountain peaks
  - **Satellite** — realistic biome colors based on Köppen climate classification and elevation (lush green rainforests, tan deserts, white ice caps, dark taiga, gray tundra), with ocean using the standard terrain palette. High elevations blend toward snow white based on climate-aware snow lines.
  - **Climate** — Köppen-Geiger classification with color swatches for all 30 climate types
  - **Heightmap** — black-to-white gradient on a fixed absolute scale (-5 km ocean floor to 6 km peaks), so the same physical height always maps to the same shade
- **View** dropdown — switch between Globe and Map (equirectangular projection)
- **Center Longitude** slider (map mode only) — shifts the map projection's central meridian to any longitude from 180°W to 180°E, scrolling the equirectangular projection so the chosen longitude is centered. Exports are unaffected (always centered on 0°).
- **Wireframe** — toggle switch to show Voronoi cell edges as a wireframe overlay
- **Show Plates** — toggle switch to color regions by plate (green shades = land, blue shades = ocean); also draws black super plate boundary lines showing tectonic super-groups
- **Auto-Rotate** — toggle switch to spin the globe continuously
- **Grid Lines** — toggle switch for latitude/longitude grid overlay on both globe and map views
- **Grid Spacing** — choose the interval between grid lines: 30°, 15°, 10°, 5°, or 2.5°

### Inspect Dropdown

The **Inspect** dropdown (in Visual Options, below the map tabs) selects a detailed visualization layer. Options are organized into groups:

- **Main views** (ungrouped at top) — Terrain, Satellite, Köppen Climate, Land Heightmap
- **Geology** — Base, Tectonic, Noise, Interior, Coastal, Ocean Floor, Hotspot, Tectonic Activity, Margins, Back-Arc, Fold Ridge, Orogenic Power, Erosion Delta (blue = eroded, red = deposited)
- **Atmosphere** — Pressure Summer/Winter (blue = low, red = high), Wind Speed Summer/Winter (with directional arrows on both globe and map)
- **Ocean** — Currents Summer/Winter (red = warm poleward, blue = cold equatorward, black = zonal; with directional current arrows)
- **Climate** — Precipitation Summer/Winter (brown = dry, green = moderate, blue = wet), Rain Shadow Summer/Winter (diverging blue = windward orographic boost, gray = neutral, red-brown = leeward rain shadow; leeward effects are seeded at downslope faces scaled by mountain height, then propagated ~1500 km downwind to show extended shadow zones like the foehn drying effect), Temperature Summer/Winter (purple-blue = cold, white = 0 C, green-yellow = warm, red = hot; fixed -45 to +45 C range), Continentality (blue = ocean, green = coast, yellow = moderate interior, orange/red = deep continental interior)
- **Elevation** — Full Heightmap (full-range B&W)

### Export

Click **Export Map** (below Visual Options) to open the export modal:

- **Type** — Color Map (terrain colors), Satellite (biome colors from Köppen classification), Climate (Köppen classification colors), Heightmap (B&W full range on fixed -5 to 6 km absolute scale), Land Heightmap (B&W on fixed 0 to 6 km absolute scale, ocean is black), or Land Mask (pure B&W — white = land, black = ocean). Satellite and Climate options are disabled when climate hasn't been computed.
- **Width** slider — 1024 to 65536 pixels (height is always width/2 for equirectangular). Large exports use tiled rendering to handle GPU texture limits.
- **Export** — downloads the selected type as an equirectangular PNG with no grid overlay
- **Export All** — downloads four maps (Satellite, Climate, Land Heightmap, Land Mask) sequentially. If climate hasn't been computed yet, it runs automatically before exporting.
- A progress overlay shows rendering and PNG encoding status during export

The same modal has an **Export Raw Data** section below the map export:

- **Contents** — raw mesh + grid (default), raw mesh only, or grid only
- **Grid** — equirectangular grid width, 64 to 4096 (height is width/2)
- **Resampling** — area-weighted mean (default) or nearest cell. Categorical
  fields are always nearest-sampled regardless of this setting.
- **Export Data (.zip)** — downloads every field the simulation computed as flat
  little-endian `.bin` files plus a `manifest.json` describing each one's dtype,
  shape, units and meaning. Format details: [tools/README.md](tools/README.md).

For scripted use, larger meshes, or NetCDF output, use the headless CLI instead:

```bash
npm install delaunator
node tools/export-planet.mjs --seed 12345 --regions 250000 --grid 1024x512 --netcdf
node tools/export-planet.mjs --list-fields
node tools/export-planet.mjs --list-basins      # closed depressions in this planet
node tools/export-maps.mjs --code <planet-code> --width 16384   # PNG maps, headless
node tools/export-maps.mjs --code <planet-code> --preserve-basins carve_list.txt \
    --name my-build                             # ...of a planet built from a carve verdict
node --test tools/test-basins.mjs               # basin regression suite
node --test tools/test-glacial.mjs              # ice placement regression suite
```

Glacial erosion carves wherever `glacial-ice.js` puts the ice. `--ice-mask FILE`
supplies that placement from outside -- one float32 per mesh region in region
order, with a `FILE.json` sidecar naming the mesh it was built on -- so the ice
can be a climatology's answer rather than a latitude threshold. Without the flag
the placement is the built-in ramp, which is Earth-calibrated and consults no
temperature; the module header says what that does and does not know.

### Sidebar & Loading

The control panel can be collapsed and expanded with the **«** toggle button in the sidebar header. On small screens (≤ 768px) the sidebar becomes a bottom sheet with a drag handle — starts collapsed, showing only the handle and header. Drag up or tap the handle to expand. A fullscreen overlay with spinner, title, and progress bar appears during every generation — fully opaque on initial load, semi-transparent on subsequent builds so the previous planet is dimmed behind it. Stage labels (shaping, plates, oceans, mountains, painting) update as the pipeline progresses.

### Tutorial & Help

A five-step tutorial modal introduces the tool on first visit (auto-shown via `localStorage`). It covers planet generation, slider controls, interactive editing, visualization, saving/sharing via planet codes, and map export. A **?** help button in the top-right corner reopens the tutorial at any time. The modal can be dismissed with the close button, backdrop click, Escape key, or the "Get Started" button on the final step.

A **What's New** modal is shown once per release to returning users (those who have already dismissed the tutorial). It highlights new features, changes, and a heads-up that saved planet codes may produce different-looking worlds due to terrain/climate reworks. The modal uses a versioned `localStorage` flag (`wo-whatsnew-seen`) — bump the `VERSION` constant in `initWhatsNew()` to trigger it again on the next release.

### Interaction

Navigation hints are shown in the sidebar panel and as a contextual tooltip when hovering the planet.

| Action | Desktop | Mobile |
|--------|---------|--------|
| Rotate globe / pan map | Drag | Drag (one finger) |
| Zoom | Scroll wheel | Pinch with two fingers |
| Highlight plate + info card | Hover | — |
| Mark plate for reshaping | Ctrl-click a plate (multi-select) | Tap the edit button (pencil), then tap plates |
| Undo pending plate | Ctrl-click the same plate again | Tap the same plate again |
| Apply pending edits | Click the Rebuild button | Tap the Rebuild button |
| Cancel all pending edits | Press Escape | — |

Hovering over a region shows an info card with plate type, elevation, coordinates, and (when climate has been computed) temperature, precipitation, and K&ouml;ppen classification. Pending plates show a colored tint (green = ocean→land, blue = land→ocean) and hover text indicates "(pending)".

### Mobile Support

World Orogen is fully usable on phones and tablets:

- **Bottom-sheet sidebar** — on screens 768px or narrower, the sidebar becomes a bottom sheet with a drag handle. Drag or tap the handle to expand/collapse. The globe stays visible above.
- **Pinch-to-zoom** — two-finger pinch zooms the globe and map, using the same smooth lerp as desktop scroll-zoom.
- **View switcher** — a dropdown in the top-right lets you switch between Terrain, Satellite, Climate, and Heightmap views without opening the bottom sheet.
- **Edit-mode toggle** — a floating pencil button (bottom-right) activates plate editing. Tap it to toggle edit mode (glows green when active), then tap plates to mark them. Tap the Rebuild button to apply all changes at once.
- **Touch-friendly targets** — buttons, checkboxes, and sliders are enlarged for comfortable finger input.
- **Performance** — detail warning thresholds are lowered on touch devices (orange at 200K, red at 500K). Export widths above 8192px are disabled on mobile.
- **Tooltips** reposition above their trigger instead of to the right, so they stay on screen.
- **Orientation** changes are handled automatically.

## How It Works

### Pipeline

1. **Fibonacci spiral** distributes N points evenly on a unit sphere with optional jitter
2. **Stereographic projection** maps the sphere points to 2D
3. **Delaunator** computes Delaunay triangulation in projected space
4. **Pole closure** connects convex hull edges to a pole point, creating a watertight mesh
5. **Coarse plate generation** on a fixed ~20,000-region reference mesh (resolution-independent), via farthest-point seed placement (with top-3 jitter for variety), round-robin flood fill with per-plate growth rates, directional bias coupled inversely to growth rate, growth-rate governor, and compactness penalty
6. **Ocean/land assignment** on the coarse mesh using farthest-point continent seeding with area budgeting
7. **Plate projection** maps coarse plate assignments onto the high-res mesh via nearest-neighbor adjacency walk, then smooths boundaries with resolution-scaled majority-vote passes
8. **Collision detection** simulates plate drift to classify convergent/divergent/transform boundaries
9. **Stress propagation** diffuses collision stress inward through continental plates via frontier BFS
10. **Elevation assignment** combines distance fields, stress-driven uplift, ocean floor profiles, rift valleys, back-arc basins, hotspot volcanism, island arcs, coastal roughening, and multi-layered noise
11. **Terrain post-processing** applies domain warping (controlled by Terrain Warp slider) using FBM simplex noise to deform the elevation field for organic coastlines and mountain ridges via greedy mesh walk, then bilateral smoothing (controlled by Smoothing slider) to blend BFS banding artefacts, two always-on detail-noise passes add domain-warped high-octave FBM relief to land cells (in physical km space, with Newton-Raphson inversion of the elev→km quartic) to break up flat continental interiors before erosion routes drainage through them — Layer 1 adds 0–100 m positive bumps for broad rolling variety, then Layer 2 adds ±50 m bipolar bumps at 2× frequency and 2× warp amplitude with magnitude biased toward the extremes for sharper finer-scale relief; both layers are dampened by 50 % at full craton/basin weight so geologically quiet regions stay characteristically subdued, and both are scaled by the orogenic-power field so noise relief tracks active mountain-building zones and quiets down over tectonically inactive crust, glacial erosion (controlled by Glacial Erosion slider) carves fjords, U-shaped valleys, and lake basins at high latitudes and altitudes, priority-flood pit resolution carves canyons through mountain saddle points to ensure all land drains to the ocean, iterative implicit stream power hydraulic erosion with sediment deposition (controlled by Hydraulic Erosion slider) carves self-reinforcing river valleys, thermal erosion (controlled by Thermal Erosion slider) softens ridges via talus-angle material transport, ridge sharpening (controlled by Ridge Sharpening slider) accentuates mountain ridgelines, and always-on soil creep gently rounds off hillslopes
12. **Wind simulation** computes a longitude-varying ITCZ by scanning for the thermal maximum at each longitude (accounting for land/sea heating differential and elevation lapse rate), builds pressure fields from Gaussian zonal bands centered on the ITCZ plus land/sea thermal modifiers and elevation barometric effects, then derives wind vectors from pressure gradients with latitude-dependent Coriolis deflection and surface friction. Computed for both NH summer and winter.
13. **Ocean currents** uses a rule-based geographic approach: classifies ocean cells by wind belt (trades, westerlies, polar easterlies) to set base zonal flow, runs three BFS passes from coastal seeds to compute distance to western and eastern coastlines (classified by coast-normal direction), deflects currents poleward near western boundaries (warm, intensified ×2) and equatorward near eastern boundaries (cold, weaker ×0.8), detects circumpolar channels at ±60° latitude for unobstructed eastward flow, smooths with 5 Laplacian passes, and classifies heat transport by meridional flow direction. Computed for both seasons.
14. **Precipitation** uses a blended dual-model approach. The complex model computes moisture advection from coasts using iterative upwind propagation driven by wind vectors, with depletion based on distance and elevation gain, plus six mechanisms: ITCZ convective uplift, frontal convergence at subpolar lows, orographic rain/rain shadow, lee cyclogenesis, polar front diffuse precipitation, and seasonal subtropical high suppression (shifts poleward in local summer to create Mediterranean dry-summer patterns). A heuristic zonal model computes smooth precipitation from ITCZ distance (with aggressive subtropical drying at 15–30°), seasonal hemisphere boost with Mediterranean subtropical suppression (up to 55% summer reduction at 25-42° latitude), continental dryness, and orographic rain shadow. The two models are blended (weight in `js/climate-config.js`, tuned against real-Earth Köppen zones) then normalized via 95th-percentile scaling. Computed for both seasons.
15. **Temperature** computes per-cell surface temperature using the ITCZ as the thermal equator (~28°C peak, warmest latitude band), with poleward cooling following a power-law curve with a tropical plateau. Modulated by a zone×latitude seasonal swing table with a slight winter-heavy asymmetry, moisture-dependent elevation lapse rate (moist to dry adiabatic, interpolated by precipitation), ocean current warmth diffused onto coastal land, and precipitation/cloud cover moderation. All constants live in `js/climate-config.js` and are auto-tuned against real-Earth Köppen zones (see `tuning/climate/`). Normalized to a fixed -45 to +45 C range. Computed for both seasons.
16. **Rendering** builds a Voronoi cell mesh with per-vertex colors and terrain displacement

### Key Algorithms

- **Seeded PRNG** — Park-Miller LCG for deterministic generation
- **3D Simplex noise** — with fBm and ridged fBm variants for terrain detail
- **Harmonic-mean distance blending** — `(1/a - 1/b) / (1/a + 1/b + 1/c)` for smooth elevation transitions
- **Domain warping** — noise-driven coordinate offsets for organic coastlines
- **Density-based subduction** — tanh mapping of density differences with undulation noise
- **BFS distance fields** — randomized frontier expansion from boundary seeds, used for elevation, coast distance, rift width, ridge profiles, and back-arc basins
- **Gaussian dome uplift** — hotspot volcanism modeled as dual-component Gaussians (thermal swell + volcanic peak) with domain-warped shape distortion, anisotropic drift elongation, summit calderas, radial rift ridges, and age-dependent texture blending

## Project Structure

```
index.html              Main page — HTML markup + import map + structured data
import.html             Import page — heightmap upload + climate visualization
styles.css              All CSS (shared by both pages)
robots.txt              Search engine crawler directives
sitemap.xml             Sitemap for search engine indexing
site.webmanifest        Web app manifest (metadata + theming)
llms.txt                AI/LLM-readable site description (AISEO)
humans.txt              Project credits
CNAME                   Custom domain config (orogen.studio)
404.html                Custom 404 page
preview.png             Social preview image (og:image / Twitter card)
js/
  main.js               Generator entry point — UI wiring, animation loop
  import-main.js        Import page entry point — file upload, import dispatch
  state.js              Shared mutable application state
  generate.js           Worker dispatcher — posts jobs, handles results
  planet-worker.js      Web Worker — hosts the pipeline off the main thread
  pipeline.js           The generation pipeline itself — DOM-free, shared by the worker and the CLI
  data-export.js        Field registry, equirectangular resampler, .bin/manifest serializer, ZIP writer
  planet-params.js      Planetary radius / gravity / orbit — non-Earth worlds
  lithology.js          Rock classification (two-layer basement/cover) + erodibility field
  basins.js             Closed-basin detection, cataloguing, rim protection, catchments, water levels
  geometry.js           Sphere geometry shared by the simulation, basins and the exporter
  min-heap.js           Priority-queue shared by the drainage conditioning and basin detection
  sha256.js             Synchronous SHA-256 for the reproducibility hashes
  planet-code.js        Planet code encode/decode (seed + sliders → base36)
  rng.js                Seeded PRNG (Park-Miller LCG)
  simplex-noise.js      3D Simplex noise with fBm and ridged fBm
  color-map.js          Elevation → RGB colour mapping + satellite biome colors
  sphere-mesh.js        Fibonacci sphere, Delaunay, SphereMesh dual-mesh
  plates.js             Tectonic plate generation (farthest-point seeding, round-robin flood fill, compactness constraints)
  coarse-plates.js      Resolution-independent plate pipeline — coarse reference grid, projection, boundary smoothing
  super-plates.js       Groups same-type plates into ~20 super plates for broad orogenic belts
  ocean-land.js         Ocean/land assignment with continent seeding
  elevation.js          Collisions, stress propagation, distance fields, elevation
  glacial-ice.js        WHERE the ice that carves the terrain sits: a supplied mask from a
                        climatology, or the Earth-calibrated latitude ramp when there is none
  terrain-post.js       Domain warping, bilateral smoothing, glacial/hydraulic/thermal erosion, ridge sharpening, soil creep
  climate-config.js     Climate simulation tunable parameters (mutable at runtime for the tuning suite)
  climate-util.js       Shared climate utilities — smoothing, ITCZ lookup, percentile selection
  wind.js               Seasonal wind simulation — pressure fields, ITCZ tracking, Coriolis wind
  ocean.js              Ocean surface currents — rule-based wind-belt gyres, coast BFS, circumpolar detection
  precipitation.js      Precipitation simulation — moisture advection, ITCZ/frontal/orographic effects, blended with heuristic
  heuristic-precip.js   Heuristic zonal precipitation model — smooth latitude/continentality/orographic patterns
  temperature.js        Temperature simulation — ITCZ thermal equator, lapse rate, continentality, ocean currents
  scene.js              Three.js scene, cameras, controls, lights
  planet-mesh.js        Voronoi mesh, map projection, hover highlight
  edit-mode.js          Ctrl-click plate multi-select + hover info
  detail-scale.js       Non-linear (power-curve) detail slider mapping
tools/
  export-planet.mjs     Headless generate + full data export CLI (see tools/README.md)
  export-maps.mjs       Headless equirectangular PNG map rendering (CPU, matches the GPU path)
  test-basins.mjs       Regression suite for endorheic basin preservation (node --test)
  test-glacial.mjs      Regression suite for ice placement and the --ice-mask path (node --test)
  test-lithology.mjs    Regression suite for lithology + differential erosion (node --test)
  test-integration.mjs  Regression suite for planet params, Gaussian grids, hydrology, planet codes
  lib/netcdf-write.mjs  Minimal NetCDF classic (CDF-1/CDF-2) writer
tuning/
  climate/              Automated climate tuning suite — scores the simulated Köppen
                        map of an imported Earth heightmap against the observed
                        Köppen-Geiger classification and optimizes js/climate-config.js
                        parameters to match (see tuning/climate/README.md).
                        Frozen in this fork — climate is handled downstream.
```

## Dependencies

Loaded via CDN import maps (no installation needed):

- [Three.js](https://threejs.org/) v0.160.0 — 3D rendering
- [Delaunator](https://github.com/mapbox/delaunator) v5.0.1 — 2D Delaunay triangulation

The headless CLI (`tools/export-planet.mjs`) needs Delaunator from npm instead:

```bash
npm install delaunator
```

## License

This project is licensed under the GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.

## Acknowledgments

Inspired by [Red Blob Games' planet generation](https://www.redblobgames.com/x/1843-planet-generation/) — Fibonacci sphere meshing, dual-mesh traversal, and distance-field elevation approach.

Additional inspiration and reference from:
- [Worldbuilding Pasta](https://worldbuildingpasta.blogspot.com/) — worldbuilding science and climate reference
- [Artifexian](https://www.youtube.com/@Artifexian) — worldbuilding tutorials and planetary science inspiration
- [Madeline James](https://www.youtube.com/@MadelineJamesWorldbuilds) ([website](https://www.madelinejameswrites.com/)) — worldbuilding methodology and climate design reference
- [Fractal Philosophy](https://www.youtube.com/watch?v=7xL0udlhnqI) — procedural terrain generation inspiration
