# CLAUDE.md

Worldbuilding project for **Vesper**, a super-Earth around a mid-K dwarf. The
geography comes from a fork of World Orogen; ExoPlaSim is one consumer of it,
and more components (biomes, hydrology, cultures) are expected alongside it.

## Layout

```
config/planet.yaml     Canonical planet/star/orbit/atmosphere parameters. Project-level.
source/                Canonical World Orogen exports. READ-ONLY — never modify or regenerate in place.
exoplasim/             The ExoPlaSim climate component (see exoplasim/README.md).
requirements.txt       Shared Python dependencies for .venv.
.venv/                 Python 3.12, already activated in this shell.
```

`config/` and `source/` are project-level and shared. Everything ExoPlaSim-specific
lives under `exoplasim/`. New components get a sibling directory and read the same
`config/planet.yaml` and `source/`.

This is **not a git repository**. There is no history to fall back on; be careful
with destructive operations, especially under `exoplasim/runs/` (20 GB of
irreplaceable model output that took days of CPU time).

## The upstream generator

The geography comes from a personal fork at
`/home/cfutro/git/planet_heightmap_generation` (World Orogen). Its
`tools/README.md` is the authoritative reference for the export format — read it
before writing anything that consumes `source/`. The fork adds, over upstream:
lithology (rock class, erodibility, scarp potential), preserved endorheic basins,
a richer export manifest, non-Earth planet parameters, and direct emission onto
Gaussian (spectral) grids.

## Reading `source/`

Four exports of the same planet, all from seed 16236323 with 2,500,001 mesh
regions, and all carrying the same `manifest.hashes.finalElevation`
(`821aa71b37a7…`). Check that hash before trusting any number quoted about the
terrain: the seed and parameters alone do not identify a build, because fixes to
the generator change the terrain under a fixed seed. This build followed a fix
to an over-erosion bug, which raised mean land elevation from 138 m to **548 m**
and the highest point from 4.5 km to 5.8 km. Anything derived from the previous
`a826bd12…` terrain is not comparable.

| Directory | Grid | Notes |
| --- | --- | --- |
| `exoplasim-T42/` | 128×64 Gaussian | the only one with `raw/` (native 2.5M-region mesh) |
| `exoplasim-T63/` | 192×96 Gaussian | grid only |
| `exoplasim-T85/` | 256×128 Gaussian | grid only |
| `grid-512x256/` | 512×256 uniform | grid only; for mapping/visualisation |

Each has `manifest.json` (the field catalogue — path, dtype, shape, units,
description, plus semantic blocks for lithology, basins, hydrology, plates),
`README.txt`, `grid/<field>.bin`, and `planet.nc`. The manifest is ~66 MB, mostly
the basin hypsometry catalogue; index into it rather than dumping it.

`planet.nc` is the easy path (105 fields, CF-1.8). Use the `.bin` files when you
need `raw/` or `gauss_weights.bin`.

### Choosing a land mask — use `surface_class`, and only `surface_class`

Two land definitions exist and they disagree by **1.9% of the planet's surface**:

| Source | Definition | Land area |
| --- | --- | --- |
| `surface_class == 1` | **authoritative**; ocean means connected to the world ocean | 43.17% |
| `land_mask` | elevation-sign test, `elevation_km > 0` | 41.26% |

The 48,092 disagreeing regions are dry closed-basin floor lying below sea level,
down to −5.6 km, which `land_mask` would flood. Preserving that terrain is the
whole point of the fork. Both fields now carry descriptions in the manifest
saying so and pointing at each other, and `manifest.landSeaMask` states the
disagreement in numbers, so this is checkable rather than folklore.

**Do not reconstruct it as `land_mask | is_endorheic`.** Only 43,554 of the
48,092 sit inside a preserved basin. The other 4,538 are smaller enclosed
depressions that `fixupTopology` kept as genuinely not-sea but that never cleared
the basin selection thresholds, so they are absent from the basin catalogue and
unflagged by `is_endorheic` (`basin_index == -1` for exactly those 4,538).
Verified: the naive union misses them and lands at 43.00% instead of 43.17% —
0.174% of the planet silently flooded. `surface_class` is the only correct
source.

`surface_class == 2` (`inland_water`) is **empty**, by design rather than
oversight. Orogen measures basin geometry but never decides water levels — that
is a precipitation-versus-evaporation balance and belongs downstream. 11.5% of
the planet's area is flagged `is_endorheic`. Filling those basins is our job,
using the hypsometry curves in `manifest.basins.preserved[]`, and it feeds back
into climate through albedo and evaporation.

### Masks in `source/maps/`

Use **`orogen-surfacemask-*.png`**: white land, grey inland water, black ocean,
so thresholding at `> 0` gives the correct binary land/sea mask. Verified against
the mesh — 0.431729 against the manifest's 0.431736. It currently holds only two
levels, because no water levels have been assigned yet; grey appears once they
are.

`orogen-landmask-*.png` is retained with its old `elevation > 0` meaning rather
than being silently redefined, so it is the *wrong* mask for land/sea. Note that
`convert_orogen.py` reads it, which means every existing ExoPlaSim boundary
condition uses the flooding convention. Whatever replaces that script should take
`surface_class` from `planet.nc`, or the surfacemask PNG if it stays
image-based.

### OPEN BUG: `elevation_km` is 10x too deep on dry basin floors

Reported upstream 2026-08-14; check whether the export predates the fix before
trusting `elevation_km` anywhere below sea level.

`elevToHeightKm` picks its land or ocean branch with the same `elevation > 0`
test that `land_mask` uses. The mask semantics were fixed; the elevation
conversion was not. So the 48,092 dry closed-basin-floor regions get the ocean
scaling of 10 km per unit instead of the land curve, and read exactly ten times
too deep. Bit-exact: `elevation_km == elevation * 10` for every one of them, with
`max |elevation_km - 10*elevation| = 0`.

The deepest basin sink reads −5.624 km when `elevation_pre_conditioning` and the
basin catalogue both say −0.562 km. `manifest.basins.preserved[].sinkElevationKm`
is **correct**; the gridded and raw `elevation_km` fields are not. 94,392 regions
(2.8e7 km², 3.8% of the planet) fall below their own basin's catalogued sink.
`manifest.landSeaMask.deepestDisagreeingElevationKm` is the same symptom and
should read −0.562.

It propagates to `orogen-heightmap-*.png`, where floors around −500 m render at
−4,300 to −4,900 m and 239,989 pixels clamp at the −5,000 m encoding floor. The
land-only heightmap is unaffected, since it only encodes positive elevations.

Two consequences worth holding onto: any hypsometry-driven water balance built on
`elevation_km` is wrong for precisely the terrain the fork exists to preserve,
and any topography handed to ExoPlaSim from `surface_class` plus `elevation_km`
would place −5 km pits in the middle of continents. Until it is fixed, derive
below-sea-level land depth from `elevation_pre_conditioning`, or from `elevation`
directly with the land branch.

### Other gotchas

- **`elevation` is not kilometres.** It is the generator's internal shaping
  parameter (nonlinear hypsometric curve; 0.5 ≈ 1.1 km, 1.0 = 6 km; ocean linear
  at 10 km/unit). Use `elevation_km` for physical orography.
- **Weighting.** `grid_cell_area` in the gridded output is the true cell area and
  sums exactly to 4πR², so it is a correct area weight on every grid.
  `grid/gauss_weights.bin` is equivalent on the Gaussian grids and `cos(lat)` on
  the uniform one; all three agree to ~2e-4. `raw/cell_area.bin` is the mesh
  region area. (Earlier exports called the gridded field `cell_area` and it was
  a mesh diagnostic, not a cell area — that trap is fixed, but any code written
  against an older export needs checking.)
- **Distance fields are in cell hops, not km.** Convert with
  `avgEdgeKm = π × 6371 / √numRegions` (`manifest.basins.resolution.avgEdgeKm`
  has it computed for this planet: 15.19 km).
- **`manifest.planet` rotation/obliquity/eccentricity are Earth defaults**
  (23.93 h, 23.44°, 0.0167), not this world's. Orogen does not consume them, so
  they were never overridden. **`config/planet.yaml` is authoritative** for
  everything except radius and gravity.
- `manifest.planetRadiusKm` at the top level says 6371; `manifest.planet.radiusKm`
  says the real 7645.2. Trust the latter.

### `source/maps/`

The four 16384×8192 equirectangular PNGs match the exports: same planet code
`01eshm059lt0b9mpgro2y83t`, 2.5M regions, and the correct radius and gravity.
Checked — the land mask agrees with the raw mesh to five decimal places
(0.412635 vs 0.412640) and the land heightmap's maximum is exactly the mesh's
5.769 km. `source/maps/manifest.json` records the code and planet parameters, and
notes that a planet code encodes sliders only, never radius or gravity.

To re-render them:

```bash
cd /home/cfutro/git/planet_heightmap_generation
node --max-old-space-size=12288 tools/export-maps.mjs \
    --code 01eshm059lt0b9mpgro2y83t \
    --radius 7645.2 --gravity 10.1989 --width 16384
```

## Planet parameters

`config/planet.yaml` is the single source of truth. Radius 1.2 R⊕, gravity
10.1989 m/s² (1.04 g⊕, implying 1.4976 M⊕), 30 h rotation, 32° obliquity,
e = 0.02, around a K2.5V star (0.80 M☉, 0.3236 L☉, 4965 K). 1 bar atmosphere at
450 ppm CO₂.

**Gravity is declared, not derived**, and Orogen's value is canonical. It scaled
this terrain's maximum relief as 1/g, so the geography in `source/` cannot be
separated from it; mass is what follows. Earlier revisions declared 1.50 M⊕ and
derived 10.2153 m/s², disagreeing with the geography by 0.16%. `derive()` now
reads `planet.gravity_m_s2` and raises if `planet.mass_earth` is inconsistent
with it, so the two cannot drift apart again.

The `model:` block in `planet.yaml` is ExoPlaSim-specific (resolution, layers,
timestep, output cadence). The rest is world-level. Do not split the file
casually: `config_sha256` in every `exoplasim/runs/*/run_manifest.json` pins its
exact contents, and `continue_exoplasim.py` refuses to resume a run whose config
hash has changed. The gravity change already broke that seal — no run under
`exoplasim/runs/` can be resumed, which is moot because all of them predate the
current geography anyway.

## The ExoPlaSim component

See `exoplasim/README.md` for the workflow and results, and
`exoplasim/notes/parameter-decisions.md` for every physical and format decision
(ExoPlaSim 3.4.2 calendar bugs, postprocessor code quirks, convergence criteria,
Köppen rate-normalisation, sign conventions). Read the notes before changing
anything about how runs are configured — most of the non-obvious choices are
already justified there.

**The existing climate results predate the current geography.** All runs under
`exoplasim/runs/` were built from `exoplasim/inputs/t42/*.sra`, which
`convert_orogen.py` derived from the *old* map PNGs — a different, coarser build
with a different land fraction. The completed flux sweep and stellar-cycle
experiments are physically valid but describe the old geography. Regenerating
boundary conditions from `source/exoplasim-T42/` invalidates them for
comparison.

`convert_orogen.py` still reads the PNGs and does its own conservative regridding
to a Gaussian grid. That entire step is now redundant: the fork emits T42/T63/T85
Gaussian grids directly off the mesh, with quadrature weights and sub-grid
orography, removing the lossy equirectangular intermediate. Replacing it with a
direct `source/exoplasim-T42/planet.nc` → SRA writer is the obvious next task.

Scripts anchor their paths in `exoplasim/scripts/_paths.py` and resolve from the
file location, not the working directory, so they can be run from anywhere:

```bash
python exoplasim/scripts/assess_convergence.py \
    exoplasim/runs/t42l10p8_s090_co20450ppm_rot30h_obl32_e020
```

## Environment

`.venv` is already active. Dependencies are pinned in `requirements.txt`; install
with `UV_CACHE_DIR=/tmp/world-uv-cache uv pip install --python .venv/bin/python -r requirements.txt`.
Building ExoPlaSim needs `gcc-fortran` and `openmpi` from the host (Arch).
Matplotlib is forced to `Agg` with its cache at `/tmp/world-matplotlib-cache`.

`exoplasim/patches/exoplasim-3.4.2-star-cycle.patch` adds a sinusoidal stellar-flux
cycle to `radmod.f90`. `build_star_cycle_exoplasim.sh` applies it, verifies the
pinned upstream SHA, rebuilds, copies the result to
`exoplasim/inputs/exoplasim_cycle_t42/`, and reverses the patch on exit — the
vendored ExoPlaSim tree in `.venv` is left clean.

## Conventions

- Prose in docs and reports uses ASCII punctuation and avoids em dashes; match it.
- Every run and analysis product records its provenance (config hash, input
  hashes, software versions) in JSON. Keep that up when adding steps.
- Claims about convergence and equilibration are stated with their exact criteria
  and are labelled honestly when they miss (the 0.85-flux case is called
  "quasi-equilibrated" for missing a threshold by 0.004 W/m²). Preserve that
  standard rather than rounding results into passes.
