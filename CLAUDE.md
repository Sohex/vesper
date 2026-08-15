# CLAUDE.md

Worldbuilding project for **Vesper**, a super-Earth around a mid-K dwarf. The
geography comes from a fork of World Orogen; ExoPlaSim is one consumer of it,
and more components (biomes, hydrology, cultures) are expected alongside it.

## Layout

```
config/planet.yaml     Canonical planet/star/orbit/atmosphere parameters. Project-level.
source/                Canonical World Orogen exports. READ-ONLY — never modify or regenerate in place.
exoplasim/             The ExoPlaSim climate component (see exoplasim/README.md).
hydrography/           Drainage, catchments, basin capacity (see hydrography/README.md).
requirements.txt       Shared Python dependencies for .venv.
.venv/                 Python 3.12, already activated in this shell.
```

`config/` and `source/` are project-level and shared. Component-specific work
lives under the component directory. New components get a sibling directory and
read the same `config/planet.yaml` and `source/`. `hydrography/scripts/orogen.py`
is a general reader for the export and is meant to be reused, not reimplemented.

Git tracks the scripts, notes, configuration, and analysis products. It does
**not** track `exoplasim/runs/` (20 GB of model output), the `source/` export
payloads, or `.venv/`. Those have no history to fall
back on, so be careful with destructive operations there; `.gitignore` says why
each is excluded.

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

### Elevation conventions, and a bug that is now fixed

`elevToHeightKm` used to branch on `elevation > 0`, the same test `land_mask`
uses, so dry closed-basin floors took the bathymetric branch and read ten times
too deep. Fixed upstream on 2026-08-14: the branch now takes its land flag from
`surface_class`, and below-sea-level land converts at 1.0 km per unit against
the ocean's 10. Verified here: dry floors are now exactly `elevation * 1.0`, the
deepest reads -562 m and matches the catalogue, ocean is untouched at -8.89 km.
An export whose `finalElevation` hash is not `821aa71b37a7...` predates this and
should not be trusted below sea level.

**The basin catalogue renamed keys in the same pass, and it is a breaking
change.** Unsuffixed keys (`sinkElevation`, `depth`, `spillElevation`) are the
generator's model parameter; `Km` and `Km3` keys are genuinely physical,
converted through the land branch. Volumes and hypsometry are integrated in
physical height. Read the suffixed keys unless you specifically want model units,
and do not assume an unsuffixed key is kilometres.

Note also that the top-level catalogue entry describes the **natural**
pre-conditioning basin, while `finalPreserved` describes the finished terrain.
They differ a lot: for the first basin, 53,968 km2 flooded at spill naturally
against 10,853 km2 on the finished surface. `hypsometry` is on the natural
terrain, which is why `hydrography/` recomputes it.

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
- **`lithology.compositionLand` is measured against `land_mask`**, numerator and
  denominator both, so it omits the dry sub-sea-level basin floors entirely. The
  distortion is not uniform: evaporite gains 1.17x on area against 1.046x for
  land overall, so its share is 20.8% rather than the published 18.6%. That is
  the expected direction, because playa fill accumulates in exactly the closed
  basins `land_mask` excludes. Compute composition from `surface_rock` and
  `surface_class` rather than quoting the table.

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

## The hydrography component

`hydrography/` resolves drainage over the native mesh and builds everything a
water balance needs short of the climate itself. See `hydrography/README.md`.

76% of the land drains to a closed basin against roughly 13% on Earth, but
**treat that as the hyper-arid limit rather than a fact about the world**. It
assumes no basin ever overflows, and a basin that overflows persistently incises
its outlet and stops being a basin. The endorheic share runs from 9% in a humid
climate to 76% in a hyper-arid one; `hydrography/README.md` has the curve. The
decision variable, `critical_aridity_index`, is pure geometry and lives in
`basins.nc`.

The open gap: the terrain was exported with drainage-enforcement carving
disabled, so rims survive that the water balance says are overtopped. Carving
belongs upstream in Orogen, which has the erodibility field, and the generator
has no per-basin hook for it yet: `--preserve-basin` only adds, `--no-basins`
is global, and the two retain about 80% and 3% of spill depth respectively.

Two things about the export that any consumer needs to know. `drain_to` is raw
steepest descent, and `drainage_terminal` is -2 for 63% of the land, which
drains into 220,649 unpreserved noise pits; integrating precipitation without
resolving that discards most of the land's water. And the catalogue's
`hypsometry` is on the natural terrain, so it overstates capacity by about 1.5x.
`build_hydrography.py` handles both. Use `data/basins.nc`, not the catalogue.

## The ExoPlaSim component

Two behaviours of ExoPlaSim worth knowing before changing anything here. Its
`configure()` clears every surface `.sra` when given a landmap, so all land
surface fields except topography and the land mask are uniform namelist defaults;
this is declared via `model.uniform_land_surface` and is deliberate, since Earth's
surface maps are tied to Earth's continents. Albedo is the exception and is
supplied from lithology by `build_surface_albedo.py`: bare rock averages 0.315
over land against ExoPlaSim's 0.22, worth about -12.5 W/m2, and 20.8% of the land
is bright evaporite because the drainage is endorheic. Note this is substrate
albedo, so the first pass is a bare-rock planet and runs cold. And its `finalize()` picks output as
the last glob match, so a run directory shared between worlds can silently emit
the wrong world's result; `run_id` therefore ends in a geography digest.

See `exoplasim/README.md` for the workflow and results,
`exoplasim/notes/lake-representation.md` for what the model can do with the
endorheic basins, and `exoplasim/notes/parameter-decisions.md` for every physical
and format decision
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
