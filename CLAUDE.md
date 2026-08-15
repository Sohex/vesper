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

Four exports of the same planet, all from seed 16236323 and identical hashes
(`manifest.hashes.finalElevation` is the same across all four):

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

### Gotchas that will silently give wrong numbers

- **`cell_area` in `grid/` is not the grid cell's area.** It is the mean area of
  the mesh regions that fell into the cell (~294 km² almost everywhere), so every
  row sums to roughly the same value. Using it as an area weight is wrong:
  it gives a land fraction of 0.486 instead of the correct **0.414**.
  - On the Gaussian grids weight by `grid/gauss_weights.bin`.
  - On `grid-512x256` weight by `cos(lat)`.
  - `raw/cell_area.bin` *is* a real area in km² and sums to 4πR².
- **`elevation` is not kilometres.** It is the generator's internal shaping
  parameter (nonlinear hypsometric curve; 0.5 ≈ 1.1 km, 1.0 = 6 km; ocean linear
  at 10 km/unit). Use `elevation_km` for physical orography. Land is
  `elevation > 0` in either field.
- **Land cells can sit below sea level.** Preserved endorheic basins are the
  point of the fork; ~5% of T42 land cells are below −434 m. Do not clamp.
- **Distance fields are in cell hops, not km.** Convert with
  `avgEdgeKm = π × 6371 / √numRegions` (`manifest.basins.resolution.avgEdgeKm`
  has it computed for this planet: 15.19 km).
- **`manifest.planet` rotation/obliquity/eccentricity are Earth defaults**
  (23.93 h, 23.44°, 0.0167), not this world's. Orogen does not consume them, so
  they were never overridden. **`config/planet.yaml` is authoritative** for
  everything except radius and gravity.
- `manifest.planetRadiusKm` at the top level says 6371; `manifest.planet.radiusKm`
  says the real 7645.2. Trust the latter.

### `source/maps/` is a different build — do not treat it as canonical

The PNGs are named for planet code `09oa7lj8kek17v5639gvhe`, which decodes to the
same seed and sliders but **510,000 mesh regions**. `source/worldorogen_seed.txt`
holds `01eshm059lt0b9mpgro2y83t`, the code for the current canonical build:
**2,500,000 regions**, basins preserved, lithology on. Region count is a
generation parameter, so the maps are a coarser, genuinely different terrain —
close in aggregate (land fraction 0.413 vs 0.414) but not the same coastlines.

Two further traps if regenerating them:

- `export-maps.mjs --code` **ignores** the code's `basinSlider` and
  `lithologyStrength` and defaults both on; `export-planet.mjs --code` **honours**
  them. The two tools disagree.
- Neither takes planet radius/gravity from the code. The maps were rendered at
  Earth defaults, so their relief scaling differs from the exports.

To render maps that actually match `source/exoplasim-*/`:

```bash
cd /home/cfutro/git/planet_heightmap_generation
node --max-old-space-size=12288 tools/export-maps.mjs \
    --code 01eshm059lt0b9mpgro2y83t --regions 2500000 \
    --radius 7645.2 --gravity 10.1989 --width 16384
```

## Planet parameters

`config/planet.yaml` is the single source of truth. Radius 1.2 R⊕, mass 1.5 M⊕,
30 h rotation, 32° obliquity, e = 0.02, around a K2.5V star (0.80 M☉, 0.3236 L☉,
4965 K). 1 bar atmosphere at 450 ppm CO₂.

Known inconsistency, unresolved: the ExoPlaSim scripts derive
`g = 1.50/1.20² × 9.80665 = 10.2153 m/s²`, while the Orogen export was generated
with `g = 10.1989 m/s²` (1.04 g⊕). A 0.16% difference — it affects the 1/g relief
scaling in the terrain. Not worth regenerating for, but do not "fix" one to match
the other without deciding which is right.

The `model:` block in `planet.yaml` is ExoPlaSim-specific (resolution, layers,
timestep, output cadence). The rest is world-level. Do not split the file
casually: `config_sha256` in every `exoplasim/runs/*/run_manifest.json` pins its
exact contents, and `continue_exoplasim.py` refuses to resume a run whose config
hash has changed.

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
