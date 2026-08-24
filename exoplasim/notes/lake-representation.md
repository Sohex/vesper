# Representing lakes in ExoPlaSim

Assessment of what ExoPlaSim 3.4.2 can and cannot do with the endorheic basins
preserved in `source/`. Read against the vendored Fortran in
`.venv/lib/python3.12/site-packages/exoplasim/plasim/src/`, not from
documentation.

## There is no lake model

`grep -in lake *.f90` over the whole PlaSim source returns nothing. Surfaces are
land, sea, sea ice, or glacier. Whatever we do with lakes is built out of those
four.

The land/sea mask is also not fractional, despite `yls` being declared
"land sea mask (frac.)". `oceanmod.f90:256` hard-binarises it:

```fortran
where (yls(:) > 0.5)
   yls(:) = 1.0
elsewhere
   yls(:) = 0.0
```

So a grid cell is entirely land or entirely slab ocean. There is no tiling and
no partial-water cell.

## Most lakes are far below the grid

Basin footprint here is the area at spill level, so it is the maximum a lake can
reach, not its likely size.

| Grid | Cell area | Basins >= 1 cell | Basins >= 4 cells | Share of total footprint in >= 1-cell basins |
| --- | ---: | ---: | ---: | ---: |
| T42 128x64 | 89,660 km2 | 118 / 3629 | 29 | 49.6% |
| T63 192x96 | 39,849 km2 | 289 / 3629 | 62 | 61.2% |
| T85 256x128 | 22,415 km2 | 572 / 3629 | 118 | 71.0% |

The median basin is 7,539 km2, about one twelfth of a T42 cell. At T42 roughly
97% of basins by count cannot be represented as mask cells at all, and the ones
that can hold 5.72% of the planet's area at full spill. Actual lakes will be
smaller than that.

The consequence worth internalising: lakes are primarily a map product, computed
on the 2.5M-region mesh where they are properly resolved. The GCM question is
narrower, namely whether the largest hundred or so perturb the climate enough to
matter.

## Four levers, in increasing cost

### 1. Surface-property perturbation, no patch required

`surfmod.f90:206-255` registers the per-gridpoint fields that can be supplied as
SRA input, written into the run directory as `N{NLAT:03d}_surf_{code:04d}.sra`
(`surfmod.f90:185`). The ones that matter here:

| Code | Array | Meaning |
| --- | --- | --- |
| 174 | `dalbcl` | background albedo |
| 175 / 176 | `dalbcl1` / `dalbcl2` | albedo below / above 0.75 um, which is what this project's two-band setup uses |
| 173 | `dz0clim` | roughness length |
| 229 | `dwmax` | maximum soil field capacity, per gridpoint |
| 232 | `dglac` | glacier mask |
| 172 | `dls` | land/sea mask |
| 129 | `doro` | orography |

This works at any lake size, because the values can be area-weighted by the lake
fraction within each cell. It is the only lever that reaches the 97% of basins
that are subgrid.

Note that ExoPlaSim's Python wrapper only ever copies two of these:
`configure()` writes `landmap` to 0172 and `topomap` to 0129 (`__init__.py:2939`)
and exposes nothing for the rest. The model reads them regardless, so supplying
them means writing the SRA files into the run directory ourselves. The run
scripts already know how to do this for 0129 and 0172.

### 2. `dwmax` is already a crude endorheic basin

This is the useful discovery. PlaSim's land surface is a bucket:

- `landmod.f90:1097-1099`: `dwatc` accumulates, `drunoff = max(0, dwatc - dwmax)`,
  and `dwatc` is then clamped to `dwmax`.
- `landmod.f90:103-104`: `drhsland = 0.25`, `drhsfull = 0.4`. The evaporation
  wetness factor `drhs` reaches 1 once soil water exceeds 40% of `dwmax`.
- `seamod.f90:145`: over sea, `drhs = 1` outright.

So a cell whose bucket is large and full evaporates at the same unlimited rate as
open water, stores water, and spills the excess to runoff. Storage, potential-rate
evaporation, and overflow are exactly the three terms an endorheic balance needs.
Setting `dwmax` per cell from basin hypsometry gets us a first-order lake without
touching the model at all.

What it does not get us: the thermal inertia of a water column (land heat
capacity comes from the scalar `soildepth` and `cpsoil`), the albedo change,
which must come from 174/175/176 separately, and ice.

### 3. Mask flip to slab ocean, no patch required

Set `dls = 0` on cells that are mostly lake. They become slab ocean with correct
albedo, unlimited evaporation, real thermal inertia, and sea ice via `icemod`.

Two limits. The binarisation at 0.5 means this is only honest for lakes covering
most of a cell, so at T42 it applies to a few dozen basins. And the mixed layer
depth is global.

### 4. Per-cell mixed layer depth, small patch

`mldepth` is a scalar namelist parameter (`oceanmod.f90:48`, wrapper at
`__init__.py:2899`), but the array it feeds is not:

```fortran
real :: ymld(NHOR,NLEV_OCE) = 0.   ! oceanmod.f90:79
...
ymld(:,jlev) = dlayer(jlev)        ! oceanmod.f90:232
```

`ymld` is already dimensioned per gridpoint and every use site indexes it that
way, including all the slab thermodynamics (`oceanmod.f90:811, 859, 930, 1057,
1101`). Only the initialisation is uniform. Making the mixed layer depth spatially
varying therefore means changing one assignment and registering a `surfcode` for
the input, with no changes to the physics.

Two properties make this cheap. `ymld` is absent from both the restart write and
the restart read (`oceanmod.f90:282-292` and `591-601`), so it is re-initialised
on every start, which means the patch applies cleanly to every segment and needs
no restart migration. And it is the same shape of change as
`patches/exoplasim-3.4.2-star-cycle.patch`, which is already maintained here, so
the build and patch-reversal tooling in `scripts/build_star_cycle_exoplasim.sh`
covers it.

This is worth doing only if a resolvable lake turns out to be shallow enough that
a 50 m slab misrepresents its seasonal cycle. Decide it after the first pass, not
before.

## What stays out of reach

- Subgrid lakes as anything other than modified land properties.
- Lake ice, unless the lake is a mask-flipped ocean cell.
- Lake temperature as a distinct prognostic from soil or slab temperature.
- Any lake-atmosphere feedback finer than one grid cell, which is most lake-breeze
  and lake-effect behaviour.

## Recommendation

Lakes are a map product first and a boundary condition second. Compute levels on
the mesh, render them at full resolution, and feed the GCM only the part it can
represent.

For the first climate pass, do not flip any mask. Carry lakes as surface-property
perturbations: area-weighted albedo in 174/175/176, roughness in 173, and
`dwmax` in 229 set from basin hypsometry. That covers every basin regardless of
size, needs no patch, and keeps the run comparable to a no-lake control by
changing only three input fields.

Mask flipping and per-cell `ymld` were held pending the first pass; CLIM-26
measured the lake moisture term at -0.9% of land precipitation, so neither was
then planned. That annual moisture response does **not** bound the seasonal
temperature, ice and evaporation-phase error from giving a lake the land soil
thermal column. HYD-21 now owns an offline/reduced seasonal bound using basin
area and hypsometric depth before any coupled implementation or run is proposed.

---

## Revisions after wiring the levers up

The recommendation above named three fields: albedo in 174/175/176, roughness in
173, and `dwmax` in 229. Albedo is done and behaves as described. The other two
need correcting, both in the direction of claiming less.

### `dwmax`: the wetness argument has the sign backwards

Lever 2 says a large full bucket evaporates at open water's rate, so set `dwmax`
from basin hypsometry. The first half is right and the prescription that follows
from it is not.

`drhs` reaches 1 once soil water exceeds **40% of `dwmax`** (`landmod.f90:103-104`).
That threshold scales with the bucket, so a *deeper* bucket needs proportionally
more water to reach the same wetness. Large is right for storage and wrong for
wetness, and on this planet the two do not point the same way: the lakes sit in
the arid cells, which are precisely where the water to fill a large bucket does
not exist. `build_surface_soil_water.py --lakes` therefore sets a **shallower**
bucket on the lake fraction, defaulting to 0.2 m against the 0.5 m land default.

### `dwmax` cannot sustain a lake at all, and this is the harder limit

A lake here is fed by its catchment. That water never reaches the evaporating
bucket:

- `dwatc` gains only from local precipitation minus evaporation
  (`landmod.f90:1097`).
- Routed river water accumulates into `driver` (`landmod.f90:1527`), a separate
  store that is advected downhill by `mkradv` and discharged at the coast.
- Nothing returns `driver` to `dwatc`.

So a cell's annual evaporation is capped by its own precipitation plus storage,
however `dwmax` is set. Lake evaporation stays underestimated and no setting of
this field fixes it. What the lever buys is the **seasonal partition**: a lake
cell that holds potential-rate evaporation while it has water, instead of going
moisture-limited on the first dry day. That is worth having and is not what the
section above promised.

This matters less than it sounds for the carve verdict, which computes lake
evaporation externally with Penman rather than reading it from the model. It
matters for the climate, as a latent-versus-sensible partition error over the
3.09% of the planet the lakes cover.

### Roughness: do not supply 173 for lakes

Assessed and rejected, rather than left undone. `dz0land` is 2.0 m and open water
is about 1.5e-4 m. At this planet's lowest-level reference height, which
`lib/lapse.py:reference_height_m` derives hypsometrically at the configured
gravity, supplying water roughness would cut turbulent exchange on those cells
by about **11x** across the span in which surface water is liquid.

Real lakes tolerate that because evaporation and a water column's heat capacity
hold the surface cool. We have neither. A cell given water's roughness, land's
heat capacity and a moisture supply capped by local rainfall would decouple from
the atmosphere and run hot, which is a larger and less physical error than the
one it replaces. Revisit only alongside a mask flip, where the water column
arrives with it.

### Separately: `dz0land = 2.0` m

Superseded by `build_surface_roughness.py`, which supplies code 173 from land
cover and subgrid relief; the excess-exchange argument that motivated it
is in that script's header and `docs/src/reference/config-rationale.md`.
