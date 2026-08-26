# Glaciers: a rough pass

*Measured 2026-08-16 on `precarve-zoned-g1281`, against the 0.968 T42
climatology `baseline_t42_s0968`.*

**This is deliberately rough.** Its only job is to confirm or refute the claim
that this world's glaciers are controlled by relief rather than latitude. The
detailed treatment belongs with the downscaling and sub-grid sampling machinery
that does not exist yet, and should not be improvised here. Nothing below is a
property of the world yet; it is a bound with its assumptions attached.

## Why the GCM cannot answer this

ExoPlaSim reports `glac = 0` everywhere, and that is not a result. A T42 cell is
roughly 300 km across and its orography is the cell-mean, so a peak inside the
cell is invisible to the model. The question is entirely sub-grid.

Measured from the 2.5M-region mesh, binned to T42 with `lib/gridding.region_cells`
and integrated by land area:

| \|lat\| | cells | cell mean | q95 | q99 | max | q99 - mean |
| --- | --- | --- | --- | --- | --- | --- |
| 0-30 | 1644 | 0.42 km | 1.15 | 1.49 | 1.71 | 1.07 |
| 30-50 | 999 | 0.28 | 0.82 | 1.10 | 1.27 | 0.82 |
| 50-60 | 769 | 0.25 | 0.73 | 0.97 | 1.08 | 0.72 |
| 60-90 | 2143 | 0.21 | 0.59 | 0.72 | 0.78 | 0.51 |

Global land mean 0.407 km; the highest single cell reaches 4.593 km, which is
4.447 km above its own cell mean. So the model is unaware of up to 4.4 km of
relief in the cells that matter most.

**A consistency check fell out of this.** 5.769 km, the peak in `source/maps/`,
times the gravity ratio 10.1989/12.81 is 4.593 km exactly. `reliefScale = 1/g` is
being applied consistently and this build is the same landscape at 12.81 m/s2.

## The join was verified before it was used

Mesh-derived land fraction against the `N064_surf_0172.sra` land mask ExoPlaSim
actually ran with: **100.00% binary agreement** at the 0.5 threshold, by index and
with no roll. Rolling 180 degrees drops cell agreement to 58%. This is checked
rather than assumed because keying on lon/lat across that boundary has silently
matched zero cells three times; see `docs/src/practice/failure-modes.md`.

## Method

For each land cell, the elevation at which the warmest month reaches freezing:

    z* = cell_mean_elevation + (T_warmest_month - 273.15) / lapse

then integrate mesh land area above `z*` within that cell. Warmest-month mean
below freezing is the standard proxy for a surface that never fully melts.

## Result

Lapse 6.5 K/km. The 0.945 rows apply a uniform -4.7 K offset to the 0.968
climatology, since the 0.945 baseline run had not finished when this was written.

| phase | dT vs 0.968 | glacier | % of land | area-weighted \|lat\| |
| --- | --- | --- | --- | --- |
| 0.968, measured | 0.0 K | 0.50 Mkm2 | 0.158% | 56.9 |
| 0.945 mean, proxy | 4.7 K | 1.85 | 0.583% | 51.2 |
| 0.945, long-cycle minimum | 8.2 K | 4.17 | 1.313% | 48.5 |
| 0.945, aligned minimum | 10.3 K | 6.90 | 2.173% | 47.1 |

At 5.5 K/km the areas roughly halve (0.23 Mkm2 at 0.968, 1.08 at the 0.945
proxy) and the picture is unchanged.

**The claim is confirmed.** Glaciers are where the relief is, not where the
latitude is. The area-weighted mean glacier latitude is 47 to 57 degrees across
every case, glaciated cells span 1 to 85 degrees, and the mean latitude *falls*
as the world cools because progressively lower-latitude peaks come in. There is
no polar ice cap: at 32 degrees obliquity the poles get too much summer
insolation, and the coldest land cell's warmest month at 0.968 is +8.9 C, at
latitude 57.

Scale: Earth's non-ice-sheet mountain glaciers are about 0.5 Mkm2, roughly 0.35%
of land. So the 0.945 mean is a few times Earth's mountain glaciation with no ice
sheets at all -- "some permanent glacier, not a lot", which is the brief.

**The cycle moves this by a factor of nearly four**, 1.85 to 6.90 Mkm2 between
cycle mean and aligned minimum. That is the mechanism behind wanting
non-commensurate periods: successive grand minima differ in depth, so advances
reach different distances and the moraine record is legible rather than
repeating. See `exoplasim/notes/parameter-decisions.md`.

## What is wrong with this and must be fixed downstream

- **The -4.7 K offset is a proxy**, not a measurement. The redo against the
  real 0.945 climatology is CLIM-34.
- **The lapse rate is no longer assumed.** It is measured from the model's own
  10-level profile by `lib/lapse.py`, and the measurement says the bracket this
  note explored sat entirely on the shallow side; the correction section below
  has the numbers and the direction. The TABLE above still shows what 6.5 K/km
  produced, as the dated record of what was computed.
- **No mass balance.** Warmest-month-below-freezing is a temperature criterion.
  A real glacier also needs accumulation, and a cold dry peak does not glaciate.
  This is the single largest reason to treat these numbers as an upper bound.
- **No lapse-rate feedback on precipitation phase, no ice albedo, no flow.**
- **Sub-grid hypsometry was computed inline here and thrown away.** It is now an
  artifact: `hydrography/scripts/build_spatial_support.py` writes the whole
  distribution of mesh elevations inside every cell of every requested ladder
  rung, over the LAND population rather than over land and seabed together, and
  `lib/gridding.py:area_fraction_above` reads a share above a threshold back out
  of it. Nothing needs to integrate it inline again. SPAT-3 closing GRID-2.

## Correction, 2026-08-19: the rate is measured, and the old bracket was one-sided

Measured on the baseline climatology of `run_8c2e1ff9ab5e` by `lib/lapse.py`:
temperature against hypsometric height over sigma 0.45 to 0.90 (about 1.1 to
3.5 km above the surface, the band `z*` lives in), per land column, land-area
weighted, output bins weighted by their record counts.

| rate | K/km |
| --- | --- |
| annual mean | 6.78 |
| warmest bin per cell, the one this note's criterion extrapolates | 7.82 |
| dry adiabat g/cp, from the configured composition | 12.75 |

Two things follow, one about the audit's expectation and one about this note.

The audit that opened PHYS-12 expected "nearer 8.5" by scaling Earth's 6.5 with
the gravity ratio. The measurement lands below that: this atmosphere is more
stably stratified relative to its own adiabats than a pure `g/cp` scaling
assumes, which is exactly why the rate had to be measured rather than derived.
The warm-season rate, the operative one here, is 7.8.

For this note's table: `z* = cell_mean + (T_warmest - 273.15) / lapse`, so at
7.8 K/km the height term shrinks to 0.83 of its 6.5 value, the freezing surface
drops, and every area in the table grows. By the note's own measured
sensitivity -- 5.5 K/km roughly halves the areas relative to 6.5 -- 7.8 roughly
DOUBLES them, putting the 0.945 cycle-mean nearer 3.5 than 1.85 Mkm2. The exact
factor is not computed here, because this note already owes a redo against the
real 0.945 climatology and the two corrections belong in one pass.

The direction is the finding. Every other caveat in the section above pushes
the same way, making the areas an upper bound; the lapse assumption pushed the
OPPOSITE way, an undercount, and the bracket 5.5 to 6.5 could never have said
so because both ends sat below the measurement. The note's areas are therefore
no longer a clean upper bound: the missing mass balance still argues they are
too large, the corrected lapse argues they are too small, and only the redo
settles which wins.

## Extension, 2026-08-24: the resolution ladder never closes this gap

The note says the question is entirely sub-grid and demonstrates it at T42. What
it does not say is whether that survives the resolution ladder, and it does --
which turns a property of the current configuration into a property of the
problem.

At 1.20 Earth radii the circumference is 48,090 km, so the ladder gives:

| rung | longitudes | cell width |
| --- | ---: | ---: |
| T21 | 64 | 751 km |
| T42 | 128 | 376 km |
| T85 | 256 | 188 km |
| T127 | 384 | 125 km |
| T170 | 512 | 94 km |

Against a 7.60 km Orogen mesh, the top of the ladder is still twelve times
coarser than the terrain the glaciers sit on. The sub-grid peak excess that
drives the whole result -- 1.065 km on land average, 2.969 km at the 90th
percentile -- stays sub-grid at every rung this project will ever run.

**So `glac = 0` is a permanent property of the model's orography rather than a
statement about the world, and no rung of the ladder converts it into one.** The
glacier answer is a diagnostic computed on the mesh, and the model's own glacier
module is carried for its orographic feedback on the cells it can resolve rather
than as the instrument that answers where the ice is. `analysis/ice_mask_freezing_height.py`
is that diagnostic and PHYS-13 owns it.

This is the same class of argument that put the ocean circulation offline: the
quantity wanted lives at a scale the atmosphere will not be run at, so it is
computed beside the atmosphere rather than inside it.

**It also settles the flow question by removing its premise.** Ice flow would be
worth porting if the model were about to resolve ice sheets; it is not, and the
areas in this note are a few times Earth's mountain glaciation with no ice sheets
at all. The port cost is separately prohibitive -- Earth geography welded into
the source, dimensional limiters sitting on g^3 and g^4 quantities that bind
about 2.2x more often at this gravity while truncating the effect being studied,
and no namelist path to any material constant. `notes/external-model-survey.md`
section 45d.

The ordering that follows is PHYS-13, then CLIM-53's accelerator verdict, and
flow is not on the board. Note for CLIM-53 that `newsnow.x` and `buildice.x`
hardcode `NLAT = 32` and `NLON = 64`, so they are T21-only and would silently
misread any higher rung -- which matters more against a ladder than against one
resolution.

## Extension, 2026-08-25: what replacing the placement would change

*Measured on `precarve-craton-10m`'s T42 export with the T21 bootstrap
climatology `bootstrap_regular_climatology.nc`, warm-season environmental lapse
rate 6.520 K/km read from `lib/lapse.py` at the time of measurement.*

`vendor/orogen/js/glacial-ice.js` now takes the ice placement from outside the
generator, and `analysis/ice_mask_freezing_height.py --write-mask` produces it.
This is the size of the difference between the two placements, evaluated on the
same ground.

**Three reasons every number here is a LIMIT rather than a state.** The build is
pre-carve, so the basins are still closed and the land area they occupy is the
uncarved one. The climatology is a BOOTSTRAP, whose numbers are not the
baseline. And the generator applies the latitude ramp to the pre-erosion
surface, while this evaluates it on the exported one, so it locates where the
two CRITERIA disagree rather than replaying either build.

| land area glaciated | latitude ramp at `glacialErosion` 0.8 | freezing-height mask |
| --- | ---: | ---: |
| any ice at all | 33.727% | 0.967% |
| index above 0.1 | 15.983% | 0.967% |
| index above 0.4 | 7.125% | 0.967% |
| index above 0.6 | 3.476% | 0.967% |

Of the 33.727% the ramp glaciates, 32.789 percentage points carry no ice under
the thermal criterion and 0.030 points are ice the ramp misses. The
area-weighted mean index over land -- the quantity `iceFlow` accumulates
downstream -- falls from 0.0678 to 0.0077, a factor of 8.8. Carving goes as
`iceFlow^0.6`, so where ice survives at all the deepening falls by roughly
`8.8^0.6`, about 3.6x, over between a sixteenth and a thirty-fifth of the area.

**The ramp's error is altitude, not latitude.** The area-weighted mean
`|latitude|` of glaciated land is 54.8 degrees under the ramp and 56.2 under the
mask, which is the same ice line to within the width of the transition. The mean
`elevation_km` of glaciated land is 0.84 km under the ramp and 2.28 km under the
mask. The ramp is putting ice on low ground at high latitude that is nowhere
near freezing in the warmest bin, and that low ground is most of the area.

That is consistent with the gate itself being right and its calibration being
Earth's: the dimensionless altitude gate is the gravity-invariant form, since
the relief ceiling and a dry-adiabatic freezing height both go as 1/g and cancel
(`notes/audits/orogen-gravity.md`). What the ramp cannot know is the surface
temperature the freezing height is measured down from, and on this planet that
is the whole of the difference above.

**No terrain has been regenerated.** Generating is loop A and the numbers above
are what a regeneration would face, not what one produced.
