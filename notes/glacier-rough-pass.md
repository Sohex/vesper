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
matched zero cells three times; see `notes/failure-modes.md`.

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

- **The -4.7 K offset is a proxy**, not a measurement. Redo against the real
  0.945 climatology.
- **The lapse rate is assumed**, 5.5 to 6.5 K/km, not taken from the model's own
  temperature profile, which is available on 10 levels.
- **No mass balance.** Warmest-month-below-freezing is a temperature criterion.
  A real glacier also needs accumulation, and a cold dry peak does not glaciate.
  This is the single largest reason to treat these numbers as an upper bound.
- **No lapse-rate feedback on precipitation phase, no ice albedo, no flow.**
- **Sub-grid hypsometry is not persisted.** It was computed inline. The reusable
  version belongs with the downscaling machinery rather than as a one-off here.
