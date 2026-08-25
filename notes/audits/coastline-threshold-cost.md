# What the binary coastline threshold costs, in both signs, on every rung

**Measured:** 2026-08-25, on `precarve-craton-10m`, terrain hash `ab0d679b`,
10,000,005 regions, at every rung of the T21/T42/T85/T127/T170 ladder.
`exoplasim/scripts/build_boundary_conditions.py:coastline_ledger` is the
measurement and it runs as part of every boundary-condition build, so the
numbers below are reproduced in `boundary_conditions_report.json` rather than
kept only here.

This is worldbuilding. Vesper is an invented super-Earth; every quantity below
is a property of that planet's modelled surface or of the gridding this
pipeline applies to it. Nothing was simulated: no World Orogen generation, no
ExoPlaSim run.

SPAT-5's first clause asks for the binary coastline threshold's area, store and
flux cost on every ladder rung. This is the area and the store; the flux
follows from the area once a climatology exists and is not derivable here.

## The threshold, and why it is not a choice about the model

`config/planet.yaml:geography_land_threshold` is 0.5 and
`docs/src/reference/config-rationale.md` records it as DETERMINED, with its
cost stated as a coastline effect and not sized. The model's mask is binary
whatever this builder writes: `oceanmod.f90:256` hard-binarises `yls` at 0.5,
so there is no tiling and no partial-water cell to write into. The threshold is
therefore a rounding, and the only open questions are how large it is and
whether 0.5 is where it should sit.

## The net hides the rounding by a factor of seven

The report already carried the net, as `land_fraction_gauss_weighted` against
`mesh_land_fraction`. The rounding is two one-signed errors that partly cancel:
land in a sub-threshold cell is dropped into the slab ocean, and ocean in a
supra-threshold cell is promoted to dry land.

| rung | land dropped | water promoted | NET | rounding / net |
| --- | ---: | ---: | ---: | ---: |
| T21 | 11.58% | 10.56% | -1.14% | 10.1 |
| T42 | 7.74% | 6.82% | -0.99% | 7.8 |
| T85 | 5.38% | 4.43% | -1.00% | 5.4 |
| T127 | 4.37% | 3.44% | -0.96% | 4.5 |
| T170 | 3.67% | 2.93% | -0.76% | 4.8 |

Dropped as a share of the mesh's own land area; promoted as a share of the land
area the model receives.

**The net is nearly flat across the ladder and the rounding is not.** The net
moves by a factor 1.5 from T21 to T170 while the land actually misassigned
falls by 3.2. A convergence claim read off the net would report the coastline
as insensitive to support when what is happening is that two errors three times
larger are cancelling. That is the reading the two halves exist to prevent.

## The terrain the fork exists to preserve is what the threshold floods

`source/README.md`'s first rule rejects `land_mask` because it puts dry closed-
basin floor below sea level under water. The threshold does the same thing to
part of the same terrain, through a different door: a cell whose land share
falls below 0.5 goes to the slab ocean whatever its elevation.

| rung | share of below-datum land dropped to ocean | land volume closure |
| --- | ---: | ---: |
| T21 | 6.44% | +2.23% |
| T42 | 3.18% | +2.49% |
| T85 | 1.32% | +1.20% |
| T127 | 0.82% | +0.43% |
| T170 | 0.47% | +0.31% |

Land volume is the extensive closure: area times elevation about the datum,
summed over the mesh's land against the same sum over what the model receives.
It carries the sign the area does not, because the preserved terrain sits below
the datum.

This is the most rung-sensitive quantity in the ledger, falling by a factor of
14 across the ladder while the area rounding falls by 3.2. It is also the one
with a rule behind it rather than a preference, which is why it is tracked
separately from the coastline area.

## How much of the grid is partial at all

| rung | cells with a partial land share | share of planet area they hold | cells the mesh does not cover |
| --- | ---: | ---: | ---: |
| T21 | 1,053 | 55.5% | 0 |
| T42 | 2,715 | 36.0% | 0 |
| T85 | 6,833 | 22.6% | 0 |
| T127 | 11,578 | 17.1% | 32 |
| T170 | 16,610 | 13.8% | 365 |

This is the size of the problem a tile path would address, and it is the number
to weigh a tile implementation against: at T42 a third of the planet's area
sits in cells that are neither wholly land nor wholly ocean.

The last column is the opposite constraint and it grows the other way. Beyond
T85 the 10M mesh stops covering every Gaussian cell and `gridding.land_weighted`
falls back to the nearest region centre, which is a point sample of a
categorical field -- the exact reduction this builder exists to avoid. 365 cells
at T170 is small, but it means a finer rung is not a free route out of the
coastline rounding: past T85 the reduction starts substituting for the mesh.

## The third class is a latent hazard and reads zero today

`surface_class` partitions the mesh into ocean, subaerial land and inland
water, and the threshold is applied to the LAND share alone. Inland water
therefore competes with land for its own cell's classification, so a cell
holding no ocean at all goes to the slab ocean once its lakes take more than
half of it.

On this build the hazard is inert: `INLAND_WATER` holds 0 regions and 0.0% of
area, so the promoted inland-water area and the count of wholly terrestrial
cells sent to the ocean are both exactly zero. Both are in the ledger so that
the first build carrying inland water reports it rather than absorbing it into
the ocean share.

## What follows

- The two halves are now in every boundary-condition report. The net stays
  beside them; it is the right number for a mass budget and the wrong one for a
  convergence claim.
- Whether 0.5 is the right threshold is a separate question from whether the
  boundary should be fractional at all, because a threshold that conserved land
  area would still not conserve which land. The below-datum row is what makes
  it a question rather than a preference.
- The flux half of SPAT-5's first clause is not derivable at this step. A cell
  flipped to slab ocean evaporates at the open-water rate and one flipped to
  land at the bucket rate, and the difference between those rates is a
  climatology.
