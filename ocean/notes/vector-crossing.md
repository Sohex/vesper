# The vector crossing onto the ocean grid, and the check that can fail

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean cells, surface velocity, sea ice and
coastlines name modelled quantities on a simulated grid, not observations of
anything.

Measured 2026-09-07. The driver is `ocean/scripts/build_ocean_grid.py`, the
operator is `lib/remap.py:Crossing.apply_vector` and the geometry is
`lib/gridding.py:GridSpec.cell_centroids` and `cell_frames`. Every number below
is in `ocean/data/<build>/ocean_grid.json` with its own tolerance beside it.

**The headline.** A tangent vector does not remap like two scalars, and the
obvious acceptance test does not say so. A rigid rotation's cartesian integral
over the sphere is analytically zero, the vector crossing returns it -- and the
componentwise remap returns it too, to round-off. The check only discriminates
when it is taken over a region that is not zonally symmetric, and there the two
separate by five orders. A second inherited claim, that the frame error vanishes
at the equator and grows with latitude, is refuted by the same measurement: what
sets the error is the destination cell's longitude span, which does not depend
on the row.

---

# 1. What the operator does

The ocean returns a surface velocity to the next ExoPlaSim baseline
(`notes/audits/ocean-and-marine-biosphere.md` section 7b; world-pt8 is the
consumer). It arrives on GOLDSTEIN's grid and has to reach the atmosphere's, at
every rung of the ladder.

"East" at one longitude is not "east" at another. Remapping the east and north
components as two independent scalars therefore averages numbers that are
components in different frames, and the result is a vector that no single frame
holds. `apply_vector` instead lifts the pair into the sphere's own three
cartesian components at the SOURCE cells' frames, remaps each of the three with
the SAME weights under the SAME semantics as an ordinary scalar -- the lift is
linear, so whatever `apply` conserves it conserves componentwise -- and projects
onto the DESTINATION cells' frames.

The frame is built at each cell's own AREA CENTROID, the mean of the position
vector over the cell, which lies inside the sphere and is not the cell centre
lifted onto it. Both integrals are closed form on a cell bounded by two sines of
latitude and two longitudes, so nothing in the geometry is quadrature. The
centroid is the reference point a cell has that is defined by the cell rather
than by a convention, and it is what makes the acceptance field exact: a rigid
rotation is linear in position, so its mean over a cell is the rotation applied
to the cell's mean position.

What the projection cannot carry is reported rather than removed. `radial` is
the part of the mapped vector that no longer lies in the destination cell's
tangent plane. It is not renormalised, because scaling the horizontal part back
up to the original magnitude would destroy the conservation the lift just
bought, and it is returned so a caller can rebuild the mapped cartesian vector
exactly.

---

# 2. The acceptance field, and why its integral is analytic

The field is a rigid rotation `omega x r`. Its integral over the sphere is
`omega` crossed into the integral of the position vector, and that integral is
the origin, so the answer is zero for any `omega` -- an analytic value rather
than a plausible one. Two properties make it usable as written:

- Its exact cell mean is `omega x c` at the cell's own area centroid `c`, so the
  source integral is zero BEFORE any remapping, to round-off, rather than to the
  order of a cell.
- `omega x c` is perpendicular to `c` and therefore lies in the tangent plane of
  the frame built on `c`. The east/north pair carries all of it, so the source
  radial component is zero and not merely small. That is asserted rather than
  assumed: 2.1e-16 against a bar of 1e-12.

The conserving remap preserves the source integral exactly, cell for cell, so
the destination integral has the same analytic answer.

---

# 3. The full-sphere integral does not discriminate, and the control is what says so

Over the whole sphere, at T42 onto a 36 x 36 equal-area GOLDSTEIN grid:

| what | residual, relative to the field's own area-weighted magnitude |
| --- | --- |
| the source quadrature | 1.7e-16 |
| the vector crossing | 1.2e-15 |
| the componentwise control | 8.7e-17 |

The control passes. The bar is 1e-12 and it misses nothing.

The reason is a symmetry rather than an accident. Both grids are uniform in
longitude, so the overlap weight between a source column and a destination
column depends only on the difference between them, while the local frames vary
sinusoidally with longitude. The frame error is then a pure phase and it cancels
around every row. A rotation about the polar axis behaves the same way, which
refutes the framing this operator arrived with: the recovered code's docstring
named the polar axis as the degenerate case and the equatorial axis as the
discriminating one, and the axis is not what separates them.

So the check has to be taken over a region that is not zonally symmetric. The
conservation law is unchanged -- `restrict` closes the columns, so all of a valid
source cell's content lands somewhere valid -- and the source's own integral over
the valid set is still closed form, `omega` crossed into a sum of analytic cell
centroids. The window used is a quarter turn placed on the rotation axis's own
longitude, so it is chosen against the field and not against either grid's
columns.

| pair | vector crossing | componentwise control |
| --- | --- | --- |
| T21 onto 36 x 36 | 5.6e-16 | 2.3e-04 |
| T42 onto 36 x 36 | 8.3e-16 | 8.8e-05 |
| T85 onto 72 x 72 | 1.1e-15 | 3.7e-05 |

The acceptance bar is 1e-12 and the control floor is 1e-09, both fixed from
`lib/remap.py`'s own round-off floors before any of this was run. The crossing
clears its bar by three orders and the control misses it by four to five.

The real ocean mask is not zonally symmetric either, which is why the wedge is
the case that matters rather than the tidy one.

---

# 4. What the lift is worth, in the units world-pt8 cares about

The per-cell gap between the two answers, over destination cells above a fifth
of the peak mapped speed, as a share of the LOCAL mapped speed:

| pair | median | p99 | max | radial part discarded, max |
| --- | --- | --- | --- | --- |
| T21 onto 36 x 36 | 0.50% | 2.07% | 2.09% | 9.6% |
| T42 onto 36 x 36 | 0.17% | 0.40% | 0.43% | 1.8% |
| T85 onto 72 x 72 | 0.070% | 0.17% | 0.20% | 0.89% |

Two things to read off it.

**It scales with the cell, not with the latitude.** Across the ocean rows of the
T42 pair the gap runs 0.36% at 76 degrees, 0.22% at 28, 0.27% at the equator and
0.36% at the other pole. The claim the operator arrived with -- that the error
vanishes at the equator and grows toward the pole -- is not what the geometry
does. The frame rotates by the destination cell's LONGITUDE SPAN, and on both of
these grids that span is the same in every row. The equatorial-versus-polar
intuition is about a fixed physical distance, and a fixed longitude span is not
one.

**The radial part is the larger number.** What the destination's tangent plane
will not hold is up to 9.6% of the local speed at the coarsest pair and 0.89% at
the finest. That is a statement about how far apart two grids' tangent planes are
over a cell, not an error, and it is why it travels in the ledger instead of
being normalised away.

---

# 5. What this does NOT establish

It is a geometry result and no ocean has run. Nothing here says what a modelled
surface velocity looks like, how much of it survives to an ice-advection term, or
whether a 0.4% velocity error matters to a sea-ice field -- world-pt8 owns that
and needs a velocity to answer it.

The GOLDSTEIN side is a candidate comparison support, as
`analysis/ocean_remap.py` already labels it. OCN-11 supplies the accepted ocean
support, and until it does no line in this work declares a cell wet.

The wedge is a geometric window and not an ocean basin. It is enough to break the
zonal symmetry, which is the only property the check needs from it; it is not a
claim about what the mask will cost.
