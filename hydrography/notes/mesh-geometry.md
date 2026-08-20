# The region mesh has a Voronoi dual, and `cell_area` is not its area

This is worldbuilding. Vesper is a fictional super-Earth and this document is
about the geometry of the mesh its terrain is exported on.

Measured 2026-08-20 on `precarve-craton`, while building the water table solver,
which is the first thing in this project to need face widths rather than only
cell areas. `hydrography/scripts/groundwater.py:Geometry` is the code.

## The mesh is Voronoi, and the dual recovers what the export does not carry

The export gives a CSR neighbour list and a `cell_area` per region. A
finite-volume scheme needs two things neither of those supplies: the width of
the face shared by two neighbours, and the distance between their generators.

Both follow from the geometry, because the regions are the cells of a spherical
Voronoi tessellation. The convex hull of the region centroids is therefore the
Delaunay triangulation, each Voronoi vertex is a triangle's circumcentre, and
the face between two neighbours is the arc joining the circumcentres of the two
triangles sharing that Delaunay edge.

The reconstruction is checkable and checks out:

| what | value |
| --- | --- |
| Delaunay edges | 7,499,997 |
| undirected pairs in the export's own adjacency | 7,499,997 |
| edges present in one and not the other | 47, each way |
| exact spherical Voronoi area, as a fraction of `4 pi R^2` | 1.00000000 |

The 47 disagreeing edges are cospherical quadrilaterals, where which diagonal
the Delaunay takes is genuinely ambiguous. The dual's own edge set is the one to
use, because it is the one that has face widths attached.

Two numbers are worth keeping for whoever next writes an operator on this mesh:
adjacent generators come as close as ten metres against a median separation of
about twenty kilometres, and the face-width-over-separation ratio runs from
0.010 to over two thousand. Great-circle distance through `arccos` of the dot
product returns exactly zero on 4,014 of those faces and takes any operator
built on it to NaN; the half-chord form through `arcsin` does not.

## `cell_area` is a different quantity from the area those faces bound

| | sum, as a fraction of `4 pi R^2` |
| --- | --- |
| exact spherical Voronoi area | 1.00000000 |
| the export's `cell_area` | 1.00068 |

Per region the two disagree by much more than that total suggests. The ratio has
a median of 1.011, a first percentile of 0.66 and a ninety-ninth of 2.71, runs
from 0.26 to 31.7, and lies within 1% on only 3.5% of the mesh.

**What this does and does not mean.** Global area-weighted means are safe: the
totals agree to 0.068%, so any quantity integrated over the whole sphere or over
a large region is unaffected. What is not safe is a per-cell area, and what is
least safe is a finite-volume divergence, which is only consistent when taken
over the area its own faces bound. Using `cell_area` there rather than the
Voronoi area put the discrete Laplace-Beltrami operator 0.57 relative RMS off
its analytic eigenvalue instead of 0.11.

So the water table solver uses both on purpose: fluxes divide by the Voronoi
area, and water volumes multiply by `cell_area`, because every other component
computes volumes that way and the solver's reduction identity against
`surface_water.py` has to be exact rather than close.

**What has not been established** is which of the two the exporter intends, or
whether any existing consumer is affected. `cell_area` is read for area
weighting across several components, and none of those uses is a divergence, so
none is wrong in the way this one would have been. Auditing them is GW-7 and it
was not done here.

## The operator's own error, which is not a geometry problem

With the true face widths and the Voronoi area, the discrete operator applied to
Legendre polynomials returns the analytic eigenvalue `-l(l+1)/R^2` to a relative
RMS of 0.022 at `l = 1` and 0.108, 0.120 and 0.124 at `l = 2, 3, 4`. The
criterion declared before the run was 0.10, so the first passes and the rest
miss.

The miss is the scheme's own truncation error and not a handful of bad faces.
The worst ten thousand cells carry under a tenth of the squared error, and
excluding every cell that touches a sliver face moves the RMS by 0.002. A
two-point flux approximation is first-order and exactly consistent only on a
mesh whose faces are perpendicular to the lines joining generators with
comparable cell sizes either side; this mesh has neither. About 12% is therefore
the mesh-scale noise floor of any field solved with this operator, and it is
reported as such rather than tuned away.

It is small against what else the water table carries: Gleeson's within-class
permeability spread is 1.5 to 2.5 orders of magnitude.

## The divide test as declared is mis-specified, and missed

Recorded because the miss is the test's fault rather than the solver's, and
because a future attempt should not re-run it in this form.

The declared criterion was that with uniform permeability and a terrain-
following water table, the groundwater catchments reproduce the surface
catchments: at least 90% of land area overall and 99% outside filled
depressions. Measured: **73.0% and 71.5%**. It misses, and the second number
being no better than the first is the tell.

The two labellings are built differently, and the difference has nothing to do
with groundwater. The groundwater trace hands each region to the neighbour
taking the largest outgoing flux, which under a terrain-following table is a
face-width-weighted steepest descent on the RAW surface. `regions.nc`'s
`terminal` comes from a priority flood on the FILLED surface, and
`hydrography/README.md` says plainly why it has to: a filled depression is flat
and has no downhill neighbour, so steepest descent would drop whole tributaries.

So the test compares a local descent rule against a global flood rule and
attributes the difference to the divides. It cannot isolate what it was declared
to isolate, and no version of the solver would pass it.

What would test the same claim honestly is a comparison against a labelling
built the same way as the trace -- the flux trace on the filled surface, against
the flood's own `receiver` pointer -- so that the only thing differing between
the two sides is whether the flow was routed above ground or below it. That is
GW-10 and it was not done here.

## The check that catches a wrong face width, and the ones that do not

Worth recording because the first attempt at this got the face width wrong and
almost every check passed anyway.

That attempt estimated each face's width from cell area and neighbour count,
treating a region as a regular n-gon and giving every face an equal share of its
perimeter. It is exact on a uniform tessellation. Here it put the operator 178x
off, with the correct answer in the MEDIAN and a spread of two orders either
side, because a face shared with a small neighbour was given the same width as
one shared with a large one.

Closure did not catch it, and could not: recharge in equals discharge out for
any symmetric face weights whatsoever, however wrong. The reduction identity did
not catch it either, because at zero permeability no face conducts and the
weights never enter. Only the eigenvalue test could, because it is the only one
with an answer fixed in advance by something outside the model.
