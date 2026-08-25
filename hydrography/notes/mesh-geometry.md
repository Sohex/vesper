# The region mesh has a Voronoi dual, and `cell_area` is not its area

This is worldbuilding. Vesper is a fictional super-Earth and this document is
about the geometry of the mesh its terrain is exported on.

Every figure below is mesh geometry and carries no climate, so the forcing does
not enter it. Where a number elsewhere in this component comes from a solved
water table, it was forced by the BOOTSTRAP climatology, not a baseline: this
build has no `baseline_climatology` and `surface_water.nc` was forced the same
way. By this project's own vocabulary a bootstrap run's numbers are not the
baseline, so those are provisional in a way these are not.

Measured 2026-08-20 on `precarve-craton`, while building the water table solver,
which is the first thing in this project to need face widths rather than only
cell areas. `hydrography/scripts/groundwater.py:Geometry` is the code. The
section on what the operator's error is an error IN was measured 2026-08-25, on
the same build and on synthetic meshes from Orogen's own generator, and says so
again where its tables are.

## The mesh is Voronoi, and the dual recovers what the export does not carry

The export gives a CSR neighbour list and a `cell_area` per region. A
finite-volume scheme needs two things neither of those supplies: the width of
the face shared by two neighbours, and the distance between their generators.

Both follow from the generator positions. The region adjacency is a Delaunay
triangulation, so the convex hull of the generators recovers it, each Voronoi
vertex is a triangle's circumcentre, and the face between two neighbours is the
arc joining the circumcentres of the two triangles sharing that Delaunay edge.

Note what this does and does not say about the export. The ADJACENCY is
Delaunay and is shared. The CELLS Orogen itself draws are not Voronoi cells:
`cell_area` is the centroidal dual, for the reason the next section gives. The
Voronoi dual is reconstructed here, in the solver, and nothing upstream is
changed by doing so.

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

## `cell_area` is the CENTROIDAL dual, and one line of Orogen says so

The two areas differ because they are duals of different point sets.

`vendor/orogen/js/sphere-mesh.js:206`, `generateTriangleCenters`, returns the
arithmetic mean of each triangle's three vertices. That is the CENTROID, not the
circumcentre, and the Voronoi vertex of a Delaunay triangle is its circumcentre.
The comment immediately above the function reads "Triangle centres (= Voronoi
vertices on the sphere)", so the code contradicts its own docstring; the centres
are also not renormalised onto the sphere.

Reconstructing the centroidal dual here reproduces `cell_area` essentially
exactly -- median ratio 1.000000, and 99.992% of regions within 0.1% -- which
settles what `cell_area` is. It is Orogen's centroidal dual, faithfully.

**The centroidal dual does not tile the sphere and the Voronoi dual does.**
Measured on this export, summing spherical excess over the per-face triangles
`(generator, corner_1, corner_2)`:

| corners | sum / `4 pi R^2` |
| --- | ---: |
| centroid | 1.00068086 |
| circumcentre | 1.00000000 |

So the 0.068% by which `cell_area` exceeds the sphere is not a mystery and not a
radius discrepancy: it is the centroidal dual's own tiling error, and the
centroid reconstruction reproduces both the per-cell values and that total.

The intuition that any one point per triangle tiles is close to a true statement
but not this one. Splitting each triangle from an interior point into three
sub-triangles does tile it, for any interior point. The dual used here is a
different decomposition -- the region of generator `i` is the union over its
faces of `(p_i, corner_1, corner_2)`, which spans BOTH triangles either side of
each face. That union tiles only when the corners are circumcentres, because
only then is it the Voronoi cell. Circumcentres lying outside their own triangle
does not break it: 36.5% of them do here, and the sum is still 1.00000000,
because a Voronoi cell is convex and contains its generator, so it is star-shaped
from `p_i` and the per-face triangles cannot fold.

## Only the Voronoi face is perpendicular, and that is why it must be used

This is the part that decides the question, and it is not about area at all.

A two-point flux approximation estimates the flux through a face as
`width * T * (h_j - h_i) / length`. That is the true flux only when the face is
PERPENDICULAR to the line joining the two generators, which is the
K-orthogonality condition. Measured over all 7,499,997 faces:

| corners | median face-to-generator angle | faces within 1 degree of perpendicular |
| --- | ---: | ---: |
| centroid | 75.746 deg | 3.8% |
| circumcentre | 90.000 deg | 100.00% |

The Voronoi face lies on the perpendicular bisector by construction. The
centroidal face does not, and the operator built on it is not a discretisation
of the Laplacian at all. Against the analytic eigenvalue:

| faces | area | l = 1 | 2 | 3 | 4 |
| --- | --- | ---: | ---: | ---: | ---: |
| centroid | centroid | 164.91 | 95.26 | 67.45 | 52.20 |
| centroid | `cell_area` | 164.91 | 95.26 | 67.45 | 52.21 |
| circumcentre | circumcentre | **0.0217** | **0.1082** | **0.1202** | **0.1243** |
| circumcentre | `cell_area` | 0.5686 | 0.5916 | 0.5954 | 0.6000 |

Centroidal faces are off by two orders of magnitude, in the same range as the
hand-estimated face width recorded below. So the choice is forced: face widths
come from the Voronoi dual because the scheme requires K-orthogonal faces, and
fluxes divide by the Voronoi area because that is the area those faces bound.

**Nothing about Orogen changes.** The Voronoi dual is reconstructed in the solver
from the export's own generator positions; `cell_area` is read and not rewritten.
It stays the denomination of water volumes, which keeps the reduction identity
against `surface_water.py` exact and leaves Orogen's own basin accounting --
`cell_area` grouped by `drainage_terminal` reproducing `finalCatchment.areaKm2`,
and `cellArea` feeding `detectBasins`, `selectBasins` and `attachHypsometry` --
untouched. No new build, so rule 7 does not bite.

## The operator's own error, and what it is an error IN

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
comparable cell sizes either side; this mesh has neither.

It is small against what else the water table carries: Gleeson's within-class
permeability spread is 1.5 to 2.5 orders of magnitude.

**That figure is a TRUNCATION error: the operator applied to a field that is
already known.** Nothing downstream reads one. What the water table hands its
consumers is the SOLUTION of `L h = f`, and the two are neither the same size
nor the same kind of quantity. The rest of this section is that difference,
measured, because reading the 12% as the noise floor of a solved field would
put every effect below it out of reach before it was built, and GW-18 and GW-24
both propose effects of a few per cent.

`groundwater.py --instrument` is the measurement. Its criteria were committed
before its first run; everything below was measured 2026-08-25.

### The truncation error does not fall as the mesh refines

Orogen's generator is a jittered Fibonacci spiral, so a mesh of this family can
be built at any resolution without a build. `groundwater.py:orogen_mesh`
reproduces `vendor/orogen/js/sphere-mesh.js:generateFibonacciSphere` line for
line and takes the jitter from the active build's own manifest, so the sweep
varies region count and nothing else.

| regions | l = 2 | l = 3 | l = 4 |
| ---: | ---: | ---: | ---: |
| 20,000 | 0.10642 | 0.11833 | 0.12325 |
| 80,000 | 0.10548 | 0.11744 | 0.12229 |
| 320,000 | 0.10552 | 0.11745 | 0.12221 |
| 1,280,000 | 0.10549 | 0.11725 | 0.12192 |
| 2,500,001, `precarve-craton` | 0.10819 | 0.12015 | 0.12427 |
| fitted order in the cell size | 0.004 | 0.004 | 0.005 |

The declared bar was an order below 0.5 and the measurement is three decimal
places from zero: **the truncation error does not converge at all.** Sixty-four
times the region count moves it by under a per cent, and twenty thousand
regions carry the same 12% as two and a half million.

That is why GW-3's Earth mesh reproduced these numbers at the same cell size,
and it says more than GW-3 could: this is a property of the jittered-Fibonacci
Voronoi family and not of a mesh, a planet or a resolution. The reason is the
geometry above. The Voronoi face is perpendicular to the segment joining its two
generators, which is what makes the two-point flux the true flux through the
point where that segment crosses it; the face is not CENTRED on that point, and
the offset is a fixed fraction of the cell whatever the cell's size. Refining
the mesh shrinks the cell and the offset together and their ratio never moves.

### The solved field's error does fall, at second order

The same meshes, handed the analytic right-hand side `-l(l+1)/R^2 P_l` and asked
for the field, scored against `P_l`. Relative RMS, with the constant a closed
sphere leaves free removed from both sides:

| regions | l = 2 | l = 3 | l = 4 |
| ---: | ---: | ---: | ---: |
| 20,000 | 3.441e-04 | 7.056e-04 | 1.213e-03 |
| 80,000 | 8.751e-05 | 1.780e-04 | 2.985e-04 |
| 320,000 | 2.183e-05 | 4.434e-05 | 7.459e-05 |
| 1,280,000 | 5.453e-06 | 1.108e-05 | 1.860e-05 |
| 2,500,001, `precarve-craton` | 2.766e-06 | 5.629e-06 | 9.437e-06 |
| fitted order in the cell size | 1.994 | 1.998 | 2.008 |

Second order, against a declared bar of first, and already three orders below
the truncation error at the coarsest mesh of the sweep. This is
supraconvergence and it is not a surprise: the inverse of an elliptic operator
damps the grid-scale part of a residual by the square of the scale ratio, and a
truncation error that this note had already measured to be distributed noise
rather than a few bad faces is exactly the part it damps.

So the operator's error and the head's error are different quantities that
behave differently under refinement, and only the second is a floor on anything
this component reports. Carried across to the active build's 10,000,005 regions
at the sweep's own fitted order, it is about 7e-07.

The declared bar on the linear solve is a relative residual of 1e-10 and it is
met everywhere in the sweep. On `precarve-craton` it is met at `l = 3` and
`l = 4` and MISSED at `l = 2`, where the residual is a difference of face
fluxes over 7,499,997 faces and has a round-off floor of its own that grows
with the mesh. That value is therefore not certified by the declared criterion.
It is not orphaned either: the sweep's own second-order trend predicts 2.80e-06
there against the 2.766e-06 measured.

### The error is a fraction of the signal, not an absolute floor

The other half of the same misreading is that a 12% relative error would swamp a
1% effect. It does not, because it is RELATIVE. Asked for the RESPONSE to a
transmissivity perturbation rather than for a field -- which is what a spatially
varying aquifer thickness is -- the operator can be scored against an analytic
answer of its own: for zonal fields on a sphere,

    div(T grad u)  =  (1/R^2) d/dz [ T(z) (1 - z^2) du/dz ]

so with `T = 1 + a P_m(z)` and `u = P_l(z)` the response per unit amplitude is
`(1/R^2) d/dz [P_m (1 - z^2) P_l']`, a polynomial evaluated exactly. On the
1,280,000-region mesh, at `l = 2`, `m = 3`:

| perturbation amplitude | error as a fraction of the response | of the field |
| ---: | ---: | ---: |
| 0.010 | 0.12193 | 8.02e-04 |
| 0.200 | 0.12193 | 1.60e-02 |

The error is the same 12%, of the response, at both amplitudes. It has to be:
the confined assembly is LINEAR in the transmissivity, so the response per unit
amplitude cannot depend on the amplitude. **A perturbation of relative size `a`
is therefore carried with an absolute error of about `0.12 a`, and a small
effect is measured with the same relative accuracy as a large one rather than
being lost in the operator's error.** GW-24's unconfined form is the case that
does NOT get this argument, because there `T` depends on the head and the
assembly is not linear in it.

### The number to quote, and what it does not cover

The 12% is the operator's and belongs in a statement about the operator. The
head's own mesh error is the second table, and it is small enough that the
discretisation is not what limits any effect this component has proposed:
GW-24's unconfined correction is 0.48% of the head at the 2 km thickness and
14% at this world's sourced median, and GW-18's thickness field is a multiplier
on the transmissivity spanning a factor of thirty. Both are orders above the
floor, and what limits them is what `notes/subgrid-water-table.md` sections 3
and 4 already measure: the model's own between-cell spread, and a run.

**This bounds the mesh error on the part of the field the mesh RESOLVES, and
nothing else.** A forcing with structure at the cell scale is not made
resolvable by a converging scheme; it is absent. That is GW-6, it is a standing
constraint rather than an error, and `notes/subgrid-water-table.md` is where it
is argued.

One note on running the arm, because it decides whether it runs at all. The
solve is over the whole mesh rather than over a free set, so the factorisation
is the cost, and `MMD_AT_PLUS_A` on its own is not enough: SuperLU pivots for
stability by default and that undoes the symmetric ordering it was just given.
On the 80,000-region mesh the factor is 6.8 million nonzeros with
`diag_pivot_thresh=0` and `SymmetricMode`, against a 2.5 GB working set and no
answer without them. The matrix is symmetric positive definite, so no pivoting
is the correct setting rather than merely the fast one.

**Do not carry that across to `solve` without measuring it there.** `solve`'s
matrix is a free set inside a Dirichlet boundary, not a whole mesh with a single
cell pinned, and the boundary coupling can leave it dominant enough that partial
pivoting picks the diagonal regardless. On a synthetic of that shape at 55,977
free cells, with a spatially coherent conductivity spanning Gleeson's 3.4 orders,
the factor is 1,042,714 nonzeros and a 0.55 GB peak IDENTICALLY with and without
the symmetric setting, so at that size it buys nothing. A four-times larger
synthetic did not finish its unsymmetric arm, so whether there is a crossover
above it is unmeasured, and `solve` is therefore left exactly as it is. What it
would be worth there is the one thing WORLD-V3UR's memory blocker has not been
tried against.

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
