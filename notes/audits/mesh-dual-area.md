# Two duals, two areas: what `cell_area` is and what it is not

Worldbuilding. Vesper is an invented planet and this note is about the
simulation of it: the mesh a terrain generator exports and the area weight
components integrate over. Every quantity named here is a modelled field.

Audited 2026-08-20 against the vendored Orogen fork and the `precarve-craton`
export, after the groundwater solver's GW-7 -- the first thing in this project
ever to need an area for a single region rather than for a sum of them.

**A first version of this note asserted that the centroidal dual does not tile
the sphere and that this explained a 0.068% excess. Both claims were tested
afterwards and both are false.** What follows is what survived the testing.

## 1. What is actually there

Orogen's regions are the cells of a spherical Voronoi tessellation of the
generators: the region adjacency is exactly the Delaunay edge set, agreeing on
all but 47 cospherical quadrilaterals out of 7,499,997 undirected pairs.

`regionCellArea` (`vendor/orogen/js/geometry.js:41`) sums the spherical excess
of each sub-triangle of the dual polygon around a region. The arithmetic is
exact. The polygon's corners come from `generateTriangleCenters`
(`vendor/orogen/js/sphere-mesh.js:206`), which returns the **centroid** of each
triangle -- the arithmetic mean of its three corners, never renormalised onto
the sphere. A Voronoi vertex is the **circumcentre**.

So `cell_area` is the area of the CENTROIDAL dual cell, and the region is the
VORONOI cell. Two different polygons around the same generator.

**The one unambiguous defect is documentation.** The comment above
`generateTriangleCenters` reads `// Triangle centres (= Voronoi vertices on the
sphere).` They are not. `regionCellArea`'s own docstring is honest about it and
calls it "the same approximation the renderer draws", so the fork knew; the
comment upstream of it is what makes the choice invisible.

## 2. Both duals tile the sphere, and that is not the problem

Tested directly, building spheres at four resolutions and summing
`regionCellArea` under each set of corners:

| regions | centroid corners | circumcentre corners |
| ---: | ---: | ---: |
| 20,001 | 1.00000000 | 1.00000001 |
| 200,001 | 1.00000000 | 1.00000028 |
| 1,000,001 | 1.00000000 | 1.00000283 |
| 2,500,001 | 1.00000000 | 1.00002906 |

The centroidal dual tiles the sphere EXACTLY, at every resolution. That is not
an accident: any choice of one point per triangle tiles, because two adjacent
regions share the two triangles either side of their common edge and therefore
share that edge of the dual.

The interesting column is the second one. A circumcentre reconstruction
degrades as the mesh refines, by three orders of magnitude between 20k and 2.5M
regions. Slivers put a circumcentre outside its own triangle, the polygon
tangles, and the `Math.abs` in the excess sum then adds a fold positively
instead of cancelling it. **A naive circumcentre dual is numerically worse at
this project's resolution than the centroidal one it would replace.**

## 3. The export's 0.068% is a radius question, not a dual question

`sum(cell_area)` over the `precarve-craton` export is 0.068% above `4 pi R^2` at
the manifest's `planet.radiusKm` of 7645.2 km, implying 7647.8 km instead. A
freshly built sphere at the same region count sums exactly, so the dual does not
explain it. Unexplained, small, and recorded as its own row rather than folded
into this one.

## 4. What the two areas cost, by aggregation

Measured per region against a circumcentre reconstruction. The errors are
essentially uncorrelated, so they cancel under aggregation:

| aggregation | centroidal against Voronoi |
| --- | --- |
| whole sphere | tiles either way |
| per basin catchment, 3,621 basins | 1st to 99th percentile 0.973 to 1.029 |
| per cell | 5th to 95th percentile 0.462 to 1.367; range 0.032 to 3.81 |

Per cell, 69% differ by more than 10% and 34% by more than 25%. Per basin, only
0.22% differ by more than 5%.

## 5. `cell_area` is internally consistent, and that is the point

The export's own metadata calls it "a MESH property -- raw output only" and says
"do not area-weight a gridded field with it; use `grid_cell_area`", which is the
field documented to sum to the planet's surface area. It also states that
grouping `cell_area` by `drainage_terminal` "reproduces `finalCatchment.areaKm2`
exactly", and `basins.drainageConsistency` publishes that check.

So `cell_area` is not a loose approximation of a Voronoi area. It is the area
Orogen's own basin catalogue, hypsometry and drainage accounting are computed
in, with a published identity tying them together. **Replacing it would break
that identity unless every consumer inside Orogen moved with it**, and
`cellArea` also feeds `detectBasins`, `selectBasins` on its `minAreaKm2`
threshold, `attachHypsometry` and lithology -- so it conditions the terrain, and
changing it is a new build under rule 7.

## 6. What to do instead

**Nothing in this repository is currently wrong.** Every existing consumer
aggregates, where the difference cancels, and the gridded weight is a different
field that sums exactly by construction.

The exposure is per-cell use, and there is exactly one: a finite-volume
divergence must be taken over the area its own faces bound. The groundwater
solver reconstructed Voronoi faces and then had to choose an area to pair them
with.

**The cheap and correct resolution is to move the solver, not the generator.**
Derive the face widths from the CENTROIDAL dual, so that areas and faces come
from the one tessellation Orogen already uses and already checks. That costs no
build, breaks no identity, and leaves the terrain untouched. Changing Orogen's
dual would mean moving its basin accounting with it, taking a new build, and
first building a circumcentre routine robust enough not to lose accuracy at 2.5M
regions -- which the naive one does.

What should change in the fork is the comment that made this invisible, and
nothing else.
