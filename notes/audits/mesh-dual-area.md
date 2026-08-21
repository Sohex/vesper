# Two duals: which one tiles, which one the operator needs, and why they differ

Worldbuilding. Vesper is an invented planet and this note is about the
simulation of it: the mesh a terrain generator exports and the geometry a
finite-volume solve on it requires. Every quantity named here is a modelled
field.

Audited 2026-08-20 against the vendored Orogen fork and the `precarve-craton`
export, after the groundwater solver became the first thing in this project to
need a single region's geometry rather than a sum over many.

**Two earlier versions of this note were wrong, in opposite directions, and the
history matters because it says how to test this.** The first asserted that a
centroidal dual cannot tile the sphere. The second retracted that after building
fresh spheres and finding the centroidal dual summing exactly -- but those were
near-regular meshes on which centroid and circumcentre nearly coincide, so the
test could not see the effect it was aimed at. **Measure on the export, not on a
synthetic mesh.**

## 1. What is there

Orogen's regions are the cells of a spherical Voronoi tessellation of the
generators: the region adjacency is exactly the Delaunay edge set, agreeing on
all but 47 cospherical quadrilaterals out of 7,499,997 undirected pairs.

`regionCellArea` (`vendor/orogen/js/geometry.js:41`) sums the spherical excess
of the dual polygon around each region, exactly. Its corners come from
`generateTriangleCenters` (`vendor/orogen/js/sphere-mesh.js:206`), which returns
each triangle's **centroid**. A Voronoi vertex is the **circumcentre**.

The comment above that function reads
`// Triangle centres (= Voronoi vertices on the sphere).` They are not, and that
is why the choice was invisible for as long as it was.

## 2. Only one of them tiles

Measured on the export, both duals against the same radius:

| corners | area sum, km2 | ratio to 4 pi R^2 |
| --- | ---: | ---: |
| circumcentre | 7.3449283955e+08 | **1.0000000000** |
| centroid (`cell_area`) | 7.3499301038e+08 | 1.0006809744 |

The circumcentre dual reproduces the sphere at the manifest's own
`planet.radiusKm` of 7645.2 km to ten significant figures. The centroidal dual
does not tile, and its 0.068% excess is its own tiling error.

**Why the intuition that "any interior point tiles" fails here.** Splitting one
triangle from an interior point does tile that triangle. But this dual gives
generator `i` the union over its faces of `(p_i, corner1, corner2)`, and each
such piece spans BOTH triangles either side of that face. The three pieces meeting
in a triangle reassemble it only when the corner is the circumcentre.

A related trap in the other direction: 36.5% of circumcentres lie outside their
own triangle on this mesh, which looks like it should fold the polygon and
inflate the sum. It does not, because a Voronoi cell is convex and contains its
generator, so no per-face sub-triangle can fold.

## 3. The area is not the reason the solver needs Voronoi faces

A two-point flux approximation is only valid on K-orthogonal faces -- the face
normal parallel to the line joining the two generators. That is a property of the
FACES, not of the areas, and it is where the two duals are not close:

| corners | median face-to-generator angle | within 1 degree of perpendicular |
| --- | ---: | ---: |
| centroid | 75.746 deg | 3.8% |
| circumcentre | **90.000 deg** | **100.00%** |

Built on centroidal faces the discrete operator misses its analytic Legendre
eigenvalue by factors of 165, 95, 67 and 52 at l = 1 to 4. On Voronoi faces it
misses by 0.11. Voronoi faces are perpendicular bisectors by construction, so
they are exactly K-orthogonal and nothing else is.

**So the solver reconstructs Voronoi faces from the export's own generators.**
It reads `cell_area` and never rewrites it.

## 4. `cell_area` must not change

The export metadata calls it "a MESH property -- raw output only" and directs
gridded weighting to `grid_cell_area` instead. It also states that grouping
`cell_area` by `drainage_terminal` reproduces `finalCatchment.areaKm2` exactly,
with `basins.drainageConsistency` publishing that check. It is the denomination
of Orogen's own basin catalogue, hypsometry and drainage accounting.

`cellArea` further feeds `detectBasins`, `selectBasins` on its `minAreaKm2`
threshold, `attachHypsometry` and lithology, so it conditions the terrain and
changing it is a new build under rule 7.

**Nothing in this repository is currently wrong.** Every existing consumer
aggregates, and per basin the two duals agree to 1.5% with only 0.22% of basins
beyond 5%. Per cell they differ by more than 10% on 69% of regions, which is the
exposure, and it has exactly one consumer.

## 5. What a change to the fork would be

The comment, and nothing else. Correcting `generateTriangleCenters`' claim that
its output is Voronoi vertices costs nothing and removes the trap. Changing the
areas would mean moving Orogen's basin accounting with them and taking a new
build, for a quantity every current consumer aggregates away.
