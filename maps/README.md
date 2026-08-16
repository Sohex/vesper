# maps/

Illustrative colour maps of Vesper in seven projections. These are pictures, not
analysis products: they exist to show what the world looks like.

```bash
python maps/build_basemap.py       # equirectangular base, ~25 s (lookup cached)
python maps/render_projections.py  # the seven projections, ~5 min at width 3000
```

| File | Projection | Property |
| --- | --- | --- |
| `vesper_winkel_tripel.png` | Winkel tripel | compromise; the usual world map |
| `vesper_equal_earth.png` | Equal Earth | equal-area pseudocylindrical |
| `vesper_lambert_azimuthal.png` | Lambert azimuthal | equal-area, east/west hemispheres |
| `vesper_lambert_azimuthal_polar.png` | Lambert azimuthal | equal-area, north/south hemispheres |
| `vesper_waterman_butterfly.png` | butterfly | octahedral, interrupted |
| `vesper_dymaxion.png` | Dymaxion | icosahedral, interrupted |
| `vesper_authagraph.png` | AuthaGraph style | tetrahedral, shape-first, rectangular |

`projections_provenance.json` records the build, the climatology, and every
orientation the search chose.

## What the colour means

Terrain is real: elevation, land mask, lithology and endorheic flags come from
the active build's native 2.5M-region mesh, one nearest region per pixel, so
coastlines and dry closed-basin floors are exact. Land and sea are taken from
`surface_class`, never from `land_mask`, so the sub-sea-level basin floors stay
land.

Colour is illustrative:

- Land tint is the biome class from `exoplasim/analysis/climatology_s096`. That
  field is T42 and was computed on the pre-carve terrain, so it is a plausible
  tint over the current geography rather than a result about it.
- Bare rock takes over above about 1.8 km, since a T42 cell cannot see a
  mountain.
- Evaporite crust and playa fill are painted from `surface_rock`, which is why
  the closed basins read pale.
- Permanent snow is where the warmest month falls below freezing after a
  6.5 K/km correction from the T42 orography to the mesh orography. There is no
  glacier model behind it; `glac` is zero everywhere in the baseline run.
- Sea ice is the annual-mean fraction.
- Relief is a hillshade of the mesh elevation, sun from the north-west, with
  vertical exaggeration.

The poles are forested rather than icy. That is not a rendering fault: at 32
degrees obliquity the polar summer reaches about +20 C in the baseline run, and
both polar caps are land.

## The graticule

Every 30 degrees, with the equator and the prime meridian drawn heavier. It is
not projected as geometry; each pixel measures its own distance to the nearest
line and shades accordingly, which is what lets the same code draw a correct
graticule on a butterfly or a tetrahedral rectangle without knowing anything
about either. Three details make it read:

- distance is divided by the local gradient of latitude or longitude, so it
  comes out in pixels and a line keeps the same weight however hard the
  projection is stretching there;
- coverage ramps across the last pixel rather than switching, so lines are
  antialiased rather than jagged;
- a dark halo is laid under the light line, so it stays legible over deep ocean
  and bright desert alike.

Lines are suppressed where latitude or longitude jumps by more than 20 degrees
in a pixel, which is a cut rather than a parallel. The earlier rule was
relative to the median gradient, and that also erased lines wherever a
projection was legitimately compressed, which is why the graticule used to
fade out in patches. Each render block is computed one row proud at each end
and trimmed, since the vertical derivative would otherwise break at every
block boundary.

## Orientation is derived, not assumed

Every projection here has a free orientation, and each one is searched against
the land mask rather than left at zero:

- central meridian, so the antimeridian crosses the least land;
- hemisphere boundary for the equatorial azimuthal pair, likewise;
- the butterfly's equatorial vertices, scored on its two cut meridians;
- the icosahedron's attitude, scored on all 30 edges, since every cut the net
  can make runs along one of them;
- the tetrahedron's attitude, on three counts at once. Its border is the map's
  only interruption, so land on it is torn; but the distortion pattern is also
  fixed in the rectangle, worst at the four side midpoints where the
  tetrahedron's vertices sit, so turning the tetrahedron decides which parts of
  the world land in it. Scoring the border alone, as the first version did,
  left exactly as much distortion over land as chance would give. The search
  now holds the border and gross size error to a standard and takes the best
  shape among what qualifies: **12% border, 4% of land area grossly
  mis-scaled, shape 1.32**. Random attitudes only resolve SO(3) to about six
  degrees, too coarse for the narrow band that meets all three, so the best
  few are refined with a shrinking step; refinement descends a penalised cost,
  because an attitude usually has to cross infeasible ground to reach the good
  region, and the winner is then picked by the standards themselves. Turning
  the tetrahedron rotates the whole map rigidly, so the samples are mapped to
  the sphere once and merely rotated per trial, which makes 40,000 attitudes
  cheap.

The Dymaxion net itself is derived too. Fuller chose his net so that Earth's
continents came out whole; the same argument applied to this world gives a
different net, so `choose_net` runs Kruskal over the face graph with each shared
edge weighted by the land along it, keeps the heavy edges inside the tree, and
scores candidates on land torn, bounding-box fill and aspect. Nets that
self-overlap are discarded, because a net has to lie flat.

With 43% land there is no clean cut anywhere: the best antimeridian available is
still 31% land.

The polar plate is the exception, and deliberately so. Its cut is the equator,
which is fixed geometry and 40% land, so nothing is searched; the only free
parameter is which meridian points down, and it reuses the central meridian.
North is on the left with that meridian at the bottom, south on the right with
its opposite at the top, so the two rims carry the same longitude where they
meet.

## The tetrahedral rectangle

The AuthaGraph-style map is built from two independent pieces, both checkable.

**The rectangle.** A flat tetrahedron is the plane folded by half turns about
the points of a triangular lattice, so its development tiles the plane, and any
fundamental domain of that symmetry is a complete world map. One such domain is
a 1 by sqrt(3) rectangle whose four side midpoints are the four vertices of the
tetrahedron: the pillowcase. It is continuous everywhere inside and cut only
along its border. `TetrahedralRectangle` develops the tiling by walking the
lattice, keeping the two shared vertices at each step and giving the new corner
the one vertex the current face does not use, and asserts that the walk closes
consistently.

**The faces.** One tetrahedron face carries a quarter of the sphere, and a
single smooth mapping of something that large has to absorb all of its
curvature at once. Snyder's equal-area polyhedral projection does it by
shearing: it cuts the face into six sub-sectors, remaps azimuth so equal
fractions of area lie either side of it, then radius so equal fractions lie
inside it. It is exactly equal-area, and it is implemented here, but measured
against its own Jacobian it puts a median anisotropy of **1.86** across the
whole map. That is not a few bad corners; it is everywhere, and it is what made
the first version of this map look smeared.

So the faces are subdivided instead, which is the idea behind AuthaGraph's 96
triangles. Each face is cut into a regular mesh of small triangles and each one
is mapped affinely onto its spherical counterpart. A small triangle is nearly
flat, so it can be nearly similar to its image; the curvature is taken up as
small kinks along the seams rather than as shear across the whole face.

Where the small triangles' corners go is then a free choice, and it is made by
optimising, over the sub-triangles,

    sum (s1/s2 + s2/s1 - 2)  +  w * sum (log(area / target))^2

shape against area. Only the multiset {i, j, k} of a mesh point is free,
because a face has the symmetry of its three corners, which keeps every face
identical and every shared edge palindromic and so keeps the four faces
agreeing where they meet. That leaves about twenty parameters. The optimiser
starts from Snyder, which gives it a sensible mesh to begin from, and the
result is cached under `build/`.

Shipped at 16 divisions a side, 256 triangles per face and 1024 in all, with
`w = 0.5`. Restarting the optimiser from perturbed meshes lands on the same
energy to four decimals, so this is the frontier and not a local minimum.

## Why this one is tuned to look right rather than to be right

Exact equal area is a bad buy on a tetrahedron, and the reason is forced rather
than incidental. A face's planar edge is longer than the spherical edge it
comes from, so the boundary has to stretch along its length, and preserving
area compresses it across by the same factor. Both numbers are fixed by the
polyhedron alone:

| | planar edge | spherical edge | stretch | forced anisotropy |
| --- | --- | --- | --- | --- |
| tetrahedron | 2.6935 | 1.9106 | 1.410 | **1.99** |
| octahedron | 1.9046 | 1.5708 | 1.213 | 1.47 |
| icosahedron | 1.2046 | 1.1071 | 1.088 | 1.18 |

That is a mean along the edge, so Jensen makes it a lower bound, not an
estimate. An exactly equal-area tetrahedral map therefore runs from about 1.0
at the four face centres to about 2.0 along every face boundary, averaging near
1.6, and measurement agrees: Snyder gives 1.86 and the best equal-area mesh
about 1.7. There is no construction that does better.

Since these maps exist to be looked at, the map steps off that constraint
deliberately. Lowering `area_weight` lets the optimiser shrink the four corner
regions, which buys shape everywhere else, and the attitude search then turns
those corners out to sea. Measured over land, with each row's attitude searched
for it:

| `area_weight` | shape (anisotropy) | land area grossly mis-scaled | border |
| --- | --- | --- | --- |
| Snyder, smooth | 1.86 | 0% (exact) | -- |
| `w = 3` | 1.61 | 1% | 7% |
| **`w = 0.5`, shipped** | **1.32** | **4%** | 12% |
| `w = 0.25` | 1.22 | 6% | 13% |

"Grossly mis-scaled" means drawn at less than 0.6 or more than 1.67 of true
size, weighted by the area a region really covers rather than the room it gets
on the page. Mild size error is not visible on a map of an invented world;
shearing is, immediately. At the shipped setting the land runs 0.63 / 1.00 /
1.18 at the 5th, 50th and 95th percentiles of scale, so sizes still read
correctly, while shape beats every other whole-world projection here: Winkel
tripel measures 1.52 over land and Equal Earth 1.45.

Pushing further is possible and was tested. `w = 0.25` looks marginally better
again but hides a fifth of the world in 5% of the page, and buying the border
down from 12% to 6% costs 15% of land area at gross scale error. Both were
rejected: they trade a visible kind of wrongness for a different visible kind.

Passing `subdivision=None` gives the exactly equal-area Snyder version back,
and any `area_weight` can be rendered by changing one argument.

For reference, the interrupted maps here are barely distorted at all -- the
Dymaxion measures 1.06 over land and the butterfly 1.20 -- because their
polyhedra are the good ones in the table above. The tetrahedral rectangle is
buying a single uninterrupted rectangle, and that is what it costs.

Both constructions were verified rather than assumed: the six Snyder sub-sector
areas sum to pi, exactly the spherical area of a tetrahedron face; cell corners
land on the tetrahedron vertices to machine precision; 600k uniform points in
the rectangle cover the sphere once, with the smooth version uniform to within
Poisson noise; and sweeping every seam in the map -- sub-triangle edges, face
edges and sector rays -- finds no step larger than a smooth map would give.

What cannot be tuned away at all is the pinch at the four side midpoints. A
tetrahedron vertex has 360 degrees of sphere around it and only 180 degrees of
flat tetrahedron, so every map of this family halves the angle there: 2.0 at
best, and about 3.0 for a piecewise-affine one, over the 1.2% of the map
nearest those four points. The attitude search puts them over water for that
reason. AuthaGraph has the same character for the same reason.

## Approximations worth knowing

The butterfly and the Dymaxion are **gnomonic on their polyhedron**, which is
exact for the face geometry: gnomonic sends great circles to straight lines, so a
spherical triangle maps onto the planar triangle its vertices span, and the
inverse is barycentric. That is the standard construction behind the Fuller
projection, but it is not Fuller's own face mapping, and the butterfly is a
Cahill-style octahedral net rather than Waterman's truncated-octahedron
construction with its clipped corners. Shapes within a face differ slightly from
the published projections; the arrangement and the interruptions do not.

The tetrahedral rectangle is **not** AuthaGraph. Narukawa's projection
subdivides into 96 triangles by a method that has never been published, so this
reconstructs the concept from parts that can be checked: the same tetrahedral
basis, the same rectangular result, the same subdivide-and-map-each-piece
approach, the same approximately-rather-than-exactly equal area that approach
forces, and the same pinch where the map has to halve the angle at a vertex.
The subdivision here is finer than 96 and its mesh is solved for rather than
given, so the two will not agree triangle for triangle.
