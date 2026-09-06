# maps/

Illustrative colour maps of Vesper. These are pictures, not analysis products:
they exist to show what the world looks like.

```bash
python maps/build_basemap.py       # equirectangular base, ~30 s (lookup cached)
python maps/render_projections.py  # the authagraph, ~50 s at width 4000
python maps/render_projections.py --all           # all seven
python maps/snapshot.py --step surface_water      # both, filed as a frame
```

`--all` is minutes rather than a minute: the iterative inverses cost far more
per pixel than the tetrahedral rectangle's does, and the Winkel tripel alone
runs several times the whole default render.

The AuthaGraph-style rectangle is the only projection drawn by default. It shows
the whole world at once with the least distortion of shape of anything here, so
it is the one to look at when the question is what the world looks like; the
other six are behind their own flags, and each carries its own orientation
search, so drawing all seven by default charged every caller for six pictures
they had not asked for.

| Flag | File | Projection | Property |
| --- | --- | --- | --- |
| default, or `--authagraph` | `vesper_authagraph.png` | AuthaGraph style | tetrahedral, shape-first, rectangular |
| `--winkel-tripel` | `vesper_winkel_tripel.png` | Winkel tripel | compromise; the usual world map |
| `--equal-earth` | `vesper_equal_earth.png` | Equal Earth | equal-area pseudocylindrical |
| `--lambert-azimuthal` | `vesper_lambert_azimuthal.png` | Lambert azimuthal | equal-area, east/west hemispheres |
| `--lambert-azimuthal-polar` | `vesper_lambert_azimuthal_polar.png` | Lambert azimuthal | equal-area, north/south hemispheres |
| `--waterman-butterfly` | `vesper_waterman_butterfly.png` | butterfly | octahedral, interrupted |
| `--dymaxion` | `vesper_dymaxion.png` | Dymaxion | icosahedral, interrupted |
| `--all` | all seven | | |

A flag replaces the default rather than adding to it, so `--dymaxion` draws the
Dymaxion and nothing else, and `--authagraph --dymaxion` draws both.

Two widths, and both are sized for the mesh underneath rather than fixed, since
a 10M-region export carries detail a 5760-wide base map cannot show:
`render_projections.py --width` is the width of a projection, and
`build_basemap.py --width` the equirectangular grid it samples from.
`snapshot.py` drives both and so calls the second `--basemap-width`.

`projections.py` is the projection layer: the maths of each one, the nets, and
the graticule. `build_basemap.py` draws the equirectangular raster,
`render_projections.py` reprojects it into a frame, `frames.py` owns the frame
manifest, and `snapshot.py` drives the first two from a point in the pipeline.

## Where the output goes

Output is namespaced by build like every other component's data, and inside
that a frame is a UUID:

```
maps/data/<build>/INDEX.json          what each frame is. TRACKED
maps/data/<build>/<frame>/*.png       the projections drawn for it
maps/data/<build>/<frame>/provenance.json   every orientation the search chose
```

**A frame is one state of the world, and it is identified by what it was drawn
from, not by when it was rendered.** The map is drawn from four things -- the
terrain, the classification, the climatology and the lake solution -- and the
frame id is a UUID with the identity of those four beside it in `INDEX.json`.
A name built from the step or the label that produced it would separate frames
only along the dimensions it encoded, which is the collision CLAUDE.md rule 6
records against ExoPlaSim runs and LPJ-GUESS runs. `maps/frames.py` carries the
argument.

`INDEX.json` is tracked and the frames themselves are not: the pixels are
regenerable while the four inputs survive, and the index is the only record of
what a UUID was once they do not. Same rule as `exoplasim/runs/`, and the same
mechanism holds it: the row is written when the frame id is taken, before
anything is drawn, and every write of the index happens under a lock on
`INDEX.json.lock` beside it. Nothing rebuilds the index from a directory
listing, so a frame whose pixels are deleted keeps its row and can be redrawn
from the inputs the row names.

## A frame at every step that can move the picture

`python scripts/pipeline.py --plan <target>` names `maps/snapshot.py --step <id>`
after every step the map is reachable from, which is most of loop A: a carve
verdict changes the next terrain and the next terrain is the picture. The point
is a series -- the coastline the carve moved, the lakes the water balance
filled, the biomes the new climatology repainted -- so that a pass can be
watched changing rather than only its end state drawn.

That is affordable because a frame is keyed on its inputs. `snapshot.py`
rebuilds the base map only when the raster on disk is not the world the config
describes, and a step that moves none of the four inputs is recorded against the
frame that already exists rather than rendering a second copy of it. So the
number of frames is the number of times the world actually changed, and the
`steps` list on each frame says which steps it stood through. `--force` renders
anyway, which is what to reach for after a change to the RENDERING rather than
to the world.

The set of steps is derived from `config/pipeline.yaml` by
`scripts/pipeline.py:map_affecting` and is never listed twice; `snapshot.py`
refuses a step the graph says cannot reach the map.

## Naming the climatology

`build_basemap.py` tints from the climatology `config/planet.yaml` names, and
that key is null whenever the active build has no baseline yet, so the bare
command fails rather than falling back. Name the climatology that does exist and
the map draws off it, tint and lapse rate from the same file:

```bash
python exoplasim/scripts/analyze_climatology.py --label bootstrap
python maps/build_basemap.py --climatology \
    exoplasim/analysis/climatology/bootstrap_regular_climatology.nc
```

The classification the tint reads is `<label>_classification.nc`, which is why
`analyze_climatology.py` has to have run under that label first.
`render_projections.py` needs no such argument: it inherits the climatology, its
caveat and the frame's whole identity from `basemap_provenance.json`, so a frame
records the world the raster shows rather than the one the config currently
names. `snapshot.py` takes `--climatology` and passes it through.

## What the colour means

Terrain is real: elevation, land mask, lithology and endorheic flags come from
the active build's native mesh, one nearest region per pixel, so
coastlines and dry closed-basin floors are exact. Land and sea are taken from
`surface_class`, never from `land_mask`, so the sub-sea-level basin floors stay
land.

Colour is illustrative:

- Land tint is the biome class from whichever climatology the active build has
  produced. Where that climatology was computed on a different terrain than the
  one being drawn -- which is the usual case early in a cycle -- it is a
  plausible tint over the geography rather than a result about it, and the map
  says so in its own provenance block.
- Bare rock takes over above about 1.8 km, since a climate-grid cell cannot see
  a mountain.
- Evaporite crust and playa fill are painted from `substrate_class`, which is why
  the closed basins read pale.
- Permanent snow is where the warmest month falls below freezing after a
  lapse-rate correction from the climatology's own orography to the mesh
  orography. The rate is `lib/lapse.py:environmental_lapse_k_per_km`,
  measured at call time from the same climatology the tint comes from and
  in its warm season, so it is this world's rate and not Earth's. There is
  no
  glacier model behind it; `glac` is zero everywhere in the baseline run.
- Sea ice is the annual-mean fraction.
- Relief is a hillshade of the mesh elevation, sun from the north-west, with
  vertical exaggeration. Standing water takes the sea's muted shading instead,
  since a lake surface is flat and the land hillshade would emboss it.

**Lakes and rivers are the exception: those are a result, not a tint.** They
come from `hydrography/data/<build>/surface_water.nc`, which solves a
closed-basin water balance against the same baseline climatology and accumulates
the same water down the drainage network. How many basins hold water, how much
of the planet they cover, and what the largest river carries are all in
`world_state.json`. Inland water is drawn
a shade greener than the sea so a lake reads as a lake rather than as a bay that
lost its connection.

Rivers are drawn at the resolution the mesh has, one region across or about
15 km, so their width carries no information. Weight does: the blend follows
discharge rather than the line getting fatter, and anything drawn at all is
drawn solidly, because a line one region wide at 20% opacity is a smudge rather
than a river.

Overflowing lakes have rivers leaving them: `basins.nc` carries the saddle
each basin spills at, so every overflowing basin's outflow is routed onto the
mesh. Those are the largest rivers on the planet, since a spilling basin
drains a catchment far larger than any single hillslope network.

One caveat carries through from the hydrography: most of the lake area drawn is
in basins pinned at their spill, which the same water balance says should have
carved their outlets by now. The lakes are real under this terrain; the terrain
is the part that has not relaxed.

If `surface_water.nc` has not been built, the map is drawn without standing
water and says so rather than failing.

The poles are forested rather than icy. That is not a rendering fault: at 32
degrees obliquity the polar summer runs well above freezing, and both polar caps
are land. The obliquity is the load-bearing number and it is a decision; the
summer temperature is a result and lives in `world_state.json`.

## The graticule

Every 30 degrees, with the equator and the prime meridian drawn heavier.
`config/planet.yaml` declares which meridian that is and which way longitude
runs, and `docs/src/reference/config-rationale.md` argues both; the heavy pass
here is on multiples of 180 degrees, so it draws the equator, the prime meridian
and the antimeridian, on the export's own longitudes rather than on any model
output's labels. It is not projected as geometry; each pixel measures its own distance to the nearest
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
