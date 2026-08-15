# maps/

Illustrative colour maps of Vesper in five projections. These are pictures, not
analysis products: they exist to show what the world looks like.

```bash
python maps/build_basemap.py       # equirectangular base, ~25 s (lookup cached)
python maps/render_projections.py  # the five projections, ~4 min at width 3000
```

| File | Projection | Property |
| --- | --- | --- |
| `vesper_winkel_tripel.png` | Winkel tripel | compromise; the usual world map |
| `vesper_equal_earth.png` | Equal Earth | equal-area pseudocylindrical |
| `vesper_lambert_azimuthal.png` | Lambert azimuthal | equal-area, east/west hemispheres |
| `vesper_lambert_azimuthal_polar.png` | Lambert azimuthal | equal-area, north/south hemispheres |
| `vesper_waterman_butterfly.png` | butterfly | octahedral, interrupted |
| `vesper_dymaxion.png` | Dymaxion | icosahedral, interrupted |

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

## Orientation is derived, not assumed

Every projection here has a free orientation, and each one is searched against
the land mask rather than left at zero:

- central meridian, so the antimeridian crosses the least land;
- hemisphere boundary for the equatorial azimuthal pair, likewise;
- the butterfly's equatorial vertices, scored on its two cut meridians;
- the icosahedron's attitude, scored on all 30 edges, since every cut the net
  can make runs along one of them.

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

## Approximations worth knowing

Both polyhedral maps are **gnomonic on their polyhedron**, which is exact for
the face geometry: gnomonic sends great circles to straight lines, so a
spherical triangle maps onto the planar triangle its vertices span, and the
inverse is barycentric. That is the standard construction behind the Fuller
projection, but it is not Fuller's own face mapping, and the butterfly is a
Cahill-style octahedral net rather than Waterman's truncated-octahedron
construction with its clipped corners. Shapes within a face differ slightly from
the published projections; the arrangement and the interruptions do not.
