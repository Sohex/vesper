# Placing the mesh on a grid the exporter does not write

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean cells, mesh regions, longitudes and
coastlines name modelled quantities on a simulated grid, not observations of
anything.

Measured 2026-09-07 on `canonical-10m-carve2`, terrain hash `f496ae9f`. The
driver is `ocean/scripts/build_ocean_grid.py` and the door is
`lib/gridding.py:spec_cells`. Every number below is in
`ocean/data/<build>/ocean_grid.json` with its own tolerance beside it.

**The headline.** There are now TWO mesh-to-grid doors and they are not
interchangeable. `region_cells` bins onto a grid the export ships and takes the
rows from that grid's own axis; `spec_cells` bins against a constructed grid's
cell BOUNDARIES, in the sine of latitude, and is the only door available for a
grid no exporter writes. On the atmosphere's own Gaussian grid the two disagree
about 37554 of ten million regions, and the disagreement is entirely in the
rows: zero columns move. The placement's rule-3 argument holds, and it holds on
a condition that has to be honoured rather than assumed.

---

# 1. Why a second door exists

GOLDSTEIN's rows and columns exist only as a constructor in `lib/gridding.py`.
There is no manifest to read an axis off, so a mesh region can be placed only
against the `GridSpec`'s cell boundaries. Rows are binned in the SINE of
latitude because the ocean's equal-area rows are uniform in the sine and
binning in the angle would put a region in the wrong row; `lib/gridding.py`'s
selftest carries the control that says so, an angle-binned restatement moving a
large share of the sphere. Either row order is accepted, because this project's
grids run north to south and GOLDSTEIN's `j` runs south to north.

# 2. The two doors disagree, and it is the rows

`region_cells` bins a region to the nearest Gaussian NODE. A node is a
quadrature abscissa and not the centre of its cell, so binning to the nearest
node and binning between the quadrature edges put a band of mesh either side of
every row boundary in different rows. Measured on the T21 grid:

| what | measured |
| --- | --- |
| regions placed differently | 37554 of 10000005 |
| mesh area placed differently | 0.376% |
| columns placed differently | 0 |

The zero is the informative entry. Both doors take the column from the same
expression, so the columns are an identity and the whole disagreement is the row
convention. If a column had ever moved, that would be rule 3 rather than this.

**The door for a grid the export ships is `region_cells`.** The export emitted
its gridded fields by that binning, so a reduction that uses the other one on
the same grid is describing a different partition from the one the artifact
beside it carries.

# 3. The placement, on the ocean grid

Ten million regions onto 36 x 36 equal-area GOLDSTEIN cells:

| what | measured | bar |
| --- | --- | --- |
| regions outside the boundaries of the cell they landed in | 0 | 0 |
| mesh area lost | 0.0 exactly | 1e-12 |
| ocean cells holding no mesh | 0 | 0 |
| mesh area in a cell against the cell's own area, worst | 0.41% | 2% |

Regions per ocean cell run 7691 to 7743 around a median of 7716, which is what a
relaxed mesh on an equal-area grid should give and is the resolution ratio
OCN-11's reductions have to work with.

The area residue is not noise to be tightened away and does not shrink like the
count. A region belongs to one cell whole, so what is left over is the net area
of the regions straddling the cell's perimeter: about `4*sqrt(N/C)` regions on a
perimeter out of `N/C` in the cell, so an imbalance of order `sqrt(C/N)` per
cell and a few times that on the worst of `C`. At these numbers that predicts a
couple of parts in a thousand typical and under half a per cent at the tail, and
the measurement is 0.41%. The bar was set from that argument before the run.

# 4. Rule 3, checked rather than inherited

`spec_cells` compares a mesh coordinate against edges built from the ocean
model's own setup. That is not the crossing rule 3 forbids -- matching one
grid's longitude LABEL against another's, where the mapping is the index -- but
it is not exempt by assertion either, because a mesh and a grid share no index
and the placement therefore HAS to go through coordinates. What makes it sound
is that both numbers are in one frame, and three checks say so:

| what it establishes | measured | bar |
| --- | --- | --- |
| `lon` is a coordinate, not a label: it reproduces from the mesh's own y-up cartesian `atan2(x, z)` | 7.6e-06 degrees | 1e-03 |
| the fold into the ocean's window is a WHOLE number of turns, a renumbering rather than a rotation | 0 | 1e-09 |
| `phi0` is a parameter of the ocean grid, not a naming of the atmosphere's: rotating it one column moves every region one column | 0 regions wrong | 0 |

The first is what separates a coordinate from a label: a label read off a file
reproduces from nothing. The second is the same shift
`lib/remap.py:_periodic_overlap` applies between the same two constructors,
which is why this placement is exactly as much or as little of a longitude
crossing as the `Crossing` operator already in the tree. The third is
`analysis/ocean_remap.py`'s origin check in the mesh's form.

**THE CONDITION.** The frame the comparison lives in is the generator's, which
is what defines this world's prime meridian, and the argument holds only while
the ocean grid's own topography is written THROUGH this placement. A `.k1`
produced any other way -- muffingen's route included -- carries a longitude frame
of its own, and the comparison is then between two frames rather than inside
one. That is a constraint on OCN-11 and not a property of this code.

# 5. What this does NOT establish

Nothing here declares a cell wet. The placement is where a region IS, not what
it is: wet area and volume, partial coasts, the per-cell depth distribution,
shelves and sills, connectivity and routed river mouths are all OCN-11's and none
of them is answered by an index. The GOLDSTEIN geometry remains a candidate
comparison support.

The 0.376% both doors disagree about on the atmosphere's grid is a measurement
of the two conventions and not an error attributed to either. It is recorded so
that the two doors are chosen rather than reached for.
