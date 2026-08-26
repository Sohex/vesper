# The atmosphere-to-ocean crossing: exact weights, and the integral that proves it

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean cells, heat flux, wet masks and coastlines
name modelled quantities on a simulated grid, not observations of anything.

Measured 2026-08-25. This is world-gb0t, the implementation half of section 5
of `notes/audits/cgenie-parallelism-and-coupling-support.md`. `vendor/cgenie`
was not edited and nothing here adopts it; what exists now is an operator and
its acceptance test, and the ocean component that would consume them is still
the candidate under OCN-3.

`analysis/ocean_remap.py` is the driver, `lib/gridding.py` holds the grid
constructors and `lib/remap.py` the operator. Every number below is in
`analysis/ocean_remap.json` with its own tolerance beside it.

**The headline.** The crossing is exact rather than approximate, the acceptance
test has an analytic answer, and the tolerance it is held to is three orders
tighter than the precedent's. The residual it actually returns is another three
orders below that, so the operator is limited by float64 and not by the method.
The part that is not settled by construction is the coastline, and that is
reported as a quantity rather than absorbed.

---

# 1. Why the geometry is cheap, and why it is not averaging

ExoPlaSim's grid is uniform in longitude with Gaussian latitudes. GOLDSTEIN's
`igrid = 0` is uniform in longitude and uniform in the SINE of latitude, so
`initialise_goldstein.F`'s `asurf(j) = rsc*rsc*ds(j)*dphi` gives every ocean cell
exactly the same area; the selftest confirms the spread across a 36 x 36 ocean's
cells is at the last bit of a double. Neither grid's row boundaries depend on
longitude and neither grid's column boundaries depend on latitude, so an overlap
area factorises into a longitude overlap times a sine-of-latitude overlap: two
one-dimensional matrices, no polygon intersection, and the planetary radius
cancels because every weight is a ratio of areas.

A mask does NOT factorise. It is genuinely two-dimensional, so the
factorisation is used for the sparsity pattern and the masked matrix is
assembled sparse over it. At a T85 atmosphere against a 72 x 72 ocean that
pattern holds 63,680 nonzeros, which is why the exact operator is cheaper than
the approximate ones it replaces rather than more expensive.

**What crosses is the flux, not the state the flux is computed from.** A surface
flux is nonlinear in the surface state, so the average of the flux is not the
flux of the average. In the architecture OCN-10 selects the flux has already
been evaluated at the atmosphere's own support, and averaging it is a
conservative linear operation. The Jensen gap `lib/gridding.py:cell_expectation`
exists for is therefore already spent on the fine grid, and what upscaling costs
is exactly one thing: the subgrid variance of the flux inside an ocean cell.

# 2. The acceptance test, and why it can fail

Two identities, and only the second catches the failure that matters.

**The crossing preserves the source's own integral.** For a field carrying a
budget, `sum(F_dst*A_dst)` must equal `sum(F_src*A_src)`. The bar is 1e-12
relative, set from the float64 round-off floor -- 2.2e-16 per operation, about
3.4e-15 for a pairwise sum over a T85 grid's 32,768 cells -- and it is three
orders tighter than the 1e-9 relative that `references/esmf/` holds itself to on
analytic fields. Measured at T21 against 36 x 36 the residual is 4.6e-18, and at
T85 against 72 x 72 it is 5.8e-18.

**The source integral is one whose answer is known in advance.** That second
identity exists because the first one cannot see the failure section 5c of the
coupling audit predicted. A crossing built on midpoint-between-nodes row
boundaries conserves perfectly well: it is a self-consistent partition of the
sphere. It just conserves against a grid the model does not use, and it does
that quietly. So the field remapped is `3*sin(lat)^2 - 1` evaluated at the
Gaussian NODES, whose integral over the sphere is exactly zero and which
Gauss-Legendre quadrature returns as exactly zero at any node count above one.
The source quadrature comes back at 2.6e-15 of the field's own scale and the
crossing returns 2.7e-15, both against the 1e-12 bar.

**The negative control fires.** The same field over midpoint row boundaries
returns 4.0e-3, seven orders above the bar and three orders above the 1e-9 floor
fixed in advance for it. So the bar discriminates, which is what makes the two
checks above tests rather than decorations.

Three more identities with right answers rather than plausible ones: both
one-dimensional overlaps are partitions of unity in both directions, which a
periodic wrap that folds instead of wrapping fails even though a crossing built
on it would still conserve; a grid crossed with ITSELF is the identity under all
three semantics, to 2.1e-16; and rotating the ocean's `phi0` by exactly one ocean
column rolls the answer by one column and changes nothing else, to 2.7e-17.
Seventeen checks in all, and the tolerances were fixed before any of them ran.

# 3. The export shipped two answers about its own cell areas, and which one is right

Every export carries both `grid/gauss_weights.bin` and `grid/grid_cell_area.bin`
and they describe the same grid. They disagreed, and world-yolx settled which of
them was wrong.

`gauss_weights.bin` reproduces the Gauss-Legendre construction to 1.4e-15 at
T21, 1.3e-15 at T42 and 3.1e-15 at T85 -- the same expression evaluated twice.
`grid_cell_area.bin` differed from those weights by 22.0 per cent in the polar
row at every rung tested, and the disagreement did not shrink with resolution:
21.99 per cent at T21, 22.10 at T42, 22.13 at T85. In the interior it shrank as
expected of an O(1/N) construction -- 0.156 per cent at T21, 0.040 at T42, 0.010
at T85 -- so the two agreed in the limit and never at the poles.

## 3a. It was a real partition of the sphere, answering a question nobody asked

`vendor/orogen/js/geometry.js:latitudeEdges` builds row boundaries midway in the
SINE of latitude between adjacent row centres, with the poles clamped. That is
the NEAREST-ROW partition, and it is exactly right for what it was written for:
`rowForLatitude` uses it to decide which grid row a mesh point falls in, and a
point should go to the row whose centre is closest to it. `gridCellArea` then
reused those edges as if they were cell boundaries.

So the file was not arithmetic nonsense. It was the exact area of the region of
the sphere that BINS into each cell, it closed to 4 pi R^2 to the last bit, and
that is why the disagreement survived: a wrong partition of the sphere passes
every conservation check a right one passes. Section 2's negative control is the
same construction and it misses the known-integral bar by seven orders.

**A Gaussian row is a quadrature abscissa and not a cell centre**, so there is no
geometric midpoint to appeal to and the partition is chosen rather than derived.
Only one choice makes a global mean of a model field equal the model's own global
mean, and it is the quadrature's: the weights ARE the sine extents of the
intervals a spectral transform integrates over, they sum to 2, and laid end to end
from the north pole they tile the sphere. So the label was the thing that was
wrong, and the file has been corrected to carry the partition the label promised
rather than the label corrected to describe the binning partition -- because
nothing in the tree wants a binning area, and four documents were recruiting
readers to use it as a weight.

The same defect had a second, independent instance: the uniform
`grid-512x256` export, whose rows ARE cell centres and whose cell boundary is
therefore unambiguously midway in LATITUDE, also took its bands from
`latitudeEdges` and was 25.0 per cent wide in its polar row, converging to
exactly 25 and never to zero.

## 3b. What the error was worth, on fields the tree already has

A 22 per cent error in one row is only worth what that row's area share is, and
on a Gaussian grid the polar row is the smallest there is. Measured on
`precarve-craton`'s own `planet.nc`, the shipped field against the quadrature:

| rung | polar rows' share of the sphere | of land | endorheic land-area fraction | global land fraction | land-mean elevation |
| --- | --- | --- | --- | --- | --- |
| T21 | 0.702 per cent | 1.438 per cent | +0.143 per cent | +0.163 per cent | -0.104 per cent |
| T42 | 0.178 per cent | 0.412 per cent | +0.010 per cent | +0.052 per cent | -0.040 per cent |
| T85 | 0.045 per cent | 0.104 per cent | -0.009 per cent | +0.013 per cent | -0.010 per cent |

Both polar rows are entirely land on this terrain, so the error is not masked
away; it is diluted, by the area share and nothing else. **That dilution is a
property of the statistic, not of the field**, so a diagnostic restricted to the
polar rows takes the full 22 per cent, and an ABSOLUTE per-cell quantity takes it
undamped in the cells it lands on.

Which is the shape of the second finding. `pedology/scripts/weathering_fluxes.py`
carried its own cell-area construction, midway in LATITUDE rather than in its
sine -- a third partition again -- 5.79 per cent wide in the polar row and again
not converging. It turns a flux per litre into an absolute silica and CO2 flux
per cell per year, and `phosphorus_budget.py` imports it for phosphorus release,
so there the polar cells were simply 5.8 per cent high with no ratio to damp it.
Both now take `lib/gridding.py:gaussian_grid`.

## 3c. What changed, and what it makes worthless

`vendor/orogen/js/geometry.js` gains `cellSinEdges`, which returns the cell
boundaries -- the quadrature intervals on a Gaussian grid, midway in latitude on
an equally spaced one -- and `gridCellArea` takes them. `latitudeEdges` is
untouched, because the binning it serves was never wrong; its docstring now says
which of the two questions it answers. Three tests in
`vendor/orogen/tools/test-integration.mjs` hold the new construction to an
identity against the weights and carry the binning partition as a negative
control that must MISS, so the check discriminates rather than decorates.

`lib/gridding.py:export_grid` constructs an export's cell boundaries from the
manifest's grid type and shape and checks them against the axis the export
ships. It never reads `grid_cell_area.bin`, which is what makes it right on both
sides of this change.

**Per rule 7, this makes a FIELD worthless and no build.** Terrain is untouched,
so no export is superseded; but every export on disk carries the nearest-row
partition in `grid/grid_cell_area.bin` and in `planet.nc`'s `grid_cell_area`,
and those are wrong until a build is regenerated. Anything derived from them is
worthless in the same narrow way: `pedology`'s weathering, phosphorus and
thermostat products, in their absolute per-cell numbers and in their polar rows.
This document's own section 3 measurements are measurements OF the superseded
field and stay true of what is on disk.

# 4. The coastline is the part construction does not settle

The ocean's wet mask comes from the Orogen mesh through OCN-11 and the
atmosphere's comes from the same mesh at a different support. They cannot agree
at the coast: a cell that is 40 per cent land at the atmosphere's support sits
inside an ocean cell that is either wet or dry.

Two rules, and both are stated rather than absorbed.

**A valid destination cell with no valid source over it is an ERROR.** It is a
mask disagreement, not sparse coverage. `lib/gridding.py:land_weighted`
backfills a grid cell that holds no mesh region from the nearest region, and
that substitution is defensible where the grid is finer than the mesh; here it
would put an invented flux into an ocean cell that integrates it for the life of
the run.

**A valid source cell with no valid destination under it goes to the nearest
valid destination cell**, by great-circle distance between cell centres, in the
CONSERVING matrix only. An intensive state in an orphaned source cell has
nothing to conserve, so it is dropped and counted. This is the single-neighbour
form of what `references/climber-x/src/geo/coast_cells.f90` does for river
discharge with an expanding stencil; the difference is that this one has to
carry heat and salt as well as water.

Measured against the export's own `surface_class`, with the ocean's wet mask
standing in as the atmosphere's ocean fraction remapped and cut at a half:

| atmosphere | ocean | orphaned source cells | source area moved | furthest move |
| --- | --- | --- | --- | --- |
| T21 | 36 x 36 | 12 | 0.567 per cent of the sphere | 12.4 degrees of arc |
| T42 | 36 x 36 | 168 | 2.086 per cent | 10.3 degrees |
| T85 | 72 x 72 | 465 | 1.344 per cent | 9.8 degrees |

Conservation holds exactly across every one of those: 4.1e-18 relative at T21,
because the rule names where the flux goes rather than letting it leak. **This
is the part of the crossing a finer ocean genuinely improves**, and it is the
argument for resolution that survives once the conservation argument is settled
by construction. It is also a demonstration and not the coupling's mask: the
real ocean wet mask arrives with OCN-11, and the threshold above is a stand-in
chosen to produce a disagreement of the right shape, not a declared cut.

# 5. What this did NOT establish

- **Nothing here is a coupling.** No ocean component reads these weights, no
  contract names which field takes which semantics, and OCN-10 is still open.
  What exists is the operator that contract will be written against.
- **The other direction is not built.** Sea surface temperature, the ice state
  and the exported surface velocity are intensive states, so the crossing back
  is an interpolation with nothing to conserve. `lib/remap.py` will carry an
  intensive field either way, but what it does coarse-to-fine is
  piecewise-constant and invents no gradients; whether that is what the next
  climate run should be forced with is section 5e's question and not settled
  here.
- **Nothing was measured in seconds.** The host is contended and the operator's
  cost is quoted as a nonzero count, which is a property of the two grids.
