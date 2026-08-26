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

# 3. The export ships two answers about its own cell areas

Every export carries both `grid/gauss_weights.bin` and `grid/grid_cell_area.bin`
and they describe the same grid. They disagree.

`gauss_weights.bin` reproduces the Gauss-Legendre construction to 1.4e-15 at
T21, 1.3e-15 at T42 and 3.1e-15 at T85 -- the same expression evaluated twice.
`grid_cell_area.bin`, whose manifest entry calls itself "exact grid geometry" and
"the correct weight for area-averaging any gridded field", differs from those
weights by 22.0 per cent in the polar row at every rung tested, and the
disagreement does not shrink with resolution: 21.99 per cent at T21, 22.10 at T42,
22.13 at T85. In the interior it shrinks as expected of an O(1/N) construction --
0.156 per cent at T21, 0.040 at T42, 0.010 at T85 -- so the two files agree in the
limit and never at the poles.

The crossing takes the weights, because the model's global budget is a
quadrature and the weights are what it is taken over. What ELSE in the tree
area-weights with `grid_cell_area` is a separate question and is world-yolx.

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
