# The water table solve: what converged, and what the exponential cost

This is worldbuilding. Vesper is a fictional super-Earth and this document is
about a numerical solver for a groundwater model of it.

Measured 2026-08-20 on `precarve-craton`, forced by the BOOTSTRAP climatology.
This build has no `baseline_climatology` and `surface_water.nc` was forced the
same way; by this project's vocabulary a bootstrap run's numbers are not the
baseline, so every FIGURE here is provisional in that sense. The convergence
findings are not, because they are properties of the discretised problem.

`hydrography/scripts/groundwater.py` is the code.

## The model that works: `T = K D`

Gleeson conductivity over a constant aquifer thickness. Transmissivity does not
depend on the head, so the matrix is fixed and the only thing left to iterate is
the active set:

    sum_j Trans_ij (h_j - h_i) + R_i A_i = S_i,   h_i <= z_i,   S_i >= 0

That is a box-constrained linear complementarity problem with a symmetric
positive definite matrix. It converges in **20 passes**, with a water balance
residual on free cells of 2.5e-18 and **closure exact at 1.7e-16** against a
declared 1e-10. No relaxation, no trust region, no damping.

`K` and `D` come from one source at one scale: Gleeson et al. (2011) put their
permeabilities at 5-100 km and say the lithology maps carrying them "represent
the shallow subsurface (on the order of 100 m)". A region is about 15 km.

## Why the exponential had to go, which is a resolution argument

Fan et al. (2007) make conductivity decay exponentially with depth at an
e-folding length `f` set by terrain slope. Two solvers were built on it and
neither converged, but the reason to drop it is upstream of both:

**`T = A exp(h/f)` has no cell-mean value at this resolution.** It is convex, so
the mean of `T` over a water table varying within a cell by `sigma` carries a
factor `exp(sigma^2 / 2 f^2)`. A region is 15 km across and holds sub-grid relief
of order 100 m; Fan's curve takes `f` down to 0.95 m on steep bedrock. That
factor is `exp(5000)`. It is not a correction to apply.

An effective `f` recovered from neighbour-elevation spread is not a way out:
that measures relief BETWEEN cells, the mesh carries nothing about relief WITHIN
one, and `f ~ sigma` is a different model wearing the same symbol.

Fan's `f` is a hillslope closure fitted over 1.25 km cells. Mixing it into a
15 km regional permeability was the original error.

The two failures are kept because they are what a future attempt would
otherwise repeat:

**Picard on the transmissivity limit-cycles.** `T` moves by a factor of e per
e-folding length, so the map is not a contraction. It needed a trust region,
release hysteresis, orphan anchoring and permanent anchors merely to stop
diverging, and then held a water balance residual flat at 1.946e-04 from pass
100 to pass 599 while the head step fell thirtyfold. **That flatness is why the
convergence bar is the residual and not the head step**: the damping that killed
the cycle would have satisfied a head-step bar on its own.

**Kirchhoff overflows.** `Phi = int T dh = K0 f^2 exp((h-z)/f)` linearises the
problem exactly, and the topography is inside the potential, so the transform
exponentiates elevation over `f`. `abs(z)/f` reaches 4588 here and overflows
float64 on 3.29% of land; even between ADJACENT cells the difference overflows on
0.62% of faces. It removes the stiffness and replaces it with a dynamic range
double precision cannot carry.

## Four defects the checks caught, and what each one was

Recorded because each was invisible to every check but one.

**Closure caught invented water.** Negative seepage on a pinned cell was clipped
to zero on the way out, which does not discard water, it creates it: 740 m3/s of
it. The free-cell residual read 3e-18 throughout. The cause was that the loop
broke immediately after a solve that had moved the head, so complementarity --
the requirement that no pinned cell is being drawn below the surface it is
pinned to -- was never re-tested on the head actually returned. Testing both
halves fixed closure to 1.7e-16.

**A count is not a mass.** Feasibility was first tested as "no pinned cell has a
negative balance". Four cells out of 2.5 million held that above zero forever at
balances of order 1e-18 m3/s, on a planet recharging at 1.3e6. The bar now
measures the mass such cells would invent, over total recharge, which is the
quantity that matters and the one round-off cannot defeat.

**A leak bar must be tighter than the closure bar it feeds.** At 1e-8 the solve
reported convergence with a leak of 2.53e-09 and then failed closure at
2.47e-09: the same number, because the clipped mass IS the closure error.

**Permanent anchors cannot be corrected.** An orphan block -- free cells with no
conducting face to anything else -- is a pure Neumann problem with non-negative
recharge and no solution, so one of its cells must be pinned. Making that pin
permanent made a bad choice unfixable. Anchors are now releasable, and a
released candidate is never re-chosen, which terminates because each release
strikes one candidate off a set that must contain a feasible member.

**Except where it must not.** The last four cells were anchors of blocks whose
`supply` is exactly zero. An isolated block taking no recharge holds no
groundwater: there is nothing to dispose of, no seepage and no water table. It
is dry and belongs out of the network. Anchoring a cell there is not merely
unnecessary -- pinning it at the surface inside a block whose terrain is not
flat drives a flux between the block's own cells, and the pinned cell absorbs
the imbalance as exactly the negative seepage that was being clipped.

## What still misses: the uniqueness identity

The matrix is symmetric positive definite, so the complementarity problem has
exactly ONE solution and any two active-set trajectories must reach it. Solving
again from every cell free rather than every cell pinned is therefore an
identity, not a comparison, and it can fail.

**It fails.** Criterion declared before the run at 1e-9 relative; measured
**2.78e-03, a maximum head difference of 12.76 m**. Both trajectories converge,
both close to 1.7e-16, and they do not agree.

The cause is diagnosed and it is the dry-block rule above. Whether a block is
ORPHANED depends on the current active set, so which cells are found to be dry
depends on the path taken to get there: the forward run marked 14, the reverse
run 80. The two trajectories therefore solve slightly different problems, and
the identity is measuring exactly that.

**So the head field is not certified.** It converges, it conserves mass to
round-off, and it is not demonstrably the unique solution of the problem it
claims to solve. The fix is to make the dry set path-independent -- determined
from the terrain and the recharge alone, before the active-set iteration starts,
rather than discovered during it -- and that is GW-12.

## The GW-4 measurement, and its standing

Reported because it is what the term was built to find, and labelled because it
rests on the two provisos above: a bootstrap forcing, and a head field that
fails its own uniqueness check.

Net groundwater exchange redistributes **0.05% of total land recharge** between
basins, with 1,836 basins gaining and 1,496 losing. The median basin shifts by
0.07% of its own recharge, but **587 basins shift by more than 10%**, which is
the distribution Fan (2019) predicts: the term is negligible in the aggregate
and concentrated in a minority of catchments.

Against `carve_verdict.py`'s own criterion, the count that would carve rises
from 1,798 to 1,881. **All 83 flips run one way, to carve, and none the other.**
That direction is worth more than the count: groundwater import raises the water
a basin actually receives above what its surface catchment delivers, and a basin
that receives more overflows more readily. Nothing here is wired into
`carve_verdict.py`; whether the term enters the criterion is a loop A decision.
