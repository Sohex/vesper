# The water table solve does not converge, and why each scheme fails

This is worldbuilding. Vesper is a fictional super-Earth and this document is
about a numerical solver for a groundwater model of it.

Measured 2026-08-20 on `precarve-craton`, forced by the BOOTSTRAP climatology.
This build has no `baseline_climatology` and `surface_water.nc` was forced the
same way; by this project's vocabulary a bootstrap run's numbers are not the
baseline, so every figure here is provisional in that sense. The convergence
findings are not, because they are properties of the discretised problem rather
than of the forcing.

`hydrography/scripts/groundwater.py` is the code. Two schemes were tried. Both
fail, for different and separable reasons, and neither failure is a coding
detail that a little more care would remove.

## The problem

Steady-state, vertically integrated, unconfined:

    sum_j Trans_ij (h_j - h_i) + R_i A_i = S_i,    h_i <= z_i,   S_i >= 0

with Fan et al. (2007) transmissivity `T(d) = K0 f exp(-d/f)` for water table
depth `d` below the surface. The complementarity is not the hard part. The
transmissivity is.

## Scheme 1: Picard on the transmissivity. Limit-cycles.

Evaluate `T` at a head, solve the linear system, repeat. `T` moves by a factor
of e per e-folding length `f`, and `f` runs from 120 m on flat regolith down to
0.95 m on steep bedrock, so the map is nowhere near a contraction and the
iterate cycles instead of converging.

What it took to get it even to stability, each item forced by a failure:

| symptom | cause | remedy |
| --- | --- | --- |
| head oscillating by 3-4 km per pass on a planet whose land spans 5 | a step much larger than `f` leaves the regime the linearisation was taken in | trust region of one e-folding length |
| active set flipping ~2,500 cells either way forever at ~385,500 free | cells exactly on the boundary release and re-pin every pass | release hysteresis at 1% of a cell's own recharge |
| matrix exactly singular after ~30 passes | a free block connected only to free cells is pure Neumann with non-negative recharge, so it has no solution | pin its shallowest cell, and never release an anchor |
| the same block re-anchored every pass | the anchor was being freed again by the release test | anchors are permanent |

With all of that it is stable and still does not converge. Run to 600 passes:

| pass | max abs head step | p99.9 | water balance residual |
| --- | ---: | ---: | ---: |
| 100 | 70.6 m | 22.33 m | 1.946e-04 |
| 158 | 70.0 m | 22.39 m | 1.946e-04 |
| 599 | 2.5 m | 0.26 m | 1.946e-04 |

**The residual is flat to four figures over five hundred passes while the head
step falls by a factor of thirty.** That is the signature of a diminishing
relaxation freezing the iterate short of the fixed point rather than finding it,
and it is the reason the convergence bar is the water balance residual and not
the head step.

**The bar was changed mid-work, and that is worth stating plainly.** It was
first declared on the head step. A diminishing step was then added to damp the
limit cycle, at which point the head-step bar became satisfiable by the damping
alone, with no convergence at all -- and the table above shows it would have
been satisfied. It was replaced by the residual, which measures the nonlinear
problem rather than the iteration chasing it. The replacement happened after
seeing the stall, so it is not a criterion declared in advance; it is a
criterion that was invalid being withdrawn.

## Scheme 2: Kirchhoff. Overflows, and the relief is why.

The Kirchhoff potential of a nonlinear diffusion is the integral of its own
coefficient, and for this transmissivity it is closed:

    Phi = int T dh = K0 f^2 exp((h - z)/f) = T f,   with   grad Phi = T grad h

so `div(T grad h) + R = 0` becomes linear. The cap maps monotonically to
`Phi <= K0 f^2`, leaving a box-constrained linear complementarity problem, which
terminates where a fixed-point iteration cycles. That is the attraction and it
is real.

**Two things defeat it here, and the first is an error in this implementation.**

**`grad Phi = T grad h` holds only where `K0`, `f` AND `z` are all constant.**
The topography is inside the potential, not only the material properties. The
first implementation used `Phi` as the unknown under a PLAIN Laplacian, which is
valid only if `A f = K0 f^2 exp(-z/f)` is spatially uniform. It is not, and the
correct form is a variable-coefficient Laplacian in `u = exp(h/f)` with
coefficient `A f`. That was an error rather than an approximation.

**The correct form is not representable in double precision on this world.**
The transform's variable is an exponential of elevation over the e-folding
length, and `f` is a regolith depth of tens of metres while the relief is
kilometres:

| quantity | median | p90 | p99 | max |
| --- | ---: | ---: | ---: | ---: |
| `abs(z)/f` over land | 1.1 | 203.7 | 1245 | 4588 |
| `abs(dz)/f` between ADJACENT cells | 0.35 | 43.9 | 401.7 | 4470 |

`exp` of the first overflows float64 on 3.29% of land and exceeds 1e30 on
15.9%. The obvious rescaling -- keep `Phi` as the unknown and put the ratio
`exp((z_j - z_i)/f)` in the face coefficient, so only the elevation DIFFERENCE
is exponentiated -- does not rescue it: that difference overflows on 0.31% of
faces taking `f` as the face mean and 0.62% taking the smaller of the two, and
exceeds 1e30 on 7.43%, because adjacent cells 15 km apart can differ by hundreds
of metres while `f` falls to a metre. `solve_kirchhoff` now measures this and
refuses rather than running, on the stricter of the two.

The run behaved exactly as that predicts. The reported residual collapsed to
7e-19, which is `Phi` underflowed to zero over most of the domain rather than a
balance being satisfied, and the matrix was exactly singular by the third pass.

**So the transform removes the stiffness in exact arithmetic and replaces it
with a dynamic range float64 cannot carry.** That is a property of this world's
terrain against this parameterisation of `f`, not of the algebra.

## What this leaves

`GW-1` is blocked on convergence and `GW-4`'s measurement is blocked behind it.
No water table is written: `build_groundwater.py` refuses to write
`water_table.nc` from an unconverged solve, because a head field in `data/`
would be read downstream as a result and nothing in the file would say it was
not one.

What still stands, and does not depend on the solve converging:

- The reduction identity passes bitwise. At zero permeability every cell
  returns its own recharge as seepage and every basin's exchange is exactly
  zero, against `surface_water.py`.
- The mesh operator is checked against an analytic eigenvalue and its error
  characterised. `notes/mesh-geometry.md`.
- The permeability table, its class mapping and the gravity term are grounded
  and independent of the solver.

Three routes remain, and the connection between the second and third is that
both attack `f` rather than the algebra:

**Newton on the full residual.** Solves the same stiff problem rather than
removing the stiffness, but converges quadratically near the solution and does
not care that the Picard map is not a contraction. The Jacobian carries the
`dT/dh` terms Picard discards, which is exactly the information that is missing.

**Rescale the exponential per cell.** Carry `Phi` in logarithms, or factor a
reference elevation out of each cell so that only local differences are
exponentiated. The measurement above says local differences are not small
enough on their own, so this needs to be paired with a floor on `f`.

**Reconsider `f` itself.** A one-metre e-folding length on a 15 km cell is the
proximate cause of every numerical difficulty here, and it comes from reading
Fan's hillslope-fitted curve with a regional dip. That is GW-9, and it is a
modelling question rather than a numerical one: what the e-folding length should
be at this mesh resolution is not obviously what it is at Fan's 1.25 km.
