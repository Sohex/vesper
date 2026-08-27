# What the periodic lake balance was getting wrong, and how it is integrated now

Vesper's closed basins, `hydrography/scripts/lake_balance.py`. Measured on
2026-08-27 against `canonical-10m-base` (8,772 basins, 8,739 with an
impoundment) forced by `bootstrap_regular_climatology.nc`, which is the first
time `surface_water.py` had been run on this build.

## What was found

`solve_periodic` refused to run at all, raising

    bin 2 did not integrate in 256 sub-steps; the step limit or the
    hypsometric curve is degenerate

and the message offered two causes. Neither was it.

**The stability bound was not what exhausted the sub-steps.** The explicit
scheme carried two limits per sub-step: an accuracy limit, no basin moving more
than 2% of its capacity, and a stability limit against the local relaxation rate
`D dA/dV`. The hypsometric slope over the live set does run a long tail --
median 704, p99 7,424 and a maximum of 133,754 km2 per km3 on basin 7522, which
holds 0.113 km3 under 159.5 km2 of surface for a mean depth of 0.7 m -- but
under this forcing the stability limit asks for at most 9 sub-steps in a bin,
and it bound 68 of the 3,985 bounded steps taken in the first cycle. The
accuracy limit bound the other 3,917.

**The cause was a small reservoir on a large river.** The accuracy limit is
`0.02 * capacity / |rate|`, and `rate` is the water passing through the basin,
not the rate its storage changes at. A basin pinned at its spill does not move
at all: everything arriving leaves over the sill. Basin 8288 holds 29.8 km3 and
takes 1,381 km3/yr from upstream once the overflow cascade has loaded, so its
limit is a ten-thousandth of a bin and every one of those sub-steps recomputes
that it is still full. The requirement is not a property of the terrain alone:
it GROWS with the cascade, cycle by cycle, as basins upstream fill and begin to
spill. The worst basin needed 50 sub-steps on the first cycle, 95 on the third
and 192 on the fifth, and crossed 256 after that. No constant bounds it, which
is why raising the limit would not have been a fix.

**A basin holding no water was evaporating a mesh cell of lake it did not
have.** Every hypsometric curve in `basins.nc` is built from mesh cells, so its
first level already carries a whole cell of surface: over the live set the area
at zero volume runs from 1.8 to 213.2 km2, median 63.7, against a mean mesh cell
of 73.5 km2. A basin whose supply cannot fill that first cell was being charged
open-water evaporation over the curve's own first area, water that is not there.
It affected 3,425 basins and destroyed 79.8 km3 per cycle against 76,502 km3/yr
passing through the catalogue, a tenth of a percent. It also suppressed the
seasonal signature of exactly the basins the product exists to describe: a playa
that empties in the dry bins was pinned at one mesh cell instead of shrinking,
so its swing read as smaller than it is.

`solve()` carries none of this. It has no time integration, so neither bound
applies to it; its cascade iteration converges in 15 passes against a limit of
200, which is the longest chain in the spill graph (14) plus one, and its answer
is bit-identical from a 1e-3 cascade tolerance to a 1e-12 one. It already read
`supply / demand` at the bottom of the curve with no mesh-cell floor, so the two
solves disagreed there by up to one mesh cell until the fix below.

## What is done about it

**The bin is integrated in closed form rather than stepped.** The forcing does
not change inside a time bin and a hypsometric curve is piecewise linear, so

    dV/dt = S - D A(V)

is a LINEAR ordinary differential equation on every segment of the curve, and
its solution is the exponential. `_integrate_bin` walks the segments a basin
crosses and solves each one exactly, computing the time to cross a segment
rather than choosing a step to take.

This is not a shortcut with a threshold, and that is the argument for it over
the two obvious alternatives. Treating a short-residence basin as instantaneous
needs a residence-time threshold that has to be derived and swept, and is only
right in the limit; integrating the stiff subset implicitly needs a criterion
for which subset and is first-order accurate everywhere else. The closed form
needs neither, is exact at every residence time, and REMOVES two constants
instead of adding one: `PERIODIC_MAX_STEP_FRACTION` at 0.02 and
`PERIODIC_STABILITY_MARGIN` at 0.5 were both numbers chosen to be small rather
than derived from anything, and there is nothing left to choose.

The work is bounded by the curve rather than by a limit. `A` is non-decreasing
in `V`, so the right-hand side is monotone in `V`, so within a bin the storage
moves in one direction only and crosses each of the curve's 128 levels at most
once. The iteration bound is the curve's own length, and reaching it means a
curve whose volume is not ascending or a step that made no progress -- both
defects, neither a bin that needed a finer integration.

Two things the closed form gets right that the stepped one did not, beyond not
raising. A basin at its spill is charged the open-water demand at the SPILL
area for the whole time it is spilling, where the stepped scheme charged the
demand at the level it started each step from. And the cycle means are now the
integrator's own exact time integrals over each bin rather than the bin-end
samples averaged, which is what a water balance needs and what a rectangle rule
on a basin that turns over inside a bin cannot give.

**A basin with no storage gets the only lake the balance allows.**
`_sub_mesh_area` returns `S / D`, the area whose evaporative demand exactly
consumes the supply, clipped to the curve's first cell. It is below one cell by
the same condition that put the basin there, and it is the same area `solve()`
returns at the bottom of the curve, so the two solves now agree at the dry end.
Nothing is fitted: with no storage to draw down, an area that evaporates more
than arrives is not a state the basin can be in.

## The evidence, against tolerances fixed before it was measured

The tolerances are declared in `lake_balance.py` as `EXACTNESS_RELATIVE`,
`BALANCE_RELATIVE` and `REFERENCE_TOLERANCE`, and each is set by what double
precision can carry over the work involved rather than by what came out.
`--selftest` prints the achieved margin on every one of them whether it passes
or not.

| Check | Declared | Achieved |
| --- | --- | --- |
| the closed form is the exponential, on one segment | 1e-13 relative | 1.3e-16 |
| a spilling basin passes its throughflow in one step, 22x its capacity in the bin | 1e-13 relative | exact |
| agreement with Radau at rtol 1e-10, over 72 stiff cases | 1e-6 of capacity | 1.3e-9 |
| every bin's water balance closes on itself | 1e-10 relative | 2.9e-18 |
| the year's water closes, per basin | 1e-10 relative | 1.1e-17 |
| the year's water closes, over the set | 1e-10 relative | 2.2e-15 |
| a basin below one mesh cell conserves water and matches `solve()` | 1e-10 and 1e-12 | exact |

The reference integrator is four orders tighter on its own tolerance than the
agreement asked of it, so what that row measures is this solver and not the
reference. The last row fails at 0.031 and 299 km2 against the defect it was
written for, and the bin-balance rows fail at 0.0178 relative on the real set
without `_sub_mesh_area`, which is how that defect was found: the balance check
was declared first and refused to pass.

On `canonical-10m-base` under the bootstrap climatology, all 8,739 live basins
close their year in 33 cycles, the worst per-bin residual is 1.09e-16 and the
worst annual residual 4.93e-16 relative, and the set closes to 2.2e-15 against
76,502 km3/yr passing through it. `surface_water.py` writes that block into
`surface_water_report.json` and refuses the result if it misses.

## What this does not settle

The cascade's timing inside a bin. An upstream basin's overflow is delivered to
the basin below it as a rate held constant across the bin, recomputed between
cycles, so a basin that begins spilling part way through a bin delivers as
though it had spilled throughout. That is unchanged, it is separate from how a
bin is integrated, and it is the same approximation `solve()`'s cascade makes.
