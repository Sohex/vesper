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

## The uniqueness identity, and what it took to pass

The matrix is symmetric positive definite, so the complementarity problem has
exactly ONE solution and any two active-set trajectories must reach it. Solving
again from every cell free rather than every cell pinned is therefore an
identity, not a comparison, and it can fail.

**It failed first, at a maximum head difference of 12.76 m against a declared
1e-9 relative.** Both trajectories converged and both closed to 1.7e-16, and
they did not agree.

The cause was that a block of cells could be found DRY during the iteration --
recharge-free and cut off from everything -- and whether a block looks cut off
depends on which cells happen to be free at that moment. The forward run marked
14 cells dry and the reverse run 80, so the two were solving slightly different
problems and the identity was measuring exactly that.

**The dry set is a static property and is now computed as one.** Water enters
this system in exactly two places: recharge falling on a conductive cell, and
the ocean boundary. A conductive cell that connects to neither, through any
chain of conducting faces, can never hold groundwater. A face conducts whenever
both its ends are in the network, whatever their pinned or free status, so that
connectivity is the same on every pass. It is one connected-components query on
the static graph, before the iteration starts. **112 regions** on this build.

With that, the identity passes and does not merely come inside its bar:

    max |h_from_all_free - h_from_all_pinned| = 0.000e+00 m

Bit-identical, with the same 112 dry regions found from either direction.

**Which model that was measured on.** 2026-08-20, before GW-15's
evapotranspiration sink and GW-17's local baselevels existed, so both
trajectories carried neither. It is an honest identity for the equation of that
date. The equation has since gained both terms and both default ON, so the
number above is evidence about the model it was run on and not about the model
the component now solves.

## The identity on the current equation, and the controls that reject it

Measured 2026-08-25 by

    python hydrography/scripts/groundwater.py --uniqueness-test

on a synthetic 2,304-cell case carrying every term the planet's solve carries:
an ocean boundary, GW-17's imposed baselevels, GW-15's spatially varying sink,
three orders of conductivity, a wet patch where the box constraint binds, and a
block enclosed by excluded lithology so the static dry-set query is reached. The
bar is `SCHEME_HEAD_RELATIVE`, unchanged, and a control is rejected by that same
criterion; there is no second bar, so nothing here was chosen after a result.

| Arm | max head difference | relative | verdict required | got |
| --- | --- | --- | --- | --- |
| identity, two trajectories | 0.000e+00 m | 0.000e+00 | pass | pass |
| a third, arbitrary active set | 5.296e-10 m | 1.250e-12 | pass | pass |
| re-solve without GW-15's sink | 5.301e+00 m | 1.251e-02 | rejected | rejected |
| re-solve without GW-17's baselevels | 6.244e+01 m | 1.474e-01 | rejected | rejected |
| `solve` mutated: pinned cells seeded below their surface | 4.467e+00 m | 1.054e-02 | rejected | rejected |

**The identity is bit-identical with both terms in.** It was not before this was
run; it was untested.

**The two dropped-term rows are the measurement of what the broken arm was
reporting.** Until world-60x0 the re-solve was given neither term, so it was
comparing two models, and the size of that comparison is the two rows above: a
tenth of the head range for the baselevels, a hundredth for the sink. The arm
could only ever have returned a miss of that order, or have been run in a
configuration where the omission did not bite. It was never evidence about
uniqueness.

**The mutated-solver row is why the identity is a test rather than a report.**
`pinned_seed` reverts one invariant of `solve`: a cell that starts below the
surface must start FREE, because the seepage accounting reads a pinned cell as
being at its surface. Both trajectories still converge and still close, and they
land 4.467 m apart. A check whose only recorded value is a pass is not evidence
that it can fail, and this is the arm that shows it going red.

The explicit anchoring path is NOT covered here and no control is offered on it:
a cell carrying GW-15's sink anchors itself, so with the sink on a free block
only reaches the anchoring code where `et_max` is zero. What a control on the
anchor release would take is tracked separately.

**The planet arm is a separate claim and is still unmeasured.** The command is

    python hydrography/scripts/build_groundwater.py --uniqueness-check --climatology <path>

at the same bar. It has no result on `precarve-craton-10m`: the reverse
trajectory starts with every land cell free, so it factors the whole 4.33
million-cell land set before a single cell is pinned, and this host has not had
the memory free to do that. The re-solve now takes the primary solve's own
argument list whole rather than a copy of it, so a term can no longer reach one
trajectory and not the other; the copy is how the sink and the baselevels came
to be missing from this arm in the first place.

## The catchment check: a real bug, and a bar that could not be met

Three attempts, and they are recorded together because two of them missed and
the reason matters more than the fact.

**Attempt 1 compared a flux trace on the RAW surface against
`regions.nc:terminal`, a priority flood on the FILLED one.** 73.0% of land area
against a declared 90%. MISS.

**Attempt 2 moved both sides onto the filled surface**, on the theory that the
surfaces were the disagreement. 71.6% against a declared 95%. MISS, and barely
moved, which killed the theory.

**That failure exposed a real defect.** `groundwater_receiver` read the module's
own sign convention backwards. A positive face flux is flow from `dst` INTO
`src`, so what leaves `src` is `-flux`; the function took the positive direction
as outgoing and therefore traced every cell to the neighbour it receives most
water FROM. It was following the water uphill. Fixing it moved the same
comparison from 71.6% to **79.3%**, which is how a genuine bug shows up in a
mis-specified test: real improvement, still a miss.

**And the bar could not have been met.** `terminal` comes from the priority
flood's discovery pointer, and hydrography's README says why the flood cannot
use steepest descent: a filled pit is flat, so descent would drop whole
tributaries. A flux trace IS a steepest-descent rule. The two routing rules
disagree on identical terrain by construction, so no solver reconciles them and
95% was declared for a comparison that cannot be an identity. **It stands as a
miss at 79.3%**, and it is reported as a measurement rather than promoted to a
pass.

**What CAN be exact is the trace machinery, and that is now the check.** Hand
`groundwater_receiver` and `trace_terminals` a flux field whose only outgoing
flux at each land cell is the face to that cell's own surface `receiver`, and
they must reproduce `terminal` on every land region. Nothing about groundwater
enters it. Declared exact before running, and measured:

| what | result |
| --- | --- |
| `terminal` reproduced, by land area | **100.0000%** |
| `receiver` itself recovered, by region count | 92.49% |

The second number is lower and is not a defect: where several faces carry equal
outgoing flux the tie is broken arbitrarily, and those cells still reach the
same terminal, which is what the catchment is. This is the check that caught the
sign bug above.

## GW-4: what the groundwater term moves

Forced by the BOOTSTRAP climatology. This build has no `baseline_climatology`
and `surface_water.nc` was forced the same way; this project's vocabulary is
explicit that a bootstrap run's numbers are not the baseline, so these are
provisional against a baseline that does not yet exist. They are not provisional
against the solver: the head field closes to round-off and the trace machinery
reproduces the surface catchments exactly. The uniqueness identity is a
separate claim and is carried by the section above, which says which equation
it was measured on.

Gleeson's within-class permeability spread is 1.5 to 2.5 orders of magnitude, so
the depth field is a bracket and every arm is reported. `--sigma` shifts each
hydrolithology by its OWN standard deviation, not by a shared one.

| | sigma -1 | sigma 0 | sigma +1 |
| --- | ---: | ---: | ---: |
| water table at the surface, share of land | 0.672 | 0.630 | 0.468 |
| exchange as a share of land recharge | 0.002% | 0.052% | 0.989% |
| basins shifted by over 10% of their own recharge | 423 | 587 | 1115 |
| basins that would carve, surface balance only | 1798 | 1798 | 1798 |
| basins that would carve, with groundwater | 1884 | 1881 | 1877 |
| verdicts flipped | 86 | 83 | 85 |
| of those, flipped towards HOLDING | 0 | 0 | 3 |

**The depth field is genuinely bracketed and the flip count is not.** The share
of land with the water table at the surface runs from 0.47 to 0.67 across the
bracket, and the exchange spans a factor of five hundred, while the number of
carve verdicts that flip moves only from 86 to 83 to 85.

**And the same basins flip.** Counts alone could hide three different
populations, so the sets were intersected:

| | value |
| --- | ---: |
| flip at all three arms | 74 |
| flip at exactly two | 10 |
| flip at exactly one | 12 |
| union over the three arms | 96 |
| pairwise Jaccard | 0.78 to 0.90 |

Against the operator-noise members the agreement is near-total: symmetric
difference from the central set of 1, 0 and 0 basins.

## Why a five-hundredfold change in exchange moves the count by three

Not because those basins sit near their thresholds. That was the first reading
and it is wrong.

**Every flip is a basin whose surface recharge is EXACTLY zero.** 86 of 86 at
sigma -1, 83 of 83 at the centre, 81 of 85 at sigma +1. Runoff is
`max(P - E, 0)` per cell, so a catchment where evaporative demand meets or
exceeds precipitation everywhere delivers nothing at all, and

    aridity index = (E_lake - P_lake) / runoff  ->  infinite

An infinite index never satisfies `index <= critical_aridity_index`, so such a
basin never carves under the surface balance, whatever its geometry. Switch
groundwater on and lateral seepage gives it a finite runoff depth and therefore
a finite index, which may fall below its threshold.

**So the flip is a switch at zero, not a response to a magnitude.** It turns on
whether any groundwater arrives, not how much, which is exactly why five hundred
times more of it changes the count by three. There are **605** such fully-arid
basins in the catalogue and 83 of them cross.

The global exchange fraction is also not the per-basin perturbation, and it was
wrong to read it that way: at sigma -1 the exchange is 0.002% of total land
recharge while **423 basins still shift by more than 10% of their own**. The
global figure is dominated by large wet basins in its denominator.

**What this means for how the number may be quoted.** "Groundwater changes 83
carve verdicts" overstates it. The defensible statement is that **83 of the 605
basins with no net surface runoff at all acquire a finite aridity index below
their threshold once groundwater is included**, and that this count is stable to
the permeability bracket and to the operator's truncation error because it is a
zero-crossing rather than a magnitude.

Whether "no surface runoff, therefore never carves" is the criterion semantics
intended by `carve_verdict.py` is a question for that criterion and not for this
component. It is raised, not answered, here.

## The flips that run the other way

Three basins flip towards HOLDING, all at sigma +1, and the mechanism is the
expected one: **all three EXPORT groundwater**, so the water reaching them falls
rather than rises.

| basin | Qg, m3/s | recharge, m3/s | Qg / recharge | catchment, km2 |
| ---: | ---: | ---: | ---: | ---: |
| 906 | -2.46 | 13.88 | -0.18 | 8,958 |
| 2277 | -0.01 | 0.01 | **-1.00** | 20,629 |
| 2784 | -32.49 | 52.08 | -0.62 | 39,209 |

Basin 2277 is the mirror image of the flips above: it loses its entire recharge
to groundwater export, its runoff depth goes to zero, its aridity index goes to
infinity, and it stops carving. The other two keep a finite index both ways and
cross the threshold downwards on the strength of losing a fifth and two thirds
of their water.

They appear only at the wet arm because that is where the exchange is large
enough, in fraction of a basin's own recharge, to take that much away.

## What GW-8's truncation error does to that answer

The mesh operator sits about 12% above its analytic eigenvalue past `l = 1`,
which is first-order truncation on an irregular mesh and is characterised in
`notes/mesh-geometry.md`. Asking what it does to the carve result means
perturbing every face coefficient by that much and re-solving.

Three members at 12% relative lognormal noise on `w/l`, against the
unperturbed 83:

| member | flips | of those, towards holding |
| --- | ---: | ---: |
| unperturbed | 83 | 0 |
| seed 1 | 84 | 0 |
| seed 2 | 83 | 0 |
| seed 3 | 83 | 0 |

**Both the count and the direction survive.** The spread is one flip in
eighty-three, an order of magnitude smaller than the spread across the
permeability bracket, which is itself small. So the operator's truncation error
is not what limits this result, and the count may be quoted rather than only the
sign.

## GW-14: the criterion's supply term, corrected

The infinite-index branch was never a statement about geometry. It was a
division by a supply term that meant SURFACE runoff, because that was the only
supply a surface-only model had. `carve_verdict.py` now takes

    r_eff = (surface_runoff * C + Qg) / C

when it is handed a water table, and is exactly the surface-only criterion when
it is not. The branch is KEPT for `r_eff <= 0`: a basin that genuinely receives
nothing, including one exporting its entire recharge underground, still never
overflows and still never incises.

`r_eff` is the SUPPLY, so it feeds the lake solver as well as the index. They
take the same quantity and giving them two would be two answers to one question.

**The reduction identity holds exactly.** Run without a groundwater field, the
report is bit-identical to the one the unmodified script produces, key for key.

Supplied with the central arm, on the bootstrap forcing:

| | surface only | with groundwater |
| --- | ---: | ---: |
| basins with no runoff at all | 1016 | **490** |
| carve list, Penman | 1742 | **1852** |
| carve list, robust (carve under BOTH evaporation estimates) | 1742 | **1742** |
| survive under both | 1016 | 897 |
| disputed between the two estimates | 863 | 982 |

110 basins are added to the Penman list and none removed. **The robust list does
not move at all.** The basins groundwater brings above the threshold are ones
that carve under Penman and not under the wet estimate, so they land in the
disputed set rather than in the bracketed one, and the carve set loop A actually
consumes is unchanged on this build.

Nothing regenerates a carve list from this. The head field is uncertified until
GW-3, and a carve list is loop A's input.

## What the solve costs, and where the cost actually was

Measured 2026-09-05. The earlier cost reasoning in `groundwater.py` was taken on
a 400,000-unknown case and concluded that the linear solve is about 30 s inside a
multi-minute run. On `canonical-10m-carve2` it is not: one confined solve ran
past forty minutes at about 16 GB and was killed without writing a water table.
A direct factorisation's work and fill are superlinear in the unknowns, so a
conclusion drawn at 400,000 does not travel to a free set that opens at
millions.

### The lever that was named, and what it measures

The proposal was to hold the evapotranspiration diagonal fixed across inner
steps so the matrix repeats and `splu` re-solves at a fraction of the cost of
factorising. Run on the uniqueness case, which carries GW-15's sink and GW-17's
imposed baselevels:

| | passes | block factorisations | reused | final residual |
| --- | ---: | ---: | ---: | ---: |
| Newton, as run | 10 | 30 | 0 | 2.25e-13 |
| diagonal held fixed | 246 | 15 | 723 | 9.58e-13 |

The reuse arrives. The matrix repeats on 241 of 246 passes, which is the 98 per
cent the lever promised, and the factorisation count halves. The iteration loses
anyway, because it needs twenty-five times the passes and every pass pays an
assembly, a partition and two water balances over the whole mesh whether or not
it factorises. The head lands 2.7e-07 m from the Newton answer rather than
bit-identical, and the frozen arm needs a pass budget of 2,000 to converge at
all against the 60 this component runs.

**The reason is the one this document already gives for a different term.** The
sink is exponential in depth at an e-folding of about a metre, so a
linearisation of `E` about the current head is good only within about a metre of
it, and a Newton step on this problem moves the head by tens of metres. Freezing
the diagonal is a fixed-point iteration on an exponential outside its own
e-folding length, which is the same disqualification Fan's exponential
transmissivity earned above. One exponential, two solvers.

The unfrozen count is worth reading beside it: the WHOLE matrix repeats on 0 of
10 passes, while its PATTERN repeats on 5. The active set settles and the
diagonal does not, so there is nothing for a factorisation to be reused on.

### Where the cost was: the free set is not one problem

GW-17's river and lake cells carry a fixed head, so they are a boundary and not
an unknown, and a face touching one contributes to a diagonal and a right-hand
side rather than to an off-diagonal. The channel network therefore CUTS the free
set: what the solver is handed is one block per interfluve plus one per island,
and the blocks exchange water only through the rivers between them, which are a
Dirichlet condition on both sides. Assembling the union and factorising it whole
pays a superlinear cost on the sum of the blocks that a direct method never has
to pay, and it makes the peak the whole rather than one block. How superlinear
decides how much that is worth, and the exponent is measured below rather than
taken from the textbook.

**Measured on `canonical-10m-carve2`, 2026-09-05**, over its 4,328,732 land
regions, 4,268,074 of them conductive under the `exclude` policy, with 363,761
river and lake cells at fixed head. The free set the solve opens on is 3,914,015
cells and it is **13,094 blocks**. The largest holds 606,808 cells, which is
15.5 per cent of it; fifty-one blocks of a thousand cells or more hold 96.3 per
cent between them and the remaining thirteen thousand are small enough that what
they cost is the per-block call and not the factorisation.

The count is a LOWER bound on the partition, because the solve also drops the
cells that reach neither recharge nor the sea and any that the sink seeds at the
surface, and dropping more cells can only cut further.

**What that is worth, and the exponent is measured rather than assumed.** The
textbook figure for a planar direct solve is `n^1.5`, which would put the sum
over the blocks at 0.280 of the union's. The exponent this file's own recorded
`splu` table actually shows, over 399,424 to 1,999,396 unknowns, is **1.235**,
and its fill grows as about `n (9.5 ln n - 53)` -- so both are nearly linear at
these sizes and the textbook exponent overstates the saving. Extrapolating that
table rather than the exponent:

| | monolithic | partitioned | ratio |
| --- | ---: | ---: | ---: |
| factorisation time | 18.7 s | 10.0 s | 0.53 |
| fill, all blocks held | 355 M | 247 M | 0.69 |
| fill, largest block alone | 355 M | 44 M | 0.12 |

That is an extrapolation from a synthetic at a fifth of the size and it is NOT a
measurement of this solve; it is what sized the design. Two things follow from
it that do not depend on its accuracy. The saving in TIME is about a factor of
two rather than the factor of three and a half the textbook exponent promises,
so the partition is worth having and is not the end of the cost. And the saving
in PEAK is only available if the factors are not all held at once, which is why
the cache is spent only after a zero-flip pass: during the passes where the free
set is largest, nothing is kept and the peak is one block's fill.

It is exact rather than an approximation to the coupled solve. There are no
entries between blocks, so there is no fill between them, and the elimination
inside a block is the sequence of operations it would have been inside the whole
matrix -- provided the cells keep their relative order, which is why the
partition sorts stably. `groundwater.py --factorisation-test` asserts that
bitwise on five interleaved lattices under a scrambled global labelling, and
carries `SymmetricMode` as the control it must reject: that ordering halves the
fill and moves the answer by 2e-14 to 6e-14 relative, which is what says the
equality is a property of this substitution and not of any two ways of solving
the same system.

**End to end against the code it replaces**, which is the bar the change is
actually held to. The partitioned solver was run against the monolithic
`spsolve` version on the uniqueness case at three sizes, 48x48, 96x96 and
128x128, each carrying the sink, the imposed baselevels, an enclosed block and a
pinned wet patch: the same pass count, an **identical converged active set**, and
a **bit-identical head field**, `max|dh| = 0.000e+00` on all three. That is the
whole of what a restructuring of this solve has to show, and it is not inferred
from the algebra arm above but measured on the solver.

`splu(A).solve(b)` replaces `spsolve(A, b)` for the same reason and at the same
bar: at one column ordering the two drive the same factorisation and the same
back-substitution, so they agree to the last bit, and the same check asserts it.

## What remains

**The external test.** GW-3, the same code on Earth topography, Earth recharge
and GLHYMPS permeability, scored against Fan et al. (2013)'s 1,603,781 well
sites. Every check that has run so far is internal -- an identity, a
conservation law, or a reduction -- and none of them can say the model is right
about a real water table. That one can, and it needs external datasets.

**A baseline climatology.** Everything above is forced by a bootstrap run.
