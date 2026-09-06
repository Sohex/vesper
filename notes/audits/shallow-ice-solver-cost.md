# What a shallow-ice flow solver costs on this project's grids

Worldbuilding frame: this prices a numerical solver against the Vesper project's
own mesh and its climate model, on this desktop. Every field named is a model
array or a mesh quantity; nothing here is about the real world.

CLIM-62 asks for a cost for each route to a shallow-ice solver and a verdict,
not a solver. The comparison is scoped to the FLOW SOLVER ALONE: neither route
supplies the surface mass balance, which is CLIM-63's object and is most of the
work, so a comparison that includes it measures the wrong thing.

## The rule and the criteria, fixed before anything was run

Declared 2026-09-05, before the prototype existed.

**The row's own rule, restated so it can be applied.** The native route wins
unless the port is cheaper by a wide margin, because it carries no grid crossing
and no external dependency. Read "wide margin" as **a factor of three or more in
implementation size**, since that is the only axis on which the two routes can
be compared before either exists: 1.2k lines of adapter around a 27k-line solver
is the measured size of the CLIMBER-X interface (`notes/external-model-survey.md`
section 4), and Route A's own operator would be smaller than the adapter alone.

**And the prior gate, which is cheaper than either cost.** If the doubly
nonlinear diffusion will not converge on this mesh then Route A is unavailable at
any cost and the comparison is over. That gate is tested against an identity with
a right answer -- the Halfar similarity solution for shallow ice at zero mass
balance, which conserves volume exactly and whose margin radius is closed form --
and never against an ice sheet that merely looks plausible.

Three thresholds on that gate, each fixed here:

- **CONVERGES** means the nonlinear residual falls monotonically to a declared
  tolerance in a bounded outer-iteration count, and does NOT mean the head step
  falls. GW-9's exponential transmissivity held a water balance residual flat at
  1.946e-04 from pass 100 to pass 599 while the head step fell thirtyfold, so a
  head-step bar would have passed the failure. The residual is the bar.
- **ACCURATE** means the margin radius lands within one cell of Halfar's
  R(t) = R0 (t/t0)^(1/18), and the profile's relative RMS error sits at or below
  the discrete operator's own noise floor. GW-8 measured that floor at about 12
  per cent relative RMS above l = 1 on this mesh family, and measured that it does
  not converge with refinement -- fitted order 0.004 to 0.005 in cell size -- so
  it is a property of the jittered-Fibonacci Voronoi family and it bounds what any
  answer here can mean. A result better than 12 per cent is not evidence of a
  better scheme; a result worse than it is evidence of a worse one.
- **CONSERVES** means the total volume is held to the linear solver's own
  tolerance over the integration, because Halfar's mass balance is zero and
  volume conservation is then an identity rather than a comparison.

**What would make the numbers not numbers.** The instrument's noise floor is
GW-8's 12 per cent, so an accuracy claim finer than that is noise however tidy.
And a timing on this host is a measurement of a machine state as much as of a
solver: every wall clock below carries the load it was taken under, and any arm
taken above a load of about 4 is reported as unusable rather than quoted.

**Rule 7 applies here too and makes the arithmetic smaller.** Nothing is
commissioned, so no route is charged for output it would invalidate. The cost of
each route is what it costs to build and to run once.

## The gate: the diffusion converges on this mesh, and what the damping is for

Measured 2026-09-05 with `analysis/shallow_ice_route_a.py` -> `analysis/shallow_ice_route_a.json`,
on synthetic meshes of 20,000 and 80,000 regions from Orogen's own generator at
the jitter the refinement sweep runs at, with a Halfar dome of 3,000 m and a
2,500 km radius. The dome is sized so 13 cells span its radius at the coarse arm
and the sphericity contamination is 1.78 per cent, below GW-8's 12 per cent
floor.

**THE LOAD, AND WHICH NUMBERS IT REACHES.** The sweep ran under
`scripts/lock_and_run` on a host whose one-minute load went from 1.41 to 18.34
across it, with the fifteen-minute figure above 21 throughout: the lock was held
but the machine had not settled from what came before. So **every wall clock in
this section is a shape and not a price**, and is marked as such. The iteration
counts, the residual traces, the profile errors, the margin positions and the
volume change are arithmetic and are unaffected by what else the machine was
doing; they are the numbers the gate is decided on.

**THE TWO-POINT SLOPE IS WHAT MAKES THE MAP NON-CONTRACTIVE, AND THAT IS THE
FINDING.** `D` depends on `|grad z_S|^(n-1)`, and a two-point flux across a face
gives only the normal component. With that approximation, undamped Picard
DIVERGES at every step size tested, in exactly GW-9's shape: the relative
residual oscillates with growing amplitude -- 1.65e-3, 1.18e-3, 1.62e-3,
9.89e-4, 2.99e-3, 1.04e-3, 4.93e-3 over the first seven passes at the explicit
diffusion limit -- and sixty passes reach neither the tolerance nor a stall. A
head-step bar would have called that convergence, which is why the bar is the
residual.

**With the full gradient reconstructed, undamped Picard CONVERGES**, in 7 passes
at 20,000 regions and 5 at 80,000, at the explicit step. That is the same
operator, the same mesh and the same dome; what changed is that `D` now sees the
whole surface slope rather than its component across each face. The transverse
component is not a refinement of the diffusivity, it is what makes the fixed
point attractive.

**Damping recovers the normal-only arm, and buys the large step for both.**

| slope | step, x the explicit limit | undamped | damped at 0.5 |
| --- | ---: | ---: | ---: |
| normal, 20k | 1 | diverges | 18 passes |
| normal, 20k | 20 | diverges | 24 |
| full, 20k | 1 | **7** | 20 |
| full, 20k | 20 | diverges | 49 |
| normal, 80k | 1 | diverges | 16 |
| normal, 80k | 20 | diverges | 21 |
| full, 80k | 1 | **5** | 17 |
| full, 80k | 20 | diverges | 26 |

Damping is not free where it is not needed: at the explicit step the full arm
costs 20 damped passes against 7 undamped. The disposition a solver would take
is a line search rather than a fixed weight, and this measurement is what says
so.

**The three thresholds fixed above, against the result.**

- **CONVERGES: PASSES.** At the explicit step the full-gradient form converges
  undamped in 5 to 7 passes and every arm converges damped, monotonically. At
  twenty times the step every arm needs damping and every damped arm converges,
  in 21 to 49 passes.
- **ACCURATE: PASSES, and it improves with refinement.** The margin radius lands
  0.35 to 0.91 cells beyond Halfar's closed form at every arm, inside the
  one-cell bar. The profile's relative RMS is 1.33 to 1.60 per cent at 20,000
  regions and **0.203 to 0.233 per cent at 80,000** -- an order better for four
  times the cells, so the SOLVED field converges even though GW-8 measured the
  operator's own truncation error not converging at all. That is the same
  supraconvergence GW-8 found for the water table, reproduced here for a
  nonlinear operator.
- **CONSERVES: PASSES exactly.** The volume change is 0.0 to the last reported
  digit at every one of the sixteen arms, and by construction rather than by
  luck: the face fluxes are antisymmetric and the assembled operator's row sums
  are zero, so the implicit step conserves `sum A H` and only the active set's
  clip could break it. It did not fire on a Halfar dome, whose support shrinks.

**WHAT THIS TEST CANNOT DISCRIMINATE, and it has to be said because the table
above invites the wrong reading.** Halfar's dome is radially symmetric, so its
surface gradient is radial and the normal component across a face between two
neighbours is very nearly the whole of it. That is why the two slope arms agree
on the profile error to within a fifth of a per cent, and it means **this
identity says nothing about which slope is more accurate on a real ice sheet**,
where the transverse component is not small. What it does discriminate, and
decisively, is the convergence behaviour above. An accuracy comparison between
the two needs a test whose surface slope is not radial, and there is not one
here.

**The gradient reconstruction is itself checked against an identity** before
either arm is quoted: on a sphere `H = R cos(theta)` has tangential gradient
magnitude exactly `sin(theta)`, and the least-squares operator reproduces it to
a median 0.30 per cent, 1.5 per cent at the ninth decile, over 19,602 cells away
from the poles. That is an order below the two-point operator's own noise floor,
so the `full` arm measures the scheme and not the reconstruction.

**SO ROUTE A IS AVAILABLE.** The prior gate the row declares is passed, and the
comparison is live.

## The cost, run twice at different loads, and what survives that

The sweep was run a second time on 2026-09-06 to get a quiet timing and the host
did not oblige: the one-minute load was 1.41 rising to 18.34 on the first pass
and 17.93 rising to 28.87 on the second. So there is no quiet number here. What
two passes at different loads DO give is a bracket and a control, and the
control is the sharper half.

**EVERY ITERATION COUNT REPRODUCED EXACTLY.** 18, 24, 7, 20, 60, 49 at 20,000
regions and 16, 21, 5, 17, 60, 26 at 80,000, identical across the two passes, as
did every profile error, every margin position and every volume change. That is
what says the gate's verdict is arithmetic and not a property of the machine.

**THE ABSOLUTE COST IS A BRACKET AND IS QUOTED AS ONE.** Per linear solve, 0.011
to 0.019 s at 20,000 regions and 0.037 to 0.052 s at 80,000, the two ends being
the two passes. The second pass is 30 to 50 per cent dearer throughout, which is
the load and not the code.

**THE SCALING IS NOT A BRACKET, because a ratio taken inside one pass carries
its own load away with it: 3.4 times the cost for 4 times the regions, in BOTH
passes.** That is close to linear in the region count for this sparsity, and it
is the number a reader should extrapolate with rather than either absolute end.

**And what it extrapolates to is the real cost item.** One sparse factorisation
per Picard pass, at 5 to 26 passes a step, is affordable at these mesh sizes and
is not the question. At the active build's 10,000,005 regions a DIRECT
factorisation is the wrong solver: GW-8 already found the LINEAR water-table
case there beyond this host's memory, at roughly 834 million factor entries and
about 10 GB, and a doubly nonlinear operator re-factorises every pass rather than
once. **A production Route A needs an iterative solver with a preconditioner,
and that is the piece this gate did not price.** The matrix is a symmetric
M-matrix at every pass, which is the condition under which a preconditioned
conjugate gradient is straightforward, so this is a known-shaped piece of work
rather than an open question -- but it is not free and it is not in the 600 to
900 lines above.

## Route B, read against the artifact, and it is not what the row assumed

Read 2026-09-05 from `references/big-mitgcm/MITgcmIS.py`. The row's open question
was the cost, the separability having been settled. Four things settle it
without a cost.

**IT IS NOT A MODULE.** The entire time loop -- 30,000 steps of assembly, solve
and isostasy -- lives inside a single string passed to
`get_ipython().run_cell_magic('time', '', ...)` at line 548. It is a notebook
export. There is no function to call and nothing to import; a port begins by
rewriting the loop, which is where all of the physics is.

**IT DOES NOT SOLVE THE FREE BOUNDARY EITHER.** `h_vecnew[h_vecnew < 0] = 0`
after the linear solve, and that is the whole of it. No active set, no
complementarity condition, no variational inequality. This is the same defect
section 45c-d of the survey found in SICOPOLIS and it is the reason Route C was
called a downgrade; Route B has it too. Route A's active set gets the margin
right by construction.

**IT DOES NOT ITERATE THE NONLINEARITY AT ALL.** `D` is evaluated once per step
from the previous step's thickness and one `cgs` solve follows. So the arms above
are not measuring the same thing MITgcmIS does: Route A converges the nonlinear
residual to 1e-8 and Route B lags `D` by a whole step and never asks. That is
what its 30,000 steps buy, at a step that is a bare module-level literal
assigned twice, `dt = 1.0` at line 223 and `dt = 1` at line 226, with no unit
stated anywhere and no stability check against the diffusivity it just built.
The unit is a year by inference from `Aglen`, and the inference is the reader's.
The `cgs` call carries no preconditioner and its failure flag is printed rather
than acted on.

**Whether that is cheaper per unit of simulated time is not answerable from
either code**, and it is worth saying so rather than estimating it. A lagged-`D`
scheme is one solve a step at a step its author chose and did not justify; a
damped Picard is 18 to 24 solves at a step this note measured against the
explicit limit. They are close enough that the comparison would turn on the
lagged scheme's actual stability limit, which is not stated and not tested in
that tree. What is not in doubt is that only one of the two has a converged
answer at the end of its step.

**AND EARTH'S GRAVITY IS INSIDE THE FLOW CONSTANT AT THE THIRD POWER.**
`rhog = 920 * 9.8`, then `Csia = 2*Aglen/(nglen+2) * rhog**nglen` at line 227.
That is GRAV-6's 2.23x, welded into a module-level literal. Beside it the surface
mass balance is an Earth-latitude equilibrium line -- 5,500 m poleward of 30
degrees, 5,000 m equatorward -- hardcoded in `smb()` and called from inside the
same loop, which is the scheme CLIM-63 refused.

**THE LICENCE IS NOT STATED.** There is no licence file in the tree.
`READ_ME.rtf` says only that "MITgcmIS is freely available in this ZENODO
repository under MITgcmIS.py", which is a statement about availability and not a
grant. `references/INDEX.md`'s own shape for this is the ArcSDM row: not stated,
check at point of use. Adoption cannot proceed on it as it stands.

## The verdict, against the rule declared above

**ROUTE A WINS, and not on a margin.** The rule was that the native route wins
unless the port is cheaper by a wide margin, read as a factor of three in
implementation size. It does not get there:

| | what would be written | the free boundary | the nonlinearity | gravity | licence |
| --- | --- | --- | --- | --- | --- |
| A, native on the region mesh | the operator, the Picard with a line search, the active set and the least-squares gradient, on `groundwater.py`'s existing `Geometry` | active set, exact by construction | converged to 1e-8 in 5 to 7 passes undamped with the full slope, 16 to 49 damped | `config/planet.yaml` | this project's |
| B, port MITgcmIS | the whole time loop, which is a notebook cell string, plus a cubed-sphere crossing under rule 3 | clipped | lagged one step, never asked | Earth's, at the third power, in a module literal | not stated |
| C, SICOPOLIS behind CLIMBER-X's adapter | ~1.2k lines of adapter around a 27k-line solver | clipped, then rebooked as a mass balance correction | its own | one parameter and two literals, with Earth-tuned limiters on `g^3` and `g^4` | GPL-3 over GPL-2-or-later |

Route B is not cheaper at all, let alone by three times, and the two things a
port is supposed to buy -- a working solver and someone else's maintenance -- are
both absent: there is no callable solver, and no licence under which to maintain
one. Route C remains what section 45c-d found, a downgrade on the free boundary,
and of the two PORT routes it is the only one with a licence and a callable
solver behind it. Route A needs neither, being this project's own code on this
project's own mesh, which is the half of the row's rule the measurement above
turns from an assumption into a fact.

**What Route A would actually cost to build.** The prototype is 400 lines
including its own test harness, and it reuses `hydrography/scripts/groundwater.py`'s
`Geometry` for the faces, widths and areas -- the piece that took a wrong first
attempt and put the operator 178x off its analytic eigenvalue. A production
version is the operator, the Picard with a line search, the active set, the
least-squares gradient and a mass balance input: call it 600 to 900 lines,
bracketed because the surface mass balance's interface is CLIM-63's and does not
exist yet. The bracket does NOT cover the linear solver: a direct factorisation
does not reach the active build's region count and an iterative solver with a
preconditioner is a separate piece of work this gate did not price.

## Three things Route A would still owe, and they are named rather than hidden

**THE TRANSVERSE SLOPE IS NOT OPTIONAL, and the gate turned that from a
suspicion into a measurement.** `D` depends on `|grad z_S|^(n-1)` and a two-point
flux gives only the normal component; with that approximation the Picard map is
not a contraction at any step size tested, and with the full gradient it
converges undamped in five to seven passes. MITgcmIS does it properly and only
because its grid is structured: `dsdx_c` and `dsdy_c` are centred east-minus-west
and north-minus-south differences, which a Voronoi mesh has no equivalent of. On
this mesh the full gradient is a per-cell least-squares solve over each cell's
neighbours, precomputable but not free, and it is now a REQUIREMENT of Route A
rather than a refinement of it.

**THE ACTIVE SET DOES NOT INHERIT THE GROUNDWATER SOLVER'S TERMINATION PROOF.**
That proof rests on the matrix being fixed and SPD, and here `D` depends on the
unknown. The project already downgraded the check from an identity to a 1e-6
measurement for GW-24's merely-Picard nonlinearity. What replaces it here is the
Halfar identity, which is stronger in one way -- it has a right answer the
solver cannot see -- and weaker in another: it certifies the answer and not the
algorithm.

**AND `saturated_thickness`'S POSITIVE FLOOR EXISTS FOR A REASON THAT AN ICE
MARGIN BREAKS.** `groundwater.py` keeps a strictly positive floor specifically so
the conductive graph does not change with the unknown, which is what makes its
static dry-set argument work. A shallow-ice margin is a changing connectivity by
construction. Route A therefore reuses the `Geometry` and the SHAPE of the active
set, and not that argument; the active set here is written fresh and is judged by
the margin radius against Halfar's closed form.
