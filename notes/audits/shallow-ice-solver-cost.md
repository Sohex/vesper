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

## The gate: the diffusion converges on this mesh, and the damping is the price

Measured 2026-09-05 with `analysis/shallow_ice_route_a.py`, on a synthetic
mesh of 20,000 regions from Orogen's own generator at the jitter the refinement
sweep runs at, with a Halfar dome of 3,000 m and a 2,500 km radius. The dome is
sized so 13 cells span its radius and the sphericity contamination is 1.78 per
cent, well below GW-8's 12 per cent floor. The host was contended throughout, so
every wall clock in this section is a shape and not a timing; the timings that
count are in the section below, taken under the lock.

**UNDAMPED PICARD DIVERGES, and it does so in the shape GW-9 recorded.** At a
step equal to the explicit diffusion limit the relative residual OSCILLATES with
growing amplitude -- 1.65e-3, 1.18e-3, 1.62e-3, 9.89e-4, 2.99e-3, 1.04e-3,
4.93e-3 over the first seven passes -- and at five times that step it climbs to
0.6 within ten. Sixty passes reach neither tolerance nor a stall; the map is not
a contraction, which is exactly what a conductance that COLLAPSES as the unknown
falls should do. `D` goes as `H^(n+2)`, so it is GW-9's sign and not GW-15's.

**DAMPED PICARD CONVERGES, and cheaply.** Under-relaxing the update at w = 0.5
the residual falls monotonically and geometrically at a ratio of about 0.49 --
1.48e-3, 6.83e-4, 3.28e-4, 1.59e-4, 7.80e-5, 3.84e-5, ... -- reaching the
declared 1e-8 in 18 passes. w = 0.25 converges in 43 and w = 0.1 does not reach
the tolerance in 60, so the damping has an optimum rather than a monotone
benefit and 0.5 is inside the basin rather than at its edge.

**And it holds as the step grows, which is the point of an implicit scheme.**

| step, as a multiple of the explicit diffusion limit | passes | profile RMS |
| --- | ---: | ---: |
| 1 | 18 | 0.0133 |
| 2 | 19 | 0.0199 |
| 5 | 21 | 0.0333 |
| 20 | 24 | 0.0553 |

Twenty times the step costs six extra passes. The accuracy degradation is time
truncation and every arm stays below the 12 per cent operator floor.

**The three thresholds fixed above, against the result.**

- **CONVERGES: PASSES at w = 0.5, FAILS undamped.** The bar was a monotone fall
  of the nonlinear residual to 1e-8 in a bounded pass count, and the residual
  was the bar precisely because GW-9's head step would have passed the failure.
- **ACCURATE: PASSES.** The margin radius lands 0.61 to 0.81 cells beyond
  Halfar's closed form at every step size, inside the one-cell bar. The profile's
  relative RMS is 1.33 per cent at the explicit step, an order below the 12 per
  cent floor, and 5.53 per cent at twenty times it, still below.
- **CONSERVES: PASSES exactly.** The volume change is 0.0 to the last reported
  digit at every arm, which is not a coincidence and is worth saying why: the
  face fluxes are antisymmetric and the assembled operator's row sums are zero,
  so the implicit step conserves `sum A H` by construction and the only thing
  that could break it is the active set's clip. It did not fire on a Halfar
  dome, whose support shrinks rather than grows.

**SO ROUTE A IS AVAILABLE.** The prior gate the row declares is passed, and the
comparison is live.

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
| A, native on the region mesh | the operator, the damped Picard and the active set, on `groundwater.py`'s existing `Geometry` | active set, exact by construction | converged to 1e-8 in 18 to 24 passes | `config/planet.yaml` | this project's |
| B, port MITgcmIS | the whole time loop, which is a notebook cell string, plus a cubed-sphere crossing under rule 3 | clipped | lagged one step, never asked | Earth's, at the third power, in a module literal | not stated |
| C, SICOPOLIS behind CLIMBER-X's adapter | ~1.2k lines of adapter around a 27k-line solver | clipped, then rebooked as a mass balance correction | its own | one parameter and two literals, with Earth-tuned limiters on `g^3` and `g^4` | GPL-3 over GPL-2-or-later |

Route B is not cheaper at all, let alone by three times, and the two things a
port is supposed to buy -- a working solver and someone else's maintenance -- are
both absent: there is no callable solver, and no licence under which to maintain
one. Route C remains what section 45c-d found, a downgrade on the free boundary,
and it is the only one of the three with a licence and a real solver behind it.

**What Route A would actually cost to build.** The prototype is 400 lines
including its own test harness, and it reuses `hydrography/scripts/groundwater.py`'s
`Geometry` for the faces, widths and areas -- the piece that took a wrong first
attempt and put the operator 178x off its analytic eigenvalue. A production
version is the operator, the damped Picard, the active set and a mass balance
input: call it 600 to 900 lines, bracketed because the surface mass balance's
interface is CLIM-63's and does not exist yet.

## Three things Route A would still owe, and they are named rather than hidden

**THE TWO-POINT FLUX SCHEME CANNOT SEE THE TRANSVERSE SLOPE.** `D` depends on
`|grad z_S|^(n-1)`, and a two-point flux across a face gives only the normal
component of that gradient. MITgcmIS does it properly and only because its grid
is structured: `dsdx_c` and `dsdy_c` are centred east-minus-west and
north-minus-south differences, which a Voronoi mesh has no equivalent of. On this
mesh the full gradient is a per-cell least-squares solve over each cell's
neighbours, precomputable but not free. Both arms are implemented and the sweep
below prices the difference.

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
