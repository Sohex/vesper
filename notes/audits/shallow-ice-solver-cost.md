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
