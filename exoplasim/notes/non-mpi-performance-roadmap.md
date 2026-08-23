# Non-MPI routes for making ExoPlaSim faster

*This is a COMPUTE roadmap for the Vesper worldbuilding project's climate
model. It excludes the OpenMP/SHTns workstream and every MPI-specific layout or
build question. Written 2026-08-21 from the measurements and source audits
already in this repository.*

## The answer in priority order

The one source-level optimization this document opened with -- store the
associated Legendre basis as P and Q plus separable factors rather than as eight
full weight matrices -- was built, measured and refused. Section 1 carries the
result, and it removes the premise the priority rested on.

What is left, and none of it is another compiler flag or a new math library:

1. measure the largest scientifically acceptable timestep at each resolution;
2. use resolution and precision staging to reduce total spin-up cost;
3. profile and reduce repeated transcendental work in radiation;
4. specialize expensive optional diagnostics and physics out of production
   builds when they are disabled.

Mixed precision for read-only transform weights was on this list and is not any
more: SHTns stores no weight matrices, so there is nothing immutable left in the
transform to demote. Section 5 says so.

The FFT, generic compiler tuning, alternate `libm` implementations, and model
output are already measured too small or actively slower on this host. They are
controls, not priority work.

## 1. Store P and Q, not eight Legendre matrices: built, and refused

Implemented in full, correct, four times smaller, and **1.56% SLOWER** on the
threaded build at T170 -- 80.57 s against 81.66 s, paired and interleaved,
faster in 0 of 4 rounds, [-2.17, -0.74]. It took the weights from 114.9 MB a die
to about 36 and from 15.1 MB a thread to 4.5, which is under the 32 MB target,
and bought nothing. The code is reverted; `verify_weight_factorisation_model.sh`
is kept, and it agreed with the eight-matrix build to 4.8e-13 at one step and
1.06e-11 at twenty against a 1e-10 bar declared before the arms ran.

**The premise was backwards.** The route was ranked first because the filter
fold -- which moved a per-mode scalar INTO the matrices to save one multiply in
three -- had returned 1.6%, read as evidence that these loops are bound by
STREAMING the weight matrices. Undoing the fold costs the same 1.6% back, and
the reading that fits both numbers is that the loops are bound by the MULTIPLY.
A second measurement from the other end agrees: at NPRO=1 one thread streams
30 MB of `pmat` and 30 MB of `qmat` per transform against 1.88 MB of each at
NPRO=16, and the two time the same, 1.45 ms against 1.48. Sixteen times the
weight traffic costs nothing measurable.

**And SHTns has since removed the matrices altogether**, computing the Legendre
functions on the fly and storing no `NCSP x NLPP` array at all, which takes
`pmat` and `qmat` out of `cache_budget.py` entirely -- 30.1 MB a die at T170,
the whole of what this route was for. The cache-boundary finding in
`rank-imbalance-and-weight-traffic.md` still stands; what no longer stands is
that this is the way to act on it.

CLIM-48, closed wontfix. Reopen only if a resolution above
T170 makes the memory headroom worth a small speed loss.

## 2. Price the timestep directly

The number of dynamics steps per simulated orbit scales inversely with
`MPSTEP`. A scientifically acceptable increase from 45 to 60 minutes would
remove 25% of the steps, a larger theoretical prize than the remaining
micro-optimizations at T42.

This is not a free performance flag. A larger step changes temporal truncation,
filter behavior, phase sampling, and possibly model stability. It therefore
needs a convergence experiment rather than a short timing bed.

### Experiment

For every production resolution, test a small declared ladder around the
current timestep:

- cold-start stability with floating-point traps enabled;
- warm-start one-orbit energy, water, and extrema checks;
- climatological differences after settling;
- wall time per simulated orbit, not per model step;
- identical output sampling in physical time.

The chosen step belongs in the resolution configuration. It must not be an
automatic rule inferred only from whether a short run crashes.

Restart conversion across different timesteps remains out of scope for the
first converter because the file stores a step count and leapfrog history, not
an elapsed-time state. Until that is implemented, timestep arms need their own
compatible starts.

## 3. Stage resolution and precision during spin-up

This is a time-to-solution optimization rather than a faster timestep. Most of
an equilibrium spin-up is spent removing large, smooth imbalances that do not
require the final horizontal resolution. A lower-resolution or single-precision
run can act as a preconditioner, after which the final configuration settles
the smaller scales and supplies the convergence diagnostics.

The prerequisite is the restart converter specified in
`restart-resolution-precision-converter.md`. Without it, savings in a cheap
spin-up cannot be transferred safely into the production state.

### Safe initial policy

- Low resolution and/or four-byte precision may be used only for exploratory
  spin-up.
- Convert into the final resolution and eight-byte build.
- Reset every diagnostic accumulation window at conversion.
- Require a declared settling interval at the final configuration.
- Apply all convergence gates only to the final eight-byte integration.

Four-byte accumulation is not suitable for the existing sub-watt convergence
thresholds. `parameter-decisions.md` bounds its annual summation error at about
0.03 to 0.2 W/m2, which is comparable to those thresholds. Promoting an
existing four-byte state to eight bytes also does not recover information
already rounded away.

The experiment must report total wall time to a converged final state. A fast
coarse run followed by a long resolution-adjustment transient has saved
nothing.

## 4. Profile the radiation transcendental functions

On the clean T85 transform profile, the largest unclassified compute was
`pow` and other `libm` work called from radiation. That is a more credible next
physics target than the FFT, but it is not yet a priced optimization.

The first step is a call-site profile that attributes transcendental cost to
individual expressions and radiation bands. Only then consider:

- hoisting quantities invariant over level, latitude, or radiation calls;
- replacing repeated powers with products when the exponent is a small exact
  integer;
- evaluating a shared expression once when several bands consume it;
- bounded lookup/interpolation for empirical absorption fits, with an error
  budget declared in radiative flux rather than in the function value.

Do not begin by substituting another math library. AMD LibM was measured 13.4%
slower than glibc on the T42 bed and changed the integration. Any approximation
or reassociation here is a numerics change and needs radiative unit tests plus
an end-to-end climate control.

## 5. Mixed precision for immutable transform data: no data left to demote

This asked whether immutable P/Q weights or selected scratch arrays could be
stored in four bytes while accumulation and prognostic state stayed in eight,
on the reading that the Legendre loops are memory-traffic limited. Both halves
of that premise are gone. The loops are multiply-bound, measured twice and
recorded in section 1, and SHTns computes the Legendre functions on the fly and
stores no weight matrices at all, so there is no immutable transform data in the
per-timestep path to demote.

What survives from this line is section 3's precision staging, where the
question is the precision of the RUN and not of a table, and CLIM-59 is where
it is tracked. The checks listed here belong to that row.

## 6. Compile or allocate optional diagnostics out of the hot path

The old T127 profile spent 5.07% in diagnostics, principally the
`(NHOR,NLEV,28)` energy array. Production disabled those diagnostics and the
bucket fell to 0.73%, while the 600-step bed fell from 120 to 73 seconds after
that and removal of runtime bounds checking.

The general rule is worth preserving:

- disabled diagnostics should allocate no large arrays;
- their per-step call sites should be absent or guarded outside hot loops;
- diagnostic builds should be distinct from production builds in provenance;
- a profiling bed must match the production namelist before its percentages
  are used to rank work.

There may be more optional modules with the same shape, but source inspection
alone is insufficient. Re-profile the future threaded production build and
look for disabled features in the sample table before changing them.

## Routes already priced too small

### Compiler flags and CPU libraries

The current `-march=znver4` gain is 2–3%. Loop unrolling and software
prefetching did not improve the production line, LTO was around 1.4%, and the
entire scalar-to-AVX-512 distance was only a few percent. AOCL LibMem was 41%
slower. Profile-guided optimization may still be screened cheaply, but these
measurements bound the likely prize and do not justify a large workstream.

Runtime bounds checking should remain a diagnostic-build feature. It cost 5%
at T42 and 17% at T127, but floating-point traps cost under 1% and remain worth
keeping because they preserve the origin of invalid arithmetic.

### FFT replacement by itself

The FFT bucket was about 4% at both T85 and T127. It scales more slowly than the
Legendre work and does not become the dominant high-resolution problem. A
standalone FFT rewrite therefore has only a few percent available before its
own overhead.

### Model output

On the measured 600-step T42 bed, restoring production output moved 12.678 s to
12.719 s while writing 230 MB, inside run-to-run noise. Reducing output can save
storage and downstream processing, but it is not a meaningful model-compute
optimization on this system.

### Accumulator loop fusion

`outaccu` runs every step, but on the corrected production profile it is around
one percent. Its statements touch separate arrays with little reusable data, so
fusion has neither a large ceiling nor an obvious locality win.

## Measurement order

1. Establish the completed OpenMP/SHTns build as the new baseline, without
   attributing any of that work here.
2. Run the timestep ladder at the resolutions that matter operationally.
3. Build the restart converter and measure staged spin-up end to end.
4. Attribute radiation `libm` time to call sites.
5. Re-profile before opening any compiler or diagnostic specialization work.

Mixed-precision weights (5) fell with the P/Q factorization: SHTns stores no
weight matrices, so there is no immutable transform data left to demote.

Every performance A/B should use interleaved, order-flipped rounds on a quiet
machine, include warm-ups, record restart differences, and compare the gain
against self-scatter. A short benchmark establishes throughput; it does not by
itself establish scientific equivalence.

