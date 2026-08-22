# Where the sixteen ranks actually spend their time

*This is a COMPUTE measurement of the Vesper worldbuilding project's climate
model on this desktop. Nothing in it is about the simulated planet. Measured
2026-08-21, on the build carrying the symmetric transforms.*

## The question

`spectral-transform-profile.md` established that the `mpi` bucket is not
addressable by making the transform faster: Open MPI busy-waits, so a rank that
reaches a collective early burns CPU that perf counts as work, and the SPREAD
across ranks is the load imbalance rather than the cost of communicating. It
did not say how large that imbalance is, or what causes it. Both are needed
before any of the remaining structural options -- threading, an uneven
decomposition, a transposed layout -- can be priced.

## The imbalance is the die boundary, and nothing else

Sixteen ranks, one per core, Open MPI's default binding, which was checked
rather than assumed: rank *r* lands on core *r*.

| | T127 | T170 |
| --- | ---: | ---: |
| ranks 0-7, CCD0, 96 MB L3 | 23.12% of samples waiting | 30.59% |
| ranks 8-15, CCD1, 32 MB L3 | **12.35%** | **20.13%** |
| gap | **10.77 points** | **10.46 points** |
| spread WITHIN CCD0 | 3.5 points | 11.4 points |
| spread WITHIN CCD1 | 0.9 points | 1.3 points |

Every rank does identical work, so the ones waiting LESS are the ones setting
the pace. **CCD1 is the critical path at both resolutions**, and the step at the
die boundary is an order of magnitude larger than the variation within either
die.

**It is the cache and not the clock.** On this processor CCD1 is the die WITHOUT
V-cache and it clocks higher; the slower half is the half with more clock and
less cache. That direction is the wrong way round for a frequency explanation
and the right way round for a working-set one.

It also reproduces `smt-rank-layout.md`'s 10.5%, which was measured from wall
clock on single jobs pinned to one die or the other. Two methods with nothing in
common, the same number.

## What balancing is worth: about five percent, and less than that in practice

The ceiling is `mean - min` across ranks, which is what perfect redistribution
would recover: **5.74% of wall at T127, 5.65% at T170.** That number bounds
every balancing scheme equally -- an uneven decomposition, OpenMP with dynamic
scheduling, anything.

The practical figure is lower, because **latitudes have to come in pairs**. The
symmetric transforms need an even count per rank, so the granularity of any
uneven split is two latitudes, not one:

| | per-latitude cost ratio CCD1/CCD0 | latitudes a CCD0 rank should hold | nearest even | worth |
| --- | ---: | ---: | ---: | ---: |
| T127, NLPP 12 | 1.140 | 12.78 | 12, which is what it already holds | **nothing** |
| T170, NLPP 16 | 1.151 | 17.12 | 18 | +2.24% of compute, about **1.8% of wall** |

At T127 the optimum is nearer the current split than the next even one, and
14/10 is 2.3% WORSE than 12/12. So an uneven decomposition has nothing to give
there at all, and 1.8% at T170.

**That is the answer on CLIM-47 and it is a refusal.** Making `NLPP` a per-rank
variable turns every `(NHOR,...)` array into a maximum rather than an exact size
and every loop bound over it into a runtime value, across the whole physics. It
cannot be worth 1.8% at one resolution and nothing at another.

The same ceiling refuses OpenMP on the balance argument. Threading would capture
at most 5.7%, it buys nothing at all on a machine whose cores are alike, and it
would mean auditing about five hundred module-level `NHOR` arrays for
thread-safety. The other argument for threading -- that SHTns needs every
latitude in one process -- is answered by a transposed mode-space decomposition
just as well, in MPI, without touching how the physics is parallelised.

## The cause is addressable directly, and that is the interesting result

The working set is the Legendre weights. `legini` builds EIGHT `NCSP x NLPP`
matrices, and every one of them is touched every timestep:

| | eight matrices, a rank | per die, eight ranks |
| --- | ---: | ---: |
| T127 | 6.05 MB | **48.4 MB** |
| T170 | 14.36 MB | **114.9 MB** |

against 96 MB of L3 on CCD0 and 32 MB on CCD1. That is the whole explanation of
the table at the top: CCD0 holds it, CCD1 does not.

**All eight are the same two arrays times separable factors.** Reading `legini`:

    qi = P * A(w)                 qc = P * B(w) * gwd(l)
    qj = Q * A(w)                 qe = Q * B(w) * gwdcsq(l)
    qu = P * A(w) * n1(w) * m     qq = P * B(w) * gwdcsq(l) * n(n+1)/2
    qv = Q * A(w) * n1(w)         qm = P * B(w) * gwdcsq(l) * m

with `A = skspgp`, `B = skgpsp`, `n1 = 1/(n(n+1))`, and `gwd`, `gwdcsq`
functions of the LATITUDE alone. So every matrix is P or Q scaled by a per-mode
factor and a per-latitude factor, and neither factor needs an `NCSP x NLPP`
array to hold it: the per-mode ones are `NCSP` vectors, the per-latitude ones
are scalars inside the latitude loop.

Storing P and Q alone would leave:

| | two matrices, a rank | per die |
| --- | ---: | ---: |
| T127 | 1.51 MB | **12.1 MB** |
| T170 | 3.59 MB | **28.7 MB** |

Both fit CCD1's 32 MB. If the die asymmetry is the working set, this removes it
rather than balancing around it, and it helps CCD0 as well.

**It is the exact inverse of the filter fold, and that is the argument for it.**
The fold moved a per-mode scalar INTO the matrices to save one multiply in
three, and measured 1.6% at T127 -- because, as `legendre-filter-fold.md` says,
those loops are bound by streaming the weight matrices rather than by the
multiplier. The same sentence read the other way says trading multiplies for a
quarter of the traffic should win, and by more.

## What would have to be established before believing that

- **The factorisation itself is checked and holds.**
  `verify_weight_factorisation.py` reproduces `legini`'s recurrence, builds the
  eight matrices the way `legini` builds them, and rebuilds each from the raw P
  and Q plus its factors: `qi` and `qj` exact, the other six within one ulp,
  which is the same products taken in a different order. It runs with a filter
  carrying an EXACT ZERO -- Cesaro at the top mode -- because that is the case
  that forbids recovering P from `qi` by division, and it is why the
  factorisation must keep P and Q raw rather than dividing the fold back out.
  Its negative control drops the zonal wavenumber from `qm`'s per-mode factor
  and is caught at 1.7.
- The order of operations changes, so the model's restart shas will move, as
  with everything else in this series. Same treatment: `compare_restarts.py` at
  a declared tolerance, against a control.
- Whether the extra per-mode vector read costs the inner loops their
  vectorisation. It is `NCSP` and shared across latitudes, so it should stay
  resident, but that is a prediction and the filter fold is the standing
  reminder that predictions about these loops have been wrong by a factor of
  five.
- Whether the gain survives at T21 and T42, where the matrices already fit
  anything. Expect nothing there, as with the symmetry work.


## CLIM-48 was built and it is a regression

*Measured 2026-08-21, T170 on sixteen, quiet machine.*

Storing P and Q and applying the per-mode and per-latitude factors, instead of
storing the eight products, takes the weights from 15.1 MB a thread to about
4.5 MB and from 114.9 MB a die to about 36, under CCD1's 32 MB of L3. It is
correct: rounding scale against the eight-matrix build at 1 and 20 steps, worst
4.8e-13 and 1.06e-11 against a 1e-10 bar declared first, with a control that
corrupts one factor and is rejected.

**And it is 1.56% SLOWER**, threaded build against threaded build, paired and
interleaved: 80.57 s against 81.66 s, spreads 1.3% and 2.3%, the two-matrix form
faster in 0 of 4 rounds, [-2.17, -0.74].

The premise was this note's own argument read backwards. The filter fold moved a
per-mode scalar INTO the matrices, saving one multiply in three, and returned
1.6%; that was read as evidence the loops are bound by STREAMING the matrices,
so cutting the matrices fourfold should pay. Undoing the fold costs 1.56%, which
is the same number back. The reading that fits both measurements is the simpler
one: **these loops are bound by the MULTIPLY, and the filter fold gained by
removing one.** Resident weight footprint is not what limits them, and a 4x cut
in it is worth nothing measurable here.

That also revises what the T170 model-code penalty was about. It was attributed
to weight pressure on CCD1; `mkdheat` then removed 106 MB of per-thread stack
and halved that penalty, and this change removes 79 MB a die of weights and does
not touch it. The pressure that mattered was the scratch, not the weights.
