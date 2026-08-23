# Is sixteen threads the right count at every resolution?

Measured on 2026-08-21 and 2026-08-22, on the SHTns build of the threaded
ExoPlaSim fork that supplies Vesper's simulated climate. The question is a
configuration one about the model executable, not about the simulated world.

CLIM-77 was opened on an inference rather than a measurement. At T42 libgomp is
a quarter of runtime, 25.5% with SHTns and 21.0% with legmod, against 10.9% at
T170; sixty-four latitudes across sixteen threads is four each, and the task
reasoned that the barriers were costing more than the work between them and
that the small rungs of the resolution ladder therefore wanted fewer threads.

**The inference is wrong, and it is wrong at exactly the resolution where it was
most likely to be right.** Sixteen is the correct count at T42, and nothing
below it comes close.

## What the decision variable was, fixed before the arms ran

Wall time, not scaling efficiency. The question is which configuration finishes
a bed first, not which uses its threads most fully; a count that leaves cores
idle is allowed to win if it wins on the clock. The 5% self-scatter floor
applies as everywhere else.

## Method

The thread count is COMPILED IN -- `!$omp parallel num_threads(NPRO)` at
`plasim.f90:59` -- so `OMP_NUM_THREADS` does nothing and each count needs its
own binary. `NPRO` must divide `NLAT`, which at T42 allows 2, 4, 8 and 16.

Each count is benched pairwise against sixteen with `bench_ab.py`, four
interleaved rounds, so every comparison is against a run seconds away on the
same machine. `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, `ulimit -s unlimited`,
one driver owning the machine. The driver is
`exoplasim/scripts/thread_count_sweep.sh`, which takes a resolution, a bed and the
counts to sweep, and builds each binary before any round runs so no compile
lands beside a measurement.

## Result: fewer threads is monotonically worse

| threads | wall | against sixteen | rounds won |
| ---: | ---: | ---: | ---: |
| 2 | 19.48 s | -308.7% [-310.3, -304.1] | 0/4 |
| 4 | 10.83 s | -123.7% [-125.8, -121.2] | 0/4 |
| 8 | 6.58 s | -36.0% [-39.2, -33.0] | 0/4 |
| 16 | 4.82 s | reference | -- |

The sixteen-thread median is 4.77, 4.85 and 4.83 s across the three arms, a 1.7%
spread and inside the floor, so the reference is stable and the comparisons are
against the same thing three times.

Two confounds both point the way that would have HELPED the low counts, and they
still lost. Fewer threads leave cores idle, so the active ones boost higher.
Eight threads under `close`/`cores` land on cores 0-7, which on this processor is
the 96 MB V-Cache die -- the best cache in the machine, to itself. Neither was
enough to buy back a third of the runtime.

## Above sixteen is already refused, on two independent grounds

Not re-measured here, because `smt-rank-layout.md` settled it: two threads per
core loses at every resolution, -113% at T42, -65% at T127 and -49.4% at T170,
and the penalty shrinks with the over-decomposition confound rather than
approaching parity. Separately, `NESP = ceil(NRSP/NPRO) * NPRO` is padded by the
thread count -- 1904 at sixteen against 1920 at thirty-two -- so a restart
written at sixteen threads runs past the end of the record at thirty-two.
Confirmed again incidentally here, as an I/O error in `get_restart_array` when
the thirty-two-thread arm tried to resume the warm T42 bed.

Eight-on-CCD0 against sixteen-across-both was also already measured at T170,
where sixteen wins by 14%. The T42 arm above is the same comparison at a second
resolution and agrees, by a wider margin: the V-Cache die buys less when the
footprint is small, so giving up half the cores for it costs more.

## Where the barrier time actually goes

Sixteen being right does not make a quarter of T42 runtime uninteresting, so the
barrier cost was attributed per thread rather than left as one number. From the
T170 frame-pointer profile, classifying a sample by its leaf: 83.8% model, 15.0%
barrier -- leaf in a stripped library with no model frame anywhere in the stack
-- and 1.1% kernel. Two runs, one taken with `--sample-cpu` and one without,
agree on all three to a tenth of a point.

Threads are pinned one to a core, verified rather than assumed: every thread
spends 100% of its samples on one core once both SMT siblings of a core are
counted as that core, which is what `OMP_PLACES=cores` promises.

| core | die | barrier |
| ---: | --- | ---: |
| 0-7 | CCD0, 96 MB | 16.7 to 18.1%, flat, mean 17.25% |
| 8-15 | CCD1, 32 MB | 8.8 to 16.5%, rising, mean 13.33% |

**The floor is 8.8% and the spread above it is imbalance.** The least-waiting
thread is on the critical path; what the others spend above it is time waiting.
Perfect rebalancing is therefore worth at most about six points of runtime, and
the rest is serial regions and barrier machinery, which no thread count and no
rebalancing reaches. This is a lower bound on waiting in both directions: the
event is cycles, and a thread that has stopped spinning and gone to sleep in a
futex burns none.

**The die accounts for 42% of the spread and the decomposition for the rest.**
T170 is NLAT 256 over sixteen threads, sixteen latitudes each and contiguous, so
the die boundary falls on the equator: threads 0-7 hold the northern hemisphere
on CCD0 and threads 8-15 the southern on CCD1. Cores 8 and 9 are both the
equatorial bands, where convection and insolation are heaviest, and they sit on
the 32 MB die. The two effects stack on the same two threads.

Which die is faster for this model is measured at T127 and NOT at T170: an
eight-rank job takes 46.19 s on CCD0 against 51.59 s on CCD1, so cache beats
clock by 10.5% there, and `smt-rank-layout.md` has no single-die arm at T170 to
say whether that carries. The barrier attribution above is T170, so the
direction of the die effect at T170 is inferred from the barrier shares
themselves rather than from a wall-time arm.

## Placement cannot collect it

That invites a free fix, since thread-to-core order is `OMP_PLACES` and not
code. Swapping the hemispheres between dies puts the critical-path bands on the
96 MB die at the cost of putting the northern bands on the 32 MB one. It also
tests the extrapolation above, since a placement swap can only pay if the die
difference measured at T127 survives to T170.

| | median |
| --- | ---: |
| hemispheres on native dies | 69.16 s |
| hemispheres swapped | 69.16 s |

**+0.006%, [-1.03, +3.66], faster in 2 of 4 rounds, and bit identical.** The
decomposition is untouched by a placement change, so the restart sha matching is
what proves the arm changed only what it claimed to.

The null is informative rather than disappointing, and it is the reason the
T127 die figure is not quoted here as a T170 one. Either the die difference
does not carry to T170, or it does and there are only eight fast cores for work
that does not fit in them, so moving the advantage relocates the bottleneck
instead of removing it. The arm does not separate those two, and does not need
to: both say placement is spent. **A static permutation cannot give the faster
die more WORK, and that is the only thing that would help.** A work queue can,
by construction: slower threads take fewer chunks and no part of the model has
to know the topology. That is an argument for the block decomposition in
CLIM-80 and the collapsed iteration space in CLIM-79, and it exists only
because the cheap shortcut was tried and failed.

The arm is `bench_ab.py --b-launch omp@<core order>`.

## Verdict

Sixteen threads at every resolution. The recommendation the task asked for is
per-resolution and it is uniform, which is the answer rather than a failure to
find one. A high barrier share is not evidence of over-threading: at T42 the
barriers cost a quarter of runtime and the work between them still parallelises
well enough that halving the team costs a third of the clock.
