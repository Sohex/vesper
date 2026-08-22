# Is sixteen threads the right count at every resolution?

Measured on 2026-08-21 and 2026-08-22, on the SHTns build of the threaded
ExoPlaSim fork that supplies Vesper's simulated climate. The question is a
configuration one about the model executable, not about the simulated world.

CLIM-65 was opened on an inference rather than a measurement. At T42 libgomp is
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
fully-unwound T170 profile, classifying a sample as barrier when its leaf is in a
stripped library and no model frame appears anywhere in the stack:

| | share of that thread's cycles |
| --- | ---: |
| least-waiting thread | 9.0% |
| most-waiting thread | 18.2% |
| team mean | 15.3% |

**The floor is 9.0% and the spread above it is imbalance.** The least-waiting
thread is on the critical path; what the others spend above it is time waiting
for work that could have been theirs. Perfect rebalancing is therefore worth at
most about six points of runtime, and the remaining nine are serial regions and
barrier machinery, which no thread count and no rebalancing reaches.

This is a lower bound on waiting in both directions: the event is cycles, and a
thread that has given up spinning and gone to sleep in a futex burns none, so
sleep is invisible to it.

The lever the number points at is the LEVEL decomposition, and that is CLIM-67's
territory rather than a second task: the SHTns wrappers are parallel over levels,
NLEV is 10 against sixteen threads, and a static schedule over ten items across
sixteen threads cannot be balanced. The wrappers already chain `nowait` so work
flows between loops, which is why the imbalance is six points rather than the
forty a single unchained loop would give.

## Verdict

Sixteen threads at every resolution. The recommendation the task asked for is
per-resolution and it is uniform, which is the answer rather than a failure to
find one. A high barrier share is not evidence of over-threading: at T42 the
barriers cost a quarter of runtime and the work between them still parallelises
well enough that halving the team costs a third of the clock.
