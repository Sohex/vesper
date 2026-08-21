# Threads instead of ranks

*This is a COMPUTE change to the Vesper worldbuilding project's climate model.
Nothing in it is about the simulated planet. Started 2026-08-21.*

## Why

SHTns has no distributed mode -- upstream offers "parallel transforms with
OpenMP" and nothing else -- so a transposed MPI layout would hand a rank the
right data and leave the library single-threaded inside it. **SHTns's
parallelism is threads**, which makes this the prerequisite for it rather than
an alternative to it.

It also has a prize that does not depend on SHTns working at all. In one
address space every gather and scatter is a memcpy, and the critical-path rank
spends 12.0% of samples in MPI at T127 and 19.7% at T170
(`rank-imbalance-and-weight-traffic.md`). That floor is transfer and
synchronisation rather than waiting.

## The shape, and why the model above it does not change

`make_plasim` already selects the MPI layer through `${MPIMOD}`, and
`mpimod_stub.f90` is an existing second implementation of the same contract. So
this is a third variant in which the ranks become threads of one process. Four
things make that fit rather than being a rewrite:

- **`mypid` appears 345 times and every one is a comparison against `NROOT`.**
  Never arithmetic, never an index. Thread 0 is root and all 345 sites stand.
- The contract is **39 routines in one file**, and the thousand-plus call sites
  across the model all go through them.
- The physics is written in **whole-array notation** -- 169 `where` statements
  against about 22 `do jhor` loops -- so `!$omp parallel do` on loops would not
  work. It does not need to. Each thread owns an `NHOR` chunk exactly as a rank
  does, and the array syntax is then correct as written.
- Each thread owning a chunk is what `!$omp threadprivate` on the module state
  gives, **with no change to any physics file**.

## Stage 0: can gfortran carry the state per thread?

The one thing that could kill the approach before any of it was written. Under
MPI the model's state is per-process static data; under threads it goes to
thread-local storage, and there was no reason to assume gfortran would take it
at this size.

**What the binary actually holds**, T170 at 16 ranks:

| | size | count |
| --- | ---: | ---: |
| module variables | **150.3 MB** | 1198 symbols |
| locals gfortran spilled out of the stack | 463.1 MB | 724 symbols |
| total static data a rank | 642.3 MB | |

The second row was the surprise and it is good news twice over. Those are
ordinary local arrays -- `zgbl`, `zfpr3`, `zu`, `zv` -- that gfortran placed in
static storage rather than on the stack, and **`-fopenmp` implies `-frecursive`,
so they become stack-allocated and therefore per-thread automatically.** They
need no declaration at all. What has to be threadprivate is the 150 MB of
module state, not the 642 MB.

The cost of that is stack rather than TLS: the largest single spilled local is
14 MB, so the real model will need `OMP_STACKSIZE` set well above the 8 MB
default. That is a launch setting, not a code change, and it has to be recorded
where the run is launched or the model will die in a way that looks like a bug.

**The spike.** `plasimmod.f90` with `!$omp threadprivate` over all 476 of its
storage variables -- the list taken from the compiled binary's symbols rather
than by reading declarations, so nothing is missed -- built at T170 dimensions
and run at 16 threads with `copyin`. Each thread writes its own identity into
`dls`, `dt` and `gu`, barriers, and reads back what only it can have written.

    threads that ran: 16 of 16
    PASS: 476 threadprivate module variables, per-thread and isolated

`copyin` matters and is not decoration: a threadprivate variable's declaration
initialiser reaches only the master's copy, where a separate process gets it for
free. `copyin` at the parallel region's entry is what restores the equivalence,
and without it every non-master thread would start from zeroed memory --
silently, and with a plausible wrong climate rather than a crash.

Measured: **97.7 MB of TLS a thread** for `pumamod` alone, 1.6 GB peak resident
across the team, on a machine with 61 GB. Not a constraint.

**So Stage 0 is a pass and the approach is memory-feasible.**

## Stage 1: `mpimod_omp.f90`, and where it stands

Written, and the build carries it: `compile.sh -j` selects
`MPIMOD=mpimod_omp`, `UTILMOD=utilities_omp` and `most_compiler_omp`, and names
the executable `..._p16_omp.x` because a threaded binary is a different
identity from an MPI one at the same resolution, not a variant of it. The build
directory carries a marker for which parallel layer it last built, because the
three share object names and a stale one links silently.

**The MPI and serial builds are provably unchanged.** `plasimmod.f90` compiled
with and without the 476 `!$omp threadprivate` lines, same filename, same
flags, gives a BYTE IDENTICAL object: the directives are comments to a compiler
not given `-fopenmp`.

**The collectives are verified.** `verify_omp_collectives.sh` checks each of the
39 routines against the answer written down in advance -- not against another
run of the model. It passes at 2 and at 16 threads.

Its first version could not have failed. Every thread contributed a CONSTANT,
so the sums came out right no matter which slice a thread read: a uniform
contribution makes any index confusion cancel. The values now depend on the
mode index, which is the only way the check can catch the mistake this code
actually invites. That flaw shipped once before in this project, in the filter
fold's indexing test, and it went in again here.

**It reproduces the MPI model.** T21, against the MPI build on the same bed:

| | records | worst |
| --- | --- | ---: |
| 2 threads, 1 step | **199 of 199 BIT IDENTICAL** | -- |
| 16 threads, 1 step | 156 identical, 43 at rounding scale, 0 beyond | 1.0e-13 |
| 16 threads, 20 steps | 143 identical, 56 at rounding scale, 0 beyond | 3.5e-11 |

Bit identity at 2 threads was not expected and is worth understanding rather
than celebrating: `mpsumsc` sums the partials in thread order 0..NPRO-1, and
Open MPI's reduce_scatter happens to sum in rank order too, so at that size
nothing is regrouped. At 16 the groupings differ and the last bits move, which
is the expected result and the one the tolerance was declared for.

**And the MPI build is unchanged, proved behaviourally rather than by
inspection.** The committed source and the working source each built a 16-rank
T21 binary; the two produced a BIT IDENTICAL restart. That is the guarantee
that matters, and it is stronger than the object-level comparison it replaces.

## Two defects found on the way, one of them upstream's

**`hurricanemod` has never had a working root guard.** It has no
`implicit none` and imported only `NHOR, NLEV, NLEP, NUGP` from `pumamod`, so
in `hurricaneini` both `mypid` and `NROOT` were UNDEFINED implicit locals and
`if (mypid==NROOT)` compared two of them. Without `-frecursive` gfortran gives
such locals static storage, so both read 0, the test was true on every rank,
and every rank has been reading the namelist and writing the hurricane
diagnostics. It survived because all 21 variables in the group are broadcast
immediately after, so the duplicated read was overwritten with the root's
values. Under `-fopenmp` they become stack locals holding whatever was there,
several threads open unit 51 at once, and the read fails. Fixed by importing
the two names; behaviour-neutral for the MPI build, which is why the fix is
safe to make here.

**Five of this file's own restart routines had the wrong record length.** They
were written from `mpimod_stub.f90`'s shape rather than from `mpimod.f90`'s
semantics, and a spectral record holds `NRSP` values, not `NESP`. At 2 threads
`NESP == NRSP` and it did not show; at 16 it wrote a restart the model could
not read back. The lesson is narrow and worth keeping: the stub is a guide to
the INTERFACE and not to the behaviour, because at one process most of the
behaviour is absent.

## The third defect, and the one that took the longest

`glacierini` computes the spectral orography, and the two threads were sharing
the array it transforms: thread 1's input was thread 0's OUTPUT. The array is

    real :: foro(NHOR) = 0.0

declared inside the subroutine. **A local variable with an initialiser has
implicit SAVE**, so it lives in static storage -- one copy for the process, and
therefore one copy shared by every thread. `-frecursive` does not move it to
the stack, because the standard says it is SAVEd.

That is a whole class, not one variable: **47 initialised locals across 9
files**, every one of them shared where a rank had its own. `icemod` has 14,
`hurricanemod` 13. Each is now `!$omp threadprivate` in its own routine.

**Why it was missed.** Stage 0 split the binary's static data into module
variables and "locals gfortran spilled out of the stack", and concluded the
second group needed nothing because `-fopenmp` implies `-frecursive`. That is
true for spilled AUTOMATICS and false for implicitly SAVEd locals, and both
land in the same `.bss` with the same `name.N` symbol shape. The distinction is
invisible in the symbol table and lives only in the declaration: an initialiser
makes it SAVE.

## What the debugging cost, and what would have found it sooner

Four hypotheses were tested and refuted before the right one -- the collectives,
uninitialised stack locals, the surfmod orography branch, and other modules with
an undeclared root guard. Three things would have shortened it:

- **The bisect that worked was N_RUN_STEPS = 0.** Comparing after zero
  timesteps separated initialisation from the timestep in one run, and the
  diagnostic output then named the exact line where the two builds first
  disagree. That should have been the first move, not the fifth.
- **A probe must not use a hard-coded unit.** The first instrumented run had
  both threads opening unit 77, which is process-global, and they clobbered
  each other -- so the probe reported that the threads held each other's
  latitudes, which was false. One unit per thread fixed it. The tool being
  debugged with was subject to the same hazard as the code under test.
- **`paste` of two files with different line counts.** A misread of that output
  produced a confident and wrong statement that the latitude chunks were
  swapped. They were identical.

## libgomp does not bind, and that has to be fixed before anything is measured

Open MPI gives rank *r* core *r*, deterministically; this project verified that
rather than assuming it. **libgomp's default is UNBOUND**, and it is worse than
merely different: the placement changes run to run, and with sixteen threads
some land on SMT siblings of the same physical core while whole cores sit idle.
On this processor it also randomises which DIE a thread gets, and
`rank-imbalance-and-weight-traffic.md` measured that the two dies are not
interchangeable at T127 and above.

Measured, three runs each:

    default                        cpus 0,1,2,13,3,10,6,14,21,7,15,16,11,12,23,29
                                   and a different set every run
    OMP_PROC_BIND=close
    OMP_PLACES=cores               cores 0..15, one thread each, every run
    Open MPI, for comparison       cores 0..15, one rank each

So the generated run script exports both, and a threaded measurement without
them is not a measurement of the threading.

## What Stage 1's measurement has to carry

The rank matrix's verdicts do NOT transfer wholesale, and one of them
specifically does not.

**SMT was refused on a mechanism that is gone.** 32 ranks lost 49% to 113%
because Open MPI busy-waits in `opal_progress`: two ranks sharing a physical
core each burn the other's cycles. That is a property of the runtime, not of
the model, and OpenMP's barrier wait is configurable -- `OMP_WAIT_POLICY` and
`GOMP_SPINCOUNT` make a waiting thread sleep instead of spin. The verdict has
to be re-measured rather than inherited.

The same fact bears on where a gain could come from at all. CCD0's ranks spend
23% of samples waiting while CCD1 paces them; under MPI that time is burned
spinning and is unrecoverable. Under OpenMP with passive waiting those cores
are idle instead. On a dedicated machine that buys nothing for LATENCY, because
there is nothing else to run -- but it is exactly the difference that should
show in the two cases where something else can use the cycles.

| arm | what it settles |
| --- | --- |
| 16 threads against 16 ranks, bound | what Stage 1 is worth |
| `OMP_WAIT_POLICY` active against passive at 16 | governs the arms below; likely null for latency, and two runs to know |
| 32 threads | whether the SMT refusal survives without the spin |
| two concurrent 8-thread jobs | the throughput incumbent, where passive waiting should matter most |

Two things to carry INTO the 32-thread arm rather than discover after it:

- **The restart incompatibility is unchanged.** A 16-way restart is unreadable
  at 32 because `NESP = NSPP * NPRO` gives 1904 against 1920. That is
  arithmetic, not an MPI property, so a favourable SMT result does not make 32
  adoptable on its own.
- **At 32 the paired layout drops out at low resolution.** `LPAIRLAT` needs
  `NPRO` to divide `NLAT/2`, and at T21 `mod(32,64)` is not zero, so 32 threads
  falls back to contiguous latitudes and silently gives up the symmetry work.
  T127 and T170 divide cleanly. A T21 arm at 32 would compare two different
  algorithms and must be stated as such or dropped.

## What Stage 1 measured: it LOSES, by 19% and 26%

Bound both sides, four interleaved rounds, order flipped, one driver owning the
machine.

| | ranks | threads | gain | rounds |
| --- | ---: | ---: | ---: | ---: |
| T127, 16 against 16 | 37.02 s | 44.01 s | **-18.9%** [-19.2, -17.9] | 0/4 |
| T170, 16 against 16 | 75.22 s | 94.71 s | **-25.7%** [-27.4, -24.5] | 0/4 |
| T127, wait policy | active | passive | -1.7% [-2.3, +0.1] | null |
| T127, 16 against 32 threads | 44.26 s | 68.74 s | **-55.2%** | 0/4 |

**The SMT refusal survives the change of runtime.** 32 ranks lost 65.1% and 32
threads lose 55.2%, so taking away Open MPI's busy-wait recovered about ten
points and left the verdict intact. Oversubscription was never mainly
`opal_progress` spinning; it is the decomposition. That question was worth
re-asking because the old answer rested on a mechanism this change removes, and
now it is measured rather than inherited.

**The wait policy is null for latency**, as it has to be on a dedicated machine:
there is nothing else to run, so a thread that sleeps at a barrier only adds
wake-up latency. It would matter for throughput, which is not measured here
because the latency result already decides the question.

Two arms were cut before they ran -- T170's wait policy and SMT. Neither could
change a decision: the policy was null at T127 and only ever fed the SMT arm,
and SMT was refused at T127 by 55 points. The T170 headline arm was allowed to
finish its fourth round rather than being stopped at three, because round four
is the second B-first round and stopping early would leave exactly the order
bias the interleaving exists to cancel.

## Where it goes, and it is not the barriers

perf on the threaded build against an MPI rank of the same bed:

| | MPI rank | threaded |
| --- | ---: | ---: |
| model code | 58.2% | 43.6% |
| **libc, all of it memcpy** | **22.6%** | **41.0%** |
| parallel runtime | 12.5% (libopen-pal + libmpi) | 10.7% (libgomp) |
| libm | 5.8% | 4.4% |

**The barriers are not the problem.** libgomp costs 10.7% against Open MPI's
12.5% -- if anything slightly cheaper. The whole deficit is memcpy, +18.5
points, which is very nearly the 18.9% measured at T127.

Every hot libc address resolves to `memcpy`, and neither dwarf nor frame-pointer
unwinding will attribute them: libc is stripped and the unwind stops there. So
the attribution was COUNTED rather than sampled, by instrumenting each routine
with the bytes it copies. Over 300 steps at T127:

| routine | staged | share |
| --- | ---: | ---: |
| `mpgallsp` | 105 GB | 56.1% |
| `mpsumsc` | 68 GB | 36.2% |
| `mpgagp` | 8.4 GB | 4.5% |
| `mpscgp` | 5.5 GB | 2.9% |

**183 GB over 300 steps: 624 MB a timestep, 4.2 GB/s of pure staging copy**, and
two routines are 92% of it.

## Why that is a verdict on the design and not a tuning problem

`mpgallsp` gathers every thread's `NSPP` slice so that every thread holds the
full `NESP` array. `mpsumsc` stages every thread's full partial so the slices
can be summed. **In one address space both are unnecessary.** The threads could
read each other's slices out of a single shared spectral array and reduce into
it in place, with no copy at all.

The copies exist because this variant keeps the MPI CONTRACT: each thread owns
a private `pp` and a private `pf`, exactly as a rank owns them. That decision is
what let the physics stay untouched -- 169 whole-array `where` statements
address a thread's chunk as written -- and it is the same decision that costs
19% at T127 and 26% at T170. **The deficit is the price of emulating private
memory inside shared memory.**

So there is no cheap fix. A threaded build that WINS has to make the spectral
state genuinely shared rather than threadprivate, which reopens the boundary the
design closed and puts the physics back in scope. That is a different and much
larger change than this one, and it should be decided on its own terms rather
than entered by momentum from here.

**And it is the plan's declared stopping condition.** SHTns needs threads, its
own ceiling is 10 to 15%, and it cannot pay off a 26% hole. Stage 2 is not
reachable from this footing.

## What is left behind, and what it is worth

The threading is CORRECT and stays: bit identical to MPI at two threads, at
rounding scale at sixteen, verified by a collectives unit test and a
model-level comparison, with the MPI build proved unchanged. It costs nothing
to keep -- `${MPIMOD}` selects it and nothing else compiles it -- and it is the
starting point for anyone who takes up the shared-state version.

Three defects it found are worth more than the performance result, and they are
all in code that predates this work:

- `hurricanemod`'s root guard has never been one, on any build.
- 47 initialised locals across 9 files are implicitly SAVE. Under MPI that is
  per-process and harmless; it is the reason a threaded build cannot simply
  reuse them, and it is a portability hazard for anyone else who threads this
  model.
- libgomp does not bind by default, and on this processor an unbound team
  randomises which die a thread lands on.

