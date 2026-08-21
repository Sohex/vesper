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
