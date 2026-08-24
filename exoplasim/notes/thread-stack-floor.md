# The per-thread stack the threaded build needs

Measured 2026-08-24 against `vendor/exoplasim/exoplasim/plasim/src` as it stands
and against two executables in `exoplasim/analysis/prerebuild_binaries_2026-08-18/`.

Worldbuilding frame: an instrument on the Vesper climate model's build. Nothing
here is about the simulated planet.

## Why the whole model is one stack frame chain

`plasim.f90` opens ONE `!$omp parallel` region and runs everything inside it:
`mpstart`, `prolog`, `master`, `epilog`, `mpstop`. It has to, because the model's
state is module data declared `threadprivate` and a thread's copy lives exactly
as long as the region does. So a thread's stack is not sized by one phase; it
has to hold the deepest call chain of the entire model.

`-fopenmp` implies `-frecursive`, which is what moves the local arrays gfortran
would otherwise place in static storage onto that stack. Those locals are the
term that matters: `exoplasim/notes/threads-instead-of-ranks.md` measured 463.1
MB of them in a T170 binary, against an OpenMP default stack of 8 MB.

That 463.1 MB is a sum over every routine and is therefore an unusable upper
bound -- the routines do not all run at once. The quantity that sizes the stack
is the heaviest CHAIN, not the heaviest total.

## The floor, and what it leaves out

`exoplasim/scripts/stack_floor.py` walks the model's own call graph out of the
parallel region's entry points and sums, over the routines on each path, the
declared local arrays whose bounds are constant at that rung. The heaviest path
is the deepest concurrent set of DECLARED locals.

At ten levels, sixteen threads and eight-byte reals, the heaviest chain is
`master -> gridpointd -> aero_main -> prepare_uvps` at every rung, and its
weight rises as NLAT squared -- exactly four times per doubling, because both
`NHOR = NLON * NLPP` and `NUGP = NLON * NLAT` carry that factor.

Run the script for the table; it is not restated here.

This is a FLOOR and nothing more. It counts declared locals. It does not count:

- the temporaries gfortran materialises for whole-array and `where` expressions,
  and the physics is written in array notation -- 169 `where` statements against
  about 22 explicit loops -- so these exist and are not small;
- saved registers, alignment padding, or frames inside SHTns and the FFT.

Nothing in the source fixes those. Only `-fstack-usage` on a real compile,
summed along the same chain, or a guard-page probe, does. That measurement is
`world-p4b` and this does not replace it.

## Why the parse can be trusted as a floor

The same declarations, summed over every routine instead of along one chain,
must reproduce what a compiled binary holds in static storage: without
`-frecursive` gfortran puts each of these locals in BSS as one symbol named
`<local>.<n>`. That is a check with a right answer, and `--validate` runs it.

- T170, sixteen ranks, eight-byte reals: the parse gives 464.4 MB against the
  463.1 MB `threads-instead-of-ranks.md` read off the binary. 0.3 per cent.
- T85, sixteen ranks, four-byte reals: 60.2 MB parsed against 57.4 MB in the
  binary's BSS over 617 symbols. 4.8 per cent.
- T42, same build: 15.1 MB parsed against 9.6 MB in the binary. 57 per cent, and
  the gap is one-signed and expected. gfortran spills to static storage above a
  size threshold, so at a small rung the smaller arrays stay in registers and on
  the stack and never appear as symbols. The parse is an upper bound on what BSS
  can show, and it converges where the arrays are large -- which is the end of
  the ladder the stack question is about.

## What this settles about `OMP_STACKSIZE`, and what it does not

`OMP_STACKSIZE` is 512M in `run_exoplasim.py` and in every gate and bench
script. The floor at the top of the ladder is a fraction of that, so 512M is
above the floor everywhere on the ladder -- which is a bound, not a derivation.
The margin between the floor and 512M covers the compiler temporaries, and
nothing here measures them, so the margin is not derived and the constant is
not either.

What has changed is that the constant is now CHECKED. `prepare_thread_stack`
computes the floor at the configured rung and refuses a stack below it, naming
the chain. At the OpenMP default of 8 MB a T170 run is refused; before this it
would have segfaulted inside the model, which reads as the model crashing.

Two things follow that are worth having on record.

**The requirement is per thread and only partly shrinks with the thread count.**
`NHOR`-shaped locals are a thread's band and divide by the team size;
`NUGP`-shaped ones are the whole globe and do not. At T170 the total of declared
locals falls by a factor of 4.7 between one thread and sixteen, not sixteen, so
adding threads is not a way out of a stack bound.

**The 32 MB per-die working set convention is not met by this model and is not
the constraint here.** That convention counts one copy per thread, so sixteen
threads have 2 MB each. The declared locals on the heaviest chain exceed that at
every rung from T42 up. This is not a regression to fix by shrinking the stack:
a stack reservation is address space, and what a chain touches in any one
routine is a few of its arrays, not all of them. The convention governs the hot
working set. The stack floor governs whether the process survives. They are
different quantities and only the second one is what `OMP_STACKSIZE` sets.

## The master thread is a separate answer

`OMP_STACKSIZE` sizes the non-master threads. The master runs on the process
stack, which is the stack rlimit, and it is a member of the team and runs the
same chain. `bench_ab.py` measured a T127 master overrunning a 16 MB limit and
segfaulting. `prepare_thread_stack` raises the soft limit to the hard one, which
on this host is unlimited; a host where the hard limit is finite and below the
floor is the case that still has no check, because the floor is known there and
the refusal is one line -- `world-p4b` again.
