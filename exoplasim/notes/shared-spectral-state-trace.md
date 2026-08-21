# Finding the shared-state race

*Appendix to `shared-spectral-state.md`. Measured 2026-08-21.*

Two races, both on arrays that Stage A made shared, found with an OpenMP-aware
race detector after hand tracing had bounded the search but not closed it.

## The races

**`dv2uv`, in `legmod.f90`.** It subtracted planetary vorticity from
`pz(1,2,v)`, ran the transform, and wrote the saved value back. That is
rank-local while every rank holds a private copy of the spectral state. It is
not thread-local once the state is shared: `pz` aliases the shared array, every
thread runs the routine over the whole of it, and subtract-transform-restore is
a write race on one element. Both threads subtract, so the mode loses planetary
vorticity twice, and the restores race with each other.

Planetary vorticity rides on the single mode `w = 2`, which is `m = 1, n = 2`,
so the transform being linear the term now goes to the RESULT instead of to the
input and the shared array is never written. In the paired branch that mode
falls in the antisymmetric loop, which gives the mirror latitude a sign the home
latitude does not have -- reading the home form for both is the mistake the
shape invites.

**The five spectral broadcasts in `prolog`.** `mpbcrn(sp,...)` and its four
siblings had every thread write the same shared words. Harmless in outcome, but
a race, and under `OMPSHARED` there is nothing to broadcast at all: root has
already written what the others are about to read. They now sit behind the same
barrier as the scatters.

## Why hand tracing did not find them

The trace bounded the race correctly and still missed it. `post-gpa` was
reproducible and `post-spa` was not, so the search went to `spectrala`; the
divergence point then moved between runs, which is characteristic of a race and
is also why each build says where it went wrong that time rather than where it
always goes wrong.

`dv2uv` was in the list of readers that had been accounted for. It reads the
spectral state and it also writes one element of it, and nothing in a reading of
the call graph distinguishes those. **A routine classified by what it is for is
not classified by what it touches**, and shared-state work needs the second.

## The detector

ThreadSanitizer needs the OpenMP runtime to publish its synchronisation, or it
reports every barrier-separated access as a race -- 94 of them here, which a
deterministic build reproduced exactly, proving them noise. Archer supplies
those annotations as an OMPT tool, and libgomp has no OMPT.

The route that works is **gfortran objects with `-fsanitize=thread`, linked
against LLVM's `libomp` instead of libgomp**. libomp carries the GNU
compatibility layer -- 256 `GOMP_*` symbols including `GOMP_barrier`, and the
Fortran-mangled `omp_get_thread_num_` -- so a gfortran build links and runs
against it unchanged, and it exports `ompt_start_tool` for Archer to attach to.

flang was the direct route and cannot do this: **flang 22 accepts no sanitizers
at all**, `-fsanitize=thread` and `-fsanitize=address` both rejected. It earned
its place anyway by rejecting a `!$omp threadprivate` directive on `ilatperm`, a
module procedure, which is invalid OpenMP that gfortran had accepted silently;
the directive was there because the list had been generated from `nm` symbols,
where a procedure reads the same as a variable. It also wants `omp_lib` declared
externally, since it ships seventeen modules and not that one.

Twenty races before, all of them on `sd`, `st`, `sz`, `sq` or `sp` and none of
them noise. None after, over 120 steps at T21 on two threads, which covers the
output path.

## The check that now stands behind it

`verify_shared_determinism.sh`, at T21 on two threads over 200 steps: three
runs, one restart hash. A different random seed gives a different hash, so the
hash tracks the computation rather than the harness.

The obvious control -- a deliberately racy build that must fail -- cannot be
built usefully. The pre-fix build does not survive five steps, because losing
planetary vorticity twice a step destroys the vorticity field and the model dies
in shortwave radiation. A control that crashes has bypassed the comparison
rather than failed it. That the pre-fix build cannot complete five steps while
the fixed one completes two hundred three times identically is the sharper
statement anyway.

## A trap that cost more than the race did

The bed. `verify_paired_decomposition.sh` was run on a cold-start bed with
`KICK = 1` and no `SEED`, and `plasim.f90` seeds that noise from the clock, so
one binary run twice did not agree with itself. It reported 94 of 199 records
wrong at ONE step against a correct rewrite. The tell was in the output and was
read past: the worst record was `seed` itself, at 4.3e+170 relative.

Seeded, the same comparison passes -- 0 records beyond rounding scale at 1 and
at 20 steps, worst 7.8e-11 against a 1e-10 tolerance -- which is what confirms
the plavor rewrite in both branches. `_bed_guard.sh` now refuses such a bed and
`docs/src/practice/failure-modes.md` class 27 is the general form.

Second trap, same afternoon: `compile.sh`'s `-j` is a FLAG and takes no
argument. `-j 2` leaves `ncpus` at its default of 4, builds a `p4` binary, and
leaves whatever `p2` was lying in `plasim/run` untouched -- so a check that
tests for the file's EXISTENCE runs somebody else's binary and reports its
crash as a result. Test freshness, not existence.
