# Sharing the spectral state

*This is a COMPUTE change to the Vesper worldbuilding project's climate model.
Nothing in it is about the simulated planet. In progress, 2026-08-21.*

## Why

`threads-instead-of-ranks.md` measured the threaded build at -18.9% at T127 and
-25.7% at T170, and found the cost: **memcpy, 22.6% to 41.0% of samples**, with
`mpgallsp` staging 105 GB and `mpsumsc` 68 GB over 300 steps. Both exist to give
a thread a private copy of something it could read in place, and both exist
because `mpimod_omp` keeps the MPI contract.

## What "reopening the physics" turned out to mean

The phrase was used loosely and hid the shape of the work. It means making the
`NHOR`-shaped GRID arrays full-globe, so that 169 `where` statements across
about 700 declarations in eight modules become slices. That is required only by
SHTns.

**It is not required here.** The staging traffic is entirely in the SPECTRAL
half, and that half is confined to the dynamical core: 217 references in
`plasim.f90`, 40 in the declarations, and **9 in all the physics modules
combined**. Zero of the 169 `where` statements touch a spectral array.

## The design

The partials stop being storage and become threadprivate POINTERS into the
shared array they are already a slice of:

    real, target  :: sd(NESP,NLEV)     ! shared
    real, pointer :: sdp(:,:)          ! threadprivate, => sd(lo:hi,:)

so all 269 reference sites are unchanged, and the gather that assembled `sd`
from every thread's `sdp` becomes a barrier. `OMPSHARED`, defined only by
`most_compiler_omp`, selects pointers over storage: under MPI a rank's partial
must be its own array, because `mpi_allgather` forbids a send buffer that
aliases its receive buffer.

## Three things this has already established

**`so` cannot be shared.** `fc2sp(foro,so)` has EVERY thread write the whole
array as its own partial before `mpsum` reduces it, so it is a per-thread
accumulator by construction. It is the only one of the seven that is a
transform output. Sharing it took the comparison from 140 to 156 of 199 records
identical.

**A partial must never cross a call boundary.** `sd(lo:hi,:)` is not
contiguous -- the stride between levels is `NESP` -- so passing it to an
explicit-shape dummy makes gfortran copy in and copy out. That copy is the
traffic this change removes, and worse, **the copy-out lands AFTER the callee's
barrier**, so a publish can release before a thread's slice is written back. The
publish therefore takes no arguments. The slice cannot be made contiguous
instead: the transform needs the full array mode-major within a level, and a
thread's modes are a sub-range of that.

**The enumeration method was wrong at first.** It covered ASSIGNMENTS to the
arrays becoming shared, but a write can happen through a CALL: `fc2sp(foro,so)`
writes all of `so` and looks like nothing in an assignment grep. The correct
enumeration is every routine taking a shared array as an output dummy, in both
the bare-name and the `ARR(1,j)` base-address form. Done properly, every
remaining call site is a read.

## Where the race is, and where it is not

It still races, and the determinism check is what catches it -- a check that was
never run against Stage 1, which turns out to BE deterministic, so the races are
new rather than inherited.

| | |
| --- | --- |
| one thread | **deterministic**, and agrees with the one-process build to 1.9e-13 |
| the pointers | correctly associated: thread 1 gets lo 254, hi 506, shape 253 x 10 |
| the publish | **sound** -- after `mpsyncsp` both threads read identical sums of sd, st, sz, sp |
| those sums | **differ run to run** |

So the shared-state logic is right, the barrier publishes what it should, and
the race is UPSTREAM of the publish: in what produces a thread's slice, not in
the publishing of it. The path to check is the one that feeds it -- the inverse
transform reads the shared arrays in `gridpointa`, the physics runs, `mktend`
produces a partial, `mpsumsc` reduces it, and `spectrala` writes the slice. A
barrier already separates `gridpointa`'s reads from `spectrala`'s writes, so
what remains is to find another reader of the shared arrays running
unsynchronised against those writes.

## Two tool facts that cost cycles

**ThreadSanitizer is unusable here.** It reported 94 races, but libgomp is not
instrumented, so it does not model `!$omp barrier` and flags every
barrier-synchronised access. The proof is that the deterministic Stage 1 build
produces the same flood.

**`compile.sh` copies sources with `cp -p`,** preserving mtimes, so an edited
file can leave a stale object and the build silently links the old one.
Deleting `plasim/bld/<file>.o` is the fix.

And one of my own: a probe must not use a hard-coded unit number. A Fortran unit
is process-global, and three separate probes in this work reported nonsense
because two threads opened the same one. One unit AND one filename per thread.
