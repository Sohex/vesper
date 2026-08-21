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

That enumeration missed one anyway, and the miss is the interesting part.
`dv2uv` takes the shared vorticity array as an INPUT and writes one element of
it, temporarily, to remove planetary vorticity before the transform. No
classification by what a routine is FOR would catch that. `shared-spectral-
state-trace.md` has both races, how they were found, and the detector that
found them.

## Status

Stage A is implemented and deterministic. At T21 on two threads over 200 steps,
three runs give one restart hash, and a different random seed gives a different
hash, so the agreement is about the model rather than about the harness.
`verify_shared_determinism.sh` is that check.

The plavor rewrite `dv2uv` needed is common to both builds, so the threaded
build's agreement with the MPI build is unaffected by it and the paired and
contiguous layouts still agree at rounding scale -- 0 records beyond it at 1 and
at 20 steps, worst 7.8e-11 against a 1e-10 tolerance.

**What has NOT been done is the measurement.** The whole case for Stage A is
that removing 105 GB of `mpgallsp` staging closes most of a 19 to 26 point
deficit, and that is a prediction until T127 and T170 are run against the
registered MPI binaries. The plan's stopping condition stands: if Stage A does
not move the memcpy share, the diagnosis is wrong and Stage B is not worth
building.

## Tool facts

**ThreadSanitizer works here, against LLVM's `libomp` rather than libgomp.**
libgomp publishes no OMPT, so TSan cannot see `!$omp barrier` and flags every
barrier-separated access -- 94 of them, proved noise by a deterministic build
reproducing the flood exactly. Archer supplies the annotations as an OMPT tool.
Build gfortran objects with `-fsanitize=thread` and link them against
`/usr/lib/libomp.so`, which carries the GNU compatibility layer, then run under
`OMP_TOOL_LIBRARIES=/usr/lib/libarcher.so`. flang cannot stand in: flang 22
accepts no sanitizers at all.

**`compile.sh`'s `-j` is a flag and takes no argument.** The thread count is
`-n`. `-j 2` leaves ncpus at the default of 4 and builds a `p4` binary, leaving
any stale `p2` in `plasim/run` in place for a careless check to pick up. Check
that a binary is NEWER than the build that was supposed to make it, never that
it exists.

**`compile.sh` copies sources with `cp -p`,** preserving mtimes, so an edited
file can leave a stale object and the build silently links the old one.
Deleting `plasim/bld/<file>.o` is the fix.

And one of my own: a probe must not use a hard-coded unit number. A Fortran unit
is process-global, and three separate probes in this work reported nonsense
because two threads opened the same one. One unit AND one filename per thread.
