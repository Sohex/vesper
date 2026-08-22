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

## Status: correct, and it recovers about half the deficit

Stage A is implemented, deterministic, and measured.

**Determinism.** At T21 on two threads over 200 steps, three runs give one
restart hash, and a different random seed gives a different hash, so the
agreement is about the model rather than about the harness.
`verify_shared_determinism.sh` is that check; `shared-spectral-state-trace.md`
has the two races it took to get there.

**Numerics.** Against the MPI build at T127 on sixteen, the threaded build
agrees at the scale of a regrouped sum, which is the standard for a change that
reassociates: worst 1.5e-13 at one step, 1.5e-11 at five, 1.4e-11 at twenty.
The two builds do NOT agree bit for bit at sixteen and are not expected to --
Open MPI's reduce-scatter orders the partials differently from `mpsumsc`'s fixed
loop. That seed amplifies, to 4.0e-10 by sixty steps and 3.3e-06 by three
hundred, which is why `bench_ab.py` reports "numerics: CHANGED" on a 300-step
bed and why that report is a statement about Lyapunov growth rather than about
correctness.

**Speed**, sixteen threads against sixteen ranks, four interleaved rounds, order
flipped, one driver owning the machine:

| | Stage 1 | Stage A | recovered |
| --- | ---: | ---: | ---: |
| T127 | -18.9% | **-11.22%** [-11.64, -8.69] | 7.7 points |
| T170 | -25.7% | **-20.66%** [-20.84, -20.17] | 5.0 points |

Self-scatter 1.6% and 0.7% on the threaded arm, against a 5% floor. Threads are
faster in 0 of 4 rounds at both resolutions -- that is, slower in every round --
so the remaining gap is real and not an artifact of taking medians.

**The prediction written down before the run was wrong.** It said the threaded
build would land level with MPI plus or minus a few percent. It did not; it
halved a deficit it was supposed to close. What was missing from the prediction
was that `mpgallsp` is a little over half of the staging and Stage B has the
rest, which the plan said plainly and the prediction did not use.

## Where the remaining time goes

A flat profile by shared object, which is attributable where a symbol-level one
is not, because libc is stripped. Shares converted to seconds against the
measured medians:

| T127 | ranks 36.28 s | threads 40.29 s | delta |
| --- | ---: | ---: | ---: |
| model code | 18.22 | 18.26 | +0.04 |
| libc, which is the staging copy | 7.01 | 15.29 | **+8.28** |
| MPI stack against libgomp | 8.51 | 4.48 | -4.03 |
| libm | 2.16 | 2.06 | -0.10 |
| | | | +4.19 against +4.01 measured |

| T170 | ranks 73.78 s | threads 88.90 s | delta |
| --- | ---: | ---: | ---: |
| model code | 37.56 | 44.13 | **+6.57** |
| libc | 11.25 | 32.83 | **+21.58** |
| MPI stack against libgomp | 20.74 | 7.77 | -12.97 |
| libm | 3.72 | 3.82 | +0.10 |
| | | | +15.28 against +15.12 measured |

Both close to within 5% of the measured difference, which is what makes the
decomposition worth reading rather than merely suggestive.

Three things follow.

**The diagnosis holds.** libc is still the largest term on both sides and it is
the staging copy. `mpsumsc` is 68 GB of the 82 GB left after Stage A, so Stage B
attacks 83% of what remains.

**Dropping MPI is worth real time**, 4.03 s at T127 and 12.97 s at T170, and
that is the term that pays if the staging goes. It is why the arithmetic can end
up FAVOURING threads rather than merely reaching parity.

**T170 has a fourth term that Stage B will not touch.** The model's own code is
6.57 s slower under threads there, 17%, while at T127 it is identical to within
0.04 s. Sixteen threads share one address space where sixteen ranks do not, and
at T170 the weight matrices are 114.9 MB a die against CCD1's 32 MB of L3 --
which is CLIM-48's finding, arrived at from the other direction. That term is a
reason to expect T170 to trail T127 after Stage B, and a reason CLIM-48 and this
work are the same problem seen twice.

**What Stage B is projected to be worth, stated before it is built:** if the
libc excess falls with the staged bytes, 83% of it goes, which is -6.9 s at T127
and -17.9 s at T170. That would put threads at roughly 33.4 s against 36.28 s
and 71.0 s against 73.78 s -- ahead at both, by about 8% and 4%. The projection
assumes the copy cost is proportional to bytes staged and that the partials
really can be written in place, and the honest form of it is that Stage B should
move the deficit by more than Stage A did and may cross zero. Given this note
already records one confident prediction that was wrong by 11 points, treat the
sign as the claim and the magnitude as a guess.
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

## mkdheat, which is where the plan did not point

The plan's Stage C was `so`, `sr`, `sak`, `sqout` and the accumulators, described
as small traffic to be done for consistency rather than for speed. That is still
true of those arrays and it is no longer the next thing.

**`mkdheat` is, and it was missed by reading the wrong guard.** It is the
frictional heating term, `spectrald` calls it at `plasim.f90:3974` under
`ndheat > 0`, and `ndheat` is 1 by default -- so it runs every timestep. The
three `mpsumsc` calls inside it were left on the old staging path in Stage B on
the grounds that `nenergy` and `nentropy` are 0; that guard governs a block
LOWER in the routine, and the reductions are above it and unconditional. A
frame-pointer profile put the routine at 12.29% of all samples, second only to
`radstep`.

It carries both halves of the problem at once:

| in `mkdheat`, per timestep | |
| --- | --- |
| `mpgallsp` | 8 calls, each staging about 2.5 MB a thread at T170 |
| `mpsumsc` | 3 calls, still copying partials into a buffer |
| stack locals | **18.7 MB a thread**, 299 MB across sixteen |

Six of those locals are `(NESP,NLEV)` FULL spectral arrays -- `zsd`, `zsz`,
`zsq`, `zsdef`, `zstf1`, `zstf2` -- one complete copy of the global field per
thread. That is the pattern Stages A and B removed everywhere else, still
standing here. One thread's 18.7 MB does not fit CCD1's 32 MB of L3 beside the
weight matrices, which are 114.9 MB a die at T170.

So the traffic argument and the cache argument name the same routine, and the
fix is the machinery that already exists: shared arrays, threadprivate slices,
`mpsumscp`, gathers that become barriers.

**Done, CLIM-50.** The six full spectral arrays became four shared ones and
three partial slots in `pumamod`; what stays on the stack is what is genuinely
per-process, the grid fields covering this process's latitudes and the `NSPP`
partials. Two collectives carry it, both build-agnostic the way `mpsumscp` is.
`mpgallspp` is `mpgallsp` with the staging removed -- the destination is one
array the team can see, so a thread writing its own slice IS the gather.
`mpzerosp` clears a replicated spectral array from a given mode upward, slice by
slice under threads and as the plain statement it replaces under MPI.

Two places the shared-array discipline bit, both caught reading the code rather
than by the detector. `zhd(:,:) = zhd(:,:)*ct*ww` in the `nenergy` block would
have been every thread writing a shared array, so the scaling moved onto the
small partial before the gather. And clearing modes 2 upward could not stay a
whole-array statement, which is what `mpzerosp` exists for.

**It is numerically INERT**, and that is the bar this change earns by
reassociating nothing: the T21 two-thread restart hashes `53b3316c9fa35f59`,
the value Stage B produced on the same bed at the same length before `mkdheat`
was touched. Bit identical against MPI at 60 steps as well, 199 of 199 records,
though that arm alone would not have shown inertness -- both builds changed
together and could have agreed with each other while both moved. Archer over
120 steps: no races. Serial, MPI and threaded all build.

**Not yet measured.** The first attempt is void: the ranks arm scattered 10.3%
against a 5% floor and its median moved from 78.04 s to 85.33 s between runs on
a binary whose only change was swapping identical implementations, because
unrelated work was running on the machine. The threaded arm in the same run was
clean at 2.8% and 84.97 s. One clean arm is not a paired result and the number
is not quotable; it is re-run when the machine is quiet.

## Two hypotheses that were wrong, and what refused them

Kept because each cost a cycle and would otherwise be re-proposed.

**The remaining `mpgallsp` calls in `spectrala` and `spectrald`.** All of them
sit behind `ndiagsp == 1` or `nenergy > 0 .or. nentropy > 0`, all off. Cold.

**mmap churn from Fortran array temporaries.** The physics creates about 250 of
them and at T170 an `(NHOR,NLEV)` temporary is 655 KB, over glibc's 128 KB
threshold, so each is an mmap/munmap pair serialising on an address space
sixteen threads share and sixteen ranks do not. Plausible and false: raising
`MALLOC_MMAP_THRESHOLD_` moved the threaded arm -1.1% and -0.5% against a 3%
bar declared before the run, and moved the MPI arm the wrong way.

Both were guesses at an attribution that a stripped libc will not give up. What
answered it was building the model with `-fno-omit-frame-pointer` and letting
perf walk out of the copy into the caller. That is the tool to reach for first
next time, not third.
