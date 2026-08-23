# AOCL has nothing to attach to, and the model's build flags are not the lever either

Measured 2026-08-20 on the host below, after `aocl-gcc` 5.3.0 was installed.

Worldbuilding frame: this file is about how the Vesper simulation's climate
model is COMPILED and how fast it integrates. Nothing here is about the
simulated world.

## The question

AMD Optimizing CPU Libraries were installed on the workstation. Do the builds
this project makes get anything from them?

## What the builds actually link

`most_plasim_t42_l10_p16.x` resolves, per `ldd`, glibc `libm` and `libmvec`,
`libgfortran`, and Open MPI. `nm -D --undefined-only` finds no `gemm`, no
`gemv`, no `getrf`, no `fftw` symbol in either the model or `burn7.x`. The
link line in `plasim/src/make_plasim` carries no `-l` at all.

That is not an accident of configuration. PlaSim carries its own spectral
transform: `fftmod.f90` for the FFT and `legmod.f90` plus `legfast32` for the
Legendre transform. There is no dense linear algebra in the dynamical core to
hand to a BLAS.

So AOCL's BLIS, libFLAME, AMD FFTW, ScaLAPACK, sparse and RNG have no
attachment point in this project. Neither has anything else here: Orogen is
JavaScript, LPJ-GUESS is a C++ vegetation model with no BLAS dependency, and
numpy/scipy carry a bundled OpenBLAS that already dispatches an AVX-512 kernel
on this host (`Core: SkylakeX` under `OPENBLAS_VERBOSE=2`) for work that in
this project is norms and a single `lstsq`.

Two AOCL components DO have an attachment point, because the model links
`libm`: AMD LibM (`libalm.so`) and AOCL-LibMem (`libaocl-libmem.so`). Both were
measured rather than argued about.

## The bed

A copy of `run_4182235e9781`'s working directory, T42 L10 p16 off its restart,
`N_RUN_STEPS` cut to 600 and `NOUTPUT`/`NSNAPSHOT` to 0 so the instrument reads
model compute and not postprocessing. Wall time of the `mpiexec` call, three
repeats, median reported. Host: Ryzen 9 7950X3D, 16 cores, 16 ranks.

The thresholds were fixed before any arm ran, and are
`score_rank_bench.py`'s rather than new numbers: under 5% median improvement is
no difference. Added to them: any arm whose `plasim_status` sha256 differs from
the baseline's is a NUMERICS change and is reported as one whatever its time.

## The measurements

| arm | median s | vs baseline | restart sha256 |
| --- | --- | --- | --- |
| A baseline, the shipped binary | 12.575 | -- | 4497a8d73f61db9e |
| CTRL same source rebuilt by hand, stock flags | 12.604 | +0.2% | 4497a8d73f61db9e |
| B A plus AMD LibM preloaded | 14.266 | 13.4% SLOWER | f9e12701ae3b45b8 |
| C `-fcheck=all` removed | 11.943 | 5.0% faster | 82bc2b3d0d49499b |
| D A plus AOCL-LibMem preloaded | 17.733 | 41.0% SLOWER | 4497a8d73f61db9e |

CTRL is what makes arm C readable, and it reproduces more than its timing. A
hand build from the same tree with the stock flags produces an executable
byte-identical to the shipped one, sha256 `1695e1aefc51c97f`, which is also
what `binary_manifest.json` records for it. **The model build is bit-reproducible
on this host**: the same sources through the same flags give the same binary, so
the sha in the manifest identifies the pair and not just the sources. The
difference arm C shows is therefore the flag, and nothing about the build
environment.

## What each arm means

**AMD LibM is a loss twice over.** It is slower than glibc's `libm` plus
`libmvec` here, and it changes the integration: the restart diverges. Whatever
the second cost is worth, it is not being bought with speed.

**AOCL-LibMem is a large loss.** Its restart is unchanged, which is what a
`memcpy` replacement should do, and it costs 41%. The model moves many small
arrays; the dispatch overhead is not amortised.

**`-fcheck=all` costs 5.0% AT T42, and the reason recorded here for keeping it
was wrong.** Removing gfortran's runtime bounds checking buys 5.0% at this
resolution, which lands on the line that says "no difference", and it changes
the restart, because the flag inhibits optimisations whose absence changes
floating-point association.

This note previously argued that the check was earning its keep because
`oceanmod.f90:236-238` compiles with `Array reference at (1) out of bounds
(2 > 1)`, "so the checking is not guarding a hypothetical". Both halves of that
are false, established 2026-08-20:

- The diagnostic is a COMPILE-TIME warning from `-O3`'s static array-bounds
  analysis, and it appears identically with and without `-fcheck=all` -- three
  warnings either way. It is not the runtime checking the flag installs, so it
  is not evidence for the flag.
- The reference is unreachable. `NLEV_OCE` is a parameter equal to 1, the block
  is guarded by `if(NLEV_OCE > 1)`, and the loop is `do jlev=1,nlem_oce` with
  `nlem_oce = NLEV_OCE - 1 = 0` and never assigned anywhere in the file. The
  body is correct for any `NLEV_OCE > 1`; gfortran simply does not use the
  guard to prune it before bounds analysis.

So it WAS a hypothetical, and the check was never guarding it. Declaring
`nlem_oce` a `parameter` -- which is what it is -- lets gfortran prove the loop
empty, and the warning goes. `plasim.f90:1629` was the same shape, a
`neqsig==5` branch guarded by `NLEV > 10` whose body indexes `sigmah(0)` when
`NLEV = 10`; `max(NLEV-10,1)` on the loop bound is the idiom the surrounding
code already uses two lines below. With both, the model builds with ZERO
warnings, and the restart sha is unchanged by either.

The decision on `-fcheck=all` is therefore reopened on its real terms and is
NOT settled here: runtime bounds checking is worth something on its own account,
and against that it costs 5% at T42 and 17% at T127
(`exoplasim/notes/spectral-transform-profile.md`). That is a cost/benefit call
with no defect on the scale to break the tie.

### Round three: the flags that survived round two, measured against the
### production line, 2026-08-20

Round two priced `-funroll-loops` and `-fprefetch-loop-arrays` at roughly 1.8%
and 1% -- but those arms ALSO dropped `-fcheck=all`, so the marginal
contribution of each flag was confounded with the removal. Measured directly
against the current production line, which already has the bounds checking out,
both are refused:

| arm | paired gain vs production | rounds won | numerics |
| --- | ---: | ---: | --- |
| `-funroll-loops` | -0.15% | 3/6 | CHANGED |
| `-fprefetch-loop-arrays` | -4.75% | 1/6 | unchanged |

Both arms carry self-scatter above the 5% floor, because a concurrent job
arrived partway through, so neither number is quotable AS A NUMBER. They are
quotable as a REFUSAL: the decision threshold is a gain above 5%, both medians
are at or below zero, and noise that is symmetric between interleaved,
order-flipped arms does not turn a real gain into a negative median. A screening
test tolerates noise the way an estimate does not.

`-fprefetch-loop-arrays` being the worse of the two is consistent with the
cache-miss profile: the Legendre routines miss least of anything in the model,
so the hardware prefetcher is already doing this job and a software hint only
adds instructions. It also comes back numerics-unchanged, which confirms round
two's "CHANGED" verdicts on these arms belonged to `-fcheck=all` and not to the
flags themselves.

So the production flag line is closed for now. What is left in codegen is
bounded by round two's other finding: the entire distance from scalar code to
AVX-512 is worth 2 to 3% here.

**Round three's refusals were re-taken on 2026-08-22 and one of them reverses.**
Everything above was measured on the MPI build with legmod as the transform;
SHTns has since taken every per-timestep transform site. `-funroll-loops` is
**+11.10%** on the threaded T170 build, 4 of 4 rounds, and is adopted --
what it unrolls is now the physics rather than the Legendre sums this round
measured. `-fprefetch-loop-arrays` and `-fstack-arrays` keep their refusals, and
`-flto` is refused again on its own timing rather than on the `zsolars` defect,
which CLIM-37 and CLIM-38 fixed. `notes/audits/model-build-flags.md` has all
four, and the reading that matters for this document: a refusal is only as good
as the profile it was taken on.

The model is compute-bound inside its own Fortran, which is where a spectral
GCM should be bound. There is no library seam to widen.

## Two build-system defects found while measuring

Both were closed the day they were found, CONS-11 and CONS-10; the evidence
below is what they were found to be, and `archive/tasks.md` carries what was
done. The dependency defect turned out to be three edges rather than the one
described here.

A CLEAN build of the model failed. `glaciermod.o`'s rule in `make_plasim` named
`plasimmod.o` but not `landmod.o`, while `outmod.o`, which does depend on
`glaciermod.o`, sits EARLIER in `OBJ` than `landmod.o` does. `compile.sh` did not
empty `plasim/bld` except when switching between the MPI and serial builds,
so a `landmod.mod` from some previous build was always lying there and the
missing edge never showed. Both halves of that are fixed: the rule names
`${LANDMOD}.o` and `compile.sh` now empties `plasim/bld` on every build. What
the fix left is measured below.

`binary_manifest.json` records the sha of every source file that goes into a
build and `exoplasim_version`, and records nothing about the compiler or the
flags. A toolchain change therefore produces a different executable under an
identical set of recorded sources. `--verify` still catches it, since it
compares the executable's own sha, but it reports it as unknown provenance and
cannot say that the toolchain is why. Installing a new compiler or a new libm
is exactly the change that lands here.

## Is `make -j` safe now, measured 2026-08-21

The graph was re-derived from the source rather than read: for every file, the
modules it `use`s, each resolved to the file defining it, checked against the
prerequisites `make_plasim` declares. Run for all three configurations
`compile.sh` can produce, since the module variables change which files are in
`OBJ`.

| configuration | OBJ members | members with no rule | missing edges |
| --- | ---: | ---: | ---: |
| serial, `most_compiler` | 35 | 0 | 0 |
| MPI, `most_compiler_mpi` | 35 | 0 | 0 |
| OpenMP, `most_compiler_omp` | 35 | 0 | **1** |

The one edge is `utilities_omp.o`, which `use`s module `mpiomp` at lines 85 and
106. `mpiomp` is defined in `mpimod_omp.f90`, and the rule names `plasimmod.o`
and `carbonmod.o` only. `make -j16 utilities_omp.o` in a clean directory builds
every declared prerequisite and then stops at `Cannot open module file
'mpiomp.mod'`.

**It does not bite, and what saves it is ordering rather than the graph.** Ten
clean `make -e -j32 plasim.x` builds of the OpenMP configuration all passed,
none failing on a missing module. `${MPIMOD}.o` is FIRST in `OBJ` and
`${UTILMOD}.o` ninth, so make launches `mpimod_omp.o` immediately and it is
finished long before anything asks for its `.mod`. The emptying of `plasim/bld`
removed the stale-`.mod` mechanism that used to hide missing edges; this one is
hidden by `OBJ` order instead.

The fix is one line, adding `${MPIMOD}.o` to the `${UTILMOD}.o` rule. With it
the clean-directory test passes and `mpimod_omp.o` is built as a prerequisite.
It is inert in the other two configurations, being an ordering constraint on an
object already in `OBJ`.

### What parallelism is worth

T21 L10 p1, OpenMP configuration, clean build each time, mean of two:

| | wall | vs serial |
| --- | ---: | ---: |
| `-j1` | 15.8 s | -- |
| `-j8` | 4.9 s | 3.2x |
| `-j16` | 4.5 s | 3.5x |
| `-j32` | 4.5 s | 3.5x |

3.5x, saturating by `-j8`. The bound is the dependency chain through
`plasimmod.o`, which nearly every other object depends on and which nothing can
start before. Rule 4 makes the unit of work `rebuild_binaries.py`'s `MATRIX` rather than one
executable, twelve on the T21/T42/T85/T127/T170 ladder, so at this rung the
operation is roughly 3 minutes against roughly 1 and the higher rungs scale it.

The win available is INSIDE each build and not across builds. Every
configuration compiles in the same `plasim/bld`, which `compile.sh` empties on
entry, so two configurations cannot build concurrently without one deleting the
other's objects. That is why `rebuild_binaries.py` is serial across
configurations and why it should stay that way until the build directory is
per-configuration.

## Round two: the optimisation flags, measured 2026-08-20

The first round refused AMD's libraries and left open whether the COMPILER has
anything to give. Same bed, same instrument, arms rebuilt from the same sources
with only `MOST_F90_OPTS` changed.

### The method correction, which came first and matters more than any arm

**Arms must be interleaved, not run as blocks.** Run in blocks, the stock
baseline measured 12.678 s in one session and 11.959 s in the next, a 6% swing
in one binary, while every other arm reproduced to within 0.5%. The moving arm
was whichever ran first after an idle stretch, which gets the boost clock; the
arms after it run warm. Blocking therefore confounds the flag with the CPU's
thermal and boost state and hands the whole artefact to one arm.

Interleaving round-robins the arms inside each round, so any drift is shared
equally, and the comparison becomes paired: each arm is measured against the
baseline that ran seconds before it. A warm-up run precedes the first timed one.
`exoplasim/notes/rank-layout-benchmark.md` requires a quiet machine, which is
necessary and is not sufficient; this is the other half.

### What the flags are worth

Eighteen interleaved rounds for the pair that decides it, six for the rest.

| arm | effect | restart |
| --- | --- | --- |
| `-march=znver4` | 2.7% median paired gain, faster in 17 of 18 rounds | changed, 1 sha over 18 runs |
| `-ffpe-trap` removed | 0.9%, inside the spread | identical |
| `-fstack-arrays` | 0.2% | identical |
| `-flto` | 1.4% | NOT REPRODUCIBLE, and see below |

`-march=znver4` is real and consistent in direction, and smaller than a blocked
measurement claimed. Isolating compute from the bed's fixed startup, by
differencing a 600-step run against a 60-step one, gives 18.173 ms against
17.797 ms per step, a 2.07% gain; startup is about 1.6 s and does not scale. So
one 5850-step orbit goes from about 106 s of model compute to about 104 s. The
bracket to quote is **2 to 3%**, the lower end being compute-only and the upper
the paired median over the whole bed. It is deterministic across all 18 runs.

`-march=znver4` and `-march=native` produce a byte-identical executable on this
host, so the explicit spelling costs nothing and is what belongs in the flag
line: it records what was targeted rather than what the build host happened to
be.

The size of the gain is the interesting part. The stock build emits **no vector
instructions at all** -- `objdump` finds zero `ymm` and zero `zmm` in the
shipped executable, against 20,371 `zmm` references with `-march=znver4`. So 2
to 3% is what the entire distance from scalar code to AVX-512 is worth here. A
spectral GCM at this resolution is bound by memory traffic and by the
collectives that close every timestep, not by arithmetic throughput, and that
number is the evidence for it. It also bounds what any further codegen work,
profile-guided optimisation included, can be expected to return.

`-ffpe-trap` costs under 1% and is kept. It unmasks the IEEE invalid, divide-by-
zero and overflow exceptions, so a NaN raises SIGFPE at the line that made it
rather than propagating. In a spectral model that distinction is stark: the
Legendre transform sums over all latitudes, so one bad gridpoint becomes every
spectral coefficient within a few timesteps and the origin is unrecoverable.
The failure it prevents is a run that keeps going, writes structurally valid
output, and reaches `assess_convergence` as a plausible number.

`-flto` is not what it looks like and is written up in
`notes/audits/zsolars-restart-overread.md`. The model integrates reproducibly
under it; the varying bytes are uninitialised padding in one restart record,
from a defect the flag exposed rather than caused. That defect is fixed, so this
was never a standing refusal -- re-measured on 2026-08-22 the flag loses every
round and changes the answer, which is the refusal it now rests on.

### Where the orbit goes

Model compute, not writing. The same 600 steps with the production `NOUTPUT`
and `NSNAPSHOT` restored take 12.719 s against 12.678 s with output off, while
producing 230 MB, so the model's own output path is inside the noise. Against a
measured 148 s per orbit and roughly 106 s of model compute, the orbit is about
five sixths compute, and the ceiling on every I/O-side lever is the rest.


## Is `make -j` safe now, measured 2026-08-21

The graph was re-derived from the source rather than read: for every file, the
modules it `use`s, each resolved to the file defining it, checked against the
prerequisites `make_plasim` declares. Run for all three configurations
`compile.sh` can produce, since the module variables change which files are in
`OBJ`.

| configuration | OBJ members | members with no rule | missing edges |
| --- | ---: | ---: | ---: |
| serial, `most_compiler` | 35 | 0 | 0 |
| MPI, `most_compiler_mpi` | 35 | 0 | 0 |
| OpenMP, `most_compiler_omp` | 35 | 0 | **1** |

The one edge is `utilities_omp.o`, which `use`s module `mpiomp` at lines 85 and
106. `mpiomp` is defined in `mpimod_omp.f90`, and the rule names `plasimmod.o`
and `carbonmod.o` only. `make -j16 utilities_omp.o` in a clean directory builds
every declared prerequisite and then stops at `Cannot open module file
'mpiomp.mod'`.

**It does not bite, and what saves it is ordering rather than the graph.** Ten
clean `make -e -j32 plasim.x` builds of the OpenMP configuration all passed,
none failing on a missing module. `${MPIMOD}.o` is FIRST in `OBJ` and
`${UTILMOD}.o` ninth, so make launches `mpimod_omp.o` immediately and it is
finished long before anything asks for its `.mod`. The emptying of `plasim/bld`
removed the stale-`.mod` mechanism that used to hide missing edges; this one is
hidden by `OBJ` order instead.

The fix is one line, adding `${MPIMOD}.o` to the `${UTILMOD}.o` rule. With it
the clean-directory test passes and `mpimod_omp.o` is built as a prerequisite.
It is inert in the other two configurations, being an ordering constraint on an
object already in `OBJ`.

### What parallelism is worth

T21 L10 p1, OpenMP configuration, clean build each time, mean of two:

| | wall | vs serial |
| --- | ---: | ---: |
| `-j1` | 15.8 s | -- |
| `-j8` | 4.9 s | 3.2x |
| `-j16` | 4.5 s | 3.5x |
| `-j32` | 4.5 s | 3.5x |

3.5x, saturating by `-j8`. The bound is the dependency chain through
`plasimmod.o`, which nearly every other object depends on and which nothing can
start before. Rule 4 makes the unit of work `rebuild_binaries.py`'s `MATRIX` rather than one
executable, twelve on the T21/T42/T85/T127/T170 ladder, so at this rung the
operation is roughly 3 minutes against roughly 1 and the higher rungs scale it.

The win available is INSIDE each build and not across builds. Every
configuration compiles in the same `plasim/bld`, which `compile.sh` empties on
entry, so two configurations cannot build concurrently without one deleting the
other's objects. That is why `rebuild_binaries.py` is serial across
configurations and why it should stay that way until the build directory is
per-configuration.
