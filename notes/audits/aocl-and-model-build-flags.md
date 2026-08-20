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

**`-fcheck=all` is not the lever it looks like.** Removing gfortran's runtime
bounds checking from the production flags buys 5.0%, which lands exactly on the
line that says "no difference", and it changes the restart, because the flag
inhibits optimisations whose absence changes floating-point association. Five
percent does not buy a numerics change, and it does not buy giving up the check:
`oceanmod.f90:236-238` compiles with `Array reference at (1) out of bounds
(2 > 1)` from gfortran itself, so the checking is not guarding a hypothetical.

The model is compute-bound inside its own Fortran, which is where a spectral
GCM should be bound. There is no library seam to widen.

## Two build-system defects found while measuring

Both were closed the day they were found, CONS-11 and CONS-10; the evidence
below is what they were found to be, and `archive/tasks.md` carries what was
done. The dependency defect turned out to be three edges rather than the one
described here.

A CLEAN build of the model fails. `glaciermod.o`'s rule in `make_plasim` names
`plasimmod.o` but not `landmod.o`, while `outmod.o`, which does depend on
`glaciermod.o`, sits EARLIER in `OBJ` than `landmod.o` does. `compile.sh` never
empties `plasim/bld` except when switching between the MPI and serial builds,
so a `landmod.mod` from some previous build is always lying there and the
missing edge never shows. The same missing edge makes `make -j` unsafe at any
level.

`binary_manifest.json` records the sha of every source file that goes into a
build and `exoplasim_version`, and records nothing about the compiler or the
flags. A toolchain change therefore produces a different executable under an
identical set of recorded sources. `--verify` still catches it, since it
compares the executable's own sha, but it reports it as unknown provenance and
cannot say that the toolchain is why. Installing a new compiler or a new libm
is exactly the change that lands here.

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
from a defect the flag exposed rather than caused.

### Where the orbit goes

Model compute, not writing. The same 600 steps with the production `NOUTPUT`
and `NSNAPSHOT` restored take 12.719 s against 12.678 s with output off, while
producing 230 MB, so the model's own output path is inside the noise. Against a
measured 148 s per orbit and roughly 106 s of model compute, the orbit is about
five sixths compute, and the ceiling on every I/O-side lever is the rest.
