# Where the grid<->spectral pipeline's time goes, and the rules declared before it was measured

*This is a COMPUTE benchmark for the Vesper worldbuilding project: it measures
how the climate model's time divides between its routines on this particular
desktop, and nothing in it is about the simulated planet. Written 2026-08-20,
before any arm was run.*

## The question, and why it is a resolution question

The standing proposal is to replace PlaSim's own FFT and Legendre transform
(`fftmod.f90`, `legmod.f90`) with SHTns. Before the mapping is worth designing,
the share has to be known: an infinitely fast transform saves exactly the
fraction of the integration the transform currently occupies, and no more.

That fraction is not a constant. The Legendre transform costs
`NLAT * NCSP * NLEV` per field, which grows as N^3 in the truncation, while the
grid physics grows as N^2. So the share rises with resolution, and this is the
whole reason the measurement is worth making at more than one:

| | NTRU | NLAT x NLON | NCSP | Legendre work vs T42 | grid work vs T42 |
| --- | ---: | --- | ---: | ---: | ---: |
| T21 | 21 | 32 x 64 | 253 | 0.13x | 0.25x |
| T42 | 42 | 64 x 128 | 946 | 1.0x | 1.0x |
| T85 | 85 | 128 x 256 | 3741 | 7.9x | 4.0x |
| T127 | 127 | 192 x 384 | 8256 | 26.2x | 9.0x |

`docs/src/reference/config-rationale.md` records that T42 was chosen against
execution time rather than against fidelity. That makes the relationship
circular in the useful direction: a transform speedup that is worthless at T42
may be what makes T85 affordable, and the affordable resolution is the thing
being bought. **The T85 arm is therefore the one that decides the answer, and
the T42 arm is the control that says how far the projection had to reach.**

## The arms

| # | resolution | ranks | start | answers |
| --- | --- | ---: | --- | --- |
| 1 | T42 L10 | 16 | warm, off `run_4182235e9781`'s restart | the share at the production configuration |
| 2 | T85 L10 | 16 | cold | the share at the resolution the proposal is really about |
| 3 | T21 L10 | 16 | cold | the third point that says whether the share moves as N^3/N^2 predicts |
| 4 | T127 L10 | 16 | cold | ADDED after arms 1-3, once their trend made T127 the resolution worth costing |

Arm 4 was not part of the original design and its existence is a result of the
first three: the share was still climbing at T85, so the question moved to where
it climbs to. It needed a T127 grid of the same build, its own boundary
conditions and its own executable, none of which existed. Nothing about the
INSTRUMENT changed to accommodate it, which is what keeps it comparable -- same
bed recipe, same 600 steps, same flags, same buckets fixed before arm 1.

Arm 3 exists because two points cannot distinguish a trend from a pair of
numbers. If the three land on the curve the table above predicts, the share at
T106 or T170 can be projected rather than measured; if they do not, the model is
doing something the flop count does not capture and the projection is void.

Arms 2 and 3 start cold because no run exists at either resolution and no
restart can be reshaped into one. That is adequate here and would not be for a
timing: which routines run, and how much work each carries, is set by the
resolution and the namelist. What a cold start changes is the physics BRANCHES a
spun-up state would take, so arm 1 is the one that carries a real climate and
the arms are not compared to each other on absolute time.

## The instrument

`exoplasim/scripts/profile_transforms.sh`, sampling every rank with
`perf record` at 997 Hz on `cpu-clock`, over a bed built by
`make_profile_bed.py`: a copy of a run directory with `N_RUN_STEPS` cut to 600
and every output stream off, which is
`notes/audits/aocl-and-model-build-flags.md`'s bed and instrument.

**The profile runs on the production executable, unmodified and unrebuilt.**
The shipped binaries carry a full symbol table, so `perf` attributes to
`mktend_`, `fc2sp_`, `mpsumsc_` and the rest with no `-g`, no flag change and no
recompilation. There is no instrumented-build confound to argue about, and the
binary is checked against `binary_manifest.json` before the bed is built.

Note which binary that is. A run directory holds the executable ITS orbits were
made with; `run_4182235e9781`'s is `1695e1ae`, which predates both
`-march=znver4` and the multi-species radiation merge. The registered one is
`02a11371`. The bed takes the registered one, because the question is what a
FUTURE run costs.

That choice also sets the baseline honestly. `-march=znver4` means the
comparison is vectorised code against a vectorised library, not against the
scalar build that `aocl-and-model-build-flags.md` measured.

**The trap that makes the obvious spelling wrong.** `perf record` inside an
mpiexec'd rank fails on this host with "Permission error mapping pages", exit
255 and a zero-byte `perf.data`. It is not a permissions problem: `perf_event_-
mlock_kb` is 516 and `RLIMIT_MEMLOCK` is 8 MB, and Open MPI's transport has
locked pages against the latter before perf asks for its ring buffer. What makes
it worth writing down is that it does NOT reproduce outside mpiexec -- perf
tested by hand works, and then every rank in the job writes nothing.
`--mmap-pages 32` is the fix, verified at 16 ranks with
`Total Lost Samples: 0`, and it is `profile_transforms.sh`'s default through
`PERF_MMAP_PAGES`. One threaded process locks nothing before perf asks, so the
refusal itself belongs to the launcher that is gone; the buffer size is kept
because it is the size this workload has been sampled at. A smaller buffer is not free: an overflow
drops whatever ran while perf was behind, which is not a random sample, so
`score_transform_profile.py` refuses to score a pass that lost any.

Two things the first passes settled, both of which were open when the arms were
designed. `02a11371` DOES read a restart written by `1695e1ae`, so the
multi-species radiation merge did not move the restart layout. And arm 2 holds
at `MPSTEP = 22.5`, three passes to an identical sha.

**A cold bed needs `SEED`.** Without it, `initrandom` (`plasim.f90:2025`) takes
the system clock, so every pass integrates different weather and the restart sha
stops being able to tell a different model from a different kick -- which is
exactly what the first T85 arm did, three passes and three shas. The cold beds
now set it and the tripwire works again. The warm arm was never affected: a
restart run has `nkits = 0` and takes no kick.

## The rules, declared here before any arm has run

**The buckets are `score_transform_profile.py`'s and are not arguments.** Two
boundaries decide the answer and both are drawn against interest:

- `spectrala`, `spectrald`, `makebm`, `minvers` and `hdiffo` are NOT transform.
  They are the semi-implicit solve and the hyperdiffusion, arithmetic on
  spectral coefficients that any transform library leaves exactly where it is.
  Counting them would inflate the addressable share.
- `mpsumsc` and the spectral collectives ARE transform. They exist only because
  the Legendre sum is split across ranks. Any scheme that changes who owns the
  latitude axis changes them, so they are counted -- and counted in their own
  bucket, so a speedup cannot be claimed over collectives it does not touch.

**The number being read is the addressable share**: Legendre plus FFT plus
transform collectives, as a percentage of the integration with startup excluded.
It is an Amdahl ceiling and is to be reported as one. A transform made
infinitely fast saves that share; a transform made twice as fast saves half of
it.

**Under 5% is no difference**, which is `rank-layout-benchmark.md`'s floor and
`aocl-and-model-build-flags.md`'s, and is not renegotiated here. Concretely: if
the addressable share at T85 comes in under 10%, the proposal is refused on the
measurement, because no achievable transform speedup can clear the floor from
there.

**A rank spread over 20% voids the arm.** The model is bulk-synchronous and
every rank does the same transform work, so ranks that disagree by much mean the
profile caught something other than the transform. Report it and re-run; do not
quote it.

**The machine must be quiet, and a warm-up must precede the first recorded
pass.** Both halves are required and neither is sufficient alone:
`rank-layout-benchmark.md` for the first, `aocl-and-model-build-flags.md` for
the second, where whichever arm ran first after an idle stretch took the boost
clock and moved 6%.

## What would mean the measurement is wrong rather than surprising

- `aerosol_tracer` or `io` above 2% of the integration. `L_AERO = 0` and both
  output streams are off in the bed, so a fat bucket there is a misconfigured
  bed, not a slow model.
- `unbucketed` heavy enough to move the answer. The buckets were written from
  `nm` on the T42 binary; a resolution or a rebuild can introduce a symbol they
  do not name, and the scorer prints the heaviest offenders for exactly this.
- The three arms disagreeing with the N^3/N^2 table by more than a factor. That
  is not a failed measurement, but it voids the projection to T106 and above,
  which is most of what the measurement is for.
- Arm 1's `plasim_status` sha differing between passes. Every pass integrates
  the same steps from the same state; a pass that disagrees with its siblings is
  a different model, not a slow one.

## What this measurement cannot settle

The share bounds the prize. It does not say what SHTns would actually deliver
against a bucket, and it says nothing about the two structural obstacles that
decide feasibility: SHTns needs every latitude in one process while
`mpimod.f90:130` gives each rank a contiguous block of `NLAT/NPRO`, and
`mktend`, `qtend`, `dv2uv` and `uv2dv` fuse physics into the Legendre weights
rather than being transforms at all. Those are read from the source, not from a
profile, and a share large enough to be worth chasing is the reason to work on
them rather than an answer about them.

One measurement this profile makes cheap and should be taken from the same
samples: the `legendre` bucket's split across `mktend`, `qtend`, `dv2uv`,
`uv2dv` and the plain `fc2sp`/`sp2fc` pair. Only the last two map onto a library
call without restructuring, so that split is the difference between a swap and a
rewrite.

---

## Results, measured 2026-08-20

Three passes an arm, 600 steps, 16 ranks, quiet machine, warm-up discarded.
Zero lost samples anywhere. The three arms were run one after another, never
concurrently, so no arm contended with another.

| arm | wall s/pass | spread | restart sha | transform % of compute | rank spread | wall ceiling |
| --- | --- | ---: | --- | --- | ---: | ---: |
| T21 L10 cold | 6.39, 6.41, 6.42 | 0.5% | agrees | 6.2% [5.1, 7.3] | 34.9% VOID | 5.0% |
| T42 L10 warm | 14.19, 14.72, 14.70 | 3.8% | agrees | 20.6% [16.7, 23.4] | 32.4% VOID | 14.3% |
| T85 L10 cold | 46.76, 46.76, 46.54 | 0.5% | agrees | **41.2% [39.1, 44.1]** | 12.1% | **34.1%** |
| T127 L10 cold | 120.77, 120.40, 121.31 | 0.7% | agrees | **52.1% [49.0, 55.7]** | 12.8% | **43.9%** |

**The T85 and T127 arms answer the question and both pass every declared
guard.** 41% of the model's compute at T85 is the grid<->spectral transform,
and 52% at T127. An infinitely fast transform takes 34% off the wall clock at
T85 and 44% at T127; halving the transform takes 17% and 22%.

The T42 and T21 arms fired the rank-spread guard and their point estimates are
NOT quotable. Both are diagnosed below and neither failure is present in the
T85 arm.

### The bucket that had to be rebuilt, and why the first T42 arm was thrown away

The first arm reported an "addressable share" of 50.4% at T42 with a 60% rank
spread. That number was wrong and the guard is what caught it.

`transform_mpi` had been declared part of the addressable share on the reasoning
that the collectives exist only because the Legendre sum is split across ranks.
That reasoning is right about WHY the collectives are there and wrong about what
the samples in them are. Open MPI polls a shared-memory queue while a rank
waits, so a rank that reaches the reduce_scatter early burns CPU in
`opal_progress`, `mca_btl_sm_poll_handle_frag` and `__vdso_gettimeofday`, and
perf counts every cycle of it as work. The same bucket measured 12.5% on the
busiest rank and 44.5% on the idlest, doing identical transform work.

**The spread was the imbalance, not the transform.** Waiting is slack: a faster
Legendre transform shortens the critical path and lets the idle ranks wait
longer. It does not save their spin. So `mpi` came out of the addressable share
and out of the denominator, the share is now taken over COMPUTE, and the wall
ceiling is computed on the critical-path rank -- the one with least spin, whose
compute is what sets the pace.

Two smaller corrections from the same arm. `outaccu` is called unguarded at
`plasim.f90:688`, every timestep whether or not anything is written, so it is
per-step model work and not the output path its name suggests; it moved out of
`io`. And the remaining column physics got a bucket of its own, so `unbucketed`
now means "a symbol the table does not know" rather than "physics".

### Why T42 and T21 still miss, and why that does not touch T85

T42's transform share of compute runs 16.7% to 23.4% across ranks. The numerator
is sound -- every rank does identical transform work -- but the DENOMINATOR
varies, because a rank holding a polar latitude band does different physics from
one holding the tropics. The floor was declared on a ratio whose denominator was
assumed constant and is not. The bracket is real and the point estimate is not;
the honest T42 reading is **17% to 23% of compute, wall ceiling about 14%**.

T21 fails worse and for an additional reason: 56% of its samples are unbucketed
and the heaviest entries are `pthread_mutex_lock`, `pthread_mutex_unlock` and
`__sched_yield`. That is synchronisation in libc rather than in the MPI
libraries, so it escaped the DSO rule and was counted as compute, inflating the
denominator. T21's 6.2% is therefore a floor, not an estimate. This is the same
fact `docs/src/reference/config-rationale.md` already records from the other side: at T21
the transposes cost more than the arithmetic they carry.

T85's spread is 12.1%, inside the floor, and its unbucketed 21.8% is `pow` and
libm called from the radiation. Nothing that voided the other two arms is
present in it.

### The share rises linearly in N, which is what N^3 against N^2 requires

Transform against non-transform compute, as a ratio: 0.066 at T21, 0.259 at T42,
0.700 at T85, 1.088 at T127.

The two clean arms fix the slope. T85 to T127 multiplies the truncation by 1.494
and the ratio by 1.554, so the exponent is **N^1.09 -- linear inside the
measurement**, exactly as the flop counts predict. T127 was projected at 51% of
compute from the T42-T85 pair before it was built, and measured 52.1%.

That the T42-to-T85 slope looked steeper, N^1.44, is the T42 arm's contamination
showing again: its denominator carries libc-side synchronisation counted as
compute, so its share reads low and the rise off it reads fast. Trust the
T85-to-T127 slope.

Extrapolating one step, **T170 lands near 59% of compute and a wall ceiling
around 53%.** That is a projection, not a measurement, and it is the last one
worth making without building the arm: past T127 the transform is most of the
model and its share stops being the interesting question.

### What shape the work is, which decides swap against rewrite

`split_legendre.py` on the same samples, as a percentage of ALL samples:

| | T21 | T42 | T85 | T127 | maps to |
| --- | ---: | ---: | ---: | ---: | --- |
| plain scalar (`sp2fc`, `fc2sp`, `sp2fcdmu`) | 1.53% | 5.43% | 15.65% | 21.68% | one library call each |
| vector (`dv2uv`, `uv2dv`) | 1.28% | 4.54% | 11.75% | 14.97% | spheroidal/toroidal, a standard primitive |
| fused (`mktend`, `qtend`) | 0.24% | 0.71% | 2.01% | 2.63% | two or three calls plus rescaling passes |
| FFT | 1.08% | 2.33% | 4.08% | 4.16% | subsumed by the library |

**94% of the transform time maps onto a standard library primitive and 6% is
the fused-physics problem, and that split holds to within a point at every
resolution measured.** The fusion obstacle is real and it is small, which is
the opposite of what the source reads like: `mktend` and `qtend` are the most
intricate routines in `legmod.f90` and they are 2% of the model.

Two single routines carry most of it, and their dominance grows with
resolution. At T85 `sp2fc` alone is 14.2% of all samples -- the heaviest
routine in the model -- and `dv2uv` is 10.9%; at T127 they are 19.7% and
14.0%, so the two of them are a third of everything the model does. Both are
INVERSE transforms, and inverse beats forward six to one, 25.1% against 4.2%,
because `sp2fc` is called for sixty-odd field-levels a timestep through `sp2fl`
while the whole forward direction is a handful of fused calls. Anything that
speeds up the spectral-to-grid direction is worth roughly six times the same
work on the grid-to-spectral direction.

The FFT row is the one that does NOT grow: 4.08% at T85 against 4.16% at T127,
flat, because it costs N^2 log N against the Legendre transform's N^3.
Replacing the FFT is worth a fixed few percent at any resolution and is not
where the return is.

### The share falls to 33.6% once the model stops doing work it no longer does

Every arm above was taken on a bed carrying `-fcheck=all` and the energy
diagnostics, because a bed inherits its source run's namelist and
`run_4182235e9781` predates both decisions. Production has since dropped
`-fcheck=all` and turned `NENERGY`/`NENER3D` off. Re-profiled at T127 on a bed
that matches what production now does:

| | as measured above | current production |
| --- | ---: | ---: |
| transform share of compute | 52.1% | **33.6%** [30.9, 36.9] |
| wall ceiling | 43.9% | **27.4%** |
| `legendre` bucket | 39.6% | 22.6% |
| `diagnostics` bucket | 5.07% | **0.73%** |
| same 600-step bed | 120 s | **73 s** |

Rank spread 17.9%, inside the floor, so the new number is quotable.

**The series above is not wrong; it measured a configuration the model no longer
runs.** Two changes moved it, and the larger one is instructive.
`-fcheck=all` inserts a check per array reference, and the Legendre loops are
almost nothing but array references -- so the bounds checking was concentrated
in exactly the code being measured, and removing it shrank the transform's
apparent share by more than a third. A profile of an instrumented build
overstates whatever the instrumentation is densest in.

The resolution TREND is unaffected, because every point in the series carries
the same instrumentation. What changes is the absolute share, and the ceiling
that follows from it.

### The diagnostics bucket was the energy accumulators, not outaccu

0.73% against 5.07%. The bucket was `adener3d`, a `(NHOR, NLEV, 28)` array read
and written every timestep -- about 10 MB at T127 -- and production turned it off
on 2026-08-19. What remains in `outaccu` is around a percent and is not worth
fusing; fusion would not help it anyway, since its statements touch distinct
arrays and there is no reuse between them to exploit.

### Which half of the transform, re-measured

`split_legendre.py` on the corrected bed, as a percentage of all samples:

| | on the old bed | current production |
| --- | ---: | ---: |
| `dv2uv` | 13.95% | **9.97%** |
| `sp2fc` | 19.72% | **5.08%** |
| `fc2sp` | 1.85% | 2.49% |
| `mktend` | 1.86% | 2.45% |
| `uv2dv` | 1.02% | 1.32% |
| `qtend` | 0.77% | 0.96% |
| inverse : forward | ~6 : 1 | **~2.1 : 1** |

`sp2fc` fell four-fold. It is the simplest and most reference-dense loop in the
set -- two array reads and a write per iteration, almost no arithmetic -- so it
was paying the most bounds-check tax and had the most to give up when the checks
went. **`dv2uv` is now the largest single transform routine**, which matters
because it is the hard case: `exoplasim/notes/symmetric-transforms.md` derives
its parity split and it needs sixteen live accumulators, with a spill risk that
has to be measured rather than predicted.

So symmetry is worth about **13.7% of wall at T127** rather than the 19% the
earlier split implied, and most of it sits in the routine least likely to give
it up easily.

### The forward-direction symmetry saving, measured

`legmod` carries a symmetry-conserving path for the forward transforms and
production never reaches it, because `mpimod.f90:130` scatters contiguous
latitude blocks and a mirror pair lands on different ranks.

Measured directly, at one rank where both paths are reachable and nothing else
differs -- the same binary built twice, once stock and once with the branch
condition forced to the parallel path:

| arm | median | spread |
| --- | ---: | ---: |
| forced-parallel | 9.27 s | 1.7% |
| symmetric | 9.04 s | 1.5% |

**+2.31%, faster in 12 of 12 rounds, range [+1.25, +3.62].** T42, L10.

Two things follow. The forward half of the symmetry work is worth roughly what
the flop count says and no more, so enabling it alone does not pay for the
decomposition change it requires. And the harness can resolve an effect of that
size, which is what this arm was for: it is the calibration that makes a larger
result on the inverse transforms believable.

### What is still not settled

The share bounds the prize; it does not say what a library would deliver against
it, and it says nothing about the decomposition. SHTns needs every latitude in
one process and `mpimod.f90:130` gives each rank a contiguous block of
`NLAT/NPRO`, so the 34% ceiling at T85 is reachable only through a transposed
layout that also rewrites the collectives -- and `mpi` is 17.2% of samples at
T85, which such a rewrite would move rather than remove.

The factor sitting in the unused north-south symmetry path is measured here only
indirectly: it acts on the 33.5% of samples the table above accounts for, and it
needs no library at all.
