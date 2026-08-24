# Is SHTns worth the restructuring it needs?

*A COMPUTE measurement on the Vesper worldbuilding project's climate model.
Nothing in it is about the simulated planet. Measured 2026-08-21.*

## Yes, and by more than the old estimate

SHTns 3.7.5, built with `-march=znver4` so it has the AVX-512 the model's own
flag line grants -- 43,384 `zmm` instructions in the library, checked rather
than assumed -- against `legmod`'s transform on the same problem: T170,
`lmax = mmax = 170`, `mres = 1`, a 256 x 512 Gaussian grid, one thread,
spectral coefficients to a full globe.

The sizes match by construction. SHTns reports `nlm = 14706`, which is exactly
`legmod`'s `NCSP`, so both transform the same coefficients onto the same points
and the comparison is of method.

| | ms per full-globe transform |
| --- | ---: |
| `legmod` `sp2fc` + `fc2gp`, NPRO=1 | 1.45 |
| `legmod` `sp2fc` + `fc2gp`, NPRO=16, x16 | 1.48 |
| SHTns `SH_to_spat`, polar optimisation on | 0.27 |
| **SHTns, polar optimisation OFF** | **0.29** |

**5.1x, in exact arithmetic.** The polar optimisation -- dropping coefficients
below a threshold near the poles, which `legmod` does not do -- is worth only 7%
of that, so the win is not an approximation `legmod` declines to make. The
comparison is against a `legmod` that already exploits north-south symmetry and
has the physics filter folded in, not against a naive transform.

## What it is worth in the model, and what that does to the deficit

The Legendre routines are 16.9% of samples in a frame-pointer profile of the
current build, and the FFT about 2 more.

    18.9% x (1 - 1/5.1) = 15.2% of total runtime

Threads are about 7% behind ranks at T170 -- 5.48% measured, plus the 1.56%
CLIM-48 costs deliberately. So SHTns lands them around **8% AHEAD**, which is
the first time the threading pays for itself rather than being paid for.

## And it removes the weights from the cache budget

SHTns computes the Legendre functions on the fly and stores no `NCSP x NLPP`
matrices at all. That takes `pmat` and `qmat` out of `cache_budget.py` entirely,
30.1 MB a die at T170, which is the whole of what CLIM-48 was for.

## A second confirmation that these loops are multiply-bound

Not what this was measuring, and the sharper evidence of the two. At NPRO=1 a
single thread streams 30 MB of `pmat` and 30 MB of `qmat` per transform; at
NPRO=16 it streams 1.88 MB of each. **The two time the same, 1.45 ms against
1.48.** Sixteen times the weight traffic costs nothing measurable.

That is CLIM-48's finding reached from the other end. There, undoing the filter
fold cost 1.56% and cutting the matrices fourfold bought nothing; here, a
sixteenfold change in weight traffic buys nothing either. Weight footprint is
not what limits this transform. It still matters for the 32 MB portability
target, which is a different constraint and is why CLIM-48 is kept.

## What is not answered

**Whether SHTns can use the thread team.** It is internally threaded over
latitude, which is the axis the model already divides between threads, and the
two cannot both own it. Synthesis is separable by latitude band and analysis is
not -- the Legendre integral needs every latitude, which is what `mpsumsc`
already reduces. That is a design question and this benchmark does not touch it.

**Whether the grid arrays can be made full-globe without rewriting the physics.**
Stage A established that a thread's slice cannot cross a call boundary: a
non-contiguous section makes gfortran copy in and out, and the copy-out lands
after the callee's barrier. Grid arrays are `(NHOR,NLEV)`, so pointing a thread
at `gfull(i0:i1,:)` hits exactly that. If a pointer works, the 169 whole-array
`where` statements never need touching, because the pointer already addresses
only this thread's band. If it does not, the eight-module restructure is back.
That probe is the next thing to run and it needs no SHTns.

## The scope estimate was wrong, and the probe says why

*Measured 2026-08-21.* The standing estimate was that SHTns needs the grid
arrays full-globe, which means 169 whole-array `where` statements across about
700 declarations in eight modules. **That is the cost of one design and not of
the requirement.**

It assumed making `NHOR` itself the globe, so every thread's arrays span every
latitude and every whole-array statement has to become a slice. The alternative
is Stage A's pattern applied to the grid: one shared full-globe array, and each
thread's name for it a POINTER to its own band. Then `NHOR` is still the band
from the thread's side, `where (dls < 0.5) gt = ...` still addresses only this
thread's latitudes, and not one of the 169 statements changes. A whole-array
expression is not a call boundary and has nothing to copy.

What DOES copy is a band passed wholesale to an explicit-shape dummy, because
the stride between levels is the globe and not the band. Measured by comparing
the address the callee sees with the caller's slice, because behaviour cannot
tell -- copy-out makes a copied argument look like a passed-through one:

| how the band is passed | copied? |
| --- | --- |
| 2-D band to explicit-shape `(NHORB,NLEV)` | **yes** |
| 2-D band to assumed-shape `(:,:)` | no |
| 1-D band, one level, to `(NHORB)` | no |
| base address `gfull(i0,1)` to `(NHORB)` | no |

Three copy-free routes, and the last is the idiom this model already uses
everywhere -- `fc2sp(gtdt(1,jlev),...)` is exactly it.

**So the work is bounded by the wholesale passes, and there are 27:**

| receiving routine | passes |
| --- | ---: |
| `gp2fc` | 12 |
| `fc2gp` | 12 |
| `uv2dv`, `fc2sp`, `calcgp` | 1 each |

Twenty-four are the FFT, which SHTns performs itself and which therefore stop
existing. Two more are Legendre routines SHTns also replaces. **`calcgp` is the
only one left**, and it has three fixes to choose from.

The change is the declarations for the transform-crossing arrays, the
association, the transform call sites, and one physics routine. It is Stage A
again, on the other half of the model.

**One risk this does not measure.** A pointer array can generate worse code than
a fixed-shape one, because the compiler must assume it may alias and cannot
assume contiguity. Stage A took that cost on the spectral side and the model
came out faster anyway, but the grid arrays are the ones the physics touches
hardest, and it is a thing to measure rather than assume. Declaring the pointers
`contiguous` where they are is the lever if it bites.

## What SHTns actually offers, surveyed before writing any call site

*3.7.5, surveyed 2026-08-21.* The API maps onto this model's transforms almost
one to one, and two of the flags remove work the integration would otherwise
have had to do by hand.

| what `plasim` does | SHTns |
| --- | --- |
| `sp2fc` then `fc2gp` | `SH_to_spat` |
| `gp2fc` then `fc2sp` | `spat_to_SH` |
| **`dv2uv`**, divergence and vorticity to wind | **`SHsphtor_to_spat`** |
| **`uv2dv`**, wind to divergence and vorticity | **`spat_to_SHsphtor`** |
| `sp2fcdmu`, the mu derivative | `SHsph_to_spat` |

`dv2uv` is the largest single transform routine at 6.64% of samples, and it has
a purpose-built vector transform rather than needing to be assembled out of
scalar ones. Every call also has `_l` and `_ml` variants -- truncated at a given
degree, and per zonal wavenumber -- and the `_ml` family is what a transpose
decomposition would be built from if the chosen one ever fails.

**`SHT_ROBERT_FORM` matches this model's wind convention.** `shtns_robert_form`
makes vector synthesis return the field multiplied by sin(theta) and analysis
divide by it first. `gu` is documented in `plasimmod.f90` as "zonal wind
(*cos(phi))", which is that form. Without the flag every vector transform would
need a cos factor applied by hand, with the pole handling that implies.

**Conventions are configuration, not arithmetic.** Normalisation is
`sht_orthonormal`, `sht_fourpi` or `sht_schmidt`; `SHT_NO_CS_PHASE` controls the
Condon-Shortley phase; `SHT_SOUTH_POLE_FIRST` controls latitude order, and is
left OFF because `inigau` hands this model north first; `SHT_PHI_CONTIGUOUS`
matches the `(NLON, NLAT)` layout the shared grid arrays already have.
`SHT_LOAD_SAVE_CFG` caches the plan, which matters because initialisation is
not cheap and every run pays it.

**THE TRAP, and it is worth stating before it is met.** SHTns's spheroidal and
toroidal potentials differ from divergence and vorticity by a factor of
l(l+1). This model already carries `1/(n(n+1))` inside `fmu` and `fmv`, so a
naive substitution applies it twice. `verify_transform_roundtrip.sh` is the
check that catches it, which is an argument for wiring the identity against
SHTns before any model code moves.

## The GPU, answered on paper because the hardware makes it answerable

The machine has an RTX 4090 and CUDA, and SHTns has `SHT_ALLOW_GPU`. It is not
attractive here, for a reason that is not the one usually given.

| | FP64 TFLOP/s |
| --- | ---: |
| RTX 4090 | 1.29 |
| 7950X3D, 16 cores, AVX-512 | 1.28 |

**Consumer NVIDIA throttles FP64 to a sixty-fourth of its FP32 rate**, so this
GPU has no double precision advantage over this CPU at all. The model is built
`-fdefault-real-8` and `config/planet.yaml` declares eight-byte precision.

**Transfer is NOT the blocker, which is the counterintuitive part.** The
spectral and grid fields crossing the transform are about 128 MB a timestep in
both directions at T170, which is 5.1 ms over PCIe 4.0 against a transform that
currently costs 48.9 ms. Ten percent overhead. The usual objection does not
apply; the arithmetic rate does.

**And GPU mode breaks the chosen parallel design.** `SHT_ALLOW_GPU` documents
that "the same plan cannot be used simultaneously by different threads anymore",
so the GPU forces one thread to drive the transform while the rest wait -- which
costs 140% of current runtime unless the device wins by a wide margin, and it
does not win at all in FP64.

**What WOULD change this is precision, not hardware.** `SHT_FP32` exists and is
GPU-only, and FP32 is 64x FP64 on this part. A mixed-precision transform -- FP32
inside the transform, FP64 state either side -- is the only version of this that
could pay, and it is a decision about how much accuracy the spectral core owes
the climate, not a performance question. `verify_transform_roundtrip.sh` would
quantify the cost exactly: the identity's residual IS the precision loss.

## The resolution ladder and the GPU are independent, and only one is worth it

CLIM-52 specifies a restart converter for progressive T21 to T170 advancement,
and its contract already covers four-byte to eight-byte real conversion in both
directions. The obvious pairing with the GPU finding above is to run the low
rungs in FP32 on the device, where it is 64x its own FP64 rate, and switch to
CPU and FP64 at the top. **The arithmetic inverts that.**

The Legendre transform is O(N^3) and the physics O(N^2), so the transform's
SHARE of runtime grows with resolution. Anchoring on the T170 measurement of
18.9%:

| rung | NLAT | transform share | ceiling on any transform-only speedup |
| --- | ---: | ---: | ---: |
| T21 | 32 | 2.8% | 2.8% |
| T42 | 64 | 5.5% | 5.5% |
| T85 | 128 | 10.4% | 10.4% |
| T127 | 192 | 14.9% | 14.9% |
| T170 | 256 | 18.9% | 18.9% |

**So the GPU is useless at both ends, for opposite reasons.** On the low rungs,
where FP32 would be tolerable because the run is being thrown away anyway, the
transform is a few percent of the work and an infinitely fast one saves nothing
worth a second code path. On the top rung, where the transform finally matters,
the precision the model is built for is the precision this GPU does not do.

The share figures below T170 are EXTRAPOLATED from one measurement through the
O(N^3) over O(N^2) argument, not measured. The conclusion survives being wrong
by a factor of two or three at the low end, which is why it is stated rather
than measured, but a T42 profile would settle it if the question is ever
reopened.

**The ladder is still worth doing, for its own reason.** Its value is not
per-step speed, it is not spinning the expensive resolution up from cold --
total orbits, not milliseconds. Nothing here weakens CLIM-52; it only detaches
it from the GPU.

**And the cap applies to the GPU generally.** Offloading the transform alone is
bounded by 18.9% even at T170. A device that actually paid would need the
physics on it too, which is a whole-model port rather than a library swap, and
would be decided on different evidence than any of this.

NOTE ON WHERE THIS BELONGS: CLIM-52 lives in main's tracker and not on this
branch, so this finding has to reach that row when the branches meet.

## FP32 on the CPU is worth 39.7%, and the GPU is not the interesting part

*Measured 2026-08-21, T170 on sixteen threads, quiet machine.*

`compile.sh -p 4` builds the model in single precision -- it is upstream's
default and this project overrides it to eight. No port, no device, one flag.

| | median | spread |
| --- | ---: | ---: |
| FP64 | 81.18 s | 0.6% |
| FP32 | 48.93 s | 0.4% |

**+39.69%** [+39.63, +39.82], faster in 4 of 4 rounds. It is the tightest
measurement in this line of work and the largest single gain: the whole
threading effort moved T170 from -25.7% to -5.5% against ranks, and one build
flag returns 39.7%.

It is stable in the narrow sense -- 300 steps at T170 completed and wrote a
restart of 213 MB against FP64's 426, which confirms the build really is
four-byte and not silently falling back. AVX-512 doubles its lanes for single,
the bytes moved halve, and the per-die working set falls from 260 MB to 130,
so the FP32 gain and the 32 MB target pull the same way rather than competing.

**THE STABILITY EVIDENCE IS 0.043% OF A SPIN-UP, and that is the open risk.**
300 steps at a 30 minute timestep is 6.2 model days. A 40-orbit spin-up is
691,200 steps. Single precision fails where accumulation is long -- a small
tendency lost against a large state, repeated -- so completing 300 steps says
almost nothing about the workload this is for. The decisive test is a full
FP32 spin-up, which costs hours and is exactly the thing the technique claims
to make cheap.

**It is a spin-up accelerator and not a production mode.** `config/planet.yaml`
declares eight-byte precision, and CLIM-52 treats a converted state as a new
initial condition that must settle at the target rather than as proof of
equilibrium. So FP32 does not replace the FP64 work; it accelerates the phase
that is thrown away, and FP64 does the final approach and the judging.

**Where the FP32 state LANDS is not established.** Its restart records are four
bytes against eight, so `compare_restarts.py` reads the two as different files
rather than different answers. That comparison is the converter's job and is
deliberately not faked here.

This reorders the ladder. FP32 on the CPU is the cheap rung and is available
now; FP32 on the GPU is the expensive rung and needs a whole-model port, since
offloading the transform alone caps at 18.9% even at T170. Try the free one
first.

## The conventions, measured rather than read

*Measured 2026-08-21 at T42, one process, by `probe_shtns_conventions.f90` and
`probe_shtns_vector_conventions.f90`.*

Reading two normalisation conventions against each other and hoping is how a
factor of sqrt(2) survives into a climate. So one spectral mode at a time is
driven through both transforms onto the same grid and the two are divided.

**The configuration that matches this model:**

    shtns_create(NTRU, NTRU, 1, SHT_ORTHONORMAL)      ! NOT SHT_NO_CS_PHASE
    shtns_set_grid(cfg, SHT_GAUSS + SHT_PHI_CONTIGUOUS, eps, NLAT, NLON)
    shtns_robert_form(cfg, 1)                          ! for the vector transforms

**The index map is the same.** SHTns's `LM(l,m)` is 0-based with m outer and l
inner, which is the order `legini` builds `lm` in, so PlaSim `lm` is SHTns
`lm-1`. Confirmed by the result rather than assumed: a wrong map gives a
mode-dependent ratio and these are flat.

**Scalar.** The model's grid is `+sqrt(2 pi)` times SHTns's, for every one of
the 946 modes. The 2 pi is the azimuthal integral, PlaSim's harmonics carrying
no `1/sqrt(2 pi)` in phi.

**Vector**, with Robert form on, against `SHsphtor_to_spat`. Feed the
potentials, read the components:

    S =  divergence / (l(l+1))            spheroidal
    T = -vorticity  / (l(l+1))            toroidal, note the sign
    gu = -sqrt(2 pi) * Vp                 Vp is SHTns's phi component
    gv = +sqrt(2 pi) * Vt                 Vt is its theta component

The `l(l+1)` is the factor this model already carries inside `fmu` and `fmv`,
and it is the one a naive substitution applies twice.

**The toroidal sign is opposite to the spheroidal, and it is not a fudge.** For
`V = grad(S) + curl(T r)` the divergence is the Laplacian of S while the
vorticity is MINUS the Laplacian of T, so the two potentials sit on opposite
sides of their sources.

That sign was missed by the per-mode probes and caught by the dense one. The
probes drove DIVERGENCE only and took the toroidal on trust, which is the
assumption a mixed field breaks; the vorticity-only arm then failed at a
relative error of exactly 2.000, which is what a pure sign flip looks like and
is why the arms are driven separately before they are driven together.

## The Condon-Shortley phase, and the mistake that hid it

**PlaSim carries the Condon-Shortley phase, so `SHT_NO_CS_PHASE` must NOT be
set.** With it set, the ratio is `(-1)^m` times the value above -- correct at
even m and sign-flipped at odd m, on both the scalar and the vector transforms.

The first version of both probes compared the MAXIMUM ABSOLUTE VALUE of the two
grids. A magnitude cannot see a phase, so it reported a clean constant
`sqrt(2 pi)` and the conclusion "the conventions agree as configured" was
recorded with the wrong flag in it. The vector probe caught it only because the
sign happened to alternate against a component that was already expected to
flip, which made the m-dependence visible where a single number was not.

**Compare signed fields, not magnitudes.** A ratio of maxima answers "are these
the same size", and the question was "are these the same field". The fix is a
least-squares scale, `sum(a*b)/sum(b*b)`, which carries the sign a maximum
throws away, and it is what both probes use now.

Had this survived, the model would have run with every odd zonal wavenumber
negated: stable, plausible, and wrong.


## The gate: SHTns and the recipe compute this model

`verify_shtns_equivalence.sh` runs the whole recipe on DENSE fields rather than
one mode at a time, because a per-mode ratio can be right on every mode tried
and still leave the recipe wrong. At T42 against `legmod`, relative:

| arm | |
| --- | ---: |
| scalar, dense spectrum | 1.5e-14 |
| u, v, divergence only | 9.1e-15, 2.0e-15 |
| u, v, vorticity only | 2.0e-15, 9.4e-15 |
| u, v, both together | 3.1e-15, 2.6e-15 |

Its control drops the Condon-Shortley phase, which is the mistake this recipe
was got wrong by once, and fails every arm. So the conversion is proved before
a single call site in `plasim.f90` moves, which was the point.

## Three more terms in the recipe, found by running the model

Measured 2026-08-22, T21, unpaired threaded build.

The gate above was passing while the model disagreed with itself by 100% at the
first step, and the reason it could is `docs/src/practice/failure-modes.md`
class 29: the driver set `plavor = 0` and never read a namelist, so two terms
that are identically one in ITS configuration and never one in the model's were
outside the comparison. Wiring `gridpointa` found all three of these.

**The planetary vorticity.** `sz` is ABSOLUTE vorticity -- `plasim.f90` writes
`plavor` into `sz(3)`, the real part of mode 2, which is l=1 m=0, the harmonic
proportional to sine of latitude -- while the wind comes from the relative part.
`dv2uv` takes it back out of its result; `sh_dv2uv` takes it off the coefficient
instead, as `zt(2) + plavor*shtinv(2)*fsp(2)`, which is legmod's `fmv(2)`
exactly. Leaving it out adds a solid-body rotation, and the model does not fail
where the wind is wrong: it runs, and dies later in the shortwave radiation.

**The spectral filter, which is part of the operator and not a setting beside
it.** `legini` folds `skspgp(n+1)` into `fsp`, `fmu` and `fmv`, so every
spectral-to-grid conversion legmod performs is filtered. The beds run
`nfilter = 2`, the exponential filter, `exp(-8 (n/NTRU)^8)`. An unfiltered
wrapper is a different operator, and its signature is a relative error that
GROWS with total wavenumber rather than announcing itself at n=1: 1.4e-5 at
n=3, 4.6e-2 at n=10, 3.4 at n=16, with the ratio between adjacent n tracking
`((n+1)/n)^8` to two figures from n=10 up. The wrappers take legmod's own `fsp`
rather than rebuilding the filter, so the two cannot drift when a filter is
added.

**The m=0 imaginary parts.** A zonal mean has no imaginary part, legmod never
reads those slots, and nothing keeps them clean --
`verify_transform_roundtrip` measured them as exactly the part a round trip
does not preserve. SHTns has no such null space and takes the coefficient it is
given, so they are zeroed on the way in. On this bed they happened to be zero
already, which is the sort of luck that decides whether a defect is found in an
afternoon or in a climatology.

`shtns_setup` also runs under `!$omp single`. `prolog` runs on the whole team
and FFTW's planner is not reentrant; four concurrent calls segfault inside
`fftw_mkplan_d`.

With all three in, `verify_shtns_model.sh` at T21 on four threads: 1.1e-13 at
one step, growing smoothly to 1.3e-10 at forty, 199 of 199 records at rounding
scale at both one step and twenty, and the SHTns path bit identical run to run.

## SHTns tunes itself, and a tuned library is not reproducible

Measured 2026-08-22, T21.

`shtns_set_grid` with `SHT_GAUSS` calls `shtns_set_grid_auto`, which BENCHMARKS
the algorithm variants available for the grid and keeps whichever wins. The
winner depends on machine timing during setup, the variants round differently,
and so the model gives a different answer from one run to the next. Eight runs
of one binary on one bed produced four distinct restart hashes, in clusters --
three identical, then a different one, then four of those.

`SHT_QUICK_INIT` selects the same grid by a fixed heuristic instead of by
timing. Eight runs, one hash, at one thread and at four.

**The diagnosis cost more than the fix, and the reason is worth keeping.** The
symptom is nondeterminism in a threaded build, which points squarely at a data
race, and there was a plausible one to find: the wrappers write the whole globe
where legmod wrote a band, so both branches were missing a barrier BEFORE the
transform. Adding it fixed five runs in six, which reads like progress toward a
race and is not. ThreadSanitizer, on the recipe in
`shared-spectral-state-trace.md`, found exactly one race in the whole model --
`mpbci` writing the shared `nshtns`, benign, now threadprivate -- and nothing
else.

What settled it was running on ONE THREAD, where the model was still
nondeterministic. No data race survives that test. The library was never
instrumented, because it is C and the sanitizer was told to ignore
non-instrumented modules, so the tool could not have found it however long it
was run.

**Two runs agreeing was what let this stand.** `verify_shtns_model.sh` asked
for bit identity twice and got it, repeatedly, while the underlying
distribution had three or four outcomes. It asks four times now. A property
that holds stochastically is not tested by a pair.

## What the whole transform is worth, and why half of it was negative

T127, sixteen threads, 300 steps, paired and interleaved, one binary with
NSHTNS flipped so no compiler difference is in the comparison. Measured
2026-08-22.

| converted | paired gain | rounds faster | self-scatter |
| --- | ---: | ---: | ---: |
| inverse only | -28.50% [-30.20, -28.24] | 0 of 4 | 0.9%, 2.4% |
| inverse + gridpointa forward | +13.17% [+10.67, +13.66] | 4 of 4 | 6.3%, 5.6% |
| all dynamical-core sites | +16.81% [+7.22, +23.00] | 4 of 4 | 7.5%, 12.2% |

The third row's scatter is well above the 5% floor and its interval overlaps the
second's almost entirely, so **the difference between +13% and +17% is not
resolved** and should not be quoted as an improvement. What the three rows do
establish is the sign, and the size of the swing.

**Half a transform is worse than none, by 41.7 points.** With only the inverse
converted the model paid for both implementations at once: legmod's weight
matrices stayed resident and its forward path still ran, while SHTns's tables
were added beside them, and every `mpsumscp` reduction still executed. The
inverse-only measurement was never evidence about SHTns and reading it as such
would have killed the work one step before it paid.

The gain also exceeds what the transform's own share predicts. `dv2uv` and its
siblings are 14.9% of T127 runtime, so a 5.1x transform caps at about 12%. The
forward conversion additionally DELETES the reductions -- SHTns integrates the
globe and returns the finished field, so there are no partials to sum -- and
that traffic was never inside the 14.9%.

What is still legmod: the `span` diagnostic, the energy and entropy block under
`nenergy`/`nentropy`, and legmod itself, which stays because `nshtns=0` is the
reference the model check compares against.

## An open discrepancy: the analysis arms miss their bar above T42

Measured 2026-08-22. `verify_shtns_equivalence.sh`, tolerance 1e-11 declared
before any of this was written:

| res | NTRU | scalar synthesis | scalar analysis | uv2dv divergence |
| --- | ---: | ---: | ---: | ---: |
| T21 | 21 | 1.0e-15 | 2.2e-14 | 1.1e-13 |
| T42 | 42 | 9.2e-15 | 2.2e-13 | 2.0e-12 |
| T85 | 85 | 2.0e-14 | 1.2e-12 | 2.2e-11 |
| T127 | 127 | 3.2e-14 | 2.3e-12 | 5.9e-11 |
| T170 | 170 | 5.6e-14 | 4.9e-13 | 9.9e-12 |

Synthesis is at rounding scale everywhere. ANALYSIS is two to three orders
worse, and above T42 the vector arms miss the bar.

**What it is not.** Not l(l+1) amplification from returning divergence as
l(l+1) times the spheroidal potential: that predicts growth like NTRU^2 and
T170 is six times BETTER than T127. Not the reference field coming from SHTns's
own synthesis rather than legmod's: swapping it changes the result in no digit,
which follows from the two syntheses agreeing to 3e-14. Not a concurrency
effect; these run on one thread.

**What the model says, which is the question that matters.** At T127, one step,
nshtns=1 against nshtns=0: every continuous field agrees at rounding scale. The
single record beyond it is `dql`, and `dql` is a floor indicator -- 737,280
cells of which 127 are nonzero and every nonzero one is exactly 1.0e-09. SHTns
has 128. One thresholded cell crossed, which a relative metric on a
floor-valued field reports as 76.

So the transform is not visibly wrong in the model, and the array check is not
meeting a bar it met at T42. Those are consistent only if the bar is measuring
something the model does not care about -- most likely the normalisation, since
`spcheck` divides by the largest spectral coefficient while the synthesis arms
divide by the largest GRID value, and the two differ by orders. That is a
hypothesis and it is not tested.

**RESOLVED: it is the Gaussian quadrature weights, and legmod's were wrong.** `probe_shtns_analysis_margin.f90` tests three mechanisms and the third
is the one.

Weights enter the ANALYSIS and not the synthesis -- legmod carries `gwd`
explicitly, SHTns applies its own -- which is the shape of the symptom. Compared
directly, `gwd` against `shtns_gauss_wts` at T127: the ratio runs from
0.999999999997 to 1.00000000001, a spread of 1.02e-10, varying with latitude
rather than a constant normalisation. Adjudicated against numpy's `leggauss`:

| latitude | inigau vs numpy | SHTns vs numpy |
| --- | ---: | ---: |
| 1, pole-most | 5.7e-11 | 4.2e-11 |
| 5 | 5.0e-12 | 1.8e-12 |
| 96, equator | 9.4e-15 | 1.2e-14 |

All three disagree near the pole and agree at the equator, and inigau and SHTns
fall on OPPOSITE sides of numpy, so their difference is the sum.

That reads as "everyone is a bit wrong", and it was wrong to leave it there.
Three double-precision implementations cannot adjudicate each other, so
`gauss_weight_reference.py` computes the same textbook algorithm in x86
longdouble -- eps about 1.1e-19, eight orders of margin on a disagreement at
1e-11 -- and the answer is not symmetric at all:

| latitude | inigau | SHTns | numpy leggauss |
| --- | ---: | ---: | ---: |
| 1, pole-most | 9.9e-11 | 2.7e-16 | 4.2e-11 |
| 2 to 5 | 1.6e-12 to 7.2e-12 | 0, bit identical | 1.7e-13 to 3.1e-12 |

**SHTns is right and `inigau` is the outlier.** numpy is also wrong at the pole,
which is why the earlier three-way comparison read as a wash. `inigau` evaluated
the Legendre polynomial as a TRIGONOMETRIC SERIES and took the weight as
`z4(1-z^2)/z5^2`, which squares that series' error into the weight.

A relative weight error of 1e-10 on products of order 0.1, summed over 192
latitudes with partial cancellation, gives an absolute coefficient error of
order 1e-11 -- which is the floor the analysis arms sit on, field-blind and
band-limit-blind, exactly as measured.

**The other two mechanisms are excluded by the same probe.** Band-limiting the
test field from n<=127 to n<=31 moves the error from 2.5e-11 to 2.0e-11, a
quarter, not the collapse that under-integration of the vector integrand would
give -- and the SCALAR analysis carries the same 2.3e-11, so it is not vector
specific at all. The FFT radix contributes nothing: an SHTns round trip at the
same nlat on nphi=384 (mixed radix, 3*2^7) against nphi=512 (pure radix-2)
gives 1.07e-14 against 1.20e-14, a ratio of 0.887 -- the mixed-radix length is
marginally BETTER.

**What followed: `inigau` was rewritten and the bar never moved.** Newton on
P_n by the three-term recurrence, then w = 2/((1-x^2) P'^2), which is the
textbook algorithm. Against the extended-precision reference at NLAT 32 to 256:
nodes to 1.1e-16 ABSOLUTE, weights to 5.7e-15 at NLAT 32 and 6.9e-13 at 256,
and the weights sum to 2 within 8.9e-16. `verify_gauss_weights.sh` is that
check and its control is the implementation this replaced, which fails at
9.9e-11.

What limits it now is the recurrence rather than the algorithm: evaluating P_n
over n terms accumulates about n*eps and the weight squares the derivative, so
the floor is about 2*n*eps, which is 1.1e-13 at NLAT 256 and is what is
measured. Converging the ANGLE instead of the node, to avoid forming 1-z^2 near
the pole where z is 0.99993, was tried on the theory that the cancellation was
the limit: it changed the weights by nothing outside noise and made the nodes
worse near the equator. The bar is 1e-12, derived from 2*n*eps, and the nodes
are judged in ABSOLUTE terms because they are cosines that pass through zero at
the equator, where a relative bar says more about which latitude sits nearest
the equator than about the algorithm.

The equivalence arms then fell by one to two orders and the 1e-11 bar, declared
before any of this existed, is met everywhere:

| res | scalar analysis | was | uv2dv divergence | was |
| --- | ---: | ---: | ---: | ---: |
| T21 | 9.9e-16 | 2.2e-14 | 9.8e-15 | 1.1e-13 |
| T42 | 7.6e-15 | 2.2e-13 | 5.3e-14 | 2.0e-12 |
| T85 | 1.4e-14 | 1.2e-12 | 1.7e-13 | 2.2e-11 |
| T127 | 9.9e-15 | 2.3e-12 | 2.2e-13 | 5.9e-11 |
| T170 | 9.9e-15 | 4.9e-13 | 6.8e-14 | 9.9e-12 |

Analysis now sits at the same order as synthesis, which is what it should
always have done.

**This is an accuracy improvement to the model and not only to the comparison.**
The weights enter every grid-to-spectral transform PlaSim performs -- `fc2sp`
carries gwd, `uv2dv` and `mktend` carry gwd/cos^2 -- and nothing in the
synthesis direction touches them. So the error was invisible to every check
that compares synthesis, which is how it survived, and it was in every forward
transform this model has ever done.

## SHT_QUICK_INIT is a determinism choice, and there is no table to trade against

`shtns_setup` asks for `SHT_QUICK_INIT` because the default `SHT_GAUSS` times
its algorithm variants at startup and keeps the winner, which is what made the
model give four distinct restart hashes in eight runs. The flag picks by a fixed
heuristic instead.

`nm` on the resulting binary shows only `_fly` kernels -- `SHsphtor_to_spat_-
fly2_m0l` and its siblings. That is not the library declining to stream stored
tables: **at the pinned revision there are no stored tables to stream.** The
matrix-based algorithms are absent from this build of SHTns, `SHT_MEM` is never
populated, and both `choose_best_sht` and `shtns_set_grid_auto` hard-code
on-the-fly as the only family available. What `SHT_QUICK_INIT` selects between
is unroll depth, and the earlier arithmetic here about what stored tables would
occupy per rung described a configuration the library cannot produce.

`shtns-algorithm-selection.md` carries the source reading, the per-rung
candidate set, the measurement showing the tuner returning a different pick
vector on all ten of ten runs, and the authorable-config route by which a
per-rung choice can be declared instead of raced for. CLIM-74 is that work.

## Where the time goes now, T42 and T170

Measured 2026-08-22, threaded build, sixteen threads, `perf record` at 997 Hz on
`cpu-clock` through `profile_transforms.sh PERF_LAUNCH=omp`. ONE pass each and
no repeat, so the wall times are indicative and the SHARES are the measurement.
The T42 bed runs `NENERGY = 1` and the T170 bed does not, so the diagnostics
bucket is not comparable between the two; nothing else here depends on it.

| | T42 legmod | T42 SHTns | T170 legmod | T170 SHTns |
| --- | ---: | ---: | ---: | ---: |
| wall | 5.80 s | 5.08 s | 96.34 s | 76.91 s |
| model code | 51.4% | 48.3% | 58.0% | 51.2% |
| libc | 16.5% | 16.4% | 25.1% | 31.2% |
| libgomp | 21.0% | 25.5% | 8.7% | 10.9% |
| libm | 7.3% | 8.4% | 4.0% | 4.9% |

**The two ends are limited by different things, which is why both were run.** At
T42 the largest non-model cost is libgomp at a quarter of runtime, sixteen
threads over sixty-four latitudes being four latitudes each. That is the price of
the team rather than evidence against it: a thread sweep at T42 finds sixteen
faster than eight by 36% and than four by 124%, and
`thread-count-by-resolution.md` has the arms. At T170 libgomp falls to a tenth
and libc rises to a third.

**SHTns's own kernels are invisible, and that is the result rather than a
measurement failure.** `SHsphtor_to_spat_fly2_l` and its siblings sample at
0.00%, the wrappers at 0.74% together. The saving is almost entirely in the
model's own code -- 55.9 s to 39.4 s at T170 -- which is legmod's Legendre
arithmetic going away. The transform is no longer a thing worth optimising.

**libc is the same 24 SECONDS in both paths at T170**, 25.10% of 96.34 against
31.21% of 76.91. The share rises only because the denominator fell. So the
wrapper temporaries -- the full-globe `zg`, `zvt`, `zvp` one per field-level --
add nothing measurable, which retires the idea that eliminating them is worth a
precision macro. It also means the largest single cost in the model is one the
transform work never touched and cannot touch.

What it most likely is: the grid arrays are full-globe with threadprivate band
POINTERS, and a band is `gd_g(lo:hi,:)`, which is not contiguous -- the stride
between levels is the globe. Every one of those passed wholesale to an
explicit-shape dummy is copied in and out again, and the physics does that
constantly. `plasimmod.f90` says so where the arrays are declared and
`probe_grid_contiguity.f90` measures it; what is new is that it is now the
biggest item. That is CLIM-75 and it is not free to fix: a layout with
contiguous bands makes the globe non-contiguous per level, so the transform
would gather instead, and the two requirements genuinely conflict.

## The legmod sweep: what is converted and what is left

Measured 2026-08-22 by grepping every caller of a legmod transform outside
legmod itself, then checking each site's guard and its default.

**Converted, all of them per-timestep.** `gridpointa` and `gridpointd` both
directions; `mkdheat`'s three `dv2uv` calls, its three analyses and its `zhe`
synthesis; `spectrald`'s forward; `mkdqtgp` in rainmod, which the sweep found --
`nprc` defaults to 1, so it synthesised the adiabatic humidity tendency through
legmod on every step of every run; and `zqout`, the output humidity, which
`nlowio = 1` and `outaccu` make per-step as well.

`zqout` was the one conversion that is not a substitution, and it is written up
under the gate below because the check mattered more than the transform.
`sqout` is THREADPRIVATE: the legmod path fills each thread's copy with the
partial from its own latitudes and `mpsum(sqout,NLEV)` turns every copy into
the whole field, so a wrapper returning the complete field would have that
reduction multiply it by the thread count, and a level-parallel write into a
threadprivate array leaves each copy holding only its own levels. It lands in
shared storage instead, the reduction is dropped, and every thread takes all of
it -- which is what `writesp` and `aasqout` need.

**Left on legmod, and each for a stated reason.**

`span` is inside `ngui > 0 .or. mod(nstep,ndiag) == 0 .or. mod(nstep,nafter) ==
0`, so it is periodic rather than per-step. The `nenergy` and `nentropy` blocks
are off by default. `glaciermod` and `surfmod` filter the orography, which is
startup and occasional. `rainmod_bm`, `rainmod_mca` and `rainmod_kuo_old` are
alternative schemes that `make_plasim` does not compile -- `RAINMOD=rainmod`.

So legmod is no longer on any per-timestep path at all, and it stays compiled
regardless: `nshtns=0` is the reference `verify_shtns_model.sh` compares
against, and deleting it would delete the comparison.

## The gate was comparing zero against zero on the output path

Found while converting `zqout`, 2026-08-22, and it is why that conversion took a
check before it took a line of code.

`verify_shtns_model.sh`'s beds set `NLOWIO = 0`. The output-humidity transform
is guarded by `nlowio > 0 .or. mod(nstep,nafter) == 0`, so at the bed's length it
never fired: `sqout` was never computed, and `aasqsp` -- the only record of it
that reaches the restart, and an accumulator rather than the field itself -- came
out all zeros in both arms. The gate reported agreement on a path it had not
run, and it would have passed a reduction that multiplied the answer by the
thread count.

The bed now sets `NLOWIO = 1`, which is also what production spin-up runs, so
the gate certifies the configuration that ships. `aasqsp` then has 4840 of 5080
values non-zero, both arms agree at rounding scale at 1 and 20 steps, and the
filter control gets sharper as a side effect, its worst going from 12 to 36.

`docs/src/practice/failure-modes.md` class 29, for the third time in this
workstream: a check whose subject is switched off is not a check, and an
accumulator is not evidence about the field it accumulates.

## It was never memcpy: a third of T170 is memset, and it is diffuse

Measured 2026-08-22, T170, sixteen threads, NSHTNS=1.

The cpu-clock profile puts about 31% of the run in four libc addresses a few
bytes apart and 11.6% in one libgomp address. Fetched by build-id from
debuginfod and resolved, they are:

  __memset_avx512_unaligned_erms    ~31%
  gomp_team_barrier_wait_end        ~11.6%

**Not memcpy.** Four adjacent addresses in libc were read as a copy loop for
several hours of work, and the model is ZEROING, not copying. The interposer
built to catch the copying found only 0.33 GiB a timestep going through
memcpy and memmove, which at any plausible bandwidth cannot be a third of the
run -- that negative result is what forced the addresses to be resolved
properly.

**The two heaviest zeroing sites are STARTUP.** `oroini_` and `roffini_` zero
13.18 GiB each, and the counts are byte for byte identical at 5 steps and at
25, which is what proves it. They vanish from a 300-step profile. A short
profile therefore OVERSTATES the libc share badly, and
`score_transform_profile.py --exclude-startup` does not save you: it drops a
startup BUCKET by symbol name, and startup work that lands in a stripped
library falls in `unbucketed` instead.

**Per step it is diffuse.** DWARF unwinding over 300 steps, as a share of the
samples whose leaf is in a stripped library:

| caller | share |
| --- | ---: |
| `gridpointd_` at +0x55, +0x68, +0x7b, +0x3984 | 10.5% |
| `mpgagp_`, `mpgallsp_`, `mpbci_`, `mpsumbcr_` | 9.5% |
| `gridpointa_+0x5627` | 3.3% |
| `sh_gp2sp`, `mkdheat_` | 2.0% |
| callchain truncated | 33.6% |

The three low offsets in `gridpointd` are the block of nine whole-array
zeroings at its head -- `gudt`, `gvdt`, `gtdt`, `gqdt`, `dudt`, `dvdt`, `dtdt`,
`dqdt`, `mmrt`. That is the largest single identifiable site and it is about
2.7% of runtime.

**So there is no low-hanging fruit here, and the reason is semantic rather than
technical.** Those arrays are tendency ACCUMULATORS: the physics adds into them
all timestep, so the zeroing is load-bearing and cannot simply be deleted.
Removing it means making the first writer of each assign rather than
accumulate, across every physics routine that touches them -- a wide,
correctness-sensitive change for a few percent. The collectives are the same
story one level down.

What this does retire is the idea in CLIM-75 that the cost is copy-in and
copy-out at explicit-shape dummies. It is not copying at all.

## Attributing the zeroing needed frame pointers, and they are free

Measured 2026-08-22, T170.

The memset attribution stalled at about a third of samples having no caller.
That looked like clipped stacks and is not: widening the DWARF dump from 8 KB to
32 KB left truncation at 44.1% and 44.4%, and 65528 bytes lost every sample to
the ring buffer. The depth histogram is the tell -- **94.4% of samples landing in
the model's own code had ONE frame**, not a long chain cut short. libdw cannot
find usable CFI at an arbitrary PC inside an -O3 vectorised loop.

Rebuilt with `-fno-omit-frame-pointer` and sampled with `--call-graph fp`, that
figure is **0.0%**. What remains at depth one is 28.7% of the stripped-library
leaves, which are samples inside libc and libgomp themselves; those are not
built with frame pointers and cannot be reached from here.

**It costs -1.17% at T170, interval [-3.18, +0.75], and the restart sha is
IDENTICAL.** So a frame-pointer profile measures the same model, which is the
part that mattered -- an instrument that perturbs the thing it measures is worth
much less. `compile.sh -g` builds it, named `_fp` apart from the registered
binaries.

That inverts what this project had been doing. DWARF was the general tool and
frame-pointer builds were ad hoc for particular investigations; on this build
DWARF fails on nineteen samples in twenty and frame pointers are free.

## And the zeroing is the physics, which is why there is nothing to pick up

With the model side fully unwound, 300 steps, 0 lost samples, 88,661 stacks, of
which 50.4% have their leaf in a stripped library. Where those come from:

| first model frame | share |
| --- | ---: |
| none -- the libgomp barrier | 30.3% |
| `radstep_` | 18.1% |
| `master_` | 16.3% |
| `gridpointd_` | 5.8% |
| `rainstep_` | 5.6% |
| `fluxstep_` | 4.9% |
| `kuo_` | 4.2% |
| `gridpointa_` | 4.1% |
| `icestep_`, `mkdheat_` | 3.4% each |

69.0% go through `master_`, 0.4% through `prolog_`. So it is per-step and it is
spread across every physics routine.

**The attribution is right and the cause was not, and the cause is one build
flag.** `MOST_F90_OPTS` carries `-finit-real=zero`, which initialises every
local real variable on entry to every routine, arrays included -- so the routine
in the table is not clearing its own work arrays, the compiler is clearing them
for it, which is why the cost is spread exactly as widely as the flag is.
Removing it is worth **+26.03%** [+22.76, +27.29] at T170, 4 of 4 rounds, with a
bit-identical restart, and the model is clean under a sNaN gate. The measurement,
the three gates and the one trap it exposed are in
`the-zeroing-is-an-init-flag.md`; the conclusion that there was no low-hanging
fruit here does not survive it.

## The model gate's 20-step tolerance is bed-dependent, and dcc is why

Measured 2026-08-24, T21, on the binaries of a full rebuild from a cleared
object cache.

`verify_shtns_model.sh` FAILS on a bed cut from `run_2b20e3324bb0`, the
85-orbit baseline, and PASSES on one cut from `run_9df7ffa14256`, a one-orbit
cold start. Same build, same thread count, same everything else.

| bed | born at 1 step | at 20 steps | verdict |
| --- | --- | --- | --- |
| 85-orbit baseline | 3.35e-12 | 3.94e-10 | FAIL against 1e-10 |
| cold start, 1 orbit | 5.88e-14 | 2.00e-11 | PASS |

Every OTHER declared criterion passes on both. Both are born inside the 1e-11
birth bound. Neither jumps more than 5.4x between adjacent samples against a
1e4 bound. The SHTns arm is bit identical over four runs on both. The
filter-dropped control is rejected on both.

### One record, and it has a gain the others do not

At 20 steps on the failing bed, `dcc` is 3.94e-10 and the next record that is
not simply its own accumulator is `xcpmea` at 8.3e-12 -- a factor of 47 below.
198 of 199 records are inside the tolerance by more than an order of magnitude.

`rainmod.f90:1977` is why:

```text
zcc = zwfac * AMAX1(0., (zrh - rcrit) / (1. - rcrit))**2
```

with `rcrit(:) = MAX(0.85, MAX(sigma, 1-sigma))` at `rainmod.f90:92`. The
derivative of cloud fraction with respect to relative humidity carries
`1/(1-rcrit)^2`, which is 44 at the eight interior levels and 400 at the top
and bottom. So `dcc` amplifies whatever difference reaches it by up to two and
a half decades, and no other record in the restart does. The prediction from
those constants -- 40 to 400 times the difference the ungained records carry,
which at 20 steps is about 4e-12 -- is 1.7e-10 to 1.7e-9. Observed: 3.94e-10.

The amplification acts only BETWEEN the two clamps: `AMAX1` at zero and
`AMIN1(...,zclmax)` at the top pass a last-bit change through as no change at
all. The failing bed holds 32.1% of its cloud field between them against the
passing bed's 21.1%, which is a factor of 1.5 of the 57x difference in seed
and not the whole of it. The rest is the flow: an equilibrated state carries
more small-scale structure for the two transforms to disagree about.

### What that makes the 1e-10 bound

Not a bed-independent criterion. It was taken unchanged from
`verify_threaded_numerics.sh`, which calibrated it against a difference born
at 3.2e-13; this one is born at 3.35e-12 on the baseline bed, ten times
larger, so the same growth crosses the same bound sooner. The bound therefore
tests seed AND growth together, while the birth bound already tests the seed
and the jump bound already tests the growth's shape.

The script's own stated physics is the bed-independent form: a last-bit
difference grows by roughly three decades every twenty steps. Both beds are
inside that -- 2.07 decades on the failing bed, 2.53 on the passing one -- and
that is the criterion the curve is printed for. Changing the gate's verdict to
rest on it is a decision to take before a run rather than after this one, so
it is not taken here.

WHAT IS NOT IMPLICATED is the transform. At one step on the failing bed every
spectral record is bit identical between the two arms and the largest
difference anywhere outside `dcc` is 1.5e-13.
