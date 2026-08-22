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

NOTE ON WHERE THIS BELONGS: CLIM-52 lives in main's TASKS.md and not on this
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
