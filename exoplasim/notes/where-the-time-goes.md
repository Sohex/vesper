# Where the time goes, after the init flag and the unroll

*Worldbuilding frame: a COMPUTE measurement of the Vesper project's climate
model on this desktop. Nothing here is about the simulated planet. Measured
2026-08-22 at b4ec421d, T170 on sixteen threads, 300 steps on
`bench/bed_t170cold`, frame-pointer build, warm-up discarded.*

This profile exists to answer one question -- is CLIM-68 still worth doing --
and it answers it no. It is kept because the question will be asked again and
because the shape of the model's cost has moved twice in one day.

## The shape

| object | share |
| --- | ---: |
| the model's own code | 57.78% |
| libgomp, which is the OpenMP barrier | **17.89%** |
| libm | **9.74%** |
| libc | 9.66% |
| libfftw3 | 3.55% |

| model symbol, self | share |
| --- | ---: |
| `swr_` | 10.31% |
| `vdiff_` | 7.63% |
| `calcgp_` | 5.08% |
| `lwr_` | 4.58% |
| `gridpointd_` | 3.46% |
| `kuo_` | 3.14% |
| `mkdheat_` | 2.57% |
| `gridpointa_` | 1.82% |
| `mkrain_` | 1.71% |

**`memset` is 0.00%.** Not reduced -- absent. It was about 31% of this same
build this morning, and it was `-finit-real=zero`.
`the-zeroing-is-an-init-flag.md` has that measurement.

## What it retires: CLIM-68

The row proposed restructuring the physics onto cache-resident column blocks,
NPROMA style. Its stated prize was the 31% memset, on the argument that blocking
does not zero fewer bytes but moves the zeroing into L1 where the physics that
overwrites those lines finds them resident. **There is no zeroing left to move.**

The second-order argument -- that `radstep`, `rainstep`, `fluxstep`, `kuo` and
`icestep` each stream the whole band independently and a resident block fuses
them -- does not survive this profile either, and two independent things say so.

**The largest model routine is arithmetic-bound, not streaming-bound.** `swr_`
is 10.31% and the call graph puts 6.69% of the whole run inside libm's `pow`,
reached only through `radstep` and then `lwr` or `swr`. A cache-resident block
does nothing for a transcendental.

**And `-funroll-loops` returned +11.10%** on this physics, measured the same day
(`notes/audits/model-build-flags.md`). Unrolling pays on loops limited by
instruction throughput and dependency chains; it is close to free on loops
limited by memory bandwidth. An 11% return is the signature of the former, which
is the opposite of the condition blocking exists to fix.

Against roughly 700 `(NHOR,...)` declarations and 169 `where` statements, that is
a refusal. It should be reopened only by a profile showing the physics
bandwidth-bound -- `profile_memory.sh`'s fill breakdown is the instrument, and
the reopening condition is a high DRAM fill share in the physics routines rather
than in the transform.

## What it promotes

**The barrier, at 17.89%**, is now the largest single item outside the model's
own code. That is CLIM-67's load imbalance seen from the other side: per-thread
attribution puts the least-waiting thread at 9.0% of its cycles and the
most-waiting at 18.2%, so roughly six points of runtime are threads waiting for
each other, and a static schedule over NLEV=10 levels across sixteen threads is
where it comes from.

**Radiation's transcendentals, at 9.74%**, have never been measured before. The
old compute roadmap listed "profile and reduce repeated transcendental work in
radiation" as route 3 and nothing had priced it. It is `pow` at 6.69% and `log`
at 1.70%, and the call graph attributes all of it to `radstep` through `lwr` and
`swr`. Two routines and one libm function is a far smaller surface than a
physics-wide restructuring, which is what makes it the better target.

That is CLIM-72, and it is deliberately framed as a question about the SCHEME
rather than about the calls: this project already holds correlated-k tables and
a harness that cross-checks the current band model against them, so "should the
radiation be replaced" has to be priced before "should these `pow` calls be
made cheaper". Optimising a scheme that is about to be swapped is wasted work.

## Filling the transform's team is worth nothing, and that was built to find out

Measured 2026-08-22. The wrappers are `!$omp do schedule(static)` over
`jlev = 1, klev`, and at T170 klev is NLEV=10 against sixteen threads. That is
not a near miss: `schedule(static)` gives threads 10 to 15 no iterations AT ALL,
and it does so independently in every one of the back-to-back calls, so `nowait`
lets a starved thread reach the next loop without giving it work there. Four
calls a step pass klev=1, where fifteen threads get nothing. `sh_sp2gp`'s own
comment already described this and offered `nowait` as the mitigation; the
mitigation does not do what it claims.

**So it was fixed, and the fix is worth nothing.** Batching the (field, level)
pairs into one iteration space -- gridpointa's four scalar fields become 40
units rather than four lots of ten -- was implemented, verified BIT IDENTICAL at
T170 over 60 steps, and measured:

| | paired gain | rounds | verdict |
| --- | ---: | ---: | --- |
| T170, 300 steps | **-0.28%** [-1.02, +0.43] | 2/4 | a clean null |
| T42, 600 steps | -2.17% | 0/4 | refused, scatter 7.8% over the floor |
| T42, 3000 steps | +1.77% [-0.66, +2.27] | 3/4 | readable, spans zero |
| T42, 3000 steps | +0.74% [-4.27, +2.78] | 6/8 | refused, scatter 6.0% |

T170 is a null with tight intervals. T42 is UNRESOLVED rather than null: three
attempts, a point estimate wandering from -2.17% to +1.77% to +0.74%, and an
interval that never clears zero. Whatever the effect is at the low rung, it is
below this machine's noise floor. The code was reverted -- it earns no place --
and the finding is kept here because the next person to read `sh_sp2gp`'s
comment will have the same idea.

**Why it is null is the part worth carrying.** The transform is under 10% of
T170: `sh_dv2uv` 0.94%, the SHTns kernels about 2.8% between them, libfftw3
3.55%. Filling six idle threads through a tenth of the run cannot return much,
and if the transform is bandwidth-bound -- CLIM-67's own measurement has the
SHTns path taking 11.17% of its demand fills from DRAM against legmod's 4.05% --
adding threads to it returns nothing at all.

**And it is the second experiment to say the imbalance model over-predicts.**
The hemisphere swap in `thread-count-by-resolution.md` CUT the measured
imbalance and bought +0.006%. This removed a real starvation and bought -0.28%.
`max - mean` over per-thread work has now failed twice as a predictor of
recoverable wall time, which is why the six points CLIM-67 quotes are an upper
bound and not a target.

**A note on beds.** The first T42 arm returned a confident -2.17% on a 2.94 s
bed against a startup of about 1.6 s, and the harness refused it for exceeding
its scatter floor. Lengthening the bed did not sharpen that number, it REVERSED
its sign. A bed too short to clear its own floor does not give a weak answer, it
gives a wrong one.
