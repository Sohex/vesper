# Where the time goes, after the init flag and the unroll

*Worldbuilding frame: a COMPUTE measurement of the Vesper project's climate
model on this desktop. Nothing here is about the simulated planet. Measured
2026-08-22 at b4ec421d, T170 on sixteen threads, 300 steps on
`bench/bed_t170cold`, frame-pointer build, warm-up discarded.*

This profile exists to answer one question -- is CLIM-80 still worth doing --
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

## What it retires: CLIM-80

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
own code. That is CLIM-79's load imbalance seen from the other side: per-thread
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

That is CLIM-84, and it is deliberately framed as a question about the SCHEME
rather than about the calls: this project already holds correlated-k tables and
a harness that cross-checks the current band model against them, so "should the
radiation be replaced" has to be priced before "should these `pow` calls be
made cheaper". Optimising a scheme that is about to be swapped is wasted work.

**Priced, and the answer is do not optimise it.**
`exoplasim/notes/radiation-scheme-price.md` counts both schemes' transcendental
work in retired instructions and finds the band-resolved candidate cheaper at
every corner, so the optimisation would be thrown away. It also says where the
`pow` share comes from: every non-integer `**` in the radiation sits inside a
`where` block, which is why none of them reaches libmvec's vector `pow`.

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
and if the transform is bandwidth-bound -- CLIM-79's own measurement has the
SHTns path taking 11.17% of its demand fills from DRAM against legmod's 4.05% --
adding threads to it returns nothing at all.

**And it is the second experiment to collect nothing by REDISTRIBUTING work.**
The hemisphere swap in `thread-count-by-resolution.md` cut the measured
imbalance and bought +0.006%; this removed a real starvation and bought -0.28%.
For a while that read as evidence that `max - mean` over per-thread work
over-predicts what is recoverable. It is not: the section below spends the slack
directly and finds it real, at very nearly the size the barrier attribution
implies. What these two experiments establish is narrower and more useful --
the slack is NOT REACHABLE BY REDISTRIBUTION. Existence and reachability are
different claims and it is worth not conflating them.

**A note on beds.** The first T42 arm returned a confident -2.17% on a 2.94 s
bed against a startup of about 1.6 s, and the harness refused it for exceeding
its scatter floor. Lengthening the bed did not sharpen that number, it REVERSED
its sign. A bed too short to clear its own floor does not give a weak answer, it
gives a wrong one.

## The slack is real, and it is about 4%, measured by spending it

Measured 2026-08-22, T170 on sixteen threads, 300 steps.

Two attempts to collect the barrier's cost by moving work around returned
nothing, which left an ambiguity that mattered: either the slack was real and
redistribution could not reach it, or there was no slack and the per-thread
spread was barrier machinery. Those two lead to opposite decisions about
CLIM-79's large shape, and no amount of further redistribution separates them.

**So it was measured from the other side: spend it.** A control patch on a
throwaway build adds a calibrated busy-wait to every thread EXCEPT the two on
the critical path -- 8 and 9, which the attribution puts at 13.7% and 14.4%
barrier wait against 20 to 24% for the rest -- at the end of the diabatic
timestep, with the delay set by an environment variable so nothing else differs
between arms. If the slack is real, spending it is free until the spend exceeds
it.

| delay a step | work added over 300 steps | wall, median | absorbed |
| ---: | ---: | ---: | ---: |
| 0 | 0.00 s | 48.35 s | -- |
| 2,500 us | 0.75 s | 47.96 s | all |
| 5,000 us | 1.50 s | 47.34 s | all |
| 10,000 us | 3.00 s | 48.85 s | about 2.0 s |
| 20,000 us | 6.00 s | 52.04 s | about 1.9 s |

Against a baseline of about 47.9 s the absorbed amount rises and then SATURATES:
1.5 s vanishes entirely, 3.0 s costs 0.95 s, 6.0 s costs 4.1 s. The saturation
point is the slack.

**About 2.0 s of a 47.9 s run, which is 4.2%**, and it sits just under the six
points the barrier attribution implies as an upper bound -- which is where a
real quantity should sit relative to its bound.

So CLIM-79's premise survives, and the two null results are re-read rather than
explained away: the slack exists at close to the predicted size, and what those
experiments showed is that a static permutation cannot reach it and that fixing
the transform does not touch it. A work queue over the physics is the mechanism
that could, and 4.2% is what it is playing for -- an upper bound still, since a
queue has its own overhead, but now an upper bound with a measured floor under
it rather than an inference.

**The first sweep was designed wrong and is kept as a caution.** It swept delays
of 0 to 800 us, which over 300 steps is 0.24 s of added work against a
run-to-run scatter of about 1 s. It produced a tidy monotone table -- 47.43 s at
zero rising to 48.33 s at 800 us -- that reads as a signal and is entirely
noise, and it would have put the knee in the wrong place. The arithmetic that
catches it is one line: 300 steps times D microseconds is 0.0003*D seconds, so
crossing a slack of seconds needs tens of thousands of microseconds a step, not
hundreds. Same failure as the 2.94 s T42 bed above: an instrument too blunt to
produce a number still produces one.
