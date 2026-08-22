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
