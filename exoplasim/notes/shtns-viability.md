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
