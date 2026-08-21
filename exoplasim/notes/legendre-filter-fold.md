# Folding the physics filter into the Legendre weights

*This is a COMPUTE change to the Vesper worldbuilding project's climate model:
it makes PlaSim's spectral transform do less arithmetic for the same answer.
Nothing in it is about the simulated planet. Measured 2026-08-20 on a Ryzen 9
7950X3D, 16 ranks, gfortran 16.2.1.*

## The change

`legmod.f90` applies a wavenumber-dependent physics filter inside the innermost
loop of every spectral transform:

```fortran
fc(1,m,l) = fc(1,m,l) + qi(w,l) * sp(1,w)*skspgp(n)
```

`skspgp` and `skgpsp` are functions of the total wavenumber alone, so they are
constant for a spectral mode and invariant in the latitude loop that encloses
them. Every use of `qi`, `qj`, `qu` and `qv` carries `skspgp`, and every use of
`qc`, `qe`, `qm` and `qq` carries `skgpsp` -- 94 call sites, no exceptions.

So the filter belongs in the weight matrices, applied once in `legini`:

```fortran
qi(lm,jlat) = qi(lm,jlat) * skspgp(n+1)
```

and deleted from all 94. That removes one multiply in three from the hottest
loops in the model. gfortran will not do it: the value is loop-invariant but
hoisting it needs a temporary array, which the compiler does not invent.

`legini` counts m and n from 0 and the transform routines count from 1, so mode
(m,n) here is `skgpsp(n+1)` there. That offset is the only place this can go
wrong, and the test below is built to catch exactly it.

## The tests, and the first one that could not fail

**Attempt one, which proved less than it claimed.** With `NGPTFILTER = 0` and
`NSPVFILTER = 0` both filters are exactly 1.0, and multiplying an IEEE double by
1.0 is the identity, so the folded model must reproduce the original bit for
bit. It does, at T21, T42 and T85 -- once both arms are built with
`-ffp-contract=off`. Without that, DELETING a multiply changes what the compiler
may contract:

    sp + qc*fc*sk  ->  t = round(qc*fc); t2 = t*sk; round(sp + t2)   two roundings
    sp + qc*fc     ->  fma(qc, fc, sp)                               ONE rounding

so the arms differ in the last bits for reasons unrelated to the fold.

**That test cannot detect a wrong mode index, and the pull request claimed it
could.** With every filter value equal to 1.0, any indexing passes: `skspgp(n)`,
`skspgp(n+1)` and anything else all multiply by 1.0. It proves the fold did not
corrupt the transform. It proves nothing about which mode gets which factor.

**Attempt two, which can fail.** Make the filter a function of the wavenumber
whose values are exact POWERS OF TWO. Scaling by 2^k is exact and rounding
commutes with it, so `(qi*sk)*sp == (qi*sp)*sk` bit for bit -- while a fold that
picked the wrong `n` applies a different power of two to most modes. Alternating
1.0 and 0.5, so an off-by-one changes the factor on EVERY mode rather than only
at a boundary, and with a deliberately wrong arm as a negative control:

| arm | result |
| --- | --- |
| unfolded | `14a1535d8495a0f74c3fae2e66b9a146` |
| folded, `skspgp(n+1)` | `14a1535d8495a0f74c3fae2e66b9a146` -- IDENTICAL |
| wrong, `skspgp(n)` | aborts in `legini` on an out-of-bounds `skspgp(0)` |

`verify_fold_indexing.sh`. The negative control failing is what makes the
positive result mean anything, and it failed harder than designed: the off-by-one
reads index 0 and `-fcheck=all` catches it before the integration starts.

T127 cannot be tested by either route. With the physics filter disabled the
ORIGINAL model reaches a NaN in `lwr_` within 20 steps and `-ffpe-trap` kills
it: the filter is load-bearing for stability at that resolution. Worth knowing
on its own account.

**Under production flags the fold IS a numerics change**, twice over: the
reassociation `(q*x)*s -> (q*s)*x`, and the FMA contraction that becomes
available once the multiply is gone. The restart sha moves at every resolution.

## What it is worth

Paired interleaved A/B, production flags, the arm order flipped every round, two
warm-up runs discarded, one bed per resolution with output off.

| resolution | original | folded | paired gain | rounds won | self-scatter |
| --- | ---: | ---: | ---: | ---: | ---: |
| T21 L10 | 6.08 s | 6.08 s | **-0.02%** | 5/10 | 2.4-3.1% |
| T42 L10 | 12.09 s | 11.96 s | **+0.94%** | 8/10 | 1.8-3.9% |
| T85 L10 | 46.68 s | 46.12 s | **+1.32%** | 8/8 | 1.4-2.2% |
| T127 L10 | 53.01 s | 52.23 s | **+1.59%** | 8/8 | 0.7-0.8% |

With `-fcheck=all` removed as well, T127 gives +1.80% over 8/8 rounds at 0.7%
self-scatter, so the fold and the bounds checking are close to independent.

**The gain is real, monotone in resolution, and small.** It tracks the
transform's share of the model -- 4% of samples at T21, 43% at T127 -- and at
T21 it is nothing at all, a 5/10 coin flip.

It was predicted at 8 to 11% and is not. The prediction assumed the Legendre
loops were flop-bound, so that removing a third of their arithmetic would return
a third of their time. It returns about a twentieth. Those loops are bound by
streaming the `q` matrices and by the accumulator dependency chain; the
multiplier has spare capacity and deleting a multiply mostly buys nothing.
`notes/audits/aocl-and-model-build-flags.md` had already measured the same fact
from the other side, where the entire distance from scalar code to AVX-512 was
worth 2 to 3%.

## The measurement failure worth keeping

The first run of this comparison reported self-scatter from 7% to 59% against a
5% floor, and gains that were indistinguishable from noise. One warm-up run was
not enough: the first TIMED round came in 30 to 40% slow in every arm, because
the machine was still climbing to its boost clock. Two warm-ups, one per arm,
took self-scatter under 1% at T127 and turned a 4/6 coin flip into 8/8.

`bench_ab.py` now guards this. It compares the GAIN against the noise rather
than looking at the noise alone -- a first cut fired on self-scatter by itself
and called a 19.4% gain measured on a 7.3%-scatter bed "not a result", which is
crying wolf on the arm that mattered most. It fired correctly on the T85 arm of
the fold comparison, where the gain was 1.3% against 15% scatter, and that arm
was re-run rather than reported.

## Upstreaming: offered as #64, declined

Offered to `alphaparrot/ExoPlaSim` as PR #64 and declined, correctly.

The maintainer's stated reason was that the fold is invalid -- that it takes
coefficients from one side of a convolution and applies them to the other. That
part is wrong, and `verify_fold_indexing.sh` above is the demonstration: `qi(w,l)`
is the kernel and carries the mode index, so each mode keeps its own factor at
every latitude and nothing crosses the summation. The pull request invited that
doubt by claiming a verification it did not have.

But the decision is right on other grounds, and they are better grounds than
speed. The fold perturbs every spectral transform in two independent ways -- the
reassociation and the FMA contraction -- to buy 1.6% at T127 and nothing at T21,
where most of ExoPlaSim's users are. Two numerics changes for a code cleanup is
not a trade a downstream user with checksum-pinned tests should be asked to take.

**The change stays local**, where the resolutions are high enough for it to
matter and the numerics change is absorbed by a re-commissioning that other work
already requires.

The separate finding next door went out as #65 and is not this patch: two static
array-bounds warnings on unreachable code, whose removal makes the build
warning-clean so a real one would be visible.

## Reproducing

The change itself is a commit under `vendor/exoplasim`; `git log` on
`plasim/src/legmod.f90` is the record. To re-check it:

    exoplasim/scripts/verify_fold_indexing.sh                 # mode indexing, with a negative control
    exoplasim/scripts/verify_fold_exactness.sh T85            # filters off, contraction off
    exoplasim/scripts/bench_ab.py --bed ... --a ... --b ...   # the A/B
