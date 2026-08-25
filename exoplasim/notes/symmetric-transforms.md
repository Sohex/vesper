# Symmetric spectral transforms: the derivation, and what it turned out to be worth

*This is a COMPUTE change to the Vesper worldbuilding project's climate model.
Nothing in it is about the simulated planet. Derived 2026-08-20, built and
measured 2026-08-21.*

**THE SYMMETRIC PATHS DESCRIBED HERE ARE DELETED.** `2508bedb` removed them
along with the decomposition that made them reachable, for the reason given in
`paired-latitude-decomposition.md`: SHTns replaced `legmod` and cannot use a
permuted layout, so a mirror pair no longer lands on one thread and there is
nothing for the algebra below to run on. `sp2fc`, `sp2fcdmu` and `dv2uv` are
plain full-latitude multiplies over `pmat` and `qmat` again. The deletion is bit
inert and that was checked: the T21 restart sha is `a7418dea572f853f` on both
sides of it.

The record is kept for the derivation, which is the part that would have to be
redone rather than reread if a mirror pair ever sits on one thread again, and
for the measured ladder at the end, which is what the two halves together were
worth and is the only number anyone should quote for them.

The decomposition that made any of this reachable in production is a separate
change with its own note, `paired-latitude-decomposition.md`. Read that one
first: without it a mirror pair sits on two different processes and none of the
algebra below can be used.

## What is on the table

A latitude and its mirror carry the same associated Legendre values up to a
sign, so one pass over the spectral modes can produce both. PlaSim uses this in
the FORWARD transforms only -- `fc2sp`, `uv2dv`, `qtend`, `mktend` -- and
production never reaches even those, because `mpimod.f90:130` scatters
contiguous latitude blocks and a mirror pair lands on different ranks.

The inverse transforms `sp2fc`, `sp2fcdmu` and `dv2uv` have no symmetric path at
all. `exoplasim/notes/spectral-transform-profile.md` puts them at **33.7% of
samples at T127**, against 4.2% for the whole forward direction. That asymmetry
is the reason this is worth doing: enabling what already exists is worth about
2%, and writing the missing half is worth an order of magnitude more.

## The parity premise, checked rather than assumed

Everything rests on two relations. With `s = (-1)^(m+n)`:

    P(-mu) =  s * P(mu)
    Q(-mu) = -s * Q(mu)        Q = the mu-derivative legini stores in zpld

Calculus gives these immediately -- the derivative of an even function is odd --
but the model does not use calculus, it uses a recurrence, and a sign convention
inside a recurrence is exactly the kind of thing that makes a plausible-looking
parity split wrong.

`verify_legendre_parity.py` reproduced `legini`'s recurrence in Python,
evaluated it at `+mu` and `-mu`, and checked both relations mode by mode. At
NTRU=10, all 66 modes, the worst relative violation was **0.00e+00** --
bit-exact, since the recurrence is symmetric in `mu` up to sign. The premise
held.

**That gate is DELETED, under `world-x3zd`**, and it is the fourth of the four
`2508bedb` left pointing at code that no longer exists;
`paired-latitude-decomposition.md` carries the other three verdicts. The
premise is a precondition of the parity split above and of nothing else in the
tree, so with the split gone it has no consumer. What it checked was also its
own second hand-transcription of `legini`'s recurrence rather than `legini`, so
a drift between the two would have passed. Anyone who needs the premise again
takes it from the compiled `legini` the way `verify_weight_factorisation.py`
does, which is a repair of a different shape from the one this script was.

So the weight matrices divide by which polynomial they carry:

| matrix | built from | parity |
| --- | --- | --- |
| `qi`, `qu` | P | `s` |
| `qj`, `qv` | Q | `-s` |

## `sp2fc` and `sp2fcdmu` are the easy case

One accumulator pair per Fourier component. Split the mode sum by parity, then

    fc(m,l) = E + O
    fc(m,k) = E - O           for qi;  E - O and E + O reversed for qj

using only `qi(w,l)`. Multiplies halve, and the weight-matrix traffic halves
with them, which matters more: at T127 each `q` is 774 KB per rank.

## `dv2uv` is not, and this is the derivation the plan asked for

`dv2uv` mixes both parities in every output, because it uses `qu` (from P) and
`qv` (from Q) together:

    pu1 += qv*pz1 + qu*pd2
    pu2 += qv*pz2 - qu*pd1
    pv1 += qu*pz2 - qv*pd1
    pv2 += -qu*pz1 - qv*pd2

For the mirror latitude `k`, substitute `qu(k) = s*qu(l)` and `qv(k) = -s*qv(l)`:

    pu1(k) += -s*qv*pz1 + s*qu*pd2  =  s * ( qu*pd2 - qv*pz1 )

**So `u` and `v` do NOT share a single parity split**, and the two terms of each
output must be accumulated separately. Writing `E`/`O` for the even/odd partial
sums of each product:

    pu1(l) = (Eqv + Oqv) + (Equ + Oqu)
    pu1(k) = (Equ - Oqu) - (Eqv - Oqv)

and likewise for the other three components. Eight distinct products --
`qv*{pz1,pz2,pd1,pd2}` and `qu*{pz1,pz2,pd1,pd2}` -- each needing an even and an
odd accumulator: **sixteen accumulators** carried through the `n` loop, reduced
once per `m`.

The saving survives that bookkeeping. Unpaired costs 8 multiplies per mode per
latitude, so 16 for the pair; paired costs 8, plus a fixed reduction per `m`
rather than per mode. Multiplies halve and each `qu`/`qv` element is read once
instead of twice.

**Sixteen live accumulators is the risk**, not the algebra. If it spills, the
saving goes to stack traffic and the change is worth nothing -- which is a thing
to measure on the built code, not to predict. `sp2fc` and `sp2fcdmu` carry no
such risk and should be built and measured first, so that a null result on
`dv2uv` is attributable.

## Not a `mod(m+n,2)` test in the inner loop

The existing forward branches test parity inside the innermost loop. That is a
branch in the hottest code in the model and it blocks vectorisation. Split the
`n` loop into two strided loops instead -- `n = m, m+2, ...` and
`n = m+1, m+3, ...` -- which needs no test at all. Worth applying to the
existing forward branches at the same time.

## Built, and checked against something that can fail

All three inverse transforms carried a symmetric path, selected by `LPAIRLAT`
and running over `NLHP` mirror pairs, with the `mod(m+n,2)` test gone from all
of them: the `n` loop was split into two strided loops instead, which needs no
test and leaves the innermost loop branchless.

The test with a right answer lifted the three loops VERBATIM out of
`legmod.f90`, compiled them against a stub module whose weight matrices it
supplied -- an associated Legendre table computed by scipy -- and ran one
spectral mode at a time. The transform is by definition the matrix-vector
product with that table, so the answer at every latitude was known in advance,
INCLUDING the mirrors the symmetric path never reads. Supplying the table also
made the check free of any normalisation or phase convention: the routine was
asked to reproduce what it was given.

66 modes at NTRU=10, 16 latitudes, both Fourier components, all four of
`dv2uv`'s outputs. Three negative controls, one per routine, each flipping
exactly the sign the derivation above warns about:

| control | rejected on |
| --- | ---: |
| `sp2fc`, parity sums swapped at the mirror | 528 values |
| `sp2fcdmu`, same | 520 values |
| `dv2uv`, the sign joining its two opposite parities | 520 values |

End to end, `verify_paired_decomposition.sh` at 8 processes: 143 of 199 restart
records bit identical, 56 at rounding scale, worst 1.9e-11 against a 1e-10
tolerance declared before the arms ran, and the wrong-permutation control out
by 6.3e15.

**What survives that gate is `verify_inverse_transform.py`**, under world-2div.
The mirror controls perturbed lines that no longer exist and the gate refused
to run, so it was re-scoped to the half that never depended on the layout: the
three routines against the matrix-vector product with the table AND the
per-mode factor they are handed, which is what the fork's surviving change to
them is. Seven controls, anchored to lines that exist.

## What it is worth

Against stock behaviour -- `LPAIRLAT` forced false gives contiguous blocks, the
forward symmetric branches unreachable and no symmetric inverse path at all --
at 16 processes, 4 interleaved rounds, order flipped each round, one driver
owning the machine:

| | contiguous | symmetric | gain | spreads | rounds |
| --- | ---: | ---: | ---: | ---: | ---: |
| T21 | 48.75 s | 48.18 s | +0.88% [-0.23, +1.27] | 2.1 / 2.7% | 7/8 |
| T42 | 38.03 s | 37.34 s | **+1.75%** [+1.46, +2.95] | 1.6 / 1.2% | 8/8 |
| T85 | 40.95 s | 39.14 s | **+4.34%** [+3.63, +5.06] | 1.4 / 0.7% | 4/4 |
| T127 | 39.94 s | 37.17 s | **+6.92%** [+6.35, +7.21] | 0.7 / 0.2% | 4/4 |
| T170 | 85.46 s | 75.55 s | **+11.65%** [+11.50, +12.21] | 1.0 / 0.3% | 4/4 |

T21 and T42 ran eight rounds rather than four because at four they were inside
the noise floor. **T21's interval still includes zero and it is reported as
what it is**: a point estimate under one percent, positive in seven rounds of
eight, on a bed where the whole transform is 6% of compute. It is not a claim
that the change is worth 0.88% at T21; it is a claim that it is worth
approximately nothing there, which is the useful thing to know for anyone
running the model at the bottom of the ladder. T42 at 8 of 8 rounds and an
interval clear of zero is a result.

Two earlier attempts at the low end are not in the table and should not be
quoted from the JSON: the first used six and nine second beds and scattered 7
to 12 percent, and the second was wrecked by gfortran runs launched alongside
round 1 -- the same contention mistake that had already destroyed one A/B in
this project, made a second time. The beds here run about forty seconds and
nothing else was on the machine.

## The accumulators did not spill

That was the one thing this note said had to be measured rather than predicted,
and the shape of the column answers it. The gain rises monotonically and
ACCELERATES with resolution: +0.9, +1.8, +4.3, +6.9, +11.7, every step a larger
increment than the last. Longer mode loops mean more work
avoided, and that is what a working `dv2uv` looks like. A spill would bend the
curve down exactly where it bends up, because the accumulators are live across
the whole `n` loop and it is the long loops that would pay the stack traffic.

The arm spreads tighten as the gain grows, which is the same fact from the other
side: at T170 the symmetric arm scatters 0.3% across four rounds.

## What is left of the estimate, and what was wrong with it

`spectral-transform-profile.md` put the whole symmetry prize at about 13.7% of
wall at T127 and this returns 6.9%, half of it. The gap is not mysterious and is
worth stating rather than rounding away: the profile priced the transform
routines' SAMPLES, and halving the arithmetic inside a routine does not halve
its wall time when part of it is streaming weight matrices and completing
collectives. The same thing was got wrong in the other direction on the filter
fold, where 8-11% was predicted and 1.6% measured.

**What the profile did get right is the ordering**, which is what it was for:
the inverse transforms were the place to work, `dv2uv` was the largest single
routine, and the payoff grows with resolution.

## Offering it upstream

This section is what the pitch WOULD have been. It is kept because the trade it
states is the trade anyone reviving this has to make again.

Unlike the filter fold, this is a capability rather than a cleanup: it removes a
restriction -- the symmetric path being unreachable on more than one process --
rather than reassociating an existing computation for a fraction of a percent.
It defaults to stock behaviour wherever `NPRO` does not divide `NLAT/2`, and at
one process it generates the code that was already there.

The honest pitch to a maintainer is the ladder, not the headline: **it is worth
approximately nothing at T21, which is where most users are**, and it grows to
a tenth of wall clock at T170. So it is for people pushing the resolution, and
it costs everyone else a decomposition change that has to be right. That is a
real trade and it should be stated as one rather than led with the 11.65%.
