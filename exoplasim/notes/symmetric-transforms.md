# Symmetric spectral transforms: the derivation, before any of it is written

*This is a COMPUTE change to the Vesper worldbuilding project's climate model.
Nothing in it is about the simulated planet. Written 2026-08-20, before the
implementation.*

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

`verify_legendre_parity.py` reproduces `legini`'s recurrence verbatim, evaluates
it at `+mu` and `-mu`, and checks both relations mode by mode. At NTRU=10, all
66 modes, the worst relative violation is **0.00e+00** -- bit-exact, since the
recurrence is symmetric in `mu` up to sign. The premise holds.

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

## What would falsify this

- `verify_legendre_parity.py` failing on a future compiler or a changed
  recurrence. It is cheap; run it before trusting any of the above.
- A single-mode synthesis disagreeing with an independently computed
  `P_{m,n}(mu)`. That is the test with a right answer and it is what catches a
  parity split that looks correct and is not.
- `dv2uv` showing no gain once built, which would mean the accumulators spilled
  rather than that the derivation is wrong. Distinguish the two before
  concluding.
