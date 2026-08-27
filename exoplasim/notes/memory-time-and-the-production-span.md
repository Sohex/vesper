# The memory time grows with the window, and the mean's interval does not shrink

Worldbuilding. Vesper is an invented planet; every number below is a property of
the simulated atmosphere's per-orbit global mean surface temperature, or of the
estimator applied to it.

Measured 2026-08-27 on `canonical-10m-base` at T21, dt 45, from
`run_893e276ee029` -- the vegetated baseline, lakes in its albedo and pedology
soil water in its column -- over **140 clean-I/O orbits**, indices 70 to 209.
One instrument throughout: no low-I/O join is inside the block, which
`segments.py` would refuse anyway.

## What was being tested

`world-k25j` recorded a hypothesis: that a lake-bearing baseline has about three
times the memory of its lake-free bootstrap. It rested on the bare-rock
baseline reading tau 5.49 where the bare-rock bootstrap read 1.82, and it said
in its own text that the magnitude was unmeasured -- an effective sample of 7.5
is where Geyer's estimator returns a number of the ordinary size without support.
140 orbits were bought to settle it.

**The hypothesis is refuted at matched window length.** Over 41 orbits, the same
length the bare-rock arms were judged on, the vegetated baseline reads tau 1.46
against the bare-rock bootstrap's 1.82. The lake-bearing arm is the LESS
autocorrelated of the two. The 5.49 that raised the question was the unsupported
reading its own row warned about.

## What is true instead, and it is a property of the estimator's input

Tau grows monotonically with the length of the window it is measured on, and
every row below is a supported estimate:

| orbits | tau | lag-1 | effective sample | 20 tau |
| --- | --- | --- | --- | --- |
| 20 | 1.00 | -0.067 | 20.0 | 20 |
| 40 | 1.35 | +0.273 | 29.7 | 27 |
| 60 | 1.90 | +0.402 | 31.5 | 38 |
| 80 | 4.54 | +0.518 | 17.6 | 91 |
| 100 | 3.75 | +0.434 | 26.7 | 75 |
| 120 | 6.15 | +0.497 | 19.5 | 123 |
| 140 | 6.30 | +0.484 | 22.2 | 126 |

The lag-1 correlation does not grow with the window; it sits between 0.27 and
0.52 throughout. Tau does, because it sums every lag, and this series has power
at frequencies a short window cannot see. An AR(1) at the measured lag-1 of
0.484 would give tau = 2.88 against the 6.30 measured on the same orbits, so
the process is not AR(1) and its autocorrelation does not decay geometrically.

**The estimator is not at fault.** Batch means use no tau at all -- the standard
error of the mean is taken from the scatter of block means -- and they agree:

| batch length | batches | SE of the mean | SE if independent | inflation |
| --- | --- | --- | --- | --- |
| 2 | 70 | 0.01238 | 0.01017 | 1.22 |
| 7 | 20 | 0.01853 | 0.01017 | 1.82 |
| 14 | 10 | 0.02302 | 0.01017 | 2.26 |
| 20 | 7 | 0.02563 | 0.01017 | 2.52 |
| 28 | 5 | 0.02532 | 0.01017 | 2.49 |

The inflation factor is the square root of tau where the batches are long enough
to be independent. It reaches 2.5 at a batch length of 20 and holds there, which
is tau 6.3 by an instrument that never computes one. The batch-35 row is dropped
from the argument: four batches cannot estimate a scatter.

## The consequence, and it is the useful half

The quantity tau exists to correct is the standard error of a mean. That
standard error **has already stopped shrinking**:

| orbits | SE of the mean, K | SE if independent, K |
| --- | --- | --- |
| 20 | 0.0201 | 0.0201 |
| 40 | 0.0185 | 0.0159 |
| 60 | 0.0191 | 0.0139 |
| 100 | 0.0231 | 0.0119 |
| 140 | 0.0255 | 0.0102 |

Between 20 orbits and 140 the independent expectation falls by half and the real
interval does not move. The effective sample size sits between 18 and 32 at
every length. **Seven times the orbits buy no better a mean**, and the bare-rock
bootstrap's control shows the same over the range it covers: effective sample
14.7 at 20 orbits and 24.0 at 44.

This is not a defect in the run. It is what a series with low-frequency power
does, and the honest reading is that this model's per-orbit global mean surface
temperature can be determined to about **0.025 K and no better**, whatever is
spent.

## What it means for `PRODUCTION_SPAN_TAU_MULTIPLE`

`lib/run_lengths.py` buys twenty tau, on a coverage criterion measured against
synthetic series whose answer is known in closed form. That derivation is sound
and it assumes what a synthetic AR series has and this one does not: a tau that
is a property of the process rather than of the window.

Applied here the rule does not close. At 60 orbits it reads tau 1.90 and asks
for 38; at 140 it reads 6.30 and asks for 126. Every span bought raises the tau
that prices the next one.

**What the criteria actually need is met long before either.** The offset
criterion discriminates at 0.15 K and requires the statistic's standard error to
be at most a third of it, which is 0.05 K. The interval above is 0.019 K at 40
orbits and never rises above 0.026. The margin is a factor of two even at the
worst reading, and it is there at 20 orbits.

So the span is not set by the mean's interval on this world. It is set by the
SLOPE's, which `assess_convergence.py` sizes separately and which does keep
tightening with the window. Pricing the span through tau prices it through the
one statistic that has already converged.
