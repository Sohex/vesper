# How long a run takes to settle, and how long a window has to be to tell

Worldbuilding frame: every number here is a diagnostic of the Vesper climate
model. Nothing in this note is about the real world.

Two different questions get confused because both are answered in orbits.

- **How long until the model has settled.** A property of the approach.
- **How long a window has to be before a criterion can SEE that it has.** A
  property of the model's variability, and independent of the first.

The second was never asked. `assess_convergence.py` used a ten-orbit window and
`compare_equilibria.py` took its standard error from the orbit-to-orbit scatter
inside one, both of which assume consecutive orbits are independent samples.
They are not, and the consequences were being read as results.

## The approach: operational experience, not a distribution

These are what this project has repeatedly seen. They are experience and are
marked as such wherever they are quoted, in `assess_convergence.py`'s
`CONVERGENCE_LENGTHS` and here. Nobody has run the ensemble that would turn
them into a measured distribution, and until somebody does they are the scale
and not the number.

| | orbits | |
| --- | ---: | --- |
| cold start at T21 | about 70 | seed- and initial-condition-dependent. `parameter-decisions.md` records one converging on all six criteria after 70 |
| reconvergence after a timestep change | 10 to 20 | the state is already at a climate |
| reconvergence after a resolution conversion | 10 to 20 | likewise; what settles is the response to a changed discretisation |

The consequence for a window is immediate. A ten-orbit window is a seventh of a
cold start's approach and the WHOLE of a reconvergence, so on the second it is
assessing orbits that are still moving by construction.

## The variability: what a window mean is worth

For a stationary series with memory, the variance of a mean over n consecutive
samples is `sigma^2 * tau / n`, where tau is the integrated autocorrelation
time: 1 for independent samples, `(1+r)/(1-r)` for an AR(1) process with lag-1
correlation r. The effective sample size is `n / tau`.
`lib/autocorrelation.py` owns that arithmetic and is the only place in the tree
that should be doing it.

**Measured 2026-08-25, on synthetic series whose answer is known in closed
form.** At r = 0.615, the lag-1 correlation measured on an 85-orbit T21 pair,
the count-based standard error is understated by a factor of 1.83 at n = 10,
1.93 at n = 20 and 2.01 at n = 70. The tau-corrected form matches the empirical
spread of independent window means to within a few per cent and errs
conservatively at short n. `scripts/smoke_test.py` runs this as a check with a
right answer.

**And tau cannot be recovered from inside the window it corrects.** Same
measurement, same true tau of 4.195: a 10-orbit span returns a median 1.46, a
20-orbit span 2.28, a 40-orbit span 2.90, an 85-orbit span 3.54. Correcting a
20-orbit window by a tau estimated from that window therefore recovers about
half the missing error and still looks rigorous. Every estimate here has to be
taken over a longer span than the window it is applied to, which is what
`compare_equilibria.py` now does and what `assess_convergence.py` reports it is
NOT doing when the run is too short to allow it.

## What the window has to be

The offset criterion is the one that bounds the answer rather than a rate. On a
settled run the exponential fit has nothing to grip, the fallback takes over,
and the statistic thresholded at 0.15 K is `|offset| + half_width`, where
`offset` is `slope * tau_expected` and `half_width` carries the same slope
error times the same `tau_expected`, about 10 orbits. So the slope's own error
is multiplied by ten and then counted twice, and the window that criterion
needs is set by

    var(slope) = sigma^2 * tau * 12 / (n (n^2 - 1))

Taking the bar as "the standard error of the statistic is at most a third of
its threshold", declared before any of this was read, and the T21 baseline's
stationary scatter of 0.07 K:

| tau | window needed | |
| ---: | ---: | --- |
| 1.0 | 21 orbits | if consecutive orbits were independent |
| 4.2 | 34 orbits | at the lag-1 correlation measured on the 85-orbit pair |
| 9.0 | 44 orbits | at r = 0.8 |

**So the window is bracketed at 21 to 44 orbits and the estimate is 34.1,
against the 10 that was in use.** `assess_convergence.py`'s default is the
ceiling of that estimate and is computed from the three inputs above rather
than typed, so it cannot be moved without moving a number that has a citation.
The bracket is what the memory costs: the independent case is the floor and
every orbit above it is bought by the autocorrelation.

Widening is close to free in orbits and buys a verdict that means what it says.
The ten-orbit window was measured to flap for twenty orbits after the state had
stopped moving, which is the same order as the extra orbits a window of this
width needs before it can return anything at all. What changes is that the
verdict at the end is a verdict. The cost lands hardest on a reconvergence,
which settles in ten to twenty orbits and would now be assessed over more
orbits than it took to settle; that is the correct answer rather than an
awkward one, because a window shorter than the approach it follows is
assessing orbits that are still moving.

`assess_convergence.py` computes this per run rather than carrying the table:
`resolving_power.window_orbits_for_offset_criterion`, from the run's OWN
residual scatter and its own autocorrelation, so it follows the resolution up
the ladder instead of being a number swept once at T21. Where tau had to be
taken from inside the window, the report says so and the required window is
labelled a lower bound.

## Neither run on disk is stationary, and the guard says so

**Measured 2026-08-25.** `stationary_enough` refuses a span whose fitted trend
moves it by more than its own scatter, a criterion fixed before it was applied.

- `run_ade7373b4c90`, 28 orbits, T21 cold start. No tail of it passes: at
  orbits 20 to 27 the fitted drift is still +0.011 K per orbit against a
  scatter of 0.065 K. A 28-orbit cold start is a third of the way through a
  70-orbit approach, so this is the expected answer and not a defect.
- The four non-overlapping 20-orbit means of the 85-orbit T21 run, 291.355,
  291.835, 291.757 and 291.651 K, are refused as a span too: the fitted trend
  moves them 0.243 K against a scatter of 0.210 K. The first of those windows
  covers orbits 0 to 19 of a cold start whose approach runs to about 70, so the
  0.210 K spread is an UPPER BOUND on the model's equilibrium variability
  contaminated by the tail of the approach, not a measurement of it.

The project therefore has no equilibrium series long enough to measure tau on,
and both instruments now say that rather than assuming a number.
`compare_equilibria.py` reports a metric indeterminate and refuses the verdict;
`assess_convergence.py` reports the window its criteria would need and labels it
a floor. What would settle it is a production span at equilibrium of at least
ten times tau, which on the estimate above is about 40 orbits AFTER the
approach has finished.
