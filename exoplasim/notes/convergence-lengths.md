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

## What the storage criterion needs, in its own units

The storage criterion tests `storage_w_m2_least_squares`, which is the
least-squares slope of the planetary heat content across the window divided by
the orbit. So it obeys the same `var(slope) = sigma^2 * tau * 12 / (n (n^2 - 1))`
as the temperature criteria, with sigma the heat content's residual scatter in
J/m2 and the threshold converted into J/m2 per orbit by multiplying by the
orbit. Nothing here passes through a temperature sensitivity, which is the
point: the threshold and the estimator are put in one unit and compared.

**Measured 2026-08-25 on the dt-30 run's 35-orbit production window.** The heat
content's residual scatter is 1.09e7 J/m2 with a lag-1 correlation of 0.80 and
an integrated autocorrelation time of 6.70, the last estimated from inside the
window and therefore a lower bound. The standard error of the storage slope:

| window | standard error, W/m2 | three times it | resolves 0.12 |
| ---: | ---: | ---: | --- |
| 10 | 0.197 | 0.590 | no |
| 20 | 0.069 | 0.208 | no |
| 28 | 0.042 | 0.125 | no |
| 35 | 0.030 | 0.090 | yes |
| 44 | 0.021 | 0.064 | yes |

**The criterion needs at least 29 orbits and it has 35, so it is a criterion
now and was a formality before.** At the ten-orbit window this project used
until the offset criterion's window was derived, the estimator's standard error
was larger than the threshold itself: a statistic that cannot see its own bar
passes or fails on its own noise. The window needed is 15.3 orbits if
consecutive orbits were independent, 28.8 at the tau measured inside the window
and 31.8 at a tau of 9, and every one of those is a floor for the same reason
the offset criterion's is. It does not bind the window: the offset criterion
needs about 38 orbits on the same run, so the window is still set by the
criterion that bounds the answer rather than a rate.

**One half of the 0.12 threshold's derivation was not what it claimed.** The
threshold was argued from two bounds that met: 0.118 W/m2 from the offset
tolerance, and "an estimator that resolves 0.11 W/m2 on a 10-orbit block", the
second taken as the difference between a run's 10-orbit and 20-orbit storage
estimates. Those two blocks are NESTED, so their difference is neither an
independent pair nor a scatter, and the figure is not the estimator's standard
error in any case; the standard error on a 10-orbit block is 0.197 W/m2, not
0.11. **The threshold is unchanged at 0.12 and rests on the temperature bound
alone**, which is where it was derived and which was fixed before any of this
was read: at the window now in use the measurable floor is 0.090 W/m2, below the
0.118 the offset tolerance implies, so the measurability half no longer binds.
`assess_convergence.py` now reports the storage row in `resolving_power` per run
and per window, which is what says whether a given verdict discriminated instead
of a sentence in a derivation saying it once.

## The derived relaxation time is not established as a ceiling

The offset criterion's fallback holds the relaxation time FIXED at what the
modelled mixed layer's heat capacity and the run's own radiative damping imply,
and holding it fixed is what makes the test cheap. That is conservative only
while the derived time is a CEILING on the true one: too long a tau inflates the
remaining offset and can only refuse a run. Three arguments were given for
expecting a ceiling and all three are arguments.

`exoplasim/scripts/check_relaxation_ceiling.py` takes the comparison the reports
have been accumulating. An artifact is evidence only when the exponential fit
was USED rather than the drift fallback, the fitted tau is finite and positive,
and it is no longer than the span it was fitted over. On an evidence artifact the ceiling is falsified when the fitted tau
exceeds the derived one by more than the fit's own standard error, and by any
amount where no such error is recorded. The unbracketed case counts AGAINST the
ceiling deliberately: the criterion's conservatism rests on the claim, so a
possible violation that cannot be dismissed is not a pass.

**Measured 2026-08-25 over all seven convergence artifacts. Three are evidence,
one falsifies, and none discriminates.**

| run | fitted | derived | evidence | verdict |
| --- | ---: | ---: | --- | --- |
| 2b20e3324bb0 | 169.13 | 9.883 | no, the drift fallback was taken | |
| 4182235e9781 | 34.78 | 9.883 | no, longer than the 26 orbits fitted | |
| 78c22fb1a1bd | -337429 | 9.883 | no, the drift fallback was taken | |
| 8044646ea7f0 | 13.68 | 9.883 | yes | ceiling falsified by 3.80 orbits, unbracketed |
| aaa95662e21a | -286755 | 9.883 | no, the drift fallback was taken | |
| ade7373b4c90 | 6.92 | 10.142 | yes | holds, unbracketed |
| ec32946bec89 | 6.27 | 10.120 | yes | holds |

**So the ceiling argument does not survive as stated, and it is not refuted
either.** The one artifact above the derived time is 38 per cent above it and
carries no error on its fit, so the excess cannot be attributed to the fit's own
noise or dismissed as it; its derived time is also the pre-correction 9.8833
built from an uncited density and specific heat pair. The two below it are below
by about the same margin as the one above is above, which is what a scatter
looks like rather than a bound. The fitted tau's own error is now recorded, and
on the one run where it exists it is 3.75 orbits against a 3.85-orbit gap: a
single run barely separates the two, which is the measurement this question
needs more of rather than a reason to prefer either answer.

What this does NOT license is treating the derived time as a floor, or widening
it to cover 13.68. The check accumulates; a run whose approach the fit can grip
adds a point, and a settled run adds none. Until it has points that
discriminate, the fallback's remaining offset should be read as an estimate
rather than the bound its comment claims, and any decision that turns on the
bound being safe needs the bracket rather than the derived number.

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
a floor. What would settle it is a production span at equilibrium of about
twenty times tau, which on the estimate above is about 80 orbits AFTER the
approach has finished; the section below measures where that number comes from
and why ten times tau is not enough.

## Fixing tau and bounding a mean are different lengths, and only one is affordable

*Measured 2026-08-25 on synthetic series whose tau is known in closed form, so
the right answer is set rather than compared against. Criteria fixed before the
run: tau is RECOVERED when the median estimate is within a tenth of the truth
and its relative RMS error is below a fifth; an interval COVERS when a nominal
95 per cent t-interval contains the true mean between 92.5 and 97.5 per cent of
the time. `lib/autocorrelation.py:integrated_time` is the estimator under test,
against an AR(1) and against a fast-plus-slow mixture that is deliberately not
AR(1).*

**Fixing tau by direct estimation is out of reach and always will be.** The
span needed scales as tau itself, so a slower process is not merely harder to
measure but harder in proportion to what makes it slow, and below that span the
estimate is biased low rather than merely noisy. The two processes have
different tau, so the same number of samples buys different multiples of it;
both columns give the estimate as a fraction of the truth, with the span in
multiples of that process's own tau in brackets:

| samples | AR(1), tau = 4.19 | mixture, tau = 10.43 |
| ---: | ---: | ---: |
| 20 | 0.56 (5) | 0.19 (2) |
| 85 | 0.85 (20) | 0.44 (8) |
| 300 | 0.93 (72) | 0.74 (29) |
| 1200 | 1.00 (286) | 0.91 (115) |
| 5000 | 1.01 (1192) | 0.97 (479) |

The AR(1) case first meets the criterion at 1200 samples and the mixture at
5000, which are about 290 and about 480 times their own tau. That agrees with
the analytic scaling for a truncated lag sum, where the span needed goes as tau
over the square of the tolerance. In orbits, at any tau this model plausibly
has, fixing tau means of order a thousand orbits AT EQUILIBRIUM. It is two
orders of magnitude beyond what this project runs and it is not a target to
design against.

**Bounding a window mean does not need tau fixed, and is affordable.**
Non-overlapping batch means take the scatter of the batch means themselves as
the error of the grand mean, which needs the batches long compared with tau but
never needs tau's value:

| batches | batch length, in tau | span, in tau | batch means covers | tau from inside the span covers |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 5 | 20 | 0.940 | 0.898 |
| 8 | 5 | 40 | 0.938 | 0.920 |
| 4 | 10 | 40 | 0.951 | 0.933 |
| 8 | 10 | 80 | 0.951 | 0.940 |

Batches of two times tau miss the criterion however many of them there are, so
the batch length and not the batch count is what has to be bought first. A span
of about twenty times tau is where an honest interval starts, and that is the
number to design a commissioning run against.

**What this means for a T21 commissioning run.** The approach is about seventy
orbits from cold, recorded above as experience. A span of twenty times tau
follows it, and tau is not measured: if it is near the value an AR(1) fitted to
a lag-1 of 0.615 implies, that span is about eighty orbits and the run is about
a hundred and fifty; if the model carries the slower component the mixture
stands in for, tau is larger and so is the span, in proportion. THE BRACKET IS
THE RESULT. What removes it is not a longer estimate of tau but a run long
enough to batch, which reports its own interval without needing tau at all.

**The twenty-orbit window is short by a factor of four even under the
optimistic branch**, which is the answer world-wdsk asked for. It is short in
the batch length, which is the term no number of extra windows repairs.
