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

**Measured 2026-08-25 on the dt-45 run's 35-orbit production window.** The heat
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

## What a criterion states: an upper bound, not a point estimate

Every convergence criterion in `assess_convergence.py` is

    |statistic| + its standard error < threshold

Five of the six tested `|statistic| < threshold` instead. That is a point
estimate, and a point estimate passes a state sitting exactly on its threshold
about half the time on the statistic's own scatter, which is a coin flip and not
a criterion. The sixth, the extrapolated offset, already tested an upper bound;
the asymmetry was defensible while nothing measured the other five statistics'
errors, and the resolving rows removed that excuse by measuring all six.

The verdict is what licenses the claim that a modelled run is equilibrated, and
every number this project publishes about the world stands on that claim. The
standing conventions all point one way: an estimate that cannot be verified is
bracketed and the bracket is reported, a convergence claim is labelled honestly
when it misses, and a gate that names a real conflict beats one that hides it.
The upper-bound form is that disposition written into the test.

**The thresholds did not move.** What changed is the form of the comparison and
not its tolerance. A threshold that looks wrong under the new form is a separate
finding and gets its own row rather than an adjustment.

**Why it was still honest to change the form after runs had been judged under
the other one.** A criterion chosen after the run it judges is not a criterion,
and adopting this one flips `run_ec32946bec89` from pass to fail. What dissolves
the objection is that nothing canonical has been run: the canonical climatology
lineage does not exist, so every build, run and climatology in the tree is
disposable whatever has consumed it, and that run is a stepping stone on a
binary already superseded. The flip therefore costs nothing that is being kept,
which is exactly the window in which a criterion's FORM can still be settled.
After the lineage is declared the same change could not be made honestly, which
is why the argument sits in the code beside the thresholds rather than in a
commit message.

**The resolving row prices the statistic the verdict was taken on.** The offset
criterion has two estimators -- the exponential fit where the series supports
one, the drift fallback where it does not -- and its row priced the fallback
unconditionally, so on a run where the fit was used the row described a
statistic the verdict had not taken and a reader could not tell which case a row
was. Both now come out of the branch that chooses the estimator, and every row
carries `statistic_error_source`; the offset row also carries
`statistic_source`, and `resolving_power.offset_statistic_source` repeats it
beside the window figures, which are the fallback form's arithmetic and price a
verdict only where the fallback was taken.

The two numbers a row reports are not the same question. `verdict_interval` is
what the criterion added to `|statistic|`; `statistic_standard_error` is what
the resolving bar is applied to. They are one number for the five direct
measurements. They differ for the offset criterion's fallback alone, where the
reported half width carries a stated 100 per cent allowance on the expected
relaxation time on top of the sampling error: the verdict must carry that
allowance, and the resolving question must not, because whether the window can
see 0.15 K is not a property of how far the run still has to travel.

**Re-assessed 2026-08-25 under the new form, every run whose output is on
disk.** One verdict moved.

| run | dt | orbits | window | before | after |
| --- | ---: | ---: | ---: | --- | --- |
| 14906cb7b914 | 30 | 35 | 35 | fails the offset | unchanged: fails the offset |
| ec32946bec89 | 45 | 50 | 35 | passes all six | fails the storage criterion |
| ade7373b4c90 | 45 | 28 | 28 | fails storage and offset, at 25 orbits and a 10-orbit window | fails all but the sea-ice slope |

`run_ec32946bec89` is the flip the decision was taken in spite of: its storage
is 0.0961 W/m2 with a standard error of 0.0303, so the upper bound is 0.1264
against 0.1200. Its offset row also changes verdict, from `resolves` false to
true, because the fit was used on that run and the row now prices the fit's own
error rather than the fallback's.

`run_ade7373b4c90` is a 28-orbit cold start, so it is shorter than the derived
window and was assessed at the longest window it supports. Five of its six rows
report that the window does not resolve their threshold, which is the honest
statement about a run a third of the way through its approach; the extra
failures are the run's length and not the change of form.

The five other artifacts in `exoplasim/analysis/convergence/` cannot be
re-assessed: their runs' output is gone, and they predate the standard errors
the new form needs, so nothing in them supports the comparison. They record
verdicts taken by earlier generations of the instrument, at a 10-orbit window,
and are superseded rather than merely old.

## The criteria against the size of the effect they have to see

The escalation route's first two rungs give the instrument something to be
measured against. T21 at dt 45 ran to a settled state and T21 at dt 30 followed
it from that restart, so any difference between their equilibria is a
timestep-dependent shift at fixed resolution, which is what alternating
resolution and timestep exists to expose.

**Measured 2026-08-26, and it is not yet a measurement of the shift.** The two
states differ by 0.103 K on the window means, 0.307 K on the fitted asymptotes
and 0.499 K on the last orbit. They differ that widely because NEITHER arm is
equilibrated by its own verdict. The dt-30 run passes the four slope criteria
and the storage criterion at a 35-orbit window and fails the extrapolated
offset, on its uncertainty rather than its magnitude. The dt-45 run passes the
other five and fails the storage criterion, by 0.006 W/m2 of a 0.12 W/m2
allowance, which is a fifth of that estimator's own standard error and
therefore a miss the window can barely see. So the shift is bracketed by roughly
0.1 to 0.5 K and the top of that bracket is contaminated by the approach both
runs have left to make.

**What the criteria can see, in each criterion's own units.** The mixed layer
drifts 0.128 K per orbit per W/m2 and the derived relaxation time is about ten
orbits, so a kelvin of remaining approach is worth about 0.79 W/m2 of storage.

| effect | as an offset, against 0.15 K | as storage, against 0.12 W/m2 | against the storage estimator's own 0.030 W/m2 |
| ---: | ---: | ---: | ---: |
| 0.10 K | 0.7x | 0.079 W/m2, 0.7x | 2.6x |
| 0.31 K | 2.0x | 0.243 W/m2, 2.0x | 8.1x |
| 0.50 K | 3.3x | 0.394 W/m2, 3.3x | 13.1x |

So both criteria can see a shift at the top of that bracket with about the
factor of three of margin the resolving bar asks for, and NEITHER can see one at
the bottom of it: a 0.1 K shift sits below both thresholds and is only 2.6 times
the storage estimator's standard error. That is the honest statement of what
this instrument settles about a timestep pair. Separating two equilibria is
`compare_equilibria.py`'s question rather than this one's, and it needs the
batch-mean interval the section below prices, not a convergence verdict on each
arm.

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

**A FITTED TAU IS ONLY A MEASUREMENT WHERE THE SERIES DETERMINES IT.** Over a
span short compared with tau, `exp(-n/tau)` is linear in n to within the fit's
own noise, so the exponential and a straight line are the same curve and tau is
whatever the optimiser drifted to. Five of the eight reports on disk carry such
a number, up to 44018 orbits and down to -337429, beside verdicts that were
correct: the criterion failed closed on the widened half width exactly as
designed while the artifact went on carrying a plausible-looking field. The
assessment now reports `relaxation_orbits_fitted` as null in that case, with
`relaxation_fit_identifiable`, the reason, and the only bound the data supports,
which is the span itself. The raw number is kept under
`relaxation_fit_raw_tau_orbits` so a degenerate fit is visible rather than
silent. No verdict moved: the raw tau still widens the asymptote's interval by
the overshoot it implies, which is how such a fit refuses through the criterion.

**The comparison is an artifact and not a table here.**
`exoplasim/analysis/convergence/relaxation_ceiling.json` carries every row, the
rule each was judged by, and the verdict over the set; the generator rewrites it
from the convergence reports and nothing else, so it costs nothing to re-take
and it moves the moment a report does. A table of rows in this note would be a
second copy of it with no way to learn that the first had changed, which is the
defect `notes/audits/frozen-derived-quantities.md` names and the reason
`lib/run_lengths.py` reads the artifact rather than declaring a bracket.

**A fitted tau is only a measurement where the series determines it.** Over a
span short compared with tau, `exp(-n/tau)` is linear in n to within the fit's
own noise, so the exponential and a straight line are the same curve and tau is
whatever the optimiser drifted to. Several reports carry such a number, up to
tens of thousands of orbits and negative on two runs, beside verdicts that were
correct: the criterion failed closed on the widened half width exactly as
designed while the artifact went on carrying a plausible-looking field. The
assessment now reports `relaxation_orbits_fitted` as null in that case, with
`relaxation_fit_identifiable`, the reason, and the only bound the data supports,
which is the span itself. The raw number is kept under
`relaxation_fit_raw_tau_orbits` so a degenerate fit is visible rather than
silent. No verdict moved: the raw tau still widens the asymptote's interval by
the overshoot it implies, which is how such a fit refuses through the criterion.

A degenerate fit is excluded twice over on some runs, and neither exclusion
depends on its size: the assessment took the drift fallback, so the tau is from
a fit it had already discarded, and the fit is degenerate by the rule above.
Read at face value the largest of them would refute the ceiling by four orders
of magnitude; read correctly it is not evidence in either direction.

**The ceiling argument does not survive as stated, and it is not refuted
either.** The artifact's verdict is what says which, and it is falsified today
by one artifact whose fit records no standard error: its excess over the derived
time cannot be attributed to the fit's own noise or dismissed as it, and its
derived time is the pre-correction figure built from an uncited density and
specific heat pair. The evidence rows sit on both sides of their derived times
by comparable margins, which is what a scatter looks like rather than a bound,
and the fitted tau's own error is comparable to the gap on every row that
records one. A single run barely separates the two, which is the measurement
this question needs more of rather than a reason to prefer either answer.

What this does NOT license is treating the derived time as a floor, or widening
it to cover the largest fit. The check accumulates; a run whose approach the fit
can grip adds a point, and a settled run adds none. Until it has points that
discriminate, the fallback's remaining offset should be read as an estimate
rather than the bound its comment claims, and any decision that turns on the
bound being safe needs the bracket rather than the derived number.

**The bracket and the ceiling are one question and are answered once.** The rows
that may test the ceiling are the rows that may bound the relaxation time --
both ask which fits are measurements of it -- so the generator emits
`fitted_bracket_orbits` over the same evidence rows, and `lib/run_lengths.py`
reads that rather than carrying a pair of its own. A fit whose standard error is
comparable to its value stays in the bracket: it cannot settle whether the
ceiling holds, which is what the row's `discriminates` field says, and it is
still a reading of how fast this model returns. Dropping it would narrow the
bracket because the reading is uncertain, which is backwards, and the bracket is
what a settling block is bought in.

## The memory sweep re-taken on a longer series, 2026-08-27

**Measured on `run_432e5e46adef` at 82 orbits**, by the rule the declaration was
taken under and unchanged: the per-orbit area-weighted mean surface temperature,
reduced by Geyer's initial monotone positive sequence over five candidate windows
starting at orbits 25, 30, 35, 40 and 45 and running to the end of the series.
The rule reproduces the declared bracket exactly on the first 70 orbits, which is
what makes the two rows below comparable rather than two different measurements.

| window starts at | 70-orbit series | | 82-orbit series | |
| ---: | ---: | ---: | ---: | ---: |
| | tau, orbits | lag-1 | tau, orbits | lag-1 |
| 25 | 2.2488 | 0.4696 | 7.8112 | 0.6670 |
| 30 | 1.8959 | 0.3767 | 8.5242 | 0.6639 |
| 35 | 1.9627 | 0.3988 | 8.3389 | 0.6782 |
| 40 | 2.0317 | 0.4241 | 7.7124 | 0.6791 |
| 45 | 2.1942 | 0.4920 | 6.9178 | 0.6862 |

**Every window of the 70-orbit series is reliable by
`autocorrelation.integrated_time`'s own span rule and none of the 82-orbit
series is**, the spans there being under ten tau at the larger value. So the
longer series does not simply replace the shorter one: it returns numbers its
own estimator declines to stand behind, and the declaration's admission rule was
that every window be reliable.

**The drift GREW rather than decayed.** Over the same five windows it went from
-0.0001 to +0.0068 K per orbit on the 70-orbit series and from +0.0032 to
+0.0087 on the 82-orbit one, and the lag-1 correlation went from about 0.4 to
about 0.67. A residual trend pushes every lag correlation up, so the second
column is an upper bound inflated by an approach that is still running, which is
the same argument the declaration carries for the first.

**The 12 added orbits are a `post_equilibrium_climatology` segment**, run at
`low_io` false where the spinup ahead of it ran at true, on
`most_plasim_t21_l10_p16.x`. They are production orbits by their declared
purpose and the assessment includes them, so this is one series and not two.

**What this settles and what it does not.** It settles that the declared bracket
does not bound this run's memory time as the run now stands: the top of it,
2.22, is exceeded by a factor of three. It does not settle what to declare
instead. The production span is twenty times whichever number is carried, so
adopting the second column would take a commissioning run at T21 from cold from
about 108 to 114 orbits to about 208 to 240, on estimates every one of which
fails the reliability rule the first column passed. What decides it is why a run
declared post-equilibrium is drifting faster at orbit 80 than at orbit 65, and
that question is prior to the number.

**Nothing is silently carrying the old pair, which is the part that is closed.**
`lib/run_lengths.py` records the observation the bracket was anchored to, and
`check_memory_bracket` refuses the moment `run_432e5e46adef`'s convergence report
is regenerated over the longer series -- both because the anchor moves and
because the report then reads a memory time above the bracket.
`assess_convergence.py`'s scatter bound refuses on the same regeneration, the
82-orbit window reading 0.0923 K against a declared 0.091.

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

## The two lengths a run is bought in, and the two taus they are bought in

Every number above is now consumed rather than only recorded.
`lib/run_lengths.py` turns the two timescales into the two lengths and keeps
them apart, because both are called tau and neither substitutes for the other:
the MEMORY time is how long the model's stationary wobble stays correlated with
itself and it prices a mean's interval; the RELAXATION time is the e-folding of
an approach and it prices a decay.

**A commissioning length is the approach plus twenty memory times.**
`lib/run_lengths.py` owns the memory-time bracket and the multiple, and this note
states neither: a span quoted here would be a second statement of a number that
moves when the sweep behind it is re-taken. The declared length is the FLOOR of
the bracket, and it is not a stopping rule: a run that reaches it and is still
refused by the convergence criteria is not finished.

**Running until the instrument stops refusing does not substitute for declaring
the length**, which is the alternative that had to be ruled in or out before the
derivation was worth building. It would be strictly better if the refusal were
the same bar. It is not: the convergence window is sized so the offset
criterion's slope error resolves its threshold, and the production span is sized
so a window mean's interval covers. The first is the smaller of the two by tens
of orbits at every point of the memory bracket, so the instruments stop refusing
well before the climatology's own interval is bought, and the refusal is a
weaker bar that can only decide whether to buy MORE.

**A settling length is `tau_relaxation * ln(A / 0.15 K)`**, where A is the
perturbation the step change made. The residual is the offset criterion's own
allowance, so what is left is smaller than the instrument that judges the next
state can see. The relaxation time is the bracket the fits support and not the
derived value, because the derived value is not established as a bound on them.
A step change worth about half a kelvin settles inside the ten to twenty orbits
this project has repeatedly seen; `scripts/smoke_test.py` fails if the
derivation stops covering that experience, which is a check on the derivation
and not a target it was fitted to.

**THE THREE NUMBERS THAT SIZE EVERY RUN ARE HELD TO THE REPORTS.** The window's
orbit scatter, its memory time and its relaxation time were declared constants
that no artifact could contradict, and two of the three were once wrong by about
a factor of two while every commissioning run on the ladder was bought in them.
Their disposition now:

| quantity | where it lives | what re-runs when the measurement moves |
| --- | --- | --- |
| the memory time | `lib/run_lengths.py`'s bracket, and nowhere else. `assess_convergence.py` sizes its window on the top of it | `check_memory_bracket` refuses when the anchor report moves or when any settled production report reads above the bracket |
| the orbit scatter | `assess_convergence.py`, with the report it was read off | `check_convergence_bounds` refuses on the same three conditions |
| the relaxation time the window is sized on | `assess_convergence.py`, bounding what each run derives | `check_convergence_bounds` refuses when any report derives more |
| the relaxation-time bracket | nowhere. It is read from `relaxation_ceiling.json` | the generator rewrites the artifact and the bracket follows |
| the fit's tail fraction | `assess_convergence.py`, once. It is a decision, not a measurement | nothing has to: there is no second copy to disagree |

All three declared values are UPPER BOUNDS taken on a series that still drifts,
so a report reading below them is what they predict and is not a defect. What
was missing is that nothing could tell that from a bound nobody had looked at
since the series moved, and the ANCHOR is what separates them: each declaration
records what the report it was taken beside reads, and a report that reads
anything else refuses the declaration whether or not the bound still holds. The
circle the audit named -- the constants size the window, the window sizes the
assessment, the assessment reports the constants, and nothing carries back -- is
closed at that refusal rather than by feeding the measurement back in, because
the measurement is a LOWER bound on the same quantity and substituting it would
shorten every span on the ladder.

A settling length establishes only that a restart is not mid-transient. It does
not establish equilibrium, a rung's climate, or an interval on any mean taken
over it.

## The 82-orbit tau is a step, not a memory

*Measured 2026-08-27 on `run_432e5e46adef`, whose 70 low-I/O orbits were
followed by a 12-orbit clean-I/O segment for the climatology.*

Fitting the whole 82-orbit series returns a memory time of 9.08 orbits, four
times the declared bound, on a fit the estimator itself marks unreliable.
Fitting each I/O regime separately returns something else entirely:

| block | drift, K/orbit | lag-1 | tau | reliable | stationary |
| --- | ---: | ---: | ---: | --- | --- |
| low-I/O, orbits 35 to 69 | +0.0037 | 0.440 | 2.00 | yes | no |
| clean-I/O, orbits 70 to 81 | -0.0049 | -0.080 | 1.00 | yes | YES |
| across the join | +0.0065 | 0.715 | 9.08 | NO | no |

**There is a step of +0.1616 K at the regime change**, and a step is a perfect
long correlation. It is what lifts the lag-1 from 0.440 to 0.715 and the memory
time from 2.00 to 9.08; the series has not grown a slow component, it has
acquired a discontinuity.

The step is not a defect in the run. A low-I/O orbit holds interval
ACCUMULATIONS and a clean-I/O orbit holds instantaneous samples, so the two
blocks report different quantities and a mean taken across them is a mean of two
measurements. `exoplasim/notes/first-output-bin.md` is where that difference is
argued.

**Within each regime the run is settled**, and the clean block is the better
evidence: it drifts NEGATIVELY, its lag-1 is indistinguishable from zero, and it
is the only block in this run that passes `autocorrelation.stationary_enough`
outright. Twelve orbits supports a memory time up to 1.2 under the
ten-tau reliability rule, so what it establishes is an upper bound of about 1.2
rather than a value.

**What this costs the ladder, which is the reason it was chased.** The production
span is twenty times the memory time. At the mixed reading it would be 180
orbits at every commissioning rung; at the low-I/O reading 40; at the clean-I/O
bound 24. The mixed reading is the one that would have been adopted by a fit over
"the whole run", and it is wrong by the width of the cost argument.

**The operational rule that follows: never fit a memory time across a change of
I/O regime.** The regimes are declared per segment and a run carries both by
design, since the cheap regime is what a spin-up is bought in and the clean one
is what a climatology is read from.
