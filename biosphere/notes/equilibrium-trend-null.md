# The equilibrium acceptance rule, and the instrument it took four measurements to find

Vesper is a simulated world. This note is about the rule that decides whether a
LPJ-GUESS run of its simulated biosphere has settled, and about four
measurements that each changed what the rule is rather than what a number in it
says. Measured on 2026-08-30 and 2026-08-31 against
`biosphere/config/equilibrium_window.yaml`, which is the declaration; this note
is the argument.

## What the rule is

A QUANTITY A CONSUMER READS is SETTLED when the upper confidence bound on its
end-to-end relative drift, taken over the whole retained record, is inside the
tolerance that quantity's own consumer owes. A run is accepted when every
assessed quantity is settled and none of their per-cell halves refuses.

WHAT THE ACCEPTANCE ARTIFACT CLAIMS IS THEREFORE NARROWER THAN "THE SIMULATED
BIOSPHERE HAS SETTLED", and it says so in those words. Every column of every
stability table keeps a drift bound on the artifact as a diagnostic, so what
the claim leaves out is visible rather than absent.

What the contract declares is `maximum_false_acceptance_rate`: the probability
that a quantity whose true drift is AT its limit is nonetheless accepted. That
is the error an acceptance gate has to control, because an acceptance gate
asserts that a run HAS settled. It is a per-quantity rate and it is also the
run-level rate, with no multiplicity correction: a run passes only when every
quantity passes, so by the intersection-union principle the probability of
accepting a run in which any quantity truly drifts at its limit is bounded by
the per-quantity level. The rate at which a SETTLED run is refused is not
declared. It is the instrument's cost, it is reported per quantity, and it is
what sizes the retained record.

Reproduce the statistic's size and cost with
`biosphere/scripts/validate_drift_statistic.py <run>`, and the timescales with
`biosphere/scripts/derive_trend_null.py <run> --timescales`.

## The first measurement: no single trending-cell fraction can be a limit

The rule began as a slope over the last 10 complete forcing cycles, with a run
refused when the spatial mean drifted more than 5 per cent end to end at over 2
standard errors, or when more than `maximum_trending_cell_fraction: 0.25` of
cells were individually flagged. The 0.25 had no derivation, and neither knob
that should have relaxed it moved it: across 100 to 1000 orbits the flagged
fraction went 0.376 to 0.365, and across npatch 5 to 20 it went 0.376 to 0.348
while the median per-cell drift fell to 0.603 of its value and the median
relative slope standard error fell to 0.658 of its value. The test is the ratio
of those two, so reducing the noise cannot satisfy it.

THE NULL COSTS NO MODEL RUN. The forcing cycle is one simulation year:
`run_manifest.forcing` records `cycle_years: 1` against a `VESPDRV8` driver of
12 intervals, and every simulated year replays byte-identical forcing. All
interannual variation in a retained record is therefore internal stochasticity
of the simulated vegetation, and a long run at fixed forcing is already the
two-seed experiment the null was thought to need. Detrend each cell's whole
record, add its mean back, cut disjoint 10-cycle windows, and apply the
contract's own test to each.

Per-window maximum over the assessed fields, 100 windows of the 1000-orbit
diagnostic `lpj_377476c097b549fcaaa07cc7644275ad`:

    median 0.392    p90 0.450    p95 0.455    p99 0.461    maximum 0.462

The same statistic on `lpj_36b0f5fc91a640ed94618b0e95420916` (npatch 20, 10
windows) gives median 0.367 and maximum 0.383, so the null is near enough
invariant to patch count as the ratio argument predicts, and the flag rate
against window position is flat (first-half median 0.394, second-half 0.388), so
the detrending left no residual drift in the estimate.

The limit therefore sat BELOW the rate at which the test flags cells with
provably no trend to find, and the null is a per-field quantity spanning three
orders of magnitude because it is set by each field's own internal variability:

| field | occupancy | null | power at 2x limit |
| --- | --- | --- | --- |
| npool.out SoilN | 1.00 | 0.001 | 0.998 |
| npool.out Total | 1.00 | 0.001 | 1.000 |
| cpool.out SoilC | 0.94 | 0.001 | 0.942 |
| cpool.out Total | 0.94 | 0.058 | 0.907 |
| tot_runoff.out Total | 1.00 | 0.133 | 0.846 |
| lai.out Total | 0.93 | 0.298 | 0.552 |
| lai.out C3G | 0.93 | 0.379 | 0.539 |
| cpool.out VegC | 0.90 | 0.387 | 0.605 |
| lai.out TeNE | 0.10 | 0.066 | 0.085 |
| lai.out C4G | 0.15 | 0.075 | 0.099 |
| lai.out IBS | 0.20 | 0.118 | 0.135 |

Two defects are visible in that table. The DENOMINATOR was every cell rather
than the cells the field occupies, so a plant functional type present on a sixth
of the grid could not reach a quarter-of-all-cells limit however hard it
drifted; nine of the thirteen types occupy under 0.37 of cells and the guard was
inert for them. And the SPREAD is irreducible by any single choice: a limit high
enough for the patch-driven grass destroys it for the slow soil pools, where it
discriminates almost perfectly.

The denominator became the occupied cells and the limit became a measurement
from the run's own detrended windows. That construction is what the per-cell
half still uses.

## The second measurement: an empirical null cannot control a family rate

A stationary field exceeds the largest of N detrended null windows with
probability 1/(N+1), which is exact PER FIELD. The gate refuses when ANY of the
64 assessed fields exceeds, and the family rate is `1-(1-1/(N+1))^F_eff`.
Measured by leave-one-out on the 1000-orbit run's 100 windows it is 0.34 against
a declared 0.05, flat across `slope_standard_errors` from 1.0 to 4.0 as the
construction requires. The helper was validated first on synthetic fields: at
100 windows it returns 0.010 for one field, 0.080 for eight and 0.460 for
sixty-four independent ones, against an analytic 0.010, 0.077 and 0.474. The
measured 0.34 corresponds to about 35 effective independent fields out of 64,
the remainder being the correlated ones -- C3G in lai, fpc and aaet move
together.

Raising N does not fix it. The finest probability an N-window empirical null can
express is 1/(N+1), so a family rate of 0.05 over about 35 effective fields
needs roughly 710 windows, 7100 retained cycles, about twenty hours of model
time, re-measured per configuration. THE REPAIR IS AN ANALYTIC NULL, which has
unlimited resolution, and an analytic null needs a standard error that is valid.

## The third measurement: the window was shorter than the thing it tested

Through `lib/autocorrelation.py`, which owns the estimator and the two verdicts
that say whether its answer means anything. On the 1000-cycle record, 400
sampled cells per table:

The integrated autocorrelation time of the spatial-mean cycle series runs 9.3 to
125.3 complete forcing cycles. The window the contract tested over was 10.
Sixty-one of the 64 assessed fields have a spatial-mean memory time longer than
the whole window, per-cell medians run 2.9 to 93.1, and the share of cells whose
own memory time exceeds the window is 0.22 to 1.00.

| table.field | spatial tau | cell tau median | cells with tau > window |
| --- | --- | --- | --- |
| lai.out BNE | 75.3 | 64.9 | 1.00 |
| lai.out TrIBE | 93.5 | 63.9 | 0.97 |
| lai.out C3G | 9.9 | 6.3 | 0.36 |
| cpool.out Total | 125.3 | 68.1 | 0.89 |
| npool.out SoilN | 60.6 | 89.0 | 0.87 |
| tot_runoff.out Total | 14.9 | 5.5 | 0.33 |

A trend fitted inside one memory time is not a trend. It is one smooth excursion
of a process that has not had time to sample its own distribution; the residuals
within the span are small because the process is smooth on that scale, so the
ordinary least-squares standard error is small and the slope looks decisive.
That is the whole of the measured 0.47 to 0.86 per-cell flag rate against a
nominal 0.081, and it is `docs/src/practice/failure-modes.md` class 34: a bed
shorter than its startup.

NO CORRECTION RESCUES THE WINDOWED FORM. Inflating the slope standard error by
the square root of the memory time moves lai.out's bare flag rate from 0.47-0.84
to 0.19-0.46: two to three times better, still two to six times nominal, still
spanning a factor of four across fields. That correction is for a MEAN over a
span and cannot make a span shorter than its memory carry a trend. A
within-window circular-shift surrogate, which keeps each cell's marginal
distribution but whitens its serial structure, returns 0.007 to 0.09 across
fields, close to the independent-samples rate: the whole excess is temporal
memory.

## The fourth measurement: a span multiple does not converge

The obvious repair is a longer record and a longer window, at the span bar
`lib/autocorrelation.py` declares: a span shorter than `RELIABLE_SPAN_MULTIPLE`
times its own memory time cannot establish that memory time. Applied to the top
of the measured bracket that asks for 1253 retained cycles, and a spin-up of
2318 cycles follows from the relaxation time. That run was bought:
`lpj_7d3c576ee4e342acb097b46fece976e0`, npatch 5, 1617 cells, 16 ranks, about 55
minutes as derived.

THE FLOOR MOVED UNDER ITSELF.

    read on the 1000-cycle diagnostic    9.3 to 125.3 cycles  ->  record floor 1253
    read on the 1253-cycle record       11.1 to 212.7 cycles  ->  record floor 2127

The top rose 70 per cent when the window it is read on grew 25 per cent. The
slowest fields are cpool.out Total 212.7, cpool.out VegC 184.1, lai.out TrBR
173.2, anpp.out TrBR 169.8, aaet.out TrBR 162.5, and the estimator marks six of
those seven NOT reliable, so 212.7 is itself a floor.

A record sized at ten times a memory time read off a shorter record is not
self-consistent, and it has now failed to converge twice. That is not an
argument for a third guess at a length. It is an argument that the span multiple
is the wrong criterion for this question: it asks whether the memory time can be
ESTABLISHED, when what the contract needs to know is whether the TOLERANCE can
be RESOLVED.

## The instrument: an upper bound on the drift

The direction is the repair. A significance test that fails to reject "no drift"
asserts nothing about equilibrium; it reports that it could not tell, and on a
record shorter than its own memory time it can never tell about anything. An
acceptance gate's own error is letting a DRIFTING run through, and the
instrument for that is an equivalence test. `lib/lpj_output.py:drift_bound`, per
field, on the spatial-mean cycle series `x[0..n-1]`:

    m = n // 2;  H1 = x[0:m], H2 = x[n-m:n]
    tau  = max over the two halves of integrated_time(H).tau
    v_i  = var(H_i, ddof=1) * tau / m
    D    = 2 * (mean(H2) - mean(H1))
    seD  = 2 * sqrt(v1 + v2)
    df   = Welch on v1, v2 at m / tau - 1 degrees of freedom each
    U    = |D| / scale + t(alpha, df) * seD / scale

and the field passes when `U <= relative_end_to_end_limit`.

Three properties follow from the direction alone, and each replaces a defect the
windowed form could not repair.

THE MULTIPLICITY DISAPPEARS, by intersection-union, as above. The Sidak factor,
the counted family size and the 0.34 family rate are all gone by a derivation
rather than by a longer record.

THE GUARD DISAPPEARS. A record too short to resolve the tolerance gives a wide
bound, the bound exceeds the limit, and the field is refused by the test itself.
There is no second span rule standing in front of the test, and being a span
rule is what made the guard self-referential.

THE RECORD FLOOR CONVERGES. `lib/lpj_output.py:cycles_for_bound` inverts the
bound: the record at which a settled field of this scatter and this memory time
closes inside the limit. The standard error falls as one over the root of the
record while the memory time grows sublinearly with it, so the requirement is
reached rather than chased.

Four details are load-bearing and each was found by measurement. The
half-to-half difference is HALF the end-to-end change of a steady drift and is
DOUBLED; left unscaled it silently doubles the tolerance it is judged against.
The span bar belongs on the half whose mean is taken, and it enters as
degrees of freedom rather than as a hard bar -- applied to the whole record
instead, a half stands on five effective samples and a family refusal rate of
0.23 was measured against a declared 0.05. Scatter and memory time are taken
on the RAW half rather than a detrended one, so a real drift inflates the error
it is judged against rather than shrinking it. And the memory time is taken at
the UPPER END of its own sampling interval rather than at the estimate, which
the measurement below forced.

## What the statistic measures, against series whose answer is known

2000 synthetic AR(1) trials per cell, 1253 cycles -- the longest record this
project has -- with a relative scatter of 0.02 and the memory time swept over
and beyond the 11.1 to 212.7 the model's own fields read. Reproduce with
`biosphere/scripts/validate_drift_statistic.py <run>`; the JSON is
`biosphere/analysis/drift_statistic_validation.json`.

THE RATE THE CONTRACT DECLARES is the share of trials that ACCEPT a field whose
true end-to-end drift is exactly the tolerance. Declared 0.05, and the bar fixed
before the sweep was 0.075 at every swept memory time:

| memory time | 1 | 10 | 20 | 40 | 80 | 125 | 175 | 213 | 300 | 400 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| at the estimated memory time | 0.000 | 0.012 | 0.020 | 0.037 | 0.051 | 0.053 | 0.068 | 0.080 | 0.099 | 0.109 |
| at its upper end, which ships | 0.000 | 0.004 | 0.005 | 0.011 | 0.018 | 0.020 | 0.027 | 0.028 | 0.031 | 0.033 |

The first row MISSES the bar from 213 upward, which is the top of the range the
model is in, and the cause is in the estimator rather than in the statistic. On
1253 samples `integrated_time` returns 1.01 for a true 1.0, 9.70 for 10.0, 34.4
for 40.0, 86.1 for 125.0 and 128.9 for 213.0. A memory time read too small makes
every standard error built on it too small, and the bias grows exactly where the
acceptance rate does. The second row is what ships, and it holds the declared
rate across the whole swept range with room to spare.

The upper end is the Madras-Sokal windowing interval, and it is validated rather
than taken on trust: it covers the memory time the series was built with 1.000,
0.775, 0.760, 0.585 and 0.545 of the time at 1, 10, 40, 125 and 213. Half to
three-quarters is NOT the 0.84 a one-standard-error upper bound would nominally
give, because the interval is symmetric about a centre that is biased low. It is
reported rather than rounded up, and what carries the declared acceptance rate
across the range is that conservatism together with the raw-half memory time and
the ignored covariance between the halves.

THE COST, which is not barred and is what sizes the record: the share of trials
that REFUSE a field with NO drift at all, at that same 0.02 relative scatter.

| memory time | 1 | 10 | 20 | 40 | 80 | 125 | 175 | 213 | 300 | 400 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| refused with no drift | 0.000 | 0.000 | 0.020 | 0.273 | 0.655 | 0.809 | 0.860 | 0.878 | 0.882 | 0.887 |

THE SIGNIFICANCE FORM, measured on the same trials for comparison, is the one
this row began with: reject "no drift" at a per-field level derived from a
declared family rate of 0.05 by Sidak over the 64 counted fields. Quoted at the
ESTIMATED memory time, which is the form as it was proposed and its most
favourable setting, its size is excellent -- 0.000 to 0.027 with no drift --
and its power collapses exactly where it matters. At a true drift AT the
tolerance it refuses 0.510, 0.446, 0.349, 0.207, 0.115, 0.089, 0.086, 0.071,
0.075, 0.068 across the same sweep, so on a field with the memory time this
model's slowest have it would PASS a run drifting at the tolerance about
ninety-three times in a hundred. That is not a defect of the multiplicity
correction; it is what a "failed to reject" verdict means on a record this
short, and it is the whole argument for the direction.

ONE REGISTERED CHECK MISSES AND THE REASON IS UNDERSTOOD. At a memory time of
one the construction was to reproduce the textbook two-sample t and accept a
drift at the tolerance 0.05 of the time. It accepts 0.000. The cause is the
raw-half memory time: an injected ramp inside a half reads as memory, so the
standard error comes out four times the independent-samples one and the bound is
far wider than nominal. That is the conservatism working as designed, and it is
recorded as a miss rather than explained away. The half of that check that does
have a clean answer -- the estimator recovering the memory time it was handed --
passes and is quoted above.

## The reported value is taken over the span the contract certifies

The contract certifies the whole retained record. It used to report a mean over
the last 10 cycles, and that mean is what `build_surface_albedo.py`,
`score_prediction.py` and `build_soil.py` consume. At a memory time of tens to
hundreds of cycles a ten-cycle mean is one effective sample and cannot be known
better than the field's marginal scatter, and ten of the assessed fields have a
marginal scatter above `relative_end_to_end_limit` outright. The gate would have
been certifying a state to 5 per cent while handing the consumer a number that
is not known to 5 per cent: class 34 again, pointed at the value rather than at
the test.

The reduced value, its temporal spread and its own relative standard error are
now taken over the whole retained record, and the standard error travels with
the value on the report. `complete_forcing_cycles` names the per-cell half's
window and nothing else.

## What a 5 per cent drift is worth to the consumer

The tolerance's own derivation, through the consumer that closes a loop:
`exoplasim/scripts/build_surface_albedo.py --mode modelled` reads `fpc.out` and
mixes tree and grass cover into the surface albedo, so the lever is the rootable
fraction times canopy albedo minus substrate albedo.

| link | value | source |
| --- | --- | --- |
| tree, grass albedo | 0.143 [0.130, 0.148], 0.209 [0.190, 0.217] | `config/planet.yaml` |
| substrate albedo | 0.25585 land mean bare rock | `exoplasim/inputs/t21/albedo_report.json` |
| rootable fraction of model land | 0.83488 | `biosphere/generated/vesper_driver_provenance.json` |
| land fraction | 0.432841 | build manifest `landFractionBySurfaceClass` |
| surface to planetary attenuation | 0.5, DECLARED, factor of two | `scripts/error_budget.py` |
| flux to kelvin | 159.7 K per flux ratio [155.3, 160.6] | `lib/sensitivity.py` |

A 5 per cent drift in both tree and grass cover is worth about +0.064 K of
surface temperature, bracketed +0.032 to +0.128 K. Against what the climate arm
resolves -- 0.15 K for the single-run convergence criterion, 0.05 K for its
resolving-power bar, 0.027 K for the standard error of a 57-orbit climatology
mean -- the derived limit spans 0.02 to 0.23 and cannot choose the number. The
bracket is dominated by one unmeasured link: there is no measured
d(planetary albedo)/d(surface albedo) in this tree, `lib/sensitivity.py`
explicitly declines to own it, `scripts/error_budget.py` declares 0.5 and says
in the same breath that the honest claim is a factor of two, and the single
measured point is 0.3045 from one paired arm. world-ckbt is open to measure it.
Until it closes, 0.05 stays a declared tolerance with this bracket recorded.

## What each quantity owes, and the twenty-seven-fold it was costing

THAT DERIVATION IS ABOUT AGGREGATE COVER, and the contract used to apply its
number to all 64 columns of the seven stability tables, each plant functional
type's own column included. The binding column was `anpp.out` TrBR at 33846
cycles -- about eight and three quarter hours of retained record before any
spin-up -- and no consumer reads it. The consumer quantities it was buying that
record for resolve on the 1253 cycles already on disk.

THE RULE, and it comes from the consumer list rather than from any measurement.
`require_lpj_acceptance` has exactly three callers, and the acceptance contract's
own closure block is a fourth reader.

  1. The assessed set is the set of consumer quantities. A quantity is a table
     plus the columns a consumer SUMS.
  2. A COLUMN THAT ENTERS A CONSUMER QUANTITY ONLY THROUGH A SUM OWES NOTHING
     BEYOND WHAT THE SUM OWES. The sum's own cycle series carries that column's
     drift at its own share, with every cancellation and every reinforcement
     already in it, and the drift bound applied to the sum bounds the sum
     exactly. Bounding each column separately at the sum's tolerance is a
     strictly stronger claim than any consumer makes.
  3. A quantity's tolerance is the TIGHTEST among the consumers that read it,
     each derived through that consumer's own chain.
  4. A column no consumer reads is not assessed. It keeps a reported drift bound
     as a diagnostic, and the artifact claims what it measured.

| consumer | quantity | tolerance | derivation |
| --- | --- | --- | --- |
| `build_surface_albedo.py --mode modelled` | `fpc.out` tree cover, grass cover | 0.05 | the albedo chain above, which is ABOUT these two |
| `build_soil.py` | `cpool.out` SoilC | 0.05 | the tightest derived in this tree, applied conservatively to an unpriced chain |
| `lpj_acceptance.yaml` closure | `cpool.out`, `npool.out`, `aaet.out`, `tot_runoff.out` Totals | 0.05 | the same, and the closure imposes no drift tolerance of its own |
| `score_prediction.py` | `anpp.out` Total, `lai.out` Total, `fpc.out` tree, grass, C4G and boreal cover | 0.2727 | tightest half-width over centre of the ten pre-registered bands |

The albedo chain's number is the one derived above and it applies to the two
quantities it was derived for. `build_soil.py`'s own chain -- soil carbon to
organic mass fraction to volumetric water capacity to the model's `dwmax` --
has not been priced, so SoilC and the closure's four pool and flux totals take
the tightest tolerance derived anywhere in this tree. That is conservative and
is recorded as conservative rather than as a derivation.

`score_prediction.py`'s tolerance IS derived, from bands fixed in
`notes/productivity-prediction.md` before the model ran a Vesper gridcell: a
drift must not carry a value at a band's centre out of the band, which is the
band's half-width over its centre. The tightest of the ten is 0.2727, set by
line 3 (1.6 to 2.8) and line 5 (0.40 to 0.70); the loosest is 1.0. Reading the
drift against the record MEAN rather than against any point in the record would
double it, and the stricter reading is taken.

ADDING A CONSUMER ADDS A ROW to `assessed.quantities`. That is the whole of the
maintenance rule, and it is what makes the assessed set a derived thing rather
than a list: a new reader that forms a sum the block does not name is a sum this
contract has not certified.

`equilibrium_window_gate.py` tests the rule against a pair that can fail either
way. One fixture's two columns drift in opposite directions and one fixture's
drift together, at the same size; judged on the sum a consumer forms the first
is accepted and the second refused, and judged column by column both are
refused. Same statistic, same tolerance, opposite verdicts.

## What is still calibrated

The PER-CELL half. It exists for the one case a spatial mean cannot see, and
that job is real: measured on the 1000-cycle record's 64 columns, a cancelling
dipole at twice the tolerance is caught by the cell half in 42 of them and by
the global half in none.

| acceptance window | cell half refuses | global half refuses |
| --- | --- | --- |
| as run | 0 of 64 | 1 of 64 |
| coherent drift at 2x the tolerance | 46 of 64 | 63 of 64 |
| cancelling regional drift at 2x the tolerance | 42 of 64 | 0 of 64 |

It runs on the assessed quantities, at each quantity's own tolerance, exactly as
the global half does. It is still a slope over `complete_forcing_cycles`, so it
still sits inside one memory time, and its size is still absorbed into a limit
measured from the run's own detrended windows rather than derived. Its family
false-refusal rate is 0.34 on 100 windows and about 0.24 on the 124 a 1253-cycle
record leaves; over eleven assessed quantities rather than 64 columns it is
smaller than either, and it is still not derived. That errs toward REFUSING, so
a pass through it is evidence and a refusal through it is not, and its
empirical-null construction needs many windows and therefore cannot adopt the
whole-record statistic. world-4hlw holds the analytic replacement.

## How long a run this world's ecology needs

Both lengths are FLOORS, both are read out of the acceptance artifacts by
`lib/run_lengths.py`, and that module states no number of its own. They are a
second timescale pair, in COMPLETE FORCING CYCLES rather than orbits, and the
module keeps them apart from the climate model's: a cycle is one simulation year
is one modelled orbit, so the two units coincide numerically only while the
driver's cycle is one year, and the run's own declared cycle length travels with
every reading.

### The retained record

The record floor is `cycles_for_bound`'s answer for the assessed quantity that
needs the most, and it is written onto the acceptance artifact whatever the
verdict -- because a run refused for not resolving its own drift is exactly the
run that says how long the next one has to be. It is a floor for one reason
only: the memory time is held at the value its own record read, and a longer
record may read a larger one. That moves the answer rather than preventing one,
which is the difference between it and the span multiple it replaced.

WHAT IT IS SET BY IS SCATTER, NOT MEMORY TIME, and that is the finding the
earlier form could not see. The record a series needs is proportional to its
variance times its memory time, so `cpool.out` SoilC, at a memory time of 483.5
and a relative scatter of 0.0030, resolves its tolerance on the 1253 cycles
already on disk, while `anpp.out` TrBR, at 324.8 and 0.1012, does not until
33846.

THE ASSESSED QUANTITIES ARE MUCH QUIETER THAN THE COLUMNS THEY ARE BUILT FROM,
and that is the finding rather than an assumption. A sum over the plant
functional types averages out the competitive shuffling between them, and the
shuffling is most of each column's variance. On the 1253-cycle record the
relative scatter of `fpc.out` tree cover about its own linear fit is 0.00807,
against 0.0220 to 0.0290 for the individual woody columns it sums; the boreal
sum is 0.0167 against 0.0220 to 0.0267 for its three. The sum is quieter than
any of its parts, which is what a competitive shift looks like from the outside.
`biosphere/analysis/drift_statistic_validation.json` prices both sets.

Measured on the same record:

| series | drift bound | needs |
| --- | --- | --- |
| fpc.out tree cover, the albedo consumer's own quantity | 0.0186 | inside 1253 |
| fpc.out grass cover, the albedo consumer's own quantity | 0.0420 | inside 1253 |
| cpool.out SoilC, the soil consumer's quantity and loop B's exit | 0.0153 | inside 1253 |
| every other assessed quantity | 0.0031 to 0.0744 against 0.05 or 0.2727 | inside 1253 |
| cpool.out Total, the carbon closure's pool | 0.0503 against 0.05 | inside 1253 |
| anpp.out TrBR, a column no consumer reads | 0.4758 | 33846 |
| aaet.out IBS, a column no consumer reads | 0.3105 | 26128 |

Every assessed quantity resolves its own tolerance on the record on disk, and
every column that does not is an individual plant functional type's. So the
record floor is the record on disk, and the eight and three quarter hours the
old assessed set asked for is not owed.

ONE ASSESSED QUANTITY IS STILL REFUSED, and it is refused for drifting rather
than for being unresolved: `cpool.out` Total bounds at 0.0503 against 0.05, on
a drift estimate of 0.0113 that a settled series of its scatter would resolve
inside the record it already has. The drift lives in `cpool.out` VegC, which
bounds at 0.1495 -- vegetation carbon is still redistributing between the woody
types while the COVER those types present has settled. That is a refusal a
longer record does not answer and a longer spin-up does, which is the
distinction the two floors exist to keep apart.

### The relaxation time, measured without its asymptote

Fitting `a + b exp(-t/tau)` needs the record to contain the turn-over, and on the
1000-cycle record the asymptote landed outside the data for 34 of 64 fields.
`lib/lpj_output.py:relaxation_time` does not need it: for that same exponential
the difference between consecutive equal blocks decays by `exp(-Q/tau)` and `a`
cancels out of the ratio, so four equal blocks give `tau = -Q / ln(r)`.

Its admissibility conditions were declared before it was run, and one of them
was added after the estimator was measured against itself: the successive
differences must share a sign, their ratio must be a contraction in (0, 1), each
difference must exceed twice its own MEMORY-CORRECTED standard error, the
contraction must sit two standard errors BELOW ONE, and the second ratio must
agree with the first inside a factor of two.

THE RESOLUTION CONDITION IS THE ONE THAT BITES. `tau = -Q / ln(r)` diverges as
`r` approaches one, so a ratio whose own uncertainty reaches one is a record with
no curvature in it and the timescale it implies is a lower bound rather than a
value. On the 1253-cycle record the model's simulated soil nitrogen returned a
contraction of 0.992 with a delta-method uncertainty of 0.344 and read out as an
e-folding time of 37647 cycles; a spin-up sized from that number would have been
about 74000 cycles against the 2500 the same record supports.

AND NO FIELD ON ANY RECORD THIS PROJECT HAS SURVIVES IT. Over 192 field-records
across the three diagnostic runs, 97 are declined for block differences that
change sign, 82 for a difference inside its own standard error, five for a ratio
that is not a contraction at all, and every one of the eight that reaches a
contraction has an uncertainty reaching one -- the best being 0.584 +/- 0.331.
The bracket of 59.7 to 846.8 cycles that the 2318-cycle spin-up was derived from
rested on ratios whose interval spanned the whole answer.

So THE RELAXATION TIME OF THIS MODEL IS NOT MEASURED. It does not have to be,
and the section below is why. The RECORD floor never depended on it; the two
are independent and one being unmeasurable was hiding the other.

### The spin-up, derived from the tolerance and nothing else

The spin-up floor takes its residual from the criterion that judges what
follows, exactly as `SETTLING_RESIDUAL_K` takes 0.15 K from the offset criterion
on the climate side. A spin-up from bare ground starts a full equilibrium level
away, so the approach remaining after `S` cycles is `exp(-S/tau)` of the level
and the drift it puts across a retained record of `L` is that times
`1 - exp(-L/tau)`. The acceptance contract refuses a record whose relative
end-to-end change exceeds `relative_end_to_end_limit`, so

    exp(-S/tau) x (1 - exp(-L/tau)) <= 0.05

and `ecological_spinup_cycles` returns the smallest `S` that satisfies it.
Starting from bare ground is the conservative reading and is stated rather than
hidden: a pool the CENTURY accelerator hands over part-grown begins closer than a
full level away and needs less.

THAT REQUIREMENT IS BOUNDED OVER ALL RELAXATION TIMES, which is what makes the
unmeasurable `tau` stop mattering. Substituting `u = L / tau`,

    S(tau) = L x ln((1 - exp(-u)) / 0.05) / u

falls away at both ends: an approach far slower than the record puts almost none
of itself into the record, and one far faster has finished. So it has a finite
maximum, and because the left-hand side of the requirement is decreasing in `S`,
any `S` at or above that maximum satisfies it at EVERY `tau`. The maximum is the
smallest such `S`, so it is the minimax answer rather than a guess at a length.

`ecological_spinup_multiple` returns it, and it is exactly proportional to the
retained record with a constant fixed by the tolerance alone: 6.86325 records at
a tolerance of 0.05, attained at a relaxation time of 7.3519 records; 17.896 at
0.02 and 3.1901 at 0.10. The proportionality is checked against
`ecological_spinup_cycles` itself by sweeping `tau` at three record lengths.

A MEASURED RELAXATION TIME COULD ONLY ASK FOR LESS, so the minimax is a ceiling
on the floor and the module prefers a measured value where one is admissible.
The gate checks the claim rather than reporting it: the returned spin-up must
hold the residual drift inside the tolerance across eight relaxation times
spanning five decades around the worst case, and one per cent shorter must fail
at the worst case.

WHAT THE SPIN-UP BUYS IS NOT CORRECTNESS. A residual transient enters the very
statistic the acceptance contract gates on, so an under-spun run is REFUSED
rather than silently accepted. What `S*` buys is that the run is not refused for
that reason, and at 6.86 records it is most of a run's cost -- which is what
makes it worth deriving rather than typing.

### What it costs

Measured throughput on the runs on record, from their own manifests, is 0.9267
seconds per simulated year at npatch 5 on 16 ranks over 1617 cells: 990.7 s for
1098 simulated years and 1824.7 s for 1998, which agree to one per cent and
imply a negligible fixed cost. Multiply the total cycle count by that; npatch 20
is about four times it.

On the record on disk that prices the next run at about 1253 retained cycles and
a spin-up of 6.86325 times it, roughly 9850 cycles in total and about two and a
half hours at npatch 5 on 16 ranks. The run this contract used to demand was
33846 retained cycles and a spin-up of 232294, about seventy hours.

The time base is correct and is not implicated in any of this. The Earth
convention that `nyear_spinup` used to carry is 500 simulation years, which at
0.5010172 Earth years per orbit is 998 of this world's, and `freenyears` 200 is
100.2 Earth years; both matched LPJ-GUESS's convention exactly. The convention is
Earth's, and this world's woody types and slow pools are not Earth's.
`build_vesper_pfts.py` reads the derived spin-up floor from whatever acceptance
artifacts the live contract has produced and records which derivation gave it;
the convention is reached only where there is no such artifact at all.

AN ARTIFACT FROM A SUPERSEDED CONTRACT IS NOT POOLED. What a resolving length
means is "the record at which THIS quantity closes on THAT tolerance", and both
halves belong to the contract it was taken under, so pooling contract 4's
artifacts would size every future run from the contract that no longer judges
it. `ecological_timescale_brackets` names the superseded ones and refuses when
none is current.

## What the contract cannot yet do, and what has never been staged

`build_surface_albedo.py --mode modelled` has never been staged, since both
`albedo_report.json` files record mode `vegetated`, and it cannot be while
`read_foliar_cover` calls `require_lpj_acceptance` and no LPJ run in the tree is
accepted. Every cover magnitude quoted from a run in this note is therefore from
a refused run and is indicative of scale only.

WHAT SEPARATES THE TREE FROM AN ACCEPTED RUN IS NOW ONE QUANTITY. Re-assessed
under contract 5, `lpj_7d3c576ee4e342acb097b46fece976e0` settles ten of its
eleven assessed quantities on the record it already has -- both the cover
quantities the albedo consumer reads among them -- and is refused by
`cpool.out` total carbon at 0.0503 against 0.05. That refusal is a drift and not
a resolution: a settled series of its scatter closes inside the record it has,
and the drift lives in vegetation carbon, which is what an under-spun run leaves
behind. The run carried a 2318-cycle spin-up against the 8600 the minimax now
derives.

A run at the derived floors -- 8600 spin-up, 1253 retained, about two and a half
hours -- is what would settle it, and it is the only thing between this tree and
its first accepted LPJ run. It is not a purchase this note makes.
