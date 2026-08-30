# The equilibrium window's trending-cell limit, and the null it was missing

Vesper is a simulated world. This note is about the acceptance rule that decides
whether a LPJ-GUESS run of its biosphere has settled, and about a number in that
rule that had no derivation. Measured on 2026-08-30 against
`biosphere/config/equilibrium_window.yaml` as it then stood, contract version
`vesper-lpj-equilibrium-window/1`.

## What the rule did

A field's cycle-mean series over the last 10 complete forcing cycles was called
trending in a gridcell when its end-to-end change exceeded 5 per cent of the cell's
window mean AND its slope exceeded 2 standard errors. The run was refused when the
spatial mean met the same two conditions, or when more than
`maximum_trending_cell_fraction: 0.25` of ALL cells were individually flagged.

The global half works. Its behaviour across run lengths is recorded on world-qcse:
TeNE and IBS cleared between 100 and 1000 orbits because their drift was real and
coherent, and the spatial-mean end-to-end change is the quantity that responded.

The 0.25 had no recorded derivation, and neither of the two knobs that should have
relaxed it moved it: across 100 to 1000 orbits the flagged fraction went 0.376 to
0.365, and across npatch 5 to 20 it went 0.376 to 0.348 while the median per-cell
drift fell to 0.603 of its value and the median relative slope standard error fell
to 0.658 of its value. The test is the ratio of those two, so reducing the noise
cannot satisfy it.

## The null costs no model run

The forcing cycle is one simulation year: `run_manifest.forcing` records
`cycle_years: 1` against a `VESPDRV8` driver of 12 intervals, and every simulated
year replays byte-identical forcing. All interannual variation in a retained record
is therefore internal stochasticity of the simulated vegetation - patch disturbance,
establishment, fire - and not forcing. A long run at fixed forcing is already the
two-seed experiment the null was thought to need.

Method. Take the 1000-orbit diagnostic `lpj_377476c097b549fcaaa07cc7644275ad`
(npatch 5, 1617 cells, 1000 cycles). Subtract a per-cell linear fit over the whole
record and add the per-cell mean back, which removes any genuine approach to
equilibrium and leaves stationary internal variability with the patch ensemble's
persistence intact. Cut 100 disjoint 10-cycle windows and apply the contract's own
per-cell test to each. Because the gate refuses a run if ANY field exceeds its
limit, the null distribution of a refusal is the per-window MAXIMUM over fields.

Reproduce with `biosphere/scripts/derive_trend_null.py <run>`.

## What the null is

Per-window maximum over the assessed fields, 100 windows, denominator as the rule
then had it (all cells):

    median 0.392    p90 0.450    p95 0.455    p99 0.461    maximum 0.462

The same statistic on `lpj_36b0f5fc91a640ed94618b0e95420916` (npatch 20, 10
windows) gives median 0.367 and maximum 0.383, so the null is near enough
invariant to patch count, as the ratio argument predicts. Flag rate against window
position is flat (first-half median 0.394, second-half 0.388), so the detrending
left no residual drift for the estimate to absorb.

So the limit sat BELOW the rate at which the test flags cells with provably no
trend to find. The observed fractions that refused the two npatch-5 runs were 0.376
and 0.365 for lai C3G against a null median of 0.379 and a null 95th percentile of
0.410, and 0.309 and 0.290 for lai Total against a null median of 0.298 and a 95th
percentile of 0.314. All four sit inside the null's 95th percentile and three of the
four sit below its median: the runs were refused for showing no more apparent
per-cell trending than a stationary surrogate of themselves.

## Why no single number could have worked

The null is a per-field quantity spanning three orders of magnitude, because it is
set by each field's own internal variability and by how much of the grid it
occupies, not by anything about equilibrium. Null median and power at twice the
contract's own drift limit, all-cells denominator:

| field | occupancy | null | power at 2x limit |
| --- | --- | --- | --- |
| npool.out SoilN | 1.00 | 0.001 | 0.998 |
| npool.out Total | 1.00 | 0.001 | 1.000 |
| cpool.out SoilC | 0.94 | 0.001 | 0.942 |
| cpool.out Total | 0.94 | 0.058 | 0.907 |
| tot_runoff.out Total | 1.00 | 0.133 | 0.846 |
| npool.out LitterN | 0.91 | 0.260 | 0.701 |
| lai.out Total | 0.93 | 0.298 | 0.552 |
| lai.out C3G | 0.93 | 0.379 | 0.539 |
| cpool.out VegC | 0.90 | 0.387 | 0.605 |
| aaet.out C3G | 0.92 | 0.379 | 0.576 |
| lai.out TeNE | 0.10 | 0.066 | 0.085 |
| lai.out C4G | 0.15 | 0.075 | 0.099 |
| lai.out IBS | 0.20 | 0.118 | 0.135 |
| lai.out TeBE | 0.24 | 0.119 | 0.184 |

Two separate defects are visible in that table.

The DENOMINATOR was every cell rather than the cells the field occupies, so a
plant functional type present on a sixth of the grid could not reach a
quarter-of-all-cells limit however hard it drifted. Nine of the thirteen types
occupy under 0.37 of cells and five occupy under 0.25; for those the measured
power at twice the contract's drift limit is 0.08 to 0.19 against a 0.25 limit.
The guard was inert for them.

The SPREAD is irreducible by any single choice. Setting the limit high enough for
the patch-driven grass and vegetation fields destroys it for the slow soil pools,
where it discriminates almost perfectly: SoilN's null is 0.001 and its power at
twice the limit is 0.998. Those pools are where per-cell disequilibrium is the
real risk, because they are what the CENTURY accelerator hands off.

## Where the excess over the textbook rate comes from

For 10 independent points a 2-standard-error slope test has 8 degrees of freedom
and a two-sided probability of 0.081. The measured null runs to 0.39. The gap is
temporal memory in the cycle means. A within-window circular-shift surrogate,
which preserves each cell's marginal distribution but whitens its serial
structure, returns 0.007 to 0.09 across fields, close to the independent-samples
rate, while the memory-preserving detrended long-run null returns 0.001 to 0.40.

That also settles a tempting shortcut: the circular-shift surrogate is trend-blind
as required, returning the same value under an injected drift as under none, but
it is not a calibrated null and using it as one would refuse everything.

The ordinary least-squares slope standard error in the reducer models no memory.
`lib/autocorrelation.py` is the module for a series whose samples carry memory,
and the reducer does not use it; repairing the standard error rather than
calibrating around it is a separate disposition and is tracked.

## What replaced it, and what replaced that

`vesper-lpj-equilibrium-window/2`. The declared fraction is gone. Each field's
limit is measured from the run being judged: the record before the acceptance
window is detrended cell by cell, cut into windows of the contract's own length,
and the field's limit is the largest share of its occupied cells that any of those
windows flags. The acceptance window itself is judged undetrended, so a drift
present throughout the record is removed from the reference and left in the
quantity being judged.

A stationary field exceeds the largest of N such windows with probability 1/(N+1),
so the retained record length, and not a chosen number, sets the rate at which the
contract wrongly refuses. The contract declares that rate, 0.05, which requires 19
null windows beside the acceptance window, so 20 windows and 200 retained cycles
in all.

Measured on the 1000-orbit run, 99 null windows, false-refusal rate 0.010, over
the 64 assessed fields of the seven stability tables:

| acceptance window | cell half refuses | global half refuses |
| --- | --- | --- |
| as run | 0 of 64 | 1 of 64 |
| coherent drift at 2x the contract's limit | 46 of 64 | 63 of 64 |
| cancelling regional drift at 2x the limit | 42 of 64 | 0 of 64 |

The last row is what the per-cell half exists for, and is the one case the global
half cannot see: opposed regional drifts that leave the spatial mean flat. The
earlier form did not separate that case from an undisturbed run on the noisy
fields, because its flagged fraction under a cancelling dipole and under no drift
differed by less than its own window-to-window scatter.

## What this leaves refused

The 1000-orbit run now passes the per-cell half on all seven stability tables and
is refused on one field only, `lai.out` TrIBE, whose spatial mean moves 0.055 over
the window at 5.65 standard errors. That is the half that works, on a plant
functional type that had not established at 100 orbits and is still expanding at
1000. The refusal is real.

The two 100-orbit runs are refused for a different and honest reason: a 100-cycle
record leaves 9 null windows, so a stationary field would face a 0.10 false-refusal
rate against a declared 0.05. Their acceptance needs a longer retained record, not
a different limit.


## The window is shorter than the thing it tests

Measured 2026-08-30 on the same 1000-orbit record, through
`lib/autocorrelation.py`, which owns the estimator and the two verdicts that say
whether its answer means anything. Reproduce with
`biosphere/scripts/derive_trend_null.py <run> --timescales`.

The integrated autocorrelation time of the spatial-mean cycle series runs 9.3 to
125.3 complete forcing cycles. The window the contract tests over is 10. Sixty-one
of the 64 assessed fields have a spatial-mean memory time longer than the whole
window, per-cell medians run 2.9 to 93.1, and the share of cells whose own memory
time exceeds the window is 0.22 to 1.00.

| table.field | spatial tau | cell tau median | cells with tau > window |
| --- | --- | --- | --- |
| lai.out BNE | 75.3 | 64.9 | 1.00 |
| lai.out TrIBE | 93.5 | 63.9 | 0.97 |
| lai.out C3G | 9.9 | 6.3 | 0.36 |
| cpool.out Total | 125.3 | 68.1 | 0.89 |
| npool.out SoilN | 60.6 | 89.0 | 0.87 |
| tot_runoff.out Total | 14.9 | 5.5 | 0.33 |

A trend fitted inside one memory time is not a trend. It is one smooth excursion of
a process that has not had time to sample its own distribution; the residuals
within the span are small because the process is smooth on that scale, so the
ordinary least-squares standard error is small and the slope looks decisive. That
is the whole of the measured 0.47 to 0.86 per-cell flag rate against a nominal
0.081, and it is failure-modes class 34: a bed shorter than its startup.

NO CORRECTION RESCUES IT. Inflating the slope standard error by the square root of
the memory time, which is the correction `lib/autocorrelation.py` owns, moves
lai.out's bare flag rate from 0.47-0.84 to 0.19-0.46: two to three times better,
still two to six times nominal, still spanning a factor of four across fields. The
correction is for a MEAN over a span and cannot make a span shorter than its memory
carry a trend.

## The record is shorter than the approach it watches

An exponential approach fitted to each field's spatial-mean series over the full
1000 cycles is INADMISSIBLE for most fields, by a criterion declared before it was
run: the asymptote must lie inside the record's own range and the e-folding time
must be shorter than the record. Thirty-four of 64 fits put the asymptote outside
the record, meaning the series has not turned over within 1000 cycles; three more
return a timescale longer than the record. Where a fit is admissible the e-folding
times are 22.9 to 775.9 cycles, mostly in the hundreds: lai.out TrIBE 685.1,
aaet.out Total 775.9, lai.out TeBS 492.4.

`lib/autocorrelation.py:stationary_enough` independently refuses the spatial-mean
series of most fields: the linear trend across 1000 cycles moves them by more than
their own standard deviation.

So the two bounds on run length are:

- On the RETAINED record, from the memory time and the module's own span bar: at
  least ten times tau, so 100 to 1250 cycles by field. The 200-cycle statistical
  bound is not the binding one for the woody types or the slow pools.
- On the SPIN-UP, from the relaxation time: several times an e-folding time of
  hundreds of cycles, so of order 2000 to 5000. `vesper_pfts.ins` sets
  nyear_spinup 998. The spin-up is short by a factor of two to five.

That is why a run at literally fixed forcing is still drifting coherently after
1000 further cycles. The time base is correct and is not implicated: 998 spin-up
cycles is 500.0 Earth years and matches the community convention exactly. The
convention is simply too short for this world's woody types and slow pools.

## Contract 2's declared rate was per field, and it refused per run

`maximum_false_refusal_rate: 0.05` was exact PER FIELD: a stationary field exceeds
the largest of N detrended null windows with probability 1/(N+1). But the gate
refuses when ANY of the 64 assessed fields exceeds, and the family rate is
1-(1-1/(N+1))^F_eff. Measured by leave-one-out on the 1000-orbit run's 100 windows,
it is 0.34, flat across `slope_standard_errors` from 1.0 to 4.0 as the construction
requires. The helper was validated first on synthetic fields: at 100 windows it
returns 0.010 for one field, 0.080 for eight and 0.460 for sixty-four independent
ones, against an analytic 0.010, 0.077 and 0.474. The measured 0.34 corresponds to
about 35 effective independent fields out of 64.

Raising the record does not fix it. The finest probability an N-window empirical
null can express is 1/(N+1), so a family rate of 0.05 needs about 710 windows, 7100
retained cycles, roughly twenty hours of model time, re-measured per configuration.
The repair is an analytic null, which needs a valid standard error, which needs a
window longer than the memory time.

## What the contract is now

`vesper-lpj-equilibrium-window/3` adds one thing and renames one thing. Both are
fully measured; nothing else changed, deliberately.

- A MEMORY-ADEQUACY GUARD, ahead of everything else. It applies
  `lib/autocorrelation.py`'s declared span bar to both spans that must carry the
  memory time: the retained record, where tau is estimated, and the acceptance
  window, where the slope is actually fitted. It fails closed on every run this
  project has: the 100-cycle runs cannot establish their fields' memory times at
  all, and the 1000-cycle run has memory times of 9 to 125 cycles against a
  10-cycle window. Every refusal now names the field, its memory time, and the
  record length that would answer it.
- `maximum_false_refusal_rate` is renamed `per_field_false_refusal_rate`, which is
  what the construction controls, with the measured family rate recorded beside it.

The trend test below the guard is unchanged and is unreachable, which is the honest
state: it is not valid at any window this model supports, and its replacement is a
piece of statistical design rather than a threshold. What that replacement has to
handle is recorded on world-ioxr, including three subtleties found and discarded in
one afternoon: a half-to-half mean difference is half the end-to-end change and
silently doubles the declared tolerance unless scaled; the span bar belongs on the
half whose mean is taken, not on the record; and guarding on an estimated tau and
then using that same estimate for the standard error selects for underestimates and
is anti-conservative near the bar. A candidate that survives those measured a
family refusal rate of 0.040 against a declared 0.05, with power 0.31 at the
contract's own drift limit and 0.92 at twice it -- promising, and not validated
enough to ship on the day it was written.

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

A 5 per cent drift in both tree and grass cover is worth about +0.064 K of surface
temperature, bracketed +0.032 to +0.128 K. Against what the climate arm resolves --
0.15 K for the single-run convergence criterion, 0.05 K for its resolving-power bar,
0.027 K for the standard error of a 57-orbit climatology mean -- the derived limit
spans 0.02 to 0.23 and cannot choose the number.

The bracket is dominated by one unmeasured link. There is no measured
d(planetary albedo)/d(surface albedo) in this tree: `lib/sensitivity.py` explicitly
declines to own it, `scripts/error_budget.py` declares 0.5 and says in the same
breath that the honest claim is a factor of two, and the single measured point is
0.3045 from one paired arm. world-ckbt is open to measure it. Until it closes, 0.05
stays a declared tolerance with this bracket recorded rather than a derivation being
invented for it.

Two things the chain also exposed: `--mode modelled` has never been staged, since
both `albedo_report.json` files record mode `vegetated`; and it cannot be, because
`read_foliar_cover` calls `require_lpj_acceptance` and every LPJ run in the tree is
refused. The cover magnitudes above are therefore from a refused run and are
indicative of scale only.

## How long a run this world's ecology needs

Both bounds are now derived rather than estimated, and both are carried by
`lib/run_lengths.py` beside the climate model's pair. They are a second pair, in
COMPLETE FORCING CYCLES rather than orbits, and the module keeps them apart: a
cycle is one simulation year is one modelled orbit, so the two units coincide
numerically only while the driver's cycle is one year, and the run's own declared
cycle length travels with every reading.

### The relaxation time, measured without its asymptote

The exponential fits above refused for 34 of 64 fields because the asymptote landed
outside the record, and a refused fit gives no number. A second estimator does not
need the asymptote: for `y = A - B exp(-t/tau)` the difference between consecutive
equal blocks decays by `exp(-Q/tau)` and `A` cancels out of the ratio, so four equal
blocks of the record give `tau = -Q / ln(r)` with `r` the ratio of successive
differences of the block means. `lib/lpj_output.py:relaxation_time` is the
implementation and its admissibility conditions were declared before it was run: the
successive differences must share a sign, their ratio must be a contraction in
(0, 1), each difference must exceed twice its own MEMORY-CORRECTED standard error,
and the second ratio must agree with the first inside a factor of two.

On the 1000-cycle record, five fields across the two assessed runs clear all four.

| field | relaxation, cycles | memory, cycles |
| --- | --- | --- |
| lai.out BNS | 464.8 | 48.7 |
| anpp.out BNS | 767.0 | 39.4 |
| npool.out Total | 846.8 | 78.3 |

The other 123 are declined, most for a block difference inside twice its own
standard error or for successive differences that change sign, which is what a
memory time of 60 to 125 cycles does to blocks of 250. So the relaxation bracket is
59.7 to 846.8 cycles with an OPEN TOP: a field the estimator declines may be slower
still, and that is why what comes out of it is a floor rather than a length.

### The two floors

The RETAINED RECORD floor is the memory time's, at the top of that bracket because a
record must serve every field it will be asked to judge, and at the span multiple
`lib/autocorrelation.py` declares for any span whose mean is taken:

    record >= 10 x 125.3 = 1253 cycles

The SPIN-UP floor is the relaxation time's, and it takes its residual from the
criterion that judges what follows, exactly as `SETTLING_RESIDUAL_K` takes 0.15 K
from the offset criterion. A spin-up from bare ground starts a full equilibrium level
away, so the approach remaining after `S` cycles is `exp(-S/tau)` of the level and the
drift it puts across a retained record of `L` is that times `1 - exp(-L/tau)`. The
acceptance contract refuses a record whose relative end-to-end change exceeds
`relative_end_to_end_limit`, so

    exp(-S/tau) x (1 - exp(-L/tau)) <= 0.05

At tau 846.8 and L 1253 that is S >= 2318 cycles. Starting from bare ground is the
conservative reading and is stated rather than hidden: a pool the CENTURY accelerator
hands over part-grown begins closer than a full level away and needs less.

### What it costs, which is the point

    spin-up   2318 cycles
    record    1253 cycles
    total     3571 cycles

Measured throughput on the runs on record, from their own manifests, is 0.9267 seconds
per simulated year at npatch 5 on 16 ranks over 1617 cells: 990.7 s for 1098 simulated
years and 1824.7 s for 1998, which agree to one per cent and imply a negligible fixed
cost. So the whole run is about 55 minutes, and about 3.7 hours at npatch 20.

THE ANSWER IS THAT IT IS CHEAP. The spin-up this world's ecology needs is 2.3 times
Earth's convention and costs under an hour, and the reason no run has ever had it is
not expense but that nobody had derived it. `vesper_pfts.ins` took the shipped 500 and
rescaled it correctly into 998 simulation years; the conversion was right and the
convention was Earth's.

Both numbers are floors and the derivation says why: the relaxation bracket's top is
open, and the memory time read off this model grows with the window it is read on, so
the next record re-reads its own floor. `build_vesper_pfts.py` reads the spin-up floor
when an assessed run exists and falls back to the rescaled convention only when none
does, recording which it used.
