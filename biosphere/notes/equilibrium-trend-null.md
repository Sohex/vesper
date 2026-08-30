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

## What replaced it

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
