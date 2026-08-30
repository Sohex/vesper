# The interval the carve verdict's open-water evaporation is evaluated over

Vesper's closed basins. Measured on 2026-08-27 against `canonical-10m-base`
(8,772 basins, `selectionSource: threshold`, so pre-carve and never carved)
forced by `bootstrap_regular_climatology.nc`, a T21 bootstrap.

Found while fixing the same defect one file over; the lake balance's version is
in `hydrography/notes/lake-balance-integration.md`.

## What was wrong

`carve_verdict.py:main()` read every climatology field through `annual_mean`
and called `reference_level_air` with no bin index, then evaluated Penman once
on that annual-mean air. `export_carve_list.py:climate_terms` did the same.
Penman is nonlinear in everything it reads, which is already the stated reason
it integrates the DIURNAL cycle rather than reading a daily mean; the seasonal
cycle is the same argument at a longer period, and
`config/land_water_ledger.yaml` holds `open_water_evaporation` at
`interval_floor: climatology_bin`. So the annual evaluation was against a
standing decision rather than a simplification anyone had chosen.

## The prediction, stated before measuring

A basin carves when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

and the right-hand side is pure geometry that does not move. Raising `E` raises
the left-hand side, so the inequality is satisfied for strictly fewer basins.
**The corrected carve set had to be a SUBSET of the old one**: basins leaving,
none entering. A basin entering would have meant the fix went in backwards, and
a mixed set would have meant something other than the evaporation had moved.

Measured: 4,772 to 3,944. **828 left the carve set and 0 entered**, subset
confirmed on the basin ids rather than on the counts. That is 9.4% of the
catalogue and 17.4% of the old carve set.

## Why the per-bin evaluation is not simply the right answer

The bin mean of the per-bin evaluation is 1.569x the single annual evaluation
over land, 1.198x over ocean. That is large enough to check before believing,
and the check does not come out one-sided.

**The magnitude is not any single input.** Taking one input per bin at a time
and leaving the rest annual, over land: temperature 1.106, humidity 1.099,
radiation 1.008, wind 0.998, diurnal range 0.999, and all of temperature,
humidity and pressure together 1.127. None of them, and no product of them,
reaches 1.569. The effect is an interaction, and the interacting parts are
Penman's two hard clamps.

**Per bin is right because of rectification.** `max(net_radiation, 0)` and
`max(es_a - e_air, 0)` are rectifiers, because a lake does not evaporate a
negative amount in a dark month. Averaging the air over a year first cancels a
real summer against a winter that never physically happens. Over this
climatology the bin mean of `max(net_radiation, 0)` is 86.4 W/m2 over land
against 77.0 at the annual mean, and 684 of 1,019 land cells go negative in
some bin while staying positive in the annual mean.

**Annual is right because Penman carries no heat storage.** It closes the
surface energy budget instant by instant, so per bin it charges the water for
bright-season energy that actually went into the water column and comes back
out in the dark season. That term integrates to zero over a periodic year and
to nothing like zero over a bin.

**The ocean measures the second effect**, because it is the one surface whose
evaporation the model already computes with the true storage in it:

| | ocean mean, mm/day | against the model |
| --- | --- | --- |
| the model's own evaporation | 1.770 | -- |
| Penman on annual-mean air | 1.685 | 0.952 |
| Penman as the bin mean | 2.019 | 1.141 |

and bin by bin the per-bin estimate runs 1.43 to 1.46 times the model in the
three brightest bins and 0.89 to 0.96 in the two dimmest, ranging 0.887 to
1.463 over the twelve. That is the signature of a missing storage term and not
a constant bias: a constant bias would have shown as the same ratio in every
bin, and would have meant the annual agreement was luck.

## What is done about it

**The interval is bracketed rather than chosen**, because the two ends are
interpretable and the spread between them is a quantity this project does not
have:

- the **annual evaluation** is the limit for a lake deep enough to hold its
  temperature through the year;
- the **bin mean** is the limit for a lake with no heat capacity at all;
- what sits between them is the water body's own heat storage, which needs a
  lake depth and a mixed-layer model. Declared, not estimated.

An ocean stores far more heat than any lake, so the 14.1% is an UPPER bound on
the bias at the bin-mean end, and a shallow playa carries almost none of it.

Both ends are bounds alongside `wet`, and the carve list is the intersection --
which is what the file already did for the two evaporation estimators and what
`export_carve_list.py` already did for the two climate arms, by taking the
larger retain. Ties break toward not carving: a basin carved in error has lost
a depression from the terrain everything else is built on, and a basin left
uncarved is still there to carve next cycle.

| | carves |
| --- | --- |
| `carve_verdict.py`, bin mean | 3,944 |
| `carve_verdict.py`, annual evaluation | 4,772 |
| `carve_verdict.py`, wet bound | 6,260 |
| **the carve list, intersection of all three** | **3,944** |
| `export_carve_list.py`, bin mean | 3,967 |
| `export_carve_list.py`, annual evaluation | 5,337 |
| **the exported list, intersection** | **3,967** |

The two scripts differ by 23 basins because they are two formulations of one
criterion -- the ratio form and the discharge form -- which the export already
documents as not nesting. The interval bracket is 25.7% of the union in the
export's discharge form.

## Two other defects found in the same files

**The cold arm was the warm arm.** `climate_terms(clim_path, args, config,
basins)` read `args.climatology` in both places it opens a climate and never
touched `clim_path`, so `climate_terms(args.endmember_climatology, ...)`
evaluated the primary climatology. The two-climate bracket was one climate
compared with itself, and it did not fail: `cut_endmember` equalled
`cut_primary`, so the reported disagreement was 0 and the bracket width 0.0% of
the union, which reads as the warm and cold ends agreeing completely. `main`
now refuses two arms whose discharge fields come out bit-identical. It cannot be
exercised on real data yet: the tree holds one regular climatology.

**The rung guard was unreachable in practice.** `require_configured_grid` was
reached only through `climatology_path()`, and both scripts call that ONLY when
`--climatology` is absent -- so every re-take, every arm and every sensitivity
skipped it. The rung appears nowhere in a climatology's name and the carve list
header takes its resolution from `config/planet.yaml`, so a sidecar could have
recorded T21 beside a T85 climatology with nothing in the tree contradicting
itself. The guard now runs on whatever climatology is in force, the endmember
arm included.

## What this costs

Nothing on the active lineage. `canonical-10m-base` carries
`selectionSource: threshold`: it is pre-carve and no carve list has ever been
applied to it, so the instrument is corrected before its first use here.

The four builds that were carved from a verdict -- `carved-zoned`,
`carved-zoned-v2`, `carved-zoned-v4` and `carved-zoned-v5`, all
`selectionSource: preserve-list` -- were taken with the annual evaluation and
with the endmember arm silently equal to the primary. Each therefore over-carved
by the direction argued above, and each preserved fewer basins than the
threshold rule does. All four are archived, payload deleted, and all predate
this lineage. Under CLAUDE.md rule 7 nothing is owed to them: a build is
disposable until the canonical climatology lineage is declared, and it is not.

# The concavity in the discharge term, and what the missing delivery phase costs

Measured on 2026-08-30 against `canonical-10m-carve1` (4,657 basins) forced by
`baseline_regular_climatology.nc`, the verdict that produced
`canonical-10m-carve2`. The evaporation interval above is one nonlinear
evaluation on an annual mean; this is the other one in the same expression, one
term over, and it is the member of the class that cannot be repaired the same
way.

## The defect, verified against the files

`export_carve_list.py:incision_retain` evaluates

    cut = coefficient * erodibility * slope**n * Q**m

with `m = 0.5` and `Q` the ANNUAL overflow at spill level. The power is
concave, so `sqrt(mean Q)` is at or above `mean(sqrt Q)` and a basin that
delivers its overflow in a season is credited with more cutting than that
season's discharge does. `Q` is `runoff * (catchment - area_at_spill) -
(E - P) * area_at_spill`, built from the catchment's annual runoff.

`hydrography/config/land_water_ledger.yaml` declares
`seasonal_phase_of_catchment_delivery` under `absences`, and the declaration is
current: the hydrography path computes catchment runoff as the annual mean of
`P - E` and cannot phase it through the year, because over one cycle at steady
state the annual mean is exactly what the cell generated while per-bin clamping
would count the wet season's supply twice. So there is no per-bin discharge to
evaluate, and the repair that fixed the evaporation interval is unavailable
here.

## Two results that do not need the phase

**The level is absorbed by the calibration, exactly.** Write the truth as
`mean(sqrt Q) = phi * sqrt(Q_annual)`, so `phi` is the whole of the concavity
error on one basin. `cut` is linear in the coefficient and linear in `Q**m`,
and `calibrate_coefficient` SOLVES the coefficient on every run from this same
discharge field against Earth's standing-basin density. A `phi` common to every
basin is therefore met by a coefficient of exactly `coefficient / phi`.
Measured: scaling every discharge to give `phi = 0.9`, `0.5` and `0.2871` and
re-solving returns `299.3884`, `538.8992` and `938.5720` against the shipped
`269.4496`, which is `coefficient / phi` to 2.4e-8, the bisection's own
resolution; with the coefficient set to `coefficient / phi` outright the retain
vector matches to 3.3e-16 and no basin changes class. **A bracket over the
common level would have zero width.** That is what settles whether the carve
list can carry a bracketed discharge: it is not two carve lists and not two
generations, because the second arm is the first one.

**The spread is bounded by the bin weights alone.** Two facts and nothing else:
an overflow is never negative, and the cycle has the climatology's own twelve
bins at `lib/climatology.py`'s weights, the lightest of which is 0.082418. An
evenly delivered discharge gives `phi = 1`; the most uneven admissible one puts
the whole year's overflow in the lightest bin and gives `phi = w_min**(1 - m) =
0.2871`. So `phi` lies in `[0.2871, 1]` whatever the phase is, the ratio
between two basins' `phi` lies within `J = 3.4833`, and with the population
level pinned by the calibration a basin's cut can move at most a factor `J`
either side of the class boundary.

Nothing in either result assumes when the water arrives.

## Against the instrument

`J = 3.4833` against `P = 2.2235`, the ratio of the Poisson bracket on Earth's
count of 15 standing basins, which is the criterion `sweep_size_floor` already
fixed: `J > P` is DECISIVE and has to be declared and carried, `J <= P` is
SUBORDINATE. **DECISIVE**, so it is carried per basin rather than dismissed. It
is not the largest lever on this coefficient: the size floor on the Earth
sample moves it by a factor of 56.3 over the span the build derives, and that
is reported beside the verdict already.

## Which basins, which is the question that matters

The verdict is categorical, so the number worth having is not how far `Q**0.5`
moves but whether a basin changes side. On the finished-depression basis, the
basis the coefficient is solved on, with `x` for the annual `cut / depth`:

| | basins |
| --- | --- |
| overflow at all | 416 |
| of those, `x >= 1`, cut | 356 |
| of those, `x < 1`, overflow but keep some rim | 60 |
| **movable: `x` inside `(1/J, J)`** | **99** |
| movable and currently cut | 67 |
| movable and currently standing | 32 |
| held cut under every admissible phase, `x >= J` | 289 |
| held standing under every admissible phase, `x <= 1/J` | 28 |

`x` over the overflowing set runs 0.23 at the 5th percentile, 6.65 at the
median and 70.8 at the 95th, so most of the population is nowhere near the
boundary and the exposure is a minority of a minority.

Against the verdict as exported, taking the union of the movable set over the
three arms the list is a maximum over:

| verdict | basins | movable |
| --- | --- | --- |
| carve | 195 | 26 |
| marginal | 41 | 24 |
| preserve | 4,421 | 49 |

**None of the 195 carved basins is in the climate bracket** -- the intersection
rule removes bracketed basins from the carve set by construction, and 0 of the
290 bracketed basins carve -- so these 26 are an exposure the two-climate
bracket does not already cover. That is why the flag is per basin: 26 of 195 is
a count, and a count cannot say which depression is at risk of being removed
from the terrain everything else is built on.

## What is done about it

`seasonal_concavity` computes the bound in the run that takes the verdict, puts
it in the sidecar's `method` block, and stamps `seasonal_concavity_movable` on
every basin record. `false` is the strong statement and is what most basins
carry: no admissible seasonal concentration of the overflow changes that
basin's class. The verdict itself is unchanged, and deliberately: the direction
of the raw concavity is one-sided, but the calibration makes what reaches the
verdict two-sided, so 32 standing basins are exposed to being cut for every 67
cut basins exposed to standing. Taking one side of that would be a preference
dressed as a bound.

`export_carve_list.py --selftest` carries the identities, and each was
demonstrated to fail on an injected defect: the Jensen direction with its
equality case, the bound being attained by the whole year in the lightest bin
and never breached, the calibration absorbing a common discount exactly while
an uncompensated one does not, the band being read off the same cut over depth
the retain is, the floor the code reports being the one the lake balance
attains, and the claim the whole disposition rests on -- that under an
admissible per-basin `phi` with the coefficient re-solved, no basin outside the
band changes class. The band and the retain read ONE arithmetic:
`cut_over_depth` is factored out of `incision_retain` precisely because the
retain is clipped and the band needs the ratio above 1.

Two more guard the disposition rather than the arithmetic. The attainment sweep
below runs `lake_balance.solve_periodic` and fails if the lake's own season
lifts the concentrated arm off the floor, with the no-storage net beside it as
the control. And `scripts/smoke_test.py` holds
`export_carve_list.LEDGER_ABSENCES_RELIED_ON` to a key under `absences:` in the
ledger: the choice to BOUND this error rather than correct it rests on the
delivery phase being unrepresentable, and closing that absence has to fail the
gate at the site counting on it rather than leave the declaration standing and
untrue.

## The bound is attained, and the lake's own season does not narrow it

The overflow is catchment delivery, whose phase is the absence, plus the lake's
own surface flux, whose phase IS representable and is carried per bin in
`surface_water.nc`. Bounding the two separately looks strictly narrower and it
is not. **Lake storage stands between the water balance and the overflow, and
it sharpens the season rather than damping it**: a lake drawn below spill has to
refill before it spills again, so a deficit season that would absorb delivery in
a no-storage model instead delays the spill and concentrates it further.

Measured through `lake_balance.solve_periodic`, the pipeline's own lake balance,
on its own synthetic set of six basins spanning a factor of a hundred in
capacity, at 12 bins whose lightest weight is 0.043478 so the floor is
`w_min**(1 - m) = 0.208514`. Two delivery phases carrying the same annual water,
against four lake seasonal swings:

| lake evaporation swing | flat delivery | whole year in the lightest bin |
| --- | --- | --- |
| 0.0 | 1.000000 | 0.208514 |
| 0.3 | 0.991854 | 0.208514 |
| 0.6 | 0.963500 | 0.208514 |
| 0.9 | 0.865745 | 0.208514 |

The concentrated arm sits ON the floor to six decimals for every basin, and the
lake's own swing does not move it by so much as a part in a million. A bound
narrowed by the lake's per-bin flux is exactly a bound that this sweep would
lift off the floor.

**The same run closes the bound's own premise**, which was stated as an
assumption and is an identity here. Under the delivery this pipeline represents,
an overflowing basin sits AT spill through the entire cycle -- cycle-mean area
over area at spill measures 1.0000 -- so its periodic annual overflow equals the
spill-level `Q` the carve verdict computes, to 1.0000. And under the adversary's
concentrated delivery the drawn-down lake evaporates over less area and passes
MORE water, which puts the worst admissible cycle at 0.03 to 2.6% ABOVE the
floor rather than below it. The bound holds and is slightly conservative.

The no-storage net is the control, and it is what made the tightening look
available: reading the same forcing as the instantaneous `delivery + (P - E)`
floored at zero gives a ratio strictly above the floor, because it lets the
lake's deficit season absorb delivery that a real lake would have stored and
spilled later. `_selftest` carries the inequality.

So `J = 3.4833` stands, the call stays DECISIVE, and the 26 carved basins are
still movable. There is no narrower bound to be had from what this pipeline
carries.

Closing the absence outright is a different decision and the ledger already
names it: it is the same decision as adopting the climate column's
`surface_runoff`, which this path rejects as an incomplete routed diagnostic.
What that would buy is not a narrower bound but a delivery phase, and with one
the correction becomes computable rather than bounded.
