# The closure checks difference two sides that do not span the same inventory

Found 2026-08-31, on run `lpj_de8a0c8bd4324639b4bd62b66867d23f`, while
disposing of the residual `notes/audits/closure-tolerance-under-written-precision.md`
left behind once the write-out quantisation was removed.

Worldbuilding. Vesper is an invented planet; every nitrogen and water quantity
below is a modelled soil, vegetation or hydrological quantity in a simulation
of it, and the tolerances are this project's own acceptance rule rather than
anything measured.

A conservation check is a right answer only when its stock side and its flux
side account for the same inventory. Two of the three closure elements in
`biosphere/config/lpj_acceptance.yaml` did not. Neither defect is a
snapshot-against-integral mismatch and neither is a tolerance question: in both
the check was differencing quantities that never spanned the same budget.

## Nitrogen: the stock omitted four of the model's six soil mineral pools

`vendor/lpj-guess/modules/ntransform.cpp` states the soil nitrogen operator's
own conservation identity over SIX pools, and enforces it every simulated day:

    n_total = NH4_mass + NO3_mass + NO2_mass + NO_mass + N2O_mass + N2_mass

`npool.out` Total carries TWO of them. `commonoutput.cpp` builds Total's
mineral term as

    availn += (NH4_mass + NO3_mass + snowpack_NH4_mass + snowpack_NO3_mass)
              * to_gridcell_average

and Total is `nmass_gridcell + nlitter_gridcell + surfsoillittern + cwdn +
centuryn + availn`. NO2, NO, N2O and N2 appear in no term of it.

They are persistent state, not a within-day intermediate. `soil.cpp` serialises
all four across the restart, and `n_gas_emission` removes only the fraction
`ftemp * (1 - wcont)` of each per day, so a remainder carries to the next day
and across the cycle boundary.

`nflux.out` NEE cannot cover the gap either, because it reports the gas when it
LEAVES those pools: `report_flux(NO_SOIL / N2O_SOIL / N2_SOIL)` is called on the
emitted increment, never on the production. Nitrogen nitrified or denitrified
out of NH4/NO3 and still resident in the other four at the cycle boundary has
therefore left the stock side and has not yet entered the flux side.

**So the residual is the change in the omitted pools.** With
`X = NO2 + NO + N2O + N2`, conservation gives
`Total[y] - Total[y-1] = -NEE[y] - (X[y] - X[y-1])`, so the check's residual is
`-(X[y] - X[y-1])`, which telescopes over any window to `-(X[end] - X[start])`.
That is bounded by the range of X and does not grow with the window, which is
exactly the signature the earlier audit measured and could not attribute: flat
in window length from ten cycles to twelve hundred.

### The test and its result

Pre-registered before measuring, with the falsification criteria fixed
(recorded on `world-w89k`): the corrected per-cycle residual
`C[y] = R[y] + (X[y] - X[y-1])` should fall inside the 1.045 kgN/ha endpoint
quantisation bound the earlier audit derived, and the four refused gridcells
should fall inside the 2.0 kgN/ha floor.

Measured over 2,024,484 gridcell-cycles and 1617 gridcells:

| statistic | uncorrected | corrected |
| --- | --- | --- |
| share of gridcell-cycles within 1.045 kgN/ha | 99.913 % | 100.000 % |
| largest per-cycle residual | 13.53 kgN/ha | 1.0019 kgN/ha |
| gridcells failing the ten-cycle window | 4 | 0 |
| largest ten-cycle window residual | 4.95 kgN/ha | 1.0000 kgN/ha |

Every one of the 2,024,484 gridcell-cycles lands inside the bound an exactly
conserving model reports at the old written precision, and the largest
corrected window residual is one endpoint quantum to four figures. The four
refused gridcells, including the two whose residual was flat in window length
and therefore real:

| lon | lat | uncorrected | corrected |
| --- | --- | --- | --- |
| -174.38 | -13.84 | +2.43 | -0.21 |
| -135.00 | -19.38 | -4.95 | -0.03 |
| -28.12 | 8.31 | +2.03 | +0.35 |
| 0.00 | -13.84 | -2.43 | +0.24 |

No gridcell inside the floor before the correction falls outside it after.

The one criterion that failed was mis-specified rather than informative. It
asked for the least-squares slope of the residual on `-(X[y] - X[y-1])` over
the gridcell-cycles selected by `|R| > 1.045`, which conditions on the
dependent variable and inflates the slope by construction; it read 1.203.
Without that selection the slope is 1.0021 over all gridcell-cycles, 1.0102
selecting on `|dX| > 1.045` and 1.0082 on `|dX| > 5.0`, and the correlation
rises to 0.998 as the omitted pools' movement grows.

NO2 is where the nitrogen sits: it accounts for 0.0377 of the 0.0400 kgN/ha
mean of X, reaches 14.58 kgN/ha, and carries a mean absolute per-cycle change
two orders above the other three. `ntransform.cpp` flushes the aerobic
fraction's nitrite into NO3 each day and leaves the anaerobic fraction's to be
drained only by the denitrification gas term, so wet soils hold it.

## Water: none of the three evaporative losses had a column the check could use

The water rule read precipitation from the pinned driver and subtracted
`aaet.out` Total and `tot_runoff.out` Total. Two of the three evaporative
losses were missing outright and the third was the wrong quantity.

**Soil evaporation and canopy interception reach no output table at all.** The
model's own annual water record adds three terms, not one:

    patch.aaet_5.add(patch.aaet + patch.aevap + patch.aintercep)

`commonoutput.cpp` accumulates `patch.aevap` and `patch.aintercep` from their
monthly arrays and sends them only to the interactive `plot()`. The two monthly
tables that do carry them, `mevap.out` and `mintercep.out`, are commented out in
every `.ins` the port ships.

**`aaet.out` Total is the transpiration of the survivors.**
`commonoutput.cpp` builds it from `indiv.aaet` summed over the individuals still
in `patch.vegetation`, and `framework.cpp` runs `vegetation_dynamics` --
establishment, mortality and disturbance by fire -- on the last day of the year,
BEFORE `output_modules.outannual`. `indiv.aaet` is reset on day zero and
accumulated daily, so an individual that transpired all year and was then killed
contributes nothing to the column. `patch.maet` is incremented by the same
`aet_total` as `patch.aaet` inside the hydrology, every day, and is unaffected;
`maet.out` is therefore the column that reports what left the soil, and it
reaches the gridcell by the same `to_gridcell_average` as `mevap`, `mintercep`
and `tot_runoff`.

So the water residual is not the implied storage change the contract calls it.
It is soil evaporation plus interception plus the transpiration of everything
that died plus the storage change, and the first two are ordinary-sized annual
fluxes.

**This has never been measured on a run.** The water block sat after carbon and
nitrogen in a loop that raised on the first failing element, so every
assessment that refused on nitrogen returned before reaching it. The refusal
recorded on `lpj_de8a0c8bd4324639b4bd62b66867d23f` names nitrogen alone for
that reason and is not a statement that water closed.

Measured now on that run, with the contract's own ten-cycle window and limits:

| statistic | value |
| --- | --- |
| gridcells failing the ten-cycle window | 1578 of 1617 |
| largest ten-cycle window residual | 3416.7 mm against a 79.5 mm limit |
| share of gridcell-cycles inside the 2.0 mm floor | 3.34 % |
| root mean square per-cycle residual | 66.0 mm |

And the residual GROWS LINEARLY with window length, which is what separates a
missing flux from a storage change. A storage term telescopes to an endpoint
difference and cannot grow; on the worst gridcell the residual runs +3417 mm at
ten cycles, +30715 at a hundred and +330709 at 1252, a steady 264 mm per cycle
on a gridcell receiving 1591 mm.

That per-cycle figure is the size the omitted losses should be. On that
gridcell the transpiration column reads 234 mm and runoff 1021 mm, so what is
missing is 17 per cent of precipitation.

### What the residual is once all four losses are subtracted

The model states the identity itself, in `soil.cpp:hydrology_lpjf`:

    initial_water_in_column + rain_melt
        = final_water_in_column + evap + aet_total + runoff

`rain_melt` is precipitation less interception, less what the snowpack takes
and plus what it releases, so over a gridcell year that identity is

    P = transpiration + soil evaporation + interception + runoff
        + change in (soil column water and ice) + change in snowpack

Four losses and nothing else, and they are exactly the ones the repaired rule
subtracts; what is left is the storage change the contract's interpretation
says it is. The check enforcing this in the model is compiled out under
`DEBUG_SOIL_WATER`, but the identity it states is the one the hydrology is
written to. So there is no fifth flux to look for -- what there was, on top of
the two with no column at all, was one of the four read off the wrong column.

Measured on `lpj_c0a9e36d42e54e7aaf245d211798e53d`, a diagnostic run of 1600
spin-up cycles and 30 retained, over 1617 gridcells:

| per-cycle residual | `aaet.out` Total and runoff only | all four losses |
| --- | --- | --- |
| population mean | 49.30 mm | 0.96 mm |
| root mean square | 65.90 mm | 6.28 mm |
| 99th percentile of magnitude | 242.6 mm | 27.4 mm |
| gridcells failing the ten-cycle window | 1584 of 1617 | 276 of 1617 |

Soil evaporation averages 36.2 mm per cycle and canopy interception 12.1,
against 96.6 of transpiration, 126.7 of runoff and 272.6 of precipitation. The
two omitted terms are 98 per cent of the mean bias.

### What the residual left after those two was

Not a store still charging. The same measurement at a quarter of the spin-up,
on `lpj_222316ba17494cb183b747c2cf0e5011` at 400 cycles, gives 306 failing
gridcells against 276 and a largest ten-cycle residual of 321 mm against 303,
so quadrupling the spin-up barely moves it. Nor is it a soil gaining carbon and
with it holding capacity under
`iforganicsoilproperties`: the per-gridcell mean residual correlates with the
soil carbon accumulation rate at +0.27, and removing that relation takes the
population root mean square from 2.52 mm to 2.43.

It is the transpiration `aaet.out` Total does not carry. At 1600 cycles the
median gridcell closes at +0.04 mm per cycle and the 564 gridcells with under
1 mm of transpiration at +0.02 mm, so the check closes wherever there is no
vegetation to lose. What remains is a strictly positive tail on the vegetated
gridcells, worth 0.7 per cent of transpiration at the median and 10.8 at the
maximum, correlating with

| against | correlation |
| --- | --- |
| burnt fraction | +0.66 |
| transpiration | +0.59 |
| leaf area index | +0.56 |
| establishment carbon flux | -0.88 |
| individual density | +0.02 |

which is turnover and not water. 1382 of 1617 gridcells carry a positive mean
residual and the most negative is -2.24 mm against a maximum of +19.14, so the
systematic part has one sign, as a column that omits a loss must.

### What the residual is once the survivors-only column is replaced too

`maet.out` in place of `aaet.out` Total closes the turnover term. Measured on
`lpj_2ffc33a5b8c749c888ca9b3a24611df4`, 1600 spin-up cycles and 30 retained
over 1617 gridcells, with the four losses read from `maet.out`, `mevap.out`,
`mintercep.out` and `tot_runoff.out`:

| statistic | value |
| --- | --- |
| gridcells failing the ten-cycle window | 81 of 1617 |
| largest ten-cycle window residual | -68.07 mm against a 17.80 mm limit |
| population mean per-cycle residual | +0.0070 mm |
| mean precipitation | 272.58 mm per cycle |
| mean transpiration, soil evaporation, interception, runoff | 97.52, 36.25, 12.10, 126.71 mm per cycle |

The correlations that identified the missing turnover are gone. Against the
per-gridcell mean residual: burnt fraction -0.028 where it was +0.66,
transpiration -0.006 where it was +0.59, leaf area index -0.008 where it was
+0.56. The sign is no longer systematic either: 901 gridcells carry a positive
mean and 716 a negative one, spanning -2.747 to +4.960 mm per cycle, where a
column omitting a loss gave 1382 positive against a most-negative -2.24.

### It is a store, and the store cannot grow

Pre-registered on `world-24xd` before the measurement: a storage term
telescopes to an endpoint difference and cannot grow with the window, and a
missing flux grows linearly. The residual was taken at every window from one to
thirty cycles, each ending on the last retained cycle, so window `w` reads
`S(end) - S(end - w)`.

| statistic | value |
| --- | --- |
| median per-gridcell slope of the residual on window length | 0.0017 mm per cycle |
| 95th percentile of that slope | 0.761 mm per cycle |
| median coefficient of determination of a linear fit in the window | 0.177 |
| mean absolute residual, one cycle to thirty | 0.625 mm to 3.249 mm |

A missing flux of rate `F` gives a mean absolute residual of `F * w` and a
coefficient of determination near one. A bounded store visited at increasing
separations grows as the square root of the window: over a thirtyfold range
that predicts a factor of 5.48 and the measurement gives 5.20, against the
factor of 30 linear growth would give. The earlier defect gave the linear
form outright, +3417 mm at ten cycles and +330709 at 1252 on the worst
gridcell, at a steady 264 mm per cycle.

The magnitudes are a store's as well. Every one of the 81 gridcells refused at
the ten-cycle window sits inside the water its own soil column can hold, at
0.453 of the available water capacity at worst and 0.112 at the median; the
capacities come from the land column property contract, per gridcell, through
the driver. Across all 48,510 gridcell-windows the residual exceeds the
column's saturation capacity exactly once, at (-5.62, 80.27), 148.79 mm against
148.0 mm at thirty cycles -- a gridcell at 80 degrees north, where the store
includes a snowpack that a soil column capacity does not cover.

### So the tolerance was judging a store, and the window decided the verdict

The check subtracted four losses from precipitation and held what was left to
`max(2.0 mm, 0.005 * window precipitation)`. What was left is the store, and the
store is a few hundred millimetres that moves between cycles on its own. Under
a repeating forcing with stochastic patch dynamics that movement is internal
variability, not disequilibrium, and it does not shrink with spin-up:
quadrupling the spin-up from 400 cycles to 1600 moved the failing gridcells
from 306 to 276 on the previous form.

The scaling gives it away without needing any of that. The limit grows in
proportion to the window and a bounded store grows as its square root, so
refusals fall off as the window lengthens:

| window, complete forcing cycles | 1 | 5 | 10 | 20 | 30 |
| --- | --- | --- | --- | --- | --- |
| gridcells refused | 113 | 110 | 81 | 48 | 22 |
| mean absolute residual, mm | 0.625 | 1.746 | 2.485 | 3.088 | 3.249 |

A conservation test whose refusal count falls as it is given more record is
answering a question about the declared window length. `failure-modes.md` class
34: the instrument was checked against the size of the effect and it is not
measuring conservation at all.

### The repair: the store is read

`vendor/lpj-guess/modules/commonoutput.cpp` writes `awater.out`. It carries the
gridcell water store at the end of each simulation year -- available soil
water as `wcont` times each layer capacity, soil ice as its volume fraction
times each layer thickness, and the snowpack, which is already a rainfall
equivalent -- in exactly the terms `soil.cpp:hydrology_lpjf` balances the
column in, and on the same `to_gridcell_average` as `tot_runoff.out`. Water
below the wilting point is absent from it for the same reason it is absent
from that balance: it is constant and cancels in the difference of two
endpoints.

The water residual is then

    (store at the end of the window - store at the end of the year before it)
        - (precipitation - transpiration - soil evaporation
           - interception - runoff)

over the same whole number of forcing cycles, which is zero when the model
conserves water. Both sides are read, neither is inferred, and the tolerance
did not move: the floor stays 2.0 mm and the relative limit 0.005.

### What the repaired water closure measures

`lpj_ed9bb44ab2b746b09de4d695651a427d`, the same 1600 spin-up cycles, 30
retained and 1617 gridcells, with `awater.out` retained and the residual taken
stock against flux. Every window from one cycle to twenty-nine:

| statistic | value |
| --- | --- |
| gridcells failing, at every window from 1 to 29 cycles | 0 of 1617 |
| largest residual at the ten-cycle window | 0.0276 mm against a 5.88 mm limit |
| written-precision resolution bound | 0.1806 mm |
| population mean residual, one cycle | -0.00003 mm |
| population mean residual, twenty-nine cycles | -0.00224 mm |

The largest residual anywhere is an order of magnitude INSIDE the residual an
exactly conserving model reports at the written precision of the columns
differenced, and the mean drifts by -7.7e-5 mm per cycle against the 0.018 mm
per cycle the thirty-six monthly columns' rounding allows. Carbon and nitrogen
close on the same run at 3.2e-5 kgC/m2 and 0.00104 kgN/ha against floors of
0.01 and 2.0, with no gridcell refused.

So the model conserves water and the previous residual was the store, entire.

### The snowpack reaches the model's own ceiling on 17 gridcells

Reading the store made this visible; it is not a closure failure and it is not
new. `soilwater.cpp:snow` caps the pack at `SNOWPACK_MAX`, 10000 mm of rainfall
equivalent, by taking `melt = -min(prec, SNOWPACK_MAX - snowpack)`. At the cap
that term is zero, the day's snowfall goes to `rain_melt` instead, and the soil
receives it as liquid water: conserving, which is why the closure above passes,
and unphysical, because the gridcell is below freezing.

Measured: 17 of 1617 gridcells sit at the cap in all 30 retained cycles, and 27
carry a pack over 1000 mm, spanning 58 degrees south to 80 degrees north. They
are gridcells whose snow accumulation exceeds its melt every cycle -- permanent
ice in everything but name -- and LPJ-GUESS has no ice sheet, so the pack
saturates and the surplus is delivered to the soil column as rain. The pack is
also what put the one gridcell-window in the previous section outside its soil
column's saturation capacity.

## Whether the repaired checks still resolve their own tolerances

Completing a stock or a loss adds every new column's written quantum to the
residual an exactly conserving model reports, so the question
`notes/audits/closure-tolerance-under-written-precision.md` settled has to be
asked again of the wider form. Measured from
`lpj_c0a9e36d42e54e7aaf245d211798e53d`, which carries the precisions
`commonoutput.cpp` writes today:

| element | columns differenced | resolution bound | floor | margin |
| --- | --- | --- | --- | --- |
| nitrogen | npool.out Total plus four soil_npool.out pools against nflux.out NEE | 0.001445 kgN/ha | 2.0 kgN/ha | 1384 |
| carbon | cpool.out Total against cflux.out NEE | 4.6e-5 kgC/m2 | 0.01 kgC/m2 | 217 |
| water | awater.out Total, thirty-six monthly loss columns and tot_runoff.out Total | 0.1806 mm | 2.0 mm | 11.07 |

All three clear the contract's required tenfold margin. Water is the tight one,
because its thirty-six columns are written at the three decimals mainline gives
a monthly table rather than at the `closure_prec_water` the annual closure
columns were raised to. The margin is measured from the artifact on every
assessment, so a precision that drops back refuses rather than passing quietly.
