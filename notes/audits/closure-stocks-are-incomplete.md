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

## Water: the losses omitted two of the three evaporation terms

The water rule reads precipitation from the pinned driver and subtracts
`aaet.out` Total and `tot_runoff.out` Total. `aaet.out` Total is TRANSPIRATION:
`commonoutput.cpp` builds it from `indiv.aaet` summed over individuals and
patches. The model's own annual water record adds three terms, not one:

    patch.aaet_5.add(patch.aaet + patch.aevap + patch.aintercep)

`patch.aevap`, bare soil evaporation, and `patch.aintercep`, canopy
interception evaporation, reach no output table at all. `commonoutput.cpp`
accumulates both from their monthly arrays and sends them only to the
interactive `plot()`. The two monthly tables that do carry them, `mevap.out`
and `mintercep.out`, are commented out in every `.ins` the port ships.

So the water residual is not the implied storage change the contract calls it.
It is soil evaporation plus interception plus the storage change, and the first
two are ordinary-sized annual fluxes.

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

That per-cycle figure is the size the two omitted fluxes should be. On that
gridcell transpiration is 234 mm and runoff is 1021 mm, so the missing term is
17 per cent of precipitation and would put the transpired share of total
evaporation at 47 per cent.

### What the residual is once all four losses are subtracted

The model states the identity itself, in `soil.cpp:hydrology_lpjf`:

    initial_water_in_column + rain_melt
        = final_water_in_column + evap + aet_total + runoff

`rain_melt` is precipitation less interception, less what the snowpack takes
and plus what it releases, so over a gridcell year that identity is

    P = transpiration + soil evaporation + interception + runoff
        + change in (soil column water and ice) + change in snowpack

The four losses are exactly the ones the repaired rule subtracts, and what is
left is the storage change the contract's interpretation says it is. There is
no third missing flux: the check enforcing this in the model is compiled out
under `DEBUG_SOIL_WATER`, but the identity it states is the one the hydrology
is written to.

Measured on `lpj_222316ba17494cb183b747c2cf0e5011`, a diagnostic run of 400
spin-up cycles and 30 retained, over 1617 gridcells:

| per-cycle residual | transpiration and runoff only | all four losses |
| --- | --- | --- |
| population mean | 48.48 mm | 1.01 mm |
| root mean square | 63.75 mm | 6.52 mm |
| 99th percentile of magnitude | 223.5 mm | 27.4 mm |
| gridcells failing the ten-cycle window | 1587 of 1617 | 306 of 1617 |

Soil evaporation averages 35.6 mm per cycle and canopy interception 11.9,
against 96.8 of transpiration, 127.3 of runoff and 272.6 of precipitation. The
two omitted terms are 98 per cent of the mean bias.

The 306 that still fail are stores that have not finished charging on a run
deliberately cut to a twentieth of the derived spin-up, and the accumulation is
not the one a soil gaining carbon under `iforganicsoilproperties` would show:
the per-gridcell mean residual correlates with the soil carbon accumulation
rate at only +0.27, and removing that relation takes the population root mean
square from 2.52 mm to 2.43. The largest single residual belongs to a gridcell
at 52.61 degrees with no soil carbon at all and a residual constant to four
figures across both halves of the record, which is a snowpack filling at a
fixed rate toward the model's 10000 mm cap.
