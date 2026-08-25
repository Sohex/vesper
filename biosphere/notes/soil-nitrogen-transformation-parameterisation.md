# The soil nitrogen transformation operator, and what it cannot see

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf, and its biosphere is LPJ-GUESS-CNP, a vegetation model written for
Earth and ported to this world's calendar. Everything below is about that
model's soil nitrogen chemistry: the constants it runs on, the response
functions those constants sit in, and what the simulated soil's gas environment
does and does not reach.

`ifntransform 1` is on in the Vesper baseline, through `global.ins` importing
`global_soiln.ins`, and it is not an emissions diagnostic. Nitrification,
denitrification and ammonia volatilisation together set how much mineral
nitrogen the simulated plants can reach, so the operator is on the path from
litter to productivity whether or not anything reads its gas fluxes.

Nothing here is verified by execution. LPJ-GUESS does not build on this tree:
`framework/vesper.h` is generated, and the chain to it runs through a baseline
run and a baseline climatology that do not exist. Every statement below is
against the source, and every bound is arithmetic on the declared forms, which
`biosphere/scripts/ntransform_gate.py` re-derives on each invocation.

## What the operator is

`modules/ntransform.cpp` is the Xu-Ri and Prentice (2008) soil nitrogen scheme
as ported into LPJ-GUESS. It runs once per simulation day, after the soil
organic matter submodel, in five steps: ammonia volatilisation, a split of every
mineral pool into an aerobic and an anaerobic share, nitrification of the
aerobic share, denitrification of the anaerobic share, and the release of NO,
N2O and N2 to the atmosphere. Six pools carry nitrogen through it -- NH4, NO3,
NO2, NO, N2O and N2 -- and its conservation identity is that what those pools
hold before the operator runs equals what they hold after plus what it emitted.

Every response function in it is a function of soil temperature, upper-soil
water content, water-filled pore space or pH, and of nothing else.

## The register

| constant | where the model reads it | value | unit | source or bracket |
| --- | --- | --- | --- | --- |
| `f_nitri_max` | `global_soiln.ins`, then `nitrification` | 0.1 | fraction of aerobic NH4 per day | Xu-Ri and Prentice (2008) table 8 eqn 1. UNVERIFIED: the paper could not be obtained |
| `f_nitri_gas_max` | `global_soiln.ins`, then `denitrification` | 0.25 | fraction of NO2 per day | the same, table 9 eqn 4. UNVERIFIED |
| `f_denitri_max` | `global_soiln.ins`, then `denitrification` | 0.33 | fraction of anaerobic NO3 per day | the same, table 9 eqn 3. UNVERIFIED |
| `f_denitri_gas_max` | `global_soiln.ins`, then `nitrification` | 0.33 | fraction of gross nitrification lost as gas | the same, table 8. UNVERIFIED. Note the crossed naming: the constant named for denitrification is read in nitrification and the one named for nitrification is read in denitrification |
| `k_N` | `global_soiln.ins`, then `denitrification` | 0.083 | kgN/m3 half-saturation | the same, table 9 eqns 3 and 4. UNVERIFIED |
| `k_C` | `global_soiln.ins`, then `denitrification` | 0.017 | kgC/m3 half-saturation | the same, table 9 eqn 2. UNVERIFIED |
| `nh3_max` | a literal in `nh3_volatilization` | 0.001 above pH 6, 0.00001 at or below | fraction of NH4 per day | the same, table 5. UNVERIFIED, and the pH branch point is a step rather than a curve |
| the volatilisation, nitrification and emission temperature responses | literals in `ntransform.cpp` | see the declaration | dimensionless | the same, tables 5, 8, 9 and 10. UNVERIFIED |
| the aerobic/anaerobic split, and the NO share of nitrification gas | literals in `substrate_partition` and `nitrification` | see the declaration | dimensionless | Pilegaard (2013), READ |
| the denitrification moisture response, and the N2O:NO and N2:N2O ratios | literals in `denitrification` | see the declaration | dimensionless | Weier et al. (1993). UNVERIFIED: the paper could not be obtained |
| `pH_soil` | nowhere | 3.5 to 8.5 declared | pH | a declared instruction parameter that no line of the model reads |

`biosphere/config/ntransform.yaml` is the machine-readable form of this table
and of every response function, and the gate checks it against
`global_soiln.ins` and `modules/ntransform.cpp` on each run. A value here that
has drifted from the source is a failure of the declaration.

## The bracket cannot be closed, and the reason is a paper

Xu-Ri and Prentice (2008), *Terrestrial nitrogen cycle simulated with a dynamic
global vegetation model*, Global Change Biology 14(8), 1745-1764,
`10.1111/j.1365-2486.2008.01625.x`, is the source of every table the operator
cites and of every constant marked UNVERIFIED above. It could not be obtained.
Weier et al. (1993), *Denitrification and the dinitrogen/nitrous oxide ratio as
affected by soil water, available carbon, and nitrate*, Soil Science Society of
America Journal 57(1), 66-72,
`10.2136/sssaj1993.03615995005700010013x`, is the source of the denitrification
moisture response and its two gas ratios, and it could not be obtained either.

So the constants cannot be checked against what they were calibrated on, and the
Earth-calibrated response-function bracket this operator needs cannot be
declared: a bracket requires knowing the range of soils and conditions each
curve was fitted over, and that is in those two papers. What is registered here
instead is the FORM of every response, its domain, and the bound each one holds
over that domain, which is arithmetic and needs no paper.

## What holds without the papers: the operator cannot create nitrogen

Every response function above multiplies a nitrogen pool mass, directly or
through a chain, so each has to lie in [0, 1] over the whole range of soil
states this world produces or the operator makes or destroys nitrogen. The gate
samples each function over its declared domain and each chain as a product of
maxima. All of them hold, and three are worth stating because they hold for a
reason rather than by construction.

**The nitrification water response is bounded by a crossing, not by a clamp.**
`nit_act` is the smaller of a rising exponential in water content and the
falling line `4 - 5w`. The falling limb alone reaches 4.0 at zero water. The
product stays at or below 1.0 only because the rising limb is below 1 everywhere
the falling limb is above it, and the two cross at exactly 1.0. Nothing in the
code enforces that, so it is a property of the two constants together and it is
what the gate re-derives.

**Volatilisation is bounded by pedology's pH range.** The ammonia multiplier
carries `exp(2 * (pH - 10))`, which grows without bound in pH, and above about
pH 14.15 it would exceed 1 and be held back only by the `min(NH4_mass, ...)`
clamp beside it. It does not get there: `pedology/config/pedogenesis.yaml` clips
its pH field, and over that range the whole multiplier reaches at most 9.2e-5,
four orders of magnitude below 1. The bound is therefore a joint property of the
operator and of the soil map that feeds it, which is why the gate declares the pH
domain as pedology's and not as chemistry's.

**These bounds are per unit pool mass, so they survive the stoichiometry
repair.** The fine-root and sapwood C:N windows were re-anchored on the tissue
mean under the plant physiology audit's finding 11, and `somdynam.cpp` builds
fine-root litter lignin-to-nitrogen from `cton_root_avr`, so root litter now
enters the simulated soil at 1.79 times its former nitrogen concentration and
mineralises faster. That moves the SIZE of every pool this operator acts on and
therefore every absolute flux it produces. It moves none of the bounds above,
because each is a bound on a coefficient. Any absolute nitrogen number from this
operator has to be re-derived on the repaired tree; the conservation and range
properties do not.

## The aerobic and anaerobic split is applied to the wrong variable

`substrate_partition` divides every mineral nitrogen pool between an aerobic and
an anaerobic share, and that split decides how much of the soil's nitrogen is
offered to nitrification and how much to denitrification. It is the single
largest control the operator has on gaseous loss.

It reads `soil.get_soil_water_upper()`, the upper layer's water as a fraction of
its AVAILABLE capacity. Its own comment describes the curve in water-filled pore
space -- denitrifier activity rising from about 50 per cent, steepest near 66 per
cent, complete by 75 per cent -- and cites Pilegaard (2013), which puts the
crossover from nitrification to denitrification at 60 per cent WFPS. Those are
two different quantities. `Soil::wfps` builds water-filled pore space as
`(wcont * gawc + gwp) / gwsats`, so it carries an offset and a scale set by the
soil's own wilting point and saturation capacity, and both are texture
dependent. `nitrification` and `denitrification` in the same file read
`soil.wfps(0)`; only this function does not.

For the medium soil, LPJ code 2, `gawc` is 75 mm, `gwp` 92 mm and `gwsats`
219.5 mm over the upper layer, so a water content of 0.5 is 0.590 in water-filled
pore space, and the 0.66 the comment names as the steepest point corresponds to a
water content of 0.705. The curve is therefore evaluated well below the point its
description places it at, by an amount that changes with texture. The curve the
comment describes is not this curve in either variable either: `richards_curve`
with a midpoint of 0.5 and a shape parameter of 7.5 reaches 0.83 at 0.75, not 1.

This is not fixed here. Substituting `soil.wfps(0)` is a one-line change that
moves the split on every gridcell, and the calibration it should be checked
against is in the Xu-Ri and Prentice (2008) paper that could not be obtained. It
is registered and filed rather than guessed at.

## Three defects the audit's finding 8 did not have

**The simulated soil's pH never reached the operator.** Finding 8 says pedology
pH reaches the code. It did not. `SoilInput::get_mineral` reads the soil map's
`ph` column into a local `SoilProperties`, and `get_lpj` sets 6.5 there, but
neither `get_soil_mineral` nor `get_soil_organic` copied it onto
`gridcell.soiltype.pH`, which stayed at the `-1.0` its constructor sets. So
`nh3_volatilization` took its no-pH branch on every gridcell.

**And that branch could not work.** It evaluated Dawson (1977)'s regression of
soil pH on annual precipitation, `3810 / (762 + climate.aprec_lastyear) + 3.5`.
Nothing in the model assigns `aprec_lastyear`; `climate.aprec`, which would feed
it, is reset at day 0 in `driver.cpp` and never accumulated. The regression
therefore evaluated at zero precipitation and returned 8.5 for every gridcell on
every day, which put every soil on the `nh3_max = 0.001` branch and fixed the pH
term at `exp(-3)`. Two further things stood between that fallback and a usable
number even with a live input: it is an Earth calibration, and its argument is
an annual precipitation sum, which on this world's shorter year is a smaller
number for the same precipitation rate and so reads as a drier, more alkaline
soil. Both are now repaired: the soil map's pH is carried onto the `Soiltype`,
and an input path that supplies none makes the operator refuse by name rather
than substitute a constant nobody chose.

`climate.aprec` remains a dead accumulator, written to `out_seasonality` as a
permanent zero. That is `driver.cpp`'s to fix and is filed separately.

**The operator's only conservation check did not exist in the binary.** The
mass-balance test at the end of `ntransform()` was an `assert`, and every
Release build compiles those out under NDEBUG. It is now a compiled-in test
against a bar relative to the nitrogen present, floored at the absolute epsilon
the assert used. The bar was fixed before any result, because there is no result
yet.

## What the operator still cannot see

Five preconditions carry the `undeclared` sentinel in
`biosphere/config/ntransform.yaml`, and `--strict` refuses while any of them
does. None of them has a default.

- **Surface pressure.** No response function takes one. The operator gives the
  same rates at this world's surface pressure as at Earth's.
- **Oxygen partial pressure.** Nitrification is obligately aerobic and
  denitrification is what happens without oxygen, so the split between them is a
  function of oxygen supply. `substrate_partition` substitutes water-filled pore
  space for it, through a curve fitted to Earth soils at Earth's oxygen partial
  pressure. `config/planet.yaml` declares this world's pO2; nothing in this
  operator reads it.
- **Soil gas diffusivity.** The water-filled-pore-space responses stand in for
  diffusion through the pore network, which depends on the gas, the pore
  geometry and the pressure. The operator carries none of the three.
- **Water table and redox state.** Every water term is the uppermost Gerten
  layer. A water table, and the reduced zone below it, has no representation.
- **The atmospheric boundary.** The NO, N2O and N2 the operator emits leave at
  the soil surface and reach no atmosphere.

The declared model boundary that follows is the same shape as the non-N/P
adequacy screen's: any nitrogen-limitation result from this operator is a result
for an Earth gas and redox environment driven by this world's water and pH, and
not a result for this world's atmosphere.

## Primary sources read

- Pilegaard (2013), *Processes regulating nitric oxide emissions from soils*,
  Phil. Trans. R. Soc. B 368(1621), 20130126, `10.1098/rstb.2013.0126`.
