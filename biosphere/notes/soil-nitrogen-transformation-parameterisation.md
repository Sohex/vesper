# The soil nitrogen transformation operator, and what it cannot see

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf, and its biosphere is LPJ-GUESS-CNP, a vegetation model written for
Earth and ported to this world's calendar. Everything below is about that
model's soil nitrogen chemistry: the constants it runs on, the response
functions those constants sit in, what each was calibrated against, and what the
simulated soil's gas environment does and does not reach.

`ifntransform 1` is on in the Vesper baseline, through `global.ins` importing
`global_soiln.ins`, and it is not an emissions diagnostic. Nitrification,
denitrification and ammonia volatilisation together set how much mineral
nitrogen the simulated plants can reach, so the operator is on the path from
litter to productivity whether or not anything reads its gas fluxes.

Nothing here is verified by execution. LPJ-GUESS does not build on this tree:
`framework/vesper.h` is generated, and the chain to it runs through a baseline
run and a baseline climatology that do not exist. Every statement below is
against the source and against the two papers, and every bound is arithmetic on
the declared forms, which `biosphere/scripts/ntransform_gate.py` re-derives on
each invocation.

## What the operator is

`modules/ntransform.cpp` is the Xu-Ri and Prentice (2008) soil nitrogen scheme,
which that paper calls DyN, as ported into LPJ-GUESS. It runs once per
simulation day, after the soil organic matter submodel, in five steps: ammonia
volatilisation, a split of every mineral pool into an aerobic and an anaerobic
share, nitrification of the aerobic share, denitrification of the anaerobic
share, and the release of NO, N2O and N2 to the atmosphere. Six pools carry
nitrogen through it -- NH4, NO3, NO2, NO, N2O and N2 -- and its conservation
identity is that what those pools hold before the operator runs equals what they
hold after plus what it emitted.

Every response function in it is a function of soil temperature, upper-soil
water content, water-filled pore space or pH, and of nothing else.

## The register: every constant against its calibration

Both source papers are now held and read. Xu-Ri and Prentice (2008) is the
source of every table the operator cites; Weier et al. (1993) is the source of
its denitrification moisture response and of its dinitrogen to nitrous oxide
ratio. Each row below gives what the paper states, what the code has, and
whether they agree.

| what | where the model reads it | value | source, and what it says | verdict |
| --- | --- | --- | --- | --- |
| `f_nitri_max` | `global_soiln.ins`, then `nitrification` | 0.1 /day | Xu-Ri table 11, Nmax, from Khalil et al. (2004), stated AT 20 C | agrees in value; the paper applies it at the peak of a curve normalised to 38 C |
| `k_N` | `global_soiln.ins`, then `denitrification` | 0.083 kgN/m3 | Xu-Ri table 11, Kn, from Li et al. (1992) | agrees, identical |
| `k_C` | `global_soiln.ins`, then `denitrification` | 0.017 kgC/m3 | Xu-Ri table 11, Kc, from Li et al. (1992) | agrees, identical |
| `f_denitri_gas_max` | `global_soiln.ins`, then `nitrification` | 0.33 | Xu-Ri table 11, RNON 0.1-4% and RN2ON under 0.1-0.2%, so 0.2-4.2% together | OUTSIDE, by 7.9x on the top of the range |
| `f_denitri_max` | `global_soiln.ins`, then `denitrification` | 0.33 | Xu-Ri table 9 eqn 3 has no maximum-rate constant, and table 11 lists none | unsourced |
| `f_nitri_gas_max` | `global_soiln.ins`, then `denitrification` | 0.25 | Xu-Ri table 9 eqn 4, the same | unsourced |
| the volatilisation temperature response | `nh3_volatilization` | -- | Xu-Ri table 5 eqn 7 | agrees, identical including the min |
| the gas emission multiplier | `n_gas_emission` | -- | Xu-Ri table 10 eqns 1 and 2 | agrees, identical |
| the nitrification temperature response | `nitrification` | -- | Xu-Ri table 8 eqn 2 | agrees; the port's min is a no-op, the function peaks at exactly 1 at 38 C |
| the denitrification temperature response | `denitrification` | -- | Xu-Ri table 9 eqn 1, which carries NO min | OUTSIDE: the port clamps it at its 22 C value |
| `nh3_max` | a literal in `nh3_volatilization` | 0.001 above pH 6, 0.00001 at or below | Xu-Ri table 5 has no such factor; eqn 6's fpH already carries the ratio | OUTSIDE: the pH dependence is applied twice |
| the aerobic/anaerobic split | `substrate_partition` | midpoint 0.5, shape 7.5 | Xu-Ri table 7 settles the VARIABLE, water-filled pore space, and gives no form | unsourced form, bracket 0.50 to 0.66 on the midpoint |
| the nitrification water response | `nitrification` | -- | Xu-Ri table 8 eqn 1 has no moisture term | unsourced, and a second moisture control |
| the NO share of nitrification gas | `nitrification` | -- | Xu-Ri table 11 RNON against RN2ON brackets it at 0.33 to 0.98 | agrees, 0.50 to 0.69 over the reachable WFPS |
| the denitrification moisture response | `denitrification` | exponent 13.036 | Weier table 2, at 60, 75 and 90% WFPS | agrees, inside the 8.09 to 14.06 that table brackets |
| the N2 share above 0.7 WFPS | `denitrification` | -- | Weier tables 4 and 5 | OUTSIDE below about 0.85 WFPS |
| the N2 branch point at 0.7 WFPS | `denitrification` | -- | Weier tables 4 and 5 report N2 at 60 and 75% WFPS | OUTSIDE: the operator produces none there |
| the N2O:NO ratio | `denitrification` | -- | attributed to Weier, which never measured NO | unsourced |
| the N2:N2O temperature sigmoid | `denitrification` | midpoint 5 C, width 10 C | in neither paper; Weier incubated at one temperature | unsourced; its sign agrees with Xu-Ri table 9 eqns 5 and 6 |
| the Michaelis-Menten divisor | `denitrification` | available-water depth | Xu-Ri table 11 gives kg m-3 of an unstated volume | unsourced |
| `pH_soil` | nowhere | 3.5 to 8.5 declared | -- | a declared instruction parameter that no line of the model reads |

`biosphere/config/ntransform.yaml` is the machine-readable form of this table
and of every response function, and the gate checks it against
`global_soiln.ins` and `modules/ntransform.cpp` on each run. A value here that
has drifted from the source is a failure of the declaration.

## The Earth-calibrated response bracket

Twelve of the twenty declared entries are what remain undeclared, and
`--strict` refuses on exactly those and on the five Vesper preconditions. Six
entries carry a bracket, three of them among the twelve; the rest carry none
because neither paper states the quantity, and no central value is fitted here
to stand in for one.

**The nitrification rate constant.** Xu-Ri table 11 gives Nmax as 0.1 per day
from Khalil et al. (2004) and states it at 20 C, while table 8 eqn 2's
temperature response is 1 at 38 C, which is where the model applies it. That
response is 0.248 at 20 C, so the operator rates 20 C nitrification at 0.0248
per day against the 0.1 per day Khalil reports. The bracket on the constant is
therefore 0.1 to 0.403 per day, the value as the paper writes it against the
value that would honour the measurement it cites. The inconsistency is Xu-Ri's,
carried faithfully into the port.

**The gas fraction of gross nitrification.** Xu-Ri table 11 gives RNON as 0.1
to 4 per cent with a mean of 2, and RN2ON as under 0.1 to 0.2 per cent. The
bracket on their sum is 0.2 to 4.2 per cent. The operator's 0.33 is outside it
by a factor of 7.9 on the top and 15 on the mean: of every kilogram of nitrogen
nitrified the model routes 330 grams to NO and N2O where the paper routes 2 to
42 grams. It lands on `NET_NITRIF`, `NO_SOIL` and `N2O_SOIL` directly, and it is
the largest single constant disagreement in the operator.

**The NO share of that gas.** RNON against RN2ON brackets the NO share at 0.33
to 0.98. The operator's curve gives 0.50 to 0.69 over the water-filled pore
space the model can reach, inside the bracket, and 0.506 at 0.60 WFPS, which is
Pilegaard's NO:N2O ratio near 1 at that crossover.

**The denitrification moisture response.** `min(1, exp(13.036 * WFPS -
11.6219))` is log-linear, so it asserts one constant factor per 15 WFPS points,
`exp(13.036 * 0.15) = 7.07`. Weier table 2 gives denitrification rates for four
soils and nine carbon-by-nitrate treatments at each of 60, 75 and 90 per cent
WFPS. Its geometric means step by 6.35 from 60 to 75 and by 8.23 from 75 to 90;
its arithmetic means step by 6.63 and 3.36. Expressed as the exponent those
extremes bracket it at 8.09 to 14.06 per unit WFPS, and 13.036 is inside. The
offset puts the function at 1 at 0.8915 WFPS, normalising it to Weier's wettest
treatment.

**The N2 share of denitrification gas.** Weier tables 4 and 5 give the N2/N2O
ratio for all four soils, nine treatments and five days at 25 C. Converted to an
N2 share by the paper's own footnote, `1 - 1/(1 + N2/N2O)`:

| WFPS | observations | N2 share, median | interquartile | Xu-Ri | the operator at 25 C |
| --- | --- | --- | --- | --- | --- |
| 60% | 109 | 0.565 | 0.333 to 0.778 | 0.978 | 0.000 |
| 75% | 129 | 0.744 | 0.545 to 0.870 | 0.978 | 0.120 |
| 90% | 163 | 0.796 | 0.630 to 0.938 | 0.978 | 0.846 |

The bracket is the union of those interquartile ranges, 0.333 to 0.938. Xu-Ri's
own tables 9 and 11 sit outside it at the high end, because eqns 5 to 7 make NO
and N2O small fixed fractions of the reduction flux and leave everything else as
N2. The operator sits outside it at the low end everywhere but the wet end: it
reproduces Weier at 0.90 WFPS and inverts the trend below it, reaching a hard
zero below 0.70 where Weier's medians are 1.3 and 2.9 and where N2 is the
majority product in 61 and 80 per cent of the observations.

**The aerobic and anaerobic split.** Xu-Ri settles the variable and gives no
form, pointing instead at PnET-DNDC's anaerobic balloon, so the bracket is on
the curve's midpoint and it is the spread between the two numbers the
operator's own comment asserts at once: Pilegaard's crossover at 0.60 WFPS and
the comment's own "steepest change around 66%". That is 0.50 to 0.66, with the
declared 0.5 at the bottom edge. The shape parameter has no source at all.

## The split is applied to water-filled pore space, and the paper says so

`substrate_partition` divides every mineral nitrogen pool between an aerobic and
an anaerobic share, and that split decides how much of the soil's nitrogen is
offered to nitrification and how much to denitrification. It is the single
largest control the operator has on gaseous loss.

Xu-Ri is unambiguous about the variable. Table 7, comparing DyN with DNDC, gives
the "soil aeration status used to allocate substrates" as "soil moisture as
water-filled pore space (WFPS)", against DNDC's redox-based anaerobic volume
fraction. The nitrification section repeats it: "Soil moisture content (WFPS) of
the top 50 cm layer is used to allocate substrates into the aerobic and
anaerobic soil fractions." `SOILDEPTH_UPPER` is 500 mm, so the layer is the
right one, and `Soil::wfps(0)` is that layer's water-filled pore space.

The function read `soil.get_soil_water_upper()`, the same layer's water as a
fraction of its AVAILABLE capacity, which is a different quantity.
`Soil::wfps` builds water-filled pore space as `(wcont * gawc + gwp) / gwsats`,
so the conversion carries an offset and a scale set by the soil's own wilting
point and saturation capacity, and both are texture dependent. `nitrification`
and `denitrification` in the same file already read `soil.wfps(0)`.

What that cost, on the current Vesper soil map's 4105 gridcells: at a fixed 0.50
WFPS, cells got an anaerobic share anywhere from 0.071 to 0.874 depending only
on their sand and clay, where the curve on WFPS gives 0.500 everywhere. The
median cell's aerobic share, which is what nitrification is offered, was 1.65
times what the curve asks for at 0.50 WFPS and 1.55 times at 0.60. The
disagreement is not a uniform bias: the curve's midpoint landed at a median of
0.606 WFPS, essentially Pilegaard's crossover, but spread from 0.372 to 0.806
across textures. The operator now reads `soil.wfps(0)`, which removes the
texture dependence. The midpoint and shape are separate and unsourced, and are
bracketed above.

**The conversion exists, but it does not give the curve its whole domain.**
`get_soil_water` refuses a `wcont` outside [0, 1], and `wcont` is a fraction of
available capacity, so `wfps(0)` is bounded below by wilting point over
saturation and above by field capacity over saturation. On the current soil map
that is 0.171 to 0.665 at the dry end and 0.573 to 0.948 at the wet end, a
median wet limit of 0.813. The consequences are structural, not marginal: 164 of
4105 gridcells can never reach 0.70 WFPS, so the operator's N2 branch is dead
on them and every denitrified atom leaves as NO or N2O; and 2032 of them are
above 0.40 WFPS even at wilting point, so denitrification never switches off.
The quantity the model does not carry is **soil water above field capacity in
the upper 50 cm**. Xu-Ri's WFPS is an aeration proxy that means something at
saturation; this model's upper layer stops at field capacity and drains the rest
to runoff within the day.

## What holds without any paper: the operator cannot create nitrogen

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

**Warm-soil denitrification is bounded by a clamp the paper does not have.**
Xu-Ri table 9 eqn 1 is an exponential equal to 1 at 22 C and rising past it:
1.21 at 25 C, 2.07 at 35 C, 3.15 at 45 C. Table 5 eqn 7 and table 10 eqn 1 both
state a `min{1, ...}`; table 9 eqn 1 does not. The port clamps it anyway, and
the conservation identity needs something to, because the coefficient multiplies
a pool mass. The bound is real and its source is the port.

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

## Three defects the audit's finding 8 did not have

**The simulated soil's pH never reached the operator.** Finding 8 says pedology
pH reaches the code. It did not. `SoilInput::get_mineral` reads the soil map's
`ph` column into a local `SoilProperties`, and `get_lpj` sets 6.5 there, but
neither `get_soil_mineral` nor `get_soil_organic` copied it onto
`gridcell.soiltype.pH`, which stayed at the `-1.0` its constructor sets. So
`nh3_volatilization` took its no-pH branch on every gridcell.

**And that branch could not work.** It evaluated Dawson (1977)'s regression of
soil pH on annual precipitation, `3810 / (762 + climate.aprec_lastyear) + 3.5`.
The port's constant was wrong as well: Xu-Ri table 5 eqn 5 gives `+ 3.8`, and
since pH enters volatilisation as `exp(2 * (pH - 10))`, 0.3 pH units are a
factor of 1.8 on the flux. Nothing in the model assigns `aprec_lastyear`;
`climate.aprec`, which would feed it, is reset at day 0 in `driver.cpp` and never
accumulated. The regression therefore evaluated at zero precipitation and
returned 8.5 for every gridcell on every day, where the paper's constant would
have given 8.8. That put every soil on the `nh3_max = 0.001` branch and fixed the
pH term at `exp(-3)`. Two further things stood between that fallback and a usable
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
does. None of them has a default, and neither paper settles any of them: both
are Earth work at Earth's surface pressure and oxygen partial pressure, and
neither carries a water table. Xu-Ri treats soil water in two sublayers of 0.5
and 1.0 m; Weier incubated 60 g of repacked soil in a sealed jar.

- **Surface pressure.** No response function takes one. The operator gives the
  same rates at this world's surface pressure as at Earth's.
- **Oxygen partial pressure.** Nitrification is obligately aerobic and
  denitrification is what happens without oxygen, so the split between them is a
  function of oxygen supply. `substrate_partition` substitutes water-filled pore
  space for it, as Xu-Ri table 7 does, through a curve fitted to Earth soils at
  Earth's oxygen partial pressure. `config/planet.yaml` declares this world's
  pO2; nothing in this operator reads it.
- **Soil gas diffusivity.** The water-filled-pore-space responses stand in for
  diffusion through the pore network, which depends on the gas, the pore
  geometry and the pressure. The operator carries none of the three.
- **Water table and redox state.** Every water term is the uppermost Gerten
  layer, which stops at field capacity. A water table, and the reduced zone
  below it, has no representation.
- **The atmospheric boundary.** The NO, N2O and N2 the operator emits leave at
  the soil surface and reach no atmosphere.

The declared model boundary that follows is the same shape as the non-N/P
adequacy screen's: any nitrogen-limitation result from this operator is a result
for an Earth gas and redox environment driven by this world's water and pH, and
not a result for this world's atmosphere.

## Primary sources read

- Xu-Ri and Prentice (2008), *Terrestrial nitrogen cycle simulated with a
  dynamic global vegetation model*, Global Change Biology 14(8), 1745-1764,
  `10.1111/j.1365-2486.2008.01625.x`. Tables 5, 8, 9, 10 and 11, and the
  nitrification, denitrification, volatilisation and gas diffusion sections.
- Weier, Doran, Power and Walters (1993), *Denitrification and the
  dinitrogen/nitrous oxide ratio as affected by soil water, available carbon,
  and nitrate*, Soil Science Society of America Journal 57(1), 66-72,
  `10.2136/sssaj1993.03615995005700010013x`. Tables 2, 4 and 5.
- Pilegaard (2013), *Processes regulating nitric oxide emissions from soils*,
  Phil. Trans. R. Soc. B 368(1621), 20130126, `10.1098/rstb.2013.0126`.
