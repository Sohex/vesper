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
| `f_nitri_gas_max` | `global_soiln.ins`, then `nitrification` | 0.022 | Xu-Ri table 11, RNON 0.1-4% mean 2% and RN2ON under 0.1-0.2%, applied by table 8 eqns 3 and 4 | agrees, inside the 0.2 to 4.2% their sum brackets |
| `f_denitri_max` | `global_soiln.ins`, then `denitrification` | 0.33 | Xu-Ri table 9 eqn 3 has no maximum-rate constant, and table 11 lists none | unsourced |
| `f_denitri_gas_max` | `global_soiln.ins`, then `denitrification` | 0.33 | Xu-Ri table 9 eqn 4, the same | unsourced |
| the volatilisation temperature response | `nh3_volatilization` | -- | Xu-Ri table 5 eqn 7 | agrees, identical including the min |
| the ammonia multiplier | `nh3_volatilization` | -- | Xu-Ri table 5 eqns 2, 3, 4 and 6 | agrees; the port's extra `nh3_max` is gone |
| the gas emission multiplier | `n_gas_emission` | -- | Xu-Ri table 10 eqns 1 and 2 | agrees, identical |
| the nitrification temperature response | `nitrification` | -- | Xu-Ri table 8 eqn 2 | agrees; the port's min is a no-op, the function peaks at exactly 1 at 38 C |
| the denitrification temperature response | `denitrification` | -- | Xu-Ri table 9 eqn 1, which carries NO min | agrees; the port's min is gone |
| the aerobic/anaerobic split, midpoint | `substrate_partition` | 0.5 | Xu-Ri table 7 settles the VARIABLE, water-filled pore space, and gives no form | unsourced, bracket 0.50 to 0.66 |
| the aerobic/anaerobic split, shape | `substrate_partition` | 7.5 | the same, and there is no bracket for it either | unsourced |
| the nitrification water response | `nitrification` | -- | Xu-Ri table 8 eqn 1 has no moisture term | unsourced, and a second moisture control |
| the NO share of nitrification gas | `nitrification` | -- | Xu-Ri table 11 RNON against RN2ON brackets it at 0.33 to 0.98 | agrees, 0.50 to 0.69 over the reachable WFPS |
| the denitrification moisture response | `denitrification` | exponent 13.036 | Weier table 2, at 60, 75 and 90% WFPS | agrees, inside the 8.09 to 14.06 that table brackets |
| the NO and N2O shares of denitrification gas | `denitrification` | 0.002 and 0.02 | Xu-Ri table 9 eqns 5, 6 and 7 with table 11's RNODN and RN2ODN | agrees, inside the 0.2 to 4.9% their sum brackets |
| the N2 share of denitrification gas | `denitrification` | -- | Weier tables 4 and 5 against Xu-Ri table 11 | OUTSIDE Weier at the high end: 0.978 against medians 0.565 to 0.796 |
| the denitrification moisture threshold | `denitrification` | 0.4 WFPS | neither paper states one | unsourced |
| the Michaelis-Menten divisor | `denitrification` | available-water depth | Xu-Ri table 11 gives kg m-3 of an unstated volume | unsourced |
| `pH_soil` | nowhere | 3.5 to 8.5 declared | -- | a declared instruction parameter that no line of the model reads |

`biosphere/config/ntransform.yaml` is the machine-readable form of this table
and of every response function, and the gate checks it against
`global_soiln.ins` and `modules/ntransform.cpp` on each run. A value here that
has drifted from the source is a failure of the declaration.

## The Earth-calibrated response bracket

Eight of the twenty declared entries are what remain undeclared, and `--strict`
refuses on exactly those and on the five Vesper preconditions. Six carry a
machine-checked bracket, one of them among the eight; a seventh, the N2 share
of denitrification gas, carries its bracket in prose because the quantity is a
disagreement between two papers rather than a constant in the model. The rest
carry none, because neither paper states the quantity and no central value is
fitted here to stand in for one.

Where a paper states a range and one number has to go into the model, the rule
is the same everywhere in the operator: the paper's stated mean where it states
one, and the top of the stated bound where it does not. That rule sets the
nitrification gas share and both denitrification gas fractions, and it is
stated once here rather than argued at each of them.

**The nitrification rate constant.** Xu-Ri table 11 gives Nmax as 0.1 per day
from Khalil et al. (2004) and states it at 20 C, while table 8 eqn 2's
temperature response is 1 at 38 C, which is where the model applies it. That
response is 0.248 at 20 C, so the operator rates 20 C nitrification at 0.0248
per day against the 0.1 per day Khalil reports. The bracket on the constant is
therefore 0.1 to 0.403 per day, the value as the paper writes it against the
value that would honour the measurement it cites. The inconsistency is Xu-Ri's,
carried faithfully into the port.

**The gas fraction of gross nitrification.** Xu-Ri table 11 gives RNON as 0.1
to 4 per cent with a mean of 2, and RN2ON as under 0.1 to 0.2 per cent, and
table 8 eqns 3 and 4 apply them to the gross nitrification flux, so the bracket
on their sum is 0.2 to 4.2 per cent. The operator runs at 0.022, RNON's stated
mean plus the top of RN2ON's bound; RN2ON can move that total by at most 0.2
percentage points, so the constant is RNON's mean to within a tenth of itself.
It ran at 0.33, outside the bracket by a factor of 7.9 on the top and 15 on the
mean, and it was read under the name of the denitrification gas constant. Of
every kilogram of nitrogen nitrified the operator now returns 978 grams to
`NO3_mass_d` where it returned 670, and routes 22 grams to `NO_SOIL` and
`N2O_SOIL` where it routed 330.

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

**The N2 share of denitrification gas, where two Earth papers disagree.** Weier
tables 4 and 5 give the N2/N2O ratio for all four soils, nine treatments and
five days at 25 C. Converted to an N2 share by the paper's own footnote,
`1 - 1/(1 + N2/N2O)`:

| WFPS | observations | N2 share, median | interquartile | Xu-Ri table 9 eqns 5 to 7 |
| --- | --- | --- | --- | --- |
| 60% | 109 | 0.565 | 0.333 to 0.778 | 0.978 |
| 75% | 129 | 0.744 | 0.545 to 0.870 | 0.978 |
| 90% | 163 | 0.796 | 0.630 to 0.938 | 0.978 |

The union of those interquartile ranges is 0.333 to 0.938, and Xu-Ri's own
tables 9 and 11 sit above all of it, because eqns 5 and 6 make NO and N2O small
fixed fractions of the reduction flux and eqn 7 leaves everything else as N2.
The operator runs at the Xu-Ri end and therefore outside Weier's bracket at the
high end. That is a declared disagreement between two Earth calibrations, and
the reason the operator sits at that end is that it IS DyN: every other equation
in it is Xu-Ri's, and this was the one place the port substituted its own. What
the port had was worse than either -- a branch at 0.7 WFPS built from three
functions neither paper states, producing no N2 at all below the branch where
Weier's medians are 1.3 and 2.9 and where N2 is the majority product in 61 and
80 per cent of the observations. Nothing here settles which Earth number this
world's soil should carry.

**The aerobic and anaerobic split.** Xu-Ri settles the variable and gives no
form, pointing instead at PnET-DNDC's anaerobic balloon, so the split's midpoint
and its shape are unsourced in two different degrees and are registered
separately. The midpoint has a bracket: it is the spread between the two numbers
the operator's own comment asserts at once, Pilegaard's crossover at 0.60 WFPS
and the comment's own "steepest change around 66%", so 0.50 to 0.66 with the
declared 0.5 at the bottom edge. The shape parameter has no source and no
bracket at all, which is why it is its own entry rather than a clause in the
midpoint's.

## What the five corrections are worth, by hand

Every number below is arithmetic on the declared forms, and none of it is
execution-verified: LPJ-GUESS does not build on this tree, so nothing here has
been run. Each is the change to a coefficient, before the Michaelis-Menten terms
and the `min()` against the pool that stand between a coefficient and a flux.

**The gas share of gross nitrification, 0.33 to 0.022.** Of every kilogram of
nitrogen nitrified, 978 grams now return to `NO3_mass_d` where 670 did, and 22
grams leave as NO and N2O where 330 did. So the mineral nitrogen a gridcell gets
from a given gross nitrification rises by a factor of 1.46, and `NET_NITRIF`,
`NO_SOIL` and `N2O_SOIL` fall by 15 in their nitrification part. This is the
largest of the five for what the simulated plants can reach.

**The two crossed reads.** The ceiling on the NO2-to-gas step moves from 0.25 to
0.33, so that coefficient of the anaerobic NO2 pool rises by 1.32. The
NO3-to-NO2 step is untouched: its constant was always read where its name says.

**The doubled pH dependence.** The daily volatilised fraction of NH4 rises by
1e3 on the gridcells above pH 6 and 1e5 on the rest. At pH 6.5 its ceiling is
2.3e-4 per day where it was 2.3e-7; at pH 8.5, 1.2e-2 where it was 1.2e-5. On an
alkaline gridcell that is a percent of the ammonium pool a day, which is a real
loss term where the port's was arithmetically absent.

**The denitrification temperature clamp.** Both denitrification steps scale
linearly in the response the clamp held at 1, so a 25 C soil denitrifies at 1.21
times the clamped rate, a 30 C soil at 1.61, a 35 C soil at 2.07 and a 45 C soil
at 3.15. Warm wet gridcells lose nitrogen faster than they did, and the effect
is monotone in soil temperature, so it is not uniform across the map.

**The gas partition.** Total denitrification does not move; only the split among
the three gases does, and all three leave the soil, so this one does not change
what the simulated plants can reach at all. It changes what is reported. At
25 C the reduced nitrogen now leaves as 97.3 per cent N2, 2.42 per cent N2O and
0.24 per cent NO at every water-filled pore space. Below 0.7 WFPS it used to
leave as 0 per cent N2 and the whole flux as NO and N2O, about 59 per cent NO
and 41 per cent N2O at 0.50 WFPS; at 0.90 WFPS it left as 84.6 per cent N2 and
15.4 per cent N2O with no NO. So `NET_DENITRIF` falls by about 37 on the drier
half of the branch and by about 5.8 at the wet end, and `N2_SOIL` gains what
they lose.

The five do not point the same way on mineral nitrogen. The gas share raises it
by 1.46 for a given gross nitrification; the crossed read, the pH correction and
the temperature clamp all lower it, by amounts that depend on the gridcell's
water, pH and temperature. Which dominates is not answerable by hand, and no
absolute nitrogen number from this operator is available until the model builds
and a baseline climatology exists to drive it.

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
median wet limit of 0.813. The consequence is structural, not marginal: 2032 of
4105 gridcells are above 0.40 WFPS even at wilting point, so on half the map the
0.4 threshold below which the operator denitrifies nothing never closes, and on
the rest it decides whether denitrification happens at all. Neither paper states
that threshold.
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
pH 10.7 it would exceed 1 and be held back only by the `min(NH4_mass, ...)`
clamp beside it. It does not get there: `pedology/config/pedogenesis.yaml` clips
its pH field, and over that range the whole multiplier reaches at most 9.2e-2,
one order of magnitude below 1. The bound is therefore a joint property of the
operator and of the soil map that feeds it, which is why the gate declares the pH
domain as pedology's and not as chemistry's. It used to be four orders below 1,
because the port multiplied a further `nh3_max` on top of `fpH` and so applied
the same pH dependence twice; that factor is gone and the margin with it, which
is what makes the pH clip load-bearing rather than decorative.

**Warm-soil denitrification is bounded by the pool, not by its coefficient.**
Xu-Ri table 9 eqn 1 is an exponential equal to 1 at 22 C and rising past it:
1.21 at 25 C, 2.07 at 35 C, 3.15 at 45 C, and 8.07 at the top of the declared
temperature domain. Table 5 eqn 7 and table 10 eqn 1 both state a `min{1, ...}`;
table 9 eqn 1 does not, and the port added one anyway, which held warm-soil
denitrification at its 22 C rate. It did not need to. Both transformations the
coefficient multiplies are already written as `min(pool, pool * coefficient)`,
so neither can take more nitrogen out of a pool than the pool holds however far
the coefficient rises. That is what `clamped_by` declares on the two
denitrification products, and the gate checks the `min()` is still in the source
rather than checking a product that no longer has to stay below 1.

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
have given 8.8, and so fixed the pH term at `exp(-3)` on every soil. Two further
things stood between that fallback and a usable number even with a live input: it is an Earth calibration, and its argument is
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
