# Where the pedogenesis constants come from, for the three blocks that had no brackets

This is a worldbuilding project. Vesper is an invented super-Earth and this
document is about the simulated soil `pedology/config/pedogenesis.yaml`
computes for it. Nothing here describes Earth except the measurements that
constants are taken from, which are named as such.

The `catena`, `water`, `andisol` and `bedrock_water` blocks of that file each
ship a bracket and an instruction not to tune inside it. The `regolith`,
`texture` and `ph` blocks did not, and those were where the values with no
derivation sat. This note carries the arithmetic behind the brackets those
three now ship, and the two findings that came out of reading the sources for
them.

## 1. The regolith asymptote is a depth-to-bedrock statistic

`regolith.maximum_depth_m` is the depth the modelled weathering front reaches
as erosion vanishes. It is a pure prefactor: every simulated cell's regolith
depth, and therefore the land-mean plant-available water capacity, scales with
it exactly.

Neither of the two papers under `erosion_coefficient_per_relief_m` bears on it,
and that is a property of the quantity rather than a gap in those two. A soil
production function integrates to a STEADY-STATE depth against an erosion rate,
which is what that key's bracket already is; it does not integrate to a maximum
reach. Heimsath's form has no asymptote at all, which is why this file replaced
it, and Portenga and Bierman measure denudation.

What does bear on it is a measured depth-to-bedrock distribution. Shangguan,
Hengl, Mendes de Jesus, Yuan and Dai (2017) compile 130,000 soil profiles and
1.6 million boreholes and define depth to bedrock as the depth from the ground
surface to contact with coherent bedrock. Their Table 1, world row, absolute
depth to bedrock over 1,590,464 observations:

| statistic | value |
| --- | --- |
| minimum | 0 cm |
| median | 670 cm |
| mean | 1,309.3 cm |
| maximum | 312,541 cm |

Their censored series, which is the same quantity truncated at 200 cm because a
soil survey stops there, has a world mean of 97.51 cm over 307,936
observations. Their predicted map has a mean absolute depth of 33.6 m and puts
85 per cent of land above 2 m.

**The measurement is not the modelled quantity, and it misses it in two
directions that the paper itself names.** Both become an end of the bracket.

*From below, the median.* The observed population is overwhelmingly eroding
land, and an eroding profile is thinner than the asymptote it would reach
without erosion. The paper adds a second reason running the same way: a soil
survey logs intact weathered rock as R horizon, so "the definition of R horizon
or hard rock is not strictly equal to bedrock", and a map built from profiles
reads as deeper than what it prints.

*From above, the mean.* Depth to bedrock counts transported fill exactly as
readily as in-situ weathering product, and this model's regolith is in-situ
product only. The mean is 1.95 times the median and the maximum is 3.1 km,
which is the size of that contamination. The observation set is dominated by
North American and European boreholes, where glacial drift and alluvium are
most of what the drill passes through.

So the bracket is 6.70 to 13.09 m and the declared value is the low end,
because every correction that separates the observation from the modelled
quantity runs upward from the median and downward from the mean, and the
depositional contamination is the larger of the two.

**This key and `erosion_coefficient_per_relief_m` are jointly constrained and
neither key's own bracket says so.** The asymptote is a prefactor and the
coefficient carries the level, and the land-mean depth the pair produces has to
land inside `regolith_depth_bracket_m`, which brackets the level. Moving one
alone can put the level outside a bracket that neither key mentions. A sweep
runs the two together.

## 2. Dry denudation is 0.85 of the reference rate, not 0.15

`regolith.dry_erosion_baseline` floors the moisture term, which is
`runoff / reference_runoff + dry_erosion_baseline`. So the key IS the ratio of
denudation at zero runoff to denudation at Earth's reference runoff, as
`b / (1 + b)`, and that ratio is something Portenga and Bierman (2011) report.

Their arid drainage basins erode at 100 +/- 17.3 m/Myr, the slowest of their
climate zones, against a global drainage-basin mean of 218 +/- 35 m/Myr. The
ratio is 0.459 and inverts to `b = 0.85`. Carrying their two stated
uncertainties through gives the bracket 0.49 to 1.79.

Two approximations set that width and they pull opposite ways. Arid basins
carry some runoff, so the true zero-runoff ratio is below theirs and the true
`b` is lower. Their global mean is drawn from steep, tectonically active,
accessible terrain, so it is high for a reference condition, which makes the
ratio too low and the true `b` higher.

The bracket reaching 1.79 -- a moisture term nearly flat in runoff -- is also
what the same paper supports directly: basin denudation has **no significant
bivariate correlation with mean annual precipitation at all**, mean basin slope
being the regressor. A weak moisture dependence is the reading that compilation
licenses. The 0.15 this key used to carry is a suppression of dry denudation to
13 per cent of the reference rate and nothing supported it.

## 3. The soil pH block, and the form its own source refutes

Slessarev, Lin, Bingham, Johnson, Dai, Schimel and Chadwick (2016) draw a
spatially random sample of 20,000 subsoil pH measurements from a compilation of
60,291 and find the distribution **bimodal**: a calcite-buffered mode near 8.2
and a gibbsite-buffered mode near 5.1, with an abrupt step between them where
mean annual precipitation begins to exceed potential evapotranspiration.
Predicting only those two values from the sign of the water balance explains 42
per cent of the observed variation, interquartile range 42 to 45 per cent.
Soils in the neutral range are uncommon, because that range is held by primary
mineral dissolution, which is slow.

### The finding: the acid end is parent-independent and this form makes it parent-dependent

In the measurement, gibbsite sets the acid end wherever leaching has run,
whatever the rock was. `soil_ph` in `build_soil.py` instead starts every parent
at its own value and takes a single `leaching_slope` down from it, so its wet
cells spread across a range where the measurement clusters on one buffer. No
single slope lands them all on 5.1.

That is a property of the FORM, not of any value in the block, and it is why
the values are bracketed rather than fitted: they cannot be sized against this
observation until the form changes. Recorded here rather than acted on, because
replacing a parameterisation is a larger change than the one this note is the
evidence for, and because a form with a parent-independent acid end is a
different model of leaching and not a retuning of this one.

### The alkaline end is derivable, and is derived at run time

`pedology/scripts/carbonate_ph.py` implements the equilibrium. `build_soil.py`
calls it on every run, `pedogenesis.yaml` states no pH for anything the
equilibrium decides, and a number written back into those keys is refused rather
than used. That is the whole disposition: the alkaline end is not a value this
project holds, it is a function of `config/planet.yaml`'s `pCO2_bar`, and the
only place it appears as a number is on `soil_report.json` under
`carbonate_system_ph`, beside the soil map it was used to build.

The equation and its constants are in that module's docstring. Two things about
it belong here rather than there, because both are findings about the source.

**The equation as PRINTED in the paper is not the equation.** Its last two terms
are set as `- H Kw K1 KH p - K1 K2 KH p`: the hydroxide and bicarbonate terms
have been run together into a product where the charge balance has a sum, and
the two charges a carbonate ion carries have been dropped. Writing the charge
balance out recovers both. The first loss is worth half a pH unit and the
implementation's test rejects it outright; the second is worth 0.0012 pH at the
paper's own pressure, which no test at the paper's published precision can see,
so the factor of two rests on the arithmetic and is recorded here as resting on
it.

**The paper's second figure is not a second test point.** Slessarev state pH 8.2
at their laboratory partial pressure of 3.45e-4 atm, and the implementation
returns 8.236 there, which is the check: both halves of it are the paper's. They
also say the expected pH is 8.3 for measurements made before 1977, and state no
pressure for that at all. An earlier version of this note supplied 3.30e-4 atm
for it, got 8.248, and recorded that both figures "round to what the paper
prints". They do not: 8.248 rounds to 8.2. What the remark does constrain, at
the one decimal it is printed to, is that the solved pH crosses 8.25 somewhere
in the ambient CO2 of the decade before 1977, and the implementation puts that
crossing at 3.281e-4 atm, 328 ppmv. That is checkable, it can fail, and it is
what the module checks. A test whose input is chosen after the answer is not a
test, and the 3.30e-4 was chosen that way.

Evaluated at Vesper's declared pCO2 rather than Earth's, the same equation gives

| condition | pH |
| --- | --- |
| water and atmospheric CO2 alone, no alkalinity | 5.586 |
| calcite equilibrium, open to the atmosphere | 8.163 |
| calcite equilibrium at ten times atmospheric, soil air | 7.498 |
| calcite equilibrium at a hundred times atmospheric, soil air | 6.832 |

*Measured on 2026-08-27, at `pCO2_bar` 0.00045. They are here as a record of what
the equilibrium returned on that day and not as a value anything reads; the
generator is the only statement of them that a consumer sees.*

The `carbonate` parent entry takes the atmospheric value, because that is the
condition a laboratory pH is measured under and the one Slessarev's 8.2 is
stated at, and this number is compared against measured pH. The soil-air arm is
the low end of its bracket: a field pH on carbonate reads below the open-air
one, and the ten-to-a-hundred-fold enrichment that says how far below is a
declared bracket on soil respiration rather than a measurement on this world.

The same two atmospheric values bracket every silicate parent, from both sides
and for a stated reason. A soil solution supplied with base cations sits above
water in equilibrium with the atmosphere and nothing else, and below the point
where calcite precipitates and takes the buffering over. What remains declared
inside that bracket is the ORDER and the SPACING of the four silicate entries;
no source here sizes that contrast. The bracket is now enforced rather than
described: `carbonate_ph.resolve` refuses a run whose declared silicate parents
have fallen outside it, which is the failure a change of `pCO2_bar` would
otherwise cause silently.

Because the derivation runs on this world's pCO2 and not on a number read off
Earth, it travels to another planet, which is the reason for deriving it rather
than adopting 8.2.

### The remaining two pH entries

`leaching_slope` is bracketed by one requirement evaluated at two runoffs: that
the term can carry the modal clastic parent down to the gibbsite buffer of 5.1,
a drop of 1.7. At Earth's reference runoff the leaching index is `ln 2` and the
slope that does it is 2.45; at five times that runoff the index is `ln 6` and
the slope is 0.95. The multiple of five is declared, not measured, and it is
what sets the bracket's width. The declared value is the low end, and the
reason is the clip rather than a preference: a steeper slope drives the felsic
and metamorphic parents onto `ph.minimum` across the wet fraction of land,
which is the railing the regolith block was rewritten to remove.

`endorheic_alkalinity_bonus` has to carry a calcite-buffered soil into the range
a closed basin reaches. Helvaci (2019) lists lake water at pH 8.5 to 11 among
the conditions the Turkish borate deposits form under, which is the one measured
closed-basin range this project holds. Helvaci's range minus the derived
carbonate pH is the bracket, so the bracket is derived too and is emitted
alongside it rather than written into the config. It is lake water in one
depositional setting and not a soil, so it bounds the entry and does not set
it.

## 4. The two texture entries, and one near-degeneracy

`clay_conversion` is the rate constant in `1 - exp(-clay_conversion * W)`. Its
reciprocal is the weathering intensity at which the conversion is `1 - 1/e`
complete, in the units W is normalised in, so the bracket is stated on that
intensity: no faster than complete by Earth's land mean, `W = 1`, and no slower
than still running at the wettest of the sixteen type localities, `W = 3.41`,
which is the wettest soil this model has been compared against. That is 0.29 to
1.00, and the declared value puts the e-folding at `W = 1.82`.

**It is near-degenerate with `clay_yield` at low weathering intensity.** The two
multiply, and where `clay_conversion * W` is small the exponential is in its
linear regime, so only the product is identified. That is the same defect the
regolith block removed by collapsing `erosion_weight` and
`erosion_reference_relief_m` into one key, and it is NOT removable the same
way: the two separate at high W, which is where the sixteen localities are and
where this world is not. They have to be swept together, and a Vesper clay
figure rests on the product rather than on either key.

`sand_to_silt_loss_ratio` is not derivable from mechanism, because the two
mechanisms point opposite ways.

*Specific surface area* goes as one over grain diameter, so silt carries more
reactive area per unit mass than sand and dissolves faster. Over the USDA class
limits the geometric mean diameters are 0.316 mm for sand and 0.010 mm for
silt, a factor of 31.6, which is a ratio of 0.032.

*The size cascade* runs the other way. A sand grain becomes silt before it
becomes clay, so sand leaves its own class faster than dissolution alone
implies while silt is replenished out of it. In the limit where silt sits in
quasi-steady state, all of the clay comes from silt and all of silt's loss is
made up from sand, so the ratio is unbounded above.

The bracket is the span those two allow, 0.032 to 31.6, and the declared value
is the one point inside it that carries a statement: loss in proportion to
abundance, no size preference. The comment beside the old value read "sand goes
first: coarse primary grains have the most surface area to lose", which has the
surface-area argument backwards, and is why the number did not follow from the
sentence printed next to it.

The entry redistributes between sand and silt at fixed clay and moves
plant-available water capacity by under one per cent, because clay's volumetric
capacity sits between sand's and silt's. A bracket that wide still moving the
capacity that little is the result, not a weakness of the bracket.

## What has not been measured

The soil step cannot run on the active build: it needs a baseline climatology
and none exists for `canonical-10m-base` yet. So the effect of these changes on
the modelled soil is stated from the law and from the sensitivities the
tuned-values audit measured, and not from a fresh soil map. The two that move
the field are the asymptote, which is exactly proportional, and the dry erosion
baseline, which raises erosion most where runoff is least. They act against each
other on the land mean: the asymptote deepens every cell and the erosion floor
thins the dry ones.
