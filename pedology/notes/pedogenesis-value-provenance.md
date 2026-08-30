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

## 3. The soil pH block, and the two ends its source sets

Slessarev, Lin, Bingham, Johnson, Dai, Schimel and Chadwick (2016) draw a
spatially random sample of 20,000 subsoil pH measurements from a compilation of
60,291 and find the distribution **bimodal**: a calcite-buffered mode near 8.2
and a gibbsite-buffered mode near 5.1, with an abrupt step between them where
mean annual precipitation begins to exceed potential evapotranspiration.
Predicting only those two values from the sign of the water balance explains 42
per cent of the observed variation, interquartile range 42 to 45 per cent.
Soils in the neutral range are uncommon, because that range is held by primary
mineral dissolution, which is slow.

### Neither end carries the rock, and the form is a crossing between them

In the measurement, gibbsite sets the acid end wherever leaching has run and
calcite sets the alkaline end wherever it has not, whatever the rock was in
either case. `soil_ph` in `build_soil.py` therefore runs the leaching index
between the two buffers rather than taking each parent down a line of its own:

    ph = G + (D - G) * exp(-leaching_slope * L / u),   L = ln(1 + q/q_ref)

with `G` the `gibbsite_buffer_ph` entry, `D` the buffer the cell's mixture
reaches at zero export, and `u` the parent's base-cation supply as a fraction of
a calcite-saturated soil's, declared in
`ph.base_cation_supply_by_category`. At `L = 0` this is `D`, which is the
calcite buffer for every rock that supplies calcium and therefore for every rock
class but the evaporite; as `L` grows every parent goes to `G`.
A line cannot express either end, because it has no asymptote and holds two
parents exactly their initial spacing apart at every slope, and an asymptote at
one end only cannot express the other.

**The fresh parent's own pH is not a state this model can be in.** That is what
the dry end turns on. The mechanism behind the alkaline mode is that calcium
released by weathering is not exported, so pedogenic calcite accumulates until
the solution saturates and calcite takes the buffering over; it needs a calcium
supply and zero export, and every rock class here supplies some calcium. The
only argument for a rock contrast at zero export is that one parent saturates
sooner than another, and Orogen has no time axis, so what a gridcell here
represents is a steady state and a claim about which soil gets there first is
one this project cannot make. `docs/src/reference/no-time-axis.md` is the rule.
The same objection retires the fresh parent as the `L = 0` endmember of the acid
relaxation, which is what it was before the dry buffer was added.

The form keeps the lithology term rather than deleting it, and that is
deliberate: at finite `L` a carbonate parent still sits above a felsic one,
which is what the paper's own wettest-quartile carbonate deviation is a
statement about. A model that collapsed either end to a single value with no
supply term would reproduce the paper's step model and destroy its second
finding, which is measurable: setting every parent to the alkaline value drives
the wettest-quartile carbonate deviation to 0.62 against the paper's 2.6.

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

`calcite_buffer_ph` takes the same atmospheric value and the same bracket, and
is a separate key because it is a separate role. The parent entry is the pH of
a solution on carbonate rock; the buffer is where a gridcell that exports
nothing ends up whatever its rock was. They coincide because the buffering
phase is the same one. Both are `derived`, and `carbonate_ph.resolve` refuses a
number written into either.

The same two atmospheric values bracket every silicate soil solution, from both
sides and for a stated reason. A soil solution supplied with base cations sits
above water in equilibrium with the atmosphere and nothing else, and below the
point where calcite precipitates and takes the buffering over. That is a pH
bound, so it is applied to a pH: bicarbonate carries the alkalinity, a supply
fraction `u` puts the fresh solution at `C + log10(u)`, and
`carbonate_ph.resolve` refuses a run whose declared supplies put it outside the
bracket. A change of `pCO2_bar` moves both ends without moving any declared
supply, which is the failure that would otherwise be silent.

Because the derivation runs on this world's pCO2 and not on a number read off
Earth, it travels to another planet, which is the reason for deriving it rather
than adopting 8.2.

### The remaining two pH entries

`gibbsite_buffer_ph` is 5.1, the mean pH of the spatially resampled NCSS
profiles with non-zero exchangeable Al. It is IMPLICIT-EARTH and not tuned: a
measured population mean is not the residual of a fit, and the derivation is in
the Methods rather than somewhere a reader cannot reach. What makes it Earth's
is the input to that derivation. The paper's eq. (7), `pH = 4.96 + 0.32
log10(CaX/AlX)`, puts 5.1 at an exchange ratio of 2.7, and that ratio is a
property of Earth's lithology, weathering and biological cycling; the Gapon
exchange chemistry travels to another planet, the population it was averaged
over does not. The repair is a cation-exchange model this project does not
have, so the bracket stands in for it: `CaX/AlX` over two decades, 0.1 to 10,
gives 4.64 to 5.28. The two decades are declared -- the paper publishes the
fit, not the range of the ratio in the field.

It is not reachable the way the alkaline end is. `carbonate_ph.py` solves an
equilibrium in `pCO2`; gibbsite solubility carries no CO2 term, so there is no
sentinel to write and no function to fill it.

`leaching_slope` is bracketed by one requirement evaluated at three points of
this world's wettest quartile: that the lithology deviation Slessarev measures
still exists in the model. The paper reports profiles in the wettest quartile
2.6 times more likely to exceed pH 6.5 where carbonate bedrock is present, so a
slope steep enough to put every carbonate-bearing cell of that quartile below
6.5 has deleted a deviation the measurement reports. In closed form,

    s = ln((carbonate - G) / (6.5 - G)) / L

which at the quartile's wettest cell gives 0.35, at its median 0.90, and at its
dry edge 1.33. The declared value is the median: the deviation survives over at
least half the quartile rather than only at its driest fringe. WHICH point of
the quartile is declared; the requirement is not.

A carbonate parent's supply fraction is 1 by construction, so that closed form
is the same expression it was under the single-buffer model and the three
readings do not move when the dry buffer is added.

**All three ends are readings of an UPPER bound, and there is no sourced lower
one.** The bound that would supply it is the measured bimodality, neutral-range
soils uncommon relative to the buffered ranges, and with two buffers in the
model that finding STOPS BEING A LOWER BOUND rather than failing. Stated as the
buffered windows holding at least the land-area density of the neutral window 6
to 7, it holds at every slope down to the smallest tried, because a shallow
slope leaves the whole land on the calcite buffer, which is one of the two
modes. The pass set is not even contiguous: it fails only in a narrow band of
slopes where the land is split between the modes and enough of it is in transit.
A quantity satisfied in the degenerate limit is not a bound.

Re-run against the sourced supply set, the derivation returns the same three
numbers and the same absence. The upper bound is unchanged because a carbonate
parent's supply is 1 by construction and its expression carries no other
lithology. The lower bound is still not there: the density test now passes at
every slope from 0.01 to 4.00 without a gap, where on the declared set it failed
in a narrow band, because a set whose silicates supply four to fourteen per cent
of a calcite-saturated soil's base cations puts the whole silicate land on one
buffer or the other at any slope. The falsification condition, `s_lo > s_hi`,
does not fire, and it does not fire for the same reason as before: there is no
`s_lo` to compare.

### The supply set, and the requirement it has to meet

Both paper findings close on one cell's supply fraction, in an inequality the
slope and the leaching index cancel out of: a carbonate parent has `u = 1` and
stays above 6.5 while `f = exp(-s L)` is at least `(6.5 - G)/(C - G)`, and the
modal silicate parent reaches a gibbsite mode one pH unit wide while `f^(1/u)`
is at most `0.5/(C - G)`. Both hold only if

    u <= ln((6.5 - G)/(C - G)) / ln(0.5/(C - G))

which at the declared buffers is 0.4319. The land-area modal silicate category
is `sedimentary_clastic`.

**The set that meets it is sourced, and the two sources are independent.**
Meybeck (1987) Table 2C gives major-ion concentrations for ten rock types in
monolithologic drainage basins under one temperate stream model; at equal runoff
a concentration ratio is a flux ratio, so the bicarbonate column relative to
sedimentary carbonate rocks is a supply fraction directly. GEM-CO2, Amiotte
Suchet and Probst (1995) as `vendor/cgenie/genie-rokgem` implements its 2003
global form, gives the slope of atmospheric and soil CO2 consumption against
runoff for 232 French monolithologic basins; converting it to an alkalinity
yield by the paper's own stoichiometry, and dividing by the carbonate row,
gives the same quantity from a different measurement. The two agree to better
than a factor of two on every category, and the requirement turns on a factor of
2.32, so they decide it.

| category | value | bracket | rejected pH reading |
| --- | --- | --- | --- |
| `igneous_mafic` | 0.1330 | [0.1330, 0.1565] | 0.8162 |
| `igneous_felsic` | 0.0401 | [0.0310, 0.0725] | 0.4244 |
| `metamorphic` | 0.0423 | [0.0310, 0.0638] | 0.4897 |
| `sedimentary_clastic` | 0.1355 | [0.0391, 0.3052] | 0.5550 |

The requirement is met at every value and at every bracket end, where the pH
encoding it replaced missed it on three of the four categories.

**The order moved, and both sources move it the same way.** The file used to
assert mafic above clastic above metamorphic above felsic. The sourced order is
clastic above mafic above metamorphic above felsic, because Earth's shales carry
carbonate cement and a bicarbonate column sees it. The clastic-over-mafic margin
is two per cent on Meybeck and is not resolvable there alone; GEM-CO2 puts it at
1.9x and settles the direction.

**A supply is not a pH, which is what the encoding got wrong.** The form uses
the parent only through `u`, and `u` sits in the denominator of an exponent
whose numerator is a leaching index, so it is a ratio of a supply flux to an
export flux. A pH is the logarithm of an activity and is proportional to no
flux, so `u = (parent - G)/(C - G)` is an affine map between quantities a
logarithm separates. It is also not a wide enough map: the silicate bracket's
floor is 5.586, so no pH inside it encodes a supply below 0.159, and three of
the four sourced supplies are below that. Reading the sources the other way
instead -- solving Slessarev eq. (6) at each lithology's measured alkalinity to
get a fresh-solution pH -- returns 6.77 to 7.30 and supplies of 0.54 to 0.72,
which fails the requirement on all four. That route is rejected on mechanism and
not on its score: it sizes a fresh parent pH, which this block says outright is
not a state the model can be in.

**What the set is worth on the modelled field**, measured on the offline
reconstruction that reproduces the shipped soil map exactly at its three
decimals, over 1,639 land cells, area-weighted, at the declared slope 0.90.
Land-mean pH falls from 7.261 to 6.717, against Slessarev's own two-value model
at 5.739; 718 cells move by more than 0.5 pH and 465 by more than 1.0. The
paper's bimodality is what moves most: the neutral 6-7 share falls from 0.255 to
0.089 and the gibbsite window rises from 0.065 to 0.399. The paper's carbonate
deviation cannot be read at the declared slope on either set -- the declared set
gives 0.90 on 6 and 14 cells, bootstrapping to [0.35, 1.88], and the sourced set
puts no non-carbonate cell of the wettest quartile above 6.5 at all, so the
ratio is undefined in every resample. What the sweep shows, as a diagnostic and
not as a bound, is that the declared set returned 0.8 to 1.0 at every slope,
which is no lithology signal at all, while the sourced set returns 12 to 18
against the paper's 2.6 in the middle of the slope bracket.

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
