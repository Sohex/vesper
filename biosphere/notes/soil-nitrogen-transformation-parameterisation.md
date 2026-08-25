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
against the source and against the papers, and every bound is arithmetic on
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

## This file is stock LPJ-GUESS 4.1.1, and what this project changed in it

`vendor/lpj-guess/modules/ntransform.cpp` arrived byte-identical to
`guess_4.1/modules/ntransform.cpp` of the 4.1.1 release apart from a stripped
MPL-2.0 licence header, and it is on by default: `data/ins/global.ins` imports
`global_soiln.ins`, which sets `ifntransform 1`, and `parameters.cpp` refuses an
instruction file that does not declare it. So every change below is a change to
default-on behaviour in a widely used community model. The CNP fork never
touched this file; this project did.

That is the reason for the review, and for the standard it was held to. "It is
mainline" is not proof the code is right, and "a paper disagrees" is not licence
to change it quietly. Both are evidence, and a divergence that survives has to
be declared where the next reader will hit it.

**Nothing below is execution-verified.** LPJ-GUESS does not build on this tree,
so every verdict is arithmetic on the declared forms and on the papers. Each
`worth` is a change to a COEFFICIENT, before the Michaelis-Menten terms and the
`min()` against the pool that stand between a coefficient and a flux.

The five `world-i2ch` re-examined as a set:

| what changed | mainline 4.1.1 | now | what settles it | verdict |
| --- | --- | --- | --- | --- |
| `nh3_max`, the ammonia multiplier's extra factor | 0.001 above pH 6, 0.00001 at or below | gone | arithmetic: `exp(2*(pH-10))` alone is 3.4e-4 at pH 6 and 1.0 at pH 10, which is the ratio the paper's text states, and mainline's own comment calls nh3_max the same quantity | KEEP AS DECLARED DIVERGENCE |
| the nitrification gas share | 0.33, read off `f_denitri_gas_max` | 0.022, read off `f_nitri_gas_max` | arithmetic: table 11 brackets RNON + RN2ON at 0.2 to 4.2 per cent, and the code applies the constant as a share and not as a ceiling | KEEP AS DECLARED DIVERGENCE |
| which constant each gas step reads | crossed against their own names | each reads the one its name says | four declarations against two uses, and the values corroborate the names | KEEP AS DECLARED DIVERGENCE |
| the denitrification temperature clamp | `min(1, .)` on table 9 eqn 1 | gone | arithmetic: unclamped the response has the Q10 of 2.04 the paper claims, clamped it has 1.68 | KEEP AS DECLARED DIVERGENCE |
| the denitrification gas partition | a branch at 0.7 WFPS built from three unstated functions | table 9 eqns 5, 6 and 7 | three citation failures against papers now held, plus a hard zero Weier refutes directly | KEEP AS DECLARED DIVERGENCE |

Six more divergences predate that set and are registered on the same terms:
`soil_ph_fallback` (mainline's Dawson regression, which no input ever reaches
because nothing assigns `aprec_lastyear`, replaced by a refusal),
`substrate_partition_argument` and `nitrification_argument` (both curves moved
from a fraction of available capacity to `Soil::wfps(0)`, neither curve
changed), `nitrification_wet_limb` (removed on Greaves and Carter),
`denitrification_gas_constant` (the other half of the crossed pair) and
`mass_balance_check` (mainline's identity, enforced in a Release build rather
than compiled out with the `assert`). `biosphere/config/ntransform.yaml` under
`mainline_divergences` carries all eleven with their arguments, and
`ntransform_gate.py` checks three things about each: that mainline's line is
recorded in the source, that it is not in the source once comments are stripped,
and that the changed line is. So a divergence can become neither a silent fork
nor a silent revert, and reverting by deleting both fails too.

### Why none of them was reverted, and what would have made one so

Two of the five turn on arithmetic that cannot be read another way. `nh3_max` is
0.001 above pH 6 and 0.00001 below, keyed on the same pH 6 the paper's text uses,
and mainline's own comment calls it the "maximum conversion ratio from NH4_mass
to NH3 gas", which is what `fpH` is; the reading that would save it, that it is a
mass-transfer coefficient rather than an equilibrium ratio, fails on the
factor-100 step at pH 6, which a transfer coefficient has no reason to have. The
nitrification gas share is a straight multiplier with no `min()` between it and
the flux, so mainline's 0.33 is what leaves as gas on every day of every
gridcell, against a paper bracket of 0.2 to 4.2 per cent.

The temperature clamp looked like typography and is not. Table 9 eqn 1 prints no
`min{1, .}` where tables 5 and 10 print one, which by itself is weak. What
refuses the clamp is that the paper states the response has a Q10 of about 2, and
only the unclamped function delivers it: 0.594 at 15 C and 1.211 at 25 C is 2.04,
while the clamped pair gives 1.68. The clamp contradicts a number the paper
states about the function it is applied to.

The crossed reads are the weakest of the five and the smallest. Their whole
behavioural content, once the gas share is settled, is the denitrification
NO2-to-gas ceiling moving from 0.25 to 0.33 -- a factor of 1.32 on a constant
that is unsourced either way. There was a cheaper option: leave the reads
crossed and move the value on `f_denitri_gas_max` instead, which is the same
model with one fewer source line changed. It was refused because it puts a
constant named for denitrification on a nitrification quantity and encodes the
crossing in the numbers, where uncrossing puts the two steps of the reduction
sequence on one ceiling, 0.33 and 0.33, as their names and four declarations
say.

The gas partition is the one that costs something. Removing mainline's branch
fixes a hard zero N2 below 0.7 WFPS that Weier tables 4 and 5 refute directly, an
N2O:NO regression attributed to a paper that never measured NO, a logistic in
temperature from a single-temperature study, and a hard zero NO above the branch.
But at 0.90 WFPS mainline's N2 share of 0.846 is inside Weier's interquartile
range and the replacement's 0.978 is above all of it. One point moves the wrong
way while the rest of the domain is repaired, and no choice of table 11 constant
reconciles the two papers: RN2ODN would have to be 25 to 45 per cent to reach
Weier's medians, against a stated 0.2 to 4.7. That cost is not hidden in the
verdict; it is `form:denitrification_n2_share`, verdict `outside`, and it is
now a DECLARED MODEL BOUNDARY rather than a residual. `--strict` no longer
refuses on it and the gate names it under its own heading instead, which is
argued below.

What would have produced a REVERT is a divergence resting on a reading of a paper
that a second reading could undo, with no arithmetic behind it and a real cost in
behaviour. None of the five is that. What a revert would have looked like is the
null result here, not a defeat: the entry leaves `mainline_divergences`, the
record leaves the source, and the gate stops asking about it.

### The comparison arm, specified and not run

The port already exists, so a stock-4.1.1 arm is a second compile and a second
run and no new code: the same build, the same climatology, the same instruction
file except for `f_nitri_gas_max`, and `ntransform.cpp` at the release in the
stock arm. The quantities to compare are `NET_NITRIF`, `NET_DENITRIF`,
`NO_SOIL`, `N2O_SOIL`, `N2_SOIL`, `NH3_SOIL` and the mineral nitrogen the
simulated plants take up.

What it would settle is bounded, and worth being exact about. It cannot say which
divergence is right, because both arms are Earth calibrations and the
disagreements above are between Earth calibrations. What it can do is bound what
they are WORTH, which none of the by-hand numbers in this note can: each of those
is a change to a coefficient, and the Michaelis-Menten terms and the `min()`
against the pool stand between a coefficient and a flux. The five also do not
point the same way on mineral nitrogen -- the gas share raises it, the crossed
read, the pH correction and the temperature clamp lower it, by amounts that
depend on a gridcell's water, pH and temperature -- so which dominates is not
answerable by hand at all. `world-9f1v` owns the arm. It cannot be run until
LPJ-GUESS builds on this tree, which needs a baseline climatology that does not
exist.

## The register: every constant against its calibration

Every source paper is now held and read. Xu-Ri and Prentice (2008) is the
source of every table the operator cites; Weier et al. (1993) is the source of
its denitrification moisture response and of its dinitrogen to nitrous oxide
ratio; Linn and Doran (1984) is where the water-filled pore space shape the
nitrification water response is drawn in comes from, and Greaves and Carter
(1920) is the measurement behind that shape. Each row below gives what the
paper states, what the code has, and whether they agree.

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
| the nitrification water response, argument | `nitrification` | -- | Linn and Doran fig. 1 draws the shape against WFPS, and p.1268 argues WFPS over water-holding capacity | agrees; it read a fraction of available capacity and now reads `Soil::wfps(0)` |
| the nitrification water response, optimum | `nitrification` | 0.6 WFPS | Greaves and Carter table 5: 20 of 22 soils form most nitrate at 50 or 60 per cent of water-holding capacity | agrees, at the top edge of the measured range |
| the nitrification water response, wet limb | `nitrification` | removed | Greaves and Carter measure 9.6 per cent of the optimum where mainline's limb is zero, and Linn and Doran fig. 1 attributes the decline to aeration, which `substrate_partition` applies | REMOVED as a double count with a contradicted zero; mainline recorded beside it |
| the nitrification water response, dry limb, its constants | `nitrification` | slope 5*log(3) per unit WFPS | Greaves and Carter table 5 supports the direction over 22 soils; no source states the slope, and the two constants are one number | unsourced |
| the nitrification water response, dry limb, what it produces | `nitrification` | -- | Stark and Firestone table on p.220 and fig. 1, placed on this model's axis by Cosby (1984) | agrees, inside the 0.333 to 0.430 their two treatments bracket |
| the NO share of nitrification gas | `nitrification` | -- | Xu-Ri table 11 RNON against RN2ON brackets it at 0.33 to 0.98 | agrees, 0.50 to 0.69 over the reachable WFPS |
| the denitrification moisture response | `denitrification` | exponent 13.036 | Weier table 2, at 60, 75 and 90% WFPS | agrees, inside the 8.09 to 14.06 that table brackets |
| the NO and N2O shares of denitrification gas | `denitrification` | 0.002 and 0.02 | Xu-Ri table 9 eqns 5, 6 and 7 with table 11's RNODN and RN2ODN | agrees, inside the 0.2 to 4.9% their sum brackets |
| the N2 share of denitrification gas | `denitrification` | -- | Weier tables 4 and 5 against Xu-Ri table 11 | OUTSIDE Weier at the high end: 0.978 against medians 0.565 to 0.796. DECLARED MODEL BOUNDARY: the operator runs Xu-Ri's partition and the Weier bracket is declared around it |
| the denitrification moisture threshold | `denitrification` | 0.4 WFPS | neither paper states one | unsourced |
| the Michaelis-Menten divisor | `denitrification` | available-water depth | Xu-Ri table 11 gives kg m-3 of an unstated volume | unsourced |
| `pH_soil` | nowhere | 3.5 to 8.5 declared | -- | a declared instruction parameter that no line of the model reads |

`biosphere/config/ntransform.yaml` is the machine-readable form of this table
and of every response function, and the gate checks it against
`global_soiln.ins` and `modules/ntransform.cpp` on each run. A value here that
has drifted from the source is a failure of the declaration.

## The Earth-calibrated response bracket

Eight of the twenty-four declared entries are what remain undeclared, and
`--strict` refuses on exactly those and on the five Vesper preconditions. Four
carry a bracket the gate re-derives against a single constant of the model --
`f_nitri_max`, `f_nitri_gas_max`, the partition's midpoint and the
denitrification moisture exponent -- of which the last two are among the eight.
Four more carry a bracket in prose, because the quantity is not any one constant
of the model: the NO share of nitrification gas, the NO plus N2O share of
denitrification gas, the N2 share of denitrification gas, which is a
disagreement between two papers rather than a number in the code, and what the
nitrification dry limb produces between two stated water potentials. The rest
carry none, because no paper states the quantity and no central value is fitted
here to stand in for one.

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
high end. What the port had was worse than either -- a branch at 0.7 WFPS built
from three functions neither paper states, producing no N2 at all below the
branch where Weier's medians are 1.3 and 2.9 and where N2 is the majority
product in 61 and 80 per cent of the observations.

THE BRACKET IS THE MODEL BOUNDARY, and the operator stays at the Xu-Ri end.
The reason is that the partition is not separable from the sequence it divides.
Table 9 eqns 1 to 4 produce the reduced NO2 flux that eqns 5, 6 and 7 split, and
eqn 1's temperature response is the same `ftemp` that appears inside eqns 5 and
6, so moving the N2 share to Weier would put one step of a four-step sequence on
a different calibration from the three above it and from the response inside it.
Weier cannot supply the rest of the step either: it is a four-soil incubation at
one temperature, measured by acetylene block, and it never measured NO, so it
has no NO limb to take. The two calibrations can be compared and cannot be mixed
term by term.

What the boundary costs is stated rather than argued away. The operator's 0.978
is above the whole 0.333 to 0.938 union, and it carries no dependence on
water-filled pore space where Weier's medians rise from 0.565 to 0.796 across 60
to 90 per cent. So no N2, N2O or NO number this model reports is a prediction of
this world's denitrification gas partition: it is the DyN scheme's partition, at
the DyN end of a two-calibration Earth disagreement, and the reported N2O:N2
ratio in particular is a floor rather than an estimate. Nothing in either paper
narrows that, and narrowing it needs a measurement on this world's soils.

`boundary` is the declaration kind this is recorded under, in
`biosphere/config/ntransform.yaml`. It is only for an entry the sources cannot
settle: the `verdict` still says what the arithmetic says, `--strict` stops
refusing on that entry alone, and `ntransform_gate.py` prints every boundary
under its own heading on every invocation, so one cannot become invisible. A
`boundary` claimed on an entry the sources do settle, or with no owner, is a
failure of the declaration and the gate has a fixture for each.

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

## The nitrification water response: the wet limb is gone, on a primary source

`nit_act` multiplies the table 8 eqn 1 rate by a curve Xu-Ri does not have. The
trail behind it ends outside this project: `ntransform.cpp` arrives complete in
the upstream fork's bulk import of LPJ-GUESS trunk SVN r8583, so no commit
explains the curve, and no held LPJ-GUESS paper adds it. What the curve looks
like says where it came from anyway. Mainline's response peaks at 0.6 and falls
to zero at 0.8, which is the water-filled pore space shape, and its origin in
the soil literature is Linn and Doran (1984).

**That paper settles the argument.** Its whole objective is to replace per cent
water-holding capacity with per cent WFPS as the index of soil aeration, because
WHC "depends on soil type, and the methods for its determination vary and are
often poorly defined" while WFPS needs only gravimetric water content and bulk
density (p.1268). That is the same objection this note makes to
`get_soil_water_upper` two sections above, argued by the paper the curve
descends from. Fig. 1 draws relative nitrification against per cent WFPS.
`nitrification` read `get_soil_water_upper`, a fraction of available capacity,
which is a third quantity again -- neither the axis the curve is drawn on nor
the index its underlying data were measured in. It now reads `Soil::wfps(0)`,
and every moisture response in the operator reads one quantity.

**It is not a nitrification measurement, and the measurement it replots is
held.** Fig. 1 carries four traces and only one of them is nitrification: the
dashed trace, attributed in the caption to Greaves and Carter (1920) and, by the
same caption, converted from an original expressed in per cent water-holding
capacity. The two solid traces are Linn and Doran's own O2 uptake and CO2
production, which are respiration. The dotted trace is denitrification, after
Nommik (1956). Their own nitrous oxide data cannot stand in for the missing
measurement: they state they "did not determine if N2O production resulted from
microbial nitrification ... or denitrification" (p.1271) and only assume
nitrification dominated because their soils sat below 70 per cent WFPS. Greaves
and Carter (1920) is the primary source, and the rest of this section is read
off it rather than off the replot.

### What Greaves and Carter measured, and on what axis

Twenty-two farm soils from Cache Valley, Utah, spanning a loose sand to an extra
tight clay and a fine sand to an organic loam. For each, 100 gm of soil with 2
gm of dried blood mixed in, in covered tumblers, incubated 21 days at 28 to 30
C with the water content restored to its treatment value every 3 days, and the
nitric nitrogen then determined by reduction and distillation; five or six
determinations per point, averaged. So the measured quantity is **net nitrate
accumulated over 21 days in an amended laboratory incubation**, not a
nitrification rate and not a field rate, and the substrate is the amendment.

The axis is per cent of each soil's **water-holding capacity**, and their
methods section defines it: the Hilgard method as modified by Briggs, in which
soil is settled by jarring into a screen-bottomed cup 5 cm across and 1 cm deep,
stood in water until it takes up its maximum, drained 30 minutes, and weighed
before and after drying. It is therefore a gravimetric water content at
drained saturation, and it runs from 31 to 78 per cent across the 22 soils. The
authors use "per cent of saturation" for it interchangeably, and their zero is
at 100 per cent of it.

That definition is what the conversion to WFPS turns on, and it is why the
conversion cannot be done from this paper. WFPS at x per cent of WHC is x/100
times the WFPS of the drained-saturated state, and that factor needs a bulk
density and a porosity per soil. Greaves and Carter report neither: table 2 is
a particle-size analysis, table 3 is WHC, and tables 8 and 9 are moisture
equivalents. Linn and Doran did not state a factor either, and comparison shows
why -- they used 1. Their fig. 1 trace read at 0.10 to 0.80 gives 0.11, 0.18,
0.32, 0.62, 0.86, 1.00, 0.40, 0.11, and Greaves and Carter's own 22-soil mean,
normalised as their fig. 2 normalises it with "the quantity produced at 60 per
cent taken as 100", is 10.9, 16.9, 30.6, 61.9, 86.0, 100.0, 39.5 and 9.6 at 10
to 80 per cent WHC. Same points, same terminal point, same normalisation: per
cent WHC was read straight onto a per cent WFPS axis, by the paper whose stated
objection to WHC is that it is an ill-defined index.

### The rising limb is supported and the falling limb is refused

**The rise is in every soil.** Twenty of the 22 form most nitrate at 50 or 60
per cent WHC -- nine at 50, eleven at 60 -- and the paper puts the optimum
"somewhere between 50 and 60 per cent of saturation" (p.373). Nineteen of the
22 rise without a reversal from 10 per cent WHC to their own optimum; the three
that do not each dip by less than 2 per cent of their own maximum at one
treatment on the way up. The two exceptions are soil 20, a fine sand of 33 per cent WHC whose
nitrate accumulation is flat and small at every treatment, and soil 22. Under
the identity relabelling that is the only conversion anyone has performed, the
measured optimum is 0.50 to 0.60 WFPS and the model's 0.6 sits at its top edge.

**There is no zero at 0.8 anywhere in it.** At 80 per cent WHC the 22-soil mean
is 9.6 per cent of the optimum and no individual soil is at zero: normalised
each to its own maximum, the values at 80 run 0.024 to 0.69. Table 5 and fig. 2
both stop at 80 and carry no point beyond. The paper's only statements past
that are prose -- nitrification "was very slight at the 90 per cent of the
water-holding capacity", and "All soils ceased to nitrify when saturated",
saturation being their 100 per cent. So the primary source puts the zero one
full decade of its own axis above where mainline puts it, and reports a
non-zero mean at the point mainline zeroes.

**The two anchors cannot come from one conversion.** Mainline's peak at 0.6
requires 100 per cent WHC to map to 1.0 WFPS, which is the identity Linn and
Doran used. Its zero at 0.8 requires 100 per cent WHC to map to 0.8 WFPS, which
would put the optimum at 0.44 to 0.48. No single linear map from the measured
axis gives both. The 0.8 that does appear in Linn and Doran is Nommik's
threshold above which significant DENITRIFICATION loss occurs, a different
process in the other direction.

**And the falling limb is a double count while the rising limb is not.** Fig. 1
is annotated across the top with two regimes, WATER LIMITING below 60 per cent
WFPS and AERATION LIMITING above it, and the abstract names the same split --
"Below 60% WFP, water limits microbial activity, but above 60%, aerobic
microbial activity decreases -- apparently the result of reduced aeration"
(p.1267). Those are two mechanisms and this operator already carries one of
them. `substrate_partition` splits every pool on WFPS and hands `nitrification`
only `NH4_mass_d`, the aerobic part, which is aeration. The rising limb is not
that: water limitation on nitrifier activity is not aeration, the partition does
not represent it, and the partition cannot, because its aerobic share GROWS as
the soil dries. Below 0.6 WFPS `nit_act` is the only dry-end limit on
nitrification anywhere in the operator.

### What changed, and what the change costs

`nit_act` is now `min(1.0, act_dry)`: the rising limb unchanged, held at 1.0
from 0.6 upwards, with the aeration decline above the optimum left to the
partition that already applies it. The bound is now an explicit clamp rather
than a property of two constants crossing.

That is a declared divergence from mainline LPJ-GUESS 4.1.1, not a fork quirk:
`vendor/lpj-guess/modules/ntransform.cpp` is byte-identical to the release apart
from its licence header, `global.ins` imports `global_soiln.ins` and that sets
`ifntransform 1`, so mainline's `max(0.0, 4.0 - 5.0*wcont)` is default-on
behaviour in a widely used model. The mainline line is recorded verbatim beside
the changed one in the operator, `biosphere/config/ntransform.yaml` declares the
divergence, and `ntransform_gate.py` checks both that the record is there and
that the code is not running it. `world-i2ch` owns the divergence set.

The cost is at the wet end, and it is worth stating rather than hiding. With
the double count gone, the whole wet-end decline is the partition's, and the
partition's is shallower than the measurement: its aerobic share at 0.7 and 0.8
WFPS is 0.63 and 0.40 of its 0.6 WFPS value, against Greaves and Carter's 0.395
and 0.096. Mainline's stacked pair went the other way, 0.32 and exactly 0 at
the same two points. The measurement sits between them. That magnitude belongs
to the partition's midpoint and shape, which `world-nga8` owns and which are
declared unsourced two sections above; it is not recovered by applying the same
control twice.

### The dry limb is one number, not two, and it has a comparand

**Two constants, one degree of freedom.** The limb is
`exp(5*log(3)*(w - 0.6))`, which is a slope of `5*log(3)` = 5.493 per unit
water-filled pore space written as two factors. Only the product enters, so any
pair with the same product is the same limb and no measurement can separate 5
from 3. A source for "base 3" and a source for "rate 5" therefore do not exist to
be found, and `world-xmiq` is one unsourced number rather than two.

**Greaves and Carter cannot settle the slope.** Against their 22-soil mean the
limb runs at 0.54 to 0.67 of the measurement over 0.1 to 0.5 WFPS, low by about a
third. Against the scatter of the same table, each soil normalised to its own
maximum, it is inside the measured range at every one of those points -- but that
range is 0.00 to 1.00 at 10 per cent WHC and 0.17 to 1.00 at 40, wide enough to
admit almost any monotone limb, so being inside it neither confirms the slope nor
refuses it. A 22-soil bracket that admits everything is not an instrument, and
saying it passed would be reporting noise.

**Stark and Firestone (1995) can, and this model can reach their axis.** They
measured nitrification in a silt loam against soil WATER POTENTIAL, which is the
axis the whole substrate-versus-dehydration question is posed on: in shaken
slurries with ammonium supplied in excess, so that only cell dehydration acts,
and in moist soil, where substrate diffusion acts too. The slurry rates fit
`k = 15.4*exp(0.58*psi)`, psi in MPa, r2 = 0.958. Relative to their -0.1 MPa
reference, the moist soil falls to 0.51 at -0.5 MPa and 0.21 at -2.7 MPa, and the
slurries to 0.79 and 0.22.

The conversion this project refused to invent for Greaves and Carter is one the
model already performs. The land column property contract builds the simulated
soil by inverting Cosby et al. (1984) eqn 1 at two declared pressures: the plant
limiting pressure of -1.5468 MPa, which is Cosby's `10^4.2` cm of water under
Earth gravity, and the drainage equilibrium `rho_w * g * L` over the declared
drainage length. And
`Soil::wfps` is `(wcont*gawc + gwp)/gwsats`, which is the volumetric water
content over the saturation capacity exactly. So Cosby eqn 1 makes potential a
power law in water-filled pore space, per gridcell, out of the sand and clay the
model reads anyway. No relabelling is performed, nothing is fitted, and the
bridge is a paper this project already holds.

Placed there, on each of the current soil map's 4105 gridcells and anchored where
Stark and Firestone anchor:

| relative to -0.1 MPa | this limb, median over the map | interquartile | moist soil | slurries |
| --- | --- | --- | --- | --- |
| at -0.5 MPa | 0.526 | 0.519 to 0.619 | 0.510 | 0.793 |
| at -1.554 MPa, this model's wilting point | 0.364 | 0.356 to 0.424 | 0.333 | 0.430 |

The two treatments bracket a question rather than a scatter: how much of the
dry-end decline `nit_act` should carry. The slurry end is dehydration alone; the
moist-soil end is dehydration plus substrate diffusion. The operator has no other
term for substrate diffusion, and `substrate_partition`'s aerobic share GROWS as
the soil dries, so the moist-soil end is the one this limb should sit near. It
does, at 0.364 against 0.333, inside the bracket by both ends.

**Three limits on citing it, and the verdict that follows.** It is ONE SOIL, so
the bracket is between two treatments of the same silt loam and not between
soils; their own discussion says the relationship "will undoubtedly change for
different soil types and microbial communities" and that diffusional limitation
should be more severe in coarse-textured soils, which is a spread they did not
measure. Their -2.7 MPa point is drier than this model's wilting point and so
outside the window a gridcell can reach, which is why only the -0.1 to -1.554 MPa
span is compared. And the measurement is a 24 h incubation of a sieved, repacked
soil at 23 C with ammonium supplied, so it is a laboratory potential response and
not a field rate.

So the slope is checkable and passes, and it is still unsourced: one soil can
refuse a slope and cannot state one. `function:nitrification_activity` therefore
stays `unsourced` and `--strict` still refuses on it, `world-xmiq` owning it,
while `form:nitrification_activity_dry_limb` carries what the slope produces and
agrees. Nothing is fitted here to stand in for a source. A curve fitted to a
22-soil mean from a 1920 amended-soil incubation, on an axis that has to be
relabelled to reach WFPS, would be an unsourced curve wearing a citation, and a
slope fitted to one silt loam would be the same thing with a better axis.

**The frame the bridge rests on, and it is now settled.** Cosby's air-entry and
wilting heads are heads in centimetres of water, and a head is a pressure only
through the local gravity. Air entry is a capillary pressure and the wilting
point is a plant pressure, so both are invariant and both of their heads scale
together when gravity changes; the closure reads them only through their ratio,
so the simulated soil's wilting point is the water content at -1.5468 MPa on
this world exactly as on Earth. Field capacity is a drainage equilibrium and is
the one state that moves, which the contract now applies. So the axis this
comparison is placed on is a pressure axis and not a column height, which is
what `world-slpa` asked and what `world-of6n` settled by adopting the contract
read literally.

## What holds without any paper: the operator cannot create nitrogen

Every response function above multiplies a nitrogen pool mass, directly or
through a chain, so each has to lie in [0, 1] over the whole range of soil
states this world produces or the operator makes or destroys nitrogen. The gate
samples each function over its declared domain and each chain as a product of
maxima. All of them hold, and two are worth stating because they hold for a
reason rather than by construction.

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
does. None of them has a default, and no source settles any of them: every one
is Earth work at Earth's surface pressure and oxygen partial pressure, and none
carries a water table. Xu-Ri treats soil water in two sublayers of 0.5 and
1.0 m; Weier incubated 60 g of repacked soil in a sealed jar; Linn and Doran
incubated 100 g subsamples compressed to a chosen bulk density in sealed jars.

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
- Linn and Doran (1984), *Effect of water-filled pore space on carbon dioxide
  and nitrous oxide production in tilled and nontilled soils*, Soil Science
  Society of America Journal 48(6), 1267-1272,
  `10.2136/sssaj1984.03615995004800060013x`. Table 1, fig. 1 and the
  objectives on p.1268. Its nitrification trace is Greaves and Carter (1920)
  replotted, not a measurement of this paper's.
- Greaves and Carter (1920), *Influence of moisture on the bacterial activities
  of the soil*, Soil Science 10(5), 361-387,
  `10.1097/00010694-192011000-00004`. The moisture-holding capacity and method
  of experimentation sections, tables 3 and 5, fig. 2, the nitrification
  results on pp.372-374 and the summary table on p.384. The primary measurement
  behind the nitrification water response, and what settles its wet limb.
- Stark and Firestone (1995), *Mechanisms for soil moisture effects on activity
  of nitrifying bacteria*, Applied and Environmental Microbiology 61(1),
  218-221, `10.1128/aem.61.1.218-221.1995`. Fig. 1, fig. 3 and the declines on
  p.220. The only held source that measures nitrification against a moisture
  axis with a stated physical meaning, and the one that brackets what the dry
  limb produces. One silt loam.
- Cosby, Hornberger, Clapp and Ginn (1984), *A statistical exploration of the
  relationships of soil moisture characteristics to the physical properties of
  soils*, Water Resources Research 20(6), 682-690, `10.1029/WR020i006p00682`.
  Not a nitrogen source: eqn 1 and table 4, which the land column property
  contract builds the simulated soil from and which are therefore the bridge
  between a water potential axis and this operator's water-filled pore space
  one.
- Pilegaard (2013), *Processes regulating nitric oxide emissions from soils*,
  Phil. Trans. R. Soc. B 368(1621), 20130126, `10.1098/rstb.2013.0126`.
