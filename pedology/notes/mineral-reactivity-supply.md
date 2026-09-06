# What pedology can supply the mineral-reactivity contract, and what it cannot

This is a worldbuilding project. Vesper is an invented super-Earth and this
document is about the simulated soil this component derives for it, and about
four quantities the vegetation model's mineral-aware arm asks that soil for.

`biosphere/config/mineral_reactivity.yaml` declares two arms for soil organic
matter protection and phosphorus sorption in the simulated soil. The texture
only arm is what runs. The mineral-aware arm refuses, because two of the four
proxies it needs still carry the `undeclared` sentinel and the contract records
that nothing in this pipeline produces them. This is the pedology-side verdict
on each:
whether this component can derive it from the World Orogen export and the
ExoPlaSim climatology it already reads, and if not, what would settle it.

Two of the four cannot be derived here and stay declared absences. The other two
are derived and emitted, as of 2026-09-05: the soil map carries an `allophane`
column and a `polyvalent` column and both are declared in the contract's
`carried_state`. THE ARM STAYS REFUSED EITHER WAY, and the last section of this
document is why: the proxies are the first of two gates and the second is on the
consuming model's side. Each section below is the evidence for one of them,
because "nothing produces it" and "the measurement exists and indexes on an axis
this project does not have" are different states and only the second says what
to go and get.

| proxy | unit | verdict | the axis that stops it, or the source that settled it |
| --- | --- | --- | --- |
| Fe-Al oxide content | kg/kg of fine earth | not derivable | the total-to-extractable step, and the one held source that measures it is single-lithology and age-indexed |
| allophane content | kg/kg of fine earth | DERIVED and emitted | Parfitt, Russell and Orbell (1983) index whole-soil allophane on rainfall and on a water balance, which is the controlling axis; the content is declared only well past the leaching threshold, because in the transition window that paper's own near-duplicate sites differ by an order of magnitude |
| aggregate capacity | dimensionless | not derivable | nothing in the pipeline resolves soil structure, and no held source transfers from what it does emit |
| polyvalent cation saturation | cmol(+)/kg | DERIVED and emitted | the emitted `cec` column and Solly et al. (2020)'s pH-indexed cation share supply both halves; a lower bound, because only one polyvalent cation is resolved at each end of the pH range |

## 1. Iron and aluminium oxide content

What the arm asks for is the oxalate- and dithionite-extractable Fe and Al: the
pedogenic, short-range-order fraction that carries the reactive surface. It is
not the parent rock's total iron and aluminium. A parent total could be declared
per rock class the way `pedogenesis.yaml` already declares a parent texture per
rock class, with the same kind of source. The step that has no source here is
total to extractable, and it is the whole of the quantity.

The best held source that measures the extractable fraction is Chadwick et al.
(2003): sixteen profiles on Kohala, sequential extraction with acid ammonium
oxalate and with Na-dithionite and Na-citrate, along an arid to humid transect,
against a leaching index defined as the annualised water balance over the
integrated porosity of the top metre. That index is a quantity this pipeline
could compute, since the climatology carries the water balance and the land
column property contract carries the porosity. Two things stop it being a
transfer function for this world.

It is ONE parent material, Hawaiian basalt, and this world's land is a mixture
of the export's rock classes, most of which are not basalt. And it is one
substrate age band: the paper's own discussion is explicit that at 2500 mm of
rainfall a soil younger than 20 ka holds its exchange properties high relative
to an older soil under the same rainfall, and its transect is at 170 ka. Orogen
has no time axis (`docs/src/reference/no-time-axis.md`), and this component's
weathering intensity deliberately folds the time integral into its Earth
reference rather than carrying an age. So the source's own independent variable
is one this project cannot supply, and a field derived from present climate
alone would be asserting a soil history the pipeline does not represent.

The andisol block already carries Dahlgren's andic criterion, Al-ox plus half
Fe-ox at or above 2 per cent, as the definition of andic development. That is a
statement about extractable oxides, and it is a floor on andic ground by
construction. Reading it back out as the proxy would be recovering the criterion
from the classification it produced, which tells a consumer nothing the andic
column does not already say.

## 2. Allophane content

DERIVED AND EMITTED as of 2026-09-05, on the axis that controls it, and the
source is the one that was sought and not reached when this section last said it
could not be. `build_soil.py`'s `andisol_properties` emits `andic`, the AREAL
fraction of a gridcell whose volcanic glass has weathered to allophane, and the
`allophane` column is that times the allophane content WITHIN andic material.

THE AXIS IS LEACHING AND NOT DEVELOPMENT. Parfitt (2009) reviews the controls on
allophane formation in tephra and Andisols and states that "the effect of time is
subordinate to these factors", the factors being the activity of silicic acid in
the soil solution, the availability of Al species, the opportunity for
co-precipitation, leaching regime, soil organic matter and pH. A distribution of
allophane content against a DEVELOPMENT measure is therefore a distribution
against the weakest of the controls, and a transfer built on it would carry the
residual variance of all the others. That is why `pedogenesis.yaml`'s `andisol`
block gates andic development on a precipitation threshold, and what was missing
was the CONTENT above the threshold rather than the threshold.

THE SOURCE. Parfitt, Russell and Orbell (1983), `10.1016/0016-7061(83)90029-0`
-- the paper Dahlgren (2004)'s 1500 mm boundary is taken from at second hand.
Four soils on one rhyolitic ash stratigraphy in the North Island, sited where
the rainfall isohyets cross the tephra isopachs at right angles so that the
stratigraphy is held constant and rainfall is the designed variant. Whole-soil
allophane is measured by acid-oxalate dissolution with the Al held in humus
complexes subtracted by pyrophosphate, so what the tables report is a
concentration in the soil and not in the clay fraction, which removes the clay
conversion this section used to require. Table VIII carries a water balance
beside the rainfall: 400, 550 and 1600 mm of winter leaching and 5, 6 and 10
months of leaching a year at 1200, 1400 and 2600 mm of mean annual rainfall.

WHAT IS DECLARED, measured on 2026-09-05 from that paper's Table III. The value
is 0.22 kg of allophane per kg of fine earth within fully andic material, the
mean over the seven sampled layers of the Mairoa Hydric Dystrandept from 0 to
104 cm: 5, 15, 26, 26, 25, 25 and 35 per cent. Mairoa is the one site in the
sequence unambiguously past the threshold the andisol block gates on -- 2600 mm,
1600 mm of winter leaching, ten months of the year, and no summer deficit in a
normal year -- so it is the site whose regime matches the ground the factor is
applied to. Its 5 per cent topsoil is not an outlier to be trimmed: pyrophosphate
recovers 2.5 of its 3.4 per cent oxalate Al as Al-humus complexes, which is the
soil organic matter control Parfitt (2009) names, operating in a profile.

THE BRACKET IS WIDE BECAUSE THE POPULATIONS DISAGREE, and it is swept rather
than narrowed: 0.02 to 0.35. Its high end is that profile's deepest sampled ash.
Its low end is Parfitt (2009)'s Icelandic field soils, 2 to 22 per cent allophane
at 11 to 42 per cent C (Sigfusson et al. 2008), a basaltic parent under a carbon
load this world's andic ground need not carry; the same review puts Atlantic
volcanic island mineral horizons under 6 per cent. Choosing among them would be
choosing a parent material for this world's arcs that nothing here supports.

WHAT THE SOURCE REFUSES TO LICENSE, and both of the papers that were sought make
the same refusal from different sides. Between roughly 1200 and 1600 mm the
content is not placeable at all. Kereone at 1200 mm of rainfall holds 0.5 to 2
per cent allophane and Tirau at 1270 mm holds 1 to 17 per cent on the same
parent material and the same stratigraphy, and Parfitt et al. read that as "only
a slight change in leaching regime or temperature is sufficient to alter the
composition of these soils". Singleton et al. (1989), `10.1071/SR9890067`, close
it from the other side: three soils 200 m apart at one site, one rainfall of
1201 mm and one vitric rhyolitic alluvium, where drainage class alone takes a
profile from allophane-dominant to no allophane at all with the switch at about
10 g/m3 of silicon in soil solution. So the scatter inside the transition window
is the size of the whole quantity, and the emitted factor is a statement about
ground the threshold has already placed well past it. The `allophanic` logistic
in `andisol_properties` carries cells through that window with a partial andic
fraction, and the content factor is not a claim about any of them individually.

## 3. Aggregate capacity

Occlusion inside aggregates is physical protection and is separate from mineral
association. The soil map carries texture, regolith depth, bulk density, pH, an
organic fraction and an andic fraction. None of those is an aggregate size
distribution or a stability measure, nothing else in the pipeline resolves soil
structure, and no held source supplies a transfer from any of the columns that
do exist to an aggregate protection capacity. This is the one of the four where
the absence is total: there is no measurement to go and index, so what would
settle it is a parameterisation keyed on quantities this component emits.

## 4. Polyvalent cation saturation

DERIVED AND EMITTED as of 2026-09-05. The arm asks for it in cmol(+)/kg, which is
an exchange capacity, so it needs both a capacity and the share of it held by
calcium, iron and aluminium. Both halves exist and the `polyvalent` column is
their product.

THE CAPACITY IS EMITTED, AND ITS LEVEL NO LONGER RESTS ON ONE REGION.
`pedogenesis.yaml`'s `exchange` block and `build_soil.py:exchange_properties`
derive cation exchange capacity additively from the clay and organic-matter
fractions the soil map carries, on coefficients that are themselves linear in the
soil map's pH: Helling, Chesters and Corey (1964),
`10.2136/sssaj1964.03615995002800040020x`, measured both contributions at six
buffered pH values over 60 Wisconsin soils with no component removed, at R2 0.86
to 0.92, and both are linear in pH within the standard errors of all six points.
The LEVEL of that pair is a statement about clay mineralogy, and it is bracketed
across nine soil orders and 37,921 pedons by Manrique, Jones and Dyke (1991),
`10.2136/sssaj1991.03615995005500030026x`, whose clay contribution runs 8.4
cmol(+)/kg of clay for Oxisols to 63 for Vertisols and which concludes that
without a parameter stratifying clay mineralogy, clay content accounts for little
of the variation in CEC. `build_soil.py` sweeps that bracket and reports what the
ANUT-8 bound is worth at each end. The exposure world-n4i0 recorded is closed by
the bracket and not by a better point value: there is no point value, because the
quantity is a mineralogy this component does not resolve.

THE SHARE IS MEASURED AGAINST A VARIABLE THIS COMPONENT EMITS. Solly et al.
(2020) partitioned effective CEC by cation over 1204 Swiss forest profiles within
nine pH classes: at pH at or above 5.5, exchangeable calcium carries 59 to 83 per
cent of it, and below pH 5 exchangeable aluminium carries 21 to 44 per cent.
Calcium and aluminium are both polyvalent, so polyvalent saturation is
substantial at both ends of that range and the quantity is recoverable from the
pH field for a reason rather than by interpolation. Iron contributes under 1 per
cent, which is why the arm's third cation does not need its own field.

IT IS A LOWER BOUND AND IS DECLARED AS ONE. Only one polyvalent cation is
resolved at each end of the pH range and magnesium sits inside the paper's "other
cations" group at both ends, so the emitted share is the resolved cation's alone
and the true polyvalent share is at least that. The direction is DOWN, which is
the opposite of the `cation_share_upper` block beside it and is deliberate: the
consumer would read this as a PROTECTION proxy, and understating protection is
the conservative direction for that, where understating a circulating nutrient
pool is not.

Chadwick et al. (2003) remains the direct measurement of the quantity and remains
unusable as a FIELD, for the three reasons this document gave before: its
effective CEC is non-monotonic in the leaching index, so a cell's value is not
recoverable from a one-sided argument about how wet it is; it is one parent
material and one substrate age band, and the paper states the age dependence in
this very quantity; and it states that the decline is irreversible under natural
conditions, so the state is a function of the climate PATH and a pipeline with no
time axis has no path to evaluate it on. What it IS used for is a check: its
Table 7 pairs pH with base saturation horizon by horizon and puts the transition
where Solly's calcium share puts it, on a different continent and a different
parent material. Two sources agreeing on where a threshold sits is what licenses
the ramp `exchange.base_saturation` declares, and the polyvalent ramp uses the
same two bounds for the same reason.

## What this means for the arm, and what it does not mean

Two of the four stay undeclared, so `mineral_reactivity_gate.py` keeps refusing
the mineral-aware arm and `active_arm` stays `texture_only`. What crosses the
interface is now five columns rather than two: `andic` and `pfixation` as before,
plus `cec`, `allophane` and `polyvalent`, all five declared in the contract's
`carried_state` and read by no equation.

One thing a later reader should not have to rediscover, and it is what stops the
two newly derived proxies being read as progress toward the arm. Closing all four
of these would still not let the arm be written. Cotrufo et al. (2013) and
Lehmann and Kleber (2015) place mineral association, aggregation and
accessibility at the centre of stable soil organic matter formation and name the
associating surfaces, and neither supplies a parameterised transfer from any of
the four to a protection coefficient. The proxies are the first of two gates, not
the only one, and the second gate is on the consuming model's side. Two of four
declared moves nothing on that second gate.
