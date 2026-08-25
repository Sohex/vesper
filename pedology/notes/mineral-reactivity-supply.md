# What pedology can supply the mineral-reactivity contract, and what it cannot

This is a worldbuilding project. Vesper is an invented super-Earth and this
document is about the simulated soil this component derives for it, and about
four quantities the vegetation model's mineral-aware arm asks that soil for.

`biosphere/config/mineral_reactivity.yaml` declares two arms for soil organic
matter protection and phosphorus sorption in the simulated soil. The texture
only arm is what runs. The mineral-aware arm refuses, because four proxies it
needs all carry the `undeclared` sentinel and the contract records that nothing
in this pipeline produces them. This is the pedology-side verdict on each:
whether this component can derive it from the World Orogen export and the
ExoPlaSim climatology it already reads, and if not, what would settle it.

The verdict is that none of the four can be derived here, so all four stay
declared absences and the arm stays refused. Each section below is the evidence
for one of them, because "nothing produces it" and "the measurement exists and
indexes on an axis this project does not have" are different states and only the
second says what to go and get.

| proxy | unit | verdict | the axis that stops it |
| --- | --- | --- | --- |
| Fe-Al oxide content | kg/kg of fine earth | not derivable | the total-to-extractable step, and the one held source that measures it is single-lithology and age-indexed |
| allophane content | kg/kg of fine earth | not derivable | the concentration inside andic material, whose held source reports the extremes and not a distribution |
| aggregate capacity | dimensionless | not derivable | nothing in the pipeline resolves soil structure, and no held source transfers from what it does emit |
| polyvalent cation saturation | cmol(+)/kg | not derivable | measured, but as a non-monotonic function of leaching index AND substrate age, and there is no age |

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

The near miss, and it stays one. `build_soil.py`'s `andisol_properties` emits
`andic`, the AREAL fraction of a gridcell whose volcanic glass has weathered to
allophane. Turning that into a concentration in the fine earth needs two
conversions. The first is arithmetic and this component already has its parts:
an areal share becomes a mass share through the bulk densities, and the andic
bulk density is declared and is already blended into the soil map's own density
column. The second is the allophane concentration WITHIN andic material, and
that is the number nothing here carries.

Parfitt (1990) measures it. Its Table 3 gives the two New Zealand soils it names
as the greatest concentrations in the country: 10 and 22 per cent allophane in
the topsoil horizons, 16 to 38 per cent below them, with up to 60 per cent
measured in deep deposits, and the acid-oxalate estimate resolving half a per
cent. So an upper end is sourced and a central tendency is not, because the
paper reports the maxima rather than a distribution over allophanic soils.

Substituting the maxima would put a factor of several into a quantity a
protection coefficient is linear in. It would also contradict the model beside
it: `pedogenesis.yaml` carries a development axis for andic material, so a
concentration that did not respond to development would say that a barely andic
cell and a strongly developed one hold the same allophane, which the andisol
model itself denies. What would settle this is a distribution of allophane
content over allophanic soils against a development measure, not a review's
extremes.

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

The arm asks for it in cmol(+)/kg, which is an exchange capacity, so it needs
both a capacity and the share of it held by calcium, iron and aluminium. This
component carries pH and texture. No held source supplies a cation exchange
pedotransfer from texture and organic carbon, so that route has no source.

The one held source that measures the quantity directly is again Chadwick et al.
(2003), and here it is measured against exactly the index this pipeline could
compute: effective cation exchange capacity rises from about 20 cmol(+)/kg in
the arid part of the transect to about 50 near 1300 mm of rainfall and falls
below 5 by 1500 mm, with base saturation near 100 per cent below the threshold
and very low above it, the turn falling at a leaching index of about 1.

Three properties of that result stop it becoming a field here, and they are
worth stating because the relation looks usable at first sight. It is
non-monotonic in the driver, so a cell's value is not recoverable from a
one-sided argument about how wet it is. It carries the same single-lithology and
single-age restriction as section 1, and the paper states the age dependence in
this very quantity. And the paper states that the decline is irreversible under
natural conditions, so the state is a function of the climate PATH and not of
the current climate; a pipeline with no time axis has no path to evaluate it on.

## What this means for the arm, and what it does not mean

All four stay undeclared, so `mineral_reactivity_gate.py` keeps refusing the
mineral-aware arm and `active_arm` stays `texture_only`. What crosses the
interface is unchanged: the soil map's `andic` and `pfixation` columns, declared
in the contract's `carried_state` and read by no equation.

One thing a later reader should not have to rediscover. Closing all four of
these would still not let the arm be written. Cotrufo et al. (2013) and Lehmann
and Kleber (2015) place mineral association, aggregation and accessibility at
the centre of stable soil organic matter formation and name the associating
surfaces, and neither supplies a parameterised transfer from any of the four to
a protection coefficient. The proxies are the first of two gates, not the only
one, and the second gate is on the consuming model's side.
