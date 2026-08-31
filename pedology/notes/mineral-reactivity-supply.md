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

Three of the four cannot be derived here and stay declared absences. The fourth
can, since this component began emitting a cation exchange complex, and it is
wiring rather than evidence that keeps it undeclared. THE ARM STAYS REFUSED
EITHER WAY, and the last section of this document is why: the proxies are the
first of two gates and the second is on the consuming model's side. Each section
below is the evidence for one of them, because "nothing produces it" and "the
measurement exists and indexes on an axis this project does not have" are
different states and only the second says what to go and get.

| proxy | unit | verdict | the axis that stops it |
| --- | --- | --- | --- |
| Fe-Al oxide content | kg/kg of fine earth | not derivable | the total-to-extractable step, and the one held source that measures it is single-lithology and age-indexed |
| allophane content | kg/kg of fine earth | not derivable | the content inside andic material against the axis that controls it, which is leaching regime and not development; the two sources that index it there cannot be reached |
| aggregate capacity | dimensionless | not derivable | nothing in the pipeline resolves soil structure, and no held source transfers from what it does emit |
| polyvalent cation saturation | cmol(+)/kg | DERIVABLE, not wired | nothing; the emitted `cec` column and the pH-indexed cation share supply both halves. world-rzyz |

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

The near miss, and what it misses turns out not to be the number that was
looked for. `build_soil.py`'s `andisol_properties` emits `andic`, the AREAL
fraction of a gridcell whose volcanic glass has weathered to allophane. Turning
that into a concentration in the fine earth needs two conversions. The first is
arithmetic and this component already has its parts: an areal share becomes a
mass share through the bulk densities, and the andic bulk density is declared
and is already blended into the soil map's own density column. The second is the
allophane concentration WITHIN andic material.

THE AXIS THIS WAS ASKED FOR IS THE SUBORDINATE ONE, and a read source says so.
Parfitt (2009) reviews the controls on allophane formation in tephra and
Andisols and states that "the effect of time is subordinate to these factors",
the factors being the activity of silicic acid in the soil solution, the
availability of Al species, the opportunity for co-precipitation, leaching
regime, soil organic matter and pH. A distribution of allophane content against
a DEVELOPMENT measure would therefore be a distribution against the weakest of
the controls, and a transfer function built on it would carry the residual
variance of all the others.

The axis that does control it is leaching regime, and this pipeline has one.
Parfitt reports the New Zealand aeolian case in exactly the form the andisol
model already runs on: halloysite and no allophane where the rainfall is 1000
mm, only allophane where it is 1600 mm, the switch turning on soil-solution
silicon crossing about 10 mg/l. `pedogenesis.yaml`'s `andisol` block already
gates andic development on a precipitation threshold for that reason, so the
control is represented; what is absent is the CONTENT above the threshold, not
the threshold.

A distribution-shaped statement of the content does exist and Parfitt (2009)
carries two. Icelandic field soils span 2 to 22 per cent allophane at 11 to 42
per cent C (Sigfusson et al. 2008), and Atlantic volcanic island mineral
horizons hold under 6 per cent. Both are ranges over real soil populations
rather than the named national maxima of Parfitt (1990) Table 3, so the factor
of several that substituting maxima would have put into a protection
coefficient is no longer the state of the evidence. What neither supplies is the
content INDEXED on anything: they are unconditional ranges, so a cell cannot be
placed inside one from what this component emits.

So the proxy stays undeclared, and the reason has changed. It is no longer that
only extremes are held. It is that the two sources that would index content on
the controlling axis were sought and could not be reached from this host:
Parfitt, Russell and Orbell (1983), `10.1016/0016-7061(83)90029-0`, which
reports allophane along a weathering sequence, and Singleton et al. (1989),
`10.1071/SR9890067`, which reports allophane and halloysite content against
soil-solution silicon. Both return 403 from the publisher, from Google Scholar
and from Sci-Hub. `references/INDEX.md` carries both rows as sought and not
reached.

WHAT WOULD SETTLE IT, restated on the right axis: allophane content over
allophanic soils against a leaching measure -- rainfall, a water balance, or
soil-solution silicon activity -- rather than against a development measure or a
review's extremes. Even then the arm stays refused, for the reason the last
section of this document gives.

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
both a capacity and the share of it held by calcium, iron and aluminium. Both
halves now exist, and this section is the one of the four whose verdict has
moved.

The capacity is emitted. `pedogenesis.yaml`'s `exchange` block and
`build_soil.py:exchange_properties` derive cation exchange capacity additively
from the clay and organic-matter fractions the soil map already carries, on the
two slopes of Sahrawat (1983), and the soil map's `cec` column is the result.
That block carries what the relation does and does not license; the short form
is that its SHAPE is checked against a second continent and holds, and its LEVEL
rests on one region's fit at R2 0.60 with n = 40 and is an open exposure that
world-n4i0 owns.

The share is measurable and is measured against a variable this component
emits. Solly et al. (2020) partitioned effective CEC by cation over 1204 Swiss
forest profiles within nine pH classes: at pH at or above 5.5, exchangeable
calcium carries 59 to 83 per cent of it, and below pH 5 exchangeable aluminium
carries 21 to 44 per cent. Calcium and aluminium are both polyvalent, so
polyvalent saturation is high at both ends of that range and the quantity is
recoverable from the pH field for a reason rather than by interpolation. Iron
contributes under 1 per cent, which is why the arm's third cation does not need
its own field.

Chadwick et al. (2003) remains the direct measurement of the quantity and
remains unusable as a FIELD, for the three reasons this document gave before:
its effective CEC is non-monotonic in the leaching index, so a cell's value is
not recoverable from a one-sided argument about how wet it is; it is one parent
material and one substrate age band, and the paper states the age dependence in
this very quantity; and it states that the decline is irreversible under natural
conditions, so the state is a function of the climate PATH and a pipeline with
no time axis has no path to evaluate it on. What it IS used for is a check:
its Table 7 pairs pH with base saturation horizon by horizon and puts the
transition where Solly's calcium share puts it, on a different continent and a
different parent material. Two sources agreeing on where a threshold sits is
what licenses the ramp `exchange.base_saturation` declares.

So this proxy is derivable and the earlier verdict that no held source supplied
a cation exchange pedotransfer from texture and organic carbon no longer holds.
What remains is wiring: the proxy is declared in
`biosphere/config/mineral_reactivity.yaml` and the emitted capacity is not yet
connected to it. That is world-rzyz, and closing it does not open the arm --
the last section of this document says why.

## What this means for the arm, and what it does not mean

All four stay undeclared, so `mineral_reactivity_gate.py` keeps refusing the
mineral-aware arm and `active_arm` stays `texture_only`. Three of them are
undeclared because nothing here can derive them; the fourth because it is not
yet wired, which world-rzyz owns. What crosses the interface is unchanged: the
soil map's `andic` and `pfixation` columns, declared in the contract's
`carried_state` and read by no equation. The `cec` column the exchange complex
adds is not in `carried_state` either, and adding it is world-rzyz's job rather
than a consequence of emitting it.

One thing a later reader should not have to rediscover, and it is what stops
world-rzyz being read as progress toward the arm. Closing all four of these
would still not let the arm be written. Cotrufo et al. (2013) and Lehmann
and Kleber (2015) place mineral association, aggregation and accessibility at
the centre of stable soil organic matter formation and name the associating
surfaces, and neither supplies a parameterised transfer from any of the four to
a protection coefficient. The proxies are the first of two gates, not the only
one, and the second gate is on the consuming model's side.
