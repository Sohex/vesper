# Audit: Orogen's lithology model

*Audited 2026-08-16, against `vendor/orogen` at `cf-fork` 3907c58. Tasks in
`TASKS.md` under `LITH`.*

Prompted by a specific question -- is the absence of a gypsum class an oversight?
-- and widened because lithology now feeds albedo, erodibility, phosphorus and
brine chemistry, so an error there propagates further than it used to.

## The gypsum question: no, and for a structural reason

`saltCrustMask()` assigns the class purely geometrically. Every cell in a
preserved basin lying below `sinkElevation + frac * depth` becomes `evaporite`;
everything above becomes `playa_clastic`. No water chemistry enters at any point.

So Orogen's `evaporite` denotes **a position in a closed basin**, not a mineral.
Adding a gypsum class would require the generator to know inflow chemistry, which
it has no basis for -- the same reason it measures basin geometry but never
decides water levels, and the same reason `inland_water` is deliberately empty.

The decision belongs downstream and is now made there:
`pedology/scripts/brine_paths.py` applies the Hardie-Eugster chemical divide to
each basin's catchment lithology. See `pedology/notes/derived-surface-classes.md`.

**But the class carries an albedo, and that number is the real exposure.**
`evaporite` is 0.50, and the module's own docstring puts it at about 1.7 W/m2 per
0.10 of error on this planet. What the literature says about it is in the albedo
section below, and it is not what was expected: the dominant control is surface
state, not mineralogy. `LITH-1`, `LITH-2`.

## Nothing in the module is cited

Twenty classes carry an erodibility, a density and an albedo: sixty numbers, and
there is **no reference anywhere in the file** -- no author, no DOI, no dataset.
The docstring justifies erodibility as "ordered by rock strength", which is an
ordering argument and a reasonable one; it is not a grounding for the values.

The ordering being defensible while the values are not is the important
distinction, because erodibility is renormalised to a land mean of 1 and so only
the ordering is consumed. **Albedo is not renormalised.** It is a boundary
condition handed to ExoPlaSim as an absolute number, so for albedo the values
matter directly and the absence of a source is a real exposure. `LITH-3`.

## Albedo: no rock has a citable broadband value, and the gypsum premise is wrong

Grounded 2026-08-16 against sixteen fetched sources; see `references/INDEX.md`.

**There is no citable "granite = 0.30" in the literature.** Nothing tabulates
broadband albedo per rock type. What exists is directional-hemispherical
reflectance spectra -- ECOSTRESS v1.0 (473 rocks, and the correct measurement
geometry for a model albedo), ASTER v2.0, USGS splib07. The defensible
construction is to solar-weight those against **this world's own K2.5V
spectrum**, which is already in the tree at
`exoplasim/inputs/stellarspectra/k25v_hr.dat`. That is a method, not a lookup,
and it is what the table should be rebuilt from. `LITH-12`.

**The gypsum-is-brighter premise is wrong, and backwards in a way that matters.**
Gypsum crust does read brighter than halite in the visible -- Slater et al.
measured White Sands at 0.49 to 0.62 across the four visible and near-infrared
Landsat bands -- but it **collapses in the shortwave infrared**, to 0.414 at
1.55-1.75 um and **0.159** at 2.08-2.35 um, on gypsum's structural-water bands.
Anhydrous halite has no such bands and stays bright right through.

Band-weighting Slater's six measurements gives White Sands gypsum **0.528 under
a K2.5V spectrum** against 0.541 under the Sun, versus clean halite pans measured
at 0.55 to 0.65. So gypsum sits at or *below* halite in broadband, and the K
dwarf widens the gap rather than closing it: **22.1% of this star's shortwave
falls beyond 1.4 um against 13.6% of the Sun's**, so 1.6x more of the incident
flux lands in the window where gypsum absorbs and halite does not.

That 0.528 is a lower bound on the effect. The six Landsat bands leave the
1.4-1.5 and 1.9-2.0 um water bands unsampled, and renormalising over them
over-weights the visible.

**Splitting evaporite into gypsum and halite with gypsum set brighter would push
the error the wrong way.** `LITH-2` is rewritten accordingly.

**Surface state dominates mineralogy anyway.** Measured halite crusts span 0.18
to 0.75, and Kampf et al. attribute the range to roughness and detrital loading
rather than composition: Atacama nucleus with detritus 0.18, rough halite 0.25,
Bonneville dry crust 0.45, transition crust 0.49, Uyuni 0.55-0.65 over twenty
years of MODIS, Pilot Valley 0.64 annual mean. Wetting is a first-order switch --
0.64 to 0.24 at Pilot Valley, 0.45 to 0.22 at Bonneville under 20-40 mm of
flooding. On a dust-loaded, periodically flooded endorheic world the low end is
not exotic, and at 1.7 W/m2 per 0.10 the measured spread is worth about
+/-4 W/m2. **A static 0.50 is wrong by up to 0.28 wherever
`hydrography/surface_water.py` puts standing water or a damp crust**, which is a
coupling the albedo build should carry and currently does not. `LITH-13`.

Two smaller corrections, both supported:

- **`playa_clastic` 0.30 sits at the bright end.** Post et al. measured 52 US
  soils with a pyranometer at 0.048-0.402, mean 0.189; 0.30 corresponds to a
  Munsell value near 6, a genuinely pale dry mud. Damp or organic-tinged playa
  falls to 0.10-0.20.
- **The basalts at a shared 0.10 hide a real spread.** Li et al. measure tephra
  below 0.05 across 350-2500 nm and fresh lava comparably dark, rising above 0.10
  only once oxidised or lichen-colonised. 0.10 is the weathered state; fresh is
  roughly 0.04, and `morb`, `oib` and `flood_basalt` carrying one value conceals
  a 0.04-0.15 range.

Henderson-Sellers and Wilson's compilation brackets the coarse picture: salt
playas and light sand deserts 0.28-0.44, semidesert and light soils 0.20-0.33,
stony deserts and common soils 0.07-0.25.

## Erodibility: the ordering mostly holds, the spread does not

Grounded 2026-08-16 against ten fetched sources; see `references/INDEX.md`.

Moosdorf et al. (2018) is directly comparable because it normalises to acid
plutonic = 1.0, which is the same form Orogen uses. Putting Orogen on that basis
by dividing through its granite value of 0.40:

| class | Orogen / granite | Moosdorf | ratio |
| --- | --- | --- | --- |
| granite, granodiorite | 1.00, 1.12 | 1.0 | 1.0, 1.1 |
| gneiss | 0.87 | 1.0 | 0.9 |
| quartzite | 0.62 | 1.0 | 0.6 |
| schist | 2.75 | 1.0 | **2.8** |
| melange | 4.00 | 1.0 | **4.0** |
| basalts | 1.62-2.12 | 1.4 | 1.2-1.5 |
| arc andesite, rift bimodal | 2.25, 2.37 | 1.1 | **2.1, 2.2** |
| carbonate | 3.25 | 1.0 | **3.3** |
| clastics | 5.50-6.50 | 1.5 | **3.7-4.3** |
| pelagic, playa | 7.50, 7.00 | 3.2 | **2.3, 2.2** |
| evaporite | 8.75 | not indexed | -- |

**The ordering is broadly supported.** Quartzite most resistant (Sklar and
Dietrich use it as their strong-rock reference), granite/gneiss next (Stock and
Montgomery put granitoids and metamorphics in one K decade; Portenga measures
igneous outcrops at 8.7 m/Myr against sedimentary 20), clastics well above
crystalline, evaporite at the top (Frumkin: bare salt denudes at 100-120 mm/yr
against ordinary rock at 0.009-0.020, a factor of thousands).

**The spread is roughly four times too wide.** Orogen runs 14x from quartzite to
salt; Moosdorf's published global index spans 3.2x, and Zondervan et al. measure
the *fluvially expressed* contrast within a mountain belt at about 4x. Those two
agree, and they are the right comparison for a landscape model: intact rock
strength does vary by four to five orders of magnitude, but channels adjust
width and slope, so the contrast a landscape actually expresses is far smaller
than strength implies. Judged against strength alone 14x looks conservative;
judged against what landscapes do, it is about 4x too generous.

Since erodibility is renormalised to a land mean of 1, this does not change how
much erosion happens -- it over-differentiates *where*, giving too much relief
contrast between resistant and weak terrain.

**Two values are wrong rather than merely wide.**

- **`carbonate` 1.30 should be near granite, roughly 0.40-0.50.** Moosdorf places
  carbonate sedimentary rocks at 1.0, in the low-erodibility group with acid
  plutonic and metamorphic. Bursztyn's limestones are among his strongest units,
  at or above granite in tensile strength. The caveat is that this is *mechanical*
  resistance and carbonate instead loses mass to dissolution, which a
  stream-power multiplier cannot represent -- but on a mostly arid, mostly
  endorheic world the mechanical value is the right one.
- **`schist` 1.10 should be roughly 0.40-0.55.** Nothing separates schist from
  granite the way 2.75x does. Bursztyn's Vishnu schist is *stronger* than
  Zoroaster granite in tension and equal in compression. Placing it modestly
  above gneiss on foliation anisotropy is defensible; 2.75x is not.

Two more are questionable rather than settled. `melange` at 4.0x relative to
Moosdorf's single metamorphic bin may be legitimate -- subduction melange really
is weak, and Moosdorf's one metamorphic class is coarse against Orogen's four.
Being more differentiated *within* a Moosdorf class is not automatically wrong;
exceeding Moosdorf's entire between-class range is. And `pelagic` at 3.00 is
ambiguous: defensible if it means unconsolidated ooze, badly wrong if it means
lithified pelagic section, since ribbon chert is among the most resistant rocks
in a melange. The intended meaning needs pinning down.

## Density: one value to change

Every Orogen density checks against Daly, Manger and Clark (1966) tables 4-1 and
4-5: granite 2.65 against 2.667 (n=155), granodiorite 2.70 against 2.716 (n=11),
gneiss 2.75 within 2.61-2.84, schist 2.80 against 2.76-2.82, basalt family 2.9
against gabbro 2.976 and diabase 2.965, quartzite 2.65 = quartz, carbonate 2.70 =
calcite 2.71.

The exception is **`evaporite` 2.2, which should be 2.1**: Frumkin states rock
salt density twice explicitly and derives diapirism from it. Minor, and nothing
consumes density yet, but it is a directly stated primary value.

## Composition against GLiM

Hartmann and Moosdorf (2012) is the Earth comparator and is already held.
Computed from `surface_rock` weighted by `cell_area` over `surface_class == 1`
-- not from `manifest.lithology.compositionLand`, which uses `land_mask` as its
denominator and omits the dry sub-sea-level floors:

| group | Vesper | GLiM Earth | ratio |
| --- | --- | --- | --- |
| metamorphic | 26.73% | 13.0% | 2.1x |
| siliciclastic | 25.62% | 30.9% | 0.8x |
| playa fill | 23.78% | part of unconsolidated, 24.6% | -- |
| plutonic | 11.20% | 6.8% | 1.6x |
| carbonate | 5.56% | 7.8% | 0.7x |
| volcanic | 3.01% | 6.2% | 0.5x |
| evaporite | 2.85% | 0.3% | 9.5x |

The evaporite excess is plausibly the endorheic setting -- Earth is about
one-fifth endorheic by catchment and this world is far more -- but Earth's
evaporite *outcrop* is 0.3% against that one-fifth, so the relationship is not
linear and 9.5x wants an argument rather than an assumption. The volcanic
shortfall and the metamorphic excess are unexplained and may point at archetype
weighting. `LITH-5`.

**GLiM's largest land class has no Orogen counterpart at all**: unconsolidated
sediment, 24.6% of Earth's land, including alluvium 4.1%, dune sand 5.3% and
loess 1.1%. This is defensible if `surface_rock` is understood as consolidated
lithology plus basin fill, with surficial cover derived downstream -- which is
exactly this project's architecture. But the field is named `surface_rock` and
described as the exposed surface, which overclaims. `LITH-6`.

Also absent, in rough order of how much they would matter: basic plutonic
(gabbro, 0.7% of Earth, and Meybeck has a gabbro row), ultramafic and
serpentinite (chemically distinctive, high Mg, and it is what Meybeck's
`misc_metamorphic` mostly is), and pyroclastics (0.6%). `melange` partially
covers the ultramafic case.

## What is good

Worth recording, because an audit that only lists faults misrepresents the thing
it audited.

- The two-layer basement/cover model with erosion stripping cover to expose
  basement is a genuinely good fit to a snapshot model with no time axis, and it
  produces exhumed orogen cores and stripped cratons without a stratigraphic
  column.
- The cover sequence is an **explicit ordered table**, first match wins, precisely
  because branch order silently decided it once before and cost a terrain
  rebuild. That is the same lesson the derived-surface-class design is built on.
- The module is candid about its own load-bearing numbers. The docstring flags
  the evaporite albedo as the one most likely to matter and explains why the
  class was zoned.
- `erodibility` is renormalised to a land mean of 1, so lithology redistributes
  erosion instead of scaling it, and the existing erosion slider keeps its
  meaning. That is a careful piece of design.
