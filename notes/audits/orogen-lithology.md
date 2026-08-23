# Audit: Orogen's lithology model

*Audited 2026-08-16, against `vendor/orogen` at `cf-fork` 3907c58. Tasks in
the `bd` issue tracker under the `lith` label.*

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

**A broadband-albedo-per-rock-type table DOES exist, and an earlier claim here
that it did not was wrong.** Hunt (1982), *Spectroscopic properties of rocks and
minerals*, chapter 3 of Carmichael's *Handbook of Physical Properties of Rocks*
volume I, Table 9, has a column headed exactly **"Albedo (0.3-2.5 um)"** against
named rock samples: rhyolite 0.69, granite (Rockport) 0.68, granite (Wisconsin)
0.91, andesite (Mt Shasta) 0.56, andesite (Colorado) 0.42, syenite 0.61,
granodiorite 0.61. The chapter also consolidates the whole Hunt and Salisbury
*Modern Geology* series, its references 10 and 13.

**The column is not an absolute albedo, and the source now held says exactly what
it is.** Graphic granite is listed at 1.13, impossible for a reflectance
fraction. Logan et al. (1973), the table's source, states the method: "The
reflectance values were then integrated over the entire 0.35- to 2.5-um
wavelength range, and **relative albedos were obtained by dividing this number by
a value obtained in similar fashion for freshly prepared MgO**." So the column is
albedo relative to an MgO white standard, and a sample can exceed it.

**And the samples are 0-74 um POWDERS**, packed in a "fairy castle" structure,
not solid rock faces. That is the single most important qualifier: powdered rock
is far brighter than the same rock as an outcrop, which is precisely the
distinction Zhuang et al. (2023) was wanted for. Logan's granite at 0.68 relative
to MgO is a fine powder, not a granite hillside.

So the table gives a real, citable, systematically ordered measurement that still
cannot be dropped into `ROCK_CLASSES`. Converting it needs two corrections
neither of which is in hand: multiply by MgO's own integrated reflectance, and
then correct powder to solid surface.

What it does establish firmly is the ORDERING and its systematic dependence on
silica, over 72 rocks spanning acidic to ultrabasic:

| rock | % SiO2 | albedo / MgO |
| --- | --- | --- |
| rhyolite | 78.6 | 0.69 |
| granite, Rockport | 75.9 | 0.68 |
| biotite granite | 67.5 | 0.81 |
| syenite | 66.4 | 0.61 |
| andesite | 66.3 | 0.56 |
| granodiorite | 61.8 | 0.61 |
| diorite | 56.6 | 0.53 |
| hornblende diorite | 54.2 | 0.41 |
| basalt | 57.9 | 0.41 |
| gabbro | 57.7 | 0.50 |
| hypersthene gabbro | 48.2 | 0.34 |

Felsic bright, mafic dark, monotonic in silica apart from anorthosite and a
couple of glassy outliers. **Orogen's table has that ordering right already**:
granite 0.30 against basalt 0.10 is the same direction, and the same rough
factor of two to three that Logan measures between granite and gabbro.

The lesson is the one this project already has a rule for. The claim "no such
table exists" was made from two compilations that classify by land cover rather
than lithology -- Henderson-Sellers and Wilson (1983) and Coakley (2003) contain
zero occurrences of granite, basalt, limestone, sandstone, quartz, shale, gneiss,
bedrock or "rock type" between them -- and generalised from a checked absence in
two sources to an unchecked absence everywhere.

What also exists, and is what the table should ultimately be built from, is
directional-hemispherical reflectance spectra -- ECOSTRESS v1.0 (473 rocks, and the correct measurement
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

### Confirmed against USGS splib07, which also settles what the library can and cannot do

The library was obtained 2026-08-16 and tested directly rather than from its
documentation. Of 3,156 spectra in `splib07a`: granite 0, sandstone 0, gneiss 0,
schist 0, quartzite 0, andesite 0, gabbro 0, basalt 8, limestone 1. **So it does
not close the crystalline-rock gap.** But it carries halite 11, gypsum 11 and
trona 6, which are the load-bearing classes, so it does close the evaporite one.

Solar-weighted over 0.35-2.5 um against `k25v_hr.dat`, matched grain-size series:

| mineral | Sun | K2.5V | change |
| --- | --- | --- | --- |
| halite | 0.832 | 0.845 | **+0.013** |
| gypsum, selenite | 0.775 | 0.746 | **-0.029** |
| trona | 0.809 | 0.782 | **-0.027** |

**The stellar response is the robust part and it confirms the Slater result from
full spectra**: anhydrous halite gains under a redder star while both hydrous
salts lose, a relative shift of 0.042 in halite's favour, because 19.1% of this
star's shortwave falls beyond 1.4 um against 11.7% of the Sun's.

**The absolute ordering is not robust.** Across samples, halite spans 0.311 to
0.877 and gypsum 0.738 to 0.858, so **the within-mineral spread is several times
the between-mineral difference**. Which of the two is brighter depends on which
specimen you pick.

And every one of these is a pure mineral separate reading 0.74 to 0.88, against
field crusts measured at 0.18 to 0.65. The gap is roughness, detrital loading and
moisture -- exactly what Kampf et al. attribute the field range to.

**The synthesis, and it strengthens rather than weakens the recommendation.** Use
the library for spectral SHAPE and stellar response, where it is authoritative;
use the field measurements for absolute LEVEL, where the library is not
applicable. And do not split evaporite by mineralogy for albedo purposes: the
within-mineral spread exceeds the between-mineral difference, and surface state
dominates both.

**Surface state dominates mineralogy anyway.** Measured halite crusts span 0.18
to 0.75, and Kampf et al. attribute the range to roughness and detrital loading
rather than composition: Atacama nucleus with detritus 0.18, rough halite 0.25,
Bonneville dry crust 0.45, transition crust 0.49, Uyuni 0.55-0.65 over twenty
years of MODIS, Pilot Valley 0.64 annual mean. Wetting is a first-order switch --
0.64 to 0.24 at Pilot Valley, 0.45 to 0.22 at Bonneville under 20-40 mm of
flooding. On a dust-loaded, periodically flooded endorheic world the low end is
not exotic, and at 1.7 W/m2 per 0.10 the measured spread is worth about
+/-4 W/m2. **A static 0.50 is wrong by up to 0.28 wherever
`hydrography/scripts/surface_water.py` puts standing water or a damp crust**, which is a
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

### Corrected 2026-08-16: the arc classes, and an inverted ordering

The two bolded arc rows above are fixed. `arc_andesite` 0.90 -> 0.50 and
`rift_bimodal` 0.95 -> 0.50 bring both to about 1.25x granite against Moosdorf's
1.1, and `arc_basalt` 0.85 -> 0.65 puts island-arc basalt with `flood_basalt`,
which is where Moosdorf indexes them both.

The ordering was also **inverted**, which is the part worth keeping. Moosdorf
puts basic volcanic rocks at 1.4 and acid at 1.1: basalt weathers faster than
andesite, because its mafic minerals are the less stable ones. Orogen had
andesite *above* basalt at 0.90 against 0.65. Note that `--lithology-strength`
could never have caught this. It compresses the whole distribution toward the
mean and preserves relative ordering by construction, so it fixes a spread that
is too wide and is blind to one that runs backwards. The two corrections are
independent and both are needed.

`melange` is deliberately **kept** at 4.0x granite. Moosdorf's single metamorphic
bin is dominated by gneiss and schist and never contained melange, which is a
sheared block-in-matrix unit and genuinely among the weakest rock masses there
is. Being more differentiated than a coarse bin is not an error.

These values had never been checked because the rules that assign them could not
fire (see the arc-terrain section), so nothing they controlled was ever exposed.

**Measured, G_melange vs H_arcerod, both at `--lithology-strength 0.682`:** the
terrain barely moves. Land mean 0.3953 -> 0.3957 km, land above 1 km 12.412 ->
12.424%, land fraction and preserved basin count identical. Composition shifts
0.045 points from `granodiorite` to `arc_andesite`, the more resistant arc rock
now surviving where it was being stripped. Through albedo that is -0.00004 on
the land mean, which is nothing. The correction buys defensibility, not a
different world. Terrain hash fb3eb4ab -> 2e06d176.

### `rift_bimodal` is 0% of land, and that one is geography

Measured on the arc-fix build: `rift_bimodal` never appears. It is worth
separating from the arc bug, because the failure looks identical from the
composition table and is not the same thing at all.

The arc rules were unreachable -- a predicate that could not be true for the
cells the rule wanted. The rift rule is reachable and simply has nothing to
match. It fires inside a continental rift, and a continental rift needs a
divergent boundary with continental crust on BOTH sides. Boundary topology is
taken from super plates deliberately, so that small-plate boundaries do not
open rifts inside a coherent continent. Each of this world's ten continents is
its own continental super plate, and no two of them are adjacent: every one of
the 60,859 boundary regions on the planet has an oceanic plate on at least one
side, with zero exceptions in any boundary class.

Two consequences follow from the same fact, and both are properties of this
world rather than defects:

- **No continental rifting.** No rift valleys, no bimodal rift volcanics, and
  `rift_bimodal` is dead weight in the class table for this configuration.
- **No continent-continent collision anywhere.** Every convergent boundary has
  an oceanic plate on one side, so all orogeny here is Andean/arc type. There is
  no Himalayan analogue and no suture belt on this planet.

If either is unwanted it is a plate-configuration question, not a lithology one:
it needs two continental super plates placed in contact, which is upstream of
everything in this file.

### Albedo: the arc classes survived the check

No change needed. `arc_andesite` 0.20 against `granite` 0.30 is a ratio of 0.67
and Logan's andesite-to-granite powder ratio is 0.69. `granodiorite` 0.28 and
`arc_basalt` 0.13 interpolate between endpoints confirmed directly.

Logan's ratios are usable across that particular pair only because andesite and
granite are both felsic to intermediate. The powder-to-slab factor is itself
lithology-dependent -- 2.2 for basalt against 4.1 for granite in Paragas -- so
powder ratios systematically compress the felsic/mafic contrast and must not be
carried across it. Applied naively, Logan would put `flood_basalt` at 0.54x
granite where three direct slab and field routes put it at 0.33x. That is a trap
for anything that tries to extend this table from powder data.

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

Two more were flagged as ambiguous and both resolve on the class definitions
rather than on further measurement.

**`pelagic` at 3.00 is correct.** The ambiguity was whether the class means
unconsolidated ooze or lithified pelagic section, and the class settles it: it is
`Pelagic ooze / abyssal clay`, assigned in the ocean domain as a cover whose
thickness grows with seafloor age. Unconsolidated, and 3.00 sits right beside
Moosdorf's unconsolidated index of 3.2. Nothing to change.

**`melange` at 1.60 is right for what the class mostly is and wrong for what its
name also claims.** It is assigned where the subduction factor exceeds a
threshold on a convergent boundary, so it is genuinely subduction melange: a
sheared, block-in-matrix unit that is one of the weakest rock masses there is,
and 4x granite is defensible for it even though it exceeds Moosdorf's whole
between-class range, because Moosdorf's single metamorphic bin is dominated by
gneiss and schist and never contained melange. Being more differentiated *within*
a Moosdorf class is not an error; the coarse bin is.

The problem is the name, `Subduction melange / blueschist`, which conflates the
weak sheared matrix with blueschist -- a competent, high-grade metamorphic rock
that behaves nothing like it. One erodibility cannot serve both. Either the
class should be renamed to the melange it actually models, or blueschist should
be split out. Renaming is the honest fix, since the assignment rule is a
subduction threshold and produces melange, not blueschist.

## Measured: what correcting erodibility would actually cost

Four builds at full 2.5M resolution, 2026-08-16, about a minute each. `A_control`
reproduces the production build exactly -- land fraction 0.4317, mean 0.4070 km, max 4.5933 km, 3629 basins, evaporite 2.847%, carbonate 5.557% -- which validates the comparison.

| variant | mask cells differing | mean land dz | land >100 m | basins | finalPreserved spillDepth median | volume median |
| --- | --- | --- | --- | --- | --- | --- |
| B, spread narrowed to 4x via `--lithology-strength 0.682` | 0 of 2,500,001 | -12.6 m | 3.49% | 3629 | 6.3% | 12.4% |
| C, carbonate 1.30->0.45 and schist 1.10->0.45 | 0 of 2,500,001 | +15.6 m | 6.52% | 3629 | 10.1% | 14.9% |
| D, both together | 0 of 2,500,001 | -3.0 m | 3.88% | 3629 | 5.7% | 8.1% |

Three things fall out, and they separate cleanly.

**The land/sea mask is bit-identical in every variant.** Not close -- zero cells
of 2,500,001 differ. Erodibility redistributes erosion above sea level and does
not move the coastline, so ExoPlaSim surface code 172 would be unchanged.

**Basin identity is fully preserved**: 3629 basins, 100% shared ids, and the
natural catalogue's `depthKm` is unchanged to under 0.05 m because it describes
the pre-conditioning basin. So an existing carve verdict remaps rather than
restarting, exactly as it does across a gravity change.

**But the finished basin geometry moves, and that is what the verdict tests.**
`finalPreserved.volumeKm3` shifts by 8.1% in the median for the combined
correction, `spillDepthKm` by 5.7%, and more than 3,000 of 3,624 basins move by
over 1%. The carve verdict is a capacity-against-water-balance test, so
hydrography and the verdict must be recomputed. The climate need not be: a 3.0 m
mean orography change is worth about 0.02 K by lapse rate.

**The corrections partially cancel, which is why the fully corrected build is the
closest to the current one rather than the furthest.** Raising carbonate and
schist resistance lifts mean elevation; narrowing the spread lowers it. Doing
only one of the two moves the terrain about five times further than doing both.

## Rock albedo, grounded at last: slabs, not powders

The missing measurement was never a table of numbers -- it was the conversion
from laboratory powders to solid rock faces. **Paragas et al. (2025), via the
POSEIDON surface albedo database, measures the same rock as slab, crushed and
powder**, which is exactly that conversion and exactly what Zhuang et al. (2023)
was being chased for.

Solar-weighted over 0.3-4 um against `k25v_hr.dat`:

| rock | slab | crushed | powder | powder / slab |
| --- | --- | --- | --- | --- |
| basalt with olivine phenocrysts | 0.100 | 0.092 | 0.220 | 2.21 |
| K1919 basalt | 0.101 | 0.047 | 0.228 | 2.26 |
| basaltic andesite | 0.073 | 0.066 | 0.201 | 2.74 |
| STM 101 andesite | 0.082 | 0.075 | 0.298 | 3.63 |
| Dalmatian granite | 0.123 | 0.284 | 0.504 | 4.08 |
| olivine gabbronorite | 0.150 | 0.132 | 0.371 | 2.47 |

Median powder/slab **2.74**, range 2.19 to 5.11. That single factor reconciles
everything: Logan's granite powder at 0.68 relative to MgO and Orogen's 0.30 for
a granite landscape were never in contradiction, they were measuring different
things.

Three results follow, and they are what the albedo section of this audit has been
missing.

**Orogen's basalt at 0.10 is essentially exact.** Paragas basalt slabs read 0.100
and 0.101. Nothing to change on `morb`, `oib` or `flood_basalt` beyond the
fresh-versus-weathered spread already noted.

**Orogen's felsic values are about 2.4x too bright.** A granite slab is 0.123
against the table's 0.30, and andesite 0.082 against 0.20.

**And the felsic-to-mafic contrast collapses on solid rock.** Granite slab 0.123
against basalt slab 0.100 is a factor of 1.2, where powders show 2 to 3 and
Orogen's table assumes 3.0. The ORDERING that Logan established survives; its
SPREAD does not. So the table is wrong in the same direction as its erodibility
was -- too much contrast between classes, not too little.

One more, in the opposite direction to the evaporites: **silicates get brighter
under this star**, by +0.006 to +0.087, because they lack the structural-water
bands that make gypsum and trona lose. The stellar correction has opposite signs
for rock and for hydrous salt, so it cannot be applied as a single factor.

A caveat that keeps this honest: a real planetary surface is neither a cut slab
nor a sieved powder. Weathering rinds, regolith and dust cover push it toward the
particulate end, which is presumably why the Hu (2012) exoplanet surface set
reads 0.22 for basaltic and 0.52 for granitoid. The slab values are a floor and
the powder values a ceiling; where a given landscape sits between them is a
question about regolith, not about rock.

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
Computed from `substrate_class` (then named `surface_rock`) weighted by
`cell_area` over `surface_class == 1`
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

### Re-measured on the arc-bearing build, 2026-08-16

The table above was measured before the arc and forearc rules could fire. With
them firing the picture changes enough that two of its three flagged divergences
are gone and a new one has appeared:

| group | ratio to GLiM, before | after |
| --- | --- | --- |
| volcanic | 0.49x | **1.28x** |
| metamorphic | 2.06x | **1.42x** |
| plutonic | 1.65x | **2.71x** |
| carbonate | 0.71x | **0.49x** |
| evaporite | 9.50x | 9.49x |
| playa fill | 0.97x | 0.97x |
| siliciclastic | 0.83x | 0.80x |

**The volcanic shortfall was the unreachable arc rules**, and it is now a mild
excess rather than a factor of two deficit. The metamorphic excess halved for
the same reason: schist lost ground to the forearc.

**Two new things to explain, and they are probably one thing.** Plutonic doubled,
because `granodiorite` is the unroofed arc root and that rule fires now too;
carbonate fell, because the arc band and forearc province are COASTAL and
overwrite exactly the shelf where carbonate is assigned. So the arc fix moved
four groups, and the two that got worse are both consequences of arc terrain
appearing where shelf terrain used to be. Whether that is right is a question
about how wide the forearc province should be at 15 km resolution, not about
archetype weighting.

### The plutonic excess is the arc root, and it is a resurfacing question

Left open when the evaporite item closed, because it is a different finding that
the arc fix created rather than revealed.

Plutonic went 1.65x GLiM to 2.71x when the arc rules began firing, and the
granite half is not what moved. `granodiorite` -- the unroofed arc root -- went
from absent to 7.29% of land, against GLiM's intermediate plutonic at 0.4% and
acid plutonic at 5.7%. However it is binned, that is an order of magnitude more
exposed arc-root batholith than Earth carries.

The mechanism is in `elevation.js`: granodiorite is the BASEMENT under the arc
belt, so every arc cell whose cover erodes away exposes batholith. On Earth the
exhumed batholiths are ancient arcs -- the Sierra Nevada, the Coast Ranges --
while an ACTIVE arc is buried under its own volcanics, because it is resurfaced
faster than it is unroofed.

**Two explanations were proposed and MEASUREMENT refuted both.** They are
recorded because the measurement is the useful part.

The first was that Orogen cannot separate an active arc from a stripped extinct
one without a resurfacing rate. The second was that resurfacing is representable
as a cover erosion cannot outpace, making `LITHO_COVER_ARC_KM` -- 1.0 km, thinner
than a LIP's 1.2 -- the culprit.

Builds at 1.0, 2.0 and 3.0 km say otherwise. Tripling the cover moves
granodiorite from 7.29% of land to 7.10%, and the stripped share of the arc belt
from 59.5% to 58.0%. The belt itself does not move at all, at 12.24%.

The direct measurement says why. Of the granodiorite cells, **94.1% were assigned
arc_andesite cover and eroded it to zero** -- maximum surviving thickness 0.0001
km. Cells that keep their cover retain 0.494 km of the 1.0 assigned. So arc-belt
erosion is bimodal: a cell either loses under a kilometre or loses more than
three, and almost nothing lies between. No plausible cover thickness sits inside
that gap, which is exactly why tripling it did nothing.

**So the arc root is exposed because arc terrain erodes deeply, which is
physical**, and there is that much arc terrain because this world was given 100
plates deliberately. Both are consequences of choices rather than defects, and
the plutonic excess against GLiM follows from them. What would change it is the
plate count or the erosion sliders, and both are worldbuilding decisions that
have already been made and liked.

### The craton archetype effectively never fires

`LITH-21`, and it is the fourth rule this audit has found producing almost
nothing.

Granite is the FALLBACK basement -- `else bm = R.granite` -- so it takes every
continental cell that is not forearc, arc belt, fold belt or craton. On the
current build that is 58.87% of land, of which 11.12% is exposed where cover has
stripped. Against GLiM's acid plutonic at 5.7% that is 1.95x.

The craton branch is what should be claiming the stable interiors, and it claims
0.75% of land. `r_t_craton = max(0, 1 - 2.5 * tecActivity) * (1 - basin)` needs
tectonic activity below 0.40 to be non-zero at all, and below about 0.22 to clear
the 0.45 threshold. Land-mean activity here is 0.668.

**The obvious explanation is wrong.** With 100 plates against Orogen's default of
24, the natural reading is a constant calibrated for a quieter world. Builds at
both plate counts say otherwise:

| | 24 plates | 100 plates | GLiM |
| --- | ---: | ---: | ---: |
| granite | 11.16% | 11.12% | 5.7% |
| gneiss | 0.44% | 0.49% | -- |
| schist | 4.16% | 11.83% | -- |
| plutonic | 2.13x | 2.71x | 6.8% |
| metamorphic | 0.67x | 1.42x | 13.0% |

Granite is 11.1% at BOTH, and gneiss is under half a percent at both. Plate count
moves schist, granodiorite and melange; it does not move the two classes this
finding is about. So the craton term is near-dead at the default plate count too,
and this is not a regime mismatch.

**What it produces is a category substitution.** Earth maps a great deal of
undifferentiated Precambrian basement as gneiss and granulite, inside the 13.0%
metamorphic class. Orogen labels the same material granite, inside plutonic. The
two errors are complementary and visible at 24 plates, where plutonic runs 2.13x
GLiM while metamorphic runs 0.67x -- an excess and a deficit of similar size,
which is what one misfiled class looks like.

**Traced, and the levers were working the whole time.** The three sweeps below
measured SURFACE gneiss, which does not move. Basement does:

| | current | mult 1.1 | mult 0.9 | half basin | no basin |
| --- | ---: | ---: | ---: | ---: | ---: |
| gneiss basement | 1.41% | 2.23% | 3.11% | 3.06% | 3.46% |
| granite basement | 46.63% | 45.81% | 44.94% | 44.98% | 44.58% |
| gneiss SURFACE | 0.49% | 0.49% | 0.49% | 0.49% | 0.49% |

So the craton branch fires, responds to its parameters, and takes basement from
the granite fallback exactly as intended. What it cannot do is reach the surface:
the cells it claims are quiet, low-relief ground, and quiet ground here carries
the thickest cover. The cover over gneiss basement is playa fill at 35.7%,
continental clastics at 15.4% and evaporite at 15.2% -- these cratons are buried
under basins.

That is the shield-versus-platform distinction, and it is correct behaviour
rather than a defect. Earth's shields are exposed because glacial scour and
epeirogenic uplift stripped their platforms, and Orogen models neither. **Visible
shields are therefore not reachable by any parameter here**; cratonic basement is.

### A granitoid-gneiss-granulite continuum: possible, cheap, and not worth it

Assessed and declined.

**It needs no new rock class.** `gneiss` already exists and is fully grounded --
erodibility 0.35, density 2.75, albedo 0.28, phosphorus 790 ppm, Meybeck gneiss.
Relabelling deeply exhumed granite as gneiss costs one rule and no new grounding.
The axis is real too: metamorphic grade is depth of formation, and `erosionDelta`
is an exhumation proxy the export already carries.

**Granulite is better covered than assumed.** Daly (1966) has "Granulite,
Lapland" at 2.93 g/cm3 hypersthene-bearing and 2.73 hypersthene-free; Moosdorf
bins it metamorphic at 1.0, Hartmann as MT at 790 ppm, Meybeck as gneiss, and
the pH parent class is metamorphic. Only albedo is genuinely absent.

**It fails on usefulness instead.** Relabelling trades one divergence for
another rather than reducing them:

| split | plutonic /GLiM | metamorphic /GLiM | sum of absolute error |
| --- | ---: | ---: | ---: |
| as built | 2.70 | 1.42 | 2.12 |
| granite > p50 -> gneiss | 1.89 | 1.84 | 1.73 |
| granite > p75 -> gneiss | 2.29 | 1.63 | 1.92 |
| granite > p90 -> gneiss | 2.54 | 1.50 | 2.04 |

Both classes are already over-represented, so moving area between them improves
one ratio by spoiling the other. The underlying fact is that this world exposes
more crystalline basement than Earth full stop -- 59% of continental basement is
the granite fallback and cover strips widely -- and a relabelling changes what
the table SAYS without changing how much basement outcrops.

If it is ever wanted for legibility rather than for the Earth comparison, p90 is
the defensible rule: the deepest-exhumed 1.10% of land, a statement about
exhumation rather than a fit to GLiM.

**The one real defect was the basin factor.** `r_t_craton` was multiplied by
`(1 - basin)`, which says a craton under sedimentary cover is not a craton. On
Earth most cratonic area IS covered -- the Russian Platform, the North American
mid-continent -- and cover is already modelled separately as `cover_rock` and
`cover_thickness`, so the factor double-counted and left the archetype claiming
under 1% of land. Removed. Cratonic basement roughly doubles to 3.46%, about a
third of Earth's cratonic share, which is the intent for a planet whose mass
implies keel-delaminating convection.

swept.

Whether cratons SHOULD exist here is a separate and easier question. This world
currently has none, and Earth's are prominent -- but a 1.881 Earth-mass planet
carries roughly twice Earth's radiogenic inventory, and vigorous convection
delaminates the depleted lithospheric keels that make a craton survive. So less
cratonic area than Earth is well supported by the mass alone, without appealing
to the planet being young. Some, but not much.

`LITH-20`, closed. The residual granite excess -- 11.12% of land against GLiM's
5.7% acid plutonic, and cratonic rather than arc -- is untouched by any of this.

### The evaporite 9.5x is a category error, not an excess

Resolved 2026-08-16, and the resolution is in the GLiM paper itself.

**GLiM's first-level `ev` is a DOMINANCE class.** A map unit is `ev` only if it
was interpreted as dominated by evaporites. Hartmann and Moosdorf say plainly
what that omits, using Asia as the worked example: evaporites there "are rarely
dominant and are therefore mapped only in a few areas as lithological class
(xx = ev: 0.28%). However, they occur subordinately in other lithological units
covering a far larger area (zz = ev: 8.52%)." Globally the same pair is **0.3%
dominant against 3.8% present** -- a factor of thirteen inside the same dataset,
between two numbers that both mean "evaporite".

Against evaporite PRESENCE, this world's 2.85% is **0.75x Earth**, not 9.5x. The
divergence was measured against the wrong column.

**And 2.85% is not the number to compare anyway.** Orogen assigns salt crust by
a geometric proxy -- the deepest quarter of a basin's relief -- and says so at
the call site: "deliberately a geometric proxy for flooding frequency and not a
water balance". `LITH-17` replaced it downstream with the ephemeral zone from the
actual water balance, which is **0.47% of land**. That is 0.12x Earth's evaporite
presence and 1.6x its dominant outcrop.

So there is no evaporite excess to explain. On the physically derived number this
world has LESS surface evaporite than Earth relative to how much evaporite Earth
actually has, which is what a dry world with low runoff should give: a salt crust
needs a sump that repeatedly floods and dries, and flooding is what this world is
short of.

The general lesson is the one this project keeps relearning. Both numbers were
called "evaporite", both came from GLiM, and picking the wrong one moved the
answer by a factor of thirteen. `LITH-5`, done.

**GLiM's largest land class has no Orogen counterpart at all**: unconsolidated
sediment, 24.6% of Earth's land, including alluvium 4.1%, dune sand 5.3% and
loess 1.1%. This is defensible if the field is understood as consolidated
lithology plus basin fill, with surficial cover derived downstream -- which is
exactly this project's architecture. The architecture was never the problem; the
NAME was. `surface_rock`, described as the exposed rock class, claimed to be the
thing this project derives downstream.

Renamed to `substrate_class`, and the description now says what it excludes as
well as what it holds. `LITH-6`, done.

One correction to the finding above while resolving it: `playa_clastic` covers
23.79% of land against GLiM's unconsolidated 24.6%, so the AREA is matched
almost exactly. The class correspondence still fails, because Earth's
unconsolidated includes dune sand and loess outside basins and this world has no
class for those -- but they are missing as a surficial-cover question, which is
`SURF-3`, not as a bedrock one.

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

## Playa clastics: the two groundings reconciled, and the class is too dark

Measured 2026-08-17, by `analysis/playa_albedo.py`. The class is about a quarter
of `precarve-craton`'s land and the largest single lever on this planet's energy
balance, and it carried two grounded values that disagreed by a factor of two.

`vendor/orogen/js/lithology.js` moved it from 0.30 to **0.19** citing Post et al.
(2000): 52 pyranometer measurements over 26 US soils, mean 0.189.
Henderson-Sellers and Wilson (1983) put salt playas and light sand deserts at
**0.28-0.44**. Both are read and both are in `references/INDEX.md`.

**They are not measuring the same thing, and the resolution is that each is
right about a different half of the question.**

Post is the right KIND of measurement -- a field pyranometer over a real surface
is exactly what a model albedo wants -- and the wrong POPULATION, since 26 US
agricultural soils are not a sample of playa mud and desert-varnished fan gravel.
Henderson-Sellers is a compilation whose "salt playa" is this project's OTHER
class: a salt-encrusted surface is `evaporite` at 0.50, and `playa_clastic` is
the clastic apron around and under it.

The ECOSTRESS library supplies the missing term, because it measures saline
desert soils as directional hemispherical reflectance -- the geometry a model
albedo wants -- but as prepared laboratory samples rather than crusted field
surfaces.

| quantity | value |
| --- | ---: |
| ECOSTRESS soils, all, sun-weighted over Post's 0.3-2.8 um band, n=69 | 0.289 |
| Post et al., field pyranometer, n=26 | 0.189 |
| laboratory over field | 1.53x |
| ECOSTRESS desert and saline soils, same weighting, n=9 | 0.330 |
| desert and saline over the whole soil population, same preparation | 1.141x |
| this star over the Sun, from the same spectra | 1.067x |
| **reconciled field albedo under a K2.5V** | **0.230** |

**The preparation offset never enters the answer**, which is the point of taking
a ratio inside one library: numerator and denominator share it and it cancels.
That 1.53x is the same trap this file already records for powders against slabs,
where the factor runs 2.2 to 4.1 and is lithology-dependent; soils are gentler
but not negligible.

What does not cancel is that Post's soils and ECOSTRESS's soils are different
soils, so the 1.141x carries whatever real difference sits between the two
populations as well as the material's own brightness. That is the residual
uncertainty here, and it is why this is reported as 0.23 rather than 0.230.

**So 0.19 is too dark, by about 20%, and 0.30 was too bright by about 30%.** On
this terrain the correction is worth about +0.0096 on land-mean albedo, near
-1.3 W/m2 and about +0.9 K, which partly offsets the darkening the lakes bring.
The mechanism for applying it without a new build is
`model.lithology_albedo_overrides`, whose `replaces` guard refuses to fire
against a class whose exported value has since moved.
