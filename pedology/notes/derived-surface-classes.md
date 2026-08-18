# Derived surface classes

*Design 2026-08-16, implemented and measured 2026-08-17 on the pre-carve build.*
`pedology/scripts/build_surface_classes.py` is the classifier and
`pedology/config/surface_classes.yaml` holds every threshold with its source.
Numbers in this file are measurements, dated, on a climatology, lake solution
and dust field that the next carve replaces.

World Orogen gives primary lithology: what the rock IS. It does not give what the
surface has BECOME under this world's climate and drainage. That second thing is
what albedo, dust emission and the phosphorus cycle actually key on, and it is
derivable here rather than upstream, because it depends on precipitation,
evaporation and lake extent that Orogen never sees.

Four candidates were proposed. One survives unchanged, one was dropped, and two
had their mechanism backwards -- pavement, which is a dust sink rather than
deflation armour, and loess, which is limited by whether a surface can trap dust
rather than by how much falls on it.

## What the literature changed about the proposal

**Desert pavement is a dust SINK, not deflation armour.** This is the correction
that matters most, because the proposed class name encoded the wrong mechanism.
Wells et al. (1995) dated pavement clasts with cosmogenic 3He and found they have
been at the surface the entire time: *"stone pavements are born at the surface."*
Deflation lag, water winnowing and shrink-swell are all rejected. Dust falls
through the clast mosaic and accumulates BENEATH it as a cumulic Av horizon, so a
pavement is an accretionary surface under a non-erodible cover. Citing it as
evidence that deflation armours a surface inverts the paper.

The consequence for us is not cosmetic. As "deflation armour" the class would
have been a *product* of dust emission and would have increased with it. As
pavement it is a *suppressor* of emission and a store of deposited dust, so it
enters the dust model with the opposite sign and enters the phosphorus budget as
a sink rather than a source.

**Silcrete gets no climatic threshold, and that is not a reason to drop it.**
Fenske et al. (2025) state that climatic calibration "will not provide duricrust
formation boundaries" for silcrete, because it forms across a broad range of
environments. The first draft of this note read that as "unplaceable, drop it".
That over-reads the claim: what it rules out is placing silcrete *by climate*.

Ullyott and Nash (2016) recognise four types, and the distinction does the work:

| type | setting | control |
| --- | --- | --- |
| pedogenic | within soil profiles on stable palaeosurfaces, >10^6 yr, laterally extensive | climate, and only coarsely: tropical to subtropical with alternating wet and dry seasons |
| groundwater | at or near a water table, or at zones of groundwater outflow | hydrological |
| drainage-line | silicification of alluvial fills in fluvial systems | hydrological |
| pan/lacustrine | silicification of sediments at the margins of ephemeral lake basins | hydrological |

The three non-pedogenic types form under climates "ranging from cold to arid",
so climate carries no information about them -- and does not need to. Their
controls are water tables, drainage lines and ephemeral lake margins, all of
which `hydrography/` already resolves. Fenske's model is itself a water-table
fluctuation model, which is the mechanism for the second row.

**Only pedogenic silcrete is dropped.** The other three are placeable here, and
on this world they are placeable unusually well, because a landscape where most
land drains to a closed basin is mostly ephemeral lake margin and interior
drainage line. Gypcrete and calcrete stay as before, on their measured bounds.

**Loess turned out to be trapping-limited, not flux-limited**, and that is the
second mechanism correction in this note. The design said to calibrate against
Muhs (2013) rather than Bettis et al. (2003), because Bettis is last-glacial
mid-continent North America and explicitly not a global typical. Muhs was read
for that on 2026-08-17 and **carries no loess accumulation rates at all**: no
table, no global range, no threshold. It is a qualitative review whose only
g/m2/yr figures are the legend bins of one redrawn model figure.

What it gives instead is better than the number it was fetched for. The Sahara
sits in Muhs' highest modern flux bin, 100-200 g/m2/yr, and produces essentially
no loess, because "arid regions rarely have the necessary conditions for trapping
loess" (Tsoar and Pye 1987, quoted at his p. 14-15). Australia is a major dust
source with little true loess. Canada had the silt supply and vegetation took the
surface before the wind could. So a rule keyed on deposition alone paints loess
exactly where Earth does not have it: on the source.

The rule therefore has two halves. A deposition threshold, which is necessary,
and a trapping surface, which is what makes it sufficient: not a barren class,
and wet enough to hold what lands on it. The threshold itself is read off the one
calibration Muhs' Fig. 44 does support, the mid-continental North American loess
belt forming Peoria Loess at 100-200 g/m2/yr in the last glacial and forming none
at 1-5 today at the same place. 50 g/m2/yr sits inside that interval and is
declared, with a 10-200 bracket around it.

## Two axes, not one chain

The `carved-zoned-v2` build lost basin fill because a lithology cover chain
decided which deposit sat on top by the order the branches happened to be
written. 203 preserved basins ended up with no fill cell anywhere. The lesson is
not "be careful with branch order" -- it is that **things which do not compete
for the same physical position must not be resolved in the same chain.**

Duricrusts are not surface covers. They are pedogenic horizons that form *within*
a profile, at or below the wetting front, and they become a surface only where
erosion has stripped what lay above. Loess and diatomite are surficial deposits
that sit ON the profile. Putting both in one precedence list forces an arbitrary
answer to a question that has none, and the arbitrary answer is whatever the code
was written in.

So the classifier emits two independent fields:

| field | question | values |
| --- | --- | --- |
| `surface_cover` | what wind and light see | `water`, `loess`, `diatomite`, `pavement`, `bare`, `soil` |
| `duricrust` | what is cementing at or below the surface | `gypcrete`, `calcrete`, `silcrete`, `none` |

A cell can be loess over calcrete, and that is a real landscape rather than a
contradiction to be resolved. Precedence exists only *within* an axis, is
declared as an ordered list, and every cell records which rule fired and which
other rules it also satisfied. A rule that never fires and a rule that always
fires are both bugs, and neither is visible without that record.

## `surface_cover`, in precedence order

Ordered by what physically sits on top, most recent deposit first.

1. **`water`** -- standing water from `surface_water.nc`. Anything beneath a lake
   is not a subaerial surface class. This is also the field that finally fills
   Orogen's deliberately empty `surface_class == 2`.
2. **`loess`** -- net aeolian deposition above 50 g/m2/yr ON A SURFACE THAT CAN
   TRAP IT: not a barren class, and above the 250 mm/yr desert boundary. Do not
   calibrate to Bettis et al. (2003), whose rates reach 17,500 g/m2/yr, because
   that is last-glacial mid-continent North America and the high end of the
   global range. Muhs (2013) has no rates to calibrate against either; what it
   has is the trapping argument above, which is the load-bearing half.
3. **`diatomite`** -- lacustrine silica from a lake that was persistent and
   productive and has since desiccated. See below; this is the interesting one.
4. **`pavement`** -- clast supply from a coarse or stony substrate, plus aridity,
   plus net dust deposition slow enough not to bury the mosaic. Suppresses
   emission and stores deposited dust.
5. **`bare`** / **`soil`** -- the fallback pair, split on whether the biosphere
   holds the surface. Already effectively present via the albedo work.

The ordering of 2 over 3 is the one place this is a genuine physical claim rather
than bookkeeping: active loess deposition buries an exposed diatomite bed. If the
dust work later shows deposition and deflation are co-located on this world, that
claim needs revisiting rather than inheriting.

## `duricrust`

Named for what it is rather than for how it forms, because the three values do
not share a mechanism. Gypcrete and calcrete are pedogenic and climate-keyed;
the silcretes this world can place are explicitly NON-pedogenic and
hydrology-keyed. Calling this axis `pedogenic_horizon`, as an earlier draft did,
would have been a name that quietly contradicted two of its own values.

Gypcrete and calcrete both need a two-sided window, and the two-sidedness is the
point: too dry and solutions never infiltrate to a wetting front, too wet and the
salt leaches straight through. Bachman and Machette (1977) give carbonate accumulation of
**0.22 to 0.51 g/cm2/kyr** and state the constraint from both directions.

| class | window | source |
| --- | --- | --- |
| `gypcrete` | below ~250 mm/yr AND potential evaporation exceeding precipitation in *every* month | Watson (1983) |
| `calcrete` | below 500-600 mm/yr, optimum 100-500; upper bound arguably 1000, lower as low as 50 | Alonso-Zarza (2003) |

Silcrete is not in that table because it has no row to put there: it is placed by
setting, below, and any climatic threshold written for it would be invented.

Gypcrete and calcrete overlap between 100 and 250 mm/yr, and Watson resolves it
explicitly: gypsum crusts occupy the driest zones, calcretes the less arid parts. So gypcrete
takes precedence where its monthly condition is met, and calcrete otherwise. Both
additionally require an ion source, which is lithological and is the point where
this meets the major-ion work that is still waiting on Meybeck's Table 5.

Note that the monthly PET-exceeds-P test needs monthly fields, not annual means.
Vesper's year is 183 days in 12 bins, and an annual mean would pass cells whose
wet season disqualifies them.

## Silcrete, placed by setting, and why it is worth placing

Silcrete is the one class here that earns its place on a cultural axis as much as
a physical one, and the two happen to align.

**The type this world can place is the type that is workable without fire.**
Ullyott and Nash state that pedogenic silcrete requires heat treatment before use
in stone tool manufacture while non-pedogenic silcrete is workable untreated --
and non-pedogenic is exactly the hydrologically-controlled set. So dropping the
pedogenic type costs nothing on this axis; it removes the variety that would have
needed a fire to be useful.

That heat treatment is real and cheap where it is needed. Schmidt et al. (2013)
put the Si-OH + HO-Si -> Si-O-Si transformation onset between **200 and 300 C**,
healing crystal defects and closing pores; effective temperatures run higher than
flint's, but the tolerated heating rate is fast enough that silcrete needs no
dedicated hearth and can be treated alongside other fire use. Brown et al. (2009)
show it was practised regularly by 72 ka and as early as 164 ka. If this world
ever grows an archaeology, that is the shape of it.

**The Earth analogue is the right one.** Nash et al. (2013) provenanced Middle
Stone Age silcrete at Tsodilo Hills to sources **220 km** away at Lake Ngami and
**295 km** on the Boteti River, in a direct line and further on the ground --
carried past adequate local quartz and quartzite. Those sources are pan-margin
and drainage-line silcretes in the Okavango, an endorheic system. A world that is
mostly interior drainage is a world with a lot of that setting, and this is
evidence that such a resource is worth carrying a long way past substitutes.

Two cautions. Silcrete distribution is **patchy rather than continuous** even
where the setting is right, so the classifier should produce sparse outcrops, not
a blanket. And silicification needs a dissolved silica source.

**The silica source is now resolved, measured 2026-08-16 on the arc-bearing
build.** `pedology/scripts/weathering_fluxes.py` computes it from the same
Meybeck Table 2C the brine paths use. Three results place it:

- **68.4% of this world's dissolved silica is delivered to closed basins**,
  measured 2026-08-17 by routing the release down the drainage network. That is
  the silcrete and diatomite supply, and it is large because the drainage is
  largely interior. Silica reaching the ocean is diluted into an enormous
  reservoir; silica reaching a closed basin has nowhere to go and concentrates
  until it saturates. Eugster and Jones' behaviour type V is exactly that, SiO2
  constant after saturation with a solid, so a closed basin fed by silica-rich
  runoff is a silica-precipitating setting by construction.

  **This number was previously recorded as 21.6% and that was the wrong
  measurement for the claim.** 21.5% is the share of silica RELEASED ON basin-floor
  cells, the `is_endorheic` weighting `weathering_fluxes.py` reports. Delivery is a
  drainage question and the drainage weighting is `terminal >= 0`, which covers
  76.2% of land against a basin floor of about a fifth of it. The two differ by a
  factor of 3.2 and only one of them is about where the silica goes. This is the
  same trap `notes/failure-modes.md` records under "one quantity, two meanings,
  three times the value", now caught a second time on a second quantity.
- **Volcanic terrain supplies 13.8% of it from 8% of land**, so it is enriched
  but not dominant. Meybeck puts volcanic rock at 200 umol/l against granite's
  150, the highest of the common classes. This is a floor rather than a central
  estimate: Meybeck's class is ordinary volcanic terrain, and Durr et al. (2011)
  note fresh unweathered ash releases far more, with local yields above
  50 t SiO2/km2/yr against a global mean of 3.3.
- **The whole-land yield is about 1.1 t SiO2/km2/yr against Earth's exorheic
  mean of 3.3.** This world is silica-poor per unit area because it is dry, and
  its total delivery is respectable only because it has twice Earth's land.

So the answer to the open question is: the silica is there, it is not primarily
volcanic, and the closed basins concentrate it. Diatomaceous basin fill is a
consumer of that silica rather than its origin, which reverses the guess this
section previously carried.

Which basins is now placed. `brine_paths.py` carries `silica_mol_per_year` per
basin, routed as concentration times local runoff over the catchment, and
`build_surface_classes.py` gates both silcrete and diatomite on it.

The gate does not bind. Every pan margin and drainage line inside an endorheic
catchment clears 100 umol/l, so **silica is not what limits silcrete on this
world; the setting is.** That follows from the 68.4% above rather than
contradicting it, and it means the interesting question for silcrete is
hydrological rather than geochemical.

## Diatomite, and why the stellar cycle matters here

Diatomite is the class that connects to everything else. Hudson-Edwards et al.
(2014) measured Bodélé diatomite at 600-610 ppm P against **40 ppm in aeolian
sand** -- but this world's land mean is 628 ppm, so diatomite is not P-rich
*material*. The Bodélé is an enormous, exceptionally deflatable source of
ordinary-P material, and its export is a mass-flux result rather than a
concentration one. The 15x contrast is against sand, not against crust.

Its formation needs a lake that is persistent enough to be productive and then
desiccates enough to be deflated, and Bristow and Moller (2018) supply the
mechanism for why it deflates so readily once exposed: saltating diatomite grains
shatter each other, so the surface manufactures its own dust.

**That is a cycle problem, not a steady-state one.** A basin permanently full
never exposes its bed; a basin permanently dry never accumulates one. The
interesting case is a basin that alternates, and this world now has a forcing
that does exactly that on a geomorphically relevant timescale: the 57-year
component of the stellar cycle is slow enough that lake level equilibrates rather
than merely breathing. Diatomite extent should therefore be assessed against lake
occupancy over the long cycle, not against a single climatology's lake mask.

The user's framing is the right test and should be preserved as one: *what
conditions allowed the Bodélé to form as it is, and can this pipeline produce
that result if the conditions are there?* Not forced -- permitted. A classifier
that hardcodes a Bodélé analogue has answered nothing.

Yu et al. (2020) dispute that the Bodélé is the primary Amazon dust source at
all. That premise should be adopted or rejected deliberately rather than
inherited, and nothing here depends on it: the deflation mechanism and the P
speciation stand whether or not the Amazon connection does.

## Closed-basin brine evolution: what the endorheic half needs

Meybeck gives the ions a lithology **releases**. In an exorheic catchment they
leave. Here mostly they do not, so a second step is required, and it turns out to
be the step that decides which evaporite class a basin gets.

Eugster and Jones (1979) sort all eight major solutes into five behaviours under
evaporative concentration:

| type | behaviour | species |
| --- | --- | --- |
| I | conserved throughout; measures the concentration factor | Na, Cl (Cl to halite saturation) |
| II | cation-anion pair precipitating a mineral; the minor one crashes at the branching point | Ca, CO3 |
| III | gradual removal across the whole range, linear but shallower | HCO3 + CO3 |
| IV | sigmoid, removed only mid-range by exchange, sorption, biogenic reduction | K, SO4 |
| V | constant once saturated with a solid phase | SiO2 |

The branching point in type II is the **chemical divide** of Hardie and Eugster
(1970). Calcite precipitation removes Ca and CO3 in equal equivalents, so
whichever is in excess at saturation dominates every later step: the first
calcite decides whether the brine becomes carbonate-rich or carbonate-poor, and
nothing downstream can undo it. Deocampo and Jones (2014) make it quantitative
with the Spencer Triangle in Ca-SO4-(HCO3+CO3).

**Those are three of Table 2C's columns**, so the release table and the fate
model compose without an intermediate step. Running our lithologies through the
calcite divide on their released Ca vs HCO3:

| lithology | Ca/HCO3 | brine path |
| --- | --- | --- |
| peridotite | 0.11 | alkaline |
| granite | 0.30 | alkaline |
| volcanic | 0.36 | alkaline |
| gneiss | 0.44 | alkaline |
| sandstone | 0.70 | alkaline |
| shale | 0.70 | alkaline |
| misc. metamorphic | 0.79 | alkaline |
| sedimentary carbonate | 0.80 | alkaline |
| halite evaporite | 1.50 | Ca-rich |
| gypsum evaporite | 3.25 | Ca-rich |

**Every silicate and carbonate lithology falls on the alkaline side, and only the
evaporites cross over.** So the prediction is that this world's closed basins are
mostly soda lakes -- trona and natron, alkaline, the Lake Magadi sort -- and that
gypsum and Ca-Cl brines occur only where the catchment already drains evaporite.
That is a testable structural claim, not a tuning knob, and it is checkable
against the salt-crust and playa-fill split the builds already carry.

It also feeds back into `duricrust`: gypcrete needs Ca and SO4 surviving to
saturation, which the alkaline path removes early. Gypcrete should therefore be
*rarer* than the climatic window alone implies, and concentrated in
evaporite-draining catchments.

**Computed per basin, weighted by discharge.** `pedology/scripts/brine_paths.py`
maps every Orogen rock class to a Meybeck lithology, weights the released
chemistry over each basin's catchment, and applies the divide. The weight is
LOCAL runoff generation, P - E clamped at zero times region area, not the
accumulated river discharge: accumulating would count every upstream region again
at every downstream one and let a handful of cells near the sink decide a basin.

Weighting by area instead was the first pass and it is retained under
`--weighting area`, because the two are not a refinement apart. Runoff on this
world spans four orders of magnitude inside one large catchment, so a dry interior
contributes area and almost no solute while a wet rim contributes the chemistry.
Measured 2026-08-17: 109 basins change side between the two weightings, and 655
basins have a catchment but generate no runoff at all under this climatology.
Those are reported UNDETERMINED rather than falling back to the area answer,
because a basin that receives no water has no brine to evolve.

The split stays lopsided either way. Under discharge weighting 2,757 of the 2,966
resolved basins are alkaline, holding 94.6% of resolved catchment area, with
Ca/HCO3 from 0.305 to 1.423 and a median of 0.766.

The check that can fail here is neither of those numbers. It is that a basin
draining exactly one rock class must return that rock's own table row, because the
weighted mean of a constant is that constant whatever the weights are. 108 basins
qualify and the worst relative error is 1.2e-15. Comparing the two weightings
against each other could only ever have told us they differ.

The 20-to-10 rock class mapping is a judgment and is written out in full rather
than defaulted, since a silent default would push every unmatched class to one
side. The largest single judgment -- whether `shelf_clastic` behaves as shale or
as sandstone -- was tested with `--clastic-as-sandstone` and moves the alkaline
area share by 0.6 percentage points, so the conclusion does not rest on it.

Two structural findings fall out. **Orogen's rock table has no gypsum
lithology**: its evaporite class is a salt crust, so a Ca-rich path can only come
from halite or from a carbonate-poor silicate mix, which independently predicts
gypcrete to be rare. And the divide is decided by *catchment* lithology, not by
what is on the basin floor, so a basin's evaporite mineralogy is set by rock it
does not itself sit on.

Two cautions. The divide is stated on total equivalents, and Mg complicates it
through Mg-calcite and dolomite, which this first pass ignores. And these ratios
are Meybeck's temperate-stream release values, not this world's -- they set the
ordering of lithologies, which is what matters here, rather than absolute
concentrations.

## What the classifier produced

Measured 2026-08-17 on the pre-carve build, at the central end of the aeolian
roughness bracket. Every one of these is regenerated after the next carve; what
survives is the classifier.

| `surface_cover` | % of land | | `duricrust` | % of land |
| --- | ---: | --- | --- | ---: |
| water | 12.59 | | none | 59.09 |
| loess | 7.58 | | gypcrete | 0.70 |
| diatomite | 4.50 | | calcrete | 37.60 |
| pavement | 5.67 | | silcrete | 2.62 |
| bare | 12.81 | | | |
| soil | 56.87 | | | |

**Gypcrete is rare for the predicted reason, and the size of the effect is the
result.** 7.30% of land meets Watson's climatic window -- under 250 mm/yr with
potential evaporation exceeding precipitation in every one of the year's twelve
bins -- against about 22% of Earth's land within the same bounds. The
carbonate-poor ion gate then cuts that to 0.70%, a factor of 10.5, because the
alkaline path removes Ca at the first calcite and only 5.2% of endorheic land
drains a carbonate-poor catchment. Gypcrete area is 19.1 times enriched in
carbonate-poor catchments over their share of the land. That was a prediction
before it was a measurement.

**Calcrete's ion gate does not bind at all**, which is also what was predicted:
every silicate and carbonate lithology in Meybeck's table supplies Ca and HCO3,
so calcrete here is climate-limited and nothing else. The gate is kept because a
lithology that failed it would have to be visible rather than silently absent.

**The silica gate does not bind either.** Every pan margin and drainage line in an
endorheic catchment clears 100 umol/l. Silcrete on this world is limited by
setting, not by supply.

**Loess spans a factor of six across the aeolian roughness bracket and the
bracket is the answer**, not the central value:

| aeolian z0 | land-mean deposition, g/m2/yr | loess, % of land | pavement, % of land |
| --- | ---: | ---: | ---: |
| 3e-6 m, the smooth end | 363.6 | 46.64 | 0.47 |
| 1e-4 m, central | 26.9 | 7.71 | 5.91 |
| 1e-3 m, the rough end | 0.002 | 0.00 | 9.57 |

Loess and pavement move in opposite directions across it, which is the continuum
falling out rather than being imposed: they are one threshold read from two sides.
The smooth end buries the clast mosaics and the rough end leaves them everywhere.

**Which silcrete setting fired is recorded, and that answers the open question
about them.** Pan margin covers 4.64% of land and drainage line 0.75%. The two are
distinguishable as SETTINGS and are kept apart in a bitmask; their products are
not distinguishable at 15.19 km, so they share one `silcrete` value. Groundwater
silcrete is left unplaced, because it sits at or near a water table and this
project models none; a proxy built out of surface elevation would be inventing
the mechanism rather than approximating it.

Two costs of the declared precedence, recorded because the ordering is a judgment.
2.73 percentage points of land carry the silcrete setting but are reported as
calcrete, which takes precedence as the narrower and two-sided window; the
bitmask recovers them. And `diatomite` takes 3.71 percentage points of barren
land ahead of `bare`, which is correct rather than a leak -- a closed basin's
strandline IS playa, and calling it diatomite is what the class is for.

## What is undetermined

- **The silica source for silcrete.** RESOLVED, and by routing rather than by
  hypothesis: 68.4% of dissolved silica reaches a closed basin, and the gate that
  tests for it never binds. Diatomaceous basin fill is a consumer of that silica
  rather than its origin.
- **Which non-pedogenic silcrete types are distinguishable here.** RESOLVED as
  above: the settings separate, the products do not, so one value plus a setting
  bitmask.
- **Whether pavement and loess can coexist.** RESOLVED as a continuum, which is
  what the design suspected. They are one deposition threshold read from two
  sides, and the bracket table above is what that looks like.
- **Diatomite over the stellar cycle.** Still open, and it is the largest
  remaining gap in this note. What is implemented is the STATIC proxy: the
  strandline band a single climatology's lake could vacate, between the solved
  lake surface and the spill level, for basins that are neither dry nor pinned at
  their spill. What the design asks for is occupancy measured over the 57-year
  component, which is the forcing slow enough for lake level to equilibrate. That
  needs the cycle run, CYC-1, and cannot be faked from one climatology.
- **Ion sources for the duricrusts.** RESOLVED. Meybeck (1987) Tables 2C and 5
  are extracted and validated in
  `pedology/data/reference/meybeck1987_tables.json`. Gypcrete needs Ca and SO4,
  calcrete needs Ca and HCO3, and Table 2C gives both by rock type.

  The apparent problem that Table 5 is **exorheic-only** dissolves on reading
  Meybeck's section 6: he assumes rock proportions are identical in the exorheic
  and endorheic areas, and excludes endorheic only because those *waters* are
  altered by evaporation and precipitation. Table 2C, the release, was never
  exorheic-specific; only Table 5's export shares are, and a closed basin has no
  export to share out.

  What a closed basin needs instead is a **fate** model, and that composes with
  Table 2C exactly. See the section below.
- **Pavement clast-supply criterion.** Still undetermined. Coarse substrate is
  necessary and the read sources give no threshold, so the classifier stands in
  "a consolidated bedrock substrate rather than basin fill" and declares it as a
  stand-in. That is the weakest rule in the file.
