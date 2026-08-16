# Derived surface classes

*Design, 2026-08-16. Nothing here is implemented yet.*

World Orogen gives primary lithology: what the rock IS. It does not give what the
surface has BECOME under this world's climate and drainage. That second thing is
what albedo, dust emission and the phosphorus cycle actually key on, and it is
derivable here rather than upstream, because it depends on precipitation,
evaporation and lake extent that Orogen never sees.

Four candidates were proposed. One survives unchanged, one was dropped, one had
its mechanism backwards, and one is blocked on machinery that does not exist yet.

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

**Loess is blocked, not undetermined.** It requires knowing where dust is emitted
and where it lands, which needs an emission scheme and a transport path. Neither
exists yet. It is listed here so the classifier has a slot for it, and it stays
empty until there is something to fill it with.

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
2. **`loess`** -- net aeolian deposition fast enough to bury clasts. BLOCKED.
   When it is unblocked, do not calibrate to Bettis et al. (2003): their mass
   accumulation rates reach 17,500 g/m2/yr and many exceed 1,500, but that is
   last-glacial mid-continent North America, the high end of the global range and
   explicitly not a global typical. Muhs (2013) is the broad compilation to
   calibrate against instead.
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
a blanket. And silicification needs a dissolved silica source; where that comes
from on this world is unresolved, though diatomaceous basin fill is an obvious
candidate worth checking rather than assuming.

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

**Computed per basin.** `pedology/scripts/brine_paths.py` maps every Orogen rock
class to a Meybeck lithology, weights the released chemistry over each basin's
catchment, and applies the divide. All 3,629 basins resolve, and the split is
lopsided: about 90% of them, holding roughly 93% of endorheic catchment area, sit
on the alkaline side, with Ca/HCO3 running from 0.41 to 1.39 and a median near
0.82. Current values are in `pedology/analysis/brine_paths.json`.

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

## What is undetermined

- **Loess**, entirely, pending dust emission and transport.
- **The silica source for silcrete.** Silicification needs dissolved silica and
  nothing here yet says where it comes from. Diatomaceous basin fill is the
  obvious candidate on this world and would tie silcrete to the diatomite class,
  but that is a hypothesis to test rather than a mechanism to assert.
- **Which non-pedogenic silcrete types are distinguishable here.** Groundwater,
  drainage-line and pan/lacustrine have different settings but the same output
  class, and it is not yet clear that this project's resolution can separate
  them. If it cannot, they should collapse to one `silcrete` value rather than
  being split on a distinction the data cannot support.
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
- **Pavement clast-supply criterion.** Coarse substrate is necessary but the
  threshold is not something the read sources give directly.
- **Whether pavement and loess can coexist.** Physically they are the same
  process at different accumulation rates, so the classifier should probably
  treat them as a continuum rather than two classes with a boundary. Deferred
  until there is a deposition rate to put on the axis.
