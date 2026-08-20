# Minerals

Ore prospectivity, as a 0-1 field per deposit type. The design and the split
between what belongs here and what belongs downstream are in
`../docs/src/reference/economic-minerals.md`; this file is what the component does and how to
run it.

```bash
python minerals/scripts/build_prospectivity.py             # tectonic and magmatic
python minerals/scripts/build_downstream_prospectivity.py  # weathering, drainage, brine
```

Writes `data/<source_build>/prospectivity.nc` and
`data/<source_build>/downstream_prospectivity.nc` on the native mesh, each with a
report beside it.

**Two artifacts, because they have different lifetimes.** The tectonic and
magmatic types are a pure function of the export and survive a re-run of the
baseline. The weathering, drainage and brine types read a climatology, a lake
solution and `pedology/scripts/brine_paths.py`, so they carry a CLIMATE in their
identity as well as a terrain and are regenerated when either moves. Keeping them
in one file would have given the durable half the disposable half's lifetime.

The downstream script needs `build_prospectivity.py` and `brine_paths.py` to have
run first, and refuses rather than falling back. Supergene copper MULTIPLIES the
tectonic porphyry field, so a climate window on its own would place copper
wherever the weather suited, which is exactly the error the genesis split exists
to prevent.

## Two constraints, and they are the point

**This is a field, not deposits.** A porphyry system is one to two kilometres
against a 15.19 km mesh cell; a vein is under a hundredth of one. An individual
body is invisible at every resolution this pipeline runs at, so placing one would
invent detail the grid cannot hold. Discrete deposits belong with the downscaling
pass -- the same place glacial overdeepening goes, for the same reason. This
field is that pass's input.

**This is never a lithology and never an erodibility modifier.**
`substrate_class` sets erodibility and therefore terrain, albedo and therefore
climate, texture and therefore the biosphere, and solutes and therefore the
weathering fluxes. Anything entering it enters the climate path, and a whole cell
would take on ore properties because of a deposit occupying a fraction of a
percent of it. Prospectivity is read only by what comes after climate.

That constraint is why this is computed here rather than inside Orogen. Every
input but one is already exported -- craton weight, fold-belt weight, stress,
substrate and basement class, cover thickness, erosion delta -- and the missing
`lipV` has its expression in the `flood_basalt` class. A value that does not
exist during generation cannot feed back into erosion, so the rule is structural
rather than a thing to remember. It also means a rule can be changed without
rebuilding the terrain.

## What the rules key on

They live in `config/prospectivity.yaml`, not in the script, so disagreeing with
one is a config edit. Each names the control it encodes, and they are first
order.

The axis worth understanding is **exhumation**, because it is what a rock map
alone cannot tell you. A porphyry ore zone is shallow -- Sillitoe (2010) puts the
parental plutons at 5-15 km with the mineralised stocks above them and the
epithermal counterparts under 1 km -- so deep erosion removes it. Orogenic gold
is the opposite, and Groves et al. (1998) sharpen what that means: it is a
crustal continuum from 15-20 km to the surface, epizonal under 6 km, mesozonal
6-12, hypozonal beyond 12. Exhumation therefore does not decide whether gold is
present, it decides which class is exposed. The same erosion that removes a
porphyry brings a mesozonal gold system to surface. Orogen tracks it through
`cover_thickness` and `erosionDelta`, which is why porphyry prospectivity here
rejects most of the arc belt rather than accepting all of it -- 59.5% of it is
stripped to its granodiorite root, and the porphyry level went with the section
above.

**Porphyry is split by arc type**, and the split is what the literature supports
rather than a weighting. Cooke et al. (2005) show the two metal associations have
different distributions: giant Cu-Mo clusters where continental crust is
thickened, while the largest gold-rich porphyries concentrate in the southwest
Pacific island arcs. Orogen already carries that distinction, since
`arc_andesite` is continental-arc and `arc_basalt` is island-arc.

On this world the split has a clear answer. Continental arc outweighs island arc
about forty to one by area, so Cu-Mo reaches high favourability across 4.65% of
land while gold-rich porphyry clears 0.5 on 0.12%. This is a copper planet rather
than a gold one, and that follows from its arcs being continental.

Kimberlite is the other entry worth noting: it needs a thick cratonic keel, and
before the craton basin double-count was removed this world had almost none. See
`../notes/audits/orogen-lithology.md`.

## What the downstream half keys on

The brine entries are the ones worth understanding, because two of them come out
of a solved chemical divide rather than out of a weighting.

**Soda ash and gypsum are DERIVED.** Hardie and Eugster's chemical divide says
calcite removes Ca and CO3 in equal equivalents, so whichever is in excess at
calcite saturation dominates every later step and the deficient one is driven
toward zero. The first calcite therefore decides irreversibly which way a brine
goes, `brine_paths.py` solves it per basin from Meybeck's release table, and
these two fields are the two sides of that answer. They are scaled by the SIZE of
the excess rather than by its sign, because a basin at Ca/HCO3 of 0.31 is a soda
lake in a way one at 0.98 is not.

Both are gated on the basin not overflowing. A basin pinned at its spill never
reaches saturation, so nothing crystallises whatever its chemistry, and the lake
solver already computed that boolean.

**The weathering and drainage rules are sourced, and sourcing moved an
invented bound by an order of magnitude.** Supergene copper's lower bound is
Reich et al. (2009): Atacama meteoric enrichment runs above 10 mm/yr and shuts
down below 1-4. There is no wet cap: Sillitoe (2005) p. 736 admits every
climate except hyperarid desert and glacial or permafrost ground, so the
exclusion that stands is the one the source names, tested on seasonal melting
rather than an annual mean. The real wet-side control is erosion outpacing
water-table descent, and because erosion answers to rainfall and slope
TOGETHER, a rainfall cap is the wrong shape of term.

**A second paper in the same anniversary volume supplied the laterites'
missing relief term**, so the two papers between them closed all three gaps. Bauxite and nickel laterite now carry a 5 degree
upper bound on regional dip, from the planation surfaces Freyssinet et al. (2005)
p. 685 measure bauxite on; the field is `lib/orogen.py:local_slope_deg`, a plane
fit through each region and its neighbours rather than a steepest drop to one,
because it is a regional TILT that the source is written about. The omission it
closes was recorded here as one-signed and is in fact two-sided -- too flat
impedes the leaching as surely as too steep strips the profile. See
`docs/src/reference/economic-minerals.md` for what is still not applied and which way it errs.

Bauxite is Price et al. (1997), whose criteria are the right shape because they
are themselves thresholds applied to gridded climate fields and validated against
observed bauxite: above 1200 mm/yr, mean annual temperature above 22 C, and 6 or
fewer months below 60 mm. The monthly half carries the mechanism, an oscillating
water table through a short dry season, so an annual mean would lose it -- and it
needed the calendar conversion, because a Vesper bin is half an Earth month, so
the COUNT transfers unchanged and the DEPTH halves.

**The placer rule's missing thresholds are quoted, not confessed.** It has no
gradient or discharge term because Slingerland and Smith (1986) state at p. 143
that the regional criteria were never established, and no distance decay because
Knight et al. (1999) show gold is progressively flattened rather than lost, so a
decay length would remove prospectivity the evidence says is still there. That is
`grounding: sourced-negative`, and it is a better outcome than a number.

**Potash is declared; lithium and borate are sourced-negative, and the config
says so per rule.** Potassium is Eugster and Jones' behaviour type IV, removed
mid-range by exchange and sorption and surviving only where concentration runs
far, so potash needs both a high K supply and extreme evaporation, and no source
fixes either threshold.

Lithium and boron are NOT among Meybeck's eight species at all, so neither can be
derived from the divide in any form, and no per-lithology release table exists
for either element to derive them from instead. What the literature does settle,
and settles for both at once, is which rock they come from: Risacher and Fritz
(2009) find Li and B in Andean salars uncorrelated with water temperature across
groundwaters that are mostly thermally influenced, and conclude both come from
ordinary alteration of volcanic rocks and especially ignimbrites. The
hydrothermal flux this pipeline cannot place therefore is not required for that
deposit type. What stays missing is the release coefficient, and Munk et al.
(2016) say plainly that the Li sources "are yet to be definitively identified and
quantified", so a number there would be invented.

**And they are ONE field, `brine_lithium_borate`, not two.** Floyd et al. (1998)
put B and Li in the same fertile ignimbrites, so the two rules take the same host
set; closure and aridity were already shared; and neither element is on the
divide. That leaves nothing this pipeline can apply that separates them, and two
names off one set of criteria is the tin-and-gem-placer error committed on
purpose. Nature does separate them and the direction is known -- borate leaves
solution as a mineral and lithium stays in the residual brine, so borate
saturates earlier and its set CONTAINS lithium's -- but placing that boundary
needs a per-basin concentration, which needs the release coefficient that does
not exist. `docs/src/reference/economic-minerals.md` carries the argument and what would split
the field again.

**Placer gold is the one rule the drainage network earns.** Orogenic gold
prospectivity is accumulated down the drainage tree exactly as runoff is, so a
river's score is the gold-bearing share of everything that drains into it. A rock
map cannot give that. What it lacks is a transport-distance term, because none of
the read sources gives one, and the consequence is one-signed: prospectivity is
carried too far downstream and the extent is an upper bound.

**Nickel laterite inherits an inference this repo has already made twice.**
Orogen's rock table has no peridotite class, and `melange` is the whole forearc
province including its serpentinite -- which is why `podiform_cr` keys on it as
obducted ocean crust. So the parent is `melange`, and the field should be read as
"this province contains the parent" rather than "this cell is peridotite".

Tin and gem placers are deliberately NOT emitted. Cassiterite needs specialised
S-type granite and gem placers need their own host suites, and the 20-class rock
table separates neither from ordinary granite. Placing them would be placing
granite twice under different names.

## Reading the numbers

Prospectivity is **relative within this world**, normalised by the land maximum
per deposit type. A 1.0 is the most favourable cell here, not a grade or a
tonnage, and comparing the number for one deposit type against another says
nothing. What is comparable is the spatial pattern within a type.

Read `grounding` before quoting any downstream number. Each rule carries one of
`derived`, `sourced`, `sourced-negative` or `declared`, and they are not equally
good: two of the brine entries fall out of a solved divide with no free choice,
`sourced-negative` means the criterion is sourced and the literature states that
the threshold it would want was never established, and the rest are judgment with
the reasoning written out. That distinction is in the config, in
the report, and as an attribute on every variable in the netCDF, so it travels
with the number rather than staying in a file someone has to remember to open.
