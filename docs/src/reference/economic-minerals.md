# Economic minerals: where each kind belongs

Design, and both halves are now built: `minerals/scripts/build_prospectivity.py`
for the tectonic and magmatic types, `build_downstream_prospectivity.py` for the
weathering, drainage and brine types. Recorded because the decision was taken in
discussion and then not written down, and because the architectural constraint
below is the part that would have been expensive to retrofit.

## The constraint that matters

**Minerals are a separate overlay field. They are never a `substrate_class`
value and never an erodibility modifier.**

`substrate_class` is load-bearing in four directions: it sets erodibility, which
sets erosion, which sets terrain; it sets bare-rock albedo, which ExoPlaSim
integrates; it sets soil texture and pH through pedology; and it maps to Meybeck
solutes for the weathering fluxes. Anything entering it enters the climate path.

An ore body is the wrong scale for that. A mesh cell is 15.19 km across; a
porphyry system is one to a few kilometres, a vein is metres. Making "copper" a
rock class would paint a whole cell with ore properties on the strength of a
deposit occupying a fraction of a percent of it, and the planet's energy balance
would move because of a mine. The same argument rules out nudging erodibility:
a change there is a change in terrain.

So the layer is **prospectivity per cell**, orthogonal to lithology, read only by
what comes after climate -- resources, settlement, culture. Nothing upstream of
it may consume it. This also means it can be rebuilt freely without invalidating
a single climate run, which is the practical benefit of keeping it out.

## The split, by genesis

The rule is simple: **a deposit belongs where the process that concentrates it is
modelled.** Orogen has tectonics, magmatism, metamorphic grade and erosion, and
no climate. Everything downstream has climate, weathering, drainage and brine
chemistry.

### Orogen: tectonic and magmatic genesis

These need only what the generator already computes, and several are the first
consumer of fields nothing else reads -- `r_t_craton` is now consumed by
`build_prospectivity.py`, while `hotspot` and `r_tectonicActivity` remain
diagnostics with no downstream reader.

| deposit type | control | field it keys on |
| --- | --- | --- |
| Porphyry Cu-Mo-Au | continental arc, upper crust | arc belt, `granodiorite` root |
| Orogenic gold | fold belt, greenschist grade, shear zones | `r_t_foldBelt`, `r_stress` |
| VMS | submarine volcanic, ridge and back-arc | `melange`, `arc_basalt`, `morb` |
| Magmatic Ni-Cu-PGE | LIP feeder systems | `flood_basalt` (`lipV` expressed through it) |
| Podiform chromite | obducted ocean crust | `melange`, forearc |
| Kimberlite / diamond | thick cold cratonic keel | `r_t_craton` |

Kimberlite is the one worth noting: it needs a craton, and until the basin
double-count was removed this world had almost none. The craton fix is what makes
diamonds placeable at all.

### Downstream: weathering, drainage and brine genesis

Cheaper than it looks, because the machinery mostly exists.

| deposit type | control | where |
| --- | --- | --- |
| Bauxite | intense tropical weathering | pedology, weathering intensity `W` |
| Laterite Ni | same, over ultramafic parent | pedology, `W` + parent class |
| Supergene Cu enrichment | weathering over a porphyry | pedology `W` over an Orogen porphyry |
| Placer Au, Sn, gem | fluvial concentration | hydrography, drainage and discharge |
| Potash, lithium-borate, soda | closed-basin brine evolution | `pedology/scripts/brine_paths.py` |

`brine_paths.py` already solves the chemical divide per basin and classifies each
as alkaline or Ca-rich. Which evaporite mineral a basin grows follows from that
path, so the brine half of the list was close to free.

Supergene enrichment is the clearest case for why the split is a split and not a
handover: it needs an Orogen deposit AND a downstream climate, and neither
component can produce it alone. The implementation makes that structural rather
than remembered -- supergene copper MULTIPLIES the tectonic file's own porphyry
field and scores zero where there is no protore, so a favourable climate on its
own can never place copper.

**The two halves are separate artifacts, and the reason is lifetime.** The
tectonic field is a pure function of the export and survives a re-run of the
baseline. The downstream field reads a climatology, a lake solution and the brine
solve, so it carries a climate in its identity as well as a terrain and is
regenerated when either moves. One file would have given the durable half the
disposable half's lifetime.

**The rules are not equally grounded and the artifact says so per rule.** Each
carries `derived`, `sourced`, `sourced-negative` or `declared`, in the config, in
the report, and as an attribute on the netCDF variable, so the distinction
travels with the number.

DERIVED: soda ash and gypsum are the two sides of Hardie and Eugster's chemical
divide, which decides irreversibly which way a brine goes, and there is no free
choice in them.

SOURCED: bauxite, nickel laterite and supergene copper's dry end, on literature
read 2026-08-17 and listed in `references/INDEX.md`. **One of the three moved by
an order of magnitude when it was sourced.** Supergene copper carried a 100 mm/yr
lower bound described as the conventional semi-arid band and chosen without a
source; Reich et al. (2009) measured the Atacama and put meteoric enrichment
above 10 mm/yr with shutdown below 1-4, so the invented bound was excluding the
best-documented enrichment province on Earth by a factor of ten.

The two Economic Geology 100th Anniversary Volume chapters that Sillitoe (2010)
delegates to arrived 2026-08-18, and between them they closed the rest. What they
settled is below.

SOURCED-NEGATIVE, which is the interesting label. The placer rule has no gradient
or discharge threshold and no transport-distance decay, and BOTH absences are
quoted rather than confessed. Slingerland and Smith (1986) say at p. 143 that the
regional criteria a rule like this wants were never established; Knight et al.
(1999) show gold is progressively flattened rather than lost with distance, so a
decay length would remove prospectivity the evidence says is still there. A
number in either place would have been invented and then quoted back as sourced.

DECLARED, what remains: **potash alone.** Lithium and borate moved to
`sourced-negative` on 2026-08-18, and became a single field in the process; the
section below is why.

### The relief term, and the supergene wet end that does not exist

Settled 2026-08-18 from Freyssinet, Butt, Morris and Piantone (2005) and Sillitoe
(2005), read off the page images because the volume's OCR is unreliable. Three
gaps that had been recorded separately turned out to be one gap.

**There is no wet bound on supergene copper.** Sillitoe p. 736 puts every climate
inside the window except hyperarid desert and glacial or permafrost ground, and
calls a tropical climate the most favourable of all for rapid upgrading. The rule
had carried a 500 mm/yr cap, declared, on the reasoning that copper is flushed
once percolation runs year round. That cap was not merely unsourced; the canonical
source contradicts it, so it is removed rather than re-derived.

**The wet-side control is erosion, and that is what unifies the three gaps.**
p. 736: the average erosion rate "must be in overall balance with the rate of
water table descent", and "Erosional efficiency depends chiefly on rainfall and
slope steepness". Rainfall reaches the wet end only jointly with slope, so a
rainfall cap was the wrong SHAPE of term, not a term with the wrong number. The
same balance is Freyssinet's criterion 2 for laterite, p. 695: the weathering
front must descend faster than the surface is lowered.

**The relief term is two-sided**, which is the correction that matters, because
both rules previously recorded the omission as one-signed. Freyssinet p. 695:
"this is favored by low relief. Conversely, relief must be sufficient to allow
leaching of chemical weathering products." Too steep strips the profile; too flat
impedes the drainage that leaching needs. Table 2 on the same page splits the ore
by exactly that: hydrous Mg silicate wants moderate relief and free drainage, clay
silicate moderate-to-low and impeded, oxide moderate-to-low and either.

**What is applied is an upper bound of 5 degrees of regional dip**, from p. 685,
where bauxite-bearing plateaus are remnants of planation surfaces "with overall
dips of 1 deg to 5 deg". The lower end of that pair is NOT applied as a floor: it
describes where those particular surfaces sat, not a limit below which bauxite
fails, and the mechanism a floor would stand for is already carried
hydrologically by `require_exorheic`. The same 5 degrees is applied to nickel
laterite on Freyssinet's own framing at p. 695, that his three constraints are
"common conditions that apply to all deep lateritic regoliths".

**Measured on `precarve-craton` 2026-08-18**, area-weighted over land, the
regional dip runs 0.12 deg at the median and 4.78 at the 99th percentile, so the
5 degree bound removes only 0.86% of land. That reads like a weak criterion and
it is not, because the deposits are not spread over land uniformly: it removes
**5.13% of melange**, the forearc province that hosts nickel laterite, six times
the rate it removes land at.

Measured on the deposits themselves, by running the generator twice at the same
commit with only the bound changed, it costs **bauxite 3.9% of its footprint and
nickel laterite 7.5%** -- nearly twice as much, which is the melange result
showing through. Every other field is bit-identical across that pair, which is
what makes the attribution good. Do not compare against reports generated before
2026-08-18: `lib/gridding.py` and `hydrography/scripts/surface_water.py` both
moved in between, and the placer field shifts across that boundary for reasons
that have nothing to do with this change.

**Removing the supergene wet cap is the large effect of the two**, and it is not
close: enrichment goes from 1.86% of land to 7.41%, four times the footprint. An
invented bound had been suppressing three quarters of the deposit.

**What is still not applied**, and its sign. The criterion is a rate against a
rate and this project models neither erosion rate nor weathering-front descent.
The anchors are in the config for whoever closes it. No slope bound is applied to
supergene copper: unlike the laterites it is not tied to planation surfaces --
Sillitoe's abstract says pediplains "are not considered to be a requirement" --
and deep profiles do form in mountainous tropics on masses isolated from drainage
incision. So on wet steep ground the supergene rule OVER-predicts, scoring
prospectivity where the real profile would be stripped before it matured.

**Tin and gem placers are not emitted**, alongside the exclusions below and for a
resolution-adjacent reason: cassiterite needs specialised S-type granite and gem
placers need their own host suites, and Orogen's 20-class table separates neither
from ordinary granite. Placing them would be placing granite twice under
different names.

### Lithium and borate: one field, and why the hydrothermal objection failed

Settled 2026-08-18 from Risacher and Fritz (2009), Risacher et al. (2003), Munk
et al. (2016), Helvaci (2019) and Huh et al. (1998), all read. Both rules had
been carried as bare geological associations on the reasoning that neither
element is among Meybeck's eight species and that boron in Earth's borate basins
is substantially hydrothermal, which this pipeline does not model. The first half
is still true. The second half is wrong for the deposit type this world can
actually hold, and it was wrong on a point that has been measured.

**The test that could have failed.** Lithium and boron are generally enriched in
thermal waters (Berthold and Baker 1976; White et al. 1976), so the expected
result was enrichment. Risacher and Fritz (2009) plot Li and B in Bolivian and
Chilean dilute springs and groundwaters against water temperature, over a sample
in which 60% of Bolivian and 78% of Chilean waters carry a thermal influence, and
find no rise with temperature and no significant difference between meteoric and
hydrothermal waters. Their section 3.5: *"hydrothermal alteration is not
specifically responsible for the high content of Li and B in waters and brines of
Andean salars. Both components originate from the alteration of volcanic rocks,
especially ignimbrites because of their rapid alteration and weathering."*
Risacher et al. (2003) section 7.7 reach the same conclusion from the Chilean
half alone, with the same test.

So the flux this pipeline cannot place is not required for the Andean type, and
the Andean type is the world's principal lithium brine province and the setting
of its Quaternary borates. What DOES require it is the Anatolian borate type:
Helvaci (2019) lists thermal springs near volcanic activity among the essential
conditions, and attributes the high boron of the Turkish source rocks to
circulating boron-rich hot water. That type is not placeable here, and naming it
as unplaceable is the honest form, the same shape as groundwater silcrete being
listed and left empty in `pedology/config/surface_classes.yaml`.

**One host set now serves both elements, and the borate rule changed.** Floyd et
al. (1998), quoted by Helvaci (2019), characterise the fertile ignimbrites of the
Turkish borate basins as *"well-evolved and fractionated ... with a high-silica
rhyolitic bulk composition, exhibit a combined high content of B, As, F, Li and
Pb"*. Helvaci puts the boron source at "andesitic to rhyolitic volcanics" and
reports 25 to 270 ppm B in the calc-alkaline andesite-to-rhyolite suite around
the deposits. The borate rule had read `{arc_andesite: 1.0, arc_basalt: 1.0}`,
which sourced boron from island-arc basalt: the wrong end of the differentiation
series for an element that concentrates in evolved melts. It now takes lithium's
host set. Plutonic classes stay at half weight because the leachable phase is
glass rather than crystal, which Hofstra et al. (2013) via Munk et al. (2016)
show for Li directly and which Wolff-Boenisch et al. (2004), already read for the
andisol work, quantifies as a dissolution-rate contrast.

**The source gate is necessary, not decorative.** Munk et al. (2016) supply the
negative case: the closed basins of the Yilgarn craton in arid Australia have the
climate and the hydrology, have held closed drainage since the Eocene, reach five
to ten times seawater salinity, and carry almost no lithium, with most brines
below 0.5 mg/L. *"Significant Li brine accumulations have not been reported from
intracratonic basins in arid regions."* A cratonic catchment with the right
climate does not make a lithium brine, which is the statement that makes the
lithology term load-bearing rather than ornamental.

**What remains sourced-negative is the release coefficient, and the canonical
review says so rather than merely omitting it.** Munk et al. (2016) on lithium
sources: *"Despite the observation that multiple potential sources of Li exist
they are yet to be definitively identified and quantified."* The direct
measurements straddle orders of magnitude on the same rocks: Price et al. (2000)
conclude that groundwater leaching of volcanic tuffs alone accounts for all the
lithium at Clayton Valley, while Jochens and Munk (2011) measured under 10 ug/L
Li released from those same rocks at ambient conditions, and Godfrey et al.
(2013) report the same at Salar del Hombre Muerto.

**The search for a release table, so nobody repeats it.** Huh et al. (1998) is
the nearest thing to a Meybeck Table 2C for lithium and is not one: it samples
the lower reaches of thirteen major rivers plus tributaries, groups them by
TERRAIN rather than by rock class, and gives shield rivers 27 to 70 nM,
tectonically active zones 160 to 580 nM, and rivers draining marine limestone and
evaporite 450 to 11,700 nM, against a flow-weighted mean of 215 nM. Its high end
is recycled evaporite lithium rather than release from fresh rock, which is why
it ranks terrains in the opposite order to the brine literature; the two are
answering different questions and only the brine question is this project's.

The compilation chapter states the absence outright. Gaillardet, Viers and Dupre
(2014), *Trace Elements in River Waters*, is where a per-lithology trace-element
table would live if one existed, and its own conclusion is that "trace-element
concentrations have been measured in a restricted number of environments, which
may not be typical", so that "extended investigations on rivers ... of variable
climatic or lithologic settings are necessary to get a better world-scale
overview of trace-element levels and controlling parameters". That is the same
shape as Slingerland and Smith on placers: the canonical review says the quantity
was never established.

The rest of the search was empty in the same direction. The boron river
literature is isotopic almost throughout, aimed at weathering regime and
palaeo-pH rather than at yield per rock type. Bluth and Kump (1994), the one
monolithologic-catchment compilation this project holds, carries bicarbonate and
silica only. A number entered in place of the missing table would be invented and
then quoted back as sourced.

**Considered and rejected: requiring the alkaline brine path for borate.**
Helvaci gives lake water at pH 8.5 to 11 among the essential conditions, which
maps onto this project's alkaline path, and the gate would have been cheap.
It would also delete the Andean borates, which are mined from salars that
Risacher et al. (2003) classify as sulfate-rich and calcium-rich rather than
alkaline. Two deposit types share one class here, so the narrower gate is wrong
for half of what the class covers.

**The correction collapsed the two rules into one, and that is the right
outcome.** Replacing borate's arc host set with lithium's -- dropping
`arc_basalt`, adding `rift_bimodal` and the half-weight plutonics -- left
`brine_lithium` and `brine_borate` identical in every key but `name`: same
host set, same `require_no_overflow`, same grounding. **Measured on `precarve-craton`
2026-08-18**, regenerating produced bit-identical fields, both at 9.62% of
land nonzero and mean 0.0145. Borate's previous report row, on its old hosts,
read 9.62% of land nonzero, 3.72% above 0.5, mean 0.0033 -- the report's
three columns: the same nonzero support, with the weights and therefore the
mean moving. Two fields that are identical by construction carry no information the
one carries, and a reader seeing two names infers two predictions. That is the
error this note already names for tin and gem placers, where the response was to
not emit them: "Placing them would be placing granite twice under different
names." So the two are emitted as one, `brine_lithium_borate`.

The precedent inside this project is exact and it is silcrete. Three settings,
one output value, "so the distinction survives without asserting that this
resolution can separate the products". Same shape: the geology distinguishes,
this pipeline does not, and the honest form is one value.

Naming it for both elements is not a compromise, because the association is what
the sources assert. Risacher and Fritz treat them in one section, "Lithium and
boron are conspicuous components of waters and brines of Andean salars", and
Bradley et al. (2013) record that lithium brines are "commonly associated with
borate mineralization in arid, closed basins". Co-occurrence is the observation.

**Nature does separate them, the direction is known, and the boundary is not
placeable here.** This is worth stating precisely rather than as a shrug. Bradley
et al. (2013): lithium "does not readily produce evaporite minerals when
concentrated by evaporation. Instead it ends up in residual brines in the shallow
subsurface", whereas borate leaves solution as a mineral. Borate therefore
saturates EARLIER in the evaporative sequence than lithium reaches an economic
brine concentration, which means the borate set CONTAINS the lithium set and the
two differ only in where the boundary between them falls. Placing that boundary
needs a per-basin Li and B concentration to test against a borate saturation
point, and a per-basin concentration needs the release coefficient this section
has just established does not exist. It is the same missing number surfacing a
second time, which is itself the useful finding: the gap is one gap, not two.

**What would split the field again**, so nobody re-derives this. Two things, and
the machinery for the second is already built.

1. A per-lithology Li and B release table. That would give inflow concentrations
   an evaporation model could carry to borate saturation -- Risacher's own
   EQL/EVP is exactly such a code -- and borate would come out the LARGER field,
   which is the prediction this argument makes and the way to check it.
2. A source giving the evaporative concentration factor at which a borate mineral
   saturates against the factor an economic lithium brine needs. Note what this
   would NOT need: hydrography already computes a graded `aridity_index` against
   `critical_aridity_index` per basin, and `require_no_overflow` throws that
   gradation away by reading only the boolean `basin_fills_to_spill`. The graded
   quantity exists and only the threshold is missing. Do not conclude from this
   rule that the pipeline lacks the field; it lacks the number.

## What this world cannot have

Worth stating so nobody looks for them.

**Banded iron formation.** Requires a specific atmospheric and ocean redox
history over deep time. Orogen has no time axis at all, and this project models
one atmospheric state, not a history.

**Anything keyed on absolute age.** Same reason. Deposit types whose defining
control is "Archean" or "Proterozoic" cannot be placed, only their tectonic
setting can. The concrete loss is komatiite-hosted Ni: Naldrett (2010) splits
magmatic sulphide deposits into a komatiite-related class and a flood-basalt
class, and puts the komatiite one at 2.7 to 1.9 Ga. Only the flood-basalt class
is available here, and `magmatic_nicu` keys on `flood_basalt` with a small
`oib` weight beside it.

**Rift-hosted deposits, on the current geography.** `rift_bimodal` is 0% of land
here, because continental rifting needs two adjacent continental super-plates and
this world has none. That is geography rather than a gap, and it removes
carbonatite-hosted rare earths along with it.

## Prospectivity, not deposits

**This layer emits a continuous 0-1 favourability field. It does not place
individual deposits.**

The reason is resolution, and it is the same reason glacial erosion is deferred
in `notes/audits/orogen-gravity.md`. A mesh cell is about 15 km across. A porphyry
system is one to two kilometres, 0.07 to 0.13 of a cell; a vein or lode is tens
to hundreds of metres, under a hundredth of one. A discrete deposit is invisible
at every resolution this pipeline currently runs at, so placing one here would be
inventing detail the grid cannot hold and then carrying it as though it were
resolved.

A field is what the model actually knows: this cell has the setting, the host and
the structural control, so it is favourable. That composes cleanly -- downstream
work can threshold it, sample it, or read it as a gradient -- and it stays honest
about the fact that favourability is all the evidence supports.

Discrete deposits belong with the downscaling pass, alongside glacial
overdeepening and the sub-grid hypsometry that `notes/glacier-rough-pass.md`
already assigns there. That pass will have the resolution to place a deposit
somewhere meaningful, and the prospectivity field is exactly the input it wants.

Two sub-grid features, resolved independently, landing in the same place is the
useful signal: **if a thing is smaller than a cell, it is not this pipeline's to
place.**
