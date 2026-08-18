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
consumer of fields nothing currently reads -- `hotspot`, `r_tectonicActivity`
and `r_t_craton` are all diagnostics with no downstream reader today.

| deposit type | control | field it keys on |
| --- | --- | --- |
| Porphyry Cu-Mo-Au | continental arc, upper crust | arc belt, `granodiorite` root |
| Orogenic gold | fold belt, greenschist grade, shear zones | `r_t_foldBelt`, `r_stress` |
| VMS | submarine volcanic, ridge and back-arc | `morb`, `backArcDist` |
| Magmatic Ni-Cu-PGE | LIP feeder systems | `flood_basalt`, `lipV` |
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
| Potash, lithium, borate, soda | closed-basin brine evolution | `pedology/scripts/brine_paths.py` |

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

DECLARED, what remains: potash, and **lithium and borate entirely -- neither
element is among Meybeck's eight species**, so unlike soda and gypsum they cannot
be derived from the divide in any form, and what stands in is a geological
association with silicic volcanic and arc volcanic catchments. That is an
association, not a mechanism this pipeline models.

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

## What this world cannot have

Worth stating so nobody looks for them.

**Banded iron formation.** Requires a specific atmospheric and ocean redox
history over deep time. Orogen has no time axis at all, and this project models
one atmospheric state, not a history.

**Anything keyed on absolute age.** Same reason. Deposit types whose defining
control is "Archean" or "Proterozoic" cannot be placed, only their tectonic
setting can. The concrete casualty is komatiite-hosted Ni: Naldrett (2010) splits
magmatic sulphide deposits into a komatiite-related class and a flood-basalt
class, and puts the komatiite one at 2.7 to 1.9 Ga. Only the flood-basalt class
is available here, and `magmatic_nicu` keys on it alone.

**Rift-hosted deposits, on the current geography.** `rift_bimodal` is 0% of land
here, because continental rifting needs two adjacent continental super-plates and
this world has none. That is geography rather than a gap, and it removes
carbonatite-hosted rare earths along with it.

## Prospectivity, not deposits

**This layer emits a continuous 0-1 favourability field. It does not place
individual deposits.**

The reason is resolution, and it is the same reason glacial erosion is deferred
in `notes/audits/orogen-gravity.md`. A mesh cell is 15.19 km across. A porphyry
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
