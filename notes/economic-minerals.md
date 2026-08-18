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
best-documented enrichment province on Earth by a factor of ten. Note also that
Sillitoe (2010), which is held and read, does NOT contain this: it delegates the
whole subject at p. 5 to Sillitoe (2005), which sits in an Economic Geology
volume that could not be fetched at all.

SOURCED-NEGATIVE, which is the interesting label. The placer rule has no gradient
or discharge threshold and no transport-distance decay, and BOTH absences are
quoted rather than confessed. Slingerland and Smith (1986) say at p. 143 that the
regional criteria a rule like this wants were never established; Knight et al.
(1999) show gold is progressively flattened rather than lost with distance, so a
decay length would remove prospectivity the evidence says is still there. A
number in either place would have been invented and then quoted back as sourced.

DECLARED, what remains: the wet end of the supergene window, any relief term for
bauxite or nickel laterite, potash, and **lithium and borate entirely -- neither
element is among Meybeck's eight species**, so unlike soda and gypsum they cannot
be derived from the divide in any form, and what stands in is a geological
association with silicic volcanic and arc volcanic catchments. That is an
association, not a mechanism this pipeline models.

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
