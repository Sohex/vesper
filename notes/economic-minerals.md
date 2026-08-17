# Economic minerals: where each kind belongs

Design, not yet built. Recorded because the decision was taken in discussion and
then not written down, and because the architectural constraint below is the part
that will be expensive to retrofit if it is got wrong.

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
path, so the brine half of the list is close to free.

Supergene enrichment is the clearest case for why the split is a split and not a
handover: it needs an Orogen deposit AND a downstream climate, and neither
component can produce it alone.

## What this world cannot have

Worth stating so nobody looks for them.

**Banded iron formation.** Requires a specific atmospheric and ocean redox
history over deep time. Orogen has no time axis at all, and this project models
one atmospheric state, not a history.

**Anything keyed on absolute age.** Same reason. Deposit types whose defining
control is "Archean" or "Proterozoic" cannot be placed, only their tectonic
setting can.

**Rift-hosted deposits, on the current geography.** `rift_bimodal` is 0% of land
here, because continental rifting needs two adjacent continental super-plates and
this world has none. That is geography rather than a gap, and it removes
carbonatite-hosted rare earths along with it.

## Open question

Prospectivity or occurrence? A continuous 0-1 favourability field is honest about
what the model knows and composes cleanly with downstream use. Discrete deposits
are what a culture actually finds and mines. The likely answer is both -- a
field, plus a sampling of it into point deposits with a declared seed -- but the
sampling rule is undecided and should not be invented in passing.
