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

**But the class carries an albedo, and that is doing mineralogy by proxy.**
`evaporite` is 0.50, a clean-halite value. The module's own docstring states the
number is worth about 1.7 W/m2 per 0.10 of error on this planet. Gypsum crust
runs far brighter, alkaline trona and natron different again, and the brine work
now predicts which a basin gets. Albedo should follow the predicted mineralogy
rather than a single constant. `LITH-1`, `LITH-2`.

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

An attempt to ground the albedo values failed and is recorded honestly: three
literature queries phrased as topic descriptions rather than titles matched
nothing above 43%. Correct DOIs are needed. `LITH-4`.

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
