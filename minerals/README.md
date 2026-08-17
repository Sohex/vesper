# Minerals

Ore prospectivity, as a 0-1 field per deposit type. The design and the split
between what belongs here and what belongs downstream are in
`../notes/economic-minerals.md`; this file is what the component does and how to
run it.

```bash
python minerals/scripts/build_prospectivity.py
```

Writes `data/<source_build>/prospectivity.nc` on the native mesh, and a
`prospectivity_report.json` beside it.

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
alone cannot tell you. A porphyry forms 1-5 km below the surface and is destroyed
by deep erosion; an orogenic gold system forms 5-15 km down and is revealed by
it. The same erosion that removes one exposes the other. Orogen tracks it through
`cover_thickness` and `erosionDelta`, which is why porphyry prospectivity here
rejects most of the arc belt rather than accepting all of it -- 59.5% of it is
stripped to its granodiorite root, and the porphyry level went with the section
above.

Kimberlite is the other entry worth noting: it needs a thick cratonic keel, and
before the craton basin double-count was removed this world had almost none. See
`../notes/audits/orogen-lithology.md`.

## Reading the numbers

Prospectivity is **relative within this world**, normalised by the land maximum
per deposit type. A 1.0 is the most favourable cell here, not a grade or a
tonnage, and comparing the number for one deposit type against another says
nothing. What is comparable is the spatial pattern within a type.
