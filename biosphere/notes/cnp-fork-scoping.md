# Scoping: LPJ-GUESS-CNP v1.0 as the biosphere model

Measured against the fork at tag `LPJ-GUESS-CNP_v1.0`
(`mateusdp/LPJ-GUESS-NTD`), cloned and diffed rather than reasoned about.

## Ordering: this does not gate the baseline re-run

The dependency runs through the soil map, not the terrain. `soilinput.cpp`
matches soilmap columns **by header name and ignores any it does not recognise**
-- our map already carries `depth`, `awc` and `bedrockfrac` that way, consumed by
our own `vesperinput.cpp`. So phosphorus columns can be added at any time without
disturbing stock LPJ-GUESS.

The real constraint is only that the P parameterisation be settled before the
*final* soil map is built, and the soil map is built late in the loop, after a
climatology exists. Missing that costs a soil-map regeneration, which is a grid
computation and not a model run. Terrain, hydrography, boundary conditions and
climate are all upstream of this and unaffected.

So: integrate in parallel with the baseline re-run, not before it.

## The port is nearly free

Our Vesper patch touches six files and the fork modifies all six, which predicted
a conflict. It does not: applied against the fork, **every hunk succeeds**, with
line offsets only and zero rejects. The CNP changes and the Vesper changes sit in
different parts of the same files.

    framework/guessmath.h     clean
    framework/guess.h         3 hunks, offsets -4 to +73
    modules/driver.cpp        3 hunks, offset +103
    modules/soil.cpp          6 hunks, offset +31
    modules/spinupdata.h      clean
    modules/CMakeLists.txt    clean

What still needs checking is `vesperinput.cpp`, which subclasses `SoilInput`.
The fork adds `kplab`, `spmax` and `pwtr` to `SoilProperties`, so the subclass
will need to carry them even if it does nothing with them.

## The input side is where the work is, and the fork is incomplete there

Three routes supply the P weathering rate `pwtr`:

1. **Soil-code table**, `data[soilcode][14]`, a constant per soil class.
2. **Texture path**, `soiltype.pwtr = soilprop.pwtr` -- which is the path we take,
   since our soil map has thirteen columns. But no soilmap column populates it,
   and the code above it reads `// FIXED FOR AMAZON FACE AT THE MOMENT`. The
   global texture path for phosphorus is **not wired to any input**. This is a
   site-scale model used at site scale; we would be the ones generalising it.
3. **Gridded**, via `file_pwtr` and the parameters `pwtr_bi`, `pwtr_pcont`,
   `pwtr_shield`, `pwtr_ea`.

Route 3 is the one the publication describes as temperature- and
humidity-dependent weathering, and it is a close fit for what this project
already computes. `pwtr_pcont` is the parent rock's phosphorus content -- exactly
the shape of the `quartz` fraction already declared per rock class in
`pedogenesis.yaml` -- and `pwtr_ea` is an Arrhenius activation energy, the same
form as the thermal term in our WHAK intensity.

So the integration is not "adopt a model", it is "supply the field the model
wants and cannot currently get", using a weathering law we already have.

## Why this world in particular

Phosphorus is rock-derived. Nitrogen arrives from the atmosphere and fixation,
which is why our nitrogen uncertainties came out second-order, but P has no such
source: its supply is weathering, and weathering is something this pipeline
already models spatially rather than assumes.

Two consequences specific to Vesper:

- **A heavily weathered world is a P-poor one.** Old, deeply weathered soils are
  P-limited on Earth, not N-limited. Our texture model over-predicts clay by
  about a factor of two, but the *pattern* validated against Earth, so the
  wet-and-old parts of this planet are real.
- **Endorheic drainage retains phosphorus.** P leached from an exorheic catchment
  reaches the sea; in a closed basin it stays. This is the mirror image of the
  carbonate-silicate thermostat result, where interior drainage *withholds*
  alkalinity from the ocean and weakens the thermostat. The same geography that
  degrades the thermostat should concentrate phosphorus, which predicts P-poor
  uplands against P-rich basin floors -- a biome-scale pattern this pipeline can
  derive rather than assert.

## Estimated work

    apply the patch to the fork, build, confirm parity      small, measured clean
    extend vesperinput.cpp for the new SoilProperties        small
    add rock-class phosphorus to pedogenesis.yaml            small, same shape as quartz
    emit pwtr/kplab/spmax per gridcell from build_soil.py    moderate
    read them in the texture path                            moderate, upstream gap
    P deposition, a declared constant as nitrogen's is       small
    re-register the productivity predictions                 required, not optional

The last one is not negotiable: the pre-registered predictions in
`notes/productivity-prediction.md` were registered against a C-N model. A C-N-P
model is a different model and needs its own registration with its own date, not
a widened band on the old one.
