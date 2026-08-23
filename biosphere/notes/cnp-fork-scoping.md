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

`VesperInput` contains a `SoilInput`; it does not subclass it. The build therefore
needs no new members in `vesperinput.cpp`. The actual defect is in the fork's
`SoilInput::get_mineral`: it copies `kplab`, `spmax` and `pwtr` from fields that
the texture reader never initializes. The vendored C-N configuration supplies
the fork's documented site defaults while `ifplim` is off. Replacing those
defaults with per-cell Vesper fields is the BIO-5 through BIO-7 chain of issues.

## The input side is where the work is, and the fork is incomplete there

Three routes supply the P weathering rate `pwtr`:

> **Unit correction, 2026-08-21:** the fork calls the texture-path value
> kgP/m2/year and divides it across one model orbit. On Vesper, “year” is
> ambiguous by a factor of about 2.02. BIO-24 must define the absolute-time
> contract before BIO-5 derives this field; the runoff-driven daily route and
> BIO-8's daily deposition input do not share that ambiguity.

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
- **Endorheic drainage retains phosphorus -- but that is the hydrological leg
  only.** P leached from an exorheic catchment reaches the sea; in a closed
  basin it stays. Aeolian transport from dry lake beds runs the other way and
  partly refills the uplands (failure-modes class 8, `notes/dust.md`), so
  expect a real but shallower upland-to-basin gradient -- a biome-scale
  pattern this pipeline can derive rather than assert.

## Integration seams tracked as issues

The source and the rock-side information are already present: the CNP fork and
Vesper port are vendored, while `pedogenesis.yaml` carries both per-class
phosphorus content and release factors and `phosphorus_budget.py` reads them.
What remains is split by contract rather than hidden inside one BIO-5 row:

    BIO-5   derive pwtr/kplab/spmax from the existing pedology fields
    BIO-6   emit the three per-cell soil-map columns
    BIO-7   consume them in the fork and remove the temporary defaults
    BIO-8   declare and propagate phosphorus deposition
    BIO-9   register the C-N-P prediction before seeing a result
    BIO-10  enable P limitation and validate parity and mass balance

BIO-9 is not negotiable: the pre-registered predictions in
`notes/productivity-prediction.md` were registered against a C-N model. A C-N-P
model is a different model and needs its own registration with its own date, not
a widened band on the old one.
