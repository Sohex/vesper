# Pedology

Turns rock into soil. Parent material sets what minerals are available, climate
sets how far they have been converted, relief and erosion set how much regolith
survives, and the biosphere sets the organic fraction. That last one is why this
is a loop and not a stage.

```bash
python pedology/scripts/build_soil.py                          # iteration 0
python pedology/scripts/build_soil.py --soil-carbon <cpool.out> --iteration 1
```

Writes `data/soilmap.txt`, which is LPJ-GUESS's own `SoilInput` format, and
`analysis/soil_report.json`.

## Why this component exists

The first version of the biosphere driver mapped World Orogen rock classes
straight onto LPJ soil texture codes. That table was parent material, not soil.
Soil texture is what a rock weathers *into under a climate*: the same granite
gives coarse grus in a cold arid place and deep kaolinitic clay in a wet tropical
one, and a lookup on rock type alone cannot tell them apart.

Every input a real pedogenesis stage needs already existed in this pipeline,
which is what made the shortcut conspicuous. Jenny's five soil-forming factors
are climate, organisms, relief, parent material and time. ExoPlaSim has the
first, LPJ-GUESS produces the second, Orogen has relief, erodibility and
lithology for the third and fourth, and the fifth is a modelling choice.

## Portability is the point

Nothing in `scripts/` hardcodes a Vesper number. Every Earth calibration lives in
`config/pedogenesis.yaml` with its source and what it is worth; the scripts read
the shared `config/planet.yaml`, a World Orogen export and an ExoPlaSim
climatology, none of which are world-specific in structure. Moving this to
another planet is a config change, not a code change.

The two places a different world would most likely need to argue with the
defaults are the weathering law's temperature e-folding, which is an Earth
silicate-weathering fit, and the fresh-regolith texture table, which assumes
Earth-like mineralogy. Both are one file away.

## What it computes

**Weathering intensity**, the Walker-Hays-Kasting form that the exoplanet
weathering literature runs on: a power of runoff times an exponential of
temperature, normalised so Earth's land mean is 1. It is an intensity rather than
a rate, because soil age is not known well enough on either world to carry
explicitly.

**Texture**, by mixing each cell's parent materials and then converting
weatherable primary minerals to clay as a saturating function of intensity.
Quartz is tracked separately and never converts. That single distinction is why
granite and basalt diverge under identical climate: both start sandy, but only
one has sand that can weather.

**Regolith depth**, from Heimsath's exponential soil production function balanced
against erosion built from Orogen erodibility, local relief and runoff.

**pH**, from parent material leached down by drainage, and pushed back up where
drainage is closed and salts concentrate instead of leaving.

**Organic fraction and bulk density**, from LPJ-GUESS's own soil carbon through
the standard reciprocal mixing rule, which is why organic soils come out light.

## The loop

```
iteration 0   build_soil.py                  no biosphere, organic carbon = 0
              LPJ-GUESS                      vegetation on mineral soil
iteration 1   build_soil.py --soil-carbon    soil carries the biosphere's carbon
              LPJ-GUESS                      vegetation on the soil it made
...           until the convergence criteria in pedogenesis.yaml are met
```

Criteria are fixed in `config/pedogenesis.yaml` before any iteration has run, in
the same way the albedo bracket's were written into a docstring before any run
finished: land-mean soil carbon stable to 2%, and fewer than 5% of land cells
moving more than 0.02 in clay content, within 6 iterations.

Verified end to end at smoke scale. Three cells run with mineral soil produce
8.1 to 9.1 kgC/m2 of soil carbon; fed back, they carry organic fractions of
0.021 to 0.024 and bulk density falls from 1500 to about 1330 kg/m3. The
biosphere is modifying the soil it grows in.

## The finding that matters most so far

**Weathering intensity is bracketed by a factor of 14.6 and the bracket is the
climate model's land hydrology, not the soil model.**

| moisture driver | land-mean W |
| --- | --- |
| runoff | 0.192 |
| precipitation | 2.808 |

Runoff is the physically correct term: weathering is limited by water passing
through the profile and carrying solutes away, not by water that falls and
evaporates. But ExoPlaSim reports 25 mm per Earth year of runoff against 892 of
precipitation, a runoff ratio of 2.8% where Earth's land manages about 35%.

That is very likely the model rather than the planet. ExoPlaSim's `configure()`
clears every surface field when handed a landmap, so soil field capacity falls
back to a uniform namelist default; `model.uniform_land_surface` in
`config/planet.yaml` already declares this. A bucket with the wrong depth
everywhere will not produce a believable runoff field.

`runoff` is the default because it is the right variable, and the bracket is
reported beside every result rather than resolved by picking the convenient one.
Under it, Vesper's soils come out markedly less weathered than Earth's: land-mean
clay 0.276 against sand 0.407, with 0.286 of the mineral fraction inert quartz.

**This is also where the loop wants to run in the other direction.** Soil water
holding capacity is exactly the field ExoPlaSim is defaulting, and this component
computes it. Feeding it back would close a second loop and probably fix the
runoff at the same time. Not built.

## Known gaps

- **Regolith depth is computed but not consumed.** LPJ-GUESS 4.1.1 has a fixed
  1.5 m profile. A thin soil over bedrock holds far less plant-available water
  than a deep one, and at present the model cannot be told. Emitted and reported
  so the gap is visible; the fix is a patch scaling available water capacity by
  depth.
- **Soil water capacity is not fed back to ExoPlaSim**, as above.
- **Time is not represented.** Weathering intensity folds the time integral into
  its normalisation, so a young volcanic surface and an ancient craton weather
  identically under the same climate. Orogen has exhumation data that could
  support a real age term.
- **Salinity is only a pH bonus.** Endorheic basins raise pH but sodicity,
  osmotic stress and the actual salt budget are not modelled, and LPJ-GUESS has
  no salinity response to receive them anyway.
- **C:N is a single declared constant**, because this world has no measured
  nitrogen cycle and the deposition rate feeding LPJ-GUESS is itself an
  assumption.
