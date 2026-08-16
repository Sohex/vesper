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
computes it. That is now built, though not enabled; see below.

## Both loops now close

**To the biosphere.** `soilmap.txt` carries a `depth` column, LPJ-GUESS's own
`SoilInput` ignores it, and `vesperinput` scales each soil layer's water capacity
by how much of that layer is really regolith rather than rock. Verified on one
cell at four thicknesses, everything else held fixed:

| profile | AET mm | LAI | NPP |
| --- | --- | --- | --- |
| 1.50 m (unscaled) | 393.2 | 2.899 | 0.375 |
| 0.60 m | 307.1 | 2.374 | 0.301 |
| 0.30 m | 264.6 | 2.068 | 0.259 |

Monotone and in the right direction. Cells with more regolith than the model's
1.5 m profile come out bit-identical to the unscaled run, which is the check that
the scaling does nothing where it should do nothing.

Two things that bit, both worth knowing before touching this again. LPJ-GUESS
carries soil water as a *fraction* of each layer's capacity, computing
`wcont = Faw_layer / soiltype.awc[layer]`, so a layer scaled to exactly zero
divides by zero and the NaN propagates through the nitrogen substrate and kills
the gridcell **silently**: zero LAI and zero evapotranspiration rather than a
crash. And because water is fractional, a thinner soil also reads as *relatively
wetter* for the same absolute water, so the net effect at any one cell can go
either way through PFT competition: in the smoke set, one cell at 0.87 m gained
LAI because temperate broadleaf evergreen displaced its grass entirely. The
controlled sweep above is the evidence, not any single cell.

### Water below the bedrock contact

Rock below the regolith is not dry, and deep-rooted woody plants use what is
there. `bedrock_water` in `pedogenesis.yaml` is how much, per unit volume,
relative to the soil above, and it varies with weathering intensity because that
is what converts impermeable rock into saprock and then saprolite:

    fraction(W) = minimum + (maximum - minimum) * (1 - exp(-shape * W))

Anchored on measurements rather than chosen. Intact bedrock between joints has
porosity of about 1% or less, so `minimum` is 0.05. Graham, Rossi and Hubbert
(GSA Today 2010) measured a Sierra Nevada Jeffrey pine site where 75 cm of soil
at 20% PAWC held 15 cm of water over 275 cm of saprock at 12% PAWC holding 33 cm,
so `shape` is set to put W = 1, the Earth land mean, near that 0.12/0.20 = 0.60.
The Swaziland Middleveld saprolites hold two to four times the available water of
their soils, which sets `maximum` at 2.0. Vesper's land mean W of 0.19 gives
0.170.

For scale: Rempe and Dietrich (PNAS 2018) measured 100-530 mm of seasonal rock
moisture, up to 27% of annual rainfall, and Lapides et al. (Biogeosciences 2024)
added a bedrock vadose zone to LPJ-GUESS itself with 180-480 mm of storage,
raising annual transpiration by a median of 100-150 mm and turning a model that
could not tell two catchments apart into one that reproduced both.

**The model cannot represent the top half of that range.** LPJ-GUESS ties a
layer's saturation capacity to its texture-derived porosity, so scaling `wsats`
above 1 drives `Frac_air` negative in `Soil::update_soil_diffusivities` and the
run dies. The fraction is therefore capped at 1.0: sub-bedrock material can match
the soil above but never exceed it. That is why Lapides et al. added a separate
bedrock layer with its own porosity rather than rescaling the existing profile.

**Measured worth, across the whole representable range**, at 0.30 m regolith with
only this parameter varying:

| fraction | | AET mm | LAI | NPP |
| --- | --- | --- | --- | --- |
| 0.050 | fresh rock | 364.2 | 2.989 | 0.368 |
| 0.170 | Vesper land mean | 354.0 | 2.912 | 0.348 |
| 0.598 | saprock, Earth mean W | 348.2 | 2.851 | 0.367 |
| 1.000 | model ceiling | 343.0 | 2.743 | 0.336 |

A span of 6% in AET and 9% in NPP, not monotone, with PFT composition stable
across all four. At three cells and `npatch 5` that is comparable to patch
stochasticity, so treat under 10% as an upper bound rather than a measurement.

This parameter entered as a hardcoded 0.02 in `vesperinput.cpp`, introduced as a
guard against the division by zero above. A first attempt to measure it, varying
that guard 40-fold, suggested it was worth 20-41% of AET. **That figure was
wrong** and should not be quoted: it changed the interpolation formula at the
same time as the value. The controlled sweep above is the number. The numerical
guard survives at 1e-4, but as a guard only, with the physics now declared,
literature-anchored and varying per cell.

**To the climate.** `exoplasim/scripts/build_surface_soil_water.py` writes
surface code 0229, `dwmax`, from the `awc` column. Land-mean capacity is 0.133 m
against ExoPlaSim's uniform 0.5 m default, and since
`drunoff = max(0, dwatc - dwmax)/deltsec`, a smaller bucket overflows sooner and
produces more runoff, which is the direction needed to fix the 2.8% ratio.

### The offline probe now validates, and says the feedback is weak

`scripts/probe_runoff_response.py` reimplements ExoPlaSim's land bucket and
checks itself before predicting anything. It used to fail that check by 5.91x,
which is what prompted looking at `mrro` and finding it river-routed rather than
local. **The probe was right and the target was wrong.** Validated against the
land water budget instead, it passes:

| | mm per Earth year |
| --- | --- |
| land budget, P - E | 167.83 |
| offline bucket at ExoPlaSim's 0.5 m | 149.96 |
| ratio | 0.89, inside the 0.50-2.00 band fixed beforehand |

Mean soil water also matches, 0.204 m against the model's 0.202, so the
prediction is usable rather than indicative:

| bucket | runoff mm/Earth-yr | ratio to P |
| --- | --- | --- |
| 0.500 m, ExoPlaSim's default | 149.96 | 0.168 |
| **0.133 m, from pedology** | **158.46** | **0.178** |

**Shrinking the bucket 3.75-fold moves runoff by 1.06x.** That is a validated
result now rather than a hint, and it settles the question the feedback was built
to answer: the soil will not on its own explain why this world's land runoff
ratio is 18.8% against Earth's 35%.

The soil water distribution shows why. The bucket sits at a median 15% of
capacity, so overflow comes from the wettest cells and seasons, which saturate at
either depth, while dry ground never fills at either. Capacity matters only in
the narrow band between. The low ratio is therefore a property of the climate,
high evaporative demand against available precipitation, and not of the bucket it
was blamed on. Implied weathering shift 0.96x, moving land-mean W from 0.50 to
0.48, which changes nothing downstream.

**It is off by default and the flag is not in `config/planet.yaml`.** It needs
`model.soil_water_source: pedology`, and that key is deliberately absent, because
adding it moves `config_sha256` and `continue_exoplasim.py` refuses to resume a
run whose config hash has changed. The code defaults to `uniform` when the key is
missing, so nothing changes for a run in flight. Add the key when re-baselining.

## Catena: physical production and topographic transport

Chemical weathering alone does not make a soil. Two terms sit on top of it, both
cheap because the inputs already existed.

**Frost shattering** breaks rock far faster than dissolution, and it needs
freeze-thaw *cycling* rather than cold: a bin whose diurnal range straddles
freezing cracks rock, one frozen solid all bin does not. The climatology carries
`mint` and `maxt` as timestep extrema, so the straddling fraction is directly
measurable. It comes out at 0.088 of the orbit as a land mean, against 34.8% of
land having a frost season at all.

**Topographic transport** moves regolith downhill, so ridges keep thin stony
lithosol and valley floors accumulate deep clay-rich fill:

    depth = depth_weathered * (1 + frost_bonus * frost_fraction)
                            / (1 + slope_transport * tan(beta))

with fines shed preferentially, moving clay to sand on slopes.

**The slope has to come from the mesh, and that is the interesting part.** Catena
is a hillslope process; a T42 cell is about 330 km across. Differencing
neighbouring cell centres gives a land-mean gradient of **0.001**, three orders
of magnitude below a real hillslope, and the term does nothing at all. Taking the
within-cell spread of mesh elevation over the 15.19 km mesh spacing instead gives
**0.031**, thirty times larger and spatially structured.

It is still an underestimate. Real catenas run at 100 m scale and gradients of
0.1 to 0.5, so even the mesh is coarse by two orders of magnitude. This term
therefore reproduces the *pattern*, mountains thin and basins deep, and not the
absolute magnitude. `slope_transport` is set against that pattern, which is a
weaker claim than a calibration and is stated as such.

### The depth model railed, and was replaced

Adding the catena made an existing weakness visible: a quarter of land sat on the
0.02 m floor and a quarter on the 5.00 m ceiling. That was
`depth = h_star * ln(production / erosion)`, Heimsath's exponential production
function at steady state, which diverges as erosion approaches zero and runs to
minus infinity as it grows.

The cause is dynamic range. Spanning 0 to 5 m through a logarithm at Heimsath's
`h_star` of 0.5 m needs production over erosion to cover a factor of **22,000**.
Nothing in these inputs covers that, so nearly every cell landed outside and was
clipped.

Two changes fixed it.

**A saturating form**, bounded at both ends by construction:

    depth = maximum_depth * P / (P + erosion_weight * E)

As erosion vanishes the profile approaches `maximum_depth`, which is the physical
statement that a weathering front cannot advance forever because water and oxygen
have to reach it through what has already accumulated. As erosion grows it thins
smoothly to bare rock. It is a parameterisation and not a derivation, unlike the
form it replaces, and that is the trade: a curve that spans the range against a
principled one that cannot be evaluated over it.

**Erosion that does not require runoff.** With erosion built purely from runoff,
it was exactly zero wherever `P - E` was, so the entire arid fraction pinned to
the ceiling: 18% of land, even after the saturating form removed the floor
problem. A slope in a desert still loses material to wind, dry ravel and creep,
so a `dry_erosion_baseline` of 0.15 now floors the moisture term.

The result:

| percentile | depth m |
| --- | --- |
| 5 | 0.09 |
| 25 | 0.29 |
| 50 | 0.66 |
| 75 | 1.46 |
| 95 | 3.24 |

Land mean 1.04 m, **0.0% of land on the floor and 0.5% on the ceiling**, against
25% on each before. `erosion_weight` is the one free scale parameter and there is
no measurement of this world's soils to fit it to, so it is calibrated against a
declared Earth-analogue target of about 1 m mean thickness. That target was
chosen before looking at what it does downstream, which is the difference between
calibrating and tuning.

The land-mean plant-available water capacity that falls out, 133 mm, sits in the
middle of Earth's typical 100-200 mm root zone. That is a check on the result
rather than an input to it.

**This changes the downstream numbers.** Water capacity was 347 mm under the
railing model and is 133 mm now, so the `dwmax` field handed to ExoPlaSim is
0.133 m against its uniform 0.5 m default rather than 0.342 m. The soil-water
feedback on runoff is therefore a larger perturbation than previously estimated,
not a smaller one.

## Where this stands, honestly

Checked over all 4,106 land cells rather than asserted.

**Sound.** Textures sum to 1 to rounding, nothing is NaN or negative, and no
field rails except `orgc` and `bulkdensity`, which are zero and uniform because
iteration 0 has no biosphere. Every relationship the model claims to represent is
present in the output with the right sign: clay rises with temperature (+0.31)
and runoff (+0.60), sand mirrors it (-0.42), pH falls with leaching (-0.57), and
water capacity tracks depth and texture as it must.

**Plausible against Earth**, where a comparison exists at all:

| | Vesper | Earth |
| --- | --- | --- |
| clay | 0.30 | 0.20-0.30 |
| sand | 0.40 | 0.35-0.45 |
| pH | 6.85 | ~6.5, higher in arid |
| plant-available water | 130 mm | 100-200 mm root zone |
| regolith depth | 1.03 m | 0.5-2 m |

### The pattern is terrain; the level is weathering

Two things are true at once and they are easy to conflate.

*Within* a configuration, depth is set by terrain and almost nothing else.
Regressing log depth on climate and soil variables gives r-squared of 0.006
against runoff, 0.027 against temperature, 0.063 against the bedrock fraction.
Relief, erodibility and slope do the work. That matches Earth, where soil depth
is a topographic story rather than a climatic one.

*Between* configurations, the weathering driver moves the whole field. Switching
the moisture variable from runoff to precipitation multiplies land-mean depth by
2.74 and water capacity by 2.83.

So the texture work, which is where most of the physics went, contributes little
to the field that actually reaches the other models: `awc = volumetric x depth`,
and depth varies 229-fold across the planet while the texture-derived volumetric
capacity varies 4.5-fold. **Depth sets the water capacity almost entirely.**
Texture matters to LPJ-GUESS's own soil physics, and to anyone reading a soil
map, but not to the number that feeds the bucket.

### A unit error was inflating the headline uncertainty

The runoff-versus-precipitation bracket was reported as 5.6x. It was comparing
this world's 892 mm of *precipitation* against Earth's 300 mm of *runoff*: not a
bracket, a mismatch of references. Each branch now normalises against its own
Earth land mean, 300 mm for runoff and 750 mm for precipitation, and the bracket
is **3.41x** on weathering intensity and 2.83x on water capacity.

That is still the largest uncertainty in this component and it is a modelling
choice rather than a defect: does chemical weathering follow water that drains
through the profile, or water that arrives on it. Runoff is the physically
correct answer and is the default. The bracket is reported beside every result.

### What is not yet earned

- `erosion_weight` is calibrated against a declared 1 m target, not measured.
  It sets the level of the field that matters most.
- The catena slope term carries the pattern and not the magnitude, because even
  the 15.19 km mesh is two orders of magnitude coarser than a hillslope.
- Nothing here has been validated against an independent product, unlike the
  Penman evaporation which was checked against the model's own ocean cells. The
  Earth comparison above is a plausibility check, not a validation.
- The soil-biosphere loop has run one iteration at smoke scale and has never
  been iterated to its convergence criteria.

## Known gaps
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
