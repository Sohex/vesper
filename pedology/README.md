# Pedology

> Figures in this document are illustrative of method, and were measured on
> builds and climates that have since moved. Current values live in
> `world_state.json`, generated from the artifacts. See the convention in
> `CLAUDE.md`.



Turns rock into soil. Parent material sets what minerals are available, climate
sets how far they have been converted, relief and erosion set how much regolith
survives, and the biosphere sets the organic fraction. That last one is why this
is a loop and not a stage.

```bash
python pedology/scripts/build_soil.py                          # iteration 0
python pedology/scripts/build_soil.py --soil-carbon <cpool.out> --iteration 1
python pedology/scripts/weathering_fluxes.py                   # CO2 and silica
python pedology/scripts/brine_paths.py                         # the chemical divide
python pedology/scripts/build_surface_classes.py               # needs brine_paths
python pedology/scripts/land_column_properties.py              # the hydraulic property contract
```

Writes `data/<source_build>/soilmap.txt`, which is LPJ-GUESS's own `SoilInput`
format, and `analysis/soil_report.json`. The soil map is per build, because
texture derives from lithology.

## Derived surface classes

`build_surface_classes.py` answers a different question from `build_soil.py`.
Soil is what the profile has become; this is what is ON it and what is CEMENTING
in it, which is what albedo, dust emission and the phosphorus cycle key on.

It emits **two independent fields**, not one:

| field | question | values |
| --- | --- | --- |
| `surface_cover` | what wind and light see | `water`, `loess`, `diatomite`, `pavement`, `bare`, `soil` |
| `duricrust` | what is cementing at or below the surface | `gypcrete`, `calcrete`, `silcrete`, `none` |

Two axes because they do not compete for the same physical position. A duricrust
is a pedogenic horizon inside a profile; loess and diatomite sit on top of one. A
cell can be loess over calcrete. The `carved-zoned-v2` build lost basin fill to a
single chain that resolved exactly this kind of non-question by branch order, and
203 preserved basins came out with no fill cell anywhere.

Precedence exists only within an axis, is declared in
`config/surface_classes.yaml`, and every region records which rule FIRED and
which others it also SATISFIED, in a bitmask. A rule that never fires and a rule
that always fires are both bugs; neither is visible without that record, and the
report states it as pass or fail rather than as a number to read.

Writes `analysis/surface_classes.nc` on the native mesh, stamped with its
source build, plus
`analysis/surface_classes_report.json`. It needs `brine_paths.py` to have run,
because the duricrust ion gates and the silcrete and diatomite silica gates read
the per-basin chemistry, and there is deliberately no fallback: guessing an ion
source would decide which evaporite a basin grows.

Three things about it are worth knowing before reading its numbers.

**It is driven by a climatology, a lake solution and a dust field, so it is the
most disposable product in this component.** The carve replaces all three.

**The loess answer is a bracket, not a value.** The aeolian roughness bracket is
worth a factor of 40 in dust emission, and loess area runs from nothing to nearly
half the land across it. Quoting the central figure alone would be quoting the
narrowest part of the widest uncertainty in the chain. Pavement moves the other
way across the same bracket, because the two are one deposition threshold read
from two sides.

**Potential evaporation is taken per output bin, not per year.** Watson's
gypcrete criterion is potential evaporation exceeding precipitation in EVERY
month, and an annual mean passes cells whose wet season disqualifies them. The
implementation is hydrography's own Penman, imported rather than rewritten,
because that one is checked against the model's evaporation where the surface
genuinely is open water and a second implementation here would have no right
answer to be checked against.

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

**Regolith depth**, from a saturating balance of weathering production against
erosion built from Orogen erodibility, local relief and runoff. Heimsath's
exponential soil production function is what that form replaced, and
`build_soil.py`'s `regolith_depth` carries the argument. The law is applied to
EACH ROCK a cell holds and the depths are then averaged, through
`lib/gridding.py`'s `cell_expectation`: the law is convex in the erodibility and
erodibility spans a factor of fourteen across the export's lithology table, so
mixing a cell's rocks into one erodibility first is a different and always
thinner answer. `notes/audits/nonlinear-spatial-reductions.md` section 1 sizes
the two orders against `regolith.minimum_depth_m`, and the soil report's
`regolith_aggregation` block carries the difference on whatever map it just
built.

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

**The moisture driver is runoff, and the choice of it is a sensitivity rather
than an uncertainty.** Weathering is limited by water passing through the profile
and carrying solutes away, not by water that falls and evaporates, and the fits
underneath the Walker-Hays-Kasting law are against runoff, so precipitation was
never a co-equal alternative. It is reported as a spread where a result depends
on it, and the spread is 3.41x on weathering intensity and 2.83x on water
capacity. It was once quoted as 14.6x, which was a mismatch of reference values
rather than a bracket; "A unit error was inflating the headline uncertainty"
below has the arithmetic.

**This world's land runoff ratio is a sixth against Earth's roughly a third**:
drier than Earth, and not extraordinarily so. It was read as 2.8% for a while and
blamed on the climate model, on the grounds that ExoPlaSim's `configure()` clears
every surface field when handed a landmap, so soil field capacity falls back to a
uniform namelist default. Both halves of that were wrong. `mrro` is river-routed
rather than locally generated, so `P - E` is the quantity that drains through the
profile, and shrinking the bucket 3.75-fold moves runoff by only 1.06x, measured
below rather than argued. The low ratio is a property of the climate, high
evaporative demand against the precipitation available, and not of the bucket it
was blamed on. See `exoplasim/notes/water-and-energy-closure.md`.

**This is also where the loop runs in the other direction.** Soil water holding
capacity is exactly the field ExoPlaSim was defaulting, and this component
computes it; `model.soil_water_source: pedology` is now set in
`config/planet.yaml`. See below for what it is worth, which is little.

**The land means hide a bimodal world, and that is the fact to carry.** Measured
on the bootstrap climatology, 2026-08-17: land runoff averages 130 mm per Earth
year and its MEDIAN is 3.8. Two thirds of the land is under 50 mm/yr. So a mean
runoff ratio of 16% describes almost none of the surface; a wet minority carries
the whole figure, and the same skew runs through everything derived from it --
plant-available water capacity averages 151 mm against a median of 84, and
regolith depth averages 1.16 m against a median of 0.66.

It is also why land-mean weathering intensity is 0.380 while the intensity of the
land-mean climate is 0.572. Weathering goes as runoff to the 0.65, which is
concave, so on a skewed distribution the mean of the intensity sits well below
the intensity of the mean. Quoting either without the other overstates by 50% or
understates by a third. Neither is wrong; they answer different questions, and
the one that feeds the solute fluxes is the mean of the intensity.

**The texture model is scored against Earth by `validate_against_earth.py`**,
over type localities whose lithology is known from published geology. It puts
`build_soil.weather_texture` itself against SoilGrids, so the only term of the
conversion the score leaves out is the catena fines loss, which needs a slope a
type locality does not carry. `analysis/earth_validation.json` holds the bias,
the correlation and the mafic-felsic divergence the model exists to reproduce,
each against its observed counterpart.

The clay bias that figure once reported was largely the scoring arm's own: it
inlined the conversion without `sand_to_silt_loss_ratio` and without either cap
on the donor pools, so it predicted more clay than the model does. Any land-mean
clay figure quoted against the pre-fix bias is comparing to a number that
measured the score, not the map.

## The land column property contract

`land_column_properties.py` reads `config/land_column_properties.yaml` and is
the single vertically explicit description of the simulated land column's
hydraulic and thermal properties: two geometries, four materials with their
vertical rules, four retention states with their ordering checked, one NAMED
retention closure with its family stated, correlated uncertainty cases, and the
aquifer boundary the contract deliberately does not own. It runs no model.

**It is also the only place those states are derived.** It writes
`analysis/land_column_properties_report.json` and
`data/<build>/land_column_states_<res>.txt`, one row per land cell carrying
saturation, field capacity, the wilting point, the closure's exponent, the
plant-available capacity of the whole declared column, and the
weathered-bedrock usable share per physical layer. ExoPlaSim's
`build_surface_soil_water.py` installs the capacity column as `dwmax`;
LPJ-GUESS takes the states and the shares through its driver file. Neither
derives a retention curve of its own any more: `soilinput.cpp`'s Cosby
inversion and `vesperinput.cpp`'s regolith rescaling are gone.

The report still computes what each of those produced, because the cost of a
removal is only legible against the number it removed. It reports pedology's
own endmember `awc` and the shipped Cosby derivation on the same cells and in
the same units, with their ratio, their absolute difference and their
correlation; it measures the region where the closure would return a field
capacity above saturation -- an analytic line in the texture simplex, not a
search -- and reports the map's margin from it rather than a count.

**The adopted states carry their own bracket, from the same table as their
mean.** Cosby et al. regressed the within-texture-class standard deviation of
each retention parameter on texture alongside the mean, and that spread is what
texture does not determine. The report evaluates the corners of the one-sigma
box in the exponent, the air-entry head and the porosity, coherently across
every layer and every state, and reports the resulting capacity envelope
against the gravity shift in the same units. It is an envelope rather than a
distribution because Cosby publishes no correlation between the three
residuals. Every consumer installs the central case.

Gravity is where the contract earns its keep, and the correction is APPLIED.
Saturation is pore geometry and does not move; the wilting point is a plant
pressure and does not move; field capacity is a drainage equilibrium over a
stated length and its matric pressure is `rho_w * g * L`, so it does. The
declared drainage length reproduces Cosby's own field-capacity suction exactly
under Earth gravity, and the check fails if it stops doing so.

The trap that makes this worth a contract: Cosby's parameters are recorded as
heads on Earth, and air entry is a capillary pressure while the wilting point is
a plant pressure, so BOTH are invariant and both of their heads scale together.
Convert one and not the other and the wilting point moves for a bookkeeping
reason. The script works in pressure so a mixed frame cannot be written, and
checks the two consistent frames against each other and against the mixed one.
`notes/land-column-property-contract.md` is the argument.

Three arms are falsifiable: the declaration check is run against ten contracts
broken in named ways and is required to catch every one; the frame check
requires two frames to agree and a third not to; and the mutation harness fails
if a mutation goes uncaught. The properties still carrying the `undeclared`
sentinel are listed by the report with the issue that owns each, and `--strict`
is the arm that refuses while any remain.

## One hydraulic description; the two columns still evolve separately

Both consumers now read the same states, so they agree on how much water the
soil holds and on how deep the column goes. They still evolve separate snow,
soil water, ice, evaporation and runoff histories from it: one shared hydraulic
STATE is a further step, and `dwmax` is a scalar bucket with no vertical liquid
structure to share. The source audit and the work required are in
`../biosphere/notes/soil-land-surface-hydraulic-consistency-audit.md` and
LSHY-1 through LSHY-7.

The checked-in `analysis/soil_report.json` predates that source correction and
still says LPJ-GUESS does not consume regolith depth. The report generator now
states the narrower truth above. The derived report was not hand-edited or
regenerated during the no-execution audit; the next authorized `soil` step
will refresh it.

**To the biosphere.** `soilmap.txt` carries a `depth` column, LPJ-GUESS's own
`SoilInput` ignores it, and the contract turns it into a usable share per soil
layer that `vesperinput` multiplies each layer's water capacity by. The
arithmetic is unchanged from when `vesperinput` derived that share itself, which
is what the sweep below measured; what moved is where the rule lives. Verified
on one cell at four thicknesses, everything else held fixed:

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
divides by zero and the NaN propagates through the nitrogen substrate and
zeroes the gridcell's vegetation **silently**: zero LAI and zero
evapotranspiration rather than a
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

The numerical guard survives at 1e-4, as a guard only, with the physics
declared, literature-anchored and varying per cell; the controlled sweep above
is the measurement.

**To the climate.** `exoplasim/scripts/build_surface_soil_water.py` writes
surface code 0229, `dwmax`, from the land column property contract's `awc_mm`
column. It is smaller than ExoPlaSim's uniform 0.5 m default, and since
`drunoff = max(0, dwatc - dwmax)/deltsec`, a smaller bucket overflows sooner and
produces more runoff, which is the direction needed to fix the 2.8% ratio.

### The offline probe validates, and says the feedback is weak

`scripts/probe_runoff_response.py` reimplements ExoPlaSim's land bucket and
checks itself against the land water budget before predicting anything. It
passes:

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
ratio is 17.8% against Earth's 35%.

The soil water distribution shows why. The bucket sits at a median 15% of
capacity, so overflow comes from the wettest cells and seasons, which saturate at
either depth, while dry ground never fills at either. Capacity matters only in
the narrow band between. The low ratio is therefore a property of the climate,
high evaporative demand against available precipitation, and not of the bucket it
was blamed on. Implied weathering shift 0.96x, moving land-mean W from 0.50 to
0.48, which changes nothing downstream.

**`model.soil_water_source: pedology` is set**, which it was not while runs were
in flight: adding the key changes the configuration `continue_exoplasim.py`
compares against the run manifest, and it refuses to resume across a changed
value, so it was held until the
baseline re-run. The code still defaults to `uniform` when the key is absent.

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
to pedology's own capacity column: it is `volumetric x column depth`, and the
depth varies far more across the planet than the texture-derived volumetric
capacity does. **Depth sets that capacity almost entirely.** The capacity the
models install is the contract's, where texture enters through the retention
closure rather than through an endmember mixture, and where the depth term is
cut at the declared column base.

### The runoff-versus-precipitation spread

Each branch normalises against its own Earth land mean, 300 mm for runoff and
750 mm for precipitation; the spread is **3.41x** on weathering intensity and
2.83x on water capacity.

**And it is a sensitivity rather than a genuine uncertainty.** Walker, Hays
and Kasting define the law on river runoff, and the 0.65 exponent traces to
Dunne (1978) and Peters (1984), both fitted against runoff, so precipitation
was never a co-equal alternative. It was a hedge against `mrro`, which looked untrustworthy
at 25 mm per Earth year, and that hedge is spent: `mrro` turned out to be
river-routed rather than local, P - E gives 168 mm, the global water budget
closes to one part in ten thousand, and an offline bucket built from ExoPlaSim's
own scheme reproduces it independently at a ratio of 0.89.

Report the spread where a result depends on it. Do not present the precipitation
branch as an equally likely world.

With that demoted, the largest remaining uncertainty in this component is
`erosion_weight`, which is calibrated rather than measured and sets the level of
the water capacity.

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
- The derived surface classes have three gaps of their own, in order of size.
  **Diatomite is placed against a static lake proxy**, the strandline band a
  single climatology's lake could vacate, where the design asks for occupancy
  measured over the 57-year stellar component; that cannot be faked from one
  climatology, because the whole class is about alternation. **The loess
  threshold has no source**: Muhs (2013) was read for it and carries no
  accumulation rates at all, so 50 g/m2/yr is read off one figure's own
  before-and-after at a single place and bracketed 10 to 200. **The pavement
  clast supply is exhaustible and this rule cannot see that**: McFadden et al.
  (1987), read 2026-08-18, gives the criterion, that pavement clasts are
  weathered bedrock from topographic highs, and states that supply stops once
  the highs are worn down and buried, which a rule keyed on present substrate
  treats as permanent. So the pavement area is an upper bound on an old surface.
  No source gives a clast-supply threshold and none is invented.

## The carbonate-silicate thermostat, and how much of it this world has

Weathering only stabilises CO2 if its alkalinity reaches the ocean. A closed
basin weathers its catchment and then precipitates the carbonate on its own
floor, so it contributes nothing to the feedback. On a world with a large
endorheic share the thermostat is weaker than the total weathering rate implies,
and `thermostat_efficiency.py` measures how much weaker:

    thermostat efficiency = weathering over exorheic land / weathering over land

```bash
python pedology/scripts/thermostat_efficiency.py --compare-builds \
    --carve-list hydrography/data/carved-zoned-v4/carve_list.json
```

Measured on `climatology_s096`, with every terrain held against that one
climate so the spread is basin geometry and not weather:

| terrain | endorheic area | decoupled weathering | efficiency |
| --- | ---: | ---: | ---: |
| `precarve-unzoned` / `precarve-zoned` | 25.6% | 19.8% | 0.802 |
| `carved-zoned-v4` (configured) | 15.8% | 12.5% | **0.875** |
| projected, after the pending carve | 11.7% | 9.2% | 0.908 |

**The endorheic share decouples less alkalinity than its area implies.**
Endorheic land here is drier -- mean weathering intensity 0.401 against 0.525
on exorheic land -- so 15.8% of the land is 12.5% of the weathering. The
correction runs in the reassuring direction and it is not large: this world's
thermostat is degraded by about an eighth, not disabled.

Two caveats attach to the last row. It is a projection of the carve list's
retain fractions onto the *current* terrain's weathering field, not a
measurement: the carved export does not exist yet. And carving removes closed
depressions, which removes their bright evaporite fill and warms the world, so
the post-carve climate is not the climate this was computed on.

`carved-zoned` reports 0.936 under the same test and should not be read as part
of a trend. Its verdict was computed with the longitude bug, so every basin was
judged by its antipode's climate; it carved 1,522 basins where the corrected
verdict carves 749. It is a different hypothesis about which basins overflow,
not an earlier point on the same curve. For the same reason there is no useful
`hydrography` comparison across builds carved under different rules, and the
three-build reading that preceded this one is withdrawn.

### Whether the prescribed CO2 is backed by an outgassing rate

`weathering_fluxes.py` computes the outgassing this world's silicate weathering
implies and compares it against what the planet could supply. That comparison had
a requirement, a supply and a margin, and no machinery: no threshold, no branch,
no exit code, and a verdict that was a constant string rather than a conclusion.
`outgassing_gate.py` is the conclusion, on the shape `bvoc_gate.py` established:

```bash
python pedology/scripts/outgassing_gate.py          # report what is undeclared
python pedology/scripts/outgassing_gate.py --strict  # refuse on any of it
```

`config/outgassing.yaml` is the declaration. The requirement half carries no
sentinels because it is computed. The supply half is declared and external:
WORLD-G60E supplied the mass-scaling exponent and the chain it claims, the
per-kilometre boundary flux bracket, the four-way source partition with its own
fluxes and brackets, the degassing efficiencies where a source separates them,
and the delivery split between ocean and atmosphere. Each is a floor or a
bracket rather than a best estimate, so a margin computed from it is a margin
that survives the range. The acceptance floor is fixed in the declaration before
any number is seen, and a floor at or below 1 is itself refused, because that is
not a margin.

The gate refuses on two different kinds of defect. A field left undeclared is
refused by name. So is a field that is declared and wrong: a bracket whose ends
are the wrong way round, a value its own bracket does not contain, fractions
that do not sum to 1, and fractions that disagree with the fluxes they were
assembled from. That second kind is arithmetic, so it has a right answer and can
be refused rather than reviewed.

Seventeen fixtures ship with it: sixteen declarations built wrong in one named
way each that must be refused by a named code, and one with every precondition
met that must be granted. A gate nothing can satisfy refuses for a reason nobody
wrote down.

### What CHOOSING a weathering scheme costs the same number

The efficiency above is a geography measurement and carries no structural
uncertainty about the weathering law itself. Three published schemes sit on disk
with their constants, they disagree about the FORM of the law rather than about a
constant inside one form, and that disagreement is a thermostat strength.
`weathering_scheme_bracket.py` measures it:

```bash
python pedology/scripts/weathering_scheme_bracket.py
```

`config/weathering_schemes.yaml` declares the three schemes, the mapping from
this world's seventeen land rock classes onto rokgem's six or seven, and the
runoff-sensitivity sweep the bracket is reported against. The mapping is part of
the bracket rather than a detail under it: under GKWM the runoff exponent is per
lithology, so reassigning a clastic class between shale and sand moves the
effective exponent with no change of scheme at all. `evaporite` maps to nothing
on purpose, since a halite crust is not a silicate weathering substrate, and the
report says how much land that leaves unmapped.

The comparison is analytic and needs no climatology, because all three schemes
separate into a power law in runoff and an exponential in temperature. It is
therefore the one thermostat number available before a baseline run exists.

## The phosphorus leg

`phosphorus_budget.py` is presently a relative geography diagnostic, not a
closed phosphorus budget. It multiplies lithological P content by a relative
release factor, groups that rank by drainage fate, and tests whether dry basin
floors geometrically concentrate delivered material. It reads the lake solution
only to exclude wet dust sources. Despite the pipeline dependency, it does
**not** read the dust-deposition field, and neither weathering nor aeolian
delivery has absolute kgP/area/time units. Nothing in the biosphere currently
consumes its output.

`biosphere/notes/abiotic-nutrient-delivery-audit.md` records the source trace.
ANUT-1 through ANUT-3 define the missing mass ledger, absolute weathering and
composition/solubility-resolved dust input; ANUT-10 owns the initial pedogenic
mineral-P state. BIO-5 through BIO-8 own the eventual flux mapping into
LPJ-GUESS-CNP. SDEC-2/SDEC-4 retain the downstream sorption and occlusion
topology.

## Known gaps
- **None of the four mineral-reactivity proxies can be derived here, and the
  reasons differ.** The vegetation model's mineral-aware arm for organic matter
  protection and phosphorus sorption asks this component for an Fe-Al oxide
  content, an allophane concentration, an aggregate capacity and a polyvalent
  cation saturation. `notes/mineral-reactivity-supply.md` is the verdict on each
  with its evidence: what the near miss actually needs, which quantities are
  measured but indexed on a substrate age this project does not carry, and the
  second gate on the consuming model's side that would still be shut with all
  four in hand.
- **The derived surface classes carry the three gaps argued above** -- the
  static lake proxy under diatomite, the unsourced loess threshold, and the
  exhaustible pavement supply the rule cannot see.
  `notes/derived-surface-classes.md` has the design: two independent axes rather
  than one cover chain, gypcrete and calcrete windows, and the reason desert
  pavement enters the dust budget with the opposite sign to the one first
  assumed.
- **The thermostat efficiency ignores the ocean's own weathering budget.**
  Seafloor weathering and carbonate burial are not represented at all, so 0.875
  bounds the continental term only.
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
