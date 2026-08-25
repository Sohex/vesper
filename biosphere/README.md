# Biosphere

> Figures in this document are illustrative of method, and were measured on
> builds, climates and model source that have since moved. Current values live
> in `world_state.json`, generated from the artifacts. See the convention in
> `CLAUDE.md`.



Replaces the assumed vegetation with a modelled one. LPJ-GUESS takes the
ExoPlaSim climatology and returns leaf area, carbon and plant functional type
composition per gridcell, which becomes the surface albedo and forest fraction
that the next climate run is forced with.

The model is mechanically integrated and has run end to end on real Vesper
cells at smoke scale. A follow-up source audit exposed ecological-calendar and
time-base correctness work below that interface; the smoke result is therefore
an integration result, not yet a meaningful biosphere result.
An interpretable full run still needs the surface-area, forcing, aggregation and
acceptance contracts tracked as BIO-11 through BIO-15 and the newly exposed
calendar, time-base and forcing corrections in BIO-21 through BIO-25.  The audit
behind the original porting choices is in `notes/lpj-guess-porting-audit.md`, the
remaining modelling gaps are evidenced in `notes/modelling-gap-audit.md`, the
implicit Earth assumptions below the port are in
`notes/implicit-earth-assumptions.md` and the absolute-day/Earth-year/orbit unit
contract they are resolved against, with the class and the reader of every
quantity that carries a year, is in `notes/time-base-unit-contract.md`, the soil decomposition, C-N-P and
pedology/groundwater seams are audited in
`notes/soil-decomposition-biogeochemistry-audit.md`, every phosphorus constant
the CNP fork runs on is registered with its source or its bracket in
`notes/phosphorus-cycle-parameterisation.md`, plant physiology and carbon
allocation are audited in `notes/plant-physiology-carbon-allocation-audit.md`,
BVOC emissions, secondary organic aerosol and atmospheric coupling are audited
in `notes/bvoc-soa-atmospheric-coupling-audit.md`, wetlands, peat and methane
are audited in `notes/wetlands-peat-methane-audit.md`, and
the ExoPlaSim-to-LPJ weather path is audited in
`notes/ecological-climate-forcing-audit.md` and the units, time base, sign,
area basis and converting side of every field crossing that seam are settled in
`notes/ecological-forcing-field-contract.md`; EFOR-1 through EFOR-8 replace the
project's artificial 12-bin forcing scaffold with a chronological,
interval-explicit contract and carry authoritative land state. Soil and
land-surface hydraulic consistency across ExoPlaSim, pedology, hydrography,
groundwater and LPJ-GUESS is audited in
`notes/soil-land-surface-hydraulic-consistency-audit.md`; LSHY-1 through LSHY-7
define one property contract, fast-state owner and water/energy ledger rather
than two independent land columns. Abiotic nutrient sources and their delivery
from rock, dust and lightning through that water ledger are audited in
`notes/abiotic-nutrient-delivery-audit.md`; ANUT-1 through ANUT-10 distinguish
total material from root-zone-available N and P and close their destinations,
and the ledger they book into, with its control volumes, its balance identity
and the fixtures that fail when it does not hold, is
`notes/abiotic-nutrient-ledger.md`.
Finally,
`notes/productivity-prediction.md` registers what the answer should be before the
model can contradict it.

## Why this component exists

The albedo bracket left the biosphere as the largest unresolved lever on this
world. Bare rock gives a land-mean albedo of 0.276 and a vegetated surface 0.179,
worth several kelvin, and the two endmembers reach the design band then in force
at non-overlapping stellar fluxes. The bracket's windows (0.977-0.994 bare against
0.952-0.970 vegetated) were measured on pre-carve terrain at T21; the flux has
since been measured directly on the current terrain and the vegetated baseline
is 0.945. The project picked vegetated, which means the current baseline is
habitable *because* the world is assumed to be vegetated. That assumption is
currently a constant in `config/planet.yaml` (`land_albedo_source: vegetated`,
0.15 on anything that can carry a canopy). This component is what turns it into
a result.

## State

| item | state |
| --- | --- |
| model obtained | `mateusdp/LPJ-GUESS-NTD`, tag `LPJ-GUESS-CNP_v1.0`, commit `b368b893`; MPL-2.0 notices restored from verified 4.1.1 |
| build | vendored CNP tree plus Vesper port; build verification awaits an available compute window |
| smoke test | stock 4.1.1 port: bundled 3-cell demo, 550 years, 73 s, expected PFTs |
| Earth-assumption audit | follow-up complete; BIO-21 through BIO-29 track cross-cutting assumptions including the latitude-use contract, PCAR-1 through PCAR-10 track plant physiology/allocation, SDEC-1 through SDEC-9 track soil decomposition/biogeochemistry, BVOC-1 through BVOC-10 track volatile carbon through chemistry/aerosol/climate, WET-1 through WET-11 track wetlands, peat and methane, LSHY-1 through LSHY-7 track the shared land-water and soil-property contract, and ANUT-1 through ANUT-10 track abiotic nutrient delivery |
| productivity prediction | registered, unscored |
| calendar and astronomy port | mechanical calendar and orbital geometry applied; natural phenology still has unreachable Earth dates under BIO-21 |
| PFT time base | one contract with four classes, in `notes/time-base-unit-contract.md`; `build_vesper_pfts.py` executes it and names each parameter's class in the file it writes and in its provenance |
| input module | `vesperinput`, runs end to end and splits across MPI ranks; its 12-bin `VESPDRV5` transport is integration scaffolding to be replaced under EFOR-1 through EFOR-8 |
| soil and water | pedology depth scales LPJ capacity and pedology AWC sets ExoPlaSim's scalar bucket at smoke scale, but the models independently derive hydraulic properties and run separate snow/soil water balances; LSHY-1 through LSHY-7 own the consistency work |
| abiotic nutrients | the ledger is DEFINED and does not CLOSE: `abiotic_nutrient_ledger.py` carries thirteen control volumes and twenty-one terms, every one of them still holding the `undeclared` sentinel with its owning issue named. Rock P and dust mass have useful relative/source artifacts but no absolute flux; `phosphorus_budget.py` does not consume dust deposition, and ANUT-2 through ANUT-6 and ANUT-10 own the terms that would close it |
| abiotic source screen | geomorphic renewal, arc tephra and marine aerosol are RETAINED against the ledger, volcanic sulfate deposition is registered and not implemented, and fire ash and lightning belong to FIRE-7 and ANUT-4. The exhumation and tephra rates come from Earth's stationary population and not from the terrain; no screen may be carried on an aerosol optical depth |
| non-N/P adequacy | screened as a critical runoff per element and per lithology, bounds and directions declared before the result. Potassium binds; iron and the trace set REFUSE for want of a release table. The declared model boundary is that any LPJ-GUESS result here is a C-N-P result and not a nutrient-limitation result |
| tissue stoichiometry | the fine-root and sapwood C:N and C:P windows are anchored on the tissue MEAN, so the proportion `canexch.cpp` applies to their nutrient demand is the one Friend et al. (1997) measured. The max-anchored form the model used to carry applied 1.79 times it on nitrogen and 2.22 on phosphorus, so every C-N number produced before the repair is worthless rather than stale and the biosphere needs re-commissioning. `notes/plant-physiology-carbon-allocation-audit.md` finding 11 |
| phosphorus parameters | every constant registered with its source or its bracket in `notes/phosphorus-cycle-parameterisation.md`; the uptake profile, the leaf C:P window, the root proportion and the labile-P saturation threshold are derived, the sapwood proportion and the litter-P saturation threshold are not, the labile-P threshold reads a pool its source did not define, and `ifplim 1` fails closed naming each. The sapwood refusal is now about the element, the proportional form and the level alone: the model applies the constant it declares |
| phosphorus sinks | leaching, fire and harvest only. Terminal occlusion is a DECLARED ABSENCE, argued in the same note, so a simulated soil that must be old carries its phosphorus depletion in its initial stocks rather than developing it |
| run harness | written; records inputs, binary and model identity in its manifest |
| albedo and forest feedback | modelled mode exists; rootable/lake and spectral corrections are BIO-17 and BIO-18 |
| aerodynamic feedback | modelled roughness is open as BIO-16 |
| BVOC/SOA feedback | LPJ source is present but off; carbon closure, PFT traits, reduced atmospheric chemistry/transport, direct optics and cloud effects are separated under BVOC-1 through BVOC-10 |
| wetlands/peat/methane | LPJ's northern-Earth peat/CH4 source is present but off; current low-latitude saturation creates water and the full extent, groundwater, peat-stock, source/sink and atmospheric closure is separated under WET-1 through WET-11 |
| full run | harness exists, but BIO-11 through BIO-15 and BIO-21 through BIO-25 must close before its output is interpreted; prior estimate ~25 min on 16 ranks |

The model is vendored at `vendor/lpj-guess/` as a git subtree from the CNP fork.
The Vesper calendar, astronomy and input-module changes live directly in that
tree, so the source being reviewed is the source being compiled and run.
Phosphorus limitation remains off until Vesper's soil-P fields and a replacement
productivity prediction land together; this is a source normalization, not a
silent change from the established C-N experiment.

```bash
cmake -S vendor/lpj-guess -B vendor/lpj-guess/build \
  -DCMAKE_BUILD_TYPE=Release -DUNIT_TESTS=OFF
cmake --build vendor/lpj-guess/build --parallel 16
```

`framework/vesper.h` is still generated and ignored; it contains the planetary
numbers. The committed port only includes it and wires the in-tree
`modules/vesperinput.*` into the build. See "Running it" below for the
generators.

## The three things that decide whether this is credible

**The calendar.** Vesper's year is a few hundred Earth days and its day is 30
hours, so nothing about LPJ-GUESS's 365 x 24 h grid survives contact. Settled and
patched: 24-hour steps, a year length derived from `lib/orbit.py` and rounded to
whole days, which keeps every per-day rate constant
calibrated against the absolute time it was calibrated against and confines the
error to daylength alone. Verified against the unpatched model on identical
forcing, where annual evapotranspiration falls to 0.492 of its former value
against an expected 0.496.

**The PAR fraction.** LPJ-GUESS assumes half of shortwave is photosynthetically
active, which multiplies straight into productivity; the wrong spectrum is a
factor of four on it (0.078 from the mislabelled `k2.dat`, the M2.5V star
K2-18, against 0.309). `FRADPAR` is derived from the BT-Settl K2.5V spectrum
at `exoplasim/inputs/stellarspectra/k25v` by `build_vesper_header.py`, which
writes the value and its derivation into `generated/vesper.h` and
`generated/vesper_provenance.json`. See
`exoplasim/notes/stellar-spectrum-audit.md`.

**The PFTs.** The shipped plant functional types are Earth's, and their
bioclimatic limits are Earth calibrations. Keeping them is defensible as an
Earth-analogue biosphere and should be declared that way rather than presented as
a prediction. Their degree-day thresholds must be rescaled to the length of a
Vesper year, and this is not optional: running the patched model on Earth's own
demo data collapses boreal needleleaf and temperate broadleaf to grass, because
Earth `gdd5min` cannot be met in a Vesper year. `build_vesper_pfts.py` does it,
deriving the factor from the configured orbit. It deliberately leaves
`phengdd5ramp` alone, which is a within-season accumulation already in absolute
time, and scaling that would be a real error.

There is a second PFT question the literature answers more sharply than expected.
Earth's 400-700 nm photosynthetic window is an accident of our star, and Lehmer
et al. 2021 predict peak pigment absorbance around a K2V at 675, 711 and 746 nm.
A 400-750 nm window on this spectrum gives 0.99 of Earth's photon flux, so a
50 nm redward shift would almost exactly cancel the dimmer, redder star. The
21% oxygen atmosphere caps how far that can go, since anything past about 800 nm
needs anoxygenic photosynthesis. See the prediction note.

## What the answer is expected to be

Registered in full, with tolerances and known biases, in
`notes/productivity-prediction.md`. The short version: a dimmer, redder star
under a clearer sky delivers 0.84 to 0.89 of Earth's PAR photons per square
metre of land, and a warmer, wetter climate more than repays it, so **per unit
area this world should be productive within about 30% of Earth**. It has 2.11
times Earth's land area, so **total NPP should land near 2.1 times Earth's**.

The line that actually decides the pipeline is not either of those. It is tree
cover, predicted at 55% of land: the climatology that produced these numbers
assumed a vegetated surface, and if LPJ-GUESS returns substantially less canopy
than that, the climate is not one that biosphere would sustain and the loop turns
again.

## Soil comes from `pedology/`

`vesperinput` takes `file_soilmap`, the map `pedology/scripts/build_soil.py`
emits, and uses LPJ-GUESS's own richer `SoilInput` path: sand, clay, silt,
organic carbon, pH, bulk density and C:N per cell rather than one of ten texture
codes. The driver file still carries a soil code per cell, derived from
lithology alone through `SOIL_CODE_BY_ROCK`, and that is the fallback when no
soil map is given. Which one was used is printed at the top of the run, because
the difference is invisible in the output otherwise.

The fallback is parent material, not soil, and the distinction matters: the same
granite gives coarse grus in a cold arid place and deep kaolinitic clay in a wet
tropical one. Prefer the soil map.

That component closes a loop with this one: soil organic matter is a product of
the biosphere, and it changes the bulk density and water-holding capacity the
biosphere then grows in. Iteration 0 runs on mineral soil; every iteration after
feeds `cpool.out` back into `build_soil.py`. See `pedology/README.md` for the
convergence criteria, which are fixed in advance, and for the finding that
the runoff-versus-precipitation choice is a 3.41x sensitivity on weathering
intensity rather than a bracket.

## Current multiyear adapter, and why it is not accepted weather forcing

The current driver file can carry several input files and `vesperinput` cycles
through them:

```bash
python biosphere/scripts/build_lpj_driver.py \
    --climatology year0.nc year1.nc ... yearN.nc
```

Each file contributes one nominal forcing year in order. This is mechanically
useful, but every file is still a 12-bin regular product and LPJ-GUESS smooths
those bins into quasi-daily weather. `build_climatology.py --per-year` preserves
differences among model orbits while discarding event order within each orbit.
It is therefore a seasonal/inter-orbit adapter, not an accepted daily weather
sequence.

EFOR-1 through EFOR-8 replace this path with consecutive, interval-explicit
ExoPlaSim forcing produced before climatological averaging. The averaged
climatology remains the right artifact for maps and equilibrium summaries. It
is not the right artifact for daily interception, snow, drought, phenology,
respiration, fire or subdaily photosynthesis.

The producer half of that replacement exists: `NECO` turns on an ecological
output stream in the climate model, written at the interval an ecological
consumer integrates over, carrying the fields
`notes/ecological-forcing-field-contract.md` declares. See
"The ecological stream" in `exoplasim/README.md`. How the resulting blocks are
replayed through a spin-up much longer than any of them, what uncertainty that
construction carries, and what outcome would count as the construction failing
are fixed in advance in `notes/forcing-replay-preregistration.md`.

The adapter's cycling mechanism was verified with three artificial years at -3,
0 and +3 K on one cell. Annual NPP locks to the supplied period and responds
strongly:

| phase | offset | NPP kgC/m2 |
| --- | --- | --- |
| 0 | -3 K | 0.508 |
| 1 | 0 K | 0.377 |
| 2 | +3 K | 0.191 |

This is an interface response test, not a meaningful weather experiment. Two
things are still worth reading off it. The response is steep, nearly halving per 3 K,
because this cell is water-limited and warming raises evaporative demand. And the
mean over the cycle, 0.359, is **4.8% below** the value at the mean climate,
0.377. That is the concavity argument made concrete: running a variable star on
its average climate over-predicts productivity, and the error is one-sided.

### Correction: the cycle does not reach survival thresholds

An earlier version of this file argued that the cycle would cross PFT
cold-survival limits on the 13.8% of land lying within 2.1 K of one, and that
thresholds do not average. **That is wrong for this model**, and the reason is
worth knowing.

`tcmin_surv` is not compared against a cold year. `vegdynam.cpp:153` tests it
against `climate.mtemp_min20`, which `driver.cpp:617-627` builds as the **mean of
the last twenty years' coldest monthly means**. That is a deliberate choice in
LPJ-GUESS: one hard winter does not extirpate a simulated species, so the model
smooths before it removes a cohort. Twenty simulation years is 9.9 Earth years, which spans about
1.2 periods of an 8-Earth-year stellar cycle, so the window averages the cycle
almost exactly.

So the cycle reaches the vegetation through **productivity**, which responds
annually and is where the 4.8% concavity effect measured above lives, and not
through survival. Biome boundaries move because growth changes, not because
plants freeze.

The 13.8% figure is still the right measure of how much land is climatically
marginal. It is not a measure of what this model will do with it.

## Resolution: a convergence ladder, not a T42/T85 choice

T42 is the current operating point because it has been a useful balance of
execution time and resolved climate, not because the biosphere or pipeline is
intrinsically limited to it. The supported climate ladder is
T21/T42/T85/T127/T170, and ongoing ExoPlaSim optimization changes the cost side
of that choice. Historical T42/T85 timings are therefore measurements of old
runs, not a permanent selection rule.

LPJ-GUESS grid cells are laterally independent, so its direct cost grows roughly
with simulated land-cell count. That does **not** make finer support cosmetic.
Each cell receives a different climate, soil and surface environment, while the
higher-resolution atmosphere changes orography, coastlines, precipitation,
extremes and feedbacks rather than merely interpolating a T42 answer. Running a
fine LPJ grid on interpolated coarse forcing would add no information, but
running it on an accepted fine climate can change both spatial pattern and
extensive totals.

The reference under every climate grid is the ~10M-region Orogen export: 7.60 km
mean edge and full sampling of the generator's ~20 km terrain-information floor.
For scale, an average global climate cell contains about 4,883/1,221/305/136/76
native regions at T21/T42/T85/T127/T170. Even T170 therefore aggregates terrain,
soil and hydrologic variation; increasing truncation does not remove the need
for conservative partial areas and ecological response units.

No universal sign can be assigned to the NPP change. A concave response to one
forcing can make `response(mean forcing)` exceed the mean local response, but
real refinement changes variance, covariance, thresholds, coastline and the
coupled atmospheric solution simultaneously. The Miami land-mean versus
per-grid-cell calculation in `productivity-prediction.md` demonstrates an
aggregation effect over that particular support; it does not predict the sign
of T42-to-T170 coupled change.

Choose the production support by the pre-registered convergence protocol in
`spatial-support-ecological-aggregation-audit.md`. Each candidate inherits a
mapped initial state through CLIM-52's restart converter where useful, then
re-equilibrates terrain/climate, hydrology, soil, vegetation and feedback loops
on its own support. Compare accepted equilibria after conservative remapping to
one common area basis. Stop at the first rung sufficient for the declared
decisions; if material results do not converge by T170, retain a resolution
bracket rather than treating the highest answer as truth.

## Spin-up is in simulation years, and that halves it

`nyear_spinup 500` reads like an absolute statement and is not. A simulation year
is about half an Earth year, so the shipped 500 gives this world roughly half the
absolute vegetation and soil development that Earth practice assumes.

`build_vesper_pfts.py` therefore scales year *counts* UP by the reciprocal of the
factor it scales annual *sums* DOWN by, and converts fractions applied once a year
as rates rather than multiplying them. Confusing those directions would be worse
than doing neither, so the class of every parameter is named in
`notes/time-base-unit-contract.md`, in that script, in the generated
`generated/vesper_pfts.ins` and in its provenance. The classes:

| parameter | class | direction | why |
| --- | --- | --- | --- |
| `nyear_spinup` | year count | up | time to reach steady state |
| `distinterval` | year count | up | disturbance return time |
| `freenyears` | year count | up | time to build an N pool before N limits |
| `longevity` | year count | up | compared against an age that counts simulation years |
| `leaflong` | year count | up | leaf lifespan, and the SLA regression converts it back to absolute months |
| `gdd5min_est` | annual sum | down | degree-days accumulated in one year |
| `greff_min` | annual sum | down | annual production per unit leaf area |
| `turnover_leaf`, `turnover_root`, `turnover_sap` | annual rate | rate conversion | a fraction, so `1-(1-r)^f` and not a multiplier |
| `estinterval`, `est_max` | seasonal cycle | none | counted in growing seasons, not absolute time |
| `phengdd5ramp` | within-season | none | already absolute time; scaling it would be a real error |

Slow soil carbon does not need integrating for all of that: `ifcentury 1` solves
the equilibrium pool sizes analytically, accumulating running means between 70%
and 80% of the spin-up and solving at the end of that window
(`guess.h:3030-3034`). What the longer spin-up buys is enough absolute time for
the *vegetation* to reach steady state before that solve happens, which slow
forest succession needs.

The cost is a doubled run, 23 minutes to about 46 on 16 ranks. Not a
consideration.

## The 30-hour day widens the diurnal range, and neither side carries it yet

A 30-hour rotation gives longer daytime heating and longer nighttime cooling than
Earth's, so the simulated diurnal temperature range is wider, and on marginal
ground a night could dip below a freezing threshold that the daily mean never
approaches.

**The forcing does not carry it.** `maxt` and `mint` are timestep extrema and do
include the full 30-hour trough, but they are extrema of `dt(:,NLEP)`, the
SURFACE temperature, and not of the near-surface air temperature. On the
bootstrap climatology they bracket `ts` in every one of 24,576 cell-bins and fail
to bracket `tas` in 15,561 of them, by as much as 28.3 K. The air-temperature
extrema this climate model does compute are `atsama` and `atsami`, output codes
201 and 202; they are written by the model and reach no product, being absent
from pyburn's `ilibrary` and from `run_exoplasim.REGULAR_CODES`. `world-j0az` owns
delivering them.

So the driver's fourth array is the surface-temperature range under its own
name, `vesperinput.cpp` does not hand it to `climate.dtr`, and a run that asks
for `ifbvoc 1` is refused rather than given a different variable under the right
name.

**LPJ-GUESS would have nowhere to put it either.** `climate.dtr` is read in
exactly one place, `bvoc.cpp`'s `daytime_temp`, for leaf temperature in the
biogenic VOC scheme, which is off. There is no daily-minimum plant mortality
anywhere in the model; every cold limit runs through `mtemp_min20`, a twenty-year
mean of monthly means. So a wider diurnal range is physically real here and
radiatively present in the climate, and the vegetation model is structurally
blind to it. Representing it would mean adding a frost plant-mortality mechanism
LPJ-GUESS does not have, which is a much larger change than this project needs.

## Cost

5.46 s of CPU per gridcell for 530 years at `npatch 5`, measured
(`notes/lpj-guess-porting-audit.md`), so 4,106 land cells at T42 are about
6.2 CPU-hours and 23 minutes of wall on 16 ranks -- the resolution table
above. The shipped demo configuration runs nearer 24 s per cell; the
difference is the patch count. The vegetation loop is not the expensive half
of this iteration; ExoPlaSim is.

## Running it

Everything flux-dependent is generated, never written down, because the year
length is a function of `orbit.baseline_flux_earth` and moves whenever the flux
does. Regenerate all three after any orbit change, and rebuild: the year length
sizes arrays at compile time, so a stale binary is silently wrong rather than
failing. The driver file records the year length it was built for and
`vesperinput` refuses a mismatch, which is the backstop for exactly that.

```bash
python biosphere/scripts/build_vesper_header.py   # vesper.h, installed into the tree
python biosphere/scripts/build_vesper_pfts.py     # degree-day limits rescaled
python biosphere/scripts/build_lpj_driver.py      # climate + soil codes + gridlist
cmake --build vendor/lpj-guess/build --parallel 16
```

### The volatile organic source is off, and off is a decision

`bvoc_gate.py` is that decision made explicit. `ifbvoc 1` reads like one switch
and is four models -- the simulated plants' production of speciated volatile
carbon, its oxidation in an atmosphere whose oxidants this project has not
derived, the aerosol that oxidation would make, and what that aerosol would do
to the model's radiation and clouds -- and only the first exists in the vendored
source. `biosphere/config/bvoc.yaml` is where each precondition is declared, the
sentinel `undeclared` refuses, and there is no default value for any field in
it. The contract is `notes/bvoc-activation-contract.md` and the cloud arm is
pre-registered in `notes/bvoc-cloud-sensitivity-preregistration.md`.

```bash
python biosphere/scripts/bvoc_gate.py                       # what is undeclared
python biosphere/scripts/bvoc_gate.py --check-run runs/<id> # accept or reject output
```

`run_lpj_guess.py` asks the gate and writes `ifbvoc` into the instruction file
from the answer, rather than inheriting the imported PFT file's zero: an
inherited zero cannot be told apart from nobody having decided. Requesting
activation in the declaration while a precondition is undeclared makes the run
exit with every unmet one named.

### One abiotic nutrient ledger, and it does not close yet

`abiotic_nutrient_ledger.py` is the graph every abiotic nutrient source books
into: thirteen control volumes from the lithosphere to a terminal coastal
export, twenty-one directed terms, and one identity, that mass held plus mass
at terminal nodes equals mass held at the start plus mass across the boundary.
It reads `config/abiotic_nutrients.yaml`, and it can FAIL: seven reduced
fixtures run on every invocation, six built to be wrong in a named way, and a
fixture that does not get the verdict it was built for is a defect in the
checker rather than in the declaration.

```bash
python biosphere/scripts/abiotic_nutrient_ledger.py            # status, exit 0
python biosphere/scripts/abiotic_nutrient_ledger.py --strict   # refuses while a term is undeclared
```

Every term carries the `undeclared` sentinel today, so the ledger is DEFINED
and does not CLOSE, and the report names the issue that owns each one. The same
run carries the source screen -- what may book in, whether its rate can come
from the terrain at all, and whether pulse timing survives averaging -- and the
non-N/P adequacy screen, whose bounds and directions are declared before any
result is seen. The contract is `notes/abiotic-nutrient-ledger.md`.

### One mineral-reactivity contract, and its second arm refuses

Every mineral control on soil organic matter in the model is a linear function of
clay, or of clay plus silt, from Parton et al. (1993). The pedology soil map
already carries an andic areal fraction and an andic phosphate-fixation share,
and until now nothing read them. Both now cross the interface into
`Soiltype::andic_frac` and `Soiltype::p_fixation_frac`, with `UNSET_SOIL_FRAC`
where an input path has no andic state, and no equation reads either -- which is
what `mineral_reactivity_gate.py` checks, alongside every declared coefficient
against `somdynam.cpp` and every declared range over the texture simplex and over
the soil map's own textures.

```bash
python biosphere/scripts/mineral_reactivity_gate.py            # status, exit 0
python biosphere/scripts/mineral_reactivity_gate.py --strict   # refuses the
                                                               # mineral-aware arm
```

The mineral-aware arm needs four proxies and all four carry the `undeclared`
sentinel: Fe-Al oxide content, allophane concentration, aggregate capacity and
polyvalent cation saturation. The pedology `andic` column is the near miss and is
registered as a partial producer, because it is an areal fraction of andic
material rather than an allophane concentration. Neither arm may be retuned to a
desired soil carbon or productivity. `notes/mineral-reactivity-contract.md` is the
contract, and it carries the two things the texture arm does that are worth
knowing: its clay control on passive SOM formation is flat above a clay fraction
of 1/3, which a quarter of this world's land is above, and its microbial
partition can leave a negative transfer fraction on sandy soils, which this
world's soil map does not reach.

### The soil nitrogen transformation operator is declared, and bounded

`ifntransform 1` is on in the baseline, and it is not an emissions diagnostic:
nitrification, denitrification and ammonia volatilisation set how much mineral
nitrogen the simulated plants can reach. `ntransform_gate.py` is the declaration
of that operator checked against the source the model reads. It fails on a
constant that has drifted from `global_soiln.ins`, on a literal no longer in
`modules/ntransform.cpp`, on a response function that leaves the range its role
allows anywhere in its declared domain, on a chain of factors whose product
would take more nitrogen out of a pool than the pool holds, and on a chain
declared to be held inside its pool by an explicit `min()` that the operator no
longer contains. A dozen reduced fixtures run on every invocation, all but one
built to be wrong in a named way.

```bash
python biosphere/scripts/ntransform_gate.py            # status, exit 0
python biosphere/scripts/ntransform_gate.py --strict   # refuses while a Vesper
                                                       # precondition is undeclared
```

`biosphere/config/ntransform.yaml` is the declaration and there is no default for
any precondition in it. All five carry the `undeclared` sentinel today: surface
pressure, oxygen partial pressure, soil gas diffusivity, water-table redox state
and the atmospheric boundary for the gases the operator emits. So the declared
model boundary is that any nitrogen-limitation result from it is a result for an
Earth gas and redox environment driven by this world's water and pH.
`notes/soil-nitrogen-transformation-parameterisation.md` is the register, and it
carries three defects the operator had: the soil map's pH never reached it, its
no-pH fallback ran on a variable nothing assigns, and its only conservation check
was an `assert` that Release compiles out.

### Respiration acclimates to a growth temperature, or not at all

`respiration_acclimated()` replaces each simulated plant functional type's
`respcoeff` with a function of the GROWTH temperature its tissue has adjusted to.
It was being handed the current day's air and 25 cm soil temperature, so the same
variable drove the acute Lloyd and Taylor response and the acclimation multiplier
and no acclimation was represented. `Climate::tacc_air` and `Soil::tacc_root`
carry the memory now: exponential running means with an e-folding time of
`acclim_resp_tau` absolute days, seeded from the first temperature each gridcell
sees and serialized so a resume keeps them.

That e-folding time has no value this project can derive, so there is no default
and `parameters.cpp` refuses `acclimated_respiration 1` without it. The baseline
takes the standard respiration path, which is what the CNP fork's own
`global_p.ins` selects and what divides `respcoeff` by the tissue C:N windows.
`acclimation_gate.py` enforces the state, reports what the memory is worth over
this world's seasonal cycle across the bracket the held literature supports, and
refuses an acclimated path with no declared memory length.

```bash
python biosphere/scripts/acclimation_gate.py            # status, exit 0
python biosphere/scripts/acclimation_gate.py --strict   # refuses an acclimated
                                                        # path with no memory
```

### The seasonal landmarks are derived, not dated

Summergreen phenology turns on two days of the simulation year: the coldest,
where the growing-degree-day sum and the annual leaf-on sum reset and chilling
detection switches off, and the warmest, where chilling detection switches back
on. LPJ-GUESS carried both as fixed ordinal dates on the Earth calendar, one pair
per hemisphere, and this world's year is shorter than half an Earth year, so the
southern date named a day the calendar does not have.

`Climate::coldest_day` and `Climate::warmest_day` replace them, read off a
running day-of-year mean of the air temperature forcing and smoothed over one of
this world's months. No hemisphere test and no thermal-lag assumption enters, and
the pair can never collapse onto one day, so each reset fires exactly once per
orbit on every gridcell including the equator. `VesperInput` hands the whole
interpolated year over at day 0, so the landmarks are the cell's own from the
first orbit and follow a driver file of several years as it cycles.

`phenology_gate.py` is the enforcement. It fails on an Earth ordinal date back in
any compiled source, on a reader no longer keyed on the derived pair, on the
chill-day count being able to leave the `Pft::gdd0` table it indexes, and on a
landmark field missing from `Climate::serialize`. Eleven fixtures run on every
invocation, four of them built to be wrong in a named way.

```bash
python biosphere/scripts/phenology_gate.py            # status, exit 0
python biosphere/scripts/phenology_gate.py --strict   # refuses on the one
                                                      # natural-vegetation event
                                                      # still on an Earth calendar
```

### Fire is GLOBFIRM, with flux and occurrence diagnostics

`run_lpj_guess.py` writes `firemodel "GLOBFIRM"` after its `import` of
`vesper_pfts.ins`, which still carries the shipped `firemodel "BLAZE"`. The later
declaration wins, so the generated instruction file is what decides this and the
PFT file is not worth editing. BLAZE wants a SimFIRE input built from Earth
observations, and it forces `weathergenerator "GWGEN"`, which wants sub-daily
statistics this world does not have; the shipped demo makes the same two
substitutions for the same reason. The mismatch is not silent: LPJ-GUESS aborts
on BLAZE with a non-GWGEN generator, so a run that starts is a run on GLOBFIRM.

GLOBFIRM writes its carbon loss in the `Fire` column of `cflux.out` and its
inferred return time plus burned fraction in `firert.out`. The harness retains
both. BLAZE has the richer daily/monthly burned-area and SIMFIRE analysis
outputs, but those are separate from GLOBFIRM's annual occurrence diagnostic.

One gate worth knowing when a short run reports no fire at all. With
`iftwolayersoil 0`, which `iforganicsoilproperties` requires,
`vegdynam.cpp:1386` accumulates fire-season days only after model year 100, from
a hard-coded constant that is not rescaled alongside `nyear_spinup` and
`distinterval`. The spin-up is far longer than that, so it never reaches a
reported year, but it does mean fire is off for the first tenth of it.

## One-off tools

- `scripts/score_prediction.py` -- one-off: scores a productivity prediction against an LPJ-GUESS run, the machinery behind BIO-2's nitrogen bracket. Registered under `one_offs` in `config/pipeline.yaml`; it generates nothing the pipeline reads.
- `scripts/check_forcing_contract.py` -- one-off: a verdict on a climate product against `notes/ecological-forcing-field-contract.md`. Runs the closure, sign and bracketing identities, and carries seven reduced fixtures, six of them wrong in a named way, so it can fail on itself. Registered under `one_offs`; it generates nothing.
