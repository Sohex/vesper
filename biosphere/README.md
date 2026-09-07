# Biosphere

> Figures in this document are illustrative of method, and were measured on
> builds, climates and model source that have since moved. Current values live
> in `world_state.json`, generated from the artifacts. See the convention in
> `CLAUDE.md`.



Replaces the assumed vegetation with a modelled one. LPJ-GUESS takes the
ExoPlaSim climatology and returns leaf area, carbon and plant functional type
composition per gridcell, which becomes the surface albedo and forest fraction
that the next climate run is forced with.

BIO-11's shared area interface is
`scripts/build_rootable_fraction.py`. It derives, per build and atmosphere rung,
the mutually exclusive rootable, solved-water and dry-barren shares of native
mesh land, reading the authoritative solved-lake mask directly so iteration-0
soil does not depend on the post-soil surface classifier. `build_lpj_driver.py` omits fully non-rootable cells; pedology's soil
carbon feedback and `scripts/score_prediction.py` apply the same fraction to
every extensive quantity and denominator. The artifact is support- and
provenance-stamped and has no nearest-cell fallback.

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
`notes/phosphorus-cycle-parameterisation.md`, the three per-cell soil phosphorus
inputs the fork reads and what this world can put in each of them are in
`notes/soil-phosphorus-input-parameterisation.md`, plant physiology and carbon
allocation are audited in `notes/plant-physiology-carbon-allocation-audit.md`,
BVOC emissions, secondary organic aerosol and atmospheric coupling are audited
in `notes/bvoc-soa-atmospheric-coupling-audit.md`, wetlands, peat and methane
are audited in `notes/wetlands-peat-methane-audit.md` and what has to be
declared and repaired before those switches move is
`notes/wetland-activation-contract.md`, and
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
| input module | `vesperinput`, runs end to end and splits across MPI ranks. Its `VESPDRV8` transport is CHRONOLOGICAL: the file carries a table of forcing intervals with explicit bounds and duration in absolute seconds and the local solar phase of each, any count, and `integrate_year` takes each absolute day's duration-weighted mean over the intervals overlapping it. So the 24-hour hydrology and biogeochemistry boundary is the consumer's rather than the format's, twelve intervals a year and one per timestep are the same code path, and the smooth curve `interp_monthly_means_conserve` manufactured between bin centres is gone because the producer states an interval mean and says nothing about the shape inside it. `config/ecological_forcing_contract.yaml` is what it has to carry; EFOR-3 and EFOR-8 own the artifact that will replace the binary |
| soil and water | pedology depth scales LPJ capacity and pedology AWC sets ExoPlaSim's scalar bucket at smoke scale, but the models independently derive hydraulic properties and run separate snow/soil water balances; LSHY-1 through LSHY-7 own the consistency work |
| abiotic nutrients | the ledger is DEFINED and does not CLOSE: `abiotic_nutrient_ledger.py` carries thirteen control volumes and twenty-one terms, every one of them still holding the `undeclared` sentinel with its owning issue named. Rock P and dust mass have useful relative/source artifacts but no absolute flux; `phosphorus_budget.py` does not consume dust deposition, and ANUT-2 through ANUT-6 and ANUT-10 own the terms that would close it |
| abiotic source screen | geomorphic renewal, arc tephra and marine aerosol are RETAINED against the ledger, volcanic sulfate deposition is registered and not implemented, and fire ash and lightning belong to FIRE-7 and ANUT-4. The exhumation and tephra rates come from Earth's stationary population and not from the terrain; no screen may be carried on an aerosol optical depth |
| non-N/P adequacy | screened as a critical runoff per element and per lithology, bounds and directions declared before the result. Potassium binds; iron and the trace set REFUSE for want of a release table. The declared model boundary is that any LPJ-GUESS result here is a C-N-P result and not a nutrient-limitation result |
| tissue stoichiometry | the fine-root and sapwood C:N and C:P windows are anchored on the tissue MEAN, so the proportion `canexch.cpp` applies to their nutrient demand is the one Friend et al. (1997) measured. The max-anchored form the model used to carry applied 1.79 times it on nitrogen and 2.22 on phosphorus, so every C-N number produced before the repair is worthless rather than stale and the biosphere needs re-commissioning. `notes/plant-physiology-carbon-allocation-audit.md` finding 11 |
| phosphorus parameters | every constant registered with its source or its bracket in `notes/phosphorus-cycle-parameterisation.md`; the uptake profile, the leaf C:P window, the root proportion, the labile-P saturation threshold and the surface microbial pool's C:P are derived, the sapwood proportion is not on either its level or its shape, and `ifplim 1` fails closed naming both halves of it so that settling either alone does not lift the refusal. Three DECLARED DIVERGENCES from the vendored CNP phosphorus values: the labile-P threshold is converted into the fork's own Hedley-labile currency, which is inert under `ifplim 0` because the P-limitation-off pin reads the same constant; the surface humus pool ramps on the slow pool's C:P line instead of holding a fixed unsourced ratio, under the identification the fork's own nitrogen ramp already makes; and the surface microbial pool, the one C:P nothing overwrites, holds its own stoichiometry, its sourced C:N of 20 times a measured decomposer N:P, instead of the active soil pool's endpoint. The sapwood refusal is about the element, the level and the shape: the model applies the constant it declares, and the constant is one scalar where the quantity varies with temperature and plant functional type |
| phosphorus sinks | leaching, fire and harvest only. Terminal occlusion is a DECLARED ABSENCE, argued in the same note, so a simulated soil that must be old carries its phosphorus depletion in its initial stocks rather than developing it |
| run harness | written; `lpj_run` records inputs, binary, model identity, and DEMO-5's stochastic root plus substream ABI, then `lpj_acceptance` runs `assess_lpj_run.py` to require exact rank/output/cell/year coverage, finite physical values, stable end windows and C/N/water closure before the run can be consumed |
| equilibrium acceptance | `vesper-lpj-equilibrium-window/6` bounds the END-TO-END RELATIVE DRIFT of THE QUANTITIES A CONSUMER READS from above over the whole retained record and passes each when the bound is inside ITS OWN tolerance, so what the contract controls is the rate at which a run that IS drifting is accepted. A run passes only when every quantity does, so by intersection-union that rate needs no multiplicity correction -- which is what the empirical per-window null could not deliver, its family rate over the 64 columns it then ran on having been measured at 0.34 against a declared 0.05. The assessed set is the quantities a consumer reads, so a reader that imposes no drift tolerance -- the closure's conservation identities -- puts none here and its columns keep reported bounds as diagnostics. A record too short to resolve a tolerance gives a wide bound and is refused by the test itself; there is no separate span guard, and the length that answers a refusal is written onto the acceptance artifact whatever the verdict. The per-cell half is the one part whose LIMIT is still calibrated rather than derived, and it errs toward refusing (`world-4hlw`); what a settled run pays through it is derived per run onto the acceptance artifact as a union bound over the assessed set rather than declared, because a family rate written beside a list that gains a row stops being true. `notes/equilibrium-trend-null.md` has the measurements; `scripts/validate_drift_statistic.py` re-takes the statistic's size and cost and `scripts/derive_trend_null.py --timescales` the timescales |
| what a consumer reads | the reduced value, its temporal spread and its own relative standard error, all over the WHOLE RETAINED RECORD -- the span the contract certifies. A ten-cycle mean of a series whose memory time is tens to hundreds of cycles is one effective sample, and ten of this model's columns have a marginal scatter above the drift tolerance outright, so reporting one would hand a consumer a number less certain than the certificate it travels with. `complete_forcing_cycles` now names the per-cell half's window alone |
| run lengths | `lib/run_lengths.py` carries a second timescale pair for this model, in complete forcing cycles and never interchangeable with the climate pair in orbits, and the two lengths are INDEPENDENT. The RETAINED RECORD is the length at which the acceptance contract's own drift bound closes on its tolerance, at the assessed quantity that needs the most. The SPIN-UP is derived WITHOUT a relaxation time: the requirement it inverts is bounded over all relaxation times, so its supremum -- a multiple of the retained record fixed by the tolerance alone -- satisfies it at every one and is the smallest length that does. The relaxation estimator declines every field-record this project has, because a contraction whose uncertainty reaches one measures the record and not the timescale, and that no longer blocks anything: a measured relaxation time could only ask for less. `notes/equilibrium-trend-null.md` has both derivations and the throughput the cost follows from |
| albedo and forest feedback | modelled mode exists; rootable/lake and spectral corrections are BIO-17 and BIO-18 |
| aerodynamic feedback | modelled roughness is open as BIO-16 |
| BVOC/SOA feedback | LPJ source is present but off; carbon closure, PFT traits, reduced atmospheric chemistry/transport, direct optics and cloud effects are separated under BVOC-1 through BVOC-10 |
| wetlands/peat/methane | LPJ's northern-Earth peat/CH4 source is present but off; current low-latitude saturation creates water and the full extent, groundwater, peat-stock, source/sink and atmospheric closure is separated under WET-1 through WET-11 |
| full run | harness and BIO-14 acceptance exist, but BIO-15 and BIO-21 through BIO-25 must close before its output is interpreted; prior estimate ~25 min on 16 ranks |

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

**The photon conversion.** A photosystem counts quanta and the radiation scheme
carries joules, so `FRADPAR` is only half of the conversion: `canexch.cpp` turns
the energy inside the window into a photon supply with `CQ`, which LPJ-GUESS
ships as the monochromatic 550 nm value for the Sun. `VESPER_CQ` is derived
beside `FRADPAR` from the same spectrum over the same window, so the pair is one
currency. `lib/stellar.py` owns both integrals and runs its controls -- the
monochromatic 550 nm identity, a solar spectrum against the shipped constant,
and the agreement of the wavelength and frequency routes -- before returning
either.

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
python biosphere/scripts/build_rootable_fraction.py # BIO-11 effective plant area
python biosphere/scripts/build_lpj_driver.py      # climate + soil codes + gridlist
qrun -p build -- \
  .venv/bin/python biosphere/scripts/build_lpj_guess.py   # compile, and record what from
```

After a run is accepted, one step turns its tables into a field other components
and the maps can read:

```bash
qrun -p light -- \
  .venv/bin/python biosphere/scripts/build_vegetation_field.py
```

It writes `biosphere/data/<build>/vegetation_<rung>.nc`: per plant functional
type foliar projective cover, leaf area index and vegetation carbon on the
atmosphere grid, each reduced by `lib/lpj_output.py:reduce_table` over the span
BIO-12's equilibrium contract certifies. A run whose retained record does not
bound its own drift emits nothing rather than a mean over a window that means
nothing, and a table absent from or changed since the acceptance artifact is
refused outright.

**The rows carry the atmosphere model's own longitude labels**, wrapped into
-180..180 by `build_lpj_driver.py` because that is the range LPJ-GUESS reads.
Wrapping a label renames a column rather than moving it, so the placement goes
through `gridding.model_label_cells`, which refuses a label off that axis
instead of rounding it to the nearest column. The export's centres would sit
half a column from every label; CLAUDE.md rule 3. `maps/build_basemap.py` is
the first consumer, and it draws the land from this rather than from the
climate classification wherever the field exists.

### Where a run's output lands, and where an instrumented build's does NOT

Three streams, and only one of them carries what a probe prints:

| | written by | carries |
| --- | --- | --- |
| `runs/<id>/run*/guess.log` | the model, per rank | its own progress and its own diagnostics |
| the terminal | `run_lpj_guess.py` | the runner's own reporting |
| **`runs/<id>/mpirun.log`** | the runner | **the model's stdout AND stderr, combined** |

`run_lpj_guess.py` launches the model under
`subprocess.run(..., capture_output=True)`, so an `fprintf(stderr, ...)` added
to the compiled model reaches `mpirun.log` and NOWHERE ELSE. The two obvious
places to look -- the rank logs and the terminal -- are silent by construction,
and nothing is broken when they are.

**A probe therefore carries a control that must fire.** Print one unconditional
line at the instrumented site, confirm it appears in `mpirun.log`, and only then
believe a zero count. Four instrumented builds on `world-n0oc` reported zero
events, all four were false negatives, and two wrong mechanisms were struck off
as measured eliminations before an unconditional control was added and did not
fire either. `docs/src/practice/failure-modes.md` class 41.

**Build the model with `build_lpj_guess.py` and not with a bare `cmake --build`.**
It runs the same two cmake commands and then writes
`vendor/lpj-guess/build/guess.provenance.json`, recording the executable's sha
beside the sha of every one of the 124 files it was compiled from -- the four
subdirectory CMakeLists' declared lists, everything a quoted `#include` reaches
from them (the generated `framework/vesper.h` among them), the build files
themselves, and the compiler and cache entries that turned that source into
those bytes. `--verify` compares the record against the tree, building nothing,
and `run_lpj_guess.py` runs the same comparison and REFUSES a binary the tree
has moved under.

This is rule 4 for the biosphere. LPJ-GUESS reads none of its source at run
time, so a binary built before an edit integrates the code it was built from and
reports nothing; the one time that was caught, it was caught only because the
new parameter reached the model through an instruction file and the parser
refused an undefined identifier. A change confined to C++ has no parser in front
of it. A hand-typed `cmake --build` leaves the record describing the previous
executable, which the next run refuses -- so the failure mode of forgetting this
script is a refusal, not a wrong number.

The build replaces the executable a run in flight has already hashed into its
manifest, so check `pgrep -x guess` before starting one.

`build_lpj_driver.py --self-test` runs the interval arithmetic and header layout
fixtures and exits: no climatology, no soil map, no model. It covers the
operator `vesperinput.cpp:integrate_year` implements and the byte layout that
module parses. LPJ-GUESS now builds and runs on the generated 183-day header;
the self-test remains the no-model interface check that can run before a
climatology or MPI launch exists.

All three land in `biosphere/generated/`, along with the gate and ledger
reports. That directory is output and is not tracked: everything in it is
re-derived by a step registered in `config/pipeline.yaml`, and what each gate
argues is in `biosphere/notes/`, which is.

`stochastic_seeds.yaml` declares the one random root. The run harness writes it
to the instruction and manifest; the model hashes it with cell coordinates,
stand, replicate patch and a named process id. Establishment, fire occurrence,
fire mortality, background mortality and disturbance therefore have separate
restart-serialized streams, and neither MPI rank nor traversal order is part of
the key. `python biosphere/scripts/stochastic_seed_gate.py` exercises fixed
vectors, rank/order reshuffling and source/serializer coverage without a model
run. `--root-seed` creates a distinct recorded replicate rather than an
unattributed change in results.

No consumer reads the final output year. `equilibrium_window.yaml` fixes one
rule for albedo feedback, prediction scoring and soil carbon: average ten
complete forcing cycles, after reducing each cycle to one mean, and refuse an
end window whose cycle means still trend or whose cell-year rows are
incomplete. `equilibrium_window_gate.py` exercises the shared reader without a
simulation. Every consumer report carries the temporal standard deviation and
labels root-seed and patch-count uncertainty `not_measured` until comparable
runs are supplied through its repeatable peer option; one realization is never
reported as zero stochastic uncertainty.

Modelled canopy is composited on the same BIO-11 partition. LPJ tree and grass
FPC are conditional on rootable ground; `build_surface_albedo.py` multiplies
them by `f_rootable` once, subtracts only the rootable substrate contribution,
and leaves solved-water and dry-barren contributions unchanged. Forest code 212
is the resulting whole-cell tree share. `rootable_albedo_gate.py` fixes ordinary,
barren, partial-lake and pure-lake answers and checks the two-band identity.

### A run saves its state, and another run continues from it

The ecological spin-up in front of a retained record is several times the
record, so most of a run's wall clock buys the approach to equilibrium rather
than the record itself. A run refused by the acceptance contract used to be
answerable only by a second run from bare ground that re-integrated the whole of
that approach; with a state file the marginal cost of more retained record, of a
second seed, or of another patch count is the record alone.

```bash
python biosphere/scripts/run_lpj_guess.py --nyear 1253 --save-state
python biosphere/scripts/run_lpj_guess.py --nyear 1253 --continue-from lpj_<parent>
```

**Every instant in a state file is a `date.year`, and `date.year` counts from
zero THROUGH the spin-up.** `--nyear` is the retained record; `nyear_spinup`, the
derived floor `build_vesper_pfts.py` writes into the PFT file, is the spin-up in
front of it. A run therefore ends at simulated year `nyear_spinup + nyear - 1`,
`vesperinput.cpp:685` stops it there and `commonoutput.cpp:785` writes no annual
row before `nyear_spinup`. Nothing that computes a save point or a resume point
is right without it, and `run_lpj_guess.py` writes `nyear_spinup` into the
instruction file from the value it read out of the PFT file, so the number it did
its arithmetic with is the number the model runs.

The continuation integrates `--nyear` MORE years, all of them retained, and
declares its parent's total as its own `nyear_spinup`. That is what makes its
output cover its own years and no others; the CENTURY accelerator window derives
from the same parameter in the Soil constructor, but it is serialized, so the
deserializer restores the parent's window over the constructor's and the
accelerator does not fire again. The continuation is a NEW run with a new id --
rule 6, unchanged by where the state came from -- and its record is its own
years, the parent's being spin-up whatever they were bought as.
`run_manifest.json` carries a `continuation` block naming the parent, the chain
back to bare ground, the simulated year it resumed at, how many years it
integrated itself and the hashes of the state files it read. `acceptance.json`
carries that block too, on a pass and on a refusal, so a consumer reading a
record can see that the spin-up in front of it was integrated by a named run
rather than assumed.

**What may be continued is guarded, not trusted.** A state file is a simulated
state and not a result: read under a different forcing, a different soil map, a
different binary or a different ecological parameter, the model integrates from
a state those inputs never produced while every number out of it is attributed
to a spin-up that did not happen -- a silent wrong answer, and a worse position
than having no state file at all. `run_lpj_guess.py:continuation_refusals`
compares the parent manifest's pinned input hashes and physical settings against
the continuing run's and refuses on any difference. The rank count is recorded
rather than refused: `PartitionedMapDeserializer` searches every file in the
state directory for a cell's coordinates, and the stochastic substreams are
derived from the cell rather than from the rank.

**A continuation is refused until the model is shown to reproduce the run it
continues.** That is a separate question from the mechanism and it has a right
answer. The first measurement said NO -- every retained table differed, over 30
retained years behind a 201-year spin-up at npatch 5 split at simulated year
211 -- and the cause was not in the serializer at all: `libraries/plib` rounded
the year-boundary sentinel `state_day -1` to 0, so every restart taken here was
an arbitrary-day restart at day 0 of `state_year` and lost the annual
accumulators outside Soil that day 0 resets. Both are repaired, and the same
measurement now finds no differing row.
`biosphere/notes/restart-state-outside-soil.md` has the finding, and every mode
of the fixture now asks the model which restart instants it PARSED rather than
restating them. `run_lpj_guess.py:continuity_verdict` reads the
fixture's report, checks it was taken against the binary about to run, and
refuses `--continue-from` while it says otherwise or is absent -- an absent
measurement is not a pass. `--save-state` is deliberately not gated: writing a
state file changes no number in the run that writes it, and only reading one
can.

`verify_lpj_restart_continuity.py` is where that question is asked. Its
model modes and what each can see are below, under what a restarted soil
column inherits; `--self-test` needs no model, holding the two integers the
runner decides against the arithmetic `framework/framework.cpp` performs on
them and exercising the parsed-instant gate on cases whose verdict is known. That test exists because the first version of `--save-state`
computed its save point from `nyear` alone and named a simulated year thousands
of years before the end of the run, so it reads the spin-up out of the PFT file
the runs import rather than carrying a copy of the number.

The fixture takes its own `--nyear-spinup`, and that is what makes it runnable:
inheriting the derived floor turns a twelve-year bed into a twelve-thousand-year
one, which is why it had never been run. It cannot go arbitrarily low. The model
refuses a spin-up at or below `freenyears`, which the generated PFT file
declares as 200 rather than the 100 the shipped instruction files carry, so the
fixture reads that value and refuses first, naming it -- three verification
attempts died one second in on preconditions only plib was checking, each
costing a host lock acquisition. A bed just above `freenyears` also puts the
split where the simulated plants are nitrogen limited, so the nitrogen pools
feeding that limitation are live state at the instant the state file is written.

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

### The wetlands, their peat and their methane are off, and off is a decision

`wetland_gate.py` is that decision made explicit, and it refuses on three
grounds rather than one. `run_peatland 1` plus `ifmethane 1` reads like two switches and
is four models -- where the simulated wetlands are, how water reaches and leaves
them, how peat carbon and redox make and consume CH4, and what an atmosphere
does with the flux. `biosphere/config/wetlands.yaml` is where each precondition
is declared and the sentinel `undeclared` refuses, on the same terms as the BVOC
declaration; the contract is `notes/wetland-activation-contract.md`.

The second ground is the vendored source itself. The gate probes it on every
invocation and refuses while the low-latitude wetland path still adds water the
simulated world never received, while runon is still a namelist scalar carrying
no catchment, while the annual water-table average is still guarded on an
ordinal `Date::next()` never produces, and while the prognostic peat hydrology
is still absent from `Soil::serialize`. A declaration claiming one of those is
closed while the probe still finds it is reported as a contradiction rather than
believed. A fixture set runs every time: mutations of the declaration built to
be wrong in a named way, the live declaration which must be refused, three
mutations of the EVIDENCE so that the contradiction branch stays checkable after
the repair that made the source agree, and a met declaration against a repaired
source which must be granted, because a gate nothing can satisfy is a wall
refusing for a reason nobody wrote down.

The third ground is what the fork can EMIT. An accepted run has to retain
`acceptance.retained_outputs` together, and each of those tables either has a
quantity behind it in the fork today or waits on a model rather than on an
output routine; `acceptance.retained_output_status` records which, and which
model each waiting table waits on. The split is not written down anywhere: the
gate re-derives it against `modules/commonoutput.cpp` on every run, prints it
and puts it on the report as `retained_output_split`, because a count beside a
list that is appended to stops being true without anything objecting.
`run_lpj_guess.py` derives one instruction row from each retained filename's
stem, so a table with no declared parameter in `modules/commonoutput.cpp` does
not produce an empty file -- it aborts the run while plib parses the instruction
file, naming one unknown parameter and nothing about why it is unknown. The gate
refuses ahead of that, naming every missing table and its owning issue, and it
checks the reverse direction as well, because a table that has become available
and is still recorded as waiting is one nobody will think to ask for.

```bash
python biosphere/scripts/wetland_gate.py                       # what is undeclared
python biosphere/scripts/wetland_gate.py --check-run runs/<id> # accept or reject output
```

`verify_lpj_restart_continuity.py` is the behavioural half of the last of those
four, where the two serializer probes are the static half. It runs the same
forcing twice and requires the resumed run to reproduce the uninterrupted one.
The default splits at a simulated year boundary and compares the output tables
from the restart point on. `--one-day` splits at an arbitrary simulated DAY and
compares the two runs' state files one day later, so state the serializer drops
appears as bytes that differ after a single day rather than as a year of
divergence. `--round-trip` writes the state again with no simulated day in
between, which names what the write and read of a state file does not carry --
the form that found `world-8yyh` in the climate model. Those two exist because
WORLD-FUJ4 gave `framework.cpp` a save point the caller can place on any
simulated day, a restart that resumes on it, and permission to do both in one
run; before that it serialized exactly once, at the end of year
`state_year - 1`. `--runner` is a fourth mode with a different subject: it
invokes `run_lpj_guess.py` for a parent that saves, a continuation from it, an
uninterrupted control twice as long, and a continuation at a patch count the
parent never ran, which has to be refused. The first three test the MODEL's
serializer and the fourth tests the RUNNER's plumbing. All four need a compiled
model and a built forcing and exit naming what is missing rather than reporting
a pass they did not earn.

`--cells` slices the driver to named lon,lat pairs. Cells are independent in
LPJ-GUESS and this project's stochastic streams are keyed by coordinate, so a
cell integrates the same trajectory in a two-cell driver as in the whole grid;
what a subset loses is reach, so its verdict is written to
`lpj_restart_continuity_subset.json` and the continuation gate in
`run_lpj_guess.py` keeps reading only the whole-grid report. `--runner` requires
it, because the runner simulates every cell in the driver it is handed.

Every mode also asks the model which restart instants it PARSED and refuses when
they differ from what the arm asked for. `state_day -1` and `save_day -1` are
the year-boundary sentinel and `libraries/plib` delivered 0 for both, so every
restart taken here was an arbitrary-day restart at day 0 of `state_year` while
the instruction file, the runner, the fixture and `--self-test` all reported a
year boundary; four documents agreeing with four copies of the same rule.
`biosphere/notes/restart-state-outside-soil.md` has the finding.

```bash
python biosphere/scripts/verify_lpj_restart_continuity.py --nyear 12 --state-year 8
python biosphere/scripts/verify_lpj_restart_continuity.py --one-day --state-year 8 --state-day 120
python biosphere/scripts/verify_lpj_restart_continuity.py --round-trip --state-year 8 --state-day 120
python biosphere/scripts/verify_lpj_restart_continuity.py --cells "11.25,19.38" --nyear 30 --nyear-spinup 201 --state-year 211
python biosphere/scripts/verify_lpj_restart_continuity.py --runner --cells "11.25,19.38;0.00,85.76" --nyear 30 --nyear-spinup 201 --state-year 211
```

### What a restarted soil column inherits

`soil_restart_state_gate.py` holds the other static half, and it covers the
whole Soil class rather than the peat hydrology alone. `Soil::serialize` is
everything a resumed column gets, and a member the class declares that the
serializer does not stream is either state the restart drops or something the
model rebuilds before it reads it -- indistinguishable in a diff, and the
difference is the whole question. `biosphere/config/soil_restart_state.yaml`
classifies every absent member as rebuilt-before-first-read, diagnostic, or
lost, each with the file:line that settles it, and the gate parses both the
class body and the serialize block and refuses on a member nothing accounts for,
on a classification the serializer now contradicts, on a name the class has
dropped, on a classification carrying no evidence, and on a non-empty `lost`
block. Six fixtures built to be wrong in a named way run every time, plus the
live declaration, which must be granted.

Writing the rebuilt ones down is the point. The sweep behind the current
classification found a handful worth carrying and sixty-odd that would have been
state-file bytes changing no result, and without the record the next sweep
re-derives all of them. A member added to the Soil class lands unclassified and
stays there until someone decides what it is.
`biosphere/notes/soil-restart-state.md` argues the classes and carries the
findings.

```bash
python biosphere/scripts/soil_restart_state_gate.py
python biosphere/scripts/soil_restart_state_gate.py --list-unclassified
```

`run_lpj_guess.py` asks the gate and writes all four switches into the
instruction file from the answer. The extent is a partition, so the gate takes
one source per class from a closed set and checks each named artifact against
the builds under `hydrography/data`; the saturated non-inundated mineral class
has none, because the saturated-area closure that would have supplied a
saturated fraction is withdrawn in
`hydrography/config/topographic_index.yaml`, and the gate reads that file rather
than carrying a copy of the verdict, so a revival refuses here too. That class
is DECLARED ABSENT under WET-12's reduced form rather than left blank, and the
gate refuses a class dropped without both the reason and the row that licensed
it; `biosphere/notes/reduced-wetland-form.md` is the reduced form and the ledger
of what taking it costs in claims. A saturated fraction is a grid-cell quantity
and the gate refuses any other support for it, which is WORLD-D9U4's constraint
made enforceable. A class whose share is a convention bracket propagates as both
arms or as neither, declared before any share is written, because `f_grad`'s
upper arm carries no terrain information and one arm of it is the depth under
another name. Peat age and depth are two-ended brackets and the gate refuses a
scalar, because a scalar is a claim to a history `no-time-axis.md` says this
world does not have.

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
longer contains. It also fails on a declared divergence from mainline
LPJ-GUESS whose mainline form the operator no longer records beside the changed
one, or which the operator has gone back to running: both halves are checked, so
the divergence can neither become a silent fork nor be silently reverted.
Reduced fixtures run on every invocation, all but one built to be wrong in a
named way.

```bash
python biosphere/scripts/ntransform_gate.py            # status, exit 0
python biosphere/scripts/ntransform_gate.py --strict   # refuses while a Vesper
                                                       # precondition is undeclared
```

WHETHER A DIVERGENCE HAS BEEN EXECUTED IS DECLARED PER ENTRY. The register's
`execution_arms` names the matched arms that have run, and every entry claims
one by name or says `none` and why; the gate refuses an entry that says nothing,
one that claims an arm the register does not carry, and an arm no entry claims.
So a divergence appended after an arm ran cannot inherit its verdict, which is
what a register-wide boolean did. `labile_carbon_microbial_share` is the entry
that claims none: it postdates the arm, and it splits the labile carbon the
operator denitrifies on into the two paths Li et al. (1992) states, so
`global_soiln.ins` carries a coefficient per path and the source half is
`labile_carbon_paths` in `somdynam.yaml`. `prepare_stock_ntransform_arm.py` copies the active
source tree, replaces only `modules/ntransform.cpp` with the exact release
object held in git, verifies its digest, and builds a provenance sidecar.
`run_lpj_guess.py --ntransform-profile stock-4.1.1 --binary <guess>` also
restores the release `f_nitri_gas_max` while keeping every other instruction
and input matched. `compare_ntransform_arms.py <vesper-run> <stock-run>` rejects
input drift and compares paired cell-years over the final ten forcing cycles.
The recorded arm changes area-weighted plant mineral-N uptake by +0.343%, but
substantially changes internal transformations and gas partitioning; the exact
runs and quantities are pinned in the arm's own block in the declaration, beside
the `standing` that says what its report is worth. Both arms fail the same
equilibrium refusal and `world-qcse` says why no run this project has made could
pass one, so this is a matched sensitivity bound between two refused runs rather
than an accepted coupled biosphere endpoint.

`biosphere/config/ntransform.yaml` is the declaration and there is no default for
any precondition in it. All five carry the `undeclared` sentinel today: surface
pressure, oxygen partial pressure, soil gas diffusivity, water-table redox state
and the atmospheric boundary for the gases the operator emits. So the declared
model boundary is that any nitrogen-limitation result from it is a result for an
Earth gas and redox environment driven by this world's water and pH.

A calibration entry can carry a `boundary` line of its own, and one does: the
denitrification N2 share, where Xu-Ri and Prentice (2008) and Weier et al.
(1993) disagree, neither is this world's, and the operator runs Xu-Ri's
partition because the partition is not separable from the reduction sequence it
divides. `--strict` does not refuse on a declared boundary, the gate names every
one under its own heading on each invocation, and a `boundary` on an entry the
sources do settle is itself a gate failure. What it declares is that no N2, N2O
or NO number this model reports is a prediction of this world's gas partition.
`notes/soil-nitrogen-transformation-parameterisation.md` is the register, and it
carries three defects the operator had: the soil map's pH never reached it, its
no-pH fallback ran on a variable nothing assigns, and its only conservation check
was an `assert` that Release compiles out.

### `modules/somdynam.cpp`'s divergences from the CNP fork are registered

`somdynam_gate.py` is the same enforcement as `ntransform_gate.py` on a
different reference point. Its register is keyed on the SOURCE FILE and not on
an element, because a gate reads a file: three of its four divergences are
phosphorus and the fourth is the labile carbon the nitrogen transformation
operator denitrifies on. The soil nitrogen transformation operator is stock
LPJ-GUESS 4.1.1, which the tree can name by Zenodo record and SVN revision.
LPJ-GUESS 4.1.1 has no phosphorus at all: `pmass_labile`, `setptoc`, the
sorption isotherm and both saturation thresholds arrived with the CNP fork, so a
divergence on this path can be measured against nothing but the commit that
subtree was imported at, and `biosphere/config/somdynam.yaml` names it.

The declaration carries the four divergences, the two saturation constants, the
five C:P ramps and the four lines the phosphorus and labile-carbon arguments
rest on. The gate
fails on a constant whose line `modules/somdynam.cpp` no longer runs or whose
value is not what that line's own initialiser evaluates to; on a ramp whose
`setptoc` call the source no longer contains, or contains a different number of
times, the call being BUILT from the declaration so neither side can move alone;
on a ramp leaving the P:C range its two endpoints allow; on an invariant the
model no longer runs, which is how the phosphorus-limitation-off pin being split
off the ramp threshold is caught even though no constant moved; and on any of
the three halves of a divergence, so it can become neither a silent fork nor a
silent revert.

Each entry also says which configuration it is live in. Two of the four are
inert under the `ifplim 0` this project runs, and that is recorded as waiting
rather than as harmless: the comparison arm that would bound what they are worth
needs two runs, and the `ifplim 1` one needs `parameters.cpp`'s refusal lifted
first. That field names the phosphorus switch; the labile-carbon entry, whose
own switch is `ifntransform`, is live either way.

```bash
python biosphere/scripts/somdynam_gate.py            # status, exit 0
python biosphere/scripts/somdynam_gate.py --strict   # refuses while a saturation
                                                     # constant has no source
```

`--strict` has nothing to refuse on. `PMASS_SAT` is Parton, Stewart and Cole
(1988) Fig. 3's own axis maximum converted into the fork's labile-P currency,
and `PCONC_SAT` is gone with the ramp it drove: a decomposer community's biomass
C:P does not vary with its resource's phosphorus content, so the surface
microbial pool holds the C:P its initialiser sets. The refusal kind stays
because the next unsourced threshold needs somewhere to land, and it is the
refusal `parameters.cpp` makes on `ifplim 1` restated where a constant is
declared rather than living only in a C++ error string.

### The CENTURY equilibrium accelerator is the daily nutrient operator

`century_acceleration_gate.py` enforces SDEC-1 on the active C-N path and on the
waiting phosphorus path. The accelerator records products of daily mineral
uptake and leaching survival fractions, takes a geometric monthly mean over the
sample years, and applies that composed operator. It samples weathering and
deposition as daily fluxes rather than summing their year-to-date diagnostics.
Accelerated phosphorus uptake, leaching and additions all call the same
`pmass_add()` labile--sorbed isotherm as daily integration, and the hidden
40,000-year solve restores every mineral-P pool afterward.

```bash
python biosphere/scripts/century_acceleration_gate.py
```

The command launches no model. Its deterministic fixtures compare daily and
monthly survival, reject the old arithmetic-mean and inverted-leaching forms,
close a reduced C-N-P ledger, preserve the P isotherm after every mutation, and
mutate the declared production source in three independently detectable ways.

### The simulated snowpack conducts what the climate model's snowpack conducts

`snow_thermal_gate.py` enforces one relation across two models. The simulated
snowpack's thermal conductivity used to be computed here from Sturm et al.
(1997) and in the climate model from Fourteau et al. (2021) Eq. (18), and at the
density the climate model declares the two differ by close to a factor of two,
so one snowfall insulated this model's soil about twice as well as that one's.
The two could never have been reconciled by matching VALUES: this model's snow
density is prognostic across the compaction ramp `modules/soil.h` declares,
where the climate model's is a single namelist key, so equal conductivities at
one density would be a coincidence at one point of two curves that diverge
everywhere else.

**The relation is declared in `lib/snow.py` and nowhere else.** Neither compiled
model can import a Python module at runtime, so each carries the adopted row as
a literal and `lib/snow.py`'s `check_restatements()` holds both to the one
table. This gate runs that check and so does `scripts/smoke_test.py`, which is
what stops the register certifying a divergence into a relation the climate
column has since moved off. A generated header was the alternative and covers
only one of the two, because the climate model's Fortran is committed source
that is read and edited by hand.

The adopted row is the FAST kinetics arm, which is the bracket's UPPER endpoint
rather than a point inside it, so the residual is one-signed: the modelled snow
may conduct LESS than the model says and cannot conduct more. `--strict` does
not refuse on it, because Fourteau's own Sect. 4.1 says which limit snow is in
is unresolved and a question the literature holds open is reported rather than
gated. `notes/audits/cryosphere-material-properties.md` argues the choice and
`analysis/ice_properties.py` evaluates all three relations across the density
range the two components span.

`biosphere/config/snow_thermal.yaml` registers the divergence. Its reference
point is a RELEASE and not a subtree commit, which is what separates it from
`somdynam_gate.py`: `update_snow_properties` arrived byte-identical to
`guess_4.1/modules/soil.cpp`'s apart from the stripped licence header, and it is
on by default, since `data/ins/global.ins` sets `iftwolayersoil 0` and `Ksnow`
becomes the conductivity of every active snow layer in the multilayer scheme's
numerical solve. The gate fails on a restatement that has drifted from
`lib/snow.py` or normalises the ice volume fraction by something other than the
density the fit was made against, on any of the three halves of the divergence,
on a compaction ramp whose ends `modules/soil.h` does not declare, and on a
coefficient of the superseded relation live anywhere in the vendored model with
comments stripped, so a second snow conductivity growing outside that one
function fails here too. `execution_arms` carries the matched arms that have run
and every entry claims one by name or says `none` and why, on the same terms as
`ntransform.yaml` and `fire.yaml`; the block is empty because no LPJ-GUESS run
this project has made has passed acceptance, and the gate refuses an entry with
no claim, a claim on an arm the register does not carry, a `none` with no
reason, and an arm no entry claims.

```bash
python biosphere/scripts/snow_thermal_gate.py            # status, exit 0
python biosphere/scripts/snow_thermal_gate.py --strict   # refuses on a divergence
                                                         # whose verdict is `gate`
```

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
and `parameters.cpp` refuses `acclimated_respiration 1` without it.
`biosphere/config/respiration_acclimation.yaml` is the declaration: the bracket,
a source for each end, and the one-factor sensitivity registered over it. Its
fast end is how fast Gifford (2003) reports respiration adjusting; its slow end
is the memory time scale QUINCY declares for the acclimation of maintenance
respiration, whose Eq. S23 is the relation this fork implements. `run_value_days`
is which arm a run sits on; its sentinel is `undeclared`, and under the path
decision below that sentinel is the settled value rather than a pending one,
because nothing is on either arm.

The sensitivity carries no pass/fail bar, deliberately: the response is
arithmetic on a first-order lag and was known before the registration was
written, so a bar here would be a criterion chosen after the result it judges.
It is carried as model-form uncertainty instead.

THE BASELINE TAKES THE STANDARD RESPIRATION PATH, and that is declared as
`path.runs` with its argument and with what would reopen it, not inherited from
whichever vendored instruction file is imported. `respiration_acclimated()`
takes no `respcoeff` argument, so switching paths replaces a coefficient
`Pft::init_cton_limits` has normalised by the tissue C:N windows with Sprugel et
al. (1996)'s two fixed reference rates, and that moves the level of sapwood and
fine-root maintenance respiration, changes what it depends on from leaf
longevity to growth temperature, and gives up the invariance under a tissue C:N
window rescaling that kept `respiration()` unmoved by the repair in
`notes/plant-physiology-carbon-allocation-audit.md` finding 11. Acclimation is
real and is not what is refused; a switch carrying two unsourced changes with it
is. `acclimation_gate.py` enforces the state, executes the registered
sensitivity on every invocation, refuses a declaration that disagrees with
itself or with the run instruction, refuses an acclimated path with no declared
memory length, and refuses a run instruction on a path the declaration did not
decide on.

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

Summergreen leaf and root litter is released over the month that coldest day
falls in, so the canopy sheds into the season its phenology restarts in.
`Date::month_of` places the day on this world's own month lengths, which are
generated with the year length; a month is otherwise an Earth calendar artifact
and this year is not twelve of anything.

`phenology_gate.py` is the enforcement. It fails on an Earth ordinal date back in
any compiled source, on a reader no longer keyed on the derived pair -- the
litter release included -- on the chill-day count being able to leave the
`Pft::gdd0` table it indexes, and on a landmark field missing from
`Climate::serialize`. Eleven fixtures run on every invocation, four of them built
to be wrong in a named way.

```bash
python biosphere/scripts/phenology_gate.py            # status, exit 0
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

### The fire operators' departures from mainline are declared and enforced

`biosphere/config/fire.yaml` registers every place this project's fire source
departs from mainline LPJ-GUESS 4.1.1, and `fire_gate.py` is the enforcement.
Its scope is the OPERATORS and not one source file, which is the structural
difference from `ntransform_gate.py` and `somdynam_gate.py`: GLOBFIRM's burned
fraction is `modules/vegdynam.cpp`'s `fire()`, its mortality is the same file's
`mortality_lpj` and `mortality_guess`, and the BLAZE effects path is
`modules/blaze.cpp`, so each entry names its own `source_file` and the gate reads
the union of them.

It carries the one divergence shape the other register gates have no case for.
`mainline_form: deleted_line` is a line the release runs that this project's
source does not run at all -- the GLOBFIRM burn-probability floor, deleted under
`world-v5j1` and argued in `notes/fire-model-audit.md` section 8. A deletion has
no changed line, so the third check the other gates make has nothing to look at,
and `absent_from_stripped` does its work instead: the fragments a reformatted
revert would reintroduce, each of which has to stay out of the code. The
recorded mainline form is searched for both verbatim and with its own trailing
comment removed, because a recorded line carrying a `//` can never survive
comment stripping and searching for it verbatim alone is a check that cannot
fail.

Execution is a property of each entry and never of the register.
`execution_arms` carries the matched arms that have run, empty and never absent,
and every entry names the arm that covers it or says `none` and why. No
LPJ-GUESS run this project has made has passed acceptance, so the block is empty
and the entry says so for itself; the gate refuses an entry claiming no arm, one
claiming an arm the register does not carry, an unexecuted entry that says
nothing about why none covers it, and an arm no entry claims.

```bash
python biosphere/scripts/fire_gate.py            # status, exit 0
python biosphere/scripts/fire_gate.py --strict   # refuses on a divergence
                                                 # whose verdict is `gate`
```

## One-off tools

- `scripts/score_prediction.py` -- one-off: scores a productivity prediction against an LPJ-GUESS run, the machinery behind BIO-2's nitrogen bracket. Registered under `one_offs` in `config/pipeline.yaml`; it generates nothing the pipeline reads.
- `scripts/check_forcing_contract.py` -- one-off: the ecological forcing contract, enforced in both directions. `config/ecological_forcing_contract.yaml` is the versioned declaration and `notes/ecological-forcing-field-contract.md` is its argument. It checks the DECLARATION against `outmod.f90:ecogp`, which is the producer that has to satisfy it -- a declared code the model does not write, a code the model writes that no row declares, a row missing a column of the schema, a process requiring a field the contract does not declare, a cadence that is neither a declared class nor the `undeclared` sentinel with an owner -- and it checks a PRODUCT against the closure, sign and bracketing identities. Fourteen declaration fixtures and eight product fixtures, all but one of each wrong in a named way, so it can fail on itself. Registered under `one_offs`; it generates nothing.
- `scripts/validate_drift_statistic.py` -- one-off: validates `lib/lpj_output.py:drift_bound`, the equilibrium contract's upper bound on a field's end-to-end drift, against series whose answer is known by construction. Four arms: the rate at which a field drifting AT the tolerance is nonetheless accepted, which is the rate the contract declares and the only one barred; the rate at which a settled field is refused, which is not barred and is what sizes the retained record; a memory time of one, where the construction must reproduce the textbook two-sample t and the estimator must recover the AR(1) memory time it was handed; and the same two rates at the per-field scatter and memory time measured on a real run's own record. The row's rejected candidate, a significance test at a Sidak-corrected level, is measured on the same trials so the two forms are compared on one set of numbers. It runs no model and reads no artifact the pipeline consumes; its JSON is evidence for a declared rate. Registered under `one_offs`.
- `scripts/derive_cover_tolerance.py` -- one-off: derives the equilibrium contract's tolerance on the two `fpc.out` cover quantities through the consumer that closes a loop with them, `build_surface_albedo.py --mode modelled`, and says whether the declared number still agrees with the chain. It DECLARES nothing: every link -- the two canopy albedos and their brackets, the substrate albedo, the per-gridcell rootable fraction, the model's own land mask, the run's cover, the measured surface-to-planetary attenuation and the one flux-to-kelvin conversion -- is read from the artifact that owns it, and the bar is `OFFSET_TOLERANCE_K / RESOLVING_FACTOR` read from `exoplasim/scripts/assess_convergence.py`. It refuses when its model-land area fraction disagrees with the one the attenuation was measured over, and when the chain stops being linear in either the cover or the attenuation. `notes/equilibrium-trend-null.md` has the argument. It generates nothing.
- `scripts/derive_trend_null.py` -- one-off: the offline instrument behind `notes/equilibrium-trend-null.md`. `--timescales` re-takes this world's ecological memory and relaxation times from a completed run, `--sigma-sweep` measures the size and power of the per-cell half against `trend.slope_standard_errors`, and the default run measures the per-cell null and its power across fields, run lengths and patch counts. It generates nothing.
