# Fire-model audit and implementation path

This note records the source and literature audit behind the `FIRE` issues in
the `bd` tracker.  It is about the fictional Vesper simulation.  It does not select
parameters by fitting a desired result, and it does not claim that an Earth
fire parameterisation becomes portable merely by setting human population to
zero.

## Scope and sources read

The code audit covered the vendored CNP fork's GlobFIRM, SIMFIRE and BLAZE
paths, the Vesper input adapter, and ExoPlaSim's convective-precipitation and
CAPE diagnostics.  The primary sources below were fetched into `references/`,
read, and recorded in `references/INDEX.md`:

- Thonicke et al. (2010), SPITFIRE;
- Li, Zeng and Levis (2012), the intermediate-complexity fire model;
- Knorr, Kaminski, Arneth and Weber (2014), SIMFIRE;
- Finney et al. (2014), cloud-ice-flux lightning;
- Romps et al. (2014), CAPE-times-precipitation lightning, including its
  supplement;
- Mangeon et al. (2016), INFERNO; and
- Rabin et al. (2017), the FireMIP phase-1 protocol and model comparison.

The findings below distinguish what those papers demonstrate from choices that
must be declared for Vesper.

## Findings

### 1. Stock SIMFIRE-BLAZE is not a population-free natural-fire model

Setting population to zero removes an inapplicable Vesper driver, but does not
remove SIMFIRE's Earth calibration.  The vendored input archive contains an
Earth satellite burned-area climatology and HYDE population history.  SIMFIRE
uses the former to distribute annual burning through the year and uses an
Earth-calibrated statistical frequency relation.  Knorr et al. (2014) supports
an ignition-saturated *control* at global Earth scale; it does not establish
that its coefficients or seasonal climatology apply to an unobserved planet.
Rabin et al. (2017) likewise classifies LPJ-GUESS-SIMFIRE-BLAZE as an empirical
burned-area model rather than a fire-spread model.

The current driver also enables BLAZE only when the input module calls itself
`cru_ncep`; this is a generator-identity test, not a capability test.  Supplying
the archive or renaming the Vesper driver would therefore conceal, rather than
solve, the missing physical inputs.

### 2. The useful division is ignition/occurrence, spread/area and effects

FireMIP's cross-model diagrams expose three separable decisions: how fires
start, how starts become burned area, and how burned area affects vegetation
and material pools.  The sources support a staged Vesper model at that seam:

1. derive natural lightning and successful ignitions from Vesper weather;
2. turn ignitions plus fuel and moisture into bounded fire occurrence and area;
3. apply audited BLAZE combustion and mortality effects.

INFERNO is a suitable reduced-complexity pattern for the first two steps, but
not a set of portable constants.  It computes flammability from weather, soil
moisture and fuel, then multiplies occurrence by fixed Earth vegetation-class
fire sizes; its lightning configuration assumes every ground strike ignites a
fire, and the evaluated version diagnosed effects without feeding them back to
vegetation.  Li et al. (2012) provides the more explicit intermediate pattern,
`fire count * average area per fire`, with wind-dependent spread.  Its fuel,
humidity, moisture, PFT spread and one-day-duration constants are nevertheless
Earth-calibrated.  SPITFIRE adds mechanistic spread and effects but also adds
four fuel classes and many Earth PFT, lightning and combustion parameters; it
is not the economical first port.

The recommended first implementation is therefore an INFERNO-like occurrence
layer with a Li-like explicit count-to-area seam, feeding BLAZE effects.  Every
Earth-derived coefficient must be exposed as a prior or bracket, not silently
accepted as a Vesper fact.

### 3. ExoPlaSim has ingredients for lightning, but its existing CAPE output is not the required local proxy

ExoPlaSim already carries convective precipitation as code 143 (`prc`) and a
CAPE diagnostic as code 322.  Vesper presently retains 143 only in compact
time-average bins and does not request 322.  More importantly, the main model
step calls `rainstep` before `hurricanestep`, where code 322 is calculated.
That makes the available CAPE post-convection.

Romps et al. (2014) proposed a flash proxy proportional to CAPE times
precipitation and found that the product of *spatially averaged* CAPE and
precipitation reproduced much of the CONUS-mean flash-rate variability.  The
paper explicitly warns that this does not validate a pointwise colocated
product: convection consumes CAPE, so a local scheme needs upstream or
pre-convective CAPE.  Its supplement also used subdaily three-dimensional
states and says daily mean thermodynamic fields are insufficient for the CAPE
calculation.  Its CONUS normalization is not a universal or exoplanet
normalization.

Consequently, merely adding code 322 to `REGULAR_CODES` and multiplying it by
143 would implement a proxy the paper did not validate.  The climate-side task
is to add a lightweight, pre-convective lightning-potential diagnostic (or an
explicit upstream construction) synchronous with convective precipitation,
retain its subdaily accumulation semantics, and benchmark its overhead before
enabling it in production.  Enabling the whole hurricane-diagnostic suite just
to obtain CAPE would couple unrelated diagnostics and may add avoidable cost.

Finney et al. (2014) is an important countercheck.  In its Earth model,
upward cloud-ice flux at 440 hPa gave a better spatial lightning distribution
than cloud-top-height, convective-precipitation and mass-flux alternatives.  It
requires cloud ice, updraft mass flux and cloud fraction at that level, which
ExoPlaSim's present physics and output contract do not supply.  It also still
requires normalization to an observed global flash rate, and the paper shows
that such scaling is model- and resolution-dependent.  Cloud-ice flux is thus
not a drop-in substitute here; it bounds the structural uncertainty in the
simpler CAPE/precipitation route.

### 4. Lightning is not ignition

All-flash rate, cloud-to-ground flashes and successful fire starts are distinct
quantities.  SPITFIRE, for example, imports an Earth lightning climatology and
then assumes a 20% cloud-to-ground fraction and 4% ignition efficiency.  Those
numbers are Earth priors, not Vesper constants.  INFERNO's assumption that each
ground strike starts a fire did not map one-for-one to burned area because
wet-weather flammability suppressed strikes in wetter environments; its varying
lightning case could reduce burning, illustrating the wet-lightning problem:
the same convection that supplies lightning often supplies rain.

The driver must therefore preserve synchronous lightning potential and
precipitation, apply an explicit cloud-to-ground and successful-ignition
bracket, and expose the wet-lightning filter.  Human ignition and suppression
are absent for Vesper rather than represented by an invented population
history.  A ubiquitous-ignition case remains valuable as a model-form control,
not as the central natural-fire case.

### 5. Fire weather is available upstream but chronology is currently lost

The ExoPlaSim archive can provide precipitation, temperature, humidity, wind,
pressure and radiation at adequate cadence.  The LPJ handoff currently drops
several of those fields and smooths or bins the temporal record.  Fire depends
on event ordering: wind and dryness before a start, rainfall coincident with a
strike, and drying after rain cannot be reconstructed from independent monthly
means.  Vesper's 30-hour solar day and 181-day orbit also make an unexplained
24-hour or 365-day convention incorrect.

The fire contract must build on BIO-13 and BIO-23 but go further: retain one
coherent chronological sequence of all fire drivers, define whether a fire
step is a local solar day or an absolute-time interval, and conserve integrated
fluxes during conversion.

### 6. BLAZE effects contain Earth assumptions and are incomplete for the CNP fork

BLAZE supplies useful fuel combustion, fire-line-intensity and mortality
machinery, but it is not already portable.  Its fire-weather and mortality
logic contains fixed Earth empirical coefficients and 365-day assumptions. In
particular, SIMFIRE uses absolute latitude 50 degrees when labelling barren,
shrub and tundra states, while BLAZE uses 30- and 50-degree bands to select
tropical/temperate/boreal mortality and litter tuning. Those latitude branches
must instead consume vegetation, fuel, weather and thermal state under the
BIO-29 latitude-use contract. The remaining effects coefficients need a named
parameter registry and Vesper calendar conversion rather than a wholesale
activation.

More seriously, `blaze.cpp` transfers and reports carbon and nitrogen during
combustion but contains no phosphorus transfers.  The older `fire()` path in
`vegdynam.cpp` has phosphorus-fire calculations and pool reductions commented
out with an unresolved note about the atmospheric fraction.  Enabling either
path under phosphorus limitation would therefore update C and N while leaving
the corresponding P in place, breaking tissue stoichiometry and the C-N-P
budget.  Fire activation must fail closed until volatilized P versus ash/soil
retention is declared, implemented for live and litter pools, and covered by
mass-balance tests.

### 7. Evaluation must diagnose mechanisms, not merely match burned area

FireMIP warns that correct burned area can result from compensating vegetation
and fire biases.  It therefore evaluates vegetation and hydrology alongside
burned area and uses fixed-driver experiments to isolate causes.  Vesper has no
observed burned-area climatology, so the analogous acceptance suite must be
comparative and budget based: no fire, GlobFIRM, ubiquitous ignition, and a
lightning-driven central/bracket case under identical accepted forcing.  It
must retain ignitions, flammability, area per fire, burned fraction, return
interval, intensity, mortality and C-N-P destinations.  Model choice should be
reported as uncertainty in downstream productivity, albedo, soils and smoke,
not collapsed into a single tuned answer.

### 8. The GlobFIRM burn-probability floor is a reporting cap applied to the state

`modules/vegdynam.cpp`'s `fire()` floored the Thonicke et al. (2001) Eqn 9
burned fraction at 0.001 per simulation year, with `// c.f. LPJF` as its whole
justification.  GlobFIRM is the live fire model, so the line was not dormant:
`run_lpj_guess.py` writes `firemodel "GLOBFIRM"` and `vesper_pfts.ins` sets
`vegmode "cohort"`, which routes `fireprob` into `mortality_guess`.  The floor is
deleted, and the deletion is registered as a declared divergence from mainline
LPJ-GUESS 4.1.1 in `biosphere/config/fire.yaml`.

**The provenance is real and it terminates in a restatement.**  "LPJF" is the
FORTRAN LPJ that `vegdynam.cpp`'s own reference block cites as Sitch et al., and
its maintained descendant carries the identical constant: LPJmL's
`src/soil/fire_prob.c`, vendored at `vendor/lpjml`, returns 0.001 for any
`fire_frac` below 0.001.  Two independent implementations state the number and
neither argues for it.  Walking the chain is the repair the OPAQUE class
prescribes, the chain has been walked, and there is no derivation at the end of
it, so the constant is not opaque.  It is not tuned either: a fit residual does
not arrive at exactly 1e-3 in two codebases, and no comment in either claims a
fit.

**What the number is, mainline states in its own source.**
`modules/commonoutput.cpp` writes `firert_gridcell += 1000.0` under
`// Set a limit of 1000 years` for any patch whose `fireprob` falls below 0.001.
0.001 per year is a reporting cap on the fire return interval -- a round cutoff
on where a diagnostic stops resolving -- and the deleted line applied that cap to
the model state.

**The cap is one per Earth year, and that is the diagnosis.**  LPJF and LPJmL
both run `NDAYYEAR 365`.  This model's year is 183 absolute days, declared as
`VESPER_YEAR_LENGTH_DAYS` in `framework/vesper.h`, so carrying the numeral across
delivered a whole Earth year of the floor every 183 days: 1.996 times the hazard
in absolute time, the same factor the VESPDRV4-to-V5 driver magic exists to catch
on nitrogen deposition.  The neighbouring quantity of the same class is converted
for exactly this reason -- `distinterval`, the generic patch-destroying
disturbance return time read 140 lines below in the same file, reaches the model
as 199.594 simulation years against an Earth calibration of 100.  The floor was
never converted because `biosphere/notes/time-base-unit-contract.md` excludes
fire rates from its registry and hands them to this review, which had not covered
the line.  The constant is therefore IMPLICIT-EARTH, which is a diagnosis and not
an endorsement: establishing that the unit is Earth's is what proves the number
wrong here.

Rescaling by 0.501 is refused rather than overlooked.  It would repair the unit
of a number that has no derivation, and the quantity it would rescale is a
diagnostic cap rather than a mechanism.

**Nothing downstream needed a nonzero burn probability**, and the source says so
rather than an assumption about it.  `fire()` already returns `fireprob = 0.0`
from its MINFUEL branch before the floor was reached, and `vegetation_dynamics`
leaves it at 0.0 wherever `has_fires()` is false, so every consumer already
handles a zero.  There is one division by `fireprob` in the tree,
`commonoutput.cpp`'s `1.0/patch.fireprob`, and it is already guarded by the same
0.001 -- the floor is what kept that guard dead for every fuel-bearing patch.
Removing the floor makes it live and does not move the reported fire return
interval: a floored patch took the else branch and reported 1/0.001 = 1000 years,
and the guard reports 1000 years.  What changes is that a patch the moisture term
says cannot carry fire stops reporting burned area it did not burn.  Every other
consumer is linear in `fireprob` and none divides.  The curve has no singularity
the floor was covering: Eqn 9's denominator has no root on the reachable domain,
with a minimum of 0.185, and `fireprob` goes to zero linearly with the fire index
at slope 0.1137.

**What the floor was worth.**  It bound wherever the annual mean daily fire
probability fell below 0.00914, which is fewer than 1.67
fully-flammable-day-equivalents out of the model year's 183 days -- the regime
where the moisture term says the patch cannot carry fire.  Where it bound, the
effect on the simulated vegetation ran through `mortality_guess`'s stochastic
draw, and the scale to judge it against is the converted generic disturbance: a
0.001 per model-year hazard of a stand-replacing fire against `distinterval`'s
0.00501, so the floor added 20 per cent to the patch-destroying disturbance rate
on every fuel-bearing gridcell the fire model says cannot burn.  Over
`nyear_spinup 2318` that is 2.32 expected fires per patch and a 90.2 per cent
chance of at least one, which is not a rounding on those gridcells' cohort age
structure.

How many gridcells that is cannot be stated from the source.  It is a property of
this world's soil moisture and fuel, and `biosphere/config/fire.yaml`'s
`comparison_arm` names the paired run that would measure it.

No LPJ-GUESS or ExoPlaSim run was performed for this audit.
