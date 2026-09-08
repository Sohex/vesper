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

Four more were read on 2026-09-07, to walk the BLAZE effects layer's own
citations to their ends rather than to add a model form:

- Noble, Gill and Bary (1980), the McArthur meters as equations, which is the
  source of every fire-danger term `blaze.cpp` runs and which settles that its
  rate-of-spread coefficient is a unit conversion rather than an empirical
  value;
- Sirakoff (1985), the two-page correction to Noble that `blaze.cpp` already
  implements and does not cite;
- Liedloff and Cook (2007), for the fuel heat yield, which that paper states as
  an assumption rather than a measurement; and
- Arora and Boer (2005), for the calm-air fire growth ratio Li et al. (2012)
  reports second-hand, which that paper also states as an assumption.

The findings below distinguish what those papers demonstrate from choices that
must be declared for Vesper.  Findings 9 through 11 were added by that second
reading; the constants themselves, with their dispositions, are in
`biosphere/config/fire_parameters.yaml` rather than restated here.

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
prescribes, and the chain is walked to its end: past the code, into the primary
description of LPJF's fire module.  Thonicke, Venevsky, Sitch and Cramer (2001)
states no minimum burned fraction anywhere.  It states the opposite -- its Eqn 5,
`A(s) = s * f(s)`, is introduced with "A is zero when fire conditions were absent
during the year" -- and its own reporting stops at "more than 900 years" for
regions unsuitable to carry fire, which is the same diagnostic ceiling and not a
floor on the state.  Sitch et al. (2003), the LPJ-DGVM description, documents the
fire terms as empirical global-Earth choices and states no floor either.  Two
implementations and two papers, and no derivation at the end of it, so the
constant is not opaque.  It is not tuned either: a fit residual does not arrive
at exactly 1e-3 in two codebases, and no comment in either claims a fit.

What Thonicke DOES fit is Eqn 8's four coefficients -- 0.45, 2.83, 2.96 and 1.04
-- by non-linear least-mean-square regression against observed area burnt in
Portugal, southern California and Kakadu National Park.  That is a published fit
with inspectable provenance and a stated Earth domain, which is a different thing
from a tuning, and it is the curve this project keeps.  It is a separate question
from the floor that sat on top of it.

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

### 9. SIMFIRE is not an alternative to BLAZE; it is BLAZE's mandatory preprocessing

Finding 1 said the driver "enables BLAZE only when the input module calls itself
`cru_ncep`", and described that as a generator-identity test. The mechanism is
both simpler and harder than that, and it changes what FIRE-4 has to do.

`modules/driver.cpp:734` reads:

```
if (firemodel == BLAZE) {
    simfire_accounting_gridcell(gridcell);
    blaze_accounting_gridcell(gridcell.climate);
}
```

So selecting BLAZE runs SIMFIRE's whole annual accounting unconditionally.
`simfire_accounting_gridcell` calls `getsimfiredata` on the first day of the
simulation, and that function opens `file_simfire` and calls `fail()` if it
cannot. The archive it opens holds Hyde 3.1 population density and a monthly
burned-area climatology. `simfire_update_pop_density` is called beside the biome
mapping every year. There is no configuration in which BLAZE runs and this does
not.

The coupling is load-bearing in the other direction too. `blaze.cpp:1012` selects
the litter combustion factor `k_tun_litter` on `gridcell.simfire_biome`, and
`blaze.cpp:403` selects which tree survival parameterisation applies on the same
field. `simfire_biome` is written only by `simfire_biome_mapping`, which is
reached only from inside `simfire_accounting_gridcell`. Removing the SIMFIRE path
without replacing that field leaves both selections reading a member that
`Gridcell`'s constructor does not initialise.

**The archive would not refuse this world; it would answer for it.** The lookup
in `getsimfiredata` keys on `floor(lon*2)/2 + 0.25` and the matching latitude, so
a Vesper gridcell carries coordinates that are perfectly valid keys into an Earth
0.5-degree archive. It would not fail to find a record. It would return the human
population history and the observed fire seasonality of the Earth location at the
same latitude and longitude, and the run would complete. This is the concrete
form of finding 1's warning that supplying the archive conceals rather than
solves the missing inputs, and it is why FIRE-4 must fail closed on a capability
rather than on a file being absent.

Consequence for the staged plan in finding 2: the ignition/occurrence seam and
the effects layer are separable as MODEL FORM but not as work. FIRE-6 cannot arm
BLAZE while FIRE-4 has removed the field BLAZE selects its coefficients with, and
FIRE-4 cannot remove SIMFIRE while BLAZE is the only thing that calls it. The
biome classification `simfire_biome` carries is the piece that has to survive the
removal, and it is derived from the model's own simulated vegetation and land
cover rather than from the archive -- `update_fire_biome` reads patch state, with
absolute latitude 50 degrees separating shrubland from tundra and barren, which is
the BIO-29 latitude branch finding 6 already names.

### 10. BLAZE's rate of spread is Noble's coefficient with a lost decimal

`blaze.cpp:227` declares `A = 3.3333e-05` under the comment `// Empirical
value`, and uses it at `:247` as `rate_of_spread = A * mcarthur_forest_fire_index
* avail_fuel`.  Noble et al. (1980) gives the Mark 5 forest rate of spread as
`R = 0.0012*F*W`, and its Table 1 fixes the units: R in km/h, W in tonnes/ha, F
dimensionless.  `blaze.cpp` computes F from Noble's own drought-factor and index
equations verbatim at `:1076` and `:1081`, so nothing rescales between them.

In the units `blaze.cpp` uses -- fuel in g/m2, spread in m/s -- Noble's
coefficient is `0.0012 / 100 / 3.6 = 3.333333e-06`.  The source compiles ten
times that, and the mantissa agreeing to five digits is what makes it a lost
decimal rather than a different calibration.  The constant is not empirical; it
is a folded unit conversion whose comment failed to say so.

Three independent readings confirm the units, so no alternative rescues it:
`MIN_FUEL = 200` is commented g/m2 and matches GlobFIRM's 200 gC/m2 threshold;
`HEAT_YIELD * avail_fuel * rate_of_spread` evaluates to kW/m only with R in m/s;
and the intensity class bounds 750, 3000 and 7000 kW/m are standard Byram
values.  Reading R as m/min, m/h or cm/s gives 2e-04, 1.2e-02 and 3.3333e-04,
none of them the declared number.

It is worth a class boundary rather than a rounding.  At F = 10 with 1000 gC/m2
of fuel, Noble gives 667 kW/m, below the lowest bound; `blaze.cpp` gives 6667,
the second-highest class.  At F = 20 with 2000 gC/m2, 5333 against 53 333, which
is 7.6 times above the top bound, so every such fire takes the maximum
combustion column and the maximum mortality.  `world-c15e`.

A separate error multiplies with it.  `available_fuel` sums `cmass`, which is
carbon, while Noble's W and Byram's w are fuel weight; no carbon-to-dry-matter
factor exists anywhere in the fire path.  `world-9zfv` carries that, filed apart
because the two have different right answers and one combined factor that made
the product look plausible would conceal both.

### 11. The combustion table's column index is wrong at both ends

`TURNOVERFRACT` is `double[13][5]`, selected by a fire-line-intensity class.

At the low end this is a defect.  `get_fireline_intensity` sets the intensity to
`-1.` when available fuel is below `MIN_FUEL` -- the documented "not enough fuel
to ignite" case -- and `get_fire_line_intensity_index` then falls through every
branch and returns `-1`.  `blaze()` passes that index straight to
`get_combustion_rates` with no guard; its only two early returns are on
flammable area and the stochastic burn draw, and neither consults fuel.  The
array is contiguous, so `[r][-1]` aliases `[r-1][4]`, the highest-intensity
column of the row above, and `[0][-1]` reads before the array.  A patch with too
little fuel to carry a fire therefore combusts all of its surface litter and half
its leaves, and is flagged burned.  The sign is backwards and the row-0 read is
undefined behaviour.  `world-voxz`.

At the high end it is dead code.  The table's header comment says columns 3 and 4
are the sprouter and seeder rates above 7000 kW/m, and the index function returns
4 for every intensity above 7000 and never 3.  `is_resprouter` exists in the file
but reaches only a survival probability, never a column.  Every intense fire
therefore takes the seeder rates, which for stems and branches to litter are 0.8
against the sprouter 0.2.

### 12. The drought index reads a 183-day rainfall total into an Earth-annual coefficient

`blaze.cpp:1069` is the Keetch-Byram drought index, and its rainfall term is
`exp(-.0441 * climate.rainfall_annual_avg / 25.4)`.  That running mean is
accumulated and reset on `date.islastday && date.islastmonth`, so it is a total
over this model's year of 183 absolute days, while the coefficient was fitted to
Earth annual rainfall over 365.

Because the term is exponential, the resulting bias depends on the climate
rather than scaling the field: holding the daily rate fixed, the drought
accumulates at 0.69 of the intended rate at an Earth-equivalent 500 mm/yr, 0.53
at 1000 and 0.46 at 1500.  A bias that varies with the cell compresses the range
of simulated fire weather and moves its spatial pattern.  It also runs the
opposite way from the guess, because rainfall enters through a negative
exponential: a shorter year makes the world read as wetter and burn less.

A SECOND INSTANCE IN THE SAME FILE, and the only literal Earth year left in the
fire path.  `blaze.cpp:1061` computes the year-end realignment of the 30-day
fire-danger ring buffer as

```
int avg_shift = AVERAGING_FFDI - (365 % AVERAGING_FFDI);
```

which is 25.  The buffer is written at `date.day % 30` and read as a rolling
maximum, so at a year boundary it has to be rotated by the year's remainder
modulo 30 or the first month of the new year mixes days from two years in the
wrong order.  For a 183-day year the shift is `30 - (183 % 30) = 27`, so the
rotation is off by two slots on every year boundary.  The repair is not 27
either: it is `date.year_length() % AVERAGING_FFDI`, so the buffer tracks the
declared year rather than a second hardcoded one.

The 30-day window is a separate question and is left open deliberately.  This
model's months are 15 days, eleven of them, and one of 18, so a "monthly" fire
danger maximum taken over 30 days spans two model months.  Whether the window is
an absolute-time memory of drying, in which case 30 days carries, or a fraction
of the seasonal cycle, in which case it does not, is a choice FIRE-6 has to make
rather than one the source settles.

**The fire path's Earth-calendar exposure is now enumerated and this is all of
it.**  Stripping comments and searching `blaze.cpp`, `simfire.cpp` and
`vegdynam.cpp` for year and month literals returns exactly: the `365` above; the
30-day window at `:1026`, `:1029` and `:1059`; three `for (i < 12)` loops over
the model's own twelve months, which are correct; and the absolute-latitude
bands of finding 6.  `vegdynam.cpp`'s `fire()` carries no calendar literal at
all -- its defect is finding 13's, which is a unit on an output rather than a
number in the source.

This is the same class as the GlobFIRM burn-probability floor of finding 8,
found the same way, and `biosphere/notes/time-base-unit-contract.md` excludes
fire rates from its registry and hands them here.  `RAINFALL_AVERAGING_SPAN = 3`
is the second half and is not obviously the same repair: a memory length may be a
count of seasonal cycles rather than a duration, and it should be decided rather
than inherited.  The rest of the block is imperial throughout and its ceiling of
800 is 203.2 mm of soil moisture deficit standing for an Earth soil profile's
water capacity.  `world-4iyz`.

### 13. GlobFIRM's burned fraction is a per-Earth-year rate on a 183-day year

Finding 8 diagnosed the burn-probability FLOOR as one-per-Earth-year and it was
deleted.  `biosphere/config/fire.yaml` records that the curve the floor sat on
top of "is a separate question".  This is that question, and it has the same
answer.

`modules/vegdynam.cpp:1431` is Thonicke et al. (2001) Eqn 9:

```
s = n / (double)date.year_length();
sm = s - 1;
fireprob = s*exp(sm/(0.45*sm*sm*sm + 2.83*sm*sm + 2.96*sm + 1.04));
```

**The argument going in is correct, and that is what hides it.**  `n` is the
summed daily fire probability over the year, so `s` is the annual MEAN DAILY
probability -- a dimensionless fraction of days.  Holding the daily fire weather
fixed, `n` scales with the year's length and `s` does not, so `s` computed over
183 days is the same physical quantity Thonicke computed over 365.  Nothing is
wrong on the way in.

**The output is a burned fraction per EARTH year.**  Eqn 9's four coefficients
were fitted against observed annual area burnt, on a world whose year is 365
days.  The model applies that fraction to a year of 183 absolute days --
`VESPER_YEAR_LENGTH_DAYS` in `framework/vesper.h`, returned by
`date.year_length()` through `Date::MAX_YEAR_LENGTH` at `guess.h:529`.  The same
fire climate therefore delivers a whole Earth year of burning every 183 days:
365.2569/183 = 1.9959 times the burned area per unit absolute time.

It is consumed as a per-model-year rate, which closes the argument.
`mortality_guess:1050` draws `randfrac(...) < fireprob` once per simulation year,
and `mortality_lpj:873` and `mortality_guess:907` set
`mort_fire = fireprob*(1-fireresist)` per simulation year.  No downstream
division by a year length absorbs it.

The same file already applies the conversion elsewhere.  `vegdynam.cpp:1005`
computes the growth-efficiency mortality as
`1.0 - pow(1.0 - KMORTGREFF_PER_EARTH_YEAR, VESPER_EARTH_YEARS_PER_ORBIT)`, which
is the rule `vesper.h` states for an Earth-calibrated fraction per Earth year,
and `distinterval` reaches the model as 199.594 simulation years against an Earth
calibration of 100.  The comment at `:1459` names `distinterval` as "the same
class of quantity read 140 lines below", converted "for exactly this reason".
The curve's own output sits two lines above that comment.

The conversion roughly halves the probability -- ratio 0.5010 at 1e-4, 0.5023 at
1e-2, 0.5142 at 0.1 -- and fire is 11.7 per cent of this world's simulated carbon
turnover, so it is a correction to a material term rather than to a diagnostic.
The equilibrium response is not simply half: a lower burn probability lets fuel
accumulate and makes each fire larger, so the sign is unambiguous and the
magnitude needs the run.  Thonicke's four coefficients are not what is wrong and
do not move.  `world-drim`.

No LPJ-GUESS or ExoPlaSim run was performed for this audit.
