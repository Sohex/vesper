# The biosphere time-base unit contract

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf, and its biosphere is LPJ-GUESS, a vegetation model written for
Earth and ported to run on this world's calendar. Everything below is about the
model's parameters and the simulated processes they drive.

LPJ-GUESS was calibrated on a planet where a year is 365.2569 days and a month
is about thirty. Vesper's orbit is roughly half that, and the port kept the
24-hour timestep so that per-day rate constants would still mean what they were
calibrated to mean. What the port did not do was say, for each quantity labelled
"per year", WHICH year. This document says it, once, for every such quantity, and
is the registry that `biosphere/scripts/build_vesper_pfts.py` executes and that
`vendor/lpj-guess` cites at the point of every conversion.

The rule that makes this necessary rather than pedantic: the same symbol was
being read in two different units in the same binary. `leaflong` was leaf
lifespan in simulation years where `growth.cpp` read it and leaf lifespan in
Earth years where `guess.h` fed it to a regression fitted in absolute months. A
unit contract that lives only at the declaration cannot catch that. Every entry
in the registry below names the file where the quantity is READ.

## The three units

**Absolute day.** Twenty-four hours. The model's timestep, and the atom of
absolute time for everything below. It is deliberately NOT Vesper's 30-hour
rotation: stepping in real Vesper days would put the rotation error into
respiration, decomposition, snowmelt and phenology alike, whereas stepping in
Earth days confines it to daylength, where it has a known sign.
`biosphere/notes/lpj-guess-porting-audit.md` argues that choice.

**Earth year.** 365.2569 absolute days. Not a unit of this world at all: it is
the unit every published calibration in the vendored source was measured in, and
it appears in the contract for exactly that reason. `VESPER_EARTH_YEAR_DAYS`.

**Orbit.** One turn of Vesper around its star, rounded to whole 24-hour steps,
which is what the model can integrate. Also called a *simulation year* in
LPJ-GUESS's own vocabulary, and the two are the same thing here. One orbit is one
seasonal cycle. `VESPER_YEAR_LENGTH_DAYS`.

The conversion between the last two is `VESPER_EARTH_YEARS_PER_ORBIT`, written
`f` below and derived from `config/planet.yaml` by
`biosphere/scripts/build_vesper_header.py`. It uses the MODELLED year, whole days,
because that is the span the model actually integrates over; the true orbital
period differs from it by the rounding and is not what any model rate acts across.
`lib/orbit.py:earth_years_per_model_year` is the one place it is computed.

## The four classes

Every quantity in the model that carries a "per year" or "in years" belongs to
exactly one of these. The class decides the conversion; nothing else does.

| class | what it is | conversion from the Earth calibration |
| --- | --- | --- |
| ABSOLUTE-RATE | a process whose speed is set by physics, chemistry or physiology and does not know the orbit: decomposition, deposition, weathering, a mortality hazard, tissue turnover | a flux `q` per Earth year becomes `q / VESPER_EARTH_YEAR_DAYS` per absolute day; a FRACTION `r` per Earth year becomes `1 - (1-r)^f` per orbit |
| ANNUAL-SUM | a threshold on something the model accumulates over one of its own years: degree-day sums, annual net production per unit leaf area | scales by `f`, linearly, because the accumulation is linear in time |
| SEASONAL-CYCLE | an event or a window that happens once per seasonal cycle by definition: leaf-out, leaf shedding, an establishment attempt, a count of years used to smooth interannual variation | unchanged. One orbit is one seasonal cycle on this world as on Earth |
| DIAGNOSTIC | a reported quantity, not a model behaviour | unchanged in the model, but its reporting interval must be named, and a comparison against an Earth number must convert |

Two traps the classes exist to keep apart.

**A fraction is not a rate.** Halving the period does not halve a fraction. The
conversion `1 - (1-r)^f` leaves 1.0 at 1.0, which is what a summergreen shedding
every leaf every year means; a linear multiplier would quietly turn "all of it"
into "half of it". Where `r` is small the two agree, which is why the error hides.

**A duration in years is the reciprocal of an annual sum.** `nyear_spinup 500`
reads like an absolute statement and is not, so it scales UP by `1/f` while
`gdd5min_est` scales DOWN by `f`. Confusing the two directions is worse than
doing neither, which is why `build_vesper_pfts.py` names the class of every
parameter it touches, in the generated file and in its provenance JSON.

## Where the boundary is

Three places, and only three, convert between units. Anything else that needs to
is a defect.

1. **`build_vesper_pfts.py`.** `vendor/lpj-guess/data/ins/global.ins` is in EARTH
   units, because it is Earth's calibration. The generated `vesper_pfts.ins` is in
   MODEL units: a year in it is an orbit, a fraction in it is per orbit, a
   degree-day sum in it is over one orbit. Everything the C++ reads from an
   instruction file is therefore already in model units.
2. **`vesper.h`**, generated by `build_vesper_header.py`. Carries
   `VESPER_EARTH_YEAR_DAYS`, `VESPER_EARTH_YEARS_PER_ORBIT` and
   `VESPER_EARTH_MONTHS_PER_ORBIT`, which is what the source uses wherever an
   Earth calibration is hardcoded in the C++ rather than declared in the
   instruction file.
3. **`build_lpj_driver.py`**, for the fluxes the driver file carries. Its
   nitrogen deposition is declared per EARTH year, because deposition is an
   atmospheric flux; `vesperinput.cpp` divides by the Earth year to reach the
   per-absolute-day rate it hands the model.

The driver file's magic is `VESPDRV5`. V5 carries the same bytes as V4 and
differs only in what the deposition field means, which is exactly the change a
magic has to catch: a V4 file read by a V5 binary looks perfectly valid and
delivers about twice the nitrogen. A later replacement of the forcing transport
under EFOR-1 through EFOR-8 is V6 or later.

## The registry

Every non-fire quantity in the vendored source that carries a year. "Read at" is
where the value is USED, not where it is declared, because a declaration-only
reading is what missed `leaflong`. Fire rates are deliberately absent: they belong
to the dedicated review in `biosphere/notes/fire-model-audit.md`.

### Plant traits, from the instruction file

| quantity | read at | class | unit in the generated file | state |
| --- | --- | --- | --- | --- |
| `gdd5min_est`, `gdd5min` | `vegdynam.cpp` establishment test against `climate.agdd5` | ANNUAL-SUM | degree-days over one orbit | scaled down, already before this contract |
| `gdd0_min`, `gdd0_max` | `vegdynam.cpp` establishment window on `climate.agdd0_20.mean()` | ANNUAL-SUM | degree-days over one orbit | in the scaled-down list; both are "no restriction" sentinels in `global.ins` and the sentinel guard leaves them alone |
| `greff_min` | `vegdynam.cpp` growth-suppression mortality, against `greff_5.mean()` | ANNUAL-SUM | kgC/m2 leaf per orbit | scaled down HERE. Left at the Earth value it condemns cohorts for having a short year |
| `phengdd5ramp` | `growth.cpp` leaf-on ramp, from the growing season's start | SEASONAL-CYCLE, within-season | degree-days, absolute | deliberately unscaled; scaling it would be a real error |
| `k_chilla`, `k_chillb`, `k_chillk` | `growth.cpp` budburst chilling requirement, against `climate.chilldays` | SEASONAL-CYCLE, within-season | absolute days and degree-days | unscaled |
| `tcmin_surv`, `tcmin_est`, `tcmax_est`, `twmin_est`, `twminusc` | `vegdynam.cpp` bioclimatic limits | not a time quantity | degrees C | unscaled; physiology carries over |
| `leaflong` | `growth.cpp` raingreen leaf replacement, `guess.h:initsla` and `init_cton_min` | duration, YEAR-COUNT | simulation years | scaled up HERE, and the two readers made to agree. See below |
| `longevity` | `vegdynam.cpp` background mortality, against `Individual::age` | duration, YEAR-COUNT | simulation years | scaled up HERE. `age` counts simulation years, so both sides are orbits once this is |
| `turnover_leaf`, `turnover_root`, `turnover_sap` | `growth.cpp:turnover_np`, once per simulation year | ABSOLUTE-RATE, fraction | fraction per orbit | converted HERE as `1-(1-r)^f` |
| `est_max` | `vegdynam.cpp` sapling establishment | SEASONAL-CYCLE | saplings/m2 per growing season | unscaled. One seed crop per growing season is the mechanism, and this world has two growing seasons per Earth year |
| `estinterval` | `vegdynam.cpp`, years between establishment events | SEASONAL-CYCLE | growing seasons | unscaled, already before this contract |
| `turnover_harv_prod` | `driver.cpp` harvested wood product pools | ABSOLUTE-RATE, fraction | per Earth year | NOT converted. Inert without land use; belongs with the cropland fail-closed guard, BIO-27 |

### Simulation controls, from the instruction file

| quantity | read at | class | state |
| --- | --- | --- | --- |
| `nyear_spinup` | `vesperinput.cpp` | duration, YEAR-COUNT | scaled up, already before this contract |
| `distinterval` | `vegdynam.cpp` disturbance probability `1/distinterval` per year | duration, YEAR-COUNT | scaled up, already before this contract |
| `freenyears` | `somdynam.cpp` and `canexch.cpp`, the N-limitation delay | duration, YEAR-COUNT | scaled up, already before this contract |

### Hardcoded constants in the C++

| quantity | read at | class | state |
| --- | --- | --- | --- |
| `KMORTGREFF` | `vegdynam.cpp`, fraction of a suppressed cohort killed per simulation year | ABSOLUTE-RATE, fraction | converted HERE against `VESPER_EARTH_YEARS_PER_ORBIT` |
| `KMORTBG_LNF`, `KMORTBG_Q` | `vegdynam.cpp` background mortality shape | dimensionless | unchanged; the equation's time enters through `longevity` and `age`, both orbits |
| `NYEARGREFF`, `NYEARAAET` | `guess.h`, lengths of the `greff_5` and `aaet_5` buffers | SEASONAL-CYCLE | unchanged. They smooth interannual variation, and one sample is one growing season on either world. What needed converting was the CONTENT of `aaet_5`, not the window |
| `NYEAR_SEASONAL` | `guess.cpp:Climate::accumulate_seasonal_cycle`, the window `dtemp_seasonal` averages the seasonal cycle of air temperature over | SEASONAL-CYCLE | unchanged, and for the same reason as the twenty-year windows below: one sample is one seasonal cycle. It asks what climate a gridcell's vegetation is adapted to, not what happened this year |
| `Climate::coldest_day`, `Climate::warmest_day` | `driver.cpp` degree-day and chilling resets, `growth.cpp` leaf-on reset, `cropsowing.cpp` and `cropphenology.cpp` | SEASONAL-CYCLE, within-season | not a converted quantity at all. They were Earth ordinal dates and are now derived per gridcell from the forcing by `guess.cpp:Climate::find_seasonal_landmarks`, which is why no year ratio enters them. `biosphere/notes/implicit-earth-assumptions.md` finding 1 |
| twenty-year windows: `mtemp_min_20`, `mtemp_max_20`, `agdd0_20`, `hmtemp_20` | `driver.cpp` | SEASONAL-CYCLE | unchanged, for the same reason, and the biosphere README's stellar-cycle argument depends on it |
| `TAU_LITTER`, `TAU_SOILFAST`, `TAU_SOILSLOW` | `somdynam.cpp:decayrates` | ABSOLUTE-RATE | Earth years. Divisor corrected HERE to the Earth year. This is the `ifcentury 0` path and is NOT the live one |
| CENTURY `K_MAX[]` | `somdynam.cpp:decayrates_century` | ABSOLUTE-RATE | already per absolute day, from Parton et al. (2010). Correct as it stands, and this IS the live path under `ifcentury 1` |
| `HALF_OMEGA`, `LAG_CONV` | `soil.cpp` analytic soil temperature | SEASONAL-CYCLE | already derived from `Date::MAX_YEAR_LENGTH`; the wave the soil column damps is the seasonal one |
| the 31-day regression buffer and its `30.0` offset | `soil.cpp` soil-temperature trend | ABSOLUTE-RATE | absolute days, and must stay absolute. The SAME buffer feeds `climate.mtemp` over `date.ndaymonth`, which is a month. One buffer, two units, both correct |
| `SECS_PER_DAY` | `guessmath.h` and the photosynthesis path | absolute day | 24 h, while Vesper's rotation is 30 h. Out of scope here and owned by PCAR-1 |
| `nfert`, `pfert` divided by `date.year_length()` | `management.cpp` | ABSOLUTE-RATE | NOT converted. Cropland, inert on the natural-vegetation baseline, BIO-27 |

### Nitrogen and phosphorus fluxes

The half of the contract BIO-24 exists for. Each of these is an input at a
boundary, so each has a declared unit and a place where it is divided down.

| quantity | read at | class | contract | state |
| --- | --- | --- | --- | --- |
| nitrogen deposition `ndep` | `vesperinput.cpp:getclimate` | ABSOLUTE-RATE | declared kgN/ha per EARTH year in the driver file, divided by `VESPER_EARTH_YEAR_DAYS` | corrected HERE. It was divided by the simulation year, delivering an Earth year's deposition every orbit |
| the NH4:NO3 split | `vesperinput.cpp` | not a time quantity | declared 50:50 in the driver file | unchanged, and unmeasured on this world. ANUT-5 owns replacing it |
| Cleveland fixation, AET term `nfix_a` | `somdynam.cpp:soilnadd` | ABSOLUTE-RATE | `aaet_5.mean()` is a sum over one orbit; the regression wants cm per Earth year, so it is divided by `f` before entering | corrected HERE |
| Cleveland fixation, intercept `nfix_b` | `somdynam.cpp:soilnadd` | ABSOLUTE-RATE | kgN/ha per Earth year, so it inherits the same per-Earth-year result and the same daily divisor | corrected HERE. It was added once per orbit |
| `anfix_calc` to daily | `somdynam.cpp:soilnadd` | ABSOLUTE-RATE | kgN/m2 per Earth year, divided by `VESPER_EARTH_YEAR_DAYS` | corrected HERE |
| phosphorus deposition `dpdep` | `driver.cpp` accumulates it; `vesperinput.cpp` never sets it | ABSOLUTE-RATE | kgP/m2 per absolute day | ZERO on this world, because the Vesper input module supplies no P deposition. Not a unit defect; a missing flux, and ANUT-1 through ANUT-10 own supplying it |
| texture-path `Soiltype::pwtr` | `somdynam.cpp:soilpadd` | ABSOLUTE-RATE | kgP/m2 per EARTH year, divided by `VESPER_EARTH_YEAR_DAYS` | corrected HERE. It was divided by the simulation year, delivering an Earth year of weathering every orbit. BIO-5 emits the field against this declaration |
| gridded-path `pwtr_bi`, `pwtr_pcont`, `pwtr_shield`, `pwtr_ea` | `somdynam.cpp:soilpadd` | ABSOLUTE-RATE | already a daily calculation driven by `patch.soil.runoff` | correct as it stands, and inactive: it needs `file_pwtr`, which this world does not supply |
| `USORB`, `USSORB` | `somdynam.cpp:somfluxes`, and `equilsom` through it | ABSOLUTE-RATE | published per Earth year, divided by `VESPER_EARTH_YEAR_DAYS` | corrected HERE. `UOCC` is gone with them: it was declared and never used, so occlusion is a declared absence and not a slow process. See `biosphere/notes/phosphorus-cycle-parameterisation.md` |
| `apwtr`, `apdep`, `anfix`, `aNH4dep`, `aNO3dep` and the `.out` columns fed from them | `commonoutput.cpp`, `miscoutput.cpp` | DIAGNOSTIC | sums over one simulation year, which is one orbit | the reporting interval and the conversion factor are named in the run manifest, under `reporting_interval` and `annual_flux_per_orbit_to_per_earth_year` |

## `leaflong`, and why it is declared in orbits

`leaflong` had two readers in two units, which is the defect this whole document
is answering. Settling it needed a modelling decision and not a conversion,
because Reich et al. (1992) regressed specific leaf area and leaf C:N on leaf
lifespan in absolute MONTHS, while `growth.cpp` divides a count of leaf-display
DAYS by it.

Both readings want an absolute lifespan; they disagreed only about which unit the
declaration was in. The decision is that the instruction file is in model units
throughout, so `leaflong` is in simulation years like every other duration in it,
and the Reich regression is fed `VESPER_EARTH_MONTHS_PER_ORBIT * leaflong` rather
than Earth's twelve.

This is exactly value-preserving for the regression. The generated file carries
`leaflong_orbits = leaflong_earth / f`, the regression multiplies by `12 f`, and
the `f` cancels: SLA, `cton_leaf_min` and the `ctop_leaf_min` derived from it are
unchanged to six figures from the shipped Earth calibration. The alternative
reading, keeping the declaration in Earth years and converting in `growth.cpp`
instead, is equally consistent but would have made the instruction file the only
place in the port still speaking Earth.

The one reader whose ANSWER changes is `growth.cpp`'s raingreen excess-leaf
allocation, and it changes by the factor the defect was worth. It compares
`aphen_raingreen`, a count of leaf-display days in the simulation year, against
`leaflong * date.year_length()`. Before, a raingreen PFT with `leaflong 0.5` had
an implied leaf lifespan of half an ORBIT; now it has half an EARTH YEAR, which
is about the whole orbit. On a cell where such a PFT holds leaves for most of the
year, the excess-leaf term therefore falls from a substantial fraction of leaf
plus root biomass shed to litter each orbit to approximately none. That is a
correction, not a tuning: a leaf that lives half an Earth year does not need
replacing twice in a year that IS half an Earth year.

## What the corrections are worth

Stated because a correct term that worsens an agreement is still information, and
because none of these has been executed: the model does not build on this tree
until a baseline climatology exists to generate `vesper.h` against, so every
figure here is arithmetic and none is a measurement.

- **Soil carbon roughly doubles at equilibrium, on the `ifcentury 0` path only.**
  Equilibrium pool size is annual litter input over the annual decay constant.
  Input per orbit is what it is; the decay constant per orbit halves; the ratio
  doubles. The live `ifcentury 1` path was already correct and does not move.
- **Nitrogen deposition roughly halves.** It was delivering a declared
  per-Earth-year amount once per orbit. Since the declared value is an assumption
  and not a measurement, the honest statement is that the model now receives what
  the driver file says it receives.
- **Cleveland fixation falls by more than half.** The AET-proportional term
  roughly halves because it was reading an orbit's evapotranspiration as an Earth
  year's, and the intercept halves because it was applied once per orbit. Both
  move the same way, so the registered `nfix_a`/`nfix_b` bracket in BIO-2 has to
  be re-derived against this contract rather than carried over.
- **Phosphorus weathering input roughly halves per orbit** and is unchanged per
  unit absolute time, which is what a rock-weathering rate should be. Labile P
  is a small, fast pool against that input, so the equilibrium labile stock is
  set by the balance of supply against uptake and leaching rather than by the
  supply alone; the effect on P limitation therefore has to be measured and is
  not read off this arithmetic.
- **Sorption to the strongly sorbed pool slows by the same factor**, and its
  equilibrium does not move at all: `USORB` equals `USSORB`, so the strongly
  sorbed pool settles at the size of the sorbed pool whatever the rate is. Only
  the approach time changes, from about 150 orbits to about 150 Earth years.
- **Growth-suppression mortality falls from 0.3 to about 0.16 per orbit**, which
  is the same hazard per unit absolute time.
- **`greff_min` halves**, which removes a spurious condemnation: growth efficiency
  is annual production per unit leaf area, so it was being measured over half a
  year and compared against a whole Earth year's threshold.
- **Tissue turnover falls to roughly the rate-converted fraction**, which is the
  same loss per unit absolute time and therefore roughly halves litter production
  per orbit while leaving it unchanged per Earth year.
- **SLA, leaf C:N and leaf C:P do not move at all.** The `leaflong` rescale and
  the months-per-orbit conversion cancel exactly. This is the behaviour the
  correction preserves, and it is worth saying so: the bug was a unit
  inconsistency, and closing it was not an excuse to change the Earth PFTs.
- **The scale factor itself moves by 0.11%**, because it is now derived from the
  modelled year of whole days rather than the true orbital period. That is below
  every threshold that reads it and is corrected only so that the factor is one
  fact rather than two.

## What is not settled here

- The phosphorus route's remaining unknowns are not units, and none of them
  needs one: `PMASS_SAT` is a mass per unit area, `PCONC_SAT` a mass fraction,
  and the `PFRAC_*` ratios are dimensionless, so no orbit enters any of them.
  `parameters.cpp` refuses `ifplim 1` while `PFRAC_LEAFTOSAP` has no derivable
  scalar, `PCONC_SAT` has no source, and `PMASS_SAT` reads a labile P pool its
  source did not define. What each one is worth, and what would settle it, is in
  `biosphere/notes/phosphorus-cycle-parameterisation.md`.
- Nothing here is verified by execution. LPJ-GUESS does not build on this tree:
  `framework/vesper.h` is generated, and the chain to it runs through a baseline
  run and a baseline climatology that do not exist. The regression checks BIO-22
  asks for are therefore no-simulation checks against the generated instruction
  file and its provenance, and belong with BIO-20.
- The Earth PFTs themselves are unchanged and remain a declared choice. Rescaling
  their limits makes them mean on this world what they meant on Earth; it does not
  make them a prediction of this world's physiology.
