# Implicit Earth assumptions below the LPJ-GUESS port

**Scope:** the Vesper adapter and the vendored LPJ-GUESS-CNP source. This was a
read-only source audit; no model build or simulation was run. BLAZE and GLOBFIRM
are deliberately excluded here so that fire can be reviewed as one mechanism
rather than scattered across the general port, except for the cross-cutting
latitude-use inventory in finding 8.

The existing porting note records the deliberate use of Earth PFTs, the
181-day/24-hour-step calendar, orbital geometry, PAR energy fraction, degree-day
rescaling, spin-up duration and the known frost limitation. Those are not
reopened. The findings below are assumptions that remained implicit underneath
those declared choices.

## 1. Natural phenology took its seasonal landmarks from Earth ordinal dates

`framework/guess.h` defined the northern coldest day as 15 January and the
southern as 15 July, aliasing each hemisphere's warmest day to the other's
coldest, and the southern date named a day the compiled calendar does not have.
Three active natural-vegetation paths read them: `modules/driver.cpp` reset the
summergreen GDD counter and switched chilling detection off and on,
`modules/growth.cpp` reset annual leaf-on accumulation, and
`modules/somdynam.cpp` released summergreen leaf litter in the first month of the
year in the north and the first month of the second half in the south.

So the southern GDD and leaf-on resets never occurred at all, and in the north
chilling detection was switched off at the northern coldest day and never
switched back on. `chilldays` was also incremented under an inclusive
upper-bound test, so losing its seasonal reset advanced it one index past the
end of `Pft::gdd0` on a cell whose monthly mean never crosses the 5 degree base
from above -- a defect the vendored source has on Earth too, which the missing
reset made reachable here.

`Climate::coldest_day` and `Climate::warmest_day` carry the landmarks now,
derived per gridcell from the temperature forcing itself: a running day-of-year
mean of the air temperature, smoothed over one of this world's months and
searched circularly for its extremum. No hemisphere test enters, and no thermal
lag, which is what scaling the two Earth dates by the year ratio would have
carried over. The two can never be equal, so each reset fires exactly once per
orbit on every gridcell, the equator included, and on a forcing of any number of
orbits: `VesperInput` hands the whole interpolated year over at day 0, so a cell
runs on its own landmarks from the first orbit and follows a multi-year driver
file as it cycles. The chill-day guard is now exclusive, which holds the count at
the last valid index; the budburst requirement it looks up is at its `k_chilla`
asymptote well before that, so the ceiling changes no result.

`biosphere/scripts/phenology_gate.py` is the enforcement, and it carries the
fixtures: northern, southern, equatorial, flat and sixteen-orbit forcings, a cell
whose air temperature never reaches the chilling base, and four cases built to be
wrong in a named way, including the two Earth ordinal dates on this calendar and
the chill-day guard the vendored source shipped.

The summergreen leaf litter release in `modules/somdynam.cpp` is on the same
pair. It sheds over the month containing `Climate::coldest_day`, placed by
`Date::month_of` on this world's own month lengths, rather than over the
January or July a hemisphere test used to select. No natural-vegetation event
now reads an Earth ordinal date, and the gate's wiring list is what holds that.

## 2. A model year is still treated as an Earth year in ecological rates

`build_vesper_pfts.py` currently rescales the two degree-day sums and three
duration controls. Many other quantities labelled per year still execute once
per 181-day orbit:

| mechanism | source examples | present interpretation |
| --- | --- | --- |
| leaf, root and sapwood turnover | PFT `turnover_*` | fraction per orbit |
| leaf and tree longevity | `leaflong`, `longevity`, individual age | model orbits |
| establishment and growth-efficiency mortality | `est_max`, `greff_min`, `NYEARGREFF` | per orbit / five orbits |
| litter and SOM turnover | `TAU_LITTER`, `TAU_SOILFAST`, `TAU_SOILSLOW` in `somdynam.cpp` | same decay fraction per orbit, on the `ifcentury 0` path only. The live CENTURY path takes its `K_MAX` from Parton et al. (2010) already on a daily basis, which is absolute time and was never affected |
| CNP P-pool kinetics | `USORB`, `USSORB` in `somdynam.cpp` | published Earth-annual rates, divided by `VESPER_EARTH_YEAR_DAYS`. `UOCC` is gone with the occlusion decision |

For an absolute-time calibration, retaining the Earth annual fraction on every
Vesper orbit makes the process about 2.02 times as fast per Earth year. But not
every annual event should be rescaled: deciduous shedding and establishment may
legitimately occur once per seasonal cycle. Fractions also require a rate
conversion such as `1 - (1-r)^f`, not a linear multiplier, where `f` is the
orbit/Earth-year ratio.

The fix is therefore an explicit registry of each annual quantity as
absolute-time, seasonal-cycle, accumulated-flux or diagnostic, followed by
parameter-specific conversions and no-simulation checks. Fire rates are outside
this audit and registry pending the dedicated fire review. This is BIO-22, and
the registry it asks for is `biosphere/notes/time-base-unit-contract.md`, which
carries the class, the unit and the READER of every quantity here and records
which conversions have been applied.

## 3. The climate handoff discards planetary pressure and longwave radiation

The ExoPlaSim climatology already carries `ps`, `rls`, `hus` and `spd`, and the
pedology component consumes them in its open-water Penman calculation. The LPJ
driver carries only `tas`, precipitation, `rss` and diurnal temperature range.

LPJ consequently reconstructs net upward longwave from air temperature and
inferred sunshine using the Earth empirical constants `A=107`, `B=0.2`,
`C=0.25`, and `D=0.5`. Its psychrometric constant is fixed near sea level.
Photosynthesis separately uses `PATMOS=100000 Pa` and `PO2=20900 Pa`; `Climate`
has no pressure field at all.

Vesper's declared sea-level composition is intentionally Earth-like, so these
values are close at sea level. The defect is spatial: mountain cells retain
sea-level CO2 and O2 partial pressures and sea-level evaporative physics despite
the GCM already resolving their surface pressure under stronger gravity.

BIO-23 carries surface pressure, net longwave, humidity and wind across the
driver, uses actual pressure in the photosynthesis and psychrometric terms, and
compares the retained LPJ EET formulation with the project's existing Penman
calculation before changing model form. The latter two fields remain part of the
artifact even if that comparison decides not to use them in baseline LPJ.

## 4. Nutrient fluxes do not have one definition of “per year”

The LPJ driver labels nitrogen deposition `kgN/ha/yr`, but `vesperinput.cpp`
divides that amount across the 181-day model year. A declared 0.5 kgN/ha per
Earth year therefore becomes 0.5 per orbit, or about 1.01 per Earth year.

The Cleveland fixation relationship also needs a narrower correction than the
existing note claims. Its AET-proportional term mostly normalises correctly when
reported per Earth year, but its non-zero intercept is applied once per orbit
and its five-year AET mean spans five orbits, about 2.47 Earth years. This changes
both the central estimate and the registered `nfix_a`/`nfix_b` bracket.

The texture-path phosphorus weathering input has the same ambiguity:
`Soiltype::pwtr` is documented as kgP/m2/year and divided across every model
orbit. BIO-5 currently uses that wording for a pedology-derived flux, while the
runoff-driven CNP route is already a daily calculation. BIO-8 deliberately uses
a daily P-deposition contract and does not have this defect.

BIO-24 defines absolute-day, Earth-year and orbit units at every nutrient input
boundary, corrects deposition and fixation, and makes BIO-5's `pwtr` contract
unambiguous before the field is emitted. The mass-balance diagnostics must name
their reporting interval as well. The definitions and the per-quantity state are
in `biosphere/notes/time-base-unit-contract.md`.

## 5. The PAR correction changes energy but not photons per joule

`build_vesper_header.py` derives `VESPER_FRADPAR` from the energy inside the
selected photosynthetic window. `modules/canexch.h` then converts that energy to
quanta using `CQ=4.6e-6 mol/J`, explicitly a 550 nm value. This is used directly
in the electron-transport-limited photosynthesis calculation.

The registered productivity prediction instead integrates spectrum-weighted
photon supply for the K star and the selected 400-750 nm window. Thus the
prediction and executable model currently use different light currencies.
BIO-25 generates a spectrum/window-weighted `VESPER_CQ` beside `FRADPAR` and
checks their combined photon supply against the registered calculation.

## 6. An active canopy scalar is tuned to Earth global totals

`ALPHAA_NLIM=0.6` is documented in `modules/canexch.h` as chosen to make global
carbon pools and fluxes agree with published estimates. This tuning is broader
than the declared Earth PFT physiology: it is an ecosystem-level scalar fitted
to Earth's realised global carbon cycle.

There is no defensible Vesper observation with which to retune it. BIO-26
therefore registers and executes a one-factor sensitivity without tuning to the
desired productivity, and propagates the response as model-form uncertainty.

## 7. Dormant modules need activation guards

These do not block the natural-vegetation, methane-off baseline, but they are
latent Earth assumptions rather than safe generic code:

- Crop sowing, harvest and test dates retain 180/364/365-day ordinals. BIO-27
  makes cropland fail closed until those dates are expressed in the Vesper
  calendar and tested.
- The methane module fixes gravity at 9.81 m/s2, atmospheric pressure at
  101325 Pa, atmospheric O2 at 209000 micro-atm and CH4 at 1.7 micro-atm; one
  production ratio is explicitly tuned to reproduce global Earth emissions.
  BIO-28 is the planetary activation guard. The subsequent dedicated audit in
  `wetlands-peat-methane-audit.md` also found a literal 40 N peat/wetland
  classifier, non-conserving wetland saturation, fixed peat structure,
  incomplete restart state and an incomplete surface/atmospheric methane
  budget; WET-1 through WET-11 own those prerequisites before `run_peatland` or
  `ifmethane` can be enabled.

## 8. Latitude is geometry, not environmental state

A follow-up inventory found several places where geographic latitude selects an
Earth ecological regime rather than contributing to geometry:

- `framework/guess.cpp` and `modules/soil.h` split detailed peat from simplified
  wetland physics at exactly 40 N. WET-2 and WET-6 own the replacement.
- `framework/externalinput.cpp` assigns missing or fixed-default crop stands to
  tropical versus temperate crop PFTs at absolute latitude 30 degrees. BIO-27
  owns this along with the crop-calendar paths.
- `modules/simfire.cpp` uses absolute latitude 50 degrees to distinguish barren,
  shrub and tundra states. `modules/blaze.cpp` uses 30- and 50-degree bands to
  choose mortality and litter-tuning functions. FIRE-6 owns these effects-layer
  substitutions; vegetation, fuel and weather state must select the mechanisms.
- Natural phenology, litter release and the snow-depth establishment check use
  the sign of latitude to choose fixed northern or southern calendar days or
  months in `driver.cpp`, `growth.cpp`, `somdynam.cpp` and `vegdynam.cpp`.
  BIO-21 already owns their replacement with forcing-derived seasonal phase.
  Crop sowing has further sign-based Earth-calendar branches under BIO-27.

The defect is not the numerical value of any one boundary. No static latitude,
absolute latitude or hemisphere flag can stand in for temperature, season,
wetness, permafrost, vegetation, productivity, fuel or disturbance regime on
Vesper. Moving the lines would preserve the hidden Earth geography.

Direct latitude remains appropriate in a small, documented allowlist: spherical
cell area and orbital/daylength geometry; coordinate lookup, indexing and
serialization; and explicitly labelled spatial aggregation or diagnostics.
Those uses describe location or geometry and do not choose ecological physics.
Where seasonal orientation is required, it must be derived from the registered
orbit and forcing rather than assumed from the sign of latitude.

BIO-29 establishes this latitude-use contract across the vendored source and
the Vesper biosphere adapters. It inventories every direct latitude branch,
requires each use to be either allowlisted geometry/coordinates/diagnostics or
replaced with forcing, process state or traits, and adds a source-level check so
new ecological latitude classifiers fail without running LPJ-GUESS. Outside the
biosphere, PHYS-13 already owns Orogen's Earth-calibrated glaciation latitude
threshold; labelled latitude-band analyses such as regional summaries are
diagnostics, not model classifiers.

## Ordering

BIO-21 and BIO-22 are correctness blockers before a natural-vegetation result is
called meaningful. BIO-23 through BIO-25 close known forcing and unit mismatches
before final interpretation. BIO-26 is an uncertainty measurement. BIO-27 and
BIO-28 are fail-closed activation guards and do not block the current baseline.
BIO-29 is the cross-cutting guard: its inventory and lint are immediate
no-simulation work, while each replacement remains with BIO-21, BIO-27, FIRE-6,
WET-2 or WET-6 as appropriate.

Fire calendar, forcing and calibration assumptions remain in the dedicated
fire-model review; finding 8 only records their latitude branches in the shared
contract.
