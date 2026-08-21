# Implicit Earth assumptions below the LPJ-GUESS port

**Scope:** the Vesper adapter and the vendored LPJ-GUESS-CNP source. This was a
read-only source audit; no model build or simulation was run. BLAZE and GLOBFIRM
are deliberately excluded here so that fire can be reviewed as one mechanism
rather than scattered across the general port.

The existing porting note records the deliberate use of Earth PFTs, the
181-day/24-hour-step calendar, orbital geometry, PAR energy fraction, degree-day
rescaling, spin-up duration and the known frost limitation. Those are not
reopened. The findings below are assumptions that remained implicit underneath
those declared choices.

## 1. Natural phenology still contains Earth ordinal dates

`framework/guess.h` defines the northern coldest day as 14 (15 January) and the
southern coldest day as 195 (15 July), and aliases each hemisphere's warmest day
to the other's coldest day. Day 195 does not exist in the compiled 181-day
calendar.

These are active natural-vegetation paths, not unused crop machinery:

- `modules/driver.cpp` uses them to reset the summergreen GDD counter and to
  disable and re-enable chilling detection.
- `modules/growth.cpp` uses them to reset annual leaf-on accumulation.
- `modules/somdynam.cpp` releases summergreen leaf litter during month 0 in the
  north and month 6 in the south, still described as January and July.

The southern GDD and leaf-on resets therefore never occur. In the north,
chilling detection is disabled on day 14 and never re-enabled on the unreachable
day 195. Because `chilldays` is incremented with an inclusive upper-bound test,
losing its seasonal reset can also advance it beyond the last valid `Pft::gdd0`
index on a sufficiently cold cell.

The replacement must derive seasonal landmarks from the forcing or its
registered orbital phase and must cover both hemispheres, the equator and a
forcing sequence with more than one orbit. Merely scaling 14 and 195 by the year
ratio would retain an Earth thermal-lag assumption. This is BIO-21.

## 2. A model year is still treated as an Earth year in ecological rates

`build_vesper_pfts.py` currently rescales the two degree-day sums and three
duration controls. Many other quantities labelled per year still execute once
per 181-day orbit:

| mechanism | source examples | present interpretation |
| --- | --- | --- |
| leaf, root and sapwood turnover | PFT `turnover_*` | fraction per orbit |
| leaf and tree longevity | `leaflong`, `longevity`, individual age | model orbits |
| establishment and growth-efficiency mortality | `est_max`, `greff_min`, `NYEARGREFF` | per orbit / five orbits |
| litter and SOM turnover | `TAU_LITTER`, `TAU_SOILFAST`, `TAU_SOILSLOW` in `somdynam.cpp` | same decay fraction per orbit |
| CNP P-pool kinetics | `USORB`, `USSORB`, `UOCC` in `somdynam.cpp` | published annual rates are divided by the model-year length, while the source comments incorrectly label them per day |

For an absolute-time calibration, retaining the Earth annual fraction on every
Vesper orbit makes the process about 2.02 times as fast per Earth year. But not
every annual event should be rescaled: deciduous shedding and establishment may
legitimately occur once per seasonal cycle. Fractions also require a rate
conversion such as `1 - (1-r)^f`, not a linear multiplier, where `f` is the
orbit/Earth-year ratio.

The fix is therefore an explicit registry of each annual quantity as
absolute-time, seasonal-cycle, accumulated-flux or diagnostic, followed by
parameter-specific conversions and no-simulation checks. Fire rates are outside
this audit and registry pending the dedicated fire review. This is BIO-22.

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
their reporting interval as well.

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
  BIO-28 wires planetary constants and local pressure and requires an explicit
  methane calibration decision before `ifmethane` can be enabled.

## Ordering

BIO-21 and BIO-22 are correctness blockers before a natural-vegetation result is
called meaningful. BIO-23 through BIO-25 close known forcing and unit mismatches
before final interpretation. BIO-26 is an uncertainty measurement. BIO-27 and
BIO-28 are fail-closed activation guards and do not block the current baseline.

Fire calendar, forcing and calibration assumptions are intentionally unresolved
here and will be assessed together in the dedicated fire-model review.
