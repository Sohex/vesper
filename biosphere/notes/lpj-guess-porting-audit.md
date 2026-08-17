# Porting LPJ-GUESS to Vesper

What LPJ-GUESS assumes about Earth, where those assumptions live in the source,
and what has to change before the model says anything about this planet. Written
after obtaining and building the model, before any Vesper-specific code exists.

## The model as obtained

LPJ-GUESS 4.1.1, released 2021-10-13, Mozilla Public Licence 2.0, from Zenodo
record 8065737. No registration is required; the licence permits modification and
redistribution of modified sources under the same terms, which matters because
this port is a fork in all but name.

```
tarball   guess_4.1.1.tar.gz   md5 519c5da9397817438b93a8990228695e  (verified)
source    /home/cfutro/git/lpj-guess/guess_4.1
build     /home/cfutro/git/lpj-guess/build
binary    /home/cfutro/git/lpj-guess/build/guess
```

It lives beside `ExoPlaSim` and `planet_heightmap_generation` in `~/git` for the
same reason those do: it is third-party source that this project modifies but
does not own, and it is not this repository's history to keep.

Build state: configures and compiles clean with CMake 4.4.2 and GCC 16.1.1, zero
errors, for a codebase last touched in 2021. MPI is found and enabled. netCDF is
**not** installed on this machine, so the CF input module is compiled out. That
is not a gap: the CF module reads CMIP-style Earth datasets, and this project
needs its own input module regardless, so netCDF is only worth installing if we
later want it for output.

Smoke test, the bundled demo (3 gridcells, 500 spinup + 50 simulation years,
cohort mode): completes in 73 s and produces the expected `lai.out`, `cmass.out`,
`anpp.out` and the rest, with boreal needleleaf dominating Lapland and temperate
broadleaf dominating Scotland. The model works as shipped.

## Where Earth is baked in

Better than feared. The astronomy is concentrated in one function, and the
calendar is a runtime member rather than a compile-time shape.

### Year length

| what | where | note |
| --- | --- | --- |
| `MAX_YEAR_LENGTH = 365` | `framework/guess.h:303` | the one constant that matters |
| `year_length()` returns it | `framework/guess.h:456` | so the whole model agrees by construction |
| `ndaymonth[12]` set to Gregorian | `framework/guess.h:352` | a member array, settable at runtime |
| `DAYS_PER_YEAR = 365` | `modules/spinupdata.h:229` | spin-up climate buffer |
| `daily_fluxes_patch[365]`, `daily_fluxes_pft[365]` | `framework/guess.h:1313,1316` | literal, must track |
| `for (d=0; d<365; d++)` | `modules/soil.cpp:189,1775` | peatland water-table averaging |

The per-day arrays that actually carry state are sized `Date::MAX_YEAR_LENGTH`
already: `Climate::qo/u/v/hh/sinehh/daylength_save/doneday`, and
`Soil::dwcontupper/dwcontlower/dthaw/wtp`. Because Vesper's year is *shorter*
than Earth's, shortening `MAX_YEAR_LENGTH` needs no reallocation anywhere. The
arrays are oversized, not undersized. This is the single most convenient fact
about the port.

Everything else in the model reads `date.ndaymonth[]` and `date.year_length()`
consistently: `driver.cpp` monthly-to-daily interpolation, `canexch.cpp` monthly
LAI, `weathergen.cpp`. Those need no change at all.

### The annual soil thermal wave

`modules/soil.cpp:364-367` carries the analytic damped-oscillation soil
temperature model with the year length written into two constants:

```
const double HALF_OMEGA = 8.607E-3;   // omega/2 = pi/365
const double LAG_CONV   = 58.09;      // 365/(2*pi), oscillation lag in days
```

This is real physics tied to the orbital period, not a convention. Both must
move with the year length.

### The astronomy

All of it is in `modules/driver.cpp:841-1000`, `daylengthinsoleet()`, which
computes daylength, insolation and equilibrium evapotranspiration from latitude
and day of year. Prentice et al 1993 equations 7 to 14. The Earth-specific
quantities:

| constant | value | Vesper |
| --- | --- | --- |
| `QOO` solar constant | 1360.0 W/m2 | 1306.56 (0.96 x 1361) |
| eccentricity | 0.01675 | 0.02 |
| obliquity | 23.4 deg | 32.0 deg |
| solstice phase offset | 10.5 days from day 0 | set by `longitude_vernal_equinox_degrees` 102.7 |
| `daylength = 24.0 * hh / PI` | 24 h rotation | **unchanged**, see below |
| `K = 13750.98708` | 12/pi x 3600, angular units to seconds per day | **unchanged**, see below |
| `FRADPAR` | 0.5 | unresolved, see below |
| `BETA` global shortwave albedo | 0.17 | avoidable |

`BETA` is avoidable because it is only applied on the `SUNSHINE` and `SWRAD`
paths. Driving with `NETSWRAD_TS`, which is what ExoPlaSim's `rss` already is,
sets `net_coeff = 1` and bypasses it entirely. That is the right choice anyway:
this project computes surface albedo from lithology and knows it far better than
a global constant does.

### PFT bioclimatic limits

`gdd5min`, `twmin`, `tcmin_surv`, `tcmax_est` and the rest, in the shipped
`.ins` files. These are Earth calibrations and they are the least mechanical part
of the port, because two different things are tangled in them: temperature
limits, which are physiology and carry over, and degree-day sums, which are
accumulated per year and therefore depend on how long a year is.

## The calendar decision

This is the load-bearing choice and everything downstream inherits it.

Vesper at the 0.96 baseline flux:

```
orbital period      180.655 Earth days = 4335.7 hours
rotation             30 hours
sidereal rotations  144.524 per orbit
solar days          143.524 per orbit, each 30.209 hours long
```

LPJ-GUESS steps one "day" at a time and every rate constant in it is per day.
There are two ways to map that onto Vesper and they fail differently.

**Option A, Earth-day steps.** A step is 24 hours, the year is 181 steps.
Absolute time stays honest: 181 x 24 = 4344 h against the true 4335.7 h, an error
of 0.19%. Every per-day rate constant in the model, maintenance respiration,
soil decomposition, snowmelt, nitrogen turnover, remains calibrated against the
absolute time it was calibrated against. The cost is that the diurnal cycle is
misrepresented, but the model only ever sees the diurnal cycle through
`climate.daylength`, so the error is confined to one variable.

**Option B, Vesper-day steps.** A step is 30.209 hours, the year is 144 steps.
The light and dark alternation is then correct. The cost is that every per-day
rate constant is now 25.9% wrong in absolute time unless each one is found and
rescaled by hand, across respiration, decomposition, hydrology, phenology and
crop sowing.

**Take Option A.** It confines the error to a single well-understood variable
instead of scattering it across every rate in the model, and it is the one that
can be checked, because the residual error has a known sign.

Under Option A, `daylength` is the lit fraction of the diurnal cycle multiplied
by 24 h, which is what `driver.cpp` already computes, so nothing has to change.

**An earlier version of this note claimed instantaneous irradiance is then 25.9%
too bright, and that was wrong.** It is exact. The error came from thinking of a
step as "a Vesper day" rather than as 24 hours of absolute time, and the same
mistake produces the opposite recommendation, that daily energy should be scaled
up by 1.25 to avoid light-starving the trees. Both are wrong and for the same
reason.

The photosynthesis routine uses the daily PAR total divided by daylength,
`canexch.cpp:718`, so what reaches it is the mean irradiance over daylight hours.
Energy and daylength both carry the same day-length factor and it cancels:

```
real   F x 30.209 x 3600 J/m2  over  f x 30.209 h   ->  F x 3600 / f
model  F x 24     x 3600 J/m2  over  f x 24     h   ->  F x 3600 / f
```

Ratio 1.000000. `F` is the same number in both, because ExoPlaSim reports a mean
flux over the timestep and `NETSWRAD_TS` takes it as exactly that. Scaling the
energy up by 1.25 without scaling daylength would over-supply light by 25% and
break the annual budget as well.

The annual total is right for the same reason: 181 steps of 86400 s is 0.19%
more absolute time than the true 180.655-day orbit, which is the year-length
rounding and nothing else.

The degree-day problem is real and separate. A Vesper year has 181 steps against
Earth's 365, so annual GDD sums are roughly halved for the same temperatures, and
Earth `gdd5min` thresholds would exclude almost every tree PFT for reasons that
have nothing to do with the climate. They should be rescaled by
180.655 / 365.2569 = 0.4946 so that a threshold means the same absolute amount of
growing season it meant on Earth. Temperature limits (`tcmin_surv`, `twmin`) are
physiology and do not scale. One second-order effect: a Vesper month is 15.1 days
rather than 30.4, so "coldest month mean temperature" is averaged over a shorter
window and will read slightly more extreme than the Earth-calibrated limits
expect.

## PAR fraction: resolved, and it found a climate bug

`FRADPAR = 0.5` is the fraction of shortwave radiation that is
photosynthetically active. It multiplies straight into GPP, so an error in it is
an error of the same size in every productivity number the model produces.

Under a K dwarf it is lower than the Sun's, because the spectrum is shifted into
the near infrared. Deriving how much lower turned up the fact that ExoPlaSim's
`k2.dat` is not a K dwarf at all: it is the star K2-18, an M2.5V at about 3450 K,
not the spectral type K2. The full finding, its evidence and its consequences for
the climate runs are in `exoplasim/notes/stellar-spectrum-audit.md`.

For this component the resolution is:

- The k2 file puts **0.078** of shortwave in 400-700 nm. That is the wrong star
  and must not be used.
- A correct K2.5V spectrum now exists as `exoplasim/inputs/stellarspectra/k25v`,
  BT-Settl interpolated to 4965 K. It puts **0.309** there against the Sun's
  0.390 computed the same way. The ratio is 0.793, so **`FRADPAR` is 0.396**.

Use **0.40**, derived from the real spectrum rather than a blackbody. For
reference a 4965 K blackbody would have given 0.43, so line blanketing is worth
about 8% of the value and was worth computing properly.

The residual caveat is that `FRADPAR` in LPJ-GUESS is a surface quantity, not a
top-of-atmosphere one, and Earth's 0.5 is higher than the Sun's 0.367 at the top
of the atmosphere precisely because atmospheric near-infrared absorption strips
NIR on the way down. Scaling by the top-of-atmosphere ratio assumes that
enrichment is the same under a redder star, which it is not exactly. Treat 0.427
as good to about 10%, not better.

## What the input module has to supply, and where it comes from

| LPJ-GUESS input | source | state |
| --- | --- | --- |
| monthly mean temperature | `tas`, 12 bins | ready |
| monthly precipitation | `pr` | ready |
| insolation as `NETSWRAD_TS` | `rss` | ready |
| cloud fraction (alternative path) | `clt` | ready |
| atmospheric CO2 | `config/planet.yaml`, 450 ppm | ready, a constant |
| soil texture class | Orogen `substrate_class` via lithology | needs a mapping, does not exist yet |
| nitrogen deposition | no source | must be declared, as the albedo constants were |
| gridlist | T42 land cells from `surface_class` | 4106 cells |
| daylength, insolation geometry | computed internally, see above | needs the patch |

Ranges over land in the existing climatology, as a sanity check that nothing is
degenerate: monthly temperature -40.15 to 48.32 C, precipitation 0 to 22.80
mm/day, net surface shortwave 0 to 412 W/m2, cloud fraction 0 to 0.988. Land
annual means 9.29 C and 2.196 mm/day.

The two entries with no source are worth naming as such rather than defaulting
quietly. Soil texture is derivable, and lithology is exactly the right thing to
derive it from, but the mapping from Orogen's rock classes to LPJ's nine soil
codes is a judgement call that should be written down. Nitrogen deposition has no
analogue at all on a world with no biosphere history and no industry, and the
honest handling is a declared constant with the sensitivity reported, in the same
style as the evaporite albedo.

## Cost

4106 land gridcells at T42. The demo runs 550 years in about 24 s per gridcell
single-threaded, so a full planet is roughly 27 hours serial, or about 1.7 hours
across the 16 physical cores. Cheap enough that the vegetation loop is not the
bottleneck; ExoPlaSim is.

## What blocks the first real run

The climatology LPJ-GUESS should consume does not exist yet. `climatology_s096`
is on `precarve-unzoned` with evaporite as one class at 0.50, and
`world_state.json` already labels it "never as Vesper's climate". The
`carved-zoned` baseline has converged at 295.18 K but has no described
climatology.

So the first LPJ-GUESS run is a shakedown on the pre-carve climatology, to prove
the input module, the patch and the calendar, and it must be labelled as such.
The first run that means anything waits on a `carved-zoned` climatology.

## The patch, written

`patches/lpj-guess-4.1.1-vesper.patch`, against the pinned 4.1.1 tarball, in the
same style as `exoplasim/patches/exoplasim-3.4.2-star-cycle.patch`.

```bash
cd /home/cfutro/git/lpj-guess/guess_4.1
patch --forward --strip=1 --directory=. < <world>/biosphere/patches/lpj-guess-4.1.1-vesper.patch
cd ../build && make -j16
```

Every planetary constant is collected in one block in `framework/guessmath.h`,
which `guess.h` and `spinupdata.h` both already include, so the assumptions are
auditable in one place rather than scattered across five files:

| constant | value | replaces |
| --- | --- | --- |
| `VESPER_YEAR_LENGTH_DAYS` | 181 | 365, in three places |
| `VESPER_STELLAR_CONSTANT` | 1306.56 W/m2 | 1360 |
| `VESPER_ECCENTRICITY` | 0.02 | 0.01675 |
| `VESPER_OBLIQUITY_DEG` | 32.0 | 23.4 |
| `VESPER_SOLSTICE_OFFSET_DAYS` | 151.2 | 10.5 |
| `VESPER_FRADPAR` | 0.396 | 0.5 |

The solstice offset is not a guess. It was fitted to the solar declination the
ExoPlaSim climatology actually reports, `zdec` over its 12 bins: rms residual
1.11 degrees, maximum 1.95. **It assumes day 0 of the LPJ-GUESS year is the
first bin of the climatology**, which the input module has to guarantee or the
seasons run out of phase with the forcing.

### A correction to the audit above

The audit said `daylength` and `K` scale with the rotation period. Writing the
patch showed they do not, and both are left alone.

`daylength = 24.0 * hh / PI` computes hh/PI, which *is* the lit fraction of the
diurnal cycle, and multiplies by 24. Under the registered calendar a step is 24
hours of real time, so that is already exactly the intended "lit fraction times
24 h". Likewise `K` converts the angular integral to seconds per 24-hour step,
which is what a step now is. Driving with `NETSWRAD_TS` supplies a genuine mean
flux in W/m2, so the energy delivered per step is flux times 86400 seconds, and
that is correct in absolute time.

So the patch is smaller than the audit implied: geometry and calendar only.
`QOO` does still matter despite the `NETSWRAD_TS` path, because `w` and
`climate.qo` are used in the equilibrium-evapotranspiration longwave term
(`driver.cpp`, Eqn 19), not only in the `SUNSHINE` branch.

### Verified

The patched model builds clean and runs. Feeding it the bundled Earth demo data
is a deliberate control: same forcing, Vesper calendar, so anything that moves is
the calendar.

- Annual evapotranspiration falls to **0.492** of the unpatched value against an
  expected 181/365 = 0.496. The year length is exactly in effect.
- Annual GPP falls further, to 0.390, because the vegetation itself changed:
  boreal needleleaf and temperate broadleaf collapse to grass. That is the
  degree-day problem, arriving on cue. Earth `gdd5min` thresholds cannot be met
  in 181 days, which is why item 2 below is not optional.

## Order of work

Done:

1. ~~Patch `Date`, `driver.cpp` and `soil.cpp` for the calendar and astronomy.~~
2. ~~Rescale the PFT degree-day limits.~~ `build_vesper_pfts.py`, factor derived
   from the configured orbit, `gdd5min_est` 500 to 247.
3. ~~Write the input module.~~ `vesperinput`, driven by one generated binary,
   with day 0 pinned to climatology bin 0 by a fitted declination phase.
4. ~~Soil.~~ Superseded by the `pedology/` component, which weathers lithology
   under the climate and closes a loop through soil carbon. Nitrogen turned out
   to be bracketed by `nfix_a`, not by deposition.
5. ~~Parallelism.~~ `vesperinput` splits cells across MPI ranks itself, strided
   so latitude bands balance. Verified 5/5/4/4 over 18 cells on 4 ranks with no
   duplication and none dropped.

Remaining, in the order they block each other:

6. ~~A run harness.~~ `run_lpj_guess.py`.

   Superseded text kept for the reasoning:

   **A run harness.** Every run so far has been hand-assembled in a scratch
   directory with a hand-written instruction file. There is no `biosphere/runs/`,
   no run manifest, no provenance record. That is against this project's own
   convention and it is the next thing to build. It should generate the `.ins`
   too: the settings are retyped each time, and parallel mode needs absolute
   import paths because each rank chdirs into its own `runN/`.
7. ~~The feedback products.~~ `build_surface_albedo.py --mode modelled`.

   **The feedback products.** Turning `lai.out` and `fpc.out` into surface albedo
   (codes 174, 175, 176) and forest fraction (212) for ExoPlaSim. This is the
   entire point of the component and it does not exist yet: the climate is still
   forced with `land_albedo_source: vegetated`, a constant.
8. ~~Scoring.~~ `score_prediction.py`, with the ten registered bands
   transcribed from the note and the simulation-year conversion taken from the
   run manifest.

Blocked on the climate side: the driver is still built from
`climatology_s096`, which is pre-carve and pre-`k25v`.

## Cost of a full run

5.46 s of CPU per gridcell for 530 years at `npatch 5`, measured on 18 cells
across 4 ranks. For 4,106 land cells that is about 6.2 hours of CPU, so roughly
25 minutes of wall clock on 16 ranks. The vegetation is not the expensive half
of this iteration.
