# The ecological forcing field contract

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf; its climate is ExoPlaSim and its biosphere is LPJ-GUESS. Everything
below is about what those two models mean by the fields they hand each other.

`biosphere/notes/ecological-climate-forcing-audit.md` establishes that the
12-bin driver is scaffolding and that the replacement must carry explicit
intervals. This document settles the other half, which the transport cannot
settle for itself: for every field crossing the seam, WHAT IT MEANS. A transport
that carries an agreed set of bytes whose meaning is not agreed is a faster way
to disagree.

The absolute day, the Earth year and the orbit are defined once, in
`biosphere/notes/time-base-unit-contract.md`, and every "per day" below is that
document's absolute day of 24 hours. This document does not restate it and does
not redefine it.

## What a row of the contract has to say

Six things, because five of them have each been the whole of a real defect
somewhere in this tree:

1. **Unit**, including the multiplier. "W/m2" and "hPa" are units; "1" is not a
   unit when the values run to 100.
2. **Time base**: which clock, and over what span. A rate per real second, a
   total per absolute day and a mean over a producer interval are three
   different numbers.
3. **Accumulated, instantaneous, or extremal.** ExoPlaSim's regular stream is
   none of the obvious one of these: see below.
4. **Area basis**: per square metre of what. Ground, leaf, land fraction, or
   whole gridcell.
5. **Sign convention.** ExoPlaSim stores upward fluxes negative and says so
   nowhere in the postprocessed metadata.
6. **Which side converts.** The rule, stated once and applied throughout: the
   PRODUCER converts into the contract unit, the CONSUMER converts out of it,
   and nothing converts twice. A field whose row does not name a converting side
   is not in the contract.

## The producer, as it actually behaves

### The regular stream is a mean, not a sample and not a total

`outmod.f90:outaccu` adds every field into an accumulator at every timestep and
increments `naccuout`; `outgp` divides by `naccuout` and writes. A record in the
regular stream is therefore the arithmetic mean over the timesteps since the last
write. It is not an instantaneous sample, and it is not a total: a precipitation
record is a mean RATE in m/s, and turning it into a depth is the consumer's
multiplication by an interval length the producer has to declare.

Three groups do not follow that rule and each is a trap of its own:

- `tempmin`, `tempmax`, `atsami`, `atsama` are extrema over the same window, not
  means. `outreset` gives the first two the sentinels `1.0e3` and `0.`, which is
  why `exoplasim/scripts/restart_schema.py` carries a `model_reset` column at
  all.
- `asndch` is never reset and accumulates across output intervals for the whole
  run (`outmod.f90`, the commented-out reset).
- The seasonal snapshot stream is genuinely instantaneous, and averaging
  snapshots across orbits gives a phase composite rather than an interval mean.

### Upward is negative

`radmod.f90` builds the shortwave upward flux as `zfu1 = -zt1*zrl1*z1mrabr` and
then forms `dswfl = dfu + dfd`, so a net flux is a SUM of a positive downward
term and a negative upward one. Every upward-signed field inherits this:
`ssru`, `stru`, `rsut`, `hfls` and `evap` are all less than or equal to zero in
the postprocessed product, while `rss` and `rst` are net downward and positive.

The consequence that matters: **incident surface shortwave is `rss - ssru`, not
`rss + ssru`.** Both are dimensionally sound and one of them is twice the other
plus a sign flip.

### Two precipitation partitions, cross-cutting, and one of them is not additive

`rainmod.f90` forms the three diagnostics as

    dprc = zprc + zprsc      convective, both phases
    dprl = zprl + zprsl      large scale, both phases
    dprs = zprsc + zprsl     solid, both processes

So `prl + prc` closes on the total, and `prsn` is a SUBSET of that same total
spanning both processes. Adding all three double counts the snow. Rain is
`pr - prsn`, and there is no separate rain diagnostic.

### `tas` is not a two-metre diagnostic

`fluxmod.f90` sets `dtsa = dt(:,NLEV) * sigma(NLEV)**zexp`: the lowest model
level's temperature brought to surface pressure along a dry adiabat. pyburn
names it `air_temperature_2m`. It is a near-surface reduction of the lowest
level, and `lib/lapse.py` is where the height of that level lives.

### The postprocessed extrema are extrema of a different variable

`outaccu` builds `tempmax`/`tempmin` (codes 321/320, `maxt`/`mint`) from
`dt(:,NLEP)`, the SURFACE temperature. It builds `atsama`/`atsami` (codes
201/202) from `dtsa`, the near-surface air temperature. Only the second pair is
an air-temperature range.

This is checkable and was checked, on
`exoplasim/analysis/climatology/bootstrap_regular_climatology.nc`: `mint <= ts <=
maxt` holds in all 24,576 cell-bins, and `tas` falls OUTSIDE `[mint, maxt]` in
15,561 of them, by as much as 28.3 K. Extrema bracket the mean of the variable
they are extrema of, so this identifies the variable without ambiguity.

Codes 201 and 202 are written by the model unconditionally, and are now in
`pyburn.ilibrary` as `tasmax`/`tasmin` and in `run_exoplasim.REGULAR_CODES`, so
a product postprocessed since `world-j0az` carries both pairs under names that
say which variable each is an extremum of. A product postprocessed before it
carries only the surface pair. What still reads the surface pair as though it
were a diurnal air range is the driver seam, below.

### Fields the contract may not name, because nothing produces them

- **`td2m`, code 168** (`world-dy7a`). Present in both `REGULAR_CODES` and
  `SNAPSHOT_CODES`.
  `outmod.f90` never writes it, pyburn does not derive it, and it is silently
  absent from the product rather than raising. Confirmed absent from the
  bootstrap climatology.
- **`uas`/`vas`, codes 165/166** (`world-1qxu`). Never written. There is no 10 m wind in this
  model's output, and no wind COMPONENT pair at all near the surface. What a
  wind consumer needs is a speed rather than a vector, and the ecological
  stream carries it: see the wind row below.

### `spd` is the speed of the mean vector, at every model level

pyburn derives code 259 as `sqrt(ua**2 + va**2)` from the divergence and
vorticity records, which in the low-I/O regime are already interval means. So
`spd` is |mean vector| and not mean speed, the two differ whenever the wind turns
within an interval, and the field carries a level axis rather than a
near-surface value. A wind consumer that takes the lowest level inherits both
properties. The gap is one-sided -- |mean vector| <= mean speed always -- so
anything driven by wind speed rather than by momentum transport reads LOW, and
no postprocessing of a mean vector recovers the difference. pyburn stamps this
onto the field as a `comment` attribute rather than leaving it to be
rediscovered.

The mean SPEED is a producer-side quantity and the ecological stream carries
it: `outmod.f90:ecoaccu` accumulates `sqrt(du(:,NLEV)**2 + dv(:,NLEV)**2)` at
every timestep into `aecowind` and `ecogp` divides it by the interval's own
step count, so what the stream writes as code 614 is the interval mean of the
speed and not the speed of the interval mean.

### `hur` is a percentage, and is not the interval mean of relative humidity

pyburn computes `rh = qq/zqsat * 100.0` and clips to `[0, 100]`. The values are
a percentage and `ilibrary` declares `%`; it declared `1`, which is what
`world-hf12` was. `zqsat` is formed on the model's own saturation branch, over
ice below `tmelt` and over liquid water at or above it, through the same
coefficients `plasimmod.f90`'s `ra1s`/`ra2s`/`ra4s` select and `rainmod.f90`
calls at every one of its saturation sites.

What remains true of the field, and is written onto it as a `comment`
attribute: it is a nonlinear function of the INTERVAL MEAN temperature and the
interval mean specific humidity, so it is not the interval mean of relative
humidity. A climatology postprocessed before the units fix carries percentage
values under a `1` attribute, so a reader checks the range rather than the
attribute.

### The bins are not the calendar

pyburn reduces the raw records to twelve bins by `np.linspace(0, ntimes,
13).astype(int)`, which distributes the remainder by truncation rather than
evenly. `lib/climatology.py` derives the counts and is the only place that
arithmetic belongs. In the clean regime the counts are `[15]*5 + [16] + [15]*5 +
[16]`, so two bins in twelve are longer than the other ten, and neither of them
is where a calendar would put a long month.

Two further properties of the producer's time axis:

- The `time` coordinate of a pyburn product is in TIMESTEPS, not seconds and not
  days.
- At the clean write interval the 182 records span 5824 of the orbit's 5850
  steps, so 0.44% of every orbit is in no bin at all. The bins tile the covered
  span, not the orbit.

### The gridcell is wholly land or wholly ocean

`lsm` takes only the values 0 and 1. A surface flux is therefore per square
metre of ground and the cell is entirely that ground; there is no land fraction
to divide by at this seam. This is not true of the Orogen mesh, where land is a
partial area per cell, and the two must not be confused: `build_lpj_driver.py`
integrates soil codes over the mesh's partial land areas and takes its land
DEFINITION from the climatology's binary mask.

## The consumer, as it actually behaves

`climate.temp` is a daily mean air temperature in degrees C. `climate.prec` is a
total depth in mm over one absolute day, all phases together;
`soilwater.cpp:snow` then sends ALL of it to the snowpack when the daily mean air
temperature is below 0 C and none of it otherwise, so a mixed or transient
interval becomes all one phase.

`climate.insol` under `NETSWRAD_TS` is a flux in W/m2 averaged over the whole
timestep, and `driver.cpp` multiplies it by a hardcoded `24 * 3600` seconds.
That constant is correct under this port precisely because the port's day is the
absolute day; it is not a claim about Vesper's 30-hour rotation, and PCAR-1 owns
the places where the rotation has to appear.

From there `climate.rad = insol * 86400` J/m2 of GROUND per absolute day,
`climate.par = rad * FRADPAR`, and `canexch.cpp` forms `apar = par * fpar *
alphaa(pft)`, still per square metre of ground. The canopy's conversion from
ground-area energy to absorbed photons is entirely LPJ's, and the producer must
not pre-apply any part of it.

`climate.dtr` is the diurnal AIR temperature range in degrees C. Its only reader
is `bvoc.cpp`'s `daytime_temp`, and `cfinput.cpp` builds it as
`dmax_temp - dmin_temp` from air temperature extrema, which fixes its meaning
beyond argument.

`gridcell.dNH4dep` and `dNO3dep` are kgN per square metre of ground per absolute
day. The driver file declares kgN/ha per Earth year and `vesperinput.cpp` does
the whole conversion, which is the pattern every row below follows.

## The contract

Contract units are what the forcing artifact carries. "Producer" means the
ExoPlaSim side, up to and including `build_lpj_driver.py`; "consumer" means
`vesperinput.cpp` and what it hands the model.

### Radiation

| field | contract unit | time base | kind | area | sign | converts |
| --- | --- | --- | --- | --- | --- | --- |
| incident surface shortwave | W/m2 | mean over the declared interval | accumulated then divided | ground | positive downward | producer, as `rss - ssru` |
| upward surface shortwave | W/m2 | mean over the declared interval | accumulated then divided | ground | positive upward | producer, as `-ssru` |
| net surface shortwave | W/m2 | mean over the declared interval | accumulated then divided | ground | positive downward | producer, `rss` unchanged |
| net surface longwave | W/m2 | mean over the declared interval | accumulated then divided | ground | positive downward | producer, `rls` unchanged |
| upward surface longwave | W/m2 | mean over the declared interval | accumulated then divided | ground | positive upward | producer, as `-stru` |
| mean cosine solar zenith | 1 | mean over the declared interval | accumulated then divided | n/a | n/a | neither; `czen` is already the night-inclusive mean of `gmu0` and is not a daytime mean |

Three shortwave terms rather than one, because they answer different questions
and one number cannot answer both. The surface energy balance and equilibrium
evapotranspiration want the NET term. The photon supply to a canopy wants the
INCIDENT term: the net term has already had the surface's own reflection removed
at the albedo ExoPlaSim was run with, which on a vegetated gridcell is largely
the canopy's own reflection, so a canopy model that then applies its own
absorption applies it to light the canopy has already been charged for. The
adapter must not name a net surface energy flux "insolation" and must not choose
between the two on the consumer's behalf. PCAR-1, PCAR-2 and BIO-25 own what the
canopy does with the incident term.

### Temperature

| field | contract unit | time base | kind | area | sign | converts |
| --- | --- | --- | --- | --- | --- | --- |
| near-surface air temperature | K | mean over the declared interval | accumulated then divided | n/a | n/a | producer supplies K, consumer subtracts 273.15 |
| near-surface air temperature maximum | K | extremum over the declared interval | extremal | n/a | n/a | producer, from `atsama` (code 201) |
| near-surface air temperature minimum | K | extremum over the declared interval | extremal | n/a | n/a | producer, from `atsami` (code 202) |
| surface temperature | K | mean over the declared interval | accumulated then divided | n/a | n/a | producer, `ts`/`tsa` |
| surface temperature maximum | K | extremum over the declared interval | extremal | n/a | n/a | producer, from `tempmax` (code 321) |
| surface temperature minimum | K | extremum over the declared interval | extremal | n/a | n/a | producer, from `tempmin` (code 320) |

Both pairs of extrema are in the contract and they are named apart. The
air-temperature pair is what `climate.dtr` means. The surface pair is what the
current driver's fourth array actually holds, and it is a legitimate field with
its own uses; it is simply not a diurnal air-temperature range.

### Water

| field | contract unit | time base | kind | area | sign | converts |
| --- | --- | --- | --- | --- | --- | --- |
| total precipitation | mm per absolute day | mean rate over the declared interval, expressed per absolute day | accumulated then divided | ground | positive downward | producer, `pr * 1000 * 86400` |
| solid precipitation | mm per absolute day | as above | accumulated then divided | ground | positive downward | producer, `prsn * 1000 * 86400` |
| liquid precipitation | mm per absolute day | as above | accumulated then divided | ground | positive downward | producer, as total minus solid |
| convective precipitation | mm per absolute day | as above | accumulated then divided | ground | positive downward | producer, `prc * 1000 * 86400` |
| large-scale precipitation | mm per absolute day | as above | accumulated then divided | ground | positive downward | producer, as total minus convective |
| evaporation | mm per absolute day | as above | accumulated then divided | ground | positive downward, so ExoPlaSim's `evap` enters negated | producer |

Both partitions travel, and the artifact carries the total plus one member of
each partition rather than all three of `prl`, `prc`, `prsn` as though they
summed. The consumer takes the phase from the producer and does not re-derive it
from a daily mean temperature at 0 C; BIO-13 owns what LPJ hydrology does with a
declared phase.

### Atmosphere

| field | contract unit | time base | kind | area | sign | converts |
| --- | --- | --- | --- | --- | --- | --- |
| surface air pressure | Pa | mean over the declared interval | accumulated then divided | n/a | n/a | producer; pyburn's `ps` is hPa and is multiplied by 100 |
| near-surface specific humidity | kg/kg | mean over the declared interval | accumulated then divided | n/a | n/a | producer, from `hus` at the lowest model level |
| near-surface wind speed | m/s | mean over the declared interval | accumulated then divided, at the lowest model level | n/a | n/a | producer, from `ecowind` (code 614) |

Specific humidity rather than relative humidity, deliberately. Specific humidity
is what the model carries and is linear, so its interval mean is a mean of the
quantity; relative humidity is a nonlinear function of temperature and pressure
whose interval mean is not the mean of the function, and pyburn's version of it
is additionally in percent and on the wrong saturation branch below freezing.
BIO-23 owns turning humidity into a vapour pressure deficit and does so from the
specific humidity, the pressure and the temperature it is given.

The wind row cannot be closed from the REGULAR product and is not a placeholder
any more. `spd` is the speed of the interval-mean vector, the mean speed is not
derivable from it, and there is no near-surface wind diagnostic in the regular
stream at all. The mean speed is a producer-side accumulation of
`sqrt(u^2 + v^2)` at the lowest model level, formed at every timestep before any
averaging, which is a quantity only the model can form; the ecological stream
forms it in `ecoaccu` and writes it as code 614. A forcing artifact takes the
wind row from there, and a wind row taken from `spd` instead is the placeholder
and must be labelled one.

### Coordinates and identity, carried once per interval

| field | contract unit | notes |
| --- | --- | --- |
| interval start, interval end | seconds of absolute time since the block's declared origin | exact, never inferred from a record number |
| interval duration | seconds | equal to end minus start; the two are both carried so a gap is detectable |
| orbital position | true anomaly, degrees | `nu`, code 50 |
| local solar phase | 1, fraction of a rotation | derived from the interval start and the rotation period, not from the calendar |
| grid identity, source build, source run | strings | the artifact is refused if they disagree between intervals |

### Static, carried once per artifact and never as weather

Cell area, the binary land mask, coordinates, and the surface identity fields.
`build_lpj_driver.py` already treats the soil code, regolith depth and bedrock
water fraction this way. Repeating a static field per interval is how a
provenance-bound quantity becomes something a consumer thinks it may interpolate.

## The identities that can fail

An artifact that satisfies the table above and none of these has not been
checked. Each of these has a right answer and can therefore be a test:

1. `liquid + solid = total` precipitation, to the artifact's stated tolerance.
2. `convective + large scale = total` precipitation, same tolerance.
3. `solid <= total` and every precipitation component is non-negative.
4. `incident - upward = net` shortwave, and `incident >= 0`, `upward >= 0`.
5. Each extremum brackets the interval mean of the variable it is an extremum
   OF, and only that one: `tasmin <= tas <= tasmax` for the air pair and
   `mint <= ts <= maxt` for the surface pair. This is the check that catches an
   extremum of the wrong variable, which is the defect it was written for.
   Neither pair is required to bracket the OTHER temperature, and requiring it
   would fail on correct data.
6. The declared intervals tile the declared block: no gap, no overlap, and the
   sum of the durations equal to the block's span.
7. The consumer's per-day series reproduces the producer's per-interval values.
   For an intensive quantity, the duration-weighted mean of the days in an
   interval returns the interval mean. For an extensive one, the sum of the days
   returns the interval total.
8. Every field the artifact declares is present, finite, and on the declared
   grid, and no field is present that the artifact does not declare.

`biosphere/scripts/check_forcing_contract.py` runs 1 to 5 and 8 against a
climatology or a forcing artifact and 6 and 7 where the artifact declares
intervals. It carries fixtures that are wrong in a named way, so a fixture that
does not get the verdict it was built for is a defect in the checker.

## What the current adapter does against this contract

Stated because the adapter is what runs today, and a contract whose violations
are not enumerated is a wish.

**Closed by this contract.** The interval-to-calendar mismatch. The driver
assigned each of the twelve bins `year_length // 12` absolute days with the
remainder on the last, which on a 183-day year is eleven months of 15 days and
one of 18, while the producer's bins span 15.08 days each except two of 16.09.
Bin 11 was therefore given 18 days of a rate measured over 16.09, inflating its
precipitation total by a factor of 1.119, and bin 5 was given 15 against 16.09,
deflating it by 0.932. The annual land-mean total moved by only 0.03% because
the errors nearly cancel, and that is exactly why it survived: the defect is a
seasonal redistribution of up to 12% within one twelfth of the year, largest at
the year boundary, which is also the seam a cyclic replay closes on.
`build_lpj_driver.py` now remaps the producer's intervals onto the model's
months conservatively, taking the interval spans from the climatology's own time
axis through `lib/climatology.py` rather than assuming them.

**Closed by this contract.** The fourth driver array is named for what it is.
`VESPDRV6` carries the same bytes as V5 and changes only the meaning of that
array, which is the change a magic exists to catch: it is the SURFACE
temperature range, and `vesperinput.cpp` no longer assigns it to `climate.dtr`,
whose only reader means an air-temperature range.

**Open, and at the driver seam.** The air-temperature extrema now reach the
regular product as `tasmax`/`tasmin`, so the range `climate.dtr` means has a
source. Carrying it needs a driver array the `VESPDRV` format does not have,
and `vesperinput.cpp` refuses `ifbvoc 1` until it does.

**Open, and owned elsewhere.** The single insolation field. The driver supplies
net surface shortwave and LPJ uses it for equilibrium evapotranspiration, where
net is right, and for PAR, where incident is right. Carrying both requires a
`Climate` member the vendored model does not have; PCAR-2 and BIO-25 own the
canopy end and EFOR-3 owns delivering the terms.

**Open, and owned elsewhere.** Precipitation phase. The artifact will carry it;
`soilwater.cpp` deciding phase from the daily mean air temperature at 0 C is
BIO-13's.

**Closed upstream.** The air-temperature extrema (`world-j0az`), the mean
near-surface wind speed (`world-1qxu`) and `hur`'s units and saturation branch
(`world-hf12`) were producer-side and are done: the first two reach a product
and the third is fixed in place. `td2m` (`world-dy7a`) is still declared and
never written. `world-ua4a` moves the BVOC refusal from the run into the gate
that is supposed to report it.

## What this document does not settle

- The transport. EFOR-1 chooses it and EFOR-3 builds it; this contract is what
  it has to carry, in whatever container.
- The cadence. EFOR-1 fixes the process intervals and EFOR-2 produces them; the
  interval rows above say what must be recorded per interval, not how long an
  interval is.
- Every process equation downstream of the seam. BIO-13, BIO-23, PCAR-1,
  FIRE-1, FIRE-2 and ANUT-4 keep their own work; what they gain here is one
  shared source sequence with agreed meanings instead of each reconstructing
  weather from a different reduction of it.
- The land column itself. Soil water, snow and the surface state are EFOR-8 and
  LSHY-4/LSHY-5, and this contract deliberately stops at the atmosphere's side
  of the surface.
