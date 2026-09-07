# The simulated polar caps are cold-filtered, then limited by water CAPTURE and not supply

Measured 2026-09-07 on `lpj_1e6a2b9ca51a4eff9592992cad96677b`, the accepted
LPJ-GUESS run, against `exoplasim/analysis/climatology/baseline_regular_climatology.nc`
on build `canonical-10m-carve2`. This is a note about a simulated world: every
quantity below is a modelled field of Vesper, and every comparison to Earth is
a distance to report rather than a target.

The 234 simulated gridcells poleward of 75 degrees hold a mean summed foliar
projective cover of 0.030. Two terms produce that, in order, and neither is a
defect:

1. **The seasonal extreme removes ten of the twelve plant functional types
   outright**, and this is the model behaving correctly.
2. **The two survivors then capture almost none of the water that arrives.**
   Supply is NOT what is short: cells elsewhere on the planet with LESS
   precipitation than the cap support seven times its cover. What is short is
   the fraction of arriving water that reaches a plant, and that is a property
   of the plants rather than of the world.

## The seasonal range is about 100 degC, and it is a hard filter

Averaged over the 234 cells, each read in its OWN warmest and coldest month
rather than in a fixed one, beside 264 tropical cells read the same way:

| | polar | tropical |
| --- | --- | --- |
| warmest-month air temperature | 31.11 degC | 23.48 degC |
| coldest-month air temperature | -68.75 degC | 17.98 degC |
| warmest-month insolation, visible + near-infrared | 390.67 W/m2 | 211.84 W/m2 |
| warmest-month precipitation | 0.212 mm/day | 3.087 mm/day |
| months above freezing | 3.39 of 12 | 12 of 12 |
| months above 5 degC | 2.87 of 12 | 12 of 12 |
| precipitation per simulation year | 77.27 mm | 632.14 mm |
| share of it falling as snow | 52.85% | 0% |

**THESE ARE THE FIGURES LPJ WAS FORCED WITH, AND THEY BLEND 31 CELLS THE
CLIMATE MODEL CALLS OCEAN.** LPJ ran on BIO-11's 1617 rootable cells, of which
234 lie poleward of 75 degrees, while the atmosphere model's own binary land
mask holds only 203 of those; the driver was built at all 234 either way, so
the table above is what the biosphere integrated and is the right support for
every biosphere number in this note. Restricted to the model's land the summer
is HARSHER, not milder: warmest month 34.46 degC rather than 31.11, insolation
398.15 W/m2 rather than 390.67, precipitation 79.36 mm rather than 77.27 and
48.9% of it snow rather than 55.7%. Nothing in the argument turns on the
difference. What the distinction DOES decide is soil: temperature and moisture
are undefined over the model's ocean, so any statement about the root zone must
use the land mask, and an earlier draft of the companion audit did not.

`tcmin_surv` is a survival limit on the coldest month's mean temperature, and
at -68.75 degC it excludes every type that declares one. The most cold-tolerant
finite limit in `biosphere/generated/vesper_pfts.ins` is -31 degC, held by BNE
and BINE; IBS is -30. The measured cover carries the filter exactly:

| poleward of 75 deg | FPC | cmass | GPP | why |
| --- | --- | --- | --- | --- |
| BNS | 0.00265 | 0.00549 | 0.00158 | `tcmin_surv` absent -- survives |
| C3G | 0.02738 | 0.01634 | 0.01475 | `tcmin_surv` absent -- survives |
| the other ten | 0 | 0 | 0 | every finite limit is above -68.75 degC |

The ten are at EXACTLY zero rather than merely small, in all three quantities.
That is the signature of a hard filter rather than of a gradient, and it is the
answer to whether the seasonal extremes are killing the vegetation: for ten of
twelve types, yes, and correctly.

## The growing season is thermally generous, so warmth is not short either

Growing degree-days above 5 degC, per simulation year, against the
`gdd5min_est` thresholds in `biosphere/generated/vesper_pfts.ins` -- BNS and
IBS 175.356, BNE and BINE 250.509, TeBS 551.119, TeNE and TeBE 1002.03, both
grasses 0:

| band | mean GDD5 | median | share below BNS's 175.4 |
| --- | --- | --- | --- |
| poleward of 75 deg | 844.4 | 752.1 | 6.8% |
| 60 to 75 deg | 658.4 | 708.8 | 19.4% |
| 45 to 60 deg | 196.6 | 135.4 | 55.3% |
| equatorward of 15 deg | 2917.5 | 3001.4 | 0% |

The cap clears BNS's threshold by nearly five times, and clears BNE's and
TeBS's as well. **The 45-to-60 band is the genuinely cold-limited one**, at a
quarter of the cap's degree-days -- which is the 32-degree obliquity's
signature, since a high obliquity moves annual insolation poleward. So at the
cap neither cold-season survival for the two survivors, nor growing-season
warmth, nor light is what holds the cover down.

## The two survivors capture almost no water, and that is measured separately

Neither survivor declares a cold limit, so within the polar band the filter is
uniform and any variation across it is something else. It is water. Sorting the
234 cells by precipitation, cover rises monotonically and by a factor of
fifteen, Spearman +0.531:

| precipitation, mm/simulation year | mean cover | cells |
| --- | --- | --- |
| 40.3 - 56.2 | 0.0040 | 59 |
| 56.2 - 70.5 | 0.0194 | 58 |
| 70.5 - 83.3 | 0.0350 | 58 |
| 83.3 - 145.9 | 0.0613 | 59 |

### The quantity is TRANSPIRATION, and that is not this project's choice

Noy-Meir (1973) states the desert water balance as `T = P - R - D - E` --
transpiration as precipitation less runoff, drainage and soil evaporation -- and
names T "the component driving the energy flow to the biotic subsystem". So
reading a cap's productivity off its capture fraction rather than off its
precipitation is the standard treatment of an arid system and not a framing
invented here. At 154 mm per Earth year the cap sits in his ARID class, between
extreme arid below 60-100 mm and semiarid above 150-250.

His pulse taxonomy also classifies it. Where the interval between rain events is
much shorter than the system's relaxation time -- "rain events clustered,
markedly seasonal rainfall, slow response, e.g. in cool, winter-rainfall
deserts" -- the input pulses accumulate into a single response, and "the total
rain of the season may then be considered a single input pulse". That is his
Figure 2b, and it is this cap: over half the annual precipitation arrives as
snow through the dark season and leaves as one melt pulse.

### The water budget closes, so the partition is a measurement

LPJ-GUESS's own annual water fluxes against the precipitation the climatology
delivered, per simulation year:

| | polar | tropical |
| --- | --- | --- |
| precipitation in (climatology) | 77.27 | 632.14 |
| transpiration through plants (`aaet.out`) | 3.66 | 314.85 |
| bare-soil evaporation (`mevap.out`) | 29.14 | 54.35 |
| canopy interception (`mintercep.out`) | 0.01 | 53.11 |
| runoff (`tot_runoff.out`) | 45.04 | 205.83 |
| total out | 77.86 | 631.72 |
| residual | +0.8% | -0.03% |

Both close, which is what makes the shares a measurement rather than an
inference: the driver delivered the climatology's own precipitation and no term
of the balance is missing at either latitude.

**Plants get 4.7% of the polar water and 50% of the tropical water.** The rest
leaves two ways, and both are consequences of the cover being low rather than
causes independent of it:

- **58% runs off.** Over half the annual precipitation arrives as snow during
  the dark season, when nothing can use it, and leaves as a melt pulse at the
  start of the warm season, before there is any leaf to transpire it.
- **38% evaporates straight off the soil**, against 8.6% in the tropics,
  because there is no canopy shading the surface. Canopy interception is
  0.01 mm against the tropics' 53.11: there is effectively no canopy.

## Supply is not the constraint. Capture is.

Transpiration as a share of the precipitation that arrived, by band, with the
warm band chosen to be AS DRY AS THE CAP:

| band | cells | precip | transpired | captured | cover | coldest month |
| --- | --- | --- | --- | --- | --- | --- |
| poleward of 75 deg | 234 | 77.3 | 3.7 | **4.7%** | 0.030 | -68.8 degC |
| under 120 mm, equatorward of 60 | 81 | **73.9** | 58.3 | **78.9%** | **0.220** | +8.7 degC |
| 120 to 250 mm, equatorward of 60 | 278 | 187.0 | 86.6 | 46.3% | 0.344 | -11.0 degC |
| equatorward of 15 deg | 264 | 632.1 | 314.9 | 49.8% | 0.839 | +18.0 degC |

**The warm dry band receives LESS water than the cap -- 73.9 mm against
77.3 -- captures seventeen times as much of it, and carries seven times the
cover.** So the cap is not short of water in the sense a mass balance would
mean. It is short of plants able to intercept water that arrives frozen, in a
dark season, and departs as a melt pulse.

That reading replaces "the cap is water-limited", which is true only if
"water" is read as delivered-to-root rather than delivered-to-ground. The three
places the polar supply goes are all traits rather than boundary conditions:
the 58% that runs off is a timing mismatch against phenology, the 38% that
evaporates off bare soil is the absence of a canopy to shade it, and the 4.7%
that gets through is what Earth-shaped roots and Earth-shaped phenology manage
on a regime Earth does not have.

## What an adapted flora could win, bounded

Schwinning and Sala (2004) give the mechanism behind the partition above:
shallow pulses "wet only the uppermost cm of the soil, where a large fraction of
soil moisture is lost by direct evaporation, due to high temperatures and low
root densities", while "the deeper the pulse depth, the larger the fraction of
precipitation leaving the soil via transpiration". Their handle on it is that
initial infiltration into dry soil runs 0.24 to 1 cm per millimetre applied,
set by a saturated volumetric water content between 0.1 for sand and 0.4 for
clay.

**Pulse depth is therefore NOT what holds this cap at 4.7%.** A melt pulse of
roughly 39 mm infiltrates 9 to 39 cm on that relation, which reaches the model's
0 to 40 cm root layer. The run's own runoff decomposition says the same thing
from the other side:

| polar runoff, mm per simulation year | |
| --- | --- |
| surface | 24.95 |
| drainage below the root zone | 20.07 |
| baseflow | 0.01 |

Twenty millimetres passes BELOW the roots, so the water reaches root depth and
leaves rather than never arriving. The surface share, 55.4% against the
tropics' 50.9%, is not a strong frozen-ground shedding signature either. What is
distinctive is the baseflow: 0.01 against the tropics' 39.19, a dead deep store
under a column sitting at -31.85 degC the year round.

So the capturable terms are the 29.14 mm a canopy would stop evaporating off
bare soil and the 20.07 mm a deeper or faster root system would intercept
before it drains. **That bounds an adapted flora at about 49 of 79 mm, roughly
62% capture**, against 4.7% today and 78.9% in the warm dry band. It is a
ceiling and not a prediction, and it is the number a sizing arm should be judged
against: not "more than now" but "how much of 62%".

### Both terms leave in season, so the ceiling survives its timing test

The bound above is a mass argument and says nothing about WHEN either term
leaves. A term departing under frozen ground is recoverable by nothing, which
would make 62% an overstatement rather than a ceiling. Measured on the driver
LPJ-GUESS read, by `biosphere/scripts/polar_water_timing.py`:

| | share |
| --- | --- |
| bare-soil evaporation in months above 0 degC | 99.6% |
| bare-soil evaporation in months above 5 degC | 96.8% |
| bare-soil evaporation strictly before the thaw month | 0.0% |
| transpiration in months above 0 degC | 99.0% |
| precipitation DELIVERED in months above 0 degC | 39.9% |

The 29.14 mm leaves while the ground is thawed and while the two survivors are
transpiring. It is genuinely on the table.

**Thaw and the growth threshold arrive together.** Over the 224 of 234 cells
reaching both, the first month above 5 degC follows the first month above
freezing by 0.018 months, which is 0.27 days. There is no interval in which
water is liquid and growth is still barred, so no phenology trigger moved
earlier recovers anything: earlier is frozen. All 77.85 mm crosses the surface
inside the 3.39-month warm window, 45.85 mm of it as the melt of the pack that
accumulated below freezing and 32.00 mm delivered warm.

`tot_runoff.out` is annual, so the 20.07 mm drainage term's month is inferred
from the surface flux rather than measured.

### What the timing does bind is the RATE

The water is not spread across the growing season. It crosses the surface in
one interval of about fifteen days:

| | share of the annual total in its single largest month |
| --- | --- |
| bare-soil evaporation | 83.5% |
| transpiration | 45.0% |
| precipitation | 30.0% |

Transpiration is spread across the season and the water is not. So the
mismatch at the cap is not between the calendar of the melt and the calendar of
leaf-out, which coincide; it is that the year's water passes the surface in
about fifteen days and the canopy that would intercept it is not deployed
within them. At LAI 0.067 there is nothing to shade the soil during the one
month that matters. The trait that binds is therefore the RATE a canopy can be
built inside that window, and the stored carbon and nitrogen that fund it
before there is any photosynthesis to pay for it -- an allocation and storage
trait rather than a water-uptake trait.

One caveat travels with the runoff figures. `landmod.f90` records LSHY-5: the
climate model's soil layers carry no water phase, so melt water there always
infiltrates whatever the soil temperature, while LPJ-GUESS carries an ice
fraction per layer and reduces available liquid under freezing. The two columns
disagree about whether water that reached the ground is liquid, and the numbers
above are LPJ's side of that disagreement.

## The cap is not inefficient at capture. It has almost no leaf

The three capture-side traits proposed for sizing -- rooting depth and
distribution, the phenology trigger, and the bare-soil evaporation path -- are
all EFFICIENCIES: water captured per unit of leaf deployed. Each is worth moving
only if the cap's efficiency is deficient, which the seven-fold gap against warm
ground was read as showing. Measured by
`biosphere/scripts/polar_capture_efficiency.py`, it is not.

Reading capture per unit leaf area as `k = -ln(1 - captured) / LAI`:

| band | cells | LAI | capture | median k | percentile of the planet's cells |
| --- | --- | --- | --- | --- | --- |
| poleward of 75 deg | 224 | 0.070 | 4.4% | 0.716 | 72 |
| 60 to 75 | 314 | 0.233 | 11.1% | 0.573 | 64 |
| 45 to 60 | 245 | 0.549 | 9.6% | 0.121 | 18 |
| 30 to 45 | 276 | 2.131 | 47.8% | 0.275 | 40 |
| 15 to 30 | 190 | 1.819 | 79.8% | 1.245 | 87 |
| equatorward of 15 | 264 | 4.411 | 60.5% | 0.209 | 32 |

**k is not a constant**: p10 0.089, median 0.376, p90 1.420 over the 1513 cells
that carry both leaf and water. So capture is NOT a function of leaf area alone,
and a Beer's law in leaf area is refused rather than assumed. What the table
supports is the ranking, and the ranking is decisive: the cap turns leaf area
into captured water at percentile 72, above the tropics at 32 and far above the
genuinely cold-limited 45-to-60 band at 18. It takes water WELL for what it
carries. It carries 0.067.

### A faster phenology displays the same nothing sooner

`lai.out` carries `growth.cpp`'s `indiv.lai`, which is `cmass_leaf * sla`: the
leaf area at FULL display. What a cell shows on a given day is
`lai_today() = lai * phen`, never more. So 0.067 is the cap's ceiling with
phenology satisfied on the first day of the season, and `phengdd5ramp` moves
WHEN that area appears rather than how much of it there is. At 0.067 the ground
is 97% bare whenever the leaves come out.

That closes the sizing question negatively for every trait proposed. The three
capture-side parameters move an efficiency already above the planetary median,
the phenology trigger has no water to recover and no area to add, and none of
the four adds leaf.

### What the cap is short of is leaf carbon, and it accumulates slowly

Standing vegetation carbon at the cap is 0.0223 kgC/m2 against an annual net
primary production of 0.0100 kgC/m2, so what stands there is about twenty
simulation years of production rather than a stock being drained by turnover.
The loop is self-limiting rather than leaky: little leaf fixes little carbon,
which buys little leaf. Nothing in the twelve types' parameter set breaks into
it, which is the finding `world-orok` records and not a defect in the run.

## What the seven-fold gap still bundles

The gap between the cap and equally dry warm ground bundles the loss of ten of
twelve types with a growing season of 3.39 months against twelve. The third
term it was thought to bundle, the melt-pulse timing, is separated above and is
not a cost: the pulse arrives inside the growing season and the two survivors
transpire through it.

Separating the remaining two is not what decides whether this world carries
plant functional types of its own. A `tcmin_surv`-relaxed arm admits ten types
still carrying Earth's phenology, establishment and allocation parameters, so a
null result from it cannot distinguish "the cap has no headroom" from "these ten
are wrong here for other reasons", and it cannot return the answer it would be
bought for.

## What this rules out

- **Not an establishment gate on the survivors.** C3G's bioclimatic limits are
  all absent: `tcmin_surv`, `tcmin_est`, `tcmax_est` and `twmin_est` are +/-1000
  and `gdd5min_est` is 0. Grass is permitted to establish in every cell on the
  planet, and does.
- **Not respiration.** Carbon use efficiency at the poles is 0.61 against 0.51
  in the tropics, so a LARGER share of what is fixed is retained there.
- **Not the photosynthetic rate.** Gross primary production per unit leaf area
  is 0.245 at the poles against 0.224 in the tropics -- per leaf, polar
  photosynthesis runs slightly faster. There is simply almost no leaf: LAI
  0.067 against 4.36.
- **Not the growing-season warmth.** Polar GDD5 is 844.4 against BNS's
  threshold of 175.356.
- **Not the light season.** The polar summer receives more light than the
  tropics ever do. At 32 degrees obliquity the summer pole stands under
  1306.56 x sin(32 deg) = 692 W/m2 at the top of the atmosphere continuously,
  against 541 W/m2 for the same geometry on Earth, and both caps are land.

## Two reinforcing loops, both properties of the obliquity

1. **The winter cold trap throttles supply.** For eight to nine months the cap
   sits below freezing and reaches -68.75 degC. Air that cold carries almost no
   water vapour, so only 77 mm arrives in a year, and over half of it arrives
   in the form that runs off.
2. **Dryness makes the summer hot, and heat costs more water.** A dry surface
   spends almost nothing on latent heat, so it reaches 31.11 degC under
   390.67 W/m2; the same surface under a closed canopy over wet soil would run
   far cooler. The heat then raises the evaporative demand on the 31 mm that
   falls while above freezing.

## What this does not settle

The climate model's polar hydrology is the weak link and this note does not
test it. 77 mm per simulation year is a small number produced at T21, where 128
columns run round a latitude circle whatever its size, and the moisture
transport into the cap is the term that sets it. A climatology at a higher rung
could move the total without anything here being wrong. What is robust to the
total is the cold filter -- ten types at exactly zero against a -68.75 degC
coldest month -- and the partition, 4.7% of arriving water reaching plants.
