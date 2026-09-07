# The simulated polar caps are cold-filtered first and water-limited second

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
2. **Water then sets what the two survivors achieve**, which is very little.

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

## The two survivors are water-limited, and that is measured separately

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

## How much each term costs

Cells equatorward of 60 degrees with under 120 mm per simulation year -- as dry
as the cap, but with a coldest month of +8.7 degC -- reach a mean cover of
**0.220**, against the polar **0.030**. So at equal water, being polar costs a
factor of about seven.

That factor is NOT attributable to the cold filter alone. It bundles three
things this comparison cannot separate: the loss of ten of twelve types, a
growing season of 3.39 months against twelve, and the melt-pulse timing that
delivers the year's water before there is leaf to use it. Separating them needs
an arm with `tcmin_surv` relaxed, which nothing here has run.

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
