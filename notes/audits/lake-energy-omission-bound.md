# Bounding the seasonal lake-energy omission, before an implementation is chosen

**Measured:** 2026-08-25, and the freezing threshold closed 2026-08-31. The
three terms are closed form throughout, from the model's own constants and this
planet's orbital period; the threshold that decides which of them dominates is
read off an existing baseline climatology and an existing solved lake set.
Nothing was simulated for either: no ExoPlaSim run, no hydrography solve.
HYD-21 asks for the bound before the choice, and this is the bound.

Worldbuilding. Vesper is an invented super-Earth. Everything below is about the
simulation of its land surface: `landmod.f90`'s soil heat solver standing in
for a lake, and what that surrogate does not carry.

## What the surrogate is

There is no lake model in PlaSim. A cell holding a lake gets `landmod.f90`'s
soil column, whose heat capacity and thermal conductivity are interpolated
between a dry and a saturated endpoint on the water that column carries. A lake
cell's store sits at or near its capacity, so the surrogate a lake gets is the
SATURATED end: a thermal inertia of about 2063 J/m2/K/s^0.5 on the build's own
median column, which `notes/audits/soil-thermal-inertia.md` measures and which
is a wet, dense soil.

The bound below was computed when the column carried one pair for the whole
planet, at an inertia of 2078. That number is within one per cent of the
saturated endpoint, so the comparison it supports is unchanged; what has changed
is that a DRY cell no longer gets a lake's inertia, which sharpens the surrogate
rather than moving this bound.

The column is 12.4 m deep over five layers. This planet's orbit is 182.80 Earth
days, so the seasonal skin depth `sqrt(2 k / (rho c omega))` is 1.94 m and the
column is 6.4 skin depths deep. It is therefore semi-infinite for the seasonal
cycle, which is what makes a surface admittance the right way to compare it
with a water column, and is worth stating because a shallower column would need
a different comparison.

For a periodic forcing at the orbital frequency the semi-infinite column
presents a surface admittance `I * sqrt(omega)` at 45 degrees, which is 1.311
W/m2/K and an effective heat capacity of 3.295e6 J/m2/K. A well-mixed water
column of depth `d` presents `i * omega * rho_w c_w d` at 90 degrees.

**The surrogate is worth 0.79 m of water.** That is the whole of what the model
gives a lake, at this planet's orbital frequency, and it is the number every
term below is measured against.

## Term 1: freezing. Absent entirely, and it dominates

`hydrography/config/land_water_ledger.yaml:soil_ice_in_the_climate_column`
declares that the five soil temperature layers carry no water and no phase, so
the surrogate has no latent heat of freeze or thaw at all.

Freezing one metre of lake ice releases `rho_i * L_f = 3.06e8` J/m2. Against
the surrogate's 3.295e6 J/m2/K, that is the heat the modelled soil column
exchanges over a seasonal temperature swing of 93 K. Ten centimetres of ice is
worth 9.3 K of it.

**So the term is a threshold, and both halves of the threshold are now
measured.** Whether a lake freezes is a question about the climate state, and
until a baseline climatology existed it could only be posed. It is closed
below. What must not be said, and was, is that 93 K exceeds any seasonal range
this planet can have: it does not.

### The measured freezing share, 2026-08-31

`hydrography/scripts/lake_mosaic_cost.py` on `canonical-10m-carve2` at T21,
from `baseline_regular_climatology.nc` (`run_67323a923013`, 57 orbits) and the
solved lake set under `support_exoplasim-T21.nc`, matched by index and never by
longitude.

A climatology bin is a fifteenth of this planet's orbit, so a bin mean
understates a seasonal minimum and the share is a BRACKET: the warm end is the
minimum over bin means of the surface temperature and the cold end the minimum
over the within-bin minima.

| population | freezes, warm end | cold end |
| --- | ---: | ---: |
| solved lake area | 55.8% | 63.7% |
| all land area | 50.9% | 56.4% |

**The freezing term is ON over most of this world's lake area at both ends of
the bracket.** The two populations agree to within six points, which is what
says the verdict does not rest on the lake set: that set is solved under the
bootstrap forcing while these temperatures are the baseline's, and a result
that holds over all land area regardless of where lakes sit survives that
mismatch.

### The equivalent ice thickness, which replaces the 93 K figure

The dominance is conditional and the single figure hid the condition. 12.1% of
solved lake area sits in cells whose seasonal surface temperature range exceeds
93 K, reaching 127 K; those cells are polar, low-lying and carry a median
annual maximum snow depth of 0.03 m, so the range is a high-latitude seasonal
cycle on a 182.80-day orbit and not a snow-surface artifact.

Inverting the comparison per cell gives the depth of ice whose latent heat
equals that cell's OWN seasonal sensible exchange, which is the number a lake
model has to clear:

| over | p5 | p50 | p95 | max |
| --- | ---: | ---: | ---: | ---: |
| all solved lake area | 0.10 m | 0.33 m | 1.23 m | 1.37 m |
| the lake area that freezes | 0.28 m | 0.58 m | 1.27 m | 1.37 m |

**Above its own equivalent thickness the omitted phase term dominates the
sensible term the depth classes are about; below it, it does not.** Over the
lake area that freezes the crossover is 0.58 m at the median, so a modelled
lake growing more than about half a metre of ice is one whose energy omission
is dominated by phase rather than by depth. The ranking in this audit survives,
and it now rests on a measured crossover rather than on a single figure that
was wrong at the tail.

## Term 2: the water column's sensible heat, bounded two ways

The seasonal surface temperature amplitude of a surface with admittance `Y` and
atmospheric damping `lambda` goes as `1 / |lambda + Y|`, so the ratio of a
lake's amplitude to the surrogate's is `|lambda + Y_soil| / |lambda + Y_lake|`.
It tends to 1 as the damping grows and to `0.79 / d` as it vanishes, so a LOWER
bound on the damping gives an UPPER bound on the omission.

The lowest defensible damping is the radiative part alone, `4 eps sigma T^3`
with the configured `land_longwave_emissivity`, which is 4.33 W/m2/K at 273.15 K
and 6.52 at 313.15 K over the liquid-water span. Turbulent and latent coupling
only add to it, so the first column below is the bound and the others show how
fast it weakens.

| lake depth | amplitude ratio at 4.33 | at 6.52 | at 15 |
| ---: | ---: | ---: | ---: |
| 0.5 m | 1.21 | 1.14 | 1.06 |
| 1 m | 1.15 | 1.12 | 1.06 |
| 2 m | 0.98 | 1.03 | 1.04 |
| 5 m | 0.57 | 0.71 | 0.93 |
| 10 m | 0.31 | 0.42 | 0.71 |
| 30 m | 0.11 | 0.15 | 0.31 |
| 100 m | 0.03 | 0.05 | 0.10 |

Ratios above one are the surrogate over-damping a lake shallower than itself.
**A lake of about 2 m is represented exactly**, and below 2 m the sign of the
error reverses.

The phase omission has a hard ceiling that needs no damping bracket at all. A
water column cannot lag the forcing by more than a quarter orbit, which here is
45.70 Earth days, and the surrogate already lags by 5.1 days at the radiative
damping. **The seasonal phase of a lake surface can be wrong by at most 40.6
Earth days**, whatever its depth, whatever its salinity, and whatever the
climate. It reaches 33 days at 10 m and saturates by about 30 m, so a depth
class finer than "shallower or deeper than about ten metres" buys nothing on
this axis.

## Term 3: the freshwater-to-brine property bracket is not a driver

HYD-21 asks for the depth and phase physics to be crossed with a property
bracket from fresh water to a closed-basin brine. It does not need crossing.
Both terms above enter through the product `rho c d`, so a property change and
a depth change are the same axis. Salinity moves `rho c` by well under a factor
1.2 between fresh water and a saturated brine, and the depth axis spans a
factor of 200 over the classes above. The property bracket is worth less than
one step of the depth bracket and cannot reorder anything.

Where salinity DOES matter is in term 1, through the freezing point, because it
decides whether the dominant term switches on at all. That is a threshold, not
a property bracket, and it belongs with the ice question.

**And that threshold has real leverage, measured 2026-08-31.** Of the solved
lake area, 55.8% has a seasonal minimum below the fresh freezing point but only
31.9% is below the NaCl eutectic, 21.1 K colder. The 23.9% in between is lake
area whose seasonal minimum lies inside the span a chloride brine's freezing
point can occupy, so whether the dominant term switches on there is decided by
the lake's salinity and by nothing else. That share is comparable to the share
already below the eutectic, not small beside it.

The eutectic is used here as a BRACKET WIDTH and never as this world's freezing
point: the lake chemistry is not declared anywhere, so what is bounded is how
much lake area the chemistry could move, which is a quarter of it. A closed-
basin brine chemistry is therefore worth deriving for the ice question, and is
still worth nothing for the property bracket in this section.

## What all three sit inside

`land_water_ledger.yaml:open_water_evaporation_from_routed_water` already
declares the larger omission: routed runoff accumulates in `driver` and
discharges to the ocean, nothing evaporates from it, and a lake cell's annual
evaporation is capped by its own precipitation however the bucket depth is set.
`open_water_evaporation` carries `flux: undeclared`, so the model supplies none
of it.

A closed basin at equilibrium evaporates its catchment's whole delivery from
its lake surface, so the term the model is missing is larger than the one it
can produce by the ratio of catchment area to lake area, which for a dry closed
basin is large. Everything measured above is a refinement of the seasonal shape
of the smaller term.

## What follows

- **Do not choose a lake implementation on depth classes.** The ranking is
  freezing, then the water column's sensible heat, then salinity as a property,
  and the first is a threshold question and the third is not a driver.
  `exoplasim/notes/lake-representation.md`'s levers 2 and 4, the per-cell
  `dwmax` and the per-cell mixed-layer depth, address only the middle term.
  Lever 3, the mask flip to slab ocean, is the only one that brings ice, and it
  is available only where a lake takes most of a cell.
- **None of the levers addresses the mass term**, which is the one the ledger
  already owns under LSHY-6 and which is larger. An implementation chosen on
  the energy bound alone would be choosing for the smaller omission.
- The bound above is the number that belongs against
  `open_water_evaporation_from_routed_water` rather than in a new ledger entry:
  the seasonal energy omission is inside the mass omission already declared
  there, not beside it.
- **The area and depth classes can be derived, and they do not settle
  anything.** `hydrography/data/canonical-10m-carve2/` carries `basins.nc` with
  the level/area/volume curves and `surface_water.nc` with a solved lake set,
  so the catalogue side is answerable on the accepted build. What it answers
  with is a population that STRADDLES the only class boundary identified above.
  Measured 2026-08-31 as mean lake depth (volume over area) at the solved
  equilibrium, area-weighted over the 1,074 basins carrying both: 45.2% of lake
  area is shallower than 10 m and 54.8% deeper, with the median at 10.7 m. The
  distribution is split almost exactly ON the boundary, which cannot rank the
  levers -- the same conclusion the three terms above reach from the other
  direction, and a stronger version of it than the pre-carve builds gave, where
  the median sat between 18 and 27 m.

  One caveat still travels with the shares and one has been retired. The
  forcing is `bootstrap_regular_climatology.nc`, a bootstrap and not a
  baseline, so the lake extents are the terrain-only-field answer; a baseline
  climatology now exists and re-solving the lake set on it costs minutes.
  The pre-carve caveat is gone: `canonical-10m-carve2` is carved, so these
  basin numbers are a state and not a limit. What survives either way is the
  straddle, because moving the climate moves lakes along the depth axis
  continuously and a split this close to the middle does not become one-sided
  under a plausible shift.

- **The accepted build still carries zero `INLAND_WATER` regions**, and the
  fact survived the move to `canonical-10m-carve2`: its T21 support file
  reports no inland-water area in any cell. That is structural rather than
  incidental -- `applyInlandWaterLevels` writes the class only for levels a
  caller supplies, and no build in `source/` was given any -- so it does not
  lift by regenerating. It lifts when a water balance feeds levels back into a
  generation. Nothing above depends on it, which is why the bound could be
  taken first: the lake area every measurement here is weighted by is the
  SOLVED lake set from `surface_water.nc`, which exists on the accepted build
  and is a hydrography product rather than a terrain class.

- **What the freezing arm still cannot say is what a lake's ice does to the
  climate, only what the model is failing to carry.** The equivalent thickness
  is the crossover a lake model has to clear, and the thickness a lake actually
  grows is not derivable from a surrogate that has no phase term at all. That
  is a property of the implementation this audit exists to precede, and the
  crossover is the number it will be judged against.
