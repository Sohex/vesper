# What two flux points say about the comfort band and its cold-extreme cap

*Measured 2026-08-27 on `canonical-10m-base`, T21 at dt 45. Two converged runs:
`run_432e5e46adef` at flux ratio 0.945 (the bootstrap, purpose spinup) and
`run_b45380e61f90` at 1.000 (the bracket, purpose diagnostic). Both pass all six
convergence criteria. This is a WORLDBUILDING project and every quantity below is
a diagnostic of a simulated planet.*

## The amplification is large, and a uniform shift cannot stand in for it

The two points differ by 0.055 in flux ratio and 8.544 K in global mean surface
air temperature. The land seasonal extremes do not move with the mean:

| band | coldest-bin dT | warmest-bin dT | coldest-bin amplification |
| --- | ---: | ---: | ---: |
| 60 to 70 S | +25.04 | -1.41 | 2.93 |
| 50 to 60 S | +20.28 | +3.80 | 2.37 |
| 50 to 60 N | +19.76 | +2.83 | 2.31 |
| 60 to 70 N | +19.60 | -2.73 | 2.29 |
| 10 S to 10 N | +6.72 to +7.23 | +6.31 to +7.05 | 0.79 to 0.85 |

High-latitude WINTER warms about three times the global mean while the same
bands' SUMMER cools slightly. A uniform shift of the land field understates the
cold tail's response by that factor, which is why `derive_design_flux.py`
projects per band from two runs rather than sliding a global sensitivity, and
why a single-point estimate of this cannot be made at any length of run.

## The cap and the comfort score pull in opposite directions

`derive_design_flux.py` scores the land-area fraction inside its comfort band
and constrains the share below the cold extreme. Both criteria, measured:

| flux | comfort score | land below -40 C | land above +38 C |
| ---: | ---: | ---: | ---: |
| 0.945 | 0.657 | 0.262 | 0.045 |
| 1.000 | 0.485 | 0.138 | 0.200 |

Warming halves the cold-extreme share and quadruples the warm-extreme share.
The comfort score falls monotonically across the measured interval, so the
unconstrained comfort optimum lies at or below 0.945.

## The trade, projected per cell between the two points

Linear per cell between the measured pair. Everything outside 0.945 to 1.000 is
EXTRAPOLATION and is marked; the doctrine this project derives fluxes under is
to bracket between two points.

| cap on land below -40 C | best comfort | flux it needs | land above +38 C there | inside the measured pair |
| ---: | ---: | ---: | ---: | --- |
| 0.025 | 0.274 | 1.100 | 0.545 | no |
| 0.050 | 0.376 | 1.060 | 0.433 | no |
| 0.100 | 0.450 | 1.020 | 0.290 | no |
| 0.150 | 0.485 | 1.000 | 0.200 | yes |
| 0.200 | 0.550 | 0.980 | 0.108 | yes |
| 0.250 | 0.609 | 0.960 | 0.059 | yes |

The declared cap of 0.05 is REACHABLE and costs 43 per cent of land above the
warm extreme for a comfort score of 0.376. That is not a habitability gain: it
exchanges a quarter of land with a brutal winter month for nearly half with a
brutal summer one.

**Three limits on the projection, each of which moves the answer.** The
interesting region lies above the measured pair, so a defensible answer near
1.06 needs a third run there rather than a projection. The per-cell projection
is linear and the cold tail's response is not, because sea-ice retreat is a
threshold process and the measured amplification is itself a function of how
much ice remains. And both runs are on BOOTSTRAP surface fields, with uniform
soil water and no lakes, which are the terms that moderate seasonal extremes.

## The criterion's instrument, not its value, is what this questions

The cap counts a cell whose COLDEST BIN mean falls below -40 C: one month in
twelve. A cell at +15 C in its warmest month and -45 C in its coldest scores
identically to one that is never habitable. At this world's 32 degrees of
obliquity the land seasonal range spans 123 K, from a warmest-bin maximum of
+47.2 C to a coldest-bin minimum of -75.8 C, so a severe but brief winter
describes a large share of otherwise ordinary land. Whether a coldest monthly
mean is the statistic that should gate the design is a separate question from
where its threshold sits, and it is the one this measurement raises.

## An unrelated finding: the flux-to-kelvin slope

The two points imply 155.4 K per unit flux ratio, on a ten-orbit tail mean of
each. `lib/sensitivity.py` declared about 202 when this was measured, from a
bracket taken before the Stephens cloud tables, the derived orographic roughness
and the 1.699x shorter hyperdiffusion landed -- and on two runs that have since
been deleted, on a build that is no longer active.

**That flag was taken up the same day and the module now declares 159.7**, from
these same two runs through its own instrument, the fitted asymptotes, rather
than through the tail means above. The 4 K between the two estimators is the
0.945 run's last twelve orbits sitting above its own asymptote, and it is inside
the spread the module declares. The measurement in this note is the corroboration
and not the declaration.

## The derivation does not run, and 0.945 stands by decision

`derive_design_flux.py` refuses on these two points: "a band's warm-season
amplification is not positive; the projection cannot be trusted". Four bands
have a negative warmest-bin response -- 70 to 80 S, 60 to 70 S, 60 to 70 N and
70 to 80 N -- so the per-band projection has nothing to scale by. The guard is
right and the refusal is the correct behaviour. A third flux point does not fix
it: the sign is a property of how this world's seasons respond and not of where
the two points sit.

So the flux was CHOSEN and not derived. `orbit.baseline_flux_earth` stays at
0.945, and the world it gives was judged acceptable on the two points measured:
a comfort score of 0.657, the best available anywhere in the projected range, a
warmest-bin land mean of 23.0 C and a coldest-bin land mean of -13.8 C. The
cold-extreme share of 0.262 is a consequence of 32 degrees of obliquity rather
than of the flux, and reaching the declared cap costs more than it buys.

A COINCIDENCE WORTH RECORDING AS ONE. The 0.945 anchor's original justification
was reverse-engineered from where the superseded physics put the tropics, and
WORLD-F9IG destroyed that justification outright. The number nevertheless sits
on the comfort plateau under the corrected optics. That is two independent
things arriving at one value, and it should not be read as the first
justification having been right after all.

The plateau is flat and the measurement cannot resolve movement across it: the
projected comfort score reads 0.6611 at 0.920, 0.6605 at 0.930, 0.6600 at 0.940
and 0.6571 at 0.945, while the score's own standard error over the settled
window is 0.0024. A move from 0.945 to 0.930 buys 0.0034, a ratio to the
standard error of 1.41 against the 2.8 a paired comparison needs. There is no
optimum to find inside that range on this evidence, and everything below 0.945
is extrapolated past the measured pair in any case.

## The soil built on this climate, and one thing the climate is missing

*Measured 2026-08-27 on `canonical-10m-base` against `precarve-craton-10m`. The
two builds differ in terrain AND in pedogenesis values, so this is not a
single-variable comparison and is not read as one.*

| | precarve-craton-10m | canonical-10m-base |
| --- | ---: | ---: |
| soil depth, median | 0.876 m | 0.203 m |
| AWC, median | 112.80 mm | 25.09 mm |
| AWC, p90 | 306.98 mm | 174.50 mm |

The land is BIMODAL rather than uniformly thin. The top quartile carries 1.454 m
of regolith and 151.66 mm of plant-available water, which is ordinary ground; the
bottom quartile carries 0.050 m and 10.07 mm, which is bare rock, and 28.9 per
cent of land sits at or below 0.10 m.

**Most of the move is a sourced change and stands.** WORLD-9CTM replaced
`dry_erosion_baseline`, an unsourced 0.15, with 0.85 derived from Portenga and
Bierman's arid-to-global denudation ratio of 0.459: the moisture term is
`runoff/reference + b`, so at zero runoff the term IS `b`, and `b/(1+b)` returns
their ratio exactly. Arid erosion therefore rose by 5.67x and depth goes as its
reciprocal. The observed median thinning of 4.3x is what that implies for a cell
carrying some runoff, and `maximum_depth_m` rising 5.0 to 6.70 pushes the other
way by 34 per cent and loses. Physics is not a knob: the term is better sourced
than what it replaced and a harsher world is the information.

**What is NOT settled is the climate the arid tail is evaluated against.**
`build_soil.py` takes its runoff from the bootstrap climatology, and a bootstrap
is terrain-only by definition, so it carries no lakes on any iteration. The
ceiling on what that omits is `area_at_spill_fraction_of_planet` at 0.0821: at
most 8.2 per cent of the planet as inland open water, sitting under the 74.4 per
cent of land that drains internally, which is where the thin soil is. The actual
figure is `surface_water.py`'s to compute.

The reason this reads as an oversight rather than a separation is inside the
graph: `soil` needs `lpj_run`, which needs `lpj_driver`, which needs
`baseline_climatology`. So one step already takes its VEGETATION from a
lake-bearing climate while taking its RUNOFF from a lake-free one.
