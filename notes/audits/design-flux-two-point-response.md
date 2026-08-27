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

The two points imply 155.4 K per unit flux ratio. `lib/sensitivity.py` declares
about 202. Both runs are converged, so this is not a spin-up artefact, and the
declared value was measured before the Stephens cloud tables, the derived
orographic roughness and the 1.699x shorter hyperdiffusion landed. The
comparison here is a ten-orbit tail mean rather than that module's own
instrument, so it is a flag for re-measurement and not a replacement value.
