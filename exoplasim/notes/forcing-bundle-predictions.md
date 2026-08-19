# Predictions for the next baseline's forcing bundle

**Context first, because this file's vocabulary invites misreading: this is
WORLDBUILDING.** Vesper is a fictional planet. Everything below is planetary
climate modelling for that fiction: offline estimates of what a handful of
namelist switches will do to a toy GCM's simulation of an invented world,
written down before the run. Nothing here refers to the real world, and
nothing here is biology. The nearest real-world analogue of this document is a
test plan.

Written 2026-08-19, measured on the baseline climatology of `run_8c2e1ff9ab5e`,
per WORKFLOW A3: every forcing change lands with a quantitative prediction of
its own effect, stated before it is run, with what result would mean "wrong".
The bundle is then checked as a SUM against the sum of these predictions, and
bisected only if the two disagree. Numbers here are recomputable:
`exoplasim/scripts/predict_ocean_terms.py` for the two ocean terms, and the
sources cited per row for the rest.

## CLIM-16: ocean horizontal heat diffusion, `nhdiff` and `hdiffk`

`oceanmod.f90` reads both keys in `oceanmod_nl`, broadcasts them, and acts on
`nhdiff > 0` at line 844: a diffusion of the SST anomaly from the ocean mean,
in flux form, with no flux across coastlines, sub-stepped `nsub` times. Two
properties of the operator do most of the predicting.

**It conserves exactly.** The global integral of the applied heating is zero by
construction, machine-zero in the offline reproduction. So this term is a
REDISTRIBUTION, not a forcing: the prediction for the global-mean flux
difference between A/B arms is 0.00 W/m2, and any appreciable global-mean
signal is an implementation flag, not a result.

**It is linear in the field and in `hdiffk`.** So one number per bracket point
scales the rest.

The offline reproduction was validated before use, on a check with a right
answer: applied to P2(sin phi), an eigenfunction of the spherical Laplacian,
it must return -6/a^2 times the field, and does to 0.92%. The check earned its
place immediately: the climatology's latitude axis runs north to south, which
sign-flips the meridional term if unnoticed, and the first draft had exactly
that error. The numbers below postdate the fix.

Applied to the baseline annual-mean SST (ocean zonal-mean contrast 41.0 K):

| hdiffk, m2/s | rms heating, W/m2 of ocean | largest single cell |
| --- | ---: | ---: |
| 300 | 0.75 | 22.5 |
| 1000 (the default) | 2.51 | 75.1 |
| 3000 | 7.53 | 225.3 |

The pattern at the default, by latitude band: -0.16 W/m2 at 0-20 degrees,
-0.27 at 20-40, +0.37 at 40-60, +0.36 at 60-90. Low latitudes lose heat to
high latitudes, which is the transport the term exists to stand in for. The
extreme cells sit at coastlines and the ice edge, where the anomaly gradient
is sharpest.

**Predicted global-mean temperature response: +0.03 to +0.12 K at the default
`hdiffk`.** The conserving term moves the mean only through the asymmetry it
creates at the ice edge: warming at 40-90 degrees removes marginal sea ice,
which darkens the planet slightly. The bracket comes from the local damping
assumption (1 to 3 W/m2/K) and the albedo-per-ice-area bracket below. Small,
and one-signed toward warming.

**What would mean wrong, in the A/B:**

- The arms' global-mean net TOA fluxes differ by more than 0.3 W/m2. The
  operator conserves, so beyond the small ice term there is nothing physical
  for a global signal to be.
- The band redistribution comes out with the opposite sign: the 0-40 degree
  surface must lose heat relative to the control and 40-90 must gain.
- The response departs grossly from linear across the 300-3000 bracket.

**A direct pointwise check now exists.** The ocean stream survives every model
call as of CLIM-12, and with `nhdiff > 0` its horizontal-diffusion code becomes
nonzero: the run itself records the flux it applied, which can be compared cell
by cell against this operator applied to that run's own SST. That is the
strongest available test and it can fail.

**A trap, flagged in advance:** `close_ocean_energy.py` asserts codes 903-906
are identically zero and refuses otherwise. On an `nhdiff > 0` run that
assertion breaks BY DESIGN. The script is right to refuse until it is taught
the new term; do not read the refusal as the run being broken.

Per the error budget's entry: bracket `hdiffk` and report the spread. A
constant diffusivity is a bound on the missing transport, not the transport,
and the argument for the term is that the transport exists, not that any value
of it improves an agreement.

## CLIM-17: the freezing point from the declared salinity

The key landed with CLIM-17's closure: `config/planet.yaml` declares
`ocean.salinity_psu: 34.7` and `run_exoplasim.py` derives `TFREEZE` from it,
271.2449 K against the compiled 271.25.

**At the declared salinity the prediction is NOTHING.** The 5.1 mK threshold
shift is three orders below the 0.23 K run-to-run spread. This is the term's
contribution to the bundle sum: 0.000 K by construction. If the next baseline
differs measurably from expectation and this key is the suspect, the key is
not the cause and something else moved.

**The bracket arm is where the prediction has content.** A lower salinity
raises the freezing point, and ocean whose coldest output bin sits between the
two thresholds becomes newly able to freeze seasonally. Counted on the baseline
climatology:

| salinity, psu | TFREEZE, K | newly freezing, fraction of planet | absorbed shortwave | global mean |
| --- | ---: | ---: | ---: | ---: |
| 30 | 271.5121 | 0.00062 | -0.05 to -0.08 W/m2 | -0.04 to -0.07 K |
| 25 | 271.7916 | 0.00190 | -0.15 to -0.24 W/m2 | -0.13 to -0.22 K |
| 20 | 272.0668 | 0.00322 | -0.26 to -0.41 W/m2 | -0.23 to -0.36 K |

The current seasonal-freezing ocean is 0.0489 of the planet and the annual-mean
ice fraction 0.016, so the fresh end adds about a fifth to the existing ice,
not a new regime. Conversion assumptions, both carried as the ranges above
rather than resolved: planetary albedo gain of 0.25 to 0.40 per unit of new
ice area (a surface contrast near 0.5, attenuated by the budget's 0.5 factor
at the low end, weakly attenuated at the high end), and kelvin via
`lib/sensitivity.py`.

Two stated biases. The coldest OUTPUT BIN averages about 480 timesteps, so
instantaneous minima are colder and the areas above are floors. And the static
count carries no ice-albedo feedback, which at the ice edge amplifies; the
high end of the range is not a hard ceiling.

**What would mean wrong, in the A/B:**

- Order of checks: AREA first, flux second. The S=20 arm should grow seasonal
  ice cover by roughly +0.3% of the planet; if no new ice forms at all, the
  bin-mean proxy failed and the flux prediction is void rather than confirmed.
- The sign: a raised freezing point must cool the simulated mean. Warming is
  wrong outright.
- Magnitude: the S=30 arm moving more than 0.2 K, or the S=20 arm beyond
  about -1 K, is outside anything the static estimate plus a generous edge
  feedback allows, and points at the implementation or the comparison.

## The bundle, summed

A3's rule: check the bundle's total against the sum of the per-term
predictions; bisect only on disagreement, and only into the terms whose
predictions were soft.

| term | change | prediction, global mean | source |
| --- | --- | ---: | --- |
| PHYS-9 | `h2oswl` 1.127 | +1.25 K (+1.4 W/m2 TOA, +7.4 W/m2 atmospheric) | scaled from the table in `exoplasim/notes/corrk-cross-check.md`; arithmetic below |
| PHYS-10 | line-list CO2 coefficients | about -0.03 K | the same table's CO2 row prices -7.9%; the landed fit is -7.2% at the planet path |
| CLIM-16 | `nhdiff = 1`, `hdiffk = 1000` | 0.00 W/m2 global by construction; +0.03 to +0.12 K via the ice edge | this note |
| CLIM-17 | `TFREEZE` from declared salinity | 0.000 K | this note |
| SPEC-1 | model reads `k25v` | 0.00 to +0.03 K | `analysis/error_budget.json`, measured on the warm state; the A/B re-measures it |
| DUST-11 | prescribed dust arm | separate note | `aeolian/notes/prescribed-dust-run.md` |

**Sum, excluding dust: +1.3 K, give or take a quarter, dominated by one term.**

The PHYS-9 arithmetic, shown because a units slip here is exactly the kind of
error this project documents: the term is linear in the PRODUCT
`h2osww * h2oswl`. The product moves from 1.3456 to 1.3456 x 1.127 = 1.5165,
a shift of +0.1709 in the same units as the corrk note's calibration row,
which prices a shift of 0.0185 at 0.80 W/m2 atmospheric, 0.15 W/m2 TOA and
0.135 K. Scaling by 9.24: +7.4 W/m2 atmospheric absorption, +1.4 W/m2 TOA,
+1.25 K. The caveat is the scaling itself: the row's sensitivities come from
`exoplasim/notes/shortwave-water-vapour.md` and are being stretched ninefold,
so this is the softest prediction in the table and the first place to bisect
if the sum misses.

**Consequence to plan for:** +1.3 K leaves the design mean, so the flux
re-derivation that WORKFLOW 6C already requires should expect to move
luminosity DOWN by about 0.6% (1.25 K over the canonical 202 K per unit flux).
That is inside the 2-3% window where the `k25v` spectrum remains valid, so no
spectrum rebuild follows. And per WORKFLOW A3, hold the flux for the A/B
itself: measure the surface first, move the flux after, on a slope measured
with the new terms in place.
