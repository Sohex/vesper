# Predictions for the next baseline's forcing bundle

**Context first, because this file's vocabulary invites misreading: this is
WORLDBUILDING.** Vesper is a fictional planet. Everything below is planetary
climate modelling for that fiction: offline estimates of what a handful of
namelist switches will do to a toy GCM's simulation of an invented world,
written down before the run. Nothing here refers to the real world, and
nothing here is biology. The nearest real-world analogue of this document is a
test plan.

Written 2026-08-19, measured on the baseline climatology of `run_8c2e1ff9ab5e`,
per `docs/src/pipeline/sequencing.md` A3: every forcing change lands with a quantitative prediction of
its own effect, stated before it is run, with what result would mean "wrong".
The bundle is then checked as a SUM against the sum of these predictions, and
bisected only if the two disagree. Numbers here are recomputable:
`exoplasim/scripts/predict_ocean_terms.py` for the two ocean terms, and the
sources cited per row for the rest.

**AMENDED 2026-08-27: EVERY KELVIN IN THIS FILE THAT CAME THROUGH
`lib/sensitivity.py` FELL BY 21 PER CENT, AND THE AMENDMENT IS VISIBLE BECAUSE
THESE ARE REGISTERED PREDICTIONS.** The slope moved from 202.0 to 159.7 K per
unit flux ratio: the bracket 202.0 was measured on is two deleted runs on a
superseded build, and 159.7 is the same measurement re-taken on
`canonical-10m-base`. No W/m2 in this file moved and no arm was re-run. What
changed is one multiplication, applied everywhere the conversion appears, so a
row's flux column and its kelvin column no longer stand in the ratio a reader
who remembers the old figure expects. **No verdict flips**: every A/B verdict in
this file is taken on W/m2 against a W/m2 resolution floor, and the two entries
whose material thresholds ARE in kelvin -- the wet soil albedo ceiling and its
`f` threshold -- carry the amendment in their own entry.
Kelvin measured directly from a run's own settled mean is NOT in this class and
did not move; the arms' own answers stand as they were.

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
| 300 | 0.72 | 21.5 |
| 1000 (the default) | 2.40 | 71.7 |
| 3000 | 7.19 | 215.1 |

The heating is a temperature tendency times the mixed layer's heat capacity, and
that capacity is read from the model rather than written down: sea water's
density and specific heat are `icemod_nl` keys, `icemod` passes them to
`oceanini`, and `lib/sea_water.py` reads the declaration and then the run's own
namelist. `predict_ocean_terms.py` carried them as literals with the specific
heat at fresh water's value, so every figure in this section is the published
one scaled by the exact ratio of the two specific heats. The operator is linear
in the capacity, which is what makes the rescale exact rather than approximate;
a regenerated table needs a declared `baseline_climatology`, and there is
deliberately no fallback.

**These are TRUE diffusivities, and the arms measured below did not run at
them.** `predict_ocean_terms.py` has always built its operator on this planet's
radius, so this table is the prediction at the coefficient named in its own
first column; `hdiffo` did not, until `world-mll`, and divided by a compiled
Earth radius instead. Arms that ran before that fix realised
(PLARAD/6.371E6)^2 = 1.4401 times the coefficient their namelist declared. The
operator is linear in `hdiffk`, so the CLIM-16 bracket, which realised
432/1440/4320, predicts 1.03/3.45/10.35 W/m2 rms rather than the three rows
above. `hdiffk` means what it says from `world-mll` on, and a new arm set to
1000 realises 1000.

The pattern at the default, by latitude band: -0.15 W/m2 at 0-20 degrees,
-0.26 at 20-40, +0.35 at 40-60, +0.34 at 60-90. Low latitudes lose heat to
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
- The response departs grossly from linear across the bracket, which spans a
  factor of ten in `hdiffk` whatever the endpoints are called.

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

## Measured: the CLIM-16 hdiffk bracket, 2026-08-20

Four arms on a bracket plus a length-matched control, all seeded from one
restart on `t42_l10_p8`, one settling orbit then fifteen labelled `diagnostic`.
**The arms realised 432, 1440 and 4320 m2/s**, which is what they are named by
below and what a new arm would have to be SET to in order to reproduce them.
Their namelists declared `HDIFFK` 300, 1000 and 3000, and under the code of the
day that is what realised the bracket above; the run record is findable by
either number and only one of them is a diffusivity. All four carry `cloud_absorption_scale` 1.192, so the set
is internally paired and is NOT comparable with the earlier bundle arms, which
ran at 1.0.

**Read PAIRED, per orbit.** Arm minus control on the same orbit index cancels
the internal variability the two share, which is what A3's one-restart rule
buys; computing each run's spread separately and combining in quadrature
double-counts it and inflates the error bars about sixfold. At five orbits and
unpaired, only the top bracket point resolved at all. Paired at fifteen, every
point does.

| arm | d net TOA, W/m2 | 0-20 deg SST | 20-40 deg | 40-60 deg | 60-90 deg | d sea ice |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 432 | +0.041 +/- 0.159 | -0.070 +/- 0.021 | -0.015 +/- 0.040 | +0.221 +/- 0.045 | +0.283 +/- 0.225 | -0.00052 +/- 0.00030 |
| 1440 | +0.009 +/- 0.166 | -0.037 +/- 0.028 | -0.020 +/- 0.044 | +0.191 +/- 0.048 | -0.104 +/- 0.196 | +0.00001 +/- 0.00025 |
| 4320 | +0.113 +/- 0.121 | -0.128 +/- 0.035 | +0.013 +/- 0.066 | +0.572 +/- 0.065 | +0.519 +/- 0.164 | -0.00124 +/- 0.00026 |

**Conservation holds, which was the primary claim.** Every arm's global net TOA
difference is consistent with zero and inside the declared 0.3 W/m2, so the
operator redistributes and does not force. The term cannot be argued for or
against on a global-mean flux, exactly as predicted.

**The poleward redistribution is confirmed at all three points.** 0-20 degrees
loses and 40-60 gains at every bracket point, each 15 of 15 orbits one-signed.
The 20-40 band sits at zero rather than the loss the offline operator's applied
heating suggested; that is the crossover latitude and the SST response is not
the applied heating, so it is a mild mismatch and not a sign failure.

**The ice-edge term is now measured, through AREA rather than flux.** At 4320
the arm loses 0.00124 of the planet in sea ice, resolved, and on CLIM-17's
conversion that is +0.079 to +0.127 K -- the prediction was +0.03 to +0.12 K at
a true 1000, so at somewhat over four times that diffusivity the term lands at
the low end of a linear expectation. At 1440 itself the ice change is +0.00001,
which is nothing: the predicted kelvin term is not resolved there and is not
claimed.

**Linearity holds, and the apparent failure was a measurement.** In this
realisation the 432 point sat 3.86x above the line through the two upper
points, which looked one-signed at 15 of 15 orbits. A replicate from a
different restart put it at -0.0359, below 1440 rather than above it and
indistinguishable from zero. The excess was an ice-state draw, not a response;
the section above has the numbers. Between 1440 and 4320, where the response
exceeds the variability, the upper over the middle is 2.99 here and 2.89 in the
replicate against a diffusivity ratio of 3.00. The ratio is unchanged by the
relabelling: 4320/1440 is 3.00 exactly as 3000/1000 was, because every arm was
scaled by the same factor.

**The spread the row asked for**, then, widened by the replicate: the 40-60
ocean warming runs +0.13 to +0.19 K at 1440 and +0.36 to +0.57 K at
4320, the two ends of each being two realisations of the SAME key, and the
sea-ice response 0 to -0.00124 of the planet. A constant diffusivity is a bound on the missing transport rather
than the transport, and this is the width of that bound.

## CLIM-19: the 300 anomaly, and the rule that settles it

Declared 2026-08-20, BEFORE the replicate ran.

The bracket's 40-60 degree ocean warming is linear between 1440 and 4320 to
0.3% and the 432 point sits 3.86x above that line, warming slightly MORE than
1440 on 3.3x less diffusivity. Both are resolved 15 of 15 orbits. The arms were
checked first for the cheap explanation and it is not there: the namelists carry
HDIFFK 300/1000/3000 with NHDIFF 1, the control NHDIFF 0, and every other key
including TSWR3 identical across the four. Those namelist numbers are the
declared coefficient and not the realised diffusivity, which is 1.4401 times
each of them for the reason given under the prediction table; the check is on
the four arms being identical but for that key, and it does not depend on which
of the two numbers names the arm.

**The test is a replicate from a DIFFERENT restart.** Internal variability in
the ice state is a property of the initial condition and repeats differently;
a response to `hdiffk` is a property of the term and repeats the same. Same
four arms, same fifteen diagnostic orbits, seeded from `MOST_REST.00035`
instead of `00039`.

**The statistic, fixed here:** paired per-orbit d(40-60 degree ocean SST)
against the arm's own control, and the ratio of the lowest point to the line
through the upper two.

- **Ratio again near 3.9, and 432 again at or above 1440:** the departure is a
  property of the response. Linearity is refuted across the bracket and the
  bound CLIM-16 reports has structure in it.
- **Ratio near 1, or 432 clearly below 1440:** the original was internal
  variability, linearity holds, and CLIM-16's spread stands as a clean bound.
- **Anything between:** report as unresolved. Do not average the two
  realisations into a verdict; two draws do not make a distribution.

**Secondary, and declared so it cannot be fitted afterwards:** the sea-ice
pattern. In the first realisation 432 and 4320 lost ice and 1440 lost none. If
the ice-albedo amplifier is the mechanism, that pattern repeats with the SST
pattern; if the SST pattern repeats while the ice pattern does not, the
mechanism is something else and this note names none.

## Measured: the CLIM-19 replicate, 2026-08-20

The declared test ran: the same four arms, fifteen diagnostic orbits, seeded
from `MOST_REST.00035` instead of `00039`.

| realised hdiffk | realisation 1, REST.00039 | realisation 2, REST.00035 |
| ---: | ---: | ---: |
| 432 | +0.2213 +/- 0.0451, 15/15 | **-0.0359 +/- 0.0763, 7/15** |
| 1440 | +0.1913 +/- 0.0479, 15/15 | +0.1252 +/- 0.0474, 14/15 |
| 4320 | +0.5720 +/- 0.0654, 15/15 | +0.3620 +/- 0.0562, 15/15 |
| ratio of 432 to the line through the upper two | 3.86x | -0.96x |

**Verdict, by the rule declared before the run: INTERNAL VARIABILITY, and
linearity holds.** The 432 point does not merely fail to reproduce its excess,
it changes sign -- one-signed positive on all fifteen orbits in the first
realisation, indistinguishable from zero and slightly negative in the second,
and clearly below 1440 rather than above it. That is the second verdict
verbatim.

**Linearity between the upper two points reproduces, which is the positive
half.** 4320 over 1440 is 2.99 in the first realisation and 2.89 in the second,
against a diffusivity ratio of 3.00. The operator is linear in `hdiffk` by
construction and the response follows it wherever the response is bigger than
the variability.

**The secondary ice check behaves as declared.** At 432 the ice went -0.00052
in the first realisation and +0.00016 in the second, flipping WITH the SST
sign; at 1440 it went +0.00001 then -0.00053. The ice and the temperature move
together and neither repeats, which is the ice-albedo amplifier working on a
variability draw rather than on the diffusivity. The 432 point is unresolvable
at this segment length for a simple reason: the response the line predicts
there is +0.04 to +0.06 K, and the ice draw is worth several times that.

**What this costs CLIM-16, and it is not nothing.** The same key gives +0.191
and +0.125 K at 1440, and +0.572 and +0.362 K at 4320: a 35 to 37% span
between realisations of one configuration. CLIM-16's reported spread was
computed inside a single realisation and is therefore too tight. The honest
bound across both is **+0.13 to +0.19 K at 1440 and +0.36 to +0.57 K at
4320**, and any future arm quoting a single realisation's error bar as the
uncertainty on a mean is understating it by about a third. Neither figure is a
bound at the DEFAULT any more: the default is a true 1000, no arm has run
there, and the nearest measured point is 1.44 times it.

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
| 30 | 271.5121 | 0.00062 | -0.05 to -0.08 W/m2 | -0.03 to -0.06 K |
| 25 | 271.7916 | 0.00190 | -0.15 to -0.24 W/m2 | -0.10 to -0.17 K |
| 20 | 272.0668 | 0.00322 | -0.26 to -0.41 W/m2 | -0.18 to -0.28 K |

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

## PHYS-11: the cloud shortwave constants, bracketed for arms

The cloud optics were stated for the Sun and none of them had ever been
re-weighted for this star. At the time of the arms that meant seven `radmod_nl`
keys: three coefficients of `swr`'s own analytic cloud fits, `acllwr`
(thermal-band mass absorption), and the `rcl1`/`rcl2`/`acl2` triplets by cloud
level. `world-f9ig` has since deleted the three coefficients and put Stephens et
al. (1984)'s tables in their place, which are stated for the Sun in the same
way, so the argument below is unchanged and now applies to the tabulated
co-albedo. The band PARTITION is star-aware through `zsolar1`/`zsolar2`; what is
solar is the physics inside each range.

**The sign is not derivable by inspection, which is the argument for arms
rather than arithmetic.** Two effects oppose. This star puts 0.618 of its flux
in range 2 against the Sun's 0.483, a naive 1.28 on range-2 cloud absorption;
but it concentrates that energy nearer 0.8-1.5 um, where liquid water absorbs
less than in the 2-3 um bands the Sun's range-2 tail reaches. So the honest
bracket on the absorption-like key is 0.78 to 1.28, spanning both signs of the
correction.

**The bracket is on the range-2 co-albedo alone, and so was the measurement.**
`acl2` is read only by the prescribed-cloud branch `nswrcl = 0` selects, and
every arm in this bundle ran `nswrcl = 1`, so the arms that were labelled as
moving both keys were moving one. The measured slope is unaffected -- what the
arms swept was the computed branch's single-scattering co-albedo, which is what
`cloudabs` now scales on the table -- and only the attribution changes.
clim-68, world-f9ig.

**Predicted magnitude bound: 0 +/- 0.5 K.** Clouds are booked at 12.0 W/m2 of
Earth's shortwave absorption beside the water vapour term; +/-28% of that is
+/-3.4 W/m2 atmospheric, and at the same atmospheric-to-TOA ratio the water
vapour calibration row exhibits (0.19), +/-0.65 W/m2 TOA, +/-0.44 K. After
PHYS-9 this is the softest entry in the bundle and the second place to bisect.

`acllwr` was in the task's list because it was untraced Earth tuning, not
because the star moves it: it is a thermal-band constant and has no stellar
dependence, so it gets no arm on this argument. That reason is closed. It is
Kiehl et al. (1998) Eq. 14's `k_l`, the equation around it is that paper's
Eqs. 12-13, and `exoplasim/notes/cloud-water-reference.md` carries the reading,
so it leaves the list. The scattering-side keys are expected near unity and are
held unless the absorption arms surprise.

**What would mean wrong, in the A/B:** an arm at scale 1.0 must be bit-identical
to the control (the keys default to the values the binary already carries); a
response beyond +/-1 K global-mean equivalent is outside what the 12 W/m2
booking allows and points at the implementation or the conversion; and any
response attributed to `acllwr` on stellar grounds is a category error by
construction.

## CLIM-32: ozone's radiative stake, `O3SCALE` arms

Added 2026-08-20. An ARM, not a bundle member: `model.ozone_scale` stays 0.794
in the baseline config and nothing here is enabled. Run at the production
resolution, T42 L10 p16, on the same binary and the same restart as the rest of
the bundle; there is no T21 leg.

`o3scale` multiplies the ozone mixing ratio directly (`radmod.f90:1802`), so it
scales the ozone PATH linearly. The absorptance is Lacis and Hansen, spelled at
`radmod.f90:2409-2411`, and is strongly saturating in path, so halving the
column removes much less than half the absorption:

    A(x) = O3VISW * 0.02118 x / (1 + 0.042 x + 0.000323 x^2)
         + O3UVW  * 1.082   x / (1 + 138.6 x)^0.805
         + O3UVW  * 0.0658  x / (1 + (103.6 x)^3)

with `O3UVW` 0.335 and `O3VISW` 0.914 from the config. The model builds its own
column at `radmod.f90:1782`, `za = a0o3 + a1o3 |sin lat| + aco3 sin lat cos(season)`;
at the compiled 0.25 and 0.11, and area-weighted `<|sin lat|> = 0.5`, the
global-mean vertical column is **0.305 cm-STP** before `o3scale`.

| slant magnification | A at 0.794 | A at 0.500 | flux moved out of the top layer |
| ---: | ---: | ---: | ---: |
| 1.5 | 0.01245 | 0.00939 | 0.99 W/m2 |
| 2.0 | 0.01506 | 0.01116 | 1.26 W/m2 |
| 2.5 | 0.01757 | 0.01285 | 1.52 W/m2 |

So **1.0 to 1.5 W/m2 of shortwave stops being absorbed aloft and continues into
the column below.** Two terms then oppose each other:

- **The top-of-atmosphere loss.** The newly penetrating beam meets the
  reflectance of everything under it, so absorbed shortwave falls by roughly
  that flux times the planetary albedo, 0.30 to 0.46 W/m2, which is **-0.20 to
  -0.31 K** on `lib/sensitivity.py`'s conversion at albedo 0.30.
- **The redistribution gain, which is why this is worth running.** `PTOP` is
  5000 Pa and the profile's `bo3` is 20 km, so the ozone maximum sits AT the
  model top and `dqo3(:,1)` carries everything above it. The absorption removed
  is therefore heating of the topmost layer, which radiates efficiently to space
  and is weakly coupled to the surface. Moving it into the troposphere and the
  surface makes it count for surface temperature where it largely did not. The
  ceiling, if all of it landed and were retained, is **+0.67 to +1.01 K**.

**Prediction: warming, +0.2 to +0.7 K in the global mean.** The sign is assigned
because the ceiling exceeds the loss across the whole magnification bracket; the
magnitude is soft, because how much of the redistributed flux is retained rather
than reflected is exactly what the ten-layer column decides and is not
calculable here.

**What would mean wrong, in the A/B:**

- Order of checks: HEATING first, temperature second. The top layer's shortwave
  heating must fall by roughly 1.0 to 1.5 W/m2. If it does not move, `O3SCALE`
  did not reach the binary and everything below is void rather than confirmed.
- Cooling beyond -0.4 K is outside the pure top-of-atmosphere loss even at the
  longest slant path, and points at the comparison rather than the physics.
- Warming beyond +1.3 K exceeds the all-absorbed ceiling and means the arm moved
  something besides ozone.
- Two assumptions carried as ranges rather than resolved: the insolation-weighted
  slant magnification, bracketed 1.5 to 2.5, and the planetary albedo used for
  the loss term, taken at 0.30.

## CLIM-33: the mixed layer depth, `mldepth` arms

Added 2026-08-20. Also an ARM: `surface.mixed_layer_depth_m` stays 50.0 and the
baseline config does not move. `oceanmod.f90:208` sets `dlayer(NLEV_OCE) =
mldepth`, so the key is a slab HEAT CAPACITY and nothing else.

That makes the leading prediction exact rather than estimated. At 50 m, with
the model's own sea water -- `oceanmod.f90`'s CRHOS 1030 kg/m3 and CPS 4180
J/kg/K, as `close_state_energy.py` reads them -- the slab holds 2.153e8 J/m2/K.
The orbital year is 1.579e7 s, so the annual forcing frequency is 3.978e-7 s-1
and `omega C` is **85.6 W/m2/K**, against a radiative damping of **1.49
W/m2/K** taken from `lib/sensitivity.py`'s canonical slope at a planetary
albedo of 0.30. The ratio is 57: this slab is
deep in the inertia-dominated regime, which the config's own note reaches
qualitatively when it says this world's half-length year damps seasonality about
twice as hard as Earth's ocean does.

Two consequences follow with no free parameters:

| depth | seasonal amplitude, relative to 50 m |
| ---: | ---: |
| 25 m | 2.00x |
| 50 m | 1.00x |
| 100 m | 0.50x |

- **Amplitude goes as 1/depth to within 1.5%**, because the damping term is 55
  times smaller than the inertia term and enters only in quadrature.
- **The phase lag is 89.0 degrees**, essentially a quarter of the orbit, and it
  barely moves with depth for the same reason.

**Prediction for the annual global mean: 0.000 K, by construction.** A heat
capacity cannot move an equilibrium. This is the term's contribution to any sum,
and it is the same statement CLIM-17 makes for its own key.

**The one path to the annual mean is SEA ICE**, and that is the part with
content. A shallower slab swings further in winter, so more ocean crosses
`TFREEZE` seasonally and the planetary albedo rises. CLIM-17's own conversion
applies unchanged: 80 to 129 W/m2 per unit planet fraction of new seasonal ice,
being S/4 times an albedo gain of 0.25 to 0.40. The area itself is not
predicted here, because the climatology it would be counted on does not exist
on this build; the A/B measures it.

**What would mean wrong, in the A/B:**

- Order of checks: AMPLITUDE first. The 25 m arm must roughly double the ocean's
  seasonal range and the 100 m arm roughly halve it. If they do not, `mldepth`
  did not take and the rest is void.
- The annual global mean moving beyond the 0.23 K run-to-run spread WITHOUT a
  matching change in ice area means something other than heat capacity moved.
  That is the falsifying result for the exact prediction above.
- A phase lag that shifts materially between the arms would contradict the
  inertia-dominated regime, and would mean the radiative damping is far larger
  than the canonical slope implies.
- Arms are 25 m and 100 m, a bracket around the declared 50 m rather than a
  search for a better value. Physics is not a knob.

## world-8h6: the cloud water reference, `clwref` arms

An ARM, not a bundle member. `clwref`, the CCM3 reference in-cloud liquid water
density that `mkclouds` gives the modelled cloud water profile, stays at CCM3's
own value in the baseline configuration and nothing here is enabled.

The bracket is a factor of two either side of that value, and the whole
argument is in `exoplasim/notes/cloud-water-reference.md`: what the two CCM3
sources do and do not settle, where the bracket's own source comes from and why
it is not the value's, and where this fork's shortwave cloud optics parts
company with the scheme the value was fixed inside. Three points, run the way
PHYS-11 is run, and the control-valued arm must be bit-identical to the
control.

**The sign IS assignable, which is what separates this term from PHYS-11.** The
modelled longwave cloud emissivity of the low layers is already saturated in
all three arms, so low cloud responds to `clwref` in the shortwave only and its
response is one-signed; the opposing longwave gain is confined to the mid and
high layers, where the modelled water path is one to two orders of magnitude
smaller. The larger arm cools the simulated global mean and the smaller one
warms it.

**Predicted magnitude: both arms beyond 2 K, plausibly 4 to 8 K**, from a
modelled low-cloud band-1 reflectance moving about 0.09 absolute against
PHYS-11's measured 3.2 W/m2 per 2.0 K. The bound is wide on purpose: the
offline estimate that priced PHYS-11 was falsified by a factor of five, in the
direction of the model responding more than the estimate allowed. If a term
this size is confirmed it is the largest lever in this model's radiation, ahead
of every entry in the bundle table below.

**What would mean wrong, and the falsifying conditions, are in the note.** So
is the decision procedure if an arm fails to equilibrate or crosses the sea-ice
transition. `clwref` reaches the model from
`config/planet.yaml`'s `model.cloud_water_reference_kg_m3`, so the arms are a
run rather than a code change.

## The bundle, summed

A3's rule: check the bundle's total against the sum of the per-term
predictions; bisect only on disagreement, and only into the terms whose
predictions were soft.

| term | change | prediction, global mean | source |
| --- | --- | ---: | --- |
| PHYS-9 | `h2oswl` 1.127 | +0.99 K (+1.4 W/m2 TOA, +7.4 W/m2 atmospheric) | scaled from the table in `exoplasim/notes/corrk-cross-check.md`; arithmetic below |
| PHYS-10 | line-list CO2 coefficients | about -0.03 K | the same table's CO2 row prices -7.5%; the landed fit is -6.7% at the planet path |
| CLIM-16 | `nhdiff = 1`, `hdiffk = 1000`, meaning 1000 m2/s since `world-mll` | 0.00 W/m2 global by construction; +0.03 to +0.12 K via the ice edge | this note |
| CLIM-17 | `TFREEZE` from declared salinity | 0.000 K | this note |
| SPEC-1 | model reads `k25v` | 0.00 to +0.02 K | `analysis/error_budget.json`, measured on the warm state; the A/B re-measures it |
| SPEC-5 | `vegetation_albedo` 0.165, bands [0.075, 0.225] | -0.51 K, bracket 0 to -0.71 | `analysis/vegetation_albedo.json`; +0.0105 on composited land mean at 0.48 K per 0.01. One-signed toward a cooler simulated mean |
| PHYS-11 | cloud absorption arms, bracket [0.78, 1.28] | 0 +/- 0.5 K, sign unassigned | this note; arms only, nothing enabled in the baseline config |
| DUST-11 | prescribed dust arm | separate note | `aeolian/notes/prescribed-dust-run.md` |
| CLIM-32 | ozone arms, `O3SCALE` 0.794 -> 0.500 | +0.2 to +0.7 K, warming | this note; arms only, `ozone_scale` unchanged in the baseline config |
| CLIM-33 | mixed layer arms, 25 m and 100 m | 0.000 K annual mean by construction | this note; arms only, `mixed_layer_depth_m` unchanged in the baseline config |
| world-ofn | `clwhsc` derived from `gascon` and `ga` | sign not assigned; see below | `notes/audits/model-earth-centrism.md` finding 5 |
| world-8h6 | cloud water reference arms, bracket a factor of two either side | beyond 2 K each arm, plausibly 5 to 10 K; larger arm cools | `exoplasim/notes/cloud-water-reference.md`; arms only, `clwref` unchanged in the baseline config |

**world-ofn is the one term in this table whose sign is not assigned, and the
reason is that it moves two things in opposite directions.** `mkclouds` carried
the CCM3 cloud-water e-folding length coefficient as an Earth literal while the
mid-layer heights it is measured against are built hypsometrically from the
model's own `gascon` and `ga`; the coefficient is now built the same way, which
shortens it in the same proportion the heights shorten. Integrated over the
T21/L10 sigma set at 25 kg/m2 of precipitable water, the pre-change liquid water
per layer ran 7.2 times the calibrated value at sigma 0.038 falling to 0.78 at
sigma 0.98, so the change thins high cloud sharply and thickens the lowest three
layers slightly, at a column total that barely moves. Thinner high cloud lowers
the longwave trap and thicker low cloud raises the shortwave reflection: both
terms point the simulated mean DOWN, which is why this is the one place the sum
above is likely to be understating a cooling rather than a warming. It is not
priced here because pricing it needs a baseline climatology and
`config/planet.yaml` declares none. Price it before the bundle is summed
against a run, not after.

**Sum, excluding dust and the unassigned cloud term: +0.6 K, spread roughly
-0.3 to +1.3, still dominated by the water vapour level but no longer
one-sided: SPEC-5 is the one term of comparable size and opposite sign.**

*Amended 2026-08-19, so the spread is reconstructable: the rows' stated
brackets sum to +0.35 to +1.37, not -0.3 to +1.3. The difference is a
downside of about -0.65 K carried implicitly by PHYS-9 -- the one row with no
stated bracket, and the prediction the table itself calls softest. The
headline numbers stand as registered; the A/B scores against them.*

The PHYS-9 arithmetic, shown because a units slip here is exactly the kind of
error this project documents: the term is linear in the PRODUCT
`h2osww * h2oswl`. The product moves from 1.3456 to 1.3456 x 1.127 = 1.5165,
a shift of +0.1709 in the same units as the corrk note's calibration row,
which prices a shift of 0.0185 at 0.80 W/m2 atmospheric, 0.15 W/m2 TOA and
0.107 K. Scaling by 9.24: +7.4 W/m2 atmospheric absorption, +1.4 W/m2 TOA,
+0.99 K. The caveat is the scaling itself: the row's sensitivities come from
`exoplasim/notes/shortwave-water-vapour.md` and are being stretched ninefold,
so this is the softest prediction in the table and the first place to bisect
if the sum misses.

**Consequence to plan for, corrected 2026-08-20:** the central sum shifts the
simulated temperature at every candidate flux by about +0.6 K, which is 0.4% of
a flux ratio on the canonical 159.7 K per unit flux. Read that as a shift in the
T(f) the re-derivation is SCORED AGAINST, not as a correction to a standing
luminosity. There is no standing value to correct: `design_flux` has not been
re-derived on this terrain, and CLIM-30 records that the 0.945 on the books is
unsupported because its cold-extreme cap was inferred from its own answer. The
derivation lands where the declared thresholds put it once the bundle is in.
That is inside the 2-3% window where the `k25v` spectrum remains valid, so no
spectrum rebuild follows. And per `docs/src/pipeline/sequencing.md` A3, hold the flux for the A/B
itself: measure the surface first, move the flux after, on a slope measured
with the new terms in place.

---

## Measured: the arm bundle, 2026-08-20

Nine arms, all seeded from `run_4182235e9781/MOST_REST.00039` on one binary
(`t42_l10_p8`, sha 44d895de), each one namelist key off the control, run 2x8
concurrent one job per die. One settling orbit labelled `spinup`, then two
orbits labelled `diagnostic`; only the diagnostic pair is read.

**Carried on every result below: the arms inherit a SUPERSEDED SURFACE.** The
donor predates the current staged fields -- code 229 did not exist when it was
prepared -- so `--superseded-surface-ok` was used deliberately and every run
manifest records `superseded_surface_override` with its reason. Valid for
differences between arms, which is all that is claimed; not valid for any
absolute climate. A recommissioning is expected regardless.

**What this A/B can resolve.** Within-run orbit-to-orbit spread on global ASR
runs 0.03 to 0.65 W/m2, median 0.19, so a two-orbit arm-minus-control
difference resolves at roughly **0.5 W/m2** taking the worst spread. Four of
the seven terms sit at or under that and are reported as unresolved rather than
as small. Kelvin uses `lib/sensitivity.py` at the measured planetary albedo
0.2534, giving 0.629 K per W/m2 of TOA forcing.

| arm | d ASR W/m2 | as K | d sea ice | verdict |
| --- | ---: | ---: | ---: | --- |
| `o3_050` CLIM-32 | -1.312 | -0.82 | -0.0001 | resolved, **prediction falsified** |
| `cld_078` PHYS-11 | -3.304 | -2.08 | +0.0004 | resolved, **prediction falsified** |
| `cld_128` PHYS-11 | +3.178 | +2.00 | -0.0015 | resolved, **prediction falsified** |
| `mld_025` CLIM-33 | -0.428 | -0.27 | +0.0053 | ice path confirmed |
| `sal_30` CLIM-35 | -0.474 | -0.30 | +0.0000 | UNRESOLVED, no ice change to carry it |
| `sal_20` CLIM-35 | -0.105 | -0.07 | +0.0016 | area confirmed, flux unresolved |
| `mld_100` CLIM-33 | -0.177 | -0.11 | -0.0040 | unresolved |
| `hdiff_on` CLIM-16 | -0.024 | -0.02 | -0.0001 | confirmed at zero |

### CLIM-16 confirmed

Predicted 0.00 W/m2 global by construction; measured -0.024, well inside the
floor. The +0.03 to +0.12 K ice-edge term the prediction also carried is below
what two orbits resolve and is neither confirmed nor refuted.

### CLIM-33 confirmed, including its exact half

The annual mean was predicted at 0.000 K by construction, with sea ice as the
only route to it. That is what happened: `mld_100` moved the ocean mean by less
than a millikelvin, and `mld_025` moved it only by growing ice, +0.53% of the
planet. CLIM-17's independently derived conversion, 80 to 129 W/m2 per unit
planet fraction of new ice, puts that at -0.42 to -0.68 W/m2; measured -0.428,
at the low end of a bracket computed for a different task. The mechanism is
confirmed by a number this prediction did not fit.

The amplitude check moved the right way and has not reached its asymptote:
ocean seasonal amplitude 1.51x at 25 m and 0.66x at 100 m, against 2.00x and
0.50x predicted. Two orbits from a state equilibrated to 50 m is a transient,
so the ratios are compressed toward 1. The falsifier was whether `mldepth` took
at all; it plainly did.

### CLIM-35 area confirmed, flux unresolved

The declared order was area first. `sal_20` grew seasonal ice by +0.164% of the
planet against +0.19% predicted, which is the check passing. Its flux, -0.105
W/m2 against -0.15 to -0.24 predicted, is under the resolution floor and is not
a measurement. `sal_30` is the useful negative: -0.474 W/m2 with NO ice change
to carry it, which is how the floor was recognised rather than assumed.

### CLIM-32 falsified, and the mechanism was in the code all along

Predicted +0.2 to +0.9 K of WARMING, with cooling beyond -0.4 K declared as
the falsifying result. Measured **-1.04 K**: wrong sign, and outside the bound.

The error was mine and it is instructive. The prediction argued that absorption
removed from the top layer would be redistributed downward and largely
retained, so the surface would gain more than the top of atmosphere lost. It
ignored the UPWARD path. `radmod.f90:2407` forms the ozone path as
`zxo3t + zmbar*zo3t`: the accumulated downward slant path PLUS a diffuse
upward traverse of the whole column. Ozone sitting above everything absorbs the
reflected beam as well as the incoming one, and that half was never going to
warm the surface -- it is planetary absorption of radiation already on its way
out. Remove the ozone and it simply escapes. The measured ASR loss, -1.31 W/m2,
is close to the whole 1.0 to 1.5 W/m2 the column stops absorbing, which is what
"almost none of it is retained" looks like.

So the arithmetic in the prediction was right about the SIZE of the absorption
and wrong about its FATE. The saturating Lacis-Hansen form and the 0.305 cm-STP
column stand; the redistribution argument does not.

### PHYS-11 falsified on magnitude, and it is much the largest term found

Predicted 0 +/- 0.6 K with +/-1 K declared as falsifying. Measured **-2.63 K at
scale 0.78 and +2.53 K at 1.28**, two and a half times outside the bound. The
two arms are near-antisymmetric about the control, which is the internal check
that the response is linear in the scale rather than an artefact.

The prediction derived its bound from the error budget's booking of clouds at
12.0 W/m2 of Earth's shortwave absorption, taking +/-28% of it through the
water-vapour row's atmospheric-to-TOA ratio of 0.19. Measured, the TOA response
is +/-3.2 W/m2 rather than +/-0.65: a factor of five. Either the 12 W/m2
booking understates cloud shortwave absorption in this model, or that 0.19
ratio does not transfer from water vapour to clouds. The note already called
this the softest entry after PHYS-9 and the second place to bisect; it is now
the first place to look, and it makes the cloud constants the largest unpriced
lever in the radiation, not a correction: +/-2.6 K from a bracket that spans
both signs of the star correction is larger than the entire bundle sum of
+0.6 K, and larger than any single term in it.

#### Which half of the prediction was wrong, measured

The bound was built from two factors: clouds booked at 12.0 W/m2 of shortwave
absorption, and water vapour's atmospheric-to-TOA ratio of 0.19. Splitting the
measurement across those same two factors says which one failed, and it is not
the one that looks soft.

| | control | scale 0.78 | scale 1.28 |
| --- | ---: | ---: | ---: |
| top of atmosphere, `rst` | 240.049 | 236.745 | 243.227 |
| surface, `rss` | 156.624 | 156.681 | 156.194 |
| atmospheric absorption | 83.425 | 80.063 | 87.033 |
| d atmospheric | -- | -3.361 | +3.608 |
| d top of atmosphere | -- | -3.304 | +3.178 |
| **d TOA / d atmospheric** | -- | **0.98** | **0.88** |

**The 12.0 W/m2 booking is right, and the arms measure it independently.** A
22% cut moves atmospheric absorption 3.361 W/m2 and a 28% rise moves it 3.608,
implying the scaled keys carry 15.3 and 12.9 W/m2 of absorption. The budget's
12.0 sits just under that, and the prediction's own +/-3.4 W/m2 atmospheric
estimate was accurate to a few percent.

**The 0.19 ratio is what fails, and it fails by a factor of five.** Measured,
the ratio is 0.88 to 0.98: nearly ALL of the atmospheric change appears at the
top of atmosphere.

The surface says why. It barely moves, +0.057 W/m2 on the low arm and -0.430 on
the high one, against atmospheric changes of 3.4 and 3.6. So the absorption
these keys add is not coming out of the surface's share. It is coming out of
what would otherwise have been SCATTERED BACK TO SPACE.

That is the general lesson and it is not about clouds. **The
atmospheric-to-TOA ratio is not a property of the atmosphere; it is a property
of what the absorber competes with.** Water vapour absorbs low in the column
and mostly takes flux the surface would have absorbed anyway, so the planet's
total absorption barely changes and 0.19 is a real number for that term. Cloud
shortwave absorption sits in a scattering medium and competes with reflection,
so a photon absorbed is a photon that does not leave, and the ratio goes to
nearly 1. Borrowing one term's ratio for the other was the whole error, and it
would repeat on any absorber placed where scattering dominates: dust and sea
salt both qualify.

**Corrected bound.** +/-3.4 W/m2 atmospheric at a ratio near 0.95 gives +/-3.2
W/m2 at the top of atmosphere, +/-2.0 K. That is what was measured. The
prediction was recoverable from its own inputs the moment the ratio was
questioned.

#### The aerosols were checked and are NOT affected

The lesson above puts dust and sea salt on notice, so both were traced rather
than assumed. They are clean, for a better reason than luck.

**Neither ever used a ratio.** The 0.19 appears nowhere under `aeolian/`,
`analysis/` or `scripts/`; both aerosols price their top-of-atmosphere forcing
with `dust_forcing.py:shortwave_forcing`, Chylek and Coakley (1974) two-stream,
and `build_sea_salt.py` calls that same function deliberately so the two are
compared through one formula. It carries the competition EXPLICITLY:

    -insolation * T^2 * ((1 - albedo)^2 * beta * tau_sca - 2 * albedo * tau_abs)

The scattering term and the absorbing term are separate, opposite in sign, and
weighted by the albedo underneath. That is the same physics the cloud arms
measured, written down properly, which is why borrowing a ratio was never
needed here.

**Sea salt is immune outright.** Its single-scattering albedo is 1 to within
1e-5 in both bands, so `tau_abs` is identically zero and the absorbing term
cannot contribute whatever albedo it is given.

**Dust absorbs, and the albedo it competes against is already bracketed.**
`dust_forcing.py` reports the forcing over ocean 0.07, vegetated land, playa
fill and salt crust 0.50 rather than one value, which is the same axis the
cloud result turns on. A dust layer over a cloud deck sees an effective
underlying albedo near or below that 0.50 bright end, so the existing bracket
already spans the cloudy case even though its top row is labelled for salt.

What is left is smaller than the error it was checked for, and is stated rather
than pursued: the clear-sky `TRANSMISSION = 0.79` above the layer is a constant
in that expression, so cloud is represented in the albedo argument and not in
the transmission. That approximation is worth what it is worth; it is not a
factor of five.

#### Derived: the weight is 1.192, and the qualitative argument was backwards

`exoplasim/scripts/cloud_band_weight.py`, step `cloud_band_weight`, 2026-08-20.
Derived the way the gas weights were: from a quantity carrying no spectrum,
liquid water's k(lambda) out of Hale and Querry (1973) Table I, through Mie to
a single-scattering co-albedo, flux-weighted over range 2 against a 5772 K Sun.

| droplet effective radius | weight, linear in co-albedo | weight, square-root |
| ---: | ---: | ---: |
| 5 um | 1.1991 | 1.1854 |
| 10 um | 1.2027 | 1.1854 |
| 15 um | 1.2053 | 1.1844 |

**1.184 to 1.205, median 1.192.** The answer barely depends on droplet size or
on which absorptance scaling is assumed, which is what a well-conditioned
derivation looks like; the bracket is reported anyway because neither was
resolved. Four checks pass before any weight is printed: the solar partition
identity against `radmod.f90:207`'s 0.517, the star's range-2 share at 0.6176
against the 0.618 recorded above, water's index endpoint against the paper, and
the Sun weighted against itself giving exactly 1.

**Two things in the prediction above are now known to be wrong, and they
cancelled partly.** First, the direction: the note argued this star
"concentrates that energy nearer 0.8-1.5 um, where liquid water absorbs less".
It is the other way round. A cooler star's range-2 flux is shifted to LONGER
wavelengths within the band, and water's k climbs steeply there -- 1e-7 near
0.8 um, 1e-3 by 2 um, 0.3 by 3 um -- so the co-albedo goes UP, not down. The
weight is above 1 for the reason the note gave for it being below 1.

Second, the naive 1.28: that was the range-2 flux SHARE ratio, 0.618 over
0.483, and the share is not the model's error. `radmod.f90` already computes
band-2 flux from this star's own `zsolar2`. What is Earth's is the FRACTION of
that flux a cloud absorbs, which depends on the composition WITHIN the band and
on nothing else. Using the share as the weight would have double-counted a
correction the model already makes.

**Adopted at 1.192** in `config/planet.yaml`. On the arms' own measured slope,
12.96 W/m2 per unit scale, that is **+2.49 W/m2 and about +1.6 K** -- larger
than the whole bundle sum, and one-signed. The flux re-derivation must follow
it rather than precede it.

**Consequence.** Two of the bundle's registered predictions are refuted and
neither refutation is a small correction. The summed prediction above is not
re-scored here: these two are ARMS, not bundle members, so the sum they do not
enter is unchanged, but the confidence the sum was quoted with should not
survive a five-fold miss on a term of this size.

## Measured: PHYS-9, 2026-08-19

The first bundle term to be run, and the largest. Predicted **+1.25 K** above.

Two bootstraps on `precarve-craton`, the second seeded from the first's final
restart so the pair is a perturbation off one state rather than two independent
spin-ups:

| run | `h2oswl` | window mean `tas` | note |
| --- | ---: | ---: | --- |
| `run_78c22fb1a1bd` | absent, so 1.0 | 292.013 K | orbits 50-59, converged on all six criteria |
| `run_4182235e9781` | 1.127 | 292.955 K | orbits 30-39, five of six |

**+0.94 K on the window means, +1.20 K on the fitted asymptote** of 293.211 K
against a prediction of +1.25 K.

Three things have to be said with that number rather than after it. The second
run is NOT converged: it misses `extrapolated_offset_lt_0.15_k` at 0.256 K and is
still warming at +0.0138 K/orbit, so +0.94 K is a lower bound and +1.20 K is a
fit rather than a measurement. The pair also carries CLIM-17's `TFREEZE` as well
as `h2oswl`, since both landed together; that term is predicted at 0.000 K by
construction at the declared salinity, so it is not a confound worth
disentangling, but the arm is a pair of keys and not one. And the relaxation was
much slower than a settled run's drift -- 34.8 fitted orbits against 9.9 expected
-- which is what a step change in forcing does and is why the asymptote is doing
work here.

**What it means for the rest of the bundle.** PHYS-9 was called the softest
prediction in the table, its sensitivities stretched ninefold from
`shortwave-water-vapour.md`, and the first place to bisect if the sum missed. It
did not miss. That does not validate the other rows, but it does remove the
term most likely to have carried the error, and the bundle sum of +0.6 K stands
on a firmer largest component than it did when it was written.

It also confirms, separately from any temperature, that the key REACHES the
radiation: `H2OSWL` was absent from every namelist before CLIM-31's neighbour
fix and is present in this run's, and the run responded.

## world-u9hq: the `h2o_sw_level` bracket, `h2oswl` arms

Registered here rather than as a separate pair, because `h2oswl`, `h2osww` and
the declared spectrum all land in the same rebuild and no single run's absolute
temperature can be attributed among them. What CAN be attributed is a difference
between two runs that differ in one key, which is what these arms are.

`config/planet.yaml` carries `h2o_sw_level: 1.163` with
`h2o_sw_level_bracket: [1.129, 1.206]`. The width is the water vapour continuum
in the near-infrared windows, whose size is not known there better than a factor
of a few: `exoplasim/notes/corrk-cross-check.md` derives the central value and
both ends from `shine2012`'s CAVIAR-minus-MT_CKD increment and `mlawer2012`'s
per-window laboratory factors, and states why neither end can be tightened from
anything this project holds. It is a bracket and not an error bar for the reason
that note's BAR 2 gives: the half-width reaches the size at which every consumer
has to run both ends.

**THE ARMS ARE THE TWO ENDS AGAINST EACH OTHER, not each against the central,
and that is a decision about the instrument rather than about the physics.**
Each end is 0.034 and 0.043 from 1.163, worth about -0.34 K and +0.43 K on the
static conversion below. The convergence criteria tolerate a 0.15 K extrapolated
offset, so a single end against the central sits about two to three times the
instrument's own declared slack. End against end is 0.077 of the key, about
0.77 K, on the same two runs. Both are branched off ONE restart so the pair is a
perturbation of a single state, which is A3's condition and is what PHYS-9 did.

**Predicted magnitude, and two estimates that do not quite agree.**

| estimate | per 0.01 of `h2oswl` | end to end, 0.077 |
| --- | ---: | ---: |
| static, `corrk-cross-check.md` | 0.10 K, bracket 0.07 to 0.15 | 0.77 K, bracket 0.54 to 1.16 K |
| PHYS-9 measured, window means | 0.074 K | 0.57 K |
| PHYS-9 measured, fitted asymptote | 0.094 K | 0.72 K |

The static estimate comes from `shortwave-water-vapour.md`'s 43.2 W/m2 of water
vapour shortwave absorption, its 0.193 top-of-atmosphere share and its 0.861 K
per W/m2. The measured pair comes from PHYS-9's own arms, a step from 1.0 to
1.127 worth +0.94 K on the window means and +1.20 K on the fitted asymptote. The
measured window slope is the low end because that run was not converged and
+0.94 K was a lower bound; the asymptote is the better of the two and it sits
just below the static estimate. **The prediction carried into the arms is
0.72 K end to end, the measured asymptote slope, with the static bracket 0.54 to
1.16 K around it.** The static number is not used as the centre because it is an
offline conversion through two ratios and the asymptote is this model's answer to
the same question.

**What would mean wrong, in the A/B:**

- Sign. More absorption at the level the key scales must warm the simulated
  mean. A cooler high arm is wrong outright, not a small result.
- Monotonicity against the central run, where one exists on the same lineage:
  1.129 below it, 1.206 above it. An inversion says the difference is scatter.
- Magnitude. Below 0.54 K or above 1.16 K end to end is outside the static
  bracket and points at the implementation or at the pairing, not at the
  continuum.
- The instrument. Each arm must meet all six convergence criteria on its
  window, and the end-to-end difference must exceed the pooled inter-orbit
  scatter of the two windows. If it does not, the result is "the bracket is
  below what this configuration can resolve" and NOT a number.
- An arm staged at 1.163 must be bit-identical to the baseline it branched from,
  which is the same control PHYS-11's scale-1.0 arm carries.

**What closing the bracket instead would take, and why it is not on this path.**
A correlated-k bundle carrying an MT_CKD or CAVIAR water vapour continuum. The
LMD Generic PCM bundle this project used has `continuum/far_wing_data` for
CO2-CO2, CO2-H2O and CO2-N2 and none for water with itself or with N2, which is
the whole reason the base ratio was one-signed and the bracket exists.

**BLOCKED on a run.** `config/planet.yaml` declares `baseline_climatology: null`
and no run exists to branch a pair from, so the arms are registered and not
measured. `h2oswl` reaches the model as `H2OSWL@radmod_namelist` through
`run_exoplasim.py`'s `SHORTWAVE_GAS_KEYS`, which writes it only when it differs
from the model default, and `verify_staged_namelists` is what makes the arm's
own namelist the record of what it integrated.

## world-trs3: the derived precipitation re-evaporation, `gamma` arms

`rainmod.f90`'s `gamma` was the fraction of a layer's sub-saturation deficit
that falling precipitation evaporates per timestep, a constant 0.01 with no
derivation on either side of the 0.007/0.01 step upstream took at a rung change.
It is now derived per cell and per level from Kessler (1969), and the constant
survives only as an override: `gamma > 0` in `rainmod_nl` restores it exactly.

**The derived form.** Kessler's table 4 rain evaporation for a Marshall-Palmer
distribution, integrated over the layer the four re-evaporation sites act on,
cancels the layer depth and leaves `gamma = 5.4395e-4 * M^0.65 * deltsec2` with
`M = P/V` and `V = 5.17*sqrt(ga/9.80665)*M^(1/8)` m/s. Substituting collapses
the exponents to

    gamma = 5.4395e-4 * (P/zvcoef)^0.577778 * deltsec2

with `P` the layer's precipitation flux in g m-2 s-1. So `gamma` is proportional
to the timestep, goes as `P^0.578`, and carries `ga^(-0.289)` through the drop
terminal speed: 0.926 of its Earth value here.

**What the constant was, in the derived form's units.** At `deltsec2` = 3600 s
and this world's gravity the derived value is 0.036 at 0.5 mm/day, 0.053 at
1 mm/day, 0.101 at 3 mm/day, 0.202 at 10 mm/day and 0.381 at 30 mm/day. The
declared 0.01 is the Kessler value at 0.055 mm/day, some fifty times below a
global-mean rate, and it is below the whole of that span. The change is
therefore a factor of 3.6 to 20 upward over the fluxes that matter, and it is
one-signed: the derived form is never below the constant at any rate the model
produces.

**The mechanism, and why the response is SUBLINEAR in that factor.** The four
sites cap each level's removal at the precipitation actually available,
`AMIN1(gamma*deficit*dsigma/deltsec2, zpr)`. At `gamma` = 0.01 the potential
removal summed over a sub-saturated column is already comparable to a typical
land precipitation flux, so the scheme sits at the transition between
deficit-limited and flux-limited. The derived value puts it firmly in the
flux-limited regime, where the answer is set by how much deficit the circulation
can maintain below cloud rather than by `gamma`. The equilibrium is reached by
MOISTENING: the sub-cloud layer wets until its deficit falls by about the factor
`gamma` rose. That is the observable, and it is what separates this change from
a simple scaling of the evaporated flux.

**Predicted, on a matched T21 L10 pair from one restart, both arms with the
energy fixer on and both segmented for a climatology rather than a diagnostic:**

- **Land P minus E falls by 10 to 50 per cent of the control's land P minus E**,
  centred on -25 per cent. The bracket is wide because it is bounded below by
  the fully deficit-limited limit, where the two arms evaporate the same mass
  and the difference is zero, and above by the no-feedback column limit, where
  at 3 mm/day and a sub-cloud relative humidity of 0.95 the share of falling
  precipitation reaching the ground goes from 0.87 to 0.14. Neither limit is
  the model.
- **Relative humidity in the lowest three levels over land rises by +0.03 to
  +0.15 absolute**, and this is the discriminating measurement rather than the
  flux: it is the mechanism by which the response saturates.
- **The reduction is largest where precipitation is lightest.** `gamma` goes as
  `P^0.578`, so the ratio of derived to constant is largest at low flux; the
  arid, lightly precipitating land where the carve criterion bites is where the
  change concentrates.
- **Runoff falls with land P minus E**, and hydrography reads that field.
- **High cloud fraction rises slightly** with the moistened column.

**What would mean wrong:**

- Sign. Land P minus E must not RISE. The derived `gamma` exceeds the constant
  at every rate the model produces, and more re-evaporation cannot deliver more
  precipitation to the ground. A higher land P minus E in the derived arm is an
  implementation fault, not a small result.
- Sub-cloud relative humidity unchanged to within the inter-orbit scatter. That
  says the code path is not reached at all, and the first thing to check is that
  the arm's own namelist did not carry a positive `gamma`.
- The change concentrated in the heaviest-precipitating cells rather than the
  lightest. That inverts the `P^0.578` dependence and points at the flux
  conversion `P = zpr*dp/ga*1000`, which is `mkrain`'s own output line read
  backwards.
- Magnitude outside 0 to 60 per cent on land P minus E. Below zero is the sign
  failure above; above the no-feedback column limit means the cap is not being
  applied.
- The instrument. Both arms must meet the convergence criteria on their windows,
  and the difference must exceed the pooled inter-orbit scatter. If it does not,
  the result is "below what this configuration resolves" and NOT a number.

**The check that can fail, and it is an identity rather than a comparison.**
An arm declaring `gamma = 0.01` in `rainmod_nl` must be BIT-IDENTICAL to the
control built before this change. The override path restores the literal
constant at all four sites with no other expression touched, so any difference
at all is a defect in the patch and not a result. That check costs one short
segment and settles the implementation independently of the physics.

**A second identity.** The derived `gamma` printed by an instrumented run at a
known layer flux and timestep must match the offline table above to rounding.
The formula has no free parameter, so a mismatch localises immediately to the
unit conversion or to `zvcoef`.

**Kessler's own caveats travel with the form** and bound how well the derived
number can do: the single-drop fit is accurate to about 40 per cent, a constant
`N0` misrepresents evaporation because the process depletes small drops
preferentially, and the rate is for standard air density with no altitude
variation. The two SNOW sites evaluate the rain expression because Kessler
supplies no snow analogue; that is declared at the value and is not a regression,
since the constant carried no flux dependence at any of the four.

**BLOCKED on a run.** `config/planet.yaml` declares `baseline_climatology: null`
and no run exists to branch a pair from, so the arms are registered and not
measured. One rebuild of every binary precedes them, per rule 4.

## world-o12h: the cloud-fraction subgrid width, `rcritwidth` across the ladder

`rainmod.f90`'s `rcrit` is the cell-mean relative humidity at which stratiform
cloud starts, and its 0.85 floor was anchored to T21 with no NLAT term anywhere.
What it encodes is the WIDTH of the subgrid humidity distribution, so the
quantity carrying a resolution dependence is `(1 - rcrit)`. Specific humidity is
a passive scalar, and Kolmogorov-Obukhov-Corrsin gives its variance across a
separation `L` as `L^(2/3)`, so the standard deviation goes as `L^(1/3)` and a
cell's width goes as `1/NLAT`:

    1 - rcrit(jlev) = (1 - rcrit_T21(jlev)) * (32/NLAT)^(1/3)

`rcritwidth` is the key, derived from NLAT unless given a positive value.
`rcritmod` and `rcritslope` keep their meanings and apply on top.

**The anchor is preserved, so T21 does not move.** At NLAT 32 the factor is
exactly 1 and the transform is skipped by construction rather than applied and
rounded. On the 0.85 floor the factor is 1.024 at T31, 1.036 at T42, 1.054 at
T63 and 1.065 at T85, taking `rcrit` to 0.869, 0.881, 0.896 and 0.906 and
`1/(1-rcrit)^2` from 44 to 58, 71, 92 and 112.

**Predicted, on a ladder sweep holding everything else fixed:**

- **T21 is BIT-IDENTICAL to the control.** This is the whole reason the anchor
  was kept at 32 rather than moved to a rung nobody runs.
- **Stratiform cloud fraction falls at every rung above T21**, monotonically in
  NLAT, because cloud starts later in a smaller cell. `rcrit` enters as
  `((rh-rcrit)/(1-rcrit))^2`, so at T42 a 0.031 rise in the floor removes cloud
  wherever the modelled relative humidity sits between 0.85 and 0.881 and
  reduces it everywhere above.
- **Planetary albedo falls and the simulated global mean warms at T42 relative
  to the same configuration before this change.** Cloud fraction is the largest
  single lever on planetary albedo, which is why this row outranks its size.
- **The T21-to-T42 cloud-fraction step SHRINKS relative to the unmodified
  ladder.** That is the point of the change and the thing the sweep measures:
  the scheme is being made less rung-dependent, not more.

**What would mean wrong:**

- T21 differing from the control at all. The factor is exactly 1 there and the
  loop is skipped, so any difference is a defect in the patch.
- Cloud fraction RISING at T42. That inverts the derivation's sign and points at
  the width being applied to `rcrit` instead of to `(1 - rcrit)`.
- `rcrit` exceeding 1 at any level and any rung. The width form is bounded below
  1 by construction; a value at or above 1 divides by zero in the cloud
  expression and means `rcritmod` has been used as the operator instead.
- The T21-to-T42 cloud step widening rather than narrowing. Then the exponent is
  wrong in magnitude even if right in sign, and 1/3 is the thing to question.

**What this does NOT claim.** The 0.85 anchor is still a fit and no paper can
supply it, because a subgrid width is a property of the mesh rather than of the
world. What changes is that the fit is now a declared function of the grid and
can be tested against the grid. The convective pair `zcca` and `zccb` does not
take this scaling and stays irreducible for a different reason, argued at the
value: their exponent is set by an unresolved convective area fraction the model
does not carry, so no exponent between the diluting and cell-filling limits is
derivable from what this scheme holds.

**BLOCKED on a run**, and on the same rebuild.

## OCN-22: the open-ocean albedo's band split, and it is not an arm

Registered here because it is a numeric change to the shortwave that changes
stored values and cannot be run in the batch that made it, which is exactly the
condition A3 attaches a prediction to. It is NOT a swept key: there is no
namelist number to move, and the change either preserves the broadband ocean
albedo exactly or it is wrong. That distinction is the whole of the check.

**What changed.** `radmod`'s upward loop replaces the open-ocean part of the
direct-beam surface albedo with a zenith fit, ECHAM-3 under `necham` or Briegleb
under `necham6`. Both fits are broadband and spectrally flat, and both were
written into `dsalb(1,:)` and `dsalb(2,:)` unchanged, so over ice-free ocean the
band index of the direct-beam surface albedo held one value. `solarini` had
already integrated `exoplasim/surfacespecs.py`'s ocean reflectance against this
star's spectrum on either side of 0.75 um to get `doceanalb(1)` and
`doceanalb(2)`; the overwrite discarded that split every radiation step, over
most of the modelled planet, for the term the shortwave weights most. The
diffuse component escaped, because `zra1s` and `zra2s` are read before the
overwrite.

The fit now carries `doceanalb`'s own band ratio across the replacement. The
shape is the spectrum's; the magnitude and the zenith dependence remain the
fit's. The shape may be applied to a directional quantity because the angular
dependence of a water surface is Fresnel's, a function of the refractive index
alone, and water's index moves by about three per cent across the 0.75 um split.
The part of `doceanalb` that is not Fresnel is water-leaving reflectance, which
is pigment and particles and which OCN-14 owns.

**Predicted magnitude, and it is small on purpose.** The band fractions are
normalised so `zsolars(1)*zofrc1 + zsolars(2)*zofrc2 = 1` identically, so the
flux-weighted broadband ocean albedo is unchanged cell by cell. The direct-beam
ocean albedo moves by +9.85 per cent in band 1 and -6.10 per cent in band 2 at
the `doceanalb` pair `solarini` computes today. Recomputing the same integral
independently on the same 965-point grid against `k25v.dat` gives +10.78 and
-6.67 per cent, and the 1.5 per cent gap in band 1 is world-a0y's `zdenom1`
mis-normalisation reached from a second direction; the prediction holds at
either value.

Because the broadband is preserved, the change reaches the top of the atmosphere
only through the difference between the two bands' two-way atmospheric
transmission, T1 - T2:

    dR_toa = f_ocean * S_toa * zsolar1 * dalpha1 * (T1 - T2)

At the ice-free ocean area fraction 0.564, a global-mean top-of-atmosphere
downward shortwave of 321.6 W/m2 (`cloud_optical_depth_bracket.json`'s band-1
incident over `zsolar1`), and a flux-weighted daily-mean ECHAM-3 ocean albedo of
0.055 to 0.080 from OCN-7's own daily means, that is **+0.02 to +0.15 W/m2 of
extra reflected shortwave, hence -0.02 to -0.13 K, COOLING**, over
|T1 - T2| = 0.05 to 0.25.

**The direction of that bracket is a claim and it is the weaker half.** T1 - T2
is taken positive because band 2 carries the water vapour bands and band 1's
attenuation is Rayleigh and ozone, which on a water-rich column is the smaller
of the two; but neither band's transmission is written out, so the sign is
argued and not measured. If T2 exceeds T1 the change warms by the same
magnitude. The MAGNITUDE bracket does not depend on the sign.

The prediction is deliberately of the same order as the 0.002 spectral bound
OCN-7 reached from the other direction, and that agreement is the point: this is
a structural correction to a physical dimension that was carrying no variation,
not a forcing. A large response would be evidence of a bug.

**What would mean wrong, in the A/B:**

- **The broadband must not move.** `zsolars(1)*dsalb(1,:) + zsolars(2)*dsalb(2,:)`
  over ice-free ocean must equal the zenith fit's own value to machine
  precision, and `dalb` over open water must be bit-identical to the control.
  This is the same identity `build_surface_albedo.py` already enforces on codes
  174, 175 and 176 over the whole grid. Any change in the broadband ocean albedo
  is a bug, not a result.
- **Two configurations must be bit-identical to the control.** `nsimplealbedo`
  sets the two `doceanalb` elements equal, so both fractions are 1; and
  `necham = necham6 = 0` drops the zenith terms entirely. Either arm differing
  from its control means the reduction is broken.
- **Magnitude.** A response beyond +/-0.5 K global mean is more than three times
  what the arithmetic above allows at its widest and points at the
  implementation or at the conversion.
- **Sign, conditionally.** Cooling is expected but is not a refutation if it
  reverses, because T1 - T2 is argued rather than measured. What WOULD refute
  the chain is a response whose magnitude requires |T1 - T2| above 1.

**Not blocked on anything.** The change is in the tree and needs no key set. It
lands with whatever pair runs next and does not need an arm of its own; it is
registered so that a bundle total which includes it is checked against a
prediction that named it.

**What OCN-7 still cannot deliver.** Its remaining half is a sourced bound on
the SPATIAL span of the ocean albedo, which needs a modelled water-leaving
reflectance. OCN-14 owns that and it does not exist, so the span is unbounded
here and is not smuggled into this prediction.

## world-jgen and world-f9ig: the band-1 cloud optics, ONE T21 pair

The two changes are in the model source together, so a control on the
pre-`world-jgen` tau relation against an arm on the current source measures the
pair. Separating them needs a third arm and nothing requires that separation.

**What the arm is, and it supersedes the definition the issue was opened with.**
`radmod` evaluated one cloud optical depth, `2*ALOG10(1.5 + W)**3.9`, for BOTH
shortwave bands. Stephens (1978) p. 2124 fits two by least squares, Eq. (10a)
over 0.3 to 0.75 um and Eq. (10b) over 0.75 to 4.0 um, and `radmod` splits its
own bands at that same 0.75 um: band 2 was getting its own equation and band 1
was getting band 2's. Each band now evaluates its own fit.

The `+1.5` offset went with it, and it was NOT replaced by zero optical depth
below 1 g/m2 as the issue proposed. Zero leaves the 1 to 10 g/m2 window on an
extrapolated fit whose implied effective radius, through Stephens Eq. (7), is
134 um at 2.15 g/m2. Below 10 g/m2, which is the bottom of the paper's Figs. 1a
and 1b, each band is continued LINEARLY in the water path, matching its own
fitted value at 10 -- Eq. (7) at the effective radius the fit implies there,
8.2 um in band 1 and 6.7 um in band 2. Continuous, zero at zero water, and no
new constant.

**The hold on `tswr1` is lifted, because there is no coefficient left to hold.**
The arm was written to hold `tswr1` at its inherited value, on the ground that
it had been tuned against the optical depth being corrected. `world-f9ig` has
since deleted `tswr1`, `tswr2` and `tswr3` and put Stephens et al. (1984) Tables
1(a) to 1(c) in `swr` in their place, so the band-1 backscatter is interpolated
from Table 1(b) at the layer's own range-1 optical depth. The control and the
arm no longer share a coefficient fitted against the quantity being corrected,
which was the whole reason the hold existed.

**Predicted magnitude: +0.73 to +2.00 W/m2 at the top of the atmosphere, +0.59
to +1.63 K, WARMING.** `exoplasim/scripts/cloud_optical_depth_bracket.py`,
re-derived on the current source. It carries the model's own per-layer cloud
water paths, its own total cloud cover with random overlap, its own zenith
geometry integrated over the day, and the band-1 surface albedo through Stephens
Eq. (12).

**The correction is not one-signed, which the issue's original hand chain got
wrong.** Over the seven layers carrying real cloud water, 7 to 90 g/m2, band-1
cloud gets LESS bright and that warms; in the top three layers, under 2.5 g/m2,
it gets brighter, because the linear continuation returns more optical depth
than the `+1.5` offset did. The second is much the smaller, and the bracket
above is the net.

**The bracket leans high** and the direction is argued rather than assumed:
band-1 gas absorption and Rayleigh scattering above and below the cloud are not
in the chain and both attenuate the change, and the column is composed by random
overlap without interlayer multiple reflection, which also attenuates it. Its
WIDTH is the spread over where in the column the cover sits, which the model
does not write out. The geometry is checked rather than trusted: the chain's
band-1 incident flux reproduces the climatology's own `rst + rsut` to better
than 1e-4, which is a quantity the other side already knows.

**What would mean wrong, in the A/B:**

- **Sign.** A correct band-1 tau is smaller than the one it replaces over every
  layer carrying real cloud water, so the arm must warm. A cooler arm is wrong
  outright, not a small result.
- **Magnitude.** Below +0.59 K or above +1.63 K is outside the bracket. Above it
  in particular points at the linear continuation, which is the only part of the
  change with no paper behind its FORM, though its endpoint is Stephens Eq. (7).
- **The top three layers.** Their contribution is the one that opposes, and it
  is small. If the arm's response is dominated by the layers under 2.5 g/m2, the
  continuation is doing more than it should be and the arm is measuring it
  rather than the band split.
- **The instrument.** Stephens (1978) p. 2127 warns that the
  cloud-over-reflecting-surface correction is unreliable above a surface albedo
  of about 0.75, which is every modelled snow and sea-ice cell, so a response
  concentrated there is outside the parameterisation's own stated domain.

**BLOCKED on the run only.** The source change is in the tree and the bracket is
current. What is missing is a T21 pair on a settled baseline, and
`config/planet.yaml` declares `baseline_climatology: null`.

## world-2esd: `th2oc`, the model's whole longwave continuum, swept

Registered as a sweep rather than a correction, because nothing this project
holds determines the value and the honest disposition is a declared bracket.

**What it is.** `radmod.f90` declares `th2oc = 0.024` and `lwr` applies it once,
as `zah2o = min(zah2o + (1 - exp(-th2oc*zsumwv)), 1)` on the pressure-weighted
water path. Sasamori (1968) carries no window absorption at all, so this single
line is the entire continuum contribution to the modelled longwave and the only
thing absorbing in the modelled window. It is the fourth member of the
per-truncation `jtune` table upstream carried; the other three are gone with
`world-f9ig`, they named themselves as tunings in their own declarations, and
this one never did.

**The ceiling is derived; the value is not.** The term adds broadband
absorptance that Sasamori leaves out, and what Sasamori leaves out is the
window, so the addition cannot exceed the share of the emitted Planck flux the
window carries. At the largest pressure-weighted path the model reaches, 6.08
g/cm2 on a T21 bootstrap climatology, and the 8-12 um window's Planck share over
250 to 300 K, that is `th2oc <= 0.038 to 0.050`.

**The inherited value is half of its own ceiling.** At the mean path, 2.23
g/cm2, the term claims 0.052 of broadband absorptance; at the maximum path
0.136, which is more than half the whole window's Planck share. This is a
first-order term carrying a residual, not a small correction, and that is the
finding rather than the arms.

**The bracket is [0.0, 0.038] and both ends are arguable.** Zero is what the
surrounding `if(th2oc > 0.)` supports and is the arm that measures the term's
entire worth; it is a bound and not a candidate, because the window continuum is
real and `shine2012` p. 536 puts it as "often the dominant cause of
wavelength-averaged absorption" in the windows. 0.038 is the tightest ceiling
across the modelled emitting temperatures. `exoplasim/scripts/lw_continuum_bracket.py`
re-derives both ends and the path they rest on.

**THE ARMS ARE THE TWO ENDS AGAINST EACH OTHER**, on the same instrument
argument `world-u9hq` makes: a single end against the inherited value sits too
close to the convergence criteria's own slack to be attributed. End to end is
0.081 of broadband absorptance at the mean path, which against a clear-sky
greenhouse trapping of order 150 W/m2 is roughly 12 W/m2 and roughly 8 K by
`lib/sensitivity.py`. **That conversion is static and leans high**: it carries no
lapse-rate response and no water vapour feedback, both of which a run has and
both of which damp it. The prediction carried into the arms is therefore an
ORDER, 4 to 9 K end to end, and its purpose is to establish that the term is
first-order rather than to predict the number.

**What would mean wrong, in the A/B:**

- **Sign.** More continuum absorption must warm the simulated mean. A cooler
  high arm is wrong outright.
- **Monotonicity** against the inherited 0.024 on the same lineage: 0.0 below it
  and 0.038 above it. An inversion says the difference is scatter.
- **Magnitude.** Below 2 K end to end contradicts the absorptance arithmetic,
  which is not an estimate but a direct evaluation of the term the model
  applies, and would point at the implementation. Above 12 K exceeds the static
  conversion, which already leans high, and would point at a feedback the
  conversion has no way to carry.
- **The instrument, and this arm is unusually safe from class 34.** The effect is
  0.081 of broadband absorptance where the convergence criteria resolve
  fractions of a W/m2; nothing here is near the noise floor. The risk runs the
  other way: an arm at 0.0 may not equilibrate on a commissioning-length window
  at all, and a pair that has not settled is not a slope.
- **An arm staged at 0.024 must be bit-identical to the control**, since that is
  the value the binary already carries.

**BLOCKED on the run, and on a source for the value.** The bracket and the
declaration are in the tree. Closing it properly, rather than sweeping it, needs
a correlated-k or line-by-line calculation with the MT_CKD continuum on this
path; `references/INDEX.md` records that Mlawer et al. (2012) does not supply an
evaluable continuum, because the coefficients ship as data with LBLRTM. That is
the same bundle `exoplasim/notes/corrk-cross-check.md` says the project does not
have and which blocks `h2o_sw_level`'s bracket. The two open together and should
be scoped together.

## The bundle splits in two, and only one half can be attributed at all

Measured against the tree on 2026-08-26, model source `a041a1e9`. A3 says bundle
the terms and check the SUM against the sum of predictions, and it attaches two
conditions to the A/B that decide which terms can enter such a sum: both arms
use the SAME BINARY and differ only by a namelist key, and both branch from ONE
restart. The five predictions batch 2 registered do not all meet the first.

**The half that is namelist-revertible on one binary.** Each of these has a
sentinel or a literal in a namelist that restores the pre-change constant
exactly, so a control arm and a bundle arm can share a binary:

| term | key | namelist | what restores the old behaviour |
| --- | --- | --- | --- |
| `world-trs3` | `gamma` | `rainmod_nl` | `gamma = 0.01`; the sentinel is `-1.0` and `> 0` restores the literal at all four sites |
| `world-o12h` | `rcritwidth` | `rainmod_nl` | `rcritwidth = 1.0`; the sentinel is `-1.0` and derives from NLAT |
| `world-2esd` | `th2oc` | `radmod_nl` | `th2oc = 0.024`, the compiled value, and the sweep moves it |
| `world-s8rv` | `clwref` | `rainmod_nl` | `clwref = 0.00021`, the compiled value |
| `world-u9hq` | `h2oswl` | `radmod_nl` | `h2oswl = 1.0`, the model default |

**The half that is compiled in, where no namelist reverts anything.**

- **`world-jgen` and `world-f9ig`, the band-1 cloud optics.** Stephens Eq. (10a)
  against Eq. (10b), the `+1.5` offset replaced by a linear continuation, and
  `tswr1` to `tswr3` deleted in favour of Stephens et al. (1984) Tables 1(a) to
  1(c). This section says it plainly: a control on the pre-`world-jgen` tau
  relation against an arm on the current source measures the pair. That is TWO
  BINARIES, which A3 admits only where the term IS compiled, as here, and the
  confound it warns about is bounded by measurement rather than by argument: the
  control built for this reverts one file and 37 of the 38 objects come out
  bit-identical to the shipped binary's.
  **A control executable was preserved all along.** `install_run_executable`
  copies the binary into the run directory, so `run_14906cb7b914` still holds
  `most_plasim_t21_l10_p16.x` at the sha its `INDEX.json` entry names, and
  `strings` finds `tswr3` in it and `cloudabs` not at all: the donor is
  pre-world-f9ig, which is what puts the cloud optics inside the drift window.
  A preserved binary wants the staging of its own date, so it controls a
  matching namelist rather than the current one, which is why the arm run here
  was built from source instead.
- **`OCN-22`, the open-ocean albedo's band split.** Its own section already says
  it is not an arm and that there is no namelist number to move. Its two
  bit-identity conditions, `nsimplealbedo` and `necham = necham6 = 0`, are
  REDUCTIONS that check the new code collapses correctly; neither isolates the
  band split from everything else those switches change, so neither is a revert.

**`world-o12h` is namelist-revertible and still contributes exactly zero to any
bundle run at T21.** `rainmod.f90` declares `RCNLATREF = 32` and derives
`rcritwidth = (RCNLATREF/NLAT)^(1/3)`, which at T21's NLAT of 32 is exactly 1,
and the applying loop is guarded by `if(rcritwidth /= 1.)`. Its own prediction
says so first: T21 is bit-identical to the control, by construction rather than
by rounding. The term has a value to measure only at T42 and above, so it cannot
enter a T21 bundle even with a route, and a T21 bundle total that named it would
be claiming a measurement of zero as an agreement.

**Three of the five have no route from `config/planet.yaml`.** `gamma`,
`rcritwidth` and `th2oc` are read from their namelists by the model and are
written there by nothing: `run_exoplasim.py` has no `declare_` function for any
of them, `expected_namelist_keys` does not cover them, and `continue_exoplasim.py`
does not reapply them per segment. `clwref` and `h2oswl` are the two that do have
one, which is why they are the two arms this batch could run. A hand-edited
namelist is not a substitute, because `verify_staged_namelists` is what makes an
arm's own namelist the record of what it integrated, and a key it does not know
about is not in that record.

**So the bundle sum cannot be checked against the sum of predictions.** Of the
five registered terms, two are unmeasurable on any pairing A3 permits, one is
identically zero at the rung a T21 bundle would run at, and three lack the route
that would let a control arm be staged with provenance. What remains is not a
bundle: it is two terms with routes, and those are worth running as their own
pairs, which is what `world-s8rv` and `world-u9hq` are.

This is not a reason to serialize the iteration. A3's argument that attribution
cannot gate anything still holds, and the terms are in the tree because they are
correct rather than because a comparison improved. What it changes is the claim
that can be made afterwards: the next baseline's total is not attributable to its
parts, and a note that says otherwise would be describing an A/B nobody can run.

**What would make the bundle checkable**, in the order that costs least:

1. A `declare_` route and an `expected_namelist_keys` row for `gamma`,
   `rcritwidth` and `th2oc`, the way `world-n1nu` built one for `clwref`. That is
   a `run_exoplasim.py` change and stales no binary.
2. Preserving one binary per rung at each source change that lands a compiled-in
   forcing term. The two compiled-in terms here are unattributable for want of a
   control executable that cost nothing to keep at the time.
3. Running `world-o12h` at T42 against a T42 control, which is the only rung
   where it has a value.

## Measured: the bundle is an order of magnitude larger than its predictions, and the other sign

2026-08-26, T21, `most_plasim_t21_l10_p8.x` at model source `a041a1e9`, five arms
branched from `run_14906cb7b914`'s `MOST_REST.00034`, host load 25 to 50 over 32
cores. That donor is 35 orbits at 292.6 K, integrated by the SUPERSEDED source
and its config, which is what makes it a usable initial condition and not a
baseline.

**Every arm relaxes to a state about 12.3 K colder than the donor.** The control
arm, which carries `config/planet.yaml`'s own values for everything except the
two energy keys:

| orbit | 0 | 4 | 9 | 14 | 19 | 24 | 31 | 59 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| global mean `ts`, K | 290.78 | 284.81 | 281.76 | 280.98 | 280.37 | 280.18 | 279.80 | 279.96 |
| sea ice, fraction | 0.0099 | 0.0340 | 0.0568 | 0.0667 | 0.0728 | 0.0737 | 0.0737 | -- |

Sea ice reaches thirteen times the donor's 0.0058 and the ice-albedo feedback is
carrying most of the amplitude. The control was later taken to sixty orbits and
settles: its last ten sit at 279.96 K, flat to 0.1 K, against the donor's last
orbit at 292.26. **The drift is -12.3 K and equilibrated**, not a lower bound.

**The sum of the registered predictions is about +1 K of WARMING.** `world-jgen`
with `world-f9ig` is +0.59 to +1.63 K warming and is the largest single term;
`OCN-22` is -0.02 to -0.13 K; `world-o12h` is identically zero at this rung;
`world-2esd` and `world-trs3` are not swept in the baseline config. Nothing in
the table predicts a term of this size and nothing predicts cooling of this size.

**What the drift is NOT.** The arms differ from the donor's configuration in four
places and only one is large enough to matter for sign:

| difference | worth |
| --- | --- |
| `energy_fixer` and `energy_diagnostics` off | +0.18 K, and WARMING, so it makes the gap larger rather than smaller. The fixer applies -0.300 W/m2 on the donor and `lib/sensitivity.py` at that run's own planetary albedo of 0.2366 gives 0.615 K per W/m2 |
| `ncpus` 16 to 8 | none. Thread count is compiled in and changes no physics |
| four `surface.cryosphere` constants, absent in the donor's config and declared now | **negligible, and it is a static check rather than an arm.** Three of the four declare the value `icemod.f90` already compiles -- `CRHOI` 920, `CPI` 2070, `CKAPI` 2.03 -- so declaring them changed nothing the model reads. The fourth, `CLFSN`, moves from the compiled 3.337E5 to the derived 333444.87, which is -0.076 per cent of the modelled snow's melting enthalpy |
| `surface.land_water_column` from a 1-layer bucket to a 2-layer scheme | unbounded here |

So the measured drift confounds the compiled forcing terms with two config
changes, and the split above says which. **It is not attributable further**, for
the reasons the previous section gives: two of the five terms are compiled in
with no namelist that reverts them, three have no route from config, and no
pre-batch-2 binary was preserved to serve as a control.

**A residual this large is the shape of one term, not of two that nearly
cancel**, which is what `docs/src/practice/failure-modes.md` class 15 would
otherwise warn about. The section "The bisect" below finds that term: it is
`world-trs3`'s derived `gamma` and finds it worth -1.89 +/- 0.08 K, which is a
seventh of the drift and not the whole of it. About -10.8 K is still
unattributed, and the terms it must live in are the two with no control binary.

**What this does not license.** Nothing here says a term should be removed.
Physics is not a knob, and a correct term that moves the simulated climate a long
way is information. What it says is that the next baseline cannot be taken on
this configuration until the two config changes above are priced, because they
are the two candidates that are cheap to price and are not yet.

**The split needs one arm, not two.** The cryosphere row above is settled by
reading `icemod.f90` against `config/planet.yaml`, and it is settled at
essentially zero, so the only config change left that could carry a response of
this size is the land water column. `run_431ecabed085` is that arm, 20 orbits of
a declared 25 before the host was taken over: the 1-layer
bucket the donor ran, against this same control on this same restart, staging
`NLANDWCOL = 0`, `DSOILWZ = 1.5` and `DSOILWF = 1` where the control stages
`NLANDWCOL = 1`, `DSOILWZ = 0.5, 1` and `DSOILWF = 0.333333, 0.666667`.

It accounts for +0.17 +/- 0.07 K of the drift and no more, so the config is not
where the drift lives and the rest is the model source's. The bisect below takes
it from there.

Removing the `cryosphere` block to build a second config arm does not work and
should not: `derive()` indexes it rather than using `.get`, deliberately, so that
a config which has lost the block fails instead of reverting the model to its
compiled values in silence.

## A3's fourth condition, and how an arm now demonstrates its own settling

A3 requires that "the segments are labelled as diagnostics, so a short A/B tail
never enters a convergence window or a climatology", and world-u9hq recorded the
condition as UNTESTABLE rather than unmet: `assess_convergence.py` read only
PRODUCTION orbits, so the label that satisfied A3 was exactly what hid an arm
from the test that would judge whether it had settled. No arm could satisfy both
halves, and the h2oswl pair above stood on its paired difference's
autocorrelation-corrected error instead.

Both halves are closed, and they had to be closed together: a guard on the label
alone would have left every arm unable to show its own settling.

`run_exoplasim.py` takes a required `--purpose`, `spinup` or `diagnostic`, with
NO DEFAULT. A default of `spinup` would leave the failure where it was and
silent in the direction that matters, because a diagnostic tail carrying the
default label is indistinguishable later from production orbits. The purpose
reaches the first segment and the manifest, and `diagnostic` sets
`canonical_lineage_eligible = false` -- which until now only `--binary` did, so
an A/B run the way A3 asks for it, as one binary differing by one namelist key,
was eligible for the canonical lineage on every arm this project has run.
`continue_exoplasim.py` already refuses a `post_equilibrium_climatology` segment
on a run carrying false, so the arm cannot acquire climatology orbits later
either.

`assess_convergence.py --assess diagnostic` takes the same six criteria over the
orbits a run declares ARE diagnostics, with `production_window`'s own
trailing-drop and interior-hole rules applied to that purpose instead. The
verdict it produces is about the EXPERIMENT and not about the planet: whether an
arm had settled far enough for its difference against another arm to mean
anything. It carries its mode in its report and plot filenames and writes
nothing into the run -- not `status`, not `equilibrium_cutoff_year_index`, not
`convergence_assessment` -- on the terms `--through` is already kept out, and
for a stronger reason: those fields are claims about where the planet has got
to, and a diagnostic segment is by declaration not evidence about that.

So an arm's settling is now demonstrable and reported, and it still cannot be
mistaken for the run's own verdict or promoted into the canonical lineage.
world-ucww.

## Measured: the `world-u9hq` h2oswl arms, 2026-08-26

T21, `most_plasim_t21_l10_p8.x` at model source `a041a1e9`, both arms branched
from `run_14906cb7b914`'s `MOST_REST.00034` and differing from the control in
`H2OSWL` alone. `run_11b43d2c56a4` carries 1.129 and `run_fc75ca5f9fd2` carries
1.206; `run_57a43e1fc3f4` is the central 1.163. The energy fixer and diagnostics
are off in all three, forced by the `epilog` defect and common to every arm.

The declared length is 11 orbits of settling, `lib/run_lengths.settling_bracket(0.75)`'s
lower end, then a 35-orbit measurement window, which is `assess_convergence.py`'s
own derived minimum. The high arm reached all 46. The low arm reached 40 before
the host was taken over, so the paired comparison runs to 40 and its window is 29
orbits. Host load 25 to 50 over 32 cores throughout, and 40 to 51 for the last
third, so no wall clock from these runs is worth keeping.

| comparison | window mean, K | window |
| --- | ---: | ---: |
| 1.206 minus 1.129, end to end | **+0.599 +/- 0.087** | 29 orbits |
| 1.129 against the central 1.163 | -0.363 +/- 0.124 | 29 orbits |
| 1.206 against the central 1.163 | +0.258 +/- 0.055 | 35 orbits |

**Sign: holds.** More absorption at the level the key scales warms the modelled
mean, and the high arm is the warm one. Every orbit of the end-to-end difference
is positive.

**Monotonicity: holds.** 1.129 sits below the central and 1.206 above it, both
resolved, with no inversion. The two one-sided differences sum to the end-to-end
one.

**Magnitude: holds.** The prediction carried in was 0.72 K end to end, PHYS-9's
fitted asymptote slope, inside a static bracket of 0.54 to 1.16 K. The measured
+0.599 +/- 0.087 K is inside that bracket and sits 1.4 standard errors below the
point prediction. It lands nearest PHYS-9's WINDOW-MEANS slope, which predicted
0.57, rather than the asymptote the prediction was taken from; the note's reason
for preferring the asymptote was that the window means came from an unconverged
run and were a lower bound, and on this pair they were not low.

**The instrument condition as registered cannot be tested, and that is
structural.** It asks that each arm meet all six convergence criteria on its
window. `assess_convergence.py` refuses both arms: `segments.production_window`
counts only PRODUCTION orbits, and the 35 orbits were added with
`--purpose diagnostic` because A3's fourth condition requires exactly that. So
A3's labelling is what makes the orbits invisible to the convergence instrument,
and no arm can satisfy both conditions at once. At 11 production orbits the arms
fail all six criteria, which is what a settling block should do.

What stands in its place is the paired difference's own scatter, corrected for
memory by `lib/autocorrelation.py`: +0.599 against a standard error of 0.087 is
6.9 to one, and the top-of-atmosphere imbalance has closed to -0.074 +/- 0.085,
so the pair is not still separating. That is an instrument and it is not the one
the prediction named.

**Which of the two conditions gives is a decision this measurement does not
make.** Either an A3 arm's segments are production and the convergence criteria
apply while A3's guard is off, or they are diagnostics and the registered
condition is replaced by the paired-difference form above. The second is what
these arms did.

## The bisect: which code path is not doing what it was written to do

A bundle missing its predicted sum by an order of magnitude and reversing sign is
the shape of a defect rather than of bad estimates, so each term below is treated
as a suspect implementation until cleared, and a check with a right answer is
preferred to a run difference wherever one exists. Verdicts are CLEARED,
DEFECTIVE or UNTESTED at this rung.

**First, the instrument, because a control that cannot reproduce itself would
make every number here worthless.** `Ct` and `Dt` do not reproduce between
identical runs. **No number in this note inherits that.** Those two columns exist
only under `nenergy > 1`, which allocates `zcnow`; every arm in this note has no
`NENERGY` key in its `plasim_namelist` at all, so `zcnow` is never allocated and
the decomposition never runs. The model state itself is bit-reproducible on this
source, shown directly by `run_0730a12ecfbd` and `run_57a43e1fc3f4` writing an
identical `MOST_REST.00000` from different cores. The instrument is clean for
every kelvin below.

**The split, measured.**

| term | worth at T21 | how it was established |
| --- | ---: | --- |
| `world-trs3`, the derived `gamma` | **-1.89 +/- 0.08 K** | A/B, `run_af3d2c9a4b05`, 25 orbits, TOA closed |
| `world-py6p`, the land water column | +0.17 +/- 0.07 K | A/B, `run_431ecabed085`, 20 orbits, TOA -0.14 +/- 0.13 |
| `surface.cryosphere` declared | about 0 | static, against `icemod.f90` |
| the energy fixer off | +0.18 K | the donor's manifest and `lib/sensitivity.py` |
| `world-o12h`, `rcritwidth` | exactly 0 | construction, `RCNLATREF` = `NLAT` |
| `world-2esd`, `th2oc` | 0 | not swept; the config value is unchanged |
| `world-f9ig` with `world-jgen`, the cloud optics | **-10.01 +/- 0.07 K** | A/B on a source-built control binary, `run_a2b1bd1e859f`, 25 orbits |
| `world-80ia`, `vdiff_lamm` | -0.15 +/- 0.06 K, not resolved | A/B, `run_5373310a7b9f`, 25 orbits |
| **residual, unattributed** | **-0.60 K** | by difference |

**The control has settled, so the total is not a lower bound.** Its last ten
orbits of sixty sit at 279.96 K against the donor's last at 292.26, flat to
0.1 K, so the drift is **-12.3 K** and equilibrated. The namelist-revertible
rows above sum to -1.89 + 0.17 + 0.23 = -1.49 K, the two config rows being
WARMING and so making the cooling smaller rather than larger, which left about
-10.8 K for the compiled-in half.

**That residual is `world-f9ig`**, and the section "The bisect finishes" below
measures it at **-10.01 +/- 0.07 K**: reverting the shortwave cloud optics to
their pre-jgen form on a binary that differs in one object file returns
+18.43 W/m2 of top-of-atmosphere shortwave from the first orbit and +10.01 K at
equilibrium, on the same restart. The bracket for it was already in the tree, in
`exoplasim/analysis/stephens_tables_vs_fits.json`, at -25.35 to -6.54 W/m2 and
-20.63 to -5.32 K. What is left over after it is -0.65 K.

### `world-trs3`, the derived Kessler `gamma`: CLEARED, and worth a seventh of the drift

Four checks with right answers, all passed:

- **The formula reproduces the registered offline table** at 0.5, 1, 3, 10 and
  30 mm/day, to the three digits the table carries, with `GAMPEXP` = 13/20 * 8/9
  and `zvcoef` = 5.17*sqrt(ga/9.80665) read from the source rather than assumed.
- **All four re-evaporation sites are structurally identical**: `zgam = gamma`,
  then the derived form on that site's OWN precipitation flux when `gamma <= 0`,
  then `AMIN1(..., zpr)` capping removal at the flux available. The two snow
  sites use `ALS` where the two rain sites use `ALV`.
- **The override is reached**, and the model says so: `run_af3d2c9a4b05`'s
  `plasim_diag` prints `precip re-evaporation: CONSTANT gamma 1.0E-002`. That is
  the first thing this note's own falsifying list says to check.
- **Every link of the stated mechanism has the predicted sign.** Over the arm's
  settled window, orbits 15 to 24, pinning `gamma` to 0.01 against the derived
  default gives +0.291 +/- 0.011 mm/day more precipitation reaching the ground,
  0.024 +/- 0.002 less cloud cover and 6.86 +/- 0.24 W/m2 less reflected
  shortwave. Reversed, that is the derived form moistening the column, growing
  cloud and reflecting more, which is the chain this note describes.

**The response saturates, exactly as this note predicted it would.** The surface
separation runs +0.36, +1.11, +1.72, +2.11 over the first seven orbits, reaches
about +2.2 by orbit nine and then flattens; the settled window gives
**+1.89 +/- 0.08 K** with the top-of-atmosphere difference closed to
-0.22 +/- 0.16 W/m2. The shortwave forcing decays with it, from 10.4 W/m2 in the
first orbit to 6.9 in the settled window, because the cloud difference itself
shrinks as the two states converge. That is the sublinearity this note argued
for: the equilibrium is reached by MOISTENING until the sub-cloud deficit falls,
not by the evaporated flux scaling with `gamma`.

**DO NOT CONVERT THE FORCING THROUGH THE STATIC SLOPE HERE.** 9.2 W/m2 at
`lib/sensitivity.py`'s 0.615 K per W/m2 is about 5.7 K, and the arm's own answer
is 1.89. The slope is local, the forcing is not constant while the state responds,
and the arm closes its own budget; the run's mean is the measurement and the
conversion is not. This is `docs/src/practice/failure-modes.md` class 34 from the
inside, and it was made once on this arm at ten orbits before the arm had settled.

**So the code is doing what it was written to do, and the prediction is very
nearly right.** It predicted the moistening, the cloud increase, the drop in
precipitation and the saturation, and all four happened. The precipitation
reduction, 13.5 per cent of the global mean, sits inside the 10 to 50 per cent it
gave for its land counterpart. What it did not carry is a number for the
radiative channel: "high cloud fraction rises slightly" is +0.024 of cover worth
6.9 W/m2 and 1.9 K, which is not slight and is the largest single effect the term
has. A term whose stated mechanism runs through cloud needs a radiative number
attached before it lands, not a qualifier.

### `world-jgen` and `world-f9ig`, the cloud optics: CLEARED, and world-f9ig is the drift

- **The Stephens coefficients are self-consistent** with the fits they name:
  10^0.2633 = 1.8336 and 1.7095*ln10 = 3.9363 for Eq. (10a), 10^0.3492 = 2.2346
  and 1.6518*ln10 = 3.8034 for Eq. (10b).
- **The linear continuation is continuous and grounded.** Below `zwfit` = 10 the
  base of the power is clamped to 10, so `log10` is exactly 1, the power is 1,
  and what remains is `ztaua_b * zlwp/10`: linear, equal to the fit at 10, zero
  at zero water. No offset and no domain error.
- **The offline chain still reproduces the model's own geometry.**
  `cloud_optical_depth_bracket.py` on `run_57a43e1fc3f4` returns a reconstructed
  band-1 incident flux of 122.97 W/m2 against the climatology's 122.97, +0.00 per
  cent, and the bracket is unmoved at +0.73 to +2.00 W/m2.
- `smoke_test.py` gates the shortwave cloud tables against the papers' tables.

**One thing the registered wording understates.** "Band 2 was getting its own
equation" reads as though band 2 did not move. It did: its prefactor went from 2
to 2.2346, its exponent from 3.9 to 3.8034, and the `+1.5` offset left with the
rest. Over the layers carrying real cloud water that is -0.8 to +4 per cent of
optical depth, but below 2.5 g/m2 both bands gain 100 to 300 per cent, which is
the continuation replacing the offset. That is the condition this note's own
falsifying list names, and it is now separated: a source-built control binary
puts the two commits together at +18.43 W/m2, and the pair already on disk
(`run_d35b554cab2a` against `run_57cecaa8391b`, both binaries pre-f9ig) puts
world-jgen's tau relation alone at 0.13 K. The optical depth is tenths of a
kelvin; the TABLES are the term. "The bisect finishes" below has both.

### `world-o12h`, `rcritwidth`: UNTESTED at this rung, not verified as zero

At T21 `RCNLATREF` and `NLAT` are both 32, the derived width is exactly 1, and
`if(rcritwidth /= 1.)` skips the transform. Measuring zero here is correct and
expected and is NOT evidence about the path: the unit factor is skipped rather
than applied, so nothing about the scaling is exercised. **The path is unverified,
and T21 cannot verify it.** A T42 pair can, and now has a route: `RCRITWIDTH`
reaches `rainmod_nl` from `config/planet.yaml` as of this batch, so an arm at
1.0 against the derived 0.7937 is a namelist key rather than a code fork.

### `OCN-22`, the ocean albedo band split: CLEARED by its own identity

Compiled in, and no namelist reverts it, but the identity it is built on is
checkable without a run and it holds. `swr` forms
`zoalbb = zsolars(1)*doceanalb(1) + zsolars(2)*doceanalb(2)` and then applies
`zofrc_b = doceanalb(b)/zoalbb` to a zenith fit that is the same in both bands,
so the flux-weighted broadband is `sum_b zsolars(b)*zofrc_b = zoalbb/zoalbb`,
exactly 1, and the split moves the band-1 and band-2 albedos apart without
moving what the two together absorb.

On the model's own numbers, out of `run_c9c24d438a94`'s `plasim_diag`: ocean
albedo 0.0755700 below 0.75 um and 0.0636681 above, overall 0.0682190. Those
three fix the band-1 flux share at 0.382372 and band 2 at 0.617628, which is the
0.6176 this note's PHYS-11 section derives independently, and the reconstruction
`z1*zofrc1 + z2*zofrc2` returns 1.000000000000. Its predicted size, -0.02 to
-0.13 K, is two orders below the drift, and it is not a candidate for it.

### `world-2esd`, `th2oc`: not in the bundle

The config value is unchanged, so the term contributes nothing to the measured
drift and there is nothing to bisect. It still has no route from config, unlike
`gamma` and `rcritwidth`, and closing that is the same one-row change.

## The bisect finishes: the residual is world-f9ig, and the tree had already priced it

2026-08-26, T21. The section above left about -10.8 K unattributed and said it
had nowhere to live except the compiled-in half of the bundle. It lives in one
term of that half: `world-f9ig` (c8debdc1), which replaced `tswr1`, `tswr2` and
`tswr3` with Stephens, Ackerman and Smith (1984) Tables 1(a) to 1(c).

**The number was in the tree before the bundle ran.**
`exoplasim/analysis/stephens_tables_vs_fits.json` prices exactly this swap at
**-25.35 to -6.54 W/m2** over the layers carrying cloud water, which it converts
to **-20.63 to -5.32 K**, and its own `settles_it` field names the experiment
that would close it: "a T21 commissioning pair on the adopted scheme".

The bundle registered **+0.59 to +1.63 K of warming** for "world-jgen and
world-f9ig" jointly. That figure comes from
`exoplasim/analysis/cloud_optical_depth_bracket.json`, which prices the OPTICAL
DEPTH swap alone -- world-jgen -- and the tables artifact names it as its own
`bracket_source`. So one of the two halves was carried into the bundle at its
full weight and the other, an order of magnitude larger and of the opposite
sign, was not carried at all. The bundle then measured -12.3 K, went looking for
the term among its five registered predictions, and could not find it there,
because the term that explains it had a bracket in the same directory and no
prediction row.

**The prediction's own falsifier fired, and it is the sign one.** It registered
"a correct band-1 tau is smaller than the one it replaces over every layer
carrying real cloud water, so the arm must warm. A cooler arm is wrong outright,
not a small result." The arm reverting the change is WARMER, which is the same
statement: on one restart the current source sits 7.6 K below the pre-jgen
control by orbit 6 and is still separating, and it reflects 18.4 W/m2 more
shortwave from the first orbit. The sign argument is sound about the band-1
optical depth and says
nothing about the backscatter fraction, which is the larger term in the same two
commits and runs the other way: this is what it costs to size a prediction
covering two changes from a bracket that prices one of them.

**Measured: the T21 pair the artifact asked for.** A control and two arms
branched from `run_14906cb7b914`'s `MOST_REST.00034`, all at T21 on `p8`, each
arm differing from the control by exactly one thing. The control is
`run_c9c24d438a94`, on `config/planet.yaml` as it stands. Orbit 0 is the same
initial state integrated by two executables, so the top-of-atmosphere difference
over it is the FORCING and not a response.

| arm | run | what it reverts | d rst | d rsut | d rss | d ts, orbit 0 |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| P | `run_a2b1bd1e859f` | world-jgen and world-f9ig in `swr` | **+18.43** | +18.43 | +12.55 | +1.56 |
| V | `run_5373310a7b9f` | `vdiff_lamm` to 160 m | -0.05 | -0.05 | -0.16 | -0.05 |

Fluxes in W/m2, arm minus control. **+18.43 W/m2 of shortwave the current cloud
optics reflect and the pre-jgen ones do not**, one-signed, and inside the
artifact's -6.54 to -25.35 W/m2 bracket. The measured drift of -12.3 K is inside
its -4.21 to -16.31 K. The bracket was right and nothing read it.

**The modelled sea ice separates the forcing from its amplifier.** Over the same
first six orbits the control's sea-ice fraction runs 0.0097, 0.0154, 0.0216,
0.0283, 0.0349, 0.0414 and reaches 0.074 by orbit 18, thirteen times the donor's
0.0058. Arm P's runs 0.0078, 0.0088, 0.0102, 0.0101, 0.0102, 0.0110 and stays
there. So the drift is not the shortwave forcing alone: the forcing cools the
modelled ocean, the ice-albedo feedback takes it the rest of the way, and
reverting the cloud optics stops both. That is why the section above reported
sea ice carrying most of the amplitude without being able to say what started it.

**This forcing PERSISTS, where `gamma`'s decayed.** Over the arm's first six
orbits the top-of-atmosphere shortwave difference runs 18.43, 17.41, 16.11,
16.50, 17.20, 16.82 W/m2 while the surface separation grows 1.56, 3.19, 4.47,
5.63, 6.46, 7.24 K. Cloud cover between the arms moves by 0.015 over the same
span. That is the signature of a change in the cloud's OPTICAL PROPERTIES rather
than in how much cloud there is: `gamma`'s forcing fell from 10.4 to 6.9 W/m2 as
its cloud difference shrank, and this one has no cloud difference to shrink. So
the sublinearity that made `gamma` worth 1.89 K instead of 7 does not apply here,
and the separation should run to the order the damping allows.

**Settled, and the budget closes.** All three arms reached 25 orbits. Over the
declared window, orbits 15 to 24:

| arm | run | window mean `ts` | TOA net |
| --- | --- | ---: | ---: |
| control | `run_c9c24d438a94` | 280.137 K | -0.525 W/m2 |
| pre-jgen, pre-f9ig | `run_a2b1bd1e859f` | 290.144 K | -0.319 W/m2 |
| `vdiff_lamm` = 160 | `run_5373310a7b9f` | 279.984 K | -0.513 W/m2 |

**Reverting the cloud optics is worth +10.01 +/- 0.07 K**, resolved at 136 to 1
against the criterion fixed before the run, with the window flat to 0.2 K and
both arms' top-of-atmosphere imbalance inside 0.53 W/m2. The prediction
registered before the arm ran was +5 to +13 K; the measurement sits in it.

So the bundle's kelvin budget against the donor's 292.26 K:

| term | K |
| --- | ---: |
| `world-f9ig` with `world-jgen`, the cloud optics | **-10.01** |
| `world-trs3`, the derived `gamma` | -1.89 |
| the land water column | +0.17 |
| the energy fixer off | +0.18 |
| `vdiff_lamm`, not resolved from zero | -0.15 |
| `world-5oyp`, the soil heat solver's thermal pair, not resolved from zero | +0.10 |
| **accounted** | **-11.60** |
| **measured drift** | **-12.30** |
| **residual** | **-0.60** |

The residual was -10.8 K before this section and is -0.70 K after it, which is
six per cent of the drift and inside what the remaining unpriced terms --
`OCN-22` at -0.02 to -0.13 K and the arms' own scatter -- can carry between them.

The soil heat solver's row is measured and is the last of the terms that were
carried here as unpriced. Its surface temperature separation does not clear its
own bar, so it enters at what the arm read rather than at zero and is labelled
that way; the term is real and is an AMPLITUDE over land, which a global mean
cannot see. `world-5oyp`, below.

**Do not turn 18.43 W/m2 into kelvin with the static slope.** At
`lib/sensitivity.py`'s 0.615 K per W/m2 it reads as 11.3 K, and the same
arithmetic on the gamma arm gave 5.7 K where the settled arm said 1.89. The
forcing decays as the two states converge, so the arm's own settled mean is the
measurement and this number is not. What 18.43 W/m2 establishes is the SIZE and
the SIGN of the term, which is what the bisect needed: it is the only term in
the bundle large enough to be the drift, and it is one-signed.

The forcing survives the instrument check the first output bin demands. Bin 0
of every orbit carries wind and humidity that do not belong with the other
eleven (`exoplasim/notes/first-output-bin.md`), so the same difference was taken
over bins 1 to 11 alone: +18.63 W/m2 against +18.43 over all twelve, a one per
cent move on a term of eighteen.

The control binary is `t21_l10_p8_production_radmod_pre_jgen_pre_f9ig`. **37 of
its 38 objects are bit-identical to the shipped binary's and only
`radmod.f90.o` differs**, so A3's objection that two binaries confound the term
with the rebuild is bounded here to one translation unit rather than argued
away, and this project's builds are reproducible enough to say so by
measurement.

**The mechanism is the backscatter fraction, not the optical depth.** At
mu0 = 0.5, at the configured co-albedo scale of 1.192:

| LWP, g/m2 | tau2 fit | tau2 table | beta2 fit | beta2 table | R2 fit | R2 table |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | 0.89 | 1.12 | 0.0407 | 0.1020 | 0.066 | 0.174 |
| 50 | 16.27 | 16.78 | 0.0300 | 0.0683 | 0.414 | 0.602 |
| 100 | 30.24 | 31.20 | 0.0256 | 0.0592 | 0.473 | 0.675 |
| 1000 | 145.27 | 145.84 | 0.0160 | 0.0291 | 0.538 | 0.706 |

The optical depths barely move: 30.2 against 31.2 at 100 g/m2. What moves is
`zb2`. The fit `tswr2*sqrt(mu0)/ln(3+0.1*tau)` returns 0.0256 where Table 1(c)
returns 0.0592, a ratio that runs 1.8 to 2.6 across the whole water-path range,
and band-2 layer reflectance rises with it from 0.47 to 0.67 while
transmissivity roughly halves. Band 2 carries 0.6176 of this star's flux.

**The artifact's own per-layer split says the same thing.** Its
`per_layer_band2` runs to -23.7 W/m2 on the layers carrying the most cloud water
where `per_layer_band1` reaches only -1.5, so more than nine tenths of the term
is in band 2. Band 1 is where the optical depth changed most and band 2 is where
the backscatter table landed, which is the same conclusion the layer optics above
reach from the coefficients rather than from the fluxes.

The pair already on disk agrees from the other side. `run_d35b554cab2a` against
`run_57cecaa8391b` is namelist-identical with both binaries pre-f9ig and differs
only in world-jgen's tau relation: 0.13 K. The tau relation is worth tenths of a
kelvin; the tables are worth the rest.

**The implementation is not defective, and that is the point.** The checks with
right answers all pass, and most of them were already written:
`stephens_tables_vs_fits.py` compares `radmod`'s tables against the papers entry
for entry -- 324 entries agree -- and shows the conservative co-albedo floor
reproducing Stephens Eq. (1) to 1.4e-5. Two more checked here:

- **The granularity is consistent to about a tenth.** Stephens tabulates against
  the CLOUD's normal optical depth and `swr` evaluates per LAYER, then adds. A
  cloud of tau_N = 30 split across 5 layers and added returns a band-2 column
  reflectance of 0.603 against 0.670 read once at tau_N = 30, a ratio of 0.90
  and in the CONSERVATIVE direction. The per-layer application is not
  manufacturing reflectance, which is the defect it was checked for.
- **`ww` reaches `fluxini` on every thread.** `readnl` sets it under
  `mypid == NROOT`; `prolog` broadcasts it at plasim.f90:387, before `fluxini`
  at 537. So the rotation-derived mixing length is not a division by a
  zero-initialised threadprivate on the seven non-root threads of a `p8` run.

So the code does what it was written to do. A tuned coefficient was carrying
about half the published backscatter and the model's previous climate rested on
it. Removing the tuning is the no-tuned-values convention working as intended,
and the cooling it exposes is information rather than grounds for putting the
tuning back.

**What this changes about how a bundle is scored.** Every registered prediction
in this note cites the artifact it came from. The one that missed by an order of
magnitude cited an artifact that prices a different term, and nothing compared
the prediction's cited source against the set of changes the prediction claimed
to cover. A prediction covering two commits needs a bracket for each of them.

## `vdiff_lamm`: a sixth bundle term, real and small

`44050c1e` made the asymptotic mixing length rotation-derived, so `fluxini`
computes `160 * (OMEGA_EARTH/ww)` = 200.5 m against the 160 m ECHAM anchors at
Earth's: +25.3% in a first-order boundary-layer parameter, over land and ocean
alike. No registered prediction covered it. `world-awm5` closed as "no numeric
and no form change", which is true of `vdiff_c` and not of the mixing-length
half delivered in the same commit.

**The prediction for it is RETROSPECTIVE and is marked so here**, because a
prediction written after the change shipped is not the same object as one
written before it and treating the two alike is worse than having neither. What
was fixed before the arm ran, and only that: the term is MATERIAL if the settled
pair separates by more than 1.0 K, and no sign was predicted, because a longer
asymptotic length raises free-tropospheric exchange and lowers the near-surface
gradient and which way that moves the global mean depends on the cloud and ice
response. The zero hypothesis under test was world-awm5's closing claim, "no
numeric and no form change".

It has a route now. `VDIFF_LAMM` reaches `fluxmod_nl` from
`model.asymptotic_mixing_length_m`, written unconditionally by
`configure_otherargs` and checked by `expected_namelist_keys`, so a control arm
is a namelist key rather than a code fork. The model says which branch it took:
`run_5373310a7b9f`'s `plasim_diag` prints `asymptotic mixing length (m)
160.0000` against the control's `200.5476`.

AND THE NUMBER IS ON THE MANIFEST, which is what a bundle needs rather than a
diag a reader has to know to open. `vdiff_lamm`, `gamma` and `rcritwidth` are
the three keys the model reads on their SIGN, deriving when it is negative, so
for those three the namelist records the selection and not the value: two runs
at different rotation rates or different rungs stage the same `-1` and
integrated different constants. `config/planet.yaml` declares the selection by
name, `derived`, and `run_exoplasim.py` reads the model's own initialisation
print back into `derived_model_constants` on the run manifest.
`check_consistency.py` refuses a run that staged a sentinel and recorded no
value, and refuses a declared literal the model did not echo back. world-et25.

**Measured, on the declared window.** `run_5373310a7b9f` against
`run_c9c24d438a94`, both 25 orbits, difference over orbits 15 to 24:
**-0.153 +/- 0.062 K**, with the standard error from `lib/autocorrelation.py`
over the paired difference. The criterion fixed before the run was
`|diff| > 2*sqrt(2)*max(SEM)`, which is 0.176 K here, so the term is **NOT
RESOLVED FROM ZERO** at this length and is reported as such rather than as a
result. The materiality threshold registered with it was 1.0 K and the arm is an
order below that either way.

Its forcing is **-0.05 W/m2** in top-of-atmosphere shortwave and -0.29 W/m2 net,
against the cloud optics' +18.43. The mixing length moves far less than its
asymptote does, because `zmixm = lambda*k*z/(lambda + k*z)` saturates: at the
lowest model level, about 331 m, it is 79.7 m against 72.4 m, a tenth, rising to
a quarter high in the column where `k*z` dominates.

## Every run directory keeps the binary it integrated with

The bisect above was declared unfinishable because "no pre-batch-2 binary was
preserved -- the only saved set is `prerebuild_binaries_2026-08-18`". That is
not so, and the refutation is one `ls`. `install_run_executable` copies the
executable into the run directory, so `run_14906cb7b914` still holds
`most_plasim_t21_l10_p16.x` at sha `b81efd05675b67d8`, which is the sha its
`INDEX.json` entry names. Every historical run carries its own control.

That binary settles what the drift window contains, from the executable rather
than from a namelist: `strings` finds `tswr3` in it three times and `cloudabs`
not at all, so the donor was integrated PRE-world-f9ig and the cloud optics
rewrite is inside the window even though c8debdc1 is an ancestor of the commit
batch 2 starts at.

A preserved binary is not drop-in runnable against today's staging -- the donor's
would refuse `CLOUDABS` and want `TSWR3` -- so it is a control for a matching
namelist rather than for the current one. The source-built control used above is
the cleaner instrument for a single term, because it differs in one object file
instead of in every change between two dates. What is wrong is only the claim
that nothing was kept.

## world-5oyp: a seventh bundle term, the soil heat solver's moisture dependence

`3aecf4ec` made the land soil heat solver read the soil water the model already
carries, interpolating conductivity and heat capacity between a dry and a
saturated endpoint instead of holding one pair. `landmod.f90` declares them:

    soildifdry = 0.2088   soildifsat = 1.4332   W/m/K
    soilcapdry = 1.1111E6 soilcapsat = 2.9689E6 J/m3/K

A factor of 6.9 in conductivity and 2.7 in capacity across the wetness range, on
every simulated land cell, and it is inside the batch-2 window: the donor
predates it. It appeared in no registered prediction, which is the second time
the sum has been checked against a total missing a live term.

**All six keys are already in `landmod_nl`**, so the model can be told to run the
old behaviour -- setting a dry endpoint equal to its saturated partner makes the
column constant and the solver bitwise what it was. The arms above staged none of
them, so every run between `3aecf4ec` and the route below integrated the compiled
pair with no artifact saying which pair that was.

It has a route now. All six reach `landmod_nl` from `surface.soil_thermal`,
written unconditionally by `configure_otherargs` at the compiled defaults and
checked by `expected_namelist_keys`, so the control arm is four namelist keys
rather than a code fork. `run_exoplasim.py` also checks the model's compiled
endpoints against
`pedology/config/land_column_properties.yaml:thermal`, which is where they are
DERIVED, at the contract's own tolerances: the six numbers exist twice and a
stale copy of a derived endpoint is a number with no live derivation.

### The registered prediction, before any arm

The magnitude comes from `notes/audits/soil-thermal-inertia.md` and from the
saturation reach declared in
`pedology/config/land_column_properties.yaml:thermal.saturation_mapping`, which
is the artifact that prices this term and nothing else. Thermal inertia in
J/m2/K/s^0.5: the retired pair is `sqrt(1.8 * 2.4e6)` = 2078, and the interval
the store can reach on this build's median column is 1342 at an empty store and
1834 at a full one. So the term is a one-signed REDUCTION in soil thermal
inertia of 12 to 35 per cent on every land cell, largest where the column is
driest.

**TOP-OF-ATMOSPHERE FORCING: ZERO, EXACTLY, AND THAT IS THE ENTRY IN THE SUM.**
The soil column is closed at its base, so its annual-mean heat flux is zero
whatever its heat capacity and conductivity are: this is a change to how a
closed reservoir stores and conducts heat, not a change to any flux crossing the
top of the atmosphere at fixed state. It enters the bundle sum at 0.00 W/m2 with
that reason, which is a different object from being absent from the list. Every
kelvin it can be worth arrives through the response to the amplitude below.

**THE DIRECT TERM, AND IT IS AN AMPLITUDE.** For a periodic surface flux on a
semi-infinite column the surface temperature amplitude goes as the inverse of
the thermal inertia, and the seasonal skin depth here, `sqrt(2*kappa/omega)` at
`kappa` near 5e-7 m2/s, is about 2 m against the column's 12.4 m, so
semi-infinite is the right limit for both the diurnal and the seasonal wave.
Inverting the inertia reduction:

- **The land seasonal peak-to-trough range of `ts` rises by 13 to 55 per cent
  per cell**, +13 at a full store and +55 at an empty one, and the land-area
  mean of that range rises by a value inside that bracket. Read off the monthly
  `ts` of the two arms over the same orbits.
- **The rise is largest where the column is driest**, because the store's own
  reach is what sets it. The per-cell change must correlate positively with the
  cell's soil water deficit.
- **The soil's thermal diffusivity falls by about a third** -- 7.5e-7 m2/s to
  5.1e-7 -- almost independently of the store, because the capacity and the
  conductivity fall together. The seasonal penetration depth therefore falls by
  about 18 per cent and the phase lag of the deep soil layers shortens with it.

**THE GLOBAL MEAN: NOT RESOLVED, AND NEGATIVE IF IT IS.** No magnitude is
registered for it, because none of the artifacts in this tree prices the
feedback that carries it: with no direct forcing, the mean can move only through
terms nonlinear in the amplitude, and the one that is one-signed is the snow
line, which a larger seasonal range pushes equatorward in winter. The zero
hypothesis under test is that the settled pair does not separate. The term is
MATERIAL if it separates by more than 1.0 K, which is the threshold `vdiff_lamm`
was registered against and is fixed here before the arm for the same reason.

**THE RESOLUTION BAR, FIXED IN ADVANCE.** A paired 25-orbit set from one
restart, one binary, differing only in the four `landmod_nl` endpoints;
difference taken over orbits 15 to 24, standard error from
`lib/autocorrelation.py` over the paired difference. The difference is RESOLVED
only where `|diff| > 2*sqrt(2)*max(SEM)`, the same bar `vdiff_lamm` was reported
under. Anything below it is reported as not resolved at this length and not as a
result.

**What would mean wrong:**

- The land seasonal range FALLS, anywhere. Lowering the thermal inertia cannot
  damp the surface more, and a cell that damps more is an implementation fault
  rather than a small result.
- The land-mean seasonal range moves by less than 13 per cent or more than 55
  per cent. Below is the store not reaching the saturation interval the contract
  declares, and the first thing to check is the arm's own `SOILSRWP` and
  `SOILSRFC`; above is the semi-infinite limit being the wrong one, which would
  mean the column's 12.4 m base is being felt at the seasonal period.
- The per-cell amplitude change does not correlate with soil wetness. That says
  `soilwtherm` is not reading the water at all, and the control arm would then
  agree with its own arm by construction.
- A resolved WARMING in the global mean. The one one-signed feedback runs the
  other way.
- A resolved global-mean separation above 1.0 K. The term is then material, the
  zero-forcing argument is not what the model is doing, and it must be bisected
  rather than folded into the bundle.

### The bar sized against the instrument, and the pair the arm is bought as

Written 2026-08-27, BEFORE either arm was launched and against
`run_c9c24d438a94`, which is the bundle's control and is not either arm.
`exoplasim/analysis/arms/land_seasonal_range.py` reports its land-area-weighted
seasonal `ts` range at **40.37 K over orbits 15 to 24, with an orbit-to-orbit
scatter of 0.433 K and an autocorrelation-corrected standard error of 0.137 K at
tau = 1.0 orbits**. The registered +13 to +55 per cent is worth 4.5 to 19 K on a
constant-column range of that order, so the effect is twelve to fifty times the
0.39 K bar that standard error implies, and the amplitude half of this
prediction is comfortably inside what the instrument resolves. The global-mean
half is registered as unresolved and nothing here changes that.

The same series carries the settling. That run's land-mean range rises from
31.99 K in its first orbit to 40.03 K by orbit 15 and is flat within its own
scatter after it, so orbits 0 to 14 are the settling block and 15 to 24 is the
window, which is the split the prediction already declared. The rise is itself
the term: `run_c9c24d438a94` branched from a donor integrated BEFORE `3aecf4ec`,
so it is an interpolating arm relaxing away from a constant-column initial
state, with no control beside it.

**BOTH ARMS ARE RUN FRESH, and that is a departure from measuring against
`run_c9c24d438a94` directly.** A3's first condition is one binary, and that
control was integrated at executable sha `37026652`, which predates both this
term's namelist route and PHYS-15's six `landmod_nl` keys; today's
`most_plasim_t21_l10_p8.x` is `98b423ce` and the two cannot be one binary.
Reusing it would confound the soil thermal pair with every source and config
change between the two dates, which is the confound the paired form exists to
remove. So `jsfm_interp` and `jsfm_const` are both new, on one binary, from
`run_14906cb7b914`'s `MOST_REST.00034` -- the restart the whole bundle branched
from -- and differing in `SOILDIFDRY`, `SOILDIFSAT`, `SOILCAPDRY` and
`SOILCAPSAT` and in nothing else. The consequence for the bundle is stated
rather than hidden: this term is measured at a later operating point than the
five terms above it, and a sum that adds them is adding two operating points.

Both arms also carry `land_albedo_source: lithology` where the base config
declares `vegetated`, because that is the mode the albedo staged in this tree
was built in and `refuse_an_albedo_field_the_config_misdeclares` stops a run
whose declaration and staged bytes disagree. Neither arm reads those fields:
`--superseded-surface-ok` discards the staged surface and `landini` takes all
three albedo bands from the restart, so the key decides what is copied into the
run directory and what the manifest records, and it is identical on both sides.
The pair still differs by four `landmod_nl` keys and nothing else, which
`expected_namelist_keys` is checked on before either arm is launched.

One further prediction follows from the settling above and is registered with
the rest: **`jsfm_const` starts at its own equilibrium and `jsfm_interp` does
not.** The donor ran the retired constant pair, so the constant arm's land-mean
range should be flat from its first orbit while the interpolating arm climbs the
way `run_c9c24d438a94` did. A constant arm that also climbs would mean the
climb is not this term's.

### The longer pair, and its window declared before it ran

Written 2026-08-27, after a first 25-orbit pair and BEFORE the 50-orbit one. The
registered window's number is not revised by anything here; it is re-read on the
longer arms and reported beside them.

The first pair's separation was still trending across the registered window --
+0.093 +/- 0.043 K per orbit -- and neither arm passed its own diagnostic
settling assessment, so 25 orbits gave a LOWER BOUND rather than a separation.
That is a length problem: the soil temperature column reaches 12.4 m at a
diffusivity near 3e-7 m2/s, so its own relaxation is of order thirty orbits and a
25-orbit arm cannot have finished it.

**The window declared before the longer pair ran is orbits 40 to 49**, on the bar
already registered: RESOLVED only where `|diff| > 2*sqrt(2)*max(SEM)` with the
standard error from `lib/autocorrelation.py` over the paired difference, plus the
same slope test, and a separation still trending at 50 orbits is reported as a
lower bound again rather than extrapolated to an asymptote.

**The pair was RE-RUN at 50 orbits rather than extended, and the reason is a
defect worth naming.** `run_exoplasim.py --ncpus 8` writes the thread count into
the loaded configuration before the manifest is stamped, so the manifest carries
8 while the configuration FILE still carries what it declares;
`continue_exoplasim.py` compares the two and refuses `model.ncpus: 8 -> 16`, and
has no flag of its own. A run prepared that way cannot be extended at all.
world-q4gh. Both arm configs now DECLARE `ncpus: 8`, which is what a pair on one
host runs at, and the 25-orbit pair was retired.

That re-run bought a check that could have failed and did not: on the longer
arms the registered window reproduces the retired pair EXACTLY -- 38.170 K
against 41.206 K for a difference of +3.036 +/- 0.146 K -- which is the
determinism this model is supposed to have at a fixed thread count, tested rather
than assumed.

## Measured: the world-5oyp soil thermal pair, 2026-08-27

`run_ae5aaf479f57` (the shipped interpolation) against `run_8f75d23c0730` (the
retired constant column at 1.8 W/m/K and 2.4e6 J/m3/K), both 50 orbits at T21 on
`most_plasim_t21_l10_p8.x` sha `98b423ce`, both branched from
`run_14906cb7b914`'s `MOST_REST.00034`, each one `diagnostic` segment at low I/O
so the whole arm is one instrument and no window can span a regime change. Their
namelists differ in `SOILDIFDRY`, `SOILDIFSAT`, `SOILCAPDRY` and `SOILCAPSAT` and
in no other key, which `expected_namelist_keys` was checked on before either arm
was launched.

**THE AMPLITUDE IS RESOLVED, IT IS THE SIGN PREDICTED, AND IT HAS SETTLED.** The
land-area-weighted seasonal range of `ts`:

| window | constant | interpolating | difference | resolved at | change |
| --- | ---: | ---: | ---: | ---: | ---: |
| 15-24, as registered | 38.170 K | 41.206 K | **+3.036 +/- 0.146 K** | 7.4 to 1 | +7.96% |
| 40-49, the extension | 38.794 K | 41.557 K | **+2.763 +/- 0.193 K** | 5.1 to 1 | +7.12% |
| 25-49, everything settled | 38.750 K | 41.588 K | **+2.838 +/- 0.103 K** | 9.8 to 1 | +7.32% |

The separation's own slope over orbits 25 to 49 is **-0.010 +/- 0.014 K per
orbit**, which is flat, so the +0.093 +/- 0.043 read over the registered window
alone was the tail of the two arms' joint settling and not this term still
growing. Neither arm passes `assess_convergence.py --assess diagnostic` even at
50 orbits, and both fail the same four criteria: they are relaxing together from
a donor that ran the scalar bucket land water column, which is a far larger
perturbation than the one under test. That is what the paired form is for, and
the difference's own flatness is the evidence that it worked.

**THE REGISTERED MAGNITUDE IS FALSIFIED, LOW, BY A FACTOR OF ABOUT FIVE.** +7.3
per cent against a registered +13 to +55. The cause the prediction named for a
low result is ruled out: both arms staged `SOILSRWP` 0.376 and `SOILSRFC` 0.7642,
the contract's own values, and the arms' median degree of saturation is 0.033,
which is the dry end of the map where the registered effect is LARGEST rather
than smallest.

What is wrong is the inversion. `amplitude ~ 1/inertia` holds for a semi-infinite
column under a PRESCRIBED periodic surface flux, and a modelled surface flux is
not prescribed: a larger swing in `ts` raises the outgoing longwave and the
turbulent fluxes that leave the surface, and those damp the swing that produced
them. A prediction sized on the column alone is an upper bound on what a coupled
surface does, and that is the standing lesson for the next amplitude prediction
in this file.

**THE MOISTURE DEPENDENCE IS DOING WHAT IT WAS WRITTEN TO DO.** The per-cell
fractional change in seasonal range correlates with the cell's soil water DEFICIT
at Pearson +0.45 and Spearman +0.49 over orbits 40 to 49. The same statistic on a
pair differing by one unrelated namelist key -- `run_5373310a7b9f` against
`run_c9c24d438a94`, `world_5oyp_null_pair.json` -- is -0.17, so the correlation
belongs to this term and not to the instrument. That pair is the null for the
other per-cell row too: its land-mean range difference is -0.8 per cent and is
NOT resolved, and 579 of its 1019 land cells fall.

**THE DEEP COLUMN CONFIRMS THE DIFFUSIVITY.** The seasonal amplitude of the
model's third soil temperature layer falls and the fourth falls further, and the
fourth layer's maximum arrives earlier. A lower diffusivity shortens the seasonal
penetration depth and the phase lag together, and both move that way.

**THE TOP-OF-ATMOSPHERE FORCING IS ZERO, MEASURED.** Over the settled window the
net imbalance separates by -0.044 +/- 0.085 W/m2 and over the first three orbits,
where the state has not yet moved, by +0.088 W/m2. Neither clears its bar. **So
the term's entry in the bundle sum stays 0.00 W/m2 and is now measured rather
than argued from a closed column.**

**AND THE REGISTERED SIGN OF THE GLOBAL MEAN IS FALSIFIED.** Over orbits 25 to 49
the near-surface air temperature separates by **+0.240 +/- 0.034 K, resolved at 7
to 1, and it is a WARMING**, against a prediction that named the snow line as the
one one-signed feedback and expected cooling. The surface temperature itself is
+0.097 +/- 0.037 K and does not clear its own bar. The mechanism is visible in
what else resolves: the modelled sea-ice fraction falls by 0.0072 +/- 0.0007,
cloud fraction rises by 0.0039 +/- 0.0007, column water vapour by 0.236 +/- 0.037
and precipitation by 0.012 +/- 0.004 mm/day. A larger land seasonal range gives
warmer summers and colder winters, and what survives into the mean is the pair of
feedbacks that are nonlinear in it; at this state the summer melt outruns the
winter growth and the sign is warming.

**IT IS NOT MATERIAL, AND THAT IS THE VERDICT THAT WAS REGISTERED FOR.** The
threshold fixed before the arm was 1.0 K and the largest resolved global-mean
separation is a quarter of it, so the term does not need bisecting out of the
bundle. What it is worth is an AMPLITUDE over land, which is what the
prediction said it would be.

**One registered falsifier fired and is not an implementation fault.** The
seasonal range falls on 157 of 1019 land cells over the extension window, against
a row that called any fall a fault. On the unrelated pair above 579 of 1019 fall,
so 15 per cent is the tail of a strongly one-signed distribution rather than the
half a null pair gives, and the row is too strict for a coupled surface where a cell's range is
set by its own energy balance and not by its column alone.

**The supplementary prediction registered with the resolution bar is FALSIFIED.**
The constant arm was predicted to start at its own equilibrium because the donor
ran the same constant pair; its land-mean range climbs from 28.80 K to about
38.8 K instead. The donor ran `NLANDWCOL = 0` at one water layer and both arms run
the declared three-layer column, so both relax from the same non-equilibrium
state. The climb belongs to the land water column, not to this term.

Numbers from `exoplasim/analysis/arms/world_5oyp_land_seasonal_range.json` and
`exoplasim/analysis/arms/world_5oyp_arms.json`.

## Measured: the star weight on its new carrier, `CLOUDABS`, 2026-08-26

`run_a61f97a32b45` against the same control `run_c9c24d438a94`, same restart,
same binary, differing in `CLOUDABS` alone: 1.0 against the configured 1.192.
Over orbits 15 to 24, **-1.373 +/- 0.033 K**, resolved at 41 to 1. So this
star's co-albedo correction is worth **+1.37 K of warming** as configured.

**This is NOT part of the drift**, and the reason matters. The donor staged
`TSWR3 = 0.006556`, which is the compiled 0.0055 times the same 1.192, so the
star correction was already in the donor and did not change across the window.
What changed is the CARRIER, and that is the `world-f9ig` term measured above.

**The prediction registered before this arm ran was 1.4 to 2.6 K of cooling, and
the measurement falls just outside it**, by 0.03 K at the near edge. It was
built from PHYS-11's own arms, which measured -2.63 K at scale 0.78 and +2.53 K
at 1.28 for a slope of 10.3 K per unit scale, and 0.192 of that is 2.0 K. The
implied slope here is 1.373/0.192 = **7.15 K per unit scale**, thirty per cent
below it.

**The weight transfers and its kelvin worth does not.** `cloud_band_weight.py`
derives 1.192 as a ratio of flux-weighted co-albedos, and a ratio is carrier
independent, which is why `world-f9ig` was right to move it from `TSWR3` to
`CLOUDABS` unchanged. What is carrier DEPENDENT is the base it multiplies:
Table 1(a)'s co-albedo is not `tswr3*mu0^2*ln(1000/tau)`, so the same 19.2 per
cent lands on a different number and buys a different number of kelvin. A slope
measured on one carrier is not a prediction for the same key on another, and
PHYS-11's 10.3 K per unit scale should be quoted against the fits it was taken
on rather than against the tables.

## PHYS-15, the first of two commits: the modelled soil albedo responds to the surface layer

Written 2026-08-26, BEFORE the arming edits and before any arm. It prices ONE
commit: the one that declares `surface.soil_albedo_moisture`, cuts
`surface.land_water_column` to three layers at 0.02, 0.48 and 1.0 m, routes
`NWETSOIL`, `SKINSRAD`, `SKINSRFC` and the three `WETSIGMA` keys into
`landmod_nl`, and stages codes 1742, 1750 and 1760 beside 174, 175 and 176. The
`dwmax` capacity half is a separate commit with its own entry below, sized from
its own artifact, because a joint entry is exactly the failure A3 records.

**THE ARTIFACT ITS MAGNITUDE COMES FROM.** `analysis/soil_albedo_wetting.json`
derives the saturated endmember per rock class and per band from Lekner and Dorf
(1988) and Twomey, Bohren and Mergenthaler (1986) applied per wavelength, and it
prices this term and nothing else. `build_surface_albedo.py` integrates it to the
land-mean pair it writes into `albedo_report.json:wetting`, quoted in
`exoplasim/notes/soil-albedo-moisture.md` and measured 2026-08-27 at T21 on the
vegetated field with lakes: land mean 0.170772 dry against 0.140616 saturated.
The saturation mapping the model reads is
`pedology/analysis/land_column_properties_report.json`'s
`surface_layer.saturation_endpoints`, 0.0 at an empty layer and 0.7642 at a full
one. Neither is a number this entry chose.

**THE FLUX-TO-KELVIN ROUTE, stated because a surface albedo is not a planetary
one.** `lib/sensitivity.py` deliberately refuses that conversion, so it goes
through the budget that owns the assumption:
`scripts/error_budget.py:albedo_to_kelvin`, land fraction 0.432841 by surface
class from the build manifest, and `DEFAULT_ATTENUATION = 0.5` with the factor of
two that file declares for itself. The planetary albedo is 0.3173, taken from
`run_c9c24d438a94`'s own annual means over orbits 15 to 24 through
`sensitivity.planetary_albedo_from_fluxes`, because `config/planet.yaml`'s
`baseline_climatology` is null and this build has no climatology to read.

**AMENDED 2026-08-27, AND THE AMENDMENT IS VISIBLE BECAUSE THE PREDICTION IS
REGISTERED.** Every kelvin below fell by 21 per cent and no threshold in the
prediction is unaffected. Nothing about this term's own artifacts moved: what
moved is `lib/sensitivity.py`'s slope, from 202.0 to 159.7 K per unit flux ratio,
because the bracket it was measured on was two deleted runs on a superseded
build. The albedo fall, the mixing curve, the land fraction and the attenuation
bracket are unchanged, so the whole amendment is one multiplication by
159.7/202.0 in the kelvin columns and its consequence for the two thresholds.
The ceiling was registered at +9.0 K bracketed +4.5 to +18.0, and the material
threshold at f above 0.0526.

**AMENDED AGAIN 2026-08-27, AND THIS ONE MOVES THE ARTIFACT.** Two things moved
under this entry and neither is the term's physics; both are visible only because
the entry states which artifact each number came from.

- **The staged field became the VEGETATED baseline with lakes composited in.**
  The entry was priced on a bare-rock field whose land mean was 0.261450 dry
  against 0.105609 saturated. The field a run now reads is 0.170772 against
  0.140616. A canopy is not a wetting surface and open water is not one either,
  so `build_surface_albedo.py` gives the covered and lake fractions of a cell the
  same value at both ends and only the bare fraction wets; the endmember swing
  falls from 0.155841 to 0.030156, a factor of 5.2, and the whole factor is
  cover.
- **`skinsrfc` went from 0.7877 to 0.7642** when the saturation mapping was
  derived from the emitted states instead of typed. Every `Sr` column below was
  computed at the old endpoint.

The table below is repriced on the field now staged, per cell through Sadeghi
Eq (13) and then land meaned; `wetting.at_full_surface_layer` in
`albedo_report.json` carries its last row so it does not have to be re-derived
again.

**THE CEILING IS NOW +1.25 K at attenuation 0.5, bracketed +0.62 to +2.49**, and
the material threshold at 1.0 K moves from f above 0.0677 to **f above 0.756**.

**THE CONSEQUENCE IS THAT THE REGISTERED BAR NO LONGER DISCRIMINATES IN `tas`,
and that is the finding rather than a caveat.** The resolution bar fixed below is
two root two times the larger standard error of a 25-orbit paired difference,
which on the held arms is 0.56 to 0.75 K in global-mean `tas` and puts the bar at
1.6 to 2.1 K. The ceiling sits below its own bar, and the ceiling is a state the
cascade cannot hold. So a paired arm reporting "not resolved" in `tas` would be
reporting its own scatter and would say nothing about this term. The `alb`
measurement below is unaffected and is now the test rather than the corroboration:
it is set by the boundary condition and the saturation the cascade hands it, not
by the circulation. What a temperature separation is still good for is the SIGN,
which is hard, and the upper bound: anything above +2.5 K is not this term.

**THE SIGN IS HARD AND IT IS ONE-SIGNED.** Every non-refused region's saturated
field is darker than its dry field, `build_surface_albedo.py` refuses to write a
field where that fails, and `wet_soil_albedo` is monotone in the saturation. So
the armed arm's land-mean surface albedo can only FALL and its absorbed shortwave
over land can only RISE. A resolved cooling is an implementation fault, not a
small result.

**THE CEILING, WHICH IS NOT THE PREDICTION.** A permanently full surface layer
maps to 0.7642 and not to 1, so the staged saturated endmember is an endpoint the
mixing approaches and never reaches. Mixed per cell and then land meaned that is
a land-mean albedo fall of 0.024633 against the 0.030156 endmember swing, worth
+1.25 K at attenuation 0.5 and +0.62 to +2.49 across the declared factor of two.
That is a state the model cannot hold and it is registered as a bound, not an
estimate. It is also far outside the regime `SLOPE_K_PER_FLUX_RATIO` was measured
in, which is itself the finding: this term cannot be settled by prediction.

**WHAT IS PRICED EXACTLY: the term per unit wetness.** `f` is the surface
layer's liquid store as a fraction of its own capacity, land-area and time
meaned; `Sr` is `skinsrad + f*(skinsrfc - skinsrad)`; the albedo is the
Kubelka-Munk mixing at `wetsigma = 1`, taken per cell and then land meaned.

| f | Sr | land-mean albedo | fall | K at 0.25 | K at 0.5 | K at 1.0 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.02 | 0.0153 | 0.170111 | 0.000660 | +0.02 | +0.03 | +0.07 |
| 0.05 | 0.0382 | 0.169141 | 0.001631 | +0.04 | +0.08 | +0.17 |
| 0.10 | 0.0764 | 0.167575 | 0.003197 | +0.08 | +0.16 | +0.32 |
| 0.15 | 0.1146 | 0.166066 | 0.004705 | +0.12 | +0.24 | +0.48 |
| 0.25 | 0.1910 | 0.163203 | 0.007569 | +0.19 | +0.38 | +0.77 |
| 0.50 | 0.3821 | 0.156775 | 0.013997 | +0.35 | +0.71 | +1.42 |
| 1.00 | 0.7642 | 0.146138 | 0.024633 | +0.62 | +1.25 | +2.49 |

The curve is strongly CONCAVE, which is the half of Sadeghi, Jones and Philpot
that matters here: a skin at a tenth of its capacity has already given up an
eighth of the swing, where a linear mix in the albedo would give a tenth.
Evaluating the mixing at the mean `f` rather than meaning it over the
distribution of `f` OVERSTATES the fall, because the fall is concave in the
saturation, so every row is an upper bound at its own `f`.

**THE ONE QUANTITY NO HELD ARTIFACT PRICES, and it is `f`.** Three reasons, and
each is a property of the tree rather than an opinion.

- `config/planet.yaml`'s `baseline_climatology` is null. This build has no
  climatology at all.
- Every run in `exoplasim/runs/INDEX.json` integrated the two-layer column at
  0.5 and 1.0 m with `dwmax` from `awc_mm`. Its `mrso` is the whole column and
  its finest reservoir is 0.5 m, twenty-five times the albedo depth, so no held
  run carries the store this term reads.
- The finest cadence any held run writes is the snapshot stream, 32 records per
  orbit. An orbit is 182.8 Earth days, so those are 4.6 Vesper days apart against
  a skin that empties under evaporation in about one. No held stream resolves the
  store's duty cycle even in principle.

So `f` is not bracketed here, and a plausible-looking central value would be a
guess wearing an artifact's name. `pedology/README.md` records the two-layer
column sitting at a median 15 per cent of capacity, and that is NOT adopted as
`f`: it prices a different reservoir, and the two corrections to it run opposite
ways, the skin filling before anything below it and emptying two orders of
magnitude faster.

**REGISTERED, AND IT IS A THRESHOLD RATHER THAN A VALUE.**

- **Top-of-atmosphere forcing: POSITIVE, and its entry in the bundle sum is
  `+1.25 K times f_eff` with `f_eff` unmeasured.** The bundle sum cannot be closed
  on this term until the arm measures it, and that is stated rather than papered
  over with a mid-range number.
- **The term reaches the 1.0 K threshold, the one `vdiff_lamm` and `world-5oyp`
  were registered against, only once `f` exceeds 0.756** -- equivalently once the
  armed arm's land-mean `alb` over snow-free land falls by more than 0.0198,
  which on this field is nearly the whole ceiling. The zero hypothesis under test
  is that the pair does not separate in `tas`, and the prediction is now that it
  HOLDS: the term is real and one-signed and its temperature consequence is
  inside the instrument's own scatter. What the arm is bought for is the `alb`
  fall and the `f` that produced it, not a temperature.
- **The realised fall is measured directly and needs no saturation diagnostic.**
  Both arms write `alb`; its land mean over the same orbits IS the fall, and the
  row of the table it lands on converts it. That measurement is what turns this
  entry from a threshold into a number.
- **It acts only where the modelled surface is snow-free, ice-free land.** Snow
  and glacier albedo override the background pair after `getalb` has mixed it, so
  the realised area is below the land fraction, one-signed and downward.

**THE RESOLUTION BAR, FIXED IN ADVANCE.** A paired 25-orbit set from one restart,
one binary, differing only in `NWETSOIL`; difference over orbits 15 to 24,
standard error from `lib/autocorrelation.py` over the paired difference, resolved
only where `|diff| > 2*sqrt(2)*max(SEM)`. The same bar `vdiff_lamm` and
`world-5oyp` were reported under.

**What would mean wrong:**

- A resolved COOLING, or a land-mean `alb` that RISES. The staged pair cannot
  brighten any region and the mixing is monotone, so either is an implementation
  fault. The first thing to check is that the arm read codes 1742, 1750 and 1760
  rather than the sentinel, which `landini` prints.
- A land-mean `alb` fall above 0.024633. That is the permanently-saturated
  ceiling and the surface layer caps at field capacity, so exceeding it means the
  fill fraction is not being clipped to one or `dsoilwfc` is not the capacity the
  cascade fills.
- `evaporite` cells moving at all. The class is staged wet-equal-to-dry, so its
  albedo must be bitwise identical between the arms at every saturation. A move
  there means the refusal was not written into the field.
- The measured fall not tracking the table at the arm's own measured `f`. The
  mixing is checked bitwise against Sadeghi Eq (13) over 400 random quadruples in
  `analysis/soil_albedo_wetting.json:model_mixing`, so a disagreement is in the
  saturation the model hands it, not in the curve.
- A separation resolved above the ceiling row, +1.25 K at attenuation 0.5. The
  attenuation is declared to a factor of two and the ceiling is a bound on the
  albedo, so anything above +2.5 K is not this term.

## PHYS-15's arm: what the paired `NWETSOIL` set will measure

Written 2026-08-27, BEFORE the arms were queued and before either had run. It
prices ONE experiment: `run_893e276ee029/MOST_REST.00209` seeded into two runs on
one binary, twenty-five orbits each, differing in `NWETSOIL` alone.

**THIS ONE IS NOT A THRESHOLD, because `f` is no longer unmeasured.** The two
entries above registered `+1.25 K times f_eff` and said no held artifact prices
`f`. That was true of the output streams and false of the RESTART:
`landmod.f90:1553` writes `dwatcl` with `mpputgp`, so every one of this run's 210
per-orbit restarts carries the surface layer's liquid store, and `dwmax` and
`dsoilwfc` are carried beside it. `dwatcl(:,1)` over `dwmax * dsoilwfc(:,1)` is
`wetalb`'s own arithmetic, so this is the model's `f` and not a reconstruction of
it. Measured over the settled block, restarts 180 to 209:

| quantity | value |
| --- | --- |
| `f`, land-area mean | 0.3283, sd 0.0146 across the thirty |
| `f`, land median | 0.0096 |
| `f`, land p90 | 0.9997 |
| share of land above `f` = 0.5 | 0.319 |

**`f` IS BIMODAL, and that is the finding rather than a detail.** A third of the
land carries a saturated skin and most of the rest is at air dry; there is very
little in between. So the mixing must be evaluated per cell and meaned, never at
the mean `f`, and the gap between the two is at its widest here: the fall at
`f` = 0.3283 is 0.0098 and the mean of the per-cell falls is 0.006745, a factor
of 1.45. Every kelvin below is on the second.

**REGISTERED: the land-mean `alb` difference will fall between 0.00378 and
0.00675, and the arms will not separate in `tas`.** The two ends are the same
per-cell calculation masked and unmasked by what covers the ground. The
background fall is 0.006745; on 0.7673 of the land, area weighted, `dalb` equals
the mixed background to 2e-4 and the term reaches the surface undiminished,
which puts the masked figure at 0.003780. Snow and glacier ice override the
background after `getalb` has mixed it, so the truth is at or below the unmasked
end, and the mask is a blend rather than a switch, so it is at or above the
masked end.

| | land-mean `alb` fall | K at 0.25 | K at 0.5 | K at 1.0 |
| --- | --- | --- | --- | --- |
| unmasked, the upper end | 0.006745 | +0.17 | +0.34 | +0.68 |
| snow and ice masked, the lower end | 0.003780 | +0.10 | +0.19 | +0.38 |

**THE INSTRUMENT, AND WHY THE TEST IS `alb` AND NOT `tas`.** A twenty-five-orbit
paired difference in global-mean `tas` carries a standard error of 0.56 to 0.75 K
on the arms this project has run, so the bar of two root two times the larger is
1.6 to 2.1 K. The prediction above is +0.19 to +0.34 K. **The pair cannot
separate in temperature and a null there means nothing**, which is fixed in
advance rather than discovered afterwards. The `alb` difference has no such
problem: its own inter-orbit scatter on the donor is 0.00039 against a signal of
0.0038, and pairing removes the common mode, so it is resolved many times over.

**What would mean wrong:**

- A land-mean `alb` difference outside 0.00378 to 0.00675. Both ends come from
  the donor's own restarts and the staged pair, so a miss means `f` moved when
  the arms diverged -- which is a real effect, the brighter arm being cooler and
  its skin wetter, and is what the width is for. Outside it by more than the
  0.00039 inter-orbit scatter is a disagreement about the mechanism.
- Any `alb` difference of the other sign, anywhere on land. The mixing is
  monotone and every non-refused region's saturated field is darker.
- `evaporite` cells differing between the arms at all. The class is staged wet
  equal to dry.
- A `tas` separation resolved above the bar. At +0.34 K predicted and a bar of
  1.6 K, a resolved temperature difference is not this term.
- The `NWETSOIL = 0` arm's land-mean `alb` differing from the staged dry field's
  0.170772 by more than the snow and ice it carries. That arm has no mixing at
  all, so its background IS the staged field.

**THE COST, stated because it is a lock.** Two twenty-five-orbit T21 runs at
`dt` 45, eight threads each, pinned to their own dies and run concurrently under
one host lock: about six minutes of wall time and twelve of model time.

### MEASURED 2026-08-27, and one half of it was wrong

`run_a1c35075747c` at `NWETSOIL = 1` against `run_598eb57c5a34` at 0, one binary
and one seed, orbits 15 to 24.

**The `alb` half is right.** The land-mean fall is 0.006397, inside the
registered 0.003780 to 0.006745 and near its upper end, resolved at 5.7 times the
arms' own bar. The wet arm's background fall, recovered from its own `dwatcl`,
is 0.006709 against the 0.006745 predicted from the donor's: half a percent.
Every per-cell claim holds -- no land cell's background brightens, worst 5.3e-16,
and the twenty-three cells staged wet equal to dry move by 4.2e-16.

The bracket was too wide at the bottom, and the reason is a fault in this entry's
reasoning rather than in the model. The lower end assumed snow and glacier ice
take the term entirely where they hide the background, costing 44 per cent. The
measured cost is 4.6 per cent, because a snow-covered cell has a frozen skin, ice
comes off the layer's capacity in the cascade, and the model already reads it as
dry -- so it was contributing nothing for the mask to take. The mask and the
wetting are ANTI-CORRELATED, and treating the mask as a switch applied to a mean
cannot see that.

**The `tas` half is wrong, and it is the more useful half.** This entry
registered that the arms could not separate in temperature and that a resolved
difference would not be this term. They separated: +0.1972 K, with `ts` agreeing
at +0.2013, resolved at 2.3 times the bar.

The bar was imported and should not have been. It came from arms seeded from a
thirty-seven-orbit control still relaxing, whose paired `tas` standard error over
twenty-five orbits is 0.56 to 0.75 K. These arms are seeded from a 210-orbit
equilibrated state and their standard errors are 0.017 and 0.031 K, twenty times
smaller at the same run length. **A paired arm's power is set by how settled its
DONOR is at least as much as by how long the arms run**, so a bar measured in one
regime is not a bar in another. Every A3 entry in this file that quotes the 1.6
to 2.1 K figure against a well-settled donor is understating its instrument the
same way.

**What the separation buys.** A fall of 0.006397 reaches +0.1972 K only at an
attenuation of 0.3045, against the 0.5 `scripts/error_budget.py` declares while
saying the honest claim is a factor of two. That is the first MEASUREMENT this
project has of that factor, it lands inside the declared bracket at the low end,
and it puts this term's forcing at 0.287 W/m2 rather than 0.471. One term is not
the budget and the attenuation is a property of the atmosphere above the surface
rather than of the surface, so it is one point; world-ckbt carries what would
make it more.

## PHYS-15, the second of two commits: `dwmax` from `evaporable_mm`

Written 2026-08-26, before the edit. It prices ONE commit: the one that makes
`build_surface_soil_water.py` install the column capacity from the states file's
`evaporable_mm` column instead of `awc_mm`, and cut the surface layer's share of
the capacity split from AIR DRY instead of from the wilting point. It is a
HYDROLOGICAL change and shares no artifact with the albedo entry above.

**THE ARTIFACT ITS MAGNITUDE COMES FROM.**
`pedology/analysis/land_column_properties_report.json`'s `surface_layer` block,
which prices this and nothing else: the surface layer's capacity, the sub-wilting
increment it adds, and `increment_over_awc` per cell. Measured on this build,
the increment over `awc_mm` runs 0.0137 at p10 to 0.0747 at p90, median 0.0225,
mean 0.0340. The identity `evaporable_mm - awc_mm` equals the surface layer's
increment cell by cell to 2.8e-14 mm, so the column gains exactly the water the
surface layer can reach and nothing else.

**TOP-OF-ATMOSPHERE FORCING: ZERO, EXACTLY, AND THAT IS THE ENTRY IN THE SUM.**
`dwmax` is a capacity, not a radiative quantity; at fixed state nothing about
this change alters a flux crossing the top of the atmosphere. Every kelvin it can
be worth arrives through the hydrological response below. This is the same object
as `world-5oyp`'s zero and a different object from being absent from the list.

**THE LAYER CUT ITSELF ENTERS AT ZERO, AND IT IS DRIVEN RATHER THAN ASSERTED.**
Splitting the top layer while holding the column's total capacity fixed cannot
move the cascade: what enters is the surface flux and what leaves is base
overflow, and neither depends on where the internal boundaries sit.
`land_column_properties.py:check_split_invariance` drives a two-layer and a
three-layer column through one flux sequence and gets agreement to 1.1e-13 mm on
both the total and the runoff over 4096 steps, against a 1e-9 mm tolerance fixed
before the first run. The soil heat solver is untouched by the cut for a separate
reason: `soilwtherm` gives each temperature layer the saturation of whichever
water layer contains its MIDPOINT, the first temperature layer spans 0 to 0.4 m,
and 0.2 m falls in the 0.02-to-0.5 m layer, so no temperature layer ever reads
the surface layer or its air-dry mapping.

**THE DIRECT TERM, AND IT IS A CAPACITY.** The column's capacity rises by 1.37 to
7.47 per cent per cell, median 2.25, on every land cell that carries a soil. It
cannot fall. Two consequences follow, both one-signed:

- **Land runoff falls.** ExoPlaSim's runoff is the overflow of this bucket, so a
  deeper bucket overflows less. `pedology/README.md` records the one measured
  point on that response: an offline bucket validated against the model's own
  soil water, 0.500 m giving 149.96 mm/Earth-yr and 0.133 m giving 158.46, a
  chord slope of -0.0416 in log capacity. Applied to the capacity bracket above,
  **land runoff falls by 0.06 to 0.30 per cent, central 0.09.** The chord spans a
  factor of 3.75 and the local slope at the operating point is not measured, so
  treat the magnitude as bounded by "well under one per cent" rather than as
  three digits; the sign is not in question.
- **Land evaporation falls slightly.** The wetness factor the limiter builds is a
  function of the store over `dwmax`, so a deeper bucket reaches any given
  wetness on more water. Same fractional scale as the capacity, same sign on
  every cell, and no artifact in this tree prices its own response, so no
  magnitude is registered for it.

**WHAT THIS MUST NOT DO, and it is checkable without a run.** The biosphere's
water does not move, and the reason is stronger than "the other column is left
alone". `biosphere/scripts/build_lpj_driver.py` reads `b`, `theta_s`, `theta_fc`
and `theta_wp` and the per-layer usable shares, and no capacity column at all:
`vesperinput.cpp` scales each layer's own available water, wilting point and
saturation by the share and rescales the aggregates to match. So this commit
touches nothing the driver reads. The check is that the driver is BYTE-IDENTICAL
across it, and a difference means the two consumers have been conflated.

**THE ZERO HYPOTHESIS AND THE MATERIALITY THRESHOLD.** The predicted global-mean
separation is ZERO, and the term is MATERIAL if a settled pair separates by more
than 1.0 K, the threshold fixed here before the arm for the reason `vdiff_lamm`
and `world-5oyp` fixed theirs. A 0.1 per cent change in land runoff is orders
below the carve criterion's own resolution, and the offline bucket experiment
already settled that this world's land runoff ratio is a property of the climate
and not of the bucket.

**What would mean wrong:**

- Land runoff RISES. A deeper bucket cannot overflow more.
- The land-mean `dwmax` does not rise by the report's own increment. The field is
  a read of `evaporable_mm` and nothing here derives a capacity, so a mismatch is
  a column read by the wrong name.
- The capacity split's three shares do not sum to one, or the three per-cell
  capacities do not sum to `evaporable_mm`. The split and the capacity are one
  arithmetic and the builder checks it cell by cell.
- The LPJ-GUESS driver moves. That is the conflation this commit exists to avoid
  and it is a failure whatever the climate does.
- A resolved global-mean separation above 1.0 K. The zero-forcing argument is
  then not what the model is doing and the term must be bisected rather than
  folded into the bundle.
