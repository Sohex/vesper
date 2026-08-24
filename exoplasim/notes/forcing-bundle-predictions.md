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

**These are TRUE diffusivities, and the arms measured below did not run at
them.** `predict_ocean_terms.py` has always built its operator on this planet's
radius, so this table is the prediction at the coefficient named in its own
first column; `hdiffo` did not, until `world-mll`, and divided by a compiled
Earth radius instead. Arms that ran before that fix realised
(PLARAD/6.371E6)^2 = 1.4401 times the coefficient their namelist declared. The
operator is linear in `hdiffk`, so the CLIM-16 bracket, which realised
432/1440/4320, predicts 1.08/3.61/10.84 W/m2 rms rather than the three rows
above. `hdiffk` means what it says from `world-mll` on, and a new arm set to
1000 realises 1000.

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

## PHYS-11: the cloud shortwave constants, bracketed for arms

Seven `radmod_nl` keys tune the cloud optics and none has ever been re-weighted
for this star: `tswr1` (cloud albedo, range 1), `tswr2` (backscatter, range 2),
`tswr3` (single-scattering albedo, range 2 -- the absorption-like one),
`acllwr` (thermal-band mass absorption), and the `rcl1`/`rcl2`/`acl2` triplets
by cloud level. The band PARTITION is star-aware through `zsolar1`/`zsolar2`;
what is Earth's is the physics inside each range.

**The sign is not derivable by inspection, which is the argument for arms
rather than arithmetic.** Two effects oppose. This star puts 0.618 of its flux
in range 2 against the Sun's 0.483, a naive 1.28 on range-2 cloud absorption;
but it concentrates that energy nearer 0.8-1.5 um, where liquid water absorbs
less than in the 2-3 um bands the Sun's range-2 tail reaches. So the honest
bracket on the absorption-like keys (`tswr3`, `acl2`) is 0.78 to 1.28, spanning
both signs of the correction.

**Predicted magnitude bound: 0 +/- 0.6 K.** Clouds are booked at 12.0 W/m2 of
Earth's shortwave absorption beside the water vapour term; +/-28% of that is
+/-3.4 W/m2 atmospheric, and at the same atmospheric-to-TOA ratio the water
vapour calibration row exhibits (0.19), +/-0.65 W/m2 TOA, +/-0.54 K. After
PHYS-9 this is the softest entry in the bundle and the second place to bisect.

`acllwr` is in the task's list because it is untraced Earth tuning, not because
the star moves it: it is a thermal-band constant and has no stellar dependence,
so it gets no arm on this argument. The scattering-side keys are expected near
unity and are held unless the absorption arms surprise.

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
  that flux times the planetary albedo, 0.30 to 0.46 W/m2, which is **-0.25 to
  -0.39 K** on `lib/sensitivity.py`'s conversion at albedo 0.30.
- **The redistribution gain, which is why this is worth running.** `PTOP` is
  5000 Pa and the profile's `bo3` is 20 km, so the ozone maximum sits AT the
  model top and `dqo3(:,1)` carries everything above it. The absorption removed
  is therefore heating of the topmost layer, which radiates efficiently to space
  and is weakly coupled to the surface. Moving it into the troposphere and the
  surface makes it count for surface temperature where it largely did not. The
  ceiling, if all of it landed and were retained, is **+0.84 to +1.29 K**.

**Prediction: warming, +0.2 to +0.9 K in the global mean.** The sign is assigned
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
and `omega C` is **85.6 W/m2/K**, against a radiative damping of **1.18
W/m2/K** taken from `lib/sensitivity.py`'s canonical slope at a planetary
albedo of 0.30. The ratio is 73: this slab is
deep in the inertia-dominated regime, which the config's own note reaches
qualitatively when it says this world's half-length year damps seasonality about
twice as hard as Earth's ocean does.

Two consequences follow with no free parameters:

| depth | seasonal amplitude, relative to 50 m |
| ---: | ---: |
| 25 m | 2.00x |
| 50 m | 1.00x |
| 100 m | 0.50x |

- **Amplitude goes as 1/depth to within 1.5%**, because the damping term is 69
  times smaller than the inertia term and enters only in quadrature.
- **The phase lag is 89.2 degrees**, essentially a quarter of the orbit, and it
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

## The bundle, summed

A3's rule: check the bundle's total against the sum of the per-term
predictions; bisect only on disagreement, and only into the terms whose
predictions were soft.

| term | change | prediction, global mean | source |
| --- | --- | ---: | --- |
| PHYS-9 | `h2oswl` 1.127 | +1.25 K (+1.4 W/m2 TOA, +7.4 W/m2 atmospheric) | scaled from the table in `exoplasim/notes/corrk-cross-check.md`; arithmetic below |
| PHYS-10 | line-list CO2 coefficients | about -0.03 K | the same table's CO2 row prices -7.5%; the landed fit is -6.7% at the planet path |
| CLIM-16 | `nhdiff = 1`, `hdiffk = 1000`, meaning 1000 m2/s since `world-mll` | 0.00 W/m2 global by construction; +0.03 to +0.12 K via the ice edge | this note |
| CLIM-17 | `TFREEZE` from declared salinity | 0.000 K | this note |
| SPEC-1 | model reads `k25v` | 0.00 to +0.03 K | `analysis/error_budget.json`, measured on the warm state; the A/B re-measures it |
| SPEC-5 | `vegetation_albedo` 0.165, bands [0.075, 0.225] | -0.64 K, bracket 0 to -0.90 | `analysis/vegetation_albedo.json`; +0.0105 on composited land mean at 0.61 K per 0.01. One-signed toward a cooler simulated mean |
| PHYS-11 | cloud absorption arms, bracket [0.78, 1.28] | 0 +/- 0.6 K, sign unassigned | this note; arms only, nothing enabled in the baseline config |
| DUST-11 | prescribed dust arm | separate note | `aeolian/notes/prescribed-dust-run.md` |
| CLIM-32 | ozone arms, `O3SCALE` 0.794 -> 0.500 | +0.2 to +0.9 K, warming | this note; arms only, `ozone_scale` unchanged in the baseline config |
| CLIM-33 | mixed layer arms, 25 m and 100 m | 0.000 K annual mean by construction | this note; arms only, `mixed_layer_depth_m` unchanged in the baseline config |
| world-ofn | `clwhsc` derived from `gascon` and `ga` | sign not assigned; see below | `notes/audits/model-earth-centrism.md` finding 5 |

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
0.135 K. Scaling by 9.24: +7.4 W/m2 atmospheric absorption, +1.4 W/m2 TOA,
+1.25 K. The caveat is the scaling itself: the row's sensitivities come from
`exoplasim/notes/shortwave-water-vapour.md` and are being stretched ninefold,
so this is the softest prediction in the table and the first place to bisect
if the sum misses.

**Consequence to plan for, corrected 2026-08-20:** the central sum shifts the
simulated temperature at every candidate flux by about +0.6 K, which is 0.3% of
a flux ratio on the canonical 202 K per unit flux. Read that as a shift in the
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
0.2534, giving 0.795 K per W/m2 of TOA forcing.

| arm | d ASR W/m2 | as K | d sea ice | verdict |
| --- | ---: | ---: | ---: | --- |
| `o3_050` CLIM-32 | -1.312 | -1.04 | -0.0001 | resolved, **prediction falsified** |
| `cld_078` PHYS-11 | -3.304 | -2.63 | +0.0004 | resolved, **prediction falsified** |
| `cld_128` PHYS-11 | +3.178 | +2.53 | -0.0015 | resolved, **prediction falsified** |
| `mld_025` CLIM-33 | -0.428 | -0.34 | +0.0053 | ice path confirmed |
| `sal_30` CLIM-35 | -0.474 | -0.38 | +0.0000 | UNRESOLVED, no ice change to carry it |
| `sal_20` CLIM-35 | -0.105 | -0.08 | +0.0016 | area confirmed, flux unresolved |
| `mld_100` CLIM-33 | -0.177 | -0.14 | -0.0040 | unresolved |
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
W/m2 at the top of atmosphere, +/-2.6 K. That is what was measured. The
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
12.96 W/m2 per unit scale, that is **+2.49 W/m2 and about +2.0 K** -- larger
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
