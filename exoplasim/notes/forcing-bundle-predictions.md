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

**Predicted magnitude bound: 0 +/- 0.6 K.** Clouds are booked at 12.0 W/m2 of
Earth's shortwave absorption beside the water vapour term; +/-28% of that is
+/-3.4 W/m2 atmospheric, and at the same atmospheric-to-TOA ratio the water
vapour calibration row exhibits (0.19), +/-0.65 W/m2 TOA, +/-0.54 K. After
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

**Predicted magnitude: both arms beyond 2 K, plausibly 5 to 10 K**, from a
modelled low-cloud band-1 reflectance moving about 0.09 absolute against
PHYS-11's measured 3.2 W/m2 per 2.6 K. The bound is wide on purpose: the
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
greenhouse trapping of order 150 W/m2 is roughly 12 W/m2 and roughly 10 K by
`lib/sensitivity.py`. **That conversion is static and leans high**: it carries no
lapse-rate response and no water vapour feedback, both of which a run has and
both of which damp it. The prediction carried into the arms is therefore an
ORDER, 5 to 12 K end to end, and its purpose is to establish that the term is
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
  BINARIES, which A3 forbids because it confounds the term with the rebuild.
  **No pre-batch-2 binary was preserved** -- the only saved set is
  `prerebuild_binaries_2026-08-18`, eight days and many source changes earlier
  -- so the control arm cannot be reconstructed from a saved executable either.
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

**Every arm relaxes to a state at least 12.5 K colder than the donor, and is
still falling.** The control arm, which carries `config/planet.yaml`'s own values
for everything except the two energy keys:

| orbit | 0 | 4 | 9 | 14 | 19 | 24 | 31 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| global mean `ts`, K | 290.78 | 284.81 | 281.76 | 280.98 | 280.37 | 280.18 | 279.80 |
| sea ice, fraction | 0.0099 | 0.0340 | 0.0568 | 0.0667 | 0.0728 | 0.0737 | 0.0737 |

Sea ice reaches thirteen times the donor's 0.0058 and the ice-albedo feedback is
carrying most of the amplitude. The donor's last orbit is 292.26 K, so the drop
is at least 12.5 K and the arm has not stopped.

**The sum of the registered predictions is about +1 K of WARMING.** `world-jgen`
with `world-f9ig` is +0.59 to +1.63 K warming and is the largest single term;
`OCN-22` is -0.02 to -0.13 K; `world-o12h` is identically zero at this rung;
`world-2esd` and `world-trs3` are not swept in the baseline config. Nothing in
the table predicts a term of this size and nothing predicts cooling of this size.

**What the drift is NOT.** The arms differ from the donor's configuration in four
places and only one is large enough to matter for sign:

| difference | worth |
| --- | --- |
| `energy_fixer` and `energy_diagnostics` off | at most +0.3 K, and WARMING, so it makes the gap larger rather than smaller. The fixer applies -0.300 W/m2 on the donor and `lib/sensitivity.py` prices that at about 0.26 K |
| `ncpus` 16 to 8 | none. Thread count is compiled in and changes no physics |
| four `surface.cryosphere` constants, absent in the donor's config and declared now | unbounded here. They set the modelled sea ice's density, heat capacity and conductivity and the snow fusion enthalpy, and sea ice is where the amplitude is |
| `surface.land_water_column` from a 1-layer bucket to a 2-layer scheme | unbounded here |

So the measured drift confounds the compiled forcing terms with two config
changes, and the split above says which. **It is not attributable further**, for
the reasons the previous section gives: two of the five terms are compiled in
with no namelist that reverts them, three have no route from config, and no
pre-batch-2 binary was preserved to serve as a control.

**The finding is the disagreement itself and it does not need the decomposition.**
A bundle whose parts were predicted to sum to about +1 K moves the modelled
global mean by at least -12.5 K at T21. That is the outcome A3 attaches the
prediction for, and per `docs/src/practice/failure-modes.md` class 15 it is also
the shape two errors that nearly cancel would NOT produce: a residual this large
is one term being wrong by a lot, not two being wrong by a little.

**What this does not license.** Nothing here says a term should be removed.
Physics is not a knob, and a correct term that moves the simulated climate a long
way is information. What it says is that the next baseline cannot be taken on
this configuration until the two config changes above are priced, because they
are the two candidates that are cheap to price and are not yet.

**The cheapest next measurement, and it needs no new tooling.** Both config
changes are declared in `config/planet.yaml` and both are therefore arm-able the
way `clwref` was: one arm with the land water column back at the 1-layer bucket
and one with the four cryosphere constants at the values the donor ran without
them, each against this same control on this same restart. That is two arms and
it splits the -12.5 K into the part the config carries and the part the source
carries, which is as far as the split can go without a control binary.

## A3's fourth condition cannot be met on an arm's first segment

A3 requires that "the segments are labelled as diagnostics, so a short A/B tail
never enters a convergence window or a climatology". `continue_exoplasim.py`
takes a required `--purpose` and will stamp `diagnostic`. `run_exoplasim.py`
takes no such flag: every run it prepares gets a first segment labelled `spinup`
and a manifest carrying `canonical_lineage_eligible: true`, whatever the run is
for. Measured on the five arms above, each of which is an A/B arm and none of
which says so.

So an arm is self-labelling only from its SECOND segment onward, and its first
one -- which for a short A/B is the whole of it -- is indistinguishable on the
manifest from a spin-up meant for the canonical chain. Nothing has been mislabelled
INTO a climatology yet, because `baseline_climatology` is null and the canonical
lineage does not exist; the gap is that the guard A3 names is not there to catch
it when one does.

`--binary` already stamps an arm's build tag and sets `canonical_lineage_eligible
= false`, so the mechanism exists and reaches only arms that carry their own
executable. A `--purpose` on `run_exoplasim.py` with the same three values
`continue_exoplasim.py` takes, defaulting to `spinup` so no existing call
changes, is what closes it.
