# Physics review: what is calibrated for the wrong world, and what is absent

Audited 2026-08-17. The question is not whether numbers are current but whether
the PHYSICS is right for this planet: 12.81 m/s2, 1.2 Earth radii, a 30-hour day,
a 182.8-day year, 1 bar, around a K2.5V star at 0.945 S-Earth.

The archetype is DUST-6. Kok's saltation constants are fitted on Earth and his
standardization corrects for air density and nothing else, so Earth's gravity sat
inside a constant used at 1.31 g. Correcting it moved emission 30%. Nothing was
wrong; an Earth-calibrated parameterisation was used unchanged.

Findings are ranked by how far they would move a number this project quotes.
Each is labelled [numeric], [inspection] or [physics] by how far it was checked.

---

## 1. The shortwave water vapour absorption is weighted for the Sun [numeric]

**This is the same defect the project already found and patched in the ozone
terms, in the same scheme, left in the larger term.**

`radmod.f90:1548` says of its shortwave code: "from Lacis & Hansen (1974) for
clear sky (H2O, O3, Rayleigh)". Those absorptances are expressed as a fraction of
TOTAL INCIDENT SOLAR FLUX, so each carries the Sun's share of flux in the band it
represents. `notes/config-rationale.md` states this exactly, for ozone:

> Each of their three terms therefore carries the Sun's share of flux in the band
> it represents, and applying them unchanged to a K dwarf puts solar band weights
> on a non-solar spectrum.

The ozone terms were re-weighted accordingly: `ozone_uv_weight: 0.335` measured
from Segura et al. (2003), `ozone_visible_weight: 0.914` computed from the k25v
spectrum. **The water vapour term was not, and neither was CO2.** The
`radmod_nl` namelist carries `o3uvw` and `o3visw` and has no equivalent for
either; `th2oc` is the LONGWAVE continuum coefficient and is a different thing.

Water vapour absorbs in the near-infrared, which is precisely where a K dwarf
puts more of its flux. The model prints its own numbers at run time
(`MOST_DIAG.00070`):

| | share above 0.75 um |
| --- | ---: |
| this star, the model's own integration | 0.5816 |
| the Sun, the scheme's own reference (`zsolar1 = 0.517`, radmod.f90:311) | 0.483 |
| **ratio** | **1.204** |

So the H2O shortwave absorptance is being applied with about 20% too little
weight on the band that does the absorbing.

**Magnitude.** Measured atmospheric shortwave absorption on the baseline is
`rst - rss = 229.46 - 169.68 = 59.78 W/m2`. Water vapour is the dominant clear-sky
shortwave absorber; taking it at 55-75% of that total gives 33-45 W/m2, and a
1.204 scaling adds **+6.5 to +9 W/m2** of atmospheric absorption with an equal
reduction at the surface. The scaling is first-order and will be sublinear
because strong bands saturate, so treat it as an upper bound -- but even half of
it is comparable to the 5.7 W/m2 that made shortwave-only dust unacceptable, and
to a third of the 21 W/m2 that spans the entire 0.85-to-0.95 stellar sweep.

**Sign.** Two effects, both pointing the same way at the top of the atmosphere
and opposite ways at the surface. Intercepting the beam higher means less reaches
a bright surface to be reflected, so the planet absorbs MORE and warms -- and
this world's closed-basin fill is bright, 0.40 to 0.50. At the surface it is
7-9 W/m2 of shortwave removed, which suppresses evaporation directly. That is
the same channel DUST-10 exists to measure, at a comparable size.

**What settles it.** Integrating the Lacis & Hansen H2O band absorptances against
the k25v spectrum the way `analysis/playa_albedo.py` and the ozone weights were
done, and adding `h2o_swweight` beside `o3uvw`/`o3visw` in the same patch. The
infrastructure exists; `radmod.f90:311` already computes a spectrum-weighted
Rayleigh cross-section and knows the solar reference share.

---

## 2. Three different values are in use for this star's band-1 flux share [numeric]

The fraction of stellar flux below 0.75 um is load-bearing: it weights the
two-band surface, snow and sea-ice albedos, and it band-weights the dust optics.
Three values are in circulation:

| value | source |
| ---: | --- |
| 0.3777 | `analysis/dust_optics.json`, integrated from `k25v.dat` |
| 0.3862 | `k25v_hr.dat` integrated the same way |
| 0.4184 | what the model itself prints at run time |

The first two differ because `dust_optics.py:64` reads `k25v.dat`, which spans
0.340-14.01 um, while the high-resolution file spans 0.200-100 um. That is a
truncation difference and it is small.

**The model's 0.4184 is not explained by either file**, and it is 10.8% from the
value the dust optics use. It is the model's number that matters most, because it
weights the snow and ice albedo -- and getting a stellar spectrum wrong in
exactly that place is what previously made snow and ice 0.10-0.17 too dark and
cost this project a baseline re-run.

**Magnitude.** For dust, negligible: the band-1 and band-2 mass extinction
efficiencies are 779.4 and 766.9 m2/kg, so a 0.04 shift in weight moves the
flux-weighted value by 0.5 m2/kg on 771.6, or 0.07%. For the two-band albedo it
is not negligible -- snow albedo is 0.745 below 0.75 um and 0.431 above, a
difference of 0.314, so a 0.041 shift in the weight moves the broadband snow
albedo by 0.013.

**What settles it.** Read `radmod.f90`'s band integration and find why it differs
from a direct integration of either file. Then make one value canonical and have
everything read it, as `config/planet.yaml` already does for gravity.

---

## 3. The carve verdict transfers Earth's basin density without asking whether density is gravity-dependent [physics]

`hydrography/notes/retain-fraction.md` fixes the incision coefficient by matching
the DENSITY of Earth's standing through-flowing impounded basins, giving
C = 161 m per (m3/s)^0.5. That is a defensible calibration and it is honest that
it is one. But the target is an Earth observable and it is being transferred to a
1.31 g planet without asking whether the density itself is a function of gravity.

Two reasons to think it is:

**Incision is faster here.** Stream power per unit bed area is `omega = rho g Q
S / W`, so at the same discharge and slope this world cuts 31% faster than Earth.
Faster incision means fewer basins survive at a given age.

**Basins are shallower here.** The export manifest states that `gravityMS2 scales
maximum relief as 1/g`, and `reliefScale` is 0.766. Shallower basins are cut
through sooner. This half is already in the calculation, because the depth term
uses this build's own `depth_at_spill` -- it is the incision-rate half that is
not.

Both push the same way: Vesper should have FEWER standing basins per unit area
than Earth, not the same number. Matching Earth's density therefore
under-estimates C, which means over-preserving.

**Magnitude.** If C should scale as g^1 through stream power, 161 becomes 210,
which from the note's own sensitivity table moves the marginal count from 98
toward the 64 at C = 218. If it scales as g^1.5 it is 240. So of order 30-40
basins, in a class of 98 -- material for a landform class, small against the
1,512 carved.

**Why it is still a finding at that size.** The same argument applies to
`Q_full`'s successor whatever it is, and the calibration currently has NO gravity
term at all, so it will be wrong by the same factor at any future C.

---

## 4. Penman has no atmospheric stability correction [physics]

`carve_verdict.py:penman_open_water` builds its transfer coefficient as
`ce = KARMAN**2 / ln(z_ref/z0)**2` -- the neutral-stability bulk formula, with no
Monin-Obukhov correction.

The setting this is applied to is the one where it is least valid. A lake at
spill level in a hot arid closed basin is COOLER than the air above it, which is
a stably stratified surface layer; stability suppresses turbulent exchange, and
the neutral coefficient therefore overstates the aerodynamic term. Over Earth's
saline desert lakes this is a factor of 2 to 5 on the exchange coefficient in
strongly stable conditions, and rarely less than 20-30% on a daily mean.

**Sign is one-way and it matters for the verdict.** Overstated evaporation means
basins are judged more likely to close, so the verdict UNDER-carves. That is the
same direction as `analysis/error_budget.json`'s existing "Penman over a dry land
column" entry, so the two compound rather than cancelling.

**Not quantified here.** Doing so needs a bulk Richardson number per cell from
`ts`, `tas` and the wind, all of which are in the climatology, plus a declared
stability function. That is an afternoon and it would convert the error budget's
one "unquantified" item into a number.

---

## 5. Penman is evaluated at the mean state, not integrated over a 30-hour day [numeric]

Saturation vapour pressure is convex in temperature, about 6.7% per kelvin at
290 K, so the mean of `e_s` over a diurnal cycle exceeds `e_s` of the daily mean.
The carve verdict evaluates Penman on 12-bin climatological means and so misses
that rectification entirely.

Measured from the model's own daily extrema (`maxt`, `mint`) on the baseline
climatology:

| | diurnal range |
| --- | ---: |
| land, area-weighted mean | 8.08 K |
| land median | 5.4 K |
| land 90th percentile | 12.3 K |
| land maximum | 18.4 K |

For a sinusoidal cycle the rectification factor is `I0(a*A)` with `a = 0.0674` per
kelvin and `A` the half-amplitude: **+1.2% at the land mean and +4.4% at the
90th percentile**. Against the basin sensitivity table that is roughly 20 to 60
basins, in the under-carving direction, compounding with findings 3 and 4.

Smaller than I expected before measuring it, and worth recording as such: the
30-hour day is a real widening of the diurnal range but the model's land range is
only 8 K, not the 15-20 K a terrestrial desert reaches, because T42 and the soil
heat capacity damp it.

`biosphere/README.md` already records the related finding that LPJ-GUESS is
structurally blind to the diurnal range. This is the evaporation half of the same
observation and it is not recorded anywhere.

---

## Knocked down, with the evidence

**The lapse rate is fine, and one constant is right by luck.** The dry adiabat
here is `g/cp = 12.75 K/km` against Earth's 9.76, a 31% difference, so a
tropospheric structure tuned for Earth would be badly wrong. Measured from the
climatology's own `ta` on model levels: the lapse rate over the bottom 3.5 km is
**6.4 to 6.7 K/km**, which is what a g-scaled moist adiabat gives and is
Earth-like because latent heat release, not gravity, sets it. `dust_forcing.py`
declares `LAPSE_K_PER_KM = 6.5` and that is correct -- but it was chosen as
Earth's moist lapse rate rather than derived, and it should read the climatology.

**Rayleigh scattering is correctly handled.** `radmod.f90:311` computes
`rcoeff` as a spectrum-weighted cross-section normalised to a 5772 K reference,
and the run prints **0.862** for this star -- correctly reduced, because a red
star Rayleigh-scatters less at lambda^-4. The same expressions carry an explicit
`(9.80665/ga)` factor for column mass per unit pressure. Both the spectrum and
the gravity are right here.

**Snow and sea-ice albedo are spectrally resolved and re-weighted.** The run
prints snow albedo 0.745 below 0.75 um and 0.431 above, combined at the star's
own band shares. This was the subject of an earlier fix and it held.

**Hadley cell width barely moves, and gravity cancels.** The Held-Hou scaling
goes as `(g H dtheta/theta)^(1/2) / (Omega a)`, and `g H = R T` is independent of
gravity, so only the rotation term changes: `(Omega a)^2` is 0.917 of Earth's,
giving a Hadley cell **4% wider**. The Rhines scale relative to planetary radius
moves 2%. The large-scale circulation is close to Earth-like and none of this is
worth a task.

---

## Proposed task rows

| PHYS-1 | Re-weight the Lacis & Hansen SHORTWAVE water vapour and CO2 absorptances for the K dwarf, as `o3uvw`/`o3visw` already do for ozone. They give absorptance as a fraction of total incident SOLAR flux, and this star puts 1.204x the Sun's share above 0.75 um where water vapour absorbs. Worth +6.5 to +9 W/m2 of atmospheric shortwave absorption and an equal reduction at the surface, against 21 W/m2 for the whole stellar sweep | `notes/audits/physics-review.md` | open, LARGEST physics finding. Same scheme, same argument and same patch as the ozone weights, left in the bigger term. Integrate the H2O band absorptances against k25v; `radmod.f90:311` already knows the solar reference share |
| PHYS-2 | Make the star's band-1 flux share canonical. Three values are in use: 0.3777 in `analysis/dust_optics.json` from `k25v.dat`, 0.3862 from `k25v_hr.dat`, and 0.4184 printed by the model itself. The model's is unexplained by either file and it weights the two-band snow and ice albedo, which is where a stellar-spectrum error previously cost a baseline re-run | `notes/audits/physics-review.md` | open. Negligible for dust (0.07% on the mass extinction efficiency), 0.013 on broadband snow albedo. Find why the model's integration differs, then have one value and let everything read it |
| PHYS-3 | Give the carve verdict's incision coefficient a gravity term. C = 161 was fixed by matching Earth's standing-basin DENSITY, which is an Earth observable transferred to a 1.31 g world; stream power goes as rho g Q S, so this world cuts 31% faster and should have fewer standing basins than Earth, not the same number | `hydrography/notes/retain-fraction.md`, `notes/audits/physics-review.md` | open. Worth 30-40 basins of the 98 marginal at g^1, more at g^1.5. Direct analogue of DUST-6, and the calibration currently has no gravity term at any C |
| PHYS-4 | Add an atmospheric stability correction to the Penman transfer coefficient. It is the neutral bulk formula, applied to a cool lake under hot arid air, which is a stable surface layer where the neutral coefficient overstates exchange by 20-30% on a daily mean and more in strong stability | `hydrography/scripts/carve_verdict.py` | open, one-signed toward under-carving and compounding with the error budget's existing dry-column Penman item. Needs a bulk Richardson number from `ts`, `tas` and wind, all already in the climatology, plus a declared stability function. Would convert the budget's only unquantified item into a number |
| PHYS-5 | Integrate Penman over the diurnal cycle rather than evaluating it at 12-bin means. Saturation vapour pressure is convex at 6.7%/K, so the daily mean of the deficit exceeds the deficit of the daily mean by +1.2% at the land-mean diurnal range of 8.08 K and +4.4% at the 90th percentile of 12.3 K | `notes/audits/physics-review.md` | open, low priority. Roughly 20-60 basins, same direction as PHYS-3 and PHYS-4. `maxt` and `mint` are already in the climatology, so no new output is needed |
