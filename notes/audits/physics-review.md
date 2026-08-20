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
represents. `docs/src/reference/config-rationale.md` states this exactly, for ozone:

> Each of their three terms therefore carries the Sun's share of flux in the band
> it represents, and applying them unchanged to a K dwarf puts solar band weights
> on a non-solar spectrum.

The ozone terms were re-weighted accordingly: `ozone_uv_weight: 0.335` measured
from Segura et al. (2003), `ozone_visible_weight: 0.914` computed from the k25v
spectrum. **The water vapour term was not, and neither was CO2.** The
`radmod_nl` namelist carries `o3uvw` and `o3visw` and has no equivalent for
either; `th2oc` is the LONGWAVE continuum coefficient and is a different thing.

Water vapour absorbs in the near-infrared, which is precisely where a K dwarf
puts more of its flux:

| | share above 0.75 um | flux-share ratio to the Sun |
| --- | ---: | ---: |
| the Sun, the scheme's own reference (`zsolar1 = 0.517`, radmod.f90:311) | 0.483 | 1.000 |
| the 4965 K blackbody every run has actually used, finding 2 | 0.5715 | 1.183 |
| k25v, the spectrum `config/planet.yaml` names | 0.6154 | 1.274 |

So the H2O shortwave absorptance is being applied with too little weight on the
band that does the absorbing, and those ratios are the flux-share estimate of how
much. The 0.5816 that `MOST_DIAG.00070` prints is the blackbody's, measured on
`radmod`'s own truncated wavelength grid rather than bolometrically, which is why
it is 1.8% above the 0.5715 here.

**Settled, and it is bigger than this estimate.** Measured 2026-08-17 by
integrating Howard, Burch and Williams' laboratory band absorptions -- which are
the data Yamamoto weighted with the solar flux to make the fit Lacis and Hansen
then fitted -- against both spectra. Derivation and evidence in
`exoplasim/notes/shortwave-water-vapour.md`, re-runnable as
`exoplasim/scripts/shortwave_band_weights.py`.

The weight is **1.346** against k25v, bracket 1.301 to 1.363, and **1.202**
against the blackbody the model is running today. Both exceed their flux-share
ratios above, because a flux-share ratio weights every band by its flux and none
by its strength, and the strong bands are the long-wavelength ones where the K
dwarf's boost is largest: 1.50 at 2.7 um against 0.99 at 0.72 um. An
absorptance-weighted average has to exceed a flux-weighted one, so the flux-share
numbers are lower bounds.

**Magnitude, measured rather than bracketed.** Water vapour absorbs **43.2 of the
59.8 W/m2** the atmosphere takes, 72% of it, against the 55-75% assumed here.
Ozone takes 4.5 and cloud the remaining 12.0. The correction is therefore
**+8.8 W/m2** of atmospheric shortwave absorption if the blackbody is kept and
**+15.0** once the spectrum is declared, with **-7.6** and **-13.0 W/m2** at the
surface. It is not sublinear in the way the estimate expected: the patch
multiplies the absorptance itself, so the extra absorption is exactly
proportional to what water vapour already absorbs.

**Sign.** Two effects, both pointing the same way at the top of the atmosphere
and opposite ways at the surface. Intercepting the beam higher means less reaches
a bright surface to be reflected, so the planet absorbs MORE and warms: **+2.9
W/m2** at the top of the atmosphere with the spectrum declared, bracket +2.0 to
+4.5 depending on how much of the extra absorption sits above cloud rather than
below it. At the surface it is 13 W/m2 of shortwave removed, which suppresses
evaporation directly. That is the same channel DUST-10 exists to measure, at more
than twice the size. In surface temperature it is **+1.45 K** for this correction
alone against the blackbody, and a further +1.04 K for the spectrum acting
through this term.

**CO2 cannot be re-weighted, because there is none.** `radmod.f90` carries CO2
only in `lwr`, from Sasamori (1968); the shortwave has ozone in band 1 and water
vapour in band 2 and nothing else, which is faithful to Lacis and Hansen rather
than a defect in the port. What its absence is worth here is **2.5 W/m2**, the
same sign as the water vapour correction, against 1.75 W/m2 solar-weighted --
and that 1.75 is inside the 1.5-2.5 W/m2 Earth's near-infrared CO2 solar
absorption is measured at, which is the check on the calculation. It is a new
absorber rather than a re-weighting, so it is a separate task.

---

## 2. The model has been running on a blackbody since the first orbit of every run [numeric]

Resolved 2026-08-17. The finding started as three irreconcilable values for the
fraction of stellar flux below ExoPlaSim's 0.75 um band edge -- a number that
weights the two-band surface, snow, sea-ice and glacier albedos, divides the
ozone transmissivity, and band-weights the dust optics:

| value | where it was |
| ---: | --- |
| 0.3777 | `analysis/dust_optics.json`, integrating `k25v.dat` |
| 0.2566 | `dust_forcing.py`, integrating `k25v.dat` a second time and a different way |
| 0.3862 | `k25v_hr.dat` integrated naively at the band edge |
| 0.4184 | what the model prints at run time |

All four are now one, and the reason the model's differed is not an integration
difference at all.

### 0.4184 is a 4965 K Planck curve, not this star

`radmod.f90:813` takes the spectrum branch only when `NSTARFILE > 0`. Otherwise
`solarini` builds a blackbody at `STARBBTEMP` and reports its partition in the
same line of `MOST_DIAG`, with nothing distinguishing the two cases. Reproducing
`solarini`'s own grid -- 1024 logarithmic points from `minwavel` to 0.75 um and
1024 more from 0.75 to 100 um, trapezoidal, band edge assigned to band 1 --
gives 0.418350 for a 4965 K Planck curve and 0.384383 for `k25v_hr.dat`. The
model prints 0.418350399 and 0.384383172. Both to nine digits.

The check that makes this an explanation rather than a coincidence:
`radmod.f90:207` states that the scheme's default partitioning of 0.517 is what
a 5772 K spectrum produces through this code. The same reproduction returns
0.517000.

### Which orbits ran on which star

`MOST_DIAG.00000` of all four runs on this build reports 0.384383; every file
from `MOST_DIAG.00001` onward reports 0.418350. The first orbit of each run used
`k25v`; every orbit after it used a blackbody.

The mechanism is a namelist that gets rebuilt. `exo.Earthlike(workdir=...)`
copies the shipped namelists over the run directory on construction, so a driver
that reconstructs the model on an existing run and calls `configure()` without
`starspec` reverts `NSTARFILE` to 0. `continue_exoplasim.py` did, and it is the
script that runs all but the first orbit of every run;
`run_stellar_cycle.py` did too. `run_exoplasim.py` was the only driver that
staged the spectrum, and it runs one orbit. Both are fixed, and
`verify_stellar_spectrum` now refuses to run when the config names a spectrum
that is not in the namelist.

The four `radmod_namelist` files still on disk carry `STARBBTEMP = 4965.0` and
`NSTARTEMP = 1` and no `STARFILE` at all, which is the same evidence read from
the other end.

### What it is worth

The model's own log prices it, because `solarini` prints the albedos it derives
from the partition:

| | k25v | 4965 K blackbody |
| --- | ---: | ---: |
| energy fraction below 0.75 um | 0.384383 | 0.418350 |
| snow albedo below 0.75 um | 0.752962 | 0.745310 |
| broadband snow albedo | 0.538198 | 0.562267 |

**+0.0241 on the albedo of every snow-covered cell**, roughly double the 0.013
estimated before the cause was known, and in the same place as the k2-to-k25v
error that already cost a baseline re-run. A photospheric model has line
blanketing in the blue that a Planck curve does not, so the blackbody is bluer
than the star, which makes snow brighter than it should be. Sea ice, glacier and
ground albedos are derived from the same partition and move with it.

Sizing it against the surface it acts on is the standing caution, and here it
cuts the other way from the k2 case: sea ice is a couple of percent of this
world, but snow is seasonal over a much larger area, so this is not the same
"negligible because there is almost none of it" argument.

For dust it is negligible, as expected: band-1 and band-2 mass extinction
efficiencies are 779.4 and 766.9 m2/kg, so the chain value moved from 771.62 to
771.70 m2/kg, 0.011%. The dust FORCING moved more, but for a different reason --
`dust_forcing.py` was weighting the bands 0.2566/0.7434, having integrated
`k25v.dat` with the bin width counted twice. Correcting it moves the global-mean
net from +0.34 to +0.55 W/m2 and leaves every sign and every threshold in
`notes/dust.md` where it was.

### And the same subroutine gets the Rayleigh normalisation wrong

Found while pricing the above, because `solarini` prints `rcoeff` on the line
after the band fractions and the two logs disagreed by more than the star does.

`radmod.f90:311` computes `rcoeff`, which multiplies the Rayleigh optical depth
at `radmod.f90:1928`, as the star's lambda^-4-weighted cross-section normalised
to a 5772 K reference. The reference `bbg1`/`bbg2` is tabulated at lines 224-226
on `solarini`'s OWN logarithmic grid. Lines 230-233 then overwrite `wv1`/`wv2`
with the spectrum file's wavelengths, and the reference integrals at lines
287-299 run afterwards -- so they pair the reference's Planck values with the
FILE's wavelengths. The two grids start in different places, 316.036 nm against
0.2 um, and lambda^-4 weighting makes the short end decisive.

The blackbody branch is unaffected, because nothing overwrites the grid there.
The defect exists only when a spectrum file is used, which is the case nothing
has run in.

| | `rcoeff` |
| --- | ---: |
| 5772 K, the identity the normalisation is built on | 1.000000 |
| 4965 K blackbody, what every run after its first orbit used | 0.862015 |
| k25v, integrated on one grid | 0.712517 |
| k25v, as `solarini` codes it | 0.209732 |

`lib/stellar.py:rayleigh_coefficient` reproduces all four; the model prints
0.862014830 and 0.209731281, and the 5772 K row is exactly 1 by construction and
comes out exactly 1.

**This inverts the priority of the fix.** Staging the spectrum into the
continuation runs, which is the obvious repair and is now done, moves `rcoeff`
from 0.8620 to 0.2097: Rayleigh scattering 3.4x weaker than the runs have had,
where the honest value is 0.7125 and is 17% weaker. So the script fix must not
reach a production run before `solarini` is patched. The patch is small -- keep
the reference wavelengths in their own arrays instead of reusing `wv1`/`wv2` --
and it belongs in the same rebuild as any other `radmod.f90` change.

### The canonical value

`lib/stellar.py` reproduces `solarini` and is the only place the partition is
computed. `analysis/dust_optics.json`, `exoplasim/data/dust/vesper_dust_aerosol.provenance.json`,
`dust_forcing.py` and `world_state.json` all resolve through it, and the model
computes the same number from the same file. `k25v.dat` is not a source for it:
it spans 0.340-14.01 um and the 0.3777 it gives is that truncation, not the star.
The naive 0.3862 differs from the model's 0.384383 for two reasons that are both
in the Fortran -- flux below `minwavel` = 316.036116751 nm is zeroed rather than
merely unresolved, and the interval spanning the band edge is added to band 1.

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

**The gravity in the Rayleigh column mass is right.** `radmod.f90:1928` carries
an explicit `(9.80665/ga)` factor for column mass per unit pressure, so the
optical depth scales correctly at 12.81 m/s2. The spectral half of the same
expression does NOT hold, and is finding 2's second half rather than a knocked-down
item: the 0.862 the runs print is a 4965 K blackbody's coefficient, and the
spectrum branch that should replace it is broken by a grid mismatch.

**Hadley cell width barely moves, and gravity cancels.** The Held-Hou scaling
goes as `(g H dtheta/theta)^(1/2) / (Omega a)`, and `g H = R T` is independent of
gravity, so only the rotation term changes: `(Omega a)^2` is 0.917 of Earth's,
giving a Hadley cell **4% wider**. The Rhines scale relative to planetary radius
moves 2%. The large-scale circulation is close to Earth-like and none of this is
worth a task.

---

## Tasks

Tracked in `TASKS.md`, not restated here. A findings document that carries its
own rows is a second copy of the tracker, and the copy nobody works from is the
one that drifts.

| finding | id |
| --- | --- |
| 1. shortwave water vapour weighted for the Sun | `PHYS-1` |
| 1b. no shortwave CO2 absorptance exists at all | `PHYS-6` |
| 2. the runs are on a blackbody, and solarini's Rayleigh reference is mis-gridded | `PHYS-2`, `SPEC-1` |
| 3. the incision coefficient has no gravity term | `PHYS-3` |
| 4. Penman has no stability correction | `PHYS-4` |
| 5. Penman evaluated at the mean state | `PHYS-5` |
