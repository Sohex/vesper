# Earth constants inherited in silence, and one decision with no artifact

Audited 2026-08-19. The question is the one `absent-and-inherited-physics.md`
asked and this document takes to different ground: where does an Earth number
reach this world without any file saying that it was chosen? That audit found
three, in the ocean, the cryosphere and the soil. This one hunts in the
RADIATION's surface and cloud terms, in the orbit block, and in the derivation
that set the flux.

**Three findings are Earth constants arriving with nothing naming them, one is a
decision that no artifact records, and three are smaller.** The first is the
largest unpriced term in the albedo budget and the direction of the error is
known; the second is the only major shortwave term that was never re-weighted for
this star while four smaller ones were.

All seven are closed. Each keeps its evidence below and each records what the
work found, because two of the seven came back with a MEASUREMENT that
contradicted the audit's own expectation, and that is the part worth keeping:
finding 3's expected lapse rate was too steep, and finding 5's codification
failed to reproduce the anchor it was codifying.

---

## 1. The vegetated land albedo is Earth's, integrated against the Sun

`analysis/rock_albedo.py` exists because of exactly this mechanism, and states it
in its own docstring:

> An albedo is a reflectance integrated against the light falling on it.
> Published rock albedos are integrated against the Sun. This star is a K2.5V at
> 4965 K, so more of its output sits at longer wavelengths, and a mineral with
> absorption features in the near infrared is darker here than its published
> value.

Every rock class in the table went through that integration, and `playa_clastic`
went through a second reconciliation in `analysis/playa_albedo.py` that ends
"then re-weight from the Sun to this star, which the spectra allow directly".

**The vegetated endmember did not.** `build_surface_albedo.py` carried
`--vegetation-albedo` at 0.15, `--tree-albedo` at 0.13 and `--grass-albedo` at
0.19, none with a source, the first justified only by a helptext reading
"0.12-0.15 covers needleleaf through broadleaf". In `--mode vegetated` the 0.15
was painted over every land region outside `barren_rock_classes`, which is about
four fifths of this planet's land, and the result is the largest single item in
`analysis/error_budget.json`: bare rock 0.25705 against vegetated 0.17907, a
land-albedo difference of 0.078 and 4.766 K.

### Why the rock argument does not cover it, and why the sign flips

Vegetation is dark in the visible and bright in the near infrared, which is the
opposite slope to a mineral with near-infrared absorption features. A redder star
therefore makes a canopy BRIGHTER where it makes those minerals darker. `k25v`
puts 0.382 of shortwave below 0.75 um where the Sun puts 0.517, and that is a
large shift of weight into the half of the spectrum a leaf reflects.

This is also why `exoplasim/notes/stellar-spectrum-audit.md` does not settle it.
That document says "ground and ocean are nearly spectrally flat over this range
and barely move", which is true of rock and soil and is false of vegetation, the
most spectrally sloped natural surface there is. It was written when
`land_albedo_source` still described a rock table; the sentence survived the
switch to `vegetated` and now covers a surface it was never measured on.

### Measured, 2026-08-19

Method, deliberately the same as `rock_albedo.py` so the two are comparable: 553
green-vegetation VSWIR spectra from `references/ecospeclib-all/`, matching
`vegetation.*vswir*spectrum.txt` and so excluding the non-photosynthetic set,
each integrated against `exoplasim/inputs/stellarspectra/k25v_hr.dat` and against
a 5772 K Planck curve, over the overlap of the reflectance range with 0.35 to
2.5 um.

| class | n | against the Sun | against k25v | ratio |
| --- | ---: | ---: | ---: | ---: |
| tree | 350 | 0.2455 | 0.2712 | 1.106 |
| shrub | 199 | 0.2548 | 0.2784 | 1.093 |
| grass | 4 | 0.3028 | 0.3329 | 1.099 |
| all | 553 | 0.2492 | 0.2742 | **1.101** |

Repeated over 0.2 to 4.0 um with the reflectance held at 0.05 outside the
measured range, which is the right order for a leaf in both tails: 0.2305 against
the Sun, 0.2629 against k25v, ratio **1.141**. The tails move it UP rather than
down, because the Sun's ultraviolet excess is also dark for a leaf. So the ratio
is 1.10 to 1.14 and the truncated figure is the conservative one.

**Leaf reflectance is not canopy albedo**, and the absolute levels above should
not be read as one: a canopy traps light between elements and shows soil through
gaps, both of which lower the near-infrared more than the visible. What transfers
is the RATIO, by exactly the argument `playa_albedo.py` uses to cancel the
laboratory-over-field offset -- numerator and denominator share whatever the
measurement does not represent. A canopy treatment would give a smaller ratio
than a leaf one, so this is an upper bound on the correction and not a value.

### What it is worth, and what was applied

**Closed as `spec-5` on 2026-08-19**, the way the row prescribed: derived like
`playa_albedo.py` derives its class, visible to the config, cited.
`analysis/vegetation_albedo.py` is the script and `vegetation_albedo` is a
pipeline step. `config/planet.yaml` carries `vegetation_albedo: 0.165` with the
bracket `[0.150, 0.171]`, the uncorrected endmember to the tails-carried ratio,
which is the bracket the paragraph above argues for: the leaf ratio is an upper
bound on the canopy correction, so the bracket runs from no correction to the
full one.

**The two other endmembers were missed and were caught later.** `--tree-albedo`
and `--grass-albedo` stayed on their argparse defaults of 0.13 and 0.19 when
0.15 moved into the config, so `--mode modelled` -- the mode
`build_surface_albedo.py`'s own docstring calls the only one that is not an
assumption -- was blending un-reweighted values while the endmember mode read a
derived one. **world-9m5** re-derived both by the same method:
`tree_albedo: 0.143` bracketed `[0.130, 0.148]` and `grass_albedo: 0.209`
bracketed `[0.190, 0.217]`, using the population leaf ratio rather than the
per-class one because the two classes differ by 0.006 in ratio and 0.001 in
albedo and the grass class has four spectra in it.

At 1.10 to 1.14 the vegetated endmember moves from 0.15 to 0.165-0.171, which
over the four fifths of land that is not barren is +0.012 to +0.017 on the land
mean. The error budget's own conversion, 4.766 K for a land-albedo difference of
0.078, gives 0.61 K per 0.01 with attenuation already in it, so **0.7 to 1.0 K of
cooling**. That is larger than the priced playa lever, larger than the priced
dust radiative forcing, and it is one-signed.

Two consequences beyond the kelvin. `baseline_flux_earth` was derived to land a
design mean at this endmember, so it moves with it. And the bare-to-vegetated gap
narrows by 15 to 20%, which is the quantity behind "the two endmembers reach a
given design mean at non-overlapping fluxes" -- the non-overlap is what makes the
orbit and the biosphere one choice, and it is measured on the pair.

### The related half: 174, 175 and 176 all carried the same field

`build_surface_albedo.py` wrote the identical array to all three codes, and
`exoplasim/notes/parameter-decisions.md` declared it honestly: "both bands get
the same value, because we have one albedo per rock class and no spectral split;
that reproduces single-band behaviour while still varying geographically". That
was a fair statement about a rock table. Under `vegetated` it asserted that a
canopy reflects equally either side of 0.75 um, which is the one thing a canopy
certainly does not do, and it was the same error as the paragraph above expressed
in bands rather than in a broadband mean. Correcting the broadband value was the
cheap fix; supplying a real 175/176 pair for the vegetated fraction was the
honest one, and the model reads them separately at `NSIMPLEALBEDO = 0`.

**The honest fix is what was done.** `config/planet.yaml` carries
`vegetation_albedo_bands`, the same integrals split at 0.75 um and anchored so
their flux-weighted combination is the broadband value, and codes 175 and 176
carry the pair over vegetated ground.

The same class then turned up inside the model, in the SNOW albedo rather than
the ground albedo, where nothing had corrected it:
`notes/audits/model-earth-centrism.md` finding 7 found `landini` masking snow
under a canopy with one spectrum-blind pair of Earth-Sun broadband ratios,
applied identically either side of 0.75 um. **world-nfh** replaced it with a
per-band mixture against `albforest(2)`, which `run_exoplasim.py` writes from
`vegetation_albedo_bands`. So the two halves of the class -- what this project
supplies and what the model computes -- now use one derivation.

## 2. The cloud shortwave constants are the one solar-weighted term never corrected

This project has now re-weighted four absorptances for a non-solar host, on the
argument that a Lacis-Hansen absorptance is a fraction of SOLAR flux and every
one of them has to move: `ozone_uv_weight`, `ozone_visible_weight`,
`h2o_sw_weight` at 1.346 and `co2_sw_weight` at 1.510. `radmod.f90` carries,
beside them:

    real :: tswr1   = 0.077  ! tuning of cloud albedo range1
    real :: tswr2   = 0.065  ! tuning of cloud back scattering c. range2
    real :: tswr3   = 0.0055 ! tuning of cloud s. scattering alb. range2
    real :: acllwr  = 0.100  ! mass absorption coefficient for clouds (lwr)
    real :: rcl1(3) = (/0.15,0.30,0.60/) ! cloud albedos spectral range 1
    real :: rcl2(3) = (/0.15,0.30,0.60/) ! cloud albedos spectral range 2
    real :: acl2(3) = (/0.05,0.10,0.20/) ! cloud absorptivities spectral range 2

all of them in `radmod_nl`, all of them Earth tunings, and none of them mentioned
anywhere in this repository. A grep for `tswr`, `acl2`, `nswrcl`, `rcl1`,
`acllwr` or `clgray` outside `vendor/` returned nothing: no note, no config
entry, no audit, no task.

**The band partition is star-aware and the within-band constants are not.**
`zsolar1` and `zsolar2` are recomputed from the spectrum, which is what the
0.382/0.618 split IS, so the weight given to each range is right. What is Earth's
is the physics inside each range. Range 2 is where this bites: the star puts
0.618 of its flux there against the Sun's 0.483, and within range 2 a K dwarf's
energy sits nearer 0.8-1.5 um, while the Sun's range-2 energy is spread further
into the 2-3 um region where liquid water absorbs much more strongly. So the
Earth tuning is being applied to a differently shaped band, and the sign is not
obvious by inspection -- which is the argument for measuring it, not for leaving
it. For scale, `exoplasim/notes/shortwave-water-vapour.md` books cloud droplets
and the multiple-scattering enhancement at 12.0 W/m2 of Earth's shortwave
absorption, beside the water vapour term that earned a 1.346 correction.

The measurement was cheap and satisfied every one of
`docs/src/pipeline/sequencing.md` A3's four conditions with no work: these are
namelist keys on the binary already built, so two arms differing by one key can
branch off one restart and be labelled diagnostic.

**Closed as `phys-11` on 2026-08-20, at 1.192, and it was never a decision.**
The argument that settled it is the one `config/planet.yaml` already made for
the four smaller absorptances: a Lacis-Hansen absorptance is a fraction of SOLAR
flux and every one of them has to be re-weighted for this host, so the arms
measure how much it MATTERS and not whether it is right, which is
physics-is-not-a-knob read the way CLIM-16 reads `nhdiff`. `cloud_band_weight.py`
is the pipeline step, and it derives the factor the way the gas weights were
derived: from a quantity carrying no spectrum, liquid water's k(lambda) out of
Hale and Querry (1973) Table I, through Mie to a co-albedo, flux-weighted over
range 2 against a 5772 K Sun. 1.184 to 1.205 across droplet effective radii 5 to
15 um and across both the weak-absorption and thick-layer scalings, so the answer
barely depends on either; 1.192 is the median.

`cloud_absorption_scale` multiplies ONE key, `cloudabs`, the scale `swr` applies
to the range-2 cloud co-albedo it interpolates from Stephens et al. (1984) Table
1(a). It used to arrive as `tswr3`, the coefficient of the analytic fit that
table replaced; the derived value did not change with the carrier, because it is
a star-over-Sun ratio of range-2 flux-weighted co-albedo and the table carries
the solar weighting in the denominator. world-f9ig. Worth about +2.0 K on the
arms' own measured slope, which is why the flux re-derivation has to follow it
and not precede it. The arms measured +/-2.6 K.

The harness applied the factor to the `acl2` triplet as well, and `acl2` has no
reader outside the `nswrcl == 0` branch of `swr`. `nswrcl` is 1 in every run
this project has made, including the PHYS-11 arms, so the whole of the measured
slope is the computed branch's co-albedo's and the `acl2` half of the write was
inert. The scaled write
is gone and the harness now stages `NSWRCL` at 1 explicitly, so a run's namelist
records which of the two cloud schemes its shortwave used. clim-68.

**What is NOT covered, and stays Earth's:** `tswr1`, `tswr2`, `acllwr`, `rcl1`
and `rcl2`. Those are albedo-like and longwave rather than shortwave
absorptances, so the Lacis-Hansen argument does not reach them, and nothing has
measured a replacement. `notes/audits/model-earth-centrism.md` adds a second
mechanism that pins `tswr1`, `tswr2`, `tswr3` and `th2oc` at T21's row whatever
the rung: `radmod`'s per-(NTRU, NLEV) `jtune` table is guarded by
`if(NDCYCLE==1) jtune=0` and `ndcycle` defaults to 1, so every branch of it is
dead. That is `world-ys9` and is OPEN.

## 3. Earth's lapse rate, in three places, on a world whose adiabat is 31% steeper

    exoplasim/scripts/dust_forcing.py    LAPSE_K_PER_KM = 6.5
    maps/build_basemap.py                t_surface = warmest_hi - 6.5 * (elev - clim_elev_hi)
    notes/glacier-rough-pass.md          "Lapse 6.5 K/km", bracketed 5.5 to 6.5

6.5 K/km is Earth's mid-latitude environmental lapse rate. None of the three said
so, none derived it, none read it from a config, and there was no `lib/` home for
it -- which is how there came to be three copies.

The dry adiabat is `g/cp`. The composition here is Earth-like so `cp` is within a
per cent of Earth's, and the ratio is therefore essentially the gravity ratio:
**12.75 K/km against Earth's 9.76**, steeper by 1.306. The moist adiabat scales
with it. So the expected environmental rate on this world is nearer 8.5 K/km than
6.5, and the bracket in `glacier-rough-pass.md` explored 5.5 to 6.5 -- entirely
on the shallow side of any defensible value.

**Closed as `phys-12` on 2026-08-19, and THE EXPECTATION IN THE PARAGRAPH ABOVE
FAILED THE MEASUREMENT.** That is the interesting half and is why the fix was to
measure rather than to declare a corrected constant. `lib/lapse.py` is the one
home: it takes the mean temperature on the 10 sigma levels from the climatology
the config names, weighted by records per output bin, builds level heights
hypsometrically with `R` and `cp` DERIVED from the configured composition rather
than declared, least-squares fits T against z per land column over sigma 0.45 to
0.90 -- about 1.1 to 3.5 km, the band every consumer extrapolates across -- and
land-area weights the result. Two seasons are exposed because the consumers want
different ones: `annual` for anything extrapolating an annual mean, `warmest`
for anything extrapolating the warmest bin. Two checks that can fail are raised
on rather than warned about: the rate must be positive and below the dry adiabat,
and the land-mean fit rms must stay under 1 K.

Measured on the baseline climatology of `run_8c2e1ff9ab5e`, 2026-08-19: annual
6.78 K/km, warmest bin per cell 7.82 K/km, dry adiabat 12.75 K/km. So the
measured rate is BELOW the 8.5 this audit predicted by scaling Earth's 6.5 with
the gravity ratio: this atmosphere is more stably stratified relative to its own
adiabats than a pure `g/cp` scaling assumes. Deliberately NOT a check in
`lib/lapse.py`: a moist-adiabat floor, because the measured annual rate sits
below the moist rate at the window-mean state and flooring there would decide the
answer instead of checking it.

That note was honest that the rate "is assumed... not taken from the model's own
temperature profile, which is available on 10 levels", and it listed the
assumption among the reasons its areas are an upper bound; it now carries the
measurement as a dated correction. The direction is what was missed: 5.5 halves
the glacier area and the measured 7.82 raises it, so this term ran opposite to
every other caveat on the same result. It is load-bearing because `docs/src/pipeline/state.md` section 5b's central geographic claim --
glaciers come from relief rather than latitude, a summit needing about 1.37 km
above its grid cell -- is `z* = cell_mean + (T_warmest - 273.15) / lapse` and
scales inversely with the rate.

`dust_forcing.py` uses it to set the temperature of the dust layer for the
longwave forcing, and `build_basemap.py` to decide where permanent ice is drawn.
The first reaches the carve verdict through `carve_verdict.py --dust-forcing`.
Both take their rate from `lib/lapse.py` now, `dust_forcing.py` the annual one
and `build_basemap.py` the warm-season one, and neither declares a number.

## 4. The longitude of perihelion is ExoPlaSim's Earth default, inside a DETERMINED block

`config/planet.yaml`, in the `orbit` block:

    longitude_vernal_equinox_degrees: 102.7

`vendor/exoplasim/exoplasim/__init__.py`, documenting the argument it feeds:

> Longitude of periapse, measured from vernal equinox, in degrees. If not set,
> defaults to Earth's (102.7 deg).

`docs/src/reference/config-rationale.md` had a section for every other key in
that block and none for this one. Nothing anywhere stated it as a choice. This
was the `TFREEZE` pattern one level up: with salinity the Earth value was
compiled into the model and invisible, and here it had been copied into the
config, which is worse in one specific way -- it read as a decision.

**What it decides is not small on this world.** It sets the phase of perihelion
against the solstices, and at `eccentricity: 0.02` the perihelion-to-aphelion
flux ratio is `((1+e)/(1-e))^2` = 1.083, so about 8% peak to peak. That is
LARGER than
`stellar_cycle.total_amplitude_flux_peak_to_peak`, which is 0.060, and the
stellar cycle is the forcing this world is built around. Perihelion phase decides
which hemisphere gets the short hot summer and which the long mild one, and
`docs/src/pipeline/state.md` section 5b establishes that summer is the term the glaciers turn on, since
summer amplification is nearly flat with latitude while winter amplification runs
three times as far. At 32 degrees obliquity the seasonal cycle it modulates is
stronger than Earth's, not weaker.

It is also already propagated: `biosphere/notes/lpj-guess-porting-audit.md`
derives a 10.5-day solstice phase offset from it, and that offset is compiled
into LPJ-GUESS through `vesper.h`.

Whether 102.7 stays was a decision and not a defect -- nothing here determines a
perihelion longitude, so it is DECLARED at best, in the sense
`config/planet.yaml` already uses for `ocean.salinity_psu`. What was a defect is
that it was indistinguishable from a value someone chose.

**Closed as `clim-23` on 2026-08-19: DECLARED, not moved.**
`docs/src/reference/config-rationale.md` now has the section, and it states the
phase is a snapshot choice of exactly the kind pCO2 is, since Orogen has no time
axis and no precession is modelled. It also states what 102.7 buys: perihelion
falls near the southern summer solstice, so the south is the amplified
hemisphere, and the amplification is hemispherically ANTISYMMETRIC where the
stellar cycle is not. Kept rather than moved because the phase is already
compiled into LPJ-GUESS through `vesper.h`'s solstice offset and carried by every
run and climatology, so an alternative costs a biosphere rebuild, a driver
regeneration and a re-run for a quantity nothing constrains. The symmetric
alternative is named, perihelion at an equinox, and choosing it is recorded as a
worldbuilding decision about the simulated seasons rather than a correction.

## 5. The derivation that set the flux has no artifact, no thresholds and no step

`baseline_flux_earth: 0.945` is marked DETERMINED "from the
habitability-by-latitude derivation rather than from a temperature target". That
derivation exists only as prose in `docs/src/pipeline/state.md` section 5b: summer and winter temperature per
band "measured on three converged runs and projected across candidate means",
with the conclusion that a lower flux moves half the land's warm-season monthly
means back inside the design comfort band at the cost of a tenth of it going
from harsh to extreme, and that 0.945 "puts the tropics near +33 C in their
warmest month rather than +36". (5b's original wording for that first clause
borrowed human-physiology vocabulary; it is restated here and rewritten there,
same content.)

There is no script, no analysis product, no row in `config/pipeline.yaml`, and no
statement of which warm-season ceiling the bands were scored against. The
three runs are not named. `compare_albedo_bracket.py` still carries the
superseded 290-293 K target and points the reader at the section.

Three of this project's own rules applied and all three failed: an artifact that
no step generates does not exist; thresholds are fixed before results are seen
and written down; record what a thing IS and where it lives.
`docs/src/pipeline/state.md` section 5b was written to replace a target whose
only justification was its own persistence, and it records that lesson in its own
text -- but the replacement was not reproducible either, and reproducibility is
what the complaint was about. It mattered operationally because section 6
requires the flux to be RE-DERIVED on every new terrain, so the next person to do
it had to reinvent the criteria rather than re-run them.

**Closed as `clim-24` on 2026-08-19, and the codification FAILED to reproduce
the anchor under its own declared rule.** That failure is the finding confirmed
by measurement rather than by argument, and it is why the fix is a script and not
a paragraph. `exoplasim/scripts/derive_design_flux.py` is a pipeline step with a
declared `writes` and a gate, and it states its thresholds at the top before it
runs. Under an unconstrained rule the winner is 0.8725, not 0.945; 0.945 comes
back only under a cap on the extreme-cold land fraction, and the value of that
cap, 0.0642, was INFERRED from the anchor rather than declared. The script says
so at the inference and the pipeline gate says so too: the next terrain's
re-derivation must declare its cap in advance, because a threshold fitted to its
own answer is what CLIM-30 exists to undo.

**Separately, and this is the bias rather than the bookkeeping, and it is not
settled by the codification:** the quantity the band criterion stands in for
couples temperature to humidity, and the criterion scores on monthly-mean
temperature alone. On a world this
arid, and one whose aridity is distributed by endorheic drainage rather than by
latitude, the two diverge and they diverge differently in each of the bands
being traded against each other. The tropics here are dry, so a temperature-only
ceiling penalises them too hard; a humid coastal margin at a lower temperature
scores better than it should. Codifying the criterion made WHICH quantity is
used checkable, and the answer is that it is still temperature alone, so the bias
stands and is now visible in a script rather than buried in prose.

## 6. A water-budget closure labelled per orbit is per Earth year

`exoplasim/scripts/build_climatology.py` computes the field named
`land_p_minus_e_minus_mrro_mm_per_orbit` with `seconds = 86400.0 *
EARTH_CALENDAR_YEAR_DAYS`. `lib/orbit.py` documents that constant as "the
Gregorian mean year. Use for annualising rates to 'per Earth year'", beside
`EARTH_SIDEREAL_YEAR_DAYS` and a paragraph explaining that the two names exist so
a reader does not have to work out which was meant.

This world's orbit is about half an Earth year, so the number was roughly twice
what its name said. Nothing outside `bootstrap_climate_series.json` and
`baseline_climate_series.json` read it, so this was a mislabel rather than a
physics error -- but it is a water-budget closure, the closure note beside it
invites the reader to interpret its magnitude, and a quantity with two meanings
is `docs/src/practice/failure-modes.md`'s own recurring shape.

**Closed as `clim-25` on 2026-08-19, by making the value true to its name rather
than by renaming it.** The closure scales by `orbital_year_seconds` read from the
run manifest, and `build_climatology.py` RAISES when the manifest does not carry
it rather than falling back on a guess. The comment at the call site names
`lib/orbit.py`'s two constants and why the wrong one was reachable.

## 7. The photosynthetic window rests on four references, all of them held

`references/INDEX.md` opens its first section by saying what turns on it: the
bare-rock-versus-vegetated item is the largest in the error budget, "so what the
biosphere does under a K dwarf sets where the planet is placed", and the
near-parity conclusion "is built entirely on the first two of these". All four
papers in that table are marked `held`, including Lehmer et al. (2021), which the
same table says the near-parity argument "turns on", and which it also records as
itself a secondhand use of Marosvolgyi and van Gorkom (2010).

`FRADPAR` is computed from the spectrum rather than taken from a paper, so the
arithmetic is sound; what was secondhand is the 400-750 nm WINDOW that the
arithmetic integrates over, and the window is the whole content of the result.
That was `docs/src/practice/failure-modes.md` class 9 on a load-bearing quantity.

**Closed as `ref-9` on 2026-08-19: all four are read.** `references/INDEX.md`
carries what each supports and what it cannot. Two things the reading changed
rather than confirmed. Lehmer et al.'s K2V is 4620 K where this host is 4965 K,
so interpolating their peaks in temperature puts the reddest optimum near 720 to
725 nm here: the 750 nm edge carries about the margin over the reddest peak that
Earth's 700 nm edge carries over chlorophyll a at 672, rather than the razor-thin
margin a true 4620 K host would give. The window stands with that margin stated.
And Kiang et al. (2007a) predates chlorophyll f, so the far-red photoacclimation
specifics in `biosphere/notes/productivity-prediction.md` are not from there and
currently cite nothing.

---

## Checked and clean

Recorded so they are not re-derived.

**The two-band snow, ice and glacier ENDMEMBERS are the model's own and are
correct.** Those albedos are derived inside `radmod` from the spectrum, which is
what the `k25v` work fixed. Finding 1 was about the surface fields this project
supplies. What that work did NOT cover, and what
`notes/audits/model-earth-centrism.md` findings 7 and 12 later found, is what
`landmod` and `seamod` apply ON TOP of those correct endmembers: a
spectrum-blind canopy mask over snow, and an Earth broadband SLOPE on the
sea-ice ramp. Both are fixed, by `world-nfh` and `world-cj4`. The lesson for
this audit's own class: deriving an endmember from the spectrum does not make
everything downstream of it star-aware.

**The shortwave band partition responds to the star.** `zsolar1` and `zsolar2`
are computed from the spectrum rather than left at the solar 0.517/0.483, so
finding 2 is about the constants inside each range and not about their weights.

**Penman's psychrometric constant uses the actual pressure**, `gamma = CP_AIR *
p_air / (0.622 * lam)` in `carve_verdict.py`, so the one place where a
gravity-dependent pressure could have entered a hardcoded Earth constant does
not.

**The economic minerals config is the model for citation state**, with an
explicit `sourced` / `sourced-negative` / `declared` taxonomy per rule and the
quoted sentence beside each threshold. Nothing in finding 1 or 3 has a
counterpart there.

**`pedogenesis.yaml`'s water and catena constants are no longer uncited**, having
been closed under `LITH-24` as DECLARED with brackets and a reported spread.

## What this audit did not cover

Stated so the coverage is not overread. I did not audit `maps/` beyond the lapse
rate it shares with finding 3, the LPJ-GUESS PFT parameter values themselves (now
carried by `world-orok`, which found the exclusion bites hardest at the polar cap, where
water CAPTURE rather than cold or supply is what holds the cover at 0.030), the
Mie code behind `analysis/dust_optics.json`, the hydrography solver internals, or
Orogen's own generation constants.

The vendored model itself was also out of scope, and this audit says so in its
opening: it hunts in this project's own configuration and analysis layer.
`notes/audits/model-earth-centrism.md` took the same question into the fork and
found the same class there twenty-nine times over, so the two documents are one
sweep in two halves rather than duplicates.

## Tasks

Tracked in the `bd` issue tracker, not restated here.

All seven are closed. The `bd` row is the authority; the column below is what
each one settled on, because two of them settled somewhere the audit did not
predict.

| finding | id | what it settled on |
| --- | --- | --- |
| 1. the vegetated albedo was never re-weighted to this star | `spec-5` | derived and bracketed; the tree and grass endmembers were missed and caught later by `world-9m5`, and the two-band half by `world-nfh` |
| 2. the cloud shortwave constants are Earth tunings | `phys-11` | the two ABSORPTION-like keys re-weighted at 1.192; `tswr1`, `tswr2`, `acllwr`, `rcl1` and `rcl2` are untouched, and `world-ys9` is open on the `jtune` table that pins them at T21 |
| 3. Earth's lapse rate in three places | `phys-12` | measured, not declared, and the audit's own 8.5 K/km expectation failed |
| 4. the longitude of perihelion is an inherited default | `clim-23` | DECLARED at 102.7, not moved, with the symmetric alternative named |
| 5. the flux derivation has no artifact | `clim-24` | codified as a pipeline step, which then failed to reproduce 0.945 under an unconstrained rule; the cap is inferred and the next re-derivation must declare one |
| 6. a closure labelled per orbit is per Earth year | `clim-25` | made true to its name, and it raises rather than guessing |
| 7. the photosynthetic window rests on held papers | `ref-9` | all four read; the window stands with a stated margin, and one downstream claim was found to cite nothing |
