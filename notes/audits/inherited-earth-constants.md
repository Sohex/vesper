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

**The vegetated endmember did not.** `build_surface_albedo.py` carries
`--vegetation-albedo` at 0.15, `--tree-albedo` at 0.13 and `--grass-albedo` at
0.19, none with a source, the first justified only by a helptext reading
"0.12-0.15 covers needleleaf through broadleaf". In `--mode vegetated` the 0.15
is painted over every land region outside `barren_rock_classes`, which is about
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

### What it is worth

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

### The related half: 174, 175 and 176 all carry the same field

`build_surface_albedo.py` writes the identical array to all three codes, and
`exoplasim/notes/parameter-decisions.md` declares it honestly: "both bands get
the same value, because we have one albedo per rock class and no spectral split;
that reproduces single-band behaviour while still varying geographically". That
was a fair statement about a rock table. Under `vegetated` it asserts that a
canopy reflects equally either side of 0.75 um, which is the one thing a canopy
certainly does not do, and it is the same error as the paragraph above expressed
in bands rather than in a broadband mean. Correcting the broadband value is the
cheap fix; supplying a real 175/176 pair for the vegetated fraction is the honest
one, and the model already reads them separately at `NSIMPLEALBEDO = 0`.

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
`acllwr` or `clgray` outside `vendor/` returns nothing: no note, no config entry,
no audit, no task.

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

The measurement is cheap and satisfies every one of WORKFLOW A3's four
conditions with no work: these are namelist keys on the binary already built, so
two arms differing by one key can branch off one restart and be labelled
diagnostic. Read it as physics-is-not-a-knob, the way `CLIM-16` reads `nhdiff`:
bracket and report the spread, do not tune within it.

## 3. Earth's lapse rate, in three places, on a world whose adiabat is 31% steeper

    exoplasim/scripts/dust_forcing.py    LAPSE_K_PER_KM = 6.5
    maps/build_basemap.py                t_surface = warmest_hi - 6.5 * (elev - clim_elev_hi)
    notes/glacier-rough-pass.md          "Lapse 6.5 K/km", bracketed 5.5 to 6.5

6.5 K/km is Earth's mid-latitude environmental lapse rate. None of the three says
so, none derives it, none reads it from a config, and there is no `lib/` home for
it -- which is how there came to be three copies.

The dry adiabat is `g/cp`. The composition here is Earth-like so `cp` is within a
per cent of Earth's, and the ratio is therefore essentially the gravity ratio:
**12.75 K/km against Earth's 9.76**, steeper by 1.306. The moist adiabat scales
with it. So the expected environmental rate on this world is nearer 8.5 K/km than
6.5, and the bracket in `glacier-rough-pass.md` explores 5.5 to 6.5 -- entirely
on the wrong side of the value.

That note is honest that the rate "is assumed... not taken from the model's own
temperature profile, which is available on 10 levels", and it lists the
assumption among the reasons its areas are an upper bound. The direction is what
was missed: 5.5 halves the glacier area, and the physically expected direction
roughly doubles it, so this term runs opposite to every other caveat on the same
result. It is load-bearing because WORKFLOW 5b's central geographic claim --
glaciers come from relief rather than latitude, a summit needing about 1.37 km
above its grid cell -- is `z* = cell_mean + (T_warmest - 273.15) / lapse` and
scales inversely with the rate.

`dust_forcing.py` uses it to set the temperature of the dust layer for the
longwave forcing, and `build_basemap.py` to decide where permanent ice is drawn.
The first reaches the carve verdict through `carve_verdict.py --dust-forcing`.

## 4. The longitude of perihelion is ExoPlaSim's Earth default, inside a DETERMINED block

`config/planet.yaml`, in the `orbit` block:

    longitude_vernal_equinox_degrees: 102.7

`vendor/exoplasim/exoplasim/__init__.py`, documenting the argument it feeds:

> Longitude of periapse, measured from vernal equinox, in degrees. If not set,
> defaults to Earth's (102.7 deg).

`notes/config-rationale.md` has a section for every other key in that block and
none for this one. Nothing anywhere states it as a choice. This is the
`TFREEZE` pattern one level up: with salinity the Earth value was compiled into
the model and invisible, and here it has been copied into the config, which is
worse in one specific way -- it now reads as a decision.

**What it decides is not small on this world.** It sets the phase of perihelion
against the solstices, and at `eccentricity: 0.02` the perihelion-to-aphelion
flux ratio is `((1+e)/(1-e))^2` = 1.083, so about 8% peak to peak. That is
LARGER than
`stellar_cycle.total_amplitude_flux_peak_to_peak`, which is 0.060, and the
stellar cycle is the forcing this world is built around. Perihelion phase decides
which hemisphere gets the short hot summer and which the long mild one, and
WORKFLOW 5b establishes that summer is the term the glaciers turn on, since
summer amplification is nearly flat with latitude while winter amplification runs
three times as far. At 32 degrees obliquity the seasonal cycle it modulates is
stronger than Earth's, not weaker.

It is also already propagated: `biosphere/notes/lpj-guess-porting-audit.md`
derives a 10.5-day solstice phase offset from it, and that offset is compiled
into LPJ-GUESS through `vesper.h`.

Whether 102.7 stays is a decision and not a defect -- nothing here determines a
perihelion longitude, so it is DECLARED at best, in the sense
`config/planet.yaml` already uses for `ocean.salinity_psu`. What is a defect is
that it is currently indistinguishable from a value someone chose.

## 5. The derivation that set the flux has no artifact, no thresholds and no step

`baseline_flux_earth: 0.945` is marked DETERMINED "from the
habitability-by-latitude derivation rather than from a temperature target". That
derivation exists only as prose in WORKFLOW 5b: summer and winter temperature per
band "measured on three converged runs and projected across candidate means",
with the conclusion that cooling "moves half the land out of sustained heat
stress at the cost of a tenth of it going from harsh to extreme" and that 0.945
"puts the tropics near +33 C in their warmest month rather than +36".

There is no script, no analysis product, no row in `config/pipeline.yaml`, and no
statement of what "sustained heat stress" is or what threshold was applied. The
three runs are not named. `compare_albedo_bracket.py` still carries the
superseded 290-293 K target and points the reader at the section.

Three of this project's own rules apply and all three fail: an artifact that no
step generates does not exist; thresholds are fixed before results are seen and
written down; record what a thing IS and where it lives. WORKFLOW 5b was written
to replace a target whose only justification was its own persistence, and it
records that lesson in its own text -- but the replacement is not reproducible
either, and reproducibility is what the complaint was about. It matters
operationally because section 6 requires the flux to be RE-DERIVED on every new
terrain, so the next person to do it has to reinvent the criteria rather than
re-run them.

**Separately, and this is the bias rather than the bookkeeping:** heat stress is
a wet-bulb quantity. The criterion as described is a dry-bulb monthly mean. On a
world this arid, and one whose aridity is distributed by endorheic drainage
rather than by latitude, the two diverge and they diverge differently in each of
the bands being traded against each other. The tropics here are dry, so a
dry-bulb threshold overstates their stress; a humid coastal margin at a lower
temperature can be worse. Nothing in the record says which quantity was used, so
this cannot be checked, which is the finding above restated.

## 6. A water-budget closure labelled per orbit is per Earth year

`exoplasim/scripts/build_climatology.py` computes the field named
`land_p_minus_e_minus_mrro_mm_per_orbit` with `seconds = 86400.0 *
EARTH_CALENDAR_YEAR_DAYS`. `lib/orbit.py` documents that constant as "the
Gregorian mean year. Use for annualising rates to 'per Earth year'", beside
`EARTH_SIDEREAL_YEAR_DAYS` and a paragraph explaining that the two names exist so
a reader does not have to work out which was meant.

This world's orbit is about half an Earth year, so the number is roughly twice
what its name says. Nothing outside `bootstrap_climate_series.json` and
`baseline_climate_series.json` reads it today, so this is a mislabel rather than
a physics error -- but it is a water-budget closure, the closure note beside it
invites the reader to interpret its magnitude, and a quantity with two meanings
is `notes/failure-modes.md`'s own recurring shape.

## 7. The photosynthetic window rests on four references, all of them held

`references/INDEX.md` opens its first section by saying what turns on it: the
bare-rock-versus-vegetated item is the largest in the error budget, "so what the
biosphere does under a K dwarf sets where the planet is placed", and the
near-parity conclusion "is built entirely on the first two of these". All four
papers in that table are marked `held`, including Lehmer et al. (2021), which the
same table says the near-parity argument "turns on", and which it also records as
itself a secondhand use of Marosvolgyi and van Gorkom (2010).

`FRADPAR` is computed from the spectrum rather than taken from a paper, so the
arithmetic is sound; what is secondhand is the 400-750 nm WINDOW that the
arithmetic integrates over, and the window is the whole content of the result.
That is `notes/failure-modes.md` class 9 on a load-bearing quantity, the file
flags it as an open exposure in its own header, and no task tracks it.

---

## Checked and clean

Recorded so they are not re-derived.

**The two-band split for snow, ice and glacier is the model's own and is
correct.** Those albedos are derived inside `radmod` from the spectrum, which is
what the `k25v` work fixed. Finding 1 is confined to the surface fields this
project supplies.

**The shortwave band partition responds to the star.** `zsolar1` and `zsolar2`
are computed from the spectrum rather than left at the solar 0.517/0.483, so
finding 2 is about the constants inside each range and not about their weights.

**Penman's psychrometric constant uses the actual pressure**, `gamma = CP_AIR *
p_air / (0.622 * lam)` in `carve_verdict.py`, so the one place where a
gravity-dependent pressure could have entered a hardcoded Earth constant does not.

**The economic minerals config is the model for citation state**, with an
explicit `sourced` / `sourced-negative` / `declared` taxonomy per rule and the
quoted sentence beside each threshold. Nothing in finding 1 or 3 has a
counterpart there.

**`pedogenesis.yaml`'s water and catena constants are no longer uncited**, having
been closed under `LITH-24` as DECLARED with brackets and a reported spread.

## What this audit did not cover

Stated so the coverage is not overread. I did not audit `maps/` beyond the lapse
rate it shares with finding 3, the LPJ-GUESS PFT parameter values themselves, the
Mie code behind `analysis/dust_optics.json`, the hydrography solver internals, or
Orogen's own generation constants. The first of those is the one I would take
next: the bioclimatic limits are Earth calibrations that
`biosphere/notes/lpj-guess-porting-audit.md` calls "the least mechanical part of
the port", and the port resolved the degree-day half while the temperature half
was carried over as physiology.

## Tasks

Tracked in `TASKS.md`, not restated here.

| finding | id |
| --- | --- |
| 1. the vegetated albedo was never re-weighted to this star | `SPEC-5` |
| 2. the cloud shortwave constants are Earth tunings | `PHYS-11` |
| 3. Earth's lapse rate in three places | `PHYS-12` |
| 4. the longitude of perihelion is an inherited default | `CLIM-23` |
| 5. the flux derivation has no artifact | `CLIM-24` |
| 6. a closure labelled per orbit is per Earth year | `CLIM-25` |
| 7. the photosynthetic window rests on held papers | `REF-9` |
