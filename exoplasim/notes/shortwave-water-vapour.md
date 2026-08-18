# The shortwave water vapour absorptance is weighted for the Sun

Measured 2026-08-17, on the `precarve-craton` baseline climatology.

This is the ozone band-weight correction repeated in the term that is ten times
larger. `notes/ozone.md` and `notes/config-rationale.md` record the ozone half.
This is the water vapour half, and it is settled by reading the papers rather
than the code.

## What the absorptance is a fraction of, from the papers

`radmod.f90`'s clear-sky shortwave is Lacis and Hansen (1974). Their water vapour
absorptance is their Eq. 21,

    A_wv(y) = 2.9 y / ((1 + 141.5 y)^0.635 + 5.925 y)

which `radmod.f90` codes verbatim at three sites. It is a fit to Yamamoto (1962),
and Yamamoto states the definition in one sentence:

> The definition of absorptivity is given by the ratio to the solar constant of
> the energy absorbed by the entire vertical air column for normal incidence.

So it is a fraction of TOTAL INCIDENT FLUX, and the Sun's spectrum is inside it,
because Yamamoto built it by weighting Howard, Burch and Williams' laboratory
band absorptivities with the solar flux and summing. Lacis and Hansen say the
same thing a second way in their Section 5a: "approximately 35% of the solar flux
is contained in the regions of significant water vapor absorption", which is
1 - p(k_1) = 1 - 0.6470 = 0.3530 from their own Table 1, and which they then use
as the literal constant 0.353 in their Eq. 39.

## Why dividing by `zsolar2` does not fix it

`radmod.f90` divides A_wv by `zsolar2`, the star's own flux share above 0.75 um,
and then applies the result to band-2 flux, which is `zsolar2` times the total.
The two cancel exactly. The absorbed flux is

    zsolar2 * S * (A_wv / zsolar2) = S * A_wv

whatever the star is. So the scheme gives a K dwarf exactly the Sun's absorbed
fraction of total flux. The `zsolar2` division is not a spectral correction at
all; it is a change of denominator that the multiplication undoes. Two
consequences worth stating because both are counter-intuitive:

- The defect is invisible in the code. Nothing in `swr` looks solar; the solar
  assumption is inside the number 2.9.
- The correction is INDEPENDENT of `PHYS-2`. Whatever the star's band-1 share
  turns out to be, it cancels out of the water vapour term. The two findings
  touch the same variable and do not interact.

## The derivation, and it runs backwards through the same three papers

`exoplasim/scripts/shortwave_band_weights.py` is the derivation; this is its
shape.

Yamamoto's construction is reversible because its ingredients separate. Howard,
Burch and Williams (1956) measured, for each near-infrared band, a TOTAL BAND
ABSORPTION in cm-1 -- an equivalent width, a property of the molecule, carrying
no incident spectrum whatever -- and fitted it as

    weak band     int A_nu d_nu = c w^(1/2) (P + p)^k
    strong band   int A_nu d_nu = C + D log10 w + K log10 (P + p)

with c, k, C, D, K and a transition value tabulated per band. Their own Eq. 11
defines the band-average fractional absorption as that divided by the band width,
and their own text states the method: "if the spectral distribution of the
radiation from a given source is known, the fraction of the total radiation
absorbed ... can be computed". So

    A_total(w) = sum_i f_i Abar_i(w)

with f_i the fraction of the incident flux in band i. Put the Sun in and this has
to reproduce Lacis and Hansen Eq. 21. Put this star in and it gives the
absorptance this star should have. The weight is the ratio.

The two spectra are BT-Settl models built through the same blend
`build_stellar_spectrum.py` uses -- 4900/5000 K interpolated to 4965 K for the
star, on the fids `build_stellar_spectrum.py` already pins, and 5700/5800 K
interpolated to 5772 K for the solar reference. Same grid, same converter, so a
model systematic divides out of the ratio rather than entering the answer.

## The checks that could have failed, and what they said

| check | expected | got |
| --- | --- | --- |
| solar flux fraction below 0.75 um | 0.517, `radmod.f90:125`, stated to hold at 5772 K | 0.5116, low by 1.1% |
| solar flux fraction above 0.9 um | 0.353, Lacis and Hansen Section 5a and Table 1 | 0.3690, high by 4.5% |
| reconstruction of A_solar against Lacis and Hansen Eq. 21, over the 0.01 to 10 cm interval they state the fit holds on | 1.00 | 1.109 to 1.319, median 1.134 |

The third is the one that matters and it partly failed, so it is worth being
exact about what survives.

**The reconstruction absorbs 11 to 16% more than Eq. 21 over the range that
matters, as a level offset rather than a shape error.** The likeliest cause is
that Yamamoto computed his curves for a real atmospheric column with a
Curtis-Godson effective pressure, while Eq. 21 is refitted as though it held at
standard pressure with the amount rescaled instead; a homogeneous 760 mm Hg path,
which is what Howard's constants describe, therefore absorbs more than Yamamoto's
column did. That is Lacis and Hansen's own approximation, not an error here.

**A level offset divides out of a ratio.** What would not divide out is the
excess being concentrated at one end of the spectrum, because the K dwarf's flux
boost is largest at the long-wavelength end. So the extreme attribution was
tested rather than argued away: charge the entire excess to the 2.7, 3.2 and
6.3 um complex and drop that complex altogether, which is Yamamoto's own curve 1.
The weight moves from 1.346 to 1.301. Dropping the 0.72 and 0.81 um bands
instead moves it to 1.363. **The weight is insensitive to the reconstruction's
level error**, and 1.301 to 1.363 is the bracket that produces.

**Both halves of that were checked afterwards against a modern line list and both
held.** `exoplasim/notes/corrk-cross-check.md` reruns this quantity from
correlated-k tables on HITRAN2020, with these same two spectra: the weight comes
out at 1.3271, inside the bracket, and the excess over Eq. 21 comes out at 1.127,
which confirms the level offset is a deficit in Eq. 21's pressure treatment rather
than an error here. The corollary is a separate defect that nothing fixes: the
model's absolute clear-sky water vapour shortwave absorption is low by about 12%,
and `h2osww` is a ratio and does not touch it.

The 0.72 and 0.81 um bands are the one assumed input here and it is worth saying
so plainly. Howard never measured them; Yamamoto estimated them from Fowle's data
and says they "cannot be neglected ... because of the large solar energy in this
region". They are carried at the 0.94 um band's shape scaled to 0.30 and 0.10,
which is a guess. It is a guess that can only make the weight SMALLER, because
those bands sit where the K dwarf's flux boost is least -- 1.05 and 0.99 -- so
omitting them would have flattered the correction, and the 1.363 endpoint is
exactly what omitting them costs.

## Which star, and it is not the one in `config/planet.yaml`

**No run on this build has been using the k25v spectrum.** `radmod.f90:813` takes
the spectrum branch only when `NSTARFILE > 0`, and every run's `radmod_namelist`
carries `STARBBTEMP = 4965.0` with `NSTARTEMP = 1` and no `STARFILE`, so
`solarini` has been building a 4965 K Planck curve. That is PHYS-2's finding and
it lands in the same rebuild as this one.

The weight is not the same number for the two, because a blackbody has no line
blanketing pushing flux out of the blue and into the near infrared, so it is less
red than the star it is standing in for:

| | flux below 0.75 um | flux above 0.75 um | H2O shortwave weight |
| --- | ---: | ---: | ---: |
| k25v, the real spectrum | 0.3846 | 0.6154 | **1.346** |
| 4965 K blackbody, what `solarini` builds today | 0.4285 | 0.5715 | **1.202** |

The 0.3846 is an independent cross-check on both sides of this: `lib/stellar.py`
makes the band-1 share canonical at 0.384383 by integrating the k25v file, and
this derivation reaches 0.3846 from the BT-Settl model the file was built from,
by a different path.

**The weight has to match whichever spectrum the run uses.** Setting 1.346 while
`solarini` is building a Planck curve would over-correct by 12%. The right answer
is to declare the spectrum and use 1.346, because a blackbody is not this star
and the choice between them is not a tuning decision.

## The weight

    h2o_sw_weight = 1.346          bracket 1.301 to 1.363

quoted at this world's own effective water path. The absorptance is evaluated at
a magnified path, so the weight is quoted there too: the baseline climatology's
area-mean column, built the way `radmod.f90:1802` builds it, is 1.680 cm, and
`zbetta = 1.66` makes the path 2.789 cm. The weight is a slowly falling function
of that path, 1.416 at 0.01 cm to 1.313 at 10 cm, so quoting it at Earth's column
or at a bone-dry one would be a different number; a constant is an approximation
worth about 7% across the range the model actually spans.

Per band, this is where it comes from:

| band | solar flux fraction | star flux fraction | ratio |
| --- | ---: | ---: | ---: |
| 6.3 um | 0.00410 | 0.00627 | 1.529 |
| 3.2 um | 0.00820 | 0.01260 | 1.537 |
| 2.7 um | 0.02410 | 0.03622 | 1.503 |
| 1.87 um | 0.04237 | 0.06425 | 1.516 |
| 1.38 um | 0.07580 | 0.10171 | 1.342 |
| 1.1 um | 0.05252 | 0.06356 | 1.210 |
| 0.94 um | 0.07569 | 0.08437 | 1.115 |
| 0.81 um | 0.04060 | 0.04260 | 1.049 |
| 0.72 um | 0.04036 | 0.04009 | 0.993 |

**Why this is 1.35 and not the 1.204 the audit estimated.** 1.204 is the ratio of
flux shares above 0.75 um, which weights every band by its flux and none by its
strength. The strong bands are the long-wavelength ones and they are exactly
where a K dwarf's boost is largest -- 1.50 at 2.7 um against 0.99 at 0.72 um --
so an absorptance-weighted average has to exceed a flux-weighted one. The audit's
number was the right shape and a lower bound, and it was described as
first-order, which it was.

## What water vapour actually absorbs here

Measured on the baseline by reimplementing `radmod.f90`'s own effective water
path against the climatology:

| | W/m2 |
| --- | ---: |
| incident at the top of the atmosphere | 321.6 |
| atmospheric shortwave absorption, `rst - rss` | 59.8 |
| of which water vapour | 43.2 |
| of which ozone | 4.5 |
| residual, cloud droplets and the multiple-scattering enhancement | 12.0 |

Water vapour is 72% of the atmosphere's shortwave absorption, inside the 55-75%
the audit assumed and now measured rather than assumed. The decomposition is
checkable in one place: the residual is cloud absorption, and 20% of atmospheric
shortwave absorption is what clouds take on Earth, so nothing is left over that
has no name.

Two things about that column integral are worth recording because both were
wrong on the first pass. `rsut` and `ssru` are signed upward-negative on this
stream, so the incident flux at the top of the atmosphere is `rst - rsut` and not
`rst + rsut`; taking the sign wrong gives a planetary albedo of -0.67, which is
loud, and a water vapour share of 30%, which is not. And the magnification
matters: at the model's own direct-beam factor `35/sqrt(1+1224 mu^2)` the water
vapour absorption comes out at 51.9 W/m2, which would leave only 3.3 for cloud.
`zbetta = 1.66` is the right representative because 56% of the sky is cloudy and
everything below cloud is on the diffuse path.

## The prediction, made before the run

Stated so the baseline re-run can falsify it. Two columns, because the run that
tests this will be the first one to see the real spectrum and the weight is a
different number under each:

| quantity | baseline | at 1.202, blackbody kept | at 1.346, spectrum declared |
| --- | ---: | ---: | ---: |
| atmospheric shortwave absorption, `rst - rss` | 59.8 W/m2 | **+8.8** | **+15.0** |
| surface net shortwave, `rss` | 169.7 W/m2 | **-7.6** | **-13.0** |
| top-of-atmosphere net shortwave, `rst` | 229.5 W/m2 | **+1.7** | **+2.9**, bracket +2.0 to +4.5 |
| mean surface temperature, `ts` | 289.71 K | **+1.45 K** | **+2.5 K**, bracket +1.6 to +3.9 |
| precipitation | 974 mm/yr | **-111 mm/yr** | **-189 mm/yr** at full compensation |
| flux ratio to restore the design mean | 0.945 | -0.0074 | -0.0127 |

**So +1.45 K is this correction and the further +1.04 K is the spectrum, acting
through this term.** The spectrum has other channels -- the two-band snow and ice
albedo, and Rayleigh scattering through SPEC-2, which is a factor of 3.4 -- and
none of those is in this number.

**The prediction is conditional on the baseline it starts from.** These are the
h2osww term alone, evaluated on a climatology produced with the blackbody, so if
the spectrum and SPEC-2 land in the same rebuild the run contains three changes
and none of them can be attributed from one result. That is worth two short
attribution segments at the new spectrum, `h2osww` off and on, before the long
converged run: cheap, and the only thing that turns this prediction into a test
rather than a story told afterwards.

**The atmospheric and surface terms are the robust ones.** They follow from the
absorptance being multiplied by a known factor and are close to arithmetic. The
+15.0 splits as +14.7 on the downward beam and +0.3 on the surface-reflected leg.

**The top-of-atmosphere term is a bracket and stays one.** The extra absorption
intercepts flux that would otherwise have gone on to be reflected, so the planet
absorbs more. How much more depends on what sits below the extra absorption: if
all of it were above every reflector, the planetary albedo 0.286 applies and the
forcing is +4.5 W/m2; if all of it were below the cloud, the surface albedo 0.119
applies and it is +2.0. Water vapour is bottom-heavy and cloud is not, so the
truth is inside, and +2.9 is a 65/35 split rather than a measurement.

**The temperature follows from the top-of-atmosphere term alone**, at 209 K per
unit flux ratio -- `WORKFLOW.md` section 5b, the converged T42 segment across
0.945-0.968, which is the direction this correction moves. That is 0.861 K per
W/m2 once the incident-to-absorbed conversion is made. The stellar sweep's 33 K
for 21 W/m2 is deliberately NOT used: it is 3.4x the local slope because it
crosses the ice transition, and `TASKS.md` BUDG-4 exists because that slope has
been reused across regimes before.

**Not included, and expected to subtract a few tenths of a kelvin.** Moving
absorption from the surface into a convecting troposphere at fixed
top-of-atmosphere flux is close to neutral for the global mean, but not exactly:
it stabilises the column and suppresses the turbulent fluxes. In the absorbing
aerosol literature the term runs -0.1 to -0.3 K per 10 W/m2 of atmospheric
absorption, which would be -0.2 to -0.5 K here. If the run comes in at the low
end of the bracket, this is the first place to look.

**The precipitation number is the atmospheric energy budget, not a guess.** The
atmosphere gains 15.0 W/m2 of radiative heating and must shed it; the only large
term available is latent heating, so precipitation falls. Full compensation is
-189 mm/yr, 19% of the total. The +2.5 K of warming pushes back at roughly
2-3% per kelvin, so the net is nearer -12 to -15%. Runoff is 15.5% of land
precipitation and amplifies a precipitation change by about 6.5x (`TASKS.md`
BUDG-1), so this is the largest single thing anything in this repository has
predicted for the carve criterion's denominator.

## What this does to the design mean, and to the flux

**The correction moves this world off its design temperature and the flux has to
be re-derived.** +2.5 K against a mean of 289.71 K is not a detail: it is
comparable to the 3.7-7.1 K that separates the vegetated and bare-rock endmembers
whose non-overlapping flux windows are why the orbit and the biosphere are one
choice. Restoring the design mean needs the flux ratio to fall by 0.0127, from
0.945 to about 0.932, using 196 K per unit flux ratio for the descending segment.
That is outside `sweep_flux_earth`'s lower point of 0.90 only in the sense that
it lands between the existing sweep points, so the bracket still spans it.

Sequence, and none of it is optional: apply, rebuild every binary, re-run the
baseline, re-derive the flux, then re-run again at the new flux. The carve list
taken on the current climatology is downstream of all of it.

## CO2: there is nothing to re-weight, and that is why it is a separate note

`TASKS.md` PHYS-1 asked for the CO2 shortwave absorptance to be re-weighted the
same way. **ExoPlaSim has no shortwave CO2 absorptance at all.** `radmod.f90`
carries CO2 only in `lwr`, from Sasamori (1968); `swr` has ozone in band 1 and
water vapour in band 2 and nothing else. Lacis and Hansen did not parameterise it
either, so this is faithful to the scheme rather than a bug in the port, and it
cannot be closed with a weight.

Adding it is a NEW ABSORBER: a namelist key, a column integral, a closed form, a
transmissivity term in band 2, and a decision about the water vapour overlap.
That work is PHYS-6, the patch is
`exoplasim/patches/exoplasim-3.4.2-co2-shortwave.patch`, and the derivation,
the overlap decision, the check it had to pass and its own prediction are in
**`exoplasim/notes/shortwave-co2.md`**. The one thing worth carrying here: the
CO2 correction has the SAME SIGN as this one and is about a sixth of its size,
so the two compound, and the flux re-derivation above has to be done with both
of them on.
