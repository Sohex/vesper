# ExoPlaSim has no shortwave CO2 absorptance, and this is what it costs

Derived 2026-08-18, on the `precarve-craton` baseline climatology.

`exoplasim/notes/shortwave-water-vapour.md` is the sibling of this note and the
place the machinery is explained; this one is CO2 and only CO2. The two
corrections have the SAME SIGN and this one is about a sixth of the size, so
they compound rather than cancel, and they are meant to land in one rebuild.

## The finding

`radmod.f90`'s `swr` carries ozone in band 1 and water vapour in band 2 and
nothing else. CO2 appears only in `lwr`, from Sasamori (1968). That is faithful
to Lacis and Hansen (1974), who did not parameterise the near-infrared CO2 bands
either, so this is an absorber that is MISSING rather than one that is wrong,
and it cannot be fixed by a weight the way ozone and water vapour were. It needs
a namelist key, a column integral, a closed form and a transmissivity term in
band 2, which is `exoplasim/patches/exoplasim-3.4.2-co2-shortwave.patch`.

It matters here and not on Earth because a K dwarf puts about 1.5x the Sun's
share of its flux into the bands CO2 absorbs in. That ratio is the whole reason
the term is worth adding to a model that has run without it for decades.

## The derivation, which is the water vapour one run on the other gas

`exoplasim/scripts/shortwave_band_weights.py` does both. Howard, Burch and
Williams (1956) measured total band absorption for CO2 in the same paper series
and the same form as for H2O, so the same integration gives the CO2 absorptance
against any spectrum: band absorptions from Howard's Table II, divided by band
width for Howard's own Eq. 11 band-average, weighted by the fraction of the
incident flux in each band, summed.

Three quantities have to be kept apart, and were not on the first pass.

**The column.** An atmos-cm is the depth the gas alone would occupy at STP, so
it is the column MASS over the STP density of CO2. The mass column of a trace
species is `vmr * (M_CO2 / M_air) * p_s / g`, because hydrostatic balance turns
TOTAL pressure into mass and not partial pressure. `radmod.f90:2395` carries
that molecular weight ratio as `zpv2pm` and `lwr` applies it. Omitting it
understates the column by 1.52, which divides straight out of a star-over-Sun
ratio and destroys every absolute quoted beside it. At Earth's gravity and this
world's mixing ratio the column is 360 atmos-cm.

**The gravity term.** At 12.81 m/s2 the same mixing ratio gives 272 atmos-cm
against Earth's 360, so this planet has about 24% less CO2 above it. That is the
same 1/g the Rayleigh term carries explicitly, and it works against the spectral
correction rather than with it.

**The pressure reduction.** `swr` does not evaluate an absorptance at the true
column. It reduces each layer's amount to an equivalent standard-pressure amount
by `sigma * ps / p0` and evaluates a fit stated at standard pressure, which is
what `radmod.f90` does for water vapour and what `lwr` does for its own CO2. For
a well-mixed gas that reduction is `sum(dsigma * sigma)`, a half on this grid.
The CO2 term uses it, or it would be the one absorber in the scheme evaluated on
a different kind of amount from every other.

## The 2.7 um overlap, which is a decision

Every CO2 band shares its interval with water vapour, and Lacis and Hansen
Eq. 21 already carries everything water vapour absorbs across the whole near
infrared, the CO2 intervals included. So CO2 must be charged only with what
water vapour leaves it, or the two gases claim the same photons twice.

Yamamoto (1962) handles this by dropping the 2.7 um CO2 band outright "because
of overlapping by the strong 2.7 um H2O band".

**The decision here is to KEEP the 2.7 um band and charge it, like every other
CO2 band, with the fraction of its interval water vapour leaves.** Three
reasons, in order of weight:

1. The overlap is COMPUTED band by band rather than argued about. Yamamoto drops
   the band because he is working by hand with a single strong overlapping
   absorber and has no cheap way to price the residue; the same integration that
   gives the CO2 absorptance gives the H2O absorptance across the same interval,
   so the residue is available.
2. `CLAUDE.md`'s physics-is-not-a-knob rule. The band absorbs, water vapour does
   not close it completely, and removing a term because it is inconvenient to
   price is exactly the move that rule forbids.
3. The decision does not decide the check. The Earth-column check below passes
   with the band in and with it out, so keeping it is not what makes the answer
   look right.

What it is worth is reported rather than hidden: dropping the band, Yamamoto's
own choice, is the low end of the bracket, and it moves the absorptance by about
a sixth and the star-over-Sun weight by less than 0.002. The weight is
insensitive to the decision; the absolute is not.

**Two things the overlap treatment does not do**, both stated because both are
places the term is known to be approximate:

- The overlap factor is evaluated at ONE water path, the baseline's own, and
  baked into the fitted constants. A drier column than the baseline gets too
  little CO2 absorption and a wetter one too much. It is the same class of
  approximation as quoting `h2osww` at one water path, and the direction is that
  the term is understated aloft, where there is less water above than the
  baseline column implies.
- The 15 um band's overlap factor is 1.0 by omission, because Howard's
  near-infrared H2O set does not include the pure rotation band that covers it.
  That band carries 3% of the total, and the direction is to overstate it.

Because the overlap is already removed inside the absorptance, the CO2
transmissivity is ADDED into `ztb2` beside water vapour's rather than multiplied
into it, matching how ozone, water vapour and the aerosol are already booked
there. Multiplying the two transmissions would remove the overlap a second time.

## The check that could have failed

Earth's near-infrared CO2 solar absorption is measured at 1.5 to 2.5 W/m2. Run
the same integration at Earth's column, Earth's gravity and Earth's mean
insolation and it has to land inside that, and nothing in the derivation was
tuned to it.

| | absorptance | W/m2 |
| --- | ---: | ---: |
| solar-weighted, Earth's column and Earth's insolation | 0.006011 | **2.05** |
| the same, dropping the 2.7 um band | 0.005137 | **1.75** |
| this star, this planet's column and insolation | 0.008359 | **2.69** |
| the same, dropping the 2.7 um band | 0.007079 | **2.28** |

It passes, both ways, which is what licenses the number for this star.

## An independent check, and the one place it disagrees

`exoplasim/notes/corrk-cross-check.md` runs the same quantity from correlated-k
tables built on HITRAN2020, at the same paths and against the same two spectra.
`co2sww` agrees to 0.03% and the absorptance to 8%, which is inside the check's
own temperature-grid uncertainty, and 1.92 W/m2 lands between the two ends of
the 2.7 um bracket below.

**The per-band attribution does not agree, and reason 1 above overstates what the
overlap treatment can do.** It is computed band by band, but at Howard's band-mean
resolution, which spreads saturation out of the water band cores and into the
wings. Correlated-k puts the 2.7 um clear fraction at 0.003 rather than 0.171,
which is Yamamoto's verdict, and the 2.0 um one at 0.718 rather than 0.479. Those
two errors have opposite signs and nearly cancel, so reason 3 survives and the
total is right for a compensating reason. Do not correct one band alone.

## The weight

    co2_sw_weight = 1.510          bracket, from the 2.7 um decision, 1.508 to 1.510

**This weight is far more robust than the water vapour one, and that is worth
knowing before anyone brackets it.** The star-over-Sun ratio runs 1.494 to 1.517
across four decades of CO2 column, and 1.508 to 1.510 across two decades of water
path. The reason is that CO2's bands are all long-wavelength, so they all sit in
the part of the spectrum where the K dwarf's boost is large and similar, whereas
water vapour's span from 0.72 to 6.3 um and the boost across them runs 0.99 to
1.53. Nothing about the CO2 weight depends on where in the column it is
evaluated.

The value belongs to the k25v spectrum. Against the 4965 K blackbody `solarini`
builds when `NSTARFILE` is 0 the weight is 1.364, because a blackbody has no line
blanketing pushing flux out of the blue and into the near infrared. Setting 1.510
while the model runs a blackbody over-corrects by a tenth, which is the same trap
`h2osww` has.

## What the model codes, and why it is a fit

`swr` cannot call an integration, so the solar-weighted absorptance is fitted to
a closed form in the absorber amount u, in atmos-cm:

    A(u) = a1 ln(1 + b1 u) + a2 ln(1 + b2 u)

with `co2sww` re-weighting the result for the host exactly as `h2osww`
re-weights Eq. 21.

**There are two fits to that form and only one of them runs.** The Howard fit is
the one `shortwave_band_weights.py` can make, because Howard's band set is the
only absorption data that script has; it is what the patch header codes and what
`corrk_cross_check.py --fit` compares against, and it is the artifact's
`co2.closed_form_fit`. PHYS-10 refitted the same form, over the same range and
by the same protocol, to HITRAN2020 through the Generic PCM correlated-k tables,
and that is what `radmod.f90` carries. At this planet's CO2 path the Howard fit
is 7.7% stronger. The artifact records the running coefficients separately, in
`co2.closed_form_fit_in_radmod`, read out of the model source rather than
restated, so the two cannot silently diverge again; the argument for the refit
is `exoplasim/notes/corrk-cross-check.md`.

Over 1 to 1e4 atmos-cm, which the model never leaves, the Howard fit is within
4% of its own integration at worst and 1.2% rms; the line-list fit is 4.6% over
the 100 to 1000 atmos-cm a T42 column occupies and worse only at the far end of
the range, where no column goes. On a 2.7 W/m2 term either residual is at most a
tenth of a W/m2, smaller than the 2.7 um bracket and far smaller than the H2O
reconstruction's own level error.

Lacis and Hansen's Eq. 21 form was tried first and fits this curve worse while
wanting a negative coefficient in its denominator, which can go singular on a
column nothing in the scheme forbids. Two logarithms are the shape Howard's own
strong-band fit has, stay positive and monotone everywhere, and cost two LOGs
per layer.

`co2sww` itself is untouched by the refit: it is a star-over-Sun RATIO and the
two derivations agree on it to 0.03%.

## The prediction, made before the run

Stated so the baseline re-run can falsify it, priced against the baseline
climatology rather than against Earth, and priced on the Howard fit `swr`
carried when it was made. PHYS-10's refit takes every row 7.2% weaker; the A3
bundle carries that as its own row rather than restating this table, because the
two changes have to be attributable separately.

| quantity | change at co2sww = 1.510 |
| --- | ---: |
| atmospheric shortwave absorption, `rst - rss` | **+2.65 W/m2** |
| surface net shortwave, `rss` | **-2.30 W/m2** |
| top-of-atmosphere net shortwave, `rst` | **+0.51 W/m2**, bracket +0.35 to +0.79 |
| mean surface temperature, `ts` | **+0.44 K** |

**The top-of-atmosphere central value is the weakest number here** and it should
be read as nearer the top of its bracket than quoted. The split between the
bracket's ends is 65/35 in favour of the below-cloud albedo, which was chosen for
water vapour because water vapour is bottom-heavy. CO2 is well mixed, so more of
this absorption sits above the cloud than water vapour's does, and the true split
is less bottom-heavy. The bracket ends are the honest statement; the central
value is inherited.

**The same sign as the water vapour correction, and about a sixth of the size.**
Both move atmospheric absorption up, surface shortwave down and the mean
temperature up, so they compound. Nothing cancels, and the flux re-derivation
that follows the water vapour correction has to be done with this term already
on, not after it.

## What has to happen next

Superseded by the A3 bundle: predictions and arms are in
`exoplasim/notes/forcing-bundle-predictions.md`; the bundle total is checked
against the summed predictions and bisected only on disagreement. What
remains:

1. The attribution segments `exoplasim/notes/shortwave-water-vapour.md` asks for,
   with `co2sww` off and then on, before the long converged run. Three physics
   changes in one rebuild and no result can be attributed to any of them.
2. Re-derive the flux with both corrections on.
