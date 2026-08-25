# A zenith term for the modelled snow albedo: form, coefficients, price

**Derived:** 2026-08-25, by `analysis/snow_albedo_zenith.py`, from
`config/planet.yaml`, the model's own star-weighted albedo constants and the
spectral snow reflectance in
`references/exocam/tools/spectral_albedos/snow100um.txt` (Shields et al. 2013,
shipped with ExoCAM). No World Orogen generation and no ExoPlaSim run: every
number here is geometry, a spectral integral, or the project's own
flux-to-kelvin conversion.

This is worldbuilding. Vesper is an invented super-Earth around a mid-K dwarf;
every quantity below is a property of that planet's simulated snow surface or
of the shortwave scheme that reflects light off it.

WORLD-SPD asked for four things and this note is them: a functional form with
its provenance, coefficients derived against the configured star rather than
inherited, the effect priced in kelvin or declared unpriceable with the reason,
and a verdict on whether the open-ocean zenith formula already in `radmod` is
the right form to reuse. It is an EXPLORATION product. Nothing in the model was
changed.

## The gap, against the code as it stands

`landmod.f90` sets the snow albedo as a linear ramp in SURFACE TEMPERATURE
between 263.16 K and `tmelt`, from `dsnowalbmn` to `dsnowalbmx`, blended toward
`albforest` by the snow-canopy mask. The broadband branch and both banded
branches carry the same ramp and none of the three carries a solar zenith
angle. `radmod.f90:swr` then overwrites the surface albedo for ICE-FREE OCEAN
with `AMIN1(0.05/(zmu0+0.15), 0.15)` under `necham`, which is a zenith
dependence and is on by default.

So inside one model, in one low-sun regime, the same physical predictor is
carried for one surface and absent for another. That is the internal
inconsistency the row is built on, and it is still true after the wave-1
changes to `radmod`.

## The form, and why it is this one rather than a fit

For a semi-infinite, strongly-scattering medium the asymptotic result of
radiative transfer is that the direct-beam albedo is the diffuse albedo raised
to the escape function of the incidence direction:

    alpha(mu0) = alpha_diffuse ** K(mu0),    K(mu0) = (3/7) * (1 + 2*mu0)

`mu0` is the cosine of the solar zenith angle. The normalisation is not a fit:
`2 * integral(K(mu) mu dmu) = 1` over the hemisphere, and the value of `mu0` at
which `K = 1` is exactly `2/3`. So the direct-beam albedo equals the diffuse
albedo at `mu0 = 2/3` and nowhere else, and A MODEL THAT CARRIES NO ZENITH
ANGLE IS IMPLICITLY QUOTING THE `mu0 = 2/3` VALUE. Adding the term is a
redistribution about the number the model already has, not a change to it.

Four properties make this the form to take rather than ClimaLand's
`alpha_0 + d_alpha*exp(-k*cos_z)`:

- **It has no coefficients to inherit.** ClimaLand's `alpha_0`, `d_alpha` and
  `k` are Earth-solar-spectrum calibrations, and this project's rule on
  undocumented values makes them a test case rather than an anchor. This form
  has one input, `alpha_diffuse`, and the model already computes it.
- **It is per band by construction.** The zenith dependence lives entirely in
  the absorption, so it is small where ice barely absorbs and large where it
  does. A broadband exponential fit cannot express that and therefore cannot be
  re-weighted for another star; this form re-weights by construction, in the
  same two bands `radmod` already carries.
- **It composes with what is there.** It acts on `alpha_diffuse`, so the
  temperature ramp and the forest blend keep their meaning: they set
  `alpha_diffuse` and the zenith term redistributes about it.
- **It is bounded by its own shape.** `alpha ** K` stays in `[0, 1]` for any
  albedo and any `mu0`, so it needs no clip. The linear co-albedo form
  `1 - alpha(mu0) = (1 - alpha_diffuse) * K(mu0)` is the same result to first
  order and is what the literature usually quotes, but it goes negative on the
  warm end of this model's own ramp, where the band-2 co-albedo reaches 0.73.
  The exponential form is preferred here for that reason and agrees with the
  linear one to within 0.002 in band 1.

Its limit is stated rather than hidden: the asymptotic result assumes a
semi-infinite pack, so it does not describe thin snow over a dark substrate,
and the amplitude it gives is an upper bound where absorption is strong.

## Coefficients, on this star, from the model's own constants

`radmod` already carries the ECOSTRESS surface blends over 965 wavelengths,
weights each by the configured spectrum and splits at the model's band
boundary; `analysis/ice_albedo.py` verifies that integration against what
`radini` printed. So the per-band diffuse albedo this form needs is not
something to derive a second time -- it is the model's own finalized
`dsnowalb*`, and taking it from there is what keeps this a composition rather
than a second opinion.

At the configured star, 38.24% of the flux falls in band 1 against the Sun's
51.70%.

| model constant | band 1 | band 2 | broadband | `mu0 = 1` | `mu0 = 0` | zenith span |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `dsnowalbmx`, cold snow | 0.9828 | 0.5853 | 0.7373 | 0.6841 | 0.8705 | **0.186** |
| `dsnowalb`, base | 0.7638 | 0.4061 | 0.5429 | 0.4643 | 0.7604 | **0.296** |
| `dsnowalbmn`, warm snow | 0.5078 | 0.2726 | 0.3626 | 0.2761 | 0.6398 | **0.364** |

**THE MISSING PREDICTOR IS WORTH AS MUCH AS THE ONE THE MODEL HAS.** The whole
range the temperature ramp can move the broadband snow albedo over is
`0.7373 - 0.3626 = 0.375`. The zenith span is half of that on cold snow and
equals it on warm snow. Two predictors of comparable weight, one of them
absent.

### The star makes it bigger, not smaller

An independent snow reflectance, weighted the same way, isolates the stellar
half of that:

| weighting of the ExoCAM snow spectrum | band 1 | band 2 | broadband | zenith span |
| --- | ---: | ---: | ---: | ---: |
| the configured star | 0.9800 | 0.5137 | 0.6920 | 0.208 |
| 5772 K Planck | 0.9830 | 0.5790 | 0.7878 | 0.150 |

5.3% of the star's flux falls outside the measured reflectance interval and the
endpoint value is held there; 3.6% for the Planck. The two datasets agree on
band 1 to 0.003 and differ by 0.07 in band 2, which is the size of the
disagreement between two spectral snow datasets and is smaller than the effect
being sized.

A redder star puts more of its flux where ice absorbs -- 61.8% above the band
boundary against the Sun's 48.3%, and further into the near infrared within
that band -- so snow under it is both darker and more zenith-sensitive. The
span is 1.39 times the Earth-weighted one. **The transfer argument runs the
opposite way from the usual caution: this is not an Earth effect that may be
smaller here, it is an Earth effect that is larger here.**

## What it does to the simulated surface, at this obliquity

Flux-weighted daily means at 32 degrees obliquity, on cold snow
(`dsnowalbmx`, broadband diffuse 0.7373), which is the polar case and the one
where the form is most defensible. `delta` is against the temperature-only
albedo the model uses now, which the form identifies as the `mu0 = 2/3` value.

| latitude | season | flux-weighted albedo | delta | daily insolation W/m2 | change in absorbed surface shortwave W/m2 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 0 | equinox | 0.7193 | -0.018 | 409.4 | +7.36 |
| 15 | winter | 0.7610 | +0.024 | 251.9 | -5.97 |
| 30 | summer | 0.7208 | -0.017 | 490.9 | +8.10 |
| 30 | winter | 0.7918 | +0.055 | 150.1 | -8.19 |
| 45 | equinox | 0.7579 | +0.021 | 289.5 | -5.97 |
| 45 | winter | 0.8309 | +0.094 | 54.2 | -5.08 |
| 60 | equinox | 0.7878 | +0.051 | 204.7 | -10.35 |
| 60 | summer | 0.7414 | +0.004 | 590.2 | -2.41 |
| 75 | equinox | 0.8257 | +0.089 | 106.0 | -9.37 |
| 75 | summer | 0.7571 | +0.020 | 658.3 | -13.05 |

At 32 degrees obliquity, latitude 60 and 75 in winter are polar night and carry
no flux to weight.

**The sign flips with latitude, which is what a global mean would destroy.**
The simulated snow absorbs MORE than the model has it within about 40 degrees
of the equator and LESS everywhere poleward, and the poleward half is where the
snow is. Annual and area-weighted over a polar cap, per unit snow-covered area:

| cap | annual mean change in absorbed surface shortwave |
| --- | ---: |
| poleward of 45 | -6.30 W/m2 |
| poleward of 60 | -7.97 W/m2 |
| poleward of 70 | -9.00 W/m2 |

## The price, bracketed rather than asserted

The global-mean price needs the snow-covered area fraction and its seasonal
distribution, which is a climatology, and the canonical climatology lineage
does not exist yet. So it is priced as a BRACKET over that fraction, and it is
an UPPER bound in two ways: it treats a surface shortwave change as a
top-of-atmosphere forcing, which ignores the atmosphere's absorption of both
beams, and the asymptotic form overstates the amplitude where absorption is
strong.

At the poleward-of-60 value of -7.97 W/m2 per unit snow area, through
`lib/sensitivity.forcing_to_kelvin`:

| snow-covered share of the planet | global forcing | equilibrium at planetary albedo 0.30 |
| ---: | ---: | ---: |
| 2% | -0.16 W/m2 | -0.14 K |
| 5% | -0.40 W/m2 | -0.34 K |
| 10% | -0.80 W/m2 | -0.68 K |
| 20% | -1.59 W/m2 | -1.35 K |

Between planetary albedos of 0.25 and 0.35 those kelvin figures move by about
8% either way, an order of magnitude less than the snow-area bracket. **The
term is not negligible at any snow cover this world could plausibly carry, and
it COOLS**, which is the direction that matters for a world whose habitability
question is a cold one. It is also a positive feedback on the snow itself:
colder means more snow means more of the planet under the term.

## Is `radmod`'s open-ocean formula the right form to reuse?

**No. Wrong form, right location.**

`AMIN1(0.05/(zmu0+0.15), 0.15)` parameterises FRESNEL reflection at a flat
dielectric interface. Three things follow, and each disqualifies it for snow:

- It REPLACES the surface albedo rather than acting on it. Reusing it for snow
  would set the simulated snow albedo to at most 0.15 at every zenith angle,
  darker than the model's bare ground.
- It saturates, and the cap is there because the expression diverges at grazing
  incidence and has to be stopped by hand. A volume-scattering surface has no
  such divergence and needs no such cap; `alpha ** K(mu0)` is bounded by its
  own shape.
- It carries no surface property at all. There is nowhere in it for a grain
  size, a density or a spectral band to enter, so it cannot be re-weighted for
  another star and cannot express the band contrast that IS the snow effect.

The two share a predictor and nothing else: for water the zenith dependence is
a reflection at the surface, for snow it is the escape probability of a photon
that has already entered the medium. Same variable, different physics.

What the open-ocean case does establish is where the term belongs.
`radmod.f90:swr` is the only place in the model where `zmu0` and the surface
albedo are both in scope, and it is already the place where one surface's
albedo is overwritten with a zenith form. `landmod` has no zenith angle at all.

Two non-local facts an implementation has to face, both checked against the
source rather than assumed:

- `radmod` uses `pumamod`, so `dsnow` is visible to it. `snowcovz`, the depth
  at which half a cell is snow-covered, is LOCAL to `landmod` and is not, so
  the snow-covered fraction has to be passed rather than recomputed. A second
  derivation of that fraction is failure class 17.
- `landmod` has already applied the forest blend and the snow-cover fraction by
  the time `radmod` runs. The zenith term acts on the diffuse albedo, so it
  composes correctly only if it is applied to the BLENDED snow value over the
  snow-covered fraction. Applying it to the unblended snow endpoint would
  double-count the canopy.

## What this does not settle

- **The density half.** The grain-size proxy is snow density, a constant in
  this model with no compaction. That is GRAV-8, and it is the prerequisite for
  a density factor, not for the zenith factor: the form above needs only
  `alpha_diffuse`, and the temperature ramp supplies it today.
- **Sea ice and glacier ice.** Both carry the same temperature-only ramp and
  the same absent predictor, `dicealbmn`/`dicealbmx` and `dglacalbmn` are
  finalized on the same footing, and
  `references/exocam/tools/spectral_albedos` ships blue marine ice and a 50/50
  mixture. The same integral runs for them unchanged. Only snow was asked for
  here.
- **The strong-absorption amplitude.** Closing the gap where the asymptotic
  form overstates needs a two-stream or delta-Eddington solve on ice optical
  constants. `lib/mie_dust.py` already carries the machinery; what is missing
  is the ice refractive index over the two bands.
- **The diffuse beam.** The form is for the direct beam. The model's `zmu00`
  is a fixed diffusivity factor of 0.5 for scattered light, and under a cloudy
  sky the surface sees mostly diffuse radiation, which has no zenith angle. An
  implementation has to weight the term by the direct fraction or it will
  overstate the effect under cloud.
