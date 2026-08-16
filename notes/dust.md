# Dust: why the model cannot help, and why it matters anyway

## What ExoPlaSim has

Aerosols are **off**: `L_AERO = 0` in every run this project has made. An
`aero_namelist` is written into each run directory and is inert.

If enabled, `l_source = 2` selects dust, and the source term is
`aerocore.f90:625`:

    mmr(:,:,NLEV,ic) = fcoeff * land

A fixed mixing ratio in the bottom layer of every land gridbox, zero over ocean,
then advected. There is no wind-speed threshold, no soil-moisture gating, and no
weighting by surface type. `l_source = 1`, the default, is photochemical haze
rather than dust.

Real dust emission is violently nonlinear in precisely the terms this omits. It
scales steeply above a threshold friction velocity, shuts off when the surface is
wet or vegetated, and concentrates in a handful of source regions: the Bodele
Depression, a dry lake bed, supplies a large share of Earth's dust by itself.

**So enabling it here would add a bias rather than remove one.** On this planet
closed-basin fill is about a quarter of the land, and a uniform-over-land source
would emit as much from forest as from salt pan -- getting wrong the exact
distinction that makes dust interesting here. Representing it properly needs a
supplied source field, which means patching `aerocore` rather than setting a
namelist value.

## Why it matters anyway, and a correction

Dust is the return leg of the phosphorus cycle, and this project has so far only
described the outbound one.

The hydrological story is that phosphorus weathered from rock is carried by
runoff, and in a closed basin it never reaches the sea -- so endorheic drainage
traps P on basin floors and starves the uplands. That was stated here as though
it were the whole account. It is not.

Aeolian transport runs the other way. Dry lake beds are the most efficient dust
sources there are, and dust is a major phosphorus input to old, deeply weathered
soils whose own apatite is exhausted -- which is exactly what the wet parts of
this planet should be. On Earth the Amazon receives a significant share of its
phosphorus from Saharan dust, and the dominant single source is a dry lake bed.

So the loop is: **water moves phosphorus into the basins, wind moves it back
out.** The net depends on the balance of the two, and this world should sit
further along that loop than Earth does, because its playa and salt-crust
fraction is roughly an order of magnitude larger.

That also weakens, in the same direction, the conclusion that this planet should
show P-poor uplands against P-rich basin floors. The gradient is probably real
but shallower than the hydrology alone implies.

## Priced before deciding

The budget's own stopping rule says to refine an input when its plausible range
exceeds the effect of the last thing refined. So, a crude bound before building
anything, using scattering-aerosol forcing
`F = -tau * S/4 * beta * (1 - a_surface)^2 * land_fraction` at an upscatter
fraction of 0.20:

| dust AOD over land | over vegetation (a=0.20) | over playa (a=0.40) |
| ---: | ---: | ---: |
| 0.03, Earth-like | -0.46 K | -0.26 K |
| 0.10 | -1.53 K | -0.86 K |
| 0.30 | -4.60 K | -2.59 K |

**Corrected: the first version of this bound used the wrong source area.** It
took closed-basin fill at 26.6% of land, which is the figure for the *uncarved*
build -- the hyper-arid limit, where no basin has had its outlet cut. That is not
this world's expected configuration.

Two corrections, both reducing it:

- **Carving removes fill.** Two carve passes took it from 26.6% to 12.4% of land.
  The uncarved figure describes a planet the drainage verdict says does not
  survive.
- **Wet fill is not a dust source.** On the build where a lake solution exists,
  4.4 points of the 16.5% sit under water, leaving **12.0% of land as dry
  closed-basin fill**. A flooded playa emits nothing.

So the realistic dry-source area is around 10 to 12% of land. Earth's
*preferential* dust sources -- topographic lows holding fine sediment, not
deserts in general -- are a few percent of land, so this world is perhaps two to
four times Earth, not ten.

| dust AOD over land | over vegetation (a=0.20) | over playa (a=0.40) |
| ---: | ---: | ---: |
| 0.03, Earth-like | -0.46 K | -0.26 K |
| 0.06, ~2x Earth source | -0.92 K | -0.52 K |
| 0.12, ~4x Earth source | -1.84 K | -1.03 K |

That is a smaller claim than the first version made. Dust is plausibly worth half
a kelvin to two kelvin, which puts it above the lake compositing and the
lithology fixes but **well below the biosphere question**, and comparable to the
structural items rather than dominating them.

**And the sign is not uniform.** Dust is darker than salt crust and brighter than
vegetation, so it warms over its own source regions and cools elsewhere. That is
the decisive argument for modelling rather than declaring: a spatially varying
sign cannot be recorded as a bias direction the way the missing ocean heat
transport can. Either it is computed or it is unknown.

## The bracket above does not span its own sign statement

Measured 2026-08-16, and it corrects the two tables above rather than extending
them.

Aerosol top-of-atmosphere forcing reverses sign at a critical surface albedo
`a_c`: below it the layer scatters and cools, above it the layer is darker than
the ground and warms. Both tables above bracket the bright case at a = 0.40 and
return cooling, so the bracket never crosses zero even though the paragraph
immediately above says the sign is not uniform. The prose was right and the
numbers did not implement it.

Computed here with a Bohren-Huffman Mie code, integrated over the Balkanski
source distribution (number-median radius 0.295 um, sigma 2.0, density 2.6) and
band-averaged on the k25v spectrum with correct bin-width weighting:

| indices | band (um) | ssa | g | a_c |
| --- | --- | --- | --- | --- |
| Di Biagio 2019, measured | 0.34-0.75 | 0.965 | 0.700 | 0.504 |
| OPAC, absorbing end | 0.34-0.75 | 0.889 | 0.716 | 0.289 |
| OPAC, absorbing end | 0.75-4.00 | 0.951 | 0.685 | 0.456 |

Salt crust is 0.40 to 0.50. **The reversal sits inside our own stated albedo
range for the surface that makes this world unusual**, and band 2 carries 61.6%
of the stellar flux (the model's own figure, not a blackbody estimate). Over
ocean at 0.07 and vegetated land at 0.18 every case cools without ambiguity.

The Mie code is validated against two independent targets before any of this was
believed: Bohren and Huffman's worked example (m = 1.55, x = 5.213, Qsca =
3.10543) to 4.5e-6, and OPAC's four published mineral components across four
wavelengths, reproducing mass extinction efficiency, single-scattering albedo
and asymmetry parameter to three or four decimals.

Rocha-Lima et al. (2018) pushes the same way and is not yet usable. Its headline
result is that fine-mode k **rises** from a minimum near 650 nm through the
shortwave infrared, which would lower `a_c` further and favour warming, and it
also means OPAC may not be an over-absorbing endmember in band 2 at all. But its
values above 950 nm exist only in a figure, its real part is an assumed constant
of 1.56 rather than a retrieval, and its data are available only on request.

**The consequence for the plan is that this does not need a GCM.** The sign is an
analytic function of ssa, g and the surface albedo map we already build, and done
offline the longwave term can be included, which ExoPlaSim's shortwave-only
aerosol block cannot do at all. What is missing is band-2 refractive indices,
not model machinery.

## Is the patch error-prone?

Less than it looks, because the parameters split cleanly.

The optics are **not** haze-locked: `l_aerorad = 1` reads extinction, scattering,
backscatter and asymmetry for both bands from a supplied `aerofile`, so mineral
dust is representable. Three of the four parameters are Earth-calibratable --
`apart` (dust effective radius, 1 to 2 microns, against a haze default of 5 nm),
`rhop` (2500 kg/m3 for mineral dust, against a default of 1000, which is water),
and the `aerofile` itself, band-averaged from published dust refractive indices
for this star's bands.

The fourth, `fcoeff`, sets the total loading and has no Vesper calibration. It
should be **bracketed rather than chosen**, anchored so that an Earth-like source
area gives Earth's dust AOD and then spanned across the playa-fraction ratio.
That is the same treatment `pedogenesis.yaml` gives every constant it cannot
verify.

Two simplifications that remain and should be stated wherever a number is
quoted: the scheme is monodisperse where real dust spans 0.1 to 20 microns, and
the source term has no wind-speed threshold even once it is weighted by surface
type.

## What to do

Not in the GCM. Model the aeolian phosphorus term the way the carve verdict and
the thermostat efficiency are modelled: outside the climate model, from its
climatology. Everything needed is already here -- the playa and salt-crust
distribution from lithology, surface winds and precipitation from the
climatology, and the basin-to-cell coupling matrix that already integrates
climate over catchments.

Two things that would need declaring rather than deriving, and should be
bracketed rather than guessed: the phosphorus content of basin fill, and an
emission efficiency. Both are Earth calibrations of the kind
`pedogenesis.yaml` already carries.

Separately unquantified: dust has a radiative effect as well as a chemical one.
Bright dust over dark ground cools and over ice warms, and a quarter of this
planet's land is a potential source. That belongs in the error budget as an
unpriced item rather than being assumed negligible.
