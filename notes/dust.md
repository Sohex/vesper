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

Computed by `exoplasim/scripts/dust_optics.py`, which is the authority for
these numbers; `analysis/dust_optics.json` is its product and carries the input
hashes. A Bohren-Huffman Mie code integrated over the Balkanski source
distribution (number-median radius 0.295 um, sigma 2.0, density 2.6) and
band-averaged on the k25v spectrum with correct bin-width weighting:

| indices | band (um) | ssa | g | a_c |
| --- | --- | --- | --- | --- |
| Di Biagio 2019, measured | 0.34-0.75 | 0.967 | 0.698 | 0.515 |
| OPAC, absorbing end | 0.34-0.75 | 0.895 | 0.714 | 0.300 |
| OPAC, absorbing end | 0.75-4.00 | 0.952 | 0.686 | 0.456 |

Salt crust is 0.40 to 0.50. **The reversal sits inside our own stated albedo
range for the surface that makes this world unusual**, and band 2 carries 61.6%
of the stellar flux (the model's own figure, not a blackbody estimate). Over
ocean at 0.07 and vegetated land at 0.18 every case cools without ambiguity.

The Mie code is validated against two independent targets before any of this was
believed: Bohren and Huffman's worked example (m = 1.55, x = 5.213, Qsca =
3.10543) to 4.5e-6, and OPAC's four published mineral components across four
wavelengths, reproducing mass extinction efficiency, single-scattering albedo
and asymmetry parameter to three or four decimals.

**Resolved, by digitising Rocha-Lima et al. (2018) Fig. 10.** Its values above
950 nm exist only as a plot, so they were extracted by colour-thresholding the
marker series and mapping through the axis frame -- reproducibly, by
`dust_optics.py --digitize`, with the extracted points committed as
`exoplasim/data/dust/rochalima2018_fig10_digitized.csv` -- then checked against
the four values the paper states in text: mixed-mode Algeria reads 0.0032 at 450 nm
against a stated 0.0030, and 0.0005 at 850 nm against 0.0005. Accuracy is about
+/-0.0005 in k.

Fine-mode k does rise from its ~650 nm minimum through the shortwave infrared,
about fivefold by 2450 nm. **That does not make band 2 more absorbing**, because
absorption efficiency goes as the size parameter times k, i.e. as k/lambda, and
lambda grows faster than k does: k/lambda is 0.0035 at 0.52 um and 0.0028 at
1.25 um. Band-2 single-scattering albedo therefore comes out *higher* than band
1, not lower.

| band 2 indices | ssa | g | a_c |
| --- | --- | --- | --- |
| Rocha-Lima Algeria fine | 0.975 | 0.684 | 0.570 |
| Rocha-Lima Mauritania fine | 0.964 | 0.687 | 0.508 |
| OPAC | 0.952 | 0.691 | 0.453 |

So `a_c` is 0.51-0.57 in band 2 against 0.50 in band 1, and salt crust tops out
at 0.50. **Dust cools over every surface on this world**, including closed-basin
fill, narrowly at the bright end. The evaporite mineralogy of the actual source
pushes the same way, being less absorbing than the Saharan silicate dust
measured here.

Two corrections this forced. OPAC is **1.0 to 2.1x** more absorbing than measured
in the near-IR, not the 8-22x inferred earlier from the ARIA Peters file --
Peters is the outlier, running 3-10x below Rocha-Lima. And the two-sided entry
this section previously argued for was an overcorrection: the sign is one-signed
cooling. What stands is the narrower point that the original bracket was
one-signed *by construction* and could not have shown a reversal had one existed.

Remaining caveats: k digitised from a figure; n taken as 1.52 (Di Biagio
measured) where Rocha-Lima assumed 1.56; k held flat above 2.45 um, which
carries about 5% of the flux; and `a_c` depends on the assumed size
distribution.

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

## What to do, revisited 2026-08-17

Reopened before the baseline run, because the first decision was taken when
nothing downstream needed a deposition field and three things now do: loess
(SURF-3), the phosphorus return leg above, and the derived-surface classifier.
A second opinion argued for fixing ExoPlaSim's scheme properly in Fortran --
Kok (2014) brittle-fragmentation emission, a subgrid wind PDF, wet scavenging,
and four to six size bins. The physics in that list is right. The conclusion
here is still not to do it in the GCM, for reasons that are specific to this
world rather than general.

### What the built-in scheme is, and why nobody should enable it

Unchanged and now agreed from both directions: `fcoeff * land` is a boundary
condition, not an emission scheme. It would emit as much from forest as from
salt pan, which is exactly the distinction that makes dust interesting here.

### Why not fix it in place

**The sign is already resolved, offline, and that was the expensive part.**
`dust_optics.py` settles it with a Mie code validated against Bohren-Huffman and
OPAC: dust cools over every surface on this world, narrowly including salt
crust, because the critical albedo is 0.51-0.57 in band 2 against a crust that
tops out at 0.50. An interactive scheme would not buy the sign. It would buy the
feedback, which is a smaller quantity.

**ExoPlaSim's radiation caps the return, and in the direction that matters.**
Its aerosol block is shortwave only, and dust's longwave warming is on the order
of a third of the net effect on Earth. So an interactive implementation would
overestimate the cooling by roughly that much, and fixing it means patching
`radmod` rather than `aerocore`. Offline, the longwave term is arithmetic.

**The emitting fraction is not the playa fraction, and this world is the case
where that distinction bites.** Kok's brittle fragmentation is calibrated on
silicate soils. A cemented salt crust is not a silicate soil and emits far less
than loose clastics; the Bodele is efficient because it is diatomite, not
because it is a salt pan. This project already separates `evaporite` from
`playa_clastic`, derives the ephemeral crust from the lake solution, and knows
which basins hold water. Those fields are the source map, and none of them exist
inside the GCM.

**The source area is at its maximum on a pre-carve terrain.** Carving removes
fill, which is why the estimate above was corrected from 26.6% to about 12% of
land. Calibrating an emission scheme now would calibrate it against a landscape
the carve verdict says does not survive.

### What to build instead

An offline dust component, on the same footing as the carve verdict and the
thermostat efficiency: computed from a climatology rather than inside the model.
Every input it needs already exists.

| term | from |
| --- | --- |
| friction velocity, and a subgrid wind distribution over the gridbox | `spd`, `ua`, `va` in the climatology; emission goes as roughly u* cubed above a threshold, so a gridbox-mean wind emits almost nothing and a Weibull over the box is not optional |
| threshold, moisture-gated | `mrso` from the climatology |
| drag partitioning | the z0 field `build_surface_roughness.py` already writes, 0.0036 to 11.1 m |
| erodible fraction | lithology, clay from pedology, erodibility per region |
| what is not a source | the lake solution, the derived ephemeral crust, and vegetation cover |
| transport | the climatology's winds |
| removal | gravitational settling, plus below-cloud scavenging from `prc` and `prl`, without which fine-mode lifetime is out by an order of magnitude |

Products: a deposition flux, which is what loess, the phosphorus return leg and
the minerals overlay all want, and an aerosol optical depth map, which with the
optics already computed gives the radiative forcing including longwave.

### When it enters the climate, and the test for whether it should

**Not in this baseline.** The order is forced rather than chosen: a dust model
reads winds, soil moisture and precipitation, so it cannot precede the
climatology it is driven by. The baseline runs without dust and the forcing is
carried as a bias with a number rather than as an unpriced item.

**Fed forward as a prescribed field, not as interactive dust.** The next
iteration's run can take a supplied aerosol distribution; that captures the
spatial pattern, which is the part the built-in scheme gets wrong, without
writing coupled emission physics in Fortran.

**The threshold for reopening, fixed here before the answer is known.** If the
offline land-mean dust optical depth exceeds **0.10**, or the forcing exceeds
**1.5 W/m2** in the global mean, the feedback is large enough that a prescribed
field is no longer defensible and the emission scheme belongs in the model. Both
numbers are chosen against the pricing above, where 0.12 was the four-times-Earth
case worth about 1.8 K over vegetation.

Below that, prescribed is the answer and the Generic PCM is not: it carries all
of this properly and costs two orders of magnitude in runtime plus a pipeline
rewrite, which is not a trade worth making for a term whose sign is settled and
whose source area is uncertain by a factor of two.
