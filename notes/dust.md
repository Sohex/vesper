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
of the stellar flux (`lib/stellar.py`, the same integration `solarini` does). Over
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

Below that, prescribed would have been the answer.

## The test fired, and this is what was decided

Measured 2026-08-17, once DUST-5 replaced a wind tail fitted to 32 snapshots with
one measured from 1,463 three-hourly samples: land-mean optical depth **0.376** at
the central roughness, against the 0.10 above. It crosses by 3.8x, and by 52x at
the smooth end of the roughness bracket.

**The emission scheme goes into the fork, and the source map stays outside it.**
Not either/or; the boundary goes where this project already puts it for every
other surface field.

`aeolian/` keeps computing the source map -- erodible fraction, clay, aeolian
roughness -- because those are pure functions of terrain, lithology, the lake
solution and the soil, exactly like the albedo, roughness and soil-water capacity
fields that already reach the model as boundary conditions. That is the half a
GCM structurally cannot know, and it is the half whose absence made the built-in
scheme useless: `fcoeff * land_mask` emits as much from forest as from salt pan.

The flux goes in the time loop, because it has to see the model's own winds. The
emission law is a threshold and then roughly a cube, and that nonlinearity is not
a detail: the same nonlinearity turned a factor of two in the wind-tail shape
into a factor of 8.8 in emission. A prescribed field cannot respond to the winds
that the dust itself changes, and that is exactly what the reopening test
objected to.

**Most of the machinery is already there**, which is what makes this a bounded
change rather than a component. ExoPlaSim carries 3-D aerosol transport in
`aerocore.f90` -- Lin and Rood flux-form semi-Lagrangian with a gravitational
settling term -- and radiative coupling in `radmod.f90` under `l_aerorad`, which
builds a per-layer optical depth in both shortwave bands from an optical-constants
file. `Model.configure` already exposes `aerosol`, `aerorad` and `aerofile`, so
none of that needs namelist surgery. And `aerocore.f90:621` already has the hook,
labelled for this:

    case(2) ! Case 2: dust
      mmr(:,:,NLEV,ic) = fcoeff*land

That one line is what gets replaced.

### What is genuinely missing, read out of the source rather than assumed

Checked 2026-08-17. The hook is one line; what sits around it is not.

**There is no longwave aerosol. At all.** `radmod.f90` puts the aerosol only in
its two SHORTWAVE bands -- `aod1`, `aod2`, `zaerr1`, `zaerr2` at line 1830
onward -- and the longwave solver contains no aerosol term of any kind. So this
model can cool with dust and cannot warm with it.

That is the same error this project already caught itself making once, and the
reason the dust question was relitigated at all: **the net cooling overestimate
was the problem.** Mineral dust absorbs and re-emits in the thermal infrared, and
for a coarse-mode-rich burden that offsets a substantial share of the shortwave
cooling. At a land-mean optical depth of 0.376, running shortwave-only is not a
small bias, and it is a bias in the direction this project is least able to
afford, because land albedo already sits inside a narrow flux window.

**So DUST-2 comes first, and it is a prerequisite rather than a follow-on.**
Until the longwave term is priced there is no way to say whether an in-model
shortwave-only dust is closer to the truth than no dust at all. If it is large,
the fork needs a longwave aerosol term as well as an emission scheme, and that is
a materially bigger change than replacing `case(2)`.

**One aerosol, one size, one density.** `plasimmod.f90:99` fixes
`parameter(NAERO = 1)` at compile time. The loops are already written as
`do jc=1,NAERO`, so raising it is structurally possible -- but `apart` and `rhop`
are scalars shared by every bin, `mmr2n` takes them as scalars to convert mass to
number density, and the optics below are single-valued too. Size-resolved dust
therefore means promoting all of those to arrays, not just changing a parameter.
The defaults are haze: 50 nm at 1000 kg/m3, against mineral dust at 0.1 to 20 um
and about 2650.

**The optics collapse to eight numbers.** `aerofile` is read by
`readdat(aerofile,1,8,aeroqs)`: one header line, then Qextinction, Qscattering,
Qbackscatter and g for band 1, then the same four for band 2. Optical depth is
then `nrho * PI * apart**2 * Qext * dz`. Everything spectral and size-resolved in
`analysis/dust_optics.json` has to be reduced to that, which is doable and has to
be done knowingly.

**And there is a sink that is documented but not physical.** The comment above
it says "put in a sink term at the bottom level to avoid infinite build-up of
haze particles", so it is a deliberate stability hack rather than an oversight --
which is worse, not better, because it means nothing physical was intended.
`aerocore.f90:869` applies
`mmr(:,:,nl,ic) = mmr(:,:,nl,ic)*10e-3` to the bottom level on every step, with
no comment: a 99% removal of the surface layer per timestep, presumably standing
in for dry deposition. Whatever emission is injected is scaled by it, so it has
to be understood before any flux is calibrated against it.

**The offline chain does not go away.** It remains the producer of the deposition
field, because loess (SURF-3) and the phosphorus return leg want size-resolved
deposition, which a single-mode in-model aerosol will not give. Two products, two
tools, and the reason stated rather than one quietly standing in for the other.

### The longwave is priced, and it changes the sign (DUST-2, 2026-08-17)

`analysis/dust_forcing.json`, from `exoplasim/scripts/dust_forcing.py`, at the
land-mean optical depth of 0.376:

| surface | shortwave | longwave | net |
| --- | ---: | ---: | ---: |
| ocean, a = 0.07 | -8.5 | +6.7 | -1.8 |
| vegetated land, a = 0.18 | -5.2 | +6.7 | +1.5 |
| playa fill, a = 0.40 | +0.4 | +6.7 | +7.1 |
| salt crust, a = 0.50 | +2.7 | +6.7 | +9.4 |

W/m2, for the fine end of the size bracket; the coarse end differs by under 16%
term by term and is in the file, though the net moves further, being a small
residual of two large ones. The longwave is a clear-sky window estimate and is an UPPER
bound, because part of the band is already opaque to water vapour and CO2.

Re-measured 2026-08-17 under PHYS-2. The two shortwave bands were being weighted
0.2566/0.7434, which was this script integrating `k25v.dat` a second time and
counting the bin width twice. The band split now comes from `lib/stellar.py`,
which reproduces what `radmod.f90:solarini` does with the hi-res spectrum, at
0.3844/0.6156. Band 1 is the more absorbing band, so under-weighting it made the
layer look more scattering than it is and every surface cooler than it is: the
shortwave term moved +0.25 to +0.5 W/m2 and the global mean net from +0.34 to
+0.55. Nothing about the sign or the threshold changed.

**The sign claim in this note's own pricing section is superseded.** Dust on this
world does not simply cool. It cools over ocean, is near neutral over vegetated
land and warms strongly over the bright closed-basin fill -- which is exactly the
surface that makes this world unusual, and exactly where the dust is.

The global mean lands at **+0.55 to +0.67 W/m2**, below the 1.5 W/m2 half of the
reopening threshold. That is not reassurance: it is ocean cooling cancelling land
warming across a net that spans 11 W/m2 by surface, -1.8 over ocean against +9.4
over salt crust, and a redistribution that large drives circulation whatever its
mean is.

**And it settles what DUST-3 must not do.** Switching ExoPlaSim's aerosol on as
shipped applies the shortwave alone, which is -3.3 to -3.8 W/m2 in the global
mean against a true +0.55 to +0.67. An error of 4.0 to 4.3 W/m2, against 21 W/m2
for the entire 0.85-to-0.95 stellar sweep that produced a 33 K range. In kelvin,
on the canonical conversion in `lib/sensitivity.py`, that is **3.3 to 3.6 K of
spurious cooling**: the second largest item in the error budget, behind only the
bare-rock-versus-vegetated question. **Shortwave-only in-model dust is far worse
than no dust**, so the fork needs a longwave aerosol term and not only an
emission scheme.

The sweep's own sensitivity would give 6.3 to 6.8 K for the same error, 1.9x
more, and it is the wrong one to use here: 33 K over 21 W/m2 is measured across
the ice transition and this world's baseline is not. Measured 2026-08-17, that
gap used to read as a factor of 3.4, and over half of it was the local
conversion being wrong rather than the regime being different. See BUDG-4.

One correction to the record while pricing it: the error budget said ExoPlaSim
could not switch aerosols on at all, because `radmod.f90` declares `aero_nl` and
never reads it. It is readable -- `aero_ini` is called at `plasim.f90:190` and
reads that namelist with `l_aerorad` and `aerofile` use-associated from `radmod`,
and `Model.configure` writes every field of it. The real limitation is narrower
and worse: the aerosol acts in the two shortwave bands only.

### Does dust change the carve verdict? The size is settled, the SIGN is not

Worth settling, because "dust is a radiation question, the carve is a water
question" is an easy assumption and it is wrong.

The global mean is negligible. +0.55 to +0.67 W/m2 is **+0.46 to +0.55 K** on the
canonical conversion, and it moves nothing.

**The size of the local effect is real.** Measured 2026-08-17 from the current
carve list, by `scripts/error_budget.py`, the number of the 1,721 overflowing
basins that stop overflowing:

| perturbation | +5% | +10% | +15% | +20% |
| --- | ---: | ---: | ---: | ---: |
| lake evaporation rises | 123 | 182 | 228 | 283 |
| catchment runoff falls | 82 | 143 | 175 | 219 |

A few percent on open-water evaporation moves a hundred basins. So dust is not
carve-neutral at any plausible magnitude, and the carve list must not be recorded
as dust-independent.

The two rows are not the same experiment and the earlier reading of this table
treated them as one. A runoff cut arrives through the criterion's denominator
alone; a lake evaporation rise arrives through its numerator, where it is worth
about 1.3x more per percent. The runoff row is also nearly independent of which
side of the water balance moved -- a 10% cut delivered by precipitation gives 143
basins and one delivered by land evaporation gives 139 -- which is what makes it
usable as a currency.

**The DIRECTION is open, and a first reading of it here was wrong.** That reading
took the top-of-atmosphere forcing -- +8.2 W/m2 over playa, +10.8 over salt crust
-- and treated it as surface warming, concluding that dust raises lake
evaporation and closes basins. It does not follow. The carve verdict is decided
by a Penman evaporation at an open water surface, which reads the SURFACE energy
balance, and for an absorbing layer that is a different quantity from the TOA
one:

| surface | TOA shortwave | surface shortwave |
| --- | ---: | ---: |
| ocean | -11.9 | -23.8 |
| vegetated land | -7.6 | -19.5 |
| playa fill | -0.1 | -12.1 |
| salt crust | +2.8 | -9.1 |

The layer absorbs about 12 W/m2 of shortwave whatever sits beneath it, so the
surface term is negative everywhere and the TOA sign over bright ground comes
from energy retained in the ATMOSPHERE, not delivered to the ground. The surface
also gains downward longwave from the dust, but that is damped by the atmosphere
already being opaque outside the window and is nothing like the TOA longwave
number.

So the two terms in Penman pull opposite ways: net radiation at the lake falls,
which suppresses evaporation, while the air above warms, which raises the vapour
pressure deficit and promotes it. For an open water surface the radiative term
usually dominates, which would mean **dust suppresses lake evaporation, keeps
basins overflowing, and makes a dust-free verdict UNDER-carve** -- the opposite
of the first reading here.

That is a preliminary indication and not a result -- and it is only half the
question, because it is only the lake.

**The catchment is the other half, and it is the larger one.** Runoff is the
denominator of the carve criterion, and it is the small residual of two large
numbers: measured 2026-08-17 on the current baseline, 126.5 mm/yr against a land
precipitation of 824.4 and an evaporation of 697.9, so 15.3% of P. Global
precipitation is constrained by atmospheric radiative cooling, and an absorbing
aerosol suppresses it through a fast adjustment that does not wait for a
temperature response. This dust absorbs about 11.9 W/m2 of shortwave where it
sits and about 7.5 W/m2 in the global mean, against a global latent heating of
77, which is a **precipitation reduction of order 10%**. What that does to runoff
depends on how much of the reduction evaporation takes with it, and the basin
counts are the criterion re-evaluated at each, from `scripts/error_budget.py`:

| E falls as fast as P | runoff | change | basins that close |
| ---: | ---: | ---: | ---: |
| 100% | 113.9 mm/yr | -10.0% | 186 |
| 80% | 99.9 | -21.0% | 261 |
| 50% | 79.0 | -37.6% | 339 |

**This limb closes basins whichever way the lake term goes**, and it is the one
an offline Penman perturbation cannot see, because it is a precipitation response
and not a surface energy balance.

How many basins that is worth was measured 2026-08-17 by
`aeolian/scripts/dust_runoff_sensitivity.py`, against 1,721 carved:

| dP | E follows | catchment runoff | basins carved |
| ---: | ---: | ---: | ---: |
| -5% | 0.80 | -8.3% | -121 |
| -10% | 1.00 | -10.0% | -51 |
| -10% | 0.80 | -16.2% | -201 |
| -10% | 0.50 | -24.1% | -309 |
| -15% | 0.80 | -23.9% | -288 |

**Do not translate that through runoff alone.** Scaling the criterion's
denominator by itself, with precipitation and evaporation left where they are,
moves 7 basins at -10% and 18 at -30%: twenty times less. The aridity index is a
ratio whose SIGN varies across the population -- the median basin here has P > E
-- so scaling the denominator makes a positive index larger and a negative one
more negative, and the two halves move toward opposite verdicts and cancel.
Reducing precipitation raises the numerator for every basin as well, so it does
not cancel. Both curves are in `aeolian/analysis/dust_runoff_sensitivity.json`.

So the answer is split in two, and DUST-10 has now settled the first half.

### DUST-10, measured: the lake term is small and it UNDER-carves

Computed 2026-08-17, and re-measured after HYD-13 corrected the grid mapping.
`analysis/dust_surface_forcing.nc` carries the per-cell dust perturbation to
`rss` and `rls`, and `carve_verdict.py --dust-forcing` adds them before Penman.
Surface forcing, not top-of-atmosphere, which is the whole point:

| surface | shortwave | longwave | net |
| --- | ---: | ---: | ---: |
| ocean, a = 0.07 | -12.3 | +4.5 | -7.8 |
| vegetated land, a = 0.18 | -10.8 | +4.5 | -6.3 |
| playa fill, a = 0.40 | -7.9 | +4.5 | -3.4 |
| salt crust, a = 0.50 | -6.6 | +4.5 | -2.1 |

Negative everywhere, including over the bright fill where the TOA term is
positive. Area-weighted the perturbation is -10.4 W/m2 shortwave and +3.8
longwave.

**Verdict with the dust surface forcing applied: 1,744 carved against 1,721
without. Twenty-three MORE basins carve.** Dust dims the lake, Penman evaporation
falls, more basins overflow. The direction is the opposite of the
top-of-atmosphere reading and the same as this note's corrected one; the earlier
claim that dust closes basins through the lake is dead.

It is small, and it survived every correction of the day almost unchanged: the
first measurement gave 24 basins on a mapping that turned out to shift every
basin half a planet, and the corrected mapping gives 23. The lake term does not
depend on the things that were wrong around it.

**So the lake term is settled, small, and in the recoverable direction.** A
dust-free verdict under-carves by about 23 basins on this limb, and an
under-carve can be corrected on a later pass while an over-carve cannot.

### What is still open is the larger limb

The catchment term is not in that number and is measured at 50 to 250 basins in
the OPPOSITE direction, over-carving, because suppressed precipitation cuts the
runoff the criterion divides by. It is two to ten times the lake term and it
decides the sign of the whole question. DUST-11 measures it and needs a
prescribed-dust climate run; nothing offline can reach it.

**That run is now specified and its parts are built.**
`aeolian/notes/prescribed-dust-run.md` carries the specification, the gates and a
prediction made before the run; `exoplasim/patches/exoplasim-3.4.2-prescribed-dust.patch`
is the model change, written and deliberately not applied; and surface code 1811
is the boundary field. The route is two changes and not six: the radiation sees a
prescribed field, and it gets a longwave term. No emission scheme, no transport,
no wet scavenging, no size bins, and `aerocore` is not on the path -- which also
sidesteps its bottom-level sink, because nothing is being advected or
accumulated.

**The Generic PCM is still out, but not for the reason given above.** That reason
was that the trade is not worth making "for a term whose sign is settled" -- and
the term is no longer small, so it no longer applies. The reason it is out now is
compute: it costs two orders of magnitude in runtime, and this project runs its
models iteratively on one workstation. A model that cannot go round the loop is
not an option however correct it is.
