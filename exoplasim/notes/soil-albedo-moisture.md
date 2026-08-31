# The modelled soil albedo and the water in the ground under it

*A worldbuilding project. Everything below is about the simulated planet Vesper
and the models that produce it: the ExoPlaSim climate column, the land water
column under it, and the albedo boundary condition this project stages into
them. Measured on 2026-08-27 against `canonical-10m-base` at T21.*

Wet ground is darker than dry ground, and until now the modelled land surface
kept one albedo through every wetting and drying cycle the land column
simulated. This note is what the term is made of: where the water state comes
from, what closes the band the model reads and no measurement covers, what the
term is worth, and the one class it refuses.

## What the modelled soil albedo does through a cycle

Two staged endmembers and a curve between them.

    codes 174 / 175 / 176     the DRY endmember, broadband and the band pair
    codes 1742 / 1750 / 1760  the SATURATED endmember, the same three surfaces

`landmod`'s `getalb` under `nwetsoil = 1` calls `wetalb`, which reads the
degree of saturation of the land column's SURFACE LAYER and mixes the two ends
through `landcolumn`'s `wet_soil_albedo`. Over a cycle:

- **Dry season.** The surface layer empties. An empty layer is a dry soil, not a
  soil at its wilting point, so the modelled albedo is exactly the staged dry
  field and not a fixed fraction of the way toward the wet one. That equality is
  bitwise: the mixing branches on the ends rather than reaching them by
  arithmetic.
- **First rain.** The surface layer fills before anything below it does, because
  the column is a cascade and the surface flux enters the top. A 0.02 m store
  fills on the first substantial rain of the season rather than over it.
- **Wet season.** The layer sits near its capacity and the modelled albedo sits
  near the saturated end. It never reaches that end: the layer's capacity stops
  at field capacity, so the saturated endmember is an endpoint the mixing
  approaches and does not touch.
- **Drydown.** The layer empties under evaporation in about a day, and the
  albedo returns to its dry value on that timescale rather than on the column's.
  This is the change that matters most and the one the previous pass identified:
  the modelled soil albedo is now right in PHASE. Craft and Horel measure a
  Bonneville salt crust recovering from 0.22 to 0.32 in eight days after
  flooding, and a 0.5 m store cannot do that.
- **Frozen ground.** The layer's store is LIQUID water; ice occupies pore space
  and comes off the capacity in the cascade, so a frozen surface reads as a dry
  one with no branch for it.

The curve between the ends is **not** linear in the albedo. Sadeghi, Jones and
Philpot derive the reflectance of a wetting soil from Kubelka-Munk two-flux
theory, and what interpolates linearly in the water content is the transformed
reflectance `r = (1 - R)^2 / (2R)`, not `R`:

    r(S) = (sigma*r_dry*(1 - S) + r_sat*S) / (sigma*(1 - S) + S)
    R(r) = 1 + r - sqrt(r*r + 2*r)

`alpha_dry*(1 - S) + alpha_wet*S`, the CLM form this row started from, has the
same two ends and the wrong shape between them, and the disagreement is largest
at low saturation, which is where a drying skin spends its time.

`sigma` is the ratio of the dry soil's scattering coefficient to the saturated
soil's. It is ONE where the water's own scattering is negligible against the
soil's, which Sadeghi verify at the strongly water-absorbing short-wave infrared
wavelengths that make up band 2. One is adopted in both bands as that limit
rather than as a fit; their visible-band fits run 0.042 to 0.528, and the
direction is known -- a value below one darkens the modelled surface sooner in a
wetting cycle without moving either end. That is the declared bracket on band 1.

**The bracket is swept OFFLINE and it does not buy a run.** Remixing the staged
pair at `wetsigma1` of 1.0, 0.528 and 0.042 with band 2 held at the limit
Sadeghi verify, land-mean fall against the skin's fill fraction:

| f | sigma1 = 1.0 | 0.528 | 0.042 | worst, in kelvin at the measured attenuation |
| --- | --- | --- | --- | --- |
| 0.05 | 0.001648 | 0.002024 | 0.005927 | +0.16 |
| 0.15 | 0.004747 | 0.005610 | 0.010028 | +0.20 |
| 0.50 | 0.014064 | 0.015305 | 0.017873 | +0.15 |
| 1.00 | 0.024668 | 0.025255 | 0.025984 | +0.05 |

The bracket is worth at most 0.20 K and it is largest where the term itself is
smallest: at a fifth of capacity it is a factor of two on a fall of half a
hundredth. It never approaches the paired-arm resolution bar, so sweeping it in
the model would buy a pair of runs that could not tell the arms apart. What it
does say is that the band-1 uncertainty is comparable to the whole term at low
wetness and negligible at high, which is the opposite of how a bracket on an
endpoint would behave and is the signature of a bracket on the SHAPE.

The one FIELD measurement of the shape this project holds sits at the adopted
arm. Idso et al. (1975) ran four pyranometer experiments over one irrigated loam
field with concurrent gravimetric sampling at eight depths, and found the albedo
a linear function of the water content of the top 0.2 cm over 0.00 to 0.18
volumetric, which on that soil is roughly the lower 0.4 of the saturation range.
Sadeghi's `sigma = 1` is the LEAST concave member of its family -- the initial
slope goes as `(r_sat - r_dry)/sigma` -- so a value below one would be more
concave than Idso measures, not less. The measurement is broadband and the
bracket is on band 1, so this does not close it; it says which end of it the
evidence points at.

## The surface layer, and how its water is conserved

The state the albedo reads is the top water layer of the layered land column,
0.02 m thick, declared in `pedology/config/land_column_properties.yaml` under
`surface_layer`. That block is the ONE top-of-column profile and both consumers
read it at their own depth: this row takes a degree of saturation over the
albedo depth, and DUST-17 takes a gravimetric water content over its emitting
depth. A near-surface saturation defined in the albedo builder would be the
second central hydrology DUST-17 exists to prevent.

0.02 m is an UPPER bound on the depth a radiation scheme sees, and the direction
is stated because it is the direction that understates how fast the modelled
surface dries. Idso et al. measure that directly: their Figure 5 puts a
single-valued albedo-against-water relation on the 0-0.2 cm layer and a near-step
function on 0-10 cm, and their own conclusion names a controlling layer "less
than 0.2 cm thick". ClimaLand's `albedo_calc_top_thickness` is 0.02 m and this
project takes the coarser of the two, because the same layer is what DUST-17's
emitting depth reads and a 0.002 m store in a tipping-bucket cascade empties
inside one model timestep.

It is TWO changes to the column and only one of them moves water.

**The split is free, and exactly.** ExoPlaSim's liquid water is a tipping-bucket
cascade: the surface flux enters the top layer, each layer fills to capacity and
passes its excess down, a withdrawal larger than the top layer holds draws on
the layers beneath, and what the base cannot take backs up and leaves as runoff.
Splitting the top layer while holding the column's total capacity fixed leaves
the total that cascade carries unchanged, because what enters is the surface
flux and what leaves is base overflow and neither depends on where the internal
boundaries are. `land_column_properties.py`'s `check_split_invariance` DRIVES
that rather than asserting it: a two-layer column and a three-layer column whose
first two capacities sum to the two-layer column's first, one flux sequence
through the same cascade, and the total and the runoff must agree. Over 4096
steps they agree to 1.1e-13 mm, which is the rounding of the two sums and not a
transfer of water; the tolerance was set at 1e-9 mm before the first run.

The 0.5 m boundary survives the split -- 0.02 + 0.48 + 1.0 -- which matters
because that is the only cut at which the climate column and the ecology column
share a boundary, and LSHY-5 needs them to share it.

**The floor is a real change and it is declared.** `available` is field capacity
minus wilting point because that is what a ROOT can remove, and evaporation from
a bare surface is not a root: the top of a drying soil goes to air dry, which is
below the wilting point and outside every state the column's capacity carried.
So the surface layer's capacity is cut from air dry instead, and the column's
capacity gains that layer's sub-wilting water. Measured on this build:

| quantity | median | p10 to p90 |
| --- | --- | --- |
| surface layer capacity | 6.29 mm | 3.27 to 7.34 |
| sub-wilting increment | 2.97 mm | 1.62 to 4.18 |
| `awc_mm`, what LPJ-GUESS reads | 48.7 mm | 22.8 to 215.6 |
| `evaporable_mm`, what ExoPlaSim installs as `dwmax` | 51.6 mm | 24.4 to 219.6 |
| increment over `awc_mm` | 0.0581 | 0.0181 to 0.0949 |

The two capacity columns are two quantities and not one number written twice: a
root cannot reach the water between air dry and the wilting point and a bare
drying surface can. `surface_layer_report` checks the identity that defines them
-- `evaporable_mm - awc_mm` equals the surface layer's increment, cell by cell,
to 2.3e-14 mm -- so the column gained exactly the water the surface layer can
now reach and not a millimetre more. The layers beneath keep their
plant-available semantics unchanged, which is why `soilsrwp` and `soilsrfc`
still describe them and the soil heat solver still reads a layer it understands.

**Air dry is declared as zero**, and the argument is that it is the state every
dry reflectance spectrum this project stages was MEASURED at, so an empty
surface layer and the staged dry albedo are the same state and the mixing has no
level shift at its dry end. That was the objection that closed this row before,
and it is answered by the column rather than by a coefficient. The cost is
bounded and one-signed: a real soil holds adsorbed water at the ambient
humidity, so the layer can give up at most `0.02 * theta_ad` more water than it
should, which is under the 2.97 mm increment. The alternative -- evaluating air
dry through this contract's own retention closure at the Kelvin potential of a
declared humidity, near -160 MPa -- is registered and not adopted: it is a
hundredfold extrapolation past the wilting point, far outside the suctions Cosby
fitted, and Clapp-Hornberger carries no adsorption branch at all. Adopting it
needs a dry-end retention form this project does not hold.

## What band 2 rests on

Band 2 carries 0.617628 of this star's flux against band 1's 0.382372, and it is
the band liquid water absorbs in. It was the open half of this row: the held
wetting measurements were a luminous reflectance table over 0.38 to 0.77 um and
two broadband in-situ pyranometer pairs, so band 1 was bracketed and band 2 was
not, and an Earth calibration cannot be carried across.

It is closed WITHOUT a band-2 measurement, and the reason is that both mechanical
theories of why wet things are darker give the wet reflectance as a closed form
in the DRY reflectance AT THE SAME WAVELENGTH:

**Lekner and Dorf (1988).** A water film over a rough surface sends diffusely
reflected light back down onto it by total internal reflection at the film's
upper face, so the surface gets more chances to absorb. Closed form and no
fitted constant: the internal-reflection probability follows from the film's
refractive index, and the change in the surface's own absorptance follows from
the drop in relative index at the substrate.

**Twomey, Bohren and Mergenthaler (1986).** In a finely divided medium,
replacing interstitial air with water lowers the particle-to-medium index ratio
and makes single scattering far more forward-peaked, so a photon needs more
scatterings to escape and is more likely to be absorbed on the way. Similarity
scaling turns that into a reduced effective single-scattering albedo and
Chandrasekhar's semi-infinite isotropic result turns that back into a
reflectance.

Applied per wavelength to the same reflectance spectra `rock_albedo_bands.py`
already integrates, and then band-integrated against this star, they produce a
band-2 wet endmember from band-2 dry data. The two bands come out wetting by
different amounts because the dry spectra differ between the bands, not because
a band-2 ratio was asserted: `granite` wets to 0.366 of dry in band 1 and 0.403
in band 2, `playa_clastic` to 0.324 and 0.373.

**The level has to be right before the map is applied**, and that is the trap
this derivation exists around. Both maps are strongly nonlinear in the dry
reflectance and both fix the endpoints, so the wet/dry RATIO is not a property
of the material alone -- it depends on where the material's dry reflectance
sits. A laboratory powder and the outcrop it came from are 2.2 to 5.1 apart on
the same rock, so each spectrum is rescaled to the class's own broadband albedo
before wetting: the shape from the library, the level from the rock table.

**The two are two arms of a bracket and never an average.** Lekner and Dorf say
plainly which medium each is for: theirs "would seem to apply best to rough
solid surfaces", the other's "to finely divided media, such as sand". This
world's land is both. The interstitial arm is what `--wetting-arm` stages by
default, because this land is regolith and playa fill far more than bare
outcrop, and the film arm is the other end.

**The bracket is a prediction that can fail, and it holds.** It must contain
every wet/dry pair this project holds:

| source | dry | wet | bracket | in |
| --- | --- | --- | --- | --- |
| Penndorf (1956) clay soil | 0.150 | 0.075 | 0.045 to 0.098 | yes |
| Penndorf (1956) sand | 0.310 | 0.180 | 0.121 to 0.197 | yes |
| Penndorf (1956) bare rich soil | 0.072 | 0.055 | 0.019 to 0.056 | yes |
| Angstrom (1925) sand, via Lekner Table III | 0.182 | 0.091 | 0.058 to 0.117 | yes |
| Angstrom (1925) black mold, via Lekner Table III | 0.141 | 0.084 | 0.042 to 0.093 | yes |
| Craft, Horel (2019) Bonneville halite crust | 0.450 | 0.220 | 0.219 to 0.301 | yes |
| Malek et al. (1990) Pilot Valley playa | above 0.75 | 0.240 | 0.558 to 0.604 | NO |

The exception is named in advance rather than tolerated. Malek's surface is a
thin halite crust over SHALLOW BRINE and its wet value is the crust with brine at
the surface, which is inundation and not a wetted skin; its dry value is stated
only as a bound. Both mechanisms describe a wetted medium and neither describes
ponded water, and the modelled surface layer caps at field capacity, so the
model cannot reach that state either.

Both implementations also reproduce their own papers' printed numbers before
they are used on anything: Lekner and Dorf's average interface reflectance
(0.0665 against the 0.0667 their own internal-reflection probability implies),
that probability (0.4749 against 0.475), Angstrom's lower estimate of it (0.4375
against 0.437) and the wet-to-dry absorptance ratio at three substrate indices
(1.068, 1.082, 1.100 against 1.07, 1.08, 1.10); and TBM's worked example carried
end to end (0.30 to 0.115 against 0.12, the residual being that the paper runs
its example through the zenith reflectance where a model albedo is a flux
albedo).

**What band 2 still does NOT rest on**, and it is one-signed. Both mechanisms
treat the liquid as non-absorbing and describe only what it does to the geometry
of scattering. Liquid water absorbs strongly in parts of band 2, so both arms are
UPPER bounds there and the bracket is open at its dark end. Turning that into a
correction needs an absolute scattering coefficient for the dry soil, which
nothing in this tree carries, or a wet-and-dry reflectance spectrum pair measured
on one sample. `references/INDEX.md` names the one that would do it.

## What the paired arm measured

Two twenty-five-orbit runs seeded from `run_893e276ee029/MOST_REST.00209`, a
210-orbit equilibrated state on this build. One binary, md5-identical in both run
directories; one seed, md5-identical; `landmod_namelist` differing in one line,
`NWETSOIL = 1` in `run_a1c35075747c` against `0` in `run_598eb57c5a34`. The wet
arm's `landini` printed the PHYS-15 banner at 0.0 to 0.7642 with both shape
sigmas at one, so it read the staged pair rather than stopping on the sentinel;
the dry arm printed nothing. Both segments are `diagnostic` at one I/O regime,
and both pass every diagnostic settling criterion. Measured over orbits 15 to 24.

| quantity | value |
| --- | --- |
| land-mean `alb`, `NWETSOIL = 1` | 0.230051 |
| land-mean `alb`, `NWETSOIL = 0` | 0.236448 |
| **the fall** | **0.006397** |
| paired scatter across the ten orbits | 0.001226 |
| the bar, two root two times the larger arm's own standard error | 0.001115 |
| resolved | yes, at 5.7 times the bar |

**The registered bracket was 0.003780 to 0.006745 and the measurement is inside
it**, close to the upper end. Both ends were computed from the donor's own
restarts before either arm ran.

**The background half of the prediction is exact.** The wet arm's own restarts
give a background fall of 0.006709, against 0.006745 predicted from the donor's:
agreement to half a percent, on a quantity derived from `dwatcl` and the staged
pair through Sadeghi Eq (13). So the mixing does in the running model exactly
what the offline derivation says.

**THE MASK COSTS ALMOST NOTHING, AND THE REASON IS PHYSICAL.** The lower end of
the bracket assumed that where snow or glacier ice hides the background the term
is lost entirely, which would have cost 44 per cent of it. The measured cost is
the difference between the background fall and the realised one, 0.006709 against
0.006397, which is 4.6 per cent. **The masking and the wetting are
anti-correlated**: a snow-covered cell has a frozen skin, ice comes off the
layer's capacity in the cascade, and the model already reads it as dry, so it was
contributing nothing to the fall for the mask to take away. Treating the mask as
a switch on the mean is what made the lower end far too low; the fault is in that
reasoning and not in the model.

**`f` barely moves between the arms**, so the term is very nearly linear in it:
0.3313 in the wet arm against 0.3257 in the dry, both with a per-orbit spread near
0.017. The albedo feedback on the skin's own wetness is inside the noise at
twenty-five orbits.

**In the units the budget consumes**, at the land fraction the build manifest
gives by surface class and the attenuation `scripts/error_budget.py` MEASURES
for itself, 0.38 with a measured span of 0.29 to 0.48:

| | value |
| --- | --- |
| top-of-atmosphere forcing at the measured attenuation | 0.358 W/m2 |
| forcing across its bracket, 0.29 to 0.48 | 0.273 to 0.452 W/m2 |
| kelvin across the same bracket | +0.19 to +0.31, centrally +0.24 |

For scale, the dust radiative item this project prices is 0.34 to 0.61 W/m2, so
the modelled soil's wetting is the same size as its dust.

**The claims that could have failed, tested where they are made.** Monotonicity
and the `evaporite` refusal are claims about the BACKGROUND the mixing produces,
not about `alb`, which carries snow that diverges chaotically between the arms.
On the wet arm's own background over the window: the worst brightening of any
land cell against the staged dry field is 5.3e-16, and the worst move on the
twenty-three cells staged wet equal to dry is 4.2e-16. Both are zero in double
precision, so no cell brightens and the refusal survives into the running model.

## The temperature separated, and the prediction said it would not

`tas` came out +0.1972 K warmer in the wet arm, against a bar of 0.0874 K taken
from these arms' own standard errors. It is resolved at 2.3 times the bar, and
`ts` agrees at +0.2013 K. **The registered prediction said the pair could not
separate in temperature, and that was wrong.**

The error is worth more than the result. The bar was imported from the arms this
project had already run, whose paired `tas` standard error over twenty-five
orbits is 0.56 to 0.75 K, giving 1.6 to 2.1 K. These arms' standard errors are
0.017 and 0.031 K, twenty times smaller, and the run length is the same. **What
differs is the donor.** The earlier arms were seeded from a thirty-seven-orbit
control still relaxing, so each arm carried a trend and the paired difference
inherited its scatter; these are seeded from a 210-orbit equilibrated state, so
each arm is stationary and the difference is nearly pure signal. The instrument's
power is set by how settled the DONOR is at least as much as by how long the arms
run, and quoting a bar measured on one regime against an experiment in another
understates the instrument by a factor of twenty.

**And the separation measures something the budget had only bracketed.** A
land-mean `alb` fall of 0.006397 reaches +0.1958 K, on the two arms' fitted
asymptotes, only at an attenuation of 0.308. That was the first point this
project had of a multiplier `scripts/error_budget.py` was declaring at 0.5 while
saying the honest claim was a factor of two.

**It is now one of four, and the budget's constant is measured rather than
declared.** `notes/audits/albedo-attenuation.md` adds three paired endmember
spin-ups already on disk and puts the constant at 0.38 with a measured span of
0.29 to 0.48. This term sits at the bottom of that span and the biosphere
endmember term at the top, which is what the attenuation being a property of the
modelled atmosphere rather than of the surface allows: this term acts where the
modelled skin is wet, which is where it rains, and the endmember term is spread
over land that includes the dry closed-basin interiors.

**And this term's own figure moves with which albedo it is denominated in.** The
0.006397 above is the REALISED fall, which is the staged mixing as the modelled
snow leaves it; the staged prediction it was registered against is 0.006745, and
back-solving on that gives 0.292. Every item the budget prices is a staged delta,
so 0.292 is this term's row there. The two differ by 5 per cent here because the
snow mask costs this term almost nothing; on the endmember term the same choice
is worth a factor of 1.45.

## The sign is disputed on salt crust, and the class is refused

The three salt-crust sources do not agree, and the disagreement is not noise.
Malek and Craft measure the surface at the time of wetting and both find it much
DARKER: above 0.75 to 0.24, and 0.45 to 0.22. A twenty-year MODIS series over
Uyuni composites annual means and finds wet or cool years BRIGHTER, 0.65 against
0.55 to 0.60.

They are not measuring the same thing. A water film darkens a crust; a wet year
also dissolves and reprecipitates the crust and keeps dust off it, which
brightens it, and a twenty-year composite cannot separate the two.

**Sadeghi's Fresnel term rules itself out as the second effect, which is new.**
It is the one physical route by which wetting brightens a surface in this
framework: a nearly continuous surface film adds a specular reflection on top of
the volume reflectance, which their Eq (19) writes as the normal-incidence
reflectance of water times the volumetric water content. That is at most 0.0204
times the porosity, under 0.01, against an interannual difference near 0.075. It
is short by close to an order of magnitude, so the Uyuni brightening is a
crust-CONDITION effect and not an optical one, and this model carries no state
for the condition of a crust.

So `evaporite` is REFUSED rather than bracketed, and the refusal is in the staged
field: the class is written with its dry pair as its wet pair, so the mixing
returns the dry albedo at every saturation. Two checks hold it there --
the wetting ratio of every refused region must be exactly 1, and a grid cell made
entirely of refused material must carry identical dry and saturated fields; at
T21 no cell is pure enough for the second, and the first covers every region.
`water` is refused on the different ground that open water is not a surface that
wets.

`evaporite` is 2.565 percent of this build's land and the brightest surface on
the simulated planet. `playa_clastic` is 25.795 percent and the largest single
class, and it is NOT refused: the dispute is about halite crust, and playa mud
and fan fill are a soil.

## The magnitude

`build_surface_albedo.py` recomputes all of these per build into
`albedo_report.json`. All are land-mean albedo deltas, the unit
`scripts/error_budget.py` consumes, except where a per-cell figure is given
beside them.

**What the staged pair is worth.** The land-mean albedo of the two fields now
written, measured on 2026-08-27 against the VEGETATED staged field with lakes
composited in, which is the field a run reads:

| quantity | value |
| --- | --- |
| land mean, dry | 0.170772 |
| land mean, saturated endmember, band 1 | 0.071698 |
| land mean, saturated endmember, band 2 | 0.183283 |
| land mean, saturated endmember, broadband | 0.140616 |
| the swing between the two staged ends | 0.030156 |

**The endmember swing is not the term, and on this field the difference is most
of it.** The staged saturated pair is the endmember at a degree of saturation of
one; the modelled surface layer caps at field capacity, which maps to 0.7642, so
the mixing approaches that end and never touches it. Evaluating the mixing there,
per cell and then land-meaned because it is concave and the two orders do not
commute:

| quantity | value |
| --- | --- |
| land mean at a full surface layer, broadband | 0.146138 |
| the fall, broadband | 0.024633 |
| the fall, band 1 | 0.017506 |
| the fall, band 2 | 0.029103 |

`albedo_report.json`'s `wetting.at_full_surface_layer` carries all four and
recomputes them per build, so nothing downstream re-derives them from the two
endmembers.

**The recombination identity holds at both ENDS and not in between, and the size
of that is worth knowing before it is found.** `build_surface_albedo.py` checks
`z1*175 + z2*176 == 174` on both staged pairs before writing, but `wetalb` mixes
all three fields independently and the mixing is nonlinear, so the broadband
field and the recombined band pair separate as soon as the skin is wet. At a full
surface layer they differ by 3.5e-5 in the land mean, 0.024633 against 0.024668.
Under `NSIMPLEALBEDO = 0` the radiation uses the PAIR and `dalbclim` is the
diagnostic, so an arm reading the modelled `alb` is reading the smaller of the
two. The gap is 0.14 percent of the term and is reported rather than corrected:
forcing the identity through the mixing would mean mixing two fields and deriving
the third, which makes the broadband diagnostic a different quantity from the
staged broadband endmember at every saturation but zero.

**Which surface the term acts on is what sets its size, and the field moved.**
The wetting maps act on SOIL. A canopy is not a wetting surface and open water is
not one either, so `build_surface_albedo.py` gives the covered fraction of a cell
and its lake fraction the same value at both ends and only the bare fraction
wets. On a bare-rock staged field the term was worth a swing of 0.155841; on the
vegetated field with lakes it is worth 0.030156, a factor of 5.2 smaller, and the
whole of that factor is cover. For scale, the bare-against-vegetated gap this
project brackets is 0.069381, and `build_surface_albedo.py` argues that gap is 15
to 19 W/m2 in absorbed flux against 21 W/m2 for the entire 0.85-to-0.95 stellar
sweep that produced a 33 K range. The wetting fall at a full surface layer is
about a third of it.

## How wet the modelled skin actually is

The ceiling above is the swing times how wet the skin is, and that fraction is
readable. `landmod.f90:1553` writes `dwatcl` into the restart with `mpputgp`, and
`dwmax` and `dsoilwfc` are carried beside it, so `dwatcl(:,1)` over
`dwmax * dsoilwfc(:,1)` is `wetalb`'s own arithmetic recovered from the file
rather than a reconstruction of it. Every orbit of a run leaves one restart, so a
run is also a sample of `f`. Measured over the settled block of a 210-orbit
`nwetsoil = 1` run on this build, restarts 180 to 209:

| quantity | value |
| --- | --- |
| `f`, land-area mean | 0.3283, sd 0.0146 across the thirty |
| `f`, land median | 0.0096 |
| `f`, land p90 | 0.9997 |
| share of land above `f` = 0.5 | 0.319 |

**`f` IS BIMODAL, and that is the finding rather than a detail.** About a third
of the modelled land carries a skin at its capacity and most of the rest is at
air dry; very little sits between. That is what a 0.02 m store in a
tipping-bucket cascade does -- it fills on the first rain and empties in about a
day -- and it is the shape a mean hides completely.

**So the mixing is evaluated per cell and then meaned, and the two orders differ
most here.** The fall at the mean `f` of 0.3283 is 0.0098; the mean of the
per-cell falls is 0.006745, a factor of 1.45 smaller. The concavity penalty is at
its largest against a distribution with its mass at the two ends, which is
exactly this one. Every figure below is the second.

**And the term reaches the surface over three quarters of the land.** Snow and
glacier ice override the background pair after `getalb` has mixed it. On 0.7673
of the land, area weighted, the model's own `dalb` equals the mixed background to
2e-4, so the term arrives undiminished there and is blended away elsewhere. That
puts the realised background fall between 0.003780 and 0.006745: the lower figure
treats the mask as a switch and the upper ignores it, and the mask is a blend, so
the truth is between them.

The same two paths agree on the level as well as the difference. The land-mean
`alb` the output stream writes over that block is 0.229733 and the land-mean
`dalb` in the restarts over the same orbits is 0.229337 -- a snapshot at one
orbital phase against a mean over the orbit, differing by less than either one's
inter-orbit scatter.

**In kelvin, and against the instrument that would have to see it.** Through
`scripts/error_budget.py:albedo_to_kelvin` at the land fraction the build
manifest gives by surface class, 0.432841, and the attenuation that file measures
for itself, across its measured span:

| f, the skin's fill fraction | fall | K at 0.29 | K at 0.38 | K at 0.48 |
| --- | --- | --- | --- | --- |
| 0.05 | 0.001631 | +0.05 | +0.06 | +0.08 |
| 0.15 | 0.004705 | +0.14 | +0.18 | +0.23 |
| 0.50 | 0.013997 | +0.41 | +0.54 | +0.68 |
| 1.00 | 0.024633 | +0.72 | +0.94 | +1.19 |

The last row is a permanently saturated skin everywhere at once, which the
cascade cannot hold, and it is +0.94 K at the measured central attenuation. **The
paired-arm instrument cannot see that in temperature.** A 25-orbit paired
difference in global-mean `tas` carries a standard error of 0.56 to 0.75 K on the
held arms, so the bar the arms are reported under, two root two times the larger
standard error, is 1.6 to 2.1 K. The ceiling sits below its own bar, and the
measured attenuation puts it further below: at 0.38 the term does not reach 1.0 K
at any fill fraction, since a permanently saturated skin everywhere at once is
worth +0.94 K and the cascade cannot hold that state anyway. The measurement that DOES resolve it is the land-mean `alb` difference,
which both arms write and which is set by the boundary condition rather than by
the circulation; a temperature separation is a consequence to look for and not
the test. `docs/src/practice/failure-modes.md` class 34 is the rule this is an
instance of, and the direction is the useful one: the term is real, one-signed
and small, and an arm that reported "not resolved" in `tas` would be reporting
its own scatter.

**The hard ceiling, both bands, all classes.** Retained from the declaration and
still checked before anything is written: a fully INUNDATED modelled surface
cannot be darker than the open water that would cover it, so the substrate's own
albedo minus open water's bounds that branch per region.

| quantity | value |
| --- | --- |
| land-mean substrate albedo | 0.259739 |
| open water albedo, from the export's rock table | 0.060000 |
| ceiling, land-mean | 0.199739 |
| ceiling, darkest land region (basalt) | 0.040000 |
| ceiling, brightest land region (evaporite) | 0.440000 |

The guard on it stands: the script refuses to write a field in which any land
region sits at or below open water's albedo, because that is the case in which
the ceiling's sign claim would be false. It bounds INUNDATION, which is a
different quantity from the saturated endmember staged in 1742, 1750 and 1760 --
the surface layer caps at field capacity and never reaches either.

**Band 1 only, from Penndorf's ratio bracket**, and **both bands on salt crust,
from the two in-situ pairs**, are both still computed per build. They are now
independent CHECKS on the derived endmembers rather than the only magnitudes
available: Penndorf's three pairs and Craft's are four of the seven the bracket
above has to contain, and the salt-crust bracket describes a class the staged
field refuses.

## What turns it on

`nwetsoil` defaults to 0, and three things have to be true before the model will
run it. The switch checks all three at initialisation rather than assuming them,
and `run_exoplasim.py` checks what it can from the configuration first, which
costs a config read instead of a launched model.

1. **`nlandwcol = 1` with the three-layer surface cut.**
   `config/planet.yaml`'s `surface.land_water_column` declares `layers: 3` at
   `layer_thickness_m: [0.02, 0.48, 1.0]`, which is the geometry `surface_layer`
   declares, and `DSOILWF` is derived from those thicknesses at the point the
   namelist is written because it is the depth share a uniform profile gives. The scalar bucket is refused
   outright: it has no surface layer, and the column it does have dries on the
   wrong timescale for an albedo. `landini` refuses the bucket; what it cannot
   check is that the first water layer IS the declared albedo depth, because
   that depth lives in the contract and not in the model, so `run_exoplasim.py`
   refuses a three-layer column cut anywhere else.
2. **The saturated pair staged.** Codes 1742, 1750 and 1760 are written by
   `build_surface_albedo.py` beside 174, 175 and 176, and the model refuses on
   the negative sentinel rather than mixing toward it. They are read on every
   start rather than carried through the restart, because they are a boundary
   condition and not a state. They join `intended_surface_codes` only under the
   switch, so at `nwetsoil = 0` they are not staged and nothing opens them. At
   `nwetsoil = 2` the third point, codes 1743, 1751 and 1761, joins them on
   exactly the same terms.
3. **`dwmax` installed from `evaporable_mm`.** The surface layer's capacity is
   cut from air dry, so the column's capacity is `evaporable_mm` and not
   `awc_mm`. `build_surface_soil_water.py` installs it, and it selects the
   column by the declared GEOMETRY rather than by a switch of its own: the extra
   water is water the top 0.02 m can give up, so a column without that layer
   gets `awc_mm` and every layer cut from the wilting point. `awc_mm` is what
   LPJ-GUESS reads and it does not move.

## The third staged point, and why a pair cannot carry a cell

The mixing is linear in the Kubelka-Munk transform and therefore concave in the
albedo, so mixing a cell's MEAN dry albedo toward its MEAN saturated one is not
the mean of the mixings of the rocks the cell holds. The staged pair reproduces
the cell's own mixture at the two ends of the saturation axis and nowhere
between them, and no other PAIR repairs it: code 174 IS the dry albedo, the
radiation reads the same field on dry ground, so it is pinned to the cell's
area-mean dry albedo by a boundary condition that has nothing to do with
wetting, and 1742 is pinned at the other end by the same argument.

`nwetsoil = 2` mixes through a THIRD staged field. Codes 1743, 1751 and 1761
carry the cell's own area mean of the per-region mixings at one declared
saturation, and `wetalb` runs the same Sadeghi curve in two segments through it,
so the composition is exact at three saturations instead of two. Every one of
the three staged fields is an exact area mean of a per-region quantity, so the
agreement at three points is a construction rather than a result. A per-cell
Sadeghi `sigma` would reach a similar place by FITTING and is refused on that
ground: it would be the residual of a fit, with no derivation to carry to
another planet.

**The declared saturation is derived from two files and chosen in neither.**
`landmod.f90`'s `drhsfull` is the fill fraction of the surface layer above which
the evaporation limiter's wetness factor reaches one, and `wetalb` maps a fill
fraction onto saturation through `skinsrad` and `skinsrfc`. The knee is the
image of one under the other. `build_surface_albedo.py` reads `drhsfull` out of
the model source and the endpoints out of `config/planet.yaml`; `wetalb`
recomputes the same expression from its own two constants, so the saturation the
field is staged at and the one it is mixed at are one fact.

**What it is worth.** `analysis/spatial_reduction_gap.py` prices the two-field
mixing at a peak of 0.147 W m-2 of global-mean absorbed shortwave against the
accepted baseline's 0.12 W m-2 state-storage tolerance, past it over a third of
the saturation range, and the three-point form at 0.042 and inside the bar
everywhere. The criterion was fixed before either was run.
`notes/audits/nonlinear-spatial-reductions.md` section 7 carries the measurement
and the builder repair it refused first.

**The order is the whole of the staging.** The mixing is applied to each REGION's
own pair and the result reduced afterwards. Reducing first and mixing after
reproduces the defect the field exists to remove, and the builder's bracket check
is what catches it: the staged knee must lie between the dry and the saturated
field cell by cell, because the mixing darkens monotonically.

`--mode scaled` refuses the third point rather than scaling it, because that mode
multiplies already-reduced endmembers by a per-cell factor and the mixing of two
scaled endmembers is not the scaling of their mixing. `landini` then stops on the
sentinel. In `--mode modelled` the third point is composited exactly as the two
ends are: a canopy's wetting ratio is one, so its mixed albedo is its dry one at
every saturation, and the blend is affine in its endmembers.

The switch itself and its five companions reach `landmod_nl` from
`surface.soil_albedo_moisture`, written unconditionally by `configure_otherargs`
and checked by `expected_namelist_keys`, so the control arm is `NWETSOIL = 0` on
the same binary rather than a code fork. The two saturation endpoints are
checked against `surface_layer.saturation_mapping` on both sides, this config
block and `landmod.f90`'s compiled defaults, at the contract's own tolerance.

The saturated fields staged are the interstitial arm, and `--wetting-arm lekner`
with a second surface build gives the other. **The bracket is the largest
uncertainty this term declares and it is bounded without that build.** The film
arm darkens every one of the 38 class-band pairs by between 0.4036 and 0.6652 of
what the interstitial arm darkens them by, median 0.5632, so scaling the staged
gap by that range brackets the film arm's field cell by cell:

| arm | endmember swing | fall at a full layer | K at the measured attenuation |
| --- | --- | --- | --- |
| interstitial, staged | 0.030156 | 0.024633 | +0.94 |
| film, at the weakest-darkening class | 0.020059 | 0.016016 | +0.61 |
| film, at the median class | 0.016984 | 0.013468 | +0.52 |
| film, at the strongest-darkening class | 0.012171 | 0.009549 | +0.37 |

The whole bracket runs +0.37 to +0.94 K, and its far end is the ceiling. Every
value in it is under the paired-arm resolution bar, so the choice of arm cannot
be settled by a run either: what the second build would buy is the exact number
inside a bracket that is already narrower than the instrument. The exact figure
is `--wetting-arm lekner --output <scratch>`, which writes nothing the model
reads.

`exoplasim/notes/forcing-bundle-predictions.md` carries what arming it is
predicted to be worth, registered before the edits and split into the two
commits that made them, per `docs/src/pipeline/sequencing.md` A3.

## What the refusal is worth in the staged field

`evaporite` is refused, and the refusal survives every step between the
derivation and the model.

`build_surface_albedo.py` assigns the class a wetting ratio of exactly
`(1.0, 1.0)` before `analysis/soil_albedo_wetting.json` is consulted at all, so
its saturated band pair is its dry band pair bitwise, and the builder then
asserts that the ratio of every refused region is exactly 1 rather than merely
close. `wet_soil_albedo` returns the dry albedo bitwise at an empty layer and
the saturated one bitwise at a full one, because both ends are branches rather
than arithmetic; in between, two equal ends give the Kubelka-Munk transform's
round trip, which at this class's staged levels is 1e-16 in double precision and
6e-9 in single. Both are orders below the 5e-6 quantum a value reaches the model
with through `.sra`, so a refused cell returns its dry albedo at every
saturation the mixing can be handed.

`playa_clastic` is NOT refused and is the largest single class on this build's
land. The dispute is about halite crust; playa mud and fan fill are a soil.
