# The modelled soil albedo and the water in the ground under it

*A worldbuilding project. Everything below is about the simulated planet Vesper
and the models that produce it: the ExoPlaSim climate column, the land water
column under it, and the albedo boundary condition this project stages into
them. Measured on 2026-08-26 against `canonical-10m-base` at T21.*

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
| surface layer capacity | 6.87 mm | 6.09 to 7.65 |
| sub-wilting increment | 3.57 mm | 2.81 to 4.52 |
| `awc_mm`, what LPJ-GUESS reads | 175.0 mm | 46.6 to 242.5 |
| `evaporable_mm`, what ExoPlaSim would install as `dwmax` | 178.5 mm | 49.9 to 246.0 |
| increment over `awc_mm` | 0.0225 | 0.0137 to 0.0747 |

The two capacity columns are two quantities and not one number written twice: a
root cannot reach the water between air dry and the wilting point and a bare
drying surface can. `surface_layer_report` checks the identity that defines them
-- `evaporable_mm - awc_mm` equals the surface layer's increment, cell by cell,
to 2.8e-14 mm -- so the column gained exactly the water the surface layer can
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
should, which is under the 3.57 mm increment. The alternative -- evaluating air
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
written:

| quantity | value |
| --- | --- |
| land mean, dry | 0.261450 |
| land mean, saturated, band 1 | 0.084953 |
| land mean, saturated, band 2 | 0.118397 |
| land mean, saturated, broadband | 0.105609 |
| the swing the mixing spans | 0.155841 |

For scale, the bare-against-vegetated gap this project brackets over the whole
simulated planet is 0.069381, and `build_surface_albedo.py` argues that gap is 15
to 19 W/m2 in absorbed flux against 21 W/m2 for the entire 0.85-to-0.95 stellar
sweep that produced a 33 K range. The wetting swing is more than twice it. That
is a CEILING on the term and not an estimate: the realised term is the swing
times how wet the modelled skin is and how much of the land it covers, and the
surface layer never reaches field capacity everywhere at once.

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
   switch, so at `nwetsoil = 0` they are not staged and nothing opens them.
3. **`dwmax` installed from `evaporable_mm`.** The surface layer's capacity is
   cut from air dry, so the column's capacity is `evaporable_mm` and not
   `awc_mm`. `build_surface_soil_water.py` installs it, and it selects the
   column by the declared GEOMETRY rather than by a switch of its own: the extra
   water is water the top 0.02 m can give up, so a column without that layer
   gets `awc_mm` and every layer cut from the wilting point. `awc_mm` is what
   LPJ-GUESS reads and it does not move.

The switch itself and its five companions reach `landmod_nl` from
`surface.soil_albedo_moisture`, written unconditionally by `configure_otherargs`
and checked by `expected_namelist_keys`, so the control arm is `NWETSOIL = 0` on
the same binary rather than a code fork. The two saturation endpoints are
checked against `surface_layer.saturation_mapping` on both sides, this config
block and `landmod.f90`'s compiled defaults, at the contract's own tolerance.

The saturated fields staged are the interstitial arm. Sweeping the bracket is
`--wetting-arm lekner` and a second surface build; the two arms differ by
roughly a factor of two in the wetting ratio, which is the largest declared
uncertainty this term carries.

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
