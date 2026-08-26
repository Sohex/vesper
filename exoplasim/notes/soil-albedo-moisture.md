# The modelled soil albedo does not know whether the soil is wet

*A worldbuilding project. Everything below is about the simulated planet Vesper
and the models that produce it: the ExoPlaSim climate column, the land water
column under it, and the albedo boundary condition this project stages into
them. Measured on 2026-08-26 against `canonical-10m-base` at T21.*

This is the declaration that closes PHYS-15. The moisture dependence of the
modelled soil albedo is ABSENT, deliberately, and this note carries the
argument, the magnitude, and the four things that would arm it.

## What is not the blocker

**Not the parameterisation.** `landmod.f90:507-508` declares `dalbclim1` for
below 0.75 um and `dalbclim2` for above, the same split
`exoplasim/scripts/build_surface_albedo.py` writes into surface codes 175 and
176, and the CLM two-band mixing is defined on exactly that pair:

    alpha_band = alpha_band_dry * (1 - S_e) + alpha_band_wet * S_e

with `S_e` the effective saturation over an albedo depth. Read from
`references/climaland/src/standalone/Soil/soil_albedo.jl`, which is the artifact
rather than a description of one; its default `albedo_calc_top_thickness` is
0.02 m and its `update_albedo!` falls back to the top layer's water content when
no layer is centred above that depth.

**Not the scalar bucket, any more.** When `notes/external-model-survey.md`
section 39 was written, ExoPlaSim's land held one store with no depth and there
was nothing to evaluate `S_e` from at all. LSHY-3 has since landed the layered
column: `config/planet.yaml`'s `land_water_column` runs `scheme: layered` at
`layers: 2` and `layer_thickness_m: [0.5, 1.0]`, `landmod` carries `dwatcl` by
layer against a per-cell capacity split, and
`pedology/config/land_column_properties.yaml`'s `saturation_mapping` is the
declared conversion from a store fill fraction to a degree of saturation. A
number of the shape `S_e` is reachable today.

**And not, mostly, the coefficients.** Section 2 is the surprise of this
declaration: several wet-against-dry albedo measurements are already held and
read in this project, including in-situ pyranometer pairs on the surface type
that covers most of this world's land.

## What the blocker is

### 1. The modelled land store cannot reach a dry soil

`saturation_mapping` states the reachable interval and its own
`what_this_bounds` states the consequence, for the thermal side, in the words
this row needs:

    sr_at_wilting_point: 0.4114
    sr_at_field_capacity: 0.7877

    "THE STORE CANNOT REACH A DRY SOIL, and that is the finding rather than a
    caveat. An empty store is the wilting point, which on this build is a degree
    of saturation near 0.41 ... it is why the dry playa and salt-crust surfaces
    cannot be given their own inertia by this route however the endpoints are
    set. Representing them needs a thin surface layer that dries below the
    wilting point, which is a change to the column and not to a constant."

`dwmax` is a plant-available capacity, so an empty store is the wilting point
and a full one is field capacity. The residual is declared zero in the same
contract, which makes that degree of saturation the effective saturation the CLM
mixing takes, so the two are the same quantity and the substitution is exact.
The reachable range of `S_e` on this modelled land is therefore

    0.4114  to  0.7877

against the 0 to 1 the mixing is defined on.

The dry endmember this project holds is at the far end of that scale: the
spectra `analysis/rock_albedo.py` and `analysis/playa_albedo.py` integrate are
prepared laboratory samples, and the field pairs in section 2 take their dry
value on a surface crust after weeks without rain. A soil at its wilting point
is neither. So arming the mixing on these ingredients would deliver, on a
modelled land cell whose store is completely empty,

    alpha = alpha_dry - 0.4114 * (alpha_dry - alpha_wet)

permanently and everywhere, swinging over only the remaining 0.376 of the
interval as the store filled and emptied. That is a LEVEL SHIFT on the modelled
land albedo carrying a seasonal term about a third of the full amplitude, and
the level shift is the larger of the two. PHYS-15 names the seasonal term; the
arithmetic delivers mostly the shift.

The size of that shift is the reason this is fatal rather than untidy. On this
world's `evaporite` class, at albedo 0.50 and the measured salt-crust wetting
ratios in section 2, the permanent shift would be **0.105 to 0.140 in albedo on
every salt-crust cell**, applied in the driest season as much as the wettest.
For scale, the entire bare-rock-against-vegetated gap this project brackets over
the whole simulated planet is 0.069.

This is the same bound that blocks the moisture-dependent soil heat capacity,
and it has the same named remedy: a thin surface layer that dries below the
wilting point. It is a change to the land column, not a coefficient.

### 2. The wet endmember is broadband where it exists, and band 2 is open

This is where the row's original framing was wrong, and the correction is worth
recording because it moves the work. `notes/external-model-survey.md` section 39
says no magnitude is quoted "because this tree ships none", meaning ClimaLand's;
what it did not check is what this project's own reference library already
holds. Four sources, all already read:

| source | surface | dry | wet | wet/dry | quantity |
| --- | --- | --- | --- | --- | --- |
| Malek et al. (1990) | thin halite crust over shallow brine, Pilot Valley playa | above 0.75 after three dry weeks; 0.64 annual mean | 0.24 | 0.32 | in-situ pyranometer pair, broadband |
| Craft, Horel (2019) | Bonneville Salt Flats halite crust | 0.45 dry summer | 0.22 under 20-40 mm of flooding | 0.489 | broadband, and it recovers to 0.32 in eight days |
| Castellani Alegria et al. (2026) | Uyuni halite pan | 0.55-0.60 in dry or warm years | 0.65 in wet or cool years | 1.08-1.18 | MODIS BRDF, twenty-year series |
| Penndorf (1956) Table 1, Sewing column | clay soil / sand / bare rich soil | 15 / 31 / 7.2 | 7.5 / 18 / 5.5 | 0.500 / 0.581 / 0.764 | luminous reflectance, 0.38-0.77 um |

Penndorf's row was read on 2026-08-26 for this note, off the render rather than
the text layer, which misreads the last wet value as 55 instead of 5.5. The
levels are Earth surfaces and do not transfer; the ratios are taken within one
column or one instrument, so the preparation and geometry cancel and the ratio
transfers, which is the discipline `analysis/playa_albedo.py` already applies to
ECOSTRESS.

So the coefficient side is much further along than the row assumed. What is
still missing on it is two things, and they are both real:

**The band split of a wet surface.** The model reads the PAIR and not the
broadband field at `NSIMPLEALBEDO = 0`, so a broadband wetting ratio cannot be
staged. Applying the dry band ratio to a wet level would assert that wetting is
spectrally flat, which `build_surface_albedo.py` already argues is false of
every material on this planet -- and it is specifically false here, because
liquid water's absorption sits above 0.75 um. This project holds the direct
evidence for that in Slater et al. (1987), read: White Sands gypsum runs 0.49 to
0.62 through the visible and near infrared and then falls to 0.414 at
1.55-1.75 um and 0.159 at 2.08-2.35 um on gypsum's structural-water bands.
Penndorf brackets band 1 and nothing brackets band 2, and under this star band 2
carries 0.617628 of the flux against band 1's 0.382372. The unmeasured band is
the larger one.

**Any pair at all for the silicate classes.** The four sources above cover salt
crust and soils. Granite, andesite and the clastics -- 60 percent of this
build's land between them -- have no wetting measurement here.

The Earth-calibrated dry and wet soil fields the CLM parameterisation ships are
not a substitute for either gap. They are a soil colour map and a
colour-to-albedo table, and only the FORM transfers, on the same rule that made
BIO-18 and PHYS-14 re-derive their albedos under this star. The sources that
would close band 2 directly were sought on 2026-08-26 and none is reachable from
this host; they are listed with their DOIs in `references/INDEX.md`.

### 3. The sign is disputed on the class it matters most for

The three salt-crust sources in section 2 do not agree, and the disagreement is
not noise. Malek and Craft measure the surface state at the time of wetting and
both find it much DARKER: 0.75 to 0.24, and 0.45 to 0.22. Castellani Alegria
composites MODIS over twenty years and finds wet or cool years BRIGHTER, 0.65
against 0.55-0.60.

They are not measuring the same thing. A water film on a crust darkens it; a wet
year also dissolves and reprecipitates the crust and keeps dust off it, which
brightens it, and a twenty-year composite of annual means cannot separate the
two. Both effects are physical and they have opposite signs on the same class.
This project has no way to separate them at present, and one saturation number
cannot carry both.

`evaporite` is 2.565 percent of this build's land and the brightest surface on
the simulated planet; `playa_clastic` is 25.795 percent and the largest single
class. Between them they are 28 percent of the land and they are the two classes
the wetting term would act on hardest. A single-signed mixing applied to all
classes would be asserting on those two something the held measurements
contradict.

### 4. The depth is 0.5 m and the term is a skin term

Layer 1 of the modelled land column is 0.5 m thick against the 0.02 m the
parameterisation names. ClimaLand's own code tolerates that by falling back to
the top layer, but what falls back is a store whose drydown timescale is the
column's rather than the skin's. Craft and Horel measure the skin's: the
Bonneville crust recovers from 0.22 to 0.32 in eight days after flooding. A
0.5 m store does not relax in eight days, so driving the albedo from it would
hold the modelled surface dark long after the simulated ground had dried, and
would miss the fast brightening entirely.

DUST-17 already owns the top emitting-layer liquid water state and forbids a
second central hydrology, so the albedo depth is a separate declaration off
DUST-17's one profile and not a number this row may define.

## The magnitude

Three figures, from loosest to tightest. All are land-mean albedo deltas, which
is the unit `scripts/error_budget.py` consumes, and all are UPPER reaches: the
realised term is these times the fraction of simulated land wet at the skin and
how wet it is, and that fraction is unknown for the reason in section 1.
`build_surface_albedo.py` recomputes the first two per build into
`albedo_report.json` under `moisture_dependence`.

**The hard ceiling, both bands, all classes.** A wetted modelled surface cannot
be darker than the open water that would cover it if the wetting went all the
way, so the substrate's own albedo minus open water's bounds the darkening
branch per region.

| quantity | value |
| --- | --- |
| land-mean substrate albedo | 0.259739 |
| open water albedo, from the export's rock table | 0.060000 |
| ceiling, land-mean | 0.199739 |
| ceiling, darkest land region (basalt) | 0.040000 |
| ceiling, brightest land region (evaporite) | 0.440000 |

**Band 1 only, all classes, from Penndorf's ratio bracket.** The broadband
albedo delta a fully wetted simulated land surface would show from band 1 alone:

| quantity | value |
| --- | --- |
| land-mean substrate band-1 albedo | 0.218108 |
| broadband delta at wet/dry 0.764 | 0.019682 |
| broadband delta at wet/dry 0.500 | 0.041699 |
| band-2 flux share carrying no bracket | 0.617628 |

**Both bands, salt crust only, from the two in-situ pairs.** Applying the
measured salt-crust wetting ratio bracket 0.32 to 0.49 to this world's
`evaporite` level of 0.50:

| quantity | value |
| --- | --- |
| per-cell albedo delta at full wetting | 0.255 to 0.340 |
| land-mean contribution, at 2.565 percent of land | 0.0065 to 0.0087 |
| permanent level shift if the mixing were armed on the current column | 0.105 to 0.140 per cell |

For scale, the two land-surface endmembers `build_surface_albedo.py` already
brackets are 0.259739 bare against 0.190358 vegetated, a land-mean 0.069381
apart, and that script argues the gap is 15 to 19 W/m2 in absorbed flux against
21 W/m2 for the entire 0.85-to-0.95 stellar sweep that produced a 33 K range.
Band 1 alone at full wetting is 28 to 60 percent of that whole gap, roughly 4 to
11 W/m2 on its scaling, with band 2 and its 0.62 flux share unbracketed on top.

This is not a term that can be dismissed on size, and section 1's level shift is
not an error that can be accepted to get it.

## What would arm it

Four preconditions, each owned somewhere else, each checkable:

1. **A modelled surface layer that dries below the wilting point**, so that
   `S_e` spans the interval the mixing is defined on rather than
   [0.4114, 0.7877]. Named by `land_column_properties.yaml`'s
   `saturation_mapping.what_this_bounds` and owned there; it is a change to the
   land column, and it is the same change the modelled soil's heat capacity
   waits on.
   The cheap alternative does not work here. Re-casting the endmember pair onto
   the model's own reachable interval -- the albedo at the wilting point and at
   field capacity, rather than at laboratory-dry and saturated -- would be a
   coefficient change rather than a column change, but no held source states a
   soil albedo at its wilting point, and section 4's eight-day drydown says the
   surface a radiation scheme sees is not at the column's water content anyway.
2. **A band-2 wetting relation**, above 0.75 um, to sit beside Penndorf's band-1
   ratio, so that a wet endmember can be derived per class from the dry spectra
   this project already integrates -- under this star, not carried across from
   Earth calibrations, in the way `analysis/vegetation_albedo.py` derives the
   canopy pair. This is the band the term is mostly made of. Blocked on the
   sources listed in `references/INDEX.md`.
3. **The sign settled for salt crust**, separating the water-film darkening that
   Malek and Craft measure from the crust-condition brightening that the Uyuni
   series composites, and a wetting pair for the silicate classes, which have
   none.
4. **An albedo depth declared off DUST-17's top emitting-layer profile**, at
   which point the depth conversion is stated once and both consumers read the
   same profile at their own depths.

Until then the field this project stages is the dry endmember and says so, and
`albedo_report.json` carries the magnitude so the omission has a size rather
than a mention.

## What the bootstrap run does

Nothing changes in what the bootstrap run integrates. Codes 174, 175 and 176
carry the same dry substrate field they carried before; the modelled soil albedo
stays constant in time and the modelled land surface keeps its dry value through
every wetting and drying cycle the land column simulates. What changes is that
the field is now declared as the dry endmember at the point where it is written,
the omission carries a per-build magnitude, and `build_surface_albedo.py`
refuses to write a field in which any land region sits at or below open water's
albedo, which is the case in which the ceiling's sign claim would be false.
