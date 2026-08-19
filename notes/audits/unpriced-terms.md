# The interior moisture source, and sea salt: two terms nothing had priced

*Worldbuilding. Vesper is an invented planet and everything below is about the
simulation of it: a toy climate model, a terrain generator and their coupling.
Every quantity named here is a modelled field of a fictional world.*

*Audited 2026-08-19, against the baseline climatology
`baseline_regular_climatology.nc` on `precarve-craton` and against the vendored
ExoPlaSim fork. Commissioned by the question "is anything else missing on the
scale of the dust coupling". Both findings were measured the same day and both
came back with the sign opposite to the one assumed: the moisture source is
about a percent of land precipitation and negative, and sea salt is the same
size as mineral dust and cools where dust warms. One thing that looked like a
first-order defect is knocked down at the end.*

Findings are tagged **[numeric]** where computed here, **[inspection]** where
read out of code or artifacts, and **[physics]** where reasoning without
checking.

---

## 1. The lake moisture source the climate model lacks is about a percent of land precipitation, and it is negative

**[inspection]** The mechanism is recorded in
`exoplasim/notes/lake-representation.md`. PlaSim has no lake. A land cell's
evaporation is capped by its own precipitation plus storage, because runoff
leaves `dwatc` and enters `driver` (`landmod.f90:1390`), a separate store that is
advected downhill and discharged at the coast, and nothing returns `driver` to
the evaporating bucket. No setting of `dwmax` fixes it.

The question this audit asked is what that is worth in the currency the carve is
decided in. The water a catchment delivers to a closed basin evaporates from the
middle of a continent in the world `hydrography/` describes and over the ocean in
the model, so what moves is where the moisture source sits, and `world_state.json`
puts the endorheic share of land at 76.2%.

### The integral

**[numeric]** Computed by `carve_verdict.py` and recorded per run in
`runoff_source.endorheic_inflow`, so these are dated readings of a generated
artifact rather than constants. Measured 2026-08-19, per Earth year:

| | km3/Earth-yr | of land P | of land runoff |
| --- | ---: | ---: | ---: |
| catchment inflow, area-weighted, unclamped | 27,459 | 10.7% | 69.1% |
| catchment inflow, clamped per basin | 28,624 | 11.2% | 72.0% |
| spill reaching the ocean | 30,960 | 12.1% | 77.9% |
| runoff generated over the lake footprint | 7,005 | 2.7% | 17.6% |
| lake net evaporation, area x (E - P) | -9,341 | -3.6% | -23.5% |
| **lake evaporation minus the model's own, over the same ground** | **-2,374** | **-0.9%** | **-6.0%** |

against land precipitation of 816.4 mm/Earth-yr and land runoff of 126.4, both
measured 2026-08-17 in `missed-couplings.md` on this climatology.

The first four lines close an identity to **0.00%**: inflow equals lake net
evaporation plus the runoff generated over the lake bed plus what spills to the
ocean, because those are the only doors out of the endorheic system at steady
state. The middle term is there because a basin's catchment includes its own
bed, so the solver's demand carries it. That identity has a right answer of zero
and is the check that makes the integral usable.

Endorheic catchments deliver 69% of the land's runoff off 76% of the land, the
endorheic share being the drier share. That is credible where the earlier figure
was not.

### Two errors, both worth recording

**The reported catchment mean is not an integral.**
`runoff_source.catchment_mean_mm_per_year` is `np.nanmean` over basins, so it
weights a 400 km2 basin the same as a 4,000,000 km2 one, and it is per Vesper
year rather than per Earth year. Multiplying it by the total catchment area gave
39,454 km3/Earth-yr, which is 97% of all land runoff arriving off 76% of the
land. The area-weighted integral is 27,459. The per-basin runoff clamp, which was
blamed for the discrepancy when this document first said so, is worth 4.2% of the
integral and is not the cause. The key now names its year and says it is
unweighted.

**Gross lake evaporation is not the moisture source.** What the atmosphere gains
by a lake existing is the lake's evaporation measured against what the model
already evaporates on that same ground, not the lake's net water demand. Those
differ in size and, here, in sign. `carve_verdict.json`'s
`penman_land_evaporation_ordering` had already measured the reason: this world's
land carries a median roughness of 0.521 m against open water's 1.5e-4, a
transfer coefficient several times larger, so the model's own evaporation exceeds
Penman on 2,484 of 4,105 land cells. The same fact retired the floor at the land
rate from the criterion. Replacing rough land with a smooth lake therefore takes
moisture out of the land atmosphere over most of this planet, and the integral
agrees: **-2,374 km3/Earth-yr, or -0.9% of land precipitation**.

### What that leaves

[physics] Through `missed-couplings.md` finding 1's amplification, `d(runoff)/
runoff = 6.5 x dP/P`, and with a recycling ratio below one, a 0.9% moisture
change is a low single-digit percentage on runoff. It is a real term and it is
not a dust-scale one, and it does not point the way the intuitive argument said
it did.

**One limitation, stated so the number is not overread.** The comparison uses
catchment-mean Penman and catchment-mean model evaporation, not values over the
lake footprint itself. Lakes sit in the low, warm, dry floor of a catchment,
where Penman is higher and the model's moisture-limited evaporation lower than
the catchment mean, so this estimator is biased against the lake and the true
term is somewhere above -0.9%. Sharpening it means the sub-grid lake fraction on
the T42 grid, which `build_surface_albedo.py --lakes` already computes for
albedo. It would have to be wrong by a factor of six to reach the size this
document first claimed.

## 2. Sea salt is the same size as mineral dust and the opposite sign

**[inspection]** As of 2026-08-19 no sea salt, sulfate or volcanic aerosol
appeared anywhere in this repository as a radiative species. Mineral dust had
been priced at +0.35 to +0.74 W/m2 and was second in the whole error budget, and
the ocean is 57% of this planet's surface.

The argument for asking was the dust episode itself. `missed-couplings.md` opens
by recording that dust's forcing had been priced carefully and its effect on the
water cycle had never been asked about, and that nobody was wrong: a seam between
two correct components went unexamined. "It lands on flux and not on the carve"
is the pre-dust reasoning, and a scattering aerosol reaches the carve by exactly
the path dust reaches it.

**Measured 2026-08-19.** The component is `aeolian/scripts/build_sea_salt.py`
with its optics in `sea_salt_optics.py` and its source function in
`sea_salt_source.py`; constants and their sources are in
`aeolian/config/sea_salt.yaml`, and the numbers below are dated readings of
`aeolian/analysis/sea_salt_baseline.json`, which regenerates.

| | all modes | spume excluded |
| --- | ---: | ---: |
| emission, Tg per Earth year | 38,684 | 7,228 |
| burden, mg/m2 ocean mean | 145.9 | 32.5 |
| optical depth, ocean mean | 0.153 | 0.068 |
| optical depth, global mean | 0.093 | 0.042 |
| **TOA shortwave forcing, W/m2 global** | **-0.675** | **-0.298** |

and across every declared bracket end the forcing runs **-0.16 to -0.89 W/m2**,
with the wet-removal lifetime carrying almost all of that spread and the scale
height and the Charnock coefficient almost none.

**Against mineral dust's +0.35 to +0.74 W/m2, this is the same magnitude with
the opposite sign.** The single-scattering albedo is 1 to within 1e-5 in both
bands, so the two-stream expression has no absorbing term at all and the sign
cannot come out positive over any surface: sea salt can only cool. Dust's
warming and sea salt's cooling are of one size, and the error budget has carried
one of them and not the other.

### Why there are two columns, and which one an Earth number is comparable with

Grythe et al.'s equation 7 has three lognormal modes and the third is centred at
a dry diameter of 30 um. Its lower tail dominates the emitted mass below the
10 um cut: modes 1 and 2 give 6.1 Pg/yr at the check's wind and the third adds
26.9. Grythe's own transport discretised the emitted spectrum into four
lognormal classes with modal radii at 80% humidity up to 8.9 um, so their
reported 8.9 Pg/yr is not the analytic integral of their own equation to 10 um.
Those drops live about five hours here and are not what a measurement network
sees. The component reports both and corrects neither.

### The checks, and that they can fail

**The Earth check is a ratio between two published source functions**, because
an absolute comparison would measure our wind treatment rather than our source
function. Monahan et al. (1986) and G13T run through the same machinery give
0.705 and 0.680 of the global production Grythe's Table 2 reports for each, and
those two ratios agree to **3.5%**. The wind treatment divides out, and what is
left is that both source functions and the mass integration are right.

**The optics carry three more.** OPAC's tabulated wet density and its growth
factor are two statements of one salt volume fraction and agree to 0.008 g/cm3,
inside the table's rounding. Inverting OPAC's volume mixing at 80% humidity
returns a refractive index of **1.3328** for the diluting medium, against water's
1.333, which appears nowhere in the calculation. And the coarse dry bin's
effective extinction efficiency in band 1 is 1.81 against the geometric-optics
limit of 2.

### Two defects found on the way, both recorded rather than worked around

**`hur` is in percent and is labelled a fraction.** Upstream ExoPlaSim's pyburn
table gives code 157 the units string "1" while PlaSim writes a percentage. No
consumer in this project had ever read `hur`, so nothing had been wrong; the
first reader of it treated the attribute as true, got 100% humidity everywhere,
and through the growth curve that is a factor of five on the optical depth.
`build_sea_salt.py` now range-checks the field instead of trusting the label.

**The shared transport solver is in advective form, not flux form.**
`build_dust.py:advect_to_steady_state` steps `u dm/dx` rather than `d(um)/dx`,
so it conserves mass only where the steering wind is non-divergent, which a
horizontal wind on a sigma surface is not. The steady-state mass residual runs
1.7 to 2.3% here for that reason and is reported as a diagnostic rather than
dressed as an identity. It still catches what it is for: an unrelaxed bin or a
sign error in the loss term shows up as tens of percent.

### What is not done

The longwave is not computed. Sea salt sits in a boundary layer whose
temperature contrast with the surface is small, and a layer at the surface
temperature has no thermal forcing to give, so the term is expected to be far
inside the shortwave bracket; expected is not measured, and this says so rather
than implying otherwise. Sulfate and volcanic aerosol are still untouched. And
the model carries ONE aerosol at a time by design -- `radmod.f90` declares
`aeroqs(8,1)`, and the fork's comment at lines 216-228 refuses the prescribed
and interactive paths at once on purpose -- so putting sea salt IN the climate
beside dust is a model change and not a second file.

---

## Knocked down, with the evidence

**The river routing is not the defect, and removing it would be worse.**
`landmod.f90:roffini` builds the routing field by iterating over the land
orography and raising every local minimum to 1 m above its lowest neighbour until
none remain (the loop from line 1216, `zoron = 1. + MIN(neighbours)`), having
first clipped below-sea-level land to zero at line 1212. Read cold, that says the
model assumes 100% of land drains to the ocean on a world where 76.2% of land
drains to a closed basin, and it looks like a first-order defect.

It is not one. Runoff has already left the evaporating bucket by the time it
becomes `drunoff`, and `driver` does not evaporate wherever it sits. A routing
that terminated in an interior cell would accumulate water in `driver` forever,
which would break the model's water conservation rather than repair its
hydrology. **The pit-fill is what keeps the budget closed.** The absent physics
is the lake, per finding 1, and the routing is downstream of that absence.

It also does not reach the carve by any other path. The project already declines
to read the model's routing: `surface_water_report.json` records
`runoff_source.used = p_minus_e` and the note that `mrro` "measures ExoPlaSim's
routing rather than inflow to our sink". The only channel from the model's
hydrology into the verdict is the P and E fields over catchment cells.

Recorded because the pit-fill loop is alarming out of context and the next reader
of `roffini` will reach for the same conclusion.

## Checked, and already tracked elsewhere

Recorded so the coverage is not overread, and so these are not re-derived. Each
is a real gap that is already a finding with a task against it: ocean heat
transport (`CLIM-16`), the cloud shortwave constants (`PHYS-11`), dust deposition
on snow and ice (`DUST-14`), the absent transpiration channel (`BIO-4`, closed,
and `WORKFLOW.md` section 4), the mixed layer depth (`notes/config-rationale.md`
and the error budget's structural items), ocean salinity through `TFREEZE`
(declared in `config/planet.yaml` 2026-08-18), T85 redistributing precipitation
into the carve verdict (`WORKFLOW.md` section 6), and the carbon cycle, which
section 4 leaves open deliberately and bounds rather than solves.

## What this audit did not cover

`minerals/`, the LPJ-GUESS PFT parameter set, the pedology weathering law's own
constants, the hydrography solver internals beyond what
`carve-criterion-terms.md` covers, and the Mie code behind
`analysis/dust_optics.json`. The last two of those were also named as uncovered
by `absent-and-inherited-physics.md`, which is now twice.

## Tasks

Tracked in `TASKS.md` and not restated here.

| finding | id |
| --- | --- |
| 1. the interior moisture source | `CLIM-26`, closed by the measurement above |
| 2. sea salt, measured | `CLIM-27`, closed by the component above |
| 2b. sulfate and volcanic aerosol, still untouched | `CLIM-28` |
