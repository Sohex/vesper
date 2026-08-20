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
size as mineral dust and cools where dust warms. Finding 2 was extended the same
day to volcanic sulfate, which is small and bounds two sulfur sources this
project cannot compute. One thing that looked like a first-order defect is
knocked down at the end.*

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

## 2. The aerosol inventory: sea salt is dust's size and cools, and everything else is small

**[inspection]** As of 2026-08-19 mineral dust was the only aerosol this project
had ever considered. It was priced at +0.35 to +0.74 W/m2 and stood second in
the whole error budget, and the ocean is 57% of this planet's surface.

The argument for asking was the dust episode itself. `missed-couplings.md` opens
by recording that dust's forcing had been priced carefully and its effect on the
water cycle had never been asked about, and that nobody was wrong: a seam between
two correct components went unexamined. "It lands on flux and not on the carve"
is the pre-dust reasoning, and a scattering aerosol reaches the carve by exactly
the path dust reaches it.

**[numeric]** Measured 2026-08-19. Global-mean top-of-atmosphere shortwave
forcing, each across every declared bracket end:

| species | W/m2 | source |
| --- | ---: | --- |
| mineral dust | +0.34 to +0.61 | wind on 43% of the planet |
| **sea salt** | **-0.16 to -0.89** | wind on 57% of the planet |
| volcanic sulfate, passive | -0.013 to -0.032 | this world's own outgassing requirement |

The components are `aeolian/scripts/build_sea_salt.py` and
`build_volcanic_sulfate.py`, with optics in `sea_salt_optics.py` and the source
function in `sea_salt_source.py`; constants and sources are in
`aeolian/config/sea_salt.yaml` and `volcanic_sulfate.yaml`, and every number
above is a dated reading of a generated artifact.

### Sea salt is the finding

**The same magnitude as mineral dust, and the opposite sign.** Its
single-scattering albedo is 1 to within 1e-5 in both bands, so the two-stream
expression has no absorbing term at all and the sign cannot come out positive
over any surface: it can only cool. Dust's warming and sea salt's cooling are of
one size, and the budget carried one of them.

| | all modes | spume excluded |
| --- | ---: | ---: |
| emission, Tg per Earth year | 38,684 | 7,228 |
| burden, mg/m2 global mean | 88.7 | 20.0 |
| optical depth, global mean | 0.093 | 0.042 |
| TOA shortwave forcing, W/m2 | -0.675 | -0.298 |

The wet-removal lifetime carries almost all of the bracket; the scale height and
the Charnock coefficient almost none.

**Two columns, and neither corrects the other.** Grythe et al.'s equation 7 has
three lognormal modes and the third is centred at a dry diameter of 30 um. Its
lower tail dominates the mass below the 10 um cut: modes 1 and 2 give 6.1 Pg/yr
at the check's wind and the third adds 26.9. Grythe's own transport discretised
the emitted spectrum into four lognormal classes with modal radii at 80%
humidity up to 8.9 um, so their reported 8.9 Pg/yr is not the analytic integral
of their own equation to 10 um. Those drops live about five hours here.

### Volcanic sulfate is a coupling, and it is small

**The source is not a scaling.** `weathering_fluxes.py` computes the CO2
outgassing this world needs at steady state to balance its own silicate
weathering, and that is a statement about magma. Carn et al.'s satellite-measured
passive volcanic SO2 flux, 23.0 +/- 2.3 Tg/yr, scaled by that requirement -- 1.10
times Earth's -- and distributed over the arc classes the fork places about the
volcanic front, gives 12.7 Tg of sulfur per Earth year. The step refuses to run
without the weathering artifact rather than falling back to a literal.

It produces a global-mean optical depth of 0.0035 and a forcing of **-0.021
W/m2**, one to two orders of magnitude below the other two. The scale height
bracket moves it not at all, which is the confirmation that settling is
irrelevant at this size and the answer is set by wet removal.

**What that measurement is worth is mostly what it bounds.** The sulfate pathway
on this world is weak: a four-day aerosol from a source far below the wind-driven
ones. Reaching sea salt's forcing would take a sulfur flux more than ten times
Earth's. That bound covers the two sulfur sources this project cannot compute at
all:

- **Explosive eruptions.** Stratospheric sulfate is where volcanic aerosol does
  its climatic work on Earth, and it is episodic. An episodic source needs a
  frequency-magnitude distribution, which needs an eruption history;
  `docs/src/reference/no-time-axis.md` says why there is none. The arc classes carry a place,
  not a rate in time.
- **Marine biogenic sulfur.** On Earth this is the LARGEST natural sulfur
  source, above volcanic, and this project has no marine biosphere at all --
  LPJ-GUESS is terrestrial. Nothing here constrains it, and the pathway bound is
  what keeps that from being an open-ended hole.

**The bound does not reach carbonaceous aerosol**, and that is the one still
open. Smoke and secondary organics ABSORB, so their forcing per unit optical
depth is larger and of the opposite sign to both scatterers, and the pathway
argument above says nothing about them. Fire is enabled in LPJ-GUESS
(`firemodel "GLOBFIRM"`) and biogenic emissions are already in the driver, so
both are estimable once a biosphere run exists on this build. `CLIM-29`.

### The checks, and that they can fail

**The sea-salt Earth check is a ratio between two published source functions**,
because an absolute comparison would measure our wind treatment rather than our
source function. Monahan et al. (1986) and G13T run through the same machinery
give 0.705 and 0.680 of the global production Grythe's Table 2 reports for each,
and those two ratios agree to **3.5%**.

**The optics carry three more.** OPAC's tabulated wet density and its growth
factor are two statements of one salt volume fraction and agree to 0.008 g/cm3,
inside the table's rounding. Inverting OPAC's volume mixing at 80% humidity
returns **1.3328** for the diluting medium against water's 1.333, a number that
appears nowhere in the calculation. And the coarse dry bin's extinction
efficiency in band 1 is 1.81 against the geometric-optics limit of 2.

**And one after the fact.** Textor et al.'s AeroCom intercomparison puts Earth's
sea-salt burden at 7.5 Tg across 16 models, or 14.7 mg/m2, with 54% diversity,
and its residence time at about half a day. This world's spume-excluded burden
is 20.0 mg/m2 with per-bin lifetimes of 0.2 to 1.2 days. Inside the model spread
on a world with warmer ocean and stronger gravity, which is where it should be.

### Two defects found on the way, both recorded rather than worked around

**`hur` is in percent and is labelled a fraction.** Upstream ExoPlaSim's pyburn
table gives code 157 the units string "1" while PlaSim writes a percentage. No
consumer in this project had ever read `hur`, so nothing had been wrong; the
first reader of it treated the attribute as true, got 100% humidity everywhere,
and through the growth curve that is a factor of five on the optical depth. Both
aerosol scripts now range-check the field instead of trusting the label.

**The shared transport solver is in advective form, not flux form.**
`build_dust.py:advect_to_steady_state` steps `u dm/dx` rather than `d(um)/dx`,
so it conserves mass only where the steering wind is non-divergent, which a
horizontal wind on a sigma surface is not. The steady-state mass residual runs
1.7 to 2.3% for that reason and is reported as a diagnostic rather than dressed
as an identity. It still catches what it is for: an unrelaxed bin or a sign
error in the loss term shows up as tens of percent.

### What is still not done

The longwave is not computed for either new species. Both are fine scatterers in
a troposphere whose temperature contrast with the surface is modest, so the term
is expected far inside the shortwave bracket; expected is not measured, and this
says so rather than implying otherwise. And the model carries ONE aerosol at a
time by design -- `radmod.f90` declares `aeroqs(8,1)`, and the fork's comment at
lines 216-228 refuses the prescribed and interactive paths at once on purpose --
so putting any of this IN the climate beside dust is a model change and not a
second file.

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
and `docs/src/pipeline/loops.md`), the mixed layer depth (`docs/src/reference/config-rationale.md`
and the error budget's structural items), ocean salinity through `TFREEZE`
(declared in `config/planet.yaml` 2026-08-18), T85 redistributing precipitation
into the carve verdict (`docs/src/pipeline/sequencing.md`), and the carbon cycle, which
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
| 2b. volcanic sulfate, measured and bounding | `CLIM-28`, closed by the component above |
| 2c. carbonaceous aerosol, which absorbs and is not bounded | `CLIM-29` |
