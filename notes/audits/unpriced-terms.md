# The interior moisture source, measured, and the aerosols that are not dust

*Worldbuilding. Vesper is an invented planet and everything below is about the
simulation of it: a toy climate model, a terrain generator and their coupling.
Every quantity named here is a modelled field of a fictional world.*

*Audited 2026-08-19, against the baseline climatology
`baseline_regular_climatology.nc` on `precarve-craton` and against the vendored
ExoPlaSim fork. Commissioned by the question "is anything else missing on the
scale of the dust coupling". Finding 1 was measured the same day and came back
about a percent, with the sign opposite to the one assumed; finding 2 is
unmeasured and remains open. One thing that looked like a first-order defect is
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

## 2. Mineral dust is the only aerosol this project has ever considered

**[inspection]** No sea salt, no sulfate, no volcanic aerosol appears anywhere in
this repository as a radiative species. `grep -ril "sea.salt\|seasalt\|sulfate"`
over the tracked tree returns `notes/economic-minerals.md`,
`references/INDEX.md` and `minerals/config/downstream_prospectivity.yaml`, all of
them about evaporite chemistry and none about an atmosphere. The ocean is 57% of
this planet's surface, and the arc classes are the one place the pipeline treats
volcanism as continuing (`WORKFLOW.md` section 3.6).

**The argument for asking is the dust episode itself.** `missed-couplings.md`
opens by recording that dust's radiative forcing had been priced carefully and
its effect on the water cycle had never been asked about, and that nobody was
wrong -- a seam between two correct components went unexamined. A scattering
aerosol over the ocean changes the top-of-atmosphere flux, which changes the
temperature, which changes P and E, which reaches the carve verdict by exactly
the path dust reaches it. "It lands on flux and not on the carve" is the
pre-dust reasoning and this project has already paid to learn that it is wrong.

**Neither the size nor the sign is known, and this document is not going to
guess one.** What can be said is what it would cost to find out, and the first
estimate of that was too cheap.

**[inspection]** The model carries ONE aerosol at a time, by design. `radmod.f90`
declares `aeroqs(8,1)`, a single set of Qext, Qsca, Qback and g per band, read
from one `aerofile` at line 1105, with one `apart` and one `rhop`. The fork's own
comment at lines 216-228 states that the prescribed-dust path and the interactive
path "are mutually exclusive and radini refuses both at once: a prescribed column
and a transported one are two aerosols, and adding their optical depths would
count one of them twice". A second species with its own optics is therefore a
model change -- a second column field, a second optical set, and a summation the
fork currently forbids on purpose -- and not a second file.

What IS reusable is the offline half, and it is most of the work:
`exoplasim/scripts/mie_dust.py`, `dust_optics.py`'s band-averaging across
ExoPlaSim's 0.75 um split on the k25v spectrum, the `aerofile` format and the
conversion identity in `dust_aerofile.py`, and the surface-field plumbing that
carries a prescribed column. An emission scheme for either species is the part
that does not exist.

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
| 2. no aerosol but mineral dust | `CLIM-27` |
