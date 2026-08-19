# Two terms nothing prices: the interior moisture source, and the aerosols that are not dust

*Worldbuilding. Vesper is an invented planet and everything below is about the
simulation of it: a toy climate model, a terrain generator and their coupling.
Every quantity named here is a modelled field of a fictional world.*

*Audited 2026-08-19, against the baseline climatology
`baseline_regular_climatology.nc` on `precarve-craton` and against the vendored
ExoPlaSim fork. Read-only; nothing here was fixed. Commissioned by the question
"is anything else missing on the scale of the dust coupling", and the honest
answer is that neither finding here is a new component. Both are terms that
already exist somewhere in the project and that nothing has ever priced.*

Findings are tagged **[numeric]** where computed here, **[inspection]** where
read out of code or artifacts, and **[physics]** where reasoning without
checking.

---

## 1. Lake evaporation never reaches the atmosphere, and it is a share of land precipitation rather than of planet area

**[inspection]** The mechanism is already recorded, in
`exoplasim/notes/lake-representation.md`. PlaSim has no lake. A land cell's
evaporation is capped by its own precipitation plus storage, because runoff
leaves `dwatc` and enters `driver` (`landmod.f90:1390`), a separate store that is
advected downhill and discharged at the coast, and nothing returns `driver` to
the evaporating bucket. No setting of `dwmax` fixes it. That note prices the
consequence as a latent-versus-sensible partition error over the fraction of the
planet the lakes cover.

**That is the wrong denominator, and it is why the item has stayed small.** What
a catchment delivers to a closed basin evaporates, in the world the hydrography
component describes, from the middle of a continent. In the model the same water
is exported to the ocean and evaporates there. Water is conserved either way and
the global budget closes either way; what moves is **where the moisture source
sits**, and the quantity that reaches the carve is that source measured against
land precipitation.

**[numeric]** Equilibrium lake area from `hydrography/analysis/carve_verdict.json`,
`bounds.penman.lake_area_km2`, is 13.7 million km2, which is 1.87% of the planet
and **4.36% of land**. The Penman rate over land implied by
`surface_water_report.json` -- global mean 3.027 mm/day, ocean mean 3.157, land
fraction 0.42817 -- is 2.85 mm/day; the same estimator's land mean is quoted at
3.4 mm/day in `carve_verdict.py`'s docstring. Across that pair:

| | km3 per Earth year | share of land P | share of land runoff |
| --- | ---: | ---: | ---: |
| lake evaporation at 2.85 mm/day | 14,289 | 5.6% | 35.2% |
| lake evaporation at 3.40 mm/day | 17,033 | 6.6% | 41.9% |
| land precipitation, 816.4 mm/Earth-yr | 256,748 | | |
| land runoff, 129.2 mm/Earth-yr | 40,635 | | |

So the absent moisture source is **5.6 to 6.6% of land precipitation**, sited in
the arid interiors rather than spread over the land. `world_state.json` puts the
endorheic share of land at 76.2%, which is why this lands here and would be a
detail on Earth.

**What it is worth is NOT computable from that alone, and the sign is not
clean.** [physics] Two steps stand between the source and the criterion. The
first is a recycling ratio -- what fraction of moisture evaporated over a
continental interior precipitates back onto the catchments -- and nothing in this
project has measured one. The second is that adding evaporation to a dry surface
converts sensible heat to latent, shallows the boundary layer and can suppress
convective precipitation while adding moisture, so the two effects oppose. An
argument that goes "more land moisture, therefore more land P, therefore more
runoff, therefore more carving" is a guess with a mechanism attached. Once past
those two steps, `missed-couplings.md` finding 1 applies as usual: runoff is a
small residual of two much larger fluxes, so `d(runoff)/runoff = 6.5 x dP/P`.

**Two cautions on the number itself.**

The equilibrium lake area inherits the criterion's per-basin runoff clamp
(`carve_verdict.py:591`, `max(runoff, 0)`, which fires on 57% of basins), and the
clamp is one-signed toward larger lakes. The table above is therefore high in
that respect.

And the obvious independent check does not work yet. Catchment inflow ought to
equal lake evaporation plus what the overflowing basins spill, and
`carve_verdict.json` reports `catchment_mean_mm_per_year.p_minus_e` of 80.682
over 241.9 million km2 of catchment. That figure is per Vesper year -- `year_s`
at `carve_verdict.py:582` is the orbital period, 180.655 days for this baseline
-- so it is 163.1 mm per Earth year, or 39,454 km3, which is 97% of the whole
land runoff over 76% of the land, and the endorheic share is the arid share. It
is not credible, and it is the same clamp reappearing. The check becomes a check
when the integral is redone unclamped and area-weighted over endorheic
catchments, which is cheap and is the first thing to do here.

**Why this is not the routing.** See the knocked-down section below.

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
| 1. the interior moisture source | `CLIM-26` |
| 2. no aerosol but mineral dust | `CLIM-27` |
