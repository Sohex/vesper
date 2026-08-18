# What physics is absent, and what is Earth's without saying so

Audited 2026-08-18. The question is content rather than process: is the physics
right for this world, are the parameters cited, and is anything significant
missing as a component.

The existing audits already cover the radiation scheme's shortwave weights, the
stellar spectrum, ozone, Orogen's gravity, the lithology and albedo table, the
four terms of the carve criterion, the grid convention, the hydrological
sensitivity and the dust pricing. This one deliberately hunts where those did
not go: the ocean, the cryosphere, the soil parameters, and the couplings that
sit between two components.

**Two findings are new physics that is absent and reachable, one is an Earth
constant inherited in silence, and two are citation gaps.** The most consequential
is the first, because the error budget records it as unfixable and the model
disagrees.

---

## 1. Ocean heat transport has a switch, and the error budget says it has none

`analysis/error_budget.json`, `other_items`:

> **no q-flux** -- "gradients too strong, ice too extensive". STRUCTURAL, and no
> cheap version exists. Declare the direction and move on.

`oceanmod.f90` carries, in the ocean namelist:

    integer :: nhdiff        = 0      ! switch for horizontal heat diffusion
    real    :: hdiffk(NLEV_OCE) = 1.E3 ! horizontal diffusion coeff. [m**2/s]

with `if(nhdiff > 0)` at line 844 doing the work, `nhdiff` broadcast at line 223
and both it and `hdiffk` in the `oceanmod_namelist` read at line 150.

**A diffusive ocean heat transport is the cheap version, and it is one namelist
key on the binary already built.** It is the standard slab-ocean treatment: it
does not reproduce gyres or overturning, but it does move heat down the
temperature gradient, which is the sign of the error the budget already declares.
So the entry is wrong in its second sentence, and the direction it tells you to
declare and move on from is testable instead.

It also happens to satisfy every one of WORKFLOW A3's four conditions for a
trustworthy A/B without any further work: same binary, differing by one namelist
key, both arms branchable off one restart, both labellable `diagnostic`. That is
the cheapest structural term in the budget to price and it is currently the only
one recorded as impossible.

Two honest qualifications, neither of which recovers the entry as written.
Choosing `hdiffk` is choosing a number, and the defensible treatment is the one
this project uses everywhere else: run the bracket and report the spread rather
than tune it. And a constant diffusivity is not the real transport, so the result
is a bound on what the missing physics is worth rather than the physics itself.
Both are reasons to bracket it, not reasons to call it unreachable.

**This is a physics-is-not-a-knob case and should be read as one.** Ocean heat
transport exists. The argument for turning it on is that it exists, not that it
improves any comparison, and if enabling it makes an agreement worse that is
information rather than a reason to switch it off again.

## 2. Dust deposition reaches the soil and never reaches the cryosphere

The aeolian component computes a deposition field. Its consumers, derived from
`needs` in `config/pipeline.yaml`, are `surface_dust`, `dust_source_fields`,
`dust_forcing`, `brine_paths` and `surface_classes`. **Nothing in the cryosphere
reads it.** The coupling runs one way only: `build_dust.py` reads snow, at line
839, purely as an emission suppressor.

That omission is one-signed and this is the wrong world to have it on. Deposition
over the band the glaciers are in, from `aeolian/analysis/dust_baseline.nc`,
area-weighted over snow-covered land between 50 and 60 degrees:

| aeolian roughness end | global | land | land 50-60 deg | snow-covered 50-60 deg |
| --- | ---: | ---: | ---: | ---: |
| smooth (`aeolian_z0` 3.0e-6 m) | 234.6 | 404.8 | 111.8 | 111.1 |
| central (1.0e-4 m) | 17.9 | 29.9 | 8.78 | 8.75 |
| rough (1.0e-3 m) | 0.002 | 0.002 | 0.001 | 0.001 |

in g/m2 per Earth year. For scale, the terrestrial dust-on-snow literature works
at end-of-season snowpack loads of about 1 to 5 g/m2 and finds broadband snow
albedo reductions of roughly 0.03 to 0.08, with radiative forcings of tens of
W/m2 during the melt season. **The central estimate here is above that range and
the smooth end is twenty times above it.**

The model has room to absorb it and no term to do so. From
`run_b014469b8091/MOST_DIAG`, under `k25v`: fresh snow overall albedo 0.538, band
1 fresh 0.745 to 0.752, band-1 aged minimum 0.4955 to 0.5006. So the snow albedo
parameterisation already carries an aging range of roughly a quarter in band 1,
driven by time since snowfall, with no dependence on what has landed on it.

**Why it matters more here than the global mean suggests.** WORKFLOW 5b's glacier
result turns entirely on summer ablation -- "cooling buys brutal winters and
barely touches the summers that control ablation" -- and darkening snow acts on
exactly that term. `notes/glacier-rough-pass.md` already labels itself a
temperature criterion with no mass balance; when that becomes a mass balance,
this belongs in it, and the field it needs already exists. Separately, 1,782 land
cells carry mean snow depth above 1 cm, so this is seasonal snow albedo across a
large fraction of land and not only a glacier question.

**And it puts a cryosphere term under the aeolian roughness bracket**, which
`dust.yaml` already flags as "not a small correction" and which spans five orders
of magnitude in the table above. That bracket was understood as controlling
emission and the direct radiative effect. It now controls a snow albedo term as
well, and nothing has said so.

## 3. Ocean salinity is nowhere declared, and Earth's arrives through the freezing point

`icemod.f90`:

    real :: TFREEZE = 271.25  ! freezing temp. for sea ice at S=34.7

`config/planet.yaml` declares no salinity. Nothing in `config/`, `notes/`,
`pedology/` or `hydrography/` declares one either, and there is no ocean salt
budget anywhere in the project. So this world's sea ice forms at Earth's
freezing point because nobody chose otherwise.

**This is the DUST-6 pattern and the world it lands on is the one where it is
least likely to be right.** The defining hydrological fact here is that the
endorheic share of land is very large, and the whole point of `brine_paths.py`
and the evaporite classes is that solutes accumulate on closed basin floors
instead of reaching the sea. A planet that routes most of its continental solute
flux into terminal basins rather than into the ocean is a planet whose ocean has
a different salt budget from Earth's, and the project computes the basin half in
detail and the ocean half not at all.

The direction is determinate even without the budget. **A fresher ocean freezes
warmer, so sea ice forms more readily and there is more of it.** Earth's 34.7 psu
gives -1.9 C; a substantially fresher ocean moves that up toward 0 C, which is a
degree or more of freezing-point shift applied to a world whose sea ice fraction
runs between 0.2% and 6.2% across the flux range and whose stellar-cycle damping
turns on the ice-albedo feedback.

The fix is not necessarily to compute the budget, which needs an ocean age and an
outgassing history this project does not have. It is to DECLARE a salinity with
its reasoning, the way `metallicity` is declared in the `star` block, and to set
`TFREEZE` from it rather than inheriting a number attached to a comment about
Earth.

## 4. The soil field that closes the loop to the climate is uncited

`pedology/config/pedogenesis.yaml`, above the `water:` block, says what the block
is for:

> Multiplied by regolith depth, this gives the millimetres of water the ground
> can actually hold ... what ExoPlaSim's dwmax bucket should be, whose overflow
> *is* its runoff. So this is the field that closes the loop back to the climate
> model.

and then gives its numbers as:

> These are standard textbook mid-range values.

`sand: 0.07`, `silt: 0.20`, `clay: 0.13`, `volumetric_capacity_organic: 0.30`.
No source. WORKFLOW section 3.6 states the convention this breaks: "Every Earth
calibration lives in `pedology/config/pedogenesis.yaml` with its source."

The values are not wrong. Volumetric available water capacity, field capacity
minus wilting point, is conventionally around 0.05 to 0.10 for sand, 0.18 to 0.22
for silt loam and 0.12 to 0.17 for clay, and the ordering with silt highest is
the standard one and is the physically right shape. What is missing is the
citation, on the one pedology parameter that reaches the climate model.

The `catena:` block below it is in the same state and matters for the same field,
because it sets the regolith depth these are multiplied by:
`frost_production_bonus: 3.0`, `slope_transport: 4.0`, `slope_fines_loss: 0.35`,
`maximum_fines_loss: 0.6`. Each carries reasoning and an order-of-magnitude
argument, which is better than nothing and is honest, but none carries a source
and none is labelled `declared` the way `dust.yaml` labels its own ungrounded
constants. `dust.yaml` is the model to copy here: every constant either cites a
paper or says it is declared and carries a bracket, and the component reports the
spread across the bracket rather than the central value alone.

## 5. CH4 and N2O are absent, and the absence is not declared

`config/planet.yaml`'s `atmosphere` block is N2, O2, Ar, CO2 and ozone. There is
no methane and no nitrous oxide, and `notes/config-rationale.md`'s `atmosphere`
section does not mention them.

This is a model limitation rather than a configuration omission: PlaSim's
longwave is Sasamori (1968) with water vapour, CO2 and ozone, so there is nothing
to set. But the config reads as a composition that was chosen, and a reader will
take the absence as a decision about the world rather than a property of the
scheme. Given that the same block already carries a careful paragraph explaining
why CO2 is prescribed and should be read as an assumption, one sentence saying
the longwave carries three absorbers and what that leaves out belongs beside it.

---

## Checked and clean

Recorded so they are not re-derived, in the spirit of
`notes/failure-modes.md`'s "things already checked and disproved".

**Gravity propagates correctly through the whole radiation scheme.** This was the
obvious place for the DUST-6 pattern to have struck the largest term in the
project, and it has not. In `radmod.f90`'s shortwave, every absorber amount
divides by `ga`: ozone at line 2183, water vapour at 2187, CO2 at 2200. Rayleigh
scattering carries an explicit `(9.80665/ga)` factor at lines 2208 and 2209. The
longwave does the same at lines 2749 and 2751. The aerosol optical depth is built
from a geometric layer thickness computed as `-dt*gascon/ga*ALOG(...)` at lines
1791 and 1866. So the 23% reduction in atmospheric column mass that follows from
1 bar at 12.81 m/s2 rather than at 9.81 is represented everywhere it should be,
and 450 ppm here really is a smaller CO2 column than 450 ppm on Earth, correctly.

**Plant-available water capacity is not a gravity case, despite looking like
one.** Field capacity and wilting point are matric-potential quantities, set by
pore geometry and surface tension, and neither depends on gravity to first order.
The operational "drained under gravity" definition does, but ExoPlaSim's land
surface is a single bucket with no water table and no drainage flux -- runoff is
overflow -- so the matric-potential reading is the right one and the Earth values
carry over. The gravity term in unsaturated flow, `K(theta)(dpsi/dz + 1)` with
`K` proportional to `g`, reaches only LPJ-GUESS's inter-layer percolation, which
does not feed the climate.

**Fire is enabled**, `firemodel "GLOBFIRM"` in `run_lpj_guess.py`, so vegetation
cover on a dry world is not being computed without a disturbance regime.

**Sea ice is thermodynamic with snow and no drift**, which is the consistent
choice beside a slab ocean, and both polar caps being land limits what ice export
would do here anyway.

**The vegetation assumption is used consistently across components.**
`build_surface_albedo.py --mode vegetated` paints everything outside
`barren_rock_classes` as vegetated, and `dust.yaml`'s source block makes the same
assumption for the same classes and flags it as generous in the same direction.
Two components sharing one wrong assumption is a smaller problem than two
components disagreeing, and they do not disagree.

## What this audit did not cover

Stated so the coverage is not overread. I did not audit `minerals/`, the
hydrography solver internals beyond what `carve-criterion-terms.md` already
covers, the LPJ-GUESS PFT parameter set, the Mie code behind
`analysis/dust_optics.json`, or the pedology weathering law's own constants
beyond their citation state. The first and last of those are the ones I would
take next, because `weathering_fluxes.py` feeds a CO2 budget that section 4 of
WORKFLOW already treats as load-bearing.
