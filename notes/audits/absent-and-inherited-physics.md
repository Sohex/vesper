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
is the first, because the error budget recorded it as unfixable and the model
disagreed.

Four are settled and one is open: `dust-14`, finding 2, dust deposition never
reaching the cryosphere. Finding 1 settled somewhere neither the audit nor the
budget entry expected, and finding 5 was overtaken twice.

---

## 1. Ocean heat transport has a switch, and the error budget says it has none

`analysis/error_budget.json`, `other_items`:

> **no q-flux** -- "gradients too strong, ice too extensive". STRUCTURAL, and no
> cheap version exists. Declare the direction and move on.

`oceanmod.f90` carries, in the ocean namelist:

    integer :: nhdiff           = 0   ! switch for horizontal heat diffusion
    real :: hdiffk(NLEV_OCE)  = 1.E3  ! horizontal diffusion coeff. [m**2/s]

at `:44` and `:53`, with `mksst` acting on `nhdiff > 0`, both broadcast at
`:274` and `:278`, and both in the `oceanmod_namelist` read at `:188`.

**A diffusive ocean heat transport is the cheap version, and it is one namelist
key on the binary already built.** It is the standard slab-ocean treatment: it
does not reproduce gyres or overturning, but it does move heat down the
temperature gradient, which is the sign of the error the budget already declares.
So the entry is wrong in its second sentence, and the direction it tells you to
declare and move on from is testable instead.

It also happens to satisfy every one of `docs/src/pipeline/sequencing.md` A3's four conditions for a
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

**The budget entry is corrected, and TWO THINGS ABOUT THIS TERM MOVED AFTERWARDS
that the audit could not have seen.**

First, the arms that ran were mislabelled. `hdiffo` divided the requested
coefficient by a compiled `parameter(PLARAD=6.371E6)`, Earth's radius, on an
angular grid where the radius is what supplies the metric, so every CLIM-16 and
CLIM-19 arm realised `(PLARAD/6.371E6)**2 = 1.4401` times the diffusivity it
declared: the 300/1000/3000 bracket really ran 432/1440/4320. **world-mll**
removed the parameter, so `hdiffo` takes the radius from `planet_nl` through
`pumamod` and `oceanini` aborts if `nhdiff > 0` reaches it with no radius
(`oceanmod.f90:290-292`). The conservation claim on that result does not depend
on the coefficient's value and is untouched; the labels are, and `world-vho` is
the relabelling.

Second, and this is the part that changes the finding's verdict rather than its
arithmetic: **the term is now STRUCTURALLY PERMANENT rather than merely
unpriced.** `nfluko` is ExoPlaSim's actual q-flux, and both of its branches
RELAX the modelled state toward a sea surface temperature and sea-ice
climatology. This world has none -- an SST field is what the model PRODUCES --
so there is nothing to relax toward, and the fields `nfluko` would have used
were built from a `-999` sentinel. **world-6fh and world-4ba** made `icemod` and
`oceanmod` refuse `nfluko` on that test rather than relax toward a constructed
field. So diffusion is the bound this term gets, and a q-flux is not available
here at any price. That is a stronger statement than the budget's original one
and it is reached by a different route: not "no cheap version exists" but "the
mechanism needs an observation this world cannot have".

## 2. Dust deposition reaches neither the soil nor the cryosphere

The aeolian component computes a deposition field. Its consumers, derived from
`needs` in `config/pipeline.yaml`, are `surface_dust`, `dust_source_fields`,
`dust_forcing`, `brine_paths` and `surface_classes`. **Nothing in the cryosphere
reads it.** The coupling runs one way only: `build_dust.py` reads snow, at line
839, purely as an emission suppressor.

That omission is one-signed and this is the wrong world to have it on. Deposition
over the band the glaciers are in, area-weighted, re-measured 2026-08-24 on the
current `aeolian/analysis/dust_baseline.nc` and the bootstrap climatology's own
`snd` and `lsm`:

| aeolian roughness end | global | land | land 50-60 deg | snow-covered 50-60 deg |
| --- | ---: | ---: | ---: | ---: |
| low | 39.50 | 67.41 | 19.38 | 19.56 |
| central | 16.81 | 28.53 | 8.53 | 8.70 |
| high | 3.76 | 6.32 | 1.97 | 2.05 |

in g/m2 per Earth year. The arms are the per-lithology roughness mosaic, one
erodible-area weighted geometric mean per end, and the artifact records the
roughness each ran at. An earlier table here read the same columns off the
single-scalar arms, whose low end sat below the drag partition's own clamp; the
central column is essentially unmoved and the low end fell by about six, so the
bracket on this term is narrower than it was and does not reach zero at the rough
end any more.

For scale, the terrestrial dust-on-snow literature works at end-of-season
snowpack LOADS of about 1 to 5 g/m2 and finds broadband snow albedo reductions of
roughly 0.03 to 0.08, with radiative forcings of tens of W/m2 during the melt
season. **Every arm here is at or above that range in the units it reports.**

**Those are not the same units, and that is the finding rather than a caveat.**
The offline chain reports a FLUX in g/m2 per Earth year; the albedo response
takes a surface-layer LOAD in g/m2. Converting one into the other needs a
reservoir: dust accumulating in the layer the light sees, diluted by fresh
snowfall, concentrated as melt removes the snow around it, and lost when the pack
goes. That reservoir is the carrier this coupling needs, and nothing in the model
holds it.

**There is nowhere to receive it, at either end.** Checked against the code
rather than assumed:

- Land snow albedo, `landmod.f90`: `zalbsnow` interpolates linearly between
  `albsmin` and `albsmax` on surface temperature over 263.16 K to `tmelt`, per
  radiation band, then blends toward the background albedo by snow depth. Two
  declared constants and a temperature ramp. There is no grain size, no snow age,
  no impurity term and no argument slot for one.
- Sea-ice albedo, `seamod.f90`: the same shape, `albicemn` to `albice` ramped by
  `dicealbdt`. Snow on sea ice has no albedo distinct from the ice.
- The model's snow state is one bulk depth `dsnowz` and one temperature
  `dsnowt`. There is no layering, so there is no surface layer for a dust load to
  live in, and no restart record that could carry one.
- The offline glacier pass, `notes/glacier-rough-pass.md`, is a warmest-month
  mean below freezing test at elevation. It has no albedo term and no melt energy
  balance, so a darkening changes nothing in it.
- The background albedo boundary field is the SNOW-FREE ground albedo, and
  `landmod` blends away from it as snow deepens, so the term cannot be smuggled
  in there either: it would wash out exactly where the snow is.

So this is not a wiring gap. The quantity is of the right kind, unlike
`world-hfsm`'s optical depths against a nutrient ledger, but it is a flux where
the receiver would need a reservoir, and the reservoir, its restart record and
the albedo function that would read it all have to be built. A peer at EMIC cost
has built exactly that: `references/climber-x/src/smb/smb_params.f90` carries
`lsnow_dust`, `w_snow_dust` and `dust_con_scale` as namelist switches, so this is
a settled design elsewhere rather than an open one.

The model has room to absorb it. From `run_b014469b8091/MOST_DIAG`, under `k25v`:
fresh snow overall albedo 0.538, band 1 fresh 0.745 to 0.752, band-1 aged minimum
0.4955 to 0.5006. So the snow albedo parameterisation already carries an aging
range of roughly a quarter in band 1, driven by time since snowfall, with no
dependence on what has landed on it.

**Why it matters more here than the global mean suggests.** `docs/src/pipeline/state.md` section 5b's glacier
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
emission and the direct radiative effect. It also controls a snow albedo term,
and nothing had said so.

**OPEN, as `dust-14`, and it is a decision rather than a defect.** Nothing in the
cryosphere reads the deposition field and the coupling still runs one way. What
is now known is that closing it is not wiring: it costs a prognostic snow-dust
reservoir, its restart record, an albedo function that reads it, and a rebuild of
every binary. The soil end the finding's title contrasts it against is not
connected either -- `pedology/scripts/phosphorus_budget.py` says in as many words
that it does not read `dust_baseline.nc`, and `ANUT-3` owns that. Deposition
currently reaches neither. Two later findings bear on the same snow albedo
and neither supplies this term: `world-nfh` fixed the canopy masking applied over
snow, and `world-e5p` made code 175 archive the albedo the radiation actually
used. Both are about the albedo the model computes from what it knows; this
finding is about a forcing the model is never told.

## 3. Ocean salinity was undeclared; its declaration remains a bracket

`icemod.f90`:

    real :: TFREEZE = 271.25  ! freezing temp. for sea ice at S=34.7

This was the state found by the audit. CLIM-17 subsequently made the inherited
Earth value visible as the explicitly DECLARED `ocean.salinity_psu: 34.7` and
made the sea-ice freezing point derive from it. There is still no ocean salt
budget, so visibility and setability do not make that value determined.

**This is the DUST-6 pattern and the world it lands on is the one where it is
least likely to be right.** The defining hydrological fact here is that the
endorheic share of land is very large, and the whole point of `brine_paths.py`
and the evaporite classes is that solutes accumulate on closed basin floors
instead of reaching the sea. A planet that routes most of its continental solute
flux into terminal basins rather than into the ocean is a planet whose ocean has
a different salt budget from Earth's, and the project computes the basin half in
detail and the ocean half not at all.

OCN-9 corrected the one-sided inference. The pre-carve exorheic delivery per
unit ocean area is about 0.51 times Earth's, but it is the zero-carve floor; the
median basin-threshold sweep is about 0.91 and the all-carve spill graph reaches
2.14 times Earth. Carving can additionally remobilise solute stored in former
terminal basins. **A fresher ocean still freezes warmer and a saltier one
colder, but drainage topology alone does not select which side of 34.7 psu this
world occupies.** The freezing-point and sea-ice mechanism remains real while
the sign becomes part of the declared salinity bracket.

Computing the budget still needs an ocean age and an outgassing history this
project does not have. The implemented fix is therefore the declared salinity
and derived `TFREEZE`; OCN-13/OCN-16 own any future salt and carbon ledgers.

**One number was not the whole of it: salinity reached the model through FOUR
compiled constants.** `notes/audits/model-earth-centrism.md` finding 10 found
`CRHOS = 1030.` labelled "at S=34.7" in both `icemod` and `oceanmod`, `CPS =
4180.` labelled sea water and carrying FRESH water's value, and `CLFI = 3.28E5`
depressed by Earth brine content, none of them reachable. **world-9hb** made all
four `icemod_nl` keys with one definition -- `oceanini` takes `prhos`, `pcps` and
`pclfi` as arguments from `icemod` rather than declaring its own copies -- and
`run_exoplasim.py` writes the first three from the declared `salinity_psu`.
`CPS` moved to sea water's 3990.34, which is a physics change. That matters for
this finding's own argument in a specific way: the snow-ice flooding threshold is
the DIFFERENCE `CRHOS - CRHOI`, so a one per cent density error was a ten per
cent threshold error, and the salinity bracket this audit says the world is least
likely to sit inside now actually moves the sea-ice feedback rather than only the
freezing point.

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
No source. `docs/src/pipeline/steps.md` section 3.6 states the convention this breaks: "Every Earth
calibration lives in `pedology/config/` with its source or an explicit
statement that it is declared."

The values are not wrong. Volumetric available water capacity, field capacity
minus wilting point, is conventionally around 0.05 to 0.10 for sand, 0.18 to 0.22
for silt loam and 0.12 to 0.17 for clay, and the ordering with silt highest is
the standard one and is the physically right shape. What was missing is the
citation, on the one pedology parameter that reaches the climate model.

**Closed under LITH-24, and the way it closed is the finding's own point made
sharper.** The block now carries Saxton and Rawls (2006) as its SOURCE and then
says what that source does and does not license: the regression is over a
CONTINUOUS mixture and explicitly excluded samples above 60 per cent clay, so a
pure endmember is outside the data it was built on, and reading endmember values
off it would be taking a number from a citation rather than from the paper --
`docs/src/practice/failure-modes.md` class 9. So the three values are DECLARED,
with the paper as the check on their ordering and magnitude rather than their
origin, and each carries a bracket: sand [0.05, 0.10], silt [0.15, 0.25], clay
[0.10, 0.18], organic [0.20, 0.40]. `build_soil.py` reports the spread across
those brackets, because this field reaches the climate and a single number would
hide what it is worth.

The `catena:` block below it was in the same state and matters for the same
field, because it sets the regolith depth these are multiplied by:
`frost_production_bonus: 3.0`, `slope_transport: 4.0`, `slope_fines_loss: 0.35`,
`maximum_fines_loss: 0.6`. Each carried reasoning and an order-of-magnitude
argument, which is better than nothing and is honest, but none carried a source
and none was labelled `declared` the way `dust.yaml` labels its own ungrounded
constants. `dust.yaml` was the model to copy, and it was copied: every one of the
catena constants is labelled DECLARED with a bracket beside it and the reasoning
kept -- `frost_production_bonus_bracket: [1.0, 10.0]` runs from no bonus at all
to the order-of-magnitude literature figure,
`slope_transport_bracket: [2.0, 8.0]` spans thinning half and twice as
aggressive -- and the component reports the spread rather than the central value
alone.

## 5. CH4 and N2O are absent, and the absence is not declared

`config/planet.yaml`'s `atmosphere` block was N2, O2, Ar, CO2 and ozone. There
was no methane and no nitrous oxide, and
`docs/src/reference/config-rationale.md`'s `atmosphere` section did not mention
them.

This was a model limitation rather than a configuration omission: PlaSim's
longwave is Sasamori (1968) with water vapour, CO2 and ozone, so there was
nothing to set. But the config read as a composition that was chosen, and a
reader would take the absence as a decision about the world rather than a
property of the scheme.

**Overtaken twice, and the second time made the absence moot.** The declaration
this finding asked for was written. Then CLIM-42 added the longwave band, so
`radmod_nl` carries `ch4` and `n2o` (`radmod.f90:190-191`, broadcast at
`:1112-1113`) and there IS something to set them to; and CLIM-43 measured what to
set, from Rugheimer et al. (2013) at the epsilon Eridani grid point, which is
this host to within 35 K. Both gases read the Sun case's values -- neither is
enhanced at this star, because the enhancement turns on below 4750 K -- so
`pCH4_bar` and `pN2O_bar` are declared at Earth's abundances and are an
ASSUMPTION in the same sense `pCO2_bar` is, with the biogenic fluxes behind them
held at modern Earth's. Priced at 0.799 W/m2 of global mean outgoing longwave
against zero, which is 40 per cent of the offline answer, and 40 per cent is
what this broadband scheme returns for CO2 as well.

**The paragraph this finding asked for was itself superseded and has been
rewritten.** `config/planet.yaml`'s preamble above the `atmosphere` key said
`radmod.f90` had no CH4 and no N2O term, four lines above the block that sets
both and prices the band. `world-wu8` rewrote it to name the CLIM-42 band and
to keep the one-signed direction argument for the absorbers that remain absent.

---

## Checked and clean

Recorded so they are not re-derived, in the spirit of
`docs/src/practice/failure-modes.md`'s "things already checked and disproved".

**Gravity propagates correctly through the whole radiation scheme.** This was the
obvious place for the DUST-6 pattern to have struck the largest term in the
project, and it has not. In `radmod.f90`'s shortwave, every absorber amount
divides by `ga`: ozone at `:2485`, water vapour at `:2491`, CO2 at `:2505`.
Rayleigh scattering carries an explicit `9.80665/ga` factor at `:2514-2515` and
`:2736`. The longwave does the same at `:3226` and `:3238`. The layer geometry
the aerosol optical depth is built on is computed as `-dt*gascon/ga*ALOG(...)`
at `:1602`, `:2004` and `:2088`. So the 23% reduction in atmospheric column mass that follows from
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
would do here anyway. What the sweep did NOT check, and
`notes/audits/model-earth-centrism.md` finding 26 later found, is that the
thermodynamics carried an Earth Arctic lead-closing scale and a 9 m thickness
clamp whose enforcement redistributed melt heat globally. Both are fixed under
`world-12c`, and `config/planet.yaml` switches the clamp off entirely at
`sea_ice_max_thickness_m: -1.0`.

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
beyond their citation state.

The vendored model's own constants were also out of scope: this audit reads the
switches and the parameters this project sets, not the Fortran behind them.
`notes/audits/model-earth-centrism.md` took the same question into the fork,
which is why several of its findings land on the same mechanisms as findings 1
and 3 here and reach further into them.
