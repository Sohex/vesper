# The wetland, peat and methane activation contract

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's dormant peatland and methane modules,
the declaration file that governs them, and the output tables an activated run
would write. No LPJ-GUESS build or simulation was performed for this document:
the model does not build on this tree until a baseline climatology exists, so
every claim here about the vendored source is checked by reading it, and the
gate re-checks each one by probe on every invocation.

The finding this rests on is `biosphere/notes/wetlands-peat-methane-audit.md`.
That document says what is true. This one says what has to be declared before
the switches may move, and `biosphere/scripts/wetland_gate.py` is the
enforcement. `biosphere/config/wetlands.yaml` is where each precondition is
either declared or carries the `undeclared` sentinel that refuses.

## 1. Four switches, four models

`run_peatland 1` plus `ifmethane 1` reads like two switches and is four models:
where wetlands are, how water arrives at them and leaves, how peat carbon and
redox make and consume CH4, and what an atmosphere does with the flux. The
vendored fork has a well-evidenced version of the third for northern Earth
peatlands, an externally prescribed stand-in for the first, a demonstrably
non-conserving version of the second, and none of the fourth.

The four switches are `run_peatland`, `ifmethane`, `ifsaturatewetlands` and
`wetland_runon`. All four are top-level `declareitem` parameters that plib
requires to be parsed, so a run instruction file that omits them fails rather
than defaults. `biosphere/scripts/run_lpj_guess.py` writes all four explicitly
into the generated instruction file, after the import, for the reason it already
writes `ifbvoc` explicitly: an inherited zero cannot be told apart from nobody
having decided.

The gate is fail-closed in one direction only. With `requested: false` it
reports what is undeclared and exits 0, because a biosphere run with peat and
methane off is a correct run. Only a request to activate can be refused.

## 2. The hydrology is not repaired, and the refusal holds until it is

Four defects below the wetland path are present in the vendored source as it
stands. Each is a statement about a named location, each has a probe in the
gate, and each probe refuses activation while it still matches. A probe that
stops matching does not by itself grant anything: the corresponding declaration
in `biosphere/config/wetlands.yaml` still has to be made, and it names the
artifact that closes the defect.

**The low-latitude wetland path creates water.** In
`vendor/lpj-guess/modules/soilwater.cpp`, the `is_true_wetland_stand()` block
computes each layer's saturation deficit `potential_layer[ly]`, clamps
`soil.rain_melt` to zero when it is smaller than the total deficit, and then
adds the FULL deficit to every layer regardless. The water added is not bounded
by the water available. `ifsaturatewetlands` does not control whether the water
is added; it controls only whether `patch.wetland_water_added_today` records it.
The recorded amount is the whole `total_potential` rather than the created part,
so even the record over-states in the case where the rain could have supplied
it. `vendor/lpj-guess/modules/soil.cpp` then tries to subtract the record from
runoff and, where runoff is smaller, carries the remainder into
`patch.awetland_water_added` -- an annual diagnostic of water this world's
simulated hydrosphere never held.

What replaces it: saturation is an OUTCOME of the water ledger, never a source
term. The wetland receives routed surface water and groundwater discharge
through the one exchange PLHY-4 requires, and it returns every loss. There is no
second withdrawal and no infiltration path that may exceed what arrived.

**`wetland_runon` is a constant, not a flux.** `soil.cpp` assigns
`soiltype.runon = wetland_runon`, a namelist scalar in mm/day applied to every
wetland stand on the planet, and only where water is already flowing. It is not
a routed quantity and carries no catchment. What replaces it is the groundwater
solver's recharge and discharge, which already resolve local baselevels,
gravity-correct hydraulic conductivity and a `sink_fraction` diagnostic on the
native mesh. `wetland_runon` stays at 0 and the declaration has to name the
exchange that supplies runon instead.

**The annual water-table average is never computed.** In `soil.cpp` the block
that averages `wtp` over the year is guarded by `date.day ==
Date::MAX_YEAR_LENGTH`. `Date::next()` in `vendor/lpj-guess/framework/guess.h`
resets `day` to 0 on the last day of the last month, so `day` never exceeds
`MAX_YEAR_LENGTH - 1` and the block is unreachable on the Earth calendar and on
this world's. `awtp` therefore holds the 0.0 it is initialised to, for the whole
run and across restarts, because it IS serialized. The consequence is not
cosmetic: `update_acrotelm_co2` computes the acrotelm CO2 concentration as a
linear interpolation in `awtp` between the atmospheric value at -300 mm and the
pore-water value at 0, so `awtp == 0` pins the simulated acrotelm to the
pore-water concentration everywhere, always. What replaces it is the ordinal
test the calendar actually reaches, and a fixture that fails when the average is
still zero after a simulated year with a non-zero water table.

**The restart state is incomplete.** `Soil::serialize` in `soil.cpp` carries
`wtp`, `awtp` and several yesterday gas stores, and does NOT carry `Wtot`,
`wtd`, `stand_water`, `mwtp`, `Frac_ice` or `rootfrac`. Four of those are
prognostic peat hydrology; `Frac_ice` is the current ice fraction, distinct from
the `Frac_ice_yesterday` that is serialized; `rootfrac` is initialised on day
zero only. An arbitrary-day restart therefore resumes with a different simulated
state than the run it continues. What replaces it is serialization of every one
of them plus an exact-continuity fixture: a run stopped and resumed mid-year has
to reproduce the uninterrupted run bit for bit, which is a check with a wrong
answer rather than a different one.

**Where these repairs go.** All four are changes to `vendor/lpj-guess`, which
this contract does not make. They are specified here so that the work is a
change with a stated target rather than a discovery, and the gate holds the
refusal until each is done and declared.

## 3. A saturated fraction is a grid-cell quantity

GW-26 settled that `f_sat_max` is a rank statistic over the terrain population
inside a cell, and the only such population this project holds is the mesh
regions inside a climate-grid cell. There is no per-region saturated fraction
and manufacturing one needs a sub-region hypsometry that GW-6 established does
not exist. `hydrography/notes/subgrid-water-table.md` section 5 carries the
argument; WORLD-D9U4 carries the decision, and it is a decision rather than an
implementation.

So the extent declaration in `biosphere/config/wetlands.yaml` carries the
support of the fraction it consumes as a field of its own, and the gate refuses
any value but `climate_grid`. A wetland extent declared on a native-mesh
saturated fraction is refused by name, not because the number would be wrong but
because there is no population it could be a fraction of.

The same section establishes that the compound topographic index transports as a
rank statistic and not as an absolute threshold: `a` carries a length, so the
whole distribution shifts with the mesh, and the shift measured between this
project's own two builds is about `ln 2`. An extent model keyed on an absolute
index threshold -- the kind CLIMBER-X calibrates on Earth at about a kilometre --
inherits that shift and returns a statement about two indices' scales. The gate
refuses a declared extent model that names an absolute index cut.

## 4. Peat accumulates, and this world has no time axis

Peat depth is a rate integrated over a duration. `docs/src/reference/no-time-axis.md`
settles that a duration is undefined here rather than unmeasured, so the terrain
cannot be asked for peatland age and no finer generation would supply it. Wania
et al. (2009b) show simulated soil carbon and heterotrophic respiration depending
strongly on spin-up duration and peatland age; Kleinen et al. (2012) reach a
present peat stock by integrating the Holocene, which is exactly the fabricated
history this project does not keep.

What an accumulation rate comes from instead is the same move
`no-time-axis.md` records: an expected-value argument over a stationary
population, or an explicitly declared bracket, and never a single number. The
declaration therefore holds `peat.age_bracket` and `peat.depth_bracket` as
TWO-ENDED brackets whose ends must differ, and the gate refuses a scalar. A
scalar peat age is a claim to a history, and the gate says so by name rather
than accepting it as a central value with an implied uncertainty.

The consequence for the model is that the active methane-producing column and
the total peat inventory are different quantities. The vendored fork fixes three
0.1 m acrotelm layers over twelve 0.1 m catotelm layers and never lets
accumulated carbon change the column depth, the hydraulic properties or the
acrotelm boundary. That fixed 1.5 m column is admissible as the active-layer
process domain under a declared bracket; it is not admissible as a prediction of
this world's peat depth, and the declaration keeps the two apart.

## 5. Wetland vegetation is a trait set, and the registry does not exist

`vendor/lpj-guess/data/ins/wetlandpfts.ins` defines Earth boreal and tundra
woody types, a moss group, and wet graminoids described by Earth genera. Its
fixed water-table thresholds, snow and degree-day limits, root distributions,
aerenchyma flag and tiller geometry determine both the simulated vegetation and
the gas transport out of it. Those are traits, and traits belong in PCAR-5's
provenance-stamped covarying registry alongside every other live-plant trait,
because the covariance is the thing that must not be broken by importing one
bundle of Earth taxa.

PCAR-5 has no registry yet. So this contract does what the BVOC contract did
when it met the same wall: it states the placement, declares the fields, and
leaves the placement itself blocked on the row that owns the registry.

The contract is that the registry must distinguish at least the bog and fen
nutrient strategies, peat moss, emergent aerenchymatous plants, wet
mineral-soil vegetation, and non-vegetated inundation, which is a state with an
area and no plants and is not a vegetation type. Roots, water-table tolerance,
exudation, litter quality and gas-transport traits couple to PLHY and SDEC
through the same registry rather than through a second copy. `wetlandpfts.ins`
is retained as ONE named transplanted-Earth bracket, kept because keeping it is
the honest central case and keeping it silently is not. Vesper traits are not
created by copying an Earth taxon's values into a new name, and not by tuning
until a methane total looks right.

## 6. What an accepted artifact retains, and the brackets fixed in advance

A plausible global methane total can arise from compensating errors in wetland
area, production ratio, oxidation and transport. WETCHIMP found participating
models spanning roughly 8.6 to 26.9 million km2 of wetland extent and 141 to 264
Tg CH4 per year while sharing climate forcing; Kallingal et al. found eleven
LPJ-GUESS methane parameters with strong equifinality, where the posterior total
flux improves while the diffusion, plant and ebullition partitions shift
sharply. Aggregate agreement does not identify the components. So acceptance is
not a total.

An accepted run retains, together: the native and aggregated surface fractions;
the daily water table and the routed and groundwater exchange that set it; the
peat depth and age bracket and the C-N-P stocks; production, oxidation,
diffusion, plant transport and ebullition as separate fluxes rather than a net
emission; the dry-soil sink; the aquatic sources; and the water, C-N-P and
atmospheric methane residuals. `biosphere/config/wetlands.yaml` lists them and
the gate's `--check-run` rejects a run directory that is missing any of them,
carries a non-finite or impossible value, covers fewer simulated cells than the
productivity output, or has no forcing hash in its manifest.

The matched cases are pre-registered before any of them can be run: a
no-methane arm, a source-only arm with no atmospheric feedback, an extent arm at
both ends of the extent bracket, a redox and production arm, a transport
partition arm, and an atmospheric-feedback arm. Registering them here is what
stops the arm being chosen after its result is seen.

**The bracket floors, and where each comes from.** These are minimum widths on
declared uncertainty, as multiplicative factors, and the gate refuses a declared
bracket narrower than its floor. They are not central values and nothing is
computed from them. A narrower bracket is a claim to know this world better than
the read Earth studies know Earth.

- Wetland extent, from WETCHIMP's 8.6 to 26.9 million km2 across models sharing
  forcing.
- Wetland methane flux, from WETCHIMP's 141 to 264 Tg per year.
- Dry-soil methane sink, from Curry's 9 to 47 Tg per year around a central 28.
- The CH4:CO2 production ratio, from the two constants the fork itself carries:
  0.085 for the detailed northern path after Wania et al. and 0.027 for the
  simplified inundated path, whose own source comment says it was changed to
  match global emissions. The spread between a mechanistic ratio fitted at seven
  northern sites and an emission factor fitted to a global total is the model-form
  spread, and it is the floor.

The gate reads each floor from the declaration and compares it against
`declared_bracket_factor`. The floor can only ever refuse a declaration; it can
never supply one.

## 7. The order

1. Repair the four hydrology defects of section 2 with peat and methane still
   off, and add the restart-continuity fixture.
2. Settle WORLD-D9U4, then derive the mutually exclusive surface fractions on
   the support it chooses, closing against BIO-11's rootable surface.
3. Place the wetland traits in PCAR-5's registry once it exists, and declare the
   peat age and depth brackets.
4. Separate production from oxidation and transport, add the dry-soil sink and
   the aquatic sources, and emit one closed surface ledger.
5. Close that ledger through BVOC-6's reduced oxidant and lifetime bracket
   before anything changes the prescribed atmospheric methane.
6. Only then run the pre-registered matched cases, with explicit permission.

Nothing in steps 2 through 6 may be attempted while step 1 stands, and the gate
is what makes that a refusal rather than an intention.
