# The wetland, peat and methane activation contract

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's dormant peatland and methane modules,
the declaration file that governs them, and the output tables an activated run
would write. No LPJ-GUESS build or simulation was performed for this document,
and no run this project has made has passed acceptance, so every claim here
about the vendored source is checked by reading it, and the gate re-checks each
one by probe on every invocation.

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

## 2. The hydrology repair, and what still refuses under it

Four defects below the wetland path were present in the vendored source. Each is
a statement about a named location, each has a probe in the gate, and each probe
refused activation while it still matched. All four are now repaired, and a
probe that stops matching grants nothing on its own: the corresponding
declaration in `biosphere/config/wetlands.yaml` still has to be made, it names
the artifact that closes the defect, and the gate reports a declaration made
while its probe still matches as a contradiction rather than believing it.

**Saturation is an outcome, never a source term.** In
`vendor/lpj-guess/modules/soilwater.cpp` the `is_true_wetland_stand()` block
computes each layer's saturation deficit `potential_layer[ly]` and infiltrates
`min(soil.rain_melt, total_potential)`, distributed in proportion to each
layer's remaining deficit, on the same terms as the peatland and upland branches
beside it. What it does NOT do is add the full deficit whether or not the rain
could supply it, which is what it did: `soil.rain_melt` was clamped to zero when
it fell short and every layer took its whole deficit anyway. `ifsaturatewetlands`
selected nothing about whether that water was added -- it selected only whether
`patch.wetland_water_added_today` recorded it, and it recorded the whole
`total_potential` rather than the created part -- so the module now fails on
`ifsaturatewetlands 1` rather than offering a switch with no referent.

`hydrography/config/land_water_ledger.yaml` books every crossing that has no
owner against its `unowned_source` boundary, and this path was one of them. The
diagnostic pair that counted it is gone with it: `patch.wetland_water_added_today`,
`patch.awetland_water_added`, the runoff clawback in
`vendor/lpj-guess/modules/soil.cpp` that charged the created water back against
a wetland stand's own runoff, and the `wetland_water_added` output table. None
of them can carry anything but zero once the path is bounded, and a table
reporting a permanent zero is worse than an absent one. What a wetland stand
does receive from outside its own column is PLHY-4's exchange, and WET-11 names
the diagnostic it gets reported through: `mwetland_exchange.out`.

**Runon is an exchange, not a namelist constant.** `soil.cpp` assigned
`soiltype.runon = wetland_runon`, a scalar in mm/day applied to every wetland
stand on the planet and only where water was already flowing, carrying no
catchment and no baselevel. Both the assignment and the read are gone. The
ledger books the two directions of the exchange that replaces it as
`capillary_rise` (groundwater to soil_liquid) and `drainage` (soil_liquid to
groundwater), PLHY-4 owns them, and their flux is `undeclared` -- so the term is
ABSENT rather than approximated. An absent term is one the ledger can name; a
constant is a second account of where the water came from. `wetland_runon` is 0
in every instruction file this project writes, granted or refused, so removing
it preserves the behaviour of every run made here.

**The annual water-table average runs on an ordinal the calendar reaches.** The
block in `soil.cpp` that averages `wtp` over the year is guarded on
`date.day == Date::MAX_YEAR_LENGTH - 1`. It was guarded on
`Date::MAX_YEAR_LENGTH`, which `Date::next()` in
`vendor/lpj-guess/framework/guess.h` never produces: `day` resets to 0 on the
last day of the last month, so it never exceeds `MAX_YEAR_LENGTH - 1` and the
block was unreachable on the Earth calendar and on this world's. `awtp`
therefore held the 0.0 it is initialised to, for the whole run and across
restarts because it IS serialized. The consequence was not cosmetic:
`update_acrotelm_co2` computes the acrotelm CO2 concentration as a linear
interpolation in `awtp` between the atmospheric value at -300 mm and the
pore-water value at 0, so `awtp == 0` pinned the simulated acrotelm to the
pore-water concentration everywhere, always. The block is kept rather than
deleted because `awtp` has a reader; what was wrong was the ordinal.

**The restart carries the prognostic peat hydrology.** `Soil::serialize` now
carries `Wtot`, `wtd`, `stand_water`, `mwtp`, `Frac_ice` and `rootfrac`
alongside the `wtp` and `awtp` it already had. `Frac_ice` is the CURRENT ice
fraction and is a different quantity from the `Frac_ice_yesterday` that was
serialized, which is why the latter was not a substitute; `rootfrac` is rebuilt
on day zero of each simulated year, so it matters for a restart taken anywhere
else. `biosphere/config/wetlands.yaml` carries the list under
`hydrology.required_restart_members` and the gate reads the serializer
statically against it.

The behavioural half is `biosphere/scripts/verify_lpj_restart_continuity.py`: a
whole run against one split at a simulated year boundary, compared table for
table and row for row over every year from the restart point, with a wrong
answer rather than a merely different one. It is a YEAR boundary because
`framework/framework.cpp` serializes exactly once, at the end of year
`state_year - 1`, and a restarted run resumes at day 0 of `state_year` and can
never reach that save point again. Neither a mid-year stop nor the zero-step
round trip that found world-8yyh in the climate model is expressible here, and
giving this model a save point that would make them expressible is WORLD-FUJ4.
The fixture has not been RUN: no compiled model and no built forcing exist until
a baseline climatology does, and it exits saying so rather than reporting a pass
it did not earn.

**What still refuses under the repair.** `hydrology.water_ledger` -- the one
daily exchange closing precipitation, snow, surface routing, PLHY uptake and ET,
groundwater recharge and discharge, storage and runoff. PLHY-4 owns it, it does
not exist, and the repairs above are what make its absence a NAMED absence
rather than a constant standing in for it.

## 3. The extent is a partition, and one of its classes has no source

The classification of section 1 is a partition of a cell's land into five
classes, so every class in it needs a field it is taken from and a class with
none is land the partition cannot account for.
`biosphere/config/wetlands.yaml` names one source per class, each from the
closed set in `biosphere/scripts/wetland_gate.py:EXTENT_SOURCES`, and the gate
checks each named artifact against the builds under `hydrography/data`. The set
is closed because a source named in prose cannot be checked against what
hydrography publishes.

Four of the five classes have somewhere to come from, and
`hydrography/scripts/build_wetness.py` resolves them together as one partition
rather than four masks. The cut is taken from the PERIODIC lake cycle and not
from the annual equilibrium, because the annual area lies between the cycle's
trough and its peak and taking open water from one solve and seasonal
inundation from the other counts the strandline twice; WORLD-T8I5 measured the
two paints against each other before choosing and
`notes/audits/wetness-partition-cut.md` carries the measurement. Open water is
a region under the solved lake in every time bin; seasonal inundation is one
under it in some bins, which is the closed-basin part of that class and no
more, because a floodplain's inundated area needs a
height-above-nearest-drainage CDF and a routing model and a seasonally
saturated soil is a water content rather than an area; playa is the depression
floor neither takes; persistent peat-forming land is the simulated peatland
stand's own state, which is a model output; dry mineral soil is the residual.
None of those is a wetland extent on its own, and what the four together are is
a partition with one class missing.

THE SATURATED NON-INUNDATED MINERAL CLASS HAS NO SOURCE AT ANY SUPPORT. It
needed a saturated fraction, and the saturated-area closure that would have
produced one is withdrawn in `hydrography/config/topographic_index.yaml`. The
withdrawal is the disposition that config declared for a missed score before any
fraction was computed, and the score missed;
`hydrography/scripts/build_topographic_index.py` implements no closure and
`hydrography/scripts/build_wetness.py` forms no saturated class and refuses to
run if that status ever says anything else.
`hydrography/notes/subgrid-water-table.md` section 7 measures why no narrower
criterion can license one: on every arm of both Earth bore sets the gain the
terrain half carries over the cell-mean depth is at or below the scatter the
support puts on it, and a narrower criterion has less support, not more.

Two routes reopen it, and each carries a decision procedure rather than a
hope. The first is a new score, and what it needs is SUPPORT rather than a
narrower criterion: enough grid cells carrying observations to put a gain of a
few thousandths of an AUC above the cells' own scatter, where two continents of
bores gave the withdrawn score 34 and 37, either from more cells or from an
areal saturation or inundation observation in place of point depths. Its
criterion is declared before those observations are in hand, and it needs a
depth field with a regime, which is GW-24's and MIN-6's question rather than
this one's. `hydrography/notes/subgrid-water-table.md` section 7 owns that
route and is where it is stated; this document does not restate its terms.
The second route is a reduced form declared under WET-12, which states what
dropping the class costs in claims before anything is built to it; the gate
refuses a class declared absent without both the reason and the row that
licensed it. THAT IS THE ROUTE THIS DECLARATION TAKES, and
`biosphere/notes/reduced-wetland-form.md` is where it is declared with its cost
ledger. The reduced extent is not a coarser version of the partition; it is a
smaller partition, because PALADYN, the published design that is the floor for
the peat and methane parts, keys its own wetland extent on the same TOPMODEL
saturated fraction. There is no floor under extent to fall back to.

`biosphere/config/wetlands.yaml` therefore declares
`extent.class_sources.saturated_mineral` absent, with the reason and the row
that licensed it, and declares `extent.convention_arm_rule` before any share
is written: a class whose share is a convention bracket propagates through
every downstream ledger as both arms or as neither. That rule is required by
both routes rather than by the reduced one, because `f_grad` is a convention
bracket whose upper arm carries no terrain information at all, and a consumer
taking one arm of it is reporting the depth under another name.

The support of a saturated fraction is settled and stays enforced beside all of
that. GW-26 settled that `f_sat_max` is a rank statistic over the terrain
population inside a cell, and the only such population this project holds is the
mesh regions inside a climate-grid cell. There is no per-region saturated
fraction and manufacturing one needs a sub-region hypsometry that GW-6
established does not exist. `hydrography/notes/subgrid-water-table.md` section 5
carries the argument; WORLD-D9U4 carries the decision.
`biosphere/config/wetlands.yaml` therefore declares `climate_grid` as the
support any revived closure would have to live at, the same reading as
`hydrography/config/wetness.yaml`'s `refusals.saturated_fraction_support`, and
the gate refuses every other value by name. Declaring it supplies no fraction:
it is a standing refusal about a revival, not a precondition of activation.

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

**What the fork can emit today, and what each remaining table waits on.** Four
of the fifteen have a quantity behind them. `mwtp.out` is `patch.soil.mwtp`, the
monthly water-table position, which the annual-average repair of section 2 also
made non-trivial. `mch4_diffusion.out`, `mch4_plant.out` and
`mch4_ebullition.out` are `Fluxes::CH4C_DIFF`, `CH4C_PLAN` and `CH4C_EBUL`, and
their table parameters are named for those files rather than for upstream's
`file_mch4diff`, `file_mch4plan` and `file_mch4ebull`, so that the retention
list and the model name one thing. The fork also writes `mch4.out`, the total,
which this list deliberately does not retain: a total is what acceptance is not.

`biosphere/config/wetlands.yaml` records that per table under
`acceptance.retained_output_status`, and the gate checks both halves of every
entry against `modules/commonoutput.cpp`. This is not bookkeeping.
`run_lpj_guess.py` derives one instruction row from each retained filename's
stem, so `mwtp.out` becomes `file_mwtp "mwtp.out"`, and plib rejects an
instruction file naming a parameter nothing declared: a retained table with no
emitter does not produce an empty file, it aborts the run while parsing, naming
one unknown parameter and nothing about why it is unknown. Activation now
refuses ahead of that, naming every missing table and the issue it waits on. The
reverse direction is checked too, because a table that has become available and
is still recorded as waiting is one nobody will think to ask for.

The other eleven have no quantity behind them, and each waits on a model rather
than on an output routine. `mch4_production.out` and `mch4_oxidation.out` are
the closest: `Soil::CH4_prod` and `Soil::CH4_oxid` are computed per layer in the
detailed peat path and never accumulated anywhere, but the simplified inundated
path has no production or oxidation term at all -- its CH4:CO2 constant IS an
emission factor with oxidation folded in -- so a table emitted now would carry a
real number over the detailed path and a zero over the simplified one, beside a
non-zero emission. That is WET-6's `process.oxidation_separate`, and it has to
land before either table means anything. `mch4_drysoil_uptake.out` and
`mch4_aquatic.out` are WET-8: the module has no aerobic dry-soil sink and no
inland-water source. `mwetland_exchange.out` is PLHY-4's exchange, the same
absence `hydrology.water_ledger` refuses on. `awetland_area.out` and
`awetland_area_native.out` are WET-2 and WORLD-D9U4: there is no mutually
exclusive wetness classification and no declared support for a saturated
fraction. `apeat_stock.out` is WET-4, which owns the depth bracket a stock is
reported against. The three `aresidual_*.out` tables are the ledgers themselves,
and a residual computed over a graph with undeclared terms is not a residual.

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
   off, and add the restart-continuity fixture. DONE. What remains of step 1 is
   `hydrology.water_ledger`: PLHY-4's one daily exchange, which the repairs turn
   from a constant standing in for it into a named absence.
2. Name a source for each of the five extent classes, and settle the saturated
   non-inundated mineral class, whose source is withdrawn. DONE, by the reduced
   form of section 3: four classes are resolved as one partition and the fifth
   is declared absent with its reason and its licence. What remains of step 2 is
   the source contract that would let the simulated stand be GIVEN those
   fractions and refuse a fallback-equal one, which is WET-1, and the closure
   against BIO-11's rootable surface, which is stated as an identity in
   `hydrography/config/wetness.yaml` and cannot be asserted until BIO-11
   exists.
3. Place the wetland traits in PCAR-5's registry once it exists, and declare the
   peat age and depth brackets.
4. Separate production from oxidation and transport, add the dry-soil sink and
   the aquatic sources, and emit one closed surface ledger.
5. Close that ledger through BVOC-6's reduced oxidant and lifetime bracket
   before anything changes the prescribed atmospheric methane.
6. Only then run the pre-registered matched cases, with explicit permission.

Nothing in steps 2 through 6 may be attempted while any part of step 1 stands,
and the gate is what makes that a refusal rather than an intention.
