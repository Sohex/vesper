# What a restarted soil column inherits

Worldbuilding. Vesper is invented, and everything below is a reading of the
vendored LPJ-GUESS-CNP source in `vendor/lpj-guess/`, not an observation of
anything. No LPJ-GUESS build or simulation was performed: the model does not
build on this tree until a baseline climatology exists, because
`framework/vesper.h` is generated from one. Every claim here is read from
source. The one thing that was executed is a `g++ -fsyntax-only` pass over the
edited translation units against a synthetic `vesper.h`, which says the changes
compile and nothing else.

## The question

`Soil::serialize` in `modules/soil.cpp` is the whole of what a resumed soil
column gets. `framework/framework.cpp` constructs a fresh `Gridcell`, runs
`landcover_init` over it, and only then calls `deserialize_gridcell`, so every
member the serializer does not stream holds whatever `Soil::init_states` gave
it. A member absent from that block is therefore one of two very different
things:

- state the restart silently drops, so the resumed run's first day is a day the
  run it continues never had; or
- something the model rebuilds from serialized state, or from a compile-time
  constant, before it reads it -- in which case adding it to the serializer is
  state-file bytes that change no result and hide the ones that matter.

The two are indistinguishable in a diff, which is why a member-by-member reading
is a check that can only differ rather than one that can fail, and why the
verdicts have to be written down once rather than re-derived. Seventy members of
the Soil class were absent. Sixty-three are absent correctly.

`biosphere/config/soil_restart_state.yaml` carries all seventy verdicts with the
file:line that settles each, and `biosphere/scripts/soil_restart_state_gate.py`
refuses on a member nothing accounts for, on a classification the serializer
contradicts, on a name the class has dropped, on a verdict carrying no evidence,
and on a non-empty `lost` block.

## Findings

### 1. Three peatland quantities are read the day before they are written

`simulate_day` runs `canopy_exchange` before `soilwater`. `hydrology_peat`,
reached from `soilwater`, is the only writer of `dmoss_wtp_limit`,
`dgraminoid_wtp_limit` and (through `update_acrotelm_co2`) `acro_co2`;
`canexch.cpp` is the only reader of all three, through `get_moss_wtp_limit`,
`get_graminoid_wtp_limit` and `get_co2`. So the simulated moss and graminoid
photosynthesis of any day runs on the water-table limits and the acrotelm CO2
that yesterday's hydrology left, and a restart resumes them at the values the
constructor set: no dessication limit, no inundation limit, and the pore-water
CO2 concentration.

The magnitudes are not marginal. `dgraminoid_wtp_limit` reaches
`graminoid_wtp_limit_lowerlimit` under total inundation, so the resumed day
silently lifts a complete shutdown of simulated graminoid photosynthesis;
`dmoss_wtp_limit` runs to 0.3. `acro_co2` was a latent difference rather than a
live one until WORLD-C4J8 moved the annual water-table average off an ordinal
`Date::next()` never produces: while `awtp` was pinned at zero,
`update_acrotelm_co2` returned `PORE_WATER_CO2` unconditionally, which is
exactly the constructed value. Unpinning `awtp` made the constructed value wrong.

All three are now serialized.

### 2. Four annual accumulators were safe only because the save point was

`dec_snowdepth`, `dthaw`, `mthaw` and `maxthawdepththisyear` are all rebuilt
from day 0 of a simulation year: `driver.cpp` resets `mthaw` and rewrites
`dec_snowdepth` from the serialized `msnowdepth` under `if (date.day == 0)`, and
`soil.cpp` resets `maxthawdepththisyear` the same way. A restart that resumes at
day 0 therefore rebuilds all four before anything reads them, and while the year
boundary was the only save point this model had, their absence was correct.

WORLD-FUJ4 made a save point on an arbitrary simulated day expressible, and the
same four stop being safe the moment one is used:

- `dec_snowdepth` is read at year end by `establish()`, and holds the
  constructed zero, so every simulated PFT with `min_snow > 0` is barred from
  establishing in the resume year.
- `dthaw` is written one day at a time with no annual reset and read at year end
  by GLOBFIRM's `fire()`, which loops over the whole year; every day before the
  restart reads as zero, so the simulated fire season is short and the burnt
  fraction is biased low.
- `mthaw` and `maxthawdepththisyear` reach only the `mald` and `maxald` outputs,
  but they report the remainder of the year as though it were the whole of it.

All four are now serialized. The classification file records that this is why:
the answer for an accumulator depends on where the save point is, and the
verdict is taken against the arbitrary-day restart because that is the stronger
question and the one the model can now be asked.

### 3. The first-call flags are not the hazard they look like

`firstTempCalc` and `firstHydrologyCalc` are set true in `Soil::init_states` and
absent from the serializer, so a resumed column does carry them into a day that
is not its first. That is not a defect, and the evidence is upstream's own code.

Every block `firstTempCalc` guards tests something beside it. The cold-start
arms are `firstTempCalc && patch.stand.first_year == date.year`; the arms next
to them are `firstTempCalc && (restart || patch.stand.clone_year == date.year)`,
commented "initialize after restarts"; and `update_from_yesterday()` is called
under `!firstTempCalc || restart || ...`. A resumed column takes the restart arm,
which recomputes exactly the deterministic quantities -- the logs of six
compile-time thermal conductivities, and `ngroundl` from the stand type -- and
leaves the prognostic ones alone. All of this is present in the vendored subtree
at its merge commit, so it is mainline behaviour and not something this project
introduced.

`firstHydrologyCalc` is simpler still: `soil.cpp:1773` is the only occurrence
outside its declaration and its initialisation, and it is a clear-on-first-use
with no test anywhere in the tree. Nothing branches on it. It is dead.

### 4. The methane model's `_yesterday` asymmetry is the design

`CH4`, `CH4_diss`, `CH4_gas` and `CO2_soil` are absent while
`CH4_yesterday`, `CH4_diss_yesterday`, `CH4_gas_yesterday` and
`CO2_soil_yesterday` are serialized, which looks like an arbitrary half of a
pair. It is not. `Soil::methane` opens each day by rebuilding today's store from
the carried copy: `CH4[ii] = CH4_yesterday[ii] + CH4_prod[ii]` and
`CO2_soil[ii] = CO2_soil_yesterday[ii] + CO2_soil_prod[ii]`, with `CH4_diss`
written on every path through `calculate_gas_ebullition`. The `_yesterday` copy
IS the carried state and today's value is derived from it, so the asymmetry is
which of the pair is prognostic, not which was overlooked.

The same argument settles `T` and `T_soil`, which come back from the serialized
`T_old` and `T_soil_yesterday`, and the layer indices `IDX`, `MIDX` and `SIDX`,
which `soil_temp_multilayer` recomputes from `ngroundl` every day before any
array access taken through them.

### 5. One derived quantity was not derived on the restart path

`init_hydrology_variables()` is the only writer of `alwhc_init`, and its call
site was guarded on `firstTempCalc && patch.stand.first_year == date.year` --
the cold-start form, with no restart arm beside it, unlike every other first-call
block in the same routine. `update_layer_water_content` reads `alwhc_init` every
day for every stand with a multilayer soil, so a resumed column derived its
liquid water fraction from zeros.

This changes no result today, and saying why is the point: `alwhc`, the only
thing `alwhc_init` feeds on that path, has no reader anywhere in the model, and
`aw_max`, which does have readers, is serialized and would recompute to the same
value from the same constants. So the finding is a latent trap rather than a
live error, and the repair is the one-line restart arm the neighbouring blocks
already have, not a new entry in the serializer. `alwhc` being written,
serialized and never read is recorded here so that a future reader of it starts
from a defined value.

### 6. Members that are never read at all

`SIDX_old` is assigned at two points in `soil_temp_multilayer` and read nowhere
in the tree. `k_CO2` and `Ceq_CO2` are computed every day by
`update_daily_gas_parameters` and never consumed, because the transport step
replaces CO2 diffusion with `co2_store`. `CH4_gas_yesterday` is written,
serialized, and read by nothing.

`CH4_oxid` and `CH4_vgc` were the only members in the gas set that
`Soil::init_states` did not touch at all, so a freshly constructed `Soil` held
indeterminate memory in both until the first `Soil::methane` write. That write
precedes every read, and only over the active layer range `[IDX, NLAYERS)`, so
a read at a shallower index would have returned whatever the allocation held.
`Soil::init_states` now sets both, and zero is derived rather than merely
defined: `CH4_vgc` is `CH4_gas_vol` over the layer volume and `CH4_gas_vol` is
initialised to zero on the line above it, and `CH4_oxid` is
`min(CH4, 0.5 * O2)` with both of those initialised to zero, so zero is what
each member's writer in `soilmethane.cpp` would produce from the state
`init_states` leaves. Neither belongs in `Soil::serialize`: within a simulated
day the write at `soilmethane.cpp` and the reads beside it share one enclosing
condition, so no read of either crosses a save point.

The vendored CNP fork leaves both uninitialised, so this is a declared
divergence from mainline and not a defect this project introduced.

## What this does not settle

The classification is read from source. It says a member is written before it is
read; it does not say the resumed run reproduces the uninterrupted one, and no
reading of source can. `biosphere/scripts/verify_lpj_restart_continuity.py` is
the half that can fail on behaviour, in three modes, and none of them has been
run: they need a compiled model and a built forcing, and neither exists until a
baseline climatology does. An absent measurement is not a pass.
