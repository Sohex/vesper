# 6. What happens next

Three nested loops, repeated as needed along a spatial-support ladder.

**A. The terrain loop, one turn of which is an ITERATION.** A generation on
the carve list, then the commissioning of the build it produces: hydrography
and boundary conditions rebuilt on the new terrain, bootstrap, derived fields,
baseline, verdict. Hydrography has to be rebuilt *after* the carve and before
the climate, because the carve changes the drainage the climate is integrated
over.

**The first turn is a special case, and it runs without ice.** `orogen` needs
an ice mask; the mask is written from a baseline climatology; the climatology
needs a build. That is a pass boundary rather than a cycle, and the way through
it is two passes:

1. Generate with `--glacial 0`. That is an honest null -- the world has no
   glacial erosion in it -- and it is chosen over the Earth-calibrated latitude
   ramp, which would be a guess wearing a number.
2. Commission that build to a baseline.
3. Write the mask from its climatology:
   `analysis/ice_mask_freezing_height.py --write-mask`.
4. Regenerate with `--ice-mask`.

So the first commissioning is PART of the regeneration rather than something
that precedes it. It cannot be collapsed into one pass, because Orogen has no
climate to place ice from before it makes the terrain. It cannot be bolted on
afterwards either: glacial, hydraulic and thermal erosion share one iteration
loop, and a mid-loop priority flood cuts outlets through the depressions
glaciation makes -- so a later glacial pass would leave those depressions
undrained and move the drainage network the climate was integrated over.

The block below is the COMMISSIONING, and it is the same block whether the
build is new or is being re-commissioned; a re-commissioning skips the
generation and starts from a build already in `source/`. Each step reads the
one before, and `check_consistency.py` reports which are still missing at any
point:

```bash
python exoplasim/scripts/rebuild_binaries.py --verify # CLAUDE.md rule 4. FIRST, always
python scripts/check_consistency.py                   # rule 8, the pre-run tier
python scripts/verify_entry_points.py                 # does every script still start
python exoplasim/scripts/verify_model_compiles.py     # does the model source compile
python scripts/smoke_test.py                          # the fast gate, per commit too

python hydrography/scripts/build_hydrography.py       # drainage, basins, coupling
python exoplasim/scripts/build_boundary_conditions.py # land mask, topography
python exoplasim/scripts/build_surface_albedo.py      # lithology albedo, no lakes yet
python exoplasim/scripts/build_surface_roughness.py   # z0
                                                      # soil_water_source: uniform
python exoplasim/scripts/run_exoplasim.py --purpose spinup   # BOOTSTRAP run; hours
python exoplasim/scripts/assess_convergence.py <run>  # NOT optional: it writes the
                                                      # run's status, and nothing may
                                                      # treat a run as settled without it
python exoplasim/scripts/build_climatology.py <run>   # then set baseline_climatology
python exoplasim/scripts/index_runs.py                # INDEX.json is the only record

python hydrography/scripts/surface_water.py           # lakes, now a climate exists
python exoplasim/scripts/build_surface_albedo.py --lakes <surface_water.nc>
python pedology/scripts/build_soil.py                 # soil, with no biosphere yet
                                                      # soil_water_source: pedology
python pedology/scripts/land_column_properties.py     # the contract's states, from that soil
python exoplasim/scripts/build_surface_soil_water.py  # 229, from those states
python aeolian/scripts/build_dust.py                  # needs the climatology AND the lakes
python exoplasim/scripts/dust_optics.py               # only if the spectrum moved
python exoplasim/scripts/dust_aerofile.py             # the model's aerofile
python exoplasim/scripts/build_surface_dust.py        # 1811, only if ndustrad = 1

python exoplasim/scripts/run_exoplasim.py --purpose spinup   # the BASELINE; hours
python exoplasim/scripts/assess_convergence.py <run>
python exoplasim/scripts/build_climatology.py <run>   # repoint baseline_climatology
python exoplasim/scripts/analyze_climatology.py <...> # Koppen, biomes, the report
python exoplasim/scripts/run_stellar_cycle.py         # A2 below; hours. The register
                                                      # orders it: carve_verdict needs it

python exoplasim/scripts/dust_forcing.py              # the surface forcing the verdict reads
python hydrography/scripts/carve_verdict.py           # the verdict, weighted per A2
python hydrography/scripts/export_carve_list.py       # the list Orogen consumes

python scripts/error_budget.py                        # re-rank against the new state
python scripts/world_state.py                         # LAST: it reads everything above
```

`python scripts/pipeline.py --plan carve_list` prints that list against the
current state, with what is already present marked skippable.

**How long the two `run_exoplasim.py` lines are, and why the length is declared
before the run starts.** A commissioning run is an APPROACH followed by a
PRODUCTION SPAN, and the span is what the climatology's error bar is bought
with: a window mean's interval covers only over a span of about twenty times the
model's memory time, made of batches of at least five, and batch length rather
than batch count is what has to be bought first. `lib/run_lengths.py` carries
that derivation and `exoplasim/notes/convergence-lengths.md` carries the
measurement behind it. The memory time is not measured and cannot be measured
here, so the answer is a BRACKET and the declared length is its floor: at T21
from cold, an approach of about seventy orbits and a span of eighty-four to two
hundred and nine, so a commissioning run is a hundred and fifty-four to two
hundred and seventy-nine orbits.

The declared length is a FLOOR and not a stopping rule. A run that reaches it
and is still refused by the convergence criteria is not finished; what the
declaration prevents is the opposite error, an orbit count chosen by whoever is
watching the wall clock, which is what the forty-orbit precision arm that
reversed sign when it was re-run at eighty-five cost.

**Declaring nothing and running until the instrument stops refusing was priced
and rejected.** It is the cheaper design and it would be strictly better if the
refusal were the same bar, but it is not: the convergence window is sized so the
offset criterion's slope error resolves its threshold, which on this model is
about thirty-eight orbits, while the climatology's own interval needs the twenty
times the memory time above. The instrument therefore stops refusing tens of
orbits before the error bar is bought, so it cannot stand in for the
declaration. Both are kept: the derived floor decides what to buy, and the
instruments' refusals decide whether to buy more.

The other half is **the carve gate: nothing outstanding may still move an
artifact the verdict is computed from.** It is not a ceremony about
irreversibility -- a wrong verdict is recoverable, because a build is replaced
wholesale -- and it is not a computation. **Read the open issues and decide.**
`bd ready` and `bd blocked` say what each one is waiting on.

Whether an issue still moves the verdict is a fact about what closing it would
CHANGE, and that lives in the issue's own prose; an issue can name an upstream
step and move nothing -- deferred, a bound, blocked on something that happens
after the carve. The `step:<id>` label records only where the work is FILED.
`conventions.md` records why no label traversal substitutes for reading the
issues.

Two imperatives in the block have no mechanical backstop, which is what earns
them their emphasis: `rebuild_binaries.py --verify` first, because nothing
else in the list notices a stale binary, and `world_state.py` last, because it
reads every artifact above it and an early run records a state that no longer
holds. The rest of the ordering is enforced by the register and by the
scripts' own refusals.

Everything below the climatology in the register -- the derived surface
classes, brine paths, the phosphorus budget, weathering fluxes, the
downstream prospectivity field (its tectonic half is terrain-only and
survives) -- is regenerated after the baseline and before anything
quotes it. They are not in the list above because they gate nothing in the
loop; they are consumers, and rule 7 governs them: when the climatology moves
they are worthless, and regenerating them is a step rather than a task.

**The graph enforces A2's ordering: `carve_verdict` needs
`stellar_cycle_run`.**

**The first climate run on a new terrain is a bootstrap, and its numbers are
not the baseline.** The surface fields that cannot be built without a
climatology are the lake compositing inside the albedo (174 to 176), soil
water capacity (229), which comes from a soil weathered under a climate, and
the dust fields (1811, 1801) when their model keys are set; the mask,
topography, roughness and the base albedo are pure functions of the terrain.
So the
loop is entered by running the model on the fields that do not need it, and
the run exists to produce the climatology the rest need.

**`model.soil_water_source` has to say `uniform` for the bootstrap and
`pedology` for the baseline, and the flip goes between them.** With it set to
`pedology`, `run_exoplasim.py` requires surface code 229 and refuses to start
without it -- the one field the bootstrap exists to make possible. Flipping
the key changes a value `continue_exoplasim.py` compares against the run's
manifest, and it will not resume across that, so the bootstrap has to be
finished -- converged, seasonal window run, climatology built -- before the
key moves.

**Lakes and soil water have to be in place before the run whose climatology
the verdict uses, not merely before the verdict.** Both change evaporation,
the numerator of the carve criterion, so a verdict taken on a climate that
lacked them is a verdict on the wrong evaporation. Their weights are not
equal -- `analysis/error_budget.json` ranks them -- so build both, and expect
the lakes to be what moves the answer.

**Leave `baseline_climatology` null until the BASELINE exists, and pass the
bootstrap's climatology explicitly with `--climatology` to each field builder
that needs it.** Every consumer raises rather than guessing while it is null,
which is the designed state and not a gap. Naming the bootstrap there instead
would hand it to every defaulting consumer -- `carve_verdict.py` among them --
and a verdict defaulted onto a climate with no lakes and no soil water is a
verdict on the wrong evaporation. `config/planet.yaml` says the same at the
key; repoint it at the baseline's REGULAR climatology once that run finishes.

Seasonal snapshots are written by default, and `analyze_climatology.py` needs
the orbital phase they carry. A segment run with `--no-seasonal-output` has to
be extended before it can produce a climatology, which has happened once.

**The flux is re-derived on a new terrain, not carried.** Terrain moves land
albedo, so the flux that lands the design mean moves with it. Take two
converged points spanning the target and interpolate; never extrapolate from
one, which this project has paid for twice, and count a point only if it
clears the remaining-offset test in `assess_convergence.py` rather than the
drift criteria alone. Which knob moves the flux is in loop C, and depends on
whether the calendar is locked.

**Between the bootstrap and the baseline, hold the flux while the surface
fields change.** If the two runs differ in flux as well, nothing measures
what the fields are worth: measure the surface first, then move the flux on a
slope measured with the surface in place. This separates SURFACE from FLUX;
it does not serialize forcing terms, which travel as one bundle per A3. The temperature target is a weaker
constraint than it looks -- it is itself a property of the terrain it was
derived on, and it cannot be held to better than the terms already
outstanding.

**The run count to a carve list.** The intersection carve in
[section 4](loops.md) needs BOTH bounding climates, and the cold one is a commissioning of its own because
its lakes, albedo compositing and soil water have to be built from a bare-rock
climatology; it cannot inherit the vegetated ones.

| | converged runs |
| --- | ---: |
| flux bracket, two points spanning more than 5 K, vegetated branch | 2 |
| vegetated baseline at the design flux | 1 |
| bare-rock arm: its own bootstrap, then its baseline | 2 |
| stellar cycle, per A2 | 1 |

Six. The bootstrap is one of the bracket points rather than a seventh, since
any converged climatology will do for fields that are themselves about to be
rebuilt. What gates what: the bracket points depend on nothing but the build;
each arm's baseline needs only its own arm's bootstrap climatology and the
fields derived from it; and the cycle run follows the vegetated baseline once
that baseline is FINAL by A2's test -- the bare arm gates nothing, being a
bound that moves nothing in the vegetated world. Runs not ordered by those
edges can go concurrently. Everything between the runs is functionally free --
nearly every other step in the graph costs minutes or seconds -- so plan in
runs and ignore the rest.

**Seed the arms rather than cold-starting them.** `run_exoplasim.py
--restart-from` changes only the spin-up path and not the equilibrium, so the
bare-rock arm starts from a converged vegetated restart -- the baseline's, or
a bracket point's if it launches before the baseline finishes -- and relaxes
across the endmember gap instead of from nothing. Expect that relaxation to be
slow:
a step change in forcing relaxes on a longer timescale than a settled run
drifts, measured at several times the drift fit on the one run that has done
it.

**Know which base you are on, because a pre-carve build restarts the count.**
The first verdict taken on one is iteration 1, whatever a previous line of
builds had already carved, and the overshoot measurement in section 4 applies
only where the terrain being re-verdicted is itself the carved one.

**A3. A forcing change costs a PREDICTION and an A/B, not an iteration.**
Landing one forcing change per iteration so each could be attributed sounds
careful and spends converged runs buying information that could never change
a decision: attribution cannot gate anything -- a correct term goes in
because it exists, per the physics-is-not-a-knob convention, and a wrong
verdict is recoverable because a build is replaced wholesale
([section 4](loops.md)). Attribution pays in exactly one place, diagnosing a
surprise, and that cost is CONTINGENT; serializing iterations pays it in full
every time against a risk that materialises occasionally.

So buy it the cheap way. A converged run measures a RESPONSE in kelvin and is
expensive; a short A/B measures a FORCING in W/m2 and is not: branch two arms
off one common restart, run a few orbits, difference the flux diagnostics, and
convert to kelvin with the canonical slope in `lib/sensitivity.py`. Every term
can be attributed this way while the whole bundle still costs ONE converged
run.

Five things make that A/B trustworthy, and none is optional:

- **Every forcing change lands with a quantitative prediction of its own
  effect**, stated before it is run, with what result would mean "wrong". That
  is the check-needs-a-right-answer convention (CLAUDE.md), applied to
  physics.
- **ONE PREDICTION PER COMMIT, AND IT NAMES THE ARTIFACT ITS MAGNITUDE CAME
  FROM.** A prediction registered jointly for two commits is sized from
  whichever of them the author was holding, and the other enters the sum at
  zero while looking accounted for. That is not hypothetical: batch 2
  registered +0.59 to +1.63 K of WARMING for `world-jgen` and `world-f9ig`
  together, sized from `cloud_optical_depth_bracket.json`, which prices the
  cloud optical depth and nothing else -- while
  `stephens_tables_vs_fits.json`, already in the tree, priced `world-f9ig`'s
  table swap at -20.63 to -5.32 K and named the experiment that would settle
  it. The measured value was -10.01 K, inside that bracket and an order of
  magnitude outside the registered one, with the sign reversed. Naming the
  source artifact is what makes the omission visible while the prediction is
  being written, because a source that prices half a change says so on its
  face.
- **Both arms use the SAME BINARY and differ only by a namelist key.** This is
  what the no-op-until-enabled convention is for. Two binaries would confound
  the term with the rebuild, and the low-I/O patch changes the restart layout,
  so arms built either side of it cannot share a restart at all.
- **Both arms branch from ONE restart, and the restart has to be SETTLED.**
  Run-to-run spread on a converged pair is larger than several of the terms
  being tested; a shared initial condition turns the comparison into a paired
  one. The second half of that is newer and was measured rather than reasoned:
  **a paired arm's power is set by how settled its donor is at least as much as
  by how long it runs.** The 1.6 to 2.1 K resolution bar this protocol quoted
  for a 25-orbit difference in global-mean temperature came from arms seeded
  from a 37-orbit control that was still relaxing. Seeded instead from a
  settled 210-orbit baseline, two arms of the same length carry standard errors
  of 0.017 and 0.031 K -- so that bar understated the instrument by more than a
  factor of twenty, and a term rejected against it was rejected against the
  donor's transient and not against the noise.
  A bar quoted from another pair therefore has to name the donor it was taken
  on, and an arm that finds itself unable to resolve a term should ask whether
  its restart or its length is the reason before it buys orbits.
- **The segments are labelled as diagnostics**, so a short A/B tail never
  enters a convergence window or a climatology.

**Bundle the terms, then check the SUM against the sum of predictions.** If
they agree, no attribution was needed and none was bought. If they disagree,
bisect, and only into the subset whose predictions were soft. What bundling
genuinely risks is `docs/src/practice/failure-modes.md` class 15, two errors that nearly
cancel; the mitigation is the per-term prediction above, not serialization,
because two cancelling errors hide just as well in a serial sequence nobody
predicted the size of.

**THE SUM CHECK IS ONLY AS COMPLETE AS THE LIST OF TERMS, and that list has
twice been short.** Batch 2 shipped `vdiff_lamm` and the soil heat solver's
conductivity as live changes that appeared in no prediction at all -- both live
by default, because the model takes a compiled value when `run_exoplasim.py`
writes no namelist key for it. A term with no namelist route is invisible to
the registry AND cannot be made into an arm, since an arm is one binary
differing by one key. So the list of terms is not assembled from the
predictions: it is assembled from the DIFF, every forcing-relevant commit in
the window, and a term that cannot be named as a namelist key is a term the
bundle is not ready to check. WORLD-SJJA and WORLD-5OYP are the two that were
found this way, and the bisect that found them cost a full paired 25-orbit set
that the sum check exists to avoid buying.

**A budgeting fact, and it is bounded rather than known.**
`build_surface_albedo.py` has a `modelled` mode that takes tree cover from an
LPJ-GUESS `fpc.out` instead of asserting a uniform vegetated endmember, and it
is the mode this project should end up in. **Adopting it requires the flux to
be re-derived.** The size of the change is bounded from above by the whole
lithology-to-vegetated span, because `modelled` leaves the deserts bare and is
therefore brighter than `vegetated` and darker than bare rock; the endmember
bracket measures that span on `canonical-10m-base` at T21, matched window and
matched I/O regime on both arms, at 4.40 W/m2 absorbed at the top of the
atmosphere and 2.924 K in global-mean surface temperature.

That is about seven steps of the design flux candidate grid, so the
re-derivation is not optional. It is NOT the largest forcing this project
varies on purpose: the stellar sweep from 0.85 to 0.95 spans 21 W/m2 absorbed
and 33 K, so the endmember span is a fifth of it in flux and a tenth in
kelvin. The cost to budget for is a re-derivation, not a re-ranking of what
dominates this world's energy balance, and it is the same cost whether or not
a carve shares the iteration. The corollary is the answer to "when does
the biosphere run": not before the carve, because its driver is built from a
climatology and a soil the carve replaces; after the re-baseline, on final
terrain, where its output can be adopted deliberately in the iteration after
that.

**A4. Dust is not carve-neutral, and the carve waits on one climate run.** Its
global-mean forcing is negligible in kelvin; the local effect is not. The
question splits in two and only one half is settled.

The LAKE half is done. Dust dims the lake, Penman evaporation falls, and more
basins overflow: a dust-free verdict UNDER-carves, modestly and in the
recoverable direction. That was settled offline by putting the dust SURFACE
forcing through `carve_verdict.py --dust-forcing`; the surface balance is the
operative one, because the top-of-atmosphere sign over bright basin fill is
energy retained in the atmosphere, not delivered to the ground.

The CATCHMENT half is larger, runs the other way, and cannot be reached
offline, because it is a precipitation response and not a surface energy
balance: reduced rainfall cuts the runoff the criterion divides by, which
OVER-carves, and over-carving is the expensive direction, a full terrain,
hydrography and boundary-condition rebuild to recover. It needs one
prescribed-dust climate run, specified with its gates and a prediction in
`aeolian/notes/prescribed-dust-run.md`. The ORDER here is a dependency and not
a serialization: the verdict divides by a runoff that dust moves, so it has to
be taken on a climatology that already has dust in it. Run it, take the
verdict on its climatology, then carve. **Do not record a carve list as
dust-independent.** `notes/dust.md` has the numbers and
`aeolian/scripts/dust_runoff_sensitivity.py` (a registered one-off) the
conversion from a precipitation change into basins.

**B. The soil and biosphere loop, at the operating support.** For a given climate:

1. `pedology/scripts/build_soil.py`, with no biosphere on the first pass.
2. `biosphere/scripts/build_lpj_driver.py`, then `run_lpj_guess.py`.
3. `build_soil.py --soil-carbon <run>/cpool.out`, then LPJ-GUESS again.
4. Repeat 3 until the criteria in `pedogenesis.yaml` are met.

Check whether step 3 moves anything before assuming it needs iterating:
LPJ-GUESS computes its own soil carbon internally, so the pedology organic
feedback may be second-order.

**A2. The stellar cycle, and where it belongs in the order.** The cycle run is
NOT a final flourish. It has to come after a baseline that is FINAL, not
merely converged -- a cycle is variance about a mean, and centring it on the
wrong one describes a world that will not exist, so anything still outstanding
that moves the mean has to land first. But it belongs BEFORE the carve
verdict, for a reason that is easy to miss.

Carving is irreversible: outlet incision does not undo when the warm phase
returns, so every cold, wet excursion that pushes a basin to overflow carves
it permanently. **The terrain therefore ratchets toward the state implied by
the cycle's wet extreme, not its mean**, and a verdict taken on the mean
climate systematically under-carves, because overflow is a threshold process
and the wet phase contributes disproportionately. Bedrock incision is slow
against both cycle periods, so it is the integral over many cycles that
matters: the verdict wants a climate between the mean and the wet extreme,
weighted by time spent overflowing.

**The two components ratchet differently, and that is the reason there are
two.** The medium one is too fast for ice to follow and too fast for an outlet
to incise, so it contributes only through the tail of its distribution. The
long one is slow enough that glaciers equilibrate and long enough that a wet
excursion is sustained, so it is the one that actually cuts. And because the
two periods are non-commensurate, the deepest minima differ in depth rather
than repeating, so the landscape ends up recording which past minima were
severe. A cycle run therefore has to be long enough to sample that: its length
is set in periods of the LONG component, not the medium one.

The cycle run also measures the damping factor, currently known only as a
factor-of-2.5 range, which converts any future amplitude choice into a climate
without another run -- measuring it once is worth more than the run that
measures it.

Sequence, then: baseline at the chosen mean, cycle run centred on it, verdict
on a cycle-informed climate, and only then the terrain loop.

**C. The vegetation-climate loop.** `build_surface_albedo.py --mode modelled`
turns the run's foliar cover into surface albedo and forest fraction, then the
climate runs again on it. Two checks belong here and neither is optional.

*The flux.* If the modelled land albedo differs much from the assumed value
the world may leave the design mean derived in section 5b. **Move the flux by
changing the star's luminosity, not the orbit.** The year length depends on
the semimajor axis, and the year is compiled into LPJ-GUESS, so moving the
orbit forces a rebuild and a driver regeneration. With the axis locked,
F = L/a^2, so luminosity is the free parameter and the calendar does not move.

The knobs swap at the lock, and the lock is what makes this a rule rather than
a preference. Before the calendar is fixed, the semimajor axis is the better
knob: it keeps the star itself untouched, so `k25v` stays valid, where
covering the same flux range with luminosity moves the effective temperature
far enough to leave the window the spectrum was interpolated in. Once
LPJ-GUESS is compiled against a year, that knob costs a rebuild and luminosity
is the only cheap one left. A bracket taken on the axis varies two things
rather than one, because the year moves with it, so the measured slope carries
a small seasonality effect through the calendar; that is second order for an
annual mean, and it is stated rather than implied.

Mind the scaling, and mind the spectrum. L scales **1:1** with F at fixed
orbit, and at fixed
radius L ~ T^4, so the effective temperature moves as F^(1/4):

| flux change | luminosity | effective temperature | band-1 fraction |
| --- | --- | --- | --- |
| 1.5% | 1.5% | 18 K | ~0.9% |
| 4% | 4% | 49 K | ~2.3% |

`k25v` is interpolated between BT-Settl models 100 K apart, and that grid step
is worth 4.7% in the fraction of flux below 0.75 um. So small adjustments are
free, but **past about 2 to 3% in luminosity, re-run
`build_stellar_spectrum.py` at the new effective temperature**, and past
roughly +35 K the target leaves the interpolation bracket entirely and the
pinned SVO grid points have to change with it. The spectrum exists to get snow
and ice albedo right, so letting it drift silently would undo the reason it
was built.

*The carve verdict.* It was taken on assumed vegetation. Real vegetation
changes the climate, which changes evaporation over catchments, which can
change the verdict. Re-run it; if basins flip, back to loop A.

**D. The declared escalation to T85.** T85 is the OPERATING SUPPORT. It is a
decision, not a rung to be discovered: the ladder below is the route to it, and
the question a comparison between rungs answers is how much the coarser ones
were wrong by, not which one to stop at.

The route escalates resolution and timestep ONE AT A TIME, never together --
and this route never has to move the timestep at all:

1. **T21 at dt 45.** Converge.
2. **Convert to T42, at dt 45.** Converge.
3. **Convert to T85, at dt 45.** Converge.
4. **Continue at NLOWIO = 0** with high-cadence orbits for dust.

**THE ROUTE IS RUN AT ONE STEP THROUGHOUT, AND THAT IS WHAT MAKES IT THREE
ENTRIES.** dt 45 is at or below the measured ceiling of every rung on it -- T21
refuses at 150, T42 at 90, T85 at 60 -- so no conversion has a step change in
front of it, and a settling block exists ONLY to make one. Both settling blocks
the route used to carry are gone rather than skipped, and the converter's
requirement is met trivially rather than by construction: every change of rung
happens at constant dt because nothing ever changes dt.

**IT IS AN ATTEMPT, AND THE FALLBACK IS NAMED.** What the ceilings license is
that each rung STARTS clean at 45; endurance at 45 is established only for T21,
over fifty orbits. If a rung will not carry 45 for a commissioning span, the
route falls back to the four-entry form -- reconverge T21 at 30, then both
conversions at 30 -- and the fallback is a step change, so it reinstates one
settling block at T21 and none above it. Running the cheap end of the ladder
first is what makes attempting it the right order: T21 is where a step that
cannot be held shows up for the least money.

**T42 at 45 has a blow-up on record and it does not bar the attempt.**
run_900548ae632e took a SIGFPE in its forty-seventh orbit at 45. The run is
gone -- not in `exoplasim/runs/`, not a stub under `archive/runs/`, not in any
`INDEX_AT_DELETION.json` -- so its source sha cannot be read and WORLD-TD3's
ninety-second reproducer cannot be re-run. What it measured is therefore
unknown and is certainly not this source, which has since taken the damping
correction, the `epilog` use-after-free and the batch-2 forcing terms. Refusing
a step on a profile the model no longer has is failure-modes class 34. So the
row is carried as evidence that is REPORTED and does not REFUSE:
`lib/rungs.py` marks it `binds: False` and `run_exoplasim.py` prints it at
launch. WORLD-TD3's reopen condition is unchanged -- if it recurs, it recurs on
a run that exists, the row binds again, and the ladder falls back to 30.

**What it saves, in model steps, which needs no clock.** Steps go as 1/dt, so
T85 at 45 against the 22.5 the route once carried is half the steps at the most
expensive rung, and both vanished settling blocks are whole spans on top of
that -- one at T21, one at T42. The cost of a settling block grows up the
ladder, which is why removing the T42 one was worth more than removing the T21
one.

This chapter DECIDES the route. `lib/rungs.py` carries it in machine-readable
form as `ESCALATION_ROUTE` and checks itself against the list above, so a
script asks the registry and the argument stays here.

**Why the reconvergence steps exist, and why they are not optional bookkeeping.**
A conversion changes the support and a timestep change moves the attractor, and
a state that changes both at once cannot say which one moved it. Reconverging at
the TARGET rung's timestep before converting means every conversion happens at
constant dt, so exactly one variable moves per step and a surprise after a
conversion is attributable to the support alone.

That is a hard requirement of the converter and not only good practice.
`convert_restart.py` copies `nstep`, and elapsed time is `nstep` times the step,
so a conversion across a change of step moves the planet in its orbit; the two
stored leapfrog levels are also a derivative over the donor's step, read by the
target as spanning its own. The converter refuses the combination rather than
taking it, and its self-test walks this route to check that every conversion the
route asks for is one it accepts. WORLD-FL9C.

**THIS ROUTE HAS NO SETTLING BLOCK, and the concept is kept because the
fallback needs it.** A SETTLING block and a COMMISSIONING verdict are named
apart so they cannot be confused. A commissioning verdict supports a claim about
this world's climate, and its window is priced so a slope's standard error
resolves its threshold. A settling block supports nothing: it exists to hand the
next conversion a restart that is not mid-transient after a step change, and any
residual drift it leaves is absorbed by the convergence that follows. Nothing
reads its climate. Applying the commissioning standard to one buys a claim
nothing consumes, at about four and a half times the orbits the purpose needs.
If the fallback to dt 30 is taken, the reconvergence of T21 at 30 is a settling
block in exactly this sense and is priced as one.

**A settling block is a LENGTH derived from the relaxation time, not a verdict
on a mean.** A perturbation decays as `exp(-n / tau)`, so the orbits needed for
a step change worth `A` kelvin to fall below a residual `r` is `tau * ln(A / r)`.
The residual is the same 0.15 K the offset criterion allows, deliberately: a
residual smaller than what the instrument judging the NEXT state can see is one
that state cannot be held responsible for. `lib/run_lengths.py` carries the
derivation. The relaxation time is bracketed by the fits that are evidence
rather than taken from the derived value, because that value is not established
as a bound on them, so a step change worth about half a kelvin settles in eight
to seventeen orbits -- which is where the ten to twenty this project has
repeatedly seen comes from.

**What a settling block does NOT establish**, and this is why it carries its own
name: that the state is equilibrated, that its climate is the rung's climate, or
that any mean taken on it carries an interval. Steps 1, 2 and 3 are all
COMMISSIONING and keep the full standard, because loop A is replayed on each of
those rungs and the carve verdict is taken on their climatologies.

**A commissioning convergence, at steps 1, 2 and 3, is one convergence window
then three-orbit increments until the criteria are met.** The first block is the
window because the criteria cannot be evaluated on fewer orbits than they are
taken over; three after, because that is the smallest increment the criteria can
judge and a longer one overshoots the exit by more than it costs to test again.
`assess_convergence.py` owns the window and derives it;
`exoplasim/notes/convergence-lengths.md` carries the arithmetic and the bracket
it lands in.

**The timestep at each rung is the highest that rung can carry, not a safety
margin.** The ceilings are now MEASURED rather than intended: WORLD-37TN's
re-taken grid has a refusing cell above each of them, at 150 for T21, 90 for
T42 and 60 for T85, and it is the first grid on which T85 was probed at all. A
ceiling still bounds the route rather than setting it, and it qualifies a step
against REFUSAL and against nothing else: 900 steps is about a seventh of an
orbit at dt 45, so a rung clean on the grid has said nothing about enduring a
commissioning span. Endurance lives in `COMMISSIONING_EVIDENCE`, and a step
that merely does not blow up is a stability floor rather than a licence.

Conversion uses CLIM-52's schema-aware converter with support-matched target
static fields; the result is an initial condition, not a continued equilibrium.
At each rung, settle loop A before trusting anything downstream, because changed
orography, coastline, precipitation and evaporation can change the carve
verdict; then rebuild hydrography coupling, groundwater, soil and ecological
forcing for that support and settle loops B and C.

There is no between-rung convergence comparison, and the reason is worth
stating so it is not reinstated by reflex. Such a comparison answers "which
rung is fine enough to stop at", and choosing T85 by decision removed that
question. It cannot be repurposed into "how wrong were the coarse rungs",
because that needs a converged finer reference to difference against and T85 is
the finest support being run. The coarse rungs are a ROUTE to T85, not
candidates being scored against it.

LPJ-GUESS gridcells are independent columns, so its direct cost grows roughly
with land-cell count, but the biosphere can still determine whether a climate
support is adequate: nonlinear ecology sees different forcing, soils and
partial areas at different support. A finer LPJ grid on interpolated coarse
forcing adds no information. A finer LPJ grid on a finer accepted climate can
change the result, and no generic NPP sign follows because variance, covariance,
thresholds, coastline and atmospheric feedback all move together. The fully
resolved fine reference for these aggregation checks is the ~10M-region Orogen
export (7.60 km mean edge), not the former 2.5M mesh.

**E. The cycle reaching the biosphere, last, on the settled world.** A2 says
when the cycle RUNS; what is left for the end is letting LPJ-GUESS see it.
Build the climatology with `--per-year` and pass the sequence to
`build_lpj_driver.py --climatology y0.nc y1.nc ...`. That matters because
productivity responds annually and saturates, so a run on the cycle mean
over-predicts it; it does *not* reach survival thresholds, which LPJ-GUESS
gates on a twenty-year mean and which therefore smooth the cycle away almost
exactly. Then the regional products, biomes, Koppen and lake maps.

Three things about the calendar and the resolution that are easy to get wrong,
each worked through in `biosphere/README.md`:

- **Spin-up is counted in simulation years**, and a simulation year is about
  half an Earth one, so the shipped default buys roughly half the absolute
  time Earth practice assumes. Year counts are scaled UP and annual sums
  scaled DOWN by the reciprocal factor; the two directions are opposite, both
  lists are named in `build_vesper_pfts.py`, and the factor is derived there
  from the configured orbit rather than written anywhere.
- **Every finer candidate can redistribute precipitation rather than merely
  resolve it.** Steeper represented relief means stronger orographic ascent and
  deeper rain shadows, so lee basins can dry and windward coasts wet. That
  reaches the carve verdict, which must be re-taken at the candidate support
  rather than carried over from T42.
- **The 30-hour day widens the real diurnal range and LPJ-GUESS cannot see
  it.** `dtr` reaches only the biogenic VOC scheme, and every cold limit runs
  through a twenty-year mean, so there is no daily-minimum plant mortality to
  trigger.
