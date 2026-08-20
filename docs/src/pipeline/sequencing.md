# 6. What happens next

Three nested loops and then a resolution change.

**A. The terrain loop, one turn of which is an ITERATION.** A generation on
the carve list, then the commissioning of the build it produces: hydrography
and boundary conditions rebuilt on the new terrain, bootstrap, derived fields,
baseline, verdict. Hydrography has to be rebuilt *after* the carve and before
the climate, because the carve changes the drainage the climate is integrated
over.

The block below is the COMMISSIONING, and it is the same block whether the
build is new or is being re-commissioned; a re-commissioning skips the
generation and starts from a build already in `source/`. Each step reads the
one before, and `check_consistency.py` reports which are still missing at any
point:

```bash
python exoplasim/scripts/rebuild_binaries.py --verify # CLAUDE.md rule 4. FIRST, always
python scripts/check_consistency.py                   # rule 8
python scripts/smoke_test.py

python hydrography/scripts/build_hydrography.py       # drainage, basins, coupling
python exoplasim/scripts/build_boundary_conditions.py # land mask, topography
python exoplasim/scripts/build_surface_albedo.py      # lithology albedo, no lakes yet
python exoplasim/scripts/build_surface_roughness.py   # z0
                                                      # soil_water_source: uniform
python exoplasim/scripts/run_exoplasim.py             # BOOTSTRAP run; hours
python exoplasim/scripts/assess_convergence.py <run>  # NOT optional: it writes the
                                                      # run's status, and nothing may
                                                      # treat a run as settled without it
python exoplasim/scripts/build_climatology.py <run>   # then set baseline_climatology
python exoplasim/scripts/index_runs.py                # INDEX.json is the only record

python hydrography/scripts/surface_water.py           # lakes, now a climate exists
python exoplasim/scripts/build_surface_albedo.py --lakes <surface_water.nc>
python pedology/scripts/build_soil.py                 # soil, with no biosphere yet
                                                      # soil_water_source: pedology
python exoplasim/scripts/build_surface_soil_water.py  # 229, from that soil
python aeolian/scripts/build_dust.py                  # needs the climatology AND the lakes
python exoplasim/scripts/dust_optics.py               # only if the spectrum moved
python exoplasim/scripts/dust_aerofile.py             # the model's aerofile
python exoplasim/scripts/build_surface_dust.py        # 1811, only if ndustrad = 1

python exoplasim/scripts/run_exoplasim.py             # the BASELINE; hours
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

The other half is **the carve gate: nothing outstanding may still move an
artifact the verdict is computed from.** It is not a ceremony about
irreversibility -- a wrong verdict is recoverable, because a build is replaced
wholesale -- and it is not a computation. **Read `TASKS.md` and decide.** Each
row in the open set says what it is waiting on.

Whether a row still moves the verdict is a fact about what closing it would
CHANGE, and that lives in the row's prose; a task can name an upstream step
and move nothing -- deferred, a bound, blocked on something that happens
after the carve. The `[step: <id>]` marker records only where the work is
FILED. `TASKS.md`'s conventions record why no marker traversal substitutes
for reading the rows.

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

Four things make that A/B trustworthy, and none is optional:

- **Every forcing change lands with a quantitative prediction of its own
  effect**, stated before it is run, with what result would mean "wrong". That
  is the check-needs-a-right-answer convention (CLAUDE.md), applied to
  physics.
- **Both arms use the SAME BINARY and differ only by a namelist key.** This is
  what the no-op-until-enabled convention is for. Two binaries would confound
  the term with the rebuild, and the low-I/O patch changes the restart layout,
  so arms built either side of it cannot share a restart at all.
- **Both arms branch from ONE restart.** Run-to-run spread on a converged pair
  is larger than several of the terms being tested; a shared initial condition
  turns the comparison into a paired one.
- **The segments are labelled as diagnostics**, so a short A/B tail never
  enters a convergence window or a climatology.

**Bundle the terms, then check the SUM against the sum of predictions.** If
they agree, no attribution was needed and none was bought. If they disagree,
bisect, and only into the subset whose predictions were soft. What bundling
genuinely risks is `docs/src/practice/failure-modes.md` class 15, two errors that nearly
cancel; the mitigation is the per-term prediction above, not serialization,
because two cancelling errors hide just as well in a serial sequence nobody
predicted the size of.

**A budgeting fact.** `build_surface_albedo.py`
has a `modelled` mode that takes tree cover from an LPJ-GUESS `fpc.out`
instead of asserting a uniform vegetated endmember, and it is the mode this
project should end up in. The two land-albedo endmembers are as far apart in
absorbed flux as the largest sweep this project varies on purpose, so adopting
it is a forcing change of the first rank and **it requires the flux to be
re-derived**. That is a cost to budget for, and it is the same cost whether or
not a carve shares the iteration. The corollary is the answer to "when does
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

**B. The soil and biosphere loop, at T42.** For a given climate:

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

**D. The resolution change.** Only after A, B and C have settled.

1. Climate at T85, with the T42 vegetation regridded as its boundary
   condition.
2. Soil and biosphere at T85 on that climatology, loop B again.
3. Climate at T85 again with T85 vegetation, unless step 2's vegetation turns
   out close to the regridded T42 field, which is a cheap comparison worth
   making first.

**The biosphere never decides the resolution.** LPJ-GUESS gridcells are
independent columns, so its cost is linear in land-cell count and trivial at
either resolution; the climate model's cost is what decides. Running it at
T85 on T42 forcing would resolve detail that is not in its input. Expect T85
to lower total NPP, and treat that as a resolution bias rather than a result:
productivity saturates with water, so averaging the forcing before the model
sees it inflates the answer.

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
- **T85 redistributes precipitation, it does not merely resolve it.** Steeper
  relief means stronger orographic ascent and deeper rain shadows, so lee
  basins dry and windward coasts wet. That reaches the carve verdict, which
  should be re-taken at T85 rather than carried over from T42.
- **The 30-hour day widens the real diurnal range and LPJ-GUESS cannot see
  it.** `dtr` reaches only the biogenic VOC scheme, and every cold limit runs
  through a twenty-year mean, so there is no daily-minimum plant mortality to
  trigger.
