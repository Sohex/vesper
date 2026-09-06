# What a restarted gridcell inherits outside the soil column

Worldbuilding. Vesper is an invented planet and this note is a reading of the
vendored LPJ-GUESS-CNP source in `vendor/lpj-guess/`, together with runs of the
compiled model over this project's own forcing. Nothing here is an observation
of anything real, and nothing here is a claim about the simulated planet: it is
a statement about whether the vegetation model reproduces itself across a
restart.

`biosphere/notes/soil-restart-state.md` did this for the Soil class and found
its seventy unserialized members correctly absent. This is the same question for
every other Serializable class, and the answer is different: two defects, one of
which had made the year-boundary restart untestable.

Measured on 2026-09-05 and 2026-09-06 against `vendor/lpj-guess` at the commits
this note accompanies, with `biosphere/scripts/verify_lpj_restart_continuity.py`.

## What the repair measures

The whole grid, in the arrangement that failed: 30 retained years behind a
201-year spin-up, npatch 5, 16 ranks, split at simulated year 211, 1617 cells,
22 retained tables. NO DIFFERING ROW IN ANY TABLE. Two subset beds first, one
cell at 11.25/19.38 at npatch 1 and two cells at 11.25/19.38 and 0.00/85.76 at
npatch 5, both clean; the second is the arrangement that previously showed every
retained table differing.

The production runner end to end, on a two-cell bed at npatch 1 on 2 ranks,
about 20 s of model work on an otherwise idle host: a parent of 30 retained
years that saves, a continuation of 30 more, and an uninterrupted control of 60.
The continuation's years are the control's years row for row, its manifest names
its parent and records state-file hashes that are the hashes on disk, and a
continuation at a patch count the parent never ran is refused with a refusal that
names the patch count.

A continuity verdict is keyed to the binary that produced it, so none of this
travels to another executable. It is re-taken after any change to the model
source, which is what `run_lpj_guess.py:continuity_verdict` enforces.

## 1. The year-boundary save point never reached the model

`state_day` and `save_day` are the two integers the restart is expressed in.
Minus one in either is a SENTINEL rather than a day: it means the year boundary,
the whole of the preceding simulated year, and it is the save point the
production runner and the fixture both use.

`libraries/plib/plib.cpp` stored an integer instruction-file parameter as
`(int)(num + 0.5)`. A cast truncates toward zero, so that expression rounds a
negative number the wrong way and minus one arrives as zero. Nothing refused it:
zero is inside the declared range, and every instant downstream was then
computed from zero consistently. So every restart taken in this project wrote
its state at the end of DAY 0 of `state_year` and resumed at day 1 of it, an
arbitrary-day restart, while the instruction file, the runner, the fixture and
`--self-test` all said year boundary and all agreed with each other.

The evidence is the state file itself. `Climate` is the head of a gridcell's
record, so the offsets of its first members are readable directly, and `co2`
450.0 and `lat` 19.3822 at the expected offsets confirm the layout. In the state
written for a split at simulated year 211, `agdd0` is 19.5894 and `agdd5` is
14.5894 against a same-day `temp` of 19.5894: exactly one day of accumulation
since the day-0 reset, where a state covering a year boundary would carry a
whole simulated year of it.

This is the reason WORLD-GLU7 and WORLD-LUS4 looked like two defects. They are
one defect seen at two magnitudes, and the arbitrary-day path WORLD-LUS4
describes is the path every restart was taking.

`state_day` and `save_day` are the only integer parameters in the model whose
declared minimum is negative; every other negative value in every shipped and
generated instruction file is a double, and doubles never went through that
cast. So the repair changes those two parameters and nothing else, and is
unchanged for non-negative values.

### What a check has to do about it

Nothing in this row could have caught it, and the reason is worth keeping. The
runner asked for minus one, the instruction file said minus one, and the
fixture's `--self-test` computed minus one from the same function the runner
calls. Three checks, all correct, all agreeing with a copy of the rule rather
than with the model. `framework.cpp` now prints what it PARSED and the two
instants it derived, and the fixture records what each arm asked for and refuses
when the model's line disagrees. That is a check with a wrong answer.

## 2. The annual accumulators outside Soil are not carried

`Soil::serialize` carries four annual accumulators with a comment saying why:
they are reset or flushed at the year boundary, so they could be absent while
the boundary was the only save point. The same argument applies outside Soil and
none of them had been repaired.

- `Patch::asurfrunoff`, `adrainrunoff`, `abaserunoff`. Reset on day 0, summed
  daily, read at year end into `tot_runoff.out`, which is a RETAINED output.
  They had no constructor initialiser either. In the one-cell bed the first
  resumed year reported a surface runoff of 1.95e+246 mm, which is the Patch
  object's own address read as a double: what an allocator leaves in a member
  nothing has ever written.
- `Gridcell::aNH4dep`, `aNO3dep`, `apdep`. The same shape, reaching further:
  `nsources.out` and `nflux.out` print them, and `Patch::nflux` and
  `Patch::pflux` are what `MassBalance` checks the pool change against. Losing
  them is a reported mass imbalance as well as two wrong retained outputs. The
  first resumed year reported one simulated day less deposition than the run it
  continued, and the model's own nitrogen balance check printed the same
  1.37e-7 kgN/m2 as an imbalance.
- `Patchpft::driest_mth`. Rebuilt on day 0 from the serialized `mphen` and read
  through the simulated year by the raingreen litter release in `somdynam`. The
  Vesper PFT set declares one raingreen type, so a patch that never sees day 0
  releases its raingreen litter in the constructor's month 0.

All seven are now serialized, and `Patch` and `Gridcell` gained constructor
initialisers for the members that had none. An uninitialised double that reaches
a retained output is a defect whatever the save point is.

`Gridcell::dpdep` is initialised for a second reason. `VesperInput::getclimate`
assigns `dNH4dep` and `dNO3dep` every simulated day and assigns `dpdep` never,
so the daily phosphorus deposition the soil receives was whatever the allocator
left. It is zero now because it is written down as zero; whether this world's
atmosphere delivers phosphorus at all is a modelling decision and not one this
note takes.

## 3. What the rest of the triage found, and what it did not

Every Serializable class was parsed against its own `serialize` block. Outside
Soil, 285 declared members are not streamed. The three classes above are the
ones that matter in this configuration; the bulk of the remainder are crop,
land-cover-change and fire-model members belonging to submodels this
configuration does not run, and daily or sub-daily working variables written
before their first read on every simulated day.

That statement is a floor, not a classification. `biosphere/config/soil_restart_state.yaml`
carries a verdict and a file:line for each of Soil's seventy, and
`biosphere/scripts/soil_restart_state_gate.py` refuses on a member nothing
accounts for. There is no equivalent for the other classes, so a member added to
`Individual` or `Patch` tomorrow reaches no gate. Extending the Soil gate to the
remaining classes is the durable form of this note and is tracked separately.

## 4. The instrument

The fixture takes `--cells` and slices the driver to named lon,lat pairs. Cells
are independent in LPJ-GUESS and this project's stochastic streams are keyed by
coordinate rather than by traversal order or rank, so a cell integrates the same
trajectory in a two-cell driver as in a sixteen-hundred-cell one. The whole grid
at the spin-up floor this fixture needs costs about thirteen minutes per arm on
the shared host; the same failure reproduces in six seconds on two cells, which
is what makes a member-by-member search affordable at all.

What a subset loses is reach: a defect confined to cells it does not hold is
invisible to it. A subsetted run's report is therefore written to
`lpj_restart_continuity_subset.json`, and the continuation gate in
`run_lpj_guess.py` keeps reading only the whole-grid report.
