# The carve overshoot: what the intersection bracket does not cover

Vesper is a generated world. This note is about loop A, the terrain loop, and
about a residual in its exit predicate. Everything here concerns a modelled
drainage verdict taken on a modelled climate of an invented planet.

Written 2026-08-24 against the loop as `config/pipeline.yaml` and
`docs/src/pipeline/loops.md` state it, and against
`hydrography/scripts/carve_verdict.py` and `export_carve_list.py` as they stand.
No overshoot has been measured, and the last section says exactly what is
missing.

## The two uncertainties, which are not the same one

Loop A's exit predicate: "Successive verdicts bracket rather than converge,
because the map is antitone: a larger carve set produces a smaller next verdict.
Exit by carving the INTERSECTION of the verdicts taken at the two bounding
climates, which makes the BRACKETED set a construction rather than a tolerance."

Two distinct things are packed into that sentence.

**The antitone oscillation.** Carving removes closed-basin fill, the brightest
lithology on this world, so the modelled land darkens, the modelled world warms,
open-water evaporation rises, and marginal basins would now stay closed. The
verdict map is order-reversing in the carve set, so successive verdicts
alternate around a fixed point rather than approaching it from one side. This is
an uncertainty about the CARVE'S OWN FEEDBACK.

**The endmember bracket.** The verdict is taken before LPJ-GUESS has ever run,
so it inherits an assumed biosphere, and the bare-rock and vegetated surfaces
reach the design mean at non-overlapping stellar fluxes. Taking the verdict at
both bounding surfaces and carving only where both cut is an uncertainty about
the VEGETATION STATE.

The exit predicate names the first and then exits on the second. The bracket is
a construction, and it is a construction over the vegetation axis; it is not a
bound on the antitone feedback, and it cannot be made into one by widening it.
Both arms run on the same PRE-CARVE terrain -- `endmember_bootstrap_run` differs
from the vegetated baseline in `model.land_albedo_source` and in nothing else --
so both are cold relative to the world their own carve produces, and both cut
more than that world would.

The set algebra makes this exact. An overshooting basin is one the applied
verdict cut; the applied verdict is the intersection, so both arms cut it; the
bracketed set is where the two arms DISAGREE; so an overshooting basin is never
in the bracketed set. Not rarely. Never. The bracket width therefore carries no
information about the overshoot, in either direction, and reporting the width as
though it were the honest uncertainty on the carve understates it by whatever
the overshoot turns out to be.

That is not an argument for abandoning the intersection. Carving the intersection
is still the smallest defensible set, and a smaller carve has a smaller feedback,
so the intersection has the smallest overshoot of the three candidate sets. It is
an argument that the exit is one measurement short, and `docs/src/pipeline/loops.md`
already says which measurement: "On iteration 2, the already-carved set should be
re-evaluated against the new climate and the number that would no longer have
carved reported."

## The measurement has a route, and the obvious one is not it

Re-running the verdict on the carved build does not answer the question. A
carved basin's rim has been breached, so the finished terrain holds no
impoundment whose catchment-to-spill-area ratio could be measured, and the basin
is absent from that build's `basins.nc` altogether. `export_carve_list.py`
already says so where it carries such an entry forward at retain 0: its verdict
is not re-decidable there. Running the verdict on the carved build measures
which of the SURVIVING basins would now carve, which is the next pass, not the
overshoot.

The route that works holds the geometry still and moves only the climate: take
the verdict again on the PRE-CARVE build's own `basins.nc` and coupling matrix,
under the CARVED build's baseline climatology. Every term except the climate is
then the one the applied verdict used, which is what isolates the feedback.
`carve_verdict.py` and `export_carve_list.py` already take `--basins`,
`--coupling` and `--climatology` separately, so no new machinery is needed to
produce the re-evaluation; `hydrography/scripts/carve_overshoot.py` is the
comparison, and it refuses a re-evaluation that does not cover every applied
carve rather than scoring an unknown as a zero.

Three properties of that route are worth writing down, because each is a way to
take the number and have it be wrong.

**It is a deliberate cross-build read.** CLAUDE.md rule 5: pointing one
component at another's output is deliberate, never defaulted. This is the
deliberate case, and it must be declared wherever it is recorded, because every
other pairing of a climatology with a basin set on this project has been a
defect.

**The background land albedo goes through one door, and the cross-build read is
declared.** `exoplasim/inputs/<rung>/orogen_<RUNG>_surf_0174.sra` is keyed by
rung alone while `surface_albedo` regenerates it per build, so it holds exactly
one build's albedo at a time and the path cannot say which. Every consumer --
`carve_verdict.py`, `export_carve_list.py`, `surface_water.py`,
`pedology/scripts/build_surface_classes.py` and
`exoplasim/scripts/dust_forcing.py` -- now resolves it through
`lib/provenance.py:staged_surface_field`, which refuses a field from another
build at the read and returns the record each of them stamps into its own
product. So a verdict says which build's albedo it used.

For the overshoot measurement the staged file is the right one, because the
carved build is what has just been commissioned; for anything re-run afterwards
on the pre-carve build it is the wrong one. That read is DELIBERATE and the two
carve scripts declare it with `--for-build`, which NAMES the build and so admits
that one and refuses every other. A flag that merely turned the check off would
have let the same silence back in.

**Both arms, or neither.** If the applied verdict was an intersection, the
re-evaluation has to be an intersection too, taken against the carved build's
own bare-rock arm. Comparing an intersection against one arm measures the
bracket and the feedback added together and attributes the sum to the feedback.
The carved build produces both arms in the course of being commissioned, so this
costs no extra climate run.

## What the number will mean when it exists

It is one-signed. The antitone argument predicts overshoot and predicts no
undershoot, so `carve_overshoot.py` reports both and a large undershoot is
evidence that something other than the carve feedback moved between the two
climates. Read it as a check on the comparison, not as a second result.

It is reported beside the bracket width in the same currency, basins, because
that is the only way to say whether the unmeasured uncertainty is small next to
the measured one. A residual of a few basins against a bracket of hundreds says
the exit predicate is sound in practice; a residual of the same order says the
loop was exited a pass early and the carved set is an upper bound that was
published as an answer.

The one prior measurement of this feedback's INPUT is on the archived
`carved-zoned` build: applying its verdict dropped closed-basin fill from 20.9%
to 12.7% of land, and that build's own note records the overshoot as not
measured. Its lineage is superseded twice over -- the verdict was computed on
antipodal climate and the lithology under it has since changed -- so the figure
is a magnitude for what a carve does to the surface, not a result to carry
forward.

## Blocked, and on what exactly

The overshoot cannot be measured now, and the blocker is not "a pre-carve base
restarts the count at iteration 1". It is more basic than that: nothing in this
lineage has taken a first verdict at all.

- `world_state.json` has `hydrography.carve_verdict` null in all three fields
  and `carve_list` null. No verdict exists for `precarve-craton`.
- `current_climate` is invalidated, with no `baseline_climate_report.json` for
  any climatology this build produced, so the input the first verdict needs is
  not there either.
- No carved build in this lineage has a payload. `archive/builds/` holds five
  `carved-zoned*` identities and every one is a stub.

So the measurement waits on a first verdict, a carve, and the carved build's
own baseline. All three are loop A actions and the carve is the user's call.
What is not blocked, and is done, is the instrument and the route: the
comparison exists with its controls, the route that isolates the feedback is
written down, and the bracket's per-basin membership is now recorded in the
carve list sidecar so a later pass can ask which basins moved rather than only
how many.
