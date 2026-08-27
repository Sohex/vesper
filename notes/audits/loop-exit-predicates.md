# The four loop exit predicates, read against what the pipeline does

Vesper is a generated world. This note is about the exit predicates of the four
pipeline loops as `config/pipeline.yaml` states them, and about the gap between
each predicate and the steps and artifacts that would let it be evaluated. Every
subject here is a modelled one: a modelled drainage verdict, a modelled soil, a
modelled vegetation cover, a chain of model runs.

Read 2026-08-27 while building `scripts/verify_joint_convergence.py`, the
finalizer `docs/src/pipeline/loops.md` declares. The finalizer re-evaluates each
loop's own exit against the final state; these are the six things that reading
turned up, each of which is about the predicate rather than about the finalizer.

## 1. The A re-take needs two hours-scale runs, not none

`docs/src/pipeline/loops.md` prices the finalizer: "Verification costs a
fraction -- no new climate run, because the operating support's baseline already
exists, and the two steps it does need, `surface_water` and `dust_forcing` at
that support, are minutes."

That is true of C and of D. It is not true of A. Loop A's exit is the
intersection of the verdicts taken at the TWO bounding climates, and the second
of those is the cold bare-rock arm, which is its own climate rather than its own
analysis of a shared one. In the graph it is `endmember_bootstrap_run` then
`endmember_baseline_run`, each declared `cost: hours`, and neither is nested
under loop D or tied to a rung. Re-taking A's intersection at the operating
support therefore buys a bootstrap and a baseline for the cold arm AT THAT
SUPPORT before any verdict can be taken.

So the finalizer's own cost is not uniform across its four rows. Three of them
read what the final state already carries. The fourth is a commissioning-scale
purchase before it can report anything, which puts loop A in the position of
being the loop whose non-closure is most expensive to detect as well as the one
whose non-closure is most expensive to repair.

## 2. Loop B's exit criteria have no consumer

`pedology/config/pedogenesis.yaml` declares them under `convergence`:
`soil_carbon_relative_tolerance`, `texture_cells_moved_tolerance` and
`maximum_iterations`, with the comment that they were "fixed before any
iteration has been run". Nothing in the tree reads any of the three. The
finalizer is the first thing to evaluate them.

That is not a criticism of the numbers, which are declared in the right place
and in the right order. It is that loop B has had a written exit and no
instrument, so "did the soil and biosphere loop exit" has been a judgement
rather than a measurement for as long as the loop has existed.

## 3. Loop B's criteria cannot be evaluated from what the pipeline keeps

Both tolerances are differences BETWEEN iterations. The `soil` step writes
`pedology/data/{build}/soilmap_{res}.txt` and `pedology/analysis/soil_report.json`,
one of each, and every iteration overwrites both. Nothing retains the previous
iteration, so at the moment the criteria could be applied the quantity they are
about has been destroyed.

`build_soil.py --iteration` is optional and defaults to `None`, so
`soil_report.json` records `"iteration": null` on the run that exists and
`maximum_iterations` has nothing to count against. The finalizer takes the
previous soil as an explicit argument for this reason and refuses when it is not
given, rather than reporting a pass it could not have earned.

The fix is a retained per-iteration soil, not a retained per-iteration report:
`texture_cells_moved_tolerance` is a share of LAND CELLS, so the land means in
the report cannot answer it however many of them are kept.

## 4. The clay-movement threshold is a number with no key

`texture_cells_moved_tolerance: 0.05` is a tolerance on the "fraction of land
cells whose clay content moves by more than 0.02". The 0.05 has a key. The 0.02
exists only in the comment above it. Two numbers decide the criterion and one of
them cannot be read by a script, which is why `verify_joint_convergence.py`
carries it as a module constant citing the comment as its source rather than
reading it.

## 5. Both carve scripts record the CONFIGURED rung, not the climatology's

`export_carve_list.py:1077` takes `resolution` from `config["model"]["resolution"]`
and `:1348` takes `land_surface` from `config["model"]["land_albedo_source"]`.
Both are written into the sidecar as a description of the climate the verdict
was taken on.

The rung guard that would keep those honest is `require_configured_grid`, and it
is reached only through `lib/paths.py:climatology_path()`, which both
`carve_verdict.py:528` and `export_carve_list.py:1074` call ONLY when
`--climatology` was not given. An explicit `--climatology` therefore bypasses
the guard, and the explicit form is exactly the one every re-take uses: the
overshoot route in `hydrography/scripts/carve_overshoot.py` passes it, and so
does any re-take at a support other than the configured one.

The result is a sidecar that can record one rung beside a climatology on
another, with nothing in the file contradicting itself. `verify_joint_convergence.py`
reads the rung off the climatology's own `lat` dimension for this reason and
never off a sidecar's `climate.resolution`.

## 6. Loop C's exit is unreachable at the configured land albedo source

Loop C's exit is "re-take the verdict on modelled vegetation". The mode that
reads modelled vegetation is `build_surface_albedo.py --mode modelled`, which
requires `--vegetation <run>/fpc.out` and raises without it; every other mode,
`vegetated` included, ignores the LPJ-GUESS output entirely. The mode comes from
`model.land_albedo_source`, and `config/planet.yaml` sets it to `vegetated`.

So the graph's `surface_albedo` -> `lpj_run` edge is correct for the mode the
loop is meant to end in and inert at the value the config carries: the pipeline
buys an LPJ-GUESS run whose fpc output the albedo step does not read, and every
baseline run downstream carries the assumed canopy. Loop C cannot exit until
`land_albedo_source` moves to `modelled`, and a verdict taken before it moves is
a verdict on the arm loop A already used rather than a re-take. The finalizer
refuses such a verdict by name.

## What no step generates

Two artifacts the finalizer's predicates need have no generator in
`config/pipeline.yaml`, and by the project's own rule an artifact no step
generates does not exist:

- a carve LIST re-taken at the operating support with both bounding climates,
  which is predicate A's object. `carve_list` writes one, but from the configured
  climatology and into the build's own directory;
- a carve VERDICT taken on a `modelled`-albedo climatology, which is predicate
  C's object. `carve_verdict` and `carve_verdict_endmember` write the two arms A
  uses and there is no third.

Both are `export_carve_list.py` and `carve_verdict.py` invoked with explicit
arguments, so the scripts exist and the rows do not. Until they land the
finalizer reports A and C as NOT EVALUABLE and names the command in each case.
