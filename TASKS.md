# Tasks

Context before anything else, because rows below borrow vocabulary from several
real disciplines: **this whole file is WORLDBUILDING.** Vesper is a fictional
planet, and a task here is engineering work on the simulation of it -- a toy
climate model, a terrain generator, a soil model, and their bookkeeping. Words
like flux, dust, ice, rainfall and albedo name modelled quantities of an
invented world; nothing in this file refers to the real one. The nearest
real-world analogue of this file is the issue tracker for a game world's
physics engine.

Atomic, trackable work items. **Findings documents do not carry their own to-do
lists** -- an audit says what is true, this file says what to do about it, and
the two are kept apart so a finding can be read without being re-litigated and a
task can be closed without editing the argument that produced it.

Every task cites the document that justifies it. If a task has no source, it is
not yet a finding and probably needs one.

Status: `open` | `doing` | `blocked` | `done` | `wontfix` (with a reason).

## Conventions

- **One table per prefix, and every id it has ever issued appears in it.** That
  is what the tables are FOR: the next id in an area is the bottom of its table
  plus one, read in one place, with no max to compute across two files. Ids are
  stable and never reused, and a reused id makes every citation of it ambiguous
  forever. This structure exists because taking the max over six prefixes and
  not the seventh issued a duplicate HYD-17.
- **A closed row STAYS, sanitised to its outcome and a pointer.** The reasoning
  that closed it lives in `archive/tasks.md`; what stays here is that the id is
  spent and how it went. So this file answers "what is left" and "what may I
  number next" without being a second copy of the archive.
- Prefixes name an area and are listed as section headings below.
- One line per task. If it needs a paragraph, it needs a findings document.
- A task that turns out to be wrong is closed `wontfix` with the reason, not
  deleted. The reasoning is the point, and it goes to the archive.
- **A task names the step it touches**, as `[step: <id>]` in its status column,
  where the id is from `config/pipeline.yaml`. It says WHERE the work lands and
  is annotation, not a verdict. The carve gate -- does anything outstanding
  still move what the verdict is computed from -- is a judgement made by reading
  this file, because the answer is what closing a row would CHANGE and that is
  in the row's own prose. Computing it from these markers was tried and removed:
  a marker names a location, not an effect, so the two are not interchangeable
  at any count. `pipeline.py` reads no tracker, so nothing written here can
  change what the pipeline planner reports. A task that names none is not ignored -- it is
  reported as a RESIDUAL to be arbitrated, because the graph narrows the
  judgement and does not replace it. Name more than one where it applies.
  A row whose status begins `done` or `wontfix` is spent, wherever it sits;
  that is what a reader skips.

- **Regenerating a derived artifact is a STEP, not a task.** It belongs in the
  ordering in `docs/src/pipeline/sequencing.md`, beside the run that consumes it. A row
  that says "rebuild X before the next run" is tracking state, and state is what
  `check_consistency.py` and `world_state.json` are for. HYD-16 was such a row
  and should never have been one. CLAUDE.md rule 7.
- **A task is code, physics, a decision or a measurement that settles a
  mechanism.** If closing it would produce only a number that the next iteration
  regenerates, it is not a task; the durable half is whatever it teaches.

---

## BIO -- biosphere

2 open of 5 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| BIO-1 | -- | -- | done, see `archive/tasks.md` |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; the span is 18.2% of NPP and it is the largest nitrogen lever | `biosphere/notes/productivity-prediction.md` | blocked on iteration 2's baseline run -- no LPJ-GUESS run exists on this build [step: lpj_run] |
| BIO-3 | -- | -- | done, see `archive/tasks.md` |
| BIO-4 | -- | -- | done, see `archive/tasks.md` |
| BIO-5 | Extend `vesperinput.cpp` for the CNP fork's `SoilProperties` fields (`kplab`, `spmax`, `pwtr`): the subclass must carry them even if it does nothing with them | `biosphere/notes/cnp-fork-scoping.md` | open, blocked on the CNP fork build itself [step: lpj_run] |

## BUDG -- the error budget and what it is denominated in

0 open of 6 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| BUDG-1 | -- | -- | done, see `archive/tasks.md` |
| BUDG-2 | -- | -- | done, see `archive/tasks.md` |
| BUDG-3 | -- | -- | done, see `archive/tasks.md` |
| BUDG-4 | -- | -- | done, see `archive/tasks.md` |
| BUDG-5 | -- | -- | done, see `archive/tasks.md` |
| BUDG-6 | -- | -- | done, see `archive/tasks.md` |

## CLIM -- climate

8 open of 36 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| CLIM-1 | -- | -- | done, see `archive/tasks.md` |
| CLIM-2 | -- | -- | done, see `archive/tasks.md` |
| CLIM-3 | -- | -- | done, see `archive/tasks.md` |
| CLIM-4 | -- | -- | wontfix, see `archive/tasks.md` |
| CLIM-5 | -- | -- | done, see `archive/tasks.md` |
| CLIM-6 | -- | -- | done, see `archive/tasks.md` |
| CLIM-7 | -- | -- | done, see `archive/tasks.md` |
| CLIM-8 | -- | -- | done, see `archive/tasks.md` |
| CLIM-9 | -- | -- | done, see `archive/tasks.md` |
| CLIM-10 | -- | -- | done, see `archive/tasks.md` |
| CLIM-11 | -- | -- | done, see `archive/tasks.md` |
| CLIM-12 | -- | -- | done, see `archive/tasks.md` |
| CLIM-13 | -- | -- | done, see `archive/tasks.md` |
| CLIM-14 | -- | -- | done, see `archive/tasks.md` |
| CLIM-15 | -- | -- | done, see `archive/tasks.md` |
| CLIM-16 | Ocean horizontal heat diffusion has a namelist switch, `nhdiff` with `hdiffk`, and the error budget records no-q-flux as STRUCTURAL with "no cheap version exists" | `notes/audits/absent-and-inherited-physics.md` finding 1, `analysis/error_budget.json` other_items | open, and it is the cheapest structural term in the budget to price rather than the impossible one. `oceanmod.f90` reads both keys in `oceanmod_namelist`, broadcasts them and acts on `nhdiff > 0` at line 844, so it is one key on the binary already built -- which satisfies every one of `docs/src/pipeline/sequencing.md` A3's four A/B conditions with no further work. Bracket `hdiffk` and report the spread rather than tuning it, the way this project handles every other ungrounded constant. Read it as physics-is-not-a-knob: ocean heat transport EXISTS, the argument for enabling it is not that it improves an agreement, and a constant diffusivity gives a BOUND on the missing transport rather than the transport. The error budget entry is CORRECTED 2026-08-18: its second sentence, that no cheap version exists, was false and the row now records that `nhdiff` is one namelist key on the built binary. What is left is the measurement -- bracket `hdiffk`, run the paired arms off one restart, report the spread. The PREDICTION is written, 2026-08-19, with its wrongness bounds: `exoplasim/notes/forcing-bundle-predictions.md`, computed by `exoplasim/scripts/predict_ocean_terms.py`, whose operator is validated against a Laplacian eigenfunction and whose applied flux the run itself now records in the ocean stream (CLIM-12), so the A/B has a pointwise right answer [step: baseline_run] |
| CLIM-17 | -- | -- | done, see `archive/tasks.md` |
| CLIM-18 | -- | -- | done, see `archive/tasks.md` |
| CLIM-19 | -- | -- | done, see `archive/tasks.md` |
| CLIM-20 | -- | -- | done, see `archive/tasks.md` |
| CLIM-21 | -- | -- | done, see `archive/tasks.md` |
| CLIM-22 | -- | -- | done, see `archive/tasks.md` |
| CLIM-23 | -- | -- | done, see `archive/tasks.md` |
| CLIM-24 | -- | -- | done, see `archive/tasks.md` |
| CLIM-25 | -- | -- | done, see `archive/tasks.md` |
| CLIM-26 | -- | -- | done, see `archive/tasks.md` |
| CLIM-27 | -- | -- | done, see `archive/tasks.md` |
| CLIM-28 | -- | -- | done, see `archive/tasks.md` |
| CLIM-29 | Carbonaceous aerosol is the one the sulfate bound does not reach: smoke and secondary organics ABSORB, so their forcing per unit optical depth is larger than either scatterer's and of the opposite sign. Fire is enabled in LPJ-GUESS as GLOBFIRM and biogenic emissions are already in the driver, so both are estimable from a biosphere run | `notes/audits/unpriced-terms.md` finding 2, `aeolian/README.md` | blocked on a biosphere run existing on this build, which is BIO-1 and BIO-2, and deliberately not worked around: the burned area and the biogenic emission are LPJ-GUESS output and inventing either would be inventing the vegetation the whole component exists to compute. The chain to reuse is complete -- OPAC carries a soot component beside the sulfate one, `sea_salt_optics.band_average` takes any distribution, and the transport and the two-stream are shared -- so this is a source term and an optics table and nothing else. Note the sign: an absorbing aerosol over this world's bright closed-basin fill warms, which is the same asymmetry `dust_forcing.py` already prices for mineral dust [step: lpj_run, sea_salt_optics] |
| CLIM-30 | DECLARE the extreme-cold land-fraction cap BEFORE the design flux is re-derived, instead of inferring it from the answer. The purged artifact inferred 0.0642, which reproduced 0.945 from an unconstrained winner of 0.8725 dry -- a threshold fitted to the answer it was meant to test, which is what the threshold convention (`docs/src/practice/conventions.md`) forbids | `docs/src/pipeline/state.md` section 5b, `exoplasim/scripts/derive_design_flux.py` | open, and a DECISION nothing in the project settles: it needs a preference no document fixes, which is why it is a task and the re-derivation is not. Re-deriving is the `design_flux` step, and the second converged point it needs is `bracket_run`; neither is tracked here. Until the cap is declared ahead of the run, the recorded flux is unsupported [step: design_flux] |
| CLIM-31 | -- | -- | done, see `archive/tasks.md` |
| CLIM-34 | Redo the glacier bound against the real 0.945 climatology with the measured 7.8 K/km lapse; the recorded -4.7 K offset is a proxy | `notes/glacier-rough-pass.md`, archive PHYS-12 | open [step: baseline_climatology] |
| CLIM-35 | Measure the salinity A/B: `TFREEZE` from declared salinity is settable (archive CLIM-17) and its global-mean worth is unmeasured | `scripts/error_budget.py` salinity row, `archive/tasks.md` CLIM-17 | open, one namelist key on the built binary, so it satisfies A3's A/B conditions [step: baseline_run] |
| CLIM-36 | Reconcile forest fraction between codes 173 and 212: `build_surface_roughness.py` reads `model.forest_fraction_assumed` (absent, so zero forest) while `build_surface_albedo.py` writes 212 at 0.5 in `vegetated` mode, so the two fields do not describe one land cover | `exoplasim/scripts/build_surface_roughness.py` note at the surface term | open, a DECISION on which source both read [step: surface_roughness, surface_albedo] |
| CLIM-32 | Price ozone's radiative effect on the simulated climate: one T21 sensitivity pair at `o3scale` 0.5 against the baseline | `exoplasim/notes/ozone.md`, `docs/src/reference/config-rationale.md` activity entry | open. Shielding of the surface is settled (0.794, measured); the radiative stake is not, and in a ten-layer model it may be small. A diagnostic pair, not a commissioning [step: baseline_run] |
| CLIM-33 | Price `mixed_layer_depth_m`: one perturbation run against the 50 m default, which sets seasonal amplitude on a half-Earth year | `docs/src/reference/config-rationale.md` ocean entry, `analysis/error_budget.json` structural items | open. The budget books it structural with no kelvin figure; a single perturbation run converts it [step: baseline_run] |

## CONS -- consistency checking

1 open of 9 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| CONS-1 | -- | -- | done, see `archive/tasks.md` |
| CONS-2 | -- | -- | done, see `archive/tasks.md` |
| CONS-3 | -- | -- | done, see `archive/tasks.md` |
| CONS-4 | -- | -- | done, see `archive/tasks.md` |
| CONS-6 | -- | -- | done, see `archive/tasks.md` |
| CONS-7 | -- | -- | done, see `archive/tasks.md` |
| CONS-8 | -- | -- | done, see `archive/tasks.md` |
| CONS-9 | Decide whether the producing binary joins `run_id`'s identity: `binary_manifest.json` guards the build, but a run's identity does not record which binary made it | `docs/src/practice/failure-modes.md` class 11 | open, and a DECISION: the UUID rule (CLAUDE.md rule 6) says identity encodes nothing, so the binary would join the manifest rather than the name [step: baseline_run] |
| CONS-9 | -- | -- | done, see `archive/tasks.md` |

## CONV -- cross-component conventions and provenance plumbing

0 open of 4 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| CONV-1 | -- | -- | done, see `archive/tasks.md` |
| CONV-2 | -- | -- | done, see `archive/tasks.md` |
| CONV-3 | -- | -- | done, see `archive/tasks.md` |
| CONV-4 | -- | -- | done, see `archive/tasks.md` |

## CYC -- the stellar cycle

1 open of 1 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| CYC-1 | Measure the cycle's DAMPING FACTOR, known only as the range 0.24 to 0.6. That one number converts any future flux amplitude into a climate without another run, so it is worth more than the run that measures it, and it survives every regeneration | `docs/src/pipeline/sequencing.md` A2, `config/planet.yaml` stellar_cycle | open. RUNNING the cycle is not this row and never should have been: `stellar_cycle_run` is a step and `carve_verdict` NEEDS it, so the graph enforces the ordering A2 argues for -- which is what this row was created to track when nothing did. What is left here is the measurement. Length is set in periods of the LONG component, 57 Earth years, not the medium one [step: stellar_cycle_run] |

## DUST -- the aeolian component

4 open of 15 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| DUST-1 | -- | -- | done, see `archive/tasks.md` |
| DUST-2 | -- | -- | done, see `archive/tasks.md` |
| DUST-3 | -- | -- | done, see `archive/tasks.md` |
| DUST-4 | -- | -- | wontfix, see `archive/tasks.md` |
| DUST-5 | -- | -- | done, see `archive/tasks.md` |
| DUST-6 | -- | -- | done, see `archive/tasks.md` |
| DUST-7 | -- | -- | done, see `archive/tasks.md` |
| DUST-8 | -- | -- | done, see `archive/tasks.md` |
| DUST-9 | -- | -- | done, see `archive/tasks.md` |
| DUST-10 | -- | -- | done, see `archive/tasks.md` |
| DUST-11 | Settle the CATCHMENT half: run one prescribed-dust climate segment and measure what it does to precipitation and runoff. Runoff is the carve criterion's denominator and only 19% of land precipitation, so it amplifies; radiatively dark mineral dust lowers the simulation's rainfall through a fast circulation adjustment that no offline calculation can reach, because it is a response of the modelled winds and not of a surface energy balance | `notes/dust.md`, `aeolian/notes/prescribed-dust-run.md` | DECIDED 2026-08-17: measure first, then decide whether to carry it. The matched pair rides inside the flux re-bracket rather than beside it, so the control is a point loop A needs anyway. The patch is resident with `ndustrad = 0` [step: surface_dust, baseline_run] |
| DUST-12 | -- | -- | done, see `archive/tasks.md` |
| DUST-13 | Decide whether the final climate carries interactive emission (`L_AERO = 1`): the reopening test in `notes/dust.md` fired, so a prescribed field is not defensible for a converged answer | `notes/dust.md`, `aeolian/notes/in-model-dust.md` | open, blocked on DUST-11's measurement -- the prescribed-dust run prices what the interactive scheme must reproduce [step: baseline_run] |
| DUST-15 | Make the measured gust samples the dust step's default input and regenerate the chain: every dust artifact on disk is 70x too thin, built from the snapshot-fitted wind tail | `notes/dust.md` last section, `config/pipeline.yaml` dust step | open, and DUST-11's run rides the wrong field until it lands [step: dust] |
| DUST-14 | Dust deposition reaches the soil and never reaches the cryosphere: nothing consumes the deposition field as a term in the model's snow and ice albedo, and `build_dust.py` reads snow only as an emission suppressor | `notes/audits/absent-and-inherited-physics.md` finding 2 | open, one-signed, and it lands on the term the glacier result turns on. Deposition over snow-covered land at 50-60 degrees is 8.75 g/m2/yr at the central aeolian roughness against a terrestrial dust-on-snow literature working at 1-5 g/m2 snowpack loads for albedo reductions of 0.03-0.08, and the smooth end of the roughness bracket is 111 g/m2/yr. The model's snow albedo carries a time-since-snowfall aging range of about a quarter in band 1 and no dependence on what has landed on it. It belongs in the glacier mass balance when `notes/glacier-rough-pass.md` stops being a temperature criterion, and the field it needs already exists. It also puts a cryosphere term under the aeolian roughness bracket, which was understood as controlling emission and the direct forcing only [step: dust, surface_albedo] |

## GRAV -- gravity

1 open of 7 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| GRAV-0 | -- | -- | wontfix, see `archive/tasks.md` |
| GRAV-1 | -- | -- | wontfix, see `archive/tasks.md` |
| GRAV-2 | -- | -- | wontfix, see `archive/tasks.md` |
| GRAV-3 | -- | -- | wontfix, see `archive/tasks.md` |
| GRAV-4 | -- | -- | done, see `archive/tasks.md` |
| GRAV-5 | -- | -- | done, see `archive/tasks.md` |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here) | `notes/audits/orogen-gravity.md` | blocked on the downscaling pass, deliberately. The audit's recommendation was Option B, carry it as a declared gap: `iceFlow` is a heuristic accumulation rather than a Glen's-law velocity, so a rigorous g^3 on it is precision theatre. Decisive on top of that, a glacial valley is 1-5 km against a 15.19 km mesh cell, so the process is SUB-GRID and cannot be resolved here at any gravity. It belongs with the downscaling machinery, which has to persist sub-grid hypsometry anyway [step: orogen] |

## GRID -- the Orogen/ExoPlaSim grid convention

0 open of 1 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| GRID-1 | -- | -- | done, see `archive/tasks.md` |

## HYD -- hydrography

2 open of 18 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| HYD-1 | -- | -- | done, see `archive/tasks.md` |
| HYD-2 | -- | -- | done, see `archive/tasks.md` |
| HYD-3 | Measure the overshoot: re-evaluate an already-carved set against the climate that carving produced, and report how many would no longer have carved | `docs/src/pipeline/loops.md` | blocked -- it applies to a verdict taken on carved terrain, and a pre-carve base restarts the count at iteration 1 [step: carve_verdict] |
| HYD-4 | -- | -- | done, see `archive/tasks.md` |
| HYD-5 | -- | -- | done, see `archive/tasks.md` |
| HYD-6 | -- | -- | done, see `archive/tasks.md` |
| HYD-7 | -- | -- | done, see `archive/tasks.md` |
| HYD-8 | -- | -- | done, see `archive/tasks.md` |
| HYD-10b | -- | -- | done, see `archive/tasks.md` |
| HYD-11 | -- | -- | done, see `archive/tasks.md` |
| HYD-12 | -- | -- | done, see `archive/tasks.md` |
| HYD-13 | -- | -- | done, see `archive/tasks.md` |
| HYD-14 | -- | -- | wontfix, see `archive/tasks.md` |
| HYD-15 | -- | -- | done, see `archive/tasks.md` |
| HYD-16 | -- | -- | done, see `archive/tasks.md` |
| HYD-17 | -- | -- | done, see `archive/tasks.md` |
| HYD-18 | -- | -- | done, see `archive/tasks.md` |
| HYD-19 | Add a clause to `vendor/orogen/tools/README.md`: erodibility is mean-normalised at generation time, not as shipped, so a consumer must not rescale against an assumed mean of 1 | `hydrography/notes/orogen-carving-request.md` | open [step: orogen] |

## LITH -- lithology

0 open of 25 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| LITH-0 | -- | -- | wontfix, see `archive/tasks.md` |
| LITH-1 | -- | -- | wontfix, see `archive/tasks.md` |
| LITH-2 | -- | -- | done, see `archive/tasks.md` |
| LITH-3 | -- | -- | done, see `archive/tasks.md` |
| LITH-4 | -- | -- | done, see `archive/tasks.md` |
| LITH-5 | -- | -- | done, see `archive/tasks.md` |
| LITH-6 | -- | -- | done, see `archive/tasks.md` |
| LITH-7 | -- | -- | done, see `archive/tasks.md` |
| LITH-8 | -- | -- | done, see `archive/tasks.md` |
| LITH-9 | -- | -- | done, see `archive/tasks.md` |
| LITH-10 | -- | -- | done, see `archive/tasks.md` |
| LITH-11 | -- | -- | done, see `archive/tasks.md` |
| LITH-12 | -- | -- | done, see `archive/tasks.md` |
| LITH-13 | -- | -- | done, see `archive/tasks.md` |
| LITH-14 | -- | -- | done, see `archive/tasks.md` |
| LITH-15 | -- | -- | done, see `archive/tasks.md` |
| LITH-16 | -- | -- | done, see `archive/tasks.md` |
| LITH-17 | -- | -- | done, see `archive/tasks.md` |
| LITH-18 | -- | -- | done, see `archive/tasks.md` |
| LITH-19 | -- | -- | done, see `archive/tasks.md` |
| LITH-20 | -- | -- | done, see `archive/tasks.md` |
| LITH-21 | -- | -- | wontfix, see `archive/tasks.md` |
| LITH-22 | -- | -- | wontfix, see `archive/tasks.md` |
| LITH-23 | -- | -- | done, see `archive/tasks.md` |
| LITH-24 | -- | -- | done, see `archive/tasks.md` |

## MIN -- economic minerals

0 open of 5 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| MIN-1 | -- | -- | done, see `archive/tasks.md` |
| MIN-2 | -- | -- | done, see `archive/tasks.md` |
| MIN-3 | -- | -- | done, see `archive/tasks.md` |
| MIN-4 | -- | -- | done, see `archive/tasks.md` |
| MIN-5 | -- | -- | done, see `archive/tasks.md` |

## PHYS -- physics calibrated for the wrong world

1 open of 10 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| PHYS-1 | -- | -- | done, see `archive/tasks.md` |
| PHYS-2 | -- | -- | done, see `archive/tasks.md` |
| PHYS-3 | -- | -- | done, see `archive/tasks.md` |
| PHYS-4 | -- | -- | done, see `archive/tasks.md` |
| PHYS-5 | -- | -- | done, see `archive/tasks.md` |
| PHYS-6 | -- | -- | done, see `archive/tasks.md` |
| PHYS-9 | -- | -- | done, see `archive/tasks.md` |
| PHYS-10 | -- | -- | done, see `archive/tasks.md` |
| PHYS-11 | The cloud shortwave and longwave constants are Earth tunings and are the only major radiation term never re-weighted for this star. `tswr1`, `tswr2`, `tswr3`, `acllwr`, `rcl1`, `rcl2` and `acl2` appear in no note, config, audit or task anywhere in this repository | `notes/audits/inherited-earth-constants.md` finding 2 | open. Four smaller absorptances were corrected on exactly this argument -- the two ozone weights, `h2o_sw_weight` 1.346 and `co2_sw_weight` 1.510 -- and clouds are booked at 12.0 W/m2 of Earth's shortwave absorption beside water vapour's term in `exoplasim/notes/shortwave-water-vapour.md`. The band PARTITION is star-aware, since `zsolar1`/`zsolar2` come from the spectrum; what is Earth's is the physics inside each range, and range 2 is where it bites, because this star puts 0.618 of its flux there and concentrates it nearer 0.8-1.5 um while the Sun's range-2 energy reaches further into the 2-3 um water bands. The sign is NOT obvious by inspection, which is the argument for measuring it. All of them are `radmod_nl` keys on the binary already built, so this satisfies every one of `docs/src/pipeline/sequencing.md` A3's four A/B conditions with no further work: bracket and report the spread, do not tune within it, exactly as CLIM-16 handles `hdiffk`. The PREDICTION and the arm design are written, 2026-08-19, in `exoplasim/notes/forcing-bundle-predictions.md`: absorption-like keys bracketed 0.78 to 1.28 spanning both signs, bound 0 +/- 0.6 K, and `acllwr` excluded from stellar arms because a thermal-band constant has no stellar dependence -- it stays in this row as untraced tuning only [step: baseline_run] |
| PHYS-12 | -- | -- | done, see `archive/tasks.md` |

## REF -- references and provenance

0 open of 9 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| REF-1 | -- | -- | done, see `archive/tasks.md` |
| REF-2 | -- | -- | wontfix, see `archive/tasks.md` |
| REF-3 | -- | -- | done, see `archive/tasks.md` |
| REF-4 | -- | -- | done, see `archive/tasks.md` |
| REF-5 | -- | -- | done, see `archive/tasks.md` |
| REF-6 | -- | -- | wontfix, see `archive/tasks.md` |
| REF-7 | -- | -- | done, see `archive/tasks.md` |
| REF-8 | -- | -- | done, see `archive/tasks.md` |
| REF-9 | -- | -- | done, see `archive/tasks.md` |

## SPEC -- the stellar spectrum and how the radiation scheme integrates it

0 open of 5 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| SPEC-1 | -- | -- | done, see `archive/tasks.md` |
| SPEC-2 | -- | -- | done, see `archive/tasks.md` |
| SPEC-3 | -- | -- | done, see `archive/tasks.md` |
| SPEC-4 | -- | -- | done, see `archive/tasks.md` |
| SPEC-5 | -- | -- | done, see `archive/tasks.md` |

## SURF -- derived surface classes

1 open of 6 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| SURF-1 | -- | -- | done, see `archive/tasks.md` |
| SURF-2 | -- | -- | done, see `archive/tasks.md` |
| SURF-3 | -- | -- | done, see `archive/tasks.md` |
| SURF-4 | Diatomite is placed against a STATIC lake proxy, not lake occupancy over the stellar cycle. The strandline band between the solved lake surface and the spill level is the ground a single climatology's lake COULD vacate; the design asks for the band it is measured to vacate over the 57-year component, which is the forcing slow enough for lake level to equilibrate | `pedology/notes/derived-surface-classes.md` | open, and blocked on CYC-1 rather than on effort. It cannot be faked from one climatology: a basin permanently full never exposes its bed and a basin permanently dry never accumulates one, so the whole class is about alternation. Re-run `build_surface_classes.py` against the cycle climatologies when they exist [step: surface_classes, stellar_cycle_run] |
| SURF-5 | -- | -- | done, see `archive/tasks.md` |
| SURF-6 | -- | -- | done, see `archive/tasks.md` |

## VOLC -- volcanism and weathering fluxes

0 open of 7 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| VOLC-1 | -- | -- | done, see `archive/tasks.md` |
| VOLC-2 | -- | -- | done, see `archive/tasks.md` |
| VOLC-3 | -- | -- | done, see `archive/tasks.md` |
| VOLC-4 | -- | -- | done, see `archive/tasks.md` |
| VOLC-5 | -- | -- | wontfix, see `archive/tasks.md` |
| VOLC-6 | -- | -- | done, see `archive/tasks.md` |
| VOLC-7 | -- | -- | done, see `archive/tasks.md` |

