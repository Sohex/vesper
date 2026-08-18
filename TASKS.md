# Tasks

Atomic, trackable work items. **Findings documents do not carry their own to-do
lists** -- an audit says what is true, this file says what to do about it, and
the two are kept apart so a finding can be read without being re-litigated and a
task can be closed without editing the argument that produced it.

Every task cites the document that justifies it. If a task has no source, it is
not yet a finding and probably needs one.

Status: `open` | `doing` | `blocked` | `done` | `wontfix` (with a reason).

## Conventions

- IDs are stable and never reused. Prefix by area: `LITH` lithology, `GRAV`
  gravity, `VOLC` volcanism and weathering fluxes, `SURF` derived surface
  classes, `CLIM` climate, `HYD` hydrography, `BIO` biosphere, `MIN` economic minerals,
  `DUST` the aeolian component, `REF` references and provenance, `CONS`
  consistency checking, `PHYS` physics calibrated for the wrong world, `BUDG`
  the error budget and what it is denominated in, `SPEC` the stellar spectrum
  and how the radiation scheme integrates it.
- One line per task. If it needs a paragraph, it needs a findings document.
- A task that turns out to be wrong is closed `wontfix` with the reason, not
  deleted. The reasoning is the point.
- A closed task moves to `archive/tasks.md` rather than staying here, so this
  file shows only what is left. Ids are never reused, so a cited id that is not
  here is there.
- **A task names the step it touches**, as `[step: <id>]` in its status column,
  where the id is from `config/pipeline.yaml`. That is what makes the carve gate
  computable: `pipeline.py --status` reports which open tasks touch a step
  upstream of `carve_list`. A task that names none is not ignored -- it is
  reported as a RESIDUAL to be arbitrated, because the graph narrows the
  judgement and does not replace it. Name more than one where it applies.

- **Regenerating a derived artifact is a STEP, not a task.** It belongs in the
  ordering in `WORKFLOW.md` section 6, beside the run that consumes it. A row
  that says "rebuild X before the next run" is tracking state, and state is what
  `check_consistency.py` and `world_state.json` are for. HYD-16 was such a row
  and should never have been one. CLAUDE.md rule 7.
- **A task is code, physics, a decision or a measurement that settles a
  mechanism.** If closing it would produce only a number that the next iteration
  regenerates, it is not a task; the durable half is whatever it teaches.

---

## Open

| id | task | source | status |
| --- | --- | --- | --- |
| SURF-6 | Nothing consumes the derived surface classes yet. Two consumers were the reason for building them and both are one field away: `aeolian/scripts/build_dust.py` should read `pavement` as an emission SUPPRESSOR and a store of deposited dust, and `pedology/scripts/phosphorus_budget.py` should read `pavement` as a P sink and `diatomite` as the deflatable source Hudson-Edwards measured | `pedology/notes/derived-surface-classes.md` | open. The pavement sign is the whole point of the Wells et al. (1995) correction: as "deflation armour" the class would have been a product of emission and increased with it, and as pavement it suppresses emission and stores what falls. Leaving the field unread keeps that correction theoretical. Do it AFTER the carve, since both consumers would be re-run on the new terrain anyway [step: surface_classes] |
| SURF-4 | Diatomite is placed against a STATIC lake proxy, not lake occupancy over the stellar cycle. The strandline band between the solved lake surface and the spill level is the ground a single climatology's lake COULD vacate; the design asks for the band it is measured to vacate over the 57-year component, which is the forcing slow enough for lake level to equilibrate | `pedology/notes/derived-surface-classes.md` | open, and blocked on CYC-1 rather than on effort. It cannot be faked from one climatology: a basin permanently full never exposes its bed and a basin permanently dry never accumulates one, so the whole class is about alternation. Re-run `build_surface_classes.py` against the cycle climatologies when they exist [step: surface_classes, stellar_cycle_run] |
| BIO-1 | Regenerate the LPJ-GUESS driver: `vesper_driver_provenance.json` is pinned to a superseded config and `check_consistency.py` has been failing on it | `scripts/check_consistency.py`, `WORKFLOW.md` A3 | deliberately DEFERRED until after iteration 2's baseline run. The driver is built from a climatology and a soil that the carve replaces, so doing it now is make-work whose only effect is silencing a warning that is currently telling the truth: no biosphere inputs exist for this build. The orbit-dependent half was durable and is done (BIO-3) [step: lpj_driver] |
| CYC-1 | Run the stellar cycle off the settled baseline. `WORKFLOW.md` section A2 requires it AFTER a converged baseline and BEFORE the carve verdict, and nothing was tracking it: carving is irreversible, so the terrain ratchets toward the cycle's wet extreme and a verdict taken on the mean climate systematically under-carves. The run also measures the damping factor, known only as 0.24 to 0.6, which converts any future flux amplitude into a climate without another run | `WORKFLOW.md` section A2, `config/planet.yaml` stellar_cycle | open, and ON THE CRITICAL PATH. The cycle binary is built and verified. Length is set in periods of the LONG component, 57 Earth years, not the medium one. The carve list of 2026-08-17 was taken on the mean climate and should be re-taken after this [step: stellar_cycle_run] |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here) | `notes/audits/orogen-gravity.md` | blocked on the downscaling pass, deliberately. The audit's recommendation was Option B, carry it as a declared gap: `iceFlow` is a heuristic accumulation rather than a Glen's-law velocity, so a rigorous g^3 on it is precision theatre. Decisive on top of that, a glacial valley is 1-5 km against a 15.19 km mesh cell, so the process is SUB-GRID and cannot be resolved here at any gravity. It belongs with the downscaling machinery, which has to persist sub-grid hypsometry anyway [step: orogen] |
| DUST-11 | Settle the CATCHMENT half: run one prescribed-dust climate segment and measure what it does to precipitation and runoff. Runoff is the carve criterion's denominator and only 19% of land precipitation, so it amplifies; absorbing aerosol suppresses precipitation through a fast adjustment that no offline calculation can reach, because it is a response of the circulation and not of a surface energy balance | `notes/dust.md`, `aeolian/notes/prescribed-dust-run.md` | DECIDED 2026-08-17: measure first, then decide whether to carry it. The matched pair rides inside the flux re-bracket rather than beside it, so the control is a point loop A needs anyway. The patch is resident with `ndustrad = 0` [step: surface_dust, baseline_run] |
| DUST-3 | Put the emission scheme in the ExoPlaSim fork, keeping the source map outside it as boundary fields. DECIDED; DESIGNED in `aeolian/notes/in-model-dust.md`, which is detailed enough that landing it is mechanical | `aeolian/notes/in-model-dust.md`, `notes/dust.md` | open. "Replace one line" was wrong: six pieces in order -- the three upstream defects, replace the bottom-level sink with a deposition velocity, wet scavenging, boundary fields plus Kok 18, the longwave aerosol term, and size bins only if needed. The longwave is the dominant RISK (inside the radiation solver) and `aerocore` running serial on NROOT is the item most likely to be UNDERESTIMATED, since nothing the emission needs is gathered. Five or six patches to a compiled model, each needing a rebuild per CLAUDE.md rule 4. SEQUENCING CHANGED 2026-08-18: A3 no longer asks for an iteration of its own, because attribution cannot gate a correct term and so cannot justify a converged run. Each piece states a prediction and is tested by a short A/B off a common restart, and the set rides one converged run. Real ordering constraints remain two: the three upstream defects come first, and the longwave term cannot be developed in parallel with another `radmod.f90` patch ITEMS 1 TO 3 AUTHORED 2026-08-18 as two patches in `PENDING_PATCHES`, each with its prediction and falsifier. Three defects became seven, and the largest is that the bottom model level was updated twice so the sedimentation flux out of the atmosphere cancelled itself -- which is what the 99%/step sink existed to hide. What remains: item 4, the three `.sra` boundary codes with their generator, `surfcode`/`mpsurfgp` registration, the gathers, WRITING `aero_namelist` FROM `dust.yaml` so the new keys are reachable at all, and Kok 18 in the source case; item 5, the longwave term, blocked on nothing now that PHYS-6 has merged; and upstream defect 1, `radmod`'s `apart` never populated, same gate. Item 6 is closed by DUST-8. `aerocore` running serial on NROOT is still the item most likely to be UNDERESTIMATED [step: rebuild_binaries] |
| HYD-3 | Measure the overshoot: re-evaluate an already-carved set against the climate that carving produced, and report how many would no longer have carved | `WORKFLOW.md` section 4 | blocked -- it applies to a verdict taken on carved terrain, and a pre-carve base restarts the count at iteration 1 [step: carve_verdict] |
| SPEC-1 | Re-run the baseline on the spectrum the config declares. Every orbit after the first of every run on this build was integrated against a 4965 K blackbody, worth +0.024 on broadband snow albedo and moving every derived surface albedo with it | `notes/audits/physics-review.md` finding 2 | DECIDED 2026-08-17: declare it. Also never a decision -- the config NAMES `k25v` and the model not reading it was a defect against the config, not an option. Half of it has already happened by itself: the driver fix means `continue_exoplasim` stages the spectrum, so orbits 77-80 ran on it. SPEC-2 has landed, so it is now safe to. What remains is the baseline re-run, which loop A does anyway [step: baseline_run] |
| CLIM-11 | The surface energy residual is entirely OCEAN, and entirely on ocean cells that never carry ice or snow, where the slab identity is exact by construction and fails by -0.29 W/m2 per unit area | `exoplasim/notes/water-and-energy-closure.md` | open, and narrowed as far as the present instruments reach. Land closes to -0.016 and that is a test that could have failed, since `landmod` integrates exactly `hfns` into exactly the column `close_state_energy.py` rebuilds, so soil capacity, the snow latent term and melt booking are all clean. On the ice-free ocean the slab has `NLEV_OCE = 1`, `MLDEPTH = 50`, no horizontal diffusion, no flux correction and no coupling-interval mismatch, and `seamod.f90` accumulates exactly `rss+rls+hfss+hfls`, so the identity should be exact. Seasonal regression of `hfns` on the slab storage gives slope 1.0946 at correlation 0.99435 over a 43 W/m2 swing, so it is an OFFSET and not a scaling. Excluded: mixed-layer depth, ocean transport, flux correction, deep-ocean flux, coupling interval, the sea-ice reservoir. NEXT INSTRUMENT is the ocean output stream, `osst` code 169 and `oheat` code 263, which this project does not postprocess -- not more model time, and not the offline radiative transfer previously proposed [step: assess_convergence] |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; the span is 18.2% of NPP and it is the largest nitrogen lever | `biosphere/notes/productivity-prediction.md` | blocked on iteration 2's baseline run -- no LPJ-GUESS run exists on this build [step: lpj_run] |


## Closed

Moved to `archive/tasks.md`, with the reasoning that closed each one. Ids are
never reused, so look there when a commit cites an id this file does not have.
