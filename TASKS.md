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
  where the id is from `config/pipeline.yaml`. That is what makes the carve gate
  computable: `pipeline.py --status` reports which open tasks touch a step
  upstream of `carve_list`. A task that names none is not ignored -- it is
  reported as a RESIDUAL to be arbitrated, because the graph narrows the
  judgement and does not replace it. Name more than one where it applies.
  `pipeline.py` reads STATUS to decide what is open, so a row whose status
  begins `done` or `wontfix` is skipped wherever it sits.

- **Regenerating a derived artifact is a STEP, not a task.** It belongs in the
  ordering in `WORKFLOW.md` section 6, beside the run that consumes it. A row
  that says "rebuild X before the next run" is tracking state, and state is what
  `check_consistency.py` and `world_state.json` are for. HYD-16 was such a row
  and should never have been one. CLAUDE.md rule 7.
- **A task is code, physics, a decision or a measurement that settles a
  mechanism.** If closing it would produce only a number that the next iteration
  regenerates, it is not a task; the durable half is whatever it teaches.

---

## BIO -- biosphere

2 open of 4 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| BIO-1 | Regenerate the LPJ-GUESS driver: `vesper_driver_provenance.json` is pinned to a superseded config and `check_consistency.py` has been failing on it | `scripts/check_consistency.py`, `WORKFLOW.md` A3 | deliberately DEFERRED until after iteration 2's baseline run. The driver is built from a climatology and a soil that the carve replaces, so doing it now is make-work whose only effect is silencing a warning that is currently telling the truth: no biosphere inputs exist for this build. The orbit-dependent half was durable and is done (BIO-3) [step: lpj_driver] |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; the span is 18.2% of NPP and it is the largest nitrogen lever | `biosphere/notes/productivity-prediction.md` | blocked on iteration 2's baseline run -- no LPJ-GUESS run exists on this build [step: lpj_run] |
| BIO-3 | -- | -- | done, see `archive/tasks.md` |
| BIO-4 | -- | -- | done, see `archive/tasks.md` |

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

4 open of 30 issued.

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
| CLIM-11 | The -0.29 W/m2 slab residual is NOT the model: read against the ocean and ice output streams the identity closes to +0.002 W/m2, and -0.158 of it is pyburn weighting a 570-timestep first bin equally with eleven 480-timestep ones | `exoplasim/notes/water-and-energy-closure.md`, `exoplasim/scripts/close_ocean_energy.py` | open, and LARGER than this row recorded, because the attribution that reduced it was refuted on 2026-08-18. The instrument this row previously named was WRONG: `osst`/`oheat` on 169/263 are pyburn's atmospheric table where `hfns` is derived arithmetic, while the ocean stream is raw service format on codes 901-906, and `ice_output` beside it -- which nobody had named -- is the better one, decomposing what the atmosphere delivers into ocean, ice, snow and surface temperature. Six identities close, each with a right answer of zero, the largest at 0.0024 W/m2 against a -25 to +16 seasonal swing; all three candidate shapes are excluded by measurement rather than argument, and `ts` is the slab temperature to 0.001 K. The bin weighting is NOT the cause: correctly derived it is +0.0734 W/m2 of planet on never-ice ocean against a -0.1517 residual, the WRONG sign, so counting it enlarges the residual instead of halving it. See CLIM-13 for why the old -0.0820 was an artefact of assuming eleven equal bins. The remainder cannot be pinned on orbits 67-76: their ocean stream was overwritten and the four estimators of the ten-orbit storage span 0.18 W/m2. NEXT INSTRUMENT is CLIM-12, which LANDED 2026-08-18: the streams now survive the whole run, so the ten-orbit storage whose four estimators spanned 0.18 W/m2 becomes an exact endpoint difference. That needs a run made after the fix; every run on this build predates it [step: bootstrap_convergence, baseline_convergence] | UPDATED 2026-08-19: the precondition is MET and a second obstacle replaced it. `run_4182235e9781` is the first run on this build made after CLIM-12, and it carries `MOST_OCEAN.NNNNN` and `MOST_ICE.NNNNN` for all 40 orbits. But `close_ocean_energy.py` REFUSES it: the ice stream does not tile the regular output at any offset, scoring 13.102 W/m2 against a worst-offset 16.657 and a seasonal scale of 31.979, so it fails both declared criteria -- beat the worst by five, and land inside 1% of the swing. It is NOT an indexing slip, which was the first hypothesis: each stream file is written about 30 s before its own `.nc`, that gap being pyburn, so index N is orbit N in both. The refusal is the check working and the cause is unknown; candidates not yet separated are the seeded start (orbit 0 follows a restart taken at different physics, h2oswl 1.0 against 1.127) and the run never having reached equilibrium. Next step is a window argument on the script so the seeded orbits can be excluded, which would test the first candidate and cannot be done by staring at it DIAGNOSED 2026-08-19 and the instrument was the defect. The refusal at 13.102 W/m2 was TILING DRIFT, not physics: `align` floored `nrec // nbin`, Vesper's year is 182.801 days so a daily stream writes 183 records against 12 bins, and 183/12 is 15.25 -- each bin a quarter-record short of the next, accumulating at about 0.42 W/m2 per orbit, which is why the score scaled linearly with window. Fractional tiling removes it and is verified to reduce to the old reshape exactly where the division is whole, 2.2e-16. WHAT IS LEFT IS THE INTERESTING PART: 0.24 to 0.32 W/m2, FLAT with window rather than falling as 1/sqrt(n), so systematic and not sampling noise; per-bin at one orbit it is +0.19 net with +/-0.3 scatter. That is the same order as the -0.29 this task exists to explain, so the instrument cannot currently resolve its own subject -- which reframes CLIM-11 from a model question to an instrument one. Leading hypothesis, being tested: at NLOWIO = 0 the regular output is INSTANTANEOUS SAMPLES while the stream is an ACCUMULATION, and a mean of samples is not a time-mean. A matched low-io/clean-io pair off one seed is running. Separately, `align`'s five-fold discrimination criterion is now VACUOUS rather than failed -- fractional tiling consumes every record so the offset scan has one candidate -- and whether a guard survives its hazard being fixed is an open judgement; relaxing it would not produce a pass either way HYPOTHESIS REFUTED 2026-08-19 by a matched pair off one seed, `run_ef3e195a7bdd` clean and `run_277e52971ec5` low-I/O, three orbits each, differing in nothing but the regime. If the residual were a statistic mismatch -- an accumulation compared against a mean of instantaneous samples -- making BOTH sides accumulations should have collapsed it. It went 0.353 to 0.861 W/m2, worse. So the disagreement is between two accumulations of the same energy flow and is not an artefact of how the regular output is sampled, which puts CLIM-11 back on the model side with 0.353 as the number to explain and the same order as the -0.29 it started from. Found on the way and fixed: the closure read its land mask from BIN 0, the one bin `align` eight lines below already excludes for straddling the restart, and under low I/O that bin carries a partial accumulation window -- the land mask, a field that cannot vary, reads 0.9534 there and exactly 1.0 in every other bin of every orbit. The script disagreed with itself and low-I/O output was unreadable to it as a result
| CLIM-12 | -- | -- | done, see `archive/tasks.md` |
| CLIM-13 | -- | -- | done, see `archive/tasks.md` |
| CLIM-14 | -- | -- | done, see `archive/tasks.md` |
| CLIM-15 | -- | -- | done, see `archive/tasks.md` |
| CLIM-16 | Ocean horizontal heat diffusion has a namelist switch, `nhdiff` with `hdiffk`, and the error budget records no-q-flux as STRUCTURAL with "no cheap version exists" | `notes/audits/absent-and-inherited-physics.md` finding 1, `analysis/error_budget.json` other_items | open, and it is the cheapest structural term in the budget to price rather than the impossible one. `oceanmod.f90` reads both keys in `oceanmod_namelist`, broadcasts them and acts on `nhdiff > 0` at line 844, so it is one key on the binary already built -- which satisfies every one of WORKFLOW A3's four A/B conditions with no further work. Bracket `hdiffk` and report the spread rather than tuning it, the way this project handles every other ungrounded constant. Read it as physics-is-not-a-knob: ocean heat transport EXISTS, the argument for enabling it is not that it improves an agreement, and a constant diffusivity gives a BOUND on the missing transport rather than the transport. The error budget entry is CORRECTED 2026-08-18: its second sentence, that no cheap version exists, was false and the row now records that `nhdiff` is one namelist key on the built binary. What is left is the measurement -- bracket `hdiffk`, run the paired arms off one restart, report the spread. The PREDICTION is written, 2026-08-19, with its wrongness bounds: `exoplasim/notes/forcing-bundle-predictions.md`, computed by `exoplasim/scripts/predict_ocean_terms.py`, whose operator is validated against a Laplacian eigenfunction and whose applied flux the run itself now records in the ocean stream (CLIM-12), so the A/B has a pointwise right answer [step: baseline_run] |
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
| CLIM-30 | Re-derive the design flux on this build's corrected physics, and DECLARE the extreme-cold land cap in advance rather than inferring it | `WORKFLOW.md` section 5b, `exoplasim/scripts/derive_design_flux.py` | blocked, and on data rather than effort: the script needs two converged sources spanning more than 5 K and only one flux exists on disk, both run payloads being 0.945. The 0.91 bracket, `run_bfa3f5269660`, is an identity stub with no data files. Recovering the purged artifact from `33949a3^` is what makes this urgent rather than routine: the unconstrained winner was 0.8725 dry and 0.875 equivalent, and 0.945 was reproduced ONLY under an extreme-cold cap of 0.0642 that the artifact itself labels `inferred` -- a threshold fitted to the answer it was meant to test, which is the thing section 7 forbids. So the recorded flux is unsupported until the cap is declared ahead of the run, and the declaration is a DECISION nothing in the project settles. Needs a second converged point at roughly 0.91 after the baseline [step: design_flux, baseline_run] |

## CONS -- consistency checking

1 open of 8 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| CONS-1 | -- | -- | done, see `archive/tasks.md` |
| CONS-2 | -- | -- | done, see `archive/tasks.md` |
| CONS-3 | -- | -- | done, see `archive/tasks.md` |
| CONS-4 | -- | -- | done, see `archive/tasks.md` |
| CONS-6 | -- | -- | done, see `archive/tasks.md` |
| CONS-7 | -- | -- | done, see `archive/tasks.md` |
| CONS-8 | -- | -- | done, see `archive/tasks.md` |
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
| CYC-1 | Run the stellar cycle off the settled baseline. `WORKFLOW.md` section A2 requires it AFTER a converged baseline and BEFORE the carve verdict, and nothing was tracking it: carving is irreversible, so the terrain ratchets toward the cycle's wet extreme and a verdict taken on the mean climate systematically under-carves. The run also measures the damping factor, known only as 0.24 to 0.6, which converts any future flux amplitude into a climate without another run | `WORKFLOW.md` section A2, `config/planet.yaml` stellar_cycle | open, and ON THE CRITICAL PATH. The cycle binary is built and verified. Length is set in periods of the LONG component, 57 Earth years, not the medium one. The carve list of 2026-08-17 was taken on the mean climate and should be re-taken after this [step: stellar_cycle_run] |

## DUST -- the aeolian component

2 open of 13 issued.

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

1 open of 17 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| HYD-1 | -- | -- | done, see `archive/tasks.md` |
| HYD-2 | -- | -- | done, see `archive/tasks.md` |
| HYD-3 | Measure the overshoot: re-evaluate an already-carved set against the climate that carving produced, and report how many would no longer have carved | `WORKFLOW.md` section 4 | blocked -- it applies to a verdict taken on carved terrain, and a pre-carve base restarts the count at iteration 1 [step: carve_verdict] |
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
| PHYS-11 | The cloud shortwave and longwave constants are Earth tunings and are the only major radiation term never re-weighted for this star. `tswr1`, `tswr2`, `tswr3`, `acllwr`, `rcl1`, `rcl2` and `acl2` appear in no note, config, audit or task anywhere in this repository | `notes/audits/inherited-earth-constants.md` finding 2 | open. Four smaller absorptances were corrected on exactly this argument -- the two ozone weights, `h2o_sw_weight` 1.346 and `co2_sw_weight` 1.510 -- and clouds are booked at 12.0 W/m2 of Earth's shortwave absorption beside water vapour's term in `exoplasim/notes/shortwave-water-vapour.md`. The band PARTITION is star-aware, since `zsolar1`/`zsolar2` come from the spectrum; what is Earth's is the physics inside each range, and range 2 is where it bites, because this star puts 0.618 of its flux there and concentrates it nearer 0.8-1.5 um while the Sun's range-2 energy reaches further into the 2-3 um water bands. The sign is NOT obvious by inspection, which is the argument for measuring it. All of them are `radmod_nl` keys on the binary already built, so this satisfies every one of WORKFLOW A3's four A/B conditions with no further work: bracket and report the spread, do not tune within it, exactly as CLIM-16 handles `hdiffk`. The PREDICTION and the arm design are written, 2026-08-19, in `exoplasim/notes/forcing-bundle-predictions.md`: absorption-like keys bracketed 0.78 to 1.28 spanning both signs, bound 0 +/- 0.6 K, and `acllwr` excluded from stellar arms because a thermal-band constant has no stellar dependence -- it stays in this row as untraced tuning only [step: baseline_run] |
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

1 open of 5 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| SPEC-1 | Re-run the baseline on the spectrum the config declares. Every orbit after the first of every run on this build was integrated against a 4965 K blackbody, worth +0.024 on broadband snow albedo and moving every derived surface albedo with it | `notes/audits/physics-review.md` finding 2 | DECIDED 2026-08-17: declare it. Also never a decision -- the config NAMES `k25v` and the model not reading it was a defect against the config, not an option. Half of it has already happened by itself: the driver fix means `continue_exoplasim` stages the spectrum, so orbits 77-80 ran on it. SPEC-2 has landed, so it is now safe to. What remains is the baseline re-run, which loop A does anyway [step: baseline_run] |
| SPEC-2 | -- | -- | done, see `archive/tasks.md` |
| SPEC-3 | -- | -- | done, see `archive/tasks.md` |
| SPEC-4 | -- | -- | done, see `archive/tasks.md` |
| SPEC-5 | -- | -- | done, see `archive/tasks.md` |

## SURF -- derived surface classes

2 open of 6 issued.

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

