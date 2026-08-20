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

6 open of 44 issued.

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
| CLIM-16 | -- | -- | done, see `archive/tasks.md` |
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
| CLIM-29 | Carbonaceous aerosol is the one the sulfate bound does not reach: smoke and secondary organics ABSORB, so their forcing per unit optical depth is larger than either scatterer's and of the opposite sign. Fire is enabled in LPJ-GUESS as GLOBFIRM and biogenic emissions are already in the driver, so both are estimable from a biosphere run. The SOA half is a switch rather than a port: LPJ-GUESS carries its own biogenic VOC scheme in `bvoc.cpp` and it is merely OFF, which is recorded in `biosphere/README.md` only in passing, as the one place `climate.dtr` is read | `notes/audits/unpriced-terms.md` finding 2, `aeolian/README.md` | blocked on a biosphere run existing on this build, which is BIO-1 and BIO-2, and deliberately not worked around: the burned area and the biogenic emission are LPJ-GUESS output and inventing either would be inventing the vegetation the whole component exists to compute. The chain to reuse is complete -- OPAC carries a soot component beside the sulfate one, `sea_salt_optics.band_average` takes any distribution, and the transport and the two-stream are shared -- so this is a source term and an optics table and nothing else. Note the sign: an absorbing aerosol over this world's bright closed-basin fill warms, which is the same asymmetry `dust_forcing.py` already prices for mineral dust [step: lpj_run, sea_salt_optics] |
| CLIM-30 | DECLARE the extreme-cold land-fraction cap BEFORE the design flux is re-derived, instead of inferring it from the answer. The purged artifact inferred 0.0642, which reproduced 0.945 from an unconstrained winner of 0.8725 dry -- a threshold fitted to the answer it was meant to test, which is what the threshold convention (`docs/src/practice/conventions.md`) forbids | `docs/src/pipeline/state.md` section 5b, `exoplasim/scripts/derive_design_flux.py` | open, and a DECISION nothing in the project settles: it needs a preference no document fixes, which is why it is a task and the re-derivation is not. Re-deriving is the `design_flux` step, and the second converged point it needs is `bracket_run`; neither is tracked here. Until the cap is declared ahead of the run, the recorded flux is unsupported [step: design_flux] |
| CLIM-31 | -- | -- | done, see `archive/tasks.md` |
| CLIM-34 | Redo the glacier bound against the real 0.945 climatology with the measured 7.8 K/km lapse; the recorded -4.7 K offset is a proxy | `notes/glacier-rough-pass.md`, archive PHYS-12 | open [step: baseline_climatology] |
| CLIM-35 | -- | -- | done, see `archive/tasks.md` |
| CLIM-36 | -- | -- | done, see `archive/tasks.md` |
| CLIM-32 | -- | -- | done, see `archive/tasks.md` |
| CLIM-33 | -- | -- | done, see `archive/tasks.md` |
| CLIM-37 | -- | -- | done, see `archive/tasks.md` |
| CLIM-38 | -- | -- | done, see `archive/tasks.md` |
| CLIM-39 | -- | -- | done, see `archive/tasks.md` |
| CLIM-40 | Route sea salt into the radiation: the `.sra` column writer on a new surface code, its aerofile column, the `run_exoplasim.py` code set and the `config/planet.yaml` key, on the pattern `build_surface_dust.py` and `dust_aerofile.py` already set for dust | `aeolian/notes/multi-species-aerosol.md` section 6 | open. CLIM-39 has LANDED, and the blocker is now data rather than machinery: `aeolian/analysis/sea_salt_baseline.nc`, the column optical depth the writer must consume, does not exist -- `build_sea_salt.py` needs a climatology and the `source/` payload, so the chain bottoms out at a climate run. The optics half is unblocked, `analysis/sea_salt_optics.json` being tracked and present. One finding to carry in: a PRESCRIBED species' aerofile column holds only RATIOS -- ssa, backscatter, band split and qlw -- with the column field carrying the magnitude, so sea salt needs four ratios per band and none of `dust_aerofile.py`'s burden-matched radius or absolute cross-section, which are specific to the interactive path where `aeroprof` actually uses `apart`. It is the species that changes an answer -- the same magnitude as dust and the opposite sign -- and the cheaper half of the pair to prescribe, being a pure scatterer over a dark surface that does not change. Volcanic sulfate rides the same machinery at one to two orders down and decides nothing, so it goes in with this or not at all [step: sea_salt]
| CLIM-41 | -- | -- | done, see `archive/tasks.md` |
| CLIM-42 | Add a CH4 and N2O band to the longwave. `config/planet.yaml` says the shape of it -- "adding either means adding a band, not a key" -- and CLIM-41 measured that the term outranks every other item in the budget | `exoplasim/notes/trace-gas-band.md` is the design and the declared tests, `exoplasim/notes/trace-gas-absorbers.md` section 5, `analysis/trace_gas_forcing.json` | open. DESIGN LANDED 2026-08-20 with the tests declared ahead of the code, on CLIM-39's precedent. The parameterisation is Donner and Ramanathan (1980), which did exactly this job and is now on disk and read: the Cess and Ramanathan (1972) band absorptance as modified by Ramanathan (1976), with Table 1's `A0` and `beta0` for CH4 1306 cm-1 and N2O 1285 and 589 cm-1. It is the right source rather than the convenient one because it is a BAND ABSORPTANCE, which is the currency Sasamori's scheme is already written in -- Byrne and Goldblatt's fits price the omission but are a top-of-atmosphere forcing, and putting a global-mean forcing inside a layer-by-layer solver is a category error however well it reproduces the mean. It also carries the water vapour overlap by multiplying the band absorptance, which is line for line what `lwr` already does to CO2 with `zth2o`. CH4's band intensity is RECOVERED and validated: `S` = 187.7 cm-1 (cm atm)-1 reproduces ten of the eleven rows of the paper's own Table 2 to 1.6%, the eleventh being a typesetting slip that reprints the cell diagonally below it. TWO THINGS LEFT, and they are different in kind. **N2O's two band intensities are not recovered** -- the paper's evidence for 1285 cm-1 is a figure and 589 cm-1 has neither, so the routes are reading Fig. 2's curve A at both pressure panels, which is two constraints on one unknown and checks itself, or obtaining McClatchey et al. (1973); until one lands there is no defensible N2O implementation and CH4 alone is the larger half. **The build and every test are BLOCKED on the host**, which has timing-sensitive work on it: rule 4 makes a rebuild of all five binaries mandatory after the source change, and the zero-abundance reduction identity needs runs. Sequencing conflict with CLIM-39 is GONE; that landed [step: rebuild_binaries] |
| CLIM-43 | -- | -- | done, see `archive/tasks.md` |
| CLIM-44 | The model is not run-to-run reproducible, and nothing owns it: the SAME binary on identical inputs gives bit-identical output at 1 timestep and differs by 16, and `plasim_status` differs between two runs even at 1 timestep, so a bitwise restart diff reports nothing | `aeolian/notes/multi-species-aerosol.md`, measured during CLIM-39 | open. Found while establishing CLIM-39's reduction identity, which is why that test is stated at one timestep. Old-vs-new differs at the SAME byte as the control pair and grows 2.4e-6 to order 1, so it is the model's own chaos rather than that change -- but it means ANY A/B here must establish its own reproducibility horizon first, and several planned ones assume bit-identity. Cause not established; 8 ranks and MPI reduction ordering is the first look and a one-rank run separates it [step: baseline_run EVIDENCE ADDED 2026-08-20 from the build-flag benchmark, which repeated runs for a different reason and therefore did not select for this: at T42 on **16 ranks**, 600 steps from one restart with `NOUTPUT = 0`, `plasim_status` came out BIT-IDENTICAL across about 48 repeats spanning four distinct binaries -- 18 runs of the stock build all at 4497a8d7, 18 of the znver4 build all at 2eb4e459, plus six each of two more. So the model is not non-reproducible in general, and whatever CLIM-39 hit is configuration-specific. That narrows this row's own first hypothesis rather than confirming it: the asymmetry runs the wrong way for a plain reduction-ordering story, since 16 ranks reduce over more partial sums than 8 and are the reproducible case here. A THIRD candidate now outranks both, found by reading rather than running: `git merge-base --is-ancestor 80f11e9 b43f17e` is FALSE, so the CLIM-39 measurement was taken on a worktree that does NOT carry the `zsolars` fix. That defect wrote a record of 8192 doubles of which 8190 were memory past the end of a 2-element array, held stable only by `-finit-real=zero`, and `notes/audits/zsolars-restart-overread.md` measured it returning a different restart sha on every one of five runs from one binary on one input. That is sufficient for the `plasim_status` half and it is already fixed; it is NOT sufficient for the gridpoint-output half, because nothing reads the record back. Also ruled out from source: the model's only random number source is `get_random_root`, gated on `nstorain` which is 0 in every run here, and a resume restores the generator state from the restart in any case -- so a clock-seeded RNG is not it. What is left needs runs. `exoplasim/scripts/reproducibility_matrix.py` is the harness: a full factorial over ranks, `NOUTPUT` and segment length, 3 repeats a cell, 24 runs, with the pattern each hypothesis would show declared in the docstring ahead of any run and `--plan` to print it. A SECOND correction, and it is to CLIM-39's control rather than to its change: a gridpoint record is written on `mod(nstep, nafter) == 0` over the ABSOLUTE step count, and on this project's production restart `nstep` is 589672 against a logged `nafter` of 32, so the next write is 24 steps away and NEITHER the 1-step nor the 16-step segment writes a record at all. An output comparison over files holding no model records cannot fail, so that control has to be re-taken at segment lengths derived from the cadence, which the harness now does. The accumulation class of `exoplasim/notes/first-output-bin.md` is NOT the cause here -- `naccuout` and the divisor are restored identically in both runs, so it is deterministic, and its signature is a corrupt FIRST record with clean ones after, against a divergence that grew with time from 2.4e-6 -- but it is why the first record after a resume is unlike its successors, and the long segment now straddles it deliberately. Blocked on permission to run the model on this host (`notes/audits/aocl-and-model-build-flags.md` has the bed). |

## CONS -- consistency checking

0 open of 11 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| CONS-1 | -- | -- | done, see `archive/tasks.md` |
| CONS-2 | -- | -- | done, see `archive/tasks.md` |
| CONS-3 | -- | -- | done, see `archive/tasks.md` |
| CONS-4 | -- | -- | done, see `archive/tasks.md` |
| CONS-5 | -- | -- | done, see `archive/tasks.md` |
| CONS-6 | -- | -- | done, see `archive/tasks.md` |
| CONS-7 | -- | -- | done, see `archive/tasks.md` |
| CONS-8 | -- | -- | done, see `archive/tasks.md` |
| CONS-9 | -- | -- | done, see `archive/tasks.md` |
| CONS-10 | -- | -- | done, see `archive/tasks.md` |
| CONS-11 | -- | -- | done, see `archive/tasks.md` |

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

5 open of 16 issued.

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
| DUST-13 | Decide whether the final climate carries interactive emission (`L_AERO = 1`): the reopening test in `notes/dust.md` fired, so a prescribed field is not defensible for a converged answer | `notes/dust.md`, `aeolian/notes/in-model-dust.md` | open, blocked on DUST-11's measurement -- the prescribed-dust run prices what the interactive scheme must reproduce. Take it AFTER CLIM-39 or the choice forecloses sea salt: `ndustrad` and `iaerint` cannot both be on, so interactive dust on today's code puts the one species of comparable magnitude permanently out of the radiation [step: baseline_run] |
| DUST-15 | Make the measured gust samples the dust step's default input and regenerate the chain: every dust artifact on disk is 70x too thin, built from the snapshot-fitted wind tail | `notes/dust.md` last section, `config/pipeline.yaml` dust step | open, and DUST-11's run rides the wrong field until it lands [step: dust] |
| DUST-14 | Dust deposition reaches the soil and never reaches the cryosphere: nothing consumes the deposition field as a term in the model's snow and ice albedo, and `build_dust.py` reads snow only as an emission suppressor | `notes/audits/absent-and-inherited-physics.md` finding 2 | open, one-signed, and it lands on the term the glacier result turns on. Deposition over snow-covered land at 50-60 degrees is 8.75 g/m2/yr at the central aeolian roughness against a terrestrial dust-on-snow literature working at 1-5 g/m2 snowpack loads for albedo reductions of 0.03-0.08, and the smooth end of the roughness bracket is 111 g/m2/yr. The model's snow albedo carries a time-since-snowfall aging range of about a quarter in band 1 and no dependence on what has landed on it. It belongs in the glacier mass balance when `notes/glacier-rough-pass.md` stops being a temperature criterion, and the field it needs already exists. It also puts a cryosphere term under the aeolian roughness bracket, which was understood as controlling emission and the direct forcing only [step: dust, surface_albedo] |
| DUST-16 | Basin deflation driven by nocturnal-inversion low-level jets is sub-resolution: the Bodele-class mechanism (`washington2005-bodele-low-level-jet.pdf`) lives in the lowest few hundred metres and ten sigma layers cannot form it, so gust sampling measures the tail the model produces and not the jet it cannot. One-signed: basin emission is understated wherever inversion jets dominate, and the shallower boundary layer at 1.31 g sharpens exactly that regime | `notes/external-review-abl-obliquity-dust-carbon.md` point 1 | open as a DECLARED gap on GW-6 and GRAV-6's precedent: the process is real, it is not representable at L10, and a parameterisation fitted to nothing would be precision theatre. The honest revisit is vertical resolution, L20 at T85 under loop D, where `close_term_energy.py` re-measures the energy terms anyway. Until then the dust artifacts carry the gust-tail emission and this row is why the basin end of it is a floor [step: dust] |

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

## GW -- groundwater

8 open of 14 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| GW-1 | -- | -- | done, see `archive/tasks.md` |
| GW-2 | Map Orogen's lithology onto Gleeson's five combined hydrolithologies and carry the per-class sigma as a bracket rather than collapsing it. The classes are the SAME Duerr set `pedology/config/pedogenesis.yaml` already maps for Hartmann phosphorus, so this is an existing mapping reused, not a new one | `hydrography/notes/groundwater-scoping.md` section 3, `references/INDEX.md` Gleeson (2011) | open. **Evaporite has no assigned permeability** -- Duerr's EV sits in Gleeson's "not assigned" row with water and ice -- and it must stay unassigned rather than take the nearest class, over a world with this much playa and salt crust. Class spread is 3.4 orders of magnitude and within-class sigma 1.5 to 2.5, so the output is a bracket [step: soil] |
| GW-3 | Calibrate against Earth: the same solver on Earth topography, recharge and GLHYMPS permeability, scored at OBSERVATION LOCATIONS against measured water table depth, threshold declared before scoring | `docs/src/reference/external-data.md`, `hydrography/data/earth_validation/` | open, and the OBSERVATIONS ARE IN HAND: 75,321 Australian sites behind 53.1 million readings, against the 83,617 Australian sites Fan had, assembled by `build_earth_wtd_sites.py` with provenance. Australia is the arid analogue and the regime this world occupies, which is why it went first. Recharge is settled and the worry was backwards -- Berghuijs et al. (2022) is 99% arid-weighted and excludes only permafrost. US wells are open via the OGC API that replaced NWIS, `parameter_code=72019`. Copernicus DEM and GLHYMPS are open. WHAT REMAINS IS THE BUILD: an Earth-side mesh at comparable cell size, the same Voronoi geometry, and a score at observation locations with the threshold declared first. Two validity conditions are already recorded: never model-everywhere against observations-everywhere, since wells sit where people live; and `T = K D` is not Fan's formulation, so a mismatch may be formulation rather than implementation and saying which is part of the result [step: hydrography] |
| GW-4 | -- | -- | done, see `archive/tasks.md` |
| GW-5 | Evaporation from a shallow water table in discharge zones is a sink ExoPlaSim structurally cannot represent: `landmod.f90` throttles evaporation on `dwatc/dwmax` and the cell dries out, so groundwater-fed playa and oasis cells under-evaporate | `hydrography/notes/groundwater-scoping.md` sections 2 and 7, `notes/audits/absent-and-inherited-physics.md` | open, one-signed, and NOT fixable by raising `dwmax` through `surface_soil_water`: a deeper bucket changes storage and timing, not the availability floor a water table sets, so it would move the answer without representing the mechanism. The real routes are a lower boundary in `landmod.f90`, which drags CLAUDE.md rule 4 and a full rebuild, or a fork of the LPJ-GUESS soil column [step: surface_soil_water] |
| GW-6 | Valley-to-ridge water table texture is sub-grid and stays sub-grid: Fan et al. (2013) put the well-articulated gradient at decameters to kilometres and find terrain signals dominate at local scales, against a 15.19 km region | `hydrography/notes/groundwater-scoping.md` section 5 | open as a DECLARED gap, on GRAV-6's precedent and for its reason: the process is real, it is not resolvable at this mesh, and a sub-grid parameterisation of it would be precision theatre. What IS resolvable is the regional recharge control and Fan's own basin-scale convergence, which is the channel GW-4 turns on. Revisit under loop D [step: hydrography] |
| GW-7 | Correct the fork comment that calls triangle centroids "Voronoi vertices on the sphere". Comment only: `cell_area` is the centroidal dual and must NOT change, being the denomination of Orogen's own basin catalogue, hypsometry and the published `basins.drainageConsistency` identity | `notes/audits/mesh-dual-area.md` | open. Measured ON THE EXPORT, the circumcentre dual reproduces `4 pi R^2` at the manifest radius to TEN significant figures and the centroidal one gives 1.00068, so only one tiles. But the reason the solver reconstructs Voronoi faces is not area: a two-point flux approximation needs K-ORTHOGONAL faces, and the median face-to-generator angle is 90.000 deg against 75.746, with the operator missing its analytic eigenvalue by 165/95/67/52 on centroidal faces against 0.11 on Voronoi. The solver reads `cell_area` and never rewrites it, so nothing upstream changes and no build is implied [step: hydrography] |
| GW-11 | -- | -- | done, see `archive/tasks.md` |
| GW-8 | The discrete elliptic operator's noise floor is about 12% relative RMS against the analytic Legendre eigenvalue above l=1, distributed truncation rather than bad faces: the worst 10,000 cells carry under a tenth of the squared error and excluding every sliver-touching cell moves it by 0.002 | `hydrography/notes/mesh-geometry.md` | open. Declared bar was 0.10 and l=1 passes at 0.0217. Reported rather than tuned away; it bounds what any mesh-scale elliptic solve on this export can claim [step: hydrography] |
| GW-9 | -- | -- | done, see `archive/tasks.md` |
| GW-10 | -- | -- | done, see `archive/tasks.md` |
| GW-12 | -- | -- | done, see `archive/tasks.md` |
| GW-13 | Sweep for other consumers that read the flux sign convention from a docstring rather than from the convention itself. `groundwater_receiver` had its sign inverted for the whole of its life -- positive flux is `dst` into `src`, so what leaves `src` is `-flux`, and it traced every cell to the neighbour it receives most water FROM, following the water uphill | `hydrography/notes/water-table-convergence.md` | open. The same convention is used by `geom_divergence`, `closure` and the coastal flux. It survived because closure holds for either sign of a symmetric exchange and the catchment comparison it broke was itself mis-specified, so nothing could see it; what caught it was an identity with a right answer, the trace reproducing `terminal` from a known flux field [step: hydrography] |
| GW-14 | Make `carve_verdict.py` key on TOTAL water delivered to a basin rather than on surface runoff alone: `r_eff = (surface_runoff * C + Qg) / C`. Keep the never-carves branch for `r_eff <= 0`, which is correct for a basin that genuinely receives nothing | `hydrography/notes/water-table-convergence.md`, archive GW-4 | open. RESOLVED IN PRINCIPLE, and the branch was never the defect: "no water therefore never overflows and never incises" is right physics, and what was wrong is that `runoff` meant SURFACE runoff because that was the only supply a surface-only model had. The infinite index for the 605 fully-arid basins was a division by a supply term missing a term, not a statement about geometry. Implement it INERT -- identical verdicts when no groundwater field is supplied, which is a bit-identical reduction test -- and do not regenerate the carve list from it, because the head field is uncertified until GW-3 and a carve list is loop A's input [step: carve_verdict] |

## HYD -- hydrography

1 open of 18 issued.

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
| HYD-19 | -- | -- | done, see `archive/tasks.md` |

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

1 open of 6 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| MIN-1 | -- | -- | done, see `archive/tasks.md` |
| MIN-2 | -- | -- | done, see `archive/tasks.md` |
| MIN-3 | -- | -- | done, see `archive/tasks.md` |
| MIN-4 | -- | -- | done, see `archive/tasks.md` |
| MIN-5 | -- | -- | done, see `archive/tasks.md` |
| MIN-6 | Replace the rainfall proxy for water table depth in the supergene rules with a modelled depth field. `minerals/config/downstream_prospectivity.yaml` says in its own comments that it substitutes rainfall because nothing better exists, while the mechanism it cites is depth: Reich and Vasconcelos (2015) put oxidation in the vadose zone and secondary sulfides below the water table | `hydrography/notes/groundwater-scoping.md` section 6, `docs/src/reference/economic-minerals.md` | open and UNBLOCKED: GW-1 closed 2026-08-20, so the mean depth field exists as the `groundwater` step's output. What is left is a judgement rather than a wait -- the head field is uncertified until GW-3 supplies an external test, and a prospectivity rule keyed on it inherits that. It supplies the DEPTH and not the descent RATE, so Sillitoe (2005)'s wet-end control -- erosion in balance with the rate of water table descent -- stays unreachable and stays declared [step: downstream_prospectivity]

## PHYS -- physics calibrated for the wrong world

0 open of 10 issued.

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
| PHYS-11 | -- | -- | done, see `archive/tasks.md` |
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

2 open of 7 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| SURF-1 | -- | -- | done, see `archive/tasks.md` |
| SURF-2 | -- | -- | done, see `archive/tasks.md` |
| SURF-3 | -- | -- | done, see `archive/tasks.md` |
| SURF-4 | Diatomite is placed against a STATIC lake proxy, not lake occupancy over the stellar cycle. The strandline band between the solved lake surface and the spill level is the ground a single climatology's lake COULD vacate; the design asks for the band it is measured to vacate over the 57-year component, which is the forcing slow enough for lake level to equilibrate | `pedology/notes/derived-surface-classes.md` | open, and blocked on CYC-1 rather than on effort. It cannot be faked from one climatology: a basin permanently full never exposes its bed and a basin permanently dry never accumulates one, so the whole class is about alternation. Re-run `build_surface_classes.py` against the cycle climatologies when they exist [step: surface_classes, stellar_cycle_run] |
| SURF-5 | -- | -- | done, see `archive/tasks.md` |
| SURF-6 | -- | -- | done, see `archive/tasks.md` |
| SURF-7 | Place groundwater silcrete and calcrete's groundwater-calcite pathway once a water table exists. `pedology/config/surface_classes.yaml` currently refuses groundwater silcrete by name, with the reason in the config: it sits at or near a water table and this project models none | `hydrography/notes/groundwater-scoping.md` section 6, `pedology/notes/derived-surface-classes.md` | open and UNBLOCKED: GW-1 closed 2026-08-20 and the `groundwater` step generates the field, so this is now work rather than a wait, carrying GW-3's uncertification with it. What becomes placeable is the MEAN depth field and a discharge mask, not Fenske et al. (2025)'s model: that hardens a layer over the RANGE of water table fluctuation at 10^5 years or longer, and both are durations `docs/src/reference/no-time-axis.md` refuses. Gypcrete is untouched, being surface and air processes rather than this mechanism [step: surface_classes]

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

