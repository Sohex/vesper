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

1 open of 7 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| BUDG-1 | -- | -- | done, see `archive/tasks.md` |
| BUDG-2 | -- | -- | done, see `archive/tasks.md` |
| BUDG-3 | -- | -- | done, see `archive/tasks.md` |
| BUDG-4 | -- | -- | done, see `archive/tasks.md` |
| BUDG-5 | -- | -- | done, see `archive/tasks.md` |
| BUDG-6 | -- | -- | done, see `archive/tasks.md` |
| BUDG-7 | Retire the "CH4 and N2O, absent from every run on disk" row from `scripts/error_budget.py` once a baseline carries the band, leaving only the shortfall row beside it | `scripts/error_budget.py`, `exoplasim/notes/trace-gas-band.md` | open, and it is a DOUBLE COUNT until it happens: the two rows price the same gas twice, one for what the scheme fails to apply and one for what the existing runs never had at all. They are both correct today, because the band is resident and no run has it. The moment a run does, the second is wrong and the budget over-counts by 1.56 to 1.98 W/m2, which would make the largest item in it larger still and misrank everything below. It is a task rather than a step because nothing regenerates it: `error_budget.py` reads no run manifest and cannot notice [step: baseline_run] |

## CLIM -- climate

5 open of 46 issued.

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
| CLIM-42 | -- | -- | done, see `archive/tasks.md` |
| CLIM-43 | -- | -- | done, see `archive/tasks.md` |
| CLIM-44 | -- | -- | done, see `archive/tasks.md` |
| CLIM-45 | Re-take CLIM-39's REDUCTION IDENTITY, which CLIM-44 showed was stated at a horizon that writes no gridpoint output: same configuration through a base-built binary and the multi-species one, at a segment length chosen against the write cadence and at a fixed rank count | `notes/audits/model-reproducibility.md`, `aeolian/notes/multi-species-aerosol.md` section 8 | open, and cheap -- seconds per run on the bed `exoplasim/scripts/reproducibility_matrix.py` uses, whose `check_binaries` guard also stops it running against a run directory's frozen executables. It needs a binary built from the pre-CLIM-39 base, which is the only part that is not free. UNVERIFIED rather than refuted: nothing suggests the identity fails, and the multi-species change is already resident and rebuilt, so this is closing a proof and not reopening a result [step: rebuild_binaries] |
| CLIM-46 | WITHDRAWN, and kept as a row because the id was issued and the mistake is the useful part. The claim was that land hypsometry is biased high by resolution and does not converge, from a two-point comparison at 2,500,001 and 10,000,005 regions. A control at 2,600,001 and 2,700,001, where the resolution is to all intents identical, spreads mean land elevation by 0.0523 km and land area above 1 km by 11.0%, against the 4x differences of 0.0355 km and 10.3%. The effect is smaller than the noise. Its cause is that `randInt(numRegions)` picks the first plate seed by INDEX, so changing the region count also changes the realisation | `notes/audits/orogen-resolution.md` | wontfix, no defect demonstrated. The standing rule it leaves: on this generator a resolution claim needs a same-resolution control before it is a claim, and two runs a few percent apart is the cheap version. The same control also withdrew a concavity result and an RMS-from-converged-terrain figure from the same note [step: orogen] |

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
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here) | `notes/audits/orogen-gravity.md` | blocked on the downscaling pass, deliberately. AMENDED: two further facts make a g^3 term precision on the wrong quantity. `iceFlow` accumulates a CELL TALLY, seeded from the glaciation index and summed downstream, where the hydraulic path was deliberately given `flowSeed = cellArea / refCellKm2` so its accumulation is a real drainage area; the two accumulators disagree about what is resolution-independent, and the same valley collects about four times the ice flow at four times the region count, carving as `iceFlow^0.6`. And per PHYS-13 the index feeding it is an Earth-calibrated latitude threshold blind to this planet's obliquity, spectrum and rotation. Fix the ice mask first. The sub-grid judgement also survives measurement: at four times the region count the mean edge is 7.59 km against a 1-5 km valley, and reaching 1.5 km edges needs about 256 million regions, 26x the run measured. See `notes/audits/orogen-resolution.md`. The audit's recommendation was Option B, carry it as a declared gap: `iceFlow` is a heuristic accumulation rather than a Glen's-law velocity, so a rigorous g^3 on it is precision theatre. Decisive on top of that, a glacial valley is 1-5 km against a 15.19 km mesh cell, so the process is SUB-GRID and cannot be resolved here at any gravity. It belongs with the downscaling machinery, which has to persist sub-grid hypsometry anyway [step: orogen] |

## GRID -- the Orogen/ExoPlaSim grid convention

0 open of 1 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| GRID-1 | -- | -- | done, see `archive/tasks.md` |

## GW -- groundwater

4 open of 19 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| GW-1 | -- | -- | done, see `archive/tasks.md` |
| GW-2 | -- | -- | done, see `archive/tasks.md` |
| GW-3 | -- | -- | done, see `archive/tasks.md` |
| GW-4 | -- | -- | done, see `archive/tasks.md` |
| GW-5 | Groundwater evapotranspiration reaches ExoPlaSim too: `landmod.f90` throttles evaporation on `dwatc/dwmax` and the cell dries out, so a groundwater-fed playa or oasis under-evaporates in the CLIMATE as well as in the water table | `hydrography/notes/groundwater-et-sink.md`, `notes/audits/absent-and-inherited-physics.md` | open, one-signed, and NARROWED by GW-15: the hydrography half of this sink now exists and is validated, so what is left is only the climate half. Still NOT fixable by raising `dwmax` through `surface_soil_water` -- a deeper bucket changes storage and timing, not the availability floor a water table sets. The routes are a lower boundary in `landmod.f90`, which drags CLAUDE.md rule 4 and a full rebuild, or a fork of the LPJ-GUESS soil column [step: surface_soil_water] |
| GW-6 | Valley-to-ridge water table texture is sub-grid and stays sub-grid: Fan et al. (2013) put the well-articulated gradient at decameters to kilometres, against a 15.19 km region | `hydrography/notes/groundwater-scoping.md` section 5, `hydrography/notes/groundwater-et-sink.md` | open as a DECLARED gap on GRAV-6's precedent, and now MEASURED TWICE. GW-3 scored the solver at 70,119 Australian bores with no skill at all. GW-15 then added the physics that was actually missing, removed the surface pinning entirely, and the skill did not appear: Pearson +0.031, 0.1% of variance, with the model spanning a factor of 1.9 in depth against the observations' 56.5. The model is a function of recharge at Spearman -0.977 and the real water table is not, at -0.156. So the variance is terrain at a scale this mesh cannot hold, and no sink, thickness or conductivity recovers it. NOW ANSWERED for the resolution route, `notes/audits/orogen-resolution.md`: a generation at four times the region count adds at most 12% to relief at any separation from 5 to 200 km, and the pre-erosion substrate is converged outright, so a finer Orogen does not hold the missing texture and cannot supply it. The generator's noise sits at fixed physical wavelengths and its erosion cuts at whatever the mesh scale is. Sub-grid hypsometry has to come from a parameterisation or external data. Revisit under loop D on that basis, not as a resolution problem [step: hydrography] |
| GW-7 | -- | -- | done, see `archive/tasks.md` |
| GW-11 | -- | -- | done, see `archive/tasks.md` |
| GW-8 | The discrete elliptic operator's noise floor is about 12% relative RMS against the analytic Legendre eigenvalue above l=1, distributed truncation rather than bad faces | `hydrography/notes/mesh-geometry.md`, `hydrography/notes/earth-calibration-criterion.md` | open. Declared bar was 0.10 and l=1 passes at 0.0217. CONFIRMED INDEPENDENTLY by GW-3: an Earth mesh built from Orogen's own generator at the same cell size reproduces 0.1062, 0.1186 and 0.1235 at l=2 to 4 against Vesper's 0.1082, 0.1202 and 0.1243, within 2%. So the floor is a property of this discretisation at this cell size rather than of one mesh, and it is what let GW-3 exclude the discretisation as an explanation for its miss [step: hydrography] |
| GW-9 | -- | -- | done, see `archive/tasks.md` |
| GW-10 | -- | -- | done, see `archive/tasks.md` |
| GW-12 | -- | -- | done, see `archive/tasks.md` |
| GW-13 | -- | -- | done, see `archive/tasks.md` |
| GW-14 | -- | -- | done, see `archive/tasks.md` |
| GW-15 | -- | -- | done, see `archive/tasks.md` |
| GW-16 | -- | -- | done, see `archive/tasks.md` |
| GW-17 | -- | -- | done, see `archive/tasks.md` |
| GW-18 | Make aquifer thickness spatially varying, from a sediment-thickness source rather than the constant 100 m that is Gleeson's map depth | `hydrography/notes/groundwater-et-sink.md` | open, and NARROWED by GW-17: a uniform 2 km with local baselevels already fixes the RANGE, taking the 95th percentile depth from 6.9 m to 52 m against an observed 42 m, so what a spatially varying thickness has left to buy is PATTERN. On the evidence it will not buy much: no resolvable field predicts observed depth above Spearman 0.27, and height above the nearest river -- the mechanism thickness would strengthen -- is the weakest at 0.064. Do it for the range being physically sourced rather than assumed, not expecting skill [step: hydrography] |
| GW-19 | -- | -- | done, see `archive/tasks.md` |

## HYD -- hydrography

2 open of 19 issued.

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
| HYD-20 | The endorheic basin catalogue is resolution-limited rather than physics-limited. `selectionCriteria` declares floors of 0.05 km depth, 1000 km2 area and 12 mesh cells, and at the build's spacing the 12-cell floor is 3526 km2, so nothing below 2646 km2 survives and the carve list is drawn from a set the mesh floored rather than the criterion selected. Preserved goes 3,621, 6,345, 9,419 across one, two and four times the region count, and at four times the smallest preserved depression is 1000 km2 exactly, the declared floor having finally taken over from the mesh. The large end does NOT settle: the biggest basin is 3,059,665 km2 at the build and 5,499,601 at four times, one south polar basin remeasured 80% larger rather than two merging, and it is not erosional, since catalogues with `--glacial 0` are byte-identical | `notes/audits/orogen-resolution.md` | open. Not a defect in the floor, which is a signal-to-noise criterion and correctly relaxes as the mesh becomes able to hold smaller depressions. It is a statement about what this build can see: the declared physical floor only starts binding near 8.8M regions. DECIDED 2026-08-20: the next generation raises the region count, and the reference export at 10,000,005 is registered as `precarve-craton-10m`. It is NOT the build to activate, because a generation meant to be commissioned should also carry the PHYS-13 ice-mask decision; the two changes belong in one generation, not two. Remaining work is that generation: pick the count, emit the full grid set T21/T42/T63/T85 and 512x256 rather than T42 alone, register the hash, and point `source_build` at it. Note what activating costs under rule 7: `precarve-craton` has been consumed by a commissioned run, so superseding it makes the climatology and everything derived from it worthless rather than stale [step: orogen] |

## LITH -- lithology

1 open of 26 issued.

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
| LITH-25 | `computeScarpPotential` gates on an absolute rise over run, the 0.002 to 0.010 smoothstep, but the gradient a mesh can measure depends on its spacing: 15.19 km cells average away drops that finer cells resolve, so the same escarpment reports a steeper gradient and crosses the smoothstep sooner at a higher region count. The 150 km influence band scales correctly through `SCARP_EDGE_KM / edgeKm`, 10 hops against 20; the gradient thresholds do not scale at all | `notes/audits/orogen-resolution.md` | open. `scarp_potential` is therefore not comparable across region counts, and a finer generation strengthens scarps for reasons partly real, the landform is better captured, and partly calibration, the thresholds were tuned at some other spacing. MEASURED at four times the region count: cells with non-zero `scarp_potential` go 12.456% to 14.999%, the mean where non-zero goes 0.1316 to 0.1986, a rise of 51%, and the 99th percentile land slope goes 5.015 to 8.427 degrees, while relief at fixed separation rose at most 12%. The gap between 51% and 12% is the calibration. Separating the two needs the thresholds restated against a known escarpment gradient rather than left at their inherited values [step: orogen] |

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

1 open of 11 issued.

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
| PHYS-13 | Orogen places the ice that carves this terrain by latitude and dimensionless elevation alone. At the build's `glacialErosion` of 0.8 the ice line falls at 58 degrees latitude, positioned by the erosion slider and nothing else, and the altitude gate is in the model's elevation parameter so the snowline tracks maximum relief, which scales as 1/g, rather than a freezing altitude set by lapse rate and surface temperature. `glacIdx` consults no temperature at all: not the K2.5V spectrum, not the 32 degree obliquity, not the 30 hour rotation, and not even Orogen's own temperature field, which `computeTemperature` produces after `runPostProcessing` has already eroded. All three leading parameters point the same way, toward less polar ice than an Earth-tuned 58 degree threshold gives | `notes/audits/orogen-resolution.md` | open. The correction needs no new physics: the ExoPlaSim climatology on this build already carries `glac`, `snd`, `sic`, `snm` and `alb`, which is this planet's ice computed under its actual obliquity, orbit and spectrum. Drive `glacIdx` from that field instead of `polarDist` and close it through loop A, which exists to feed a commissioned climate into the next generation; the terrain-before-climate ordering that blocks it inside one pass is not a constraint across iterations. SHAPE OF THE FIX, decided 2026-08-20: an ice mask consumed AT GENERATION like the carve list, `--ice-mask FILE`, not a post-hoc terrain alteration. Glacial, hydraulic and thermal erosion share one iteration loop and a mid-loop priority flood at `GLACIAL_MID_FLOOD_FRAC` exists to cut outlets through the depressions glaciation makes, so a separate later pass would leave them undrained and move the drainage network and the basin catalogue. Bootstrap with `--glacial 0`, an honest null, then commission and regenerate with the mask. Two constraints on it: mapping `glac` from the ExoPlaSim grid back to the mesh is a rule 3 crossing and must go by index or a shared coordinate source, and per `docs/src/reference/no-time-axis.md` a current ice field can only say WHERE ice was, never how long, so the strength slider stays a declared choice rather than a measurement. Subsumes the mechanism half of GRAV-6 [step: orogen] |

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

