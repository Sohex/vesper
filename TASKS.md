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

25 open of 28 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| BIO-1 | -- | -- | done, see `archive/tasks.md` |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; recompute the bracket after BIO-24 fixes the Cleveland fit's time basis rather than carrying forward the smoke-run span as though five model orbits were five Earth years | `biosphere/notes/productivity-prediction.md`, `biosphere/notes/implicit-earth-assumptions.md` finding 4 | blocked on BIO-24 and iteration 2's baseline run -- no accepted LPJ-GUESS run exists on this build [step: lpj_run] |
| BIO-3 | -- | -- | done, see `archive/tasks.md` |
| BIO-4 | -- | -- | done, see `archive/tasks.md` |
| BIO-5 | Derive the Vesper parameterisation for LPJ-GUESS-CNP's three soil-P inputs: convert the existing rock-class P release and per-cell weathering into `pwtr` in the explicit absolute-time unit settled by BIO-24, and ground `kplab` and `spmax` in the pedology fields that control sorption | `biosphere/notes/cnp-fork-scoping.md`, `biosphere/notes/implicit-earth-assumptions.md` finding 4 | open; this is the scientific mapping the current site defaults stand in for, while its flux-unit contract is blocked on BIO-24 [step: soil] |
| BIO-6 | Emit `pwtr`, `kplab` and `spmax` as named, provenance-stamped columns in the per-cell soil map, with validation that every land cell receives finite non-negative values | `biosphere/notes/cnp-fork-scoping.md` | open, blocked on BIO-5 settling what the columns mean [step: soil] |
| BIO-7 | Extend the vendored `SoilInput` mineral-texture path to read the three soil-P columns, remove the temporary site defaults, and fail if phosphorus limitation is enabled without complete inputs | `biosphere/notes/cnp-fork-scoping.md` | open, blocked on BIO-6 defining the file contract [step: lpj_run] |
| BIO-8 | Declare a Vesper phosphorus-deposition rate and carry it through the LPJ driver and `vesperinput` in kgP/m2/day; do not substitute the fork's Earth observation file | `biosphere/notes/cnp-fork-scoping.md` | open; independent of BIO-5 through BIO-7 but required before activation [step: lpj_driver, lpj_run] |
| BIO-9 | Register a dated C-N-P productivity prediction, uncertainty bracket and scoring rule before the first phosphorus-limited run; the existing C-N prediction must remain unchanged | `biosphere/notes/cnp-fork-scoping.md` | open, blocked on BIO-5, BIO-8 and BIO-24 settling the model inputs and their time basis [step: lpj_run] |
| BIO-10 | Enable `ifplim`, retain the P-pool/source/flux diagnostics needed to close the budget, verify that `ifplim 0` preserves the C-N baseline, and require P mass-balance checks before a full run is interpreted | `biosphere/notes/cnp-fork-scoping.md` | open, blocked on BIO-6 through BIO-9 and BIO-14's general run-acceptance contract [step: lpj_run] |
| BIO-11 | Define and emit one provenance-stamped rootable surface fraction from the derived surface classes and solved lakes, then use it to weight LPJ extensive outputs, soil-carbon feedback and productivity denominators so salt, playa and open water do not produce biomass | `biosphere/notes/modelling-gap-audit.md` finding 1 | open; prerequisite to interpreting either the C-N or C-N-P result, and to BIO-16 and BIO-17 consuming the same surface [step: lpj_driver, lpj_run, soil] |
| BIO-12 | Replace every last-output-year read with one shared equilibrium-window reducer: reject a trending end window, average a fixed number of complete forcing cycles, and report temporal plus patch/seed uncertainty for feedback, scoring and soil carbon | `biosphere/notes/modelling-gap-audit.md` finding 2 | open; prerequisite to downstream use of an LPJ run and to BIO-15 through BIO-17 [step: lpj_run, soil, surface_albedo] |
| BIO-13 | Carry Vesper-native rainfall event information in the driver -- a daily sequence or monthly wet-day counts -- and make `VesperInput` honour `ifrainonwetdaysonly` while conserving every monthly total; measure smooth versus event-distributed precipitation before adoption and do not import Earth GWGEN observations | `biosphere/notes/modelling-gap-audit.md` finding 3 | open; the present adapter silently supplies smooth daily rain despite the wet-days declaration [step: lpj_driver, lpj_run] |
| BIO-14 | Make LPJ completion an assessed artifact: require every rank, output, cell and year, reject non-finite or invalid values, retain the stock runoff and N pool/flux diagnostics, and report C, N and water closure before albedo, soil or scoring may consume a run | `biosphere/notes/modelling-gap-audit.md` finding 4 | open; independent of phosphorus and the base acceptance layer BIO-10 must extend [step: lpj_run] |
| BIO-15 | Implement the soil-biosphere convergence assessor declared in `pedogenesis.yaml`: compare successive compatible iterations, apply both fixed tolerances, write the residuals and fail at the maximum iteration rather than leaving convergence as a manual judgement | `biosphere/notes/modelling-gap-audit.md` finding 5 | open, blocked on BIO-12's equilibrium statistic and BIO-14's accepted-run artifact [step: soil, lpj_run] |
| BIO-16 | Close the aerodynamic feedback by deriving code 173 from the same aggregated, rootable tree/grass cover as modelled albedo, retaining lake and orographic terms, and recording the `lpj_run` back-edge in the pipeline | `biosphere/notes/modelling-gap-audit.md` finding 6 | open, blocked on BIO-11 and BIO-12 defining the cover being consumed [step: surface_roughness, lpj_run] |
| BIO-17 | Correct modelled albedo and forest compositing so persistent lakes remain water with zero forest and partial barren/lake cells use BIO-11's rootable fraction; cover ordinary, barren and partial-lake cells with a regression fixture | `biosphere/notes/modelling-gap-audit.md` finding 7 | open, blocked on BIO-11 and BIO-12; the current modelled blend runs after lakes and repaints them [step: surface_albedo] |
| BIO-18 | Derive K-star broadband and two-band albedos for the modelled tree and grass endmembers, write distinct codes 175 and 176, and report their stellar-flux-weighted recombination residual instead of copying the broadband field into all three codes | `biosphere/notes/modelling-gap-audit.md` finding 8, archived SPEC-5 residual | open; the scientific derivation is independent of an LPJ run and should precede adoption of modelled albedo [step: vegetation_albedo, surface_albedo] |
| BIO-19 | Make bootstrap and iterative baseline soil states unambiguous in the pipeline, and reject any soil map whose recorded climatology hash does not match the climatology used by the LPJ driver | `biosphere/notes/modelling-gap-audit.md` finding 9 | open; the current graph names bootstrap climate while the builder defaults to baseline climate [step: soil, lpj_driver] |
| BIO-20 | Add a no-simulation regression suite for the Vesper adapter and consumers: binary round-trip, coordinate bijection, year cycling, precipitation conservation, incomplete-output rejection, equilibrium reduction, rootable/lake/two-band compositing, and the seasonal-phase, time-base, pressure/radiation, nutrient-unit and photon-supply contracts added by BIO-21 through BIO-25 | `biosphere/notes/modelling-gap-audit.md` finding 10, `biosphere/notes/implicit-earth-assumptions.md` | open; keep the interface testable without building or launching LPJ-GUESS, adding fixtures alongside BIO-11 through BIO-25 as their contracts land [step: lpj_driver, lpj_run, surface_albedo, surface_roughness] |
| BIO-21 | Replace the natural-vegetation January/July and day-14/day-195 phenology constants with forcing-derived seasonal landmarks, including GDD, chilling, leaf-on and summergreen litter resets; cover both hemispheres, the equator and multi-orbit forcing, and prove `chilldays` cannot exceed its lookup table | `biosphere/notes/implicit-earth-assumptions.md` finding 1 | open; correctness blocker -- day 195 is unreachable in the 181-day calendar and the southern annual resets never occur [step: lpj_run] |
| BIO-22 | Classify every non-fire LPJ ecological quantity expressed per year as absolute-time, seasonal-cycle, accumulated-flux or diagnostic, then convert PFT turnover/longevity/mortality/establishment windows, litter and SOM decay, and CNP P-pool kinetics according to that registry with no-simulation regression checks; do not apply one blanket multiplier to fractions | `biosphere/notes/implicit-earth-assumptions.md` finding 2 | open; correctness blocker -- absolute-time rates presently execute once per 181-day orbit and can run about 2.02x fast, while genuinely seasonal events must remain once per orbit [step: lpj_run] |
| BIO-23 | Extend the Vesper driver with ExoPlaSim surface pressure, net longwave radiation, humidity and wind; use local pressure for CO2/O2 partial pressures and the psychrometric term, and compare the corrected LPJ EET against the existing Vesper Penman calculation before choosing whether to change model form | `biosphere/notes/implicit-earth-assumptions.md` finding 3 | open; the source climatology already carries `ps`, `rls`, `hus` and `spd`, but LPJ currently discards them, reconstructs Earth longwave and gives every elevation sea-level photosynthesis [step: lpj_driver, lpj_run] |
| BIO-24 | Define an explicit absolute-day/Earth-year/orbit unit contract for every N and P input and diagnostic; correct N deposition, the Cleveland fixation slope/intercept and averaging window, and texture-path `pwtr`, then update BIO-2 and BIO-5 through BIO-10 to consume the same convention | `biosphere/notes/implicit-earth-assumptions.md` finding 4 | open; prerequisite to interpreting the C-N baseline or defining BIO-5 -- the declared N deposition is currently delivered once per orbit and annual P weathering would inherit the same approximately 2.02x error [step: soil, lpj_driver, lpj_run] |
| BIO-25 | Generate a stellar-spectrum- and photosystem-window-weighted photon conversion constant beside `VESPER_FRADPAR`, use it in LPJ photosynthesis, and require the combined `FRADPAR * CQ` photon supply to reproduce the registered productivity calculation | `biosphere/notes/implicit-earth-assumptions.md` finding 5 | open; LPJ currently converts the K-star's PAR energy with a fixed 550 nm Earth coefficient, so the executable and prediction use different light currencies [step: lpj_run] |
| BIO-26 | Register and run a one-factor sensitivity for the active Earth-global `ALPHAA_NLIM` canopy scalar, without retuning it to the desired Vesper productivity, and propagate the response as model-form uncertainty in productivity and feedback products | `biosphere/notes/implicit-earth-assumptions.md` finding 6 | open; measurement task, blocked on an accepted baseline under BIO-12 and BIO-14 and the correctness fixes BIO-21 through BIO-25 [step: lpj_run, surface_albedo, soil] |
| BIO-27 | Make cropland fail closed until its sowing, harvest and test dates no longer use 180/364/365-day Earth ordinals; express the calendar contract in Vesper seasonal coordinates and cover both hemispheres before permitting non-zero cropland | `biosphere/notes/implicit-earth-assumptions.md` finding 7 | open; dormant activation guard and not a natural-vegetation baseline blocker [step: lpj_run] |
| BIO-28 | Make methane fail closed until the module consumes Vesper gravity, local surface pressure and declared atmospheric O2/CH4, and record a decision or sensitivity for its Earth-global emissions tuning before `ifmethane` can be enabled | `biosphere/notes/implicit-earth-assumptions.md` finding 7 | open; dormant activation guard while `ifmethane 0`, independent of the baseline C-N/C-N-P run [step: lpj_run] |

## BUDG -- the error budget and what it is denominated in

0 open of 7 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| BUDG-1 | -- | -- | done, see `archive/tasks.md` |
| BUDG-2 | -- | -- | done, see `archive/tasks.md` |
| BUDG-3 | -- | -- | done, see `archive/tasks.md` |
| BUDG-4 | -- | -- | done, see `archive/tasks.md` |
| BUDG-5 | -- | -- | done, see `archive/tasks.md` |
| BUDG-6 | -- | -- | done, see `archive/tasks.md` |
| BUDG-7 | -- | -- | wontfix, see `archive/tasks.md` |

## CLIM -- climate

5 open of 48 issued.

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
| CLIM-45 | -- | -- | done, see `archive/tasks.md` |
| CLIM-46 | WITHDRAWN, and kept as a row because the id was issued and the mistake is the useful part. The claim was that land hypsometry is biased high by resolution and does not converge, from a two-point comparison at 2,500,001 and 10,000,005 regions. A control at 2,600,001 and 2,700,001, where the resolution is to all intents identical, spreads mean land elevation by 0.0523 km and land area above 1 km by 11.0%, against the 4x differences of 0.0355 km and 10.3%. The effect is smaller than the noise. Its cause is that `randInt(numRegions)` picks the first plate seed by INDEX, so changing the region count also changes the realisation | `notes/audits/orogen-resolution.md` | wontfix, no defect demonstrated. The standing rule it leaves: on this generator a resolution claim needs a same-resolution control before it is a claim, and two runs a few percent apart is the cheap version. The same control also withdrew a concavity result and an RMS-from-converged-terrain figure from the same note [step: orogen] |
| CLIM-47 | -- | -- | wontfix, see `archive/tasks.md` |
| CLIM-48 | Store the Legendre weights as P and Q alone. `legini` builds EIGHT `NCSP x NLPP` matrices and every one is touched every timestep, which is 48.4 MB a die at T127 and 114.9 MB at T170 against CCD0's 96 MB of L3 and CCD1's 32 MB -- and that footprint is the measured cause of a 10.5 point load imbalance that makes CCD1 the critical path at both resolutions. All eight are P or Q times a per-MODE factor and a per-LATITUDE factor, both separable and neither needing an `NCSP x NLPP` array, so storing two matrices instead of eight leaves 12.1 MB and 28.7 MB a die and both fit CCD1. It is the exact inverse of the filter fold, which moved a per-mode scalar INTO the matrices to save one multiply in three and returned 1.6% because those loops are bound by streaming the matrices rather than by the multiplier; the same sentence read backwards is the argument for this | `exoplasim/notes/rank-imbalance-and-weight-traffic.md`, `exoplasim/notes/legendre-filter-fold.md` | open. It is a reassociation, so restart shas move and it needs `compare_restarts.py` at a declared tolerance with a control, as the symmetry work did. Two things to establish before believing the prize: whether the extra per-mode vector read costs the inner loops their vectorisation, and whether anything survives at T21 and T42 where the matrices already fit [step: rebuild_binaries, baseline_run] |

## CONS -- consistency checking

0 open of 12 issued.

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
| CONS-12 | WITHDRAWN, raised in error, and kept because the misreading is the reusable part. The claim was that `world_state.json` carried a wrong `endorheic_fraction_of_land`, 0.7630 against an artifact value of 0.2676. Both numbers are right and they are different quantities: `is_endorheic` marks cells INSIDE a preserved basin, which is 0.2676 of land area, while the generator reads `basins.drainageConsistency.fractionOfLand`, which is `totalByFieldKm2 / landAreaKm2` and is the land DRAINING INTO those basins, their catchments, at 0.7630. Catchments are far larger than the depressions they feed. Regenerating with `scripts/world_state.py` reproduces 0.7630 exactly, as it should | `world_state.json`, `docs/src/reference/builds.md` | wontfix, no defect. Two things worth carrying out of it. `world_state.json` is GENERATED and `docs/src/reference/builds.md` says so in as many words, never edited and re-run after anything that changes a value, so the correct response to a suspect entry is to re-run the generator and read the manifest field it names, not to edit the file. And the key name is the trap: `endorheic_fraction_of_land` reads as basin membership and means catchment share, which is what caused this. It is not worth a schema change on its own, but if the key is ever touched it should say catchment [step: any] |

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

## FIRE -- fire modelling

9 open of 9 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| FIRE-1 | Add a lightweight ExoPlaSim lightning-potential diagnostic from pre-convective or explicitly upstream CAPE and synchronous convective precipitation, with defined substep accumulation, units and provenance; do not multiply the existing post-`rainstep` code 322 by code 143, and benchmark diagnostic overhead before production activation | `biosphere/notes/fire-model-audit.md` finding 3 | open; the source change and no-simulation checks can proceed independently, but any timing or climate run requires explicit permission [step: baseline_run] |
| FIRE-2 | Preserve one coherent chronological fire-weather sequence through climatology and the LPJ driver -- lightning potential, convective and total precipitation, temperature, humidity, wind, pressure and radiation -- define the fire step against Vesper's 30-hour solar day and 181-day orbit, and conserve integrated fluxes rather than independently smoothing the fields | `biosphere/notes/fire-model-audit.md` finding 5 | open, builds on BIO-13's rainfall events and BIO-23's missing weather fields; do not infer event order from monthly means [step: baseline_climatology, lpj_driver] |
| FIRE-3 | Convert lightning potential to successful natural ignitions through explicit all-flash-to-ground and ignition-efficiency priors, a synchronous wet-lightning filter and uncertainty brackets; human ignition and suppression are absent, while ubiquitous ignition remains a labelled model-form control rather than the central case | `biosphere/notes/fire-model-audit.md` finding 4 | open, blocked on FIRE-1 and FIRE-2 defining the lightning and weather contract [step: lpj_driver, lpj_run] |
| FIRE-4 | Replace SIMFIRE's archive and `cru_ncep` generator-name gate with an explicit fire-driver capability interface that accepts Vesper-native ignition and weather fields and fails closed when any required capability is missing; do not ingest GFED burned-area seasonality or HYDE population | `biosphere/notes/fire-model-audit.md` findings 1 and 2 | open; independent interface work, required before the reduced occurrence model can reach BLAZE [step: lpj_driver, lpj_run] |
| FIRE-5 | Implement a reduced-complexity occurrence-to-area layer at the FIRE-4 seam: INFERNO-like fuel/moisture flammability and a Li-like explicit `fire count * bounded area per fire`, with wind-dependent spread as a declared option and every Earth fuel, humidity, duration, PFT-spread and fire-size coefficient exposed as a prior or bracket rather than a Vesper calibration | `biosphere/notes/fire-model-audit.md` finding 2 | open, blocked on FIRE-2 through FIRE-4; SPITFIRE complexity is deliberately deferred unless this model-form bracket proves inadequate [step: lpj_run] |
| FIRE-6 | Audit and port the retained BLAZE effects layer: replace latitude/biome tuning and 365-day assumptions with Vesper seasonal and absolute-time contracts, register combustion and mortality coefficients with provenance and brackets, and cover fuel consumption, intensity and cohort mortality with no-simulation fixtures | `biosphere/notes/fire-model-audit.md` finding 6 | open; can proceed independently of the ignition model but must close before BLAZE activation [step: lpj_run] |
| FIRE-7 | Implement phosphorus-conserving fire effects for the CNP fork across live, transitional-litter and soil-litter pools, explicitly partitioning combusted P between atmospheric loss and retained ash/mineral soil; add `P_FIRE` diagnostics, stoichiometric invariants and a fail-closed guard so neither BLAZE nor GlobFIRM fire can run with P limitation while P routing is absent | `biosphere/notes/fire-model-audit.md` finding 6 | open; correctness blocker for any C-N-P fire run and an extension of BIO-10's P mass-balance gate [step: lpj_run] |
| FIRE-8 | Make fire an assessed artifact: retain ignition potential, ground flashes, successful starts, flammability, area per fire, burned fraction, return interval, fire-line intensity, mortality and C-N-P destinations, reject impossible ranges, and close C-N-P mass before productivity, soil, albedo or smoke consumes the result | `biosphere/notes/fire-model-audit.md` finding 7 | open, blocked on FIRE-3 through FIRE-7 and extends BIO-14's general accepted-run contract [step: lpj_run, soil, surface_albedo] |
| FIRE-9 | Pre-register and execute a matched FireMIP-style model-form comparison under identical accepted forcing: no fire, retained GlobFIRM, ubiquitous ignition, and lightning-driven central plus bracket cases; score fire mechanisms together with vegetation and hydrology and propagate the spread into productivity, albedo, soil and CLIM-29 smoke rather than selecting the case that best fits a desired outcome | `biosphere/notes/fire-model-audit.md` finding 7 | open, blocked on FIRE-8 and an accepted baseline under BIO-12 and BIO-14; execution requires explicit permission because these are LPJ-GUESS runs [step: lpj_run, soil, surface_albedo] |

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

1 open of 2 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| GRID-1 | -- | -- | done, see `archive/tasks.md` |
| GRID-2 | Sub-grid hypsometry is computed three times and persisted nowhere, and four threads are parked behind the machinery that would hold it. The quantity is the distribution of MESH elevations inside each model grid cell, 15.19 km regions under a T42 cell about 300 km across. `notes/glacier-rough-pass.md` integrated it inline and threw it away, `analysis/ice_mask_freezing_height.py` rebuilds it, and `maps/build_basemap.py` applies the same lapse correction against the render raster instead, so three implementations of one criterion sit on a field no artifact carries. It is measurable today and it is large: land peak excess `orog_max - orog_mean` averages 1.065 km and reaches 2.969 km at the 90th percentile, so the model evaluates its own high ground about 7.8 K too warm on average and 21.9 K in the top tenth, and land below freezing in the warmest month goes from 0.002% on the model's orography to 1.657% on the 15.19 km mesh. Those are `analysis/ice_mask_freezing_height.py`'s numbers, restated from the climatology product the pipeline consumes; an earlier draft of this row carried 1.355 km, 9.9 K and 1.395% from a raw run file. Waiting on it: PHYS-13's ice mask, the glacier mass balance that note says it cannot have without it, the saturated-area FRACTION that would replace SURF-7's all-or-nothing at-surface flag, and GRAV-6 | `notes/glacier-rough-pass.md`, `notes/audits/orogen-gravity.md`, `hydrography/notes/groundwater-et-sink.md`, `docs/src/reference/economic-minerals.md` | open, and what is missing is the GENERATOR rather than a stale artifact: no step in `config/pipeline.yaml` writes this field, so under the graph convention it does not exist, and the row that declares it belongs in the same commit as the script. SHAPE: one artifact per (build, grid), on the precedent of the coupling matrices the `hydrography` step already writes, which are the interface every mesh-to-grid consumer integrates over, with `lib/gridding.py` owning the region-to-cell mapping as it already does. Cheap, minutes, and no climate run: the export and a grid definition are the whole input. Rule 3 binds -- map by index or a shared coordinate source, never by longitude -- and land comes from `surface_class` per rule 1, or the closed-basin floors below sea level drop out of the very distribution this exists to hold. SCOPE, and one artifact must not be read as both halves. This delivers the MESH-UNDER-CELL half, which the export holds. The BELOW-MESH half is a different quantity and is not deliverable at any region count: `notes/audits/orogen-resolution.md` finds the generator designs nothing below about 20 km, so relief inside a mesh region has to be parameterised or come from external data, and GW-23's Earth measurement is its only calibration -- 14.8 m of within-cell standard deviation and 25.9 m from mean to minimum at 15.19 km, against a three-number target of 3.96 to 36.18 m of model spread with the observations at 28.19 m. That half stays GW-6's. Nor is this a decision to build the downscaling pass: discrete mineral deposits and glacial overdeepening are parked there too, and this row is only the one input all of them named [step: hydrography] |

## GW -- groundwater

5 open of 25 issued.

| id | task | source | status |
| --- | --- | --- | --- |
| GW-1 | -- | -- | done, see `archive/tasks.md` |
| GW-2 | -- | -- | done, see `archive/tasks.md` |
| GW-3 | -- | -- | done, see `archive/tasks.md` |
| GW-4 | -- | -- | done, see `archive/tasks.md` |
| GW-5 | Groundwater evapotranspiration reaches ExoPlaSim too: `landmod.f90` throttles evaporation on `dwatc/dwmax` and the cell dries out, so a groundwater-fed playa or oasis under-evaporates in the CLIMATE as well as in the water table | `hydrography/notes/groundwater-et-sink.md`, `notes/audits/absent-and-inherited-physics.md` | open, one-signed, and NARROWED by GW-15: the hydrography half of this sink now exists and is validated, so what is left is only the climate half. Still NOT fixable by raising `dwmax` through `surface_soil_water` -- a deeper bucket changes storage and timing, not the availability floor a water table sets. The routes are a lower boundary in `landmod.f90`, which drags CLAUDE.md rule 4 and a full rebuild, or a fork of the LPJ-GUESS soil column [step: surface_soil_water] |
| GW-6 | Valley-to-ridge water table texture is sub-grid and stays sub-grid: Fan et al. (2013) put the well-articulated gradient at decameters to kilometres, against a 15.19 km region | `hydrography/notes/groundwater-scoping.md` section 5, `hydrography/notes/groundwater-et-sink.md` | open as a DECLARED gap on GRAV-6's precedent, and now MEASURED TWICE. GW-3 scored the solver at 70,119 Australian bores with no skill at all. GW-15 then added the physics that was actually missing, removed the surface pinning entirely, and the skill did not appear: Pearson +0.031, 0.1% of variance, with the model spanning a factor of 1.9 in depth against the observations' 56.5. The model is a function of recharge at Spearman -0.977 and the real water table is not, at -0.156. So the variance is terrain at a scale this mesh cannot hold, and no sink, thickness or conductivity recovers it. NOW ANSWERED for the resolution route, `notes/audits/orogen-resolution.md`: a generation at four times the region count adds at most 12% to relief at any separation from 5 to 200 km, and the pre-erosion substrate is converged outright, so a finer Orogen does not hold the missing texture and cannot supply it. The generator's noise sits at fixed physical wavelengths and its erosion cuts at whatever the mesh scale is. Sub-grid hypsometry has to come from a parameterisation or external data. CONFIRMED INDEPENDENTLY ON EARTH 2026-08-20, which is the stronger test because Earth's topography has real content at every scale where Vesper's generator has none. The same case at 15.19, 10.74 and 7.60 km raises the CEILING from 0.6335 to 0.7408 and more than doubles the model's own spread, 1.91 to 4.33 m, while Pearson does not move at all: +0.0703, +0.0665, +0.0711. The dependence structure says why. The model is a monotone function of recharge at every cell size, Spearman -0.90 to -0.93, where the truth is -0.27; and against elevation the two diverge as the mesh refines, the observations' terrain signal strengthening +0.125 to +0.161 to +0.199 while the model's decays to -0.003. So GW-6 is NOT a resolution gap and was never going to be closed by a finer mesh on either planet. GW-21 then took the formulation apart and found the limit is neither: on RANK the model reaches Spearman +0.28, which is precisely its recharge field's own +0.27, so it recovers what its forcing carried and adds nothing, and the dependence structure can be tuned onto the observed values by transmissivity alone without agreement following. The binding constraint is the input fields IN THAT REGION, and GW-22 then showed the qualifier matters: scored on United States bores labelled unconfined the same model reaches Pearson +0.2599 against Australia's +0.0703, with a spread of 20.72 m against 28.94 where Australia gave 1.91 against 18.36. So this gap is regional rather than universal. Where the evapotranspiration sink takes the recharge locally the depth is a recharge map and the sub-grid texture this row is about is the least of its problems; where lateral flow carries it the solve does real work. `water_table.nc` ships `sink_fraction` to tell the two apart per cell. Which regime VESPER is in is not yet established, because the current forcing is a bootstrap climatology on pre-carve terrain. SHARPENED 2026-08-20: the 0.28 those rows measured against was never computed by anything, and is corrected to 0.2098, what a gradient boosting fit over every resolvable field reaches held out with the model EXCLUDED, from `earth_calibration.py --stage benchmark`. Both halves hold at once: the fields support only 0.2098 against a between-cell ceiling of 0.6335, and the model reaches +0.0193 alone and ADDS to that fit, 0.2098 to 0.2381. An earlier draft of this row called it information-negative from the -0.0467 and 0.1435 against 0.1411 measured on bores that fail the consistency filter; on the clean subset it is not, and that reading is withdrawn. See GW-21 [step: hydrography] |
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
| GW-20 | The project's ONLY external test cannot currently be re-run. GW-3 scored the solver at 70,119 Australian bores and its inputs survive as registered work, `hydrography/scripts/build_earth_wtd_sites.py` plus the fetched ETOPO, GLHYMPS and recharge products, but the harness that turned them into a result does not: the Earth mesh construction, the DEM and permeability sampling onto it, the solve and the scoring were written inline and are gone. So the one number that judges this component, R^2 = 0.017 against a bar of 0.07 and a ceiling of 0.857, is recorded but not reproducible, and cannot be recomputed at any other cell size | `hydrography/notes/earth-calibration-criterion.md`, `notes/audits/orogen-resolution.md` | done, rebuilt 2026-08-20 as `hydrography/scripts/earth_calibration.py`, registered under `one_offs` and named in `hydrography/README.md`, with `--edge-km`/`--regions` setting the mesh and per-stage caching under `data/earth_validation_cache/edge<E>km/`. The code was recovered from the session transcript rather than rewritten, and the mesh, sampled fields, permeability, solve and score were found intact in ANOTHER session's scratchpad on /tmp and preserved. REBUILDING IT CAUGHT A STALE RESULT: the recovered harness is the ORIGINAL GW-3 configuration, predating GW-15's sink and GW-17's baselevels, and it pins 93.8% of cells at the surface with a median depth of 0.00 m, which is the pathology the pipeline gate names. So the archived R^2 was measured on a configuration the model no longer uses. With the current physics restored the baseline behaves, 0.0% pinned and a 5.12 m median depth. The immediate use is the confound that note DECLARES in its own words, that a miss confounds formulation with implementation: Fan's figures are at about 1 km and this ran at 15.19 km, 15 times coarser in the direction her own finding says dominates. The fetched ETOPO is 1 arc-minute, 1.855 km, so the same case can be scored at 7.6 and 4 km without new data and without ExoPlaSim, the recharge being an external product. Either outcome closes something: skill rising with cell size makes resolution the binding constraint, and since Orogen's designed terrain has no content below about 20 km that also proves Vesper can never reach it by generating finer; skill flat at 4 km indicts the `T = K D` formulation and puts GW-9's abandonment of Fan's exponential decay back on the table [step: hydrography] |
| GW-21 | The water table's dependence structure is wrong, and it is not resolution and not the arithmetic. Measured on Earth at three cell sizes: the model is a monotone function of recharge, Spearman -0.90 to -0.93 where the truth is -0.27, and carries no elevation signal at all, -0.047 to -0.003, while the observations' strengthens as cells shrink, +0.125 to +0.199. The suspect is the sink rather than the transmissivity: with `E(h) = et_max exp(-(z-h)/lambda)` the supply-and-sink balance has the closed form `d = lambda ln(et_max A / supply)`, which `et_balance_depth` computes directly and which mentions no neighbour, so at `lambda` = 1 m it is a local clamp that can express a few metres of range from recharge alone and leaves nothing for the lateral term to organise | `hydrography/notes/earth-calibration-criterion.md`, `hydrography/notes/water-table-convergence.md` | done 2026-08-20, and the answer is neither of the two things this row first blamed. CORRECTION, and it is this row's own: an earlier draft proposed re-running with Fan's depth-decaying transmissivity and said GW-9's argument against it "was never measured". Both halves were wrong. The argument IS quantitative and sits in `groundwater.py:transmissivity`: `exp(h/f)` is convex, so a cell mean over a water table varying within the cell by `sigma` carries `exp(sigma^2/2f^2)`, and 100 m of sub-grid relief against the 0.95 m `f` Fan reaches on steep bedrock makes that `exp(5000)` -- not a correction to apply but a statement that the parameterisation has no value at this cell size, and refining to 7.60 km does not rescue it. And the proposed experiment is the one that already failed: `water-table-convergence.md` records that Picard on the transmissivity LIMIT-CYCLES, because `T` moves by a factor of e per e-folding length so the map is not a contraction, and the damping that stopped the oscillation held the water balance residual flat instead of converging. Two solvers died there. Do not re-run it. MEASURED INSTEAD. The sink does set depth locally: median depth is linear in `lambda`, 5.12/10.23/25.51/50.86/126.35 m at 1/2/5/10/25 m, and rho(model, recharge) is pinned at -0.90 throughout. But raising `D` from 100 m to 20 km moves rho(model, elevation) -0.047 to +0.466 and rho(model, recharge) -0.902 to -0.245, bracketing the observed +0.125 and -0.272, so the dependence structure IS reachable with the formulation in hand, and agreement still does not follow: Pearson stays between +0.058 and +0.101 across `D` over 200-fold, `lambda` over 25-fold, and the sink on or off. SECOND CORRECTION: this row and GW-6 quoted Pearson and said no skill. On RANK the model has Spearman +0.276, +0.291 at best, which is not nothing and the claim is withdrawn. What it is, though, is exactly the recharge field's own rank agreement, -0.272, so the flow solve recovers what its forcing already carried and adds under 0.02. R^2 stays negative and worsens with `D`, -0.365 to -7.025, because the spread overshoots to 45 and 75 m against the observed 28.5. So the limit is the INPUTS, not the formulation and not the mesh: recharge, GLHYMPS at 15 km and a cell-mean elevation do not determine this water table, and no rewriting of a flow equation adds information its forcing lacks. What would reopen it is a better forcing, not a better solver: a recharge product that is not itself a climatology-derived proxy, or a permeability field with real aquifer structure rather than surface lithology. The saturated-thickness form `T = K (h - z_bottom)` remains the one untried formulation worth a look, being linear in head and so free of both the Jensen term and the factor of e that killed the exponential, but on this evidence it would be tidying rather than a fix [step: hydrography] |
| GW-22 | The Earth calibration has only ever been run on Australia, and `docs/src/reference/external-data.md` says take Australia FIRST, not only. The second leg was never run and was never tracked. It matters more now than when that was written: GW-21 concluded the limit is the INPUT FIELDS rather than the formulation or the mesh, and that conclusion rests on one region; and the confound flagged against it, that Australia contains the Great Artesian Basin and a bore screened below the water table measures a potentiometric head in a confined aquifer, is Australia-specific | `docs/src/reference/external-data.md`, `hydrography/notes/earth-calibration-criterion.md` | done 2026-08-20, and it OVERTURNS the Australian verdict rather than confirming it. On 73,451 bores the USGS labels unconfined -- a water table, which is what the solver computes, rather than the potentiometric head a confined bore reads -- the cell-mean correlation is +0.2599 against Australia's +0.0703, and the model's spread is 20.72 m against an observed 28.94 where in Australia it was 1.91 against 18.36. Fitting a two-parameter bias correction on half the cells and scoring the other half over twenty splits gives held-out R2 = 0.0646 with a standard deviation of 0.0216, clearing the declared 0.07 bar on 9 of them: it STRADDLES the bar rather than passing it, and one favourable split alone read +0.0748. Repeated on the wide domain, modelled from 15 to 60 north and 130 to 60 west so neither land border is a fake coast, the answer moves in the third decimal, which retires that worry and says `score_bbox` was already doing the work. THE FINDING IS A REGIME. Where the sink takes the recharge locally the steady state collapses to `lambda ln(et_max A / supply)` and the flow solve contributes no variance; where lateral flow carries it the solve does real work, which is why the same model behaves differently on two continents. That is measurable per cell, so `water_table.nc` now ships `sink_fraction` beside the depth and a consumer can ask which of three cases a cell is in: a recharge map near 1, a pinned and seeping cell when `at_surface` is also set, or a flow solution. So the depth field's restriction is a regime statement rather than a blanket failure. What this does NOT establish is which regime Vesper is in: the sink takes a median 97.0% of recharge on this build, but that is a bootstrap climatology on pre-carve terrain, and applying a carve verdict once took endorheic land from 60.10% to 43.06%, so the premise is a property of the open loop. Re-read `sink_fraction` after a baseline. Canada stays unfetched and stays lowest priority for the reason already recorded. Superseded plan: cheaper than it looks, because the expensive parts are already built and region-independent: the Fibonacci mesh and its `Geometry` are GLOBAL and cached, GLHYMPS is local at 1.3 GB and needs only a different bbox, and recharge is a global Zenodo product read by byte-range window, so the US needs different row and column indices and nothing else. What is new is a US DEM window and the wells. The USGS route in `external-data.md` is live and needs no key, PROBED 2026-08-20: `api.waterdata.usgs.gov/ogcapi/v0/collections/field-measurements/items?parameter_code=72019` returns 200 with `value`, `unit_of_measure` of ft, `monitoring_location_id`, `vertical_datum` and a `qualifier` carrying Static, which is the pumping filter `external-data.md` says Fan's four columns cannot support. `numberMatched` is absent so the volume is unknown and it will need pagination. Fan's own `US_obs_wtd.txt.gz` at 7 MB is the fallback and is on the recovered AquaKnow mirror, but it lacks the well-use metadata that makes the filter possible. Canada stays lowest priority for the reason already recorded, that it is provincial rather than national and repeats the US regime. What would settle GW-21: if US rank skill is again its recharge field's own, the input-limited finding generalises; if it is materially better, Australia's confined systems were the problem and the Australian verdict is too harsh [step: hydrography] |
| GW-23 | Groundwater leaves an unconfined aquifer through unmapped first and second order gullies at hundreds of metres, not by crossing a 15.19 km cell to a mapped river, and the solver knows only the mapped network. So recharge has no local outlet and piles up until a REGIONAL gradient can carry it, which is the shape GW-21 measured: without the evapotranspiration sink 91.4% of cells pin at the surface, and with it the depth is a local recharge balance. The discarded input is real and is the right size. Binning the 1.855 km ETOPO into the 15.19 km Australian mesh, at a median 95 pixels a cell, the within-cell relief the mean throws away is 14.8 m of standard deviation and 25.9 m from mean to minimum at the median, against an observed water table depth whose median is 7.6 to 9.8 m. The relief being averaged out is about three times the signal being predicted | `hydrography/notes/earth-calibration-criterion.md`, `docs/src/reference/external-data.md` | done, tested at its bounding limit 2026-08-20, and it does not hold, but the investigation found a real defect and that is the part to keep. `C -> infinity` is a fixed head at the cell's own valley floor, which `fixed_head_m` supplies for free, so the family is bounded without a solver term: 22,284 of 26,372 Australian land cells have 10 m or more of relief to vent into, at a median of 30 m, and draining them INVERTS the correlation below the cell mean, +0.0703 to -0.0575, and costs below the bore, +0.1341 to +0.0492, over-draining to a median 25.5 m against an observed 12.9. Both ends of the family are without skill, so a finite `C` would have to beat both, which nothing measured suggests. THE DEFECT, which is the useful half: the model's depth hung from the DEM CELL MEAN and the bore's from the ground at the bore, two different surfaces, and bores sit a mean 8.4 m below their cell mean because wells are in valleys. Hanging both from the DEM takes Pearson +0.0703 to +0.1341 with NO change to the model, the largest single improvement anywhere in this work, and it means every earlier score including GW-3's and GW-21's was measured against the wrong surface. It is still not skill: R2 -3.05, spread 36.18 m against 28.19 observed, because a flat cell head plus each bore's own ground gives the depth all of the relief. Do NOT mix the observation's own surveyed elevation in: it disagrees with the DEM at sd 51.97 m and that lands whole on the residual, R2 -18.0. What remains open is the middle of the bracket rather than the drain: the true within-cell behaviour sits between a model spread of 3.96 m and 36.18 m with the observations at 28.19, which is a three-number target for any sub-grid parameterisation. Superseded approach kept for the record: a head-dependent sub-grid drain per cell rather than a finer mesh, `Q = C max(0, h - z_valley)`, with `z_valley` from the within-cell DEM minimum or `mean - 2 sigma` and `C` scaled by a drainage-density proxy. That is a Robin condition and the solver has no term for it; the existing `fixed_head_m` gives its `C -> infinity` limit for free and BOUNDS the effect, so run that first. TWO THINGS TO WATCH, both from measurements already taken. A direct correlation of within-cell relief against observed depth is only -0.10 Spearman against recharge's -0.27, but that test is confounded by the bias `external-data.md` records from Fan, that wells sit in valleys, so the bores already report valley-floor depths and relief cannot predict them; the real mismatch is that a CELL-MEAN head is being scored against VALLEY-SITED bores, which is what this term would fix. And a sub-grid drain is another LOCAL sink, so it risks reproducing the evapotranspiration sink's pathology; the difference that matters is that its balance is set by terrain, `h = z_valley + Q/C`, where the sink's is set by recharge, `d = lambda ln(et_max/supply)`, which is the dependence GW-21 found missing. Falsifier: Pearson rising off 0.07 and rho(model, elevation) moving toward the observations' +0.20. NOTE THE ASYMMETRY: this is measurable on Earth, where a 1.855 km DEM sits under a 15.19 km cell, and NOT on Vesper, where `notes/audits/orogen-resolution.md` finds the generator designs nothing below about 20 km, so `z_valley` there would have to be parameterised rather than measured. That is GW-6 with a job to do [step: hydrography] |
| GW-24 | Depth-dependent transmissivity is the governor the model lacks, and it is not the same lever as transmissivity MAGNITUDE. GW-21 raised constant `D` from 100 m to 20 km and did move the dependence structure onto the observed values, rho(model, elevation) -0.047 to +0.466 and rho(model, recharge) -0.902 to -0.245, but that raises `T` everywhere at once. A depth-dependent `T` raises it only where the table is SHALLOW, which is a self-limiting local governor rather than a global flattening, and the uniform sweep cannot stand in for it | `hydrography/notes/water-table-convergence.md`, `hydrography/notes/earth-calibration-criterion.md` | open, and NOT via Fan's exponential, which is excluded twice over and measured both times: `exp(h/f)` is convex so a cell mean carries `exp(sigma^2/2f^2)`, which at 100 m of sub-grid relief against the 0.95 m `f` Fan reaches on steep bedrock is `exp(5000)`, and Picard on it LIMIT-CYCLES because `T` moves by a factor of e per e-folding length so the map is not a contraction. Two solvers died there; do not re-run it. The form to try is saturated thickness, `T = K (h - z_bottom)`, the textbook unconfined case the constant `D` currently approximates as confined. It is depth-selective, so it carries the governor; it is LINEAR in head, so cell-averaging is exact and there is no Jensen term; and it is Boussinesq, which under the Kirchhoff substitution `u = (h - z_bottom)^2 / 2` is exactly linear for uniform `K` and `z_bottom`, so the fixed point that killed the exponential does not arise. The `h <= z` box constraint survives the substitution monotonically, as `u <= (z - z_bottom)^2 / 2`, so the complementarity structure is kept [step: hydrography] |
| GW-25 | Every Earth score this project has recorded hangs the model's depth from the DEM CELL MEAN and the observation's from the ground at the bore. Those are different surfaces: measured on 48,552 Australian bores, the DEM at a bore sits a mean 8.4 m below its cell mean, because wells are in valleys, which is the bias Fan states and `external-data.md` records. Hanging both from the DEM takes Pearson +0.0703 to +0.1341 with no change to the model at all, which is the largest single improvement found in this work and makes it a measurement defect rather than a modelling one | `hydrography/notes/earth-calibration-criterion.md` | done 2026-08-20, and the answer is that the default STAYS `cell-mean`. Treating the two surfaces as ends of a subdued-replica family, `h(x) = h_cell + alpha (z(x) - z_cell)`, and sweeping alpha shows the three measures disagreeing about where to stand: Pearson peaks near 0.8, Spearman and R2 are both best at 1.0, and the rule declared before looking at correlations, match the observed 28.49 m spread, picks 0.36 and is worse than either end on two of three. The attribution is why: `delta`, where a bore sits within its own cell, predicts observed depth at Pearson +0.1290 ALONE, against the model's +0.0703, so the near-doubling is a geometric covariate folded into the model's score rather than model skill recovered. A joint least squares on model and `delta` reaches R2 +0.0208 against a statistical fit's 0.2098. `--surface dem-at-bore` stays available and labelled, and this sharpens GW-21: the model does not merely fail to add to its recharge forcing, it is out-predicted by pure geometry carrying no groundwater physics at all. Two things to settle before flipping it. The corrected scoring is not simply better: Pearson nearly doubles while R2 goes -0.365 to -3.05 because the spread goes to 36.18 m against 28.19 observed, since a flat cell head plus each bore's own ground gives the depth ALL of the within-cell relief. Neither surface is right and GW-23 brackets the truth between them, so the honest default may be to report both rather than pick. And restating the archived numbers matters: GW-3's R2 = 0.017 against a bar of 0.07, and GW-21's conclusion that the inputs limit the score, were both measured against the wrong surface. Neither verdict flips, since +0.1341 is still far under the 0.2098 a statistical fit reaches, but the numbers quoted should say which surface they used. Do NOT reach for the observation's own surveyed elevation as the third option: it disagrees with the DEM at sd 51.97 m, median -0.28 m, so it is scatter rather than offset and it lands whole on the residual, giving R2 -18.0 [step: hydrography] |

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
| HYD-20 | The endorheic basin catalogue is resolution-limited rather than physics-limited. `selectionCriteria` declares floors of 0.05 km depth, 1000 km2 area and 12 mesh cells, and at the build's spacing the 12-cell floor is 3526 km2, so nothing below 2646 km2 survives and the carve list is drawn from a set the mesh floored rather than the criterion selected. Preserved goes 3,621, 6,345, 9,419 across one, two and four times the region count, and at four times the smallest preserved depression is 1000 km2 exactly, the declared floor having finally taken over from the mesh. The large end does NOT settle: the biggest basin is 3,059,665 km2 at the build and 5,499,601 at four times, one south polar basin remeasured 80% larger rather than two merging, and it is not erosional, since catalogues with `--glacial 0` are byte-identical | `notes/audits/orogen-resolution.md` | open. Not a defect in the floor, which is a signal-to-noise criterion and correctly relaxes as the mesh becomes able to hold smaller depressions. It is a statement about what this build can see: the declared physical floor only starts binding near 8.8M regions. DECIDED 2026-08-20: the next generation raises the region count, and the reference export at 10,000,005 is registered as `precarve-craton-10m`. It is NOT the build to activate, because a generation meant to be commissioned should also carry the PHYS-13 ice-mask decision; the two changes belong in one generation, not two. Remaining work is that generation: pick the count, emit the full grid set, register the hash, and point `source_build` at it. THE LADDER IS NOW T21/T42/T85/T127/T170, T63 is dropped, and `precarve-craton-10m` carries all five with every grid verified to the same terrain hash; `rebuild_binaries.py`'s MATRIX gains the two missing rungs so no truncation has an export without an executable. What is left of this row is the DECISION to activate, which rule 7 makes consequential and which should carry PHYS-13's ice-mask change in the same generation. CEILING ON THE COUNT, measured: Orogen's terrain noise sits at fixed physical wavelengths and the finest it DESIGNS is about 20 km, so `elevation_pre_erosion` gains only 8% of semivariance per lag doubling below 10 km against 26% at 50 km. The 7.60 km mesh at 10,000,005 regions already oversamples that floor by 2.6x. MEASURED AT 25,000,000 and the catalogue SATURATES there: preserved goes 3,621, 6,345, 9,419, 9,649 across 2.5M, 5M, 10M and 25M, so four times again the region count adds 2.4% against the 48% the previous step added, while detected depressions keep climbing 286,534 to 646,725 because every new one is below the declared 1000 km2 floor. That is the crossover working as predicted: once `minCells` stops binding the physical criterion governs and refining further buys nothing. So 10,000,005 is not adequate-with-margin, it is where this converges, and it is the recommendation; 25M costs 3,875 s and 21.1 GB against 10M's 1,093 s and 9.0 GB, superlinear in time, for 230 more basins; 50,000,000 is 3.40 km rather than the sub-km it looks like, since edge goes as the square root of the count, and a 1 km edge needs about 577 million regions, roughly 17 hours and 517 GB. Below the floor a finer mesh manufactures erosion texture at the mesh scale rather than resolving terrain. Sub-km relief is a synthesis problem, not a region-count one. Note what activating costs under rule 7: `precarve-craton` has been consumed by a commissioned run, so superseding it makes the climatology and everything derived from it worthless rather than stale [step: orogen] |

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
| MIN-6 | Replace the rainfall proxy for water table depth in the supergene rules with a modelled depth field. `minerals/config/downstream_prospectivity.yaml` says in its own comments that it substitutes rainfall because nothing better exists, while the mechanism it cites is depth: Reich and Vasconcelos (2015) put oxidation in the vadose zone and secondary sulfides below the water table | `hydrography/notes/groundwater-scoping.md` section 6, `docs/src/reference/economic-minerals.md` | open and UNBLOCKED: GW-1 closed 2026-08-20, so the mean depth field exists as the `groundwater` step's output. What is left is a judgement rather than a wait, and GW-3 HAS now supplied the external test, so the judgement can be made on evidence. The verdict is regime-dependent: on Australian bores the model has no per-cell skill, Pearson +0.0703 with a spread of 1.91 m against an observed 18.36, while on United States bores labelled unconfined it reaches +0.2599 with a matched spread. `config/pipeline.yaml`'s gate says plainly that MIN-6 cannot use this field, and that stands for a per-cell depth. What GW-22 adds is the per-cell test rather than a blanket one: `water_table.nc` ships `sink_fraction`, and a prospectivity rule could in principle key on cells where the sink did NOT set the depth. Whether enough of this planet qualifies is not yet knowable, because the current forcing is a bootstrap climatology on pre-carve terrain. It supplies the DEPTH and not the descent RATE, so Sillitoe (2005)'s wet-end control -- erosion in balance with the rate of water table descent -- stays unreachable and stays declared [step: downstream_prospectivity]

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
| PHYS-13 | Orogen places the ice that carves this terrain by latitude and dimensionless elevation alone. At the build's `glacialErosion` of 0.8 the ice line falls at 58 degrees latitude, positioned by the erosion slider and nothing else, and the altitude gate is in the model's elevation parameter so the snowline tracks maximum relief, which scales as 1/g, rather than a freezing altitude set by lapse rate and surface temperature. `glacIdx` consults no temperature at all: not the K2.5V spectrum, not the 32 degree obliquity, not the 30 hour rotation, and not even Orogen's own temperature field, which `computeTemperature` produces after `runPostProcessing` has already eroded. All three leading parameters point the same way, toward less polar ice than an Earth-tuned 58 degree threshold gives | `notes/audits/orogen-resolution.md` | open. The correction needs no new physics, but NOT the field first proposed here. `glac` is identically zero across all twelve bins and all 8,192 cells while the run manifest records `"glaciers": true`, and `snd` peaks at 0.29 m, which is seasonal snow. That zero is a fact about the grid, not the planet: the sub-grid peak excess `orog_max - orog_mean` on land averages 1.355 km area-weighted and 2.969 km at the 90th percentile, so at the measured warm-season lapse rate of 7.372 K/km the model evaluates its own high ground about 7.8 K too warm on average and 21.9 K too warm in the top tenth, and no ice survives that. USE INSTEAD the warmest-month surface temperature, which is smooth and genuinely resolved at T42, corrected to each mesh region's own elevation through `lib/lapse.py`'s freezing-height criterion. `maps/build_basemap.py` lines 379-385 already implement exactly this for the basemap's ice shading. Measured that way, land below freezing in the warmest month goes from 0.002% on the model's orography to 1.657% on the 15.19 km mesh, and `analysis/ice_mask_freezing_height.py` is the registered computation. Raising the model truncation is not the alternative: T85 is about 150 km against a 7.59 km mesh, so the peak excess stays sub-grid at any resolution the model can afford. Drive `glacIdx` from that mask instead of `polarDist` and close it through loop A, which exists to feed a commissioned climate into the next generation; the terrain-before-climate ordering that blocks it inside one pass is not a constraint across iterations. SHAPE OF THE FIX, decided 2026-08-20: an ice mask consumed AT GENERATION like the carve list, `--ice-mask FILE`, not a post-hoc terrain alteration. Glacial, hydraulic and thermal erosion share one iteration loop and a mid-loop priority flood at `GLACIAL_MID_FLOOD_FRAC` exists to cut outlets through the depressions glaciation makes, so a separate later pass would leave them undrained and move the drainage network and the basin catalogue. Bootstrap with `--glacial 0`, an honest null, then commission and regenerate with the mask. Two constraints on it: mapping `glac` from the ExoPlaSim grid back to the mesh is a rule 3 crossing and must go by index or a shared coordinate source, and per `docs/src/reference/no-time-axis.md` a current ice field can only say WHERE ice was, never how long, so the strength slider stays a declared choice rather than a measurement. Subsumes the mechanism half of GRAV-6 [step: orogen] |

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
| SURF-7 | Place groundwater silcrete and calcrete's groundwater-calcite pathway once a water table exists. `pedology/config/surface_classes.yaml` currently refuses groundwater silcrete by name, with the reason in the config: it sits at or near a water table and this project models none | `hydrography/notes/groundwater-scoping.md` section 6, `pedology/notes/derived-surface-classes.md` | open and UNBLOCKED: GW-1 closed 2026-08-20 and the `groundwater` step generates the field, so this is now work rather than a wait. GW-3 HAS since supplied the external test and the answer is regime-dependent rather than a blanket uncertification: no per-cell skill in the arid Australian case, Pearson +0.2599 on United States unconfined bores, and `sink_fraction` in `water_table.nc` separates the two per cell. The gate's measured AUC of 0.52 to 0.57 for this mask was taken over ALL cells, so re-measuring it restricted to cells where the sink did not set the depth is the cheap thing to try before treating the mask as dead. What becomes placeable is the MEAN depth field and a discharge mask, not Fenske et al. (2025)'s model: that hardens a layer over the RANGE of water table fluctuation at 10^5 years or longer, and both are durations `docs/src/reference/no-time-axis.md` refuses. Gypcrete is untouched, being surface and air processes rather than this mechanism [step: surface_classes]

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
