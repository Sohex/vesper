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
  `DUST` the aeolian component, `REF` references and provenance.
- One line per task. If it needs a paragraph, it needs a findings document.
- A task that turns out to be wrong is closed `wontfix` with the reason, not
  deleted. The reasoning is the point.
- A closed task moves to `archive/tasks.md` rather than staying here, so this
  file shows only what is left. Ids are never reused, so a cited id that is not
  here is there.

---

## Open

| id | task | source | status |
| --- | --- | --- | --- |
| MIN-2 | Place the weathering, drainage and brine deposit types downstream: bauxite, laterite Ni, supergene Cu, placers, and the brine evaporite minerals | `notes/economic-minerals.md` | unblocked: the baseline climatology exists. The brine half is close to free: `brine_paths.py` already solves the chemical divide per basin, which is what decides the mineral |
| VOLC-3 | Place silcrete and diatomite against per-basin silica supply. A large share of silica delivery goes to closed basins, which is the setting for both, but per-basin placement needs hydrography run on the arc-bearing build | `pedology/notes/derived-surface-classes.md` | unblocked: hydrography and a climatology both exist on this build |
| BIO-1 | Regenerate the LPJ-GUESS driver: `vesper_driver_provenance.json` is pinned to a superseded config and `check_consistency.py` has been failing on it | `scripts/check_consistency.py`, `WORKFLOW.md` A3 | deliberately DEFERRED until after the iteration-2 re-baseline. The driver is built from a climatology and a soil that the carve replaces, so doing it now is make-work whose only effect is silencing a warning that is currently telling the truth: no biosphere inputs exist for this build. The orbit-dependent half was durable and is done (BIO-3) |
| CYC-1 | Run the stellar cycle off the settled baseline. `WORKFLOW.md` section A2 requires it AFTER a converged baseline and BEFORE the carve verdict, and nothing was tracking it: carving is irreversible, so the terrain ratchets toward the cycle's wet extreme and a verdict taken on the mean climate systematically under-carves. The run also measures the damping factor, known only as 0.24 to 0.6, which converts any future flux amplitude into a climate without another run | `WORKFLOW.md` section A2, `config/planet.yaml` stellar_cycle | open, and ON THE CRITICAL PATH. The cycle binary is built and verified. Length is set in periods of the LONG component, 57 Earth years, not the medium one. The carve list of 2026-08-17 was taken on the mean climate and should be re-taken after this |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here) | `notes/audits/orogen-gravity.md` | blocked on the downscaling pass, deliberately. The audit's recommendation was Option B, carry it as a declared gap: `iceFlow` is a heuristic accumulation rather than a Glen's-law velocity, so a rigorous g^3 on it is precision theatre. Decisive on top of that, a glacial valley is 1-5 km against a 15.19 km mesh cell, so the process is SUB-GRID and cannot be resolved here at any gravity. It belongs with the downscaling machinery, which has to persist sub-grid hypsometry anyway |
| SURF-1 | Implement the derived-surface classifier against the two-axis design | `pedology/notes/derived-surface-classes.md` | unblocked: `baseline_climatology` now names the baseline |
| SURF-2 | Re-run `brine_paths.py` weighted by discharge rather than catchment area | `pedology/notes/derived-surface-classes.md` | unblocked, same reason as SURF-1 |
| SURF-3 | Loess: needs the deposition field DUST-1 produces | `pedology/notes/derived-surface-classes.md` | blocked on DUST-1 |
| DUST-1 | Offline dust component at `aeolian/`. SETTLED and the answer is that this is a dusty world: land-mean optical depth 0.740 at the central roughness against a reopening threshold of 0.10, 7.098 at the smooth end, and only the roughest end below it at 0.0005 | `aeolian/README.md`, `notes/dust.md` | done as an offline component. What remains is not DUST-1's: the pre-committed threshold has fired by 7.4x, so DUST-2 and DUST-3 inherit it |
| CLIM-2 | The first output bin of every orbit has corrupted WIND and HUMIDITY only -- ua 13.0x, spd 7.57x, va 0.72x, hus 0.71x, worsening downward from 1.02 at the model top. Every scalar is clean to better than 2%, so the flux bracket, soil, code 229, convergence assessments and the baseline run all stand | `notes/failure-modes.md` class 14 | open -- find the cause in the first output interval's accumulation. Two parts of this are now closed: the consumers are corrected (CLIM-3), and code 259 is settled, being exactly `sqrt(ua^2+va^2)` per sample in the snapshot product |
| DUST-2 | Price the dust radiative forcing from DUST-1's optical depth and the optics already computed, longwave included, and put it in the error budget as a number rather than an unpriced item | `notes/dust.md`, `analysis/dust_optics.json` | unblocked: DUST-1 is settled. The optical depth to price is 0.740 at the central roughness, so the forcing is very unlikely to come in under the 1.5 W/m2 half of the reopening test either |
| DUST-3 | Put the emission scheme in the ExoPlaSim fork, keeping the source map outside it as a boundary field. DECIDED, see `notes/dust.md`. Transport and radiative coupling already exist -- `aerocore.f90` is Lin-Rood with gravitational settling, `radmod.f90` builds per-layer AOD under `l_aerorad`, and `configure` exposes aerosol, aerorad and aerofile -- and `aerocore.f90:621` already has a `case(2) ! dust` hook whose whole body is `mmr(:,:,NLEV,ic) = fcoeff*land`. That line is the change | `notes/dust.md` | open. Four parts: erodibility/clay/z0 as .sra boundary fields; Kok emission in the source case against the model's own u*, soil moisture and snow; an `aerofile` from `analysis/dust_optics.json`; and `apart`/`rhop` off haze values, which needs NAERO checked for whether size bins survive. Patch and rebuild per CLAUDE.md rule 4. SEQUENCING: this moves land-surface forcing, so `WORKFLOW.md` A3 applies -- do not land it in the same iteration as the carve |
| HYD-3 | Measure the overshoot: re-evaluate an already-carved set against the climate that carving produced, and report how many would no longer have carved | `WORKFLOW.md` section 4 | blocked -- it applies to a verdict taken on carved terrain, and a pre-carve base restarts the count at iteration 1 |
| CLIM-5 | Patch PlaSim's low-I/O output path, which writes a corrupt first record per model call. Runs now set `NLOWIO = 0` to avoid it at 25x the output volume, so this buys back disk rather than correctness | `exoplasim/notes/first-output-bin.md` | open -- do NOT write it on the obvious candidate. `naccuout` persists across runs through the restart while the accumulators do not (`plasim.f90:767`), but a counter error scales a field uniformly and this defect changes its vertical structure, 1.04 at the model top against 11.2 at the surface. The mechanism is unpinned |
| CLIM-6 | Give `denergy` an accumulator in `outmod.f90` beside `ashfl`, `alhfl` and the rest, so the 28 energy terms are time means like every other field on the stream rather than instantaneous values at the output step | `exoplasim/notes/water-and-energy-closure.md` | open, low priority -- NLOWIO = 0 already makes the terms usable, so this buys back the 25x output volume rather than correctness. Same standing caution as CLIM-5: it is a patch to a compiled model and every binary needs rebuilding |
| CLIM-1 | Name the residual between the top-of-atmosphere and surface budgets. `ntr - (rss+rls+hfss+hfls)` is -0.776 W/m2 on the settled baseline, of which the snowmelt fusion booking in `hfns` accounts for +0.422, leaving a gap of -0.523 over the climatology and -0.354 on the one clean orbit | `exoplasim/notes/water-and-energy-closure.md` | open, and much narrower. Three leads are CLOSED: the 11+15 non-cancellation is the missing latent heat of fusion in `mklsp`, predicted and measured to 0.7%; the radiation diagnostics are consistent to 0.03 W/m2; and the gap is NOT the fixed -0.455 this was opened on. The reason the 28-term instrument had not answered is that `denergy` is written unaccumulated while every flux it is compared against is a time mean, so under NLOWIO = 1 the terms are snapshots -- 7.14% error on a pair that is algebraically an identity, 0.15% under NLOWIO = 0. Next step is several more NLOWIO = 0 orbits, to separate the annual residual from the +/-7 W/m2 seasonal swing it hides in |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; the span is 18.2% of NPP and it is the largest nitrogen lever | `biosphere/notes/productivity-prediction.md` | blocked on the re-baseline -- no LPJ-GUESS run exists on this build |


## Closed

Moved to `archive/tasks.md`, with the reasoning that closed each one. Ids are
never reused, so look there when a commit cites an id this file does not have.
