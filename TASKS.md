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
| BIO-1 | Regenerate the LPJ-GUESS driver: `vesper_driver_provenance.json` is pinned to a superseded config and `check_consistency.py` has been failing on it. Its provenance also records no generator, so the check cannot name the remedy | `scripts/check_consistency.py` | unblocked: the baseline climatology exists, so the driver can be rebuilt against it |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here) | `notes/audits/orogen-gravity.md` | blocked on the downscaling pass, deliberately. The audit's recommendation was Option B, carry it as a declared gap: `iceFlow` is a heuristic accumulation rather than a Glen's-law velocity, so a rigorous g^3 on it is precision theatre. Decisive on top of that, a glacial valley is 1-5 km against a 15.19 km mesh cell, so the process is SUB-GRID and cannot be resolved here at any gravity. It belongs with the downscaling machinery, which has to persist sub-grid hypsometry anyway |
| SURF-1 | Implement the derived-surface classifier against the two-axis design | `pedology/notes/derived-surface-classes.md` | unblocked: `baseline_climatology` now names the baseline |
| SURF-2 | Re-run `brine_paths.py` weighted by discharge rather than catchment area | `pedology/notes/derived-surface-classes.md` | unblocked, same reason as SURF-1 |
| SURF-3 | Loess: needs the deposition field DUST-1 produces | `pedology/notes/derived-surface-classes.md` | blocked on DUST-1 |
| DUST-1 | Offline dust component at `aeolian/`. Chain is sound and converged; the answer spans below to well above the reopening threshold across two unmeasured parameters, the wind tail and the aeolian roughness | `aeolian/README.md`, `notes/dust.md` | open -- needs the high-cadence segment in DUST-5 to decide it. Everything else is done |
| DUST-5 | Run one orbit of high-cadence output off the settled baseline, sampled every four timesteps, and fit the gust distribution from it instead of from 32 snapshots. That is the parameter the reopening test turns on: at the measured k = 3.98 dust is below threshold, at an Earth-like k = 2.0 it is well above | `aeolian/README.md` | unblocked: the baseline has settled |
| CLIM-2 | The first output bin of every orbit has corrupted WIND and HUMIDITY only -- ua 13.0x, spd 7.57x, va 0.72x, hus 0.71x, worsening downward from 1.02 at the model top. Every scalar is clean to better than 2%, so the flux bracket, soil, code 229, convergence assessments and the baseline run all stand | `notes/failure-modes.md` class 14 | open -- find the cause in the first output interval's accumulation. The Penman consumers are NOT affected in the way first thought; see CLIM-3 |
| CLIM-3 | Point `carve_verdict.py`, `surface_water.py` and `export_carve_list.py` at the SNAPSHOT climatology's `spd` for the Penman wind, which is the mean of instantaneous speed at 7.501 m/s, rather than the binned product's speed of the time-mean vector at 5.14. The binned-vs-snapshot gap is `build_climatology.py` averaging ua and va as vectors while spd averages as a scalar. Current usage is 7.732, accidentally +3.1%, because the corrupt bin's inflation and vector averaging nearly cancel. Do NOT fix by dropping bin 0 alone: that gives 4.997 and turns a 3% error into 33% | `notes/failure-modes.md` class 15 | open -- a 3.1% change, so no rebuild of the lake solution is required on its account |
| DUST-2 | Price the dust radiative forcing from DUST-1's optical depth and the optics already computed, longwave included, and put it in the error budget as a number rather than an unpriced item | `notes/dust.md`, `analysis/dust_optics.json` | blocked on DUST-1 |
| DUST-3 | Feed the dust field forward as a PRESCRIBED aerosol distribution in the next iteration's run, and reopen the interactive question only if DUST-1 exceeds a land-mean optical depth of 0.10 or 1.5 W/m2 of global forcing | `notes/dust.md` | blocked on DUST-1 and DUST-2 |
| HYD-3 | Measure the overshoot: re-evaluate an already-carved set against the climate that carving produced, and report how many would no longer have carved | `WORKFLOW.md` section 4 | blocked -- it applies to a verdict taken on carved terrain, and a pre-carve base restarts the count at iteration 1 |
| HYD-6 | Check the first real verdict against the retain change: how many basins sit between 0 and 1 for want of discharge rather than for uncertainty, and whether the lake-fed spillers are a handful or a population. The second decides whether `Q_full` is worth calibrating | `hydrography/notes/retain-fraction.md` | blocked on the carve verdict, which can now be run |
| CLIM-5 | Patch PlaSim's low-I/O output path, which writes a corrupt first record per model call. Runs now set `NLOWIO = 0` to avoid it at 25x the output volume, so this buys back disk rather than correctness | `exoplasim/notes/first-output-bin.md` | open -- do NOT write it on the obvious candidate. `naccuout` persists across runs through the restart while the accumulators do not (`plasim.f90:767`), but a counter error scales a field uniformly and this defect changes its vertical structure, 1.04 at the model top against 11.2 at the surface. The mechanism is unpinned |
| CLIM-1 | Name the constant -0.455 W/m2 gap between the top-of-atmosphere and surface budgets from the 28-term decomposition on codes 360-387. Needs a settled run: a segment taken off a restart sits several W/m2 out of balance | `exoplasim/notes/water-and-energy-closure.md` | unblocked: the baseline is settled and carries codes 360-387 and 460-487 |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; the span is 18.2% of NPP and it is the largest nitrogen lever | `biosphere/notes/productivity-prediction.md` | blocked on the re-baseline -- no LPJ-GUESS run exists on this build |


## Closed

Moved to `archive/tasks.md`, with the reasoning that closed each one. Ids are
never reused, so look there when a commit cites an id this file does not have.
