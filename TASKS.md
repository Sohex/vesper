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
| MIN-2 | Place the weathering, drainage and brine deposit types downstream: bauxite, laterite Ni, supergene Cu, placers, and the brine evaporite minerals | `notes/economic-minerals.md` | blocked on the re-baseline -- all of these need a climate. The brine half is close to free: `brine_paths.py` already solves the chemical divide per basin, which is what decides the mineral |
| VOLC-3 | Place silcrete and diatomite against per-basin silica supply. A large share of silica delivery goes to closed basins, which is the setting for both, but per-basin placement needs hydrography run on the arc-bearing build | `pedology/notes/derived-surface-classes.md` | blocked on the re-baseline |
| BIO-1 | Regenerate the LPJ-GUESS driver: `vesper_driver_provenance.json` is pinned to a superseded config and `check_consistency.py` has been failing on it. Its provenance also records no generator, so the check cannot name the remedy | `scripts/check_consistency.py` | blocked on the re-baseline -- the driver is built from a climatology, so rebuilding it now would only have to be redone |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here) | `notes/audits/orogen-gravity.md` | blocked on the downscaling pass, deliberately. The audit's recommendation was Option B, carry it as a declared gap: `iceFlow` is a heuristic accumulation rather than a Glen's-law velocity, so a rigorous g^3 on it is precision theatre. Decisive on top of that, a glacial valley is 1-5 km against a 15.19 km mesh cell, so the process is SUB-GRID and cannot be resolved here at any gravity. It belongs with the downscaling machinery, which has to persist sub-grid hypsometry anyway |
| SURF-1 | Implement the derived-surface classifier against the two-axis design | `pedology/notes/derived-surface-classes.md` | blocked on the re-baseline -- there is no climatology at all now, and `baseline_climatology` is null until one exists |
| SURF-2 | Re-run `brine_paths.py` weighted by discharge rather than catchment area | `pedology/notes/derived-surface-classes.md` | blocked on the re-baseline, same reason as SURF-1 |
| SURF-3 | Loess: needs the deposition field DUST-1 produces | `pedology/notes/derived-surface-classes.md` | blocked on DUST-1 |
| DUST-1 | Finish the offline dust component: the machinery, source map and grounded physics are built at `aeolian/`, but emission comes out 3-4 orders above Earth's 2000 Tg/yr and transport does not reach steady state, so no number from it is usable | `aeolian/README.md`, `notes/dust.md` | open -- two known issues, the unconverged advection and an undiagnosed factor of ~50 in the emission chain. The wind-shape error is already found and fixed |
| DUST-4 | Decide whether a 12-bin climatology can drive dust at all. Mean friction velocity over the source cells sits BELOW threshold, so emission is entirely a tail property and the driver has averaged the tail away. This is a different failure mode from the one `notes/dust.md` reopens on, which is feedback magnitude | `aeolian/README.md` | open -- a pipeline decision rather than a component one. Options are higher-cadence output from the baseline run, a fitted subgrid distribution, or conceding the in-model route for a reason the note did not anticipate |
| DUST-2 | Price the dust radiative forcing from DUST-1's optical depth and the optics already computed, longwave included, and put it in the error budget as a number rather than an unpriced item | `notes/dust.md`, `analysis/dust_optics.json` | blocked on DUST-1 |
| DUST-3 | Feed the dust field forward as a PRESCRIBED aerosol distribution in the next iteration's run, and reopen the interactive question only if DUST-1 exceeds a land-mean optical depth of 0.10 or 1.5 W/m2 of global forcing | `notes/dust.md` | blocked on DUST-1 and DUST-2 |
| HYD-3 | Measure the overshoot: re-evaluate an already-carved set against the climate that carving produced, and report how many would no longer have carved | `WORKFLOW.md` section 4 | blocked -- it applies to a verdict taken on carved terrain, and a pre-carve base restarts the count at iteration 1 |
| HYD-6 | Check the first real verdict against the retain change: how many basins sit between 0 and 1 for want of discharge rather than for uncertainty, and whether the lake-fed spillers are a handful or a population. The second decides whether `Q_full` is worth calibrating | `hydrography/notes/retain-fraction.md` | blocked on the re-baseline |
| CLIM-1 | Name the constant -0.455 W/m2 gap between the top-of-atmosphere and surface budgets from the 28-term decomposition on codes 360-387. Needs a settled run: a segment taken off a restart sits several W/m2 out of balance | `exoplasim/notes/water-and-energy-closure.md` | blocked on the re-baseline |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; the span is 18.2% of NPP and it is the largest nitrogen lever | `biosphere/notes/productivity-prediction.md` | blocked on the re-baseline -- no LPJ-GUESS run exists on this build |


## Closed

Moved to `archive/tasks.md`, with the reasoning that closed each one. Ids are
never reused, so look there when a commit cites an id this file does not have.
