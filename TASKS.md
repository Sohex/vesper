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
  classes, `CLIM` climate, `HYD` hydrography, `BIO` biosphere, `MIN` economic minerals, `REF`
  references and provenance.
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
| SURF-3 | Loess: needs a dust emission scheme and a transport path, neither of which exists | `pedology/notes/derived-surface-classes.md` | blocked |
| HYD-1 | Route the zero-catchment-runoff basins whose lake surface has P > E through the marginal band, rather than the blanket preserve the guard at `export_carve_list.py:180` gives them | `hydrography/notes/carve-verdict-open-items.md` | blocked on the re-baseline -- the set has to be re-measured against a climate before it can be handled |
| HYD-2 | Make the retain fraction a function of overflow discharge rather than of distance from the evaporation threshold | `hydrography/notes/carve-verdict-open-items.md` | open -- the discharge is already computed, so the mapping can be written now, but it cannot be checked until a verdict exists |
| HYD-3 | Measure the iteration-1 overshoot: re-evaluate the already-carved set against the new climate and report how many would no longer have carved | `WORKFLOW.md` section 4 | blocked on the re-baseline |
| HYD-4 | Validate the lake equilibrium solver against something. It is the least-checked product in the pipeline, and it is the only thing that could break the crust-extent against crust-albedo degeneracy | `exoplasim/notes/parameter-decisions.md` | blocked on the re-baseline |
| CLIM-1 | Name the constant -0.455 W/m2 gap between the top-of-atmosphere and surface budgets from the 28-term decomposition on codes 360-387. Needs a settled run: a segment taken off a restart sits several W/m2 out of balance | `exoplasim/notes/water-and-energy-closure.md` | blocked on the re-baseline |
| BIO-2 | Quote productivity with the `nfix_a` bracket 0.102-0.367 carried through rather than the central value alone; the span is 18.2% of NPP and it is the largest nitrogen lever | `biosphere/notes/productivity-prediction.md` | blocked on the re-baseline -- no LPJ-GUESS run exists on this build |


## Closed

Moved to `archive/tasks.md`, with the reasoning that closed each one. Ids are
never reused, so look there when a commit cites an id this file does not have.
