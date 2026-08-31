# Handoff

Written 2026-08-30 at the close of a long session (76 commits, ~45 rows
closed). This file is a MAP for the next session: state, the one next action,
and the context that would otherwise be re-derived. Everything here cites a
row, a note or a file; nothing here is the record itself. Delete or rewrite
this file when its content is consumed -- it describes one moment.

## State at handoff

- Tree clean at `87f484376`. Gates: `smoke_test.py` 70/0,
  `verify_entry_points.py` 160/160, `verify_model_compiles.py` passes.
- `check_consistency.py`: 61 checks, 6 failures, ALL pre-existing
  carve2-transition staleness -- old runs on `canonical-10m-base` (the
  flux-to-kelvin slope pair, cold-start temperature, staged-namelist and
  staged-surface drift on superseded runs) plus one ocean gate report that
  re-stales whenever `continue_exoplasim.py` changes. None names current work.
- No model process running; the host lock is free.
- `bd ready` is the work list. This file only annotates rows whose context is
  expensive to rebuild.

## The single highest-value next action

**Run LPJ-GUESS at the derived spin-up.** ~55 minutes at npatch 5, command on
`world-qcse`. `biosphere/generated/vesper_pfts.ins` already carries
`nyear_spinup 2318` (derived, not Earth's convention; the derivation is in
`lib/run_lengths.py` and `biosphere/notes/equilibrium-trend-null.md`).

What it unblocks, in order:
1. The first LPJ run that can PASS acceptance (every prior run failed; the
   why is `world-qcse`'s comments -- spin-up short by 2.32x, and the old trend
   statistic was invalid inside one memory time; contract /3 fixed both).
2. `require_lpj_acceptance` consumers, notably `build_surface_albedo.py` --
   which unblocks `world-cp27`'s last two initial-state fields
   (`land_background_albedo`, `land_vegetation_cover` in
   `config/partial_surface.yaml`), which unblocks the SPAT-5 seam chain.
3. Re-reading the ecological timescales on the longer record
   (`assess_lpj_run.py` records them whatever the verdict); watch whether the
   1253-cycle record floor moves under itself -- the memory time grows with
   the window it is read on.
4. `world-rnns` (the replacement trend statistic): its author declined to ship
   it twice, and the condition that lifts the decline is exactly this run's
   record length. Design and three known subtleties are on the row.
5. `world-gfut` (the labile-carbon bracket sweep) needs an accepted arm pair.

## The SPAT-5 chain (the session's main line of work)

Order is enforced by `analysis/partial_surface_decision_gate.py`, which names
the next actionable seam. Chain: `world-cp27` (initial state; 6 of 8 fields
available, 2 blocked on LPJ acceptance above) -> `world-uiq7` (tile state:
both tiles advance on their own fraction; landmod is already fraction-shaped,
seamod's `dls < 0.5` must become `dlf < 1.0`) -> `world-yqd7` (dual exchange)
and `world-9do0` (tile diagnostics) -> `spat-5` closes -> `ocn-11` unblocks.
The evidence note is `notes/audits/partial-cell-tile-state.md`.

## Open P1s with context

- `world-zq9v` -- gate `biosphere/config/fire.yaml`. The register exists (born
  in the per-entry execution shape) but NOTHING READS IT, so the fire
  divergence can silently revert. The row carries the four assertions the gate
  must make, verified by hand, including the `deleted_line` shape no existing
  gate handles.
- `world-rnns` -- see next-action item 4 above. Do not ship it early; the
  author's declines were correct both times and the row says why.

## Rows whose context is on the row (read comments before starting)

- `world-hkcj` (BIOGEM threading): a measurement series was planned, queued
  and deliberately cut at session close; the full plan, the measured 36x36 and
  72x72 curves, the futex mechanism, and the unbuilt-with-trigger NOWAIT
  change are all handed over on the row and in
  `notes/audits/cgenie-parallelism-and-coupling-support.md` section 7g.
- `world-90y6` (converted arm's ice remembers its donor, 16 sigma): needs a
  third arm from a different donor. The controlled-pair result is in
  `notes/audits/resolution-ladder.md`.
- `world-3y2v`: T85 is now the only rung with no endurance row.
- `world-mqzk` is blocked on `world-ckbt` (surface-to-TOA attenuation), which
  is unstarted and is the one unmeasured link in the drift-tolerance chain.
- `world-6247` (`SCRIPT_DIRS` omits four components): still open; every
  file-based smoke check remains blind to `aeolian/`, `ocean/`, `minerals/`,
  `analysis/`. Two checks work around it with their own file lists.

## Standing intent recorded this session (not yet formalised)

The author stated 2026-08-30, deciding `world-dnrr`: the carve2 lineage "quite
literally is canonical" and re-running the commissioning is not wanted. The
FORMAL canonical-climatology-lineage declaration (CLAUDE.md vocabulary) has
still not been made, so rule 7 still classes every run as disposable. A future
session proposing anything that would invalidate the carve2 commissioning
should raise that contradiction explicitly first; declaring the lineage may
now be the author's preference and would change what rule 7 protects.

Also settled: design flux CHOSEN at 0.945 on the choosable ground
`cold_extreme_cap_admits_no_candidate` (record in
`exoplasim/analysis/design_flux.json`, decision in
`notes/audits/design-flux-two-point-response.md`). The cap stands as a
recorded unmet preference; do not re-litigate without new terrain or a
re-declared cap.

## Stale derived artifacts (resting state, NOT a work list)

Regenerate only when a step needs them (working agreement):
- Pedology: `soil_report.json`, `soilmap_T21.txt`, `land_column_states_T21.txt`
  carry FIVE stacked pH changes (two-buffer form, sourced supplies, playa
  remap, sump buffer) and need one `build_soil.py --state baseline
  --iteration 0` run.
- Aeolian: `dust_baseline.{json,nc}`, `dust_intensity_levers.json` need one
  `build_dust`-class run after the weighting fixes (commit e39de02e9 records
  which fields moved). The sea-salt and sulfate deposition carriers have never
  been generated; they materialise when their steps next run.

## Traps found this session worth knowing before they bite again

- Every T21 restart on disk predates the partial-cell tile records, so
  `convert_restart` REFUSES all of them; a ladder conversion needs a fresh
  donor cut first (cost: it turned a ~3 h series into ~4.5 h). Fresh donors
  now exist for T21/T42 at dt 45 -- see `exoplasim/runs/INDEX.json` and
  `exoplasim/analysis/ladder/*.conversion.json`.
- `bd update --notes` REPLACES; use `bd comment`. `bd create --id X --force`
  silently overwrites. (Both already in bd memories.)
- Multi-agent staging: stage by explicit path AND verify the file belongs to
  one change -- two HEAD-breaking collisions this session came from shared
  files (`smoke_test.py`, `verify_high_cadence_rescue.py`) carrying two
  agents' work.
- The delegation agreement added this session (CLAUDE.md working agreements,
  argued in `docs/src/practice/working-agreements.md`): an agent works the
  rows it files; the delegator draws file boundaries wide enough that a thread
  does not immediately cross them.
