# Handoff

Written 2026-08-31 at the close of a long session (232 commits). This file is a
MAP for the next session: state, the decisions waiting, and the context that
would otherwise be re-derived. Everything here cites a row, a note or a file;
nothing here is the record itself. Delete or rewrite this file when its content
is consumed -- it describes one moment.

## State at handoff

- Tree clean at `681f3265c`. `smoke_test.py` 74/0, `verify_entry_points.py`
  219/219, all ten ExoPlaSim binaries and the LPJ-GUESS binary verify against
  the source they were built from.
- `check_consistency.py`: 62 checks, 10 failures. They are NOT one class and
  the count is up from 6 because three of them are new information, not new
  damage. Read them as four groups:
  - **carve2-transition staleness, pre-existing** -- namelist drift and staged
    surface drift on superseded runs, and one ocean gate report that re-stales
    whenever `continue_exoplasim.py` changes, exactly as the previous handoff
    predicted.
  - **recordless runs, newly legible** -- three failures name
    `run_432e5e46adef` as "in no run index at all". That wording is new because
    `exoplasim/runs/INDEX.json` became a ledger this session (`world-ww6z`); the
    runs were already gone. `world-qc1k` carries the enumeration.
  - **toolchain drift** -- the host's cmake moved 4.4.2 to 4.4.3 mid-session, so
    the ExoPlaSim binaries report their build toolchain as moved. The LPJ binary
    was rebuilt and is current; the ExoPlaSim ten were not. One
    `rebuild_binaries.py` clears it and nothing is blocked meanwhile.
  - **generated inputs stale under this session's own config changes** -- the
    three-point wet-soil albedo key moved under `vesper_provenance.json` and
    the driver provenance. Resting state under the working agreement:
    regenerate when a step needs them.
- `bd ready` is the work list. This file annotates only rows whose context is
  expensive to rebuild.

## The two decisions waiting, and neither is mine to take

**1. Buy the 2.5-hour LPJ acceptance run, with `--save-state`?**

The equilibrium problem is SOLVED. `lpj_de8a0c8b` -- the derived 8600-cycle
spin-up with a 1253-cycle record -- is the first LPJ run this project has made
that passes every equilibrium check. `assess_lpj_run.py` evaluates stability
before closure, so reaching a closure refusal means every assessed quantity
settled.

What stands between here and a first ACCEPTED run is conservation, and both
defects are diagnosed with their fixes committed. The run that would attempt
acceptance has not been bought.

If it is bought, pass `--save-state`. It is deliberately ungated -- writing a
state changes no number in the run that writes it -- costs about 300 MB, and if
`world-glu7` is later repaired that run becomes continuable rather than
re-bought. Buying those 2.5 hours WITHOUT it is the only irreversible choice
available. `world-gqyp` carries the argument.

**2. Finish the water closure's second clause first?**

`world-24xd` pre-registered two conditions. The first is met: the failing
gridcells went 1578 to 81 of 1617 once soil evaporation, canopy interception and
`maet.out` replaced the survivors-only transpiration. The second -- that what
remains is bounded and does not grow with window length -- is UNMEASURED. It is
a read on tables already on disk in `biosphere/runs/lpj_2ffc33a5`, about twenty
minutes, no new run. Worst remaining residual 68.07 mm against a 17.80 mm limit,
and the three worst cells are tropical, which is where the burnt-fraction
correlation would live if the remainder is still turnover.

Doing this before the acceptance run is cheap insurance: if the pairing is
wrong, the acceptance run refuses on water and the whole spin-up is re-bought.

## What the biosphere's instruments now are, because they were all replaced

Do not read older notes on these as current; all four were rebuilt this session.

- **The equilibrium test** is `vesper-lpj-equilibrium-window/6`: an upper
  confidence bound on each field's end-to-end drift over the WHOLE record,
  passed when the bound is inside the tolerance. Intersection-union makes the
  run-level false-acceptance rate equal the per-field one, so there is no
  window, no memory-adequacy guard and no multiplicity correction left to be
  self-referential. `biosphere/notes/equilibrium-trend-null.md`.
- **The tolerance is scoped to what a consumer reads**, and a reader that
  imposes no drift tolerance is not a consumer for that purpose: the closure's
  conservation identities hold at any state of drift, so their four pool and
  flux totals left the assessed set and keep reported bounds as diagnostics.
  Contract 6 is that change, and it supersedes every acceptance artifact on
  disk. `world-mxmr`, `world-cqoc`.
- **The cover tolerance is derived through its consumer**, by
  `biosphere/scripts/derive_cover_tolerance.py`, which declares nothing and
  reads every link from the artifact that owns it. The bar is the climate arm's
  own offset tolerance over its own resolving factor. `world-mqzk`.
- **The spin-up is derived and needs no relaxation time**: the requirement is
  bounded over all tau, so 6.86325 times the record satisfies it at every tau.
  8600 cycles. `world-nhhm`.
- **The closure stocks were incomplete.** `npool.out` Total is not the nitrogen
  stock -- four of six soil mineral pools are missing from it and their gas is
  reported only when it leaves, so nitrogen resident in them at a cycle boundary
  is on neither side. `notes/audits/closure-stocks-are-incomplete.md`.

## Open P1s with context

- `world-glu7` -- a run resumed from an LPJ state file is not the run it
  continues, at the 1e-4 level, everywhere, from the first resumed year. No pool
  is lost; something is not carried exactly. `Soil` is fully triaged at 0 lost,
  and every OTHER Serializable class is untriaged -- `Individual` streams 100 of
  139, `Patch` 50 of 67 -- which is where it most likely lives. `--one-day` is
  the instrument that localises it. Blocks `world-gqyp`.
- `world-24xd` -- see decision 2 above.
- `world-ccx6` -- `build_soil.py` takes no `--grid` and its climatology
  declaration carries no rung. Blocks the escalation route above T21, and
  `world-512r` settled the policy question in front of it: the soil WAITS for a
  climate at its rung, because a remapped climatology is not a less determined
  climate, it is another world's. The rung is not a stage.
- `world-qcse` -- its original question is answered and its title's premise was
  corrected twice. Keep it open only until an accepted run exists.
- `fire-7`, `pcar-1` -- untouched this session.

## Rows whose context is on the row (read comments before starting)

- `world-90y6` closed GRADED, and its finding is not what the row asked. North
  of 40 degrees two converted arms are THE SAME ARM from initial fields 0.0326
  apart; south of 40 they separate ordered with the initial ice. The global
  verdict is the southern margin alone, and two mechanisms of opposite character
  were being averaged. `notes/audits/resolution-ladder.md`.
- `world-3y2v` -- T85 still has no endurance row, and it is blocked on
  `world-j1po` (no T85 staged surface family, no restart template) rather than
  on the arm's price. Buying only CONVERTED arms there repeats the two-point
  ambiguity world-90y6 had to buy a third arm to escape; a cold arm needs the
  soil chain first.
- `world-l11x` closed: `assess_lpj_run.py` is 8.6x faster with a bit-identical
  report. Parallelism was measured and DECLINED with the numbers on the row, so
  the decision can be re-taken if run sizes grow.
- `world-bga1` and `world-ov4u` -- the two directions of the worktree link
  hazard, both closed. A per-file linked file is a symlink INTO the main
  checkout: `references/pdf/` now exists so papers land in one wholly-ignored
  directory, and `link_worktree.py` keeps a ledger that reports a linked file
  whose target moved. `notes/audits/worktree-write-through.md` enumerates 594
  per-file links across 11.6 GB.

## Traps found this session worth knowing before they bite again

- **A check that can only agree with itself.** Four separate agents found
  defects in their OWN work by trying to buy a measurement rather than by
  reading, and the recurring shape is a test holding its own copy of the number
  it is testing. The sharpest: a save-point self-test that restated the rule
  instead of calling the function, and so passed over a state file written
  thousands of simulated years early.
- **`date.year` counts through the spin-up.** A run of `--nyear 1253` behind
  `nyear_spinup 8600` ends at simulated year 9852. Anything computing a year
  from `--nyear` alone is wrong by the spin-up, silently.
- **Rule 3 in the form that does not announce itself.** The land column states
  file keys rows by a CLIMATOLOGY's longitude convention while the export's
  `lon.bin` calls the same column something else. Matching them does not match
  zero cells -- it matches ALL of them half a planet away, agreeing on 0.5955.
  `lib/gridding.py` now owns that label axis and refuses a label not on it.
- **A gate must not RUN what it only needs to START.** 27 of 212 entry points
  build no argparse parser, so `--help` was an argument nothing read and
  starting them once rewrote nine tracked artifacts.
  `docs/src/practice/failure-modes.md` class 38.
- **`pgrep -f <pattern>` matches your own command line.** It cost this session
  twelve minutes of believing a finished profile was still running. Use
  `pgrep -x`, or wait on an artifact.
- **Check a row's existing labels before adding one.** Three rows were
  double-labelled today and the gate caught each; `bd export` also lags the
  live database, so a batch-label failure is sometimes only a stale export.
