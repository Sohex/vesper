# Handoff

Written 2026-09-07 at the close of a fan-out session: 231 commits, 24 merges,
55 rows closed. This file is a MAP for the next session: state, the decisions
waiting, and the context that would otherwise be re-derived. Everything here
cites a row, a note or a file; nothing here is the record itself. Delete or
rewrite this file when its content is consumed -- it describes one moment.

## State at handoff

- Tree clean. `smoke_test.py` 83/0, `verify_entry_points.py` 222/222, model
  source compiles, all ten ExoPlaSim binaries and the LPJ-GUESS binary verify
  against the source they were built from.
- `check_consistency.py`: 65 checks, 4 failures, and none of the four is a
  defect anybody needs to chase:
  - **two are carve2-transition staleness** on runs rule 7 already calls
    disposable -- namelist float round-trip and staged surface fields on
    superseded runs.
  - **two are the flux bracket refusing to certify itself.** The arms exist and
    the slope measured cleanly; the declaration is deliberately not moved onto
    them. `world-jejw` is why, and it is the first thing to read below.
- `bd ready` is the work list, 135 rows. This file annotates only what is
  expensive to rebuild.

## The one defect that matters most, and it is an instrument

**`world-jejw` (P1): the offset criterion is decided by a fit whose own tau is
not determined.** `assess_convergence.py` selects its statistic on
`fit_usable`, which never consults the fit's own tau standard error even though
it computes and records it. The consequence is inverted: a fit that fails
COMPLETELY routes to the drift fallback and PASSES, while one that half
succeeds FAILS.

It is demonstrated rather than suspected. Two arms differing only in flux ratio
were bought on the configured build -- `run_c919391cf715` at 0.945 and
`run_9d5dbf9bd3d9` at 1.000, 130 orbits each, matched in build, geography,
executable, window and I/O regime -- and **there is no window on which both
pass**: cold passes at 15-69 and fails at 75-129, warm the reverse, every other
criterion passing on all four assessments. Two runs that differ only in flux
cannot have a settling verdict that flips with the window.

The slope itself is fine: **173.8 K per unit flux ratio**, asymptotes 173.78,
window means 173.22, last-ten 173.05, agreeing to 0.73, and two of the three
estimators use no fit at all. What is missing is a criterion that can certify
it. Fixing this changes verdicts tree-wide and needs a declared sigma
threshold, which is why the agent that found it did not take it: the repair
unblocks its own row.

Also touching this: `world-nmtp`, and the note in `lib/rungs.py` that
`run_14906cb7b914`'s fitted relaxation is degenerate and supports nothing, are
the same estimator and the same blind spot.

## The acceptance run, which is specified and not yet bought

**Nothing blocks it now.** `world-24xd` is resolved -- water closes stock
against flux like carbon and nitrogen, via a new retained table `awater.out`,
and a fresh short run refused 0 of 1617 gridcells at every window from 1 to 29
cycles. `world-glu7` is resolved and the restart is continuous on main's own
binary, 23 tables, `continuous: true`. The binary carries every merged change.

`world-qcse` carries the specification: 8600 spin-up, 1253 retained, npatch 5,
16 ranks, `--save-state`, about 2.3 hours. Pass `--save-state`: it is ungated,
changes no number in the run that writes it, and `world-gqyp` is now resolved,
so a state file is continuable rather than dead weight.

## What the fan-out changed under you

Read this before trusting any pre-2026-09-06 measurement of these:

- **The water closure was holding a store to a conservation tolerance.** The
  residual telescopes -- 5.20x over a thirtyfold window against 5.48 predicted
  for a store and 30 for a missing flux -- and the burnt-fraction correlation
  went +0.66 to -0.028. Refusals FELL as the window lengthened.
- **Every restart this project ever took was an arbitrary-day restart.**
  `plib.cpp` stored an integer parameter as `(int)(num + 0.5)`, which rounds
  negatives toward zero, so the `-1` year-boundary sentinel arrived as 0.
  Failure-modes class 40.
- **The photon currency is this star's**: `VESPER_CQ` 4.864175e-06 mol/J,
  +5.7 per cent, and `FRADPAR` 0.4624 -> 0.4913 -- the second because the solar
  reference was integrated over a truncated spectrum.
- **Soil stoichiometry moved**: `NMASS_SAT` 0.05 -> 0.002 kgN/m2, the C:N ramps
  read NO3 + NH4 rather than NH4 alone, microbial C:P 80 -> 63.31.
- **The area weight is the Gauss-Legendre quadrature**, 28 sites. Nothing moved
  past its own scatter, but every generated JSON in that path changes in its
  fifth or sixth significant figure.
- **The acceptance contract is `vesper-lpj-equilibrium-window/6`**: four
  quantities left the assessed set, which removes the only refusal standing on
  the record on disk.
- **The run index's payload-gone flags were false.** `register()` read a row's
  absence as a deletion, so registering from a worktree marked 33 of 34 runs
  gone. `upsert()` replaced it.

## Decisions taken, so they are not re-argued

- **The broadband radiation scheme stays.** The band-resolved candidate
  measured 8.6x per column on the corner this model runs, trebling a T85
  commissioning from 5.2-10.6 hours to 16-33.
  `docs/src/reference/design-intent.md` carries it. The three holes it would
  have closed are now rows in their own right -- `world-2kp9`, `world-skiq`,
  `world-6fjn` -- and the route to them is gap fill inside the present scheme.
  A future proposal to swap has to beat the measurement rather than restate the
  holes; `world-0603` carries the ratio it would need.
- **`surf_ocn_sic` stays; `config2xml.py` and `dzu` are deleted.** The
  distinction that decided all three is whether the deadness is mainline cGENIE
  or Vesper-created.

## Rows whose context is expensive to rebuild

- `world-ccx6` closed: the soil takes `--grid` and the climatology declaration
  answers per rung, refusing a rung with nothing declared rather than serving
  another. `world-j1po` is now blocked only on there being no climatology above
  T21 -- not on a carrier.
- `world-s6wz` (P2): a worktree that satisfies rule 4 commits a manifest
  describing binaries only it has. Cleared here by rebuilding in main, but the
  gap is real and the row weighs whether the manifest should be tracked at all.
- `world-giqf` (P2): 285 declared members outside `Soil` are unstreamed and
  ungated. `Soil` has a gate; no other class does.
- `world-0xpg` (P2): 17 gridcells sit at the 10000 mm snowpack cap in all 30
  retained cycles, where the melt term is zero and snowfall reaches the soil as
  liquid water below freezing. Conserving, which is why the closure passes, and
  not what the surface is.
- `world-x0u5`: the canopy snow capacity is 19 to 31 times smaller than
  Hedstrom and Pomeroy measure. The row's question changed on reading the
  paper -- `forcap` is being asked to be an optical saturation scale AND a mass
  capacity, and those are two quantities. Free to fix properly: the store is
  inert until bio-33 and `forintc` ships at zero.
- `world-g5xe`: 195 tracked JSON files carry absolute paths and nine name
  worktree directories that no longer exist.
- `wet-12` declared the reduced wetland form, and the recorded one was
  circular: it named the withdrawn saturated fraction as the part that forfeits
  least. Extent has no floor, and the reduced form is a smaller partition rather
  than a coarser one.

## Traps this session paid for

- **A check that can only agree with itself**, again. `stability_check` wrote an
  abort message then a bare Fortran `stop`, which exits 0 -- so aborted and
  finished were the same status, and a 200-step bed was measuring its own
  length. Failure-modes class 38, in a new place.
- **A count of one kind of operation is not a cost.** An instruction count said
  the band-resolved scheme was cheaper at every corner; measured, it is 8.6x.
  The excluded terms were most of the cost. Second time on that row.
- **A timing with no rung attached** cost a 7x mispricing. Every per-orbit
  figure in the tree now carries its rung, its rank count and its date.
  `native_runtime_seconds` is CPU-seconds, not wall: about 14x the wall span on
  T21 at 16 threads.
- **An ordering keyword, not an algorithm**, was the groundwater cost:
  `permc_spec="MMD_AT_PLUS_A"` took 154 s on one block and 11 s on a larger one.
  Removing it took a 35-minute, 29 GB bound to 301 s at 16 GB.
- **A generated input is upstream of the binary that compiles against it.**
  Rule 4's rebuild is not the first step of a close-out: the LPJ build failed
  because `vesper.h` had not been regenerated ahead of it.
- **Held material is invisible unless something indexes it.** Four pedology
  PDFs sat on disk with no `INDEX.md` row while three rows recorded them as
  unreachable; two source trees the same; and a worktree's per-file-linked
  directory nearly took the flux bracket's convergence reports and did take the
  T21 wetness field. Run `link_worktree.py --check` before a worktree dies.
- **Zenodo 403s this host** on IP reputation, on open-access records. A browser
  or different egress gets past it.
