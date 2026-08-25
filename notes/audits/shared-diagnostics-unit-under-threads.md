# The diagnostics unit is one unit for the whole thread team, and the guard convention had thirteen holes

Established 2026-08-25 by static analysis of
`vendor/exoplasim/exoplasim/plasim/src`, against the tree at the time. No model
run was involved and none could settle it: the question is which statements a
thread other than 0 can reach, and that is a property of the source.

Worldbuilding frame: this is about the Vesper climate model's own diagnostics
file, `plasim_diag`. Nothing here concerns the simulated planet.

## The mechanism

`opendiag` (`plasim.f90`) opens `nud` under `if (mypid == NROOT)`, and `nud` is
6 (`plasimmod.f90`). A Fortran unit belongs to the PROCESS, so that one `open`
makes unit 6 the diagnostics file for every thread, and `write(*,...)` lands
there too because `*` is unit 6 as well.

The whole run is one parallel region: `plasim.f90` opens
`!$omp parallel num_threads(NPRO)` around `mpstart` through `mpstop`. So every
statement in the model's call tree is executed by NPRO threads unless something
on the path tests `mypid == NROOT`. Under the deleted MPI arm a non-root rank's
`write(nud,...)` went to that rank's own stdout and out of the record, so an
unguarded site was invisible; under threads it lands in `plasim_diag` beside
root's.

libgfortran locks a unit for the duration of one data transfer statement, so a
single unguarded `write` is reordered rather than torn. That lock does not span
statements, and it does not make `open`/`write`/`close` on one unit safe from
several threads.

## Method

Two passes, both purely lexical over the 24 sources that write to `nud`:

1. For each of the 928 `write(nud,...)` sites, walk the enclosing `if` / `else`
   / `select` nesting and record whether any enclosing branch tests
   `mypid == NROOT`. 214 sites were not dominated by one, in 36 subroutines.
2. Build the call graph from `call` statements with the same guard walk, and
   take the fixed point of "every call site is either under a root guard or in
   a subroutine that is itself root-only". That reduced 36 subroutines to 13.

The pass over-reports rather than under-reports: `end where` is
indistinguishable from `end subroutine` to it and clears the guard stack, which
turns a guarded site into a reported one and never the reverse. Every reported
site was then read in context.

## What the 13 were, and the three classes they fall into

**Root-only by caller, and correct as they stood.** `print_planet`, `initpm`,
`initsi`, `readnl`, `setzt`, `noise`, `printseed`, `wrorb`, `wrspam`, `wrzs`,
`dayseccheck`, `nlreplace`, `tpcore`, `ytp`, `aerocore`, and all of
`restartmod`'s record readers and `surfmod`'s `get_surf_array`,
`check_surf_header` and `surface_ini`. These are reached only from inside a
`if (mypid == NROOT)` block, so the guard is on the path rather than at the
statement.

**Guarded by `npro == 1`, and correct as they stood.** The orography summaries
in `surfini` and `glacierini`. `NPRO` is the thread count, so that test admits
exactly the single-thread case and there is no second writer to race.

**Genuinely reachable by every thread.** Thirteen sites, in three kinds, and the
kind decides the fix:

| Site | Kind | Fix |
| --- | --- | --- |
| `gridpointd` "Semi-Lagrangian q running" and its `flush(nud)` | global fact, every timestep | `mypid == NROOT` |
| `spectrald` "### Damping with" | global fact, every timestep until the countdown ends | `mypid == NROOT` |
| `master` "HC OUTPUT step" | global fact, every high-cadence step | `mypid == NROOT` |
| `master` "HC STORM CAPTURE step" | global fact | `mypid == NROOT` |
| `carbonstep` CPU-time line under `ntime > 0` | per-thread number, reported like a global one | `mypid == NROOT`, matching `fluxstop`, `rainstop`, `radstop` and `miscstop` |
| `surflx` three `zbz0` cell reports | the THREAD's own gridcell | `!$omp critical`, and the thread id put in the record |
| `surflx` wrong-`ntsa` refusal | uniform refusal | root prints, every thread stops |
| `glacierprep` `glacpersist` refusal | uniform refusal | root prints, every thread stops |
| `landini` `tau_veg`/`tau_soil` refusal | uniform refusal | root prints, every thread stops |
| `vegini` NVEG/NSTARTEMP refusal | uniform refusal | root prints, every thread stops |
| `calini` 365-day calendar refusal | uniform refusal, in a module with no thread identity | `!$omp critical` |
| `fftini` unsupported-resolution refusal, both FFT variants | same, and `trigs` is threadprivate so each thread initialises its own | `!$omp critical` |
| `mpabort` `Abort_Message` block and its `nud` lines | `open`/`write`/`close` of unit 44 from any thread that aborts | `!$omp critical`, thread id in the record |

`mpabort` is the one whose failure was not merely a reordering. Unit 44 is
opened, written three times and closed; one thread's `close` landing between
another's `open` and its writes disconnects the unit under it. `Abort_Message`
is load-bearing outside the model -- `exoplasim/__init__.py` and the profiling
shells test for its existence -- so the file must still be written by whichever
thread aborts, which is why this is a critical section and not a root guard.

## The rule a new site follows

Three questions, in order. Is the quantity the same on every thread? Then
`mypid == NROOT`. Is it the thread's own chunk? Then `!$omp critical` with
`mypid` in the record, because a root guard would silence every thread but 0 and
those are the cells worth seeing. Is the module a leaf with no `pumamod`
dependency, so it has no `mypid` to test? Then `!$omp critical`, rather than
giving the module a dependency on the parallel layer.

A refusal takes the shape `landini` and the LSHY-3 checks already use: the
message under `if (mypid == NROOT)`, the `stop` outside it, so every thread
stops and one copy of the message is written.
