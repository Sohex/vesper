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
`dayseccheck`, `nlsetkey`, `tpcore`, `ytp`, `aerocore`, and all of
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

## The fourteenth site, which the first pass could not see

Re-deriving the population under a gate found one more, and the reason it was
missed is a property of the search rather than of the site.
`hurricanemod`'s `entropy_deficit` ends on

    if (xhi<0) print *,smr,R(NLEP),sms-sme,smo-smb,smrt,R(i600)

A `print` statement writes to unit 6, which is `nud`, so it reaches
`plasim_diag` exactly as `write(*,...)` does. The first pass enumerated
`write(nud,...)` and `write(*,...)` and no `print`, so this site was never in
the 928. `entropy_deficit` is called from `hurricanestep`'s `do jhor=1,NHOR`
over the thread's own chunk, and `hurricanestep` is called from `plasim.f90`
with nothing between it and the parallel region, so every thread reached it.

It falls in the second class above: the numbers are the thread's own gridcell,
so a root guard would silence exactly the cells worth seeing. It is serialised
with `!$omp critical (nudwrite)` and the thread id goes in the record, matching
`surflx`'s per-cell reports.

Two more sites were guarded correctly and unreadably.

`mpstart` wrote its thread-count mismatch message under `if (mypid == 0)`.
`NROOT` is `parameter(NROOT = 0)`, so the branch taken was already the right
one; the spelling was the only site in `plasim/src` that did not say `NROOT`,
and it now does. The branch is unchanged.

`calini` writes its calendar summary under `if (kpid == 0)`, where `kpid` is a
dummy argument. It is not a thread test at all: `prolog` passes -1 and `readnl`
passes 0, so the block is a print-on-the-second-call flag and is reached from
`readnl` alone. `readnl` is root-only by caller, so one thread reaches it. That
is the argument, and it rests on the two call sites rather than on anything
visible in `calmod`.

## What holds the convention now

`exoplasim/scripts/lint_diag_writes.py` is the two passes above written down,
and `scripts/smoke_test.py` runs it. It permits four answers and nothing else:
an enclosing `mypid == NROOT`, an enclosing `npro == 1`, an enclosing
`!$omp critical`, and an enclosing procedure the call graph shows is reached
only through those. It reads the boolean structure of a condition at depth zero
rather than searching its text, because `.and.` binds tighter than `.or.` in
Fortran and the two need opposite treatment, and it narrows the caller
population when a branch is pinned to an integer literal through a dummy
argument, which is how `calini`'s summary is cleared without a table row.

Its exception table is empty. Every write in `plasim/src` today is answered by
one of the four rules, which is the state worth holding: a table row is an
argument a future reader has to re-check, and a rule is one the parser checks
for them.

The pass answers twenty reduced fixtures on every invocation before it reports
on the tree, because a gate whose only evidence is that the tree is clean
cannot distinguish a clean tree from a pass that reports nothing. Removing the
`.and. mypid == NROOT` from `landini`'s LSHY-3 message reports its three
writes; removing the `!$omp critical` from `surflx` reports its four; an
unguarded `write(nud,...)` dropped into `surflx` is reported where it stands.
