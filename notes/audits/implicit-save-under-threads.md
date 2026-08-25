# Every implicitly SAVEd local in the model is written before it is read, and two were dead

Established 2026-08-25 by static analysis of
`vendor/exoplasim/exoplasim/plasim/src`, against the tree at the time. No model
run was involved and none could settle it: the question is which storage class a
declaration acquires and which reference to a variable comes first, and both are
properties of the source. The behavioural evidence that the transformation this
note audits preserved results is already in the tree, at commit `5f7febfb`: T21
against the MPI build on one bed, 199 of 199 records bit identical at two
threads.

Worldbuilding frame: this is about how the Vesper climate model's Fortran
declares its variables and what that means under threads. Nothing here concerns
the simulated planet.

## The mechanism

A local declared with an initialiser inside a procedure acquires the SAVE
attribute. SAVE means static storage, and under `-fopenmp` static storage is ONE
copy for the whole thread team. `plasim.f90` opens a single parallel region
around `mpstart` through `mpstop`, so every procedure in the model's call tree
executes on NPRO threads and every SAVEd local in it is shared.

`notes/audits/uninitialised-reads-and-implicit-save.md` measured the race
directly -- eight threads writing their own id into a SAVEd scratch array and
reading it back, seven of the eight reading a value another thread wrote -- and
established that gfortran reports nothing about the class at any warning level
tried. `-frecursive` does not move such a variable to the stack, because the
standard says it is SAVEd.

Under the MPI arm each rank was a process and its SAVEd locals were private by
construction, so the identical source was correct. world-38b removed that arm.
There is no configuration left in which this class is harmless.

## Method, and the direction it errs in

`exoplasim/scripts/lint_implicit_save.py` walks procedure nesting over the 38
translation units and reports every declaration inside a procedure body that
acquires static storage: a type declaration carrying `::` and an initialiser, an
explicit `save` attribute or statement, and a `data` statement. `parameter`
declarations and module-level variables are excluded, because their static
storage is the point.

**It over-reports.** A declaration it cannot attribute to a procedure is
reported rather than dropped; a `program` unit is treated as a procedure body
even though its variables are static by definition; and an `!$omp
threadprivate` directive answers a site only when it is in the SAME procedure,
so a directive in a module specification part does not clear a local of the same
name. The one shape it does not look for is a declaration without `::` that
carries an initialiser, because Fortran has no such shape.

The population is 66 sites. The verdict for each lives in the script's
`CLASSIFIED` table and in the `threadprivate` directives in the source, not in
this document, so that a declaration which CHANGES falls out of its answer and
is reported again.

## The classification

**Forty are threadprivate, and that was the thread port's blanket answer rather
than a per-site one.** `5f7febfb` transformed every initialised local it found,
which was the right move for a port that had to reproduce the MPI build. It left
open the question this note answers: whether `threadprivate` is the correct
answer at each site, or whether it silently split a variable that was meant to
persist globally across calls.

It did not. **Every one of the 66 has a WRITE as its first reference in the
procedure, or is never written at all.** That was checked by walking each
procedure body from its declaration and finding the first statement that
mentions the name. The consequence is that no site carries a value from one call
to the next in a way any caller observes, so `threadprivate` cannot have changed
a result -- each thread's copy is overwritten before it is read, exactly as the
shared copy was. Under the original single-process build the same property held
for a different reason, which is why the source was correct under ranks.

The forty fall into two shapes, and the shape does not change the verdict:
per-cell scratch sized `NHOR`, `NUGP` or `NLEV` (`icemod`'s snow and conductive
flux temporaries, `oceanmod`'s tridiagonal coefficients, `seamod`'s ice fraction
and roughness, `carbonmod`'s and `glaciermod`'s gather buffers, `hurricanemod`'s
column profiles), and search indices that a level loop sets before it uses
(`i850`, `i600`, `i200`).

**Nine data-statement tables and two more are read-only constants, and one
shared copy is what every thread wants.** `radmod`'s `orb_params` carries the
nine Berger orbital series -- amplitudes, rates and phases for obliquity,
eccentricity and the moving vernal equinox -- each written once by its `data`
statement and read nowhere but the three series sums. `calmod`'s month
abbreviations and `fft991mod`'s factor list are the same shape. A shared
read-only table cannot race, and converting nine long coefficient tables to
`parameter` would be transcription risk for no correctness gain, so they stand
as they are.

**Thirteen belong to two standalone programs.** `buildice.f90` and
`newsnow.f90` are `program` units, not in `plasim/CMakeLists.txt`'s source list
and not in `config/pipeline.yaml`; `notes/audits/dormant-exoplasim-modules.md`
holds the facility and CLIM-53 owns the verdict on whether they get wired. A
main program's variables have static storage whatever their declaration, and one
of these is one thread, so the class does not reach them.

**Two are safe because of their CALLER, and are the only two like that.**
`trc_routines`'s `filns` caches a grid-derived constant behind a `first` flag,
which is the shape that races under threads: one thread clears the flag before
another tests it, and the second reads the cache before it is filled. It does
not race here because `tracermod`'s `tracer_main` calls `tpcore` inside
`if (mypid == NROOT)`, so one thread reaches it. That matches the root-only-by-
caller class `notes/audits/shared-diagnostics-unit-under-threads.md` established
for the same call chain. It is worth naming separately because it is the only
site in the population whose safety cannot be read off the procedure it is in.

**Two were dead, and deleting them was the fix.** `icemod`'s `addfci` declared
`zsum(2)`, referenced only from commented-out lines, and `hurricanemod`'s `gpot`
declared an `i850` it never used. Both had been given `threadprivate` directives
by the blanket transformation, which is a thread-local copy of nothing.

**Four were constants declared as variables, which is what put them in the class
at all.** `landmod`'s `getalb` held `aa = 5.2` and `yy = 4.0`, `icemod`'s
`make_ice_thickness` a fourteen-element seasonal factor table, and `surfmod`'s
`surface_ini` a format string, each never assigned and each carrying a
`threadprivate` directive. They are `parameter` now, which removes the storage
question rather than answering it and drops them out of the population.

## What this does not say

Zero defects found is not zero risk. The property that makes all 66 safe --
first reference is a write -- is a property of the current code and not of the
declaration, so it is re-established by the pass and not by the source. A new
initialised local added tomorrow would be reported by this pass and by nothing
else, which is why it is registered as a check rather than run once. The
`CLASSIFIED` table is keyed on the variable name, so a rename or a new site is
reported; it is not keyed on the body, so a site that stops being
write-before-read keeps its row. That is the gap, and closing it needs a
dataflow pass rather than a lexical one.
