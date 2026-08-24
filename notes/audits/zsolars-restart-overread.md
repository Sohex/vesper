# The restart writes 8190 elements of buffer for a 2-element array

Measured 2026-08-20, on the T42 L10 p16 build, found by a compile-flag benchmark
rather than looked for.

Worldbuilding frame: this is about the Vesper simulation's climate model and the
files it checkpoints itself into. Nothing here is about the simulated world.

## What is true

`radstop`, the routine that writes radiation state into the restart, saved the
two solar constants with

    call mpputgp('zsolars',zsolars,2,1)

`zsolars` is declared `real :: zsolars(2)` at `radmod.f90:354`. `mpputgp`
(now `mpimod_omp.f90:680`) declares its dummy `p(kdim,klev)`, allocates a LOCAL
`z(NUGP,klev)`, gathers into it with `mpgagp`, and writes all of `z`:

    real :: z(NUGP,klev)
    call mpgagp(z,p,klev)
    if (mypid == NROOT) call put_restart_array(yn,z,NUGP,NUGP,klev)

`mpgagp` moves `NHOR` elements per rank, so at T42 p16 each of the sixteen ranks
contributes 512 elements read from an array that holds 2. The record that
reaches the file is `NUGP` = 8192 doubles, of which two are solar constants and
8190 are whatever the gather buffer and the memory past `zsolars` happened to
hold.

That is `docs/src/practice/failure-modes.md` class 18, and it is the FIRST
INSTANCE IN THE OPPOSITE DIRECTION. The documented one moves too FEW elements,
so a field comes back short and the counter beside it survives intact. This one
moves too MANY, so nothing comes back wrong and the surplus is written out.

## What it costs, which is not the physics

Nothing reads the RECORD. `get_restart_array("zsolars",...)` is commented out at
`radmod.f90:1177`, and `solarini` recomputes the pair and broadcasts it
(`radmod.f90:865`). So the integration is untouched, and no run on disk is
wrong because of this.

The pair ITSELF is not a diagnostic and was not one then: `swr` forms the
band-weighted surface albedo from it at `nstartemp = 1`
(`radmod.f90:2918`), so `zsolars` is live model state every timestep. What
nothing reads is the restart copy.

What it damages is RESTART IDENTITY. A restart's bytes are the restart's name in
this project: `continue_exoplasim.py` records `input_restart_sha256` on every
segment and CLIM-31's `seed_copy_sha256` names a donor by content. Those hashes
are stable today only because `-finit-real=zero` is in the shipped flags and
holds the gather buffer at zero. They are stable by a compiler option, not by
construction, and the option is one that a performance change would plausibly
touch.

## How it surfaced

The build-flag benchmark's `-flto` arm returned a DIFFERENT restart sha on every
run, five runs and five shas, from one binary on one input. That looked like
non-determinism in the model and was not.

Comparing the restarts record by record: 397 of 398 records were bit-identical,
both between two `-flto` runs and between `-flto` and the stock build. The
integration is reproducible under `-flto`. The single differing record was
`zsolars`, and within its 8192 elements exactly one differed:

| | element 11 |
| --- | --- |
| stock build | 0 |
| `-flto` run 1 | 0.168103 |
| `-flto` run 2 | 0.255179 |

`-flto` lets the compiler stop zeroing a local buffer it can prove is never
read, so the padding stops being a stable zero and starts being memory.

## The sweep

Every `mpputgp` and `mpgetgp` call site in `plasim/src` was resolved to its
argument's declaration and compared against the `(NHOR,klev)` shape those
routines assume. **459 call sites, one genuine mismatch, and it is this one.**
The 24 sites in `simba.f90` that a shape match flags are `allocatable` and are
allocated `NHOR` at `simba.f90:171`, so they are correct.

## The fix was already in the same file

`zsolars` is a global pair, not a distributed gridpoint field, so it does not
belong in a gather at all. `solarini` wrote it correctly eleven hundred lines
earlier, with `put_restart_array("zsolars",zsolars,2,2,1)`. One call was the
right idiom and the other was not, in one file, for one array.

## The fix, and what was checked before and after it

Applied at `radmod.f90:1782`, replacing the gather with
`if (mypid == NROOT) call put_restart_array('zsolars',zsolars,2,2,1)`.

Checked BEFORE changing it, because the low-I/O change made restart layout
something a resume depends on: `reseek` (`restartmod.f90`) scans the file
reading a name record and then doing a bare `read` to skip the data record, so
records are located by NAME and skipped WHOLE. A record's size cannot dislodge
any other record, and nothing reads `zsolars` back in any case.

Checked after, on the same bed and the same restart:

- the record is 16 bytes where it was 65536, and still carries 0.38237199 and
  0.61762801, the two solar constants
- the restart still holds 398 records, and 397 of them are BYTE-IDENTICAL to the
  stock build's. The only record that moved is the one that was meant to. The
  fix is physics-neutral by measurement, not by argument
- restarts cross the change in both directions: the fixed binary resumed a
  restart written by the stock one, and the stock binary resumed a restart
  written by the fixed one, both cleanly and with no missing-record message

## A second defect in the same nine lines

`nwriunit` is unit 34, and `solarini`'s nine `put_restart_array` calls run at
INITIALISATION, when no restart file is open on that unit. Fortran therefore
connects `fort.34` and writes them there. Every run directory carried one:
`run_4182235e9781/fort.34` is 432 bytes and holds exactly those nine name and
data pairs, `zsolars` through `doceanalb`.

**This is a correctness defect and not a design question**, which is the
opposite of what it looks like at first. Whether those nine arrays SHOULD be
checkpointed is a design question, and the source already answers it: the
matching READ block is commented out in full at `radmod.f90:990-1032` -- nine
`get_restart_array` calls, their broadcasts, and the `nrestart > 0` branch
around them -- and `solarini` was made unconditional in its place, so the arrays
are recomputed from the namelist at every start. Someone decided these are not
restart state. What was left behind is the WRITE half of the round trip they
disabled, pointed at a unit nothing opens.

Both defects are present in `alphaparrot/ExoPlaSim` master verbatim, checked
before anything was offered, and both went upstream as PR #63.

So there is no reading under which the nine calls were doing their job. Removed,
with the reason recorded at the site. Verified in isolation, at stock flags so
the `-march=znver4` change could not hide it: `fort.34` is no longer created,
the restart still holds 398 records, and it is BYTE-IDENTICAL to the build
carrying only the `zsolars` fix. The removal is exactly neutral.

## What the surviving write is FOR, and what the converter does with it

world-5rq settled the question this record left implicit. The write is a
CONFIGURATION FINGERPRINT and not checkpointed state, and `radstop` says so at
`radmod.f90:1764-1781`: it is the only place a restart handed on without its run
directory records which two-band split of the stellar constant the orbits were
integrated with, and a restart whose `zsolars` disagrees with the configuration
it is resumed under is a different star.

`exoplasim/scripts/restart_schema.py` carries the matching conversion policy and
gives the same reason: `zsolars` is TARGET and deliberately not REQUIRE_EQUAL,
because a converter that refused over it would be refusing on a value the model
is about to discard and recompute.
