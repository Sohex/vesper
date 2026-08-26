# `epilog` deallocates `adenergy` and then writes it to the restart

*Worldbuilding frame: this is a defect in the Fortran of the Vesper project's
climate model. Nothing here is about the simulated planet.*

**Established 2026-08-26 on model source `a041a1e9`, against the binaries
`exoplasim/binary_manifest.json` registers.** Every run this tree prepares
carries `NENERGY = 1`, and every such run takes SIGSEGV at the end and leaves a
truncated `plasim_status`.

## The defect, in two lines of one subroutine

`vendor/exoplasim/exoplasim/plasim/src/plasim.f90`, `subroutine epilog`, which
begins at line 922:

    934:  if(nenergy   > 0) deallocate(denergy)
    935:  if(nenergy   > 0) deallocate(adenergy)
    ...
    1148: if (nenergy > 0) call mpputgp('adenergy',adenergy,NHOR,28)

The array is freed at 935 and read at 1148, under the same guard, so the
use-after-free fires exactly when the energy diagnostics are on and never
otherwise. Line 1148 is the newer of the two: it belongs to the change that put
the energy diagnostics' partial window into the restart, and it was added below
a deallocation that already stood.

## What it does

`mpputgp` hands the freed array to `mpgagp`, which reads it. The fault is a
segmentation fault rather than a wrong number, so nothing silently continues:

    #3  mpgagp_
    #4  mpputgp_
    #5  epilog_
    #6  MAIN__._omp_fn.0

The run has already finished integrating when this happens. What is lost is the
restart: `plasim_status` is written up to the point of the fault and ends one
record marker short, so `compare_restarts.py` and `restart_format.py` refuse it
as truncated and the run cannot be continued.

## The evidence, one variable at a time

Every arm below is a cold 200-step run except where noted, on the registry's own
executables, and the beds differ from each other in the named key alone.

| arm | rung | threads | flags | `NENERGY` | outcome |
| --- | --- | ---: | --- | ---: | --- |
| prepared run directory, copied whole | T21 | 16 | -O2 | 1 | **SIGSEGV in `epilog`** |
| staged bed | T21 | 16 | -O2 | 1 | **SIGSEGV in `epilog`** |
| staged bed | T21 | 8 | -O2 | 1 | **SIGSEGV in `epilog`** |
| staged bed | T42 | 16 | -O2 | 1 | **SIGSEGV in `epilog`** |
| staged bed, one full orbit | T42 | 16 | -O3 | 1 | **SIGSEGV in `epilog`** |
| the same T42 bed, one key changed | T42 | 16 | -O2 | **0** | completes, returncode 0 |

So it is not resolution-dependent, not thread-count-dependent, not
optimisation-dependent, and not a property of a hand-built bed: the control is
the last row, which is the identical directory with `NENERGY` and `NENER3D`
taken to zero and nothing else touched.

**Why no probe caught it.** `exoplasim/scripts/stability_probe.py` writes
`NENERGY = 0` and `NENER3D = 0` into every bed it builds, so the whole stability
grid runs the branch that does not fault. The probe is not wrong to do that --
it is switching off a diagnostic it does not read -- but it means the grid says
nothing about this path.

## The fix, and what it costs

Move the two deallocations below the restart write, or delete them. The program
is ending and the process is about to exit, so nothing depends on the memory
being returned; the deallocations exist for tidiness and the restart write is
the load-bearing statement. `denergy` at 934 is not read afterwards and is safe
either way, but the pair should move together so the next reader does not have
to work out which of them mattered.

Rule 4 applies: this is a change under `vendor/exoplasim`, so every binary is
rebuilt after it and `rebuild_binaries.py --verify` is what says the registry
agrees.

## What is now worthless, and what is not

Every run on this source that reached `epilog` with `NENERGY = 1` has a
truncated restart and cannot be continued. Its OUTPUT is unaffected: the fault
is after the integration and after the output stream is closed, so records
already written are sound, and the ocean stream `oceanmod` writes during the run
is likewise complete up to the last interval.
