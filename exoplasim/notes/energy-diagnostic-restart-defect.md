# The energy diagnostics are deallocated before the restart writes them

**Context first, because this document can be read standalone: this is
worldbuilding for Vesper, a fictional super-Earth, and everything below is a
property of the ExoPlaSim source that simulates its climate. Measured
2026-08-26 on model source `a041a1e9`, the source all ten registered binaries
were compiled from.**

## What is true

`plasim.f90`'s `subroutine epilog` deallocates the energy diagnostic arrays and
then writes one of them into the restart 213 lines later:

    922    subroutine epilog
    ...
    934      if(nenergy   > 0) deallocate(denergy)
    935      if(nenergy   > 0) deallocate(adenergy)
    936      if(nener3d   > 0) deallocate(dener3d)
    937      if(nener3d   > 0) deallocate(adener3d)
    ...
    1148     if (nenergy > 0) call mpputgp('adenergy',adenergy,NHOR,28)
    1149     if (nener3d > 0) call mpputgp('adener3d',adener3d,NHOR,NLEV*28)

Both statements are inside `epilog` and both are reached on the same pass. The
allocation at `plasim.f90:476` and the write agree on shape, `(NHOR,28)`, so the
dimensions are not what is wrong: the array the write dereferences has already
been returned to the heap.

The restart write was added by `8f6a591b`, which carried the diagnostics across
a restart so that the first output record after a resume stopped dividing a
partial accumulation by a whole window. The deallocate above it was already
there and was not moved.

## What it costs

**Every run with `nenergy > 0` fails after its last timestep**, in `epilog`,
gathering `adenergy`:

    #3  mpgagp_   at ./mpimod_omp.f90:218
    #4  mpputgp_  at ./mpimod_omp.f90:740
    #5  epilog_   at ./plasim.f90:1148
    #6  MAIN__._omp_fn.0 at ./plasim.f90:69

`nenergy > 0` is not an unusual setting. `config/planet.yaml` declares
`energy_diagnostics: true`, and it has to: `energy_fixer: true` is declared
beside it and `run_exoplasim.py:declare_energy_fixer` raises if the diagnostics
are absent, because the fixer is driven by `denergy26` and `denergy27`. **The
standing production configuration therefore cannot complete a single orbit.**

ExoPlaSim invokes the executable once per orbit and `epilog` runs on every
invocation, so this is an orbit boundary and not a run boundary: a multi-orbit
run dies at the end of orbit 1, and no restart is written for the next.

**It is a use-after-free and not a clean abort, so the failure mode is not
guaranteed to be a crash.** Whether the freed page is still mapped decides
between a segfault and a silent write of whatever now occupies it. A restart
carrying a corrupt `adenergy` would be read back by `plasim.f90:1590` and
divided as a partial window, which is the quantity `8f6a591b` existed to get
right. The crash is the benign outcome.

## The measurement

Six one-orbit T21 runs, `most_plasim_t21_l10_p8.x`, each branched from
`run_14906cb7b914`'s `MOST_REST.00034`, each pinned to its own eight cores, host
load 5.0 to 7.4 over 32 cores. Everything is held except the two energy keys.

| run | configuration | `nenergy` | outcome |
| --- | --- | ---: | --- |
| `run_d3606ce0d265` | dry adiabatic, fixer off, diagnostics 2 | 2 | SIGSEGV in `epilog`, 8 threads |
| `run_c24776f33d25` | the same again | 2 | SIGSEGV in `epilog`, 8 threads |
| `run_383b871dcfad` | the same a third time | 2 | SIGSEGV in `epilog`, 8 threads |
| `run_614579b8b44a` | dry adiabatic, fixer ON, diagnostics 2 | 2 | SIGSEGV in `epilog`, 8 threads |
| `run_2260d97ffc27` | dry adiabatic, fixer off, diagnostics off | 0 | **completed, exit 0** |
| `run_352f6a4180e5` | full physics at `config/planet.yaml`'s own settings | 1 | SIGSEGV in `epilog`, 8 threads |

Every `nenergy > 0` row fails and the one `nenergy = 0` row does not.

`run_2260d97ffc27` is the control and is what makes this an isolation rather
than an observation: it differs from `run_d3606ce0d265` in `energy_diagnostics`
and in nothing else, and it is the only row that finishes.
`run_352f6a4180e5` is the standing production configuration, full physics with
the fixer on, and it fails identically at the same source line at `nenergy = 1`.

## What the diagnostics are still good for

The failure is at the restart write, after the last timestep and after
`plasim_diag` is complete. **A crashed run's `plasim_diag` carries the whole
orbit**, `CONVDECOMP` and `ENERGY FIXER` lines included, so a one-orbit
diagnostic arm still returns its measurement. What is lost is the restart, which
is what makes multi-orbit work impossible rather than merely awkward.
`exoplasim/analysis/arms/read_conversion_decomposition.py` reads those lines out
of a crashed run for that reason.

## What would settle it

Moving the four deallocate statements below the restart write, or deleting them
and letting the program exit reclaim the memory. The arrays are freed at the end
of the run and nothing between the two points allocates against them, so the
deallocate buys nothing that the process exit does not.

Any such change is a `vendor/exoplasim` change and stales all ten registered
binaries by rule 4, so it lands with a full rebuild and before any arm that has
to survive an orbit boundary.

## A second defect in the same control: `Ct` and `Dt` did not reproduce

`nenergy = 2` adds the conversion decomposition `world-0ov` is diagnosed from.
Three runs of one configuration on one restart with one binary agreed bit for
bit on `d02`, `d26`, `d27`, `cimp`, `cvadv`, `ctm` and `ctp`, and disagreed on
`ct` and `dt`. Those two terms are the only ones evaluated on `zcnow`, which
each thread copied whole from the shared divergence while every thread advanced
its own slice of that divergence later in the same routine with nothing between
the two. `world-tqh4`; the mechanism, the barrier that closes the copy against
the advance and the reproduction evidence are in
`exoplasim/notes/convdecomp-reproducibility.md`.

The two defects are independent and only the epilog one is in this note's
subject. The decomposition race never reached the state -- the racy runs and the
clean ones write the same `plasim_status` -- so it is invisible to a gate that
compares restarts, where the epilog defect kills the run outright.

## The gate sees it too

`scripts/check_consistency.py` reports the five crashed runs under "runs that
staged no namelists": they have a manifest and no `plasim_namelist`, because the
segfault comes before the wrapper moves the staged files into place. So a run
killed by this defect also loses the artifact that
`verify_staged_namelists` exists to leave behind, and there is nothing to read
back to find out what it integrated. The `plasim_diag` in the `_crashed`
directory is the only record, and it is not the one the gate reads.
