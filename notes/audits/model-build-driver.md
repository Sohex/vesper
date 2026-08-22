# The model build driver fails silently, and in the direction that costs most

*Worldbuilding frame: this is about the BUILD SYSTEM of the Vesper project's
climate model, a vendored fork of ExoPlaSim. Nothing here concerns the simulated
planet. Measured 2026-08-22 against `vendor/exoplasim/exoplasim/compile.sh` and
`configure.sh` at 335de8c7.*

Every defect below shares one shape: an input the driver does not recognise
produces a DIFFERENT BUILD rather than an error, and the difference is invisible
downstream. `archive/tasks.md` CLIM-22 is this shape having already cost a
session -- every 16-rank binary silently built single precision -- and the
mechanism that allowed it is untouched.

## An unrecognised `-r` builds T21 and says nothing

`compile.sh` parses the resolution with a `case` whose arms are the seven names
`T21`..`T170` and the seven latitude counts `32`..`256`. The final arm is
`\?)`, which in `bash` `case` matches the literal character `?`. It is not a
wildcard, and there is no `*)` arm, so the message it guards --
"INVALID RESOLUTION PASSED! Reverting to T21." -- is unreachable for every value
that could trigger it. An unmatched value falls through leaving the defaults set
at the top of the file: `resolution="t21"`, `latitudes=32`, `longitudes=64`.

MEASURED: `./compile.sh -j -p 8 -n 16 -r 170 -v 10` exits 0, announces
`most_plasim_t21_l10_p16_omp.x`, and writes `parameter(NLAT_ATM = 32)` into
`plasim/bld/resmod.f90`. `170` is not one of the accepted spellings; `T170` and
`256` are.

`-p` has the same structure and a worse default. Its arms are `4`, `8`,
`single`, `double`, then the same unreachable `\?)`, and the default at the top
of the file is `prec=4`. An unrecognised precision therefore does not fail, it
builds SINGLE precision -- which is CLIM-22 exactly. `-fdefault-real-8` changes
the width of `real` in every declaration and every interface, so a real*4 object
linked against real*8 ones does not fail to link. It computes.

## A defaulted build leaves the caller's name pointing at a stale binary

The executable is named from the resolution the parse arrived at, and
`compile.sh` removes `plasim/run/$executable` only for the name it is ABOUT to
write. When the resolution silently falls back, the name falls back with it, and
any earlier binary under the name the caller expects survives the build
untouched.

MEASURED, and it is how the two defects surfaced: three arms intended to differ
only in `-finit-real`, built with `-r 170`, came out byte for byte identical.
All three were copies of a `most_plasim_t170_l10_p16_omp.x` dated 06:09 that no
build in the session had written. The real arms, built with `-r T170`, differ:
1193 memset call sites against 623.

A caller cannot defend against this with the exit status, and checking that the
file exists is not enough either. The check that works is to assert that
`plasim/bld/resmod.f90` holds the latitude count that was asked for and that the
executable is newer than the build that claimed to write it.

## The declared flag line reaches the MPI build only

`docs/src/reference/environment.md` says the flags are declared in
`config/planet.yaml` under `model.compile_flags.f90_opts` and nowhere else, and
`rebuild_binaries.py` writes that line in full before every build -- into
`most_compiler_mpi`.

`compile.sh -j` reads `most_compiler_omp`, and nothing derives that file from
`config/planet.yaml`. It is a hand-edited file in the vendored tree, and the two
agree today only because they were kept in step by hand.

This is not a corner. **All twelve binaries in `exoplasim/binary_manifest.json`
are MPI**, and every threaded binary -- what `verify_shtns_model.sh`,
`verify_shared_determinism.sh`, `thread_count_sweep.sh` and `bench_ab.py` build
and measure, which is this entire optimisation workstream -- is ad hoc and
unregistered. `binary_manifest.json` records `most_compiler_mpi`'s contents as
the toolchain, so a threaded binary that did reach the registry would carry a
provenance record of flags it was not built with.

## What is NOT wrong: the object tree

`compile.sh` does `rm -rf *` inside `plasim/bld` on every run, so every build is
a full build and no flag change can leave a stale object behind. A T170 threaded
build costs about 21 s on this host. The "stale objects after a flag change"
worry is answered, and answered in the safe direction.

## configure.sh

380 lines, of which the majority is commented-out X11 detection for a GUI this
project builds as `guimod_stub`. What it does that matters is regenerate
`most_compiler*`, so it can reintroduce defaults into files the project treats
as generated -- `-fcheck=all` is the one that costs, 5% at T42 and 17% at T127.
`rebuild_binaries.py` already works around it by writing the whole
`MOST_F90_OPTS` line rather than appending through `compile.sh -O`, and says so
where it does it.

Nothing in this project calls `configure.sh`. The upstream install path in
`vendor/exoplasim/exoplasim/__init__.py` does, and this project does not use it:
the package is installed editable and built through `rebuild_binaries.py`.
