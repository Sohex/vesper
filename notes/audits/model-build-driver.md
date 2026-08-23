# The model build driver fails silently, and in the direction that costs most

*Worldbuilding frame: this is about the BUILD SYSTEM of the Vesper project's
climate model, a vendored fork of ExoPlaSim. Nothing here concerns the simulated
planet. Measured 2026-08-22 against `vendor/exoplasim/exoplasim/compile.sh` and
`configure.sh` at 335de8c7.*

Every defect below shares one shape: an input the driver does not recognise
produces a DIFFERENT BUILD rather than an error, and the difference is invisible
downstream. Closed CLIM-22 is this shape having already cost a
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
agree today only because they were kept in step by hand -- CLIM-81 had to edit
both.

There is a third file and a third shape. `most_compiler`, the serial one, keeps
`MOST_F90_OPTS` on **line 2** where the other two have it on line 3, and it
carries `-fcheck=all` and no `-march=znver4`, so it is the checked profile's
flags under the production profile's name. Line 3 matters: `compile.sh`'s `-O`
hook is `sed -i '3s/$/ '$optimization'/'`, and `rebuild_binaries.py` asserts
`F90_OPTS_LINE = 3` for the file it writes. On the serial file line 3 is
`MPIMOD=mpimod_stub`, so `-O` appends the flag to the MPIMOD line and the build
takes a module name with a compiler flag stuck to it.

This is not a corner. **All twelve binaries in `exoplasim/binary_manifest.json`
are MPI**, and every threaded binary -- what `verify_shtns_model.sh`,
`verify_shared_determinism.sh`, `thread_count_sweep.sh` and `bench_ab.py` build
and measure, which is this entire optimisation workstream -- is ad hoc and
unregistered. `binary_manifest.json` records `most_compiler_mpi`'s contents as
the toolchain, so a threaded binary that did reach the registry would carry a
provenance record of flags it was not built with.

## The object tree is safe, and it is what serialises the matrix

`compile.sh` does `rm -rf *` inside `plasim/bld` on every run, so every build is
a full build and no flag change can leave a stale object behind. A T170 threaded
build costs about 21 s on this host. The "stale objects after a flag change"
worry is answered, and answered in the safe direction.

The cost is that there is ONE build directory for every configuration, so two
builds cannot run at once without one deleting the other's objects. That is the
reason `rebuild_binaries.py` is serial across the twelve executables of its
`MATRIX`, and the reason the marker files exist at all: `MPI`, `OMP`, `PREC4`,
`PREC8` and `FPTR` are there to notice that the shared directory now holds
objects from a different configuration. A per-configuration build directory
would retire all five and make the matrix parallel as well as each build.
`aocl-and-model-build-flags.md`, "Is `make -j` safe now", has what parallelism
inside one build is worth -- 3.5x, saturating by `-j8` -- and CONS-13 is the one
missing dependency edge that holds it.

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

## Replaced, 2026-08-22

`compile.sh`, `configure.sh`, `make_plasim`, the three `most_compiler*` files,
the two `most_precision_options*` files and the four `most_{ice,snow}_build*`
scripts are deleted. `vendor/exoplasim/exoplasim/plasim/CMakeLists.txt` and
`exoplasim/scripts/build_model.py` replace them.

**What the equivalence rests on.** The bar was fixed before the first build:
byte-identical, or else an identical restart plus an identified reason for the
difference. Byte-identity FAILED and the reason is known. CMake compiles out of
tree, so gfortran embeds absolute paths in the runtime-diagnostic strings it
carries -- `In file '<134 characters>', around line 590` -- where the old
in-tree build embedded a bare filename. The same 1640 strings appear in both,
each about 100 characters longer, `.rodata` grows 35,776 bytes, and every
address downstream of it shifts. `-ffile-prefix-map` does not reach those
strings, because they arrive through `#line` directives in the file CMake
preprocesses for its module scanner.

What is identical, and this is the stronger statement:

| | T170 omp | T21 p8 mpi |
| --- | --- | --- |
| `.text` size | 1,917,198 both | 1,900,942 both |
| defined symbols, name and size | 2187, no differences | 2194, no differences |

and the restart, T170 over 60 steps on the SHTns path with `NLOWIO=1`, is
`017a35afbdf1f5f6c621c734` from both builds. Two successive registry rebuilds
also produced identical shas for all twelve executables, so the new build is
reproducible run to run.

**And the gate agrees, which is the check that did not need interpreting.**
`verify_shtns_model.sh` at T21 on four threads passes whole through the new
build: both arms at rounding scale, the control rejected, and the four-run
bit-identity arm giving `fcf46ebbb1302d3a` -- THE SAME HASH the same gate
produced before the build system was replaced. So the binary the new system
produces computes bit for bit what the old one's did, which is a statement about
the model rather than about the executable's layout.

**What each defect became.** The silent `-r` and `-p` defaults are gone: every
argument is checked against a list and an unrecognised value exits non-zero
having written nothing. The stale-binary trap is gone: the output name is
removed first and written last, only on success. The three divergent compiler
files are one declaration. The hand-declared dependency edges are gone --
Ninja scans the `use` statements, so CONS-13's missing edge is not a thing that
can be missing. And `plasim/bld` is one directory per configuration, so the
registry builds concurrently.

**What parallelism turned out to be worth.** A clean rebuild of all twelve
registered executables: **19 s**, against about 300 s for the serial
`compile.sh` loop it replaced. That is 3.5x inside each build, which
`aocl-and-model-build-flags.md` measured, multiplied by four builds at once.

**The upstream Python API is fenced rather than followed.** `exoplasim/__init__.py`
configured and compiled on demand -- `sysconfigure()` ran `configure.sh`, and
`Model.__init__` ran `compile.sh` whenever the executable it wanted was absent.
Both are removed. A `Model` that finds no executable now raises, naming
`build_model.py`, because building one silently is how the registry stopped
describing the binaries that existed. Note the editable install resolves to the
main checkout rather than to a worktree, so this takes effect on merge.
