# Environment

`.venv` is already active. Dependencies are pinned in `requirements.txt`; install
with `UV_CACHE_DIR=/tmp/world-uv-cache uv pip install --python .venv/bin/python -r requirements.txt`.
Building ExoPlaSim needs `gcc-fortran` and `openmpi` from the host (Arch).
Matplotlib is forced to `Agg` with its cache at `/tmp/world-matplotlib-cache`.

**ExoPlaSim is installed EDITABLE from its subtree, not from PyPI**, so a build
compiles in place and `requirements.txt` deliberately does not name it:

    uv pip install --python .venv/bin/python -e vendor/exoplasim

What the fork contains and where it came from is in [the vendored upstreams](vendored-upstreams.md).

LPJ-GUESS is likewise compiled from its subtree. Its generated `vesper.h` must
exist first because the configured orbital year sizes arrays at compile time:

    python biosphere/scripts/build_vesper_header.py
    cmake -S vendor/lpj-guess -B vendor/lpj-guess/build \
      -DCMAKE_BUILD_TYPE=Release -DUNIT_TESTS=OFF
    cmake --build vendor/lpj-guess/build --parallel 16

The build directory and generated header are ignored. Do not substitute an
external checkout: run manifests point at the vendored binary and hash it.

Two things about the model are worth knowing before you touch it. Several
changes are NO-OPS until a namelist key turns them on -- `h2osww` defaults to
1.0, `ndustrad` to 0, `nsolcycle` to 0 -- so a rebuilt binary reproduces the runs
that exist, and enabling one is a configuration decision that moves the mean.
And the low-I/O change ALTERS THE RESTART LAYOUT, so a run started before it
cannot be resumed by a binary built after it.

The stellar cycle is one of those switches rather than a separate build. There
is no cycle executable and no cycle tree: `nsolcycle` defaults to 0 and both
amplitudes to 0.0, which reduces the guard to `gsolinst = gsol0`, so the
ordinary binaries are bit-exact identical to an unpatched model until a cycle is
configured on.

## The toolchain: what is declared, and what optimised libraries do not buy

The model links glibc's `libm` and `libmvec`, `libgfortran` and Open MPI, and
nothing else. PlaSim carries its own FFT and Legendre transform, so there is no
BLAS, LAPACK or FFTW seam anywhere in the build, and an optimised replacement
for one of those has nothing to attach to. Installing such a library is not a
reason to rebuild.

`config/planet.yaml` declares `model.optimization_flag`, which
`rebuild_binaries.py` passes to `compile.sh`'s `-O` hook and which lands on
`MOST_F90_OPTS`. **Declare it there and nowhere else.** Editing the generated
`most_compiler_mpi` looks equivalent and is not: `configure.sh` rewrites that
file, and `compile.sh` appends `-O` to a COPY of it under `plasim/bld`, so an
edit there is both temporary and invisible to anything reading the original.

The two seams that do exist, `libm` and `memcpy`, were measured against AMD's
implementations of both, and the optimisation flags were measured against each
other. What that settled, with the numbers, is in
`notes/audits/aocl-and-model-build-flags.md`: AMD's libraries are refused,
`-march=znver4` is adopted and is worth 2 to 3%, `-ffpe-trap` costs under 1% and
is kept for what it catches, and `-flto` is refused. The model is compute-bound
inside its own Fortran, and 2 to 3% is what the whole distance from scalar code
to AVX-512 is worth here, which is also the bound on what any further codegen
work can return.

Benchmarking it again needs one rule beyond a quiet machine: **interleave the
arms**. Run as blocks, whichever arm goes first after an idle stretch gets the
boost clock and the comparison measures the CPU's thermal state instead of the
flag. That was worth 6% on a stock baseline against itself.
