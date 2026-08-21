# Environment

`.venv` is already active. Dependencies are pinned in `requirements.txt`; install
with `UV_CACHE_DIR=/tmp/world-uv-cache uv pip install --python .venv/bin/python -r requirements.txt`.
Building ExoPlaSim needs `gcc-fortran` and `openmpi` from the host (Arch).
Matplotlib is forced to `Agg` with its cache at `/tmp/world-matplotlib-cache`.

**ExoPlaSim is installed EDITABLE from its subtree, not from PyPI**, so a build
compiles in place and `requirements.txt` deliberately does not name it:

    uv pip install --python .venv/bin/python -e vendor/exoplasim

What the fork contains and where it came from is in [the vendored upstreams](vendored-upstreams.md).

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

## Reading the artifacts from the command line

The host carries the netCDF and NCO command-line tools, and they are the first
reach for inspecting a run or a climatology: `ncdump -h` for structure and
global attributes, `ncks` to subset, `ncdiff` to difference two files, `ncwa` to
collapse dimensions, `ncatted` for metadata. Model output and analysis products
are all `NETCDF4`, which is HDF5 underneath, so `h5diff` compares two files at
the value level with a tolerance, `-d` absolute and `-p` relative. That is the
tool for asking whether one run reproduces another.

Beside them: `yq` for `config/pipeline.yaml` and `config/planet.yaml`, the GDAL
command-line tools for the reference shapefiles and the Copernicus DEM COGs,
`dot` for rendering a graph, `valgrind` for the Fortran, and `ncdu` for the run
tree. In the venv, `dask`, `flox` and `bottleneck` back xarray over anything
run-sized; see [large data](large-data.md), which they assist and do not
replace.

### Three traps, because these tools assume Earth

**NCO does not know this planet's grid, and `ncwa -a lat,lon` is an UNWEIGHTED
mean.** The files carry no Gaussian weight variable, so there is nothing for
`-w` to find and nothing warns. Over a Gaussian latitude grid that produces a
global mean which is wrong and entirely plausible, which is the worst shape a
number can have. Area weights come from the project's own grid convention in
`lib/gridding.py`.

**The time axis is not a calendar.** It is `units = timesteps` with no
`calendar` attribute and raw counts for values. Any operator that assumes an
Earth calendar either refuses or silently imposes 365 days on a planet whose
orbital period is not that. Time-bin weights come from `lib/climatology.py`, and
the orbital period from `lib/orbit.py`.

**`ncra` needs a record dimension and these files have none.** `time` is a fixed
dimension, so a time mean wants `ncks --mk_rec_dmn time` first.

### What was refused, so the search is not run twice

**CDO**, whose distinguishing operators are calendar climatologies and Earth
remapping. On a `units = timesteps` axis the first family is unusable and the
second is a trap, and NCO covers the reductions that remain. **nccmp**, because
every file here is `NETCDF4` and `h5diff` already compares values with a
tolerance. **cartopy**, because `maps/projections.py` is the projection layer
and cartopy's value is Earth coastlines and features. **ccache**, because its
handling of Fortran `.mod` outputs is not reliable enough to trust against
rule 4; the build-time lever that was measured instead is `make -j`, in
`notes/audits/aocl-and-model-build-flags.md`.
