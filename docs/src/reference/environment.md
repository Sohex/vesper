# Environment

`.venv` is already active. Dependencies are pinned in `requirements.txt`; install
with `UV_CACHE_DIR=/tmp/world-uv-cache uv pip install --python .venv/bin/python -r requirements.txt`.
Building ExoPlaSim needs `gcc-fortran` and `openmpi` from the host (Arch).
Matplotlib is forced to `Agg` with its cache at `/tmp/world-matplotlib-cache`.

**ExoPlaSim is installed EDITABLE from its subtree, not from PyPI**, so a build
compiles in place and `requirements.txt` deliberately does not name it:

    uv pip install --python .venv/bin/python -e vendor/exoplasim

What the fork contains and where it came from is in [the vendored upstreams](vendored-upstreams.md).

**A NEW INTERPRETER MEANS A PYFFT REBUILD.** `pyburn` imports `exoplasim.pyfft`
to build the Gaussian grid, so without it no raw model output can be read and
every run fails at postprocessing. It and `pyfft991` are f2py extensions tagged
with a CPython ABI, and git tracks only their `.f90` sources -- so rebuilding
`.venv` on a new Python leaves the sources looking untouched and the extensions
unloadable. Rebuild them, which needs `meson`, `ninja` and `gfortran` on the
host because numpy refuses the distutils backend on Python >= 3.12:

    python exoplasim/scripts/build_pyfft.py            # build both and verify
    python exoplasim/scripts/build_pyfft.py --check    # do they load? exit 1 if not

`scripts/smoke_test.py` carries the same check, so a stale ABI is caught before
a run rather than after one. That matters because the failure does not look like
itself: `exoplasim/__init__.py` turns any postprocessing exception into
`_crash()`, which moves the run to `<run>_crashed/` and reports "ExoPlaSim has
crashed or begun producing garbage" for a model that integrated perfectly.

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

## Working in a git worktree

A worktree carries the tracked tree and nothing else, so everything this project
deliberately keeps out of history is absent from it: the Orogen export payloads,
the reference PDFs and bulk datasets, the climate run output, the Earth
validation caches, `node_modules`, `.venv`. A script that runs in the main
checkout dies on a missing file there.

    python scripts/link_worktree.py                    # from inside the worktree
    python scripts/link_worktree.py --worktree PATH    # from the main checkout
    python scripts/link_worktree.py --check            # report only, exit 1 if incomplete

It symlinks that set back to the main checkout. The set is DERIVED from the
ignore rules on every run rather than listed, so a new build or a new cache is
picked up without editing it, and a wholly-ignored directory is linked as a unit
while a directory holding tracked content is linked entry by entry. The second
half of that is what puts the links at `exoplasim/runs/<id>` rather than over
`runs/` itself, whose `INDEX.json` is tracked. Re-running it is how a worktree
picks up a build added since: a correct link is left alone, a stale one is
repaired, and one whose target has been archived away is removed.

Entry-by-entry linking has a consequence worth knowing before writing anything
in a worktree. A wholly-ignored directory is one symlink, so a write inside it
lands in the main checkout. A directory holding tracked content beside ignored
payload CANNOT be one -- `references/INDEX.md`, `exoplasim/runs/INDEX.json` and
the `source/` READMEs are all tracked -- so what the worktree gets is a real
directory of one link per EXISTING file, and a file created there afterwards is
real in the worktree alone. Being ignored, it is never committed; when the
worktree is removed it is gone, while a tracked row describing it survives and
outlives its own artifact. `references/` lost thirteen PDFs that way, with their
`INDEX.md` rows still standing. `--check` now fails on payload in that state and
names it, so the fix is to copy it to the main checkout before the worktree
goes. `notes/audits/worktree-stranded-payload.md` carries the shapes and the
evidence.

Two things are held back. Everything compiled from tracked source that a
worktree may have edited -- `vendor/exoplasim` and the LPJ-GUESS build -- is not
linked, because a link both hides the worktree's own edit behind the main
checkout's binary and lets a rebuild in the worktree overwrite that binary,
which is CLAUDE.md rule 4 with the safety off. Build them in the worktree, or
pass `--model-binaries` when the worktree does not touch the model. Regenerable
output is not linked either, for the narrower reason that the worktree's build
would land in the main checkout.

`exoplasim/bench` IS linked, as a wholly-ignored directory, and the transform
gates work inside it. Each one opens by deleting its own subdirectory there, so
two worktrees running the same gate delete each other's work. Every gate that
does this takes an environment variable naming somewhere else instead --
`BANDED_WORK`, `SHTNSMODEL_WORK`, `TNUMERICS_WORK`, `WFCHECK_WORK`,
`GAUSSWEIGHTS_WORK` -- and a worktree sets the one it needs.

`.venv` IS linked, and ExoPlaSim is installed editable from the main checkout's
`vendor/exoplasim`. So `import exoplasim` in a worktree reads the main
checkout's model source whichever tree the interpreter was invoked from. That is
a property of the editable install rather than of the link, and it is not
fixable from the worktree side: model work in a worktree needs its own editable
install and its own compile.

A symlink is a file and not the directory it stands in for, so an ignore rule
ending in `/` does not cover the link that replaces it. `.gitignore` carries a
second set of patterns for exactly these paths, and the script's last act is to
check `git status` in the worktree and name any link that block still fails to
cover. That is what keeps the block complete: it is tested rather than
remembered.

## The git hooks, and the beads export they keep honest

`core.hooksPath` points at `.beads/hooks`, so the hooks are TRACKED and apply in
every checkout and worktree without an install step. Each is a thin shim around
`bd hooks run <name>`, inside `BEGIN/END BEADS INTEGRATION` markers that
`bd hooks install` regenerates. Anything added inside those markers is lost on
the next beads upgrade; this project's own additions go after the END marker.

`pre-commit` carries one such addition, and the reason is a trap worth knowing
independently of the hook. **`export.auto` is debounced, not per-write.** A burst
of `bd` commands writes `.beads/issues.jsonl` once, near the start, and the rest
of the burst never reaches it -- so the export routinely lags the database by an
arbitrary amount, and `git status` showing it clean is not evidence that it
agrees with the tracker. Measured 2026-08-24: a `bd update` exported
immediately, and a following run of two `bd supersede` calls and five
`bd dep add` calls changed the file not at all.

The hook therefore regenerates the export from the database and re-stages it,
but only when it is ALREADY STAGED -- that is, only when you have decided to
record tracker changes in this commit. It does nothing otherwise. An unstaged
export that trails the database is the resting state of a derived artifact, not
a work list, and an earlier version of this hook that tested mtimes to warn
about it fired on every commit, because the beads block above it runs `bd` and
touches the database first.

Two consequences to expect rather than be surprised by. The export is a
whole-database snapshot, so regenerating it picks up any concurrent session's
bead writes as well: there is no per-session subset of it, and a shared snapshot
is the correct outcome where a stale one is not. And during a rebase, merge or
cherry-pick the hook stands down, because the export belongs to the commit being
replayed rather than to the database as it stands now. `WORLD_SKIP_BEADS_EXPORT=1`
bypasses it.

## The toolchain: what is declared, and what optimised libraries do not buy

The model links glibc's `libm` and `libmvec`, `libgfortran` and Open MPI, and
nothing else. PlaSim carries its own FFT and Legendre transform, so there is no
BLAS, LAPACK or FFTW seam anywhere in the build, and an optimised replacement
for one of those has nothing to attach to. Installing such a library is not a
reason to rebuild.

`config/planet.yaml` declares `model.compile_flags`, and
`exoplasim/scripts/build_model.py` reads it. **Declare flags there and nowhere
else**, which is now enforceable rather than advisory: there is no generated
compiler-options file to edit instead. There used to be three of them, they
disagreed, and only the MPI one was ever written from this declaration -- so
the threaded build, which is what this project measures, took whatever a
hand-edited file happened to hold. `notes/audits/model-build-driver.md`.

The two seams that do exist, `libm` and `memcpy`, were measured against AMD's
implementations of both, and the optimisation flags were measured against each
other. What that settled, with the numbers, is in
`notes/audits/aocl-and-model-build-flags.md`: AMD's libraries are refused,
`-march=znver4` is adopted and is worth 2 to 3%, `-ffpe-trap` costs under 1% and
is kept for what it catches, and `-flto` is refused.

**`-finit-real=zero` is not in the production line, and the reason is the
largest single figure in this file.** It initialises every local real on entry
to every routine, arrays included, and it cost 26.03% of T170 [+22.76, +27.29].
The model does not need it: the restart is bit identical without it, and
identical again under `-finit-real=snan` with the trap masked, so nothing
uninitialised reached a stored value at the length the bench runs. The poisoned
initialisation lives in the `poisoned` profile as `-finit-real=snan` alongside
`-Og`, and not in `checked`: at the optimisation levels this project ships
gfortran folds the signalling NaN to a quiet one at compile time, so a
deliberate uninitialised read in the model exits 0 rather than trapping.
`notes/audits/uninitialised-reads-and-implicit-save.md` finding 1 has that
control. Both profiles are run with `NSHTNS=0`: every build is threaded, and on
the SHTns path the model faults inside the library at the first timestep before
it can say anything about itself.
`exoplasim/notes/the-zeroing-is-an-init-flag.md` has the three gates and
the one thing they exposed, which is in SHTns rather than in the model. The model is compute-bound
inside its own Fortran, and 2 to 3% is what the whole distance from scalar code
to AVX-512 is worth here, which is also the bound on what any further codegen
work can return.

Benchmarking it again needs one rule beyond a quiet machine: **interleave the
arms**. Run as blocks, whichever arm goes first after an idle stretch gets the
boost clock and the comparison measures the CPU's thermal state instead of the
flag. That was worth 6% on a stock baseline against itself.

## One host, many agents: the lock, and what it does not do

This project fans work out across many agents on ONE machine, and several of
them run model integrations, builds and profiles. Anything that uses the CPU
for more than a moment runs under the wrapper:

    scripts/lock_and_run -m "what you are doing" python exoplasim/scripts/run_exoplasim.py ...

It waits for the lock, runs the command, and releases. The command's stdin,
stdout and stderr are inherited untouched, so a command under the wrapper
behaves exactly as it does without one and can still be piped on either side.

**ONE STEP, BECAUSE EVERY STEP A CALLER CAN SKIP HAS BEEN SKIPPED HERE.** The
lock used to be spelled out at each call site -- take a directory, write a
description into it, do the work, remove the directory -- and each of the three
ways to get that wrong cost real work. Leaving the directory behind held the
host against every other agent until a human noticed. Writing the description
WITHOUT taking the directory put three integrations on thirty-two cores at a
load of fifty-one and made every wall clock from that session unusable. Taking
it twice without releasing deadlocked a session against itself, silently,
because the waiter sleeps and prints nothing; the tell was a `who` file
describing a phase that had obviously finished. A wrapper has one step, so
there is nothing left to skip, and nesting is free rather than fatal: the
command runs with `WORLD_LOCK_HELD` set and a wrapper that sees it runs
straight through.

**THE LOCK THAT DECIDES IS AN `flock`, AND THAT IS WHY THERE IS NO STALE-LOCK
RECOVERY.** `/tmp/world.lock.flock` is a file that is never deleted; the lock
is the kernel's advisory lock on an open descriptor, and the kernel drops it
when the last descriptor closes. Releasing is therefore not an action anyone
has to remember, or survive long enough to perform. `kill -9` releases it. A
crash releases it. A stale flock does not exist, so there is no procedure for
clearing one and no judgement call about whether a holder is a corpse.

The descriptor is deliberately INHERITED by the command. Kill the wrapper and
the integration it started keeps running -- orphaned, but still using the
machine -- and it keeps holding the lock. That is the right answer and it is
the one a pid-file scheme gets wrong, because a pid file names the wrapper and
the wrapper is not the thing using the CPU.

**`/tmp/world.lock` REMAINS, AS THE VISIBLE CLAIM.** The wrapper creates that
directory while it holds the lock and removes it on release, with the same
`who` file as before, so `cat /tmp/world.lock/who` still answers "who has the
machine, and since when". It also keeps the two protocols interoperable in both
directions: anything still written as `until mkdir /tmp/world.lock` blocks
against the wrapper and is blocked by it. `mkdir` rather than `touch` is why
that half works at all -- it tests and takes in ONE atomic step, where the
obvious spelling, look for the file and then create it, lets two agents that
check in the same moment both find it free.

A wrapper that is killed leaves that directory behind, and its corpse is
identified rather than guessed at: `who` carries a marker line only the wrapper
writes, and a marked directory can only be a dead wrapper, because a live one
would still hold the flock the reader has just acquired. The wrapper clears it
and says so. An UNMARKED directory is a claim taken by hand, and the wrapper
waits for it instead of deciding it is dead.

The lock lives outside the repository because each fan-out agent works in its
own worktree, so an in-tree path is a different file for every agent and
coordinates nothing.

`scripts/lock_and_run --self-test` exercises all of this against a private lock
path: that the three streams and the exit status pass through, that a second
command cannot start until the first has finished, that a nested wrapper runs
instead of deadlocking, that a dead wrapper's claim is cleared and a
hand-rolled one is not. `scripts/verify_entry_points.py` runs it.

**HOLDING THE LOCK DOES NOT PARTITION THE HOST.** It keeps other AGENTS off the
machine. It says nothing about how you divide the machine between your own
processes, and reading it as a reservation is how two 16-thread integrations
came to be co-scheduled on 32 logical cores under a single claim. Two paired
arms at once are `p8` binaries pinned to their own cores, under ONE
`lock_and_run`; two `p16` at once is never right. An unpinned pair fights over
the same CCDs, which breaks the 32 MB per-die target below -- that target is
stated per thread TEAM, and two teams on one die exceed it silently. Worse, the
pair then measures the scheduler as much as the model, which is the one thing a
paired experiment exists to avoid.

**A timing is a measurement of a machine state as much as of a model.** Record
the load beside any timing worth keeping: one without the machine state it was
taken under cannot be compared against a later one, and a number taken while
someone else is integrating is not a slow number, it is a number of a different
experiment. Where the choice exists, price work in something the scheduler
cannot move -- retired instructions under `OMP_WAIT_POLICY=passive` survive
contention that wall clock does not, though once a run is threaded that count
needs its own correction, because a thread spinning at a barrier retires
instructions in proportion to how long it waits.

**What this replaced, and why.** A 172-line `scripts/machine.py` carried a claim
file, a process-table scan, a load threshold and a worker-count API. The parts
beyond "do not start heavy work while someone else is running" were not the job,
and the vocabulary actively misled: a "claim" reads as though it reserves the
host, which is the reading that put two 16-thread integrations on one box. What
it became was four lines an agent could follow without opening a script, and
what those four lines then cost was a caller executing three of them. The
wrapper is the same policy at one call: the mechanism may be a script again,
but the API it exposes is a command prefix and not a set of primitives to
sequence correctly.

## The per-die working set targets 32 MB

**A thread team's working set on one die targets 32 MB, and a change that takes
it above that is a regression even when this machine gets faster.** Declared
2026-08-21 as a standing constraint on compute work.

The number is not this desktop's. The 7950X3D has 96 MB of L3 on CCD0 and 32 MB
on CCD1, so half its threads have a cushion the other half does not, and a
working set between the two is fast on one die and slow on the other -- which is
exactly the load imbalance `exoplasim/notes/rank-imbalance-and-weight-traffic.md`
measured. 32 MB is the smaller of the two and is what a part with more cores and
no stacked cache is likely to have or less. Targeting it is targeting the
hardware this model might be moved to, not the hardware it is on.

The rule carries an assumption, stated so it can be argued with rather than
discovered: **above 32 MB is taken to be slower than below it, always.** That is
a design target and not a measured law -- a 33 MB set with good locality can beat
a 31 MB set that is streamed -- but the exception needs the argument, not the
rule.

**It applies per DIE and per TEAM, so count the copies.** A threadprivate array
is one copy a thread, and eight threads on a die means eight. That is how the
Legendre weight factorisation first missed: two matrices and six per-mode
vectors came to 35.76 MB a die because the vectors, identical on every thread,
were threadprivate. Shared, the same change is 30.82 MB. The 4.94 MB of
redundant copies was the whole margin, and no arithmetic about the matrices
would have found it.

**A change that moves toward the ceiling is worth a small loss on this
machine.** The weight factorisation costs 1.56% here and is kept, because the
CCD that never had a cache problem is hiding the case the change exists for.

**The target is the TOTAL, and the total is not there yet.** Counting only the
component you just changed is how this rule gets quoted as satisfied while
nothing fits. At T170 the transform loops cycle about 74 MB a die: 30.8 MB of
weights, 31.5 MB in the six Fourier fields `mktend` holds, and 12.0 MB of
spectral state. The factorisation took the weights from 120 MB, four times the
die, to about one -- a large win on what had been the dominant term, and not a
fit. The next terms are named by that arithmetic rather than guessed at.

**Do not set a headroom margin below 32 until the total is near it.** Interrupts,
the OS, page tables and set-associativity are all real and together are worth
something like a tenth of the capacity, which is a correction to make when the
budget is close. Applied to a working set 2.3x over, a tighter number is false
precision and invites treating the trimmed component as the whole.

### Counting it: `exoplasim/scripts/cache_budget.py`

The budget is a script rather than a table here, because a number recomputed by
hand every time a resolution or a thread count moves is a number that gets
quoted stale, and because a component of it gets quoted as the whole. Both
happened while this rule was being written.

    python exoplasim/scripts/cache_budget.py --res T170 --threads 16 --per-die 8

**Three storage classes, and the third is the one that catches people.**
`private` is one copy a thread, so eight on a die. `shared` is one copy that
every thread reads all of. `sliced` is one copy of which a thread only ever
touches its own `NSPP` rows -- the reduction partials are this, since thread t
writes slot t and reads only its own rows of every slot. Counting a sliced array
as fully resident overstates it by the thread ratio, which at sixteen threads on
an eight-thread die is a factor of two.

At T170 on sixteen the total is about 260 MB a die against the 32 MB target,
eight times over, and the largest terms are the reduction partials rather than
anything the transforms hold. The physics modules' grid arrays and their
compiler temporaries are not counted, so it is a floor.

### A second target, not yet a priority: the hot loops in L2

Zen 4 has 1 MB of L2 a core, private to it. **The innermost transform loops
should aim to fit their working chunk in that**, which is a different and
stricter goal than the die-level one and needs a different technique: blocking
the latitude and level loops so the chunk in flight is a megabyte rather than
streaming a whole `NCSP x NLPP` matrix past. One latitude's column of `pmat` is
118 KB at T170 and fits easily; the loop as written walks 1.88 MB of it.

This is recorded so it is not rediscovered, and it is NOT the current priority.
The die-level total is eight times its target, and blocking for L2 under that is
optimising the inner loop of a problem whose outer numbers are wrong.

## Benchmark in FP32 while optimising, confirm in FP64 before believing

A four-byte build runs T170 39.7% faster than an eight-byte one, so the
iteration loop during optimisation work is roughly a third shorter for nothing.
Use it. Precision is declared in `config/planet.yaml`, so an FP32 arm means
changing `model.precision_bytes` for the arm rather than passing a flag.
Declared 2026-08-21.

**The hazard is that FP32 is a different cache regime, not just a faster one.**
It halves the per-die working set, 260 MB to 130 at T170. So a change whose
benefit is about FOOTPRINT can measure well in single precision and vanish in
double, and the reverse: the Legendre weight factorisation is exactly that shape
of change, and it is kept for a cache argument that FP32 would have muddied.

The rule that follows: **iterate in FP32, confirm in FP64, and quote only the
FP64 number.** Anything cache- or bandwidth-bound is confirmed rather than
extrapolated, and a result that appears only in one precision is a finding about
the precision rather than about the change.

Production is unaffected. `config/planet.yaml` declares eight-byte precision and
FP32 is a benching and spin-up tool, never a run anything is read from --
CLIM-59 and main's CLIM-52 carry that.
## Reading the artifacts from the command line

The host carries the netCDF and NCO command-line tools, and they are the first
reach for inspecting a run or a climatology: `ncdump -h` for structure and
global attributes, `ncks` to subset, `ncdiff` to difference two files, `ncwa` to
collapse dimensions, `ncatted` for metadata. Model output and analysis products
are all `NETCDF4`, which is HDF5 underneath, so `h5diff` compares two files at
the value level with a tolerance, `-d` absolute and `-p` relative. That is the
tool for asking whether one run reproduces another.

Beside them: `yq` for `config/pipeline.yaml` and `config/planet.yaml`, the GDAL
command-line tools `gdalinfo` and `ogrinfo` for the reference shapefiles and the
Copernicus DEM COGs,
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
