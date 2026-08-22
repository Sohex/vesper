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
uninitialised reaches a stored value. It now lives in the `checked` profile,
where a read of an uninitialised local traps instead of quietly returning a
zero. `exoplasim/notes/the-zeroing-is-an-init-flag.md` has the three gates and
the one thing they exposed, which is in SHTns rather than in the model. The model is compute-bound
inside its own Fortran, and 2 to 3% is what the whole distance from scalar code
to AVX-512 is worth here, which is also the bound on what any further codegen
work can return.

Benchmarking it again needs one rule beyond a quiet machine: **interleave the
arms**. Run as blocks, whichever arm goes first after an idle stretch gets the
boost clock and the comparison measures the CPU's thermal state instead of the
flag. That was worth 6% on a stock baseline against itself.

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
