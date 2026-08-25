# The production flags cannot see an uninitialised read, and nothing sees an implicit SAVE

Measured 2026-08-21 on a Ryzen 9 7950X3D, gfortran 16.2.1 20260810.

Worldbuilding frame: this file is about how the Vesper simulation's climate
model is COMPILED and which programming defects the build can and cannot
detect. Nothing here is about the simulated world.

## The question

The Stage 1 thread port found an implicit SAVE. Two things follow that the
build flags are the natural place to ask about: what would have caught that
defect earlier, and what else of its class is sitting in `plasim/src`.

## The production flag line, as measured

    MOST_F90_OPTS=-O3 -cpp -fopenmp -ffixed-line-length-132
                  -ffpe-trap=invalid,zero,overflow -ffpe-summary=none
                  -finit-real=zero -march=znver4 -fno-omit-frame-pointer
    MOST_PREC=-fdefault-real-8

**Both findings below were acted on and the flag line moved.** The declaration
is `config/planet.yaml`'s `model.compile_flags.f90_opts` and there is no
`most_compiler*` file any more. `-finit-real=zero` has LEFT the production
profile, so the poisoned initialisation is where a check can be run rather than
where it is paid for on every production orbit.

**Where it went is the `poisoned` profile, and it had to carry `-Og` with it.**
A profile that is `f90_opts` plus `-finit-real=snan` compiles at the
optimisation level this finding measures the fold at, and the control below
confirms the fold in the model rather than in a test case. `checked` is
therefore `-fcheck=all`, which is what it certainly delivers, and `poisoned` is
`-fcheck=all -finit-real=snan -Og`, which is the arm the control faults in.
Both are run with `NSHTNS=0`: every build is threaded since world-38b, and on
the SHTns path the model faults inside the library at the first timestep.

**The check that works at production's optimisation level is the PROPAGATING
one** -- `-finit-real=snan` with `FE_INVALID` masked, restart compared bit for
bit against production on the same bed -- because it needs no trap to arm and a
quiet NaN in a stored record is as visible as a signalling one. It is
`exoplasim/scripts/verify_uninitialised_reads.sh`, a `checks` row in
`config/pipeline.yaml`, and its positive control is a deliberate uninitialised
read reaching the temperature tendency. The flag traded for that gate,
`-finit-real=zero`, was worth 26.03% of T170.
`exoplasim/notes/the-zeroing-is-an-init-flag.md` carries the measurements on
both sides of the trade.
`notes/audits/model-build-flags.md` is the re-run of the whole line on the
threaded build and is what a reader wanting today's flags should read.

`notes/audits/aocl-and-model-build-flags.md` keeps `-ffpe-trap` on the argument
that a NaN then raises SIGFPE at the line that made it, rather than propagating
through the Legendre transform into every spectral coefficient. That argument
is sound for a NaN the arithmetic produces. It does not reach an uninitialised
read, and the reason is the flag sitting four positions to its right.

## Finding 1: `-finit-real=zero` makes the uninitialised case untrappable

`-finit-real=zero` fills uninitialised reals with 0.0. Zero is a valid operand,
so nothing raises IEEE invalid and `-ffpe-trap` never fires. The read succeeds
and its result is indistinguishable from a deliberate zero.

A four-element array is declared, never assigned, and summed. Same source, same
flag line, only `-finit-real` varied:

| build | result | exit |
| --- | --- | --- |
| production line, `-finit-real=zero` | 0.0 | 0 |
| production line, `-finit-real=snan` | NaN | 0 |

The signalling NaN does not trap either, and that is the second half of the
finding. Holding `-ffpe-trap=invalid,zero,overflow -finit-real=snan` fixed and
varying only the optimisation level:

| level | outcome |
| --- | --- |
| `-O0` | SIGFPE, trapped at the arithmetic |
| `-Og` | SIGFPE, trapped at the arithmetic |
| `-O1` | quiet NaN, exit 0 |
| `-O2` | quiet NaN, exit 0 |
| `-O3` | quiet NaN, exit 0 |

At `-O1` and above gfortran folds the signalling NaN to a quiet one before any
instruction executes, so the trap is never armed. `-fno-tree-vectorize` at `-O3`
does not change this, which rules out vectorisation as the mechanism.

The consequence for any debug configuration is that it is NOT the production
line with one flag swapped. Poisoned initialisation only survives to the
arithmetic at `-Og` and below, so a debug build is a different optimisation
level, and its numerics differ from production for that reason alone rather
than by choice. `-Og` rather than `-O0` is the useful ceiling: it traps, and it
is far cheaper than `-O0`.

### The same fold, measured inside the model

*Measured 2026-08-24 at 743c67a9, T21 on sixteen threads, one timestep on a cold
T21 bed with `NSHTNS=0`, gfortran 16.2.1 20260810.*

The table above is a four-element local declared, never assigned and summed in
one scope, which the compiler can see whole. That leaves open whether the fold
also reaches a local written in one branch and read in another, which is the
shape a real defect takes. It does.

A positive control was placed at the head of `gridpointd` in `plasim.f90`: two
never-assigned four-element locals, one summed straight and one whose only
assignment sits under a condition the compiler cannot prove false, both summed
into a `volatile` scalar so neither store can be deleted or sunk. Each
optimisation level was built twice from the `checked` flag line, once with the
control and once without, and both binaries were run on the same bed.

| level | with the control | without it |
| --- | --- | --- |
| `-Og` | SIGFPE in `gridpointd` at the first timestep | exit 0, restart written |
| `-O2` | exit 0, restart written | exit 0, restart written |
| `-O3` | exit 0, restart written | exit 0, restart written |

The control is the only difference between the columns, so the `-Og` fault is
the control's and the clean exits at `-O2` and `-O3` are the trap not arming.

The mechanism is visible in the instruction stream. At `-Og`, `gridpointd`
opens with two loops that store the signalling NaN word into the stack slots
and then a `vaddsd` chain that reads them, so the pattern reaches the FPU. At
`-O2` and `-O3` there is no `vaddsd` chain at all: the sums are folded at
compile time and only the `volatile` stores survive. Both binaries carry two
NaN constants in `.rodata` -- `0x7ff4000000000000`, the signalling pattern the
initialisation stores, and `0x7ff8000000000000`, the quiet one the fold
produced -- which is the fold on disk. The branch-merged site folds too,
because the compiler splits the path and folds each side separately.

So the fold is a property of the optimisation level and not of how visible the
variable's life is, and `-Og` is where a poisoned initialisation is a trap.
`config/planet.yaml` declares that arm as the `poisoned` profile;
`checked` is `-fcheck=all` and does not see this class.

## Finding 2: implicit SAVE is invisible to the compiler and is a race under threads

A local declared with an initialiser inside a procedure acquires the SAVE
attribute and persists across calls. gfortran reports nothing, at any warning
level tried: `-Wall -Wextra`, plus `-Wsurprising`, plus
`-Wcharacter-truncation -Wconversion-extra`. No diagnostic mentions SAVE.

Under `-fopenmp`, a SAVEd local is SHARED by every thread. Eight threads each
write their own id into a SAVEd scratch array, spin, then read it back; **seven
of the eight read a value another thread wrote.**

This is why the defect belongs to the thread port specifically. Under MPI each
rank was a separate process, so its SAVEd locals were private by construction
and the same source was correct. The identical source under OpenMP is a data
race. Ranks and threads did not merely differ in performance here; they differed
in which defects are reachable, and a source that had only ever run under ranks
had never exercised this class.

**That asymmetry is now the whole picture rather than half of it.** world-38b
removed the MPI and serial build paths, so every build this project makes is
`-fopenmp` and every one of the sites below is in the reachable half. There is no
longer a configuration in which an implicit SAVE is harmless by construction.

## What is on disk

Counted 2026-08-21: 103 procedure-body declarations in
`vendor/exoplasim/exoplasim/plasim/src` carried an initialiser and so were
implicitly SAVEd, across `plasim_dummy.f90` 18, `icemod.f90` 14,
`icemod_template.f90` 14, `hurricanemod.f90` 13, `mpimod.f90` 8,
`mpimod_multi.f90` 7, `carbonmod.f90` 5, `mpimod_stub.f90` 5, `cpl.f90` 4,
`glaciermod.f90` 4, and 11 in the remainder.

**Re-counted 2026-08-24 by the same method: 46 sites across nine files.**

| file | sites |
| --- | ---: |
| hurricanemod.f90 | 13 |
| icemod.f90 | 11 |
| carbonmod.f90 | 5 |
| glaciermod.f90 | 4 |
| landmod.f90 | 4 |
| oceanmod.f90 | 3 |
| seamod.f90 | 3 |
| mpimod_omp.f90 | 2 |
| surfmod.f90 | 1 |

Most of the fall is deletion rather than repair: `plasim_dummy.f90`,
`icemod_template.f90` and `cpl.f90` went under world-cmz, and `mpimod.f90`,
`mpimod_multi.f90` and `mpimod_stub.f90` under world-38b, which is 56 sites in
files no configuration compiled. `icemod.f90` lost three to the salinity, lead
and cold-start work and to world-ro6's removal of unreferenced procedures. The
files that remain are the files that are compiled, so the population is now
entirely live.

Counted by walking each file and tracking procedure nesting, matching type
declarations that carry `::` and an initialiser, and excluding `parameter`
declarations and module-level variables, whose static storage is intended. The
method's limits are its own: continuation lines are not joined, and a genuinely
intended persistent counter is indistinguishable from a scratch array by shape
alone.

So the count is a POPULATION, not a defect count. A site that is written before it is
read on every call is harmless whatever its storage class, and some of these
are deliberate. What the number establishes is that the class is enumerable,
that it is larger than the one site the thread port tripped over, and that
nothing in the build or the test path distinguishes the harmless members from
the rest.

## What this note does not claim

Poisoned initialisation does not catch implicit SAVE. An implicitly SAVEd
variable IS initialised, exactly once, so no `-finit-real` value ever reaches
it. The two findings share a chapter because they share a remedy discussion
about build flags, not because either flag addresses the other's defect. The
remedy for finding 1 is a build configuration; the remedy for finding 2 is
reading the sites the population count leaves standing.
