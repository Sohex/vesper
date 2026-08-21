# The production flags cannot see an uninitialised read, and nothing sees an implicit SAVE

Measured 2026-08-21 on a Ryzen 9 7950X3D, gfortran 16.2.1 20260810.

Worldbuilding frame: this file is about how the Vesper simulation's climate
model is COMPILED and which programming defects the build can and cannot
detect. Nothing here is about the simulated world.

## The question

The Stage 1 thread port found an implicit SAVE. Two things follow that the
build flags are the natural place to ask about: what would have caught that
defect earlier, and what else of its class is sitting in `plasim/src`.

## The production flag line

    MOST_F90_OPTS=-O3 -cpp -fopenmp -ffixed-line-length-132
                  -ffpe-trap=invalid,zero,overflow -ffpe-summary=none
                  -finit-real=zero -march=znver4 -fno-omit-frame-pointer
    MOST_PREC=-fdefault-real-8

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

## Finding 2: implicit SAVE is invisible to the compiler and is a race under threads

A local declared with an initialiser inside a procedure acquires the SAVE
attribute and persists across calls. gfortran reports nothing, at any warning
level tried: `-Wall -Wextra`, plus `-Wsurprising`, plus
`-Wcharacter-truncation -Wconversion-extra`. No diagnostic mentions SAVE.

Under `-fopenmp`, which the production line carries, a SAVEd local is SHARED by
every thread. Eight threads each write their own id into a SAVEd scratch array,
spin, then read it back; **seven of the eight read a value another thread
wrote.**

This is why the defect belongs to the thread port specifically. Under MPI each
rank is a separate process, so its SAVEd locals are private by construction and
the same source is correct. The identical source under OpenMP is a data race.
Ranks and threads therefore do not merely differ in performance here; they
differ in which defects are reachable, and a source that has only ever run
under ranks has never exercised this class.

## What is on disk

103 procedure-body declarations in `vendor/exoplasim/exoplasim/plasim/src`
carry an initialiser and so are implicitly SAVEd:

| file | sites |
| --- | ---: |
| plasim_dummy.f90 | 18 |
| icemod.f90 | 14 |
| icemod_template.f90 | 14 |
| hurricanemod.f90 | 13 |
| mpimod.f90 | 8 |
| mpimod_multi.f90 | 7 |
| carbonmod.f90 | 5 |
| mpimod_stub.f90 | 5 |
| cpl.f90 | 4 |
| glaciermod.f90 | 4 |
| remainder | 11 |

Counted by walking each file and tracking procedure nesting, matching type
declarations that carry `::` and an initialiser, and excluding `parameter`
declarations and module-level variables, whose static storage is intended. The
method's limits are its own: continuation lines are not joined, and a genuinely
intended persistent counter is indistinguishable from a scratch array by shape
alone.

So 103 is a POPULATION, not a defect count. A site that is written before it is
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
reading the 103 sites.
