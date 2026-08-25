# Pairing mirror latitudes on a process

*This is a COMPUTE change to the Vesper worldbuilding project's climate model.
Nothing in it is about the simulated planet. Measured 2026-08-21.*

**THE DECOMPOSITION DESCRIBED HERE IS DELETED.** `2508bedb` removed it, 515
lines: seven paired branches in `legmod`, four in `mpimod`, three blocks each in
`utilities` and `utilities_omp`, the un-permute in `plasim.f90`, `ilatperm`
itself, the `LPAIRLAT` parameter and `shtnsmod`'s refusal that could no longer
fire. The reason is the one the last section below anticipates from the other
side: SHTns replaced `legmod` and cannot use a permuted layout, so the
precondition this change existed to establish had nothing left to unblock. The
deletion is bit inert and that was checked -- the T21 restart sha is
`a7418dea572f853f` on both sides of it. `NLHP = NLPP / 2` survives in
`plasimmod.f90:108` as a plain parameter.

world-38b has since removed the MPI build entirely, so `mpimod.f90` and the four
permutation sites named below no longer exist either.

The record is kept for two things it establishes independently of the layout:
what a survey of "routes through the primitives" misses, and what the two halves
of the symmetric transform were measured to be worth.

Four gates were left pointing at the deleted code, and the four verdicts are
not the same verdict:

| gate | what became of it |
| --- | --- |
| `verify_latitude_pairing.py` | DELETED with its pipeline row under `world-h8o`. It lifted `ilatperm`; there was no function left to point it at, and repointing it would have made it pass vacuously |
| `verify_symmetric_transform.py` | RENAMED and re-scoped to `verify_inverse_transform.py` under `world-2div`. Its mirror controls perturbed lines that no longer exist, but the half that never depended on the layout survives and is worth more than it was: the fork's per-mode factor hoist is live in those loops |
| `verify_omp_collectives.f90` | REPAIRED under `world-8y1e`. It called `ilatperm` and so failed in the compiler rather than with a message. Only the placement block was broken, and the scatter is contiguous again, so that block asserts the contiguous placement instead |
| `verify_legendre_parity.py` | DELETED with its `config/pipeline.yaml` and `exoplasim/README.md` rows under `world-x3zd`. It checked the parity premise, whose only consumer is the deleted parity split; and it checked a Python transcription of `legini` rather than `legini`, so it could pass on a recurrence the model no longer had. `symmetric-transforms.md` keeps what it measured |

## The change

A latitude and its mirror carry the same associated Legendre values up to a
sign, so a transform that holds both on one process does half the multiplies
and reads the weight matrix once instead of twice. PlaSim already has code for
this in its four FORWARD transforms and production never reaches it, because
`mpimod` scatters a CONTIGUOUS block of latitudes to each process and a mirror
pair therefore lands on two different ones. The existing branch is guarded by
`NLPP < NLAT`, which is only false on one process.

So the scatter is permuted instead. Process `r` gets the northern block
`r*NLHP+1 .. (r+1)*NLHP` followed by those latitudes' mirrors in reverse, which
makes local `l` and `NLPP+1-l` a mirror pair on every process. The gather undoes
it, so nothing above `mpimod` can tell which layout it is running under.

Two new compile-time parameters carried it, both in `plasimmod.f90`:

    NLHP     = NLPP / 2
    LPAIRLAT = (mod(NLAT,2*NPRO) == 0)

`LPAIRLAT` is a `parameter` rather than a runtime switch on purpose: the
alternative branch is then dead code the compiler deletes, and legmod's
selection costs nothing in the hottest loop in the model. Where `NPRO` does not
divide `NLAT/2` the permutation is the identity and the stock contiguous layout
stands -- T21 on 32 processes is the case that reaches it. **At `NPRO == 1` the
permutation is ALSO the identity**, and legmod's rewritten branch generates
exactly the code that was there before, which is what makes the single-process
executable the trivial case rather than a special one.

## Where the permutation lives, including the two sites the plan missed

Every grid-space transfer routes through five primitives, and the survey that
said so was right about those five and incomplete about the rest.

| site | what it does |
| --- | --- |
| `mpimod.f90:mpscgp` | permutes the global array before `mpi_scatter` |
| `mpimod.f90:mpgagp` | un-permutes after `mpi_gather` |
| `mpimod.f90:mpgallgp` | un-permutes after `mpi_allgather`, on every process |
| `mpimod.f90:mpgacs` | un-permutes the NLAT axis of a cross section |
| `plasim.f90:prolog` | permutes `sid`, `gwd`, `csq`, `rcs` before their scatters |
| `utilities.f90:makeareas` | **missed by the plan.** `mpgarn` is a raw gather, and the cell bounds are built from ADJACENT gathered latitudes, which only means anything in global order. Un-permutes the gather and permutes `lat1`/`lat2` back before their scatter |
| `utilities.f90:finishlat` | **missed by the plan.** Another raw gather, and this one writes a latitude record straight to a file |

`mpreadgp`, `mpwritegp`, `mpwritegph`, `mpgetgp`, `mpputgp` and `mpsurfgp` all
delegate to `mpscgp`/`mpgagp` and needed nothing, as did the semi-Lagrangian
tracer scheme, which gathers to a global array before it does anything
meridional. `deglat`, `cola` and `rcsq` derive locally from the four scattered
vectors and follow automatically; so do legini's weight matrices, which is the
whole reason the four vectors are the right place to put the permutation.

**The lesson from the two misses is that "routes through the primitives" was
checked and "is a collective at all" was not.** The survey that found them is
`grep` for `mpi_gather|mpi_scatter|mpi_allgather` outside `mpimod.f90`, which
takes a second and should have been the first thing run.

## What is deliberately NOT permuted

`oceanmod.f90` scatters four latitude vectors -- `cphi`, `cphih`, `dphi`,
`dmue` -- and they stay in global order. `mpscrn` scatters in place, so the
root's first chunk is written back over itself and its copy of the global array
survives; `hdiffo` is the only reader of any of them, runs on the root alone,
and indexes them by GLOBAL latitude beside the fields it diffuses. The
scattered values on the other processes are never read at all. Permuting these
would leave `hdiffo` reading a reordered array through global indices, which is
silently wrong rather than an error, and `nhdiff = 1` is a live setting in this
project (CLIM-16). The four calls carry a comment saying so.

## `ilatperm` belongs in `pumamod`

It was written into `mpimod.f90`, which left the single-process executable with
an unresolved reference: `most_compiler` selects `mpimod_stub` and
`utilities_stub`, and neither has it. As a module procedure of `pumamod` it is
compiled by both builds and needs no external declaration at any call site.
This is worth recording because the failure is invisible on the configurations
this project runs -- every registered binary is 8 processes or more -- and
would have surfaced only for someone building the model at one.

## The tests, and what each of them can fail

**`verify_latitude_pairing.py`** lifted `ilatperm` verbatim out of
`plasimmod.f90`, compiled it inside a parameter module for one `(NLAT, NPRO)`,
and checked the map it produced. It is DELETED, with the function it tested.
What it checked: bijection, mirror pairing, the northern block contiguous and
northernmost, identity at one process, identity where the divisibility fails.
Seventeen cases across the resolution ladder.

It carries TWO negative controls, because the two ways the function goes wrong
fail different checks:

| control | what it does | what catches it |
| --- | --- | --- |
| stride | northern block advances by `NLPP` instead of `NLHP` | latitudes collide: not a bijection |
| shift | southern blocks handed round the processes by one | still a bijection, still `NLPP` latitudes each -- only the PAIRING is wrong |

The shift control is the one that matters. A bijection check alone waves it
through, and it is the mistake this function actually invites.

**`verify_paired_decomposition.sh`** built the paired and contiguous layouts
from the same source and ran them on the same bed, at the same process count.
It is deleted with the layout it compared. It could not be a checksum
comparison and it is worth saying why, because the next comparison across a
decomposition change meets the same thing: each process sums the transform over
the latitudes it holds and `mpsumsc` completes the partial sums, so changing
WHICH latitudes a process holds regroups a floating point sum and the last bits
move legitimately. What separates a regrouped sum from a different computation
is size, not exactness.

Two things in it failed hard, and both are the shape to reuse.

- **The round-tripped fields must be bit identical.** `dls`, `doro` and `darea`
  are scattered in and gathered back with no arithmetic between. A permutation
  applied on the way in and not on the way out shows up there and nowhere else,
  and `darea` is the sharpest of the three because it is the Gaussian weight of
  a process's OWN latitude and only survives if the `sid`/`gwd` scatter is
  permuted the same way the grid fields are.
- **The shift control must be rejected by the same comparison.** It is, and it
  is rejected over one step, because over twenty it reaches a floating point
  exception in `fluxstep` instead of producing a restart at all. Running the
  control short is what keeps the comparison exercised rather than bypassed by
  a crash.

### What the runtime test did NOT cover

Two of the five primitives were correct by construction and unexercised. The
permutation they would have carried is gone, so the gap is closed by deletion
rather than by a check, but the two facts underneath it are still true of the
threaded layer and still decide whether anything new needs its own gate.

- `mpgallgp` **has no call sites at all** in the model as it stands. Its
  permutation is written to match `mpgagp` and nothing runs it.
- `mpgacs` runs on the diagnostic interval, but a cross section is only ever
  handed to `guiput`, which is a stub without the GUI. A wrong permutation
  there cannot reach a restart, so this test would not see it.

Anything that starts using either one needs its own check first.

## The result

T21, 8 processes, 20 steps, tolerance 1e-10 declared before the arms ran.

| | records | worst relative difference |
| --- | ---: | --- |
| bit identical | 143 of 199 | -- |
| at rounding scale | 56 | `dcc` at 2.2e-11 |
| beyond it | 0 | -- |

and the control, over one step:

| | records | worst relative difference |
| --- | ---: | --- |
| bit identical | 155 of 199 | -- |
| at rounding scale | 3 | -- |
| beyond it | 41 | `sqm` at **21.5** |

**Twelve orders of magnitude separate the passing arm from the failing one**,
which is what makes the pass mean something: the verdict does not depend on
where in that gap the tolerance was put. The tolerance itself is a 20-step
bound and nothing more -- two integrations that differ in the last bits diverge
on their own from there, and a longer arm would legitimately exceed it.

## What it is worth

Measured on its own, before any of Phase 3 was written, so that a later result
is attributable. T127, 16 processes, 4 interleaved rounds, the same bed and
binary pair throughout:

| | contiguous | paired | gain |
| --- | ---: | ---: | ---: |
| T127, 300 steps | 39.61 s | 39.11 s | **+1.43%** [+0.87, +2.07] |

Faster in 4 of 4 rounds, self-scatter 1.0% and 0.7%. That is the forward
symmetric branches becoming reachable in production, net of the permutation
copies, and it lands where the 2% ceiling implied by
`spectral-transform-profile.md`'s forward-direction share put it.

At T21 and T42 this half was never resolved on its own and is not quoted: six
and nine second beds gave 7 to 12 percent self-scatter, over the 5% floor. What
IS resolved at the low end is the two halves together, in
`symmetric-transforms.md`, which is the number that matters for the decision
anyway.

## What this half is FOR

+1.43% would not justify a decomposition change on its own, and the plan said
as much before any of it was written. Its value is that it is the precondition:
with mirror pairs on one process the inverse transforms can use the same
symmetry, and those are 15% of samples at T127 against the forward direction's
7%. Together the two halves are worth **+6.92% at T127 and +11.65% at T170** --
see `symmetric-transforms.md` for that table and for what the sixteen
accumulators in `dv2uv` turned out to cost.

## What this unblocks

`NLPP` is now the number of latitudes a process holds AND an even number of
them arranged in mirror pairs, which is the precondition for the symmetric
INVERSE transforms. Those are `sp2fc`, `sp2fcdmu` and `dv2uv`, **33.7% of
samples at T127** against 4.2% for the whole forward direction, and they have
no symmetric path at all. The derivation is in `symmetric-transforms.md`.

The layout matrix in `smt-rank-layout.md` settled production at 16 processes,
and 16 divides `NLAT/2` at every resolution on the ladder, so the constraint
this change imposes costs nothing where the model actually runs.
