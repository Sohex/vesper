# The model is run-to-run reproducible; the rank count is what is not free

Measured 2026-08-20 for CLIM-44, on `exoplasim/scripts/reproducibility_matrix.py`.
The artifact is `exoplasim/analysis/reproducibility_matrix.json`.

Worldbuilding frame: this is about the Vesper simulation's climate model and
whether it produces the same bytes twice. Nothing here is a claim about the
simulated world.

## What is true

**At a fixed rank count the model is bit-reproducible.** Eight configurations --
8 and 16 ranks, output off and at the production setting, segments that cross no
write and that cross two -- three repeats each, from one restart through one
binary. `plasim_status`, `plasim_output`, `plasim_snapshot` and the model half of
`plasim_diag` are identical across every repeat of every cell.

**The rank count is not free.** 8 ranks and 16 ranks integrate the same restart
to different answers. After 56 timesteps, 83 of the restart's 199 records differ,
all of them at round-off:

| record | max abs difference / max abs value |
| --- | ---: |
| `dcc` | 3.737e-10 |
| `acc` | 5.636e-11 |
| `asthr` | 2.708e-11 |
| `atsolu` | 6.928e-12 |
| `assol` | 6.866e-12 |

That is a reduction-ordering signature -- the decomposition changes which partial
sums are formed in which order -- and it is benign per step and unbounded per
run, because this model is chaotic. **So an A/B must hold the rank count fixed.**
It is not enough to hold the binary fixed: ExoPlaSim compiles one executable per
(resolution, layers, ranks) triple anyway, so a rank change is already a binary
change, and this is the reason it is also a numerics change.

**Writing output does not perturb the integration.** At each rank count the
`NOUTPUT = 0` and `NOUTPUT = 1` cells produce the same restart, byte for byte.

## What CLIM-39 hit, and what it did not

CLIM-39 reported gridpoint output bit-identical at 1 timestep and different by
16, and `plasim_status` different between two runs of one binary even at 1
timestep. **Neither reproduces here**, on either tree:

| tree | `zsolars` record | reproducible at fixed rank count |
| --- | ---: | --- |
| current, `aeafe3c` | 16 bytes | yes, all 8 cells |
| the pre-fix binaries a run directory keeps | 65536 bytes | yes, all 8 cells |

The second row is the one that matters, and it refutes the obvious explanation.
The `zsolars` overread of `notes/audits/zsolars-restart-overread.md` was live in
CLIM-39's tree -- `git merge-base --is-ancestor 80f11e9 b43f17e` is false -- and
that record really did carry 8190 doubles of memory past a 2-element array. But
under the shipped `-finit-real=zero` it is **deterministic**: measured here, its
surplus is stable across three repeats at a fixed rank count and differs only
between rank counts. So the overread is a source of RANK-dependence, not of
run-to-run difference, and "it was the overread" is refuted as the cause even
though the overread is real.

What remains is that CLIM-39's tree was a worktree build. The same audit records
`-flto` making exactly those bytes vary run to run, five runs and five hashes
from one binary. A build difference is therefore the surviving candidate, and it
is not pursued further: the binaries that tree produced are gone, no current
configuration shows the behaviour, and the fix that removes the mechanism
entirely is already resident.

**Its control could not have shown what it was read as showing.** A gridpoint
record is written on `mod(nstep, nafter) == 0` over the ABSOLUTE step count. This
project's production restart carries `nstep` = 589672 and the model logs `nafter`
= 32, so the next write is 24 steps away -- and neither a 1-step nor a 16-step
segment reaches it. Measured on that restart, a 23-step segment writes a
`plasim_output` of 32816 bytes, one record of one constant field; a 56-step
segment writes 13,413,736. A comparison over the former cannot fail.

This does not establish what CLIM-39's own bed did, whose namelist is not
recoverable. It establishes that the recorded control does not distinguish "the
output agreed" from "there was no output", and the reduction identity that
CLIM-39 states at one timestep rests on the same comparison.

## Method, and what could have failed

Three hypotheses with their patterns were declared in the script before it ran,
and the matrix could have shown any of them: output-writing, rank reduction
ordering, or everything reproducing. The third is what showed, and it is
reported as a PATTERN rather than a mechanism for the reason above.

Two traps were hit and are now closed in the harness rather than in prose:

**`plasim_diag` is not a comparison target as it stands.** It ends with the run's
own resource accounting, so two identical integrations differ in it by
construction. Filtering it by reading a diff is not enough: every line of that
block is conditional on its counter being non-zero, so a filter built from one
pair of runs passed three cells and failed three more when "Page faults"
appeared. The filter is now taken from the block that writes it,
`plasim.f90:1011-1042`. With it, the difference is exactly zero.

**A bed copied out of `exoplasim/runs/` carries the binaries that run STARTED
with.** The first pass of this matrix measured a tree three commits behind the
one `rebuild_binaries.py --verify` calls current, and nothing about the result
said so -- it ran, it was self-consistent, and it was about the wrong code.
`--verify` answers "are the installed binaries current" and cannot answer "is the
bed running them". The harness now hashes the bed's executables against
`exoplasim/binary_manifest.json` and refuses to start otherwise.

That accident is the reason the pre-fix row in the table above exists at all, and
it is the measurement that refutes the comfortable answer. It is recorded because
the useful half of it was not planned.

## What this changes

The standing instruction CLIM-39 left -- that any A/B expecting bit-identity must
first establish its own reproducibility horizon with a same-binary control -- is
**narrowed rather than withdrawn**. Bit-identity is available. What it requires
is the same rank count, and a same-binary control is one cheap run rather than a
bound carried on every claim.

Not measured here: rank counts other than 8 and 16, resolutions other than T42,
and segments longer than 56 timesteps. The build-flag benchmark covers the last
of those in one configuration, at 600 steps on 16 ranks with output off, where
`plasim_status` came out bit-identical across about 48 repeats spanning four
binaries (`notes/audits/aocl-and-model-build-flags.md`).
