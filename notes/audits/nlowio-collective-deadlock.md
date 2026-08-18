# NLOWIO = 0 deadlocks a patched model, and it was read as a compiler regression

Measured 2026-08-18.

## What is true

`nlowio` is read from `plasim_nl` inside `if (mypid == NROOT)` in `prolog`
(`plasim.f90`, `call readnl`) and is **never broadcast**. Its compiled default is
1 (`plasimmod.f90:151`). So when the namelist says `NLOWIO = 0`, rank 0 holds 0
and every other rank holds 1. `nafter`, `noutput`, `nsnapshot`, `nstps` and the
high-cadence keys beside it are all broadcast; `nlowio` and `nstpw` are the two
that are not.

That divergence is harmless on its own. It becomes a deadlock when a collective
is placed behind a condition that reads `nlowio`, and
`exoplasim-3.4.2-lowio-first-record.patch` does exactly that. It rewrites, in
`gridpointd`:

    -      if (mod(nstep,nafter)==0 .and. nqspec == 1) then
    +      if ((nlowio > 0 .or. mod(nstep,nafter)==0) .and. nqspec == 1) then
            do jlev=1,NLEV
             ... gp2fc / fc2sp ...
            enddo
            call mpsum(sqout,NLEV)          <-- COLLECTIVE
           endif

`nafter` is broadcast, so the original condition was the same on every rank. The
patched one is not: at `NLOWIO = 0` the non-root ranks take it EVERY timestep and
NROOT takes it every `nafter` steps. The ranks desynchronise on the first step,
and from then on the non-root `mpsum` (a Reduce) meets NROOT's later `mpsumbcr`
(an Allreduce) in `fixer`. The run never advances.

The change itself is correct and is the point of the patch -- under low I/O,
`outaccu` sums `sqout` every timestep, so `sqout` has to be current every
timestep. What was missing was the broadcast that makes `nlowio` mean the same
thing on every rank.

## Evidence

`gdb` on the live ranks of a stuck run, all inside `gridpointd_`, which is where
the patched condition sits:

    ranks A: MPI_Allreduce <- mpsumbcr_ <- fixer_ <- miscstep_ <- gridpointd_
    ranks B: MPI_Reduce    <- mpsum_    <- gridpointd_

`/proc/<pid>/stat` on a stuck rank: utime 19011 ticks against stime 9. Pure user
time, because OpenMPI busy-waits in `opal_progress`. A deadlock here presents as
sixteen pinned cores, so "the machine is busy" is not evidence that work is
happening. `plasim_output` stops at 40 bytes, the header, and never grows.

The decisive control is the binary `238af7e8`, which ran the successful
`NLOWIO = 0` orbits and does NOT carry the low-I/O patch. It completes
`NLOWIO = 0` normally. Adding the patch is what breaks it; adding the broadcast
under the patch fixes it again.

## Why it did not bite sooner

Two halves, introduced seven hours apart, neither harmful alone:

| when | what |
| --- | --- |
| 2026-08-17 15:21 | `c57537f` adds `disable_low_io()`, so runs start setting `NLOWIO = 0`. Nothing before this ever set it, so `nlowio` matched the default on every rank and the divergence did not exist. |
| 2026-08-17 22:50 | `e3ef2f3` makes `lowio-first-record` resident, which is what puts a collective behind `nlowio`. |
| 2026-08-17 23:56 | The post-equilibrium climatology runs `NLOWIO = 0` on binary `238af7e8`, built BEFORE the patch. It completes: 90 s of model time per orbit. |
| 2026-08-18 01:59 | System update, gcc-fortran 16.1.1 to 16.2.1. Irrelevant, and this is the coincidence that misdirected the diagnosis. |
| 2026-08-18 onward | Every run uses a binary carrying both halves. Every one of them deadlocks. |

So the collision needed a binary built after 22:50 AND a namelist written after
15:21, and the first run to have both was on 2026-08-18.

## Timings, T42 p16, one full orbit of 5850 steps

Identical namelists and inputs, sequential, no competing load, 16 ranks bound one
per physical core.

| binary | patches | NLOWIO | wall | CPU |
| --- | --- | --- | --- | --- |
| pre-update, gcc 16.1.1 | with low-I/O | 1 | 77.7 s | 1567% |
| post-update, gcc 16.2.1 | with low-I/O | 1 | 78.2 s | 1567% |
| private-prefix, gcc 16.1.1 | with low-I/O | 1 | 78.0 s | 1567% |
| pre-update, gcc 16.1.1 | with low-I/O | 0 | deadlock, killed at 300 s | busy-wait |
| post-update, gcc 16.2.1 | with low-I/O | 0 | deadlock, killed at 300 s | busy-wait |
| `238af7e8`, gcc 16.1.1 | WITHOUT low-I/O | 0 | runs normally | -- |
| current, gcc 16.2.1 | + broadcast, no low-I/O | 0 | 79.3 s | 1567% |
| current, gcc 16.2.1 | + broadcast + low-I/O | 0 | 89.7 s | 1568% |
| current, gcc 16.2.1 | + broadcast + low-I/O | 1 | 84.0 s | 1567% |

89.7 s reproduces the 90 s per orbit that `run_manifest.json` recorded for the
2026-08-17 `NLOWIO = 0` segment, which is the check that the fix restores the
behaviour rather than merely stopping the hang.

The low-I/O patch costs 84.0 s against 78.0 s at `NLOWIO = 1`, about 8%, because
the humidity transform now runs every timestep instead of every `nafter`. That is
what the patch is for and it is the honest price of it.

## What was believed instead, and why it was wrong

It was recorded that the low-I/O patch "roughly triples model runtime" under
gcc-fortran 16.2.1 by adding `use restartmod` to `plasim.f90` and pessimising the
translation unit, and a private gcc 16.1.1 prefix was built from the pacman cache
to escape it. Every part of that is wrong. The three compilers give 77.7, 78.2
and 78.0 s on the same orbit. The patch costs 8%, not 200%, and it costs it
through a per-timestep spectral transform that is visible in the diff. The
pre-update binary deadlocks exactly like the post-update one.

What was almost certainly measured is the wall-clock jump from `NLOWIO = 1` to
`NLOWIO = 0` -- 104 s to 387 s per orbit, recorded independently in
`continue_exoplasim.py` -- which is pyburn chewing 2.4 GB an orbit instead of
96 MB, single-threaded, with model time unchanged at 88 to 90 s. That is a real
cost, it is not a regression, and it is not the compiler either.

The test that separated all of this takes 78 seconds: run the OLD binary under
the NEW conditions. Two changes landed within a day of each other and the one
that was easier to name got the blame.

## The second defect: silent corruption of every NLOWIO = 0 field

CONFIRMED 2026-08-18 by A/B experiment, having first been deduced from source.

`outgp` and `outsc` are entered by every rank together, since the guards above
them (`mod(nhcstp,nafter)`, `noutput`) are broadcast. Inside, the output fields
sit in `if (nlowio .eq. 0) ... else ... endif` pairs and BOTH branches call the
collective `writegp`: the first writes `dmld`, `dt`, `dwatc`, `dsnow`, the second
writes `aadmld`, `aadt`, `aadwatc`, `aadsnow`. With `nlowio` divergent the
collective COUNTS still match, which is why this never hung -- but NROOT
contributes its INSTANTANEOUS slice to the same `mpgagp` that the other fifteen
ranks feed ACCUMULATED slices into. The assembled global field is part sample and
part interval mean.

### The experiment

Two binaries from one source, differing in the broadcast patch and nothing else:
A carries it, B does not; neither carries the low-I/O patch, so neither
deadlocks. Both run 320 steps at `NLOWIO = 0` from the same restart, 16 ranks.

    A run twice   -> byte-identical output. The comparison is deterministic.
    A against B   -> 474 differing data records out of 8902. Every HEADER
                     identical, so the record structure is untouched and only
                     values move.

### The signature is exactly the predicted one

For code 139, surface temperature, the latitude rows identical between A and B
are rows 0, 1, 2 and 3, and no others. `NLAT/NPRO = 64/16 = 4`, and `mpgagp`
assembles rank r's slice at global offset `r*NHOR`, so rows 0-3 ARE rank 0's
share. Rank 0 is NROOT, which reads the namelist and holds `nlowio = 0` in both
builds, so its rows cannot differ. The other 60 rows are ranks 1-15, which held
the default 1 without the broadcast and wrote accumulated values there.

**93.75% of the grid was wrong in every affected field.**

### What is affected, and by how much

Only the fields written inside those if/else pairs. Fields written unconditionally
beside them -- `aprl` 142, `aprc` 143, `ashfl` 146, `alhfl` 147, `aroff` 160,
`acc` 164, `atsa` 167, `ats0` 169 -- are byte-identical in A and B, which is the
internal control on the result.

| code | field | max difference |
| --- | --- | --- |
| 139 | surface temperature | 14.93 K |
| 183, 207 | | 7.19, 3.26 |
| 162 | cloud cover, per level | 0.99 (i.e. clear against overcast) |
| 208, 210 | | 0.94, 0.97 |
| 174, 175 | albedo | 0.57, 0.46 |
| 268, 269 | | 1.60 |
| 159 | u-star cubed | 0.52 |
| 140, 141 | soil wetness, snow | 0.048, 0.073 |
| 161, 265, 298 | liquid water and others | 1e-4 and below |

Codes 210 and 211 differ on only 7 and 15 rows because they are zero over most of
the domain; where they are non-zero they differ like the rest.

### What it means for data on disk

Every `NLOWIO = 0` field written by a binary without the broadcast is affected.
That is orbits 66 to 76 of `run_8c2e1ff9ab5e` -- the high-cadence diagnostic and
the whole post-equilibrium climatology -- and anything derived from them. The
convergence scalars come from codes 167 and 169, which are outside the if/else
and are clean; the surface, cloud, albedo and roughness fields are not. Tracked
as CLIM-20.

Note the direction, because it is the opposite of what was intended: the point of
`NLOWIO = 0` was to obtain instantaneous samples rather than accumulations, and
what it produced was accumulations over 60 of 64 latitude rows.

### What is NOT affected: CLIM-1 and DUST-5 are both clear

Checked structurally and confirmed on the same A/B output.

**Only the regular stream can carry this.** `nlowio` appears in `outsc`, `outsp`,
`outgp`, `outdiag`, `outreset` and `outaccu`, and nowhere else in `outmod.f90`.
`snapshotgp`, `snapshotsc`, `snapshotdiag`, `hcadencegp` and `hcadencesc` do not
reference it at all, so the snapshot and high-cadence streams are structurally
immune whatever `nlowio` says.

**CLIM-1, the energy closure, is clear.** `outdiag` has no `nlowio` branch: the 28
energy terms accumulate unconditionally in `outaccu` and are divided by
`naccuout`, which is incremented at outmod.f90:2650 OUTSIDE the `if (nlowio > 0)`
guard and zeroed at 2501 outside it too. So `naccuout` is the same on every rank
and the terms scale identically. On the A/B output all 28 terms (codes 360-387),
all 28 three-dimensional terms (codes 460-487) and every comparand the closure
compares them against (318, 320, 321, 142, 143, 146, 147, 176-179) are present and
BYTE-IDENTICAL. The independent check is code 142 `aprl`, which divides by
`naccuout` and does not differ -- had `naccuout` diverged, it would have.

**DUST-5, the gust distribution, is clear.** Codes 131, 132 and 259 do not appear
in the raw regular stream at all; `pyburn` derives the winds from vorticity 138
and divergence 155, neither of which differs, and the gust fit reads the
high-cadence stream, which has no `nlowio` reference in its writer.

The 26 codes that DO differ are 139, 140, 141, 159, 161, 162, 170, 173, 174, 175,
183, 184, 207, 208, 209, 210, 211, 265, 267, 268, 269, 298, 322, 323, 324, 327 --
surface state, cloud, albedo, roughness and the entropy diagnostics.

## The fix

`exoplasim-3.4.2-nlowio-broadcast.patch` adds `call mpbci(nlowio)` and
`call mpbci(nstpw)` to `prolog`, beside the `mpbci(nafter)` already there.
`nstpw` is used off-NROOT nowhere today -- its one use, in `write_atmos_restart`,
is inside `if (mypid == NROOT)` -- so broadcasting it changes nothing now and
stops the same trap being reset later.

At `NLOWIO = 1` the patch is a no-op by construction: the non-root ranks were
already holding 1, which is what they now receive.
