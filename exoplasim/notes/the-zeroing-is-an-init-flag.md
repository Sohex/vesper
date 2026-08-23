# A third of T170 was `-finit-real=zero`, and the model does not need it

*Worldbuilding frame: a COMPUTE measurement of the Vesper project's climate
model on this desktop. Nothing here is about the simulated planet. Measured
2026-08-22 at e1733979, T170, sixteen threads, on `bench/bed_t170cold`.*

`shtns-viability.md` measured `__memset_avx512_unaligned_erms` at about 31% of
T170 and attributed it, with frame pointers, to `radstep`, `rainstep`,
`fluxstep`, `kuo`, `icestep` and `gridpointd` -- every physics routine clearing
its own work arrays, per step, with no single site worth attacking. The
attribution was right and the CAUSE was not. The model is not clearing those
arrays. **The compiler is**, because `MOST_F90_OPTS` carries
`-finit-real=zero`, which initialises every local real variable on entry to
every routine, arrays included.

## What the flag is worth

Two T170 threaded binaries from the same source, differing in that flag alone,
paired and interleaved on the 300-step bed:

| arm | median | spread | restart |
| --- | ---: | ---: | --- |
| `-finit-real=zero` | 70.55 s | 2.1% | `9808bc54507629ad` |
| flag removed | 51.97 s | 6.1% | `9808bc54507629ad` |

**+26.03% [+22.76, +27.29], faster in 4 of 4 rounds**, and the restart is
IDENTICAL. The 6.1% self-scatter is above the 5% floor, and the gain is four
times it.

It is also visible statically: 1193 memset call sites in the binary with the
flag, 623 without.

## The two zeroings are different, and only one is the flag's

This matters because the earlier reading -- the tendency arrays are
ACCUMULATORS, so their zeroing is load-bearing and cannot be deleted -- is
correct and is about the other one.

`gudt(:,:)=0.` and its eight siblings at the head of `gridpointd` are
STATEMENTS IN THE SOURCE. The physics adds into them all timestep, the values
are live, and no compiler flag removes them. They are 0.655 MB each at T170.

What the flag adds to that routine is exactly three arrays:

| `gridpointd_`, static memset bytes | |
| --- | ---: |
| with `-finit-real=zero` | 35.062 MB |
| without | 3.604 MB |
| difference | 31.457 MB, over exactly 3 calls |

31.457 MB is three `(NLON,NLAT,NLEV)` arrays at 10.486 MB each: `zgq`, `zmmr`
and `znrho`, full-globe scratch locals, one copy per thread. `zgq` is not read
at all in this configuration -- its branch needs `nqspec == 0`, and both the
default and the bed are 1. So the flag zeroes 503 MB a timestep across sixteen
threads in this routine alone, and a third of it is for an array nothing looks
at.

The same split runs through the physics. Static bytes the flag ADDS, per call,
for the routines that run every step:

| routine | added |
| --- | ---: |
| `swr_` | 36.4 MB over 115 sites |
| `gridpointd_` | 31.5 MB |
| `lwr_` | 14.0 MB |
| `vdiff_` | 11.7 MB |
| `kuo_` | 6.7 MB |
| `mkshallow_`, `mkdca_`, `mkradv_`, `mkclouds_`, `hdiffo_`, `calcgp_` | 2 to 7 MB each |

Static, so it is bytes per CALL at each site rather than per step, and sites
inside branches that do not run are counted. It is the shape that matters: the
cost is spread across every physics routine because the flag is, which is
exactly the profile that read as "diffuse, nothing to attack".

## Whether the model needs it: three gates, and it does not

**The trapping gate, `-finit-real=snan` with `-ffpe-trap=invalid`.** A read of
an uninitialised local then traps at the instruction that uses it. On the legmod
path the model is CLEAN: 60 steps at T21 and at T170, at `NLOWIO=0` and
`NLOWIO=1`, no trap.

On the SHTns path it traps at step 1, and NOT in the model -- see below. That
means this gate cannot certify the SHTns path, because the library trips before
the gate can say anything about the model.

**The propagating gate, which answers what the trapping one cannot.**
`-finit-real=snan` with FE_INVALID left MASKED: an uninitialised read then
produces a quiet NaN and propagates into the restart instead of aborting. A NaN
that reached any prognostic field would spread and change the file. T170, 60
steps, `NSHTNS=1`, `NLOWIO=1`:

| arm | restart |
| --- | --- |
| `-finit-real=zero` | `017a35afbdf1f5f6c621c734` |
| `-finit-real=snan`, INVALID masked | `017a35afbdf1f5f6c621c734` |
| no init flag | `017a35afbdf1f5f6c621c734` |

Bit identical. Nothing uninitialised reaches a stored value.

**The 300-step restart**, above, agrees at the length the bench runs.

## The trap is SHTns reading its own scratch

Located rather than inferred. Under `-finit-real=snan` the fault is in
`_an12_l`, called from `spat_to_SH_fly2_l`, called from `sh_gp2sp` at
`shtnsmod.f90:367`, called from `gridpointd` at `plasim.f90:3717` -- the
temperature-tendency analysis, at the first timestep. It is not the `zqout`
path: it traps identically at `NLOWIO=0`.

The faulting instruction is the tail of a horizontal reduction:

    vshuff64x2 $0x1,%ymm0,%ymm2,%ymm1
    vblendpd   $0x3,%ymm2,%ymm0,%ymm0
    vaddpd     %ymm0,%ymm1,%ymm0
    vaddpd     -0x20(%rsi),%ymm0,%ymm0     <-- faults

and the address it reads, `0x1554b3f2f2e0`, is below both of the wrapper's
locals -- `zlm` at `0x1554b3f347c0` and `zg` at `0x1554b3f357c0` -- so it is
inside `_an12_l`'s own frame, some 21 KB deeper in the stack. `zg` holds real
data at the fault and `zlm` is half written.

Three candidates were excluded by measurement rather than by argument:

**The model does not hand SHTns an undersized buffer.** `shtns.h` documents the
spatial argument as "a double array of size shtns->nspat", and
`probe_shtns_nspat.c` asks the library: with `SHT_PHI_CONTIGUOUS`, nspat is
exactly NUGP at T21, T42, T127 and T170. No overrun.

**It is not an unwritten m=0 imaginary part.** In the direction where `zlm` is
an INPUT, `sh_sp2gp` assigns every element from 1 to NCSP and sets the m=0
imaginary part to an explicit zero. In the direction that traps, `zlm` is
SHTns's output.

**It is not the accumulator argument in disguise.** The fault address is in no
model array.

What remains is SHTns reading a lane of its own scratch that it has not written
-- a vector tail over-read, on the reading the disassembly supports. **The init
flag neither causes it nor cures it**: by the time SHTns runs, the stack below
the Fortran frame holds leftovers from previous frames in either build, and the
flag only zeroes locals at entry to routines that then overwrite them. The
exposure is that `-ffpe-trap=invalid` sets MXCSR for the whole process, so
SHTns runs under the model's trap mask whether or not its authors assumed
masked exceptions, and its padding read is safe only while the stack garbage
under it is not a NaN.

## Adopted, and what the adoption verified

`config/planet.yaml` no longer declares `-finit-real=zero`, and `checked`
carries `-finit-real=snan`. The threaded build compiles byte for byte to
`c8b9c41046fd93ca`, which IS the arm the 26.03% was measured on, so the
declaration delivers the measured thing rather than something believed
equivalent to it.

Determinism was the one thing this could have cost, archive CLIM-44 being this
model failing to reproduce and having been paid for once. It survives:
`verify_shtns_model.sh` passes whole on the new build, both arms at rounding
scale, the control rejected, and the four-run bit-identity arm giving one hash,
`fcf46ebbb1302d3a`.

`checked` is usable because `rebuild_binaries.py` builds MPI binaries, which
cannot run SHTns and therefore never meet the trap below. On a threaded build
snan still trips on the library first, which is CLIM-82.

## What this retires

The conclusion that there is no low-hanging fruit in the zeroing, and the
reading that only a wholesale restructuring of the physics could reach it. The
largest single cost in the model was one flag, it is worth 26% at T170, and the
model computes the same answer without it.

It does not retire the case for cache-resident column blocks. It resizes it:
blocking's first-order claim was that it moves this zeroing into L1, and most of
this zeroing should not happen at all. What survives is the second-order claim
-- that `radstep`, `rainstep`, `fluxstep`, `kuo` and `icestep` each stream the
whole band independently and a resident block fuses them -- and the work queue
it would supply.
