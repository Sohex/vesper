# The EMBM-free surface-flux path, and what the ocean's threading cost to do

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean, sea ice, surface flux and wind stress
name modelled quantities and numerical properties of the adopted offline ocean
host, not observations of anything.

Established 2026-08-25 on this machine. This is the implementation half of
OCN-19 and the precondition half of OCN-10, and it follows
`notes/audits/cgenie-parallelism-and-coupling-support.md`, which profiled the
component and proposed three changes without making any of them. That note's
numbers are taken as given here; what is new is what happened when the changes
were made.

`analysis/cgenie_omp.py` is the driver. It exports the vendored tree out of the
repository and builds it there, because a worktree's ignored build products are
symlinks into the shared checkout.

**The acceptance test, fixed before any arm was built.** Each of these changes
is correctness-preserving by construction -- an array temporary removed, a
storage class changed, a loop divided over threads that carries no dependence --
so the bar is not a tolerance:

> Every float variable in a shipped regression case's reference netCDF must be
> BIT-FOR-BIT identical to what the tree at commit `4038c12e` writes.

Two cases carry it, and each is 20 model years -- `koverall_total` 10000 at
`kocn_loop` 5 is 2000 ocean timesteps at the shipped `nyear` of 100 -- so this
is a short-run regression surface and not an equilibrated comparison.
`eb_go_gs` is EMBM, GOLDSTEIN and the sea ice, compared on the 28 variables of
GOLDSTEIN's year-20 annual average. `eb_go_gs_ac_bg` adds ATCHEM and BIOGEM,
compared on the 88 variables of the BIOGEM three-dimensional fields, and it is
in the set for one specific reason rather than for coverage: its per-column
arrays are what gfortran's own `-fmax-stack-var-size` warning is about.

The reference is the tree's own output at `4038c12e` rather than the shipped
`genie-knowngood/` netCDF, because the shipped file was written by a different
compiler on a different host and the difference between it and the unedited tree
is not a property of any arm. Both are reported. Against the shipped reference
the unedited tree differs in two cells of `uvel` at 7e-27 of the field range and
in two BIOGEM fields at 3e-45, and every arm below reproduces those same
differences and no others.

---

# 1. Two flux paths, both gated on an atmosphere, and a third routine nothing calls

The question this settles is whether an EMBM-free configuration needs a new
surface-flux path or only a switch. It needs a new path, and less of one than
the shape of `genie.F` suggests.

## 1a. What the recipe actually gates

`genie.F` calls a surface-flux routine in exactly two places and each is behind
an atmosphere flag:

- `flag_ebatmos .and. flag_goldsteinocean` calls `surflux_wrapper`, which is
  EMBM's `surflux` in `genie-embm`.
- `flag_plasimatmos` calls `plasim_surflux_wrapper` inside the sea-ice block,
  which is `surflux_goldstein_seaice` in **`genie-goldsteinseaice`**.

`initialise_genie.F` defaults both flags false, and with both false three things
do not happen: no flux routine runs, so `latent_ocn`, `sensible_ocn`,
`netsolar_ocn` and `netlong_ocn` are never filled; no wind stress is produced,
because `ocean_stressx2_ocn` and its three companions are outputs of `embm` on
one path and of PLASIM on the other; and `istep_ocn` is never incremented,
because both of its increments sit inside those same two gates.

**The EMBM-free recipe runs.** A configuration naming only `goldstein` and
`goldsteinseaice` initialises both components, integrates 20 model years and
reaches the shutdown banner. So nothing in the build or the initialisation
requires an atmosphere; what is missing is the forcing and the ocean's own step
counter, which is a gate to add rather than a dependency to break.

`MODULE_NAMES` in `genie-main/makefile` is fixed, so `genie-embm` is compiled
and linked whatever the recipe says. Removing EMBM is therefore a RUNTIME
question and not a build one: with `flag_ebatmos` false, `embm`, `tstipa` and
`surflux` are never entered and their 12.5 to 19.9 per cent of the instructions
is never retired.

## 1b. The routine the third path should call already exists

`surflux_goldstein_seaice` is the right one and it carries no EMBM dependency:
it lives in `genie-goldsteinseaice`, and its arguments are FLUXES and transfer
COEFFICIENTS from outside -- incoming solar and longwave, a net heat flux over
open ocean, evaporation, and latent and sensible coefficients "for consistency
in energy" with whatever computed them -- rather than a state to recompute them
from. That is the contract the offline architecture wants, for the reason
section 5b of the parallelism note gives: a surface flux is nonlinear in the
state, so evaluating it on the atmosphere's fine grid and averaging is
conservative where recomputing it from a coarse mean is not.

It also carries the sea-ice surface temperature solve, which is why the host's
sea ice cannot be switched off on a driven path: the ice is the host's and the
surface energy balance over it is part of this routine.

## 1c. A fourth routine is compiled, wrapped, and reachable from nothing

`surf_ocn_sic` in `genie-goldstein` describes itself as "a surflux routine for
models using only c-GOLDSTEIN's ocean and sea-ice modules". It is in
`genie-goldstein/src/fortran/makefile`, it has `surf_ocn_sic_wrapper` in
`genie_loop_wrappers.f90`, and no line of `genie.F` calls that wrapper.

It is not the path to use, and the reason is the one that chose the
architecture. It takes atmospheric STATE -- lowest-level temperature, humidity,
pressure, height and winds -- and computes the turbulent fluxes from bulk
formulae on the OCEAN's grid, and it derives the wind stress from the
lowest-level winds itself. That is the coarse-grid flux evaluation the offline
exchange exists to avoid. It is recorded because it is a working, compiled
routine whose name makes it look like the answer.

## 1d. What the third path has to carry

The gate is small: advance `istep_ocn`, fill the coupling fields from the
supplied forcing, call `surflux_goldstein_seaice`. The fields are exactly the
ones `genie.F`'s PLASIM block accumulates and averages, and they are the input
to OCN-10's contract rather than something this path invents:

| field | what it is | where it goes |
| --- | --- | --- |
| `insolar_sic`, `inlong_sic` | incoming shortwave and longwave | the sea-ice surface balance |
| `netheat_sic` | net heat flux over open ocean | the open-ocean balance |
| `surft_atm_sic`, `surfq_atm_sic`, `surfp_atm_sic` | lowest-level temperature, humidity, pressure | the ice surface solve |
| `surf_windspeed_sic` | lowest-level wind speed | the transfer coefficients |
| `latent_coeff_atm`, `sensible_coeff_atm` | the atmosphere's own transfer coefficients | so the two sides agree on the flux |
| `latent_ocn`, `sensible_ocn`, `netsolar_ocn`, `netlong_ocn` | the open-ocean fluxes | in and out of the routine |
| `evap_ocn`, `precip_ocn`, `runoff_ocn` | the freshwater forcing | GOLDSTEIN's salinity |
| `ocean_stressx2_ocn`, `ocean_stressy2_ocn`, `ocean_stressx3_ocn`, `ocean_stressy3_ocn` | wind stress at the u and v points | GOLDSTEIN's momentum |
| `go_solfor`, `go_fxsw` | insolation and shortwave at the surface | BIOGEM |

The wind stress rows are where OCN-17's declared `scf` bracket enters. `scf`
scales the stress linearly in `goldstein.F`, so the ocean's stability ceiling
inherits that 1-to-3 bracket linearly and what comes out of any stability sweep
on this path is a bracket rather than a number.

---

# 2. The heap traffic is gone, and the disassembly is the evidence

`eosd` declared its two levels as `real t(2)` and `real s(2)`, and its only
caller holds them in a four-dimensional array, so `tstepo_flux.F` passed the
strided sections `ts1(1,i,j,k:k+1)` and `ts1(2,i,j,k:k+1)`. gfortran packs each
into a contiguous heap temporary before the call and frees it after, once per
wet cell per level per ocean timestep, which the parallelism note measured at
9.26 per cent of retired instructions at 36 x 36 x 16, 7.87 at 72 x 72 x 16 and
11.17 in the BIOGEM configuration.

`eosd` now takes `t1, t2, s1, s2` as scalars and the call site passes four
element references. The arithmetic below the argument list is character for
character what it was.

**Checked at the instruction level rather than inferred.** `objdump` over the
built executable finds zero `malloc@plt` and zero `free@plt` anywhere in
`tstepo_flux_`, against the two and three the parallelism note disassembled in
the basic block around the one call, and the one `call eosd_` is still there.

Two things went with it in the same pass, for the same reason.

`tstepo.F` filled `ts_t1`, `ts1_t1`, `rho_t1`, `ts_t2`, `ts1_t2` and `rho_t2`
every ocean timestep and nothing read them: their only consumers were the
commented-out OpenMP SECTIONS beside them. That fill loop is the whole of
`tstepo`'s own 0.56 per cent at the shipped grid.

`tstepo_flux_t` was a second copy of the whole tracer-transport routine, called
only from those same commented-out sections. It is deleted rather than carried,
because a duplicate of the routine the threading was about to restructure is a
guarantee that the two diverge silently.

---

# 3. `-fno-automatic` came out, and the enumeration did not have to happen

The parallelism note's section 3e sized the fallback: 390 declared names and 81
arrays across seventeen files, of which `surflux.F` alone carries 141 and 20.
None of that was needed.

`makefile.arc`'s gfortran SHIP block now compiles `-frecursive` instead of
`-fno-automatic`, and both shipped regression cases reproduce the tree at
`4038c12e` bit-for-bit. So nothing in the exercised paths depended on the
implicit SAVE that static storage supplies, which is what the note's enumeration
of explicit `DATA` and `save` statements predicted and this is the test of.

The TEST block and the Win32 `gfc` block go with it. `gfc` is gfortran under
another name and the same argument applies unchanged; this host has no `gfc`, so
that one is reasoned rather than tested.

## 3a. A stack raise is required, not precautionary

`-frecursive` implies an unlimited `-fmax-stack-var-size`, which is the point --
the same locals go on the stack where each thread has its own. With the default
8 MB limit the `eb_go_gs_ac_bg` regression case **segfaults before its first
timestep**, and completes with the limit raised. BIOGEM's per-column arrays are
exactly what gfortran's own `-fmax-stack-var-size` warning is about, which is
why that case is in the acceptance set.

`genie.job` now raises the shell's stack limit and exports `OMP_STACKSIZE`
before running the executable, because `OMP_STACKSIZE` is the limit for every
thread other than the master and its default is smaller still. A caller that has
already chosen one keeps it.

## 3b. Bit-for-bit is not the whole of what the flag changes

A static local sits in `.bss` and the loader zeroes it, so a routine reading one
of its own locals before writing it reads zero on the first call. An automatic
local is whatever was on the stack. That is a behaviour change no regression
comparison is guaranteed to reach, because the value it changes is one nothing
was supposed to read.

The `initpoison` arm compiles `-finit-real=snan -finit-integer=-2147483647`, so
an uninitialised local cannot pass for a plausible answer and a read before
write propagates instead of hiding. **It reproduces the same output
bit-for-bit** on both regression cases, which is the evidence that nothing in
the exercised paths reads a local it has not written.

That test is not a proof over the whole tree. It covers what 20 model years of
the two shipped configurations execute, which is EMBM, GOLDSTEIN, the sea ice,
ATCHEM and BIOGEM on their ordinary paths and none of their error or restart
branches.

---

# 4. The threading, and the one dependence that had to be removed

OpenMP over the ocean's `j` index in `genie-goldstein`'s `tstepo_flux`, `co`,
`coshuffle` and `velc`, with `-fopenmp` in the gfortran SHIP block. Every other
`!$omp` string in `vendor/cgenie` sits inside a comment -- `!` in free form,
`c$$$` in fixed form -- so enabling the flag globally makes directives of these
three files and of nothing else.

## 4a. Why the parallel axis is j

`tstepo_flux`'s nest is `k` outside, then `j`, then `i`, and each index carries
a different kind of dependence:

| index | what carries | what it costs a decomposition |
| --- | --- | --- |
| `k` | `fb(l,i,j)`, the flux through the bottom face, filled at level k for level k+1 | k stays sequential; the levels are separated by the barrier at the end of the `!$OMP DO` |
| `i` | `fw(l)`, the flux through the western face, filled by the previous cell in the row | nothing: it is contained within one `j` |
| `j` | `fs(l,i)`, the flux through the southern face, filled by the previous ROW | the one dependence that has to be removed |

It is removed by recomputation. The northern flux of row `j-1` is a pure
function of `u`, `ts1` and `k1` at that row, and **`ts1` is not written anywhere
in this routine** -- `tstepo` updates it after this returns -- so each row
computes its own southern face from the same expression, at `j-1`, that the
recycling filled `fs` from. That is exact rather than approximate, it costs one
extra meridional flux per cell, and it makes the store `fs(l,i) = fn(l)` dead,
which is deleted.

The two reductions are order-independent: `limps` counts flux-limited points and
`dmax` takes a maximum. With `SCHEDULE(STATIC)` beside them the result is
reproducible across thread counts and not merely correct at each, which is what
the acceptance test below actually demonstrates.

`velc` runs three worksharing loops where it had one loop and two tails, and
they have to stay separate: the periodic copy takes `u(1,imax,j,k)` into
`u(1,0,j,k)` and the vertical velocity reads `u(2,i,j-1,k)`, so both read a
NEIGHBOURING row of what the first loop wrote. `co` and `coshuffle` are column
algorithms and needed only their scratch made private.

## 4b. Two scalars had to be initialised inside the region, and that is the class

`diffv` in `tstepo_flux` and `icond` in `co` are assigned per cell only under a
condition, and elsewhere keep the value they had -- for `diffv` with `iediff` at
0 that is `diff(2)` from before the loop, which is never reassigned at all. A
private copy that inherited nothing would be whatever was on the stack. Both are
now set inside the parallel region.

That is the same failure the `-fno-automatic` demonstration is about, arriving
from the other side: the flag made locals wrongly SHARED, and a `private` clause
makes them wrongly UNINITIALISED. The `omppoison` arm below is what tests for
the second.

`velc`'s `dzu` needed the opposite treatment. It is declared `dzu(2,maxk)` in
`ocean.cmn`, so it is one array for the whole process, but it is per-column
scratch: `velc` is its only reader and writer anywhere in the tree and
`initialise_goldstein.F` only zeroes it. It is now a local. The COMMON member is
left in place because removing it changes the block's layout in every file that
includes the header.

## 4c. What the acceptance test actually ran

Every arm below reproduces the tree at `4038c12e` BIT-FOR-BIT, on all 28
variables of `eb_go_gs`'s GOLDSTEIN year-20 average and all 88 of
`eb_go_gs_ac_bg`'s BIOGEM three-dimensional fields.

| arm | what it isolates | threads |
| --- | --- | --- |
| `serial` | the source changes, compiled with the `!$OMP` sentinels as comments | 1 |
| `initpoison` | a local read before it is written | 1 |
| `omp` | the threads | 1, 2, 4, 8, 16 |
| `omppoison` | a private copy that needed an inherited value | 16 |
| `omp`, BIOGEM case | the threads under fourteen tracers and a raised stack | 1, 16 |

---

# 5. What it cost, measured

The cost case is the parallelism note's own: the shipped `worjh2` topography at
36 x 36 x 16, EMBM and GOLDSTEIN and the sea ice, 100 model years at `nyear`
100, `-mcmodel=medium`. Retired instructions COUNTED with `perf stat -e
instructions` rather than sampled, which on a single-threaded process is exact
and does not move with what else the host is doing.

**The instrument agrees with the sibling's.** The unedited tree at `4038c12e`
retires 302.68 G instructions over that run. The parallelism note's section 2c,
which SAMPLED the same case at a period of 2e7, reports 302.4 G. Two different
instruments, 0.09 per cent apart.

## 5a. Serial, one thread

| arm | source | storage class | G instructions | against `4038c12e` |
| --- | --- | --- | ---: | ---: |
| `4038c12e` | unedited | `-fno-automatic` | 302.68 | -- |
| `serialstatic` | this branch | `-fno-automatic` | 247.88 | -18.11% |
| `serial` | this branch | `-frecursive` | 246.64 | **-18.51%** |
| `omp`, 1 thread | this branch | `-frecursive -fopenmp` | 246.64 | -18.52% |

Three readings.

**The whole serial saving is 18.5 per cent**, and the ocean's tracer transport
allocates nothing.

**The OpenMP runtime is free at one thread.** 246.638 G against `serial`'s
246.640 G is four significant figures of agreement, so a threaded binary run on
one thread is not paying for the option.

**And the storage class is nearly free here, which it was not before the last
commit.** With the earlier form of the threading, where every row recomputed its
own southern face, `-fno-automatic` cost 288.99 G against `-frecursive`'s
247.34 -- 16.6 per cent -- because the recomputation loop writes `fs` and reads
`pec` and `ups`, none of which a compiler can hold in a register when they live
in static storage. With one recomputation per thread block the gap is 0.5 per
cent. The lesson is not about either flag: **an instruction count taken under
one storage class does not price a source change under the other.**

## 5b. Threaded, and the bound is not what binds

Taken on this host at a one-minute load average between 3.3 and 6.0, which
`scripts/machine.py` calls marginal at the top of that range. The internal
control says it does not matter here: the `serial` arm at load 16.7 and the
`omp` arm at one thread at load 6.0 differ by 0.6 per cent in the wall clock and
by 0.001 per cent in instructions, and every arm's `task-clock` is within 0.3
per cent of its wall clock, so the process never waited for a core.

| threads | wall, s | speedup | G instructions | instructions added |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 14.47 | 1.00 | 246.64 | -- |
| 2 | 11.84 | 1.22 | 248.44 | +0.7% |
| 4 | 11.09 | **1.30** | 251.87 | +2.1% |
| 8 | 11.78 | 1.23 | 258.96 | +5.0% |
| 16 | 13.27 | 1.09 | 275.69 | +11.8% |

**The speedup peaks at 1.30 on four threads and falls away above it**, against
the parallelism note's sixteen-thread Amdahl bound of 12.47. Amdahl's law is not
what stops it, and neither is correctness: every one of these arms reproduces
`4038c12e` bit-for-bit.

What stops it is **synchronisation granularity, and the instruction column is
the evidence**. The added instructions are the OpenMP runtime's own -- entering
a region, dividing a loop, and waiting at a barrier -- and they grow from 0.7
per cent at two threads to 11.8 at sixteen. At this grid `tstepo_flux` enters
one parallel region per ocean timestep and takes a barrier at each of 16 levels,
and the loop it divides is 36 rows. On sixteen threads that is two or three rows
of 36 columns between barriers, which is a few hundred cells of work per
synchronisation.

Two things follow and both are measurements rather than readings.

**A bigger grid moves this and a smaller one does not.** The work between two
barriers scales as the cell count while the barrier cost does not, so 72 x 72
puts four times the work between the same barriers and 144 x 144 sixteen times.
The recommended ocean grid is therefore also the grid at which this
decomposition starts to pay, and the shipped 36 x 36 is the worst case for it.

**And the routines still serial matter more than the bound suggests.** EMBM is
12.5 per cent of this case and is not threaded, which is defensible because the
offline recipe does not run it -- but that means this measurement understates
what an EMBM-free configuration would get, and there is no EMBM-free
configuration to measure yet.

The honest summary is that the threading is correct, reproducible, and worth
1.3 times at the shipped grid, and that what it is worth at the grid this
project would actually run has not been measured.
