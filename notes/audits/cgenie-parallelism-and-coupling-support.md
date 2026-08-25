# Can the EMIC ocean be parallelised, and what grid should it be coupled on

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean, wind stress, sea ice and tracer transport
name modelled quantities and numerical properties of a candidate component, not
observations of anything.

Measured 2026-08-25 on this machine. This is OCN-19's profile, the parallelism
half of OCN-20's precondition, and the support and regridding half of OCN-11 and
OCN-10. `analysis/cgenie_profile.py` is the driver and
`analysis/cgenie_profile.json` is what it wrote; the toolchain, the make
overrides and the host load are stamped in its provenance block rather than
repeated here.

It follows `notes/audits/cgenie-build-cost-and-grid-ceiling.md`, which prices a
configuration and says nothing about which routine spends it. Everything here
takes that document's costs as given and answers the question it left open.

**The vendored tree was not edited.** It was exported out of the repository with
`git archive` and built there, because a worktree's ignored build products are
symlinks into the shared checkout.

---

# 0. The decision rule, fixed before the profile was read

Written into `analysis/cgenie_profile.py` and committed before any profile
existed, because a threshold chosen after the run it judges is not a threshold.

A routine counts as **parallel over cells or tracers** when its loop nest runs
over the grid indices and any loop-carried dependence in it is a reuse of a face
flux the previous cell already computed -- recoverable by recomputing one row,
one column or one level per thread -- with no global reduction other than a sum
or a maximum. It counts as **serial** when it carries a dependence over the whole
domain that no such recomputation removes; the banded triangular solve in
`ubarsolv.f` is the case that motivates the definition. Anything that is neither
is **unclassified** and is counted against the parallel fraction, not for it.

The verdict on question one then follows from the SERIAL share at 36 x 36 x 16:

| serial share | verdict |
| --- | --- |
| below 5 per cent | parallelisation is possible and Amdahl is not what limits it at this grid |
| 5 to 20 per cent | possible, with a stated ceiling that must be quoted with every speedup claim |
| above 20 per cent | the algorithmic route of OCN-20 comes before the threading route of OCN-19 |

The bound at one grid is not the bound at another, because the serial part grows
faster with resolution than the parallel part does. Both grids are therefore
reported and the trend between them is part of the answer.

---

# 1. What is serial today, re-verified rather than inherited

`notes/audits/ocean-and-marine-biosphere.md` section 9h and the cost note's
section 1a both say the component is serial. Re-checked here over the whole of
`vendor/cgenie`, because the claim is load-bearing for everything below:

- **No live OpenMP sentinel exists anywhere in the tree.** Every `$omp` string in
  it sits in exactly three files -- `genie-goldstein/src/fortran/tstepo.F`,
  `genie-biogem/src/fortran/biogem.f90` and `genie-main/genie.F` -- and every one
  is commented out in column one, so none is a directive.
- **Four `use omp_lib` statements are live** and no call stands against any of
  them: `biogem.f90`, `tstepo.F`, `genie.F` and `genie_loop_wrappers.f90`. The
  earlier count of three missed the last.
- **`-fopenmp` is commented out** in `genie-main/makefile.arc`'s gfortran block,
  and enabling it would parallelise nothing.
- **MPI exists in one directory**, `genie-plasim/src/fortran`, which is the
  vendored Planet Simulator this project does not use.

So the question is not what to re-enable. It is what a decomposition would have
to be, and what it would buy.

---

# 2. The instrument, and two ways it lies

Both are recorded before the numbers, because both were found by checking the
instrument against the size of the effect and both would have produced
ordinary-looking answers.

## 2a. The kernel throttles the sampler and does not say so

`perf record -e instructions -c <period>` takes one sample per `period` retired
instructions, which is what makes the shares a property of the binary and its
input rather than of what else this contended host is doing. But
`kernel.perf_event_max_sample_rate` is **3000 samples per second** here, and
above it the kernel stops delivering samples for the rest of the window and
resumes afterwards. `perf report` still says `Total Lost Samples: 0`, because no
record was lost. What is lost is the assumption the profile rests on -- that
every retired instruction is equally likely to be sampled.

Measured: at a period of 2e7 the shipped grid asks for about 750 samples per
second and delivers exactly the number the period demands, 1409 for 28.18e9
instructions. At 2e6 the same run asks for about 7500 and delivered **7133 of
the 14090 the period demands**, and the per-symbol shares moved by up to a factor
of 1.6 -- eight standard errors on the counts involved. So the tighter sampling
was not a better measurement of the same thing. It was a measurement of a
different, throttled process.

**More samples come from a longer run at a period the kernel will honour, never
from a shorter period.** `analysis/cgenie_profile.py:sample_rate_ok` now computes
the rate a period asks for, compares it with the cap, and stamps a `throttled`
line in the artifact when it is over. Every number below is from a period of 2e7,
under the cap, over 100 model years.

## 2b. The profile is a function of the model state, and a short run is a transient

The same binary and the same input at 36 x 36 x 16, profiled over the first 10
model years, the first 100 and the first 1000:

| | 10 years | 100 years | 1000 years |
| --- | ---: | ---: | ---: |
| `tstepo_flux` | 44.14% | 53.28% | 52.14% |
| `co` | 6.25% | 3.31% | 2.61% |
| `eos` | 2.91% | 1.41% | 1.19% |
| `ubarsolv` | 2.06% | 1.83% | 1.73% |
| **parallel share** | 87.57% | **88.85%** | **89.06%** |
| **serial share** | 2.06% | **1.84%** | **1.75%** |

The moves between 10 and 100 years are what a cold start predicts. The run begins
from a uniform 10 degree ocean, so convective adjustment and the
equation-of-state calls inside it do their heaviest work at the start and fall
away as the model stratifies, and `tstepo_flux`'s isoneutral diffusion branch is
entered only where `dzrho < -1e-12`, so it turns ON as stratification develops. A
profile of the first ten years of a spin-up is a profile of a transient.

**Between 100 and 1000 years the decomposition has settled**, which is the
statement that matters: the parallel share moves by 0.2 points and the serial
share by 0.09, against 1.3 and 0.2 between 10 and 100. The sixteen-thread bound
with the heap traffic removed is 12.47 at 100 years and 12.51 at 1000. Individual
routines are still drifting -- `co` is still falling -- so this is not a claim
that the model has equilibrated, only that what the verdict rests on has stopped
moving.

It is recorded because the difference between the shortest run and the longest at
one grid is LARGER than the difference between the two grids, so a reader
comparing grids across run lengths would be comparing states. Every number in
section 2c is at 100 model years for both grids.

## 2c. Where the instructions actually go

100 model years, `nyear = 100`, GOLDSTEIN's `debug_loop` off, both grids compiled
`-mcmodel=medium` so they differ only in the grid. 15,121 samples at the shipped
grid and 68,236 at the doubled one; the sampled totals, 302.4 G and 1364.7 G
instructions, agree with the cost note's independently fitted per-year slopes to
1 and 1 per cent, which is the check that the sampler saw the whole run.

| routine | what it is | 36 x 36 x 16 | 72 x 72 x 16 | divides? |
| --- | --- | ---: | ---: | --- |
| `tstepo_flux` | ocean tracer advection and diffusion | 53.28% | 46.60% | yes |
| `tstipa` | EMBM's implicit temperature and humidity step | 8.48% | 14.41% | yes |
| `velc` | ocean velocity from the density field | 6.20% | 5.55% | yes |
| `co` | convective adjustment | 3.31% | 2.82% | yes |
| `goldstein` | the ocean driver's own per-cell loops | 3.06% | 3.02% | yes |
| `surflux` | EMBM's surface fluxes and the sea-ice surface solve | 2.49% | 2.82% | yes |
| `eosd` | vertical density gradient, called per wet cell | 1.96% | 3.92% | leaf |
| `eos` | density, called per wet cell | 1.41% | 1.75% | leaf |
| **`ubarsolv`** | **the barotropic streamfunction solve** | **1.83%** | **2.57%** | **no** |
| `embm` | the EMBM driver's own per-cell loops | 1.55% | 2.66% | yes |
| `jbar` | pressure torque, integrated down each column | 1.50% | 0.97% | yes |
| `tstepo` | the tracer step's wrapper, and its dead array copies | 0.56% | 0.38% | yes |
| `tstepsic` | sea-ice transport | 0.56% | 0.53% | yes |
| libm `pow`, `log`, `exp` | leaves of the per-cell nests | 3.91% | 3.44% | leaf |
| libc `malloc`, `free`, `memmove` | see 2d | 9.26% | 7.87% | see 2d |

Three readings, and none of them is what the audit expected before it was
measured.

**The tracer transport is the run.** `tstepo_flux` alone is half the
instructions at both grids, and everything around it is a per-cell or per-column
loop. Section 9h of the ocean audit predicted that the tracer loops are where the
work worth dividing sits; this is the measurement, and it is more concentrated
than the prediction.

**The barotropic solve is small and it is the only thing that is genuinely
serial.** 1.83 per cent at the shipped grid. Section 4a of the cost note reached
the same conclusion without a profiler, from two shipped grids with nearly equal
cell counts and a 1.69-fold difference in barotropic work, and returned 0.956
against a decision rule fixed beforehand. Two independent instruments, the same
answer.

**EMBM is an eighth of the work at the shipped grid and a fifth at the doubled
one**, 12.52 per cent against 19.89. That is a component the offline architecture
does not run, and section 4 is where it is spent.

## 2d. A tenth of the run is heap traffic for a two-line defect

The libc row above is not the model doing anything. `tstepo_flux.F:82` reads

    call eosd(ec,ts1(1,i,j,k:k+1),ts1(2,i,j,k:k+1),zw(k),
   1     rdza(k),ieos,dzrho,tec)

and `eosd` declares its arguments as `real t(2)`. Those two actual arguments are
STRIDED sections of a four-dimensional array, so gfortran packs each into a
contiguous temporary before the call and frees it after. Disassembling the built
executable confirms it at the instruction level: the basic block around the one
`call eosd_` in `tstepo_flux_` contains two `call malloc@plt` and three
`call free@plt`.

That happens once per wet cell per level per ocean timestep -- of order 11,500
times a step at 36 x 36 x 16, of order 11.5 million times over ten model years --
and it accounts for **9.26 per cent of the instructions at the shipped grid and
7.87 per cent at the doubled one**.

It is not a decomposition question. It is work that should not exist, it is
removable by giving `eosd` two scalars or copying into a two-element local, and
the change is correctness-preserving, so `genie-knowngood/` and `nccompare`'s
relative tolerance in units in the last place are already the right test for it.
**It is a larger, cheaper saving than anything threading offers per unit of
effort**, and it should be taken before any parallelisation is designed, because
it also removes an allocator call from the middle of what would become the
innermost parallel loop.

A second and much smaller instance of the same class sits beside it. `tstepo.F`
fills `ts_t1`, `ts1_t1`, `rho_t1`, `ts_t2`, `ts1_t2` and `rho_t2` on every ocean
timestep and nothing reads them: their only consumers are the commented-out
OpenMP sections. Section 9h of the ocean audit found that statically and could
not say what it cost. It costs about 0.56 per cent of the run at 36 x 36 x 16 and
0.38 at 72 x 72 x 16, which is the whole of `tstepo`'s own share.

## 2e. The split, and what Amdahl's law then bounds

Applying section 0's classification, which was fixed before any of this was seen:

| | 36 x 36 x 16 | 72 x 72 x 16 |
| --- | ---: | ---: |
| parallel over cells, columns or tracers | 88.85% | 89.34% |
| serial (`ubarsolv`) | 1.84% | 2.60% |
| heap traffic (section 2d) | 9.26% | 7.87% |
| unclassified | 0.14% | 0.03% |

**The serial share at the shipped grid is 1.84 per cent, which is below section
0's 5 per cent line, so the verdict is that parallelisation is possible and
Amdahl's law is not what limits it at this grid.**

Two bounds follow, and the difference between them is the section 2d defect
rather than any property of the decomposition:

| threads | heap traffic left in place | heap traffic removed first |
| ---: | ---: | ---: |
| 2 | 1.80 | 1.96 |
| 4 | 3.00 | 3.79 |
| 8 | 4.49 | 7.07 |
| 16 | 5.99 | 12.47 |
| 32 | 7.18 | 20.18 |

The right-hand column is the one a budget should be built on, because the heap
traffic is a defect with a two-line fix rather than a part of the calculation.
The left-hand column is what a thread team would actually see if someone
parallelised the loops and left the allocator calls inside them, and it is
recorded so that outcome is recognisable if it happens. Neither column is a
prediction of measured speedup: Amdahl's law is an upper bound and says nothing
about memory bandwidth, synchronisation or load imbalance across a land-sea mask.

**The bound at 36 x 36 is not the bound at 144 x 144, and the trend is measured
rather than assumed.** Two corrections have to be made before the two grids can
be compared, and both are stated rather than folded in.

The two profiles ran at different `ndta`, 5 and 10, because that is what EMBM's
stability demanded at each grid, and `tstipa` and the `embm` driver run once per
EMBM step where everything else runs once per ocean step. So EMBM is
double-weighted at 72 x 72 and every ocean share there is diluted. Rescaling the
doubled grid to `ndta = 5` puts `ubarsolv` at 2.81 per cent rather than 2.57,
and the share therefore grows by **1.53 per doubling** rather than the 1.41 the
raw numbers give.

In absolute terms `ubarsolv` went from 5.53e5 to 3.51e6 instructions per ocean
timestep, a factor of 6.34 where its loop bounds predict 7.79. Everything else
scaled as cell count predicts -- the rescaled total ratio is 4.13 against 4.00 --
so the whole of the deviation is in `ubarsolv` itself, and it is what a short
inner loop costs: the shipped grid's band is 37 wide against 73, so it amortises
its per-element overhead over half as much work. That overhead does not grow, so
the asymptotic behaviour is the r^3 the loop bounds give, and the share ratio
approaches 2.0 rather than staying at 1.53.

Both readings, applied to the 2.81 per cent the doubled grid corrects to:

| grid | serial share, at the measured 1.53 | at the asymptotic 2.0 | sixteen-thread bound, heap removed |
| --- | ---: | ---: | --- |
| 72 x 72 x 16 | 2.8% | 2.8% | 11.3 |
| 144 x 144 x 16 | 4.3% | 5.6% | 8.7 to 9.7 |
| 288 x 288 x 16 | 6.6% | 11.2% | 6.0 to 8.0 |

**The serial fraction grows, and it does not become the binding constraint
anywhere below muffingen's declared 72 ceiling or at twice it.** What binds first
is the r^4 in the timestep, which threads do not touch, and that is OCN-20's
territory rather than OCN-19's.

---

# 3. What the parallelisation would cost to do

The loop structure is not the obstacle. A build flag is.

## 3a. `-fno-automatic` makes every procedure-body local SHARED, silently

`makefile.arc`'s gfortran SHIP block compiles with
`-O2 -O3 -funroll-loops -msse -fno-automatic`. That last flag moves every
procedure-body local into static storage, and a variable in static storage is one
variable for the whole process. Inside an OpenMP region it is therefore shared
between the threads unless it is named in a `private` clause, and the compiler
will not tell you which ones you missed.

Demonstrated rather than argued, with a loop whose right answer is known:

    subroutine work(n, out)
      real :: scratch(512)
      !$omp parallel do private(j)          ! scratch NOT named
      do i = 1, n
         do j = 1, 512
            scratch(j) = real(i) + real(j)
         end do
         out(i) = scratch(1) + scratch(512) ! must be 2*i + 513
      end do

Built `gfortran -O2 -fno-automatic -fopenmp`, run over 200,000 iterations:

| threads | private clause | cells wrong |
| --- | --- | ---: |
| 1 | none | 0 of 200,000 |
| 8 | none | 177,007 of 200,000 |
| 8 | `private(scratch)` | 0 of 200,000 |

Three things follow, and each is a property of the work rather than of this
example. The wrong answer is a real answer -- finite, plausible, and dependent on
the thread count, which is the failure class `docs/src/practice/failure-modes.md`
already tracks. `private` on such a variable is ACCEPTED by gfortran rather than
rejected, so there is no compile-time gate to lean on. And the whole combination
draws exactly one diagnostic for the entire translation unit,
`Flag '-fno-automatic' overwrites '-frecursive' implied by '-fopenmp'`, which
says nothing about which variables are affected.

**So the cost of the threading route is an enumeration, not a directive.** Every
procedure-body local in every routine reachable from a parallel region has to be
classified and named. That is the same class as CLIM-51's 103 implicitly SAVEd
declarations in `plasim/src`, arriving here through a compiler flag instead of
through the language default.

## 3b. The flag looks removable, and the test for it already ships

The cheaper route is to stop compiling that way. `-frecursive` gives an unlimited
`-fmax-stack-var-size`, which puts the same locals on the stack where each thread
has its own. It costs stack rather than correctness, and `makefile.arc` already
carries the commented line `F90FLAGS += -frecursive` beside gfortran's own
warning about exactly this.

Whether the tree DEPENDS on the implicit SAVE that `-fno-automatic` supplies is
checkable, and the evidence says it very likely does not. Across
`genie-goldstein`, `genie-embm` and `genie-goldsteinseaice` there are ten `DATA`
statements and eight `save` statements, every one of them EXPLICIT and every one
in a driver, a diagnostic or a netCDF writer rather than in a per-cell kernel:
`goldstein.F:133,137` and `gold_seaice.F:65,69` (the first-call flags and the
initial energy and water inventories), `diag2.f:22`, `diagosc.F:56-58`,
`diagend.F:48` and `diagend_embm.F:47` (diagnostic accumulators), and the month
name tables in the three `netcdf*.F` writers. A routine that declares its
persistence keeps it under either flag.

**And the acceptance test needs no invention.** `genie-knowngood/` ships
per-component netCDF output for four configurations and `genie-main/src/c/compare.cpp`
builds `nccompare`, which takes a relative tolerance in units in the last place.
Dropping `-fno-automatic` is a change that should move nothing at all, so its
right answer is bit-for-bit or within a declared ULP count -- a test that can
fail, in the sense `CLAUDE.md` requires, rather than a comparison that can only
differ. That check is the first step of any parallelisation work here and it
needs no threads.

## 3c. What a decomposition would look like, read off the loop structure

Stated because the profile below only says where the time is, not whether it can
be divided. Each of these is a reading of the source, not a measurement.

- **`tstepo_flux.F`, the ocean tracer transport.** A `k, j, i, l` nest with a
  flux-recycling optimisation: `fw(l) = fe(l)`, `fs(l,i) = fn(l)` and
  `fb(l,i,j) = fa(l)` carry the face flux forward so it is computed once per
  face. That is a loop-carried dependence in all three space indices, and it is
  the recoverable kind -- each of those is a pure function of `ts1` and `u` at
  the face, so a thread taking a block of rows recomputes one row of meridional
  flux at its southern edge and owns the rest. The code already does exactly this
  at `i = 1` for the western doorway. The two genuine reductions, `limps` and
  `dmax`, are a count and a maximum.
- **`co.F` convection, `krausturner.F`, `jbar.f`, `velc.f`.** Column algorithms:
  a vertical recurrence inside each `(i,j)`, no horizontal dependence. Parallel
  over columns, which is thousands of them at any grid considered here.
- **`eos.f` is pure.** It declares no locals and is called from the innermost
  tracer loop, so it needs nothing.
- **EMBM's `tstipa.f` is a Jacobi iteration**, not a line solve: `tq2` holds the
  whole previous iterate and the update reads only `tq2` and `tq1`. Every one of
  its `nii` sweeps is data-parallel over cells with a barrier between them. It is
  what runs: `embm.F:195-199` selects it under `dimpa`, which this build defines,
  and `tstepa.f` -- the explicit alternative, and the same flux-recycling shape as
  `tstepo_flux.F` in two dimensions -- takes no samples at either grid. The ocean
  is the other way round: `goldstein.F` selects `tstipo` under `dimpo`, which is
  NOT defined, so `tstepo` runs.
- **`surflux.F` is one `(i,j)` loop** running from line 683 to about 1600,
  including the sea-ice surface Newton iteration and the land column. Parallel
  over cells, and the routine with the most locals to classify by a wide margin.
- **`ubarsolv.f` is the exception.** Its forward elimination and back
  substitution over the banded factors are sequential along the band by
  construction, and no recomputation removes it. This is the routine the decision
  rule in section 0 is about.

## 3d. The working set is cache-resident at every grid considered

`CLAUDE.md`'s thread-team target is 32 MB on one die, and the two figures that
bear on it are three orders of magnitude apart.

The executables themselves are enormous: `size -A` gives 1.39 GB of static
storage at 36 x 36 x 16 and 4.42 GB at 72 x 72 x 16, of which 4.38 GB at the
doubled grid is `.lbss` -- the large-data section `-mcmodel=medium` exists to
address. Every field cGENIE holds is in a named COMMON sized from the grid macros
and `-fno-automatic` makes the locals static too, so all of it is allocated
whether or not it is touched.

The arrays a GOLDSTEIN timestep actually touches are a different quantity
entirely: `ts`, `ts1`, `u` and `rho` come to roughly 2 MB at 36 x 36 x 16 in
double precision, roughly 8 MB at 72 x 72 x 16, and roughly 32 MB at
144 x 144 x 16. That is arithmetic from the declared shapes rather than a
measurement of a cache miss rate, and it says a thread team on the tracer
transport would be dividing an arithmetic-bound problem rather than a
bandwidth-bound one, up to about the grid at which the horizontal refinement
stops being affordable anyway.

The gap between the two is `embm.cmn`'s seasonal arrays and `ocean.cmn`'s
`maxnyr`-dimensioned storage, which the timestep does not read.

## 3e. How big the enumeration is

The seventeen files a parallel region over the ocean, the sea ice and the surface
fluxes would contain are `tstepo`, `tstepo_flux`, `co`, `krausturner`, `velc`,
`jbar`, `wind`, `eos`, `ediff`, `get_hosing`, `tstepa`, `tstipa`, `surflux`,
`radfor`, `ocean_alb`, `tstepsic` and `tstipsic`. Counting the distinct names on
their type-declaration lines gives **390 declared names, 81 of them arrays**.

That is an UPPER bound and is meant to size the work rather than to be the list.
It does not separate dummy arguments, which need no private clause, from genuine
locals, and it counts loop indices. Two things about its shape matter more than
its total. `surflux.F` alone carries 141 of the 390 and 20 of the 81, over a
third of the work in one routine, which is the routine with the sea-ice Newton
iteration and the land column inside it. And the ocean's own tracer path --
`tstepo_flux`, `co`, `eos`, `ediff`, `velc`, `jbar` -- comes to 90 names and 30
arrays, which is a day's careful reading rather than a project.

Set beside section 3b, that is the argument for trying `-frecursive` first: if
dropping `-fno-automatic` reproduces `genie-knowngood/` within a declared ULP
count, all 390 become automatic and thread-private by the language's own default,
and the enumeration does not have to happen at all.

---

# 4. The coupling resolution

## 4a. The constraint that was binding is not binding any more

`notes/audits/ocean-and-marine-biosphere.md` section 6 chose offline coupling on
one property, that the ocean's cost must not climb the atmosphere's ladder, and
the cost note measured that the ocean's price is a function of its own grid
alone. What was left was that a finer ocean costs r^4 on one serial core, so the
grid was an affordability question with connectivity on the other side of it.

Section 2 removes most of that. Three reductions compose, and each is measured or
read off a measurement rather than assumed:

**EMBM is 12.5 to 19.9 per cent of the instructions, and the offline
architecture does not run it.** The exchange OCN-10 selects has ExoPlaSim
computing the fluxes, so EMBM's prognostic temperature and humidity step,
`tstipa` plus the `embm` driver's own loops, has nothing to do.
`surflux`'s 2.8 per cent stays, because something equivalent -- the
`surflux_goldstein_seaice` path the cost note's section 6 identifies -- still has
to turn a supplied climatology into ocean fluxes.

**EMBM also sets the timestep, and a configuration that does not integrate it
does not inherit its limit.** The cost note's section 3 proves this with pairs at
matched `dtatm` and ocean timesteps a factor of two apart failing in the same
cell, and its section 3c puts the ocean's OWN limit at `nyear` about 13 at
36 x 36 x 16 and about 52 at 72 x 72 x 16, against the 100 both profiles above
used. So the ocean could take between two and eight times the timestep EMBM
permits, and everything left after EMBM is removed runs once per ocean step.

**A tenth of what remains is the section 2d heap defect**, which is a two-line
fix with an acceptance test that already ships.

Composed at 72 x 72 x 16, from the measured 13.65 G instructions per model year:
remove EMBM's 17.1 per cent, take the ocean's own timestep instead of EMBM's, and
remove the heap traffic, and the physics costs about **5.3 G instructions per
model year**. A 20,000-year spin-up is then about 1.1e14 instructions, roughly
two hours on one core of this machine, against about one hour for the shipped
36 x 36 x 16 grid as it stands today. Threading is a further bound of about 11 at
sixteen threads on top of that.

**This is a composition of measured pieces under one assumption that has not been
demonstrated**, and the assumption is named in section 6: no EMBM-free
configuration has been built or run here, and `genie.F` offers only two
surface-flux paths, each gated on an atmosphere module being in the recipe. Read
it as the budget the resolution decision should be argued against, not as a
measurement of a configuration that exists.

## 4b. So the grid is chosen from connectivity, and the recommendation is 72 x 72 x 16

With the cost constraint no longer binding below muffingen's ceiling, only one
constraint is left, and it is the one section 8b of the ocean audit and section
5c of the cost note both arrive at: straits, sills and partial coasts. At 10
degrees of longitude they are not represented at all. Section 5d above adds that
the coastal mask disagreement between the atmosphere and the ocean is the part of
the regridding contract that resolution genuinely improves, where conservation is
settled by construction.

**Recommendation: 72 x 72 x 16 with `igrid = 0`**, the equal-area grid.

- It is the top of muffingen's declared `[1-72]` range for `par_max_i` and
  `par_max_j`, so it is the finest grid the generator claims, and the cost note
  demonstrated that the MODEL runs there.
- 5.0 degrees of longitude is a factor of two better than the shipped grid on the
  only axis that is still constraining.
- `igrid = 0` rather than 1 or 2 for a reason the model states itself:
  `tstipa.f:37-46` raises EMBM's implicit iteration count from 4 to 16 on the
  constant-latitude grid, and `initialise_goldstein.F:496` makes the equal-area
  grid's cell area exactly constant, which is what section 5b's factorised
  conservative remap wants.
- The ocean's grid does not depend on the atmosphere's rung, so this choice
  survives the ladder moving.

**36 x 36 x 16 is rejected**, and not on cost -- it is the cheapest thing here.
It is rejected because it is the resolution at which the connectivity this
project needs is provably absent, and the affordability argument that was the
only reason to accept that is no longer true.

**64 x 64 x 16 is the tiebreak candidate, and there is a decision procedure
rather than a preference.** 64 ocean longitudes divide the longitude count of
every rung on the ladder except T31's 96: T21's 64, T42's 128, T63's 192, T85's
256, T106's 320, T127's 384 and T170's 512 are all whole multiples. That makes
the longitude half of the conservative remap a whole-number block sum with no
partially covered columns anywhere, which removes a class of bookkeeping rather
than a class of error -- section 5b's operator is exact either way. It costs 79
per cent of 72 x 72's cells and gives 5.625 degrees instead of 5.0. So:

> OCN-11's bathymetry contract reports, for each candidate ocean grid, the
> straits and sills the Orogen mesh says exist and how many of them the grid
> resolves. Take the COARSEST candidate that resolves the connections the mesh
> carries; among candidates that pass, prefer one whose longitude count divides
> the accepted atmosphere rung's.

That is a criterion fixed before the inventory is seen, which is the point.

## 4c. What has to be measured before this is committed to

Three things, none of which this document does, and each of which can fail.

1. **A stability sweep at the chosen grid**, by the cost note's method: GOLDSTEIN's
   own `Cn` from `diag.f` with GOLDSTEIN's `debug_loop` on and EMBM's off, over
   `nyear`, with the criterion `Cn < 1` at every diagnostic step. A grid that
   RUNS is not a grid that is STABLE, and past the limit the model completes and
   reports success. The ocean's limit also inherits OCN-17's declared `scf`
   bracket of 1 to 3 linearly, so what comes out is a bracket rather than a
   number.
2. **muffingen actually producing `.k1`, `.paths` and `.psiles` at that grid**
   from the Orogen bathymetry, and REPORTING the island count it resolves --
   `GOLDSTEINMAXISLES` is a compile-time bound. Whether muffingen runs outside
   MATLAB at all is world-crky and is not settled.
3. **Whether an EMBM-free configuration exists.** Section 4a's budget rests on
   it, and nothing here or in the cost note has built one.

---

# 5. The regridding contract

Two directions, and they are not one problem. Upscaling the atmosphere's fluxes
onto the ocean is a conservative reduction with an exact answer. Downscaling the
ocean's surface state onto the atmosphere invents structure the ocean never
resolved. They cost different things and the contract states them separately.

## 5a. The precedent chose matching, and the reason does not transfer intact

Holden et al. (2016) is the only existing PlaSim-to-GOLDSTEIN coupling and it
restricted itself to T21 against a matched 64 x 32 GOLDSTEIN grid, in its own
words "in order to avoid the need for interpolation", with 10 atmospheric levels
against 32 ocean levels. That coupling was SYNCHRONOUS: fluxes crossed every
coupling step, so an interpolation error entered the prognostic state, was
integrated, and could feed back through the atmosphere's own response to the sea
surface temperature it had just perturbed.

This project's coupling is offline and passes climatologies, so the feedback path
is gone. What is NOT gone, and is worse in one specific way, is the residual: a
climatology is a fixed field, so a non-conservative regrid of it is a permanent,
spatially structured spurious source term in an ocean that is then integrated for
of order 10,000 model years. A synchronous coupling's interpolation error is
transient and partly self-cancelling; an offline one's is a bias with unlimited
time to express itself.

**So the precedent's requirement changes shape rather than relaxing.** What the
offline architecture needs is not "avoid interpolation" but "the crossing
conserves exactly and the residual is reported", and that is a much more
tractable requirement than Holden's, because it is a property of a weight matrix
computed once rather than of an error budget accumulated every coupling step.
Grid matching is one way to get it and not the only one, and section 59c of
`notes/external-model-survey.md` already records that declining to match means
writing the regridder the offline architecture needs anyway.

## 5b. Upscaling the fluxes: exact, and cheaper than it looks

Both grids are separable in longitude and latitude. ExoPlaSim's is a Gaussian
grid: uniform in longitude, Gaussian latitudes. GOLDSTEIN's shipped `igrid = 0`
is uniform in longitude and uniform in the SINE of latitude, and
`initialise_goldstein.F:496` sets `asurf(j) = rsc*rsc*ds(j)*dphi` -- so every
ocean cell has the same area exactly, by construction.

Two consequences follow that make this direction cheap.

**The conservative weight matrix factorises.** Because neither grid's cell
boundaries depend on the other coordinate, the overlap area of an atmosphere cell
and an ocean cell is the product of a longitude overlap and a sine-of-latitude
overlap. The operator is therefore two one-dimensional sparse matrices rather
than a polygon intersection, it is exact in floating point, and the planetary
radius cancels out of it because every weight is a ratio of areas.

**The nonlinearity has already been resolved on the fine grid.** A surface heat
flux is not a linear function of the surface state, so the average of the flux is
not the flux of the average. In the architecture OCN-10 selects -- ExoPlaSim
computes the fluxes and there is no EMBM -- what crosses the boundary IS the
flux, already evaluated at 1.4 degrees or wherever the accepted rung sits, and
averaging it is a conservative linear operation. In the eliminated offline-regrid
architecture EMBM would have recomputed the fluxes from coarse state, putting the
nonlinearity back at the ocean's resolution. **That is an argument for the chosen
architecture that the record did not carry**, and it sits beside the tuning
argument in section 59b rather than replacing it.

What upscaling costs in fidelity is therefore exactly one thing: the subgrid
variance of the flux inside an ocean cell is discarded. `lib/gridding.py` already
owns the vocabulary for saying what that means field by field -- extensive,
intensive, categorical, moments, expectation -- and its `cell_moments` and
`cell_expectation` are the operators for the cases where a mean is not enough.
Which of them a field takes is a property of the field and belongs in the
contract, not in the call.

## 5c. Five rules the upscaling contract has to state

Read off `references/esmf/` (survey section 55), which is the only conservative
regridder in this project's reading with a stated conservation identity and a
test suite against it.

1. **Normalisation is a property of the field.** A destination value formed as
   intersection area over DESTINATION area is what makes an extensive total
   close; a destination value formed as intersection area over COVERED area is
   what makes a density or a rate come out right. A field of ones must remap to
   ones under the second and does not under the first. Net heat flux, freshwater
   flux and wind stress are rates per unit area and take the second; an integrated
   total takes the first.
2. **The coverage fraction travels with the result.** A fraction that is
   discarded is a conservation identity that can no longer be evaluated.
3. **An unmapped destination cell is an error, not a backfill.** `lib/gridding.py`
   currently backfills a grid cell with no mesh region from the nearest region
   centre and returns the count; the precedent refuses by default and makes the
   fallback the opt-in. For the ocean crossing the refusal is the right default,
   because an ocean cell with no atmosphere over it is a mask disagreement rather
   than a sparse-coverage artifact.
4. **The closure residual is registered in advance and stamped in the artifact.**
   `sum(F_src * A_src * frac_src)` against `sum(F_dst * A_dst)` has a right answer
   and the precedent holds itself to about 1e-9 relative on analytic fields.
   `lib/gridding.py:transfer_ledger` already emits the shape of this; what it has
   never been checked against is a field whose integral is known.
5. **The source cell edges are the ones the model's own quadrature implies, not
   midpoints between centres.** The ocean side is unambiguous: `igrid = 0` puts
   the row boundaries at exact equally spaced values of the sine of latitude and
   `asurf(j) = rsc*rsc*ds(j)*dphi` is the area the model integrates over. A
   Gaussian grid has no cell boundaries at all -- it has nodes and quadrature
   weights -- so the edges have to be CONSTRUCTED, and the only construction that
   makes the remap conservative against the atmosphere's own budget is the one
   whose cell areas reproduce the Gaussian weights. Taking midpoints between
   Gaussian latitudes instead gives areas that differ from the weights, so the
   crossing would then conserve against a grid the model does not use, and it
   would do it quietly. This is where `CLAUDE.md` rule 3 bites in its area form:
   the two sides must share one coordinate source, and the areas are part of the
   coordinate.

## 5d. The mask disagreement is the part that does not have a clean answer

The ocean's wet mask comes from the Orogen mesh through OCN-11, and the
atmosphere's comes from the same mesh at a different resolution. They will not
agree at the coast, and they cannot: a cell that is 40 per cent land at 1.4
degrees is inside a 5.6-degree ocean cell that is either wet or dry.

So there are atmosphere cells carrying an ocean flux that lie over ocean-grid
land, and ocean cells with no atmosphere-ocean area over them at all. Every joule
and every kilogram in the first class has to go somewhere named, or the global
heat and freshwater budgets do not close and the residual is not noise but a
coastline-shaped field. This is the same problem `references/climber-x/`'s
`coast_cells.f90` solves for river discharge, cited in OCN-11's notes, and the
contract's answer has to cover heat and salt as well as water.

**This is the part of the crossing that a finer ocean grid genuinely improves**,
and it is the argument for resolution that survives when the conservation
argument is settled by construction.

## 5e. Downscaling the ocean's surface state: what it invents

Three fields cross the other way, and none of them is a flux: the sea surface
temperature the next climate run is forced with, the ice state world-pt8 wants,
and the surface velocity `goldstein.F:780-781` exports as `ustar_ocn` and
`vstar_ocn`.

These are intensive state variables, so the crossing is an interpolation and
there is nothing to conserve. What it costs in fidelity is that the atmosphere
receives a field with no structure below the ocean's cell size, while its own
dynamics respond to gradients at its own. At an ocean with 5.6-degree longitudes
against a T85 atmosphere at 1.4, every sea surface temperature gradient the
atmosphere sees is a 5.6-degree gradient smeared over four cells: no fronts, no
western boundary currents, and no straits at all, because `notes/audits/ocean-and-marine-biosphere.md`
section 8b already records that straits and sills are not represented at these
resolutions in the first place.

Two rules follow, and both are about honesty rather than accuracy.

**Interpolate the state, never a flux, in this direction.** A bilinear or
higher-order interpolation of an intensive field is the right operator and a
conservative one is not, because there is no integral to preserve. The
temptation to route a heat transport this way -- an extensive quantity dressed as
a state -- is what the two words in `lib/gridding.py` exist to separate.

**The support mismatch is declared, not smoothed.** An adequacy claim about the
coupled system is a claim about a field carried across a factor of four in
longitude at T85 and a factor of eight at T170, and the honest form of it names
that factor. Section 5c of the cost note reaches the same place from the price
side: a finer ocean is affordable and a matched one is not, so the mismatch is a
fact to state rather than one to avoid.

---

# 6. What this did NOT establish

- **Nothing here is an implementation and nothing here is an adoption.** No
  OpenMP directive was written, no file under `vendor/cgenie` was edited, and the
  tree still has no consumer, no pipeline row and no step. The one compiler flag
  that moves is the code model, applied to both grids equally so it cancels out
  of the comparison.
- **Amdahl's law is an upper bound and no threaded build exists to measure
  against it.** Load imbalance across a land-sea mask, memory bandwidth, false
  sharing on the COMMON blocks, and the cost of the halo recomputation section 3c
  proposes are all unpriced. The bound says what is not achievable; it does not
  say what is.
- **The EMBM-free configuration section 4a's budget rests on has not been built
  or run.** `genie.F` offers exactly two surface-flux paths for GOLDSTEIN and each
  is gated on an atmosphere module being in the recipe; the second is the right
  shape and that is a reading of the recipe, not a demonstration. Until it exists,
  section 4a is arithmetic over measured pieces rather than a measurement.
- **The ocean's own stability limit has not been measured without EMBM.** The
  cost note's sweep had EMBM live throughout, and its statement that the ocean
  could carry two to three times the timestep EMBM permits was read off
  GOLDSTEIN's `Cn` in configurations EMBM was still integrating. The limit also
  inherits OCN-17's declared `scf` bracket of 1 to 3 linearly, so it is a bracket.
- **The profile is a function of the model state and neither run reaches
  equilibrium.** Section 2b shows the shares moving between 10 and 100 model
  years by more than they move between the two grids. A 20,000-year spin-up
  spends nearly all of its time in a state that no run here samples, and the
  direction of the drift -- the parallel share rising, the serial share falling --
  is what makes the verdict safe rather than what makes it precise.
- **The 72 x 72 probe's circulation means nothing.** Its pair of advective
  wind-speed fields were made by replicating the 36 x 36 fields into 2 x 2
  blocks, so the instruction counts and the loop structure are a 72 x 72 x 16
  model's and the ocean state is not evidence about anything.
- **Every number is at Earth's parameters.** `rsc` and `const_rEarth` are
  hardcoded, the topographies are Earth's, and
  `notes/audits/ocean-tier-implicit-earth.md` findings 1 and 4 stand exactly as
  they were.
- **The regridding contract in section 5 is a specification, not an
  implementation.** No weight matrix has been built, no conservation identity has
  been evaluated, and `lib/gridding.py` still has never been checked against a
  field whose integral is known -- which is the check survey section 55a asks for
  and is cheap and separable from everything else here.
- **The resolution recommendation depends on an inventory that does not exist.**
  Section 4b's decision procedure needs OCN-11 to report which straits and sills
  the Orogen mesh carries and how many each candidate grid resolves. Until that
  exists the recommendation is the coarsest defensible default rather than the
  answer the procedure returns.
- **The host was not quiet and no wall-clock number is quoted.** A sibling held a
  wall-clock claim on this machine throughout and the one-minute load ran between
  7 and 63 on 32 logical cores. cGENIE is one serial process, so it never waited
  for a core; retired instructions are a property of the binary and its input and
  every quantitative claim above rests on those. The profiling added one serial
  core to a contended machine and that is recorded rather than hidden.
