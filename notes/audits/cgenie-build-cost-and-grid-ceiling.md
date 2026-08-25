# cGENIE built and run here: what an ocean-year costs, and where the grid stops

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean, sea ice, wind stress and Courant number
name modelled quantities and numerical properties of a candidate component,
not observations of anything.

Measured 2026-08-25 on this machine. This is the measurement half of OCN-3 and
the whole of OCN-18. `vendor/cgenie` is a CANDIDATE under OCN-3: nothing reads
it, it has no pipeline row, and nothing here adopts it. The vendored tree was
not edited; everything below is build flags, namelist values and generated
configurations.

`analysis/cgenie_cost.py` is the driver and `analysis/cgenie_cost.json` and
`analysis/cgenie_stability.json` are what it wrote. The toolchain, the exact
make overrides, the code model and the host load are stamped in their
provenance blocks rather than repeated here.

**The headline.** It builds, it runs, and it reproduces the reference the tree
ships for it. The shipped 36 x 36 x 16 grid is not near any limit: the model's
own maximum Courant number there is far below one at the default timestep. What
caps the grid is not the ocean at all. Two independent limits sit in two
different components, they are controlled by two different namelist integers,
and at every grid tested the one that binds first is the ATMOSPHERE's.

---

# 1. What it took to build, and why one executable is not enough

Five things stood between the vendored tree and a running executable, the fifth
only at large grids. None of them is a source defect and none needed the tree
changed. The netCDF FORTRAN bindings, absent on this host until 2026-08-25 and
a separate package from the C ones, are the sixth and are recorded in
`docs/src/reference/vendored-upstreams.md` because they are a host fact rather
than a property of the tree.

**The tree expects to be somewhere else.** `genie-main/user.mak` sets
`GENIE_ROOT = $(HOME)/cgenie.muffin` and `RUNTIME_ROOT = ../../cgenie.muffin`,
and looks for netCDF under `/usr/local`. Three make overrides settle it, and
they are recorded in the artifact's `make_overrides`.

**It compiles one executable per GRID, exactly as ExoPlaSim compiles one per
(resolution, layers, ranks).** `genie-main/genie_control.f90` fixes the
atmosphere as `ilon1_atm = GENIENX, ilat1_atm = GENIENY` and
`genie-goldstein/src/fortran/ocean.cmn` fixes the ocean and sea ice as
`GOLDSTEINNLONS/NLATS/NLEVS`. All are cpp macros with `#ifndef` fallbacks, and
the fallbacks DISAGREE: the atmosphere falls back to the 64 x 32 IGCM grid and
the ocean to 36 x 36, so a bare `make` fails in `genie.F` where a sea-ice array
accumulates an atmosphere field, on shapes that are not conformable. That reads
like a source defect and is a missing configuration. The macros cannot be passed
on the make line either: `makefile.arc` takes them from `GENIE_FPPFLAGS`, and
only `genie.job` composes that, from the `<build>` block of a config file. So
the build has to run THROUGH a config, and CLAUDE.md rule 4 applies here too.

**The make target has to be `genie.exe`.** The default target also builds
`nccompare`, whose `src/c/compare.cpp` includes `netcdf.hh` from the legacy
netCDF C++ interface, which this host does not have.

**The build has to be serial.** The dependency generator `genie-main/finc.py`
is Python 2 and fails under this host's `python`, so no `.d` files are written
and a parallel make races on `.mod` files.

Above roughly a gigabyte of static COMMON the LINK fails with
`relocation truncated to fit: R_X86_64_PC32 ... defined in COMMON section`.
Every field cGENIE holds is in a named COMMON sized from the grid macros and
`user.mak` compiles with `-fno-automatic`, so all of it is static;
`genie-embm/src/fortran/embm.cmn` alone dimensions its seasonal arrays
`(maxi, maxj, maxnyr)` with `maxnyr = 400`. The 72 x 72 x 16 grid is past that
line and needs `-mcmodel=medium`. That is a build flag, not a source change, and
this document prices it separately, with a same-grid control in section 4, so it
does not sit inside a grid cost.

## 1a. There is no thread-count axis, and that is a property worth having

The component is serial, and more completely than the build flags suggest.
`makefile.arc`'s gfortran section carries `-fopenmp` commented out, but enabling
it would parallelise nothing: every `$omp` string in `genie-goldstein`,
`genie-embm`, `genie-biogem` and `genie-main` is commented out in column one, so
none is a sentinel and none is a directive. Three live `use omp_lib` statements
have no call against them. MPI appears only in
`genie-plasim/src/fortran/plasimmod.f90`, the vendored PlaSim this project does
not use. `notes/audits/ocean-and-marine-biosphere.md` section 9h reads the two
commented blocks and finds them abandoned plumbing experiments rather than a
decomposition of the physics; this adds only that there is nothing to enable.

So cost is one core against simulated years, cores cannot absorb a grid
refinement, and parallelism is across EXPERIMENTS rather than within a run.
**That is why an offline ocean can run beside a commissioning instead of
queueing behind it**: it wants one core, and the atmosphere's thread team wants
the rest.

## 1b. It reproduces the reference the tree ships for it

`make testebgogs` cannot answer this on its own, for two reasons that are both
about its own configuration. Its comparison tool is `src/c/nccompare.exe`, which
needs the legacy netCDF C++ header; and its reference is a GOLDSTEIN
annual-average netCDF that `goldstein.F:732` writes only when the namelist
`debug_loop` is true, which defaults false. The shipped test therefore runs,
writes no reference-shaped output, and reports a failure that is about the test.
Setting `debug_loop` changes no physics: every use of it in `goldstein.F` guards
a print, a diagnostic dump or the averaging call, never the state.

With that set, the shipped `eb_go_gs` case runs and its output is compared field
by field against `genie-knowngood/`. The result is in the artifact's
`knowngood` block.

---

# 2. The instrument, and the criterion fixed before the runs

`initialise_goldstein.F:515` sets the tracer timestep as
`tv = sodaylen*yearlen/(nyear*tsc)` and assigns it to every level. `nyear` is a
plain namelist integer defaulting to 100, and there is no Courant test, no CFL
test and no stability test anywhere in `genie-goldstein` or `genie-embm`. So a
criterion has to be applied from outside, and it was fixed in the driver's
docstring before any run: **the advective Courant number must stay below one and
the diffusive numbers below one half**, evaluated at a nominal surface current.

The model turns out to carry a better instrument than that estimate.
`genie-goldstein/src/fortran/diag.f` forms

    tv = max(|u|*rc(j)*rdphi, |w|*rdz(k)) * dt(k)
    if (j < jmax) tv = max(tv, |v|*cv(j)*rdsv(j)*dt(k))
    cnmax = max(cnmax, tv)

over every wet cell and prints it as `Cn`, from the flow GOLDSTEIN actually
diagnosed rather than from a nominal speed. **Nothing acts on it**, and `diag`
is called only behind `debug_loop`, so at the shipped settings the number is
never computed. That is the trap this section exists to name: a run past the
limit prints nothing, completes, and reports success. **A grid that RUNS is not
a grid that is STABLE**, and the difference is invisible unless the diagnostic
is switched on.

The criterion applied below is the same one in the same form, read off the
better instrument: **a configuration is STABLE when the model's own `Cn` stays
below one at every diagnostic step and the tracer field stays finite and
bounded.** The bound on temperature is a blow-up detector rather than a physical
range, so a marginal but real circulation is never called unstable on it alone.

## 2a. There IS one guard, it is in the atmosphere, and diagnostics disarm it

`genie-embm/src/fortran/surflux.F:900-954` solves the sea-ice surface
temperature by Newton iteration. When it does not converge it prints
`warning sea-ice iteration failed at`, dumps the state, and executes a bare
Fortran `stop` -- **unless EMBM's own `debug_loop` is set**, in which case it
clamps `tice` and carries on.

This is the only hard stability guard in the physics, it is not in the ocean,
and turning EMBM diagnostics on converts it into a silent clamp. Every sweep
below therefore sets `debug_loop` on GOLDSTEIN ONLY: the ocean's `diag` runs and
the atmosphere's guard stays live.

---

# 3. Two limits, two components, two namelist integers

`nyear` is the ocean's tracer timestep. It is NOT the atmosphere's.
`initialise_embm.F:578` sets `dtatm = tv/ndta`, and `ndta` is a separate
namelist integer in `ini_embm_nml` defaulting to 5, so EMBM takes `ndta`
sub-steps per ocean step. Whether a refinement forces the OCEAN's step down or
only EMBM's is the entire ceiling question, because they cost differently.

## 3a. Sweeping nyear alone finds a boundary, and it is the atmosphere's

At the default `ndta = 5`, with the ocean's `diag` on and EMBM's guard live, the
boundary is sharp on both grids and every failure below it is the same failure:
EMBM's sea-ice surface solve does not converge and `surflux.F` stops the model.
The ocean's own `Cn` at the last stable point is nowhere near one.

| grid | last unstable nyear | ocean dt there | first stable nyear | ocean dt there | model's own Cn at the first stable point |
| --- | ---: | ---: | ---: | ---: | ---: |
| 36 x 36 x 16 | 30 | 12.2 d | 35 | 10.4 d | 0.39 |
| 72 x 72 x 16 | 140 | 2.61 d | 160 | 2.28 d | 0.33 |

Two things fall out. The ocean has a factor of about three of unused Courant
headroom at the point where the model stops, on BOTH grids, so the ocean is not
what stops it. And the boundary moves by a factor of about 4.6 for a factor of 2
in resolution: close to quadratic, and not the linear scaling an advective limit
would give.

A control was run for the obvious objection, that switching GOLDSTEIN's
diagnostic on is itself doing something. It is not: 72 x 72 x 16 at
`nyear = 100, ndta = 5` stops at the same EMBM step in the same cell,
`(36, 43)`, with `debug_loop` false and with it true. The instability is a
property of the configuration and the diagnostic only makes it visible.

## 3b. The boundary is a property of dtatm alone, and it was predicted before it was run

If the boundary belongs to EMBM it should depend only on `dtatm = dt_ocean/ndta`
and not on the ocean's step at all. Read off the 72 x 72 rows above, EMBM's
boundary sits at `dtatm` between 0.457 d and 0.522 d. That predicts, before any
run: at `nyear = 100` the model is unstable for `ndta <= 7` and stable for
`ndta >= 8`, and at `nyear = 50` unstable for `ndta <= 14`.

Every one of those held. The sharpest evidence is not the boundary but a
coincidence the hypothesis requires: `nyear = 100, ndta = 5` and
`nyear = 50, ndta = 10` are the same `dtatm` at ocean timesteps a factor of two
apart, and both stop **at the same EMBM step, in the same cell**. The ocean's
step is not in it.

So a finer ocean grid does not force the ocean's timestep down. It forces
EMBM's. And EMBM is a one-layer atmosphere on the same horizontal grid, so
sub-stepping it is cheap: over `ndta` from 8 to 20 at 72 x 72 the wall clock of a
ten-year integration does not trend, which says the ocean's three-dimensional
tracer work is the cost and EMBM is not.

## 3c. Take the atmosphere out of the way and the ocean's own limit appears

The prediction written down beforehand had a fourth clause, and that one FAILED,
which is how the second limit was found. It said that at `nyear = 50` on
72 x 72 the model would become stable once `ndta` reached 16. Raising `ndta` to
16 does exactly half of that: EMBM's solve stops failing, as predicted. The run
is still not stable, and now for the other reason. The model's own `Cn` reaches 1.61, over the criterion,
and the same configuration at `ndta = 20` gives 1.64.

That is the ocean's own ceiling, and it behaves as an advective Courant number
should. At 72 x 72 with `ndta` converged, `Cn` is 0.81 at `nyear = 100` and 1.64
at `nyear = 50`: a ratio of 2.01 for a factor of two in timestep, linear to
within a percent. Solving for `Cn = 1` puts the ocean's own largest usable
timestep at 72 x 72 x 16 at about 4.5 days, `nyear` about 81.

The same construction at 36 x 36 x 16, where `Cn` is 0.134 at the shipped
`nyear = 100`, puts the ocean's own limit there at about 27 days. **The shipped
grid runs at a seventh of the timestep its own ocean could carry**, and what
holds it there is EMBM.

## 3d. What the connector says about the same ceiling

`references/muffingen` is the generator for the `.k1`, `.paths` and `.psiles`
files a new geography needs, and it is where OCN-18's second half sits.

- **Its declared range is `[1-72]`** for `par_max_i` and `par_max_j`, annotated
  on every example configuration. So the grid probed here is exactly the top of
  what the connector claims, and 131 of the shipped cGENIE configurations sit at
  36 x 36 with nothing above it. A declared range is not a demonstrated one, and
  this document demonstrates only that the MODEL runs at 72 x 72.
- **Islands are found automatically, not drawn.** `muffingen.m:1033` calls
  `find_grid_islands` and `:1054` `find_grid_islands_update`; the `ginput` calls
  in `source/fun_grid_edit_*.m` are an interactive editor for corrections, not
  the primary path. The manual records that the ORIGINAL generator required
  hand-drawn island paths and that getting them wrong could leave the
  circulation unsolved somewhere without any complaint, which is why the
  automatic version exists.
- **It is GPL-3**, where cGENIE is MIT. Consuming its output is unencumbered;
  vendoring the generator would bring a different licence into this tree.
- **Whether it runs under Octave is still not demonstrated.** The manual says
  MATLAB throughout and states no position either way. The only interface that
  would decide it is the low-level `netcdf.*` family, used in about 170 places,
  which Octave provides through a separate package; nothing here has run it.
  That is the one part of OCN-18 this document does not close.

The island count itself is a compile-time bound, `GOLDSTEINMAXISLES`, and the
sweeps here compile with it raised well past the shipped default, so it is a
build parameter rather than a ceiling.

## 3e. How the criterion fixed in advance actually did

Worth recording, because the point of fixing a criterion before the run is that
it can then be wrong in a way you can see. The a-priori form -- advective
Courant below one at a nominal 0.5 m/s surface current -- puts the boundary at
`nyear` about 61 on 36 x 36 and about 171 on 72 x 72. Measured, it is between 30
and 35, and between 140 and 160.

So it was conservative by about 1.8x on the shipped grid and by about 1.1x on
the doubled one, and it named the right order of magnitude both times. It did
that for the wrong reason: a nominal-speed advective criterion is LINEAR in
resolution and the limit that actually binds is quadratic, so the two agree only
where they happen to cross, which is near 72 x 72. The criterion was a usable
predictor and is not the explanation, and the difference is exactly what the
model's own `Cn` and EMBM's guard were needed to tell apart.
