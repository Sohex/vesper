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
`initialise_embm.F:578` sets EMBM's step as `dtatm = dt_ocean/ndta`, and `ndta`
is a separate namelist integer in `ini_embm_nml` defaulting to 5. Whether a
refinement forces the OCEAN's step down or only EMBM's is the entire ceiling
question, because the two cost differently.

**`ndta` is not a sub-step knob on its own, and getting that wrong produces
runs that are stable and meaningless.** How often EMBM is CALLED is `katm_loop`
at `genie.F:300`, which defaults to 1, so EMBM steps once per genie step. With
`kocn_loop` genie steps per ocean step, EMBM advances `kocn_loop/ndta` times the
ocean's time over the same interval, and that is one only when
`kocn_loop == ndta`. The shipped default pairs `ndta = 5` with `kocn_loop = 5`
for exactly that reason. Raising `ndta` alone slows the atmosphere's clock
relative to the ocean's rather than sub-stepping it, and nothing in the model
complains. Every configuration below therefore moves `kocn_loop` with `ndta`.

## 3a. Both grids stop, both times in the atmosphere, and the ocean is nowhere near its limit

`analysis/cgenie_stability.json` is `nyear` crossed with `ndta` on both grids,
with the ocean's `diag` on and EMBM's guard live. Two things are true of every
cell of it.

The ocean's own `Cn` depends on `nyear` and on nothing else: at 36 x 36 x 16 it
is 0.13 at `nyear = 100` whether `ndta` is 2, 5, 10 or 20, and 0.52 at 72 x 72
under the same four. That is what an advective Courant number has to do, and it
is the check that the sweep is measuring what it claims to.

And every failure is EMBM's sea-ice surface solve stopping the model. At the
shipped `ndta = 5`:

| grid | last unstable nyear | ocean dt there | first stable nyear | ocean dt there | ocean Cn at the first stable point |
| --- | ---: | ---: | ---: | ---: | ---: |
| 36 x 36 x 16 | 30 | 12.2 d | 35 | 10.4 d | 0.39 |
| 72 x 72 x 16 | 140 | 2.61 d | 160 | 2.28 d | 0.32 |

So the ocean has a factor of about three of unused Courant headroom at the point
where the model stops, on BOTH grids. It is not what stops it.

A control was run for the obvious objection, that switching GOLDSTEIN's
diagnostic on is itself doing something. It is not: 72 x 72 x 16 at
`nyear = 100, ndta = 5` stops at the same EMBM step in the same cell,
`(36, 43)`, with `debug_loop` false and with it true. The diagnostic makes the
instability visible; it does not cause it.

## 3b. The atmosphere's boundary is a property of its own timestep, and only that

Across the whole grid, EMBM fails exactly when `dtatm` is above a threshold that
does not depend on the ocean's step:

| grid | largest dtatm that failed | smallest dtatm that survived |
| --- | ---: | ---: |
| 36 x 36 x 16 | 2.435 d | 2.087 d |
| 72 x 72 x 16 | 0.522 d | 0.457 d |

Every one of the 64 cells is on the correct side of that. The sharpest evidence
is not the boundary itself but a coincidence the hypothesis requires: pairs with
the SAME `dtatm` and ocean timesteps a factor of two apart fail in the same
cell. At 72 x 72, `nyear = 50, ndta = 10` and `nyear = 100, ndta = 5` are both
`dtatm = 0.7305 d` and both stop in cell `(36, 43)`; `nyear = 25, ndta = 10` and
`nyear = 50, ndta = 5` are both `1.4610 d` and both stop in cell `(1, 19)`. The
ocean's step is not in it.

The threshold falls by a factor of between 4.0 and 5.3 for a factor of 2 in
resolution. That is quadratic, not the linear scaling an advective limit would
give, and EMBM's step is semi-implicit with a FIXED four iterations
(`tstipa.f:44`, `nii = 4`, `cimp = 0.5`) against a temperature eddy diffusivity
whose shipped amplitude is 5.0e6 m2/s. The iteration count was chosen for the
shipped grid and does not move with it.

## 3c. Sub-step the atmosphere properly and the ocean's own limit appears

Because EMBM's limit is on `dtatm`, raising `ndta` should relieve it and expose
whatever is behind it. It does. At 72 x 72, `nyear = 100` is unstable at
`ndta = 5` and stable at `ndta = 10` with the ocean's `Cn` at 0.52. Pushing
further, `nyear = 50, ndta = 20` puts `dtatm` at 0.365 d, comfortably inside
EMBM's threshold, and EMBM's solve does not fail. The run is still not stable,
and now for the other reason: the model's own `Cn` reaches 1.04, over the
criterion.

That is the ocean's own ceiling. It behaves as an advective Courant number
should: at 72 x 72 `Cn` is 0.52 at `nyear = 100` and 1.04 at `nyear = 50`, a
ratio of 2.00 for a factor of two in timestep. So the ocean's own largest usable
timestep at 72 x 72 x 16 is about 7.0 days, `nyear` about 52.

The same construction at 36 x 36 x 16, where `Cn` is 0.13 at `nyear = 100`,
puts the ocean's own limit there at `nyear` about 13. **At both grids the ocean
could carry a timestep about three times the one EMBM permits**, and the two
limits scale together: the ocean's `Cn` at fixed `nyear` is 4.00 times larger at
72 x 72 than at 36 x 36, the same quadratic as EMBM's.

## 3d. What that makes the cost of a refinement

Both limits are quadratic in resolution and both components do work proportional
to cell count, which is also quadratic. Refining the horizontal grid by a factor
r therefore costs r^2 in cells and r^2 in timesteps in EACH component:
**cost per model year goes as r^4**, whichever limit is being respected.

`notes/audits/ocean-and-marine-biosphere.md` section 9g estimated r^4 by reading
`ubarsolv`'s loop structure and attributing the growth to the barotropic solve.
The exponent is right and the attribution is not: the growth is the TIMESTEP, in
the atmosphere first and the ocean second, and section 4 finds the barotropic
solve too small to see at these grids.

## 3e. How the criterion fixed in advance actually did

Worth recording, because the point of fixing a criterion before the run is that
it can then be wrong in a way you can see. The a-priori form -- advective
Courant below one at a nominal 0.5 m/s surface current -- puts the boundary at
`nyear` about 61 on 36 x 36 and about 171 on 72 x 72. Measured at the shipped
`ndta = 5`, the boundary is between 30 and 35, and between 140 and 160.

So it was conservative by about 1.8x on the shipped grid and by about 1.1x on
the doubled one, and it named the right order of magnitude both times. It did
that for the wrong reason twice over: the limit it describes is the ocean's,
which is a factor of three further away than where the model actually stops, and
a nominal-speed advective criterion is LINEAR in resolution where both real
limits are quadratic. The two agree only where they happen to cross, near
72 x 72. The criterion was a usable predictor and was not the explanation, and
telling those apart is exactly what the model's own `Cn` and EMBM's guard were
needed for.

## 3f. What the connector says about the same ceiling

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

