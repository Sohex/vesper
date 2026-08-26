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

**On the instrument.** This host runs several agents at once, and for most of
the window in which the timings below were taken a sibling was integrating an
atmosphere on sixteen threads at a load of 16 to 28 on 32 logical cores. Wall
clock under that is a measurement of a machine state, so the durable price of an
ocean-year is quoted in RETIRED INSTRUCTIONS, which are a property of the binary
and its input and are the same on a busy host and an idle one. cGENIE is
strictly serial, so there is not even a spin-wait to correct for. Seconds are
quoted beside them for whoever has to wait, with the load they were taken under
and marked contaminated where they were.

**The headline.** It builds, it runs, and it reproduces the reference the tree
ships for it to the last bit that matters. The shipped 36 x 36 x 16 grid is
nowhere near the ocean's own limit. What caps the grid is not the ocean at all:
two independent limits sit in two different components, controlled by two
different namelist integers, and at every grid tested the one that binds first is
the ATMOSPHERE's. Both are quadratic in resolution, so a refinement costs r^4 --
which is the exponent the audit already estimated, for the wrong reason. A
doubled 72 x 72 x 16 ocean runs and is affordable; matching a T85 atmosphere is
not, and the constraint that survives is therefore the support mismatch rather
than the price.

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
by field against `genie-knowngood/`. All 28 floating-point variables present in
both files agree, and the worst difference over all of them, relative to that
field's own range, is 4.3e-21. Nothing is present in the reference and missing
from the run. **The build here reproduces cGENIE's own reference output**, which
is the strongest single statement this document can make about whether it works,
and it is the acceptance criterion OCN-19 will need.

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

## 3a. Both grids stop in the atmosphere, and the ocean is nowhere near its limit

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
could carry a timestep between two and a half and three times the one EMBM
permits**, and the two limits scale together: the ocean's `Cn` at fixed `nyear`
is 4.00 times larger at 72 x 72 than at 36 x 36, the same quadratic as EMBM's.

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
files a new geography needs, and it is where OCN-18's second half sat. Run here
2026-08-25 under GNU Octave 11.3.0 for world-crky; the driver, the settings
files and the sixteen run logs are the reproduction recipe below.

- **It runs outside MATLAB, headless, with no netCDF package installed.**
  Sixteen configurations of muffingen v0.9.26 completed under `octave-cli`
  11.3.0 with `opt_user=false`, `opt_plots=false` and no display. The Octave
  install carried no packages at all, which settles the netCDF question by
  demonstration on the arm this project would use, and settles the arm's
  graphics need as well.
- **THE KNOWN-GOOD ANSWER REPRODUCES EXACTLY.** `genie-paleo/wworld` is a
  36 x 36 world generated under MATLAB by muffingen v0.64, and it ships its own
  generator configuration and run log. Re-run from that configuration, Octave
  wrote a `.k1`, a `.paths` and a `.psiles` BYTE-IDENTICAL to the shipped three,
  and its log reports the same island result: both pole islands found and
  uncounted, `total # true islands = 1`. The wind-stress and albedo files differ,
  and that is the v0.64 to v0.9.26 generator change rather than the interpreter:
  the newer version splits planetary from cloud albedo, renames the output, and
  labels its winds "idealized". Reproducing a known-good answer before attempting
  a new one is what makes the higher-resolution runs interpretable.
- **THREE INCOMPATIBILITIES, NONE OF THEM netCDF AND NONE OF THEM A TOOLBOX.**
  `muffingen.m` queries `get(0,'Diary')`, which is a MATLAB root-object property
  Octave's root does not carry; the call raises before anything else runs, and
  `diary off` unconditionally is equivalent in both. The other two are
  unguarded plots: `source/make_grid_winds_zonal.m` emits two zonal wind-stress
  profile figures with no `opt_plots` test -- the flag is not even a parameter of
  that function, so no settings file can suppress them -- and `muffingen.m`
  emits the generic zonal-mean planetary albedo and cloud albedo profiles the
  same way. Every other figure in the tree is behind `opt_plots`. On a host with
  a working graphics toolkit those three would cost a `.ps` file each and
  nothing else; this host has none, because Octave 11.3.0 is installed without
  `fltk` and without `gnuplot`, so `available_graphics_toolkits` is empty and
  `figure` raises. Installing either package would let muffingen run unpatched
  except for the `Diary` line.
- **`par_max_i` AND `par_max_j` ARE IGNORED ON THE ARM THIS PROJECT WOULD USE,
  and this is the finding that matters most for the connector.** On `par_gcm` of
  `mask`, `k1` or `k2`, `muffingen.m` calls `fun_read_k1`, which returns
  `[jmax,imax] = size(gk1)` from the input file, and then re-generates the GENIE
  grid from those. The code says so: "imax and jmax are deduced from the file".
  Twelve runs that raised `par_max_i` and `par_max_j` to 48 x 40 and 72 x 72 over
  the shipped 36 x 36 example masks all produced 36 x 36 output. The connector's
  resolution is set by the MASK OROGEN WRITES, not by a muffingen parameter, and
  the `[1-72]` annotated on those two parameters governs only the GCM arms.
  `par_max_k` is honoured everywhere.
- **Given genuinely high-resolution input it produces a complete, well-formed
  world at 48 x 40 and at 72 x 72.** A 48 x 40 mask gives a `.k1` of 50 by 42,
  exactly the shape of the shipped `fm0000bb.k1`, with the same border rows, the
  same duplicated zonal wrap columns, land runoff codes in 91 to 94, and a
  `.psiles` of 48 by 41. A 72 x 72 mask gives 74 by 74 and 72 by 73. Nothing
  degrades and nothing warns.
- **The automatic island machinery gets the right answer at 72 x 72.** A
  synthetic 72 x 72 world built with six separate landmasses -- a meridional
  barrier clear of both poles, four isolated islands and a north polar cap --
  was resolved into exactly six: six land masses found, six true islands, six
  paths built as `#2` through `#7` with `#1` the ignored border, the south pole
  correctly left as the uncounted open-ocean reference, and `.psiles` indices
  running 0 to 7. This is the check the manual's warning asks for, at the
  resolution that was only declared before.
- **The `k1` arm round-trips bathymetry exactly, and the `mask` arm has none.**
  A 72 x 72 `.k1` carrying thirteen distinct depth levels came back with the
  mask preserved cell for cell and every level unchanged, and the island machinery
  ran on it with the corner handling engaged. The `mask` arm cannot do this: it
  sets a uniform ocean depth, so its `.k1` is a flat ocean at `par_min_k`. So the
  Orogen connector writes a `.k1`, not a mask `.dat`, and the land runoff codes
  it carries pass through unchanged rather than being recomputed.
- **What this did NOT establish.** No generated world was integrated by cGENIE.
  Usability is claimed structurally -- shape, border convention, wrap columns,
  land codes, level range, island count against `.psiles` and `.paths` -- against
  the shipped `fm0000bb` and `wworld` as the format reference, and not by a run.
- **It is GPL-3**, where cGENIE is MIT. Consuming its output is unencumbered;
  vendoring the generator would bring a different licence into this tree.

The island count itself is a compile-time bound, `GOLDSTEINMAXISLES`, and the
sweeps here compile with it raised well past the shipped default, so it is a
build parameter rather than a ceiling.

**To reproduce.** muffingen is GPL-3 and is not vendored, so the recipe is
stated rather than a script kept in the tree. Copy `references/muffingen`'s
`muffingen.m`, `source/`, `DATA/` and the input mask outside the repository;
replace the `get(0,'Diary')` line with `diary off;`; comment out the two
figure blocks in `source/make_grid_winds_zonal.m` and the two albedo-profile
figure blocks in `muffingen.m`, or install a graphics toolkit instead; write a
settings file with `par_gcm='k1'`, an eight-character `par_wor_name`,
`opt_user=false` and `opt_plots=false`; and run
`octave-cli --no-gui --eval "muffingen('<settings>')"`. The input file's own
dimensions set the output resolution.

---

# 4. What an ocean-year costs

Every case is a topography that ships COMPLETE at its stated grid: `.k1`,
`.paths`, `.psiles`, the four wind-stress components, the two wind-speed fields
and the restoring climatologies. The level count is a property of the topography
rather than a free parameter, because a `.k1` entry is the index of the level a
column's floor sits on, so the highest wet value in the file IS the grid's level
count and land is coded above it. The vertical arm therefore changes topography
with the levels.

Each grid is built once and then priced twice, at 20 and at 100 model years, so
that the fixed initialisation subtracts out of the difference. The price is
retired instructions, for the reason in the opening: this host was busy and
seconds measure that.

| topography | grid | wet cells | code model | ndta | Ginstr per model year | initialisation, in model-years of work |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| `worbe2` | 36 x 36 x 8 | 6210 | small | 5 | 1.776 | -0.66 |
| `worjh2` | 36 x 36 x 16 | 12511 | small | 5 | 2.973 | -1.09 |
| `worri4` | 36 x 36 x 32 | 25042 | small | 5 | 5.330 | -1.32 |
| `g3660l` | 36 x 60 x 8 | 9108 | small | 5 | 2.805 | -0.64 |
| `igcmv3` | 64 x 32 x 8 | 8688 | small | 5 | 2.681 | -0.63 |
| `worjh2` | 36 x 36 x 16 | 12511 | medium | 5 | 3.061 | -1.12 |
| `dan_72` | 72 x 72 x 16 | 49680 | medium | 10 | 13.781 | -0.98 |

**Initialisation is not a cost.** Every intercept is smaller than one model
year's worth of work and every one is NEGATIVE, which says the fixed cost is
below the resolution of this fit and that the per-year cost rises very slightly
with run length instead. The two-length design was there to separate the two and
its answer is that there is nothing to separate: a pass pays for its years and
almost nothing for starting.

**The code model costs 3 per cent.** The same topography, grid and settings
compiled `-mcmodel=medium` instead of the shipped `small` is 1.030 times the
instructions. That is the whole of the 72 x 72 build flag's price, and it is why
the doubled grid's number is not contaminated by it.

## 4a. The barotropic solve is not what runs at these grids

This arm was designed before it ran, and it is the one place where the answer
could have gone either way. `ubarsolv.f` eliminates over `do i=1,n*m-1` with an
inner width of at most `n+1`, where `n = imax`: work per call is about
`imax^2 * jmax`, QUADRATIC in longitudes and only LINEAR in latitudes. Tracer
work in `tstepo` is `imax*jmax*kmax*lmax`, symmetric in the two. Two shipped
topographies have nearly the same cell count and very different barotropic work:

    36 x 60 :  imax^2 * jmax =  77760
    64 x 32 :  imax^2 * jmax = 131072    ratio 1.686

The decision rule, fixed in the driver's docstring before the runs: above 1.30
reads barotropic-dominated, below 1.10 reads tracer-dominated, and between them
leaves the question open.

**Measured, 2.681 / 2.805 = 0.956.** That is below 1.10, so the arm returns
TRACER-DOMINATED, and it lands within 0.6 per cent of the 0.95 that hypothesis
predicted. A 1.69-fold increase in barotropic work costs nothing measurable:
at these grids `ubarsolv` is invisible.

That settles the runtime half of
`notes/audits/ocean-and-marine-biosphere.md` section 9h without a profiler, and
it is why section 3d relocates the r^4 growth onto the timestep. It does NOT
settle where the time goes WITHIN the tracer work, which is OCN-19's.

## 4b. Levels are cheaper than cells, and the doubled grid is dearer

Doubling the level count costs 1.674 (8 to 16) and 1.793 (16 to 32), both below
the 2.0 that cell count alone would give, because the two-dimensional work --
the surface fluxes, the barotropic solve, EMBM -- is shared across levels.

Doubling both horizontal dimensions costs 4.502 at the same `nyear`, against a
wet-cell ratio of 3.971: 1.13 times per cell. Combined with section 3's result
that the usable timestep falls as r^2, a doubling of the horizontal grid costs
about 18 times per model year, against the r^4 = 16 that both scalings predict.

**The instructions do convert to seconds at a rate that does not change with the
grid**, which is worth stating because a first pass suggested otherwise. Taken
from the same `perf` runs, throughput is 11.4 to 18.3 G instructions per second
across every physics case INCLUDING 72 x 72 x 16, which sits at 15.5 in the
middle of that band. An earlier pass put the same case at 6.3, and that was the
contended window rather than a property of the grid: a grid needing
`-mcmodel=medium` looked like a candidate for its state leaving cache, and
re-measuring at lower load says it is not. So r^4 in instructions is r^4 in
seconds for the physics.

## 4bb. The instrument, checked against the size of the effect

The same seven cases were measured twice, hours apart, at one-minute load
averages spanning 8 to 42. **The instruction counts agree to within 0.04 per
cent**, against effects here of 1.67 to 4.50. That is three orders of magnitude
below the smallest thing being claimed, and it is why the numbers above are
quoted in instructions and the seconds are not.

## 4c. Storage is nothing, and the biogeochemistry is where the published costs are

A physics-only run writes 564 kB whatever its length, because these
configurations suppress the periodic output; storage is a function of what is
asked for, not of the integration.

The shipped BIOGEM regression case, `configs/eb_go_gs_ac_bg_test.xml` run
unaltered except for its length, is the same 36 x 36 x 8 grid as the cheapest
physics case with fourteen GOLDSTEIN tracers, ATCHEM and BIOGEM on. It writes
9.8 MB at ten model years and 23.8 MB at fifty, so about 350 kB per model year,
and that is the storage figure a spin-up should be planned against rather than
the physics-only one.

Priced the same way, that configuration costs **4.764 Ginstr per model year**
against the physics-only 1.776 at the same grid: the ocean biogeochemistry is
**2.68 times the physics**, in instructions.

In seconds it is nearer five and a half times, and that gap is real rather than
an artifact of the host. BIOGEM retires 5.94 G instructions per second where
every physics case is between 11.4 and 18.3, and the two independent
measurements of it, hours apart at different loads, gave 5.96 and 5.94. So this
is a property of the geochemistry -- long-latency arithmetic in the carbonate
system, and a much larger working set for fourteen tracers -- and not of the
machine. **The biogeochemistry costs 2.7 times the instructions and 5.4 times
the seconds**, and a spin-up should be planned on the seconds.

## 4d. This reproduces the published EMIC costs, within a hardware generation

`notes/audits/ocean-and-marine-biosphere.md` section 9a quotes Capirala and
Olson (2026) at 24 to 48 hours on a single core for 36 x 36 x 16 to physical and
biogeochemical steady state at 20,000 model years, and Edwards and Marsh (2005)
at 20 kyr in about a day on a PC.

Taking the BIOGEM price above and scaling it to sixteen levels by the physics
case's own level ratio -- an ASSUMPTION, since the biogeochemistry was not run
at sixteen levels here -- gives about 8.0 Ginstr per model year, so 20,000 years
is about 1.6e14 instructions, and at the 5.94 Ginstr per second the
biogeochemistry runs at, about 7.5 hours on one core of this machine.

That is a factor of three to seven faster than the published figure, on a part
about twenty years newer than Edwards and Marsh's and several generations newer
than Capirala and Olson's. **The published costs are reproduced rather than
merely quoted**, which is what OCN-3 asked for, and the physics-only figure that
sits underneath them is smaller again: 2.973 Ginstr per model year at
36 x 36 x 16, about 1.7 hours for the same 20,000 years.

---

# 5. What this means for a T85 atmosphere

`notes/audits/ocean-and-marine-biosphere.md` section 6 chose offline coupling on
one property: the ocean's cost must not climb the T21/T42/T85/T127/T170 ladder
with the atmosphere. This is the measurement that decision rests on, so it owes
two statements rather than one.

## 5a. The ocean's grid is independent of the atmosphere's rung, and now measured

Nothing in section 4's table depends on an atmosphere resolution. EMBM runs on
the OCEAN's grid, not the atmosphere's, and the chosen coupling hands the ocean
a climatology rather than a synchronous partner. So the ocean's price is a
function of its own grid and its own timestep, both fixed here, and moving the
atmosphere from T21 to T170 does not touch a number in it.

That is exactly why the published synchronous coupling was eliminated: it gets
its affordability from matching the three horizontal grids, and matching at T170
puts a serial ocean at 512 x 256.

## 5b. Matching the atmosphere is not affordable, and this says by how much

T85 is 128 latitudes by 256 longitudes (`lib/rungs.py`). Against 36 x 36 that is
25.3 times the surface cells, an equivalent uniform refinement of about 5.0.
Section 3d's r^4 -- r^2 in cells and r^2 in timestep, in both components -- puts
a T85-matched ocean at about 640 times the shipped grid's price per model year,
so around 1900 Ginstr per model year against 2.973, and a 20,000-year physics
spin-up at about 3.8e16 instructions. At the 11 to 18 Ginstr per second the
physics runs at, that is between three and six weeks on one core, and carrying
biogeochemistry alongside multiplies it again by about five in seconds.

**That figure is an extrapolation across a factor of five in resolution from a
single measured doubling, and it should be read as an order of magnitude.** What
it is enough to settle is the decision it was taken for: matching is not
affordable, offline coupling on an ocean grid chosen for the ocean is, and the
eliminated architecture's 512 x 256 at T170 is worse again by another factor of
sixteen.

## 5c. The constraint that survives is the support mismatch, not the cost

Because the ocean does not have to match, its cost is not the binding
constraint. What remains is that any adequacy claim about coupled behaviour is a
claim about fields carried across a large jump in support.

At T85 the atmosphere resolves 1.4 degrees of longitude. The shipped ocean grid
resolves 10, and it is equal-area in the sine of latitude, so its rows are about
3.2 degrees near the equator and much wider in latitude towards the poles. Every
ocean field crossing to the atmosphere is therefore a 10-degree field being read
on a 1.4-degree grid, a factor of about seven in longitude, and every atmosphere
field crossing the other way is an average over about fifty atmosphere cells.

The doubled grid halves that to a factor of about three and a half, and section
3 prices the halving at roughly eighteen times per model year, or more once the
throughput effect in 4b is settled. **So a finer ocean is affordable and a
MATCHED one is not, and the honest form of an adequacy claim is a declared
support mismatch rather than a resolution that was chosen and then not
mentioned.** Section 8b's straits, sills and partial coasts ride entirely on
this, and at 10 degrees they are not represented at all.

## 5d. EMBM is a second atmosphere AND the thing that caps the ocean's grid

Section 9f already records that adopting this host means the modelled land and
the modelled ocean read different atmospheres, and that OCN-10 has to name which
terms EMBM owns. This measurement adds a harder reason to care. EMBM is not only
a second atmosphere whose ownership is unstated; it is the component that sets
the ocean's timestep and therefore prices every grid refinement. At both grids
tested the ocean could carry a timestep between two and a half and three times
the one EMBM permits.

That has a consequence specific to the offline architecture. In offline
full-flux coupling the ocean is driven by an ExoPlaSim climatology, so what
EMBM's prognostic atmosphere is still FOR has to be answered rather than
inherited -- and the answer changes the ceiling, because a configuration that
does not integrate EMBM does not inherit its stability limit and gets the
factor of three back. Nothing here establishes that such a configuration exists
or is correct. Section 6 is the reading of what the recipe would have to do.

---

# 6. Three viability findings this measurement walked into

Configuring the runs meant reading how `genie.F` assembles a recipe and what
GOLDSTEIN hands back. Three things there bear on the coupling architecture
rather than on cost, and each answers a question OCN-3 already asks.

**There are exactly two surface-flux paths for GOLDSTEIN, and each is gated on
an atmosphere module being in the recipe.**

- `genie.F:279` calls `surflux_wrapper` under `flag_ebatmos .and.
  flag_goldsteinocean`. That is EMBM's, and it computes the fluxes from EMBM's
  own prognostic state.
- `genie.F:402` calls `plasim_surflux_wrapper` under `flag_plasimatmos`, and it
  sits INSIDE the `flag_goldsteinseaice` block at `:397`.

There is no third path, and with neither atmosphere flag set GOLDSTEIN receives
no surface forcing at all. Every one of the shipped configurations carries an
atmosphere for this reason: of the ones in `genie-main/configs`, all include
either `embm` or `plasim`, and all include `goldsteinseaice`.

**The second path is nevertheless the shape offline full-flux coupling wants.**
`genie_loop_wrappers.f90:154-175` shows `surflux_goldstein_seaice` taking
surface temperature, humidity and pressure, insolation, downward longwave, net
heat, wind speed and the latent and sensible transfer coefficients as INPUT, and
returning the ocean's latent, sensible, net solar and net longwave fluxes. It
does not need PlaSim integrated; it needs those fields to exist. So a driver
supplying an ExoPlaSim climatology has an interface to aim at rather than a
surface to invent.

**Two things about that interface have to be settled before it is called a
route, and they belong to OCN-10 and to section 7b rather than here.** It also
returns `dhght_sic`, `dfrac_sic`, `temp_sic` and `albd_sic`: the path is the
SEA-ICE one, and it computes an ice state, which is exactly the quantity section
7b makes ExoPlaSim's `icemod` authoritative for. And it is reachable only with
`flag_goldsteinseaice` set, so "run the host's sea ice diagnostically or not at
all", which OCN-3's own notes ask about, is not a namelist choice on this path:
the sea-ice module is what carries it.

## 6a. The ocean's own ceiling moves with a declared bracket

`scf` scales the wind stress, and it enters linearly: `goldstein.F:199-205` sets
`dztau = scf * stressxu_ocn` and the three companions, and the same four lines
appear in both surface-flux paths. GOLDSTEIN's flow is frictional-geostrophic
and linear in the stress, so the diagnosed velocity, and with it the model's own
`Cn`, is very nearly proportional to `scf`.

OCN-17's verdict is that `scf` is a DECLARED BRACKET, conventional range 1 to 3,
not derivable from a stress-product ratio. Its shipped default is 2.00 and every
run here used it. So the ocean's own timestep limit in section 3c inherits that
bracket: at the bottom of it the ocean would tolerate about twice the timestep
reported, and at the top about two thirds of it. **The ceiling is a bracket for
the same reason the forcing is**, and the published tuned values, all between
1.18 and 1.67, are Earth fits that do not transfer.

## 6b. The ocean surface velocity IS exported, and world-pt8's premise holds

OCN-3 asks whether the host exports an ocean surface velocity in a form
world-pt8 can consume, and notes that the published in-tree coupling never
needed it because it hands ice the other way. It does export one, as a
first-class coupling field rather than a diagnostic.

`genie_loop_wrappers.f90:322-323` documents `ustar_ocn` and `vstar_ocn` as
OUTPUTS of `goldstein_wrapper`, `genie_global.f90:278-279` dimensions them
`(ilon1_ocn, ilat1_ocn)`, and `goldstein.F:780-781` fills them every ocean step
as `u(1,i,j,kmax)` and `u(2,i,j,kmax)`, the top-level zonal and meridional
velocity. `gold_seaice_wrapper` already consumes them at `:297`, so the field is
not merely present, it is a field the model already couples on.

**Two things a consumer has to know.** The field is NON-DIMENSIONAL: `u` is
scaled by `usc`, which `initialise_goldstein.F:376` sets to 0.05 m/s and `:2104`
exports as `go_usc`, so the scale travels with the field but the multiplication
is the consumer's. And `usc` is one more hardcoded scale of the kind
`notes/audits/ocean-tier-implicit-earth.md` catalogues, so what it means on this
planet is that audit's question and not this one's.

---

# 7. What this did NOT establish

- **Nothing here is an adoption.** `vendor/cgenie` still has no consumer, no
  pipeline row and no step. OCN-3's remaining halves -- the Vesper-parameter
  port, and the forcing contract OCN-10 owns -- are untouched. Section 6
  answers the sea-ice-configuration and surface-velocity questions by reading
  the recipe, which is a reading and not a demonstration: no configuration
  without EMBM has been built or run here.
- **The 72 x 72 probe's circulation means nothing.** `dan_72` ships a real
  topography and real wind stress, but not the pair of advective wind-speed
  fields EMBM reads unconditionally, so those were made by replicating the
  36 x 36 fields into 2 x 2 blocks. Replication changes no arithmetic the
  timestep performs, so the wall clock and the Courant numbers are the wall
  clock and the Courant numbers of a 72 x 72 x 16 model. The ocean state it
  produces is not read and is not evidence about anything.
- **The scaling exponents are two-point estimates from ONE doubling**, between
  two different topographies. They are enough to say the growth is quadratic
  rather than linear in the timestep, and they are not enough to carry a
  five-fold extrapolation without the bracket section 5b states.
- **This is not a profile.** It says what a configuration costs, not which
  routine spends it. OCN-19 owns that, and the barotropic-versus-tracer arm here
  answers only the one question it was designed for.
- **Whether the frictional-geostrophic closure survives this planet's radius and
  gravity is untouched.** `notes/audits/ocean-tier-implicit-earth.md` findings 1
  and 4 stand exactly as they were: every number here is a number at Earth's
  hardcoded `rsc` and `const_rEarth`, on Earth topographies, and the whole sweep
  would have to be redone after those move.
- **The host was not quiet, and no wall-clock number here is clean.** Other
  agents ran throughout at one-minute loads of 8 to 42. cGENIE is one serial
  process on a 32-thread part so it never waited for a core, but it did not have
  the cache or the memory bandwidth to itself, and the per-repeat wall fits vary
  by up to a factor of two. Every one is labelled with the load it was taken
  under and marked contaminated in the artifact rather than deleted. The
  instruction counts do not have this problem and every quantitative claim above
  rests on those; where seconds appear they come from the `perf` runs' own CPU
  time and are upper bounds. world-ap7w owns the quiet-host re-take.
- **No generated world was integrated by cGENIE.** Section 3f establishes the
  connector's output at 48 x 40 and 72 x 72 structurally, against the format the
  shipped worlds use. Whether GOLDSTEIN solves the barotropic streamfunction on
  a 72 x 72 `.paths` is a separate question and is not answered here.
