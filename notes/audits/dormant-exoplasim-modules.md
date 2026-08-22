# What the climate model ships and this configuration never runs

Audited 2026-08-21, against the namelists ExoPlaSim itself wrote in
`exoplasim/runs/run_066cd40559fc/` rather than against the scripts that intend
them, and against the build path in `vendor/exoplasim/exoplasim/compile.sh`.

Worldbuilding frame: this file is about which parts of the Vesper simulation's
climate model source are compiled into the binary and never execute. Nothing
here is about the real world.

## The question

Rule 7 asks what a change makes worthless. This asks the neighbouring question,
which nothing had asked: of the modules `plasim/src` compiles into the
executable this project runs, which ones never do any work, and which of those
absences are decisions somebody made rather than defaults nobody has read.

Nine modules are dormant. Three are settled decisions, two are already owned by
open rows, one is a diagnostic that is cheaper to enable than anybody has
assumed, and one is a facility this project is half-using without knowing it.

## 1. The inventory, and the switch that gates each one

| module | gate, as the run namelists carry it | where the value comes from |
| --- | --- | --- |
| `simba.f90` | `NVEG = 0` | `configure(vegetation=)`, never passed, defaults False |
| `carbonmod.f90` | `NCARBON = 0`, `NCO2EVOLVE = 0` | `configure(co2weathering=, evolveco2=)`, both default False |
| `hurricanemod.f90` | `NSTORMDIAG = 0`, `HC_CAPTURE = 0` | `configure(stormclim=)`, never passed |
| `aeromod.f90`, `aerocore.f90` | `L_AERO = 0` | `model.dust_emission` absent from `config/planet.yaml` |
| `tpcore.f90`, `trc_routines.f90` | `NQSPEC = 1` | model default; nothing writes it |
| `lsgmod.f90` | `NLSG = 0`, and `OCEANCOUP=cpl_stub` at compile time | `compile.sh:305` |
| `guimod.f90`, `pumax.c` | `NGUI = 0`, and both are `_stub` variants | `configure.sh:307-310` |

Sub-features off inside modules that do run: `NHDIFF = 0` (deliberate, CLIM-16),
`NFLUKO = 0`, `NDUSTRAD = 0`, `NENERGY` and `NENER3D` unset, `NHCADENCE = 0`,
`newsurf = 0` in `landmod`, and every one of `ndiaggp*`, `ndiagsp*` and `ndiagcf`
at 0, which is why `outmod`'s `outdiag` never fires.

The aerosol pair is not re-argued here. `L_AERO` and `NDUSTRAD` are owned by
DUST-11, DUST-13 and DUST-15, in that order, and this audit adds nothing to
them.

## 2. Three of the dormant modules are settled and need no row

**`simba.f90`.** Vegetation reaches the model as prescribed surface fields,
through `land_albedo_source: vegetated` and the `.sra` files
`build_surface_albedo.py` writes. Turning `NVEG` on would put a second and much
poorer vegetation model inside the climate model, competing with the one that is
already the downstream component. The albedo and roughness feedbacks SIMBA
provides are the ones this project derives explicitly and can cite.

**`carbonmod.f90`.** `pCO2_bar` is prescribed, and the block in
`config/planet.yaml` says why. Worth recording for whoever reaches for the
switch: the namelist it would activate carries `KACT = 0.09`, `KRUN = 0.045`,
`PEARTH = 79.0` and `WMAX = 80.9`, which are Earth weathering constants, so
enabling it is a PHYS-class inheritance and not a free mechanism. The weathering
this project does believe in is `pedology/scripts/weathering_fluxes.py`, offline
and citable.

**`lsgmod.f90`.** Doubly off, and the second half is the one that matters:
`compile.sh:305` exports `OCEANCOUP=cpl_stub`, so `cpl.f90` is never compiled
and every coupling call in `oceanmod` reaches an empty return. `lsgmod.o` links
as dead weight. The package also ships the full LSG source in `lsg/src/`, which
nothing builds. This is recorded so that nobody rediscovers a 3-D ocean sitting
in the tree and reads it as an available option: OCN-3 already took the resolved
tier off the board on spin-up cost, and a coarse geostrophic ocean does not
escape that argument.

## 3. `hurricanemod` is two instruments, and only one of them is resolution-bound

The module resolves nothing and feeds nothing back. Its index half computes
Emanuel potential intensity (`POTI`), an entropy deficit, a ventilation index,
lower-atmosphere absolute vorticity and a genesis potential index, then masks
cells against thresholds taken from Komacek et al. (2020): `GPITHRESH = 0.37`,
`VITHRESH = 0.145`, `VMXTHRESH = 33.0`, `LAVTHRESH = 1.2e-5`,
`MINSURFTEMP = 298.15`. Those are ENVIRONMENTAL indices. They exist precisely
because a model of this class cannot resolve a storm, and evaluating them on
coarse fields is their intended use rather than a compromise.

The module is more planet-aware than most of what this project has audited:
`hurricaneini` takes `RD = gascon` from the model rather than assuming Earth
air, and `absvorticity` uses the model's own `omega`. Its thermodynamic
constants are Earth's, which is right for a 1 bar N2/O2 atmosphere, and every
threshold above is a namelist key rather than a compiled constant.

The capture half, `HC_CAPTURE` with the high-cadence stream, is what the
resolution intuition is about. It hunts a RESOLVED vortex: at least 33 m/s in
the lower atmosphere with 20.5 m/s at the surface, over at least 30 gridcells,
sustained at least 256 timesteps. The bar on this world is higher than on Earth,
not lower, and by two compounding factors.

Grid spacing at each rung of the ladder, from `config/planet.yaml`'s radius:

| | Vesper | Earth |
| --- | --- | --- |
| T42 | 375 km | 313 km |
| T85 | 188 km | 156 km |
| T127 | 125 km | 104 km |
| T170 | 94 km | 78 km |

A given spectral truncation buys 20% coarser spacing in kilometres, because the
planet is 1.2 radii. And the storms are not correspondingly larger: the 30-hour
rotation puts the rotation rate at 0.800 of Earth's, which raises the
deformation radius, while 12.81 m/s2 puts the scale height at 0.766 of Earth's,
which lowers it, for a net 0.957. Earth GCMs want roughly 50 km before
cyclone-like vortices reach realistic intensity, which on this planet is about
T320. **The capture half is off the resolution ladder this project has, and
converting up does not put it on.**

**The index half is available at T42 for one namelist key and has never been
asked for.** There is no prognostic hurricane state, so it is valid on a
continuation segment as well as a cold start, unlike every surface field, which
`landmod`'s `landini` takes from the restart. The decision procedure, if it is
ever wanted: it is a diagnostic, so it goes on a segment declared `diagnostic`,
`stormclim=True` has to be reapplied per segment for the CLIM-17 reason, and
what it returns is whether this climate's environments admit cyclogenesis at
all, not how many storms there are.

Cost when off is not quite zero and is near enough: `hurricanestep` is called
unconditionally at `plasim.f90:3485` and returns after zeroing seven `NHOR`
masks and computing one surface wind speed.

## 4. The hurricane accumulators are restart records and belong to the reset set

`plasim.f90:958` writes eight accumulated index fields into the restart --
`agpi`, `aventi`, `alaav`, `ampoti`, `avrmpi`, `acapen`, `alnb`, `achim` -- and
`outmod.f90:2490` zeroes all eight at every output interval boundary. That is
exactly the criterion `reset_restart_accumulators.py` states for a record
belonging to its set, and they were not in it. They are present in real
restarts: records 76 to 83 of `MOST_REST.00002` in the run cited at the top.

They are in the set as of this audit. Nothing on disk is affected, because the
records only ever hold zero while `NSTORMDIAG = 0`; the point of the fix is that
the first seeded run with the diagnostics on would otherwise have opened
mid-window on somebody else's partial sum, which is the CLIM-31 failure the
script exists to prevent.

## 5. The glacier accelerator is half-wired, and this project has been feeding it

Glaciers are ON: `glaciermod.f90` is compiled, `NGLACIER = 1`, `GLACELIM = 2.0`,
`ICESHEETH = -1.0`, from `surface.glaciers.enabled`. What is not wired is the
offline half of the facility.

`configure.sh:99-108` builds two standalone programs from `plasim/src`:

- `buildice.f90`, 29 lines. Reads `newdsnow` and `lsm`, adds 400 m of ice to
  every land cell capped at 3 km, writes `newdsnow` back. A one-shot ice-sheet
  initialiser.
- `newsnow.f90`, 112 lines. Reads five consecutive years of snow tendency plus
  `newdsnow`, `newxsnow`, `lsm` and `persist`, and extrapolates the five-year
  mean forward by a `deltyrs` argument, but only where snow already exists and is
  either shrinking or has persisted a full year. Land is capped at 3 km, sea ice
  at 10 m.

`newsnow` is a TIME ACCELERATOR for land ice, the same device `NCVEG` is for
vegetation carbon: it lets an ice sheet reach equilibrium in a few hundred
orbits instead of the tens of thousands the model would have to integrate.

**`glaciermod` already writes its inputs.** `glaciermod.f90:610-619` calls
`finishup` on `newdsnow`, `newxsnow`, `lsm` and `persist` at the end of every
model call, and all four are sitting in the run directory cited at the top. So
every glacier-enabled run this project has made has produced the accelerator's
feedstock and nothing has ever consumed it. Nothing in ExoPlaSim's own Python
invokes either program either -- only `configure.sh`, which builds them -- and
`newsnow`, `buildice` and `newdsnow` appear nowhere in this repository outside
`vendor/`.

Two things stand between the facility and use, and neither is fatal. Both
programs hardcode `NLAT = 32, NLON = 64`, so they are T21 and would silently
misread a T42 field. And `buildice`'s uniform 400 m over land is an
initialisation for a planet whose ice sheets are already known, which this one's
are not.

This matters because it is upstream's answer to the question
`notes/glacier-rough-pass.md` leaves open. That note refutes latitude control
and defers the real treatment to sub-grid machinery that does not exist; it is
also true that 3 km of ice does not accumulate inside an affordable run, and
this is the device that exists for exactly that. CLIM-53 owns the verdict.

## 6. What is not compiled at all

Recorded so it is not re-derived. Selected against at build time:
`fft991mod.f90` (only T63 and T106, or an explicit `-f`, choose it),
`rainmod_bm.f90`, `rainmod_mca.f90` and `rainmod_kuo_old.f90` (`RAINMOD=rainmod`
selects Kuo), `p_mars.f90` and `p_exo.f90` (`PLAMOD=p_earth`; `p_exo.f90` is
referenced by no build script in the tree and its own header says it is the
Earth module), `mpimod_stub.f90`, `mpimod_omp.f90`, `mpimod_multi.f90`,
`utilities_stub`, `cpl.f90`, `guimod.f90`, `pumax.c`. Dead in the tree rather
than deselected: `icemod_template.f90`, `plasim_dummy.f90`, and `outdiag.f90`,
which duplicates the subroutine that lives at `outmod.f90:1070`.

The package's `glacier/` directory is not a model and not related to
`glaciermod.f90`. It is a PDF, a readme and `N032_surf_0129.sra`, an Earth T21
surface file with Antarctica and Greenland lowered a few kilometres so that
GPlaSim can grow ice back onto them. There is nothing in it for this world.
