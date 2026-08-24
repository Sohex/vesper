# What the climate model ships and this configuration never runs

Audited 2026-08-21, against the namelists ExoPlaSim itself wrote in
`exoplasim/runs/run_066cd40559fc/` rather than against the scripts that intend
them, and against the build declaration in
`vendor/exoplasim/exoplasim/plasim/CMakeLists.txt`.

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

The inventory below is what the run cited above carried, and the gates are
unchanged except where a later finding is named. What did move is what a run
now DECLARES about them: `NENERGY`, `NENER3D`, `NENERGYFIX`, `L_AERO` and
`NDUSTRAD` are written by `run_exoplasim.py` and checked by
`expected_namelist_keys`, so a segment cannot silently drop back to a compiled
default, and the dormant Earth constants inside each module say so at the switch
under world-9d1.

## 1. The inventory, and the switch that gates each one

| module | gate, as the run namelists carry it | where the value comes from |
| --- | --- | --- |
| `simba.f90` | `NVEG = 0` | `configure(vegetation=)`, never passed, defaults False |
| `carbonmod.f90` | `NCARBON = 0`, `NCO2EVOLVE = 0` | `configure(co2weathering=, evolveco2=)`, both default False |
| `hurricanemod.f90` | `NSTORMDIAG = 0`, `HC_CAPTURE = 0` | `configure(stormclim=)`, never passed |
| `aeromod.f90`, `aerocore.f90` | `L_AERO = 0` | `model.dust_emission` absent from `config/planet.yaml` |
| `tpcore.f90`, `trc_routines.f90` | `NQSPEC = 1` | model default; nothing writes it |
| `lsgmod.f90` | `NLSG = 0`, and `cpl_stub.f90` is compiled in its place | `CMakeLists.txt:184` |
| `guimod.f90`, `pumax.c` | `NGUI = 0`, and both are `_stub` variants | `CMakeLists.txt:152,182` |

Sub-features off inside modules that do run: `NHDIFF = 0` (deliberate, CLIM-16),
`NFLUKO = 0`, `NDUSTRAD = 0`, `NENERGY` and `NENER3D` unset, `NHCADENCE = 0`,
`newsurf = 0` in `landmod`, and every one of `ndiaggp*`, `ndiagsp*` and `ndiagcf`
at 0, which is why `outmod`'s `outdiag` never fired. `outdiag.f90` has since been
deleted under world-cmz: it duplicated the subroutine in `outmod.f90` and no
configuration compiled it.

`NFLUKO` is no longer merely off. world-4ba made it refuse: with no code-169
climatology supplied, `xclsst` holds its `-999` sentinel, and relaxing this
world's SST toward a sentinel is what `nfluko` would have done, so `icemod` and
`oceanmod` abort by name rather than run it.

The aerosol pair is not re-argued here. `L_AERO` and `NDUSTRAD` are owned by
DUST-11, DUST-13 and DUST-15, in that order, and this audit adds nothing to
them. Two later findings do touch them and are recorded so they are not
re-derived: world-906 found that `APART` and `RHOP`, the dust grain radius and
density in `aero_nl`, were never written by anything, so in-model dust would
have settled about 1790 times too slowly; and world-24v gave `aeroqlw` a writer
and turned the `l_aerorad` pin into a switch.

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
`CMakeLists.txt:184` names `cpl_stub.f90`, so every coupling call in `oceanmod`
reaches an empty return, and `plasim/src/lsgmod.f90` is a four-line stub that
links as dead weight. The package used to ship the full LSG source in `lsg/src/`
and the real `cpl.f90` beside it; world-6ak and world-cmz deleted both as trees
nothing builds, reads or ships. That does not change the verdict, it records it:
nobody can now rediscover a 3-D ocean sitting in the tree and read it as an
available option. OCN-3 already took the resolved tier off the board on spin-up
cost, a coarse geostrophic ocean does not escape that argument, and
`vendor/cgenie` is the candidate this project actually carries for that role.
`oceanmod.f90:437` also refuses `nlsg > 0` unless `n_days_per_year` is exactly
360, which is a count of SIDEREAL DAYS PER ORBIT and is this planet's 146, not a
setting.

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
converting up does not put it on.** CLIM-55 owns that as a declared gap.

**The index half is available at T42 for one namelist key and has never been
asked for.** There is no prognostic hurricane state, so it is valid on a
continuation segment as well as a cold start, unlike every surface field, which
`landmod`'s `landini` takes from the restart. It is not quite free: the eight
indices are output codes 322 to 329, which are not in `REGULAR_CODES`, and
`stormclim=True` has to be reapplied per segment for the CLIM-17 reason. What it
returns is whether this climate's environments admit cyclogenesis at all, not
how many storms there are. CLIM-54 owns the decision.

Cost when off is not quite zero and is near enough: `hurricanestep` is called
unconditionally from `plasim`'s main loop and returns after zeroing seven `NHOR`
masks and computing one surface wind speed.

Two later findings sit on this module and neither changes the verdict above.
world-khn declared the four thresholds that are in GRID units rather than
physical ones, `SIZETHRESH` among them, at `hurricanemod.f90:70-82`; world-9d1
put the trigger and the consequence there too. The positive form both settled
on: this diagnostic is not meaningful on this world as parameterised, it stays
off, and if the index half is ever wanted its thresholds are a re-derivation and
not a namelist tweak. `LAVTHRESH = 1.2e-5` per s is an absolute vorticity
threshold against a planetary vorticity 0.8 of Earth's at a 30-hour day, so 25
per cent too strict, and `CPD = 1005.7` at `:51` is hardcoded rather than tied
to `acpd`.

## 4. The hurricane accumulators are restart records and belong to the reset set

`epilog` writes eight accumulated index fields into the restart --
`agpi`, `aventi`, `alaav`, `ampoti`, `avrmpi`, `acapen`, `alnb`, `achim` -- and
`outmod` zeroes all eight at every output interval boundary. That is
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

The deleted `configure.sh` built two standalone programs from `plasim/src`,
and the CMake build that replaced it has no rule for either:

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

**`glaciermod` already writes its inputs.** `glacierstop`
(`glaciermod.f90:702-709`) calls `finishup` on `newdsnow`, `newxsnow`, `lsm` and
`persist` at the end of every model call, and all four are sitting in the run
directory cited at the top. So
every glacier-enabled run this project has made has produced the accelerator's
feedstock and nothing has ever consumed it. Nothing in ExoPlaSim's own Python
invokes either program either -- only the deleted `configure.sh` ever built
them, and nothing builds them now -- and
`newsnow`, `buildice` and `newdsnow` appear nowhere in this repository outside
`vendor/`.

Two things stood between the facility and use, and neither was fatal. Both
programs hardcoded `NLAT = 32, NLON = 64`, so they were T21 and would silently
misread a T42 field; **world-zsa fixed that**, and both now take the grid from
the generated `resmod.f90` like every other unit in the tree
(`newsnow.f90:4`, `buildice.f90:4`, `:21-23`), so a mismatch is a build error.
What remains is `buildice`'s uniform 400 m over land, which is an initialisation
for a planet whose ice sheets are already known, and this one's are not.

This matters because it is upstream's answer to the question
`notes/glacier-rough-pass.md` leaves open. That note refutes latitude control
and defers the real treatment to sub-grid machinery that does not exist; it is
also true that 3 km of ice does not accumulate inside an affordable run, and
this is the device that exists for exactly that. CLIM-53 owns the verdict.

## 6. What is not compiled at all

Recorded so it is not re-derived, and stated against the `_sources` list in
`plasim/CMakeLists.txt:148-185`, which is the whole of the build declaration.

**The parallel layer was a variant slot, and is not one any more.**
`CMakeLists.txt` selected `mpimod_omp` with `utilities_omp` for `PARMODE=omp`,
`mpimod` with `utilities` for `mpi`, and `mpimod_stub` with `utilities_stub` for
`serial`, and all three were built: omp was what every binary in
`exoplasim/binary_manifest.json` was, and the other two were the references the
threaded build was checked against. **world-38b removed the MPI and serial build
paths and the parmode axis with them.** `mpimod_omp.f90` and `utilities_omp.f90`
are literals in `_sources` (`CMakeLists.txt:145-158`), `mpimod.f90`,
`mpimod_stub.f90`, `mpimod_multi.f90`, `utilities.f90` and `utilities_stub.f90`
are gone from the tree, and there is one decomposition. `NPRO` is the OpenMP
thread count.

Selected against by a variant slot: `fft991mod.f90`, which only T63 and T106
choose and neither is in the registry matrix. `p_mars.f90` was the other
example; `PLASIM_PLANET` and both alternative planet modules went under
world-58v and world-cmz, so `p_earth.f90` is a literal too.

Compiled by no configuration and reachable by no flag, because the build named
them as literals with no variant slot at all: `rainmod_bm.f90`,
`rainmod_mca.f90` and `rainmod_kuo_old.f90` (`CMakeLists.txt` hardcodes
`rainmod.f90`), `cpl.f90`, `guimod.f90`, `pumax.c`, and `p_exo.f90`, which
`CMakeLists.txt` did not admit as a legal value.

Dead in the tree rather than deselected: `icemod_template.f90`,
`plasim_dummy.f90`, `readdat.f90`, `resmod_def.f90`, and `outdiag.f90`, which
duplicated the subroutine in `outmod.f90`.

All of the above were deleted: the thirteen sources under world-cmz, the
unreferenced procedures inside the surviving modules under world-ro6, and
`mpimod_multi.f90` with the rest of the MPI path under world-38b. What survives
as a `_stub` and is compiled on purpose: `guimod_stub.f90`, `pumax_stub.c`,
`cpl_stub.f90` and the four-line `lsgmod.f90`. `newsnow.f90` and `buildice.f90`
survive uncompiled and on purpose, for `clim-53`.

`notes/audits/dead-code-and-unreachable-paths.md` carries the line counts, the
two build options that were accepted and could not compile, and the consequence
that every file in this list was hashed into the binary registry -- which
world-cmz also fixed, by reading the compiled set out of `CMakeLists.txt`
instead of globbing `SRC/*.f90`, so keeping a file on purpose costs no
provenance.

The package's `glacier/` directory was not a model and not related to
`glaciermod.f90`. It was a PDF, a readme and `N032_surf_0129.sra`, an Earth T21
surface file with Antarctica and Greenland lowered a few kilometres so that
GPlaSim can grow ice back onto them. There was nothing in it for this world, and
world-6ak deleted it.
