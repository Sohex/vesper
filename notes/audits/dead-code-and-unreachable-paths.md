# Dead code and unreachable paths in the climate model fork

Audited 2026-08-24 against `vendor/exoplasim/` at 3ca1ab37, the five binaries in
`exoplasim/binary_manifest.json`, and the namelists the model itself wrote in
`exoplasim/runs/run_2b20e3324bb0/`.

Worldbuilding frame: this file is about the source of the Vesper project's
climate model, a hard fork of ExoPlaSim. Every quantity named here is a property
of that source or of a run of it. Nothing here is about the real world.

## The question

`notes/audits/dormant-exoplasim-modules.md` asked which compiled modules never do
any work. This asks the wider question in the form the hard fork makes
answerable: what is in the tree that NO configuration can reach, what is
compiled that no run executes, and what is reachable only through an interface
that would abort if anyone used it.

The hard fork declaration is what makes this worth doing. Upstream
compatibility stopped being a constraint on 2026-08-21, so "upstream might want
it" no longer defends anything in this tree, and the question of whether a file
should exist can be asked of every file.

## What "dead" means here, and the three instruments

Kept apart throughout, because they carry different weight.

**Unreachable by construction.** No configuration reaches it: the file is in no
build, or the guard is a `parameter` that is constant-false. Only a source edit
falsifies these.

**Unreachable by configuration.** Never executes in the configurations this
project builds and runs. Falsified by one named namelist key. Every such claim
below was checked against the model's own namelist echo in
`run_2b20e3324bb0/MOST_DIAG.00003` rather than against source defaults or
harness intent, because three keys (`NSWRCL`, `NEWRSC`, `NMOMENT`) appear in no
harness file at all and are invisible any other way.

**Unreferenced.** No caller anywhere, established by a comment-stripped sweep of
every procedure name across all 58 files in `plasim/src`, not merely the
compiled set, with `external` statements, procedure pointers and C-side symbols
checked separately. The only such reachers in the tree are
`mpimod_omp.f90:554` for the OpenMP intrinsics and `nresources` in
`pumax_stub.c`.

The compiled set is the `_sources` list in `plasim/CMakeLists.txt:147-184` with
the four variant slots resolved. All five registry binaries are `parmode: omp`,
l10, p16, at T21 through T170; MPI and serial builds exist only as the ad hoc
arms of the verification scripts.

## 1. Three defects that are not dead code

Found by this sweep, and they outrank everything below it.

### 1a. Continuation segments run on undeclared damping and the wrong filter

`exoplasim/scripts/run_exoplasim.py:1898-1900` passes `physicsfilter`,
`filterkappa` and `filterpower`. `exoplasim/scripts/continue_exoplasim.py:463`
passes only `physicsfilter`, so the other two take `configure()`'s defaults.

`config/planet.yaml:388` declares `filter_power: 16`. The staged
`plasim_namelist` of `run_2b20e3324bb0` records `NFILTEREXP = 8`.

`NDEL`, `NHDIFF`, `TDISSD`, `TDISSZ`, `TDISST` and `TDISSQ` are absent from that
namelist entirely, because `continue_exoplasim.py` never calls
`declare_hyperdiffusion` and `configure()` does not write them either. A
continuation therefore integrates on `plasim.f90`'s hard-coded T21/T42 damping
rather than the values this project derived.

`expected_namelist_keys` at `run_exoplasim.py:1462` checks neither set, which is
why `verify_staged_namelists` passes. `run_2b20e3324bb0` reached
`climatology_complete` at 85 orbits on this namelist.

### 1b. `configure(orography=)` wrote a key no namelist declared

`OROSCALE` went into `landmod_namelist` from both the `configure` and the
`modify` path, and `oroscale` appeared in no `namelist /.../` statement in any
of the 57 source files. `landmod.f90`'s read of `landmod_nl` is bare, with no
`iostat`, so an unrecognised name aborts the run rather than being ignored: the
documented API would have stopped the model, and was latent only because this
harness never passes `orography=`.

`oroscale` now lives in `planet_nl`, in every planet module, and the API writes
it to `planet_namelist`. It could not go into `landmod_nl`: every reader of it
-- `surfmod`'s scaling of `doro` and `glaciermod`'s three -- runs inside
`surface_ini`, and `planet_ini` is the only namelist read that happens before
that. `p_mars.f90` already assigned it in the planet module, so `planet_nl` is
also where it was already treated as belonging.

`NDESERT` had the same shape on one branch: the `desertplanet=True` arm wrote it
to `plasim_namelist`, which declares it, and the `desertplanet=False` arm wrote
it to `landmod_namelist`, which does not -- so turning a desert planet back off
aborted the run and left `NDESERT = 1` standing. Both arms now write
`plasim_namelist`.

### 1c. The LSG coupler stub has drifted from its call sites

`cpl_stub.f90` is compiled (`CMakeLists.txt:183`); `cpl.f90` is not.

- `cpl_stub.f90:5` declares `clsgini(kdatim,ktspd,kaomod,pslm,kxa,kya)`, six arguments.
- `cpl.f90:136` declares it with eight.
- `oceanmod.f90:345` calls it with eight: `clsgini(ndatim,ntspd,naomod,nlsg,ngui,zls,nlon,nlat)`.

The count is not the whole of it. The call passes the integers `nlsg` and `ngui`
into positions four and five, where the stub declares `pslm(kxa,kya)`, a real
array. `clsgstep` is the same: ten passed at `oceanmod.f90:466`, eight declared
at `cpl_stub.f90:19`.

These are non-module externals with no explicit interface, so nothing diagnoses
it and the link succeeds. `NLSG = 1` on this build corrupts the stack rather
than coupling an ocean.

## 2. Unreachable by construction

### 2a. Fifteen files, about 14,560 lines, are in no build configuration

`cpl.f90` 1998, `icemod_template.f90` 1934, `rainmod_bm.f90` 1424,
`rainmod_kuo_old.f90` 1324, `rainmod_mca.f90` 1015, `plasim_dummy.f90` 817,
`pumax.c` 4721, `mpimod_multi.f90` 463, `guimod.f90` 422, `p_mars.f90` 175,
`p_exo.f90` 170, `outdiag.f90` 72, `readdat.f90` 34, `resmod_def.f90` 30, and
the pair `newsnow.f90` 113 with `buildice.f90` 30.

The pair is the exception and is not dead: it is the offline glacier
accelerator, CLIM-53 owns the verdict, and
`notes/audits/dormant-exoplasim-modules.md` section 5 argues it is the only
device in the tree for reaching an equilibrium ice sheet inside an affordable
run. The CMake build has no rule for either program, so adopting it means
writing one.

Three of the others need sharpening beyond "not compiled".

**The three alternative convection schemes have no build knob at all.**
`CMakeLists.txt:152` names `rainmod.f90` as a literal. Unlike the FFT, parmode
and planet slots there is no cache variable, so selecting Betts-Miller or moist
convective adjustment is not a supported operation, and each of those files
declares its own incompatible `rainmod_nl` whose keys would abort the current
build.

**`p_mars.f90` is an accepted build value that cannot compile.**
`CMakeLists.txt:61` admits `p_earth` and `p_mars`. `p_mars.f90:13-14` has a
doubled comma inside the `planet_nl` continuation, so `-DPLASIM_PLANET=p_mars`
is a validated, documented option that fails the build.

**`p_exo.f90` is not a legal value, and it is the one that matters.**
`CMakeLists.txt:61` does not list it, and `exoplasim/scripts/build_model.py`
never passes `-DPLASIM_PLANET` at all, so the CMake cache default `p_earth`
always wins. `p_exo.f90` is the only planet module exposing `alv`, `als` and
`tmelt` in `planet_nl`: the latent heat of vaporisation, the latent heat of
sublimation, and the melting point. On this fork the module that would make the
modelled water thermodynamics settable is the one the build refuses.

**All fifteen are hashed into the binary registry.** `model_sources()` globs
`SRC/*.f90`, so editing any file in this list invalidates all five binaries'
provenance while changing nothing they contain.

### 2b. Whole trees, about 120 MB

`plasim/bld/` 34 MB of staged sources and objects from the deleted build,
already described as non-authoritative at
`exoplasim/scripts/restart_schema.py:66`; `lsg/` 25 MB, whose real
`lsgmod.f90` is 9313 lines with no build rule anywhere while the compiled
`plasim/src/lsgmod.f90` is a 4-line empty return; `plasim/dat/` 17 MB of Earth
and Mars boundary sets that the staging path never reads; and the sibling models
and utilities `puma/`, `sam/`, `cat/`, `octave/`, `cat_tools/`, `images/`,
`tools/`, `basespecs/`, `hazeconstants/`, `exoplasim_testers/`.

`postprocessor/` belongs with them, and it ships: every `burn7` call site in
`__init__.py` is commented out, yet `.gitignore:39` re-includes
`plasim/run/burn7.x` by name and `__init__.py:604` copies the whole of
`plasim/run/` into every model working directory.

The LaTeX and PDF manuals, about 66 MB, are vendored documentation and are a
separate judgement.

### 2c. The Python package: about 12,500 of 17,670 lines

The whole world repo reaches this package through six names: `exo.Earthlike`
with `configure`/`run`/`_edit_namelist`/`_add_postcodes`/`cfgpostprocessor`/`exportcfg`,
`exo.Model` as a type annotation, `exo.__file__`, `exoplasim.compile_pyfft`,
five `pyburn` entry points, and `exoplasim.makestellarspec.convert`. Dynamic
access was ruled out first: `getattr`, `setattr`, `eval`, `exec`, `globals`,
`__dict__` and `vars` return zero hits across all eleven top-level modules, and
the one string dispatch, `configure(otherargs=)`, routes into namelist keys
rather than Python attributes.

Whole modules with no caller, four of which cannot be imported at all:

- `hpc/`, 2180 lines over 11 files. Seven fail to import on Python 2 implicit relative imports, and two of the modules they import, `identity.py` and `crawldefs.py`, exist nowhere in the tree. `launcher.py:8` imports `sets`, removed in Python 3.0. There is no cluster anywhere in this project's configuration or documentation.
- `pRT.py`, 2091 lines. `petitRADTRANS` is absent from the environment and from `requirements.txt`. `__init__.py:19-21` wraps the import in a bare `try/except: pass`, so the failure is invisible.
- `surfacespecs.py`, 1471 lines and zero definitions. Its only referrer is `pRT.py`, yet `__init__.py:16` imports it unconditionally and `surfacespecs.py:1469` globs and loads every `basespecs/*.npy` at import time. Every `import exoplasim` in this project pays that to populate data only the dead bridge reads.
- `randomcontinents.py` 577, `colormatch.py` 361, `quickprocess.py` 39. The last hardcodes another planet's radius and gravity and a path under someone else's home directory.

Inside `__init__.py`, about 2150 lines are dead, split between the residue of
the deleted build system and an upstream read-and-plot API nothing calls:
`modify` 769, `loadconfig` 238, `runtobalance` 234, `image` 169, `inspect` 123,
`transit` 104, `finalize` 92, and the tidally-locked presets. In `pyburn.py`,
`advancedDataset` 838 lines is reached only when `variables` is a dict and every
call site passes a list; the non-grid branches of `_transformvar` and
`_transformvectorvar` total 486 lines and every consumer passes `mode="grid"`;
the csv, hdf5 and npz writers are selected only by output types
`config/planet.yaml:476` does not use. In `gcmt.py`, about 1671 of 1947 lines
are unreachable from the one live entry point, `load`, and `gcmt.py:1740-1947`
is 208 lines of which 191 are commented out.

### 2d. Build-system residue

`docs/src/reference/vendored-upstreams.md` records the CMake rework as the
fork's largest deliberate divergence and lists what it deleted. Two tracked
files that belong on that list survived it: `make_most`, which builds `most.x`
from `most.c`, and `most.c` itself, 125 KB of X11 launcher. Their only remaining
reference in the tree is a line of vendored LaTeX at `Puma_UG_17/install.tex:58`
reading "Used by configure script", and that script is gone.

`sysconfigure()` at `__init__.py:139-214` is now a `chdir` round trip and
nothing else; every real statement is commented out. Its `except Exception as e`
handler prints `result.stderr`, and `result` is never bound in the function, so
the handler raises `NameError` over whatever it caught. `Model.__init__:464`
still calls it and then touches a `firstrun` marker into the installed package
directory. The marker is present and untracked here, so the block is dormant in
this tree and fires on a fresh clone.

`config/planet.yaml:289` states that `rebuild_binaries.py` writes line 3 of
`most_compiler_mpi` before every build. That script contains no reference to
`most_compiler`; the flag line reaches CMake through `build_model.py`.

`vendor/exoplasim/.gitignore` explains its every rule in terms of what
`configure.sh` and `compile.sh` write into the source tree, and names
`most_compiler*`, `firstrun`, `F90_INTEGER`, `F90_REAL` and `most_info.txt`,
which are untracked leftovers sitting in the tree from before the rework.

`setup.py` never runs, because `pyproject.toml` declares hatchling as the build
backend, and the two disagree on both version and dependencies. `MANIFEST.in` is
a setuptools input hatchling ignores. `pyproject.toml:36` reads
`[projection.optional-dependencies]`, a typo for `project`, so hatchling
silently ignores the whole extras table and `pip install exoplasim[pRT]`
installs nothing.

## 3. Unreachable by configuration

About 9000 lines of compiled physics never execute in a production run. The
module-level cases are already inventoried in
`notes/audits/dormant-exoplasim-modules.md` section 1 and are not re-argued.
What that audit does not carry:

**Entropy diagnostics, about 452 lines across nine modules.** Three separately
declared switches named `nentropy`, in `plasim_nl` (`plasim.f90:1412`),
`icemod_nl` (`icemod.f90:288`) and `oceanmod_nl` (`oceanmod.f90:160`), plus
`nentro3d`. Nothing writes any of them, and none appears in any run namelist.
381 lines over 42 sites are gated on `nentropy > 0` and 71 lines over 31 sites
on `nentro3d > 0`. This is the largest namelist-gated block in the model.

**`tpcore.f90` and `trc_routines.f90`, 2221 lines, are dead twice over.** The
outer gate is `NQSPEC = 1`, which advects moisture spectrally: `tracer_main`'s
only call is `plasim.f90:3914`, under `nqspec == 0`. The inner one needs no
namelist at all. `iord`, `jord`, `kord`, `fill`, `mfct`, `deform` and `cnst` are
`parameter` declarations in `tracermod.f90:20-31` and `aeromod.f90:21-37`, so
about 780 further lines cannot run at any setting: `FCT3D` 276 lines behind
`mfct=.FALSE.`, `qckxyz` with `filns` and `filew` 257 lines behind
`fill=.FALSE.`, `fxppm` 46 behind `iord=2`, both `lmtppm` arms behind an
always-zero `LMT`. `tracermod`'s own initialisation does run, so the model
builds a tracer grid nothing then advects.

**Two terms are multiplied by zero rather than branched around, and are
therefore computed on every step.** `newrsc` defaults 0 (`radmod.f90:902`), so
`radmod.f90:2553-2554` evaluates a full layer-resolved Rayleigh expression and
multiplies it by zero for every level and every column, while `:2763-2767`
carries the formulation that runs. The eight `*iaeron` aerosol products at
`radmod.f90:2788` and after behave the same way.

**`vdiffo` cannot run at any namelist setting.** `oceanmod.f90:15` is
`parameter(NLEV_OCE = 1)`, so its only call at `:926` sits behind a
constant-false `if(NLEV_OCE > 1)`, and inside the routine `nlem_oce` is 0, which
makes both its loops zero-trip. The multi-layer ocean vertical diffusion solver
is unreachable by construction.

**The GUI costs real work on every diagnostic step.** `guimod_stub.f90` is
always the compiled variant, so all 23 `call gui*` sites reach empty returns.
The argument preparation is not behind `ngui`. `subroutine energy`
(`plasim.f90:2395-2412`) exists only to build `ziso(6)` for `guiput`, `diag`
calls it unconditionally at `:2329`, and its five inputs `umax`, `t2mean`,
`precip`, `evap` and `olr` have no other consumer in the model: producing them
costs five `mpgagp` full-globe gathers at `plasim.f90:2821-2829`.

**`nprint` gates 902 lines** of instrumentation across six modules, and reads
back 0. This is diagnostic code designed to be off; it is recorded as a size,
not as a defect.

## 4. The transform migration is less complete than the tracker records

CLIM-64 closed on the finding that `mkdheat` was the last hot legmod caller. It
is not.

`mkdheat` calls legmod's `sp2fl` at `plasim.f90:5092` and `:5106`, inside
`if(nenergy > 0)` and outside every `nshtns` branch. The comment above that
block describes it as off by default. `config/planet.yaml:502` declares
`energy_diagnostics: true`, `run_exoplasim.py` maps that to `NENERGY = 1`, and
`run_2b20e3324bb0/plasim_namelist:3` records `NENERGY = 1`. The block runs every
timestep of every production run, and legmod's Legendre loops run with it,
alongside `dv2uv` at `plasim.f90:3578` and `:3587` and the `sp2fc` that `sp2fl`
calls.

The rest of legmod is correctly retained and must not be read as dead. The axis
is `nshtns`, a runtime namelist key, not the parmode: `verify_shtns_model.sh`
builds one `--parmode omp` binary and runs `NSHTNS=1` against `NSHTNS=0` on it,
so the threaded build carries both transforms and switches at runtime.
`invlega`, `invlegd`, `sp2fcdmu`, `mktend`, `qtend` and `uv2dv` are the
reference arm. Only `sp3fc` (`legmod.f90:465-476`, 12 lines) is unreachable in
every configuration.

The mirror holds on the other side: all four `use shtnsmod` sites that import
the wrappers are inside `#ifdef OMPSHARED`, so about 398 lines of
`shtnsmod.f90` are compiled with zero references in any MPI or serial build.
That is exactly why `CMakeLists.txt` has to link SHTns in every configuration.
`shtns_teardown` (`shtnsmod.f90:563-573`) has no caller in any configuration;
the process exits holding its plans.

`dealias_gp` (`plasim.f90:2876-2945`) and the legmod pair it alone calls,
`fc2sp_t` and `sp2fc_t`, total 128 lines reachable only at `NSHTNS=0` with
`NDEALIAS=1`. `plasim.f90:694-697` stops the model outright if `ndealias > 0`
and `nshtns == 1`, `ndealias` defaults 0, `config/planet.yaml` has no key for
it, and commit 4f08c119 records the remedy as refuted.

## 5. Namelist keys that gate nothing

Of 368 namelist variables across the 16 groups the compiled model declares, 167
have no writer anywhere in this project. Most of those are ordinary unused
knobs. These are different: they are declared, documented, sometimes written,
and read by no code at all.

**`nfixer` was the one that mattered, and it is gone.** It was declared in
`miscmod_nl`, defaulted to 1 and was commented "switch for negative humidity fix
(1/0 : on/off)", and it was tested nowhere: `miscstep` calls `fixer`
unconditionally, so setting NFIXER=0 was silently ignored. Not fork damage --
the independently vendored PlaSim under
`vendor/cgenie/genie-plasim/src/fortran/miscmod.f90` has the same unconditional
call. The switch was removed rather than wired up, because negative specific
humidity is not a state this model may integrate and "off" was never a
configuration.

What the fixer does is now recorded where it runs: it borrows moisture between
columns, first within a column, then along a latitude row, then globally, and
the global step conserves the column integral except in the branch where the
global deficit exceeds the global surplus, where it takes the field to zero and
is a sink. So a global moisture or latent-energy closure carries no term for it
and a regional one does. `exoplasim/scripts/close_state_energy.py` names it
against its vapour reservoir.

**`aeroqlw` is fork-added and has no writer.** Declared `radmod.f90:906`,
defined `:305`, default 0.0, and described by the fork's own comment at
`:286-292` as the thermal-IR absorption ratio for the interactive aerosol, the
twin of `dustqlw`. `dustqlw` is written from the dust provenance file at
`run_exoplasim.py:1143`; `aeroqlw` is written by nothing. The emitted-dust path
therefore has no longwave term available, which is the remaining live half of
why `run_exoplasim.py:1242` pins `l_aerorad = 0`. The other half of that pin's
stated rationale, that `aero_ini` never populates `apart` from the namelist, was
fixed at `aeromod.f90:277` and the comment has not caught up.

Read by nothing, and now retired: `nflux`, `npackgp`, `npacksp`, `nsurf`,
`zeta` (whose own declaration said "not used", and which the shipped
`carbonmod_namelist` template wrote into every run) and `ngpitrigger` are gone
from their groups, their declarations, their broadcasts and the shipped
templates. Of the seven, only `nguidbg` and `zeta` were ever staged: `npackgp`,
`npacksp` and `nsurf` appear in no template and no run, contrary to the earlier
reading of this row.

`nguidbg` and `sellon` are out of `plasim_nl` and out of the template but keep
their module declarations, because their only reader is the uncompiled
`guimod.f90` and the compiled build passes `sellon` to bodyless column routines
in `guimod_stub.f90`. `nguidbg` also had a second declaration inside the private
`plasim_nl` that `carbonmod`'s `psurfupdate` WRITES; that group emits a new
`plasim_namelist`, and `plasim_nl` no longer declares the name, so it had to go
from there too or the next read would abort.

## 6. What this makes of the earlier audit

`notes/audits/dormant-exoplasim-modules.md` was written on 2026-08-21, one day
before commit 137c55a2 replaced the build system. Its mechanism claims all
survive: which module is dormant, and which switch gates it, still read
correctly against `plasim/CMakeLists.txt`.

Its build-side citations do not. `compile.sh:305` and `configure.sh:99-108` and
`:307-310` name files the rework deleted. One claim is inverted rather than
merely stale: its section 6 lists `mpimod_omp.f90` and `utilities_stub` among
what is "not compiled at all", when `mpimod_omp.f90` is the parallel layer of
every binary in the registry and the stub pair is the serial reference the
threaded build is checked against.

That section has been rewritten against the CMake source list.

## 7. What is deliberately retained, so it is not re-derived

The MPI build is the correctness reference the threaded build is compared
against, at bit-identity at T21 and at rounding scale at T170, and
`docs/src/reference/vendored-upstreams.md` records that it stays until the
serial build replaces it in that role. `mpimod.f90`, `utilities.f90`,
`mpimod_stub.f90` and `utilities_stub.f90` are not dead code and nothing in this
audit argues otherwise.

`fft991mod.f90` is reachable: `build_model.py:57-58` selects it for T63 and
T106. Neither rung is in the registry matrix and neither is the declared
resolution, so no binary contains it today, but a `--res T63` build compiles it.

Within the omp parallel layer, 205 lines have no caller in the compiled set:
`mpgallgp`, `mpsumr`, `mpgathersp`, the four `mpread*`/`mpwrite*` routines
superseded by `restartmod`, `mrsum`, `mrbci`, `mpgadn`, and four `finish*`
variants. `mpsumscp`, `mpgallspp` and `mpgallsp` all survived the SHTns
migration and are live; this is not transform fallout.

## 8. Elsewhere, unreferenced and small

`earthveg` (`specblock.f90:1426-1747`) is 322 lines of spectral data referenced
by no file in the tree; the other thirteen arrays in that module are read by
`radmod.f90`. `hcadencediag` (`outmod.f90:2320-2427`) has one call site and it
is commented out at `plasim.f90:807-809`. `hcadencesc` (`outmod.f90:2178-2201`)
has none at all, which is a gap rather than only waste: the ordinary and
snapshot output paths each call a scalar writer beside the gridpoint one, and
the high-cadence path does not, so that stream writes no scalar record.
`outdiag.f90` is an orphan file duplicating the live `outdiag` at
`outmod.f90:1096`, and adding it to the build would be a duplicate-symbol error.

Zero-caller procedures across the compiled set, about 485 lines in total:
`fyppm`, `cosa`, `cosc`, `n2mmr`, `orb_print`, `mkicecf`, `surfname`,
`check_equality`, `dtodat`, `nweekday`, `ndayofyear`, `writecolumn`,
`guistep_puma`, `guihorlsg`, `readvariablecode`, `_getknownwordlength`,
`_getmarkerlength`, `f2py_compile`.

`zsolars` is written to the restart at `radmod.f90:1799` and its only read is
commented out at `:1144`. That is deliberate residue of the CLIM-38 repair
recorded in `notes/audits/zsolars-restart-overread.md`, but the record still
carries a `REQUIRE_EQUAL` policy at `restart_schema.py:311`, so a converter can
refuse a restart over two values nothing consumes.
