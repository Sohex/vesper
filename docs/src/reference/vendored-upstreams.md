# The vendored upstreams

Four of this project's components are other people's code, vendored as git
subtrees so that a change to the component and the change to whatever consumes
it land together, and so provenance is a commit in this repository rather than
the state of a directory outside it. Each is a maintained fork rather than a
pinned dependency; `docs/src/reference/design-intent.md` says why that is the
expected end state for an external model here and how to weigh a candidate
against it.

## World Orogen

The geography comes from a personal fork of World Orogen, vendored into this
repo at `vendor/orogen/` as a git subtree from the `cf-fork` branch of
`raguilar011095/planet_heightmap_generation`. Pull upstream with `git subtree pull --prefix vendor/orogen
orogen-fork cf-fork --squash`. Its generated output stays untracked: the
subtree's own `.gitignore` excludes `out/`, which runs to GB.

Its `tools/README.md` is the authoritative reference for the export format -- read it
before writing anything that consumes `source/`. The fork adds, over upstream:
lithology (rock class, erodibility, scarp potential), preserved endorheic basins,
a richer export manifest, non-Earth planet parameters, and direct emission onto
Gaussian (spectral) grids.

## ExoPlaSim

### It is a HARD fork as of 2026-08-21

**Upstream compatibility is no longer a constraint on this component.** The
threaded build is the line this model develops along, and the changes it wants
next -- full-globe grid arrays for SHTns above all -- cannot be made while the
MPI decomposition also has to keep working. Declared deliberately, so that the
next person does not preserve a contract nobody is holding.

What that licenses: removing the MPI path, restructuring the parallelism,
changing array shapes and storage classes across the physics, and taking any
change that is faster without asking whether it could be contributed back.
Pulling upstream is not expected to be possible again, and a change worth
sending upstream now has to be written for upstream separately.

**The MPI path is gone, and with it the parmode axis.** `mpimod.f90`,
`mpimod_multi.f90` and `mpimod_stub.f90` are deleted, `plasim/CMakeLists.txt`
names `mpimod_omp.f90` and `utilities_omp.f90` as literals, and
`build_model.py` has no `--parmode`: the model has one parallel layer, threads
over a shared address space. Removing it is not a speedup -- the threaded build
never linked it -- and what it costs is the second reference. Every correctness
check in this component used to compare a threaded arm against an MPI one, and
that comparison caught the dv2uv planetary vorticity race, the weight
pre-scaling bug and the mkdheat inertness result. **A gate that needs an
independent reference now has to build one**, from a control patch or from a
standalone driver, rather than reaching for a second registered runtime.

The executable's name follows: `most_plasim_<res>_l<levels>_p<ranks>.x` carries
no parallel-mode suffix, because there is one parallel mode and an axis with one
value separates no two builds.

### Names upstream's Python loads and never binds

Worth reporting to `Sohex/ExoPlaSim` and to `alphaparrot/ExoPlaSim` under it,
because none of it is fork-introduced: every one of these is present in the
subtree-import commit and so has never executed upstream either. They are
repaired here.

`pyburn.py`, five names in the two transform helpers. `_transformvar` reads
`vairable` for its own `variable` argument, and `nlats` for `nlat`;
`_transformvectorvar` reads `variable` for `uvar`/`vvar`, `gpvuar` for `gpuvar`,
`rottlgridvar` for `rottlgriduvar`, and `nlats` again. Four of the six sites
sit under `mode` values of `spectral`, `synchronous` or `syncfourier`; the
`variable` pair is under `mode='grid'`, in the arm that takes a spectral vector
pair with no level axis. Three of the names are MISSPELLINGS of names bound in
the same scope, which is the evidence that those branches have never run
anywhere.

`_transformvectorvar` also gave each of `rottlgriduvar` and `rottlgridvvar` one
half of the joined meridian instead of both, where the scalar loop in
`_transformvar` writes both halves of its one array. Repairing the names alone
would have left a branch that returns zeros for half of each component instead
of raising, which is worse.

`gcmt.py`, four more. `_csvData.__setitem__` names the file it writes in
neither branch (`fname`, and `lf.filestem` for `self.filestem`); `eqstream`
takes `file` and its body reads `dataset`; `orthographic` reads `zlon` for its
`lon` argument. `streamfxn`, the deprecated shim over `eqstream`, forwards a
`time` argument `eqstream` does not accept, and now refuses it rather than
raising `TypeError` from the callee.

`scripts/smoke_test.py` holds the vendored package to this from now on: see the
check named there. It is the instrument the `world-ro6` sweep should have used.

### Before the fork

The climate model is a personal fork vendored at `vendor/exoplasim/`, a git
subtree from the `master` branch of `Sohex/ExoPlaSim`. Pull upstream with
`git subtree pull --prefix vendor/exoplasim exoplasim-fork master --squash`.

The branch name differs from Orogen's deliberately. There, the fork's master
tracks upstream and `cf-fork` is a separate line of work; here the fork's master
IS the integrated line, so a second branch would be a copy with no owner. The
consequence to know: GitHub defaults the head of a new pull request from a fork
to that fork's default branch, so a PR against upstream must name its head
branch explicitly or it will offer the whole stack.

It is installed EDITABLE, so the source you read is the source that compiles and
the source that runs. Its build artifacts stay untracked: everything the build
writes lands under `vendor/exoplasim/build/`, one directory per configuration,
and in `plasim/run`, which is what keeps a twelve-binary rebuild from leaving
the working tree dirty.

`vendor/exoplasim/build/patched/` is the one exception to "one directory per
configuration", and it exists because a configuration is not the whole identity
of a binary. A verification arm that corrupts `plasim/src` in place writes a
control patch marker first; `build_model.py` reads that marker, builds under
`patched/` with a hash of what was patched in the directory name, and refuses
to publish at all. Nothing under `patched/` is a registry entry, and a build
directly under `build/` is always from the committed model source.

**The fork does not use upstream's build system and no longer carries it.**
`compile.sh`, `configure.sh`, `make_plasim`, the `most_compiler*` files,
`make_most` and the X11 launcher `most.c` it built are deleted, along with
`setup.py` and `MANIFEST.in`, which were setuptools inputs the declared
hatchling backend never read and which disagreed with `pyproject.toml` on both
version and dependencies; `plasim/CMakeLists.txt` and
`exoplasim/scripts/build_model.py` replace the build, `pyproject.toml` is the
only packaging declaration, and `exoplasim/__init__.py` raises rather than
compiling on demand -- its `sysconfigure()` and `printsysconfig()`, which ran
and read the generated configuration, are gone with the files they served. That
is the largest single divergence from upstream and it is deliberate: pulling
upstream will conflict there, and the resolution is always to keep this side.
`notes/audits/model-build-driver.md` says what the old one did wrong.

Upstream's sphinx documentation goes with it. `docs/` and `.readthedocs.yaml`
built the readthedocs site for the pip-installable package: an install path
this fork replaced, a feature list advertising modules it deleted, and an API
page for functions that are gone. Nothing here built it and a fork whose build,
install path and module set all differ has no upstream doc to keep current.
What the aerosol pages carried is in `aeromod.f90` and `radmod.f90`, which is
where a claim about the model is checked anyway; the PlaSim reference manual
and user guide under `plasim/doc` are kept and are a different tree. world-chv.

**The sibling models, the LSG ocean and the Earth boundary sets are deleted
too.** Upstream ships PUMA, SAM, CAT and their manuals, the octave and cat_tools
utility trees, an `images` gallery, the burn7 postprocessor, the full LSG source
and `plasim/dat`'s Earth and Mars surface sets alongside the model. Nothing here
built, read or shipped any of it: burn7 is superseded by `pyburn`, and LSG was
closed as an option by `notes/audits/ocean-and-marine-biosphere.md` finding 3
and `notes/external-model-survey.md` section 6. What is deliberately KEPT and
why: `plasim/src/specs` and `basespecs`, which are the ECOSTRESS spectra
`surfacespecs.py` was derived from and `analysis/ice_albedo.py` reads;
`stellarspectra`, a live search root in `run_exoplasim.py`; `hazeconstants`,
which `Model.configure` stages an aerosol file from; `tools`, because
`exoplasim/notes/restart-resolution-precision-converter.md` names `pusca.f90`'s
`convert_restart_array` as the base for the converter that note designs; and the
PlaSim reference manual, the PlaSim user guide and the PUMA user guide, which
document the model and the spectral core that run. The CAT manuals and the
climate report went with the models they describe.

The fork carries, over upstream: the low-I/O restart and broadcast repairs, the
pyburn reader fix, the shortwave weights for a non-solar host, the dust and
aerosol stack, and the stellar cycle. `exoplasim/patches/README.md` says which
are upstream pull requests and which are ours to keep.

## LPJ-GUESS CNP

The terrestrial biosphere model is vendored at `vendor/lpj-guess/` from
`mateusdp/LPJ-GUESS-NTD`, pinned to tag `LPJ-GUESS-CNP_v1.0` and commit
`b368b893c4324840b43c56866e901c6916858afb`. Pull that exact upstream state
with:

    git subtree pull --prefix vendor/lpj-guess \
      https://github.com/mateusdp/LPJ-GUESS-NTD.git \
      LPJ-GUESS-CNP_v1.0 --squash

The tagged fork omits the MPL-2.0 licence files present in the LPJ-GUESS 4.1.1
distribution. The subtree restores those files verbatim from the md5-verified
Zenodo 8065737 tarball; an upstream pull must not delete them.

The subtree also carries the Vesper calendar and astronomy port and the
`vesperinput` module directly. There is no second source copy and no patch stack:
the source under `vendor/lpj-guess/` is what is compiled. `framework/vesper.h`
is the exception only because it is generated from `config/planet.yaml` before
each build and therefore ignored.

### There is nothing upstream to integrate

The pull command above is documented so the pinned state can be reproduced, not
because it is outstanding work. The fork is at or ahead of LPJ-GUESS 4.1.1
everywhere, so a pull would bring back nothing but the licence headers the fork
deliberately omits.

Compared at the subtree import against the official `guess_4.1` release, 80
files differ, and 36 of those differ only by the stripped MPL-2.0 header. The 44
with real content differences are led by exactly the C-N-P files -- `somdynam.cpp`,
`canexch.cpp`, `guess.cpp`, `guess.h` -- which is the fork's own phosphorus work
rather than upstream drift; `somdynam.cpp` alone is 695 lines only in the fork
against 162 only in the release. The decisive test is the SVN `$Date` header per
file across `modules/` and `framework/`: the fork is older than 4.1.1 on ZERO
files, identical on 26, and newer on 8, those eight carrying 2022 dates against
the release's 2021. They are `canexch`, `commonoutput`, `driver`, `growth`,
`landcover`, `somdynam`, `vegdynam` and `framework/guess.cpp`.

**The method, so a future pull can be re-checked cheaply rather than
re-investigated.** Fetch the official tarball from Zenodo record 8065737
("LPJ-GUESS Release v4.1.1 model code", SVN r10118, MPL-2.0, published
2021-10-13; files `guess_4.1.1.tar.gz`, `.zip` and `releasenotes_4.1.1.txt`),
compare the `$Date` header of each file under `modules/` and `framework/`
against the vendored tree, and treat a diff that is only the licence header as
no diff at all. What matters is whether any file's release date is NEWER than
the fork's, and no file's is.

### `modules/ntransform.cpp` is stock 4.1.1 with declared divergences

Know this before opening that file. As it arrived, the soil nitrogen
transformation operator was byte-identical to `guess_4.1/modules/ntransform.cpp`
apart from the stripped licence header, and it is ON BY DEFAULT: `global.ins`
imports `global_soiln.ins`, which sets `ifntransform 1`. The CNP fork never
touched it. So every change this project has made there is a divergence from
default-on behaviour in a widely used community model, not a fork quirk it
inherited, and the same is true of `data/ins/global_soiln.ins`, where one
constant is off its release value.

Every one of them is declared. Each carries mainline's own line verbatim beside
the changed one in the source, under a `DECLARED DIVERGENCE FROM MAINLINE`
heading; `biosphere/config/ntransform.yaml` holds the register under
`mainline_divergences` with the verdict, what settles it and what it is worth;
and `biosphere/scripts/ntransform_gate.py` checks that mainline's form is
recorded, that it is not what the model runs, and that the changed line still
is, so a divergence can become neither a silent fork nor a silent revert. None
of it is execution-verified, because LPJ-GUESS does not build on this tree. The
argument for each is
`biosphere/notes/soil-nitrogen-transformation-parameterisation.md`.

### The phosphorus path diverges from THIS FORK, not from a release

`modules/somdynam.cpp` and `modules/soil.cpp` carry three declared divergences.
`PMASS_SAT`, the labile-P saturation threshold of the soil organic C:P ramp, is
converted out of Parton, Stewart and Cole (1988)'s resin-extractable currency
into the CNP fork's own Hedley-labile one; the surface humus pool, which the fork
left with no phosphorus ramp and a fixed C:P, ramps on the slow pool's line off
labile P, and its initialisation moves to that line's phosphorus-poor end. Each
records the fork's own line verbatim beside the changed one under the same
`DECLARED DIVERGENCE FROM MAINLINE` heading.

The reference point is what differs from `ntransform.cpp`. Release 4.1.1 has no
phosphorus at all, so there is no release form to compare against and no Zenodo
record or SVN revision to name: the whole path arrived with the CNP fork, and
the only fixed point is the subtree commit at the top of this section.
`biosphere/config/somdynam.yaml` holds the register keyed on that commit and
`biosphere/scripts/somdynam_gate.py` checks it, on the same three conditions as
the nitrogen gate plus the ramps those constants drive and the two lines the
phosphorus argument rests on. It also records which configuration each
divergence is live in, because two of the three are inert under the `ifplim 0`
this project runs and are waiting rather than harmless. The argument is
`biosphere/notes/phosphorus-cycle-parameterisation.md`.

### The snowpack's conductivity diverges from the RELEASE, and the relation is not declared in this subtree

`modules/soil.cpp`'s `update_snow_properties` carries one more declared
divergence, and its reference point is the release rather than the fork.
The function arrived byte-identical to `guess_4.1/modules/soil.cpp`'s apart from
the stripped licence header -- the CNP fork never touched it -- and it is ON BY
DEFAULT: `data/ins/global.ins` sets `iftwolayersoil 0`, which selects the
multilayer soil temperature scheme, and `Ksnow` becomes the conductivity of
every active snow layer in that scheme's numerical solve whether
`ifmultilayersnow` is 1 or 0. So it sits with `ntransform.cpp` rather than with
the phosphorus path: a change to default-on behaviour in a widely used community
model.

What changed is WHICH RELATION turns the simulated snowpack's density into a
conductivity. The release computes it from Sturm et al. (1997) and the climate
model computes it from Fourteau et al. (2021) Eq. (18), and at the snow density
`landmod` declares the two differ by close to a factor of two, so one snowfall
insulated the vegetation model's soil about twice as well as the climate
model's. `biosphere/config/snow_thermal.yaml` holds the register with mainline's
own lines recorded verbatim, and `biosphere/scripts/snow_thermal_gate.py` checks
it on the same three conditions the other two gates use.

**The relation itself is NOT declared in this subtree, and that is the point.**
`lib/snow.py` states it once. A Fortran model and a C++ model cannot import a
Python module at runtime, so `landmod.f90` and `soil.cpp` each carry the adopted
row as a literal and `snow.check_restatements()` holds both to that one table;
`scripts/smoke_test.py` runs it as well as the gate. Generating a header for the
vegetation model was the alternative, and the tree already generates
`framework/vesper.h` that way, but it covers only one of the two consumers,
because the climate model's Fortran is committed source that is read and edited
by hand. A route that fixes one side and leaves the other restating is two
declarations again with the drift moved.

The two cannot be reconciled by matching VALUES in any case. The vegetation
model's snow density is prognostic across the compaction ramp `modules/soil.h`
declares, where the climate model's is one namelist key inside that span, so
equal conductivities at one density would be a coincidence at one point of two
curves that diverge everywhere else. `notes/audits/cryosphere-material-properties.md`
argues the choice, the bracket it sits at the upper endpoint of, and what the
change is worth.

No other file under `vendor/lpj-guess/modules/` carries a register of this kind.
Where one does, it belongs beside these three.

Vendoring the CNP source does not itself enable phosphorus limitation.
`data/ins/global.ins` and the run harness keep `ifplim 0` until the gridded
weathering, sorption, deposition, and replacement productivity prediction are
ready to land as one scientific change.

## cGENIE

The candidate offline ocean is vendored at `vendor/cgenie/`, a git subtree from
the `master` branch of `derpycode/cgenie.muffin`, MIT licensed. Pull upstream
with `git subtree pull --prefix vendor/cgenie cgenie-fork master --squash`.

The remote is `cgenie-fork`, pointing at `Sohex/cgenie.muffin`, which is a fork
of `derpycode/cgenie.muffin`. A fork is expected rather than optional here,
because OCN-19 and OCN-20 are fork-shaped by construction: one deletes dead
duplication and threads BIOGEM's tracer loops, the other replaces the barotropic
streamfunction solve. As with ExoPlaSim, the fork's `master` is the integrated
line, so upstream is pulled into the fork first and this subtree pulls from the
fork.

That is also why it is a subtree and not an extraction under `references/`.
External source this project reads and will not edit is held there instead, and
`references/INDEX.md` records which trees those are and why.

**It builds, links and runs here.** Getting there took the five things below.
They are independent, and none of them is a source defect.

`genie-main/user.mak` sets `GENIE_ROOT = $(HOME)/cgenie.muffin` and
`RUNTIME_ROOT = ../../cgenie.muffin`, so the tree expects to sit at
`~/cgenie.muffin`, and `NETCDF_DIR=/usr/local` where this host has netCDF at
`/usr`. Overriding all three on the make line is enough to get past them.

It needs the **netCDF FORTRAN bindings**, which are a separate package from the
C ones and were absent here until 2026-08-25. Nothing else in this project
noticed, because every other netCDF consumer goes through Python's `netCDF4`;
`nf-config`, `libnetcdff` and `netcdf.mod` are what cGENIE wants and none of
them ships with `netcdf`. With `netcdf-fortran` installed the build gets past
every netCDF-dependent module.

**cGENIE compiles one executable per GRID, the same way ExoPlaSim does.**
`genie_control.f90` fixes the atmosphere as `ilon1_atm = GENIENX,
ilat1_atm = GENIENY` and `genie-goldstein/src/fortran/ocean.cmn` fixes the ocean
and sea ice as `GOLDSTEINNLONS/NLATS/NLEVS`. A bare `make` defines none of them,
and the `#ifndef` fallbacks in the two files DO NOT AGREE: the atmosphere falls
back to the 64 x 32 IGCM grid and the ocean to 36 x 36, so `genie.F`'s sea-ice
accumulation of an atmosphere field fails to compile with shapes that are not
conformable. That reads like a source defect and is a missing configuration.
The build therefore has to run THROUGH a config, because `makefile.arc` takes
the macros from `GENIE_FPPFLAGS` and only `genie.job` composes that, from the
`<build>` block of a config file. That is CLAUDE.md rule 4's shape on a second
model.

Two smaller things shape the command. The default make target also builds
`nccompare`, whose `src/c/compare.cpp` includes `netcdf.hh` from the legacy
netCDF C++ interface, which this host does not have, so the target has to be
`genie.exe`. And the build must be SERIAL: the dependency generator
`genie-main/finc.py` is Python 2 and fails under this host's `python`, so no
`.d` files are written and a parallel make races on `.mod` files.

A fourth thing appears only at larger grids. Every field cGENIE holds is in a
named COMMON block sized from the grid macros, so all of it is static; past
roughly a gigabyte of it the default `-mcmodel=small` cannot reach it and the
LINK fails with `relocation truncated to fit: R_X86_64_PC32 ... defined in
COMMON section`. A 72 x 72 x 16 ocean is past that line and needs
`-mcmodel=medium`. Like the grid macros, that is a per-configuration build flag
rather than a source change.

**The gfortran block compiles `-frecursive -fopenmp`, and the storage class is
a correctness setting rather than a preference.** Upstream's `-fno-automatic`
puts every procedure-body local in static storage, which is one variable for
the whole process and therefore SHARED between threads unless a `private`
clause names it; gfortran accepts `private` on such a variable rather than
rejecting it, and the only diagnostic the combination draws names no variable.
`-frecursive` puts the same locals on the stack, where each thread has its own,
so `genie.job` raises the shell's stack limit and sets `OMP_STACKSIZE` before
running the executable. That is required and not precautionary: with the
default limit the shipped `eb_go_gs_ac_bg` regression case segfaults before its
first timestep. The threaded routines are in `genie-goldstein`, and the thread
count is `OMP_NUM_THREADS` as usual.

`analysis/cgenie_cost.py` is the driver that puts those together; it records the
exact command and the toolchain in its provenance block, builds through a
generated config per grid, and times the result.
`notes/audits/cgenie-build-cost-and-grid-ceiling.md` is what it found, including
what the shipped configurations actually cover and where the usable resolution
stops.

`analysis/cgenie_profile.py` is the second driver, built on the first. It samples
retired instructions per symbol at two grids and maps each symbol to the file and
component that defines it, so it says which routine spends the cost rather than
what a configuration costs, and it carries the guard against this host's
`kernel.perf_event_max_sample_rate`, which throttles a sampler without reporting
a lost record. `notes/audits/cgenie-parallelism-and-coupling-support.md` is what
it found: what fraction of the work a thread team could divide, what Amdahl's law
then bounds, and the coupling resolution and regridding contract that budget
supports. It builds from a `git archive` export outside the repository, because a
worktree's ignored build products are symlinks into the main checkout and a build
in place would write there.

`analysis/cgenie_omp.py` is the third driver and the one that CHANGES the tree
rather than measuring it as it stands. It builds named ARMS -- the tree serial,
the tree threaded, the tree compiled with uninitialised locals poisoned -- and
holds every one to the same acceptance test: every float variable in a shipped
regression case's reference netCDF, bit-for-bit against what a named git
revision writes. It also counts retired instructions with `perf stat` on the
audit's own cost case, which is exact on a single-threaded process and is not a
cost on a threaded one, since a thread waiting at a barrier retires
instructions in proportion to how long it waits.
`notes/audits/cgenie-embm-free-path-and-threading.md` is what it found.

Run output goes to `OUT_DIR = $(HOME)/cgenie_output`, outside this repository,
which is why the ignore rules here cover only what a build leaves in the tree.

**What it carries that open rows already name.** `genie-goldstein` is the
frictional-geostrophic ocean OCN-19 and OCN-20 are about, and its `invert.f`
ordering is what OCN-20 quotes. `genie-knowngood` ships the four reference
configurations OCN-19 validates against. `genie-ecogem` is the ecosystem tier
OCN-4 names beside MARBL. And `genie-plasim` is not a coupling module but a
whole PlaSim: a strict SUBSET of `vendor/exoplasim`'s module set -- no `p_exo`,
no aerosol core, no MPI, no `hurricanemod` or `glaciermod` -- plus one file,
`geniemod.f90`, 71 lines of pure declaration. So this subtree carries the Earth
model ExoPlaSim's exoplanet capability was added to, and the entire
atmosphere-ocean coupling surface is one array-declaration module. It bears on
OCN-17's evaluation of the published ExoPlaSim-to-cGENIE path, and on OCN-10,
whose forcing contract can be read off it rather than derived. Note that the two
published paths are different arrangements: the offline regrid route moves wind
stress, winds and albedo, while this in-tree one hands over seventeen fields.
`notes/external-model-survey.md` section 30.

**IT IS NON-DIMENSIONALISED AGAINST EARTH'S RADIUS AND GRAVITY, and half of the
scales are configurable so the other half is easy to miss.**
`genie-goldstein/src/fortran/initialise_goldstein.F:376-391` hardcodes
`rsc = 6.37e6` and `gsc = 9.81` and exposes neither, while `sodaylen`,
`sidaylen`, `yearlen`, `nyear` and the depth scale `par_dsc` are all in
`ini_gold_nml`. Four derived scales then carry the hardcoded pair everywhere:
`tsc = rsc/usc` sets every non-dimensional time, `rhosc` carries both, `opsisc`
scales the reported overturning streamfunction and `rfluxsc` the heat flux.

At this planet's 1.20 Earth radii and 1.306 Earth gravity, three of those are
wrong by a factor of 1.20 and `rhosc` by about eight percent. **Nothing will
warn.** The model is entirely non-dimensional, so an overturning reported
through `opsisc` comes back twenty percent low, in the right units, looking
exactly like an answer.

Rotation is the one planetary constant the source does parameterise, and it
FAILS OPEN. `fsc` takes `4*pi/sidaylen` only when the solar and sidereal day
lengths differ by more than 0.001, and otherwise reverts to Earth's
`2*7.2921e-5` for backwards compatibility. A configuration that sets the two
equal gets Earth's Coriolis scaling silently. Set them to this planet's values
and check `fsc` took the intended branch before believing any run.

None of this is an argument against the model, and OCN-12 owns the full
inventory. It is here because it is the thing most likely to produce a
plausible wrong number for somebody who assumed a namelist covered the planet.
`notes/external-model-survey.md` section 10e has the derivation and the
precedent, including a scale that was already wrong once and corrected outside
the code.

**`genie-paleo`'s plots are deleted from this checkout, and the rest is kept.**
That directory arrived at 1.4 GB, of which 678 PostScript plots of Earth
paleogeographic reconstructions were 1290 MB. They are pictures, nothing here
will read them, and removing them took the subtree from 1.6 GB to 289 MB.

What was KEPT is the other 49 MB and the reason is specific: 486 base
configurations under `genie-main/configs` point into those directories through
`ea_1`, `go_1`, `gs_1` and `bg_par_pindir_name`, and the `.k1`, `.paths`,
`.psiles` and `.dat` files they resolve to are the format reference for what a
cGENIE bathymetry configuration IS. OCN-11 has to produce one from the Orogen
export, and OCN-18 asks specifically what muffingen's `.k1`, `.paths` and
`.psiles` generation does above 36 x 36. Deleting the whole directory would have
cost that and saved 49 MB.

The deletion is IN THE FORK, not only in this checkout, so a `git subtree pull`
will not reintroduce the plots. `git subtree split --prefix=vendor/cgenie`
reconstructs real upstream lineage rather than rooting at the squash commit, so
the split branch was a descendant of `008fd490` and pushed to the fork's master
as a fast-forward carrying only the deletion. That is the general recipe for
sending a change up: split, push the split branch to `cgenie-fork master`,
delete the split branch.

One consequence remains. The history here still carries the deleted blobs, so
`.git` did not shrink and will not; the saving is working-tree size, which is
what `grep` and `find` pay.

Unrelated and worth knowing before someone blames the deletion: four directories
that configs reference, `fkl_np10`, `fkl_pp01_DH`, `fm0450ab` and `wppcont1`,
are absent upstream and were never in the vendored commit.

## LPJmL

LPJmL is vendored at `vendor/lpjml/`, a git subtree squashed from the `master`
branch of `PIK-LPJmL/LPJmL` at commit `572e2b906ac2c55b2ee6661a93e4633b126254e4`.
Pull upstream with:

    git subtree pull --prefix vendor/lpjml \
      https://github.com/PIK-LPJmL/LPJmL.git master --squash

**It is a subtree and not a `references/` extraction for the usual reason**: the
work intended on it is fork-shaped and it is planned into the pipeline. It was
first pulled down as read-only comparison material for SPITFIRE, its
process-based fire scheme, because FIRE has ten open rows of ten issued and the
only scheme consulted was BLAZE, which arrived with LPJ-GUESS. That reading is
still worth doing and is now incidental to why it is here.

**It is NOT the biosphere, and it does not supersede LPJ-GUESS.** This project's
terrestrial biosphere is the LPJ-GUESS CNP fork above, and nothing in that
changes. LPJmL is a different model in the same LPJ family -- same lineage,
different code base, different scope -- and having both vendored is deliberate
rather than duplication. Read `vendor/lpj-guess/` when the question is about this
project's biosphere; read `vendor/lpjml/` when the question is about LPJmL.

**Nothing reads it yet**, exactly as with cGENIE: it is vendored source with no
pipeline step, so `config/pipeline.yaml` has no row for it and rule 7's "what is
now worthless" question does not reach it. When a step does consume it, that step
and its artifact go into the pipeline graph in the same commit.

**It ships its own `.gitignore`**, covering objects, the binaries `configure.sh`
generates into `bin/`, the `Makefile.inc` that script writes, and the `output/`
and `restart/` directories. Nested ignore files are honoured, so unlike cGENIE --
whose own ignore file carries only `.DS_Store` -- this subtree needs no build
rules added to the repository root.

**The licence is AGPL-3.0**, which is stronger copyleft than the other three
subtrees carry: MIT for cGENIE, MPL-2.0 for LPJ-GUESS, GPL for ExoPlaSim. That
is not a problem for reading or for local modification, and it is worth knowing
before any of it is copied into code that leaves this repository.
