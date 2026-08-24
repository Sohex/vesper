# Where the climate model's behaviour changes with resolution, and where it changes without saying so

Audited 2026-08-24, against the ExoPlaSim fork's source at
`vendor/exoplasim`, its build system, the project-side drivers under
`exoplasim/scripts` and `lib/`, and the data actually staged on disk. Claims
about the current tree state were checked against the artifacts rather than
against the scripts that intend them.

Worldbuilding frame: this file is about the Vesper simulation's climate model
and the code that drives it. Every quantity named here is a modelled quantity
of an invented world, and nothing here is about the real one.

## The question

`exoplasim/notes/resolution-tuned-parameters.md` asked what "tuned for coarse
resolution" means and found four places the model source branches on
truncation. This asks the wider question that leaves open: across the whole
fork, and not only its physics, where does a change of rung change behaviour,
and where does it change behaviour SILENTLY?

The distinction that organises everything below is not coarse against fine. It
is whether a rung the code was not prepared for produces an ERROR or a NUMBER.
The build path refused cleanly almost everywhere. The model's own data path
refused nowhere at all, and that asymmetry is what most of the fixes below
closed: the surface reader, the restart reader and the thread decomposition each
refuse by name now.

Two findings were already owned and are not restated: the hyperdiffusion branch
at `plasim.f90:1544` and the unreachable `jtune` radiation table at
`radmod.f90:968`. The first was a `NTRU==42` preset; world-rt1 and world-8bs
between them made the project write every hyperdiffusion key on every segment,
so the preset is overwritten rather than relied on. The second is `world-ys9`
and is **still open**.

## 1. A filter-stripped control binary occupied the registry's build directory

The strongest finding was about the tree as it stood rather than about a rung.

`build_model.py` composes a build directory tag from resolution, layers, ranks,
parmode, profile, frame pointers and the flag delta. The STATE OF THE MODEL
SOURCE is not in that tag. `verify_shtns_model.sh` builds its control arm by
patching `shtnsmod.f90` in place to remove the spectral filter, so the arm and
the registry entry claim the same tag and the same directory.

Measured on 2026-08-24:

| artifact | sha256 head | filter sites |
| --- | --- | --- |
| `build/t21_l10_p16_omp_production/plasim.x` | `6523babd32fbb827` | 0 |
| `plasim/run/most_plasim_t21_l10_p16_omp.x` | `e0830f776ecb9a4f` | 11 |

The preprocessed source beside the first, at
`build/t21_l10_p16_omp_production/CMakeFiles/plasim.x.dir/src/shtnsmod.f90-pp.f90`, carries the `CONTROL PATCH IN PROGRESS` marker and no filter sites at all.

`--no-publish` fixed the half of this that overwrote the shipped binary. The
`drop+extra` hash in the tag existed to keep arms apart, and a source patch was
an arm the tag could not see.

**Closed by world-70k.** `build_model.py` scans `plasim/src` for the control
markers `CONTROL PATCH IN PROGRESS` and `! CONTROL:` before it composes a tag
(`patched_sources()`). When any file carries one, the build goes under
`build/patched/` with `source_state()`, a hash of the patched files, in the
directory name, so two control arms cannot collide with each other or with the
registry's directory; and publishing is refused outright with a message naming
the files and telling the caller to pass `--no-publish` or restore the source.
`--print-path` still hands the caller exactly the path it asked for, which is
now the patched root.

## 2. The restart carried its resolution and nothing read it

`epilog` writes `nlat`, `nlon`, `nlev` and `nrsp` into the restart
(`plasim.f90:927-930`). No `get_restart_integer` call anywhere in `plasim/src`
asked for any of the four; `read_atmos_restart` took `nstep` and `naccuout` and
stopped there.

The read itself is `read (nreaunit) pa(1:k1,:)` at `restartmod.f90:155`, with no
`iostat`. A Fortran unformatted sequential read that consumes fewer values than
the record holds is legal and advances to the next record. So a donor written at
a HIGHER truncation is silently reinterpreted: the target takes the leading
`k1*k3` values of a longer record, which is not a truncated field but a flat
reslicing that runs out of the donor's level one into its tail. A donor at a
LOWER truncation hits end of record and the program dies without an `iostat` to
name why.

`restart_schema.py` does validate all four, but that is the converter, and an
unconverted restart from another rung passed every gate the run path had.

**Closed.** `read_atmos_restart` reads all four back
(`plasim.f90:1227-1230`) and refuses a restart whose grid is not this model's,
with the comment at `:1200` recording that `epilog` had been writing them for
nothing. The converter is what a caller reaches for once the refusal names the
mismatch.

## 3. A missing surface field is still not an error, and the default planet is all land

`get_surf_array` handles an absent file by writing
`* Init <name> [code = NNNN] internally *` to the diagnostic unit, leaving
`kread` at 0, and returning. `mpsurfgp` then leaves the array untouched, so the
module default stands. **This half is unchanged**, and is `world-qml`.

The default that matters is `plasimmod.f90`'s `dls = 1.`: the land/sea mask is
all land. With `doro` at 0.0 the model integrates a flat, entirely continental
planet, and says so only in a line indistinguishable from the thirteen other
codes it initialises internally on a healthy run.

**The wrong-grid half is closed by world-fuh.** The same routine read the
eight-word SRA header and never compared words 5 and 6 against `NLON` and
`NLAT`; it now refuses by name, `Surface file on the wrong grid`, printing both
the header's dimensions and the model's (`surfmod.f90:170-187`). A named refusal
rather than a warning, because the alternative is a run that produces numbers.
`surface_ini` still re-splits and relabels any `surface.sra` under the RUNNING
model's NLAT, so a file arriving by that route is renamed rather than checked --
but the per-code reader downstream of it now catches the mismatch.

The project's own `read_sra` compares the entire header and refuses, and
`surface_field_report` raises when an intended code is absent. Those remain the
defences that keep the absent-field half latent rather than live.

## 4. The postprocessor selected its transform on the wrong variable

`pyburn` set `nlat = min(header[4], header[5])`, the latitude count, and six
sites then branched on `if nlat in [192,320]` to choose the FFT extension. 192
and 320 are the LONGITUDE counts of T63 and T106.

`pyfft.f90:791` accepted `8,16,32,48,64,96,128,256,384,512,1024,2048,4096`.
Neither 192 nor 320 was in it, and `pyfft.f90:1195-1203` ends in a bare Fortran
`stop`, which from inside an f2py extension terminates the Python process with
status 0: no exception, no traceback, no output.

| rung | nlat | nlon | selected | outcome |
| --- | --- | --- | --- | --- |
| T63 | 96 | 192 | `pyfft` | process death |
| T106 | 160 | 320 | `pyfft` | process death |
| T127 | 192 | 384 | `pyfft991` | correct, for the wrong reason |

The model binary was unaffected: `build_model.py` carries
`NEEDS_FFT991 = {"T63", "T106"}` and routed correctly.

**Closed by world-i38, and the selector is DERIVED rather than listed.**
`_fftmodule(nlon)` at `pyburn.py:33-71` decides on the longitude count and on
the factorisation `8 * 4**k * r` with `r` in 1, 2, 3, which is what the radix-4
module can transform, so a rung nobody has added yet routes correctly without a
new list entry. Six literal `nlat in [192,320]` tests became one call.

## 5. The FFT module T63 and T106 select was not thread safe

`fftmod.f90:32` carries `!$omp threadprivate(lastn,nallowed,trigs)`, added by
this fork so a thread owns what a rank owned. `fft991mod.f90` carried no
`threadprivate` directive at all, while holding the same mutable planner state:
`lastn`, `ifax(10)` and the allocatable `trigs`, with a lazy init that
deallocated and reallocated `trigs` and reran `set99` whenever `n` changed.
`build_model.py` selects `fft991mod` for exactly T63 and T106, and this project
builds `omp`, so those two rungs, and only those two, put every thread into a
race over the planner while other threads read it.

It was never reachable, because neither rung is in the build matrix; it was
armed by the act of adding one. **Closed by world-988**:
`fft991mod.f90:54` carries `!$omp threadprivate(lastn,nallowed,ifax,trigs)`, and
the block above it records that gfortran gives every thread its own copy and
what the per-thread cost is.

## 6. The two FFT tables disagreed with each other and with the ladder

As audited:

| rung | NLON | in `fftmod`'s `nallowed` | in `fft991mod`'s |
| --- | --- | --- | --- |
| T21 | 64 | yes | yes |
| T31 | 96 | yes | yes |
| T42 | 128 | yes | yes |
| T63 | 192 | no | yes |
| T85 | 256 | yes | yes |
| T106 | 320 | no | yes |
| T127 | 384 | yes | yes |
| T170 | 512 | yes | no |

`fft991mod` documented `T170 - N512 : 8-4-4-4` while its own `nallowed` omitted
512, and `set99` factorises 512 without complaint: the refusal was the table,
not the algorithm. `fftmod`'s table also accepted `NLON = 48`, which its radix
loop cannot do -- `la = 6` satisfies `la >= 4`, one `dfft4` pass sets `la = 1`,
and neither the radix-3 nor the radix-2 tail fires, so 32 of the required 48 is
applied, and the inverse overshoots instead. No error, wrong transform.
Unreachable then and now, because NLAT 24 is not a rung.

**Closed by world-78e: each table now holds what its own module can transform.**
`fftmod.f90:7` reads `8,16,32,64,96,128,256,384,512,1024,2048,4096` -- 48 gone,
96 present -- and `fft991mod.f90:12` reads `4,64,96,128,192,256,320,384,512`,
with 512 added and the reason recorded at `:22`. `pyfft.f90:791` carries the
same correction as `fftmod`. The tables are the two modules' capability, not a
restatement of the ladder, which is why they differ from each other and from
`lib/rungs.py` on purpose.

## 7. Silent resolution dependence in the physics

The four branches the earlier note found are the places the model KNOWS about
resolution. These are places it does not.

**`landhoskn0 = 15.0`** at `plasimmod.f90:455` is an absolute wavenumber, and
its own comment states the assumption: the default gives 0.1 at n=N for T21. The
Lander-Hoskins physics filter therefore confined the tail at T21 and removed
most of the resolved spectrum at T85 and above, while the neighbouring
exponential filter normalises by NTRU and does not have that problem, which is
what made the omission legible. **Closed by world-cjk**: `legmod.f90:193` sets
`zhoskn0 = landhoskn0 * real(NTRU) / zhoskref`, so the critical wavenumber
scales with the truncation and the filter damps the same FRACTION of every
rung's spectrum; the reference at `:176-192` records that the declared value is
T21's and reproduces it there exactly.

**Global means are area-weighted, and the four in `gridpointa` had no consumer.**
The mean 2m temperature, precipitation, evaporation and outgoing longwave were
formed as `sum(zgp) / (NLON * NLAT)` with no Gaussian weight, which on an
unequally spaced grid over-weights high latitudes by an amount that is a
function of NLAT. Measured on the way to fixing it: the only reader of
`t2mean`, `precip`, `evap` and `olr` is `energy`, which passes them to
`guiput`, and `plasim/CMakeLists.txt:156` compiles `guimod_stub.f90`, whose
`guiput` is a bodyless `return`. So no published number moved, and the
quadrature error had never reached one. The prints that a person does read did
move: the land and sea percentages from `ncountsea` were a fraction of the
CELLS, and the "Mean:" topography and icesheet-height lines in `surfmod` and
`glaciermod` were arithmetic means over the grid.

The weighting now goes through two helpers in `plasim.f90`, `gpareamean` for a
distributed field and `ugpmean` for a gathered one, so there is one place that
turns a gridpoint field into a global mean rather than a convention that half
the call sites follow. `icemod.f90` and `utilities.f90` already weighted
correctly.

**Absolute per-cell thresholds, now declared as anchored.** Each of these
compares a gridcell mean, or a raw cell count, against a constant chosen for one
cell size, with no NLAT term. Each now says so at its declaration, and
`exoplasim/notes/resolution-tuned-parameters.md` section 4 carries the list a
convergence experiment has to account for. `snowcovz` was a bare 0.01 at seven
sites and is now one named key in `landmod_nl`.

| what | where | why it moves |
| --- | --- | --- |
| storm size, 30 flagged gridpoints | `hurricanemod.f90:70-82` | an unweighted count, 1.46% of the globe at T21 and 0.37% at T42 |
| snow cover fraction, `dsnow/(dsnow+snowcovz)` | `landmod.f90:17-25` | the canonical subgrid snow parameter, 1 cm scale |
| critical relative humidity, 0.85 | `rainmod.f90:132-138` | encodes an assumed within-cell humidity distribution |
| convective cloud fit, 0.245 and 0.125 | `rainmod.f90:1926-1931` | a fit to gridcell mean convective rain rate |
| land or sea at mask > 0.5 | `icemod.f90` and eleven more | which cells count as ocean is rung dependent |
| glacier flag at snow depth > 30.0 | `glaciermod.f90:266` | a partially glaciated coarse cell and a full fine cell differ |

**`gamma` is a fraction per timestep, not a rate.** `rainmod.f90:2255` and three
siblings divide by `deltsec2`, so the mass evaporated per step is
`gamma*(qsat-q)*dsigma` independent of step length. Upstream had a `NTRU==42
.and. NLEV==10` override of 0.007 against 0.01 elsewhere, which reads as a
partial correction for exactly this applied at one rung. **world-j6v removed
it**: what `gamma` is worth depends on the STEP LENGTH and not on the
truncation, this project gives T21, T42 and T85 the same 45-minute step, so the
override could not have been compensating one, and a constant that steps 43 per
cent on a rung change makes every cross-rung comparison carry a
precipitation-physics change it did not ask for. `gamma` is one number at every
rung and is a `rainmod_nl` key. The timestep dependence itself is untouched and
is `world-khn`'s class, not this one's.

**`newsnow.f90` and `buildice.f90`** declared `NLAT = 32`, `NLON = 64` as their
own parameters and never referenced `resmod`, so the offline snowpack
accumulator could not be right at any other rung and did not know the grid it
was handed. **Closed by world-zsa**: both `use resmod` and take `NLAT` from
`NLAT_ATM`, so a mismatch is a build error. Whether either program should exist
is `clim-53` and is open.

## 8. `NLPP` truncated, and the polar band went unowned

`plasimmod.f90` has `parameter(NLPP = NLAT / NPRO)`, integer division, restated
verbatim in four other modules. `mpstart` checked that the runtime thread count
equalled the compiled `NPRO`; nothing in the Fortran checked that `NPRO` divides
`NLAT`.

When it does not, `NPRO*NLPP < NLAT` and the consequences are all silent: the
scatters and gathers cover only part of the globe, `makeareas` at
`utilities.f90:24-36` reads cell areas from latitude slots the gather never
filled, and under `OMPSHARED` the pointer association at `plasimmod.f90:1099`
leaves the top of the globe owned by no thread, so the forward transform
analyses a field whose polar band is whatever the array was initialised to.

The build guard at `plasim/CMakeLists.txt:71-79`, added by this fork, refuses
such a build, and that is the right place for it. It was the ONLY place, so a
hand-edited `resmod.f90` or a direct `cmake` call bypassed it. **Closed by
world-ljj**: `mpstart` in `mpimod_omp.f90` refuses a thread count that does not
divide `NLAT` at run time, with the consequences recorded at `:570-575`, so the
two guards now cover the two ways in.

## 9. Only one rung is provisioned, and the rest have binaries without data

The model's surface files are named `N<NLAT>_surf_<code>.sra` at
`surfmod.f90:203`, from the compile-time NLAT. Provisioning as audited, with the
vendored Earth and Mars boundary sets in `plasim/dat/` counted:

| rung | binary in the matrix | `plasim/dat/` | `exoplasim/inputs/` |
| --- | --- | --- | --- |
| T21 | yes | 15 codes | 8 codes |
| T31 | no, the build refused | 15 codes | none |
| T42 | yes | 15 codes | 7 codes, no 0229 |
| T63 | no | none | none |
| T85 | yes | none | none |
| T106 | no | none | none |
| T127 | yes | none | none |
| T170 | yes | none | none |

**`plasim/dat/` no longer exists.** world-6ak deleted it, 17 MB of Earth and
Mars boundary sets, as a vendored tree nothing builds, reads or ships. So the
`plasim/dat/` column is now empty at every rung and `exoplasim/inputs/` is the
whole of the provisioning: T21 and T42 only.

T85, T127 and T170 have executables and no surface data anywhere. Combined with
finding 3's unchanged half that is a flat all-land planet under T21's damping.
The project-side `declare_hyperdiffusion` and `surface_field_report` both refuse
before the model is reached, which is why this is a hazard rather than an
incident.

T31 was the inverse case: a complete vendored data set, an `mpstep` branch, a
`radmod` tuning branch, and a build that could not happen because
`plasim/CMakeLists.txt` omitted 48 from its allowed latitude counts while
`lib/rungs.py` included T31. **Closed by world-ajx**: 48 is in the allowed list
at `plasim/CMakeLists.txt:60` with the reason at `:55-59`, and the two
restatements that disagreed with the ladder were brought onto it. T31 can build;
it has no `exoplasim/inputs/` data, which puts it with T63 and above.

## 10. The ladder was written down nine more times

`lib/rungs.py` is the declared single source and is self-checking. The
restatements as audited, and how each differed:

| where | content | divergence | state |
| --- | --- | --- | --- |
| `exoplasim/scripts/cache_budget.py` | seven rungs with NTRU, NLAT, NLON | invented `T213`, which is on no ladder, and offered it through `choices`; dropped T63 and T106 | closed, world-qjq: takes `rungs.geometry` and `sorted(rungs.RUNGS)` |
| `vendor/exoplasim/exoplasim/plasim/CMakeLists.txt` | allowed NLAT | omitted 48, so T31 could not build | closed, world-ajx |
| `vendor/exoplasim/exoplasim/__init__.py` | full dispatch chain | no T31; `force991` assigned and never read; the T106 arm tested `"T106"` twice where `"t106"` was meant | closed, world-ajx: see the block at `__init__.py:350-356` |
| `exoplasim/scripts/verify_transform_roundtrip.sh` | a `case` over six rungs | refused T63 and T106 as unknown | closed, world-ltc: resolves through `rungs.geometry` |
| `exoplasim/scripts/stability_probe.py` | five rungs | also globbed the binary across layers and ranks | closed, world-ltc |
| `exoplasim/scripts/filter_timestep_matrix.py` | five rungs | restated SPAT-2's reasoning in the comment above the copy | closed, world-ltc: imports `rungs` |
| `exoplasim/scripts/verify_gauss_weights.sh` | five latitude counts | 48, 96 and 160 never checked | closed, world-ltc: the sweep is `rungs.RUNGS` |
| `exoplasim/scripts/verify_latitude_pairing.py` | a case list | comment said it covered the ladder; three rungs absent | still a list, and correctly so: `CASES` at `:43-48` pairs NLAT with NPRO to exercise both sides of the divisibility fallback, which is a test matrix and not a ladder restatement |
| `scripts/archive_builds.py:56-57` | five expected grids | it gates archiving at `:150`, and the ACTIVE build has neither `exoplasim-T63` nor `grid-512x256`, so that build can never be archived | **UNCHANGED**; `EXPECTED_GRIDS` is still a literal list, and no lint reaches it. `world-txz` |

`world-ltc` also widened the SPAT-2 lint that had missed most of them:
`smoke_test.py` globbed `*.py` only, so every shell copy was unpoliced, and its
path test recognised `.sra`, `.nc`, `.json` and `.rest` but not `.x`. The lints
read shell scripts, know `.x`, and catch a table.

## 11. Artifacts that survived a rung change

`lib/paths.py`'s `climatology_path` resolved the baseline climatology from a
config key with a clean-IO guard and no grid guard, so changing
`model.resolution` neither renamed nor refused it, and seven consumers then took
the grid from the file. `scripts/error_budget.py` mixed rungs inside one
product: its albedo half resolves through `rungs.model_grid(config)` while its
hydrology half read whatever the climatology carried.
`check_consistency.py`'s `land_sea_mask` built the mask the same way, so a T21
coupling could pass the convention check for a T42 run.

**Closed by world-876.** The grid is checked in `climatology_path` itself
(`lib/paths.py:71-75`), so the one resolver refuses a climatology that is not on
the configured rung and every consumer inherits the guard rather than repeating
it. `error_budget.py:473-477` records what that buys: before
`require_configured_grid` existed a budget could combine a T21 hydrology with a
T42 albedo and say nothing.

`lib/builds.py`'s `resolution_of` returned any `T<digits>` tag without
consulting the ladder, so `exoplasim-T99` yielded `"T99"` and three
surface-field generators wrote `inputs/t99/`; `grid_export` in the same module
called `rungs.geometry` and refused, which was the pattern this one was missing.
**Closed by world-yop**: `resolution_of` calls `rungs.geometry(tag)`
(`lib/builds.py:197`) and refuses anything off the ladder.

`config/pipeline.yaml` declared the hydrography coupling by glob rather than by
rung, so a stale coupling read as present for a rung that never generated it and
`--purge` would not name it. **Closed by world-yop**: the `hydrography` step's
`writes` names `coupling_exoplasim-{res}.nc` explicitly, with the glob kept
beside it so an older rung's coupling is still enumerable.

## 12. One Earth constant found on the way, which is not about resolution

`oceanmod.f90` carried `parameter(PLARAD=6.371E6)`. The module used only
`resmod`, not `pumamod`, so this was the ocean's ONLY definition of the planet
radius and the configured value never reached it. The grid there is angular,
`dlam` being `2*pi/NLON` with `cphi` and `dphi` in radians, so
`zfac = hdiffk/plarad/plarad` is what supplies the metric.

At the declared 1.20 Earth radii, the simulated ocean's horizontal heat
diffusion therefore ran 1.44 times stronger than the coefficient it was given:
an arm labelled 1000 m2/s behaved as 1440.

**Fixed by world-mll**: `hdiffo` takes the radius from `planet_nl` through
`pumamod` (`oceanmod.f90:1376`, applied at `:1394`) and `oceanini` refuses
`nhdiff > 0` with no radius (`:291-296`).

`nhdiff` defaults to 0 in `oceanmod` and `config/planet.yaml` sets
`horizontal_diffusion: false`, so no current run reaches it. `clim-16` was closed
on a MEASURED 300/1000/3000 bracket, and those arms did reach the old code, so
the diffusivity labels on that result are the ones affected -- 432/1440/4320 is
what ran. Its conservation claim does not depend on the coefficient's value and
is untouched. world-vho is the relabelling.

This belongs to the family in `notes/audits/inherited-earth-constants.md`, whose
clean list does not cover `oceanmod`.

## Checked and clean

Recorded so they are not re-derived.

**The rewritten `inigau` is the model for how a rung-dependent numerical property
should be handled.** `gaussmod.f90:29-36` states a floor of roughly `2*n*eps`,
about 1.1e-13 at NLAT 256, and measures it. It flattened the SHTns against
legmod divergence from two orders of growth across the ladder to within a factor
of about 20, and the analysis arm now passes a 1e-11 bar at every rung.

**`resmod.f90` is generated into the build directory**, precisely so a build for
one resolution cannot be read as another's. The `resmod_def.f90` that sat beside
it in the source tree was a dead template and was deleted under world-cmz.

**Thread divisibility is enforced at build time**, at
`plasim/CMakeLists.txt:71-79` and in `build_model.py`, and, since world-ljj, at
run time in `mpstart`. The rank wording this entry used to carry is gone with
the MPI path: world-38b removed the MPI and serial builds and the parmode axis,
so `mpimod_omp.f90` is the only decomposition and `NPRO` is the thread count.

**`verify_gauss_weights.sh`'s tolerance is set at the hardest rung** and applied
everywhere with the derivation stated, which is the correct direction.

**Output records are self-describing.** `outmod.f90` and the per-record headers
carry NLON, NLAT, NLEV and NTRU, so a reader that honours the header cannot
mis-size a record. That was the contrast that made finding 2 legible, and the
restart reader now honours its own header too.

**The postprocessor's decoding is otherwise fully generic.** `pyburn` takes every
dimension from the record headers and derives the truncation as `(nlon-1)//3`;
`gcmt.py` carries no grid literal at all.

**The climatological ozone path is latent, not live.** `radini` zeroes
`dqo3cl` and reads it from a `_0237.sra` that exists at no rung, but `no3`
defaults to 1 and nothing in this project sets 2, so the synthetic distribution
is what runs. world-ayx placed that synthetic profile on a pressure-equivalent
basis; see `notes/audits/model-earth-centrism.md` finding 11.

## What this audit did not cover

Stated so the coverage is not overread. I did not audit the SHTns library source
under `vendor/shtns-src`, the LSG ocean or the coupled `cpl` path, the CAT, SAM
and PUMA trees beside `plasim`, the `burn7` C postprocessor, or the vertical
discretisation beyond the NLEV conditionals that sit next to the horizontal ones.
All but the first and the last are moot now: world-6ak and world-cmz deleted the
LSG tree, `cpl.f90`, the CAT, SAM and PUMA trees and `burn7.x`, none of which
any build compiled or any run reached.
Resolution dependence in the consumers downstream of the climate model, in
hydrography, pedology and the biosphere, is the subject of the SPAT rows and is
not touched here.

## Tasks

Tracked in the `bd` issue tracker under the `resdiv` label, not restated here.
