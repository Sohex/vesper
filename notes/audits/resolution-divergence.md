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
The build path refuses cleanly almost everywhere. The model's own data path
never refuses at all.

Two findings are already owned and are not restated: the hyperdiffusion branch
at `plasim.f90:1443` and the unreachable `jtune` radiation table at
`radmod.f90:947`, both in `world-ys9` and in the note above.

## 1. A filter-stripped control binary occupies the registry's build directory

The strongest finding is about the tree as it stands rather than about a rung.

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
build directory half is untouched, and `--print-path` hands its caller exactly
that path. The `drop+extra` hash in the tag exists to keep arms apart; a source
patch is an arm the tag cannot see.

## 2. The restart carries its resolution and nothing reads it

`plasim.f90:904-907` writes `nlat`, `nlon`, `nlev` and `nrsp` into the restart.
No `get_restart_integer` call anywhere in `plasim/src` asks for any of the four;
`read_atmos_restart` takes `nstep` and `naccuout` and stops there.

The read itself is `read (nreaunit) pa(1:k1,:)` at `restartmod.f90:155`, with no
`iostat`. A Fortran unformatted sequential read that consumes fewer values than
the record holds is legal and advances to the next record. So a donor written at
a HIGHER truncation is silently reinterpreted: the target takes the leading
`k1*k3` values of a longer record, which is not a truncated field but a flat
reslicing that runs out of the donor's level one into its tail. A donor at a
LOWER truncation hits end of record and the program dies without an `iostat` to
name why.

`restart_schema.py:776-789` does validate all four, but that is the converter.
An unconverted restart from another rung passes every gate the run path has.

## 3. A missing surface field is not an error, and the default planet is all land

`get_surf_array` at `surfmod.f90:174-176` handles an absent file by writing
`* Init <name> [code = NNNN] internally *` to the diagnostic unit, leaving
`kread` at 0, and returning. `mpsurfgp` then leaves the array untouched, so the
module default stands.

The default that matters is `plasimmod.f90:710`, `dls = 1.`: the land/sea mask
is all land. With `doro` at 0.0 the model integrates a flat, entirely
continental planet, and says so only in a line indistinguishable from the
thirteen other codes it initialises internally on a healthy run.

The same routine reads the eight-word SRA header into `ih(:)` at `:134` and
`:141` and never compares `ih(5)`, `ih(6)` against `NLON`, `NLAT`. The header
carrying the grid is read and discarded. `surface_ini` at `:303-362` compounds
it: any `surface.sra` is re-split and relabelled under the RUNNING model's NLAT,
so a T21 field becomes `N064_surf_*.sra` with no check.

The project's own `read_sra` compares the entire header and refuses, and
`surface_field_report` raises when an intended code is absent. Those are the
defences that make this latent rather than live.

## 4. The postprocessor selects its transform on the wrong variable

`pyburn.py:649` sets `nlat = min(header[4], header[5])`, which is the latitude
count. Six sites then branch on `if nlat in [192,320]` to choose the FFT
extension. 192 and 320 are the LONGITUDE counts of T63 and T106.

`pyfft.f90:791` accepts `8,16,32,48,64,96,128,256,384,512,1024,2048,4096`.
Neither 192 nor 320 is in it, and `pyfft.f90:1193-1199` ends in a bare Fortran
`stop`, which from inside an f2py extension terminates the Python process with
status 0: no exception, no traceback, no output.

| rung | nlat | nlon | selects | outcome |
| --- | --- | --- | --- | --- |
| T63 | 96 | 192 | `pyfft` | process death |
| T106 | 160 | 320 | `pyfft` | process death |
| T127 | 192 | 384 | `pyfft991` | correct, for the wrong reason |

The model binary is unaffected: `build_model.py:60` carries
`NEEDS_FFT991 = {"T63", "T106"}` and routes correctly. Only the postprocessor's
selector is wrong, and the fix is `nlat` to `nlon` in six places.

## 5. The FFT module T63 and T106 select is not thread safe

`fftmod.f90:26` carries `!$omp threadprivate(lastn,nallowed,trigs)`, added by
this fork so a thread owns what a rank owned. `fft991mod.f90` carries no
`threadprivate` directive at all, while holding the same mutable planner state:
`lastn`, `ifax(10)` and the allocatable `trigs`.

`fft991mod.f90:38-43` does its lazy init on that shared state, deallocating and
reallocating `trigs` and rerunning `set99` whenever `n` changes. `build_model.py:168`
selects `fft991mod` for exactly T63 and T106, and this project builds `omp`. So
those two rungs, and only those two, put every thread into a race over the
planner while other threads read it.

Not currently reachable, because neither rung is in the build matrix. It is
armed by the act of adding one.

## 6. The two FFT tables disagree with each other and with the ladder

| rung | NLON | in `fftmod.f90:7` | in `fft991mod.f90:12` |
| --- | --- | --- | --- |
| T21 | 64 | yes | yes |
| T31 | 96 | yes | yes |
| T42 | 128 | yes | yes |
| T63 | 192 | no | yes |
| T85 | 256 | yes | yes |
| T106 | 320 | no | yes |
| T127 | 384 | yes | yes |
| T170 | 512 | yes | no |

`fft991mod.f90:19` documents `T170 - N512 : 8-4-4-4` while its own `nallowed`
omits 512, and `set99` factorises 512 without complaint. The refusal is the
table, not the algorithm.

`fftmod`'s table also accepts `NLON = 48`, which its radix loop cannot do:
`la = 6` satisfies `la >= 4`, one `dfft4` pass sets `la = 1`, and neither the
radix-3 nor the radix-2 tail fires, so 32 of the required 48 is applied. The
inverse at `:80-90` overshoots instead. No error, wrong transform. Unreachable
today because NLAT 24 is not a rung; it is a trap for anyone adding T15.

## 7. Silent resolution dependence in the physics

The four branches the earlier note found are the places the model KNOWS about
resolution. These are places it does not.

**`landhoskn0 = 15.0`** at `plasimmod.f90:452` is an absolute wavenumber, and
its own comment states the assumption: the default gives 0.1 at n=N for T21. The
Lander-Hoskins physics filter at `legmod.f90:184-186` therefore confines the
tail at T21 and removes most of the resolved spectrum at T85 and above. The
neighbouring exponential filter normalises by NTRU and does not have this
problem, which is what makes the omission legible.

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
| storm size, 30 flagged gridpoints | `hurricanemod.f90:71` | an unweighted count, 1.46% of the globe at T21 and 0.37% at T42 |
| snow cover fraction, `dsnow/(dsnow+0.01)` | `landmod.f90:406`, `simba.f90:483` | the canonical subgrid snow parameter, 1 cm scale |
| critical relative humidity, 0.85 | `rainmod.f90:92` | encodes an assumed within-cell humidity distribution |
| convective cloud fit, 0.245 and 0.125 | `rainmod.f90:1904` | a fit to gridcell mean convective rain rate |
| land or sea at mask > 0.5 | `icemod.f90:213` and eleven more | which cells count as ocean is rung dependent |
| glacier flag at snow depth > 30.0 | `glaciermod.f90:211` | a partially glaciated coarse cell and a full fine cell differ |

**`gamma` is a fraction per timestep, not a rate.** `rainmod.f90:2238` and three
siblings divide by `deltsec2`, so the mass evaporated per step is
`gamma*(qsat-q)*dsigma` independent of step length. The T42 override at
`rainmod.f90:88` is a partial correction for exactly this, applied at one rung.

**`newsnow.f90:5-6`** declares `NLAT = 32`, `NLON = 64` as its own parameters and
never references `resmod`. The offline snowpack accumulator cannot be right at
any other rung and does not know the grid it was handed.

## 8. `NLPP` truncates, and the polar band goes unowned

`plasimmod.f90:102` is `parameter(NLPP = NLAT / NPRO)`, integer division,
restated verbatim in four other modules. `mpstart` checks that the runtime rank
or thread count equals the compiled `NPRO`; nothing in the Fortran checks that
`NPRO` divides `NLAT`.

When it does not, `NPRO*NLPP < NLAT` and the consequences are all silent: the
scatters and gathers cover only part of the globe, `makeareas` at
`utilities.f90:24-36` reads cell areas from latitude slots the gather never
filled, and under `OMPSHARED` the pointer association at `plasimmod.f90:1099`
leaves the top of the globe owned by no thread, so the forward transform
analyses a field whose polar band is whatever the array was initialised to.

The only guard is `CMakeLists.txt:72-77`, added by this fork, which refuses the
build. That is the right place for it. It is worth knowing it is the ONLY place,
because a hand-edited `resmod.f90` or a direct `cmake` call bypasses it.

## 9. Only two rungs are provisioned, and three have binaries without data

The model's surface files are named `N<NLAT>_surf_<code>.sra` at
`surfmod.f90:190`, from the compile-time NLAT. Provisioning as it stands:

| rung | binary in the matrix | `plasim/dat/` | `exoplasim/inputs/` |
| --- | --- | --- | --- |
| T21 | yes | 15 codes | 8 codes |
| T31 | no, the build refuses | 15 codes | none |
| T42 | yes | 15 codes | 7 codes, no 0229 |
| T63 | no | none | none |
| T85 | yes | none | none |
| T106 | no | none | none |
| T127 | yes | none | none |
| T170 | yes | none | none |

T85, T127 and T170 have executables and no surface data anywhere. Combined with
finding 3 that is a flat all-land planet under T21's damping. The project-side
`declare_hyperdiffusion` and `surface_field_report` both refuse before the model
is reached, which is why this is a hazard rather than an incident.

T31 is the inverse case: a complete data set, an `mpstep` branch, a `radmod`
tuning branch, and a build that cannot happen because `CMakeLists.txt:57` omits
48 from its allowed latitude counts while `lib/rungs.py` includes T31.

## 10. The ladder is written down nine more times

`lib/rungs.py` is the declared single source and is self-checking. Surviving
restatements, and how each differs:

| where | content | divergence |
| --- | --- | --- |
| `exoplasim/scripts/cache_budget.py:33-39` | seven rungs with NTRU, NLAT, NLON | invents `T213`, which is not a rung, and offers it through `choices`; drops T63 and T106; its comment claims it mirrors `build_model.py` |
| `vendor/exoplasim/exoplasim/plasim/CMakeLists.txt:57` | allowed NLAT | omits 48, so T31 cannot build |
| `vendor/exoplasim/exoplasim/__init__.py:545-571` | full dispatch chain | no T31; `force991` is assigned and never read; the T106 arm tests `"T106"` twice where `"t106"` was meant |
| `exoplasim/scripts/verify_transform_roundtrip.sh:29-33` | a `case` over six rungs | refuses T63 and T106 as unknown |
| `exoplasim/scripts/stability_probe.py:63` | five rungs | also globs the binary across layers and ranks |
| `exoplasim/scripts/filter_timestep_matrix.py:103` | five rungs | restates SPAT-2's reasoning in the comment above the copy |
| `exoplasim/scripts/verify_gauss_weights.sh:54` | five latitude counts | 48, 96 and 160 never checked |
| `exoplasim/scripts/verify_latitude_pairing.py:44-48` | a case list | comment says it covers the ladder; three rungs absent |
| `scripts/archive_builds.py:56-57` | five expected grids | disagrees with both builds on disk, so the newer one can never be archived |

The SPAT-2 lint does not reach most of these. `smoke_test.py:995` globs `*.py`
only, so every shell copy is unpoliced, and the path test at `:867` recognises
`.sra`, `.nc`, `.json` and `.rest` but not `.x`.

## 11. Artifacts that survive a rung change

`lib/paths.py:44` resolves the baseline climatology from a config key with a
clean-IO guard and no grid guard, so changing `model.resolution` neither renames
nor refuses it. Seven consumers then take the grid from the file. `scripts/error_budget.py`
mixes rungs inside one product: its albedo half resolves through
`rungs.model_grid(config)` at `:186` while its hydrology half reads whatever the
climatology carries at `:459`. `check_consistency.py:99` builds the land/sea mask
the same way, so a T21 coupling can pass the convention check for a T42 run.

`lib/builds.py:177-194` returns any `T<digits>` tag without consulting the
ladder, so `exoplasim-T99` yields `"T99"` and three surface-field generators
write `inputs/t99/`. `grid_export` in the same module does call `rungs.geometry`
and refuses, which is the pattern this one is missing.

`config/pipeline.yaml:85` declares the hydrography coupling by glob rather than
by rung, so a stale coupling reads as present for a rung that never generated it
and `--purge` will not name it.

## 12. One Earth constant found on the way, which is not about resolution

`oceanmod.f90:24` is `parameter(PLARAD=6.371E6)`. The module uses only `resmod`,
not `pumamod`, so this is the ocean's ONLY definition of the planet radius and
the configured value never reaches it. The grid there is angular, `dlam` being
`2*pi/NLON` with `cphi` and `dphi` in radians, so `zfac = hdiffk/plarad/plarad`
at `:1344` is what supplies the metric.

At the declared 1.20 Earth radii, the simulated ocean's horizontal heat
diffusion therefore runs 1.44 times stronger than the coefficient it is given:
an arm labelled 1000 m2/s behaves as 1440.

`nhdiff` defaults to 0 in `oceanmod` and `config/planet.yaml` sets
`horizontal_diffusion: false`, so no current run reaches it. `clim-16` was closed
on a MEASURED 300/1000/3000 bracket, and those arms did reach it, so the
diffusivity labels on that result are the ones affected. Its conservation claim
does not depend on the coefficient's value and is untouched.

This belongs to the family in `notes/audits/inherited-earth-constants.md`, whose
clean list does not cover `oceanmod`.

## Checked and clean

Recorded so they are not re-derived.

**The rewritten `inigau` is the model for how a rung-dependent numerical property
should be handled.** `gaussmod.f90:29-36` states a floor of roughly `2*n*eps`,
about 1.1e-13 at NLAT 256, and measures it. It flattened the SHTns against
legmod divergence from two orders of growth across the ladder to within a factor
of about 20, and the analysis arm now passes a 1e-11 bar at every rung.

**`resmod.f90` is generated into the build directory**, at
`CMakeLists.txt:134-141`, precisely so a build for one resolution cannot be read
as another's. The `resmod_def.f90` in the source tree is a dead template.

**Rank divisibility is enforced twice at build time**, at `CMakeLists.txt:72-77`
and `build_model.py:138-142`, and once at runtime in each parallel module.

**`verify_gauss_weights.sh`'s tolerance is set at the hardest rung** and applied
everywhere with the derivation stated, which is the correct direction.

**Output records are self-describing.** `outmod.f90:15-22` and the per-record
headers carry NLON, NLAT, NLEV and NTRU, so a reader that honours the header
cannot mis-size a record. This is the contrast that makes finding 2 legible.

**The postprocessor's decoding is otherwise fully generic.** `pyburn` takes every
dimension from the record headers and derives the truncation as `(nlon-1)//3`;
`gcmt.py` carries no grid literal at all.

**The climatological ozone path is latent, not live.** `radmod.f90:1227` zeroes
`dqo3cl` and reads it from a `_0237.sra` that exists at no rung, but `no3`
defaults to 1 at `:173` and nothing in this project sets 2, so the synthetic
distribution is what runs.

## What this audit did not cover

Stated so the coverage is not overread. I did not audit the SHTns library source
under `vendor/shtns-src`, the LSG ocean or the coupled `cpl` path, the CAT, SAM
and PUMA trees beside `plasim`, the `burn7` C postprocessor, or the vertical
discretisation beyond the NLEV conditionals that sit next to the horizontal ones.
Resolution dependence in the consumers downstream of the climate model, in
hydrography, pedology and the biosphere, is the subject of the SPAT rows and is
not touched here.

## Tasks

Tracked in the `bd` issue tracker under the `resdiv` label, not restated here.
