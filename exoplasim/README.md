# Superhabitable-world climate workflow

This component converts the canonical World Orogen geography in `../source/` to
ExoPlaSim boundary conditions and runs reproducible climate experiments. All
commands below are run from the project root.

It is described here by what it does and how to run it, and it carries no
results. `world_state.json` holds current values, `runs/INDEX.json` records what
each run physically was, and the dated measurements live in `notes/`.

## Reproduce geography conversion

On Arch Linux the required host tools are `gcc-fortran` and `openmpi`. Activate
`.venv`, install `requirements.txt`, then run:

```bash
python exoplasim/scripts/build_boundary_conditions.py   # land mask + topography
python exoplasim/scripts/build_surface_albedo.py        # background albedo
python exoplasim/scripts/run_exoplasim.py --run-years 1
```

`build_boundary_conditions.py` integrates the mask and topography from the
native mesh. Do not go back to remapping the equirectangular PNGs: that is a
lossy intermediate now the fork emits Gaussian grids directly, and the PNG mask
uses the `elevation > 0` convention, which floods every dry closed-basin floor.

The first command writes generated model inputs to `exoplasim/inputs/t42/` and
geography diagnostics to `exoplasim/analysis/geography/`. The second prepares the 0.90-S-Earth
experiment, records a complete manifest, and runs one smoke orbit. It refuses
to overwrite an existing run with climate output. See
`exoplasim/notes/parameter-decisions.md` for physical and format assumptions.

### Output mode: low-I/O to spin up, clean to finish

**Spin-up orbits run `NLOWIO = 1`; any orbit a climatology will be built from
runs `NLOWIO = 0`.** The two regimes write different things and the difference is
not a quality setting, it is what the output MEANS.

    # spin up cheap
    python exoplasim/scripts/continue_exoplasim.py --run <id> --orbits 40 \
        --purpose spinup --low-io

    # then the orbits you will actually read
    python exoplasim/scripts/continue_exoplasim.py --run <id> --orbits 10 \
        --purpose post_equilibrium_climatology

**Why the split.** Under `NLOWIO = 1` the model accumulates over the output
interval and writes 12 bins; under `NLOWIO = 0` it writes instantaneous records
which `pyburn` then averages to the same 12. Nothing downstream changes shape.
But an accumulation cannot be undone, and a sample set is strictly more
information, since the mean can be recovered from it. The accumulation is also
not the average you would compute yourself: binned `spd` under `NLOWIO = 1` sits
between the speed of the time-mean vector and the mean of instantaneous speeds.
So anything reading variance, extremes or single records -- DUST-5's gust
distribution above all -- needs `NLOWIO = 0`. Scalars and slowly-varying fields
do not care, which is why spin-up does not.

**A separate defect, now fixed.** `NLOWIO = 1` used to corrupt the first output
record of every model call, because `naccuout` survived a restart while the
accumulators did not. `exoplasim-3.4.2-lowio-first-record.patch` repairs that and
was VERIFIED on 2026-08-18 against a binary with it reversed. It does not touch
the accumulation semantics above, so it does not merge the two regimes. Every run
currently in `exoplasim/runs/` was written before it and still needs the
corrections in `exoplasim/notes/first-output-bin.md`.

**Cost.** The clean regime's postprocessing overhead has not been measured
since the reader fixes and is probably near 1.0x; re-measure with one segment
in each regime on a run made after both fixes before quoting a ratio.
`notes/audits/pyburn-postprocessing-cost.md`.

What is NOT in doubt is the raw volume: about 2.4 GB an orbit at `NLOWIO = 0`
against about 96 MB at `NLOWIO = 1`, so the clean regime writes roughly 25x the
bytes even though they are deleted after postprocessing. Model time differs by a
few percent at most. Delete run directories once their climatologies are
extracted.

**So the regime is chosen on INFORMATION, not on speed.** Low I/O writes interval
accumulations and the clean regime writes instantaneous samples; a mean can be
recovered from samples and an accumulation cannot be undone. Buy samples for the
orbits something will read and not for the orbits nothing will.

**What CHOOSES it.** The regime follows the declared purpose
rather than a flag anyone has to remember. `run_exoplasim.py` prepares a run and
integrates it toward equilibrium, so its block is low I/O by default and
`--clean-io` overrides; it now also registers that block as a `spinup` segment,
so the orbits it writes carry a purpose and a `low_io` like any other.
`continue_exoplasim.py` derives the default from `--purpose`: `spinup` gets the
cheap regime, `post_equilibrium_climatology` and `diagnostic` get samples, and
`--low-io`/`--clean-io` force either.

**What enforces it.** `--purpose post_equilibrium_climatology` with `--low-io` is
refused outright. Every segment records `low_io`, `build_climatology.py` refuses a
tainted window unless given `--allow-low-io`, and a segment with no `low_io` key
counts as tainted, because the honest default for a segment nobody labelled is to
assume the cheap regime. `close_term_energy.py` and `close_state_energy.py`
require clean windows for the same reason. See `segments.py` and CLIM-9.

### Which resume path is valid

`continue_exoplasim.py` continues the SAME run: same config, same surface fields,
more orbits. That is the normal spin-up path and it has no caveat.

`run_exoplasim.py --restart-from` seeds a NEW run from another's restart. It is
valid only when the surface fields are identical and the FORCING differs -- flux
or CO2, which are namelist parameters rather than restart state. Use it for a
flux bracket, where it saves the whole cold spin-up.

It is refused, by a guard, whenever a surface field differs, and the refusal is
correct rather than conservative: `landmod`'s `landini` takes `dwmax`, `dz0clim`
and all three `dalbcl` bands from the restart when `restart > 0`, so the `.sra`
files are read only on a cold start. A seeded run with a new albedo or a new soil
water field would discard it and reproduce its parent while every manifest
recorded the new values. Every step of loop A adds a surface field, so every step
of loop A needs a cold start. See `docs/src/practice/failure-modes.md` class 13.

Continue a validated experiment from its latest restart (seasonal snapshots are
written by default; `--no-seasonal-output` skips them for a segment, and a run
without them has to be extended before it can produce a climatology):

```bash
python exoplasim/scripts/continue_exoplasim.py --run <run_id> --orbits 5 \
  --purpose spinup
```

`--purpose` is required and is not inferred from the other flags. It says what
the orbits are FOR: `spinup` while the run is approaching equilibrium,
`post_equilibrium_climatology` for orbits meant to be read as this world's
climate, and `diagnostic` for orbits run to measure the model rather than the
planet -- an I/O verification, a high-cadence wind sample for DUST-5, a block on
a differently patched binary. A diagnostic segment does not move the run's
status, is skipped by the convergence window, and is refused as climatology
input. The vocabulary lives in `segments.py`.

A resume also compares the STELLAR SPECTRUM BY CONTENT, not by name. The config
names `k25v`, `build_stellar_spectrum.py` writes `k25v.dat` from the `star`
block -- spectral type, effective temperature and metallicity -- and it rewrites
it in place, so the whole radiative input can change while every recorded name
stays the same. The run manifest carries a sha256 of both spectrum files and a resume
refuses across a change to either. A run prepared before that existed is stamped
on its first resume and says so.

A resume also compares the SURFACE FIELDS BY CONTENT, and for the same reason.
`surface_field_report` checks that the `.sra` files are present in the run
directory, which is the whole of what a cold start needs; a resume takes those
fields out of the restart instead, so presence says nothing about what the model
will read. `restart_surface.py` reads the restart's own copy of every code this
project stages and requires it to equal, exactly, the field the run was built
from. The mapping from surface code to restart record is checked against
`surfmod.f90`'s `surfcode` table so a renumbering upstream fails loudly, and
`--self-test` builds the restart a substitution would have written and requires
the check to refuse it.

Two things are asserted and they have different references. First, the restart
still carries the surface the run ADOPTED: the staged `.sra` for a cold-started
run, `MOST_REST.seed` for one seeded with `--restart-from`. A departure means the
model re-derived a boundary field between segments, and that is never
overridable. Second, the staged `.sra` in the run directory IS what the model
reads -- true after a cold start, and deliberately false for a run prepared with
`--superseded-surface-ok`, which adopts a donor's surface and discards the staged
one. `stage_surface_extras` then rewrites the current `.sra` into that directory
on every resume, so the directory and the manifest come to assert a surface the
model will not read. That is refused on resume and restated with
`--superseded-surface-ok` on `continue_exoplasim.py`, which stamps the segment
with what it is actually integrating.

Codes 129 and, under `NVEG = 2`, 212 and 229 are excluded and reported by name
rather than dropped: the model transforms the orography before it reaches the
restart, and coupled vegetation owns forest cover and field capacity once a run
starts, so in neither case is the staged file the right answer.

Assess a spin-up and create a separate five-orbit seasonal climatology only
after it passes:

```bash
# Run ids are UUIDs and carry no meaning, so start from the index.
python exoplasim/scripts/index_runs.py

python exoplasim/scripts/assess_convergence.py exoplasim/runs/<run_id>
# The TOA criterion above reads a diagnostic that is offset from the planet's
# actual energy tendency; this says by how much on the run in front of you.
python exoplasim/scripts/close_state_energy.py exoplasim/runs/<run_id> \
  --first 67 --last 76
# --run is required: a continuation cannot recompute a name, and being handed
# one cannot silently resolve to a different run.
python exoplasim/scripts/continue_exoplasim.py --run <run_id> --orbits 5 \
  --purpose post_equilibrium_climatology
python exoplasim/scripts/build_climatology.py exoplasim/runs/<run_id> \
  --start-year 46 --end-year 50
python exoplasim/scripts/analyze_climatology.py
```

Current results are in `world_state.json`, never here. What each completed run
was is in `exoplasim/runs/INDEX.json`; run ids are UUIDs, so that index is the
only thing that maps one to its physics.

Principal products are in `exoplasim/analysis/climatology/`:

- `baseline_regular_climatology.nc`: all retained variables in 12 time bins.
- `baseline_snapshot_climatology.nc`: 32 orbital-phase samples.
- `baseline_classification.nc`: annual fields plus categorical class indices.
- `baseline_climate_report.json`: global diagnostics and area fractions.
- PNG maps for temperature, precipitation, circulation, snow, sea ice,
  hydrology, Köppen--Geiger classes, and broad ecological interpretation.

The Köppen map is explicitly rate-normalized: precipitation rates are
annualized to 365.2425 days before applying empirical Earth thresholds. That
avoids classifying this world's 180.7-day orbital year as artificially dry.
It is a worldbuilding interpretation, not a dynamic vegetation simulation.

## Building the model

`exoplasim/scripts/build_model.py` is the only way the model gets built. It
reads the flag line from `config/planet.yaml`, checks every argument against a
list, and builds through CMake and Ninja in a directory of its own per
configuration.

    python exoplasim/scripts/build_model.py --res T170 --ranks 16
    python exoplasim/scripts/build_model.py --res T21  --ranks 8

**An argument it does not recognise is an error, and that is the point.** What
this replaced defaulted silently three times over -- an unrecognised `-r` built
T21, an unrecognised `-p` built SINGLE precision under a declaration of eight,
and the executable was named for the resolution the parse arrived at, so a
caller that copied the name it asked for could get an untouched binary from a
previous session. `notes/audits/model-build-driver.md` has the measurements.

| option | what it selects | name |
| --- | --- | --- |
| `--ranks <n>` | the OpenMP thread count; NLAT must divide by it. There is one parallel layer, `mpimod_omp` over a shared address space, so there is no parmode to choose and no suffix to carry | `most_plasim_<res>_l10_p<n>.x` |
| `--profile` | a flag set from `config/planet.yaml`; `production` by default, `checked` adds `-fcheck=all`, and `poisoned` adds `-finit-real=snan` at `-Og`, which is the one level the poison survives to the arithmetic at | no change to the name |
| `--frame-pointers` | a PROFILING build. DWARF cannot unwind the -O3 code -- 94% of model samples get no caller -- and this costs -1.17% with an identical restart sha, so the profile measures the same model | `..._fp.x` |
| `--extra-flag=`, `--drop-flag=` | for a verification arm that varies a flag ON PURPOSE. Both are part of the build directory's identity, and a `--drop-flag` naming a flag the declaration does not carry is an error rather than a no-op | no change to the name |

Precision is not an option: `config/planet.yaml` declares it. It is the worst of
these to get wrong, because `-fdefault-real-8` changes the width of `real` in
every declaration and a real*4 object linked against real*8 ones does not fail
to link -- it computes.

The paired latitude decomposition is retired. It existed to let legmod fold a
mirror pair together, SHTns replaces legmod and cannot use the permuted layout
anyway, and carrying two grid layouts was a cost with nothing left to buy.

`NSHTNS` is a NAMELIST switch and deliberately not a build flag, so the two
transforms can be compared inside one binary with no compiler difference in the
comparison -- which is what `verify_shtns_model.sh` relies on. It defaults to 1
on the threaded build and 0 on the MPI build, which cannot run SHTns at all;
legmod remains reachable with `NSHTNS = 0` and is the reference that gate
compares against.

## The ecological stream

A second output stream, off by default, written for the biosphere rather than
for a climate diagnostic. `NECO = 1` in `plasim_nl` turns it on and `NECOSTEP`
sets its interval in timesteps; left at zero, `NECOSTEP` becomes `mtspd`, which
is the number of timesteps in one absolute 24-hour day exactly, since prolog
recomputes `mpstep` so that `mtspd * mpstep * 60 = day_24hr`. It writes
`plasim_eco` on unit 143.

WHY IT EXISTS. The two alternatives for driving a nonlinear ecological model
are both bad: accepting the twelve-bin seasonal reduction as weather, or writing
the model's whole raw payload every timestep, which is about 15 GB per T42
orbit. This writes nineteen surface fields plus its own interval bounds, reduced
at the producer with the operator each field's meaning calls for.
`biosphere/notes/ecological-forcing-field-contract.md` is what it carries and
why; `outmod.f90`'s `ecoaccu`, `ecogp` and `ecoreset` are the whole of the
implementation.

Three properties worth knowing:

- **It cannot change a result.** The restarts written with `NECO = 1` and
  `NECO = 0` are byte-identical: it reads state that already exists and writes
  to a unit of its own.
- **Its interval survives a restart.** The nineteen accumulators, the four
  extrema and the counter `naccueco` are all serialized, behind an `ecovers`
  marker on the same terms as `accuvers`, so an interval that straddles a
  segment boundary is closed at the right step and covers the span it declares.
  `compare_eco_streams.py` is the check, and it separates a structural failure
  from a bounds failure from a payload failure because they have different
  owners.
- **Its VALUES across a boundary are only as exact as the model's state.** They
  are not: `world-8yyh`, measured in
  `notes/audits/ecological-stream-restart-continuity.md`.

The codes are 600 to 602 for the interval and 610 to 628 for the fields;
`compare_eco_streams.py` carries the name of each. Nothing teaches them to
pyburn yet, because EFOR-3 reads the stream directly rather than through it.

## If you add or remove a restart record

A `put_restart_integer`, `put_restart_real`, `put_restart_seed`,
`put_restart_array`, `mpputgp` or `mpputsp` call site is not a local change.
Three tools read the restart by NAME and none of them can guess at a record
nobody has classified: `convert_restart.py` treats an unknown name as a hard
error rather than passing it through, `reset_restart_accumulators.py` decides
from the same table whether a record is a partial accumulation to clear before
a seeded run, and `restart_surface.py` gates every resume on named records.

**One file to update: `exoplasim/scripts/restart_schema.py`.** Add a `POLICY`
entry saying what the record IS -- its semantic class, what a conversion does
with it, and for a gridpoint field whether remapping it averages an intensity
or moves an inventory, on which surface class, within what bounds. The
mechanical half -- which names exist, from which module, with what shape and
which variable behind them -- is parsed out of the call sites themselves and
out of `plasim/CMakeLists.txt`, so nothing there needs touching.

**The gate is `scripts/smoke_test.py`,** check "the restart schema covers every
record the model writes". It goes red for a call site with no policy, for a
policy entry no call site writes, and for an accumulator whose reset in
`outreset` no longer matches what the schema records. Its own control shows it
failing in all three directions.

**An accumulator's clean value is not always zero.** `tempmin` resets to 1.0e3
and `atsami` to 1.0e10, both running minima; `asndch` and `aanrho` are reset by
the model nowhere. The schema carries all four and the gate holds them against
the source, which is how `atsami` was found after a by-eye reading had missed
it.

Nothing else needs changing. `convert_restart.py --self-test` and
`reset_restart_accumulators.py --self-test` both read the schema, so they cover
the new record the moment it has a policy.

## The land liquid water scheme is a selection

`nlandwcol` in `landmod_nl` chooses between two schemes, and the default is the
one that has always run.

- **`nlandwcol = 0`, the scalar bucket.** One store `dwatc`, one capacity
  `dwmax`, and runoff is the store overflowing. Bit-identical to what this model
  did before the selection existed.
- **`nlandwcol = 1`, the layered column.** `nlsoilw` layers, shaped by
  `dsoilwf`, with the lower boundary `nlandwdrain` selects: impermeable, or free
  drainage that produces a `ddrain` flux out of the base. The surface flux
  enters the top layer, a withdrawal deeper than the top layer holds draws on
  the layers below, each layer fills to capacity and passes its excess down, and
  whatever the base cannot take backs up and leaves as surface runoff.

**The bucket is the layered column at one layer with an impermeable base, and
that is checked rather than claimed.** The arithmetic lives in
`plasim/src/landcolumn.f90` as two kernels that depend on nothing -- no
plasimmod, no resmod, no grid -- so `verify_land_column_reduction.sh` compiles
them alone and drives both with one flux sequence, requiring bitwise equality.
That file is separate for exactly this reason: a kernel that can only be run by
running the model can only be compared by running the model.

Bit-exactness is a design constraint there and not an aspiration.
`column_step` takes a layer to its capacity with `AMIN1(cap, w)` rather than the
algebraically identical `w - AMAX1(0., w - cap)`, because the second is not the
same number in floating point and a reduction right to a tolerance cannot tell a
refactor from a physics change.

**`dwatc` remains the column total under either scheme.** `fluxmod`'s
evaporation limiter, `simba`'s water stress factor, `aeromod` and `outmod` all
read it and none of them knows about layers, so the layered scheme sums its
layers into it every step. At one layer the sum is over one element and is
exact.

**The evaporation limiter is three axes rather than one shape.** `drhsfull`,
`drhslow` and `nrhsexp` give
`beta = min(1, max((theta - drhslow)/(drhsfull - drhslow), 0)) ** nrhsexp`,
and the defaults are this model's active form exactly. The default path is taken
by BRANCH and not by algebra: `x**1.0` with a real exponent is not bitwise `x`,
which is also why `nrhsexp` is an integer. The point of the general form is that
the other implementation this project has read, cGENIE's ENTS, sits on the same
three axes at a different point -- exponent four with the knee at a full store,
which at equal fractional fill differs from this model by a factor of 39 at 0.4
of capacity and 2 at 0.8 -- so the disagreement becomes a runtime bracket rather
than a code fork. Neither limiter is right; that is what a bracket is for.

### Soil phase, and the two columns that share no boundary

`nlandwphase = 1` freezes and thaws the liquid in each water layer against the
soil temperature at its depth. It needs `nlandwcol = 1`, because phase is a
property of a layer and the scalar bucket has no layer to freeze, and it is off
by default: this model's five soil temperature layers carry no water and no
phase, so melt water always infiltrates whatever the soil temperature.

**Ice occupies pore space, and that is the whole infiltration impedance.** A
layer's capacity in `column_step` is its capacity less its ice, so a frozen
layer fills and overflows sooner. No conductivity is involved, which is just as
well: the column has none, and `pedology/config/land_column_properties.yaml`
carries `flow.saturated_conductivity` as undeclared. `frozen_impedance` in
`landcolumn.f90` is the conductivity form, declared for the gradient-driven
hypothesis that would need it and used by nothing.

**Water and energy close for the same reason.** The mass that changes phase and
the temperature change are the same number read two ways: a layer at
`tmelt - dT` freezes at most the water whose latent heat of fusion would raise
it back to `tmelt`, and that latent heat is exactly what the layer's
temperature is then moved by. The latent heat is `als - alv`, because this model
declares vaporisation and sublimation and derives fusion as the difference, and
that is the constant `landmod.f90`'s own snowmelt already uses.

`verify_land_column_reduction.sh` checks both identities across a sweep in
temperature and in how much water is present, so the energy-limited and
water-limited branches are visited in both directions. **The tolerance is
DERIVED and not chosen**: the energy identity is read through a temperature
difference near 273 K, so its relative precision is degraded by the ratio of the
absolute temperature to the change, and each case carries its own bound computed
from machine epsilon and that ratio. The reported number is the worst residual
over its own bound, which comes out the same at default real and at
`-fdefault-real-8` -- the confirmation that the bound is the right scaling
rather than a number that happened to pass. Control 4 books the exchange at the
latent heat of vaporisation instead of fusion and must break the energy identity
while leaving the water identity intact, which is exactly what a wrong latent
heat does.

**`dsoilwz` is the water column's layer thicknesses, and the model has never had
them.** `dwmax` is a capacity in metres OF WATER and says nothing about depth;
`dsoilwf` is a share of that capacity. Phase needs a depth, because the
temperature that decides it lives on the `dsoilz` soil temperature layers and
those share no boundary with the water layers. The mapping is declared and is by
midpoint: a water layer takes the temperature of whichever temperature layer
contains its centre. The default is one water layer of 1.5 m, which is the
property contract's column base and LPJ-GUESS's physical profile.

**The snow half is not here.** ExoPlaSim's snow density is the constant
`rhosnow = 330` kg/m3 with no compaction while LPJ-GUESS ages its snow from 275
to 500, and snow DEPTH is what sets the insulating thickness over the soil, so
the two columns insulate differently from the same snowfall. GRAV-8 owns that,
and it is the one snow term that carries gravity.

**What is NOT here.** The third registered hypothesis is gradient-driven flow,
and it needs an unsaturated conductivity and a matric potential.
`pedology/config/land_column_properties.yaml` carries both as undeclared, so
registering two schemes and naming the third is the honest state rather than a
gap. Neither scheme has an infiltration capacity either, so saturation excess is
the only runoff mechanism in the land column under both.

## Every script here

The workflow above uses a few of these. The rest are tools you will not find
unless told they exist.

| script | what it does |
| --- | --- |
| `compare_eco_streams.py` | does a run taken in pieces write the same ecological stream as the same run taken whole? Reports structure, interval bounds, payload and headers separately, because they have different owners |
| `build_boundary_conditions.py` | land mask and topography, integrated from the Orogen mesh |
| `build_surface_albedo.py` | background land albedo from lithology, optionally composited with solved lakes |
| `build_surface_roughness.py` | aerodynamic roughness length per cell, surface code 0173 |
| `build_surface_soil_water.py` | feeds pedology's soil water capacity back as `dwmax` |
| `build_stellar_spectrum.py` | this star's spectrum from BT-Settl, checked against the blend at source resolution |
| `sra.py` | writes ExoPlaSim's `.sra` surface format; imported by the builders above |
| `run_exoplasim.py` | prepare, validate and run an experiment |
| `continue_exoplasim.py` | resume a prepared run from its latest restart; `--purpose` says what the segment is for |
| `restart_surface.py` | checks a restart's own copy of every staged surface field against the field it was built from; gates every resume, and `--self-test` proves it can fail |
| `convert_restart.py` | converts a restart across horizontal resolution and real precision onto a template the target executable wrote, and reports what the conversion cost; the result is a NEW INITIAL CONDITION, never a continuation. `--inspect`, `--check-template`, `--dry-run`, and `--self-test` proves every transform against a control that fails |
| `compare_equilibria.py` | asks whether two runs sit at the same equilibrium, against their own year-to-year scatter rather than a picked tolerance: two independent ten-orbit means of one climate differ by about sqrt(2) times the standard error, and the bound is two sigma on that. Reports the pattern difference alongside without thresholding it. For the resolution ladder, SPAT-8 |
| `build_restart_template.py` | cuts a restart TEMPLATE from a short run of the target build: the target's own record set, precision, seed shape and staged static fields, with a clean accumulation window. `--audit` says whether an existing template is clean. `convert_restart.py` refuses a template that is not |
| `restart_format.py` | the one restart parser and writer: names are found by position, not by looking like text. Imported, not run |
| `restart_schema.py` | what every restart record IS. Its inventory is parsed out of the model's own `put_restart_*` call sites; its semantic policy is reviewed and checked in, and `check_policy_covers_source` fails if the model grows a record the policy has not been told about |
| `restart_transforms.py` | the two ways a field crosses a resolution change: triangular spectral projection by `(m,n)`, and a conservative Gaussian remap whose longitude is a fraction of a turn from the cell index and never a degree |
| `restart_convert_selftest.py` | the converter's proof; reached through `convert_restart.py --self-test` |
| `segments.py` | what a run's segments were for and which orbits that makes usable; the one reader of the manifest's `segments` list |
| `finalize_existing_segment.py` | record a completed segment after post-run bookkeeping failed; takes the same `--purpose` |
| `run_stellar_cycle.py` | run or resume a superposed-sinusoid stellar-flux experiment |
| `rebuild_binaries.py` | rebuild every executable and record which patches each contains. `--verify` builds nothing and answers whether the executables on disk are CURRENT -- rule 4's passive half -- and `verify()` is the same answer as data, which `scripts/check_consistency.py` calls rather than restating |
| `verify_model_compiles.py` | does the model SOURCE compile: `gfortran -fsyntax-only` over every translation unit `plasim/CMakeLists.txt` names, sorted into `use` order and run under the flag line `config/planet.yaml` declares, `-fdefault-real-8` included. It builds no binary and answers nothing about whether one is current, which is rule 4's separate question. Run it after a Fortran edit, before paying for a build; `scripts/smoke_test.py` runs it too, and `--skip-compile` opts out |
| `make_profile_bed.py` | short, output-free copy of a run directory to measure the model on, and the one thing in the tree that produces a bed a gate can use. Every run directory on disk was written by an older model, so a verbatim copy carries namelist keys the current `namelist` statements no longer declare and lacks the ones the current model requires; this PRUNES the first, reading what the model declares out of `plasim/src` rather than listing it, and FORCES the second from `config/planet.yaml` -- the `icemod_nl` cold-start profile among them, without which `ice_cold_start` refuses to start. Both are reported in the bed manifest. `--verify-namelist-map` prints the namelist-file-to-group pairing it derives and builds nothing |
| `profile_transforms.sh` | sample every rank with `perf` over a bed; `perf_rank.sh` is its per-rank wrapper |
| `score_transform_profile.py` | fold per-rank `perf` samples into the transform's share of compute |
| `split_legendre.py` | split that share by routine: what maps to a library call and what does not |
| `sweep_compiler_flags.py` | interleaved A/B of compiler flag sets, timed on a bed |
| `bench_ab.py` | interleaved A/B of two executables on one bed, paired per round |
| `stack_floor.py` | the per-thread stack floor at a rung: the heaviest call chain out of the parallel region, summed over its declared local arrays, with `--validate` against a built binary's BSS |
| `verify_filter_fold.sh` | bit-identity check on the Legendre filter fold, filters off |
| `verify_land_column_reduction.sh` | LSHY-3's reduction: the scalar soil water bucket must be the EXACT reduction of the layered land column, not a similar answer. Compiles `landcolumn.f90` -- the model's own kernels, dependency-free and therefore compilable alone -- against a driver that runs both over twenty thousand steps and requires BITWISE equality of store and runoff, at default real and at `-fdefault-real-8` because that is what the model builds with. Three negative controls must break the equality: two layers, free drainage, and the wetness limiter at exponent two |
| `verify_land_column_reduction.f90` | its driver, and the flux sequence is fixed rather than random so a failure is reproducible; it visits overfilling, emptying past zero, exact-capacity arrivals and long dry spells |
| `verify_legendre_parity.py` | checks P and its mu-derivative have opposite parity in legini's own recurrence |
| `compare_restarts.py` | compares two restarts record by record and separates a regrouped sum from a different computation |
| `verify_inverse_transform.py` | lifts `legmod.f90`'s own `sp2fc`, `sp2fcdmu` and `dv2uv` and checks each against the matrix-vector product with the table and the per-mode factor it was handed, one spectral mode at a time, planetary vorticity included as an offset on the vorticity coefficient; the per-mode factor is not unity, so a factor carried to the wrong mode cannot hide, and seven negative controls must each be rejected |
| `verify_weight_factorisation.py` | lifts `legini` out of `legmod.f90`, compiles it, and checks each of the eight weights the transform sites need against the same product formed from scipy's associated Legendre functions and the filter's declared formula; which table pairs with which factor is read back out of the transform routines rather than assumed, and five negative controls perturbing the extracted Fortran must each be caught |
| `verify_weight_factorisation_model.sh` | the other half of that: runs the two-matrix form against the eight-matrix one on one bed at a declared tolerance, with a control that corrupts one per-mode factor and must fail |
| `cache_budget.py` | what the dynamical core cycles through one die's L3, by component and storage class, against the 32 MB target; the class to get right is `sliced`, one copy of which a thread touches only its own rows |
| `probe_grid_contiguity.f90` | decides by ADDRESS whether gfortran copies a thread's band of a shared full-globe grid array at a call boundary, which is what bounds the SHTns change; behaviour cannot tell, because copy-out makes a copied argument look passed-through |
| `verify_threaded_numerics.sh` | the threaded build's grid bands, in two parts: it runs `verify_banded_transform.sh` as the independent side, then runs the model itself against a build whose bands overlap by one latitude row, which must disagree. world-38b left one build, so the growth curve that needed two model arms lives in `verify_shtns_model.sh` instead |
| `verify_banded_transform.sh` | the model's own scatter, weight matrices, banded partial sums and cross-thread reduction against a reference computed outside the model, at bounds declared before the arms run, with three controls: the quadrature weight dropped from the analysis, the bands slid one row, and the scatter handing every thread the same latitudes |
| `verify_banded_transform.f90` | the driver the above compiles against the model's own modules at NPRO threads with OMPSHARED, so the bands are bands; its header says what is subject, what is reference, and what a pass does not cover |
| `banded_transform_reference.py` | the reference that driver reads: associated Legendre functions from scipy, a different formulation from `legini`'s recurrence, on a Gauss-Legendre rule pinned by two identities that need no method -- the weights sum to 2 and every Legendre moment up to degree 2*NLAT-1 integrates to zero |
| `verify_transform_roundtrip.sh` | analysis composed with synthesis must be the identity, which is a right answer rather than a comparison; the one transform check that survives replacing legmod with SHTns, with a control that drops the quadrature weight |
| `verify_transform_roundtrip.f90` | the driver the above compiles against the model's own modules, not a lifted copy; it asserts nfilter is 0, because a filtered round trip returns the coefficient times the filter and is not the identity |
| `probe_shtns_concurrency.f90` | can several threads call SHTns at once on ONE config, which is the assumption the parallel-over-field-and-level design rests on; each thread transforms a spectrum only it knows and every result must be bit identical to the serial answer. Re-run it whenever SHTns is upgraded |
| `probe_shtns_conventions.f90` | drives one spectral mode through legmod and SHTns onto the same grid and divides, SIGNED, to measure the scalar normalisation instead of reading it; a magnitude comparison cannot see the Condon-Shortley phase and the first version of this reported the wrong flag |
| `probe_shtns_vector_conventions.f90` | the same for `SHsphtor_to_spat` against `dv2uv`, which is where the l(l+1) between spheroidal/toroidal potentials and divergence/vorticity lives; zeroes plavor first, since dv2uv would otherwise add planetary vorticity to the mode under test |
| `probe_shtns_forward_conventions.f90` | the ANALYSIS direction against `fc2sp` and `uv2dv`, which is not the inverse's recipe read backwards: the forward filter is `skgpsp` and not `skspgp`, legmod's quadrature weight is explicit where SHTns applies its own, and the l(l+1) between potentials and divergence/vorticity appears again |
| `probe_shtns_analysis_margin.f90` | why the analysis arms miss their bar above T42, tested three ways: band-limiting the field, the FFT length at fixed nlat, and the Gaussian weights themselves. The weights are the answer -- legmod's and SHTns's differ by 1e-10 near the poles, and weights enter analysis and not synthesis |
| `gauss_weight_reference.py` | Gauss-Legendre nodes and weights in extended precision, as an adjudicator. Three double-precision implementations disagreed about the polar weights -- inigau, SHTns and numpy's leggauss -- and double precision cannot settle which is right |
| `verify_gauss_weights.sh` | the model's `inigau` against that reference, with the implementation it replaced as the control. Weights enter the FORWARD transform and nothing else, so an error here is invisible to every check that compares synthesis |
| `profile_memory.sh` | where the model's data comes from: this core's L2, this die's L3, the OTHER die's cache, or DRAM. Zen-specific events, because the generic LLC ones report not-supported on AMD and no uncore PMU is exposed. Two passes, so nothing multiplexes. It is what makes the 32 MB per-die target checkable rather than assumed |
| `count_memops.c` | an LD_PRELOAD interposer counting memcpy, memmove and memset by CALL SITE, for when a profile says a library address is hot and cannot say who reached it. It only sees calls through the PLT, and it is too slow for memset -- a hash probe on every call turns a ten second run into minutes -- so for that one use `perf --call-graph dwarf` instead. Its negative result is what proved the hot addresses were not memcpy |
| `build_model.py` | builds ONE executable and refuses any argument it does not understand. The flag line comes from `config/planet.yaml`, the build directory is per configuration so two configurations cannot delete each other's objects, and Ninja derives the Fortran module dependencies instead of the rule set declaring them by hand |
| `probe_shtns_nspat.c` | how many reals SHTns requires a spatial field to hold, against the NUGP the model allocates. `shtns.h` documents the spatial argument as being `shtns->nspat` long and says the library uses it as its own scratch, so an undersized buffer is a plausible cause of a fault inside a transform. It is not the cause here -- nspat is exactly NUGP at every resolution -- and the probe is kept because that is worth being able to re-establish in one command after any SHTns upgrade or layout change |
| `run_shtns_probe.sh` | builds and runs any of the four probe drivers, in the model's own configuration -- double precision, OMPSHARED, NOPAIRLAT. They were each built by hand once, which made them findable but not runnable |
| `verify_shtns_equivalence.sh` | the gate before any call site moves: the whole conversion recipe on DENSE fields against legmod, scalar and vector, divergence and vorticity separately and together, with a control that drops the Condon-Shortley phase |
| `verify_shtns_equivalence.f90` | its driver; needs SHTns, so pass a prefix or set SHTNS_PREFIX |
| `verify_shtns_model.sh` | the gate the array comparison cannot give: one binary, NSHTNS=0 against NSHTNS=1 on one bed, with the growth curve and a control that drops the spectral filter. The model reaches the transform through a namelist, a build and a call site that comparing arrays never touches |
| `build_shtns.sh` | builds SHTns from a PINNED revision into `vendor/shtns-install`, untracked as `.venv` is; not a subtree because it is used unmodified, not a system install because that puts a build dependency where no manifest records it. Every convention the model reconciles against is a property of the revision |
| `verify_omp_collectives.sh` | checks each of `mpimod_omp.f90`'s 39 routines against the answer written down in advance, with mode-dependent data so an index error cannot cancel |
| `verify_uninitialised_reads.sh` | does anything uninitialised reach a stored value, at the optimisation level this project ships? The trapping form of that question does not arm at `-O2` -- gfortran folds the signalling NaN before the arithmetic -- and `-finit-real=zero` was removed because it hid the answer rather than reporting it. This is the PROPAGATING form, which needs no trap: the declared flag line plus `-finit-real=snan` with FE_INVALID masked, on one bed, restart compared bit for bit against the shipped build's, at NSHTNS=1 and at NSHTNS=0 so a difference the library's own scratch over-read causes is separated from one the model causes. The control is a deliberate uninitialised read reaching the temperature tendency and must break the agreement |
| `verify_snow_mask.sh` | the canopy/snow mask against the one-scalar mask it replaced. The reduction is BITWISE and is the check that can fail: at the shipped values `snowmaskmod` must return the canopy cover itself, so `landmod` evaluates the expression it evaluated before. The other answers are written down too -- a treeless cell and a buried canopy hide nothing, a half-buried one hides half as much, a fully loaded canopy is the snow endmember -- and it reports the plant area index upstream's two masking fractions imply under a gap-fraction reading |
| `verify_snow_mask.f90` | its driver. `snowmaskmod` reads no model state, so this compiles it alone: no grid, no namelist, no build tree. At `-fdefault-real-8`, because the reduction is bitwise and a different working precision is a different check |
| `verify_shared_determinism.sh` | runs the threaded build repeatedly on one bed and requires one answer, since a race on shared state shows up as a different answer each time rather than a wrong one; the control is the pre-fix source and must NOT reproduce itself |
| `_bed_guard.sh` | sourced, not run: refuses a comparison bed that starts cold with a kick and no fixed SEED, because the model seeds that noise from the clock and every comparison on such a bed measures the clock |
| `bench_rank_layout.py` | rank layouts incl. oversubscribed, on latency AND throughput, with pinning verified |
| `thread_count_sweep.sh` | sweeps the OpenMP thread count at one resolution, pairwise against sixteen. The count is compiled in via `num_threads(NPRO)`, so it builds one binary per count first and `OMP_NUM_THREADS` does nothing; NPRO must divide NLAT. `notes/thread-count-by-resolution.md` has the verdict |
| `attribute_barrier_wait.sh` | splits the OpenMP barrier cost per thread and maps each thread to the core and die it ran on, which is what separates "the team is too big" from "the work is unevenly spread". Needs a frame-pointer build (`compile.sh -g`) and records `--sample-cpu`, so the die column is measured rather than inferred from thread creation order |
| `index_runs.py` | index every run by what it is, since a UUID says nothing |
| `assess_convergence.py` | spin-up convergence against the predeclared criteria, over the last `--window` PRODUCTION orbits |
| `close_state_energy.py` | closes the energy budget against the PROGNOSTIC STATE, which is the check the flux diagnostics cannot fail |
| `close_term_energy.py` | closes the model's INTERNAL budget against PlaSim's 28 terms; three identities, and it refuses a low-I/O window |
| `score_rank_bench.py` | scores the rank-layout arms from restart mtimes, against the thresholds `notes/rank-layout-benchmark.md` declared before they ran |
| `reproducibility_matrix.py` | CLIM-44: is the model run-to-run reproducible, and what decides it. A full factorial over ranks, `NOUTPUT` and segment length, three repeats a cell, compared by sha256 with no tolerance to choose. The three hypotheses and the pattern each would show are declared in the module docstring, ahead of any run; `--plan` prints the matrix and runs nothing |
| `diff_restarts.py` | compares two restarts record by record and ranks the differences, because a sha256 answers "are these the same" and nothing else. Handles records that changed LENGTH between builds, which `zsolars` did. The scale to read it against is in the docstring: 1e-14 is arithmetic reordering, 1e-10 after tens of steps is that seed amplified by the model's chaos, 1e-3 is a physics change |
| `verify_restart_continuity.py` | is a run taken in two model calls the same experiment as the same run taken whole? Compares the final restarts of an N-step run and of a `--split` + rest pair, record by record and with no tolerance. `--round-trip` instead takes N steps and then restarts taking NONE, so that every record which differs is state the restart does not carry, named directly rather than inferred from a diverged integration. Refuses a bed that sets `KICK` without `SEED`, because the initial perturbation is then drawn from the system clock and the two cold starts are two different experiments. This is a different property from `reproducibility_matrix.py`'s, which repeats one model call; world-8yyh was reproducible per call and discontinuous at every boundary |
| `co2_overlap_589.py` | CLIM-42's CO2 overlap for the N2O 589 cm-1 band, which sits inside CO2's 15 um band and cannot be carried without one. Adopts a correlated-k band-mean transmissivity from the LMD bundle, corrected onto the window the N2O band occupies, and cross-checks it against Ramanathan (1976) Appendix A over Dickinson's seven 15 um bands. Compares ABSORPTANCES, not transmissivities: the two conventions average over different spectral widths and comparing transmissivities disagrees by 6x for pure bookkeeping reasons. Writes `analysis/co2_overlap_589.json` |
| `trace_gas_band_model.py` | CLIM-42's band model: Donner and Ramanathan (1980) Eq. (1)-(3) with its three band parameter triples, checked against the tables they come from. The CH4 intensity is RECOVERED from the paper's own Table 2 rather than taken on trust, and the per-molecule conversion is checked on the one band both sources carry. `--fortran` prints the `radmod.f90` parameter block so the constants are not retyped; writes `analysis/trace_gas_band_model.json` |
| `reset_restart_accumulators.py` | zeroes the ACCUMULATOR records in a copy of a restart, leaving prognostic state byte-identical. `run_exoplasim.py --restart-from` calls it, so a seeded run does not inherit the donor's partial accumulation window; CLIM-31 |
| `close_ocean_energy.py` | closes the SURFACE budget against the ocean's and the ice module's own output streams, which nothing else here reads; six identities. Reads the per-orbit `MOST_OCEAN.NNNNN`/`MOST_ICE.NNNNN` files the wrapper now preserves; on a run from before that fix it reaches only the final orbit, and the report records which layout it read |
| `predict_ocean_terms.py` | predictions for the two ocean namelist terms of the next baseline bundle, BEFORE the run: the `nhdiff` redistribution and the salinity-derived `TFREEZE` bracket. Validates its diffusion operator against a Laplacian eigenfunction and writes nothing; the write-up is `notes/forcing-bundle-predictions.md` |
| `verify_ocean_flux_channel.py` | verifies the prescribed ocean heat-transport channel, surface code 903 under `oceanmod_nl` `nfluko = 1`, end to end on a field whose answer is stated before the run. `stage` writes the field and the prediction into a run directory; `check` answers each criterion against that run's `ocean_output` stream and exits non-zero if one fails |
| `build_climatology.py` | average an equilibrated segment into climatologies; refuses orbits declared `diagnostic`, and refuses to mix I/O regimes |
| `analyze_climatology.py` | diagnostics, maps and a rate-normalised Koppen interpretation |
| `derive_design_flux.py` | codifies `docs/src/pipeline/state.md` section 5b's flux choice: declared comfort-band thresholds scored per candidate flux from two converged points, with the humidity-coupled variant beside the dry score; writes `analysis/design_flux.json` |
| `analyze_smoke.py` | audit and plot a one-orbit smoke run |
| `bench_pyburn_read.py` | times and VERIFIES pyburn's raw reader on a synthetic output file of one orbit's geometry, because no real raw file survives postprocessing to re-measure against. `--verify` requires every variable to match a reference reader in value, shape and dtype; `notes/audits/pyburn-postprocessing-cost.md` has what it measured |
| `analyze_stellar_cycles.py` | phase-folded response of completed cycle runs |
| `compare_flux_sweep.py` | compare equilibrated reports across a flux sweep |
| `run_albedo_bracket.sh` | drives the four bracket cases end to end, threading each run's announced id into its own continuation; `--self-test` checks that threading without a model |
| `compare_albedo_bracket.py` | compare albedo endmembers and decide whether the bracket resolved |
| `dust_indices.py` | the one place a refractive-index dataset NAME becomes n and k; `aeolian/config/dust.yaml` declares which |
| `dust_optics.py` | two-band mineral dust optics, and the sign of its forcing |
| `mie_dust.py` | Bohren and Huffman Mie code with lognormal size integration |
| `dust_aerofile.py` | writes the model's `aerofile` from the dust optics, and validates the round trip |
| `dust_forcing.py` | prices the dust radiative forcing per surface, shortwave and longwave |
| `build_surface_dust.py` | the prescribed dust column as surface code 1811, read only at `ndustrad = 1` |
| `shortwave_band_weights.py` | integrates the H2O and CO2 band absorptances against this star, for `h2osww` and `co2sww`, and fits the CO2 closed form to Howard's bands. That fit is NOT the one `radmod.f90` runs: PHYS-10 refitted it to the line list, and the script reads the running coefficients out of the model source so the artifact carries both |
| `cloud_band_weight.py` | the CLOUD half of the same correction: Mie over liquid water's k(lambda) to a co-albedo, flux-weighted over range 2 against a 5772 K Sun, for `model.cloud_absorption_scale` which scales `tswr3` and `acl2`. PHYS-11 |
| `spectral_tail.py` | asks whether the damping is absorbing the cascade at the truncation or reaching down into the resolved scales: fits the flow's own inertial range in zonal wavenumber and reports where the spectrum leaves it. Criteria are fixed in the script. Writes `analysis/spectral_tail.json` |
| `transform_exactness.py` | asks whether the model's Gaussian quadrature is exact for the terms `calcgp` forms, by building each one twice from the same coefficients -- once on the model's grid and once on a grid three times finer -- and comparing the global mean. A quadratic and a cubic are carried as controls that can fail. Writes `analysis/transform_exactness.json` |
| `dry_energy_order.py` | fits the dry adiabatic energy sink against the flow amplitude across arms built by `scale_restart.py`, and reports the exponent. Each candidate term enters the energy tendency at a different power of that amplitude, so the exponent names the term; the attribution bands and the fitting window are fixed in the script. Writes `analysis/dry_energy_order.json` |
| `scale_restart.py` | multiplies the atmospheric state in a restart by a factor, holding the global mean and the orography fixed, so a sink can be re-measured at several flow amplitudes. Each candidate term carries a different exponent in that factor, which is what makes the ladder a test rather than a comparison. Writes the scaled restart and a JSON provenance record |
| `stability_probe.py` | finds where the model REFUSES to start and prices a timestep, in minutes rather than orbits: output off, a few hundred steps, cost extrapolated by the measured linearity in step count. Cannot see a LATE blow-up. Writes `analysis/stability_probe.json` |
| `filter_timestep_matrix.py` | sweeps filter strength against timestep on a wall-clock budget, pricing every job before it starts and logging what the deadline refused. Writes `analysis/filter_timestep_matrix.json`. See `notes/physics-filter-stability.md` |
| `build_pyfft.py` | builds the `pyfft` and `pyfft991` f2py extensions `pyburn` needs, for the running interpreter, and verifies them by importing what it produced. `--check` reports whether they load. Rebuild after any Python change; `docs/src/reference/environment.md` says why |
| `shtns_variant_sweep.py` | measures which SHTns on-the-fly variant is fastest per transform type at each rung, forcing each through an authored `shtns_cfg` and confirming what the library loaded. Writes `analysis/shtns_variant_sweep.json`; the adoption threshold is declared in the file. See `notes/shtns-algorithm-selection.md` |
| `corrk_cross_check.py` | checks those two absorptances against correlated-k tables from a modern line list; `--checks` runs the falsifiable checks, `--bands` the per-band CO2 comparison. Needs the LMD Generic PCM bundle, and writes nothing |

## What the flux sweeps established

The measurements are in `notes/parameter-decisions.md` with the terrain and
spectrum each was made on. The one that everything converts through lives in
`lib/sensitivity.py`, which is the only flux-to-kelvin sensitivity this project
has; `scripts/check_consistency.py` recomputes it from the runs it was measured
on. The rest of the error budget is an analysis product rather than a stored
value: regenerate it with `scripts/error_budget.py` once a current climatology
exists. Two results from the sweeps are methodological and outlive any particular
terrain.

**The response is strongly nonlinear across the ice-albedo transition**, so a
flux-to-temperature slope measured on one side of it does not transfer to the
other. Bracket the target between two converged points that span it rather than
extrapolating from a sensitivity measured elsewhere. Doing the latter once
predicted 291.9 K for a run that converged at 287.47 K.

**Equilibration is reported against its predeclared criterion, including when it
misses.** One cold case ran 76 spin-up orbits to a stationary temperature and ice
trend but a mean TOA imbalance of -0.504 W/m2, missing the strict limit by
0.004 W/m2. It is labelled quasi-equilibrated rather than having the rule quietly
relaxed around it, and anything derived from it is a qualitative endpoint rather
than a result of equal confidence. Preserve that standard.

## Two ExoPlaSim behaviours to know before changing anything

Two behaviours of ExoPlaSim worth knowing before changing anything here. Its
`configure()` clears every surface `.sra` when given a landmap, so all land
surface fields except topography and the land mask are uniform namelist defaults;
this is declared via `model.uniform_land_surface` and is deliberate, since Earth's
surface maps are tied to Earth's continents. Albedo is the exception and is
supplied from lithology by `build_surface_albedo.py`, because this world's bare
rock is markedly brighter than ExoPlaSim's uniform default and a large minority
of the land is bright closed-basin fill. **Measure it on the surface the model
sees, not from the rock table**: a lithology change confined to basin fill moves
the VEGETATED land albedo roughly twice as far as the bare one, because
vegetation masks bare-rock variation but not the barren classes. Note this is
substrate albedo -- real vegetation arrives from LPJ-GUESS in loop C of
`docs/src/pipeline/sequencing.md`. Current values are in `world_state.json`.

Its `finalize()` picks output as the last glob match, so a run directory shared
between worlds can silently emit the wrong world's result. **`run_id` is a
UUID** (CLAUDE.md rule 6): a derived identifier separates runs only along the
dimensions it encodes, and a UUID collides with nothing, including along
dimensions nothing here models.

What a run *was* lives in `run_manifest.json`, which gains a `physical` block, and
in `exoplasim/runs/INDEX.json`, generated from those manifests by
`index_runs.py`. That index is tracked even though `runs/` is not, because it is
the only record that survives deleting the output. A continuation must now be
given `--run`; it cannot recompute a name, which is the safer direction.

`model.energy_diagnostics` adds PlaSim's 28-term energy decomposition on codes
360-387. The postprocessor ships 119 codes and none of those, so
`run_exoplasim.py` registers them at run time; patching the vendored tree would
be undone silently by any reinstall of the untracked `.venv`.

Setting it to `2` adds the conversion decomposition, a control for `world-0ov`:
the adiabatic conversion's reference half as the semi-implicit scheme applies it
and as it stands at time t, printed at `ndiag` cadence in `denergy02`'s own
units. Four extra spectral transforms a timestep, for a diagnostic arm.

`model.conversion_time_level` takes that conversion's divergence half back to
time t, so the two halves meet. It CHANGES WHAT THE MODEL INTEGRATES, and it
takes that half out of the semi-implicit treatment in the temperature equation,
so the timestep it is stable at is the explicit gravity-wave one -- 9.5 minutes
at T42 on this planet against a configured 22.5 -- and the model refuses the
setting above that limit rather than integrating something that is not a
solution. `world-0ov`.

`model.robert_filter` sets PNU, the leapfrog time filter's coefficient. Absent
leaves the value ExoPlaSim's own namelist carries. It is here because the filter
is what damps the leapfrog computational mode, and that mode is a candidate
carrier of the adiabatic sink: it alternates sign every step, so it contributes
to the second time difference the sink is built from whatever the timestep is.

**Patched source and per-configuration binaries.** ExoPlaSim compiles a separate
executable for every (resolution, layers, ranks) configuration, so
patching the source
and running rebuilds *only the configuration you are running*. Every other binary
keeps the old code until something asks for it. This is failure class 11 and it
fired three times in a single day.

**The `p` in `most_plasim_t42_l10_p16.x` is THE PARALLEL WIDTH, and the `-p`
flag to the build this replaced was PRECISION IN BYTES.** One letter, two
meanings, and they sat one call apart: `__init__.py` builds the name as
`"most_plasim_t%d_l%d_p%d.x" % (nsp, layers, ncpus)`, while `compile.sh` matched
`-p` against `4`, `8`, `single` and `double` with a default of `prec=4`.

Read the name as threads. `_p16` is the sixteen-thread build, `_p8` the
eight-thread one, and neither says anything about precision -- `config/planet.yaml`'s
`precision_bytes` does, and `rebuild_binaries.py` is what carries it to the
flag. The one case where the two readings agree is `_p8`, because `-p 8` happens
to be a precision `compile.sh` recognises; that coincidence is why the bug in
`notes/audits/compiled-precision.md` hid in the 16-rank builds only, where
`-p 16` matched no case and silently fell back to single precision while the
config declared double.

So, two rules:

- **After any change under `vendor/exoplasim`, rebuild everything**: `python
  exoplasim/scripts/rebuild_binaries.py`. It deletes every executable, rebuilds
  the matrix, and writes `exoplasim/binary_manifest.json` recording the sha of
  every source file each executable was built from. There is no separate
  star-cycle build: the cycle is a namelist switch on these same binaries.
- **`rebuild_binaries.py --verify` is the cheap check** that no executable is
  stale against the source on disk, which is failure class 11 and has fired
  repeatedly. It builds nothing and links nothing: a manifest read and a sha per
  source. It hashes the files rather than reading the subtree commit, so an
  uncommitted edit under `vendor/` is caught rather than waved through, and it
  reports four separate findings -- an executable of unknown provenance, one
  built from source that has moved, a MATRIX row with no binary at all, and a
  toolchain that has moved under unchanged sources. `check_consistency.py` CALLS
  it, through `rebuild_binaries.verify()`, and reports every one of those as a
  failure once anything is built; a tree with no executables at all is unbuilt
  rather than inconsistent and gets a warning. Nothing restates the comparison
  a second time, which is how the copy that used to live in `check_consistency`
  came to see only the first two findings and only in `plasim/run`.

The ozone patch is *resident* in the source -- it must be applied for any build
to be correct. The star-cycle patch is not: it is applied and reversed around its
own build, because a cycle binary and a steady binary are different things and
only one can be in the tree at a time.

**The shortwave gas absorptances are all weighted for the Sun, and all three
corrections live in namelist keys.** Lacis and Hansen give ozone, water vapour
and (in this fork) CO2 as fractions of SOLAR flux, so a K dwarf needs each one
re-weighted: `o3uvw`/`o3visw`, `h2osww`, and `co2sww`. The last is a new absorber
rather than a weight, because `swr` has no shortwave CO2 at all, so its default
is 0.0 and not 1.0. Arguments in `exoplasim/notes/ozone.md`,
`exoplasim/notes/shortwave-water-vapour.md` and
`exoplasim/notes/shortwave-co2.md`; values in `config/planet.yaml`; derivation in
`exoplasim/scripts/shortwave_band_weights.py`. Both were derived from Howard's
1950s band data and both have been checked against correlated-k tables built from
a modern line list, in `exoplasim/notes/corrk-cross-check.md`; the CO2 closed
form was then REFITTED to that line list, so `radmod.f90`'s four coefficients no
longer come from Howard and the patch header records that its own do.
**Two cloud constants and two surface emissivities are namelist keys now, and
one of them is derived rather than set.** `clwhsc` and `clwref` in `rainmod_nl`
are the CCM3 cloud-water e-folding length coefficient and the reference
in-cloud liquid density that `mkclouds` carried as bare literals; `clwhsc`
below zero, which is the default, means DERIVE it from this planet's own
`gascon` and `ga`, because it is a length and the heights it is measured
against are already built that way. `elwland` and `elwsea` in `radmod_nl` are
the surface longwave emissivities `lwr` carried as one literal; they are
declared in `config/planet.yaml` at the values that literal had. `dql` is not
diagnostic: it sets shortwave cloud optical depth and longwave cloud
emissivity, so `clwhsc` moves both. `exoplasim/notes/cloud-water-reference.md`
is what the two CCM3 sources say about `clwref`, the bracket that gets run as
arms around it, and where this fork's shortwave cloud optics parts company with
the scheme those sources describe.

The model source is `vendor/exoplasim`, a git subtree from the `master` branch
of the personal fork, installed editable so the source you read is the source
that compiles. There is no patch stack to keep applied: a model change is a
commit in that directory, and `rebuild_binaries.py` builds whatever is there.
`exoplasim/patches/` holds the changes as they were AUTHORED, with the argument
and evidence for each in its header; its README says which are open as upstream
pull requests.

See `exoplasim/notes/lake-representation.md` for what the model can do with the
endorheic basins, `exoplasim/notes/baseline-equilibration.md` for why the TOA
convergence criterion is failing on an instrument rather than on a state, and
`exoplasim/notes/parameter-decisions.md` for every physical
and format decision
(ExoPlaSim 3.4.2 calendar bugs, postprocessor code quirks, convergence criteria,
Köppen rate-normalisation, sign conventions). Read the notes before changing
anything about how runs are configured — most of the non-obvious choices are
already justified there.

**The completed runs under `exoplasim/runs/` span several eras**, and which era a
run belongs to is not visible in its id, because ids are UUIDs. `INDEX.json` is
the only thing that knows: its `physical` block records the geography, spectrum
and surface albedo each run actually used, and that triple is what decides
whether two runs are comparable. Ask the index; do not infer an era from a name.

What survives a baseline re-run is decided by what a result depends on, not by how
old it is. A flux-versus-temperature slope measured on superseded terrain stays
the best measurement of that slope, because it turns on sea ice and the Planck
response rather than on which basins are bright. A mean surface temperature from
the same run does not survive at all. Judge each quoted number by which of those
it is.

`k2.dat` is the star K2-18, an M2.5V, not spectral type K2. On this nearly
ice-free world the corrected spectrum moves absorbed shortwave by only
0.04 W/m2, because it acts on snow and ice; it is not null on the cold branch,
so it matters for the stellar cycle.

**And no run has yet used `k25v` past its first orbit.** `continue_exoplasim.py`
rebuilt the namelists without `starspec`, so `solarini` fell back to a 4965 K
blackbody: 0.4184 of the flux below 0.75 um against the spectrum's own share, and
+0.024 on broadband snow albedo. Fixed in the drivers, gated on a `radmod.f90`
patch, and tracked as `SPEC-1`/`SPEC-2`. Take band shares from `lib/stellar.py`,
not from a run log.

**Convert a surface albedo change through the atmosphere, or better, measure
the model's own planetary albedo.** Multiplying a surface-albedo delta by full
top-of-atmosphere insolation ignores everything above the surface and overstates
the forcing; roughly half of a surface change reaches the top of the atmosphere,
before any feedback. Predicting a lithology fix that way was a kelvin out.
