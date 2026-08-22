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

## Every script here

The workflow above uses a few of these. The rest are tools you will not find
unless told they exist.

| script | what it does |
| --- | --- |
| `build_boundary_conditions.py` | land mask and topography, integrated from the Orogen mesh |
| `build_surface_albedo.py` | background land albedo from lithology, optionally composited with solved lakes |
| `build_surface_roughness.py` | aerodynamic roughness length per cell, surface code 0173 |
| `build_surface_soil_water.py` | feeds pedology's soil water capacity back as `dwmax` |
| `build_stellar_spectrum.py` | this star's spectrum from BT-Settl, checked against the blend at source resolution |
| `sra.py` | writes ExoPlaSim's `.sra` surface format; imported by the builders above |
| `run_exoplasim.py` | prepare, validate and run an experiment |
| `continue_exoplasim.py` | resume a prepared run from its latest restart; `--purpose` says what the segment is for |
| `segments.py` | what a run's segments were for and which orbits that makes usable; the one reader of the manifest's `segments` list |
| `finalize_existing_segment.py` | record a completed segment after post-run bookkeeping failed; takes the same `--purpose` |
| `run_stellar_cycle.py` | run or resume a superposed-sinusoid stellar-flux experiment |
| `rebuild_binaries.py` | rebuild every executable and record which patches each contains |
| `make_profile_bed.py` | short, output-free copy of a run directory to measure the model on |
| `profile_transforms.sh` | sample every rank with `perf` over a bed; `perf_rank.sh` is its per-rank wrapper |
| `score_transform_profile.py` | fold per-rank `perf` samples into the transform's share of compute |
| `split_legendre.py` | split that share by routine: what maps to a library call and what does not |
| `sweep_compiler_flags.py` | interleaved A/B of compiler flag sets, timed on a bed |
| `bench_ab.py` | interleaved A/B of two executables on one bed, paired per round |
| `verify_filter_fold.sh` | bit-identity check on the Legendre filter fold, filters off |
| `verify_fold_exactness.sh` | the same with FMA contraction disabled, which is what makes it exact |
| `verify_fold_indexing.sh` | proves each spectral mode gets its OWN filter value, with a negative control |
| `verify_legendre_parity.py` | checks P and its mu-derivative have opposite parity in legini's own recurrence |
| `verify_latitude_pairing.py` | compiles `plasimmod.f90`'s own `ilatperm` and checks the paired decomposition gives each process a mirror pair per slot pair, with two negative controls |
| `verify_paired_decomposition.sh` | runs the paired and contiguous layouts against each other on one bed, plus a wrong-permutation control that must fail |
| `compare_restarts.py` | compares two restarts record by record and separates a regrouped sum from a different computation |
| `verify_symmetric_transform.py` | runs `legmod.f90`'s own `sp2fc` and `sp2fcdmu` against an associated Legendre table computed by scipy, one spectral mode at a time, with a negative control for each |
| `verify_weight_factorisation.py` | checks that `legini`'s eight weight matrices are all P or Q times a per-mode and a per-latitude factor, which is what CLIM-48 rests on, with a negative control |
| `verify_weight_factorisation_model.sh` | the other half of that: runs the two-matrix form against the eight-matrix one on one bed at a declared tolerance, with a control that corrupts one per-mode factor and must fail |
| `cache_budget.py` | what the dynamical core cycles through one die's L3, by component and storage class, against the 32 MB target; the class to get right is `sliced`, one copy of which a thread touches only its own rows |
| `verify_omp_collectives.sh` | checks each of `mpimod_omp.f90`'s 39 routines against the answer written down in advance, with mode-dependent data so an index error cannot cancel |
| `verify_shared_determinism.sh` | runs the threaded build repeatedly on one bed and requires one answer, since a race on shared state shows up as a different answer each time rather than a wrong one; the control is the pre-fix source and must NOT reproduce itself |
| `_bed_guard.sh` | sourced, not run: refuses a comparison bed that starts cold with a kick and no fixed SEED, because the model seeds that noise from the clock and every comparison on such a bed measures the clock |
| `bench_rank_layout.py` | rank layouts incl. oversubscribed, on latency AND throughput, with pinning verified |
| `index_runs.py` | index every run by what it is, since a UUID says nothing |
| `assess_convergence.py` | spin-up convergence against the predeclared criteria, over the last `--window` PRODUCTION orbits |
| `close_state_energy.py` | closes the energy budget against the PROGNOSTIC STATE, which is the check the flux diagnostics cannot fail |
| `close_term_energy.py` | closes the model's INTERNAL budget against PlaSim's 28 terms; three identities, and it refuses a low-I/O window |
| `score_rank_bench.py` | scores the rank-layout arms from restart mtimes, against the thresholds `notes/rank-layout-benchmark.md` declared before they ran |
| `reproducibility_matrix.py` | CLIM-44: is the model run-to-run reproducible, and what decides it. A full factorial over ranks, `NOUTPUT` and segment length, three repeats a cell, compared by sha256 with no tolerance to choose. The three hypotheses and the pattern each would show are declared in the module docstring, ahead of any run; `--plan` prints the matrix and runs nothing |
| `diff_restarts.py` | compares two restarts record by record and ranks the differences, because a sha256 answers "are these the same" and nothing else. Handles records that changed LENGTH between builds, which `zsolars` did. The scale to read it against is in the docstring: 1e-14 is arithmetic reordering, 1e-10 after tens of steps is that seed amplified by the model's chaos, 1e-3 is a physics change |
| `co2_overlap_589.py` | CLIM-42's CO2 overlap for the N2O 589 cm-1 band, which sits inside CO2's 15 um band and cannot be carried without one. Adopts a correlated-k band-mean transmissivity from the LMD bundle, corrected onto the window the N2O band occupies, and cross-checks it against Ramanathan (1976) Appendix A over Dickinson's seven 15 um bands. Compares ABSORPTANCES, not transmissivities: the two conventions average over different spectral widths and comparing transmissivities disagrees by 6x for pure bookkeeping reasons. Writes `analysis/co2_overlap_589.json` |
| `trace_gas_band_model.py` | CLIM-42's band model: Donner and Ramanathan (1980) Eq. (1)-(3) with its three band parameter triples, checked against the tables they come from. The CH4 intensity is RECOVERED from the paper's own Table 2 rather than taken on trust, and the per-molecule conversion is checked on the one band both sources carry. `--fortran` prints the `radmod.f90` parameter block so the constants are not retyped; writes `analysis/trace_gas_band_model.json` |
| `reset_restart_accumulators.py` | zeroes the ACCUMULATOR records in a copy of a restart, leaving prognostic state byte-identical. `run_exoplasim.py --restart-from` calls it, so a seeded run does not inherit the donor's partial accumulation window; CLIM-31 |
| `close_ocean_energy.py` | closes the SURFACE budget against the ocean's and the ice module's own output streams, which nothing else here reads; six identities. Reads the per-orbit `MOST_OCEAN.NNNNN`/`MOST_ICE.NNNNN` files the wrapper now preserves; on a run from before that fix it reaches only the final orbit, and the report records which layout it read |
| `predict_ocean_terms.py` | predictions for the two ocean namelist terms of the next baseline bundle, BEFORE the run: the `nhdiff` redistribution and the salinity-derived `TFREEZE` bracket. Validates its diffusion operator against a Laplacian eigenfunction and writes nothing; the write-up is `notes/forcing-bundle-predictions.md` |
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
| `shortwave_band_weights.py` | integrates the H2O and CO2 band absorptances against this star, for `h2osww` and `co2sww`, and fits the closed form the CO2 patch codes |
| `cloud_band_weight.py` | the CLOUD half of the same correction: Mie over liquid water's k(lambda) to a co-albedo, flux-weighted over range 2 against a 5772 K Sun, for `model.cloud_absorption_scale` which scales `tswr3` and `acl2`. PHYS-11 |
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

**Patched source and per-configuration binaries.** ExoPlaSim compiles a separate
executable for every (resolution, layers, ranks) triple, so patching the source
and running rebuilds *only the configuration you are running*. Every other binary
keeps the old code until something asks for it. This is failure class 11 and it
fired three times in a single day.

**The `p` in `most_plasim_t42_l10_p16.x` is RANKS, and the `-p` flag to
`compile.sh` is PRECISION IN BYTES.** One letter, two meanings, and they sit one
call apart: `__init__.py` builds the name as
`"most_plasim_t%d_l%d_p%d.x" % (nsp, layers, ncpus)`, while `compile.sh` matches
`-p` against `4`, `8`, `single` and `double` with a default of `prec=4`.

Read the name as ranks. `_p16` is the sixteen-rank build, `_p8` the eight-rank
one, and neither says anything about precision -- `config/planet.yaml`'s
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
  repeatedly. `check_consistency.py` runs the same check. It hashes the files
  rather than reading the subtree commit, so an uncommitted edit under
  `vendor/` is caught rather than waved through.

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
a modern line list, in `exoplasim/notes/corrk-cross-check.md`.
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
