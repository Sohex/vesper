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

Continue a validated experiment from its latest restart (seasonal snapshots are
written by default; `--no-seasonal-output` skips them for a segment, and a run
without them has to be extended before it can produce a climatology):

```bash
python exoplasim/scripts/continue_exoplasim.py --orbits 5
```

Assess a spin-up and create a separate five-orbit seasonal climatology only
after it passes:

```bash
# Run ids are UUIDs and carry no meaning, so start from the index.
python exoplasim/scripts/index_runs.py

python exoplasim/scripts/assess_convergence.py exoplasim/runs/<run_id>
# --run is required: a continuation cannot recompute a name, and being handed
# one cannot silently resolve to a different run.
python exoplasim/scripts/continue_exoplasim.py --run <run_id> --orbits 5
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
| `build_stellar_spectrum.py` | this star's spectrum from BT-Settl |
| `sra.py` | writes ExoPlaSim's `.sra` surface format; imported by the builders above |
| `run_exoplasim.py` | prepare, validate and run an experiment |
| `continue_exoplasim.py` | resume a prepared run from its latest restart |
| `finalize_existing_segment.py` | record a completed segment after post-run bookkeeping failed |
| `run_stellar_cycle.py` | run or resume a superposed-sinusoid stellar-flux experiment |
| `rebuild_binaries.py` | rebuild every executable and record which patches each contains |
| `index_runs.py` | index every run by what it is, since a UUID says nothing |
| `assess_convergence.py` | spin-up convergence against the predeclared criteria |
| `build_climatology.py` | average an equilibrated segment into climatologies |
| `analyze_climatology.py` | diagnostics, maps and a rate-normalised Koppen interpretation |
| `analyze_smoke.py` | audit and plot a one-orbit smoke run |
| `analyze_stellar_cycles.py` | phase-folded response of completed cycle runs |
| `compare_flux_sweep.py` | compare equilibrated reports across a flux sweep |
| `compare_albedo_bracket.py` | compare albedo endmembers and decide whether the bracket resolved |
| `dust_optics.py` | two-band mineral dust optics, and the sign of its forcing |
| `mie_dust.py` | Bohren and Huffman Mie code with lognormal size integration |

## What the flux sweeps established

The measurements are in `notes/parameter-decisions.md` with the terrain and
spectrum each was made on. The flux calibration is an analysis product rather
than a stored value: regenerate it with `scripts/error_budget.py` once a current
climatology exists. Two results from them are methodological and
outlive any particular terrain.

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
`WORKFLOW.md`. Current values are in `world_state.json`.

Its `finalize()` picks output as the last glob match, so a run directory shared
between worlds can silently emit the wrong world's result. `run_id` used to name
everything physical to prevent that -- geography digest, spectrum, flux at
thousandths -- and each of those was added after a near-miss.

**`run_id` is now a UUID, and that is the fix rather than a retreat from one.** A
derived identifier separates runs only along the dimensions it encodes, and the
encoded set is just a list of everything someone has thought of so far. The ozone
band-weight patch changed the physics and moved nothing in it, so the pre-patch
and post-patch runs at the same flux computed the same name and shared a
directory. A UUID collides with nothing, including along dimensions nothing here
models.

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

So, two rules:

- **After any patch, rebuild everything**: `python
  exoplasim/scripts/rebuild_binaries.py`. It deletes every executable, rebuilds
  the matrix, and writes `exoplasim/patches/binary_manifest.json` recording which
  patches are compiled into which sha256. The star-cycle tree is separate and
  needs `build_star_cycle_exoplasim.sh` afterwards.
- **After any `.venv` reinstall, do the same.** `.venv` is untracked and
  reinstallable, and a reinstall restores pristine ExoPlaSim and discards every
  applied patch with no warning at all. `rebuild_binaries.py --verify` is the
  cheap check that tells you whether that has happened; `check_consistency.py`
  runs the same check.

The ozone patch is *resident* in the source -- it must be applied for any build
to be correct. The star-cycle patch is not: it is applied and reversed around its
own build, because a cycle binary and a steady binary are different things and
only one can be in the tree at a time.

See `exoplasim/README.md` for the workflow and results,
`exoplasim/notes/lake-representation.md` for what the model can do with the
endorheic basins, and `exoplasim/notes/parameter-decisions.md` for every physical
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

What survives a re-baseline is decided by what a result depends on, not by how
old it is. A flux-versus-temperature slope measured on superseded terrain stays
the best measurement of that slope, because it turns on sea ice and the Planck
response rather than on which basins are bright. A mean surface temperature from
the same run does not survive at all. Judge each quoted number by which of those
it is.

The stellar spectrum was wrong for three eras -- `k2.dat` is the star K2-18, an
M2.5V, not a K dwarf -- but fixing it changed absorbed shortwave by 0.04 W/m2 on
this nearly ice-free world, because it acts on snow and ice and there is almost
none of either. It is not null on the cold branch, so it still matters for the
stellar cycle. That asymmetry is the general lesson: a correction's size depends
on how much of the surface it acts on, so estimate it against the state you are
in rather than the state it was first measured on.

**Convert a surface albedo change through the atmosphere, or better, measure
the model's own planetary albedo.** Multiplying a surface-albedo delta by full
top-of-atmosphere insolation ignores everything above the surface and overstates
the forcing; roughly half of a surface change reaches the top of the atmosphere,
before any feedback. Predicting a lithology fix that way was a kelvin out.
