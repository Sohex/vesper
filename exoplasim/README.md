# Superhabitable-world climate workflow

This component converts the canonical World Orogen geography in `../source/` to
ExoPlaSim boundary conditions and runs reproducible climate experiments. All
commands below are run from the project root.

It is described here by what it does and how to run it, and it carries no
results. `world_state.json` holds current values, `runs/INDEX.json` records what
each run physically was, and the dated measurements live in `notes/`.

`convert_orogen.py` is superseded by `build_boundary_conditions.py`, which
integrates the mask and topography from the native mesh rather than from the
equirectangular map PNGs.

## Reproduce geography conversion

On Arch Linux the required host tools are `gcc-fortran` and `openmpi`. Activate
`.venv`, install `requirements.txt`, then run:

```bash
python exoplasim/scripts/build_boundary_conditions.py   # land mask + topography
python exoplasim/scripts/build_surface_albedo.py        # background albedo
python exoplasim/scripts/run_exoplasim.py --run-years 1
```

`convert_orogen.py` is superseded. It remapped the equirectangular PNGs, which is
a lossy intermediate now that the fork emits Gaussian grids directly, and its
mask used the `elevation > 0` convention that floods dry closed-basin floors.

The first command writes generated model inputs to `exoplasim/inputs/t42/` and
geography diagnostics to `exoplasim/analysis/geography/`. The second prepares the 0.90-S-Earth
experiment, records a complete manifest, and runs one smoke orbit. It refuses
to overwrite an existing run with climate output. See
`exoplasim/notes/parameter-decisions.md` for physical and format assumptions.

Continue a validated experiment from its latest restart (seasonal snapshots
are omitted during spin-up unless `--seasonal-output` is supplied):

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
python exoplasim/scripts/continue_exoplasim.py --run <run_id> \
  --orbits 5 --seasonal-output
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

## What the flux sweeps established

The measurements are in `notes/parameter-decisions.md` with the terrain and
spectrum each was made on, and the current calibration is in
`analysis/error_budget.json`. Two results from them are methodological and
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
