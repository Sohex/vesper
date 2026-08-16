# Superhabitable-world climate workflow

> Figures in this document are illustrative of method, and were measured on
> builds and climates that have since moved. Current values live in
> `world_state.json`, generated from the artifacts. See the convention in
> `CLAUDE.md`.



This component converts the canonical World Orogen geography in `../source/` to
ExoPlaSim boundary conditions and runs reproducible climate experiments. All
commands below are run from the project root.

> **The results in the sections below predate the current geography and are kept
> as a record, not as a description of this world.** The current baseline is
> 0.945 S-Earth, vegetated, on terrain `5bed5549` (`carved-zoned-v4`), and is
> running rather than settled. The flux comes from three converged points on the
> superseded `carved-zoned` terrain -- 0.92/287.47 K, 0.94/291.29 K,
> 0.96/295.15 K, a slope of 192.2 K per unit flux ratio -- corrected by -0.81 K
> for the v4 lithology fix. See `notes/parameter-decisions.md`.
>
> Two further things invalidate numbers quoted here. The stellar spectrum was
> `k2.dat`, which is the star K2-18, an M2.5V, and not a K dwarf at all; it is
> now `k25v`, built for this star's 4965 K. That correction turned out to be
> radiatively null on this warm world, at 0.04 W/m2 of absorbed shortwave, but it
> is not null on the cold branch and matters for the stellar cycle. And every
> carve verdict before the longitude fix integrated each basin's climate from its
> antipode. See `notes/stellar-spectrum-audit.md` and
> `../hydrography/README.md`.
>
> The 0.95 recommendation below belongs to a configuration older than any of
> that: it assumed a uniform 0.22 albedo and a blackbody star, and the albedo
> bracket has since shown the habitable flux depends on the biosphere.
>
> Every run in those sections was built from
> the equirectangular map PNGs of an earlier, 510k-region World Orogen build.
> `../source/` now holds full data exports of a 2.5M-region build with preserved
> endorheic basins and lithology-modulated erosion, including T42/T63/T85
> Gaussian grids emitted directly off the mesh. The results below remain
> physically valid for the geography they used; they are not a description of
> the world as it now stands, and `convert_orogen.py` is superseded by reading
> `../source/exoplasim-T42/planet.nc` directly.

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

## Stellar-flux sweep result

**The T42 flux sweep below is the FIRST era and is superseded on every axis
that matters.** It used the old 510k-region map PNGs, a uniform 0.22 land albedo,
and a blackbody star rather than a measured spectrum; its output has since been
deleted, with the derived products kept in `archive/runs/`. It is retained here
because the *shape* of the response is still instructive -- flux moves surface
temperature and sea ice together, steeply, across this range -- and for nothing
else. Do not quote its temperatures.

| Flux | Climate state | Surface T | Precipitation | Planetary sea ice |
| --- | --- | ---: | ---: | ---: |
| 0.85 S-Earth | quasi-equilibrated cold/ice-rich | 258.80 K | 1.11 mm/day | 32.54% |
| 0.90 S-Earth | strictly equilibrated cool | 280.90 K | 2.15 mm/day | 6.33% |
| 0.95 S-Earth | strictly equilibrated warm | 292.08 K | 2.90 mm/day | 0.264% |

Measured on the superseded pre-carve terrain under a blackbody star. The current
flux calibration is 150.2 K per unit flux ratio, from the T21 bracket, and is
carried in `analysis/error_budget.json`.

At the time, the 0.95 case was the strongest candidate for the intended habitable
world: it reaches the desired 290--293 K range without CO2 tuning. Its broad
land-area interpretation is about 26.8% continental mixed/temperate forest,
10.5% temperate/subtropical forest, 8.8% tropical forest, 8.3% seasonal
tropical woodland/savanna, and 0.7% tundra; however, dry climates remain
substantial (about 16.1% desert and 16.3% steppe/semidesert).

The 0.85 case ran 76 spin-up orbits. Its last ten-orbit temperature and ice
trends were effectively stationary, but its mean TOA imbalance was
-0.504 W/m2, missing the predeclared strict limit by 0.004 W/m2. It is labeled
quasi-equilibrated rather than silently relaxing the rule. Its subsequent
five-orbit climatology has a -0.573 W/m2 mean TOA residual and should be used
as a qualitative cold sensitivity endpoint, not at the same confidence as
0.90 and 0.95.

Endpoint products are under `exoplasim/analysis/climatology/s085/` and
`exoplasim/analysis/climatology/s095/`; the combined report and plot are under
`exoplasim/analysis/sweep/`. Because the response is strongly nonlinear and crosses a
large ice-albedo transition, these three points should not be linearly
interpolated as a precise flux-to-temperature calibration.
