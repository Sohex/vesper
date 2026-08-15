# Superhabitable-world climate workflow

This component converts the canonical World Orogen geography in `../source/` to
ExoPlaSim boundary conditions and runs reproducible climate experiments. All
commands below are run from the project root.

> **These results predate the current geography.** Every run here was built from
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
python exoplasim/scripts/convert_orogen.py
python exoplasim/scripts/run_exoplasim.py --run-years 1
```

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
python exoplasim/scripts/assess_convergence.py \
  exoplasim/runs/t42l10p8_s090_co20450ppm_rot30h_obl32_e020
python exoplasim/scripts/continue_exoplasim.py --orbits 5 --seasonal-output
python exoplasim/scripts/build_climatology.py \
  exoplasim/runs/t42l10p8_s090_co20450ppm_rot30h_obl32_e020 \
  --start-year 46 --end-year 50
python exoplasim/scripts/analyze_climatology.py
```

The completed 0.90-S-Earth baseline equilibrated after output indices 0--45;
indices 46--50 form the independent climatology. Its five-orbit global means
are 280.90 K surface temperature, 2.15 mm/day precipitation, -0.42 W/m2 net
TOA radiation, -0.15 W/m2 downward surface heat flux, and 0.0633 planetary
sea-ice fraction. This is a cool baseline, not yet the intended warmer regime.

Principal products are in `exoplasim/analysis/climatology/`:

- `baseline_regular_climatology.nc`: all retained variables in 12 time bins.
- `baseline_snapshot_climatology.nc`: 32 orbital-phase samples.
- `baseline_classification.nc`: annual fields plus categorical class indices.
- `baseline_climate_report.json`: global diagnostics and area fractions.
- PNG maps for temperature, precipitation, circulation, snow, sea ice,
  hydrology, Köppen--Geiger classes, and broad ecological interpretation.

The Köppen map is explicitly rate-normalized: precipitation rates are
annualized to 365.2425 days before applying empirical Earth thresholds. That
avoids classifying this world's 189.6-day orbital year as artificially dry.
It is a worldbuilding interpretation, not a dynamic vegetation simulation.

## Stellar-flux sweep result

The complete T42 sweep held geography, atmospheric composition, rotation,
obliquity, eccentricity, and model physics fixed:

| Flux | Climate state | Surface T | Precipitation | Planetary sea ice |
| --- | --- | ---: | ---: | ---: |
| 0.85 S-Earth | quasi-equilibrated cold/ice-rich | 258.80 K | 1.11 mm/day | 32.54% |
| 0.90 S-Earth | strictly equilibrated cool | 280.90 K | 2.15 mm/day | 6.33% |
| 0.95 S-Earth | strictly equilibrated warm | 292.08 K | 2.90 mm/day | 0.264% |

The 0.95 case is the strongest candidate for the intended unusually habitable
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
