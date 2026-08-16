# Eight-Earth-year stellar-cycle experiments

> **Superseded as a description of this world, kept as the only measurement of
> how hard the ocean damps a cycle.** These ran on the 510k-region terrain, at a
> 0.90 baseline, with a blackbody partition rather than a measured spectrum, and
> before the carve. The world now sits near 0.945 on `carved-zoned-v4` with the
> `k25v` spectrum.
>
> The number worth carrying forward is the damping: the 50 m slab reduced the
> swing to **0.239 of the static span**. Treat even that as provisional, because
> damping goes as 1/sqrt(1 + (w*tau)^2) and tau = C/lambda depends on the
> feedback strength, which has moved. The current lambda of 1.31 W/m2/K gives
> tau near 4 Earth years and predicts about 0.30 at an 8-year period. The two
> disagree by more than rounding, and the measured value is the one taken on a
> colder world with far more sea ice, where the feedback is stronger and tau
> longer. Re-measure on the settled world before quoting either.
>
> Both cases below use amplitudes centred on 0.90. The active configuration is a
> single case centred on the current baseline; the control case was dropped.

Two T42 ExoPlaSim integrations were started from the equilibrated 0.90-S-Earth
baseline and run for 62 local orbits (4.024 stellar cycles). The physical orbit
was held fixed at 0.599624 AU and 189.6145 Earth days. Irradiance varies as a
smooth sinusoid at every radiation timestep; both runs use the same fixed
4965 K stellar spectral partition.

| Case | Flux range | Mean T | Local-orbit mean T range | Sea-ice range | Precipitation range |
| --- | --- | ---: | ---: | ---: | ---: |
| Extreme | 0.85--0.95 S-Earth | 280.73 K | 276.68--284.69 K | 3.36--10.58% of planet | 1.915--2.387 mm/day |
| Control | 0.88--0.92 S-Earth | 280.96 K | 279.33--282.64 K | 4.87--7.80% of planet | 2.050--2.236 mm/day |

These are statistics over the final two stellar cycles after averaging each
local orbit, so local seasons do not dominate the ranges. Sea ice is reported
as fraction of the whole planetary surface, as in the model output, rather
than fraction of ocean area.

The seasonally adjusted harmonic response is 7.96 K peak-to-peak in the
extreme experiment and 3.24 K in the control. In both cases the fitted global
temperature maximum occurs about 1.44 Earth years after maximum stellar flux.
This thermal lag and the finite cycle duration keep the planet far from the
separate equilibrium climates found in the static 0.85 and 0.95 experiments.
The extreme run therefore oscillates substantially without alternately
entering the static ice-rich and nearly ice-free endpoint states.

Bright-phase minus dim-phase composites from cycles three and four give an
area-weighted temperature difference of +2.80 K for the extreme run and
+1.22 K for the control. The corresponding planetary sea-ice differences are
-2.01 and -0.83 percentage points, and precipitation differences are +0.129
and +0.045 mm/day. The temperature response is strongest at high latitudes,
consistent with sea-ice feedback. Regional precipitation changes are
heterogeneous: the warmer phase intensifies or shifts circulation rather than
uniformly making every grid cell wetter.

## Periodic-settling assessment

Because 8 Earth years equals 15.409 local orbits, the orbital seasons and
stellar-cycle phase do not repeat together after one stellar cycle. The
appropriate state is quasiperiodic, not an identical year-by-year loop.
Comparison of local-orbit means interpolated to equal stellar phase suppresses
that seasonal aliasing.

For the extreme run, cycle four differs from cycle three by +0.027 K in mean
temperature with a phase-curve RMSE of 0.069 K. The control differs by
-0.130 K with an RMSE of 0.153 K. Mean TOA net radiation over the final two
cycles is -0.463 W/m2 and -0.543 W/m2, respectively. These residuals imply a
small continuing cool drift, especially in the control, but the cycle-to-cycle
changes are small compared with the forced excursions. Four cycles are
sufficient for physically informed worldbuilding interpretation; a longer run
would be appropriate before claiming a strict periodic radiative equilibrium.

## Files

- `stellar_cycle_timeseries.png`: all four cycles, including resolved local
  seasonality.
- `stellar_cycle_phase_response.png`: cycles three and four using local-orbit
  means at equal stellar-cycle phase.
- `stellar_cycle_bright_minus_dim_maps.png`: regional composites around flux
  maximum minus flux minimum.
- `stellar_cycle_report.json`: exact machine-readable statistics.

The full NetCDF fields and restart checkpoints remain in the two corresponding
directories under `runs/`. They occupy about 3.8 GiB per experiment.
