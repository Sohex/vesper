# Aeolian

Mineral dust: emission, transport and deposition, computed offline from a
climatology rather than inside the climate model. `notes/dust.md` carries the
decision and why it is not in the GCM; this component is the implementation.

```bash
python aeolian/scripts/build_dust.py                          # the baseline source map
python aeolian/scripts/build_dust.py --variant arid_bare_ground   # the bracket on vegetation
```

Products are `analysis/dust_<variant>.json` and a matching `.nc` carrying annual
mean optical depth, deposition flux, emission and the erodible fraction.

## STATUS: the chain is sound, the answer is undetermined by two parameters

The 1700x emission excess is found, and it was not in this component.

**Every orbit of the source run carries a corrupted first output bin.** In
`spd`, `ua` and `va` the lower troposphere is inflated, worsening downward from
1.02x at the model top to **7.5x at the bottom level**. It is systematic rather
than a restart shock: orbits 86, 87, 88 and 90 of `run_b014469b8091` give 7.50,
7.57, 7.54 and 7.58. The snapshot climatology built from the same run is clean,
so it is the regular 12-bin output's first interval specifically.

Emission goes as roughly u* cubed above a threshold, so that one bin was
**100.00% of the annual total**. Excluding it, and with the transport fixed
below, emission falls from 3.4e6 Tg per Earth year to a range that brackets
Earth's ~2000. Nothing in the physics or the unit conversions was wrong; the
driver was.

`flag_anomalous_bins` now detects it and the run excludes it loudly rather than
consuming it. That is a guard against a known defect, not defensive habit.

**Transport now converges, and the reason it did not is worth recording.** It
was neither a CFL violation nor a cycle nor a missing sink. The explicit step is
CFL-limited by the convergence of the longitude grid at the poles, which gave a
14 km cell and a 597 s step, so 4000 iterations covered 27.6 days -- against a
relaxation time of 1351 days for the 0.5 um bin where nothing rains. The
integration was being stopped at 2% of the way. A declared polar `cos(lat)`
floor of 0.2 relaxes the step to 6029 s and the iteration cap is now 40000; all
cases report converged.

### What it says now, and why that is still not an answer

| | z0 = 3e-6 m | z0 = 1e-4 m | z0 = 1e-3 m |
| --- | ---: | ---: | ---: |
| emission, Tg per Earth year | 3381 | 19.7 | 0.0 |
| land-mean optical depth | 0.0597 | 0.0002 | 0.0000 |
| deposition, g/m2 per Earth year | 4.80 | 0.018 | 0.000 |

Earth for scale: about 2000 Tg per year and a land-mean dust optical depth near
0.03. So the smooth end of the roughness bracket is Earth-like and the rough end
is nothing, across a bracket that spans plausible playa surfaces.

The wind tail moves it as much again. Emission in Tg per Earth year against the
Weibull shape:

| shape k | z0 = 3e-6 m | z0 = 1e-4 m |
| ---: | ---: | ---: |
| 1.50 | 77165 | 7235 |
| 2.00 | 29394 | 1465 |
| 2.50 | 14056 | 398 |
| 3.00 | 7885 | 131 |
| 3.98 (measured) | 3469 | 20 |

**The measured 3.98 is an upper bound and the table above is therefore a lower
bound at every row.** It comes from snapshots 5.7 days apart, and averaging
removes variance, which biases the shape parameter high. Earth's near-surface
winds sit at 1.5 to 2.5. So the physically expected range is the top of that
table, not the bottom.

### The reopening test, stated without tuning

`notes/dust.md` reopens the in-model question at a land-mean optical depth above
0.10. At the measured k = 3.98 the answer is **0.0000 to 0.0597, below the
threshold**. At an Earth-like k of 2.0 the same chain gives roughly eight times
the emission at the smooth end, which puts the optical depth around 0.5 and
**well above it**.

So the test is not decidable on what is currently known, and the two parameters
that decide it are the wind tail and the aeolian roughness. Neither is a
free knob to be set by what answer is wanted.

### What would decide it, and what to ask for

The wind tail can be measured rather than fitted. ExoPlaSim writes high-cadence
output and every run directory already contains a `highcadence.nl`. A short
high-cadence segment off the settled baseline would give a real gust
distribution instead of a Weibull fitted to 32 samples.

**Specification: one orbit, sampled every four timesteps.** At the 45-minute
timestep that is 3-hourly, 1464 samples per cell over 183 days, which resolves
the diurnal cycle on a 30-hour day and gives enough independent samples to fit
the tail rather than the body. Bottom-level wind alone is sufficient. That is
roughly 190 MB per variable at T42, and about 25 minutes of model time at the
rate the baseline is running.

It cannot be run until the baseline settles, and the request is recorded rather
than acted on.

## What the component does get right

**The source map, which is the thing the built-in GCM scheme gets wrong.**
`fcoeff * land_mask` emits as much from forest as from salt pan. This one starts
from `substrate_class`, which is consolidated lithology plus closed-basin fill,
so the only unconsolidated material is the fill; weights the two barren classes
separately, because a cemented salt crust is not a silicate soil and Kok's
fragmentation theory does not describe halite cement; removes standing water
from the solved lake extent; and removes snow. It gives **16.4% of land** as
bare erodible ground against a playa fraction of 23.9%, and the difference is
lakes and crust rather than an assumption.

**The physics is grounded rather than recited.** Every constant traces to a
fetched primary source, listed in `config/dust.yaml` with what it is worth:
Kok et al. (2014) for emission, Kok (2011) for the emitted size distribution,
Fecan et al. (1999) for the soil-moisture threshold, Marticorena and Bergametti
(1995) for the drag partition, Sportisse (2007) for below-cloud scavenging.

**Settling uses this world's gravity.** At 12.81 m/s2 a given particle falls
1.31x faster than terrestrial intuition, so dust lifetime here is shorter than
Earth analogues suggest. That is in `settling_velocity` rather than in a comment.

## What is declared rather than derived, and where the brackets are

`config/dust.yaml` carries all of it. The widest are, in order:

1. **The subgrid wind shape**, now measured but only to an upper bound.
2. **The aeolian roughness of the erodible surface**, bracketed 3e-6 to 1e-3 m,
   worth a factor of 40 in emission across that range. The grid-cell roughness
   field is deliberately NOT used: its median over source cells is 0.49 m, and
   feeding a 15 km orographic variance to a scheme built for centimetre-scale
   roughness elements returns zero emission everywhere. Sheltering of a patch by
   the terrain around it is therefore not represented, which biases emission up.
3. **The evaporite erodible weight**, 0.1 with a bracket of 0.0 to 0.3, for
   crust cementation. Declared suppression, not measured efficiency.
4. **Vegetation cover**, which is not modelled at all. Non-barren land is
   assumed to carry a canopy and not emit, following the project's existing
   declared position, which is generous on a world whose median land runoff is a
   few mm per Earth year. `--variant arid_bare_ground` brackets it.

## What it does not do

No dust-climate feedback: the climatology is an input and does not respond.
No vertical structure: a well-mixed column of declared scale height advected by
a single steering wind. No inter-bin microphysics, which is correct for mineral
dust because it neither coagulates nor grows appreciably.
