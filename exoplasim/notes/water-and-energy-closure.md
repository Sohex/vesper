# What closes, what does not, and what was misread

Written after a claim in this project's own state file that ExoPlaSim's runoff
field was defective. It is not. The field is doing exactly what the code says;
it was being read as something it is not. This note records the check, the
correction, and what the same method found about the energy budget.

## The method

Close a budget from the model's own outputs. It costs no model time, and it is
the check that separates "this number is surprising" from "this number is wrong".

## Water: everything conserves

Globally, over the 0.96 climatology, area-weighted:

| | mm per Earth year |
| --- | --- |
| precipitation | 1077.04 |
| evaporation | 1076.95 |
| **P - E** | **0.10** |

One part in ten thousand. And `hfls` matches `evap` through the latent heat of
vaporisation to 0.15%. The water fields are sound.

## But mrro is not local runoff, and the name says it is

Over land, P - E = 892.14 - 724.31 = **167.83** mm per Earth year, which in
steady state has to leave as runoff. The reported `mrro` land mean is **25.35**,
a factor of 6.6 short. Nothing is stored as ice: this run has **no glaciers on
land at all** and peak snow depth is 0.267 m against the 2 m glacier threshold.

The water is not missing. Integrated globally, `mrro` recovers **99.1%** of what
the land budget says was generated:

| | mm per Earth year |
| --- | --- |
| mrro, global integral | 71.257 |
| land runoff x land fraction, 167.83 x 0.4284 | 71.898 |
| recovered | 0.991 |

It has been **moved**. The ocean-only mean of `mrro` is 105.66 mm per Earth
year, and 726 ocean cells carry more than 1 mm per year, which no soil bucket
can produce because there is no soil under the ocean.

`landmod.f90` explains it exactly. `roffstep` calls `mkradv`, which advects
runoff downhill and **modifies its argument in place**:

```
zroff(jlon,:)   = zroff(jlon,:)   - zuroff*zrivin/zarea(jlon)
zroff(jlon+1,:) = zroff(jlon+1,:) + zuroff*zrivin/zarea(jlon+1)
```

So after routing, `drunoff` is local generation minus river outflow plus river
inflow: a net water divergence. Every observation follows from that and none of
it is a defect.

- Negative values, which `AMAX1(0., dwatc-dwmax)/deltsec` cannot produce, are
  cells exporting more river water than they generate.
- Nonzero ocean values are river discharge arriving at the mouths.
- The land shortfall correlating 0.60 with precipitation is wet cells exporting
  the most.
- Global conservation to 99.1% is the routing conserving water, as it should.

**The CF metadata calls it `surface_runoff`, and that is what misled the
reading.** Treat it as river-routed net water flux. For anything that needs water
draining through the *local* profile, weathering above all, use `P - E`.
`pedology/config/pedogenesis.yaml` does, via `weathering.runoff_source`.

The land runoff ratio is therefore **18.8%** against Earth's roughly 35%: drier
than Earth, not extraordinarily so.

## Energy: part of the gap was arithmetic

The same method, applied to the surface energy budget, finds a gap that is not a
gap.

| | W/m2 |
| --- | --- |
| TOA net radiation | -0.356 |
| rss + rls + hfss + hfls | +0.384 |
| reported `hfns` | +0.090 |

The naive sum and `hfns` differ by 0.294. Global snowmelt is 27.81 kg/m2/yr,
whose latent heat of fusion is **0.294 W/m2**. To three digits. So `hfns`
correctly books the energy consumed melting snow and the naive sum does not;
`hfns` is the field to use, and that part of the apparent discrepancy was a
mistake in the reading rather than in the model.

## What genuinely does not close

TOA net of -0.356 against `hfns` of +0.090 leaves **0.446 W/m2**, which matches
the 0.27-0.45 already recorded in `world_state.json` under "energy budget
closure". Fusion does not explain it: snowfall's fusion equivalent is 0.403 and
snowmelt's 0.294, and the difference of 0.109 is the wrong size and would need
snow to be accumulating somewhere it demonstrably is not on land.

### It is structural, not climatic

The decisive observation is that the gap barely moves across runs that differ
enormously:

| run | mean TOA | mean surface | gap |
| --- | --- | --- | --- |
| `convergence`, 0.90 flux, old terrain | -0.4960 | -0.0431 | **-0.4529** |
| `convergence_carved`, 295.18 K | -0.4635 | +0.0021 | **-0.4656** |
| `convergence_s096`, 292.97 K | -0.3374 | +0.1081 | **-0.4455** |

TOA moves by 0.16 W/m2 and the surface by 0.15, but the gap between them holds at
**-0.455 +/- 0.010**, a 2% spread across different fluxes, terrains and mean
temperatures.

That rules out every state-dependent candidate. Sea-ice mass change, snow
accumulation and latent heat of fusion all scale strongly with climate, and a
term that varies by 2% while the climate varies by tens of kelvin is not one of
them. This is a fixed offset.

### It is not the classic spectral leak either

The usual suspect in a spectral model is kinetic energy removed by hyperdiffusion
and never returned as heat. That is closed here. `plasim.f90` calls `mkdheat`
whenever `ndheat > 0`, which is the default, and that routine recomputes the wind
field before and after *both* Rayleigh friction and biharmonic diffusion,
converts the kinetic energy difference to a temperature tendency, and adds it.
Surface-friction dissipation is handled separately in `fluxmod.f90` under the
same switch.

### What it most likely is

An atmosphere genuinely losing 0.45 W/m2 would cool about 1.4 K per Earth year,
given a column mass near 1e4 kg/m2. These runs do not cool. So the energy is not
actually leaving: something that heats the atmosphere is missing from the sum of
`rst`, `rlut`, `rss`, `rls`, `hfss` and `hfls`, or one of those diagnostics is
offset from the others by a constant.

Distinguishing those two needs the model to say, not us.

## What to do about it, and one trap on the way

**PlaSim already has the instrument.** `denergy(NHOR,28)` is a 28-term energy
decomposition, written to output codes 360-387 when `nenergy > 0`, with a 3D
version on 460-487 under `nener3d`. Turning it on for a short segment would name
the term carrying the 0.455 directly. It is not a spin-up; it is a diagnostic
run.

`nenergy` is not exposed by ExoPlaSim's Python API, so it needs a namelist edit,
exactly as this project already does for `STARFILE` in
`run_exoplasim.py:stage_stellar_spectrum`.

**The trap.** `rainmod.f90:524-527` assigns the latent heat constants the wrong
way round:

```
if(zt(jhor) > TMELT) then
  zzal=als      ! sublimation, above the melting point
else
  zzal=alv      ! vaporisation, below it
endif
```

`zt` is the updated layer temperature and `als` is sublimation, so that is
inverted. It sits inside `if(nenergy > 0)` and feeds only `denergy(:,15)`, so it
does not touch the physics and no completed run is affected. But it means that
**the moment anyone turns these diagnostics on to chase this residual, term 15
will be wrong**, by `als - alv` on every condensing gridpoint. Fix it before
reading term 15, or read the other 27.

## Status

Unresolved, but narrowed from "energy does not close" to "a constant -0.455 W/m2
offset between the top-of-atmosphere and surface diagnostics, not climate
dependent, not the spectral dissipation leak, resolvable with a diagnostic run".

It remains larger than the |mean TOA| < 0.5 W/m2 convergence criterion it sits
inside, which is why it is worth resolving rather than merely recording.

## What this cost, and what it saved

An afternoon of arithmetic on existing outputs. It removed a phantom
postprocessing defect from the project's state file, corrected a land runoff
ratio by a factor of 6.6, moved land-mean weathering intensity from 0.19 to 0.50
and narrowed a bracket from 14.6x to 5.6x, and showed that a third of the
recorded energy closure gap was never real.

None of it needed a model run. The general lesson is the cheap one: before
calling a surprising number a defect, close the budget it belongs to.
