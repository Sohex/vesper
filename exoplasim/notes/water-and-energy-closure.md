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

### It is the residue of a much larger seasonal swing

Adding seasonal resolution to the climate series changed the picture again. The
gap is not a flat offset within an orbit; per time bin it oscillates hard:

```
orbit 70  +2.23 +6.66 +6.98 +0.95 -7.52 -6.85 -0.79 +0.23 +1.99 -1.39 -5.14 -3.89
orbit 71  +4.70 +7.53 +5.88 -0.30 -7.04 -4.97 -1.82 -0.23 +0.56 -0.50 -4.87 -4.22
orbit 72  +4.30 +6.71 +6.06 +0.34 -6.70 -4.64 -1.46 +1.07 +0.71 -1.55 -6.42 -1.52
```

A swing of about 14 W/m2 peak to peak, repeating closely from orbit to orbit,
against an annual mean of -0.55, -0.44, -0.26. **The seasonal signal is 34 times
the annual residual.**

That is what it should look like. Heat goes into the slab ocean and the
atmosphere through one half of the orbit and comes back out through the other,
while the top of the atmosphere lags, so TOA and surface disagree strongly within
a year and should cancel across it. The interesting quantity is therefore not a
constant leak but **the part that fails to cancel**, which is a very small
difference between two large numbers.

That reframes the search. A 3% asymmetry in a 14 W/m2 seasonal storage term
produces the whole residual, so the candidates are things that would bias one
half of the orbit against the other, rather than things that lose energy
uniformly. Note the bins are equal-length time averages, so their plain mean is
the true annual mean and the residual is not a sampling artifact of the 12-bin
output.

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

`nenergy` is not exposed by ExoPlaSim's Python API, so `run_exoplasim.py` now
edits `plasim_namelist` directly, exactly as it already does for `STARFILE`. Set
`model.energy_diagnostics: true` in `config/planet.yaml`, optionally with
`model.energy_diagnostics_3d`, and the 28 codes are added to the regular output
and recorded in the run manifest.

Both keys default to false when absent and neither is in `planet.yaml`, so
nothing changes for a run in flight and `config_sha256` does not move until
someone chooses.

**The package needs rebuilding for the term-15 fix to take effect**, since the
`.x` binaries are compiled. That is safe to do at any time: each run directory
holds its own copy of the executable, so a rebuild cannot disturb a job already
running, and with `nenergy` at its default of 0 the rebuilt binary is
behaviourally identical to the current one.

**The trap, now fixed.** `rainmod.f90:524-527` assigned the latent heat
constants the wrong way round:

```
if(zt(jhor) > TMELT) then
  zzal=als      ! sublimation, above the melting point
else
  zzal=alv      ! vaporisation, below it
endif
```

`zt` is the updated layer temperature and `als` is sublimation, so that was
inverted. Confirmed on three independent counts before changing anything:

1. **Physics.** Sublimation is the vapour-to-ice transition and applies below the
   melting point; vaporisation applies above it.
2. **The same file disagrees with itself.** The prognostic code at
   `rainmod.f90:731`, `850`, `944` and `1047` all read
   `if(ztnew < TMELT) then zlcp=ALS else zlcp=ALV`, under the comment "update
   constants (ice/water phase)". That is the correct convention, roughly 200
   lines below the diagnostic that inverts it.
3. **So does the surface scheme.** `fluxmod.f90:687` uses
   `where(dt(:,NLEP) > TMELT .or. dls(:) < 0.5)` to select `ALV`, and `ALS`
   otherwise: vaporisation over warm or ocean surfaces, sublimation over frozen
   ones.

`zzal` has exactly four references in the file and all of them are the
diagnostic, so nothing else could have depended on the inverted sense. The fix
is `patches/exoplasim-3.4.2-energy-diagnostics.patch`, a single comparison
operator plus the reasoning as a comment. It compiles clean under the project's
`-fdefault-real-8` promotion.

It sits inside `if(nenergy > 0)`, which defaults to 0, so **no completed run is
affected and the prognostic physics never was.** It would only ever have been
wrong for someone enabling these diagnostics to chase this residual, which is
precisely what is now recommended.

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
