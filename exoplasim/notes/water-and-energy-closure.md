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

So that one stands, and now stands better characterised: the surface budget is
internally consistent once fusion is counted, so the residual sits between the
surface and the top of the atmosphere rather than within the surface terms.
Candidates not yet eliminated are sea-ice mass change, sublimation partitioning
between `prsn`, `snm` and `evap`, and whether dissipated kinetic energy is
returned as heat.

It remains larger than the |mean TOA| < 0.5 W/m2 convergence criterion it sits
inside, which is the reason it is worth resolving rather than merely recording.

## What this cost, and what it saved

An afternoon of arithmetic on existing outputs. It removed a phantom
postprocessing defect from the project's state file, corrected a land runoff
ratio by a factor of 6.6, moved land-mean weathering intensity from 0.19 to 0.50
and narrowed a bracket from 14.6x to 5.6x, and showed that a third of the
recorded energy closure gap was never real.

None of it needed a model run. The general lesson is the cheap one: before
calling a surprising number a defect, close the budget it belongs to.
