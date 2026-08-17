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

Those three readings look like a fixed offset, and that reading did not survive
a settled run: see below, where the baseline gives -0.5230 over its climatology
and -0.3539 on a single clean orbit. All three above were made under
`NLOWIO = 1`, whose corrupt first output record contaminates `hfss` and `hfls`
and therefore `hfns`, so their agreement with each other is partly a shared
defect. What the three do still rule out is a strongly state-dependent term: sea
ice, snow accumulation and fusion all scale with climate far more steeply than
anything here moves.

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

Both keys default to false when absent, so nothing changed for the runs that were
in flight when this was written. Both are now set in `planet.yaml`, which means
every run made under it carries the decomposition rather than only a diagnostic
segment.

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

## The instrument was reading a snapshot, and that is why it had not answered

Measured 2026-08-17 on the settled baseline `run_8c2e1ff9ab5e`.

**`denergy` has no accumulator.** Every other field on the regular output stream
is summed across the output interval and divided at write time -- `ashfl`,
`alhfl`, `assol`, `asthr`, `atsol`, `atthr` and the rest, at `outmod.f90:2519`
onward. The 28 energy terms are not: `outmod.f90:1162` writes
`denergy(1,jdiag)` straight from the live array. So under `NLOWIO = 1` the terms
are **instantaneous values at the output timestep** while every flux they would
be compared against is a **time mean over the interval**. Comparing them is
comparing two different averages of two different things.

The size of that is measurable, because one pair is an algebraic identity.
`denergy(:,7)`, the atmosphere's heating from the sensible flux, and `dshfl`,
the surface sensible flux written as `hfss`, reduce to the same expression:
substituting `ztn` from `fluxmod.f90:534` into both gives `term 7 = -dshfl`
exactly. They should agree to rounding.

| | -hfss | term 7 | difference |
| --- | ---: | ---: | ---: |
| orbit 65, `NLOWIO = 1` | 22.4054 | 20.8062 | +1.5992, 7.14% |
| orbit 66, `NLOWIO = 0` | 22.5691 | 22.5351 | +0.0340, 0.15% |

**A factor of 47, on a quantity that is an identity.** Under `NLOWIO = 0` the
model writes instantaneous records and pyburn averages them, so both fields get
the same average and the identity shows through. Under `NLOWIO = 1` it does not.

So the "first look" reading below, taken on `run_b014469b8091` under
`NLOWIO = 1`, is snapshot noise at the several-percent level and its individual
term values should not be quoted. Runs from 2026-08-17 set `NLOWIO = 0` for
unrelated reasons -- the corrupt first output record, `first-output-bin.md` --
and that incidentally makes this instrument work. The alternative fix, an
accumulator for `denergy` in `outmod.f90` beside the others, is CLIM-6 and buys
back the disk rather than the correctness.

## The large-scale condensation lead is closed: it is the missing latent heat of fusion

The standing lead was that terms 11 and 15 fail to cancel while 12 and 16 cancel
to machine precision. That is real, it survives the averaging change, and it is
now named.

**`mklsp` has no ice phase.** `rainmod.f90:362` sets
`zlcpe(:) = ALV/(acpd*(1.+ADV*dq(:,jlev)))` unconditionally, at every level and
every temperature, and the scheme's heating is `ztn = ztn - zdqdt*zlcpe`. So
large-scale condensation always releases the latent heat of vaporisation, even at
210 K where the condensate is ice. Term 15, after the energy-diagnostics patch,
books `ALS` below `TMELT`. The two therefore differ by `(ALS - ALV)/ALV` wherever
condensation happens below freezing, and that is **13.344%**.

That is exactly what the vertical structure shows. The residual as a fraction of
term 11, by model level:

| sigma | T, K | 11 | 15 | sum | sum / 11 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.094 | 202.9 | +0.0411 | -0.0465 | -0.0055 | 13.4% |
| 0.188 | 219.7 | +0.6969 | -0.7898 | -0.0930 | 13.3% |
| 0.297 | 237.6 | +1.7281 | -1.9587 | -0.2306 | 13.3% |
| 0.421 | 252.8 | +2.1544 | -2.4419 | -0.2875 | 13.3% |
| 0.554 | 263.8 | +2.5868 | -2.9200 | -0.3332 | 12.9% |
| 0.691 | 272.2 | +2.6334 | -2.8791 | -0.2456 | 9.3% |
| 0.818 | 279.0 | +2.5851 | -2.8026 | -0.2175 | 8.4% |
| 0.922 | 283.7 | +3.7487 | -3.8854 | -0.1367 | 3.6% |
| 0.983 | 287.0 | +0.9854 | -1.0073 | -0.0219 | 2.2% |

A flat 13.3% at every sub-freezing level, falling away exactly as the level
crosses the melting point. And the prediction closes: 69% of term 11 acts below
273.15 K, so the expected residual is `-0.13344 * 11.8533 = -1.5817` against a
measured **-1.5714**, 0.7% apart.

**It is not the gap and cannot be.** It is nearly constant at -1.55 to -1.59
while the gap swings +/-7 W/m2 through the orbit, the correlation across bins is
0.32, and it is a physics deficiency rather than an accounting error: the model
is self-consistent in using `ALV` throughout, so it loses no energy. What it
loses is the fusion, which it never had.

### What that means for the patch, which is kept

`patches/exoplasim-3.4.2-energy-diagnostics.patch` made term 15 book `ALS` below
freezing. Its justification, recorded above, was that "the same file disagrees
with itself" at `rainmod.f90:731`, `850`, `944` and `1047`. Those lines are real
and they do switch phase -- but they are in `mkrain` and `kuo`, not in `mklsp`.
`mklsp` has no phase switch to be inconsistent with, so the patch did not restore
agreement, it created a deliberate disagreement between the diagnostic and the
scheme it describes.

The patch is kept anyway, because the disagreement now has a **predicted value**
rather than being noise: 11 + 15 should equal `-0.13344` times the part of term
11 acting below `TMELT`, and it does, to 0.7%. So the pair still works as an
audit -- the test is agreement with that prediction rather than cancellation to
zero -- and it additionally measures a real shortcoming of the model, worth
**1.57 W/m2** of atmospheric heating that Earth would have and Vesper's
atmosphere does not. Reverting would hide that and buy nothing.

## What the gap is not, and where it now stands

On the settled baseline, all twelve bins, area-weighted:

| | W/m2 |
| --- | ---: |
| TOA net, `ntr` | -0.6188 |
| `rst + rlut` | -0.6188 |
| surface `hfns` | -0.0958 |
| `rss + rls + hfss + hfls` | +0.3139 |
| **gap, `ntr - hfns`** | **-0.5230** |

**It is not a radiation-diagnostic offset.** `ntr` equals `rst + rlut` to
0.0000. The atmosphere's radiative heating computed inside `radstep` from the
flux profiles, terms 9 + 10, agrees with the same quantity taken from the flux
diagnostics, `ntr - (rss + rls)`, to +0.115 W/m2 on the climatology and +0.028 on
the clean orbit, with per-bin residuals that change sign. Terms 17 to 20 are
nested sub-splits of term 9 -- 17 + 18 = 9 and 19 + 20 = 18, both to five
digits -- and must not be summed with it.

**It is not the fixed constant this note claimed.** The -0.455 +/- 0.010 across
three earlier runs does not hold: the settled baseline gives -0.5230 over the
climatology and -0.3539 on orbit 66. Some of that spread is the corrupt first
output record, which contaminates `hfss` and `hfls` and therefore `hfns` on every
`NLOWIO = 1` run, and the earlier three were all of that kind.

**What survives.** The melt booking accounts for -0.4221 of it:
`hfns - (rss + rls + hfss + hfls)` is the latent heat of fusion consumed by
snowmelt, which `hfns` correctly includes. That leaves
`ntr - (rss + rls + hfss + hfls) = -0.7761` as the quantity still unexplained,
against a gridpointd physics sum of +0.4763.

## Status

Open, and narrowed twice more on 2026-08-17 against the settled baseline.

**Closed.** The large-scale condensation lead, which was the standing candidate:
it is the missing latent heat of fusion in `mklsp`, predicted and measured to
0.7%, worth 1.57 W/m2, and structurally incapable of being the gap. The radiation
diagnostics, which are internally consistent to 0.03 W/m2 on a clean orbit. And
the claim that the gap is a fixed -0.455.

**Found.** The instrument was mis-deployed. `denergy` is written unaccumulated
while everything it would be compared against is a time mean, so under
`NLOWIO = 1` -- which is every run made before 2026-08-17 -- the 28 terms are
snapshots and disagree with the fluxes by several percent. That is why turning
the decomposition on did not answer the question. Runs now set `NLOWIO = 0` and
the terms are usable; exactly one clean orbit exists, orbit 66 of the baseline.

**Open.** `ntr - (rss + rls + hfss + hfls) = -0.7761` on the climatology,
-0.7761 against a gridpointd physics sum of +0.4763. The next step is not a new
instrument but more of the clean one: several `NLOWIO = 0` orbits, so the annual
residual can be separated from the +/-7 W/m2 seasonal storage swing it hides
inside. That is CLIM-1 still.

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

## The 28 terms, first look, and why its numbers are not quotable

Measured 2026-08-17 on `run_b014469b8091`, the 0.945 bootstrap on
`precarve-craton`, at orbits 65, 75 and 85. This was the first run made with
`nenergy` on, and it was made under `NLOWIO = 1`, so every term in it is an
instantaneous snapshot at the output step rather than a mean over it. The
sensible-heat identity above puts that error at 7% of the term. **Individual
values from this reading are superseded and should not be quoted**; what
survives is the structural observation it produced, that 11 and 15 fail to
cancel while 12 and 16 cancel exactly, and that is now closed above.

Two of the 28 are not fluxes and must never be summed with the rest:
`denergy01` is an absolute column enthalpy, of order 1.9e9, and `denergy28` is
the mass-weighted column temperature in kelvin, which `radmod.f90:1290` builds as
`sum(dt * dsigma)`.
