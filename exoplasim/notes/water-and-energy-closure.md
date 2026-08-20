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

### It is not the classic spectral leak, though it is in the spectral core

The usual suspect in a spectral model is kinetic energy removed by hyperdiffusion
and never returned as heat. That one is closed. `plasim.f90` calls `mkdheat`
whenever `ndheat > 0`, which is the default, and that routine recomputes the wind
field before and after *both* Rayleigh friction and biharmonic diffusion,
converts the kinetic energy difference to a temperature tendency, and adds it.
Surface-friction dissipation is handled separately in `fluxmod.f90` under the
same switch. Measured, the dissipation returned as heat matches the adiabatic
generation of kinetic energy to -0.0073 W/m2 out of 2.21, so nothing is lost on
that path.

The leak is on the other side of the same step, in the conversion itself rather
than in the friction that closes it, and it runs the other way: the adiabatic
dynamics CREATES energy. "The offset is not a diagnostic" below has the
measurement.

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

Distinguishing those two needs the model to say, not us. It did, on 2026-08-18,
and the first is right: see "The offset is not a diagnostic" below.

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

**What survives.** The melt booking accounts for most of it. On the clean
orbit 66, `hfns - (rss + rls + hfss + hfls)` -- the latent heat of fusion
consumed by snowmelt, which `hfns` correctly includes -- is -0.4222, leaving
`ntr - (rss + rls + hfss + hfls) = -0.7761` unexplained on that orbit; a
gridpoint physics sum of +0.4763 was recorded against it. On the climatology
in the table above the same identities give -0.4097 and -0.9327, but that
climatology's `hfss` and `hfls` carry the corrupt first record, so the clean
ten-orbit block below is the measurement to quote.

## The budget on ten clean orbits, and what carries the gap

Measured 2026-08-17 on orbits 67-76 of `run_8c2e1ff9ab5e`, the first block run
with `NLOWIO = 0` throughout, so the 28 terms and the flux diagnostics are for
once averaged the same way. Mean over ten orbits, with the spread across them.

| | W/m2 |
| --- | ---: |
| TOA net | -0.6039 +/- 0.1576 |
| surface `hfns` | -0.1678 +/- 0.1490 |
| **gap** | **-0.4361 +/- 0.0737** |

So the gap is real and it is near the -0.455 this note was opened on -- but its
spread across orbits is 0.074, not the +/- 0.010 that three earlier runs
suggested. It was never that constant; the earlier agreement was three samples of
a noisy quantity, all taken under the same defective output path.

**Of the four channels -- radiative, sensible, latent, and the melt booking
-- three close; the latent one does not, and its failure is the answer.**

| channel | from the terms | from the fluxes | difference |
| --- | ---: | ---: | ---: |
| radiative, 9 + 10 vs `ntr - (rss+rls)` | -100.657 | -100.634 | -0.023 +/- 0.031 |
| sensible, term 7 vs `-hfss` | +22.415 | +22.411 | **+0.004 +/- 0.026** |
| latent, 11 + 12 + 14 vs `-hfls` | +77.781 | +77.366 | **+0.415 +/- 0.071** |

The sensible line is worth pausing on: under `NLOWIO = 1` that same identity was
out by 1.6 W/m2, and it is now out by 0.004. That is the snapshot artifact
disappearing, measured, and it is why the earlier 28-term reading could not have
answered anything.

**The latent channel does not close, and its partner is the melt term.** The
atmosphere releases 0.415 W/m2 more condensation heating than the surface
supplies as latent flux. Over the same ten orbits,
`hfns - (rss + rls + hfss + hfls)` is **-0.4211 +/- 0.0085** -- the snowmelt
booking. Same magnitude, opposite sign, and far steadier than either the gap or
the latent residual on its own.

That pairing is the mechanism. Water leaves the surface as vapour, is condensed
in the atmosphere releasing latent heat, falls as snow, and consumes the heat of
fusion again when it melts at the surface. Both books are right. What is wrong is
the habit of writing the surface budget as `rss + rls + hfss + hfls`, which omits
the melt term that `hfns` correctly includes -- which is why the four-flux route
gives -0.857 while the gap against `hfns` is -0.436.

**And the whole budget closes once dissipation is counted.** The atmosphere is
also heated by about 1.105 W/m2 of frictional and diffusive dissipation, terms 6,
8, 13, 21, 23 and 25, which has no surface-flux counterpart because that energy
came from the atmosphere's own kinetic energy rather than from below. Summing:

    flux input       -0.857
    dissipation      +1.105
    latent asymmetry +0.415
    -------------------------
                     +0.663   against a term sum of +0.643

Closed to 0.02 W/m2. **The 28-term decomposition is self-consistent.** It was
never going to show a leak, because there is no leak in the decomposition.

## The residual is in the reported TOA net, and the state says so

Measured 2026-08-17 by `exoplasim/scripts/close_state_energy.py`, which is the
first check in this file that could have failed.

Every closure above compares one flux diagnostic with another, and the note kept
running out of road for the same reason each time: `ntr` equals `rst + rlut` by
construction, terms 9 and 10 are the divergence of the very flux profile that
`rlut` and `rls` are read off, and `rst` and `rlut` sum to `ntr` and so cannot be
cross-checked against each other. That family of checks can only ever show the
diagnostics agreeing with themselves. `docs/src/practice/failure-modes.md` class 17.

The model has a quantity the radiation code does not supply: its own prognostic
state. Conservation gives an identity rather than a tolerance,

    d/dt (planetary heat content)  =  mean net TOA radiation

with the heat content assembled from the ocean mixed layer, the sea-ice and snow
mass as latent heat, the twelve-metre soil column, the atmosphere's enthalpy and
its column vapour. If a slab or sea-ice storage term were carrying the reported
imbalance, this would show it at the full size of the imbalance.

**It does not.** On the ten clean `NLOWIO = 0` orbits of the baseline:

| | W/m2 |
| --- | ---: |
| mean TOA net, `ntr` | -0.5995 |
| d(planetary heat content)/dt | -0.0547 |
| **residual** | **-0.5448** |

The most sensitive single reservoir settles it on its own. A sustained 0.6 W/m2
loss would cool the 50 m mixed layer by 0.0769 K per orbit. The measured drift is
**0.0038 K per orbit**, twenty times smaller.

**The offset is structural, and this time that claim is measured across states
rather than across three readings of one.** Each run on its own convergence
window -- the trailing ten orbits, which is what the criterion reads -- plus the
baseline's longer settled block:

| run | window | storage, W/m2 | mean TOA, W/m2 | residual, W/m2 |
| --- | --- | ---: | ---: | ---: |
| `run_8c2e1ff9ab5e` baseline, 289.71 K | 67-76, `NLOWIO = 0` | -0.055 | -0.600 | **-0.545** |
| `run_8c2e1ff9ab5e` baseline | 45-65, `NLOWIO = 1` | +0.057 | -0.500 | **-0.556** |
| `run_b014469b8091` bootstrap, 289.15 K | 81-90 | +0.112 | -0.515 | **-0.626** |
| `run_524fbed77a9a`, 289.00 K | 36-45 | +0.084 | -0.483 | **-0.567** |
| `run_bfa3f5269660`, 281.92 K | 36-45 | +0.010 | -0.552 | **-0.562** |

Three further twenty-orbit windows on the same runs -- 26-45 on both short runs
and 71-90 on the bootstrap -- give -0.628, -0.533 and -0.569. Across all eight
the residual is **-0.573 +/- 0.035 W/m2**, and not one falls outside -0.53 to
-0.63.

The storage column is what makes this a test and not a tautology: it moves by
0.17 W/m2 across the five, correctly reporting a run still climbing as still
climbing, while the residual does not move with it. Widening the two short runs'
windows from ten orbits to twenty raises their storage to +0.279 and +0.125,
because those windows reach back into the approach -- and their residuals do not
follow, staying at -0.628 and -0.533. Each window is reported with two
estimators, a least-squares slope and an endpoint difference, and they agree
throughout: -0.055 against -0.063 on the baseline's clean block, +0.084 against
+0.109 on `run_524fbed77a9a`.

**And the criterion's verdicts do not track the physical imbalance at all.**
`run_524fbed77a9a` is recorded `equilibrated_for_worldbuilding`; it passes
`|mean TOA| < 0.5` at -0.483, and its heat content is rising at +0.084 W/m2. The
baseline fails at -0.600, and its heat content is flat to -0.055. The criterion
passed the run further from equilibrium and failed the one closer to it. Both
storages are small enough that this is not a claim about which run is better
spun up; it is a claim that the quantity being thresholded is not the quantity
the threshold is for.

**And the incoming shortwave is not the offset.** Over a full orbit the time mean
of the inverse-square distance factor is exactly `1/sqrt(1 - e^2)`, so
`<rst - rsut>` has to equal `GSOL0 / (4 sqrt(1 - e^2))` -- a number that comes
from the run's `planet_namelist` and the geometry of an ellipse, and that the
radiation code never sees. Measured on the clean block: **321.5887 against
321.6006, a residual of -0.0119 W/m2.** That tests the disc average, the
eccentricity weighting, the zenith-angle integration and the radiation call
frequency together, and a 0.2% error in any of them would have been the whole
gap. The same identity on this run's `NLOWIO = 1` block is out by -0.0633,
five times worse, which is the size a corrupt first record of twelve would give;
but the 0.91 run's `NLOWIO = 1` window closes to -0.0100, so that is a
suggestion rather than a demonstration and it is not offered as one.

**Where in the column it sits.** Splitting the same measurement at the surface:

| layer | flux route says | state stores | residual |
| --- | ---: | ---: | ---: |
| surface, `hfns` | -0.168 | -0.046 | -0.122 |
| atmosphere, `ntr - hfns` | -0.432 | -0.009 | -0.423 |
| planet, `ntr` | -0.600 | -0.055 | -0.545 |

Four fifths of it is in the atmospheric column, between the two flux
diagnostics, rather than at the surface.

So the residual is named as a quantity: **-0.573 +/- 0.035 W/m2 that the
reported top-of-atmosphere net radiation carries and the planet's heat content
does not experience, four fifths of it in the atmospheric column.** It is not
storage, it is not the insolation, it is not a seasonal sampling artifact, and
the next section says what it is.

## The offset is not a diagnostic: it is energy the model makes

Measured 2026-08-18 by `exoplasim/scripts/close_term_energy.py` on
`run_8c2e1ff9ab5e` orbits 67-76, the only `NLOWIO = 0` block on disk and
therefore the only window on which the 28 terms and the fluxes are averaged the
same way.

### Two of the three candidates were never possible

The section above left `rsut`, `rlut`, or an unbooked atmospheric heating. The
first two are excluded by reading what the accumulators are, before any
arithmetic.

`outmod.f90:2562-2568` builds `assol`, `asthr`, `atsol`, `atthr` and `atsolu`
by adding `dswfl(:,NLEP)`, `dlwfl(:,NLEP)`, `dswfl(:,1)`, `dlwfl(:,1)` and
`dfu(:,1)` every timestep and dividing by `naccuout` at write time. Those are
the arrays the radiation itself wrote. `radmod.f90:1211` forms
`dflux = dlwfl + dswfl`, `landmod.f90:704` drives the surface with
`dflux(:,NLEP)`, and `radmod.f90:1231-1233` builds the atmospheric heating as
the divergence of the same profile. **The reported top-of-atmosphere net is the
number that heats the model, not a rendering of it**, so there is no
reported-against-actual split for a flux diagnostic to be offset in.

That also disposes of the softer version of the candidate. A radiation scheme
that is physically wrong but internally consistent does not produce this
signature; it produces a different climate, in balance. A settled state that
reports a persistent negative net at the top requires an energy SOURCE inside
the model of the same size, and nothing else.

`rsut` is pinned twice over: it drives nothing, and `<rst - rsut>` matches
`GSOL0 / (4 sqrt(1 - e^2))` to -0.0119 W/m2, so it cannot carry half a watt on
its own.

### Three identities, and the one that fails

Each has a right answer and each could have gone the other way.

| identity | right answer | measured, W/m2 |
| --- | --- | ---: |
| gridpoint physics: `denergy 6+7+8+9+10+11+12+13+14+21+22` = `denergy04` | 0 | **+0.0264 +/- 0.0002** |
| adiabatic dynamics: `denergy26 - denergy27` = 0 | 0 | **+0.3446 +/- 0.0129** |
| kinetic energy in steady state: `-denergy27` = `denergy 21+22+23+25` | 0 | **-0.0073 +/- 0.0067** |

The first says the total gridpoint temperature tendency is accounted for by the
parameterisations that book their own heating: there is no large unbooked
heating in the physics half. Term 22 belongs in that sum and is easy to miss,
because it is written as a kinetic-energy expression while `fluxmod.f90:916`
adds exactly it to `dtdt` -- it is the vertical-diffusion frictional heating.

The second is the answer. `denergy26` is the column enthalpy change across
`spectrala`, the adiabatic spectral step, and `denergy27` is minus the kinetic
energy change across the same step (`plasim.f90:3125-3145`, and the loop order
at `plasim.f90:622-661` is what makes that step adiabatic: `gridpointa`,
`spectrala`, `gridpointd`, `spectrald`). Adiabatic dynamics moves energy between
enthalpy and kinetic energy and creates none, so the two must cancel. They do
not. **The spectral core is a net energy source of +0.3446 +/- 0.0129 W/m2**,
which is 15.6% of the conversion it is performing.

The third is what makes the second readable rather than a mis-pairing. In a
settled run the kinetic energy is flat, so the adiabatic generation must equal
the frictional dissipation returned as heat. It does, to -0.0073 W/m2 out of
2.21. Had that failed, `denergy27` would not have been the adiabatic generation
and the second identity would have meant nothing.

One thing the second identity does NOT separate, and the distinction matters to
whoever takes this further. Both terms compare the raw new state against the
Robert-Asselin-filtered old one, because `atm`, `adm` and `azm` are `stm`, `sdm`
and `szm` as the previous step's filter left them (`plasim.f90:2873-2875` and
`2979-2984`). So what is measured is the step AS THE MODEL EXECUTES IT --
semi-implicit leapfrog with the filter's first part folded into the reference
state -- and not the continuous adiabatic equations on their own. The filter
enters both terms the same way, so it cannot be read out of their difference
here. The number is right for the model that produced these runs, which is what
the residual needed; attributing it between the time scheme and the semi-implicit
conversion is a separate question and needs a different experiment.

### The atmospheric budget closes once it is counted

| | W/m2 |
| --- | ---: |
| radiative convergence, `ntr - (rss + rls)` | -100.6349 +/- 0.1362 |
| sensible from the surface, `-hfss` | +22.4110 +/- 0.2228 |
| latent from the surface, `-hfls` | +77.3712 +/- 0.2881 |
| **flux route total** | **-0.8527 +/- 0.0785** |
| latent asymmetry, `(11+12+14) - (-hfls)` | +0.4143 +/- 0.0750 |
| adiabatic non-conservation, `26 - 27` | +0.3446 +/- 0.0129 |
| spectral diffusion of heat, `denergy24` | -0.0099 |
| **total** | **-0.1038 +/- 0.1071** |

against a measured atmospheric enthalpy tendency of **-0.0074 W/m2**, taken as
the trend of `denergy01`, which is the model's own column enthalpy and needs no
reconstruction. The flux route alone was out by -0.85; with the two transfers
counted it is out by -0.10, inside its own scatter.

The latent asymmetry is real physics rather than an error, and it pairs with
the melt term inside `hfns`: convection books the latent heat of sublimation
below freezing, so the atmosphere receives the heat of fusion that the surface
pays back when the snow melts. Its partner, `hfns - (rss+rls+hfss+hfls)`, is
-0.4210 +/- 0.0089 over the same orbits and equals `-ALF * rho * snm` to six
digits.

### What is left

| | W/m2 |
| --- | ---: |
| reported `ntr` | -0.5995 |
| planetary heat content trend | -0.0547 |
| residual | -0.5448 |
| adiabatic source, measured | +0.3446 |
| **still unattributed** | **-0.2002** |

Of that remainder, the atmospheric column carries -0.10 +/- 0.11, which is
consistent with zero, and the surface carries -0.12, which is the second
residual and is treated below.

The result is stable across sub-windows of the block. Trimming the first
orbits, which follow a restart, gives an adiabatic source of 0.3446, 0.3426,
0.3419 and 0.3426 on 67-76, 68-76, 69-76 and 70-76, and a planetary residual of
-0.5448, -0.5475, -0.5614 and -0.5201.

### Why the four runs could not have answered this by themselves

The obvious lever is that the four runs span two flux ratios, so a residual
that is a fixed fraction of a flux should move between them. It cannot be read
that way, and the arithmetic says so in advance rather than afterwards. Taking
the baseline residual as the anchor and rescaling to the 0.910 run:

| if the residual were a fixed fraction of | predicted shift, W/m2 |
| --- | ---: |
| absorbed shortwave | +0.0328 |
| outgoing longwave | +0.0328 |
| incoming shortwave | +0.0202 |
| reflected shortwave | -0.0112 |

against a run-to-run scatter of 0.035 on the residual itself. **Every
hypothesis predicts a shift smaller than the noise**, so the flux-ratio
leverage is exhausted before it starts and no reading of it is evidence either
way. The residual as a fraction of absorbed shortwave runs -0.2374%, -0.2484%,
-0.2606% and -0.2742% across the four; as a fraction of reflected shortwave,
-0.5914%, -0.6096%, -0.5978% and -0.6731%. Neither ordering separates.

The same three identities computed on the three `NLOWIO = 1` windows give an
adiabatic residual of 0.354, 0.439 and 0.343. Those are snapshot readings of an
unaccumulated array and their numbers are not quotable; they are recorded
because they are the same size, and as a suggestion rather than a
demonstration. `close_term_energy.py` refuses such a window unless
`--allow-low-io` is passed.

## The surface residual is in the ocean, and the land is clean

Measured 2026-08-18 on the same block. `hfns` is -0.1678 +/- 0.1571 W/m2 and the
reconstructed surface heat content falls at -0.0461, leaving **-0.1217**. Split
by surface type, with each half compared against `hfns` over the same cells:

| | flux, W/m2 of planet | storage | residual |
| --- | ---: | ---: | ---: |
| land | -0.0273 +/- 0.0582 | -0.0113 | **-0.0160** |
| ocean | -0.1406 +/- 0.1413 | -0.0349 | **-0.1057** |

**The land closes.** That is a test of the reconstruction rather than of the
model, and it could have failed: `landmod` integrates exactly `hfns` into
exactly the soil column that `close_state_energy.py` rebuilds from `SOILCAP`
and the five layer thicknesses, with `ts` standing in for the layer the model
does not write. It did not fail, so the soil capacity, the snow latent term and
the melt booking are all right, and none of them is carrying the residual. The
land figure is -0.016 to -0.027 across the four sub-windows.

**The ocean does not, and the failure is on the cells that never see ice.**
Splitting the ocean by whether a cell carried sea ice or snow at any point in
the window:

| | area fraction | flux, W/m2 of planet | storage | residual |
| --- | ---: | ---: | ---: | ---: |
| never any ice or snow | 0.5201 | -0.1589 | -0.0072 | **-0.1517** |
| ice or snow at some point | 0.0517 | +0.0183 | -0.0503 | **+0.0687** |

Per unit area of the ice-free cells that is -0.2917 W/m2 delivered and
-0.0138 stored. The flux column is exact and sums to the ocean row; the
storage column does not (-0.0575 summed against -0.0349), because storage
here comes from the ten-orbit trend whose four estimators span 0.18 W/m2 (the
next section) -- read it as indicative, not as a partition.

On those cells the slab identity is exact by the model's own construction, and
this is what makes the failure a finding rather than a mismatch of two
reconstructions. `oceanmod.f90:15` sets `NLEV_OCE = 1`, the run's namelist sets
`MLDEPTH = 50`, `nhdiff` defaults to 0 so there is no horizontal ocean
transport, `nfluko` to 0 so there is no flux correction, `NLSG = 0` so there is
no deep-ocean flux, and `ncpl_atmos_ice` and `ncpl_ice_ocean` are both 1 so no
accumulation interval can be mis-divided. `seamod.f90:172` accumulates
`dshfl + dswfl(:,NLEP) + dlwfl(:,NLEP) + dlhfl` and `oceanmod.f90:810` adds it
to a single 50 m layer. So `CRHOS * CPS * mld * d(SST)/dt = hfns`, locally, with
no transport term to hide in.

Two checks say the failure is an offset rather than a scaling. Regressing `hfns`
on `CRHOS * CPS * mld * d(ts)/dt` bin by bin over those cells, across a seasonal
swing of -25.7 to +17.0 W/m2, gives a slope of 1.0946 and a correlation of
0.99435; a centred difference on 12 bins understates the derivative of the
fundamental by 4.5%, which is most of the slope excess, so the heat capacity and
the flux pairing are right. And `hfls` equals `-ALV * rho * evap` on those cells
to 0.000000 W/m2, which is what `fluxmod.f90:687` promises over an ocean
surface and which would have caught a phase-constant error in the latent flux.

The residual is not uniform. Area-weighted by latitude band over the ice-free
ocean it runs -0.52 poleward of 60S, -0.04 between 60S and 30S, -0.64 between
30S and the equator, -0.41 from the equator to 30N, -0.09 between 30N and 60N,
and +3.26 on the 0.0065 of the planet that is ice-free ocean poleward of 60N.
Its per-cell correlation with `hfls` is +0.40 and with `ts` is -0.40.

**So the surface residual is narrowed here, and the section below settles which
side of it carries the gap.** Excluded at this point: the land reconstruction,
the soil heat capacity, the snow and melt booking, the mixed-layer depth, ocean
horizontal transport, a flux correction, a deep-ocean flux, a coupling-interval
mismatch, the sea-ice reservoir (the residual is largest where there is no ice),
and a scaling error in either the capacity or the latent flux. What was still
open was which side of `CRHOS * CPS * mld * d(SST)/dt = hfns` carries the -0.29,
since both were read from the same output stream and nothing in the regular
output reports the slab temperature independently of `ts`. The ocean and ice
streams do, and they are read below.

## The ocean's own books close exactly, and the surface residual is in the annual mean

Measured 2026-08-18 by `exoplasim/scripts/close_ocean_energy.py` on all four
runs. This is the first reading of the ocean and ice streams in this project.

### The instrument exists, and it reaches one orbit per run

`oceanmod.f90:oceanout` writes `ocean_output` and `icemod.f90:iceout` writes
`ice_output`: raw service-format files, an 8-integer header record and one
`NLON*NLAT` float32 record per field, with their own accumulator and their own
interval of `nout` timesteps. Nothing in this project postprocesses them, and
the codes are not the ones the atmospheric table uses. `ocean_output` carries
901 to 906, 910, 939, 972 and 990; `ice_output` carries 701 to 714, 739, 741,
769, 772 and 790 to 796. `osst` and `oheat` appear on 169 and 263 only in
`plasim_dummy.f90:673`, and in pyburn's table those two numbers are the
atmosphere's `tsa` and `hfns`.

**`ice_output` is the better of the two and nothing had named it.** It carries
the whole decomposition of what the atmosphere delivers: 701 the flux as the
ice module received it, 703 the part spent changing the surface temperature,
704 melting snow, 705 melting ice, 706 what is passed on to the ocean, 702 what
comes back from the ocean.

**Both files are truncated at every model call.** `oceanmod.f90:331` opens with
a bare `open(unit,file=...,form='unformatted')`, and so does its counterpart in
`icemod`; the wrapper calls the model once per orbit. So what survives on disk
is the run's LAST ORBIT and nothing else. Orbits 67-76, which every other
measurement in this note uses, were overwritten ten times over.
`close_ocean_energy.py` therefore reads the final orbit, refuses a window whose
records are not there, and refuses again if the stream does not tile the regular
output.

### Six identities, each with a right answer of zero, and all six close

On ocean cells carrying neither ice nor snow in any stream record -- a mask
taken at the stream's own 32-timestep resolution, which is stricter than the
same mask taken from twelve binned records. Per unit area of that mask, against
a seasonal swing of -25 to +16 W/m2 on the baseline:

| | right answer | baseline, orbit 80 | 0.945 bootstrap, orbit 90 | 0.945 short, orbit 45 | 0.910, orbit 45 |
| --- | --- | ---: | ---: | ---: | ---: |
| D1 `xheat` = `rss+rls+hfss+hfls` | 0 | -0.0024 | -0.0093 | +0.0065 | -0.1276 |
| D2 `yheat` = `xcflux` | 0 | 0.0 | 0.0 | 0.0 | 0.0 |
| D3 `xheat - xcflux` = `xsmelt` | 0 | -1.4e-09 | +8.2e-10 | -1.8e-05 | -1.5e-05 |
| D4 `CRHOS*CPS*mld*d(SST)/dt` = `yheat` | 0 | -7.1e-05 | -3.3e-04 | -1.3e-04 | +2.3e-05 |
| D5 `hfns - (rss+rls+hfss+hfls)` = `-ALF*rho*snm` | 0 | +8.5e-08 | -1.2e-08 | +1.2e-09 | -4.3e-08 |
| **D6 `hfns` = `CRHOS*CPS*mld*d(SST)/dt`** | **0** | **+0.0024** | **+0.0096** | **-0.0063** | **+0.1277** |

Stated before the numbers were looked at: any one of these coming out at
-0.29 W/m2 per unit area would have located the residual, and each could have.
D6 is the surface residual itself, evaluated over a window the two streams
share exactly. **It is zero.**

Three things make that a measurement rather than a coincidence of definitions.
The streams' land mask is the regular stream's array, index for index, and the
script refuses if it is not -- `CLAUDE.md` rule 3. The four ocean fields that
must vanish in this configuration do vanish identically: 903 the flux
correction, 904 the vertical diffusion, 905 the horizontal diffusion, 906 the
deep-ocean flux. And `ts` really is the slab temperature: `xts` equals `xsst`
to 0.000 K on these cells, and the snapshot `ts` matches the ocean stream's own
`ysst` to 0.001 K, which is the interpolation error across a 32-timestep gap.

D3 is worth naming because it is the one term the ice module withholds from the
ocean: snow falling into open water is melted on the spot and its heat of fusion
charged against the flux (`icemod.f90:1253`). The atmosphere books the same
quantity as `snm`, D5 confirms `hfns` includes it, and the two agree to five
decimals. So the chain is complete: the atmosphere sends, the ice module passes
on everything but the snow fusion, the ocean receives exactly that, and the slab
integrates exactly what it receives.

**All three of the shapes the residual could have had are excluded.** A term
applied to the slab that the atmosphere does not book: D1 and D3 close. An
ordering or timestep offset between accumulation points: D4 closes, and it is
the integration itself. A unit or area-weighting difference between the two
grids: the masks are the same array.

### What does NOT carry it: the bin weighting, re-measured

This section previously identified the culprit as the first output bin, held to
cover 570 timesteps against 480 for the other eleven, and attributed 54% of the
ice-free ocean's residual to it. **Both the mechanism and the attribution were
wrong, and the error was in this project's own instrument.** Re-measured
2026-08-18; the correction is recorded rather than the original claim, because
a reader grepping for the number should land on the right one.

**Where the 570 came from.** `close_ocean_energy.py` derived the first bin's
span as `total_steps - 11 * ordinary_bin_steps`, having assumed every bin is
tiled by the same 15 stream records. That assumption is what was under test, so
the derivation could only ever return the whole discrepancy as bin 0's excess.

**What pyburn actually does.** It reduces the raw stream to twelve bins with
`np.linspace(0, ntimes, 13).astype(int)` and divides each bin by its own count,
so every bin mean is correct and the remainder of `ntimes / 12` is spread across
the bins by truncation. The run's `burnout` log gives the record counts
directly, and they differ by regime:

| regime | raw records per orbit | records per bin |
| --- | ---: | --- |
| `NLOWIO = 1` | 36 | 3 in every bin, exactly equal |
| `NLOWIO = 0` | 182 | `[15, 15, 15, 15, 15, 16, 15, 15, 15, 15, 15, 16]` |

So the CHEAP spin-up regime bins exactly, and the CLEAN regime a climatology is
required to use is the one whose bins are uneven -- the opposite of what the
original mechanism asserted. Bin 0 holds fifteen records, the FEWEST.

**The bin centres confirm it to the timestep.** At `NLOWIO = 0` the write
interval is 32 steps, so the predicted gaps between centres are
`(c_i + c_i+1)/2 * 32` = `[480, 480, 480, 480, 496, 496, 480, 480, 480, 480,
496]`, and that is exactly what orbits 66 to 76 carry. Orbits 60 to 65 and 77 to
80, all `NLOWIO = 1`, carry a flat 480. `lib/climatology.py` recovers the record
count from the centres this way and raises if no count reproduces them.

**What the correct weights are worth**, on orbits 67-76, against the surface
residual measured over the same orbits. Masks here are taken from the binned
file rather than from the ice stream at its own resolution, so they are the
weaker version of the same masks:

| mask | area | weighting error, W/m2 of planet | surface residual, W/m2 of planet |
| --- | ---: | ---: | ---: |
| land | 0.4282 | -0.0089 | -0.0160 |
| never any ice or snow | 0.5201 | +0.0734 | -0.1517 |
| ice or snow at some point | 0.0517 | -0.0192 | +0.0687 |
| planet | 1.0000 | +0.0452 | -0.1217 |

**The sign is the finding.** On the ice-free ocean the weighting error is
POSITIVE where the residual is negative, so correcting the weights makes the
residual slightly worse rather than explaining half of it. The old table had
this the other way round on every ocean row, and that agreement was the main
evidence for the attribution. It does not survive.

Two smaller things follow from the same arithmetic. The claim that the
coefficient is a property of the call length, with the 0.910 run at 2.8 times
the baseline's, was derived from the same residual and goes with it: what
differs between runs is the record count, and the coefficient follows from
`ntimes mod 12`. And the observation that `NLOWIO = 0` does not remove the
residual is still true and is now unsurprising, since the weighting was never
what caused it.

**A separate truncation that no weighting fixes.** At `NLOWIO = 0` the 182
records cover `182 * 32 = 5824` of the orbit's 5850 timesteps, so the last 26 --
0.44% of every orbit -- are never written at all. An annual mean from a binned
file is a mean over 99.56% of the orbit however it is weighted.

### What is left, and it is the whole of it

Nothing of the -0.292 W/m2 per unit area is attributed: the weighting error
on the ice-free ocean is +0.141 per unit area, the wrong sign to explain any
of it (the table above). Nor can the residual be pinned on the block it was
measured on, because the storage there has to come from a ten-orbit trend and
that trend's estimators disagree by more than half the quantity: the same
10-orbit heat content on orbits 67-76 gives -0.0138 W/m2 as a least-squares
slope of the ten annual means, -0.0372 as a least-squares slope of all 120
bins, -0.0548 as the difference of the first and last annual means, and
-0.1958 as the difference of the first and last bins. A spread of 0.18 on a
quantity being asked for 0.29.

What DOES close is the slab identity itself: D6 says the model's ice-free
ocean conserves energy exactly wherever the question is put to it over a
window the model itself defines. So the residual is unattributed rather than
explained -- a fact about annual-mean bookkeeping against a trend estimator,
not about the slab -- and the Status section says what it waits on.

**The same weighting acts on every annual mean this project takes from a 12-bin
file.** On orbits 67-76 it is worth -0.0255 W/m2 on `hfns` and -0.0653 on
`ntr`, so the planetary residual of -0.545 keeps -0.48 of it and the
top-of-atmosphere half of this note is not overturned. Anything else that
quotes an annual mean of a strongly seasonal field off a binned climatology
carries an error of the same shape, and it is largest exactly where the seasonal
cycle is largest.

**The next instrument is to stop destroying the streams.** `ocean_output` and
`ice_output` should be moved out of the run directory at the end of each model
call, the way `MOST.NNNNN.nc` already is, so that the next climatology block
carries them. Then this closure runs on the block the verdicts are read from
rather than on the one orbit that happened to survive, the ten-orbit storage
trend is replaced by an exact endpoint difference, and the remaining residual
either is there or is not. That costs 200 MB an orbit and no model time.

## What this changes for the convergence criterion

`|mean TOA| < 0.5 W/m2` is applied to a diagnostic that sits 0.57 below the
planet's actual heat tendency, and every recent miss on this baseline is smaller
than that difference. The full argument,
including where the "missing it by 0.0014" came from, is
`exoplasim/notes/baseline-equilibration.md`. Two things from it belong here.

The miss is **not** rounded into a pass and the criterion is **not** amended on
the strength of this measurement: choosing a new quantity for a criterion
immediately after measuring that the new quantity passes is the move
`docs/src/practice/conventions.md` forbids, whatever the physics says.

But the baseline's failure to converge is not a spin-up fault and more orbits
cannot fix it: the criterion reads a quantity that a non-conserving model holds
away from zero at equilibrium. That is what A2 needed to know.

## Status

Narrowed on 2026-08-17 against the settled baseline; the top-of-atmosphere half
named on 2026-08-18, and the surface half on the same day against the ocean's
and the ice module's own output streams.

**Closed.** The large-scale condensation lead, which was the standing candidate:
it is the missing latent heat of fusion in `mklsp`, predicted and measured to
0.7%, worth 1.57 W/m2, and structurally incapable of being the gap. The radiation
diagnostics, which are internally consistent to 0.03 W/m2 on a clean orbit. The
claim that the gap is a fixed -0.455. **A storage term in the slab ocean or the
sea ice**, eliminated against the prognostic state, on eight windows across four
runs. **The incoming shortwave**, closed against the declared solar constant and
the orbit to 0.012 W/m2. **`rsut` and `rlut` as diagnostic offsets**: the
accumulators hold the same arrays that heat the model, so there is nothing for
an offset to be relative to.

**Found.** The instrument was mis-deployed. `denergy` is written unaccumulated
while everything it would be compared against is a time mean, so under
`NLOWIO = 1` -- which is every run made before 2026-08-17 -- the 28 terms are
snapshots and disagree with the fluxes by several percent. That is why turning
the decomposition on did not answer the question. Runs now set `NLOWIO = 0` and
the terms are usable.

**Found, and it is the third candidate.** The adiabatic spectral step does not
conserve enthalpy plus kinetic energy: `denergy26 - denergy27` is +0.3446 +/-
0.0129 W/m2 where it must be zero, and the kinetic-energy identity that would
have refuted the pairing passes at -0.0073. So the reported top-of-atmosphere
net is a true measurement of an atmosphere that is genuinely radiating away
energy the numerics create. With that counted, the atmospheric budget closes to
-0.10 +/- 0.11 W/m2 against the model's own column enthalpy.

**Found, and it is the surface half.** The slab identity does not fail. Read
against the ocean's and the ice module's own output streams, on a window the two
share exactly, `hfns` equals `CRHOS * CPS * mld * d(SST)/dt` to 0.002 W/m2 on
the ice-free ocean of the baseline and to 0.01 or better on two of the other
three runs, against a seasonal swing of 41 W/m2; the four intermediate identities
between the atmosphere's fluxes and the slab's temperature close to five
decimals or better. So the apparent -0.29 W/m2 is not the slab.

**What it is NOT is the bin weighting**, which this note attributed it to until
2026-08-18 and which the section above now refutes with the measurement. The
weighting error is +0.0452 W/m2 of planet and POSITIVE on the ice-free ocean
where the residual is negative, so counting it correctly makes the residual
slightly larger rather than removing half of it.

**Still open, and now the whole of it.** The residual on the ice-free ocean is
unattributed, and it cannot be pinned on orbits 67-76 because the ocean stream
for those orbits was overwritten and the storage there has to come from a
ten-orbit trend whose four estimators already span 0.18 W/m2. That truncation is
fixed as of CLIM-12 -- each model call's streams are moved aside rather than
overwritten -- so the next climatology block carries the whole series and the
trend is replaced by an exact endpoint difference. No run on this build predates
that fix, so the measurement waits on the next one.
That is what is left of CLIM-11: not more model time, and not the offline
radiative transfer this note once proposed.

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
