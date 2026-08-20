# The first output bin of every orbit, and why only the wind blows up

Measured 2026-08-17. The defect is real output, not a misreading, and the
mechanism is now pinned to specific lines with a chain of predictions that could
each have failed and did not.

## What it looks like

In every `MOST.NNNNN.nc` this project produced with `NLOWIO = 1`, the first of
the twelve time bins carries wind and humidity that do not belong with the other
eleven. Bottom-level wind speed ratio between bin 0 and the rest: 7.50, 7.57,
7.54, 7.58 on `run_b014469b8091` orbits 86 to 90, and 7.97 on
`run_bfa3f5269660`. Systematic, not a restart shock.

Every scalar is clean to better than 2%: `ts`, `tas`, `pr`, `evap`, `mrso` and
`ps` all sit within 0.4 to 1.8% of the other bins. Results that do not read a
wind are unaffected, and that covers the flux bracket, the soil, surface code
229 and the convergence assessment.

Mean absolute eastward wind by model level, bin 0 against bins 1 to 11, on
`run_b014469b8091` orbit 86:

| level | bin 0 | bins 1-11 | ratio |
| ---: | ---: | ---: | ---: |
| 0 | 18.86 | 18.21 | 1.04 |
| 1 | 21.08 | 20.19 | 1.04 |
| 2 | 19.60 | 18.64 | 1.05 |
| 3 | 21.65 | 14.74 | 1.47 |
| 4 | 25.26 | 10.77 | 2.35 |
| 5 | 28.27 | 7.69 | 3.68 |
| 6 | 30.82 | 5.56 | 5.54 |
| 7 | 32.89 | 4.32 | 7.62 |
| 8 | 34.32 | 3.77 | 9.10 |
| 9 | 34.44 | 3.09 | 11.16 |

The three top levels agree to 5% and everything below level 2 is wrong by a
factor that grows downward. That break between level 2 and level 3 is the whole
clue, and it is not a boundary-layer boundary. It is a **restart transfer
boundary**, and where it falls is set by the MPI rank count.

## The mechanism

Three defects in PlaSim's low-I/O path (`NLOWIO = 1`, the default set at
`plasimmod.f90:151`). The first causes almost all of it.

### 1. The spectral output accumulators are saved through the gridpoint restart path

`plasim.f90:857` and `plasim.f90:1086` pass `aaso`, `aasp`, `aast`, `aasqout`,
`aasd` and `aasz` to `mpputgp` / `mpgetgp`. Those six are spectral arrays,
`(NESP,NLEV)`. `mpgagp` and `mpscgp` (`mpimod.f90:130` and `:149`) declare their
dummy argument `pp(NHOR,klev)`, so they move `NHOR` elements per level out of an
array whose levels are `NESP` long.

For T42, 10 layers, 16 ranks: `NHOR = NLON*NLAT/NPRO = 512`,
`NESP = ceil(NRSP/NPRO)*NPRO = 1904`. The transfer covers the first
`NHOR*NLEV = 5120` elements of a `NESP*NLEV = 19040`-element array:

- model levels 1 and 2 whole,
- level 3 through its 1312th coefficient, which is zonal wavenumber m = 18,
- levels 4 to 10 not at all. They come back zero.

Confirmed directly in `run_b014469b8091/plasim_restart`: the 8192 values stored
per level for `aasd`, `aasz`, `aast` and `aasp` are sixteen identical copies of
one 512-element block, which is what gathering a replicated array through
`mpi_gather` produces. The planetary-vorticity coefficient of `aasz`, which sits
at element 3 of every level, appears in the saved data at flat offsets 2, 1906
and 3810: stride 1904, three times only. Its values are 152.50, 155.33, 155.24,
which are 90 times the live `sz(3,jlev)` of 1.6941, 1.7264, 1.7248 in the same
file.

### 2. `naccuout` survives whole, so the divisor does not match the accumulators

`naccuout` is restored correctly (`plasim.f90:767` and `:1002`), and it is not
zero. A model call here is `N_RUN_STEPS = 5850` steps with output every
`NSTPW = 160`, so the last write is at step 5760 and 90 further steps accumulate
without ever being written. `mod(90,160) /= 0`, so 90 goes into the restart, and
`plasim_restart` holds exactly `naccuout = 90`.

The next call reaches its first write at step 160 with `naccuout = 90 + 159 =
249`, because `outsp` runs before that step's `outaccu`. Levels 1 and 2 hold 249
samples and are right. Levels 4 to 10 hold 159 and are divided by 249, so every
field on them is scaled by

    delta = 159 / 249 = 0.638554

**The counter is not the bug on its own.** It is restored correctly and it is
the right count for the accumulators it was saved with. The gridpoint
accumulators -- `aprl`, `ashfl`, `aevap`, `aadt` and the rest -- go through the
same `naccuout` and through `mpputgp`, which is the right routine for them, and
they are clean in the first record. Same counter, same divisor, correct answer.
That is the discriminator: what fails is the six arrays whose type does not
match the routine.

### 3. The planetary vorticity turns a 36% deficit into a 160 m/s jet

`outmod.f90:394` subtracts the whole planetary vorticity from the divided
accumulator:

    aasz(3,jlev) = aasz(3,jlev) - plavor

`sz` carries absolute vorticity and `plavor = sqrt(8/3)` is its (n=1,m=0) mode.
On levels 4 to 10 the accumulator only holds `delta*plavor`, so the written
relative vorticity keeps a residual of `-(1-delta)*plavor` in that one
coefficient. That coefficient is exactly solid-body rotation, so the record
carries a retrograde zonal wind

    u = -(1 - delta) * Omega * a * cos(phi) = -160.8 * cos(phi) m/s

identical at every affected level, zonal only, with no vertical structure of its
own. This is why:

- only the wind blows up, and only `ua`: `plavor` is in the vorticity and
  nowhere else, and solid-body rotation has no meridional component;
- the ratio worsens downward: the spurious wind is the same at every affected
  level while the real wind it is added to falls from 15 to 3 m/s;
- the profile looks like a missing boundary layer. It is not. Nothing about
  surface drag is involved.

Temperature escapes for a different reason: `aast` holds the deviation from the
reference profile `t0`, and `outsp` adds `t0(jlev)*ct` back after the division,
so `delta` only scales the deviation. `aasp` likewise carries `ln(ps)` with
`log(psurf)` added at write time.

### 4. Separately: `sqout` is stale, and its phase walks

`plasim.f90:3218` recomputes `sqout` only when `mod(nstep,nafter)==0`, but
`outmod.f90:2556` accumulates it on every timestep, and `sqout` is not in the
restart. So at the start of each model call it is zero until the first step whose
absolute `nstep` is a multiple of 160, and `aasqout` sums zeros until then. A
call is 5850 steps and `5850 mod 160 = 90`, so the phase advances 90 each orbit
and repeats every 16 orbits. This is on top of the `delta` scaling, and it is why
humidity is worse than `delta` alone predicts.

## What was measured, and what would have refuted it

Each of these was stated before it was looked at.

**Orography is constant in time, so bin 0 must differ from bins 1-11 in `sg` if
and only if the accumulator is divided by a divisor that does not match it.**
Bins 1 to 11 of `sg` are bit-identical to each other in all four runs. Bin 0 is
not.

**The `sg` anomaly must be zero for the spectral coefficients that survive the
restart and exactly `(delta-1)/3` for those that do not** -- one third because
pyburn averages three model records into each bin. The 512-element cut falls
between m = 5 (486 reals) and m = 6 (560), so the prediction is: nothing at
m <= 5, a mixture at m = 6, a single exact number at m >= 7.

| zonal wavenumber | measured diff/ref |
| --- | ---: |
| m = 0 to 5 | 7.9e-10 of the anomaly variance, i.e. float32 zero |
| m = 6 | -0.0209, and not a uniform scaling within the block |
| m >= 7 | -0.120482, relative residual 3e-7 at every latitude |

`3*(-0.120482) + 1 = 0.638554 = 159/249`, against `naccuout = 90` read out of the
restart file. Six significant figures, from a divisor nothing else could have set.

**A run with a different `N_RUN_STEPS` must give a different number.**
`run_bfa3f5269660` is 6018 steps, so 37 records, 98 left over.
`159/(98+159) = 0.618677`, so the ratio must be -0.127108 rather than -0.120482.
Its restart holds `naccuout = 98` and its measured ratio is -0.127108.

**The first call of a run has nothing left over to restore, so its bin 0 must be
clean.** Orbit 0 of all four runs has an `sg` anomaly of exactly zero.

**A call restarted from a differently-sized call must carry the other call's
number.** `run_524fbed77a9a` orbit 1 measures -0.127108, the 6018-step value,
while orbits 2 to 45 measure -0.120482. The ratio in orbit N is set by the length
of the call that wrote the restart, which is what the mechanism says and nothing
else does.

**The wind anomaly must be a level-independent `cos(phi)` solid-body rotation of
amplitude `-(1-delta)*Omega*a/3 = -53.59 m/s` on levels 3 to 9, and absent on
levels 0 to 2.** Measured on `run_8c2e1ff9ab5e`, bin 0 of orbits 60-65
(`NLOWIO = 1`) against bin 0 of orbits 66-71 (`NLOWIO = 0`, same window), after
removing the `-(1-delta)/3` deflation of the real wind:

| level | fitted cos(phi) amplitude, m/s | residual rms, m/s |
| ---: | ---: | ---: |
| 0 | +2.50 | 2.16 |
| 1 | +3.88 | 2.74 |
| 2 | +3.71 | 2.29 |
| 3 | -53.44 | 1.27 |
| 4 | -53.46 | 1.03 |
| 5 | -53.52 | 0.79 |
| 6 | -53.56 | 0.64 |
| 7 | -53.54 | 0.60 |
| 8 | -53.55 | 0.69 |
| 9 | -53.45 | 0.60 |

The amplitude is from `plavor = sqrt(8/3)`, the rotation rate and the radius,
with no fitted constant. A wrong sign, a wrong magnitude, a shape that was not
`cos(phi)`, or a level dependence below level 2 would each have killed it.

**Humidity must carry a further error whose size repeats every 16 orbits,**
because `5850 mod 160 = 90` and `gcd(90,160) = 10`. Predicted bin-0 ratio
`((160-k+1)/249 + 2)/3` with `k` fixed by `nstep` at the start of each call,
against measured bottom-level `hus`, over orbits 1 to 90 of `run_b014469b8091`:
correlation 0.946, mean absolute difference 0.028. The measured values sit
systematically about 0.025 low, which is the same staleness acting on records 1
and 2.

## What was eliminated

- **`naccuout` persisting across the restart, on its own.** It is restored
  correctly and matches what it counted. The gridpoint accumulators share it and
  are clean. It becomes harmful only because six of its partners are truncated.
- **A short first output interval.** The first interval is a full 160 steps like
  every other one. What is short is the sample count on levels 4 to 10, and the
  cause of that is the restart transfer, not the interval.
- **pyburn and the 36-to-12 binning.** Bins 1 to 11 of a constant field are
  bit-identical and bin 0 is not, and the anomaly is confined to spectral
  coefficients beyond the 512th. Binning has no way to produce that.
- **The vector transform from divergence and vorticity.** Same evidence: it
  cannot act on a time-constant orography field.
- **`nkits`, the kick-start.** It is forced to 0 on any restart
  (`plasim.f90:179`), and every affected call is a restart.
- **Any `cos(phi)` factor in `ua` and `va` on sigma levels.** Refuted earlier by
  pointwise measurement, see `notes/failure-modes.md` class 15. The `cos(phi)`
  in this note is a different thing: an additive solid-body wind in one record,
  not a multiplicative factor on the field.

## Still wrong, and not the subject of the fix

- `outsp` runs before that timestep's `outaccu` while `outgp` runs after, so the
  spectral stream is a mean over 159 of the 160 steps in its interval and the
  gridpoint stream is a mean over all 160. Each is self-consistent; the two
  cover slightly different windows. A constant 1-in-160 effect in every record.
- `outmod.f90:274` divides `arasc` by `naccuout` and then writes `rasc`, the
  instantaneous value. Code 54 under low I/O is not the accumulated quantity it
  claims to be.

## Untested prediction

The level at which the profile breaks is `NHOR*NLEV/NESP` counting from the model
top, and `NHOR` is `NLON*NLAT/NPRO`. At T42 with 8 ranks rather than 16,
`NHOR = 1024` and `NESP = 1896`, so about 5.4 levels would be restored and the
break would sit between levels 5 and 6 instead of between 2 and 3. Nothing in
this project has run that configuration with `NLOWIO = 1`, so it is a prediction
and not a measurement.

## What to do about it

**Read every existing run in `exoplasim/runs/` with the corrections above.**
All of them were written by binaries without the patch, so their `NLOWIO = 1`
orbits carry the bin-0 defect in full and their `NLOWIO = 0` orbits carry the
CLIM-20 rank-divergence corruption; see
`notes/audits/nlowio-collective-deadlock.md`.

**For new runs, use `NLOWIO = 1` for spin-up and `NLOWIO = 0` for any orbit a
climatology will be built from.** The bin-0 defect is fixed, so low-I/O spin-up
output is now trustworthy beyond the scalars. What is NOT fixed is the
accumulation itself: `NLOWIO = 0` writes instantaneous samples and `NLOWIO = 1`
writes interval accumulations, and a sample set is strictly more information,
since a mean can be recovered from it and an accumulation cannot be undone.
Anything reading variance, extremes or single records -- DUST-5's gust
distribution above all -- needs the samples. The cost of the clean regime is now
small: measured 2026-08-18 with the postprocessor fixed, an orbit is 86.8 s at
`NLOWIO = 1` against 110.4 s at `NLOWIO = 0`, a ratio of 1.27 where it used to be
3.7 (`notes/audits/pyburn-postprocessing-cost.md`).

The measurements below are of UNPATCHED behaviour, on `run_8c2e1ff9ab5e`:

| | first-bin wind ratio | binned `spd` / snapshot mean-of-speed |
| --- | ---: | ---: |
| orbit 65, `NLOWIO = 1` | 9.00 | 1.141 |
| orbit 70, `NLOWIO = 0` | 1.03 | 1.002 |
| orbit 76, `NLOWIO = 0` | 0.99 | 0.999 |

That second column corrects an earlier conclusion. The binned `spd` under
`NLOWIO = 1` sits between the speed of the time-mean vector and the mean of
instantaneous speeds, and that was written up here as a property of the binning
that would outlive the switch. **It does not.** The cancellation happens in the
model's output accumulation, and with the accumulation off the binned field is
the mean of the speed to 0.2%.

So the corrections written against it are gone rather than improved:

- `carve_verdict.py:turbulent_forcing` read wind from the snapshot product and
  dropped the corrupt record from humidity. It now reads both from the binned
  climatology, which is strictly better at 182 samples an orbit against 32.
- `build_dust.py:speed_bias_correction` scaled the binned wind up by a per-cell
  factor with a median of 1.554. On clean output that would be a 55% error.

**Both corrections are wrong once the defect is gone, so the guard refuses
instead.** `build_climatology.py` stamps `low_io` on every product it writes and
`lib/paths.py:require_clean_io` raises on a tainted one. A product with no stamp
counts as tainted, because every climatology built before 2026-08-17 was.

**Do not drop the first bin from anything else.** Every scalar is within 2% and
dropping a bin discards a twelfth of the year for nothing.

## The patch

`exoplasim/patches/exoplasim-3.4.2-lowio-first-record.patch`, **resident since
2026-08-18 and compiled into all five binaries.** It:

- saves and restores the six spectral accumulators with the root-only array
  routines the model already uses for `sz`, `sd`, `st`, `sp` and `so`, under new
  record names so a pre-patch restart is recognised rather than silently
  misread;
- adds the five accumulated orbital scalars and the eight accumulated hurricane
  indices, which were not in the restart at all and carry the same `delta`;
- writes a sentinel record, and on a restart that lacks it discards the partial
  interval and zeroes `naccuout` so the counter cannot again outlive the
  accumulators;
- recomputes `sqout` every timestep under low I/O, at the cost of `NLEV`
  transform pairs per step;
- (the code 54 `arasc`/`rasc` fix was SPLIT OUT of this patch on 2026-08-18 and now lives in `exoplasim/patches/exoplasim-3.4.2-arasc-output.patch`, which is resident alongside it; the parent carried a duplicate copy until then and the two could not both apply).

It is in `RESIDENT_PATCHES` and `rebuild_binaries.py --verify` unwinds it with
the rest of the stack.

### It is verified, as of 2026-08-18

Until this date the patch had never been run: it was authored, reviewed, applied
and rebuilt, and every run in `exoplasim/runs/` predated it. "It fixes bin 0" was
a design claim, not a measurement. It has now been measured.

The test uses this note's own probe. `sg` is surface geopotential and is CONSTANT
in time, so any bin-to-bin spread in it is pure artifact and the right answer is
zero. Two orbits per arm at `NLOWIO = 1`, because the defect needs a RESTART and
a first call has nothing to restore; two arms, one binary with the patch and one
with it reversed and recompiled, so the test can fail.

| arm | orbit | max bin-to-bin spread in `sg` | bottom-level `ua`, bin 0 |
| --- | --- | ---: | ---: |
| patch reversed | 0, cold start | 0.000e+00 | +0.018 m/s |
| patch reversed | 1, restarted | **9.209e-02** | **-34.68 m/s** |
| patch resident | 0, cold start | 0.000e+00 | +0.020 m/s |
| patch resident | 1, restarted | **8.651e-06** | **-0.58 m/s** |

The artifact falls by a factor of about 10,600, and the spurious solid-body
rotation predicted above -- which shows here as bin 0 carrying -34.68 m/s of zonal
wind against +0.13 m/s in the other eleven bins -- is gone.

Two things make this a test rather than a demonstration. Orbit 0 is EXACTLY
identical in both arms, which is what this note predicts for a first call and
what shows the probe is measuring the defect rather than inventing one. And the
reversed arm reproduces the defect, so a clean result from the patched arm cannot
be an insensitive metric.

The residual 8.651e-06 is not a defect. `sg` is accumulated over about 480
timesteps in single precision and then divided, so the rounding floor is between
2.6e-06 and 5.7e-05. The reversed arm's 9.2e-02 is four orders of magnitude above
any such floor.

### What the patch does NOT fix

It repairs the RESTART RESTORATION of the accumulators. It does not change WHAT
is accumulated, and the second defect on this page is about exactly that: under
`NLOWIO = 1` the model accumulates in a way that puts binned `spd` at 1.141 of
the snapshot mean-of-speed, between the speed of the time-mean vector and the
mean of instantaneous speeds. That number is a property of accumulating inside
the model and it survives this patch untouched.

That is why the two I/O regimes are still used for different things, and why the
guidance below is scoped by regime rather than withdrawn.

It changes the restart file layout: a restart written by a patched binary cannot
be read by an unpatched one, and the reverse, so a run started before 2026-08-18
cannot be resumed by a current binary. It buys back the 25x output volume that
`NLOWIO = 0` costs; it does not change any result, because no current result
comes from the low-I/O path.


## A third instance, 2026-08-19: `--restart-from` imports the donor's accumulators

Found while trying to use a low-I/O run as a control for CLIM-11. Under
`NLOWIO = 1` the land mask -- a field that CANNOT VARY -- reads **0.95341** in
orbit 0, bin 0, and exactly 1.0 in every other bin of every orbit. A constant
field is the cleanest possible probe for a normalisation error, and that is a
4.7% one.

**It needs BOTH seeding and low I/O**, which is what identifies the cause:

| run | seeded | NLOWIO | orbit 0 bin 0 |
| --- | --- | ---: | ---: |
| `run_b572b2e503b5` | no, cold start | 1 | **1.0** |
| `run_277e52971ec5` | `--restart-from` | 1 | **0.95341** |
| `run_ef3e195a7bdd` | `--restart-from` | 0 | 1.0 (nothing accumulates) |

`run_exoplasim.py --restart-from` hands the donor file to
`model.configure(restartfile=...)` whole, and PlaSim reads it whole -- including
`naccua`, `naccuout` and the accumulator arrays, which `seamod.f90:121-126`
restores by name. So a seeded run opens with another run's partial accumulation
already in the bucket and divides the sum by its own count.

**This is not the defect the patch above fixed.** That one was per model call and
is gone: orbits 1 and 2 are exactly 1.0 here. This is the first record of a RUN,
and only of a seeded one.

**Scope, which is small but not nothing.** One bin of one orbit. `align` and
`_bin_mean` in `close_ocean_energy.py` both drop bin 0 already, and a climatology
window is taken from the end of a run rather than its first orbit. What it did
cost is a measurement: the low-I/O arm of the CLIM-11 control pair is
contaminated, and a cold-started replacement is useless because a cold start
begins at 268.8 K with no ice-free ocean at all, so there is nothing to compare.

**It matters more now than it would have yesterday**, because low I/O became the
default for a prepared run the same day, and every seeded spin-up in
`WORKFLOW.md` section 6 -- the bracket points, the endmember arm -- is
`--restart-from`. CLIM-31 is the fix, and it is DONE as of the same day.
`exoplasim/scripts/reset_restart_accumulators.py` zeroes the accumulator records
in a copy of the restart before the seeded run reads it, and a seeded low-I/O run
now reads 1.0 in orbit 0 bin 0, matching the cold-start control.

**So all three instances of this class are now closed**, and they were three
different mechanisms wearing one symptom: `naccuout` surviving a restart while
the accumulators did not, fixed by the patch; the collective placed behind an
unbroadcast `nlowio`, fixed by the broadcast; and a seeded run inheriting a
donor's partial window, fixed by zeroing the copy. The symptom each time was a
first output record that did not mean what the records after it meant.
