# The first output bin of every orbit has no boundary layer

Measured 2026-08-17. It is genuinely bad output rather than a misreading, and
the mechanism is specific enough to say what it is and what it is not.

## What it looks like

In every `MOST.NNNNN.nc` this project has produced, the first of the twelve time
bins carries wind and humidity fields that do not belong with the other eleven.
Checked on `run_b014469b8091` orbits 86 to 90 and on `run_bfa3f5269660`: the
bottom-level wind speed ratio between bin 0 and the rest is 7.50, 7.57, 7.54,
7.58 and 7.97. Systematic, not a restart shock.

**Only the wind and humidity fields are affected.** Every scalar is clean to
better than 2%: `ts`, `tas`, `pr`, `evap`, `mrso` and `ps` all sit within 0.4 to
1.8% of the other bins. So results that do not read a wind are unaffected, and
that covers the flux bracket, the soil, surface code 229 and the convergence
assessment.

## What it is

Not a short average, not an instantaneous state, and not an averaging mistake of
ours: it is in the raw model output before anything of this project's touches it.

The vertical profile of bin 0 against bins 1 to 11 is what identifies it. Mean
absolute eastward wind by model level, top to bottom:

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

**The free troposphere agrees to within 5% and the boundary layer is missing.**
Bin 0's wind is nearly constant with height, 19 m/s at the top and 34 m/s at the
surface. The other eleven bins decay downward from 18 to 3 m/s, which is what
surface friction does. Humidity is 27% low in bin 0 and temperature 1.6% low,
which is the same story: the surface fluxes have not mixed moisture upward, and
temperature barely moves because a fractional change on 281 K is small.

So bin 0 is a state in which the surface drag and boundary-layer mixing have not
been applied. Each orbit is a separate model call resuming from a restart, which
is why it recurs every orbit rather than only at the start of a run.

## Why it cannot be an average of physical states

A time-mean cannot exceed the extremes of what it averages. Bin 0's mean
absolute eastward wind is 34.4 m/s against 4.94 for the 32 instantaneous
snapshots of the same orbit, a factor of 7, while its maximum exceeds the
snapshot maximum by only 1.54. A mean seven times larger than the mean of the
instantaneous states it supposedly averages is not an average of them.

## Where it comes from: the model's low-I/O path, not the postprocessor

Tested 2026-08-17 by running one orbit from the bootstrap's final restart with
`NLOWIO = 0` instead of the default 1, in an isolated copy, and putting the
result through the same `pyburn` averaging the production files get.

| | bin 0 over bins 1-11, bottom-level specific humidity |
| --- | ---: |
| `NLOWIO = 1`, production | 0.714 |
| `NLOWIO = 0`, this test | 0.930 |

0.930 is inside the ordinary seasonal spread: the twelve bins run 0.00706 to
0.00797 and bin 0 is simply the lowest of them. The defect is gone.

**Only the model's output path changed.** `pyburn` did the same 12-bin averaging
in both cases -- it reported "going from 182 timestamps to 12" for the test -- so
the averaging is not at fault and neither is the vector transform that derives
`ua`, `va` and `spd` from divergence and vorticity.

`NLOWIO = 1` is PlaSim's default, set in `plasimmod.f90:151`. It is the path that
accumulates fields over the output interval and divides in place at write time,
`outmod.f90:257` onward, using `naccuout`. That is where the bug is.

**What is not yet pinned is the line.** A wrong divisor would scale a field
uniformly, and this does not: the ratio runs 1.04 at the model top and 11.2 at
the surface, so bin 0 has a different vertical structure rather than a scaled
one. Whatever is wrong involves the vertical weighting of the accumulation, not
just its count. Note also that each orbit is a separate model call and the first
output interval of a call is short -- the records are stamped 319 then every 480 --
so a partial first interval is the most likely trigger.

## What to do about it

**Wind is read from the snapshot product, humidity from the binned one with the
first record dropped.** `hydrography/scripts/carve_verdict.py:turbulent_forcing`
is the single place that does it, and the other two Penman consumers call it.

The two fields need different treatment because only one of them has a second
problem. Measured on the baseline climatology, bottom level, area-weighted:

| | m/s |
| --- | ---: |
| binned `spd`, all 12 bins | 8.411 |
| binned `spd`, bins 1-11 | 5.033 |
| `sqrt(ua^2+va^2)` from the binned components, bins 1-11 | 3.922 |
| snapshot `spd`, 32 instantaneous samples an orbit | 7.356 |
| `sqrt(ua^2+va^2)` from the snapshots, per sample then averaged | 7.356 |

**In the snapshots `spd` is exactly `sqrt(ua^2+va^2)`.** Ratio 1.000, correlation
1.0000, at every model level and on all five orbits of the climatology. So code
259 is the instantaneous wind speed and there is nothing to resolve about what it
means. An earlier reading of this table had the two disagreeing by 1.5x in the
snapshots; they do not.

**The binned `spd` is the one that is not what it looks like.** It sits between
the speed of the time-mean vector, 3.922, and the mean of instantaneous speeds,
7.356 -- so the model's output accumulation cancels the reversing component over
each output interval, and pyburn then averages twelve of those. That is a
partially vector-averaged speed, and a bulk aerodynamic formula wants the mean of
the speed. The snapshot product gives it, at 32 samples an orbit: coarse, but
unbiased, which the binned product is not.

So for wind, dropping bin 0 was never the fix. It takes 8.411 to 5.033 and the
answer is 7.356.

**For humidity there is no snapshot to switch to** -- `hus` is absent from the
snapshot output despite code 133 being requested -- and none is needed, because
the binned mean is the right kind of average and only the one record is bad.
Dropping it costs a twelfth of the seasonal cycle. That is worth 0.08%:
reconstructing the record by cyclic interpolation from its neighbours instead
gives an annual mean 0.08% from simply dropping it, against the 2.1% the
uncorrected mean is out by.

**Do not drop it from anything else.** Every scalar is within 2% and dropping a
bin there discards a twelfth of the year for nothing.

The repair is gated on detecting the defect, by the first bin's global-mean
bottom-level wind exceeding three times the median of the others. The observed
factor is 7 to 9 and no seasonal cycle approaches 3. A climatology from an
`NLOWIO = 0` run passes through untouched.

## What this cost before it was found

The uncorrected Penman wind of 8.411 m/s was 14% high, not the 3% a comparison
against a different run had suggested, and it fed the carve verdict, the lake
solution and the soil. The near-cancellation of two errors is what made it look
harmless: the corrupt record inflates and the vector averaging deflates, and how
completely they cancel depends on the run. Do not carry a cancellation between
runs.
