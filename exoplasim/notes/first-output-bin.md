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

## What to do about it: refuse the product, do not correct it

**Run with `NLOWIO = 0` and read the binned climatology normally.** Both defects
are artifacts of the low-I/O path and they vanish together. Measured on this run:

| | first-bin wind ratio | binned `spd` / snapshot mean-of-speed |
| --- | ---: | ---: |
| orbit 65, `NLOWIO = 1` | 9.00 | 1.141 |
| orbit 70, `NLOWIO = 0` | 1.03 | 1.002 |
| orbit 76, `NLOWIO = 0` | 0.99 | 0.999 |

That second column is the important one and it corrects an earlier conclusion in
this note. The binned `spd` under `NLOWIO = 1` sits between the speed of the
time-mean vector and the mean of instantaneous speeds, and that was written up
here as a property of the binning -- as something that "outlives" the switch,
because a bulk flux wants a mean of the speed and a bin mean of a vector cannot
be one. **It does not outlive it.** The cancellation happens in the model's
output accumulation, not in the binning, and with the accumulation off the binned
field is the mean of the speed to 0.2%.

So the corrections that were written against it are gone rather than improved,
and that is deliberate:

- `carve_verdict.py:turbulent_forcing` read wind from the snapshot product and
  dropped the corrupt record from humidity. It now reads both from the binned
  climatology, which is strictly better at 182 samples an orbit against 32.
- `build_dust.py:speed_bias_correction` scaled the binned wind up by a per-cell
  factor with a median of 1.554. On clean output that would be a 55% error.

**Both corrections are wrong once the defect is gone, so the guard refuses
instead.** `build_climatology.py` stamps `low_io` on every product it writes,
and `lib/paths.py:require_clean_io` raises on a tainted one. A product with no
stamp counts as tainted, because every climatology built before 2026-08-17 was.

**Do not drop the first bin from anything else.** Every scalar is within 2% and
dropping a bin discards a twelfth of the year for nothing.

