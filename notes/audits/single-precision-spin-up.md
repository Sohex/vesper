# What a single-precision spin-up would change, and the one bound the source cannot tighten

*Worldbuilding frame: a NUMERICS audit of the Vesper project's climate model.
Nothing here is about the simulated planet. Derived 2026-08-25 at 4e14c4a9 from
the model source and `config/planet.yaml`; the two compiled checks are stated
where they are used. No timing was measured and none is quoted.*

`compile.sh -p 4` builds the model at four bytes per real, which is upstream's
default and what this project overrides to eight. The question clim-59 asks is
whether that is sound for a SPIN-UP, given that single precision fails exactly
where accumulation is long, and a spin-up is the longest accumulation the
project runs. The answer is available from the source and one bound is not.

## What four bytes actually reaches

`plasim/CMakeLists.txt` appends `-fdefault-real-8` only when `PLASIM_PRECISION`
is 8. At 4 the flag is simply absent, so DEFAULT-kind reals demote and
explicitly kinded ones do not. That splits the model in two:

| stays eight bytes | demotes to four |
| --- | --- |
| `gaussmod`'s PI | the whole prognostic state in `plasimmod` |
| `legmod`'s Gaussian-latitude and Legendre SETUP scratch (`amsq`, `z1`..`z3`, `zpli`, `zpld`) | `legmod`'s transform TABLES: `pmat`, `qmat`, `fsp`, `fgp`, `fmu`, `fmv`, `fmm`, `fmq`, and the filter vectors `skgpsp`, `skspgp` |
| `shtnsmod`'s transform buffers (`shtinv`, `shtl1`, `zlm`, `zg`) and the SHTns library itself | every gridpoint field, every surface reservoir, every restart record |

COMPILED, not read: `build_model.py --res T21 --ranks 16 --precision 4
--no-publish` builds and links, and the generated `build.ninja` carries no
`-fdefault-real-8`. The mixed explicitly-kinded and default-real code compiles
at four bytes without a kind mismatch.

## Where single precision fails, in this model's arithmetic

The failure mode is STAGNATION: an additive update whose increment is smaller
than half an ulp of the value it is added to rounds back to that value, and the
increment is lost entirely. Every reservoir in this model is updated that way
and holds an absolute temperature in Kelvin, so the ulp is fixed by the
temperature and the stagnation threshold falls out of the heat capacity.

`oceanmod` is the binding case: `zsst = zsst + F * dtmix/(CRHOS*CPS*ymld)`. At
the declared mixed layer depth and the declared sea water constants the slab
holds 2.055e8 J/(m2 K), one timestep of 2700 s moves it 1.314e-5 K per W/m2,
and the float32 ulp at 290 K is 3.052e-5 K. So:

| reservoir | net flux below which one step changes nothing |
| --- | --- |
| soil layer 1 (0.4 m) | 0.005 W/m2 |
| soil layer 5 (6.4 m) | 0.087 W/m2 |
| **ocean slab, 50 m** | **1.16 W/m2** |
| an ocean layer at 200 m | 4.65 W/m2 |
| an ocean layer at 500 m | 11.6 W/m2 |
| an ocean layer at 2000 m | 46.5 W/m2 |

The slab's threshold is two orders of magnitude below the instantaneous net
surface flux at an ocean gridpoint, which the diurnal cycle carries to O(100)
W/m2 either side of zero. So no increment that matters is lost, and what
single precision leaves behind is unbiased rounding noise of half an ulp per
step. That noise is confined by the slab's own restoring: at feedbacks from 1 to
5.53 W/(m2 K) the restoring timescale is 13,800 to 76,100 steps and the
stationary spread is 7.3e-4 to 1.7e-3 K.

**The instrument against the size of the effect.** The convergence criterion is
0.05 W/m2 per orbit, which at those same feedbacks is 0.009 to 0.050 K of
remaining offset. The single-precision noise floor is below it by at least a
factor of five in the least favourable pairing. Rounding noise is not what
would make a single-precision spin-up unsound.

## The verdict is conditional on there being no deep reservoir

`NLEV_OCE` is 1, `nlsg` defaults to 0 and `nco2evolve` to 0, so the
longest-memory state in this model is the 50 m slab and there is nothing behind
it. The table above is why that matters: the stagnation threshold scales with
layer depth, and no gridpoint's net surface flux reaches 11.6 W/m2 on a time
mean, let alone 46.5. **Adopting any deep ocean makes single precision unsound
at once**, and that is not a remote possibility here -- the ocean tier is an
open design question and cGENIE is vendored as a candidate under OCN-3. The
condition belongs with the verdict, not under it.

## The transform: the exposure is in the reference arm, not the production one

`shtnsmod` declares its buffers `real(kind=8)` and `complex(kind=8)` and
converts the state across explicitly, and SHTns is built double. So at four
bytes the production transform still COMPUTES in double and what changes is a
four-byte round trip of the state at each call, 6e-8 relative.

`legmod` is the opposite: its tables and its accumulation are default reals and
demote with everything else. The accumulations run over O(NLAT) and O(NTRU), so
sqrt(N) eps puts the per-transform relative error near 1e-6, roughly twenty
times the round trip alone. The consequence is not that legmod is unusable but
that a four-byte build cannot be checked against its own legmod reference on
equal terms -- and legmod is the arm `verify_shtns_model.sh` compares against
and the arm the `poisoned` profile has to be run with.

## What four bytes breaks that is not physics

The `time4*` counters in `miscmod`, `rainmod`, `fluxmod` and `carbonmod`
accumulate elapsed seconds with sub-millisecond increments. Once such a counter
reaches 1e4 accumulated seconds its float32 ulp is 9.8e-4 s and an increment
below half of that is lost outright, so the counter stops advancing while the
run continues. Nothing in this project reads them, but a cost figure taken off a
four-byte build is not a cost figure.

`glaciermod`'s `persistt` is the one clock in the model STATE: seconds of
continuous snow cover per gridpoint, counted up to `glacpersec`. At the declared
one-orbit persistence the target is 1.97e7 s, where the float32 ulp is 2 s
against a timestep-sized increment, so it advances. Its headroom shrinks as the
persistence target grows, and that is worth knowing before the target is raised.

## The one bound the source cannot tighten: the deadband

The 1.16 W/m2 threshold is not only a stagnation floor, it is a DEADBAND. A step
whose net surface flux falls inside it moves the slab not at all, and the
fraction of steps that do is set by how much time the flux spends near zero. For
a diurnal cycle of amplitude A the fraction is (2/pi) asin(1.16/A), and every
dropped increment is at most half an ulp:

| diurnal amplitude | steps inside the deadband | one-sided worst case over a 691,200-step spin-up |
| --- | --- | --- |
| 50 W/m2 | 1.48% | 0.156 K |
| 100 W/m2 | 0.74% | 0.078 K |
| 300 W/m2 | 0.25% | 0.026 K |

The right-hand column assumes every dropped increment has the same sign and is
maximal, which it will not be: the drops happen at the zero crossings, where the
flux is changing sign, so they very largely cancel. But that is an argument, not
a bound. The strict worst case is not below the 0.009 to 0.050 K the criterion
corresponds to, and it cannot be tightened from the source, because tightening
it needs the distribution of the net surface flux near zero and that is a model
output. **This is the whole of what is not established.**

## The comparison has to be climatological, and cannot be state-by-state

clim-59 records the FP32 state's landing place as unknown "because four-byte
restart records against eight are different files to `compare_restarts.py`
rather than different answers". Converting the file is necessary and not
sufficient. The state handed across each transform is truncated to seven digits,
the flow is chaotic, and two trajectories from one initial state decorrelate
within model days whatever the converter does. A four-byte and an eight-byte run
are not the same trajectory and no field-by-field restart comparison of them
means anything. What can be compared is the CLIMATE: time means over many
orbits, judged against the eight-byte run's own orbit-to-orbit spread.

## The experiment that settles it, with its criterion fixed in advance

Two arms at T21 on sixteen threads, identical configuration and
`cold_start_seed`, one at the declared precision and one at four bytes, 40
orbits each with the last 20 as the comparison window. FIXED BEFORE ANY RESULT
IS SEEN: the four-byte climate is indistinguishable if the difference in
global-mean surface temperature over the window is inside two standard errors of
the eight-byte window's own orbit-to-orbit spread, and likewise for global-mean
TOA balance. A miss is a bias; either way the bound is what gets reported, and
it is the bound the deadband row above cannot supply.

The four-byte arm can now be built:

    python exoplasim/scripts/build_model.py --res T21 --ranks 16 --precision 4 --no-publish --print-path

It cannot be RUN. `run_exoplasim.py` and `continue_exoplasim.py` compose the
registry's executable name and take the copy in the run directory; neither has a
flag naming a binary, and an arm at a precision other than the declared one is
refused publication because the registry's name carries no precision. That is
the blocker on this half, and it is a missing route rather than a busy machine.
