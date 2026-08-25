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

## The empirical half, taken 2026-08-25

*Measured at 5dcdf64b. The four-byte arm is `build_model.py --res T21 --ranks 16
--precision 4 --no-publish`, run through `run_exoplasim.py --binary`, which is
the route world-u5pf opened. Both halves of the experiment below were run on a
host shared with other work; which numbers that touches is stated where it
matters.*

### The declared 40-orbit experiment gave an answer that was not stable

The experiment fixed above -- two arms, 40 orbits each, the last 20 as the
comparison window -- was run. Over orbits 20-39 the four-byte arm was 0.077 K
COOLER, and `compare_equilibria.py` failed it at 2.04 sigma on surface
temperature with the other five metrics inside.

That result does not survive a longer run. Repeated at 85 orbits, the four-byte
arm is 0.187 K WARMER over orbits 65-84, at a nominal 7.24 sigma and failing four
metrics of six. **The sign reversed.**

The 40-orbit reading carried its own warning and the warning was right. Forty
orbits against a 0.10 K orbit-to-orbit scatter give a slope standard error of
1.40 mK/orbit, which over the 118.2 orbits a 691,200-step spin-up takes at dt 45
is plus or minus 0.33 K -- wider than the 0.026 to 0.156 K the experiment existed
to test. The lever arm could not carry the claim, that was written down before
the 85-orbit arm was run, and the failure it predicted is the one that arrived.

### What is actually there, and why the sigma figures are not it

Both 85-orbit runs reproduce their 40-orbit predecessors BIT FOR BIT over orbits
0-39, in both arms, maximum absolute difference 0.0 K. So everything separating
the two arms is the compiled precision and nothing is a different draw.

The difference series carries a lag-1 autocorrelation of 0.615. Eighty-five
orbits therefore hold about twenty independent samples, and the whole-run mean
difference in global-mean surface temperature is

    +0.0139 +- 0.0352 K

which is four tenths of a standard error. The window figures that look decisive
are smaller than the eight-byte arm's OWN variability: its four non-overlapping
20-orbit means are 291.355, 291.835, 291.757 and 291.651 K, a spread of 0.210 K
and a range of 0.480 K, against window differences of -0.068, -0.105, +0.062 and
+0.187 K. `compare_equilibria.py` takes its standard error on the raw orbit count
and so understates it by about a factor of two; that is world-yj9o, and it is not
specific to this comparison.

**The deadband failure mode is not detected.** Its signature is growth with step
count, and over 497,250 steps -- 72 per cent of a 691,200-step spin-up -- no
systematic offset survives the model's own low-frequency variability at the 0.035
K the run can resolve. That is not the same as a bound at the 0.009 to 0.050 K
the criterion corresponds to, and reaching one needs a window long enough to hold
seventeen independent samples, which at this autocorrelation is about seventy
orbits rather than twenty. `exoplasim/notes/convergence-lengths.md` carries that
arithmetic and the convergence lengths the experiment should have been designed
against.

### The cost, and it is large at the rung that matters

`bench_ab.py`, paired and interleaved with the arm order flipped every round, on
cold beds. Cold because a warm bed cannot be shared: see the next section.

| rung | bed steps | fp64 median | fp32 median | gain | range | rounds fp32 faster |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| T21 | 600 | 1.667 s | 1.162 s | +27.5% | +18.8 to +33.4 | 6 of 6 |
| T21 | 6000 | 13.151 s | 9.491 s | +29.4% | +26.0 to +40.0 | 6 of 6 |
| T85 | 600 | 26.777 s | 15.765 s | +39.2% | +29.6 to +46.0 | 6 of 6 |
| T85 | 1500 | 69.245 s | 39.933 s | +40.5% | +15.4 to +51.3 | 6 of 6 |

**Three of those four sets are REFUSED by `bench_ab.py`'s own guard**, and the
fourth passes only as readable-and-noisy: round-to-round self-scatter runs 24.5
to 275.7 per cent against the 5 per cent floor this project declares as "no
difference". The host carried other work throughout and no single set clears the
declared criterion. That is stated rather than rounded into a pass.

What carries the claim instead is two things the scatter cannot reach. The
direction is 24 of 24 rounds, and interleaving flips the arm order every round,
so contention cannot systematically favour one arm. And the T85 median reproduces
`exoplasim/notes/shtns-viability.md`'s independent 39.69 per cent [+39.63,
+39.82], measured at T170 on sixteen threads from a different bed on a different
day, to within a point.

The gain grows with the rung because four bytes buys memory bandwidth, and it is
largest where a spin-up is actually expensive.

### The four-byte arm cannot read an eight-byte restart

Both warm beds aborted inside the model at `plasim.f90:67` on the four-byte
binary, because `make_profile_bed.py` had copied in the eight-byte run's
`plasim_restart`. This is not a Python check refusing a file: it is the model
reading four-byte records out of an eight-byte file. `restart_surface.py:302`
carries the same assumption one level up, which is why a four-byte run cannot be
continued at all (world-73sn) and why the runs above are single
`run_exoplasim.py --run-years N` calls.

So the technique's own handoff -- four bytes spins up, eight bytes does the final
approach and the judging -- **requires CLIM-52's converter in both directions and
is not an optional convenience.** Nothing in the tree converts a state between
precisions today.

### Where that leaves the decision

Measured: the cost gain is real and about 40 per cent at T85, and no climate
difference survives the model's own variability over 85 T21 orbits at the 0.035 K
those orbits can resolve. Not measured: a bound at the criterion's own 0.009 to
0.050 K, which needs a window of about seventy orbits rather than twenty; and
anything at all about a four-byte arm at T85 or above, where the state is larger
and the transform tables demote further.

Not physics but plumbing, and it is the whole of what stands in the way: a
precision converter both ways, a restart reader that takes its record width from
the executable that wrote the file, and a route by which a four-byte run can be
continued rather than taken in one call.
