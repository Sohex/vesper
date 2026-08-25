# The resolution ladder: what a rung costs, and what the convergence test costs on top

*Worldbuilding frame: measurements on the Vesper project's climate model, an
ExoPlaSim fork, run on this desktop. Nothing here is about the simulated planet.
Measured 2026-08-24 at b30e314e, sixteen threads, `-O2` production profile.*

The goal this serves is stated as a wall clock: reach a converged state at the
highest resolution in the least real time. That makes two things worth
measuring separately -- how long a rung takes to converge, and how many of
those orbits the convergence TEST is responsible for rather than the physics.

## The T21 baseline, re-run

`run_8044646ea7f0`, T21, 85 orbits from cold in one segment, 13.2 s per orbit
under light contention.

It replaces `run_2b20e3324bb0`, which is not a baseline: under world-1nz
`continue_exoplasim` never wrote the hyperdiffusion namelist, so 84 of that
run's 85 orbits ran the model's compiled defaults -- `ndel` 2 rather than 4,
grad^4 rather than grad^8, humidity damped about 7.4 times too hard. The two
runs differ by **3.0 K** in equilibrium temperature, 294.46 K against 297.48 K.
That is the size of the defect, and it is why every number the old run produced
goes with it.

## The criterion decides convergence later than the physics does

Assessed at every truncation from 20 to 85 orbits with `--through`:

| orbits | verdict |
| --- | --- |
| 20-57 | fail |
| 58-64 | pass |
| 65-66 | fail |
| 67-71 | pass |
| 72-74 | fail |
| 75 | pass |
| 76-77 | fail |
| 78-85 | pass |

The temperature is stationary from orbit 58: 297.42 K with a spread of 0.07 K
over 58-85, and no trend inside it. The verdict is not. It flips seven times
after the run has settled, and the criterion that flips is always the same one,
`extrapolated_offset_lt_0.15_k`, which tests `|offset| + half_width`. That
statistic wanders between 0.031 and 0.255 on a state that is not moving.

So **the run converges at orbit 58 and the test agrees at orbit 78**. Twenty
orbits separate the two, and they are not physics. At T21 they cost four
minutes. The same twenty orbits at T170, where one orbit is about 128 times a
T21 orbit, cost about eight hours.

The mechanism was already known and is now priced: at a ten-orbit window the
standard error of the fitted offset is 0.113 K against an allowance of 0.075 K,
so the test is noise-dominated by half again. What is new is that the noise has
a wall-clock price, and that the price is paid at the most expensive rung.

## A twenty-orbit window declares convergence thirteen orbits sooner

The obvious objection to widening the window is that a wider one needs more
orbits before it can return any verdict at all, so the two effects trade. They
do, and the trade has an interior optimum. Swept over `run_8044646ea7f0`,
window against every truncation from 35 to 85 orbits:

| window | first pass | stable from | verdict flips |
| ---: | ---: | ---: | ---: |
| 10 | 58 | 78 | 7 |
| 15 | 62 | 78 | 3 |
| **20** | **65** | **65** | 1 |
| 25 | 68 | 68 | 1 |
| 30 | 71 | 71 | 1 |
| 35 | 76 | 76 | 1 |

"Stable from N" is defined as: every truncation from N to 85 passes. It was
fixed before the sweep was run, because a criterion chosen after the run it
judges is not a criterion.

From window 20 upward the test stops flickering entirely -- first pass and
stable point coincide -- and from there each extra five orbits of window costs
about three orbits before the verdict arrives. Below 20 the noise costs more
than the width saves: window 10 first passes at 58 and then takes another
twenty orbits to stop changing its mind.

So **window 20 declares convergence at orbit 65 where window 10 declares it at
78**. Widening the window does not make the test slower to satisfy; at the
current width the test is paying thirteen orbits for its own scatter. At T170,
where the note's own measurement puts an orbit at 1530 s at dt 22.5, thirteen
orbits is five and a half hours.

This is one run. It is not yet the default, and what would settle it is the
same sweep on the T42 arms -- a window chosen on a single rung is a window
fitted to one dataset.

## The T42 rung does not run at the T21 timestep

`timestep_minutes` is one scalar in `config/planet.yaml` and does not move with
the rung, so the first T42 arm ran at T21's 45 minutes. It integrated 46 orbits
of ordinary climate and took a SIGFPE inside the 47th, on a gridpoint at -12.81
K at the second level from the top. Nothing was building towards it: level 2's
coldest cell held 190 K with no trend for the whole run. That is world-td3, and
`exoplasim/notes/physics-filter-stability.md` carries the evidence.

What it costs the ladder is the reason it belongs here. That note's stability
grid qualifies a step on 400-step probes and one full orbit, and it already
warns that a probe cannot see a late blow-up. **An orbit cannot see one
either.** A step qualified for a two-orbit diagnostic is not thereby qualified
for a commissioning run of eighty-five, and the only instrument that settles it
is a run of the length actually intended.

**That argument is superseded and the reason it fell is worth keeping.** It ran:
the ladder therefore runs T42 at dt 22.5, which that note measures as clean at
kappa 8 and which is also the step T170 needs, so the rungs are compared at one
step rather than each at its own margin -- which matters, because the same note
measures the spectral core's energy residual falling from 0.31 to 0.12 W/m2 when
the step is halved. Every clause of that is still true. What went is its
CONSUMER: comparing rungs at one step serves a between-rung convergence
comparison, and SPAT-8, which was that comparison, closed when T85 was declared
the operating support rather than a rung to be discovered. Nothing now needs the
rungs to share a step, so nothing pays for running T42 at 22.5.

The ladder runs the ESCALATION ROUTE instead -- T21 at 45, T42 at 30, T85 at
22.5, in `docs/src/pipeline/sequencing.md` section D -- where the step changes
so that resolution and timestep never move together and every conversion happens
at constant dt. The residual measurement above still bears on it: a coarser step
carries a larger energy residual, so T21 at 45 carries the largest of the three,
and that is a cost of the route to be reported rather than a reason to change
it. WORLD-F997. A ladder whose rungs each ran at their own
largest stable step would be comparing equilibria that differ by their
truncation error as well as by their resolution.

## The timestep moves the equilibrium further than the resolution does

Three converged runs, 85 orbits each from cold, same build and same staged
surface, differing one factor at a time:

| run | rung | dt, min | mean T, K |
| --- | --- | ---: | ---: |
| `run_8044646ea7f0` | T21 | 45.0 | 297.458 |
| `run_aaa95662e21a` | T21 | 22.5 | 295.867 |
| `run_4fe5e6050df5` | T42 | 22.5 | 294.567 |

| factor | change |
| --- | ---: |
| timestep, 45 to 22.5 at T21 | **-1.591 K** |
| resolution, T21 to T42 at dt 22.5 | **-1.300 K** |
| both | -2.891 K |

**Halving the step moves the equilibrium by more than doubling the resolution
does.** That is not a small correction to a resolution study; it is the larger
of the two effects, and a ladder whose rungs each ran at their own largest
stable step would have reported it as resolution sensitivity.

It is also consistent in sign and rough size with what
`exoplasim/notes/physics-filter-stability.md` measures independently: the
spectral core's adiabatic non-conservation falls from 0.31 to 0.12 W/m2 when
the step is halved at T42. A step that conserves worse warms the model, and
halving it removes about 0.19 W/m2 of spurious heating -- which against this
model's sensitivity is the right order for 1.6 K.

The practical consequence for the ladder is the one already taken for stability
reasons: **every rung runs at one step**, dt 22.5, which is what T170 needs.
Otherwise the ladder measures resolution plus truncation error and reports the
sum as resolution.

## A conversion may double the truncation and no more

*Measured 2026-08-24 at 220b7142. 900-step arms, each in a bed built from the
target rung's own run directory, each compared against a control that differs
only in which restart is copied in.*

| conversion | jump | target step | result |
| --- | --- | --- | --- |
| T21 -> T42 | 2x | dt 22.5 | integrates, and relaxed for 16 orbits |
| T42 -> T85 | 2x | dt 22.5 | completes |
| T85 -> T170 | 2x | dt 15 | completes |
| **T42 -> T170** | **4x** | **dt 15** | **traps** |

The controls are what make this readable. At each target, that rung's OWN
restart was run through the same bed at the same step: T85's completes at dt
22.5 and T170's completes at dt 15. So the bed, the namelist, the binary and
the step are all cleared, and the only thing left is the conversion.

The converter treats the two jumps identically -- the action counts are the
same to the record, 40 remapped, 10 projected, 114 reset, 27 from the template
-- so this is a property of the DATA, not of the policy. Doubling the
truncation produces a state the target integrates; quadrupling it does not.

**So the ladder must step, and the step is a factor of two.** T21 to T170
directly is a factor of eight and is not available; neither is T42 to T170. The
rungs that exist -- 21, 42, 85, 127, 170 -- give 21 -> 42 -> 85 -> 170 as three
doublings, which is the path.

What this does NOT establish is a mechanism, or where between 2x and 4x the
boundary lies. T42 -> T127 is 3x and untested. And 900 steps is a floor: a
converted state that starts has not been shown to survive commissioning, which
is the same caveat the stability grid carries.

## The refusal at T170 is not a startup transient

The cheapest hope for the ladder was that a rung refuses a coarse step only
because a cold start is the most violent thing it ever integrates, and that a
spun-up state handed to it would sail through. It does not.

T170 at 900 steps:

| initial state | dt 15 | dt 30 |
| --- | --- | --- |
| cold | runs | refuses |
| T170's own restart | completes | traps |
| converted from T85 | completes | traps |

dt 30 fails from every initial condition available, including the rung's own
balanced restart, and the dt 15 column is the control showing those same states
are otherwise fine. The refusal is a property of the configuration rather than
of the transient, so buying a coarser step by arriving at it gently is not
available and the ladder pays T170's dt 15 in full.

## A converted arm relaxes four times faster and lands somewhere else

*Measured 2026-08-25 on the corrected damping. Both arms T42 at dt 30, same
build, same staged surface, differing only in where they start.*

| arm | start | converged at |
| --- | --- | ---: |
| A, `run_1d39fef9bfc2` | cold | 85 orbits asked, converged |
| B, `run_42aaf441b10b` | converted from converged T21 | **19 orbits**, held at 22 |

Arm B is inside a tenth of a kelvin of arm A's equilibrium after TEN orbits and
meets every convergence criterion at nineteen. That is the ladder's whole
promise and it is real: the relaxation from a converted state is a small
multiple of the convergence window, not a fresh spin-up.

**And the two arms are not at the same equilibrium.** `compare_equilibria.py`,
whose criterion was fixed before either arm ran, against the runs' own
ten-orbit scatter:

| metric | A | B | B - A | bound | sigma |
| --- | ---: | ---: | ---: | ---: | ---: |
| surface temperature, K | 294.549 | 294.655 | +0.106 | 0.042 | 5.1 |
| air temperature 2 m, K | 294.078 | 294.166 | +0.088 | 0.041 | 4.3 |
| precipitation, mm/day | 2.796 | 2.742 | -0.054 | 0.016 | 7.0 |
| TOA shortwave up, W/m2 | 76.175 | 78.667 | **+2.492** | 0.158 | 31.5 |
| TOA longwave up, W/m2 | 245.596 | 243.093 | **-2.503** | 0.191 | 26.2 |

It does not close with time. Checked at 22, 35, 45 and 47 orbits the surface
offset sits at +0.100, +0.111, +0.096 and +0.106 K, and the two TOA terms hold
near +2.5 and -2.5 throughout. A run still relaxing narrows; this does not.

**The two shortwave and longwave terms cancel.** Net TOA differs by about 0.01
W/m2, so both arms are in energy balance -- they are balanced at DIFFERENT
CLOUD STATES. Arm B reflects 2.5 W/m2 more and emits 2.5 W/m2 less, which is
more cloud, and it is warmer underneath by a tenth of a kelvin and drier by
five hundredths of a mm/day. The spatial correlation is 0.990 at 1.81 K RMS, so
this is one climate with a systematic offset rather than two different worlds.

What that means for the ladder is a judgement rather than a measurement. The
offset is 0.03% of the surface temperature and is detectable only because the
criterion is strict -- the bound is the runs' own scatter, about 0.04 K. Whether
a tenth of a kelvin and a cloud partitioning of 2.5 W/m2 is acceptable for a
rung's output depends on what that rung's output is for, and the honest
statement is that a converted arm is NOT a substitute for a cold arm at the
precision this criterion can see.

Not established: the mechanism. The donor's cloud and humidity structure is
remapped into the target and the model may simply keep it, which would make
this an imprint rather than a second equilibrium; a third arm converted from a
DIFFERENT donor would separate those. Nor is it known whether the offset shrinks
as the jump shrinks, which matters because the ladder's later hops are T42 to
T85 and T85 to T170 rather than T21 to T42.

## What is not yet measured

Whether a wider window removes the flicker without moving the orbit at which a
run first genuinely settles, and what the T42 rung converges to -- both from
cold and from a converted T21 state. Those are the arms this note is waiting
on.
