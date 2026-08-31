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

`timestep_minutes` was one scalar in `config/planet.yaml` with no dependence on
`model.resolution`, so the first T42 arm ran at T21's 45 minutes. It integrated
46 orbits of ordinary climate and took a SIGFPE inside the 47th, on a gridpoint at -12.81
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

The ladder runs the ESCALATION ROUTE instead, decided in
`docs/src/pipeline/sequencing.md` section D and carried in machine-readable form
by `lib/rungs.py:ESCALATION_ROUTE`, which is where its per-rung steps are read
and not from here. The route's invariant is that resolution and timestep never
move together, so every conversion happens at constant dt. This paragraph
carried a copy of the route's three steps and had gone stale against it, which
is the whole reason a restatement is a pointer now.

**Three quantities, not three answers.** What each rung CAN take, what the route
RUNS it at, and what a commissioning-length run has SHOWN about a pair are
different facts, and each is now declared once in `lib/rungs.py`:
`STABILITY_CEILING_MINUTES`, `ESCALATION_ROUTE` and `COMMISSIONING_EVIDENCE`.
The ceiling is checked against the probe grid it is read out of and the route
against the chapter that decides it. The one-step-for-the-whole-ladder quantity
is deliberately not among them, for the reason above: it has no consumer. The residual measurement above still bears on it: a coarser step
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

The practical consequence for the ladder is a constraint on COMPARISON rather
than a step: two rungs compared at different steps measure resolution plus
truncation error and report the sum as resolution, so a comparison across rungs
is readable only where the step is held. Which step the route actually runs is
`lib/rungs.py:ESCALATION_ROUTE`.

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

## A converted arm relaxes fast and lands on the same climate, with more ice

*Measured 2026-08-30 at 2896f62ef. Both arms T42 at dt 45, `canonical-10m-carve2`,
the same staged surface family, the same binary, `NFILTEREXP` 16 and
`FILTERKAPPA` 8.0 in both. They differ in their initial state and in nothing
else. Two p8 arms, one pinned per die of this 7950X3D, under one host lock;
94 s an orbit at loads of 11 to 28.*

| arm | run | start | orbits |
| --- | --- | --- | ---: |
| A, cold | `run_88ed6f9d34ae` | cold | 144 |
| B, converted | `run_8d0ae7e2d02c` | `run_0d41aa82c287` at T21, converted at constant dt | 89 |

Both meet all six convergence criteria. So does the donor, at 108 orbits. Orbit
counts are `lib/run_lengths.py`'s commissioning bracket rather than round
numbers: the bottom of it was bought first, the comparison came back
indeterminate on one metric, and the declared response was to buy the top.

**The donor then moved the bracket it was bought under, which is the loop
working rather than a defect in the purchase.** `run_0d41aa82c287` is the
longest settled T21 production window this project has, and it reads a memory
time of 3.8616 orbits and an orbit scatter of 0.1235 K against declared upper
bounds of 3.69 and 0.117. Both are raised to cover it. Neither reading is about
a different world -- it is the same build and the same rung as the readings they
replaced, over a longer settled window, and this model's memory time grows with
the window it is measured on. The consequence for these arms is arithmetic and
was not acted on: the a-priori span at the new top is 77.4 orbits rather than
73.8, so a purchase made today would ask for 147 and 92 rather than 144 and 89.
Both arms converged at what they bought, and `production_span_from_report` is
the operative rule wherever a run exists.

**THE PAIR THIS REPLACES WAS NOT CONTROLLED, and that is why it was re-run
rather than re-measured.** `run_1d39fef9bfc2` and `run_42aaf441b10b` differed in
the physics filter exponent as well as in where they started -- `NFILTEREXP` 16
against 8, at the same `FILTERKAPPA`, `MPSTEP`, `NDEL`, `NHDIFF` and `TDISSZ` --
which is the knob `exoplasim/analysis/filter_gamma_pair.json` exists to price
and which moves the spectral bite fraction from 0.714 to 0.688 and the depth
from -0.609 to -1.391 dex between those same two runs. Their bounds were also
taken before WORLD-YJ9O, from the orbit-to-orbit scatter over the root of the
count. A better instrument on that pair would still have been comparing two
things that differ in two ways, so the arms were cut again with the filter held
at what `config/planet.yaml` declares. Reconstructed identity for both old arms
is under `archive/runs/`; WORLD-WW6Z carries how their records were lost.

### The comparison

`compare_equilibria.py` on the corrected instrument: sigma and the integrated
autocorrelation time over the longest stationary tail of each run, the bound
`2 sqrt(2) max(SEM_A, SEM_B)`, and a twenty-orbit window declared before the
arms ran.

| metric | A | B | B - A | bound | sigma |
| --- | ---: | ---: | ---: | ---: | ---: |
| surface temperature, K | 283.825 | 283.808 | -0.017 | 0.129 | 0.3 |
| air temperature 2 m, K | 283.554 | 283.499 | -0.055 | 0.127 | 0.9 |
| precipitation, mm/day | 2.0563 | 2.0532 | -0.0031 | 0.0104 | 0.6 |
| **sea ice fraction** | **0.0626** | **0.0767** | **+0.0141** | **0.0018** | **16.1** |
| TOA shortwave up, W/m2 | 93.270 | 93.250 | -0.020 | 0.582 | 0.1 |
| TOA longwave up, W/m2 | 228.504 | 228.507 | +0.002 | 0.405 | 0.0 |

Surface temperature agrees to 1.256 K RMS at a spatial correlation of 0.9972.

**THE TWO CLOUD TERMS ARE GONE, and they are refuted rather than merely
unconfirmed.** The old table read +2.492 and -2.503 W/m2 and called the arms
balanced at different cloud states. Controlled, they read -0.020 and +0.002
against bounds of 0.582 and 0.405, so this comparison would have seen a
partitioning a fifth the size of the one claimed and sees nothing. The 2.5 W/m2
belonged to the filter exponent, not to the conversion.

**The three thermal and hydrological rows agree, and one of them agrees at a
precision that cannot speak to the old number.** The smallest surface
temperature difference this comparison can resolve is 0.129 K and the old table's
offset was +0.106 K, which is below it. So that row is an unasked question here
rather than a refutation: what can be said is that the arms agree to within
0.017 K and that nothing of the old size is visible.

**WHAT SURVIVES IS A DIFFERENT FINDING, in a different variable.** The converted
arm carries 22.5 per cent more sea ice, 0.0767 against 0.0626, at sixteen sigma
on a bound of 0.0018. It is polar and in both hemispheres -- +0.076 in fraction
north of 60 degrees, +0.042 south of -60, +0.037 and +0.019 in the two
mid-latitude bands, and identically zero equatorward of 40 -- and 292 of 8192
cells differ by more than 0.05 in fraction. Both arms' ice series are stationary
over their tails, drifting by 3e-5 in fraction an orbit or less against a scatter
of 0.0017, so this is two settled ice margins and not one arm still moving.

**So the two arms are not at the same equilibrium, and the criterion is
all-or-nothing, but every row the previous verdict rested on has collapsed.**
The thermal state, the hydrological cycle and both radiative terms are
indistinguishable at this precision. What differs is the ice, which the
conversion remaps out of the donor and which has hysteresis at its margin, and
which shows no detectable radiative or thermal consequence at the precision that
sees it. That is the imprint the previous section could only guess at, now
located in a specific field.

Not established: whether the ice offset is the donor's imprint or a second ice
margin the target rung would reach on its own, which a third arm converted from
a DIFFERENT donor would separate; and whether it shrinks as the jump shrinks,
which matters because the ladder's next hop is T42 to T85.

### The relaxation is real and it is what the ladder is for

The converted arm meets every convergence criterion at 89 orbits against the
cold arm's 144, and its donor cost 108 at the rung below where a T42 cold start
costs about four times an orbit. The escalation route's promise -- that a
conversion buys the approach rather than a fresh spin-up -- holds on this pair.

### T42 endures dt 45

Neither arm failed. 144 orbits from cold and 89 from a converted state, on this
source, at the step `lib/rungs.py:ESCALATION_ROUTE` runs T42 at. That supersedes
`run_900548ae632e`'s blow-up in its forty-seventh orbit, which was taken on
source batch 2 has replaced and which `docs/src/pipeline/sequencing.md` section D
carries as reported rather than binding. `COMMISSIONING_EVIDENCE` is where the
row lives.


## What is not yet measured

Whether a wider window removes the flicker without moving the orbit at which a
run first genuinely settles. The sweep above is one run at one rung, and the T42
arms this note now carries are the dataset it should be re-taken on.

Whether the converted arm's ice offset is the donor's imprint or a second margin
the target rung reaches on its own; a third arm converted from a different donor
separates those. And whether it shrinks as the jump shrinks, which is what the
T42-to-T85 hop needs to know.
