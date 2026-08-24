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

The ladder therefore runs T42 at dt 22.5, which that note measures as clean at
kappa 8 and which is also the step T170 needs. The rungs are then compared at
one step rather than each at its own margin -- which matters, because the same
note measures the spectral core's energy residual falling from 0.31 to 0.12
W/m2 when the step is halved. A ladder whose rungs each ran at their own
largest stable step would be comparing equilibria that differ by their
truncation error as well as by their resolution.

## What is not yet measured

Whether a wider window removes the flicker without moving the orbit at which a
run first genuinely settles, and what the T42 rung converges to -- both from
cold and from a converted T21 state. Those are the arms this note is waiting
on.
