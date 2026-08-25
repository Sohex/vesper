# How long the model takes to converge, and what a window can see

*Worldbuilding frame: the Vesper project's climate model, an ExoPlaSim fork.
Nothing here is about the simulated planet.*

Two different things are written down here and they have different standing. The
first is OPERATIONAL EXPERIENCE with no artifact behind it, recorded because it
existed only in one person's head and it is what says how long a commissioning
run has to be. The second is measured, with the run named.

## Operational experience: how many orbits

*Recorded 2026-08-25 from the project owner's own experience of running this
model. NOT a measurement, and there is no artifact to check it against. Treat it
as the number to design an experiment around and not as one to quote.*

- **A cold start at T21 converges in about 70 orbits**, give or take, and the
  spread is real: it depends on the initial condition and on the cold-start
  seed.
- **A reconvergence is much shorter, generally 10 to 20 orbits.** That is what a
  run needs after a timestep change or after a resolution conversion, because
  the state handed in is already near the attractor and only the support or the
  truncation error has moved.

The escalation route in `docs/src/pipeline/sequencing.md` section D alternates
resolution and timestep and reconverges at each change, so it pays the long
number once and the short one at every step after.

## Measured: what a twenty-orbit window can and cannot see

*Measured 2026-08-25 at 5dcdf64b on `run_4769c35837df` -- T21, dt 45, ten levels,
sixteen threads, eight-byte, 85 orbits from cold in one call. One run; no
comparison is involved in anything in this section.*

The four non-overlapping 20-orbit means of global-mean surface temperature:

| orbits | mean, K |
| --- | ---: |
| 5-24 | 291.355 |
| 25-44 | 291.835 |
| 45-64 | 291.757 |
| 65-84 | 291.651 |

Spread 0.210 K, range 0.480 K. The run is past the cold-start length above, so
this is not the climb: it is the model's own low-frequency variability, and the
per-orbit series carries a lag-1 autocorrelation of 0.615, which puts about five
independent samples in a 20-orbit window.

**So a criterion evaluated on one 20-orbit window is measuring variability as
much as state.** A window's standard error goes as `sqrt((1+r1)/(1-r1)/n)`, so
reaching 0.05 K against a 0.10 K orbit-to-orbit scatter at that autocorrelation
needs about seventeen independent samples -- roughly seventy orbits, not twenty.
world-wdsk carries the question of what the window length should be, world-omn
what the window's noise costs a convergence verdict, and world-yj9o the same
fact reaching `compare_equilibria.py`, which takes its standard error on the raw
orbit count and so understates it by about a factor of two.

## Why this note exists rather than a number in a criterion

An experiment is designed against a convergence length, and until this was
written down there was none to design against: the four-byte precision arm of
clim-59 was specified at 40 orbits with the last 20 as its window, which is
inside the cold-start climb, and the result reversed sign when it was re-run at
85. `notes/audits/single-precision-spin-up.md` records what that cost.
