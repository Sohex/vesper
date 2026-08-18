# Is the baseline equilibrated? The window, and the instrument the criterion reads

*Measured 2026-08-17 on `run_8c2e1ff9ab5e`, the 0.945 T42 baseline on
`precarve-craton`, 77 annual outputs on disk (orbit indices 0 to 76). Opened as
CLIM-7: the run fails `|mean TOA| < 0.5 W/m2` while an earlier assessment had
recorded it as missing by 0.0014 W/m2, and those two statements cannot both
describe the same window.*

The question matters beyond bookkeeping. `WORKFLOW.md` section 6 A2 requires the
stellar cycle run to come after a baseline that is FINAL rather than merely
converged, so whether this run can honestly be called quasi-equilibrated is on
the project's critical path.

Two things were suspected and only one of them is true. The assessment does not
choose a favourable window. But the criterion it applies reads a diagnostic that
is offset from the physical quantity by more than the criterion's own threshold,
and that is the finding.

## Where the 0.0014 was evaluated

On **orbits 51 to 60**, when the run was 61 orbits long. The convergence report
committed at `02c916e~1` records `mean_toa_balance_w_m2 = -0.5014128`,
`completed_orbits = 61`, `window_orbits = 10`: a miss of 0.0014128.

It was stale before it was written down. The commit that quotes it, `02c916e`,
also *contains* the reassessment at 66 orbits, which reads -0.6010 on orbits 56
to 65 -- a miss of 0.1010, seventy times larger. The prose was written from the
number the run had before the commit's own work re-measured it, and the sentence
then propagated into `TASKS.md` and into `water-and-energy-closure.md`. Nothing
chose the window; a number outlived its measurement.

## The assessment does not shop for a window

`assess_convergence.py` takes the mean over the **last `--window` orbits**,
default 10, and there is no other window it could take: `output_files` requires
the annual outputs to be contiguous from year zero, the slopes and the means all
read `series[-w:]`, and the window is a command-line argument fixed before the
data is read. There is no search over start points and no selection among
candidate windows. On that count the suspicion is unfounded and should be
dropped.

**The free variable is the run length, not the window.** Re-assessing the same
run at every length it has passed through gives:

| run length | window | mean TOA, W/m2 | verdict |
| --- | --- | ---: | --- |
| 41 to 44 | trailing 10 | -0.4251 to -0.4946 | passes |
| 45 to 46 | trailing 10 | -0.5010, -0.5003 | fails by 0.0010, 0.0003 |
| 47 to 59 | trailing 10 | -0.4429 to -0.4957 | passes |
| 60 | orbits 50-59 | -0.5051 | fails by 0.0051 |
| 61 | orbits 51-60 | -0.5014 | fails by 0.0014 |
| 66 | orbits 56-65 | -0.6010 | fails by 0.1010 |
| 77 | orbits 67-76 | -0.5995 | fails by 0.0995 |

The criterion has been failed at **every** run length from 60 orbits onward, and
the 0.0014 is the smallest miss in that failing sequence. The three sign changes
before it are not the run approaching equilibrium and then leaving it. Per-orbit
`ntr` over orbits 55 to 76 has a standard deviation of 0.1697 W/m2, so the
standard error on a ten-orbit mean is 0.0537 -- and the criterion is being
applied at a distance from its threshold of 0.001 to 0.1. **A verdict decided
inside one standard error of the estimator is not a verdict about the run**, and
between run lengths 44 and 60 that is all it was.

Over the whole settled stretch, orbits 55 to 76, mean TOA is -0.5683 W/m2 while
surface temperature sits between 289.61 and 289.86 K with no trend. That is the
pair of statements the rest of this note is about.

## The criterion reads an instrument with an offset larger than its own threshold

`|mean TOA| < 0.5 W/m2` is a proxy. What it is trying to test is that the planet
is neither gaining nor losing energy, and the reported top-of-atmosphere net
radiation is the model's own statement about that, never checked against the
planet.

It can be checked, from the prognostic state, and
`exoplasim/scripts/close_state_energy.py` does it: the ocean mixed layer, the
sea-ice and snow mass as latent heat, the soil column, the atmosphere's enthalpy
and its column vapour, differenced in time. Conservation says that derivative
equals the mean TOA net. Measured on this run, orbits 67 to 76, the block written
with `NLOWIO = 0` throughout:

| | W/m2 |
| --- | ---: |
| mean TOA net, `ntr` | -0.5995 |
| d(planetary heat content)/dt | -0.0547 |
| **residual** | **-0.5448** |

The planet is not losing 0.6 W/m2. It is losing 0.05, which is a ninth of the
criterion's threshold. `water-and-energy-closure.md` carries the same measurement
on four runs spanning 7.8 K of mean temperature and finds the residual at
-0.573 +/- 0.035 W/m2 -- a structural offset in the diagnostic, not a property of
this run.

So the criterion is failing by 0.10 on a number that is 0.57 too negative. Every
recent miss is smaller than the offset, and the sign of the offset is the sign of
every miss.

## What this does and does not license

**It does not turn the miss into a pass, and the label does not move.**
`sufficiently_equilibrated_for_worldbuilding` is all-or-nothing, and it stays
`quasi_equilibrated` with `abs_mean_toa_lt_0.5_w_m2` named. The recorded verdict
is -0.6010 on orbits 56 to 65, a miss of 0.1010, taken when the run was 66
orbits long; re-assessed at the 77 orbits now on disk it would read -0.5995 on
orbits 67 to 76, a miss of 0.0995. Both fail, and neither is close. This project
calls a run quasi-equilibrated for missing by 0.004 W/m2; a run missing by 0.1 is
not entitled to better treatment because a separate measurement suggests the
threshold is being applied to the wrong quantity.

The assessment was deliberately **not** re-run here. It rewrites the run
manifest and would leave `INDEX.json` and `world_state.json` needing
regeneration, and the verdict it would produce is the same verdict on the same
criterion. The eleven-orbit staleness is worth recording rather than quietly
refreshing, because it is the same shape as the defect this note opened on.

**It does not license amending the criterion here.** The threshold and the
quantity a criterion tests are fixed before results are seen
(`WORKFLOW.md` section 7). Choosing a new quantity -- state storage rather than
reported TOA -- immediately after measuring that the new quantity passes and the
old one fails is the exact move that convention forbids, whatever the physics
says. If the criterion should read storage, that is a decision to take
deliberately, with its threshold argued from the estimator's noise rather than
from this run's answer, and it belongs to CLIM-7's remaining work rather than to
this note.

**What it does settle, for A2.** The baseline's failure to converge is an
instrument fault and not a spin-up fault. There is no missing model time here:
running further will not move -0.5995 to -0.4999, because the number it would
have to move is not the number the planet has. Whatever else "final" requires --
and A2 is explicit that it is stricter than "converged", and that anything
outstanding which moves the mean has to land first -- **equilibration is not what
is holding the cycle run back.** That was worth knowing before spending a
hundred more orbits on it.

## What was checked and is not the answer

- **The window rule.** Fixed, singular, and applied without choice. See above.
- **A missing storage reservoir.** The model has one ocean layer
  (`oceanmod.f90:15`, `NLEV_OCE = 1`), `glac` is identically zero everywhere on
  this run, and every other reservoir is in the closure above. Nothing is
  quietly absorbing the difference.
- **The seasonal storage swing.** The gap oscillates by about +/-7 W/m2 within an
  orbit, so any single-bin reading is noise. Every number here is a full-orbit
  mean over ten or more orbits, and the residual per orbit runs -0.453 to -0.595
  across the ten clean ones.
