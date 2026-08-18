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

## What was decided, and what it changed

**The criterion now tests state storage. Decided 2026-08-17, CLIM-7.**

The quantity moved because the old one was measured against the thing it
proxies, not because a run failed on it. `|mean TOA| < 0.5 W/m2` was thresholding
a diagnostic carrying a structural offset larger than the threshold itself, and
`0.5` was picked early, never revisited, and applied to two different quantities.

**The threshold was derived and committed before it was applied to anything**,
which is the only way this could be taken without being the retrofit
`WORKFLOW.md` section 7 forbids. Two bounds, agreeing to 6%:

| bound | value |
| --- | ---: |
| what matters: 0.15 K offset tolerance / (0.128 K per orbit per W/m2 * 9.9 orbits) | 0.118 W/m2 |
| what is measurable: 10-orbit against 20-orbit windows of the same run, three runs | 0.111, 0.115, 0.116 W/m2 |

Below about 0.11 the estimator is measuring its own sampling noise; above 0.118
the criterion admits more drift than the temperature tolerance already forbids.
**0.12 W/m2** is where they meet, and it is four times tighter than what it
replaces. The physical half is derived FROM the temperature tolerance rather than
invented beside it, so the two criteria now bound one thing in two units.

**The outcome.** Re-assessed on the 77 orbits on disk, the baseline stores
-0.0547 W/m2 and passes. Three of the four runs pass;
`run_b014469b8091` does not, and it fails on the temperature-offset criterion
rather than on energy.

The ordering anomaly is gone, and that is the check on the change rather than a
by-product of it. Under the old criterion `run_524fbed77a9a` passed at
TOA -0.483 while its heat content ROSE at +0.084 W/m2, and the baseline failed at
-0.600 while storing -0.055. The criterion was passing the run further from
equilibrium. On storage the two are ranked the right way round: the baseline at
-0.055 sits nearer zero than 524fbed at +0.084, and 524fbed is now the one closer
to its threshold.

**Reported TOA is still recorded, and so is its difference from storage**, as
`reported_toa_minus_storage_w_m2` in every assessment. That gap is CLIM-1 and is
open: the candidates are `rsut`, `rlut`, or an atmospheric heating neither flux
diagnostic books. A criterion that stopped reading the number would also stop
anyone noticing when it changed.

**What it settles for A2.** The baseline's failure was an instrument fault, not a
spin-up fault, and no amount of further model time would have moved it. What
"final" still requires is the physics that moves the mean -- PHYS-1, the spectrum
declaration, dust -- and none of that is equilibration.

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
