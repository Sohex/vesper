# What the physics filter costs the spectrum, and why its price carries no timestep

*Worldbuilding frame: a numerical experiment on the Vesper project's climate
model. Nothing here is about the simulated planet. Measured 2026-08-25 on runs
already on disk, by `exoplasim/scripts/filter_spectral_cost.py`, whose criteria
were committed before any arm was read.*

`model.filter_kappa` was chosen on the trap boundary of
`physics-filter-stability.md` -- the longest stable step at every rung measured
-- and nothing had priced it against the spectrum it damps. This note prices
it, and in doing so refutes the rate law the price was expected to follow.

## The claim under test

`resolution-tuned-parameters.md` and `config/planet.yaml` both carried the same
arithmetic: the filter multiplies each coefficient by `f(n) = exp(-kappa
(n/NTRU)^gamma)` at both transform directions every timestep, so a field is
multiplied by `f^2` per step and the implied damping rate is

    r(x) = 2 kappa x^gamma / dt,   x = n/NTRU

That rate is what `scripts/check_consistency.py` solves against the flow's own
cascade rate `x/tau_vorticity` to get the confinement crossover
`x = (dt / (2 kappa tau))^(1/(gamma-1))`, and it is what makes the filter two to
three orders of magnitude stronger than the derived hyperdiffusion.

**It has a timestep in it, and the spectrum does not.**

## The refutation: halving the timestep changes nothing

Two pairs, each differing in `MPSTEP` alone. Same binary by sha, same cold start
with the same declared seed, same surface fields, same damping, same filter.
Depth of the kinetic-energy spectrum below its own extrapolated inertial range
at 0.8 of the truncation, in dex, with the scatter across disjoint two-orbit
windows:

| rung | gamma | dt | measured depth | dt pair moves it by | the rate law predicts |
| --- | ---: | ---: | ---: | ---: | ---: |
| T21 | 16 | 45 | -1.418 +/- 0.049 | | |
| T21 | 16 | 22.5 | -1.409 +/- 0.056 | **+0.009, 0.1 sigma** | -0.296 |
| T42 | 16 | 45 | -0.602 +/- 0.020 | | |
| T42 | 16 | 30 | -0.609 +/- 0.026 | **-0.007, 0.2 sigma** | -0.169 |

The bite point does not move either: 0.567 against 0.552 at T21 (0.5 sigma) and
0.714 against 0.714 at T42. Both predicted moves are several times the
three-sigma bar the instrument reports, and neither arrives. At T21 the measured
sign is opposite to the predicted one.

`run_aedeb37bc7f2` is the third dt pair and is NOT counted here: its window is
orbits 10 to 23 of a cold start where the others are orbits 30 upward, so it is
a spin-up window and not comparable.

## Why: the filter is folded into the transform, not applied to the state

`legmod.f90` builds `skgpsp(n)` and `skspgp(n)` from `filterkappa` and
`nfilterexp` and folds them into every per-mode transform weight -- `fsp`,
`fgp`, `fmu`, `fmv`, `fmm`, `fmq` -- so no transform weight is ever used
unfiltered. What passes through those weights is the NONLINEAR TENDENCY:
`gp2fc` and `fc2sp` deliver into `spt`, `stt`, `sdt`, `szt`, `sqt`, and the
prognostic spectral state `sp, st, sd, sz, sq` is advanced from those in
spectral space. **No filtered value is ever stored back into the prognostic
state.**

So the filter attenuates the SOURCE that fills a wavenumber, by a fixed
multiplicative factor, rather than damping the state at a rate. A source scaled
by `f^2` against an unchanged sink equilibrates at

    spec(x) / law(x) = f(x)^2 = exp(-2 kappa x^gamma)

which has no `dt` in it at all -- and the tendency's own `dt` factor cancels
exactly, because a tendency applied for a shorter step is applied more often.
That is the null above, and it is a property of where the filter is applied
rather than of how strong it is.

## The gamma pair confirms the replacement law

Two arms at T42 and dt 30, differing in `NFILTEREXP` alone, all-level derived
damping, 12 two-orbit windows each:

| | gamma 8 | gamma 16 | moves by |
| --- | ---: | ---: | ---: |
| depth at 0.8N, dex | -1.391 +/- 0.029 | -0.602 +/- 0.020 | **+0.789, 22 sigma** |
| `exp(-2 kappa x^gamma)` predicts | -1.166 | -0.196 | +0.970 |
| the rate law predicts | -2.256 | -1.493 | +0.763 |
| bite point | 0.688 +/- 0.007 | 0.714 +/- 0.000 | +0.026, 3.8 sigma |
| eddy KE removed | 2.64% +/- 0.19% | 1.87% +/- 0.20% | -0.8%, 3.2 sigma |

Both laws get the gamma response roughly right, because gamma sits inside the
exponent in both. Only one of them survives the timestep test. The measured
response is 0.81 of what `f^2` alone predicts, and the shortfall is the same at
both gammas -- the residual depth of 0.22 dex at gamma 8 and 0.40 dex at gamma
16 is what is left after the filter, which is the hyperdiffusion and the
spectrum's own steepening away from a fit taken over the middle band.

The bite point at gamma 8 is the sharper test of the two laws, because they
disagree by a factor of nearly two there: `f^2` puts it at
`(ln2 / 2 kappa)^(1/gamma)` = 0.675 and the rate law puts it at 0.381. It
measures 0.688.

## The price of kappa, and what can and cannot see it

Under `f^2` the filter's contribution to the depth at a fraction `x` of the
truncation is `2 kappa x^gamma / ln 10` dex, and the bite point is
`(ln 2 / 2 kappa)^(1/gamma)`. At gamma 16 those give

- **The confinement criterion cannot resolve kappa.** The bite point goes as
  `kappa^(-1/16)`, so a factor of two in kappa moves it by 0.030 of the
  truncation. That is below the instrument's quantisation at T21 (1/21 = 0.048)
  and at the quantisation at T42 (1/42 = 0.024) against a measured three-sigma
  bar of 0.000 to 0.021. The criterion in `world-6on`'s acceptance -- a kappa
  chosen against the spectral confinement criterion -- is not achievable at the
  gamma the model runs, and this is measured rather than argued.
- **The depth at 0.8N can, barely.** A factor of two in kappa is worth 0.136 dex
  before the 0.81 calibration and 0.110 after it, against a three-sigma bar of
  0.062 to 0.077 dex at T42 and 0.148 to 0.168 dex at T21. So a factor of two is
  about 1.5 bars at T42 on eighteen two-orbit windows and is BELOW the bar at
  T21.
- **The whole of kappa 8 is worth 0.196 dex at 0.8N** against a total measured
  depletion of 0.602 dex at T42 and 0.969 to 1.418 dex at T21. The filter is
  therefore a third of the depletion at T42 and a fifth at T21, and the
  hyperdiffusion and the natural steepening carry the rest.

## What the filter is not

**It is not what confines this model's spectrum at gamma 16.** The claim it
displaced -- that the filter is two to three orders of magnitude stronger than
the hyperdiffusion at every scale -- was computed from the rate law, at gamma 8,
and it is a gamma 8 statement even on its own terms: at `x = 0.5` the same
arithmetic gives 12.0 hours at gamma 8 and 3072 hours at gamma 16, against the
derived hyperdiffusion's 19551, so the ratio at half the truncation falls from
1632 to 6.4. Under the law that survives the timestep test the comparison is not
between two rates at all, because only one of the two mechanisms is a rate on
the state.

The direct evidence is in `resolution-tuned-parameters.md`, section "The
hyperdiffusion reaches every level now, and it moves the spectrum": a
hyperdiffusion change that reaches all ten levels moves this spectrum by 0.45
dex at ten sigma, which is more than twice what the whole of kappa 8 is worth.

## The top of every rung is dead, and it is not the filter that kills it

The last wavenumbers of every run measured sit at floating-point roundoff. The
last one standing 100x clear of the measured roundoff floor is m=20 of 21 at
T21 and m=41 of 42 at T42, and at T42 the wavenumbers below that fall away
across four decades in three steps. A triangular truncation gives `m = NTRU`
exactly one meridional mode against `NTRU` at `m = 1`, so the mode count alone
empties the top whatever the damping does. This is why the depth is read at 0.8
of the truncation: read at the truncation it returns the roundoff floor, which
is what the first version of the instrument did.

## What this leaves undecided

The weakest kappa that still RUNS. Nothing here touches it: the trap boundary is
a stability question, the filter buys timestep by a mechanism this note does not
address, and `physics-filter-stability.md`'s grid is the only instrument that
discriminates kappa at all. Since the spectrum cannot choose kappa and the trap
can, kappa 8 stands on the criterion that can see it -- which is where it was.
