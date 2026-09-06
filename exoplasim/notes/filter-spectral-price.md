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

That rate is what `scripts/check_consistency.py` solved against the flow's own
cascade rate `x/tau_vorticity` to get the confinement crossover
`x = (dt / (2 kappa tau))^(1/(gamma-1))`, and it is what made the filter look two
to three orders of magnitude stronger than the derived hyperdiffusion.

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
a spin-up window and not comparable. That run has no row in
`exoplasim/runs/INDEX.json` and no payload (`archive/runs/RECORDLESS.json`). Its
arm survives in `exoplasim/analysis/filter_spectral_cost_arms.json`, which
carries seven windows against the twenty-two and eighteen of the counted arms
and so corroborates the exclusion; the orbit range itself, and the build and
executable the arm ran on, are not recoverable from any artifact.

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

- **The confinement criterion cannot resolve kappa, and the gamma pair is what
  says so.** The FILTER's own bite point goes as `kappa^(-1/16)`, so a factor of
  two in kappa moves it by 0.036 of the truncation weaker and 0.035 stronger.
  That is below the quantisation at T21, 1/21 = 0.048, and ABOVE it at T42,
  1/42 = 0.024, so on the filter's law alone the criterion would read a factor
  of two at T42 and at every rung above it. It does not, because what a run
  measures is not the filter's bite but the TOTAL departure from the inertial
  range, and the hyperdiffusion sets the rest of it. The gamma pair is the one
  arm that moves the filter's bite: 0.675 to 0.822 predicted, 0.688 to 0.714
  measured, a transfer of 0.18. Applied to a kappa doubling that is 0.0065 of
  the truncation, below the quantisation of every rung through T127 (1/127 =
  0.0079) and 1.1 wavenumbers at T170. So the criterion in `world-6on`'s
  acceptance -- a kappa chosen against the spectral confinement criterion -- is
  not achievable at the gamma the model runs. The transfer is ONE calibration
  point and the conclusion holds for any transfer below 0.65.
- **The depth at 0.8N can, and only at T42.** The depth is LINEAR in kappa, so
  the two directions are not the same size: halving kappa is worth 0.098 dex and
  doubling it 0.196, and 0.81 of each after the calibration below. Against a
  three-sigma bar of 0.062 to 0.077 dex at T42 that clears in both directions;
  against 0.148 to 0.168 dex at T21 neither direction does. Eighteen two-orbit
  windows at T42 and twenty-two at T21, so this is the rung and not the sample.
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

## The top of every rung is dead, and the filter is most of why

The last wavenumbers of every run measured sit at floating-point roundoff. The
last one standing 100x clear of the measured roundoff floor is m=20 of 21 at
T21 and m=41 of 42 at T42. This is why the depth is read at 0.8 of the
truncation: read at the truncation it returns the roundoff floor, which is what
the first version of the instrument did.

**The cause is the filter, not the mode count.** Between m=35 and m=42 at T42,
which is 0.83 of the truncation to the truncation, the eddy kinetic energy falls
10.3 dex. The four contributions, three of which are arithmetic:

| | dex over m=35 to m=42 |
| --- | ---: |
| the filter, `f^2 = exp(-2 kappa x^gamma)` at kappa 8, gamma 16 | **-6.57** |
| the hyperdiffusion, by difference | -2.68 |
| the mode count, at most | -0.90 |
| the run's own inertial slope, -2.04 | -0.16 |
| measured | **-10.32** |

A triangular truncation holds `n = m .. NTRU`, so zonal wavenumber `m` carries
`NTRU - m + 1` meridional modes and `m = NTRU` carries one against `NTRU` at
`m = 1`. That is a factor of 42 at T42, which is 1.6 dex across the whole
spectrum and 0.9 dex across this band. **It cannot put anything at roundoff**,
and the claim it could -- that the mode count alone empties the top whatever the
damping does -- is wrong by six decades. The mode count is the smallest of the
three mechanisms here.

The filter's share is scale-free, so the shape repeats at every rung: over the
same fraction of the truncation, 0.81N to N, `f^2` falls 6.7 dex at T21.

## How many zonal wavenumbers a rung carries

CARRIED means the filter leaves at least half of what the flow's own inertial
range would put there, which is the same factor two `spectral_tail.py` calls a
bite and is not a new criterion. Under `f^2` that is
`m <= (ln 2 / (2 kappa))^(1/gamma) * NTRU`, a fixed FRACTION of every rung: 0.822
at kappa 8 and gamma 16. `exoplasim/scripts/filter_spectral_cost.py --rungs`
writes the table to `exoplasim/analysis/rung_spectral_budget.json`, reading no
run -- both halves are arithmetic on the truncation and on `config/planet.yaml`,
which is the only way to answer it at a rung nothing has been integrated at.

It is an UPPER BOUND, and the reason is the spectrum's coordinate: the filter is
a function of total wavenumber `n` and the count is in zonal wavenumber `m`,
which sums over `n >= m`. `f` decreases in `n`, so a wavenumber this says is not
carried is certainly not carried, and one it says is carried may still be
depleted by the modes above the diagonal.

So T21 carries 18 zonal wavenumbers of 22 and T42 carries 35 of 43.

**What that costs is not its share of the wavenumbers, and it does not distort a
rung comparison.** Legendre work at zonal wavenumber `m` is `NTRU - m + 1` modes
over every Gaussian latitude, so the uncarried band is the small end of a
triangle: 18 per cent of the wavenumbers at every rung and 3.3 to 4.0 per cent
of the Legendre work. The spread across the whole ladder is under a percentage
point, because the filter is scale-free in `n/NTRU` -- so the dead top cancels
out of a rung-to-rung ratio rather than biasing one, and no cost accounting has
to carry it. The roundoff-dead top is smaller again: one wavenumber at each of
T21 and T42, which is 0.40 and 0.11 per cent of the Legendre work.

**The accounting that IS wrong is comparing rungs by their truncation.** The
resolved degrees of freedom are `(NTRU+1)(NTRU+2)/2` -- 253 at T21 and 946 at
T42, a factor of 3.74 where the truncations differ by 2. SPAT-11's own recorded
Legendre ratios already use modes times latitudes and are right: this table
reproduces 7.9x at T85 and 26.2x at T127 against T42 exactly. It is the
"42 usable zonal wavenumbers against 21" framing that does not survive, not the
cost model underneath it.

**And it bounds what a spectral diagnostic can be read at.** This instrument
reads its depth at 0.8 of the truncation and REFUSES unless the spectrum there
stands 100x clear of the measured roundoff floor. `spectral_tail.py` now keeps
the band above the model's truncation for the same reason -- it is the only
place the floor is visible -- and ends its tail fit at the last wavenumber
standing that same margin clear of it, reporting the band it fitted over beside
the slope. The two criteria never needed it and are unchanged: the bite point is
the FIRST departure below the inertial range and the pile-up test looks for
excess above it, so roundoff at the top cannot move either. The tail slope did,
and a tail slope fitted to `m = NTRU` is a measurement of the distance from the
run's inertial range down to the floating-point floor rather than of the
damping. Both scripts take the 100x margin from one definition, in
`spectral_tail.py`; `smoke_test.py` holds the fit to a synthetic tail of known
slope, which the old band cannot return.

## Every arm here is attributable, which is what makes it usable

Each pair above is same-binary BY SHA, taken from the runs' own manifests: the
two T21 arms share `e7b54a87` and the two T42 arms share `68ad59e5`, and each
pair also shares its cold start, its seed, its surface fields by sha and its
stellar spectrum. That is stated because the neighbouring instrument does not
have it: the 22 entries in `analysis/stability_probe.json` name no executable,
no sha and no template, and the file's single `generated` field is rewritten on
every write, so they cannot be tied to a build at all (world-qnue). Nothing in
this note calibrates against them.

## What the confinement gate asserts now, and what it stopped claiming

`scripts/check_consistency.py`'s filter confinement row is a CONFIG gate: it
reads `config/planet.yaml` and has no run to look at. Under the refuted law it
solved a crossover with `dt`, `tau_vorticity` and the rung in it, and it
disagreed with the measured bite in both directions at once -- 0.625 against
0.561 at T21 and 0.655 against 0.714 at T42.

It now asserts the filter's own half-attenuation point,
`(ln 2 / (2 kappa))^(1/gamma)`, which is what the surviving law gives and which
carries no timestep and no rung. **The floor is unchanged at 0.60.** The quantity
moved; the bar did not. Re-choosing the bar to suit a corrected quantity would be
choosing a criterion after seeing the result it judges.

**It no longer claims to predict a run's bite point, and it must not be read as
doing so.** What a run measures is the total departure from its own inertial
range, which carries the hyperdiffusion and the spectrum's own steepening as
well, and the hyperdiffusion is the larger of the two levers. That criterion is a
property of a RUN and it has an instrument: `spectral_tail.py` fails a run whose
bite point sits below 0.6 of the truncation. A config gate that asserts a
measured quantity from arithmetic is a gate with no right answer available to
it.

**What the gate can see is gamma.** `x_f` goes as `kappa^(-1/gamma)`, so at gamma
16 the 0.60 floor holds up to kappa 1229 and the row is in practice a guard on
the exponent -- which is the degree of freedom that moves it, and which was 8
until it was derived. The gate reports that ceiling beside its verdict so the
margin is a number, and reports the filter's own depletion at 0.8 of the
truncation beside it because that one IS linear in kappa.

## What this leaves undecided, and it is more than it was

The weakest kappa that still RUNS. Nothing here touches it: the trap boundary is
a stability question, the filter buys timestep by a mechanism this note does not
address, and `physics-filter-stability.md`'s grid is the only instrument that
discriminates kappa at all.

So the spectrum cannot choose kappa and the trap can -- and the trap's current
evidence cannot be attributed to a build. `filter_kappa` 8.0 therefore rests on
a boundary measured under an unknown executable, on damping that reached one
model level in ten, and at gamma 8 where the model runs gamma 16. **kappa 8 is
not settled; it is determined by the only criterion that can see it, and that
criterion has to be re-taken.** What this note removes is the possibility that
the spectrum could settle it instead. world-37tn carries the re-measurement and
world-qnue the attribution.
