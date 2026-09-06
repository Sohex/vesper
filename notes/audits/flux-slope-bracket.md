# The flux-to-kelvin bracket on the configured build, and the criterion that could not judge it

*Worldbuilding frame: measurements of the Vesper project's simulated climate,
made with an ExoPlaSim fork. Nothing here is about the real world. Measured
2026-09-06 on build `canonical-10m-carve2`. world-fciz.*

`lib/sensitivity.py` declares one flux-to-kelvin sensitivity for the whole
project, dT/df near the baseline flux, measured as a secant between two
converged runs. Its declared pair was on `canonical-10m-base` and its currency
checks had been failing on that since the build moved. This note is the pair
bought to replace it, the measurement they give, and the reason the declaration
does not yet rest on them.

## Both arms had to be bought, and that is not what the row assumed

The row was written expecting one purchase: a warm arm at flux ratio 1.000 to go
with `run_5994d1f9624e`, which is on the configured build. That pairing is not
available. `_verify_bracket` refuses two arms whose `geography` digests differ,
and the digest covers every staged surface field:

| | `run_5994d1f9624e` | staged now |
| --- | --- | --- |
| geography digest | `007497a0` | `75475656` |
| surface code 1720, the SPAT-5 subaerial area fraction | not staged | staged, and required of every run |
| executable | `fdf01605` | `fe21b1b6` |

Nine model sources changed between them, over SPAT-5's coastal tiling and the
three-point wet-soil mixing. A bracket taking one arm from each side of that
measures the surface and the model as well as the flux. Both arms are therefore
bought together, and `exoplasim/scripts/run_flux_bracket.sh` buys them together
so the pairing cannot be assembled out of whatever happens to be lying around.

## The pair

Identical in structure and in everything but the flux they were given: orbits
0-69 at low I/O declared `spinup`, orbits 70-129 at clean I/O declared
`post_equilibrium_climatology`, assessed over the last 55 orbits of the clean
block. 55 is `assess_convergence.py`'s derived default window, the shortest at
which the offset criterion is meant to discriminate rather than flip; the clean
block is 60 so the window clears the I/O join that `segments.py` refuses a window
to span. **The window and the acceptance criterion were fixed before any number
was seen**: all six convergence criteria on each arm over orbits 75-129.

| | cold | warm |
| --- | --- | --- |
| run | `run_c919391cf715` | `run_9d5dbf9bd3d9` |
| flux ratio | 0.945 | 1.000 |
| assessed window | 75-129, clean I/O | 75-129, clean I/O |
| geography | `75475656` | `75475656` |
| executable | `fe21b1b6` | `fe21b1b6` |
| fitted asymptote | 280.1690 +/- 0.1887 K | 289.7269 +/- 0.0929 K |
| settled-window mean | 280.2544 K | 289.7813 K |
| last ten of the window | 280.2420 K | 289.7597 K |
| sea-ice mean fraction | 0.0753 | 0.0190 |
| six-criteria verdict | fails `extrapolated_offset_lt_0.15_k` | passes |

## The measurement, which three estimators agree on

| K per unit flux ratio | estimator |
| ---: | --- |
| 173.78 | the two fitted asymptotes |
| 173.22 | the two settled-window means |
| 173.05 | the last ten orbits of each assessed window |

The chord is 9.5579 K over 0.055 in flux ratio. The three span 0.73 K per unit
flux ratio, and two of them use no fit at all: they are direct averages over 55
orbits of settled data on each arm. **The slope is far better determined than
either arm's asymptote is**, because the two arms' fits err in the same direction
and the difference does not carry their full errors. Propagating the two
asymptote half-widths gives 3.82 K per unit flux ratio, five times the span the
estimators actually show.

The modelled sea-ice mean fraction falls from 0.0753 to 0.0190 across the chord,
so the ice-albedo feedback over that retreat is inside the number rather than
outside it.

## Why the declaration does not yet rest on this pair

`_verify_bracket` requires each arm's convergence report to carry
`sufficiently_equilibrated_for_worldbuilding`. The cold arm's does not. It fails
one criterion of six, and the failure is not drift:

| | cold | warm | threshold |
| ---: | ---: | ---: | ---: |
| `abs_temperature_slope` standard error | 0.0012 | | 0.05 |
| `abs_state_storage` standard error | 0.0101 | | 0.12 |
| offset criterion, exponential-fit form | 0.2741 | 0.1473 | 0.15 |
| offset criterion, drift form | 0.0418 | 0.0519 | 0.15 |

Every direct measure of settling passes with an order of magnitude to spare. The
offset criterion's statistic is `|offset| + its standard error`, and on the cold
arm the offset is -0.0854 K, inside the threshold; what fails is the standard
error, 0.1887 K, which is the fitted asymptote's own.

**Neither arm's offset criterion resolves, and both reports say so.** The
resolving bar is a standard error of at most a third of the threshold, 0.05 K.
The cold arm's is 0.1887 and the warm arm's is 0.0929, so both rows read
`resolves: false`. The warm arm's pass is by 1.8 per cent of the threshold. The
two arms are not on opposite sides of a settling question; they are on opposite
sides of a coin flip in an estimator that the report itself declines to trust at
this window.

## What decides which estimator is used, and why it is the wrong test

`assess_convergence.py` chooses between the two forms above on `fit_usable`,
which asks only that the asymptote and its half-width are finite, that the offset
is under 20 K and that the half-width is under 5 K. A half-width of 0.1887 K
clears that by a factor of 26, so the fit is used and its error becomes the
criterion's error.

The same function computes two things that would have decided it correctly and
consults neither:

- `relaxation_fit_identifiable`, whose test is `0 < tau <= fit_span`. On the cold
  arm it reads true, with the verdict "tau = 39.988 orbits is shorter than the 84
  orbits fitted, so the series carries the curvature that fixes it".
- `relaxation_orbits_fitted_standard_error`, which on the same fit is **40.35
  orbits against a tau of 39.99**. The fitted time constant is not
  distinguishable from zero at one sigma. The curvature does not fix it.

So the identifiability verdict is taken on a comparison that cannot fail for any
tau shorter than the span, however badly determined, and the branch that selects
the estimator does not look at the identifiability verdict at all.

**The consequence is that the criterion is easier to pass the worse the fit is.**
The two runs that pass on comparable clean windows do so because their fits
collapsed outright: `run_67323a923013` over 57 clean orbits and
`run_5994d1f9624e` over 46 both report `relaxation_orbits_fitted: null`, which
routes them to the drift form and its much tighter error. A fit that fails
completely yields a pass; a fit that half-succeeds yields a failure. On this cold
arm the drift form passes by a factor of 3.6 on exactly the same data.

This is the shape `docs/src/practice/failure-modes.md` class 34 describes: the
instrument's own scatter is larger than the effect it is asked to resolve, so the
number it returns is noise however tidy it looks.

## What was not done, and why

The cold arm was not extended. The report will not price what would close the
criterion in its exponential-fit form -- `window_orbits_for_offset_criterion`
prices the drift form only, and says so -- so no extension length can be fixed in
advance, and buying orbits until a verdict flips is not a criterion. Extending
is also the wrong direction on the mechanism above: a more settled run has less
curvature, which drives the fit toward collapse and toward a pass by the fallback
rather than by the run being better.

The bracket was not re-declared on the arms' low-I/O blocks either, and it could
not have been. **There is no window on which both arms pass**, and the four
assessments the purchase took say why:

| window | regime | cold | warm |
| --- | --- | --- | --- |
| 15-69 | low I/O | passes all six | **fails** `extrapolated_offset` |
| 75-129 | clean I/O | **fails** `extrapolated_offset` | passes all six |

The two arms fail on opposite windows, on the same criterion, and on nothing
else: every other criterion passes on every one of the four. The arms are
identical in build, staged surface, executable, structure and window indices and
differ only in the flux they were given, so a verdict that flips between them and
between windows with no pattern is the estimator and not the settling. On the
low-I/O window the warm arm reads an offset of +0.0864 K with a half-width of
0.1090, from a fit whose tau is 9.88 +/- 17.41 orbits -- again a time constant
indistinguishable from zero, again admitted by `fit_usable`.

That is the whole of the case for world-jejw, and it is why no choice of window
rescues this pair.

## The warm arm is the noisiest settled run this project has, and that moves the window

`assess_convergence.py` sizes its default verdict window on a declared orbit
scatter that has to bound every settled production report, and
`check_convergence_bounds` refuses when a report exceeds it. The warm arm reads
0.144917 K over its 55-orbit clean window against a declared 0.124, so the bound
is raised to 0.145 and the anchor moves to it from `run_67323a923013`.

The window grows as the two-thirds power of the scatter, so **the default verdict
window goes from 55 orbits to 61**. The arms were bought with a 60-orbit clean
block, sized to hold the 55 that was the default when the purchase was specified.
They cannot hold 61. Re-buying this pair means a clean block of at least 66
orbits, not 60.

The cold arm reads 0.103598 over the same window indices. The 40 per cent
difference between the two arms is the flux and nothing else: same build, same
staged surface, same executable, same window, same regime.

## Cost, with the machine state it was taken under

*Measured 2026-09-06, T21 at 16 threads, one host lock held throughout.* 260
orbits across the two arms in 64 minutes wall clock, 23:22:04 to 00:26:05 MDT:
14.8 s per orbit averaged over both blocks, 15.4 s in the low-I/O spin-up and
about 16 s in the clean-I/O block. Host load sampled every 15 s throughout, 270
samples: mean 18.9, minimum 9.2, maximum 32.5 on 32 logical cores. The two arms
ran sequentially inside one lock, so the load above the model's own sixteen
threads is other agents.
