# Audit: what a kelvin does to this world's water, and why the runoff answer is a residual

*Measured 2026-08-18 on `precarve-craton`, from the raw per-orbit output of the
two converged runs that share a surface: `run_bfa3f5269660` at flux 0.91 and
`run_524fbed77a9a` at flux 0.945, both carrying geography digest `3a17c498`.
Neither has a climatology built, so the fields are averaged from `MOST.*.nc`
directly. This closes BUDG-1 and supplies the one input BUDG-2 needed.*

The error budget prices every item in kelvin. The decision it has to inform is
the carve list, which is decided in millimetres of runoff. This is the
conversion between them, and the useful part of it is not the number but the
shape: **the two large terms very nearly cancel, so the answer is a small
residual and is known to about a factor of two.**

Findings are tagged **[numeric]** where computed here and **[inspection]** where
read out of code.

---

## 1. The secant, and the window it is taken over

**[numeric]** Land means are taken with Gaussian quadrature weights over cells
where `lsm > 0.5`, annualised to Earth years. That is deliberately the same land
mean `scripts/error_budget.py:land_water_balance` uses, because the response
measured here is multiplied by the amplification measured there and two
different land means do not compose.

| window | run | T (K) | P | E | P - E |
| --- | --- | ---: | ---: | ---: | ---: |
| last 5 orbits | run_bfa3f5269660 | 281.941 | 690.45 | 576.17 | 114.28 |
| last 5 orbits | run_524fbed77a9a | 288.996 | 799.52 | 680.08 | 119.44 |
| last 10 orbits | run_bfa3f5269660 | 281.922 | 690.88 | 576.76 | 114.12 |
| last 10 orbits | run_524fbed77a9a | 288.996 | 798.10 | 675.66 | 122.45 |
| last 20 orbits | run_bfa3f5269660 | 281.846 | 692.38 | 577.39 | 114.99 |
| last 20 orbits | run_524fbed77a9a | 288.811 | 797.81 | 673.80 | 124.00 |

All depths are mm per Earth year. The secants that follow from them:

| window | dT (K) | dP (%/K) | dE_land (%/K) | d(P - E) (%/K) |
| --- | ---: | ---: | ---: | ---: |
| last 5 orbits | 7.055 | +2.239 | +2.556 | +0.640 |
| last 10 orbits | 7.075 | +2.194 | +2.424 | +1.031 |
| last 20 orbits | 6.966 | +2.186 | +2.397 | +1.125 |

**The ten-orbit window is what `HYDROLOGICAL_RESPONSE_PER_KELVIN` carries**, at
+2.19 %/K on land precipitation and +2.42 %/K on land evaporation. Ten orbits
because that is the clean window this project uses everywhere else it averages a
converged run.

## 2. The finding: the runoff response is a residual, and it is the thing that is uncertain

**[numeric]** Read the two tables together. Land P and land E are stable across
the window: they move by 2.4% and 6.6% of themselves between the 5-orbit and
20-orbit answers. Their difference is not: **+0.64 to +1.13 %/K, a factor of
1.8, from windows that agree to within a few percent on both inputs.**

That is arithmetic rather than noise in the climate. Runoff is 15% of
precipitation here, so each amplified channel is about twenty times the residual
between them, and a 1% wobble in either input is a 13% wobble in the answer.
Which is also why the runoff response cannot be inherited from Earth intuition:
it is a property of how dry this world's land is, not of its climate
sensitivity.

The consequence for the budget is bounded and stated rather than fixed: the
kelvin-to-runoff conversion is good to roughly a factor of two. `error_budget.py`
declares itself a factor-of-two instrument for ORDERING, so this is inside its
tolerance and outside the tolerance of anything that wanted a number.

## 3. The trap: the amplification is not the sensitivity

**[inspection]** `error_budget.py` reports that runoff amplifies a fractional
precipitation change by `P/R` and a land evaporation change by `E/R`, currently
6.5 and 5.5. It is tempting to multiply the precipitation response by the first
of those and call the product the temperature sensitivity of runoff. It is not,
and the error is an order of magnitude.

    precipitation-only   +2.19 %/K * 6.51  =  +14.3 %/K      WRONG for a kelvin
    a kelvin             +2.19 * 6.51 - 2.42 * 5.51 = +0.92 %/K

A kelvin raises land evaporation slightly FASTER than precipitation on this
world, so the two amplified terms nearly annihilate. The amplification is the
right factor for a perturbation that changes the hydrological cycle without
changing the temperature -- which is what a dust suppression or a circulation
change can be, and what
`aeolian/analysis/dust_runoff_sensitivity.json` deliberately measures -- and it
is the wrong factor for every item in the error budget, because every one of
them is priced in kelvin.

The `+0.92 %/K` above uses the baseline state's `P/R` and `E/R` rather than the
cold run's, which is the project's rule that a correction is estimated against
the state you are in. Taking the secant directly on the 0.945 run's own P and E
gives +1.03 %/K, and the difference between the two is smaller than the window
spread in finding 2.

## 4. What is NOT measured: the lake limb

**[inspection]** The carve criterion has three water terms and this measures
two. Precipitation and land evaporation set the runoff the criterion divides by;
open-water Penman evaporation sets the numerator, and no measurement of its
response to temperature exists on this world. Both converged runs would support
one, but Penman needs the lowest-level air, the wind and the surface radiation
from a built climatology, and neither run has one.

`HYDROLOGICAL_RESPONSE_PER_KELVIN` therefore carries a bracket for that channel,
with declared ends:

* **low, the land rate.** A lake responding no faster than the moisture-limited
  ground beside it is the weakest response that is physically arguable.
* **high, 6.7 %/K**, the convexity of saturation vapour pressure this project
  already measured under PHYS-5. At fixed relative humidity the deficit driving
  Penman's aerodynamic term grows at that rate; the radiative term does not grow
  at all, so the Penman response at fixed radiation and wind is strictly inside
  this end.

The bracket is worth tens of basins per kelvin, which is material against a
population of a couple of thousand, and it reaches the basin count only. Runoff
is `P - E_land` and contains no lake evaporation, so the runoff column carries
no bracket at all. That asymmetry is why `error_budget.py` returns a single
runoff number and a RANGE of basins.

## 5. Direction, which is the part that survives a re-measurement

**[numeric]** Warming closes basins. Both evaporation channels rise faster than
precipitation, so a warmer world sends less water to a lake surface that
evaporates more, and the overflowing population shrinks. Every budget item that
warms the world therefore carries a negative basin count and every item that
cools it a positive one.

That is the same sign as the antitone verdict map in `WORKFLOW.md` section 4:
carving darkens the land, the world warms, and basins that were marginal stay
closed. The two arguments are independent and agree, which is the only thing
about this measurement that does not depend on the window it was taken over.

## 6. Reproducing this

Average `pr`, `evap` and `ts` over the last N `MOST.*.nc` of each run, take the
land mean with `numpy.polynomial.legendre.leggauss(nlat)[1][::-1]` weights over
`lsm > 0.5`, scale m/s to mm per Earth year, and difference the two runs against
their global-mean `ts`. The runs are named at the top of this note; both are
`equilibrated_for_worldbuilding` in `exoplasim/runs/INDEX.json` and both sit on
the same surface fields, which is what makes the pair a flux secant rather than
a comparison of two different worlds.

**Do not pair either of them with `run_8c2e1ff9ab5e`.** That is the baseline,
and it carries a different geography digest because the lakes are composited
into its albedo. A secant across it would vary the surface and the flux
together.
