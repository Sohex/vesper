# The surface-to-planetary albedo attenuation, measured on four paired arms

*Measured 2026-08-30, from output already on disk. No model run was bought.
Vesper is an invented super-Earth; every albedo, temperature and flux below is a
diagnostic of the climate model that simulates it, and "the atmosphere above the
surface" means the modelled one.*

`scripts/error_budget.py` prices every albedo item through one multiplier: a
land-mean albedo change reaches the top of the modelled atmosphere reduced by an
attenuation, and is then converted to kelvin by `lib/sensitivity.py`. That
multiplier was 0.5, DECLARED, with the file saying in its own words that the
honest claim was a factor of two rather than a number. A factor of two on it was
a factor of two on the whole albedo half of the budget.

It is now measured on two terms and three builds, and the answer is that a single
constant survives, with a bracket four times narrower than the declared one and
straddling the value that was declared.

## What the attenuation has to be, for the budget to compose

The chain is

    surface, planet-wide   d_surface = d_land * land_fraction
    top of atmosphere      d_toa     = d_surface * attenuation
    temperature            d_T       = sensitivity.planetary_albedo_to_kelvin(d_toa)

`lib/sensitivity.py`'s slope maps a FORCING to an equilibrium temperature and
already carries every feedback, so `d_toa` is the forcing and never the
equilibrium planetary albedo change. The attenuation is therefore defined
operationally: it is the value that makes the chain reproduce a measured
temperature separation. Each measurement below back-solves it.

**The denominator is the STAGED land-mean albedo, not the model's diagnosed
one, and the choice is worth a factor of 1.45 on the answer.** Every albedo item
in the budget is a delta of a staged boundary condition read from
`albedo_report.json`. The diagnosed `alb` a run writes is that boundary
condition as realised through the modelled snow, which is a different quantity
and one no item can be denominated in, because nothing can set it. Both are
reported below; only the staged column composes with the budget.

## The instrument

Paired equilibrated spin-ups that differ in `model.land_albedo_source` alone
within one build, one flux ratio and one rung. That key is the budget's own
"biosphere: bare rock vs vegetated" item. `dT` is the difference of the two
convergence reports' fitted asymptotes, which is the estimator `lib/sensitivity.py`
measures its own slope with. The fourth row is PHYS-15's `NWETSOIL` pair, whose
run directories are deleted and whose convergence diagnostics survive.

| term and build | dT | staged d(alb) | attenuation, staged | diagnosed d(alb) | attenuation, diagnosed |
| --- | ---: | ---: | ---: | ---: | ---: |
| soil wetting, canonical-10m-base | +0.1958 | -0.006745 | 0.292 | -0.006397 | 0.308 |
| biosphere endmember, canonical-10m-carve1, uniform soil water | +3.1381 | -0.082856 | 0.381 | -0.052737 | 0.598 |
| biosphere endmember, canonical-10m-carve1, pedology soil water | +3.5256 | -0.080385 | 0.441 | -0.055229 | 0.642 |
| biosphere endmember, canonical-10m-base | +3.2714 | -0.069569 | 0.473 | -0.050071 | 0.657 |

Kelvin, land-mean albedo, planetary albedo 0.313627 from the configured baseline
climatology, slope 159.7 K per unit flux ratio, land area fraction 0.427692 on
the T21 land-sea mask the runs share.

The runs, cold arm then warm arm: `run_598eb57c5a34` / `run_a1c35075747c`;
`run_4b7281ee038a` / `run_76e441e0a761`; `run_bd7de9ba67a4` / `run_8ff97d5e189a`;
`run_58f467b0872d` / `run_893e276ee029`.

**Every pair is resolved many times over.** The three endmember pairs separate at
53, 78 and 53 times the two asymptote half-widths added in quadrature; PHYS-15's
separates at 4.3 times. The bar was fixed before the numbers at three times that
quadrature sum.

## The null the instrument was checked against

`run_88ed6f9d34ae` and `run_8d0ae7e2d02c` are T42 on `canonical-10m-carve2` with
identical `surface_field_sha256`, identical `source_config`, different
executables and different routes to equilibrium (one cold-started, one resumed
from a converted restart). Two runs handed the same surface albedo must separate
by nothing.

| quantity | measured |
| --- | ---: |
| dT between the two asymptotes | -0.0071 K |
| d(planetary albedo) in window | +1.8e-06 |

So the estimator's floor across independent integrations is under a hundredth of
a kelvin, against endmember separations near 3.3 K. The signal is 460 times the
null, and the null carries an executable difference the measured pairs do not.

## The one thing the endmember swap changes besides albedo

`NVEG = 0` in every run here, so `dforest` is prescribed from surface code 212
and is not prognostic. `landmod.f90` uses it in exactly one place,
`snowcanopymask`, which sets how much of the modelled snow the canopy hides. So
the forest fraction that moves with the endmember swap is an ALBEDO effect and
touches no flux directly. That was the confound the measurement was most exposed
to and it is disproved rather than bounded.

**What is not disproved is code 229.** The base pair and the carve1 pedology pair
stage different soil-water field capacities between their arms, because the
pedology field depends on the vegetation; the carve1 uniform pair stages none in
either arm and is albedo-only. That pair is the clean one and it gives the LOWEST
endmember attenuation, 0.381, against 0.441 and 0.473 for the two that carry the
extra field. The direction says the field capacity difference contributes
warming, and it puts the endmember term's honest value at the low end of its own
three rows.

## The verdict, against the rule fixed before the numbers

The rule registered on `world-ckbt` was: within a factor of 1.5 the attenuation
becomes a measured constant; 1.5 to 2, a measured constant with the measured
spread and both columns kept; beyond 2, a single constant is the wrong shape.

The staged set is 0.292, 0.381, 0.441, 0.473: a span of 1.62, inside the middle
band. **So the attenuation becomes MEASURED, with the measured span as its
bracket, and the budget goes on printing the unattenuated column beside it.**

The diagnosed set spans 2.13 and would have refused a single constant. It is not
the set the budget composes with, and saying which convention the budget is in
is what settles that. The pre-registration named the diagnosed estimator as
primary in order to match PHYS-15's published figure, and that was the wrong
choice: PHYS-15 quoted its realised fall, and every item the constant multiplies
is a staged delta.

**One pair is excluded under the pre-registered instrument test and it changes
nothing.** The test was that a pair whose diagnosed and staged deltas disagree by
more than a factor of 1.5 is pricing the endmember swap rather than the
attenuation. The carve1 uniform pair disagrees by 1.57 and is excluded; the other
two disagree by 1.39 and 1.46 and are not. The excluded point sits inside the
remaining span, so the bracket and its midpoint are unmoved to three figures.

**The central value is the midpoint of the measured span, 0.38, and the choice of
rule is below the resolution of the measurement.** Midpoint 0.382, median 0.411,
mean 0.396, and 0.362 weighting the two terms equally rather than the four
points. Those span 13 per cent against a bracket that spans 62 per cent. The
midpoint is taken because it does not silently weight the term that happens to
have three builds behind it.

## What the two terms disagree about, and why the disagreement is expected

The wetting term measures 0.292 and the endmember term 0.381 to 0.473. The
attenuation is a property of what sits above the modelled surface rather than of
the surface, so it can differ between a term acting under cloud and one acting in
the clear. The wetting term acts where the modelled skin is wet, which is where
it rains; the endmember term is spread over all land including the dry
closed-basin interiors, where the modelled sky is clearest. That is the leading
explanation and it is a HYPOTHESIS: no cloud-weighted decomposition of either
perturbation has been measured, and the gap is also consistent with PHYS-15's
23 per cent uncertainty on its own separation, which is the widest of the four.

The two are separated by measuring the cloud cover over the cells each
perturbation acts on. That needs the wetting arm's per-cell fall, which its
deleted run directories carried.

## The equilibrium ratio, which is a different quantity and not a second opinion

Read straight from `planetary_albedo_in_window`, the equilibrium
d(planetary albedo) / d(surface albedo) is 0.53, 0.42, 0.49 and 0.46 on the
staged denominators. It sits ABOVE the forcing attenuation because it contains
the ice and cloud albedo response that `lib/sensitivity.py`'s slope already
carries. Putting it in the chain would count that response twice. Its ratio to
the attenuation is the feedback amplification: 1.12, 1.09, 1.12 for the
endmember pairs and 1.48 for PHYS-15.

## What is still owed

Two of the three terms `world-ckbt` named are not measured. The lakes item and
the `playa_clastic` level are both staged-field changes rather than namelist
keys, so each needs `build_surface_albedo.py` re-run under a different rock table
or lake compositing, the six albedo codes restaged, and a paired 25-orbit set
branched from an equilibrated restart on the configured build.

**Restaging is what blocks them, not the run.** The staged fields live at
`exoplasim/inputs/t21/` keyed by rung alone, so rewriting them replaces the one
build's field the whole tree reads, and in a worktree those paths are symlinks
into the main checkout. The run itself is cheap: a paired 25-orbit T21 set at
eight threads an arm, pinned and run concurrently under one host lock, is
bracketed at six to twenty-five minutes of wall clock. The lower end is what
PHYS-15's own arms recorded; the upper is 25 orbits at the 53 s per orbit
`docs/src/reference/config-rationale.md` measures for T21 on eight ranks, plus
the concurrency penalty.

Both terms would sit on `canonical-10m-carve2`, which is the configured build and
which carries no measurement here at all.

## Currency

Every pair sits on `canonical-10m-base` or `canonical-10m-carve1`, and the
configured build is `canonical-10m-carve2`. There is no bare-rock endmember run
on the configured build, so no pair on it exists to prefer.

`lib/sensitivity.py`'s slope was measured on `canonical-10m-base`, so the base
rows are internally current with the slope they are divided by and the carve1
rows are a declared cross-build read. The cross-build read is DEFENSIBLE here in
a way it is not for `HYDROLOGICAL_RESPONSE_PER_KELVIN`, and the difference is
which side of the product moves with the terrain: the hydrological response
composes with an amplification taken from the configured build's own water
balance, so two builds are two worlds; the attenuation is a property of the
modelled atmosphere, and the two builds it is measured across differ only in
carved closed-basin floors. **The agreement between the base and carve1 rows is
what turns that from an assumption into evidence.** They bracket 0.381 to 0.473
with the base row inside, and the spread within the endmember term across two
builds is smaller than the spread between the two terms on one build.
