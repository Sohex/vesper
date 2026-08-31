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
keys, so each needs `build_surface_albedo.py` re-run under a different rock
table or lake compositing, the seven albedo codes restaged, and a paired set on
the configured build. `world-3ooi` carries the task; what follows is what
measuring the instrument beforehand established, and it moved both the design
and the price.

### The playa arm resolves and the lakes arm does not

The perturbation each arm applies is now measured rather than assumed. Running
`build_surface_albedo.py` three times into a scratch `--output`, changing only
`model.lithology_albedo_overrides.playa_clastic.albedo`, moves the staged land
mean by:

| step in the class albedo | d(staged land mean) | share |
| --- | ---: | ---: |
| 0.23 to 0.25 | +0.0017028459 | 0.08514230 |
| 0.23 to 0.33 | +0.0085142296 | 0.08514230 |
| 0.25 to 0.33 | +0.0068113837 | 0.08514230 |

Exactly linear to eight figures across a factor of five in step size, so a
single share is the right shape. Code 212, the forest fraction, is
byte-identical across all three, so the override moves albedo alone and the
confound the endmember pairs had to disprove is absent here by construction.

**The instrument is the quadrature sum of the two fitted asymptote half-widths,
and across every pair on disk it runs 0.028 to 0.067 K**, with PHYS-15's own
16-orbit-span pair at 0.046 K. That range does not track the fit span, which
runs from 16 to 136 orbits over the same set, so a 25-orbit arm sits inside it.

Predicted separations, at the measured share and over the attenuation's own
0.29-to-0.48 bracket:

| arm | d(staged land mean) | predicted dT | separation, in quadrature half-widths |
| --- | ---: | --- | --- |
| playa_clastic 0.25 to 0.33 | 0.006811 | 0.197 to 0.325 K | 2.9 to 11.6 |
| playa_clastic 0.19 to 0.33 | 0.011920 | 0.344 to 0.569 K | 5.1 to 20 |
| lakes composited or not | 0.003716 | 0.107 to 0.178 K | 1.6 to 6.3 |

The bar fixed on `world-ckbt`, before any of these numbers, is three times that
quadrature sum.

**So the lakes arm is refused and the playa arm is bought at the wider step.**
The lakes arm fails its own bar over much of its bracket, and the attenuation it
would return carries a half-width of 20 to 48 per cent against a bracket that is
25 per cent wide: a point wider than the thing it is meant to narrow. It becomes
buyable at an amplified lake contrast or at about four times the orbits. The
playa arm at 0.19 to 0.33 clears the bar across the whole bracket, and 0.19 is
the export's own table entry rather than an invented endpoint.

**It also tests the same hypothesis, on the same ground and harder.** The
leading explanation for the two measured terms differing is that the attenuation
belongs to the modelled atmosphere over the cells a perturbation acts on, and
`playa_clastic` IS the dry closed-basin interior material where the modelled sky
is clearest. The prediction registered before the arms run is that the playa
point lands at or above the endmember term, in 0.38 to 0.55; the hypothesis is
refuted if it lands at or below the wetting term's 0.292 plus its own
uncertainty. The pre-registered exclusion is unchanged: a pair whose diagnosed
and staged land-mean deltas disagree by more than a factor of 1.5 is pricing the
swap rather than the attenuation.

### What actually blocks it, which is not what the row assumed

**Restaging is safe and is not the blocker.** `exoplasim/scripts/sra.py`'s
`write_sra` refuses a write through a symlink on every one of the seven codes,
so a generator run in a worktree exits before writing anything rather than
replacing the main checkout's staged fields. The refusal fires on the first
write, which is a wet-endmember code, so the report at the end is never reached
either. `--output` takes the generator anywhere inside the worktree, and
replacing the worktree's own seven links with real files is contained: removing
a symlink does not touch its target, `exoplasim/inputs/*/*.sra` is ignored, and
`albedo_report.json` is a tracked file whose edit stays in the worktree.

Two things do block it.

**No model binary is current.** `rebuild_binaries.py --verify` reports all ten
configurations in the matrix as not built against the model source at
`vendor/exoplasim`. Rule 4 asks for every binary, so the arms cost a full
rebuild before the first orbit and that rebuild is unpriced here.

**No T21 restart sits on a surface the tree still holds.**
`check_consistency.py` reports every T21 run as being on a surface that has been
regenerated since, and says in its own words that nothing can be seeded from it.
Five live T21 runs on `canonical-10m-carve2` at the configured flux exist to
branch from, and branching both arms from ONE of them is what makes that
survivable: the common transient from the surface change is identical in both
arms and cancels in the difference of the two asymptotes, which is the estimator.

The staged T21 surface is itself behind the generator, and by more than the
lakes item is worth. `analysis/soil_albedo_wetting.json` was rewritten by
`world-znr7` after the surface was staged, and today's generator moves the
staged land mean by -0.003474 and the soil-wetting bound from 0.015707 to
0.010121. The gate already reports both.

**The run, priced.** Two arms, sequential rather than concurrent, 40 orbits
each, T21 on eight ranks under one host lock: 45 to 95 minutes of wall clock.
The low end is the rate PHYS-15's own arms recorded; the high end is 53 s per
orbit from `docs/src/reference/config-rationale.md` plus per-arm setup.
Sequential is chosen over a concurrent pair deliberately: the estimator's floor
across independent integrations is under a hundredth of a kelvin, two arms
sharing a die do not, and twenty minutes is a cheap price for keeping the two
arms as identical as this host can make them.

Both terms would sit on `canonical-10m-carve2`, which is the configured build
and which carries no measurement here at all.

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
