# Pre-registered aerosol-number and cloud structural sensitivity

Worldbuilding. Vesper is an invented planet and everything here is about the
simulation of it: a climate model's radiation and cloud fractions, an offline
aerosol chain, and the modelled organic and marine particles that would feed
them. Nothing below is an observation.

Registered 2026-08-24, before any arm exists and before any number from any of
them has been seen. Its subject is a MODEL-FORM question -- whether the
simulation gets an aerosol-to-droplet pathway at all -- and a model-form
question settled after the results are in is not settled, it is chosen. The
decision rule, the metrics and the threshold below are fixed as of this
document. If an arm's result makes a rule look badly posed, the rule is
re-registered in this file with the reason and the old one kept visible in the
history; it is not quietly widened.

The finding this rests on is
[`bvoc-soa-atmospheric-coupling-audit.md`](bvoc-soa-atmospheric-coupling-audit.md)
finding 9 and `notes/audits/ocean-and-marine-biosphere.md` finding 5c. The
activation contract that carries the declaration fields is
[`bvoc-activation-contract.md`](bvoc-activation-contract.md).

## Why this is registered rather than decided

Mass and optical depth do not determine particle number. The climate model has
no prognostic aerosol activation and no droplet-number response: its cloud
optical properties do not read any aerosol field, and its shortwave aerosol
effect acts in the cloud-free fraction of each layer only. So the pathway is
absent rather than approximated, and its absence cannot be presented as a
direct-optics approximation.

The Earth literature says the omission can be large: a chamber study
demonstrating nucleation from pure biogenic oxidation products including an
ion-induced path with no sulfuric acid, a microphysics study putting that
mechanism into a global model and finding a large change in preindustrial
aerosol-cloud forcing, and an uncertainty study attributing a substantial share
of industrial-period indirect forcing variance to natural emissions. **Those
numbers are not this world's priors and none of them is a target.** What they
establish is only that the question cannot be dismissed, which is why it gets a
registered decision rule instead of a judgement call at the end.

## The three arms, fixed

Exactly three, shared across every aerosol source. There is no fourth, and in
particular no arm in which a cloud effect is represented by scaling a column
optical depth.

**Arm 1, no cloud effect.** The registered default and the null. Direct optics
only, in the radiation's clear-sky fraction, from whatever prescribed organic
aerosol climatology BVOC-8 produces. Costs nothing beyond the direct-optics
work and is what the project keeps if the rule below does not fire.

**Arm 2, an offline particle-number bound.** Not a model of clouds: a bound on
cloud condensation nuclei computed offline from the same archived meteorology
the chemistry-transport operator runs on. It carries, explicitly and separately:

- the terrestrial precursor source, from the vegetation model's speciated
  volatile flux under the activation contract;
- the marine precursor source, which is OCN-8's to specify and nobody else's;
- the sea-salt seed population, which is CLIM-40's;
- the highly-oxidized-molecule and sulfur nucleation yields as separate
  branches;
- ion production, particle size, hygroscopicity, the pre-existing seed aerosol
  surface, and the cloud supersaturation the activation would occur at.

Those six are the registered contents of the bound, and they are registered
because each is a quantity an optical-depth multiplier silently averages over.
An arm 2 that cannot report all six is not arm 2.

**Arm 3, a minimal activation and cloud-albedo port.** A droplet-number
response in the climate model. Built ONLY on the verdict below, and it is one
shared port: no source-specific cloud module, no separate marine and terrestrial
pathways. If it is built, every source feeds the same activation.

## The decision rule

Arm 3 is built if and only if

    (T1 or T2) and T3

evaluated between arm 1 and arm 2 on the same terrain, the same rung and the
same spectrum.

**T1, the terrain test.** Arm 1 and arm 2 disagree about the carve verdict on a
basin that BOTH bounding climates agree about. The qualification is the whole
of the test: loop A already exits by carving the intersection of the verdicts
taken at the two bounding climates, so a basin the bounding climates disagree
about is already handled as bracketed, and an arm that only moves such a basin
has changed nothing the project does. A basin the bracket agrees on and the
arms do not is outside what the existing construction covers. The carve list
leaves the project and the terrain cannot be un-carved, so the count that fires
this test is one.

**T2, the vegetation test.** Arm 1 and arm 2 disagree about the modelled
vegetation verdict that loop C exits on, over more than 2 percent of simulated
land cells. Two percent is registered here as a fraction rather than derived,
and it is registered because loop C's exit is a re-taken verdict rather than a
tolerance: below a fraction this size the verdict is being moved by cells that
sit on a classification boundary, and the same cells move under a re-run of
either arm.

**T3, the instrument test.** The two arms are distinguishable at all.
`exoplasim/scripts/compare_equilibria.py` reports them as NOT at the same
equilibrium, on the criterion that script already declares -- two sigma on a
difference of two means, against each run's own year-to-year scatter -- for at
least one of its registered metrics. Reflected top-of-atmosphere shortwave is
the metric where a cloud change appears, and it is already in that script's
list.

T3 is a necessary condition and never a sufficient one. Its job is the
instrument check: an arm difference smaller than the model's own scatter is
noise however tidy the verdict table looks.

### The three outcomes, all registered

- **T3 fails.** The arms are indistinguishable. Arm 3 is not built and the
  verdict is `not_material`. If T1 or T2 fired while T3 failed, the verdict is
  `below_resolution` and the response is a longer production window on both
  arms, not a new model: a verdict difference under the noise floor is a
  statement about the run length.
- **T3 holds and neither T1 nor T2 fires.** The arms differ measurably and no
  decision the project makes moves. Verdict `not_material`; arm 3 is not built.
  The spread is carried into `analysis/error_budget.json` as an item like any
  other, which is where a real but non-deciding uncertainty belongs.
- **T3 holds and T1 or T2 fires.** Verdict `material`. Arm 3 is authorised,
  once, as one shared port.

The verdict is written to `clouds.materiality_verdict` in
`biosphere/config/bvoc.yaml`, with the scoring artifact that produced it. The
gate refuses arm 3 without it (`BVOC-CLOUD-VERDICT-MISSING`), and refuses any
arm carried on the prescribed aerosol optical-depth field
(`BVOC-CLOUD-OPTICAL-PROXY`).

## What is not registered here

- The direct optical effect of the organic aerosol. That is BVOC-8 and it is a
  separate bracket with a separate result; this document is only about particle
  number and clouds. The two are kept apart because the audit's finding is that
  they have been conflated, and because the direct effect is mostly scattering
  while the number effect need not share its sign.
- The marine and sea-salt SOURCES. OCN-8 and CLIM-40 own those and feed arm 2;
  neither builds a cloud pathway of its own. That is the point of registering
  one shared sensitivity rather than three.
- Any Earth forcing magnitude as a prior or a target for any arm.

## What voids this registration

The registration is against a fixed comparison. It is void, and re-registered
before the arms are re-run, if any of these changes between registration and
scoring: the carve criterion, the stellar spectrum, the resolution rung, the
terrain build, or the source bracket the arms are driven by. A registration
scored across one of those changes is measuring the change.

## Where the numbers go when they exist

Nowhere in this file. Arm results are dated records with a "measured on" under
`notes/`; the current values live in `world_state.json`; the spread's price in
the project's own currency goes through `scripts/error_budget.py`, which is
what makes it rankable against everything else the project is unsure about.
This document holds the rule and the threshold and keeps holding them
unchanged.
