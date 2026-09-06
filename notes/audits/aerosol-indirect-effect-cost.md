# What the aerosol indirect effect would cost, and what it is worth

Worldbuilding frame: this prices a process that is ABSENT from the Vesper
project's climate model. Every field named is a model array or an offline
product of this project; nothing here is about the real world.

CLIM-64 asks for a price. The chain is aerosol -> condensation nuclei -> droplet
number -> cloud effective radius -> cloud albedo, and in this model nothing
connects any link of it: the aerosol products are offline and mass-based, the
clouds are diagnostic, and the cloud optical depth is a fitted function of liquid
water path alone.

**The term is absent, not wrong, and that decides how the price is read.** A
process belongs in the model because it exists, not because including it improves
a comparison. So the price below is one input to a decision and the other input
is the size of the term; "it would move the answer away from something familiar"
is not on the list.

## The criteria, fixed before the size of the term was computed

Declared 2026-09-05, before any forcing was evaluated.

**The threshold is the project's own and it is not invented here.**
`notes/dust.md` fixes a reopening test at a global-mean forcing of 1.5 W/m2, and
`aeolian/config/dust.yaml` carries it as a value nothing may be tuned toward or
away from. On `lib/sensitivity.py`'s conversion at the baseline planetary albedo,
1.5 W/m2 is 1.03 K.

- **WORTH BUYING if the term's global-mean shortwave forcing crosses 1.5 W/m2**
  anywhere inside the bracket that the aerosol number and the activation route
  can be established to, since a term that size is already the scale at which
  this project reopens an in-model question.
- **NOT WORTH BUYING NOW if the whole bracket sits below 1.5 W/m2**, in which
  case the finding is that it is absent and small, recorded with its bracket, and
  a row says what would move it.
- **REPORTED AS UNRESOLVED if the bracket spans the threshold**, which is a
  statement about the aerosol number and not about the parameterisation.

The size is computed as a SENSITIVITY rather than as a forcing, because a forcing
needs a droplet number this project cannot yet supply: the quantity evaluated is
what the model's own cloud albedo does per factor of two in droplet number,
through the model's own optical depth and two-stream, and the bracket is then the
range of droplet-number factors the three aerosol products can support.

**The instrument against the effect.** The conversion from a top-of-atmosphere
flux to a temperature is `lib/sensitivity.py` and nothing else; the slope it
carries is local and measured across two converged points, so a forcing quoted
here is converted once and never re-derived. The model's cloud optical depth is a
fit whose stated range bottoms out at a liquid water path of 10 g/m2, and the
sensitivity is evaluated on the baseline climatology's own cloud fields rather
than on a chosen column, so the answer is weighted by where the model actually
puts cloud.

**Rule 7.** Nothing is commissioned. Adding this term costs what the next cycle
runs, and no existing run is charged against it.

## What the term is worth: measured 2026-09-05

`analysis/aerosol_indirect_sensitivity.py` -> `analysis/aerosol_indirect_sensitivity.json`,
on `run_67323a923013`'s baseline climatology on build `canonical-10m-carve2`. The
computation is offline arithmetic over a climatology and takes seconds, so no
host load is quoted: nothing here is a timing.

**A doubling of droplet number is worth 3.6 to 4.6 W/m2 in the global mean, 2.5
to 3.2 K.** A quadrupling is 7.0 to 9.1 W/m2 and 4.8 to 6.2 K. The bracket is
the product of two arms, both reported because neither is a refinement of the
other:

| vertical shape of the cover | surface underneath | x2 in N | x4 in N |
| --- | --- | ---: | ---: |
| the liquid water profile | no | 4.27 W/m2 | 8.21 W/m2 |
| the liquid water profile | yes | 3.62 W/m2 | 7.01 W/m2 |
| uniform over levels | no | 4.64 W/m2 | 9.10 W/m2 |
| uniform over levels | yes | 3.82 W/m2 | 7.54 W/m2 |

The surface arm is the one that matters for the size: a cloud over a bright
surface adds less at the top of the atmosphere than the same cloud over a dark
one, because the surface was already reflecting what the cloud now reflects. It
takes about 15 per cent off. The vertical shape is worth about 9 per cent and is
the assumption, since the model emits only the column total `clt` and its
per-level cover cannot be rebuilt from a climatology.

**The verdict against the criterion fixed above: WORTH BUYING, by a factor of
2.4 to 3.1.** Every corner of the bracket clears 1.5 W/m2, so the answer does not
turn on either arm.

**And the inversion is the number to carry away. A change of 25 to 33 per cent
in droplet number already reaches the reopening threshold.** That is what makes
this decidable without a droplet number: the threshold is not at some large
aerosol perturbation this world might or might not produce, it is at a third of
one doubling. Any indirect effect this model could have that is not
approximately zero is material.

## Why the instrument can carry that, and where it cannot

**The chain is the model's own at every step.** The cloud water is `rainmod`'s
CCM3 diagnostic reconstructed from the climatology's `prw`, `ta` and `ps`; the
optical depth is `radmod`'s two Stephens (1978) fits; the reflectance is the
Stephens, Ackerman and Smith (1984) tables that `world-f9ig` put into `swr`, read
back out of `exoplasim/scripts/stephens_tables_vs_fits.py`, which already checks
that transcription against the model source entry for entry. Nothing here is a
textbook parameterisation standing in for the model.

**The cloud fraction is an identity rather than a reconstruction, and the
reconstruction is what refused it.** Distributing the emitted `clt` over levels
as `dcc_k = 1 - (1-clt)^s_k` reproduces the model's own random-overlap total
exactly for any shape `s`, and it does: the maximum error over the whole field
is 3.0e-08, which is `clt`'s own float32 spacing. The alternative was to rebuild
the per-level cover from relative humidity, and that arm was run and REFUSED --
the median ratio of the reconstructed total to the emitted `clt` came out at 0.0
over 24,575 cells carrying cloud, because `mkclouds`'s convective branch is keyed
on `icclev`, which is not an output. The stratiform half alone is not this
model's cloud.

**The closed form and the tables agree to the term the closed form omits.**
`radmod`'s band-1 branch is `A = x/(1+x)`, whose derivative in `ln(tau)` is
`A(1-A)` -- Twomey's sensitivity exactly. Over 114,371 layers with optical depth
above 1 and a zenith cosine above 0.05, the tabled difference and
`A(1-A) ln2 / 3` differ by a median of 22 per cent and 39 per cent at the ninth
decile. That gap is `d ln beta / d ln tau`: the backscatter fraction is itself
tabulated against optical depth and the closed form holds it fixed. The tabled
path is the one quoted, because it is the one the model runs.

**Three things this measurement does not include, and their signs.** It stops at
the SHORTWAVE, so the longwave cloud response -- opposite in sign and smaller,
since `radmod`'s longwave cloud absorption is grey in liquid water path -- is
absent. It uses monthly-mean zenith cosines in a function that is not linear in
them. And it takes no account of what the atmosphere above the cloud absorbs. The
first is the only one that could move the answer toward the threshold, and it
would have to remove three quarters of the term to reach it.

**What it deliberately does not do is compute a droplet number.** The factor of
two is a unit of the derivative and not a prediction. Whether this world's
aerosol can move droplet number by a third is the aerosol side's question, and
the answer to it does not exist yet.
