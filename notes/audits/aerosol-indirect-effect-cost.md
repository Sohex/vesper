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
this decidable without a droplet number: the threshold sits at a third of one
doubling rather than at some large aerosol perturbation this world might or might
not produce. Below it the term still is not zero -- a tenth of a doubling is
about 0.4 W/m2, or 0.3 K -- so what the bracket has to establish is not whether
the effect exists but whether this world's aerosol moves droplet number by more
than a quarter, which is a question about the aerosol products and not about the
parameterisation.

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

## What it would cost to build, and the chain has four breaks rather than one

Read 2026-09-05 from `vendor/exoplasim`, `aeolian/`, `references/pyrcel/` and
`references/cloudmicrophysics/`. The row names the updraft as the blocker. It is
one of four, it is correctly identified, and it is stated slightly wrong.

**1. Hygroscopicity, and it is the only genuinely missing aerosol input.**
`notes/external-model-survey.md` section 50a establishes that a condensation
solver reads four numbers per mode -- number concentration, median dry radius,
geometric standard deviation, and kappa -- and that composition collapses into
the last before the run starts. Three of the four already exist in this
project's products in some form. The dust optics carry a number median radius
and a geometric standard deviation; the sea salt source function IS a number
flux in three declared lognormal modes before it is converted to mass; the
volcanic sulfate carries an OPAC size. **Kappa exists nowhere**, and grep for it
across `aeolian/` returns von Karman's constant and two hygroscopic-GROWTH
headers. The repair is three declared literals from a paper that is already on
disk: `references/pdf/petters2007-kappa-hygroscopicity.pdf`, currently recorded
as `held` and declined for the OPTICS use because pairing a different growth
curve with OPAC's refractive indices would break a consistency the optics script
checks. That reason does not reach the CCN-activity half of the same paper. This
is the cheapest link in the whole chain.

**2. The parameterisation, and it is already acquired.** Section 50b's
instruction was to price acquiring an ARG implementation rather than PySDM's
tree. Two are now on disk and neither is recorded.
`references/pyrcel/pyrcel/activation/_arg2000.py` is 112 lines of which about 60
are arithmetic, BSD 3-clause, in this project's own language, and carries four
schemes side by side in its legacy module so the row gets a second opinion free.
`references/cloudmicrophysics/src/AerosolActivation.jl` is 435 lines, Apache-2.0,
Julia, and is the better-retargeted of the two: gravity, the water and ice
densities, the surface tension and the six ARG fit coefficients are all fields of
one runtime parameter struct rather than module globals, so a Vesper instance is
a constructor call. It also carries a Korolev and Mazin (2003) correction for
pre-existing liquid and ice that a diagnostic-cloud host actually needs, since a
model that re-diagnoses `dcc` and `dql` every step would otherwise activate a
full aerosol population into an already-cloudy cell every timestep. **Gravity
enters both at FIRST POWER and in one place**, inside the coefficient that
multiplies the updraft, so on this world gravity and the updraft are degenerate
in the ARG algebra: a 1.31x gravity is exactly a 1.31x updraft. That is worth
knowing because it puts the gravity retarget two orders of magnitude below the
uncertainty in the next item.

**3. The updraft, which is the blocker, and the row states it one step too
early.** It is not that nothing supplies a vertical velocity.
`plasimmod.f90` declares `dw(NHOR,NLEV)`, the grid-scale pressure velocity,
`calcgp` writes it at every level every timestep, and exactly one branch reads
it -- `rainmod`'s `clwfac` cloud suppression, which is dead by default because
`clwcrit1 < clwcrit2` makes `rainini` set `clwfac = -1` and the guard that reads
it tests `clwfac > 0`. So a cell-mean omega is there and free. What is missing is
a SUB-GRID updraft, and the model carries neither of the two quantities every
GCM-scale ARG implementation builds one from. **There is no TKE anywhere in the
model source**: `vdiff` is first-order Louis K-theory with no prognostic
turbulence energy and no boundary-layer depth. There is no convective mass flux
either: `kuo` is moisture-accession and iterates a parcel temperature, never a
parcel velocity, and `mkshallow` is implemented as an enhanced diffusivity. The
one convective velocity scale that exists is `fluxmod`'s `freeconv`, which the
source itself identifies as `w*` with the gravity and the inversion height folded
into a fitted coefficient -- ocean-only, not exposed as a velocity, and with the
inversion height absorbed so it cannot be recovered. `dust3` is `u*^3`, a surface
shear scale with no closure relating it to a cloud base. **So this is a new
boundary-layer diagnostic and not a plumbing change**, and it is the item that
carries the cost.

**4. The radiation side, which no row has named and which section 50 does not
mention.** Even given a droplet number there is nothing in `radmod` to hand it
to. The shortwave optical depth is `ztaua * log10(W)^ztaup` on the liquid water
path alone, and the fork's own comment records that below the fit's range the
effective radius is frozen at what the fit implies -- 8.2 um in band 1 and 6.7 um
in band 2. Above it the radius is whatever a 1978 Mie fit to eight standard cloud
models implies, and it is not a variable anywhere. The longwave is worse: grey
mass absorption at `acllwr`, taken as Kiehl's liquid-only coefficient
deliberately, because this model carries one condensate and no phase split, where
the same paper's ice coefficient does carry an effective radius. So both bands
and the longwave would each need re-deriving. **The one piece of good news is
structural and already banked**: `world-f9ig` moved the two-stream to table
lookups keyed on optical depth and zenith cosine, so those tables are agnostic
about how the optical depth was obtained and `tau` can be re-sourced without
touching them.

## The price, bracketed

**The parameterisation is free and the aerosol side is nearly free.** ARG is
acquired, permissively licensed, and about 60 lines of arithmetic plus a Kohler
critical-supersaturation function and five thermodynamic functions -- under 300
lines self-contained. Kappa is three literals from a paper on disk.

**The updraft closure and the radiation re-derivation are the price, and neither
can be quoted as a point.** Bracketed by what the tree already shows: a
boundary-layer velocity scale means a new diagnostic in `fluxmod` with its own
derivation and its own gravity retarget, and re-deriving `tau(W, r_e)` in two
bands plus a longwave coefficient carrying an effective radius is the same shape
of work as `world-f9ig`, which is a landed and inspectable precedent for the
radiation half.

**And a run cost that is near zero.** ARG is closed form -- no root find in the
2000 scheme -- so it is a handful of transcendentals per mode per cloudy column,
against `lwr`'s 1,045 per column per call. Nothing about this term is expensive
to integrate; what is expensive is deriving what feeds it.

**Rule 7 makes the whole of that the price.** Nothing is commissioned, so no run
is charged against adding the term.

## What this changes about where the term sits

`notes/audits/absent-and-inherited-physics.md` is the enumeration of absent
physics and does not carry this term: grep it for droplet, activation,
supersaturation or condensation nuclei and it returns nothing. Its stated scope
was the ocean, the cryosphere, the soil parameters and the couplings between
components, so the omission is a scope boundary rather than an oversight. The
term is now priced and it is the largest single absent shortwave term this
project has costed.
