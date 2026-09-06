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
