# The land column property contract

**Recorded:** 2026-08-24
**Enforced by:** `pedology/scripts/land_column_properties.py`
**Declared in:** `pedology/config/land_column_properties.yaml`
**Finding it rests on:** `biosphere/notes/soil-land-surface-hydraulic-consistency-audit.md`, finding 1
**Books into:** `hydrography/config/land_water_ledger.yaml`

Worldbuilding. Vesper is an invented planet. Everything below is about the
simulation of it: a declared column of soil, the states and flow properties it
is described by, and what each consumer currently derives for itself instead.
No model was run to produce any of it.

This is a CONTRACT, not a finding. What is true about the repository's soil
handling is in the audit; what to do about it is in `bd`. This document says
what the description IS, where gravity enters it, and what fails when a
consumer disagrees with it.

## What is being fixed, and what is deliberately not

One pedology artifact becomes two incompatible vadose-zone soils. Pedology
writes texture, organic carbon, bulk density, regolith depth, a plant-available
water capacity, weathered-bedrock fraction, andic fraction and phosphate
fixation. ExoPlaSim's surface builder reads exactly one of those columns,
converts millimetres to metres, and installs it as a bucket depth. LPJ-GUESS
reads none of it: it re-derives saturation, matric potential, wilting point and
field capacity from the Cosby texture regressions, mixes in an ideal organic
soil down a fixed Earth profile, and then scales the result by regolith depth
and weathered-bedrock fraction.

This contract does not yet REPLACE either path, and doing so is not the next
step. The audit's own recommended order puts the comparison first for a reason:
removing a consumer-side pedotransfer without knowing what it currently
produces would change every downstream result and leave nothing to attribute
the change to. So the contract is the description, and the enforcing script is
the comparison, and neither runs a model.

## What the column is

**Two geometries, and they are different cuts.** The physical column is fifteen
uniform 100 mm layers to 1.5 m: where water is and where it moves. The rootable
column stops at 1.0 m: where uptake may draw from. The physical base is where
`soil_liquid` stops and the aquifer begins, which is the contact the land water
ledger's `drainage` and `capillary_rise` terms cross.

The 1.5 m is not invented. It is LPJ-GUESS's own physical profile, and it is
the only column in this pipeline that HAS a depth. ExoPlaSim's `dwmax` is a
capacity in metres OF WATER and carries no thickness at all; pedology's regolith
depth scales LPJ's per-layer capacities without moving a layer boundary or a
root, so a shallow-regolith cell has the same rooting depth as a deep one with
less water in it.

**Four materials, with their vertical rules stated.** Mineral texture, uniform
down the column because pedology emits one texture per cell and both consumers
apply it to every layer. Organic, on a fixed Earth profile under
Lawrence and Slater in LPJ-GUESS and as an additive bulk term with no profile at
all in pedology -- two rules for one soil, currently agreeing because the
checked-in soil map carries zero soil carbon everywhere, and not once the
biosphere loop has turned. Andic, which pedology gives a hydraulic effect and
LPJ-GUESS does not have at all, so the two derivations diverge most exactly
where the soil is most volcanic. Weathered bedrock below the regolith contact,
capped so sub-bedrock material can at most match the soil above it, which is a
limitation of the model recorded so it is visible in results.

**Four states, and the ordering between them is checked.** Residual, wilting
point, field capacity, saturation. Available capacity is the difference of two
of them and is the only one either consumer receives -- which is why nothing
downstream of `dwmax` can recover saturation, and why water between field
capacity and saturation has nowhere to sit on either side of the boundary.

**One named retention closure.** Clapp-Hornberger, of the Brooks-Corey family,
on Cosby Table 4 parameters. The project already sits entirely in that family
and had not recorded that a family was chosen. The alternative is van Genuchten
and the two differ most in their tails, which matters here specifically because
wetting and drying through near-saturation is a common state on this land rather
than an edge case: a large share of it is playa clastic, which LITH-27
establishes is derived from endorheism and therefore moves with the carve list.
Naming the closure is not an argument for switching it. It is what makes a later
comparison have something to vary.

## Where gravity enters, and where it does not

This is the one place a super-Earth changes a soil property with no free
parameter and no fitted term, so the contract declares it before deriving
anything from it.

Total hydraulic potential per unit weight is `H = psi / (rho_w g) + z`. The
three states respond differently:

- **Saturation** is the geometry of the pore space. Gravity does not enter. The
  same number on any planet.
- **The wilting point** is a PLANT pressure: the suction a root can generate
  before it can no longer extract. Gravity does not enter that either, so the
  Earth-derived value crosses over as a pressure unchanged.
- **Field capacity is an operational drainage equilibrium**: the water a column
  retains against gravity once free drainage has stopped, over a stated drainage
  length. The length is a length; the matric pressure that balances it is
  `rho_w * g * L`, and THAT scales with gravity.

The drainage length is declared as 1 m, and the reason is not roundness.
Cosby's own field-capacity suction, as `soilinput.cpp` inverts at it, is 100 cm
of water, which is exactly `rho_w * g_earth * 1 m`. That coincidence is what
licenses stating field capacity as a drainage length at all, and the contract
check fails if the declared length stops reproducing it -- because at that point
the derivation would have silently changed meaning.

**The consequence is one-directional.** At this planet's gravity the same soil
drains to a drier state than it would on Earth. Field capacity falls, the
wilting point does not move, and plant-available water therefore falls by more
in relative terms than field capacity does. The size is set entirely by the
retention exponent, so it is small and tightly bracketed, and the report carries
it per cell.

It is reported and NOT applied. Two reasons, and both are in the report: it
holds the Clapp-Hornberger exponent fixed while gravity moves, which is an upper
bound on the effect rather than a correction factor; and the within-texture-class
variance Cosby reports is larger than the shift, so applying it as a point
correction on top of an unquantified spread would be a false precision. What
makes it worth carrying anyway is that it moves every cell the same way, so
unlike a scatter it does not average out over the map.

Neither consumer applies it today. Pedology's declared endmember capacities are
Earth field observations, and LPJ-GUESS's Cosby inversion evaluates at the Earth
suction.

## Uncertainty is CORRELATED, and that is a declaration

Three named cases -- low, central, high -- applied coherently across every layer
and every property at once. Never sampled per layer.

The reason is physical. A soil that holds more water at one depth holds more at
the next, because the same texture, the same organic content and the same
endmember uncertainty produced both. Drawing each of fifteen layers
independently averages the spread away by roughly the square root of fifteen and
reports a precision the description does not have. The contract check refuses a
declaration that says otherwise.

The declared endmember brackets come from `pedology/config/pedogenesis.yaml` and
carry their own argument: they are DECLARED values with Saxton and Rawls (2006)
as the check on their ordering and magnitude rather than their origin, because
that regression parameterises a continuous mixture and excluded samples above 60
percent clay, so a pure endmember is outside the data it was fitted on.

**The reported spread is a lower bound and is labelled as one.** Cosby et al.
report the residual variance within texture classes and it is large -- larger
than everything else in this contract's uncertainty put together. No artifact
here carries it. Until one does, the bracket the report quotes is what the
declared endmembers alone are worth.

## The check that would fail

**The declaration check**, on the contract itself: the physical layers must
reach the declared column base, the rootable base must lie above it, the vadose
base and the aquifer contact must be one number, every state named in the
ordering must exist and must say whether it depends on gravity, the wilting
point must carry a suction's sign, the drainage length must reproduce Cosby's
own field-capacity suction under Earth gravity, the closure must name a family
that has an implementation, the uncertainty cases must be declared correlated,
and the aquifer property must be declared separately scaled.

That check found a real omission in this contract on its first run: `available`
was declared silent on gravity while both states it is built from were not.

**Ten declaration mutations**, because a check returning nothing proves nothing:
a contract that is correct and a check that cannot fail look identical from
there. Each mutation breaks the declaration in one named way -- a layer count
that no longer reaches the base, a root below the column, a contact depth stated
twice and differing, a state in the ordering that does not exist, a state silent
on gravity, a wilting point with the wrong sign, a drainage length that is not
Cosby's, a closure family with no implementation, uncertainty drawn per layer,
and an aquifer property shared with the vadose zone -- and the check is required
to catch every one.

**The Cosby inversion's air-entry clamp.** Cosby's equation holds only below the
air-entry value; at and above it the pore space is full and the water content
is saturation. Unclamped, the inversion returns a field capacity above
saturation wherever the air-entry suction falls below the field-capacity
suction, and substituting the regressions and `silt = 1 - sand - clay` gives
that region as a straight line in the texture simplex,
`1.58 * sand + 0.63 * clay < 0.17`, which needs no search. WORLD-NGA10 clamped
both heads at the air-entry value in `get_mineral`, so the CODE is what keeps
this right and the map has never had to. The check transcribes the clamp,
because a transcription that kept the unclamped form would be comparing against
a model that no longer exists, and it still measures the soil map against the
line as a MARGIN and not as a count -- a count of zero says nothing about how
close the nearest cell is, and a map that moved inside would have its field
capacity set by the clamp rather than by its texture.

## What one central case would move

The removal this contract was built to make attributable is WORLD-OF6N, and it
needs a decision that is about the world rather than about the code: what the
central case IS. The report costs the candidates so that decision is a
one-liner, and it runs no model.

Two of the three are readings of this contract as it already stands. The
closure is NAMED -- Clapp-Hornberger on Cosby Table 4 -- and field capacity is
declared as a drainage equilibrium at `rho_w * g * L`, so taken literally the
states follow from texture and the only open part is which gravity the suction
is evaluated at. Adopting the second of those is the contract read literally
and is the only candidate that moves BOTH consumers. The third would replace
the named closure with pedology's endmember mixture, which is not a retention
curve at all: it carries no saturation, no air entry and no exponent, so
saturation and matric potential would have to come from somewhere else. That
makes it a change to the declaration rather than a choice between consumers,
which is what the report says and why it is not offered as a symmetric option.

The allophane term is not what separates the two derivations, which is worth
saying because the contract's own materials block reads as though it is. It
contributes a median of nothing and a p90 of a few per cent of pedology's
capacity, against a median absolute gap of tens of millimetres, and it is zero
on more than half the cells. What separates them is that the declared endmember
capacities and the Cosby retention inversion are different numbers for the same
texture.

Two things the report cannot supply, and says so rather than estimating. The
change in land runoff ratio and in E over R needs a paired baseline, and
`probe_runoff_response.py` needs a climatology to back-solve potential
evaporation from; `config/planet.yaml` declares `baseline_climatology` null. And
nothing on the LPJ-GUESS side is execution-verified: it does not build on this
tree until a baseline climatology exists, so every LPJ number here is a
transcription of the vendored source evaluated offline.

## What is still undeclared, and who owns it

Eight properties carry the sentinel. Each names its issue, and `--strict` is the
arm that refuses.

| what is missing | owner |
| --- | --- |
| the residual water content, which both active evaporation limiters pin at zero without saying so | LSHY-3 |
| saturated hydraulic conductivity; no artifact here carries one, and LPJ's `perc_base` is a dimensionless drainage exponent that cannot be converted into one | LSHY-3 |
| unsaturated conductivity: the family is declared, the number is not | LSHY-3 |
| infiltration capacity, absent from both columns, which is why saturation excess is the only runoff mechanism | LSHY-3 |
| the frozen-pore impedance EXPONENT; the form is declared and is used by nothing, because the scheme that would need it has no conductivity to reduce | LSHY-5 |
| composition-dependent thermal properties, so that one material property is not split into an organic profile for water and a constant for heat. The PHASE half is now declared: freeze and thaw exchange the latent heat of fusion with the soil temperature layer that contains the water layer's midpoint, and the sensible heat capacity is still a constant | LSHY-5 |
| the organic fraction and one vertical rule for it | SDEC-10 |
| Cosby's within-texture-class variance, the dominant uncertainty term | LSHY-1 |

## What this contract does not own

Aquifer transmissivity. A vadose-zone unsaturated conductivity is not an aquifer
permeability: they live at different material and spatial scales, they are
measured on different things, and copying one into the other is a category error
with plausible units. The groundwater solver owns transmissivity as `K` times
saturated thickness, scales it separately, and is joined to this column at the
declared contact through the ledger's `drainage` and `capillary_rise` terms and
through nothing else. The contract check refuses a declaration that says the two
are scaled together.
