# The land column property contract

**Recorded:** 2026-08-24
**Enforced by:** `pedology/scripts/land_column_properties.py`
**Declared in:** `pedology/config/land_column_properties.yaml`
**Finding it rests on:** `biosphere/notes/soil-land-surface-hydraulic-consistency-audit.md`, finding 1
**Books into:** `hydrography/config/land_water_ledger.yaml`

Worldbuilding. Vesper is an invented planet. Everything below is about the
simulation of it: a declared column of soil, the states and flow properties it
is described by, and the one derivation every consumer of it reads. No model was
run to produce any of it.

This is a CONTRACT, not a finding. What is true about the repository's soil
handling is in the audit; what to do about it is in `bd`. This document says
what the description IS, where gravity enters it, and what fails when a
consumer disagrees with it.

## What this replaced

One pedology artifact became two incompatible vadose-zone soils. Pedology
writes texture, organic carbon, bulk density, regolith depth, a plant-available
water capacity, weathered-bedrock fraction, andic fraction and phosphate
fixation. ExoPlaSim's surface builder read exactly one of those columns,
converted millimetres to metres, and installed it as a bucket depth. LPJ-GUESS
read none of it: it re-derived saturation, matric potential, wilting point and
field capacity from the Cosby texture regressions, mixed in an ideal organic
soil down a fixed Earth profile, and then scaled the result by regolith depth
and weathered-bedrock fraction. They were different numbers for the same soil.

The contract is now the only derivation. The enforcing script evaluates the
closure declared below and emits per-cell states; ExoPlaSim installs the
capacity column as `dwmax` and LPJ-GUESS takes the states and the per-layer
weathered-bedrock shares through its driver file. `soilinput.cpp`'s Cosby
inversion and `vesperinput.cpp`'s regolith rescaling are gone, and the
enforcing script still computes what each of them produced, because the cost of
a removal is only legible against the number it removed. Nothing here runs a
model.

The comparison came first for a reason, and the order was not optional:
removing a consumer-side pedotransfer without knowing what it produced would
have changed every downstream result and left nothing to attribute the change
to.

## What the column is

**Two geometries, and they are different cuts.** The physical column is fifteen
uniform 100 mm layers to 1.5 m: where water is and where it moves. The rootable
column stops at 1.0 m: where uptake may draw from. The physical base is where
`soil_liquid` stops and the aquifer begins, which is the contact the land water
ledger's `drainage` and `capillary_rise` terms cross.

The 1.5 m is not invented. It is LPJ-GUESS's own physical profile, and it is
the only column in this pipeline that HAS a depth. ExoPlaSim's `dwmax` is a
capacity in metres OF WATER and carries no thickness at all; regolith depth
scales the per-layer capacities without moving a layer boundary or a root, so a
shallow-regolith cell has the same rooting depth as a deep one with less water
in it.

**The base is a cut, not a limit reached by accident.** Regolith runs deeper
than it on part of this map, and every capacity in this contract stops there,
because nothing in the pipeline draws water from below it: the physical column
ends there, the rootable column ends higher, and the root distribution spans
exactly the physical layers. A plant-available capacity integrated past the base
would be counting water no root and no evaporating surface can reach.

**Four materials, with their vertical rules stated.** Mineral texture, uniform
down the column because the soil map carries one texture per cell and every
layer takes it. Organic, on a fixed Earth profile under Lawrence and Slater in
LPJ-GUESS and as an additive bulk term with no profile at all in pedology's own
capacity -- two rules for one soil, currently agreeing because the checked-in
soil map carries zero soil carbon everywhere, and not once the biosphere loop
has turned. Andic, which the adopted closure has no term for at all, so the
hydraulic effect of allophane is a declared absence. Weathered bedrock below the
regolith contact, capped so sub-bedrock material can at most match the soil
above it, which is a limitation of the model recorded so it is visible in
results, and emitted as one usable share per physical layer.

**Four states, and the ordering between them is checked.** Residual, wilting
point, field capacity, saturation. All four are derived here; the states file
carries three of them plus the closure's exponent, and the residual is
undeclared. Available capacity is the difference of two of them and is the only
one ExoPlaSim receives -- which is why nothing downstream of `dwmax` can recover
saturation, and why water between field capacity and saturation has nowhere to
sit on that side of the boundary. LPJ-GUESS receives all three, so the ceiling
its water-filled pore space carries is the ratio of field capacity to
saturation rather than the absence of a number.

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

It is APPLIED. The emitted states evaluate field capacity at `rho_w * g * L`
with this planet's gravity, and both consumers install them. Leaving it out was
never the neutral choice: a contract that declares field capacity as a drainage
equilibrium and then evaluates it at Earth's suction is claiming that this
world's columns drain to Earth's equilibrium, which is a claim and not an
abstention.

**The adopted value is an UPPER BOUND, and that is a declared limitation of
it.** The derivation holds the Clapp-Hornberger exponent fixed while gravity
moves, so the whole of the shift is carried by the defining pressure and none
of it by the shape of the curve. A retention curve whose exponent also responded
to gravity would fall by less. Two things make the bound worth adopting anyway:
the shift moves every cell the same way, so unlike a scatter it does not average
out over the map; and its direction is not in doubt. What it is NOT is a
precision. Cosby's within-texture-class variance is larger than the shift and no
artifact here carries it, so the adopted states have no declared bracket at all.

**The trap this section exists to close.** Cosby's parameters are recorded as
heads in centimetres of water, and a head is a pressure only through the local
gravity. Air entry is a CAPILLARY pressure and the wilting point is a PLANT
pressure, so both are invariant and both of their heads scale by the same factor
when gravity changes. The closure reads them only through their ratio, so done
consistently the wilting point does not move at all: field capacity is the only
state whose defining pressure carries gravity. Converting one of the two
invariant heads and not the other mixes two frames and moves the wilting point
for a bookkeeping reason. The enforcing script evaluates the closure in PRESSURE
so that a mixed frame cannot be written, and checks the two consistent frames
against each other and against the mixed one, which must differ.

## Uncertainty is CORRELATED, and that is a declaration

Three named cases -- low, central, high -- applied coherently across every layer
and every property at once. Never sampled per layer.

The reason is physical. A soil that holds more water at one depth holds more at
the next, because the same texture, the same organic content and the same
endmember uncertainty produced both. Drawing each of fifteen layers
independently averages the spread away by roughly the square root of fifteen and
reports a precision the description does not have. The contract check refuses a
declaration that says otherwise.

**The adopted central case has no declared bracket, and saying so is the
point.** Its uncertainty is Cosby's within-texture-class residual variance,
which Cosby et al. report as large -- larger than everything else in this
contract's uncertainty put together -- and which no artifact here carries. LSHY-1
owns it.

The declared endmember brackets in `pedology/config/pedogenesis.yaml` bracket
something else. They carry their own argument: DECLARED values with Saxton and
Rawls (2006) as the check on their ordering and magnitude rather than their
origin, because that regression parameterises a continuous mixture and excluded
samples above 60 percent clay, so a pure endmember is outside the data it was
fitted on. What they parameterise is pedology's endmember mixture, which is the
derivation that was NOT adopted. The soil map still carries that capacity as
pedology's own field, so the bracket is still what THAT number is worth; it says
nothing about the states a consumer installs.

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

**The frame check**, and it is the one this contract most needed. Cosby's
parameters are heads recorded on Earth, and two of the three defining pressures
are invariant under a change of gravity while one is not. The check evaluates
the states three ways: in an Earth-equivalent frame, with every head in
centimetres of water on Earth and the field-capacity head raised so it still
names the pressure a column of the declared length exerts here; in a Vesper
frame, with every Earth-recorded head divided by the gravity ratio so it names
the same pressure as a head of water on this world; and in the frame-mixed
variant that converts the plant pressure and leaves the capillary one alone.
The first two are required to agree exactly, because the closure reads the two
heads only through their ratio and a frame change multiplies both by the same
factor. The third is required NOT to agree, because a check that cannot fail
proves nothing, and the size of its disagreement is what a consumer would have
shipped had it converted one head and forgotten the other. The check also
requires the adopted pressure derivation to reproduce the head derivation, since
they are one closure written twice.

**The air-entry clamp.** Cosby's equation holds only below the air-entry value;
at and above it the pore space is full and the water content is saturation.
Unclamped, the closure returns a field capacity above saturation wherever the
air-entry suction falls below the field-capacity suction, and substituting the
regressions and `silt = 1 - sand - clay` gives that region as a straight line in
the texture simplex, `1.58 * sand + 0.63 * clay < 0.17`, which needs no search.
WORLD-NGA10 put the clamp in `get_mineral`; WORLD-OF6N moved the derivation here
and the clamp with it, so both defining pressures are clamped at the air-entry
pressure and the CODE is what keeps this right rather than the map. The line is
derived at Cosby's own field-capacity suction and this world's is larger, so the
region where the clamp binds is strictly smaller than the line and the margin
reported is the conservative one. The map is measured against it as a MARGIN and
not as a count -- a count of zero says nothing about how close the nearest cell
is, and a map that moved inside would have its field capacity set by the clamp
rather than by its texture.

## The central case, and the two column extents that had to agree first

The removal this contract was built to make attributable is WORLD-OF6N, and it
needed a decision about the world rather than about the code: what the central
case IS. It is the contract read literally -- the named closure,
Clapp-Hornberger on Cosby Table 4, evaluated at the suction the potential
convention declares, `rho_w * g * L` at this planet's gravity. Both the closure
and the pressure were already declared here, so nothing about the adopted case
is a new choice; what was open was only which gravity the suction was evaluated
at, and a contract that declares field capacity as a drainage equilibrium has
already answered that.

Two alternatives were costed and neither was adopted. The same closure at
Cosby's own suction reads the contract's closure and ignores the contract's
potential convention. Pedology's endmember mixture is not a retention curve at
all -- it carries no saturation, no air entry and no exponent -- so adopting it
would have replaced the named closure and left saturation and matric potential
with nowhere to come from. That makes it a change to the declaration rather
than a choice between consumers, which is why it was never a symmetric option.

**The gravity question and the column-extent question are one decision, and the
extent one was the larger.** WORLD-SLPA raised the gravity conversion from the
LPJ-GUESS side and is subsumed here, because landing it beside the adopted case
would have applied the conversion once and a frame-mixing artifact once.
WORLD-7702 is the other half: pedology's capacity followed regolith depth past
the declared column base while LPJ-GUESS's stopped at it, and that disagreement
was a factor of three across the map against the gravity shift's few per cent.
It is settled by cutting pedology's capacity at the base, which is the
definition of a plant-available capacity rather than a truncation of one, and
by the adopted capacity being an integral over the declared layers. Both
consumers now cover the same column by construction.

The water that cut drops is not lost physics. It belongs to the aquifer term
this column already declares it stops at, and
`hydrography/config/land_water_ledger.yaml` registers
`transient_saturated_storage` as an UNREPRESENTABLE owned by PLHY-4 because the
groundwater solver is steady-state and carries no storage change. Deepening the
declared column instead -- 26 layers of 100 mm to cover the regolith p90 and 43
for the p99, against `NSOILLAYER` appearing across six LPJ-GUESS modules and a
`rootdist` in the instruction file that must carry one component per layer --
becomes worth doing when that solver gains storage, and not before. That is the
successor to this decision, not an alternative to it.

The allophane term is not what separated the two derivations, which is worth
saying because the materials block above reads as though it is. It contributed
a median of nothing and a p90 of a few per cent of pedology's capacity, against
a median absolute gap of tens of millimetres, and it is zero on more than half
the cells. What separated them is that the declared endmember capacities and
the Cosby retention inversion are different numbers for the same texture. The
adopted closure carries no allophane term either, so the andic hydraulic effect
is now a declared absence rather than a disagreement between consumers.

Two things the report cannot supply, and says so rather than estimating. The
change in land runoff ratio and in E over R needs a paired baseline, and
`probe_runoff_response.py` needs a climatology to back-solve potential
evaporation from; `config/planet.yaml` declares `baseline_climatology` null. And
nothing on the LPJ-GUESS side is execution-verified: it does not build on this
tree until a baseline climatology exists, so every LPJ number here is a
transcription of the vendored source evaluated offline, and the changes to
`soilinput.cpp` and `vesperinput.cpp` are checked for syntax and read for
meaning and nothing more.

## What is still undeclared, and who owns it

Eight properties carry the sentinel in the declaration, and the table below adds what adopting one closure left undeclared. Each names its issue, and `--strict` is the
arm that refuses.

| what is missing | owner |
| --- | --- |
| the residual water content, which both active evaporation limiters pin at zero without saying so | LSHY-3 |
| saturated hydraulic conductivity; no artifact here carries one, and LPJ-GUESS's `perc_base` is a dimensionless drainage exponent that cannot be converted into one. That is why `perc_base` stayed in `soilinput.cpp` when the retention states left it: what this contract supplies is the exponent the fit reads, not a conductivity | LSHY-3 |
| unsaturated conductivity: the family is declared, the number is not | LSHY-3 |
| infiltration capacity, absent from both columns, which is why saturation excess is the only runoff mechanism | LSHY-3 |
| the frozen-pore impedance EXPONENT; the form is declared and is used by nothing, because the scheme that would need it has no conductivity to reduce | LSHY-5 |
| composition-dependent thermal properties, so that one material property is not split into an organic profile for water and a constant for heat. The PHASE half is now declared: freeze and thaw exchange the latent heat of fusion with the soil temperature layer that contains the water layer's midpoint, and the sensible heat capacity is still a constant | LSHY-5 |
| the organic fraction and one vertical rule for it | SDEC-10 |
| Cosby's within-texture-class variance, which is the uncertainty of the ADOPTED central case and the dominant term. Until it exists the adopted states carry no bracket at all | LSHY-1 |
| the andic hydraulic effect. The adopted closure has no allophane term, so the material is declared, its fraction is read, and nothing in the retention curve responds to it | LSHY-1 |

## What this contract does not own

Aquifer transmissivity. A vadose-zone unsaturated conductivity is not an aquifer
permeability: they live at different material and spatial scales, they are
measured on different things, and copying one into the other is a category error
with plausible units. The groundwater solver owns transmissivity as `K` times
saturated thickness, scales it separately, and is joined to this column at the
declared contact through the ledger's `drainage` and `capillary_rise` terms and
through nothing else. The contract check refuses a declaration that says the two
are scaled together.
