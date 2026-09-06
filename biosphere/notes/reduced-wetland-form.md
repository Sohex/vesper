# The reduced wetland, peat and methane form, and what each part costs in claims

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's dormant peatland and methane modules,
the hydrography artifacts that would say where its simulated wetlands are, and
the output tables an activated run would write. Nothing here has been run. No
LPJ-GUESS simulation this project has made has passed acceptance.

WET-1 through WET-11 specify a TARGET: vertically resolved substrate,
temperature, saturation and redox; production and oxidation exposed as separate
fluxes; stock-conserving peat transition semantics; a mutually exclusive surface
partition; and a closed atmospheric methane budget. That chain rests on SDEC-3's
vertical soil coordinate, PLHY-4's water exchange, PLHY-5's crossing, PCAR-5's
trait registry and BVOC-6's oxidant bracket, all of which are open. A standard
that does not ship is worth less than a design that does, and PALADYN is a
design that works, is published and runs at this project's throughput. So this
document declares, per part, what is taken INSTEAD when the target's dependency
is not there, what taking it forbids the project to claim, and the condition
under which a part falls back. It is written before anything is built to it, so
the comparison is made while neither option is under schedule pressure.

`biosphere/notes/wetland-activation-contract.md` says what has to be declared
before the switches may move. This document says what the declaration is allowed
to be when the target cannot supply it, and what that costs. The enforcement is
`biosphere/scripts/wetland_gate.py`; the declaration is
`biosphere/config/wetlands.yaml`; the finding underneath both is
`biosphere/notes/wetlands-peat-methane-audit.md`.

## 1. The three parts are not equally viable, and extent has no floor at all

The three parts of the reduced form are EXTENT, PEAT STOCK and METHANE. Two of
them fall back to a published design at a nameable cost. The third does not fall
back, and that is the substance of this document.

PALADYN's wetland extent IS a TOPMODEL saturated area fraction, taken wherever
the simulated surface is snow free, with potential peatland the fraction wet for
at least three months of a simulated year. This project has no such fraction.
The saturated-area closure that would have produced one is withdrawn in
`hydrography/config/topographic_index.yaml`, on the verdict of a score that
config declared before any fraction was computed;
`hydrography/scripts/build_topographic_index.py` implements no closure and
`hydrography/scripts/build_wetness.py` forms no saturated class.
`hydrography/notes/subgrid-water-table.md` section 7 measures why no narrower
criterion could license one: on every arm of both Earth bore sets the gain the
terrain half carries over the cell-mean depth is at or below the scatter the
support puts on it, and a narrower criterion has less support rather than more.

So the floor and the target rest on the same missing quantity. The reduced
extent is not a coarser version of WET-2's partition; it is a SMALLER
PARTITION, with one of its five classes absent and no design anywhere in this
project's survey able to supply it. That is why the recorded reasoning had to be
redone: extent was recorded as the part that forfeits least, forfeiting only
mutual exclusivity, and mutual exclusivity is the one thing the reduced extent
actually keeps.

Two routes reopen the class and the gate admits only these two: a new score
whose SUPPORT can resolve a gain of a few thousandths of an area under the
curve, whose terms belong to `hydrography/notes/subgrid-water-table.md` section
7; or this reduced form, which is what `world-4fr6` selects between.

## 2. Extent: the reduced form

WHAT IS TAKEN. The four classes the wetness classification resolves on the
native mesh, crossed to the climate grid as area shares by
`lib/gridding.py`'s categorical reduction:

- persistent peat-forming land, from the simulated peatland stand's own state;
- seasonal inundation, from the periodic lake cycle's per-region wet-bin count,
  which is the closed-basin third of that quantity and no more;
- open lake and playa water, from the same cycle's always-wet set;
- dry mineral soil, the residual.

WHAT IS DROPPED. The saturated non-inundated mineral class, at every support.

THE CONDITION. There is none to wait for. The other two parts fall back when a
consumer needs a field and a dependency is missing; this part has no full form
to fall back FROM, because the only published design for it needs the same
withdrawn fraction. It is taken now and it is reopened only by a revived closure
under the route above.

THE COST LEDGER.

**Dropping the saturated non-inundated mineral class costs the claim that this
world's simulated land carries any wetland that is not standing water.** Every
wetland the reduced extent can express is a lake surface or its seasonal margin.
A fen or a bog on a groundwater-fed slope with no open water has no area at any
support, so the modelled wetland extent is a LOWER BOUND on wetland extent
rather than an estimate of it, and the extent arms `wetlands.yaml` registers can
only be occupied on one side: `extent_low_end` is reachable and `extent_high_end`
is not a bracket end this form can reach by widening a number.

**It costs the claim that the dry mineral share is dry mineral soil.** The
partition is closed, so the land the absent class would have taken stays in the
residual. A per-cell `dry_mineral` share is therefore an UPPER BOUND, and any
flux computed as a dry-surface rate times that share overstates the area it acts
on by exactly the saturated area. This bites hardest on the near-surface aerobic
methane sink, whose whole subject is the dry surface.

**It costs the two-support structure and the constraint that went with it.**
`hydrography/config/wetness.yaml` states exclusivity at two supports because the
saturated share was the one class that could only arrive at the climate grid,
taken OUT of the resolved classes rather than added beside them. With that class
absent, every class in the partition is resolved on the native mesh, nothing
crosses from a coarser support to a finer one, and `extent.downscaling_rule` has
nothing to rule on. The rule and its refusal are kept, because a revived closure
brings the second support back with it.

**Dropping the floodplain and saturated-soil thirds of seasonal inundation costs
the claim that a seasonal inundated area is a planetary total.**
`hydrography/notes/land-water-ledger.md` splits the quantity in three and only
the closed-basin third is published. A floodplain's inundated area needs a
height-above-nearest-drainage distribution and a routing model, neither of which
exists here; a seasonally saturated soil is a water content rather than an area,
and forming one from a water content is what
`hydrography/config/wetness.yaml`'s `area_proxy_from_cell_mean` refuses. So the
seasonal share is a lower bound too, and it may not be compared with an Earth
seasonal inundated area as like with like.

**Dropping river water costs the claim that the open-water class is open
water.** `surface_water.nc` carries a discharge and no channel width, so a river
area cannot be formed from it and discharge stays a rank and a mask. A lake
evaporation taken over the open-water area therefore omits channel evaporation,
and a methane budget that assigns an aquatic source per unit open water assigns
none to rivers.

## 3. Peat stock: the reduced form

WHAT IS TAKEN. PALADYN's structure. Acrotelm and catotelm confined to the top
soil layer; transfer from one to the other at a critical acrotelm carbon;
catotelm carbon shifted to lower layers as density fills one; peat carbon
vertical diffusivity zero; areal change limited to a fixed fraction of a
simulated year with a minimum fraction seeding every cell.

WHAT IS DROPPED. WET-4's vertical resolution of the peat column.

THE CONDITION. This part falls back when a consumer needs a peat carbon stock
and SDEC-3's vertical soil coordinate does not exist. The consumer that reaches
it first is the C-N-P closure `aresidual_cnp.out` reports, through
`apeat_stock.out`.

THE COST LEDGER.

**Dropping the vertically resolved peat column costs every depth-resolved redox
claim.** Nothing the reduced form emits can say at what depth the simulated peat
is anoxic, where in the column production happens, how far a methane molecule
travels through unsaturated peat before it is oxidised, or how the profile
responds to a water table that moves within the column. A peat stock is one
number per simulated cell reported against a declared depth bracket, and the
bracket is a declaration rather than a measurement.

**Setting the peat carbon vertical diffusivity to zero costs the claim that the
simulated peat column has a memory of a water table that fell.** Carbon does not
move between the acrotelm and the catotelm except through the transfer rule, so
a catotelm that is re-aerated by a falling table cannot lose carbon by any path
the model contains, and a drying trend appears in the simulated stock only
through the areal term.

**Capping areal change costs the claim that the simulated peatland area responds
to climate on any timescale shorter than the cap allows.** With a minimum
fraction seeding every cell, the simulated peatland is present everywhere at a
small area and the cap is what sets how fast it grows, so a modelled expansion
rate is the cap and not a result. The cap and the seed fraction are numerical
devices, and this document does not present either as physics.

**Taking this part imports three constants whose disposition has to be recorded
before it is taken.** The critical acrotelm carbon at which transfer happens is
sourced from Earth northern peatlands and is IMPLICIT-EARTH by this project's
vocabulary: the derivation is sound and is for the wrong planet, so Earth's
figure is a comparison to report the distance from and never a target to solve
onto. The areal-change cap and the seeding fraction have no derivation at all
and are irreducible numerical devices, which is the fourth disposition
`CLAUDE.md` allows and the one that has to be argued rather than assumed.
`notes/audits/tuned-values.md` is where all three are enumerated when this part
is taken; taking it without that enumeration is what makes a fallback a silent
gap.

## 4. Methane: the reduced form

WHAT IS TAKEN. A constant fraction of simulated heterotrophic respiration
wherever that respiration is anaerobic, with one fraction per surface class of
the partition in section 2.

WHAT IS DROPPED. WET-6's separation of production from oxidation, and WET-7's
diffusion, plant-conduit and ebullition pathways.

THE CONDITION. This part falls back when a consumer needs a surface methane flux
and WET-6's substrate interface does not exist. It is BRACKETED and never a
measurement, and the bracket is the one `wetlands.yaml` floors at the spread
between the two CH4:CO2 constants the vendored fork itself carries: a
mechanistic ratio fitted at northern Earth sites against an emission factor
whose own source comment says it was changed to match a global Earth total.

THE COST LEDGER.

**Dropping the separation of production from oxidation costs the claim that any
number this form emits is a production rate or an oxidation rate.** The constant
IS an emission factor with oxidation folded into it, so what the model can
report is a net surface flux and nothing else. The specific test WET-6 exists to
make possible -- telling a plausible total apart from compensating production
and oxidation errors -- cannot be run at all under this form, so an agreement
between the modelled total and any other total is uninformative about mechanism.
`mch4_production.out` and `mch4_oxidation.out` therefore stay unemitted under
this form rather than being filled with the halves of a constant.

**Dropping transport costs the claim that the modelled seasonal methane cycle
has a shape of its own.** With no diffusion, no plant conduit and no ebullition,
the emitted flux is the seasonal shape of simulated heterotrophic respiration
multiplied by a constant. There is no winter storage and release, no thaw pulse,
no dead-tiller venting and no pressure- or wind-driven bubble release, so the
`transport_partition` arm has nothing to partition and the three pathway tables
carry no quantity this form produces.

**Dropping the dry-soil sink and the aquatic sources costs the claim that the
ledger is complete.** The reduced form emits a source over anaerobic surfaces
and nothing over the rest, so the modelled net surface exchange is a gross
source, and the distance between a gross source and a net exchange is the sink
this form does not have.

**THE PROHIBITION THAT SURVIVES INTACT.** This form yields a surface FLUX under
a transplanted-Earth assumption and never an abundance. An abundance needs the
oxidant and lifetime calculation WET-10 owns and this form does not contain, so
a flux computed this way MAY NOT move `atmosphere.pCH4_bar` in
`config/planet.yaml`, whatever its magnitude. `wetlands.yaml`'s
`atmosphere.trace_gas_state_closed` stays false under this form and the gate
refuses the combination by name. WET-10 has a reduced form of its own and it is
WET-10 that names it.

## 5. What a reduced part is labelled with, and what a consumer may do with it

A fallback is labelled ON THE ARTIFACT, so nothing downstream can read a reduced
field as a full one. Two requirements, and the gate holds both.

**A class declared absent carries `absent_because` and `licensed_by`.** The
reason is what the class needed and what happened to it; the licence is the row
that decided the reduced form. A class dropped with only one of the two is
refused, because a reason without a licence is an omission with an excuse
attached and a licence without a reason is a decision nobody can check.

**A bracketed class enters a downstream ledger as BOTH arms or as neither.**
`biosphere/config/wetlands.yaml` carries it as `extent.convention_arm_rule` and
the gate refuses a rule that admits a single arm. It is stated before any share
is written, because that is the constraint `world-4fr6` puts on either route. The saturated-area closure's `f_grad` is a
CONVENTION bracket rather than a spread around a central value: at one end of it
the closure's ranking is the depth's ranking exactly, so that arm contains no
terrain information at all. A consumer that takes one arm of a convention
bracket is reporting a single number that has no central case behind it. The
rule therefore binds a revived closure as well as this reduced form: a class
whose share is a convention bracket propagates as two arms through every ledger
that reads it, or it does not propagate.

## 6. What the whole reduced form cannot claim, in one place

- That the simulated wetland extent is an estimate. It is a lower bound whose
  missing part is every saturated surface that is not under standing water.
- That the simulated dry mineral area is dry. It is an upper bound by the same
  quantity.
- That any emitted methane number is a production, an oxidation or a transport
  flux. It is a net surface flux carrying a bracket at least as wide as the
  spread between the two constants the fork already holds.
- That the simulated peat has a depth structure, an age, or a response to a
  falling water table other than through its area.
- That this world's atmospheric methane is its own. It stays prescribed.
- That any total emitted under this form can be compared with an Earth total as
  like with like. Extent, seasonal area and the ledger are each incomplete in a
  named direction, and the comparison would be between a bounded quantity and a
  measured one.

Each line above is a claim the project may not make while the reduced form is
taken, and each is recoverable by exactly the row named beside its part. None of
them is recoverable by widening a bracket.
