# The land water ledger

**Recorded:** 2026-08-24
**Enforced by:** `hydrography/scripts/land_water_ledger.py`
**Declared in:** `hydrography/config/land_water_ledger.yaml`
**Finding it rests on:** `biosphere/notes/soil-land-surface-hydraulic-consistency-audit.md`, findings 2 through 4

Worldbuilding. Vesper is an invented planet. Everything below is about the
simulation of it: a declared graph of control volumes, the water transfers
between them, and the ownership rules that decide which component may advance
each store. No model was run to produce any of it.

This is a CONTRACT, not a finding. What is true about the repository's land
water handling is in the audit; what to do about it is in `bd`. This document
says what the ledger IS, what it conserves, and what fails when it does not.

## Why a ledger and not a water balance

The pipeline already has a water balance. `P - E` closes at equilibrium, the
groundwater solver conserves its own recharge against its own seepage to a
declared tolerance, and the basin code adds only signed exchange to catchment
supply. Every one of those is correct on its own terms, and together they do
not show that the coupled system conserves water.

The reason is that a balance answers "does the total come out", and the question
the coupled column needs answered is "which store was debited, by which
component, over which interval". Those are different questions, and finding 3 of
the audit is what happens when the first is mistaken for the second: one
interval's positive `P - E` over one cell is routed as ExoPlaSim's bucket
overflow, taken again by the surface-water builder as catchment runoff, and
taken a third time by the groundwater builder as recharge. Each consumer is
internally consistent. The quantity has been spent three times.

A ledger cannot do that, because every unit of water has one source, one
destination and one place it currently sits, and drawing on a store harder than
it holds is refused rather than clipped.

## What conserves, over what control volume

**Substance.** Water mass, always. Never a depth and never a volume. A depth is
a mass only once an area is named; a volume is a mass only once a density and a
phase are named, and the simulated snow at 330 kg/m3 in ExoPlaSim and at 275 to
500 kg/m3 in LPJ-GUESS is the same water in different volumes. Phase is a
property of a term, never a different quantity.

**Unit.** Kilograms of water, absolute, per declared interval. A per-area depth
is converted using the area of the node it is stated on, and that area is an
argument rather than an assumption. Three denominators appear in the graph --
the land cell, the open-water fraction of a cell, and a coastline segment --
and a script that reconstructs a mass from another script's depth with its own
area changes the mass silently.

**Time.** Earth days, on the contract in
`biosphere/notes/time-base-unit-contract.md`. Closure is required at EVERY
interval down to the model timestep, not only at the reporting interval: a
ledger that closes annually and not daily is the finding rather than the fix.
Each term therefore declares an `interval_floor`, the shortest interval its
producing artifact can supply it on, and a term reaching a store that empties
and refills inside the reporting interval cannot be produced only as an annual
mean. An annual total is not a partition.

**Control volumes.** Ten nodes, declared with their kind, their area basis and
their owner. Four kinds, and the kind decides what the closure check does:

| kind | what it means |
| --- | --- |
| `boundary` | water may ENTER the domain here; its stock is not tracked and what crosses it counts as external input |
| `reservoir` | an ordinary control volume: inflow, outflow, a stock |
| `terminal` | water stops here FOR THIS LEDGER; no term may take it as a source |
| `interface` | owned by another component; the ledger books the crossing and not the pools behind it |

The vertical cut is `column_base_m`: where `soil_liquid` stops and
`groundwater` begins. It is a declaration and not a measurement, and it is
taken from the only column in the pipeline that has a depth at all: LPJ-GUESS's
fifteen 100 mm layers, five upper and ten lower, reaching 1.5 m. ExoPlaSim's
`dwmax` is a capacity in metres OF WATER and carries no column thickness;
pedology's regolith depth scales LPJ's per-layer capacities without moving a
layer boundary or a root. It is deeper than the 1.0 m root zone the abiotic
nutrient ledger cuts its solute pools at, because a root zone is where uptake
happens and a vadose column is where water moves. The lateral cut is the climate-grid
land cell, and land comes from `surface_class`.

The nodes that carry the argument:

- **`canopy`** has no owner. ExoPlaSim's land column has no interception store
  at all, so precipitation reaches its bucket undiminished; LPJ-GUESS removes
  `patch.intercep` before its own snow partition. Declaring the node with no
  owner is what makes the absence countable. Leaving it out of the graph would
  make interception someone's evaporation without saying whose.
- **`snowpack`** and **`soil_liquid`** each carry a shadow copy. A shadow copy
  is not a second store; it is the same water counted twice, and it is the
  mechanism by which two vegetation-capable models evaporate one millimetre.
- **`open_water`** is the node whose area basis differs. Lake evaporation is a
  depth over the open-water fraction and catchment runoff is a depth over the
  land cell, and averaging them into a cell mean is finding 7.
- **`vegetation`** is an interface, and transpiration crosses it twice: once as
  `root_uptake` from the soil into the plant, once as `transpiration` from the
  plant to the atmosphere. Two terms rather than one so that a transpiration
  flux with no matching uptake leaves a residual at the interface instead of
  being absorbed.
- **`unowned_source`** is a boundary that should carry nothing, and carries
  nothing: the one crossing ever booked here is now a declared absence with a
  measured bound. It stays because a ledger that cannot name an unowned
  crossing cannot report one, and because it is the destination one of the
  graph mutations is built on.
- **`ocean_export`** is terminal, and terminal means this ledger stops. It
  carries water across the coastline and claims nothing about what happens on
  the other side.

**The identity.** For every node `n` over one interval, with everything in kg of
water,

    S_n(t+dt) - S_n(t) = sum of terms into n  -  sum of terms out of n

and over the whole domain,

    (water held in reservoirs) + (water at terminal nodes)
        = (water held at the start) + (water across boundary nodes)

Every term has exactly one source and one destination inside the declared node
set, so water is neither created nor destroyed by a transfer, and the domain sum
is a real invariant rather than a bookkeeping convention.

**The second invariant, on the same terms.** Each evaporative component --
canopy evaporation, soil evaporation, transpiration, snow sublimation,
open-water evaporation, groundwater uptake -- is debited by exactly ONE term
from exactly ONE store. This is not implied by the mass identity: two models
withdrawing the same millimetre from the same store conserve mass perfectly and
are still wrong, so ownership is checked separately and can fail separately.

**The energy axis, carried only as far as it must be.** A phase change costs
energy, so every term declares a transform and the latent heat that transform
consumes or releases, and the pairing is checked. ExoPlaSim declares a latent
heat of vaporisation and of sublimation and never declares one of fusion:
`landmod.f90`'s snowmelt uses the difference directly. Fusion is therefore not
an independent constant here, and the ledger holds the identity
`sublimation = vaporisation + fusion` as a check, so a fourth constant appearing
somewhere that disagrees with the three that exist is caught rather than
averaged.

## The check that would fail

A balance nothing can violate is not a check, and neither is a structural check
that has never been shown to catch anything. Both arms are exercised on every
invocation.

**Eleven closure fixtures**, ten of them built to be wrong in a named way. A
fixture that does not get the verdict it was built for is a defect in the
checker and exits non-zero even when the declaration is clean.

| fixture | what it does | required verdict |
| --- | --- | --- |
| `balanced` | rain and snow in, through canopy, pack, soil and channel, out at the coast, every store left with a stock | closes |
| `shadow_store_evaporation` | one model's soil evaporation debits the store; another's transpiration is credited to the atmosphere and debits a store the domain does not have | residual |
| `double_debited_component` | one evaporative component withdrawn by two terms from two stores | refused |
| `pme_triple_ownership` | one interval's positive `P - E` spent as bucket overflow, as catchment runoff and as recharge | refused |
| `area_basis_mismatch` | lake evaporation debited on the open-water area and credited on the land cell | residual |
| `interval_mismatch` | an instantaneous rate integrated over two different durations across a handoff | residual |
| `unnamed_phase_change` | melt booked as a move: ice in, liquid out, no process named | refused |
| `wrong_latent_heat` | ice taken straight to vapour carrying the latent heat of vaporisation | refused |
| `terminal_source` | water drawn back out of `ocean_export` | refused |
| `store_below_zero` | a withdrawal larger than the store holds | refused |
| `negative_transfer` | a transfer with a negative mass, standing in for a reversed sign convention | refused |

`shadow_store_evaporation` is the one worth reading, because it is finding 4
written down. Both withdrawals are physically real and both models are
internally consistent. What the ledger sees is that the domain was credited two
evaporative fluxes and debited one, because the second store is not in the
domain, and the residual is exactly the water the coupled system does not have.

`store_below_zero` is the second, because it is the shape `landmod.f90`'s
floor at zero would have if the climate column could reach it. It cannot. The
sweep in `bucket_floor_bound` is the third arm of this check and it closes that
crossing rather than bounding it loosely: three properties of the model make
the store land at exactly zero and never past it, so what the floor can create
is the rounding of a difference that cancels. The fixture stays because the
shape is general -- any coupling that draws on a store harder than it holds has
it -- and refusing it is what the ledger does instead of clipping.

**The unowned crossing, measured.** `bucket_floor_creation` was the one term
booked from the `unowned_source` boundary, and it is now a declared absence
with a bound rather than a term with a sentinel. The three properties that
close it are `fluxmod.f90`'s AMAX1 against the level humidity, which stops the
surface evaporation ever adding water; the cap on that same withdrawal at the
store divided by the timestep, applied on every land cell whatever the snowpack
holds; and `rainmod.f90` building the snowfall rate as a subset of the total
precipitation rate, so precipitation minus snowfall is rain and cannot be
negative. The snow-exhaustion path this was first attributed to is the wrong
suspect: the cap is unconditional, so a snowpack absorbs part of a withdrawal
already limited against the soil store and leaves the store fuller than bare
ground would. The report carries the residue per land cell and timestep, the
ceiling over one orbit if every land cell clipped at every step, and that
ceiling against the land soil water stock. Two arms remove one closing property
each and are required to produce a macroscopic crossing, which is what makes
the bound a measurement rather than a null result; change any of the three
properties and this becomes a term again.

**Twelve graph mutations.** `check_graph` returning nothing on the real
declaration proves nothing on its own: a check that cannot fail and a
declaration that is correct look identical from there. So the declaration is
loaded twelve more times, broken in one named way each, and the check is
required to catch every one -- a dangling destination, a term drawing on the
terminal node, a phase change with no process, a melt term carrying the latent
heat of vaporisation, a transform whose phase pair is reversed, one evaporative
component on two terms, a reservoir with no outflow, a reservoir with no
inflow, a store with no owner, a node with no area basis, a fusion constant
that disagrees with the difference of the other two, and a boundary receiving
from inside the domain. A mutation the check passes is a defect in the check.

## What refuses, and until when

Every term carries the `undeclared` sentinel today. The ledger is
DEFINED and does not CLOSE, and the two are different states. `--strict` refuses
while any of three conditions stands:

**Terms with no flux.** Each names the issue that owns it.

| what is missing | owner |
| --- | --- |
| the climate-side land column that produces infiltration, storage, surface runoff and drainage as separate terms | LSHY-3 |
| the authoritative snow, melt, soil temperature and liquid/ice chronology | LSHY-5 |
| vegetation demand and fulfilled uptake by source, without a second withdrawal | LSHY-4 |
| the signed groundwater lower-boundary exchange | PLHY-4, GW-5 |
| open-water and lake fluxes at a sub-annual interval on the open-water denominator | LSHY-6 |
| routed dissolved and particulate transport on this ledger | ANUT-6 |

**Stores with a shadow copy.** Four nodes have one. `snowpack` is held as
ExoPlaSim's `dsnowz` and again as LPJ-GUESS's `soil.snowpack`; `soil_liquid` as
`dwatc` and again as `soil.wcont`; `soil_ice` has no owner at all in the climate
column and is held only by LPJ-GUESS's `Frac_ice`; `open_water` is hydrography's
lake storage and also ExoPlaSim's `driver`, which accumulates routed runoff and
discharges it to the ocean without ever evaporating.

**Terms whose producer cannot resolve their destination.** Three today, and all
three are the same shape: an annual-mean artifact feeding a store that changes
inside the year. `lake_precipitation` and `baseflow` reach `open_water`, and
`capillary_rise` reaches `soil_liquid`, from a steady-state groundwater solution
and an annual surface-water forcing. The fix is not a better annual mean.

## Named absences

An absence is not a term with an undeclared flux. A term with `undeclared` has a
place in the graph and no number; an absence has no place in the graph, because
the state it would move does not exist in either model. Eight are declared, each
with what is missing, what it makes unrepresentable, and its owner. The two that
constrain other components hardest:

- **Soil water above field capacity.** LPJ-GUESS's `wcont` is a fraction of
  AVAILABLE capacity and is refused outside `[0, 1]`, so the wettest state the
  column can represent is field capacity and there is nowhere for water between
  field capacity and saturation to sit. Water-filled pore space derived from it
  therefore has a ceiling that is a property of the texture rather than of the
  weather, and near-saturated soil is a state the column cannot enter. On this
  world that is not an edge case: a quarter of the land is playa clastic, so
  wetting and drying through near-saturation is a common condition. LSHY-1 owns
  it, and any process keyed on water-filled pore space inherits the ceiling.
- **Open-water evaporation from routed water.** Routed runoff accumulates in
  `driver` and discharges to the ocean; `dwatc` gains only from local
  precipitation minus evaporation. A cell's annual evaporation is therefore
  capped by its own precipitation however the bucket depth is set, and a
  terminal lake -- whose entire balance is evaporation of water its catchment
  delivered -- is not representable. Raising `dwmax` on the lake fraction buys
  the seasonal partition and not the annual total, and the builder says so.

## What would settle the parts that are open

Nothing here needs a climate run or an LPJ-GUESS run. The declaration, the
fixtures and the graph mutations are deterministic, and the ledger is the
artifact the other rows book into rather than a result about the world.

What closes it, in the order the audit recommends: LSHY-1 supplies one
vertically explicit property artifact so both consumers read the same column;
LSHY-3 replaces the scalar bucket behind a selectable reduction and produces
infiltration, storage, surface runoff and drainage as separate terms with
separate owners; LSHY-5 makes snow and soil phase one state rather than two;
LSHY-4 closes vegetation demand and fulfilled uptake so `root_uptake` and
`transpiration` are the same water; PLHY-4 supplies the signed lower-boundary
exchange that turns `drainage` and `capillary_rise` into a contact rather than
an outlet; and LSHY-6 supplies the subgrid response that gives `open_water` its
own denominator back.

Until then the ledger's value is that it says, of any proposed coupling, which
store is being debited and by whom -- and refuses the answer "all of them".
