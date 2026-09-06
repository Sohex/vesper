# The abiotic nutrient ledger

**Recorded:** 2026-08-24
**Enforced by:** `biosphere/scripts/abiotic_nutrient_ledger.py`
**Declared in:** `biosphere/config/abiotic_nutrients.yaml`
**Finding it rests on:** `biosphere/notes/abiotic-nutrient-delivery-audit.md`

Worldbuilding. Vesper is an invented planet. Everything below is about the
simulation of it: a declared graph of control volumes, the mass transfers
between them, and the screens that decide which sources may book into it. No
model was run to produce any of it.

This is a CONTRACT, not a finding. What is true about the repository's abiotic
nutrient handling is in the audit; what to do about it is in `bd`. This document
says what the ledger IS, what it conserves, and what fails when it does not.

## Why a ledger and not a table

A table of nutrient sources answers "what is there". It cannot answer the
question the biosphere actually needs answered, which is "how much of element X
crossed into the pool a root can take from, and where did the rest of it go".
The difference is that a ledger BALANCES: every unit of mass has one source, one
destination and one place it currently sits, and the sum is checkable.

The audit's finding 1 is what a table looks like when it is mistaken for a
budget. `pedology/scripts/phosphorus_budget.py` multiplies a phosphorus
concentration in ppm by a relative lithology release factor, groups the result
by drainage and concentrates it geometrically over basin floors. Every step is
defensible as a geomorphic hypothesis and none of them carries mass, so the
output cannot be added to anything, subtracted from anything, or checked against
anything.

## What conserves, over what control volume

**Species.** Element mass, always. Never a compound mass: a sulfate mass and a
sulfur mass differ by 2.996 and nothing in a column name says which is meant. A
term's phase says what compound the element is in; the ledger's arithmetic is on
the element.

**Unit.** Kilograms of element, absolute, per declared interval, on ABSOLUTE
time. A per-area rate is presentation and is converted to absolute mass before
anything is summed, using the receiving area the term declares. This is the
single design decision that makes an area-basis error fail rather than rescale:
a dust flux quoted per square metre of receiving land, an erosion flux quoted
per square metre of source hillslope and a basin-floor enrichment quoted per
square metre of floor are three denominators, and a script that reconstructs a
mass from someone else's rate with its own area changes the mass silently.

**Time.** Earth days, on the contract in `notes/time-base-unit-contract.md`.
Never orbital phase; a term stated per orbit changes meaning when the flux
moves, and BIO-24 already owns that correction where LPJ-GUESS meets it.

**Control volumes.** Thirteen nodes, declared in the config with their kind and
their area basis. Four kinds, and the kind decides what the closure check does
with the node:

| kind | what it means |
| --- | --- |
| `boundary` | mass may ENTER the domain here; its stock is not tracked and what crosses it counts as external input |
| `reservoir` | an ordinary control volume: inflow, outflow, a stock |
| `terminal` | mass stops here FOR THIS LEDGER; no term may take it as a source |
| `interface` | owned by another component; the ledger books the crossing and not the pools behind it |

The nodes that carry the argument:

- **`parent_material`** is a finite stock, not a supply. Drawing on it decrements
  it. That is what makes depletion a result the ledger can produce rather than
  an outcome it assumes, and it is what Porder et al. (2007) need in order for
  both fast advection and slow renewal to lower root-zone phosphorus.
- **`regolith_mineral`** stands between every particulate input and the root
  zone. A refractory particle is not a nutrient until something dissolves it,
  and the ledger refuses a term that puts one straight into solution. The
  vendored fork would not: `snow_pinput()` stores deposited phosphorus in snow
  and `pmass_add()` then places it into the labile and sorbed mineral system, so
  a bulk dust phosphorus passed there makes every deposited apatite grain
  immediately available. The repository's own Bodele measurement is 2 to 4
  percent water-soluble.
- **`soil_solution`** is the only node a root takes from, and the only node the
  ledger allows a term to deliver dissolved material to.
- **`endorheic_store`** is a reservoir and not an accumulator. It has two
  outflows -- deflation back to the atmosphere, and burial -- because a closed
  basin that only receives is the shape that makes retention look like
  fertility. Material there may stay dissolved in lake water, precipitate in
  salt, bind to sediment, become biogenic material, or lie beneath an
  inaccessible floor.
- **`coastal_export`** is terminal, and terminal means this ledger stops. It
  carries species and phase across the coastline and claims nothing about what
  they are worth on the other side. "Reached the ocean" is a boundary condition
  and not marine nutrition; OCN-13 begins exactly here.
- **`lithosphere`** exists because the graph check demanded it. Without an
  inflow, `parent_material` is a stock that can only fall, and the check that
  says "a reservoir with no inflow cannot be filled" fired on the first run of
  this ledger and produced the `exhumation` term below.

**The identity.** For every node `n` over one interval, with everything in kg of
one element,

    S_n(t+dt) - S_n(t) = sum of terms into n  -  sum of terms out of n

and over the whole domain,

    (mass held in reservoirs) + (mass at terminal nodes)
        = (mass held at the start) + (mass across boundary nodes)

Every term has exactly one source and one destination inside the declared node
set, so mass is neither created nor destroyed by a transfer, and the domain sum
is a real invariant rather than a bookkeeping convention.

## The check that would fail

A balance nothing can violate is not a check. Seven reduced fixtures run on
every invocation, six of them built to be wrong in a named way, and a fixture
that does not get the verdict it was built for is a defect in the checker and
exits non-zero even when the declaration itself is clean.

| fixture | what it does | required verdict |
| --- | --- | --- |
| `balanced` | import, weather, take up, return, leach, export | closes to 1e-9 kg |
| `area_basis_mismatch` | the same delivery debited on the source area and credited on the receiving area | residual |
| `dangling_destination` | a term delivering to a node the ledger does not have | refused |
| `terminal_source` | a term drawing on `coastal_export` | refused |
| `element_change` | a transfer carrying nitrogen through a phosphorus ledger | refused |
| `stock_below_zero` | a weathering flux drawing harder than `parent_material` holds | refused |
| `negative_transfer` | a transfer with a negative mass, standing in for a reversed sign | refused |

`area_basis_mismatch` is the one worth reading. It is not a conservation error
inside one script; it is what a cross-component handoff looks like when one side
writes a per-area rate and the other reads it and multiplies by an area of its
own. The ledger books the debit and the credit as given, which is why the
difference appears rather than being absorbed.

The structural checks fail for the same class of reason: a transform declared as
`none` that changes phase, a phase change with no process named, a reservoir
with no outflow, a per-area rate at a node with no declared area basis.

## What refuses, and until when

Every one of the twenty-one terms carries the `undeclared` sentinel today. The
ledger is DEFINED and does not CLOSE, and the two are different states. Each
term names the issue that owns its flux:

| what is missing | owner |
| --- | --- |
| absolute rock-to-solution weathering, and the material boundary that says how much rock reacts per area and time | ANUT-2 |
| source-, size- and phase-resolved particulate deposition, with a separate soluble fraction | ANUT-3 |
| composition-aware total-flash nitrogen fixation and its deposition | ANUT-4 |
| a spatial, chronological, species- and wet/dry-resolved nitrogen interface | ANUT-5 |
| dissolved and particulate transport on the accepted land-water and sediment ledger | ANUT-6 |
| exhumation and regolith production rates | ANUT-7, and see below |
| the pedogenic initial mineral state the fork currently initialises to zero | ANUT-10 |
| uptake and return read back across the LPJ interface | ANUT-9 |

The ledger is what those issues deliver INTO. None of them can be checked
against anything until the graph they book into exists, which is why this comes
first.

## Where a rate has to come from, given the terrain has no time axis

Two of the retained sources are rates, and `docs/src/reference/no-time-axis.md`
settles where a rate may come from here: World Orogen has no time axis at all,
so a duration is UNDEFINED rather than unmeasured, and
`params.hydraulicErosion`, `thermalErosion` and `glacialErosion` are
dimensionless intensity sliders that can ORDER cells and cannot SCALE them.

**Exhumation and regolith production** therefore take their absolute scale from
Earth's cosmogenic denudation compilation, on the expected-value argument that
document sets out: Earth is also one randomly chosen moment of a stationary
population, and so is this terrain. Portenga and Bierman (2011) give igneous
outcrops 8.7 +/- 1.0 m/Myr and sedimentary 20 +/- 2.0 m/Myr. Heimsath et al.
(1997) supply the depth dependence rather than the scale, `-(de/dt) = (77 +/- 9)
exp(-(0.023 +/- 0.003) h)` um/yr with `h` the soil thickness in cm, so the rate
is a function of the soil depth pedology already carries. The number that
results is calibrated to an outcome and is labelled that way; it is not a
measurement of how old anything is.

**Volcanic ash** takes its rate from the same argument, already made in
`pedology/config/pedogenesis.yaml` for the andic soils and carried here
unchanged: a flood basalt emplaces most of its volume in about a million years
and then persists for hundreds, so a randomly chosen moment finds a given
province erupting with probability of order one percent. The tephra term is
therefore live on ARC provinces only and zero on `flood_basalt` and `oib` -- by
the expected-value argument and not by a measurement.

## The ANUT-7 register: what may book into the ledger

The materiality test is in the config and was fixed before any magnitude was
looked at. A candidate is retained if EITHER its plausible upper-bound delivery
of some element reaches a tenth of the incumbent term for that element at the
same node, OR it is the only term supplying that element at that node. The tenth
is not a taste: the incumbent terms are bracketed by more than a factor of two
-- the lightning source by about five, the soluble fraction of deposited dust
phosphorus between 2 and 4 percent, foliar C:P by a coefficient of variation of
79 percent -- and a term at a tenth of an incumbent that uncertain cannot move
the incumbent's bracket. The second rule exists because "smaller than the
incumbent" is not a reason to drop a source when the incumbent is zero.

| candidate | verdict | rule | pulse timing |
| --- | --- | --- | --- |
| geomorphic renewal (uplift, erosion, exhumation) | retain | only supplier | not material |
| fresh volcanic ash and tephra | retain, arc only | magnitude | MATERIAL |
| marine aerosol | retain for Ca, Mg, K, S; refused for P | magnitude | not material |
| volcanic sulfate deposition | register only | magnitude | eruptive only |
| fire ash | out of scope | -- | FIRE-7 owns it |
| lightning fixation | out of scope | -- | ANUT-4 owns it, and it is already a term |

**Tephra keeps its pulse timing, and that is the one case where it matters.**
Gislason and Oelkers (2003) find glass is the first phase altered in terrestrial
weathering, from a faster intrinsic rate, higher solubility, greater surface area
and better access to undersaturated water. Wolff-Boenisch et al. (2004) put a
1 mm basaltic glass grain at a 500 year lifetime against 4500 for rhyolitic,
with field rates on andesitic ash about twenty times the laboratory ones. So a
fall is a step change in `regolith_mineral` followed by decades to centuries of
elevated release. Dahlgren et al. (2004) then supply the reason a time-average
gets it wrong at both ends: phosphorus availability is HIGH in young volcanic
material and FALLS with development as allophane and Al-humus fix it. An
averaged ash flux gets the sign of the phosphorus effect wrong early and late.

**Marine aerosol is retained for base cations and refused for phosphorus.**
Chadwick et al. (1999) measured atmospheric input directly at a Hawaiian site
and report 600 +/- 400 mg Ca/m2/yr, mostly marine aerosol, and find the
atmosphere becoming the dominant source of calcium, magnesium and potassium on
substrates older than about 1e5 years, with potassium crossing over sooner
because it is more mobile than the alkaline earths. The same paper is why the
phosphorus arm is refused rather than scaled down: the seawater contribution of
phosphorus to marine aerosol is vanishingly small, because marine uptake holds
surface-water phosphorus down. Booking a phosphorus flux there would invent one.

**Volcanic sulfate is registered and not implemented.** Andres and Kasgnoc
(1998) put time-averaged non-eruptive subaerial emission near 12 Tg SO2 per
year, about 6 TgS, which spread over the whole planetary surface is of order
10 mg S/m2/yr before any land-ocean split. Meybeck (1987) Table 5 gives global
riverine sulfate release at 280e12 g/yr, about 93 TgS over exorheic land, two
orders larger per unit area. That is below the materiality fraction, so it is
registered rather than absent -- and the second retention rule would pull it
back in on any surface where the weathering sulfur term is zero, which the
adequacy screen below is what identifies.

**No screen may be carried on an optical quantity.** The checker refuses, by
name, any candidate whose carrier names an aerosol optical depth, an extinction,
or the repository's dust, sea-salt or volcanic-sulfate optics products. An
optical depth folds in refractive index, size distribution and water uptake, and
no elemental mass comes back out of it. The two live traps: the volcanic sulfate
product is built from a DEGASSING inventory and describes stratospheric aerosol,
while tephra is fragmentation and falls proximally -- two materials, two source
processes, two geographies; and the sea-salt product's soil-facing output today
is optics rather than deposition mass.

## The ANUT-8 screen: conservative, and able to come back adequate

The vendored CNP fork can limit the simulated plants by nitrogen and phosphorus
only. A run that is not N- or P-limited has not thereby been shown to be
unlimited, so the question is whether any other element could be co-limiting.

### The bound, its direction, and the failing case, all fixed first

The screen asks whether abiotic supply can build the ecosystem's entire standing
circulating pool of an element, from nothing, within one soil residence time.

    T_acc(X) = Q_high(X) / F_low(X)     Earth years

    ADEQUATE       T_acc <= residence time
    NOT SETTLED    otherwise. A BRACKET, never a verdict of inadequate
    NOT ADEQUATE   only if the OPTIMISTIC supply is below the PESSIMISTIC pool

It deliberately does NOT compare supply against gross uptake. In steady state
most of gross uptake is met by litter return, so a supply-against-gross test
fails for every element including the ones Earth ecosystems plainly have enough
of, and a threshold nothing can pass is not a screen. Running that version
first is what showed it: the critical runoff came out near 100 m per year.

Directions, each chosen so that it runs against the conclusion:

- `Q_high`, the standing pool, is an UPPER bound. Vitousek and Sanford (1986)
  Table 2, the maximum above-ground nutrient content over their moist tropical
  forest compilation rather than a mean; sulfur from the maximum biomass in the
  same table at CENTURY's MINIMUM plant C:S ratio, 190 by mass (Parton, Stewart
  and Cole 1988), the minimum ratio being the maximum sulfur content. A larger
  pool takes longer to build, so up is conservative.
- `F_low`, the supply, is a LOWER bound: dissolved release from rock weathering
  only, with marine aerosol, dust, ash, recycling and every other term set to
  ZERO. All of those add.
- The residence time is `root_zone_depth / denudation rate`, so 1 m of soil at
  Portenga and Bierman's ordinary-rock rates turns over in 5.0e4 to 1.15e5 Earth
  years. The screen runs on the SHORT end, because a shorter residence time is a
  stricter test. Chadwick et al. independently find the atmosphere taking over
  as the dominant base-cation source below 1e5 years, which is the same order
  and is not what set this number.
- The multiplier taking the above-ground pool to the whole circulating pool --
  roots, and the exchangeable pool on soil colloids -- is a BOUND READ OFF A
  FIELD. Pedology emits an exchange complex: `build_soil.py` derives cation
  exchange capacity from the soil map's clay and organic-matter fractions, base
  saturation from its pH, and the exchangeable pool per element over this
  ledger's own `root_zone_depth_m`, which it reads from here so both sides count
  the same column. `pedology/config/pedogenesis.yaml`'s `exchange` block carries
  every source and what each does and does not license.
- THAT BOUND IS LARGER THAN THE BRACKET IT REPLACED, AND THE OLD UPPER END WAS
  NOT AN UPPER BOUND. The 5 came from Chadwick et al. (1999) Fig. 2a, whose
  Hawaiian exchange cations reach about 6 meq/100 g over the top metre -- a
  deeply weathered, base-depleted basaltic profile. This world's simulated soils
  sit above the base-saturation transition that Chadwick et al. (2003) and Solly
  et al. (2020) independently place near pH 5.5 to 6.5, so most of their
  capacity holds bases rather than aluminium. The screen had therefore been
  reporting a critical runoff BELOW what its own construction supports, which is
  a screen permissive where it advertises conservatism. The repair is the larger
  number and not a narrower one.
- THE FOUR ELEMENTS DIFFER BY MORE THAN THE SCALAR CAN CARRY. Sulfate is an
  ANION and the cation exchange complex holds none of it, so sulfur's soil term
  is not the cations' at all: it is the anion exchange capacity of
  variable-charge andic material, small over most of this world's land and worth
  several times the above-ground sulfur pool where the andic fraction is high.
  `belowground_and_exchangeable_by_element` records each element's bound; the
  scalar this script reads is the maximum over them and is conservative for
  every element until world-vf8j moves the read onto the map.
- AND IT CANNOT GO STALE. `build_soil.py` refuses when the declaration, or any
  entry in the map beside it, falls below what the emitted field implies. Soil
  carbon grows the capacity's organic term through loop A, so the check fires on
  the iteration that outgrows the declaration.

Because `F_low(X) = c_min(X) * runoff`, with `c_min` from Meybeck (1987) Table
2C -- the per-lithology representative stream analysis already extracted and
arithmetically validated in `pedology/data/reference/meybeck1987_tables.json` --
the screen's durable output is a CRITICAL RUNOFF per element and per lithology,

    R_crit(X) = Q_high(X) / (c_min(X) * residence time)

which a cell either clears or does not. Concentration times runoff is Meybeck's
own model form, not a construction added here.

### What the screen returns

Critical runoff in mm per Earth year, from the report the script writes. A cell
whose local runoff exceeds its lithology's value is supplied by rock weathering
alone, with every other source counted as zero:

| lithology | K | Ca | Mg | S |
| --- | --- | --- | --- | --- |
| granite | 96.6 | 49.9 | 14.1 | 27.0 |
| peridotite | 193.1 | 38.9 | 0.9 | 7.0 |
| gneiss | 77.2 | 32.4 | 7.7 | 15.0 |
| misc metamorphic | 85.8 | 1.4 | 1.2 | 7.0 |
| volcanic rocks | 55.2 | 12.6 | 2.7 | 83.7 |
| sandstone | 36.8 | 22.1 | 6.9 | 8.8 |
| shale | 38.6 | 4.8 | 1.8 | 5.9 |
| sedimentary carbonate | 59.4 | 0.8 | 0.7 | 9.8 |
| gypsum evaporite | 22.1 | 0.3 | 0.2 | 0.1 |
| halite evaporite | 7.7 | 0.6 | 0.1 | 0.3 |

Read off it directly: **potassium is the binding element**, needing tens to
nearly two hundred millimetres of runoff a year on crystalline lithologies where
calcium and magnesium need a few. Sulfur is easy except on volcanic rock, whose
sulfate release is the lowest in the table. Magnesium is adequate almost
everywhere, and on peridotite it is adequate by three orders.

**What is not settled, and the bracket on it.** The land fraction below these
thresholds needs the runoff field, which this screen did not read. The
registered land-structure table in `notes/productivity-prediction.md` puts 14.2%
of land hyper-arid below 250 mm of precipitation and a further 15.7% semi-arid
between 250 and 500 mm, and runoff never exceeds precipitation, so **between 0
and 29.9% of land could fall below the potassium threshold.** That bracket is
the honest width; narrowing it is an evaluation, not a judgement, and the recipe
is in the next section.

**And where it is not settled, the retained sources are what close it.**
Marine aerosol at Chadwick's measured Earth rate of 600 mg Ca/m2/yr would build
the upper-bound calcium pool in about 3,200 years, two orders inside the
residence time, with no rock contribution at all. That is the reason the arid
fraction is a bracket rather than a finding of inadequacy, and it is the reason
ANUT-7 retains marine aerosol.

### What the screen refuses, and why a refusal is a result

- **Iron** is not in Meybeck's major-ion set, because dissolved iron in stream
  water is low and its supply is particulate and reductive rather than
  congruent. No artifact here carries an iron phase, so there is no lower bound
  on supply to be conservative with. The screen refuses rather than assuming.
- **The trace set** -- Mn, Zn, Cu, B, Mo, Ni, Co, Cl -- refuses for a sharper
  reason. Gaillardet, Viers and Dupre (2014) is where a per-lithology
  trace-element release table would live if one existed, and its own conclusion
  is that trace-element concentrations have been measured in a restricted number
  of environments which may not be typical. There is no table to be conservative
  WITH. Chloride is the exception in principle, being in Meybeck's set, but its
  supply over land is marine aerosol rather than rock and that deposition mass
  is undeclared.

A refusal says the element cannot be screened with what exists. That is a
different statement from adequacy and the report keeps them apart.

### The model boundary this declares

Recorded so it cannot be discovered later as an excuse: every element in this
section is outside the vendored model. Any LPJ-GUESS result on this world is a
C-N-P result and not a nutrient-limitation result while iron and the trace set
are refused and while potassium is unsettled over an arid fraction that has not
been measured. Promoting any element here to an active limiter is a separately
justified implementation task and not a consequence of this screen.

## What would settle the parts that are open

- **The per-cell potassium evaluation.** Local runoff generated at each mesh
  region is already computed, as `P - E` clamped at zero, by
  `pedology/scripts/solute_routing.py:region_runoff_m3_yr` -- and not by the
  model's `mrro`, which is river-routed net divergence and understates land
  runoff by 6.6x. The lithology per region comes from the export, and the join
  to the climate is `lib/gridding.py:climatology_cells` and nothing else. What
  it needs that this document did not use is the export payload and a
  climatology, which is a minutes-scale read.
- **The exchangeable-cation bracket.** Closed: pedology emits the field and the
  multiplier is a bound read off it, refused when it falls below what the soil
  implies. The LEVEL of the capacity relation no longer rests on one region's
  fit: it is bracketed across nine soil orders and 37,921 pedons, and
  `build_soil.py` re-evaluates this multiplier at the bracket's ends. It still
  moves this number by more than the declared bracket beside it does, so read
  the sweep in `soil_report.json`'s `level_probe` rather than the single value.
- **Everything the ledger holds open.** Twenty-one undeclared terms, each named
  above with its owner. The ledger closes when they are declared, and not before.
