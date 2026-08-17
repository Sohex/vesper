# Validating the lake solver against Earth

HYD-4. Registered 2026-08-17, **before any data was fetched**, in the same way
`pedology/scripts/validate_against_earth.py` and the albedo bracket's criteria
were. What follows above the line is the test; what follows below is what it
found.

## Why this test and not another

The lake solver is the least-checked product in this pipeline, and it is load
bearing twice over. It sets how much of the land is under water, which is worth
-0.0139 on land-mean albedo and therefore reaches the climate; and it is the only
thing that could break the degeneracy between salt crust extent and salt crust
albedo, which `exoplasim/notes/parameter-decisions.md` records as unresolvable
from the climate side alone.

Penman was validated the same way and that is the standard being matched here:
applied to ocean cells, which already are open water, it reproduced the model's
own evaporation to 1.7%. The lake solver has had nothing equivalent.

## What is under test, and what is not

The solver has three parts and only one of them is testable on Earth.

**Testable: the equilibrium relation.** A terminal lake grows until net
evaporation over its surface consumes what its catchment delivers:

    R * (A_catch - A_lake) = (E_lake - P_lake) * A_lake

which rearranges to a prediction that needs no hypsometry at all:

    A_lake / A_catch = R / (R + E_lake - P_lake)

**Not testable here: the hypsometric inversion**, which converts a solved area to
a level and a volume through curves measured on this world's terrain. Earth
basins have their own hypsometry and testing it would test the data, not the
code.

**Not testable here: the overflow cascade**, for the same reason. What can be
said about it is already known -- it reports whether it reached a fixed point,
and a non-converged solution is refused rather than used.

## The prediction, and the thresholds

Stated before fetching anything.

The relation above is scale free, so the test is of the *ratio* rather than of an
area. Earth's terminal lakes span four orders of magnitude in size, so a relation
that is right in form should hold across all of them or fail visibly.

- **PASS** if the median absolute log-ratio of predicted to observed lake area is
  below 0.30, i.e. predictions typically within a factor of 2.
- **MARGINAL** between 0.30 and 0.60, a factor of 2 to 4.
- **FAIL** above 0.60, or if the error correlates with basin size, aridity or
  lake area, since any of those would mean the relation is missing a term rather
  than being noisy.

A factor of 2 sounds loose and is not. The quantity that reaches the climate is
lake *area* as a fraction of land, and the albedo consequence of getting it wrong
by a factor of 2 is smaller than the playa albedo uncertainty already carried at
1.99 W/m2. What would matter is an order of magnitude, or a bias.

## What will make this noisy, stated in advance

**Earth's terminal lakes are mostly not in equilibrium**, and the ones people
have measured best are the ones humans have most altered. The Aral Sea and Lake
Chad are the two most cited endorheic basins on Earth and both are dominated by
twentieth-century irrigation withdrawal; their present areas measure water
management, not climate. Where a pre-development balance is published, that is
the number to use, and the choice has to be made per basin and stated.

**Lake area is not a smooth function of climate.** A basin's hypsometry can be
flat, so a shallow terminal lake swings enormously in area for a small change in
volume. Lake Eyre is the extreme case, dry in most years and 9,500 km2 in a wet
one. A single-year area is the wrong observable; a long-term mean is the right
one and is not always available.

**Runoff over the catchment is the least certain input**, because for most of
these basins it is itself modelled rather than gauged.

None of that is a reason to skip the test. It is a reason to report the scatter
and to say which basins carry it.

## Method

1. Choose basins on whether a defensible water balance exists, before looking at
   whether the prediction works for them.
2. Take catchment area, lake area, catchment runoff, lake evaporation and lake
   precipitation from the literature, one source per basin where possible.
3. Compute the predicted ratio from the relation above and compare.
4. Report every basin, including the ones that miss.

---

## Results

**Not run. The test is well posed and the data it needs is not in the sources
that describe these basins.** Recorded on 2026-08-17, because the reason is a
property of the literature rather than a gap in effort, and the next attempt
should start from the right kind of source instead of repeating this one.

### What the test needs, restated after trying

    A_lake / A_catch = R / (R + E_lake - P_lake)

`R` has to be **runoff depth over the catchment, measured independently of the
lake**. That independence is the whole test. Everything else can be had.

### What the reviews actually carry

Two were fetched and read.

Yapiyev et al. (2017) is the right kind of paper and gives the terms for two
basins. For Great Salt Lake: catchment 55,000 km2, historically more than 89,000;
lake 2,470-5,490 km2; lake precipitation 370 mm/yr and lake evaporation 1,000
mm/yr; inflow composed of 66% river discharge and 31% direct precipitation. For
Issyk-Kul: basin about 40,000 km2, and a balance of roughly 300 mm/yr each from
direct precipitation, surface runoff and groundwater against 800 mm/yr of lake
evaporation and under 100 mm/yr of irrigation withdrawal.

**Every one of those terms is expressed per unit LAKE area, and that is what
makes them useless here.** A balance written that way closes on itself: Issyk-Kul
gives 300 + 300 = 600 in against 500 + 100 out, which is exact, and confirms
nothing except that the authors' numbers are consistent. To predict an area, the
supply term has to be a depth over the *catchment*, which none of these state.

Great Salt Lake shows the circularity plainly. Its river inflow can be recovered
from the balance -- (E - P) x A_lake, or about 2.77 km3/yr at the mid-range area
-- and that is the same equation being tested, so it can only ever agree with
itself. The paper's own numbers are also mutually inconsistent at the 19% level:
the 31% precipitation share implies a lake evaporation of 1,193 mm/yr against the
1,000 stated, which is what a century of disequilibrium looks like in a table.

Wurtsbaugh et al. (2017) supplementary gives annual volume series for the Aral
Sea, Great Salt Lake, Owens, Urmia, Walker and the Dead Sea. It is the best
available evidence for the disequilibrium warned about above -- the Aral falls
from 100% to 4.5% of maximum volume within the record -- and carries no water
balance at all.

### What would actually work

Gauged discharge per basin, in km3/yr, against a stated catchment area. That
exists, in per-basin hydrology papers and in discharge archives, but it is one
source per basin rather than one source for the set, so the assembly is the work
and the definitions have to be reconciled by hand: basin against catchment
against active catchment, and lake area at which epoch.

The alternative is gridded: basin polygons from HydroBASINS crossed with a
gridded runoff product, area-averaged properly. That gives `R` over the catchment
by construction and would let the test run over dozens of basins instead of a
handful. `pedology/scripts/validate_against_earth.py` already establishes the
pattern for reaching a live dataset from inside this project, so the machinery is
not novel; the data volume is the cost.

**Do not accept a two-basin version of this test.** With N of 2, definitional
ambiguity in the catchment area alone moves the answer by more than the pass
threshold, and a result that cannot distinguish the model from the bookkeeping is
worse than none.
