# Which solve decides open water: the two lake paints, measured

Worldbuilding. Vesper is an invented planet; everything below is about a solved
lake surface on its native mesh, the depression catalogue that mesh was built
from, and the partition of a simulated cell's land those two produce.

WORLD-T8I5 asked which of two solves decides the open-water class of the wetness
partition, and said what would settle it: the two paints measured against each
other. This is that measurement.

`hydrography/scripts/surface_water.py` solves the closed-basin lake balance
twice. `paint_lakes` paints the ANNUAL equilibrium onto mesh regions and writes
`surface_water.nc:lake`. `paint_lake_cycle` runs the same fill once per time bin
of the climatology and writes `lake_cycle_bins_wet`, `lake_cycle_fraction` and
`lake_cycle_decided`. The annual area sits between the cycle's trough and its
peak, so the two paints cut the same ground differently, and the class boundary
they disagree over should be the strandline. Should be, and now is measured.

## What was measured, and on what

Measured 2026-09-05 on `canonical-10m-carve2`, from `surface_water.nc`,
`regions.nc` and `basins.nc` on that build, with land taken from
`surface_class` per CLAUDE.md rule 1. 4,328,732 land regions over 318,136,715
km2. Twelve time bins.

| quantity | share of land area | regions |
| --- | ---: | ---: |
| annual equilibrium paint, `lake` | 2.1189% | 90,058 |
| cycle, wet in every bin | 2.0673% | 87,847 |
| cycle, wet in some bins and not all | 0.1576% | 6,925 |
| depression footprint | 10.9131% | -- |
| basin cycle decided | 45.2959% | 1,960,758 |

The two paints, region by region:

| | share of land area | regions |
| --- | ---: | ---: |
| annual lake AND wet in every bin | 2.0589% | 87,489 |
| annual lake, NOT wet in every bin | 0.0600% | 2,569 |
| of which wet in some bins | 0.0599% | 2,566 |
| of which the basin never closed its year | 0.0001% | 3 |
| wet in every bin, NOT annual lake | 0.0085% | 358 |
| wet in some bins, NOT annual lake | 0.0976% | 4,359 |

## The three findings

**The disagreement is the strandline, and nothing else.** The symmetric
difference between the two paints is 0.0685% of land area, which is 0.032 of the
annual lake area, and every square kilometre of it lies INSIDE the depression
footprint. Zero regions wet in some bin lie outside the footprint, and zero
regions wet in every bin do. So the two solves disagree about where a basin's
shoreline sits within its own depression and about nothing else. That is what
licenses taking the cut from either one; it is not what chooses between them.

**The choice is forced by exclusivity, not by the sizes.** Taking open water
from the annual paint and seasonal inundation from the cycle claims the same
2,566 regions twice: they are under the annual equilibrium lake and dry in some
bin. `build_wetness.py:assign_exclusive` refuses that by construction, so the
partition has to be cut against one solve. The cycle is the one that carries the
seasonal class at all, and it is the more determined of the two states, so the
cut is taken from the cycle: wet in every bin is open water, wet in some bins is
seasonal inundation, and the depression floor neither takes is playa.

**The residual does not move, to the digit.** Under the annual cut the land
shares are 2.1189% open water, 8.7942% playa and 89.0869% dry mineral. Under the
cycle cut they are 2.0673% open water, 0.1576% seasonal inundation, 8.6882%
playa and 89.0869% dry mineral. The re-cut moves area only among the three
classes inside the depression footprint, which follows from the first finding
rather than being a second one. A consumer of `dry_mineral` is unaffected by the
choice; a consumer of `open_water` gains a seasonal class and loses 3.2% of the
area it had.

## The rule for a basin with no cycle

`paint_lake_cycle` paints nothing for a basin whose year did not close, so those
regions carry a `bins_wet` of zero that cannot be told apart from dry all year.
For them the annual equilibrium paint is the only statement that exists, and it
decides. On this build that rule moves 3 regions, and 54.7% of land regions
drain to a basin with no decided cycle while carrying no annual lake either, so
they are dry mineral soil under both cuts. The count is written into
`hydrography/analysis/wetness_report.json` as
`two_lake_paints.undecided_and_annual_lake_regions` rather than absorbed,
because it is the one place two solves meet in one partition and a rule that
moved a large area would need re-deciding rather than reporting.

## What re-measures this

`build_wetness.py` recomputes every figure in the second table on each run and
writes them into `wetness_report.json` under `two_lake_paints`, including the
share of the disagreement lying outside the depression footprint. That number is
the one to watch: it is zero here, and a non-zero value means the two solves
have started disagreeing about something other than a shoreline, which is a
different finding from this one and would reopen the choice.
