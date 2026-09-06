# What a revived saturated-area closure needs, and which of it exists

Worldbuilding. Vesper is an invented planet. What follows is about a withdrawn
closure in this project's hydrography component, about the Earth observation sets
its score was run against, and about what a second attempt would need before it
could be run.

`hydrography/config/topographic_index.yaml` carries `closure.status: withdrawn`.
`hydrography/notes/subgrid-water-table.md` section 7 measured why: the score
reaches its verdict through 34 and 37 grid cells, the gain the terrain half
carries over the cell-mean depth is at or below the scatter those cells put on
it on every one of eight arms, and a criterion restricted to fewer cells
resolves it less. The route back is stated there as SUPPORT rather than as
instrument: enough grid cells carrying observations, from more cells or from an
areal saturation or inundation observation in place of point depths.

WORLD-4FR6 put the two routes as alternatives: acquire the support, or take
WET-12's reduced form. This document is the evidence behind taking the reduced
form now while leaving the first route open, and it is written so the next
attempt starts from a list rather than from a search.

## The support exists, and it is not the thing that blocks the route

Checked 2026-09-05. `docs/src/reference/external-data.md` carries the access
routes; the papers are in `references/pdf/` and their rows are in
`references/INDEX.md`, marked held.

**The AWS Registry of Open Data does not carry this class of data.** Swept over
all 1199 dataset records and the rendered index: no occurrence of GLWD, WAD2M,
GIEMS, SWAMPS, wetland fraction, inundated area fraction or saturated fraction.
What it has is open water with the wrong footprint or the wrong record, plus one
catalogue that is land-surface-model output rather than observation. The standing
rule to check the registry first is discharged.

**Four gridded products elsewhere would supply the support.** GIEMS-
MethaneCentric v1.1 is monthly at 0.25 degree from 1992 to 2020 and is the only
one carrying inundated-and-saturated wetlands as a layer distinct from a static
permanent-open-water layer. WAD2M is monthly at 0.25 degree from 2000 to 2020.
Tootchi et al. (2019) is static at 15 arcsec and is defined as persistent
near-saturated soil from flooding OR SHALLOW GROUNDWATER, which is the union a
saturated fraction is about. GLWD v2 is a static maximum-extent inventory at 15
arcsec as percent of cell.

**The support gain is between two and three orders of magnitude, and it comes
from cells the harness already solves.** The withdrawn score reached 34 and 37
cells because a cell counts only where bores fall in it, and bores are clustered.
An areal product carries a value in every land cell, so within the SAME cached
Earth windows the harness already holds -- `data/earth_validation_cache/`
`etopo_aus.nc` and `etopo_us.nc` at 15.19 km -- a 0.25 degree product supplies of
order 10^4 land cells over the Australian window alone. That is arithmetic on a
continent of roughly 7.7e6 km2 and a 0.25 degree cell of roughly 7e2 km2 at those
latitudes, so it is a bracket of 10^3 to 10^4 rather than a count, and either end
is far above 34. No new Earth solve domain is needed.

So the DATA is not the blocker, and recording that is the point of this
document: a future reader must not conclude from the withdrawal that the support
cannot be got.

## Four things do block it, and none of them is a download

**1. The criterion has to be declared before the observations are in hand, and
nothing has been fetched.** That is what made the first withdrawal a clean
outcome rather than an argument, and it is WORLD-4FR6's condition on this route.
Fetching any of the four products before a criterion exists spends the one thing
that would license the next score. Nothing has been fetched. The papers have
been, because reading them is how a criterion gets written.

**2. An areal observation is a DIFFERENT INSTRUMENT, not more support in the old
one.** The withdrawn score is a per-bore discrimination: it asks whether a
cell-scale predictor ranks a shallow bore above a deep one, and reports an area
under the curve. An areal product replaces that with a per-cell predicted
fraction against a per-cell observed fraction, which is a skill statistic on a
continuous quantity and has no area under a curve in it. So the declared bar,
the comparator and the attribution identity all have to be written again for the
new statistic. A criterion copied across from the bore score would be a
criterion for a question the new data cannot be asked.

**3. The predictor half does not exist on any build.** The closure multiplies a
terrain rank statistic by a function of the water table depth, and no build under
`hydrography/data/` carries `water_table.nc`. The `groundwater` step is
`cost: minutes` and needs `surface_water` and a baseline climatology first.

**4. The depth field needs a regime, which is GW-24's and MIN-6's question.** It
scores below the 0.5 of no discrimination on one of the two Earth bore sets and
above it on the other, so a cell mixing the regime where lateral flow sets the
head with the regime where the local sink does has no verdict attached to it.
`groundwater_access_<grid>.nc` already reports the `sink_fraction` regime share
per cell, including its two kinds of no-answer.

## What the reduced form does with this

The reduced form is taken now and does not close this route.
`biosphere/notes/reduced-wetland-form.md` declares the saturated non-inundated
mineral class absent with its reason and its licence, and states the cost in
claims; `biosphere/config/wetlands.yaml` carries the declaration and
`biosphere/scripts/wetland_gate.py` reads hydrography's closure status on every
invocation, so a revival refuses in the consumer until what it may key has been
decided again.

One constraint binds both routes and is declared before any share is written:
`extent.convention_arm_rule`. A class whose share is a convention bracket
propagates through every downstream ledger as both arms or as neither. The case
it exists for is `f_grad`, whose upper arm makes the closure's ranking the
depth's ranking exactly, so a consumer taking that arm alone is reporting the
depth under another name.
