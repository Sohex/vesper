# Iteration 1 — superseded drainage verdict, corrected lithology

**Do not consume this export's basin verdict as a finished answer.** The terrain
is a valid Orogen run and internally consistent, and its lithology and albedo are
current. The *drainage hypothesis* it was built from is known to be wrong, and
wrong in two opposing ways at once.

- terrain hash: `3899a0c57d1eee2f…` (`hashes.finalElevation`)
- basin catalogue hash: `2d1f8e57c26b60c6…` — **unchanged** from `out/vesper/`
  and from every earlier export, so every per-basin artifact keyed by id still
  resolves against this one
- verdict applied: `/home/cfutro/docs/world/hydrography/data/carve_list.txt`
  (3,629 entries: 1,872 preserve, 235 partial, 1,522 carve)

## What is wrong with the verdict

**Under-carved, ~228 basins.** The verdict preserved 1,248 basins on the grounds
of zero catchment runoff. Zero *catchment* runoff says nothing about
precipitation falling directly on the lake surface: where P > E over open water
at the basin floor, a lake gains water regardless of catchment supply and will
fill and spill. Of the 1,248, 1,020 are genuinely dry pans (E > P) and **228
should have carved** — 6.3% of the catalogue, ~15% of what this export preserves
outright. Robust to sampling: evaluating the forcing at each sink cell rather
than over the catchment gives 235, with only 9 of the set flipping.

Those 228 should probably not carve at retain 0 either. They spill on lake-surface
P − E alone with zero catchment inflow, which makes them the lowest-discharge
spillers on the planet, and incision scales with discharge.

**Over-carved, the marginal band.** The verdict came from a climate run on the
*pre-carve* planet. Carving removes closed depressions, which removes their fill:
this export drops closed-basin fill from 20.9% to 12.7% of land. Darker surface,
warmer world, higher evaporation — and marginal basins this verdict carved would
have stayed closed under the resulting climate. That overshoot is not yet measured.

## What has been corrected since the verdict

The basin surface is now **zoned**. It used to be a single class, "Evaporite /
playa fill", carrying clean-halite albedo 0.50 across every cell of every closed
basin. It is now salt crust (`evaporite`, 0.50) in the sump — the lowest quarter
of each basin's relief, 10.6% of basin fill by area — and `playa_clastic` (0.30)
over the rest, for an **effective 0.321** where the verdict's climate assumed 0.50
everywhere.

That is a large cooling correction and it partly offsets the carving feedback
above. Both belong in iteration 2. It also changes erodibility (playa 2.80 against
salt 3.50), which is why the terrain hash moved from `407cdd8586d43f4b`.

`playa_clastic` is a **new rock class, id 19**. Readers that hardcode the class
table rather than reading `manifest.lithology.rockClasses` will not see it.

## Why the verdict was not re-run

Fixing the 228 in isolation would produce an export still wrong in the other
direction, from the same superseded climate. Iteration 2 re-runs the baseline on
*this* surface — which is the point, since the albedo it now carries is the input
the next verdict needs — and fixes both errors together.

## What is good here

The terrain, mesh, lithology, albedo, routing and every field in `raw/` and the
grids are valid outputs of the stated hypothesis. `manifest.basins` records that
hypothesis (`selectionSource: preserve-list`, `drainageHypothesis.entries`,
`carvedByRetainZero: 1522`). Use this export as the climate input for iteration 2.
Do not use its basin verdict as the answer.

`out/vesper/` is the pre-carve planet (terrain `26fc76914da14289…`), also carrying
the zoned lithology. Its basins come from the threshold default rather than a
climate verdict, so it carries no drainage error — only the same albedo correction.
