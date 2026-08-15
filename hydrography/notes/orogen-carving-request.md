# Request to the World Orogen fork: per-basin carving control

Two items from building the hydrography component against the 2026-08 export
(`finalElevation` `821aa71b37a7...`). The first is a feature we are blocked on;
the second is a correctness issue in an existing product.

## 1. We need a way to say "carve these specific basins"

### Why

The fork's position is that Orogen measures basin geometry and the water balance
decides endorheism. We can now compute that verdict, but there is nowhere to put
it.

A basin that overflows year on year incises its outlet. Over the 1e4 to 1e6 years
a landscape takes to relax, that cuts the sill, drains the lake, and the
depression stops existing. So the terrain as exported is not in equilibrium with
any climate for those basins: hydraulic, thermal and glacial erosion all acted on
it, but drainage-enforcement carving was suppressed globally so basins would
survive, leaving rims standing that the outflow crossing them would have cut.

The overflow test is pure geometry up to one climate number. A basin overflows
exactly when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

Across the 3,629 preserved basins the right-hand side has a median of 2.92, and
the resulting endorheic land share runs from 9% of land in a humid climate to 76%
in a hyper-arid one. The figure is meaningless without the climate, and the
terrain should differ between those cases.

### What we need

A per-basin subtractive control: basins named should be dropped from the
preserved set and carved normally, as `--no-basins` would do to all of them.

**It must take a file, not repeated flags.** The list is 10^2 to 10^3 entries:
2,583 basins carve at a ratio of 2, 857 at 5. Something like

    --carve-basins FILE      # newline-separated basin ids, carve these
    --preserve-basins FILE   # or the inverse: preserve only these

Either shape works. The inverse form may be cleaner since it makes the preserved
set explicit and reproducible rather than a function of the size floors plus a
subtraction.

### The hard part: id stability across the iteration

This is the design question we cannot answer from outside.

The loop is: export -> we compute the verdict -> you regenerate with the carve
list -> we re-derive -> repeat. It has to converge, and carving is monotone
(carving only removes basins), so it should. But it only works if the ids we send
back mean the same thing in the next export.

Ids are `orogen-basin-r<sink>-<hash>`. Carving changes the terrain, so surviving
basins may re-measure differently and take new ids, which would break the loop on
the first iteration. The README already warns ids are not stable across region
counts, for the same underlying reason.

We need either:

- ids stable under carving for basins that were not carved, or
- a mapping in the manifest from each input id to the basin it became, or to
  nothing if it was carved away.

The second is probably more honest, since carving a basin can genuinely merge or
split its neighbours.

### Please also record what was done

`manifest.basins` should echo the carve list it was given and note which entries
matched, so an export says on its face which drainage hypothesis produced it. Two
exports of the same seed will now differ in terrain, and the hash alone will not
say why.

### Optional refinement, not a requirement

Binary carve-or-preserve is coarse. Physically, incision depth should scale with
overflow discharge and inversely with sill erodibility, both of which you already
have. A basin that barely overflows should end up a through-flowing valley with a
residual lake, not a fully cut trench. If accepting a per-basin overflow discharge
and letting the erosion model do proportional incision is not much more work than
the binary switch, it would be substantially more physical. If it is more work,
the binary switch is enough for now.

## 2. `finalCatchment.areaKm2` understates catchments by about 2x

`manifest.basins.preserved[].finalCatchment.areaKm2` is documented as "the area to
integrate precipitation over". It is computed on raw steepest descent, so it
counts only regions whose unrouted flow path reaches the sink.

On this planet `drainage_terminal` is -2 for 63% of the land, which drains into
220,649 unpreserved pits, mostly single mesh cells of noise. All of that land is
excluded from `finalCatchment`. Totals: 7.12e7 km2 across all basins by the
catalogue, against 2.417e8 km2 once the pits are filled and drainage resolved,
with a median per-basin ratio of 1.955.

We reproduce the catalogue's figures exactly by selecting `drainage_terminal ==
sink`, so this is a definition difference rather than a discrepancy. But anyone
following the manifest's own advice integrates precipitation over roughly half
the real catchment.

The catchment-to-floor ratio is the clearest symptom. On the catalogue's numbers
the first basin has a catchment 0.13x its floor area, which would make a lake
impossible. With drainage resolved the median ratio across all basins is 7.6x,
which is the concentration that lets lakes exist at all.

Suggested fix: either document that the field is measured on unrouted flow and is
a lower bound, or fill depressions before measuring it. We do the filling
ourselves in `hydrography/scripts/drainage.py` (priority flood, about 4 seconds
over the 2.5M-region mesh), so we do not need the field changed, but it is
currently a trap for anyone who trusts the note.

## What we will send you

Once the hook exists: a newline-separated list of basin ids to carve, plus a JSON
sidecar giving each one its `critical_aridity_index`, catchment area, area at
spill, and the modelled `(E - P) / runoff` that produced the verdict, so the
decision is auditable rather than a bare list.

---

# Follow-up on the 2026-08 regenerated export (`27b7479a...`)

## Verified working

- **Routing coverage.** `drainage_terminal == -2` is 0 on land, from 63%.
- **Catalogue stability.** `hashes.basinCatalogue` is present, counts are
  unchanged (81,904 detected, 3,629 preserved, 922/6,880 below sea level), and
  `selectionSource`, `drainageHypothesis` and `incision` are all populated.
  This export carries an empty hypothesis, so it is the all-preserved endmember.
- **`finalCatchment.areaKm2` is fixed.** It now totals 2.412e8 km2 against an
  independent priority flood's 2.409e8, with a per-basin ratio of median 1.0000,
  p5 0.999, p95 1.000. Median catchment/floor is 7.46x against our 7.46x. Two
  independent implementations agreeing to 0.1% per basin is about as good a
  cross-check as this gets.

## One new bug: `drainage_terminal` contradicts `finalCatchment`

The two products in the same manifest no longer agree.

| Measure | Endorheic land |
| --- | ---: |
| implied by summing `drainage_terminal` per sink | 1.846e8 km2, 58.19% |
| declared by `finalCatchment.areaKm2` | 2.412e8 km2, 76.01% |

878 of 3,629 basins differ by more than 1%. For **every one of those 878**,
`drainage_terminal` routes *nothing* to the basin at all: implied catchment under
1 km2, while `finalCatchment` declares a full one. Together they account for
5.656e7 km2, 17.82% of the planet's land.

### The cause is single and clean

For exactly those 878 basins, the sink's own `drainage_terminal` is not itself.

- 2,751 of 3,629 sinks are self-rooted.
- The 878 that are not are precisely the disagreeing set, no overlap either way.
- Of the 878, 874 point at ocean-draining land and 4 at another basin sink.

So those sinks were never roots of the spanning tree. The traversal reached them
from outside and gave them a parent, after which nothing can accumulate in them.

The affected set has a systematic signature: median sink elevation 0.026 km
against 0.001 km for all basins, median spill 0.114 km against 0.062 km. These
are the shallower, higher-sitting basins, the ones an outward-growing front from
the ocean reaches most easily.

### Suggested fix

Either make basin sinks roots unconditionally, seeded before traversal and never
eligible to receive a parent, or make the traversal a priority flood keyed by
level rather than a plain spanning tree.

The level-keyed version is what keeps the two fronts apart without special
casing. Seed ocean-margin land and basin sinks each at their own elevation and
pop lowest-first: a basin sink sits below its rim, so the basin front claims the
interior long before the ocean front can cross the rim from outside, and the
meeting point falls out as the spill. Our `hydrography/scripts/drainage.py` does
this in about 4 seconds over the 2.5M-region mesh if the implementation is
useful as a reference.

`finalCatchment` is the one that is right. It matches an independent
implementation to 0.1%; `drainage_terminal` does not.

## Please confirm: the terrain changed with no carve applied

`finalElevation` went `821aa71b` to `27b7479a` while `drainageHypothesis` is
empty, and mean land elevation fell from 548 m to 521 m. We assume that is the
routing fix feeding hydraulic erosion, which would be a genuine improvement since
erosion was previously routing over unfilled pits. Worth confirming it is
intended, because with an empty hypothesis nothing else on the artifact signals
that the terrain moved, and every number computed against the previous build is
now stale.
