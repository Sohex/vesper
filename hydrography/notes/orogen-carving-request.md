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

---

# Follow-up on the re-regenerated export (back to `821aa71b`)

## Verified

- **The revert is clean.** Every hydrography product reproduces the earlier
  `821aa71b` run to the digit: endorheic share 76.1%, 99,085 pits filled at a
  median 8.5 m, capacity 6,838,339 km3, critical aridity index median 2.921. The
  terrain really is bit-identical.
- **`drainage_terminal` and `finalCatchment` now agree exactly.** Independent
  check: max per-basin `|ratio - 1|` is 0.00e+00 and both total
  2.418487e+08 km2 to every digit. The `drainageConsistency` block is doing what
  it says.
- **Area conservation is exact.** Land draining to sinks, plus land draining to
  ocean, plus non-land, sums to 7.349930e+08 km2 against a planet surface area of
  7.349930e+08, relative error 0.00e+00.
- **River statistics reproduce.** Sink inflow median 3.189e4 km2 against the
  reported 3.19e4; inflow-to-floor ratio 7.64x against 7.64x.
- **Our two routings now agree to 99.888% of land**, 1,207 regions of 1,079,435,
  with a per-basin catchment ratio of median 1.0000 and p5/p95 both 1.000. The
  residual is all one direction, ours assigning to ocean where the export assigns
  to a basin, which is our documented tie-break: we seed the heap so that a
  saddle reached at equal level goes to the ocean. Neither is wrong; it is a
  choice at exact ties.

We keep our own flood rather than reading `drainage_terminal`, because it also
produces the filled surface the hypsometry is built from, and because an
independent implementation is what surfaced the last two routing bugs. The
comparison is now a standing regression check recorded in
`hydrography/data/<build>/hydrography_report.json` under
`cross_check_vs_export_routing`.

## One correction, on the reported percentage rather than the routing

`drainageConsistency` reports 79.74% of land endorheic. That figure is normalised
by `land_mask`, not by `surface_class`.

    catchment total          2.418487e+08 km2
    / land by surface_class  3.173228e+08 km2  ->  76.22%
    / land by land_mask      3.032877e+08 km2  ->  79.74%

`land_mask` excludes the dry sub-sea-level basin floors, 1.91% of the planet,
which is the trap the manifest itself now documents and which `surface_class`
exists to avoid. Using it here inflates the endorheic share by about 3.5 points.

The routing is right and the areas are right; only the normalisation is wrong.
Worth fixing because a headline percentage is exactly the number that gets
quoted downstream.

## And a correction against us

Our guess that the previous build's terrain move came from the routing fix
feeding hydraulic erosion was wrong, and wrong in a way we could have checked:
the ordering in `pipeline.js` settles it, with post-processing ending at 462 and
all basin and river work at 504 and beyond, so routing cannot feed erosion. We
offered a plausible mechanism without checking whether it was reachable, and
labelled it a probable improvement. The actual cause, a pre-erosion absolute
protection floor that also served as its own assertion baseline, was a
regression. Noted so the reasoning error is on the record alongside the fix.

## A published fraction with the same denominator problem

Following up on the endorheic percentage: the artifact itself is clean there.
`drainageConsistency` publishes absolute areas only, so the 79.74% exists in
prose rather than in the manifest and nothing downstream can pick it up.

`lithology.compositionLand` is a different matter, because it is in the artifact
and it looks like data to build on. It is measured against `land_mask` in both
numerator and denominator: its areas sum to 3.032877e+08 km2, which is the
`land_mask` land area exactly, against 3.173228e+08 for `surface_class`.

The distortion is not uniform, and it lands where it hurts most:

| class | published | vs surface_class | area gain |
| --- | ---: | ---: | ---: |
| evaporite | 0.1860 | 0.2083 | 1.17x |
| oib | 0.0129 | 0.0133 | 1.08x |
| pelagic | 0.0133 | 0.0136 | 1.07x |
| (land overall) | | | 1.046x |

Evaporite is the most affected class, which is the physically expected direction:
playa and salt-pan fill accumulates in exactly the closed basins below sea level
that `land_mask` drops. So the table systematically under-reports the lithology
most characteristic of a planet with preserved endorheic drainage.

Nothing we have built consumes it. Our albedo boundary condition is computed from
the `rock_albedo` field cell by cell, not from the table, and every land fraction
in `hydrography/` is computed against `surface_class`. The measurable consequence
is small: land-mean bare-rock albedo is 0.3083 by the table's denominator against
0.3135 by `surface_class`, about 0.7 W/m2.

But we did quote 18.6% from it in two places before checking, and a biome or soil
model taking parent-material fractions from this table would inherit the same
bias. Suggest measuring it against `surface_class`, or naming the denominator in
the block the way `landSeaMask` already does.

---

# Follow-up: denominator corrections verified

All confirmed against the regenerated export.

- `elevHash` 821aa71b and `catalogueHash` 2d1f8e57 both unchanged, and every
  hydrography product reproduces to the digit: endorheic 76.1%, 99,085 pits at a
  median 8.5 m, capacity 6,838,339 km3. Manifest numbers moved, terrain did not.
- `drainageConsistency.fractionOfLand` is 0.7621534408214234, which matches what
  the export's own `drainage_terminal` field gives to every digit, and
  `landAreaKm2` matches `surface_class` land area exactly. `build_hydrography.py`
  now asserts the published fraction against the field rather than recomputing
  it, so a future denominator change on either side fails the build.
- `compositionLand` reproduces our figures exactly: evaporite 0.2083 on
  6.608997e+07 km2, total 3.173228e+08, fractions summing to 1.000000.
- `computeScarpPotential`: verified zero. Not one of the 48,092 dry-floor regions
  on the native mesh has nonzero scarp potential, which is a stronger statement
  than the subset measured upstream. Agreed it is a closed trap rather than a
  finding.

## One measurement on erodibility, offered as a note

The reasoning for keeping `buildErodibility` on `elevation > 0` is sound and we
are not disputing it: when it runs, no basin has been preserved yet, so nothing
below sea level is land and the two definitions coincide at that moment. The
divergence only exists downstream.

The measurement is still worth recording, because the delivered field no longer
has the property the docs claim for it. `tools/README.md` describes erodibility as
"mean-normalised to 1 over land". On the finished export:

    unweighted mean over land_mask (elevation > 0)   0.9735
    unweighted mean over surface_class (subaerial)   0.9995

Neither is 1, which is expected once cover stripping swaps cells to their basement
erodibility after the normalisation was applied. So the normalisation is a
property of the field at the time it was computed, not of the field as shipped.
Worth a clause in the docs, since a downstream model that assumes mean 1 and
rescales against it would introduce a few percent of error for no reason.

## On the two implementations

Your framing is the right one and worth keeping: a consistency check between two
products derived from the same root cannot find a bug in the root. That is
structural, not a matter of care. It is the argument for keeping an independent
implementation even after both agree, which is why the cross-check now runs every
build rather than having been a one-off comparison.

---

# First carve verdict, from the 0.96 S-Earth baseline

Climate: T42, vegetated land surface, glaciers enabled, K2 spectrum, converged on
all six criteria at 292.88 K. Five-orbit climatology, integrated over each
basin's catchment through the T42 coupling matrix.

## Verdict

| estimate | carve | survive | dry | with lake | lake, % of planet |
| --- | ---: | ---: | ---: | ---: | ---: |
| Penman (primary) | 1522 | 2107 | 1248 | 962 | 0.945% |
| land-evaporation sensitivity | 2119 | 1510 | 1248 | 967 | 1.479% |

**1,522 of 3,629 basins carve** under the primary estimate, 1,510 survive under
both, and 597 (16.5%) are sensitive to which evaporation estimate is used.

Note that "carve under both estimates" is no longer independent evidence. Penman
is floored at the land rate, so it exceeds it by construction and the Penman
carve set is a strict subset of the land-evaporation one. The two numbers are
identically 1,522. A genuinely independent second estimate would be needed to
recover that check.

## How open-water evaporation was estimated, and why it can be trusted

The model has no lake in these basins, so it does not report lake evaporation.
Land evaporation is moisture-limited: ExoPlaSim scales it by a wetness factor
that only reaches 1 above 40% of field capacity, so it understates what a lake
would evaporate and over-carves.

The primary estimate is the Penman combination equation, evaluated with water's
albedo (0.06 rather than the 0.15 to 0.50 substrate) and water's roughness
length. Net radiation is recomputed by backing shortwave out of the model's `rss`
using the albedo the run was actually given, then re-absorbing it at water's
albedo.

**It is validated against the model itself.** Applied to ocean cells, which are
already open water and which the model gives water's own roughness, Penman
reproduces the model's own evaporation to within a few percent -- the figure is
computed on every run and written into `carve_verdict.json` rather than quoted,
having twice been quoted from a hardcoded constant. Reproducing open-water
evaporation over open water from surface fields alone is what makes the same
calculation trustworthy over a lake the model does not have.

A third approach, dividing land evaporation by the reconstructed wetness factor,
was tried and rejected: it gives a land mean of 20.5 mm/day, because the division
is unstable wherever soil is dry, which is precisely where endorheic basins live.
Recorded so it is not attempted again.

## The dry basins are the interesting result

1,248 basins have **zero catchment runoff**. That is 34% of all basins and it agrees
independently with the lithology: evaporite is 20.8% of this planet's land, and
evaporite forms in closed basins. Two separate models, terrain chemistry and
climate, arriving at the same picture.

The two classes of surviving endorheic terrain overlap rather than split:
1,248 basins with zero catchment runoff, 962 holding a real lake above 1 km2
(together about 1% of the planet's surface), and about a hundred basins in
both -- zero-runoff floors whose lake is fed by rain on its own surface. So
zero catchment runoff does not by itself rule out standing water or overflow.

## What we are not sending yet

The 597 estimate-sensitive basins are one reason to hold. The larger one is that
this verdict inherits an assumed biosphere: the vegetated land surface was chosen,
not modelled, and the albedo bracket showed that choice is worth 3.7 to 7.1 K,
which flows straight into precipitation, evaporation and runoff here.

The 1,522 basins that carve under both estimates are robust to the evaporation
question, though not to the vegetation one. They are a defensible first list if
the terrain iteration should start now, with the remainder settled on the second
pass once LPJ-GUESS has replaced the assumption.
