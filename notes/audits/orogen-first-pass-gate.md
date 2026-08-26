# What is outstanding against the first Orogen export

*Swept 2026-08-26 against 224 open `bd` rows, the three Orogen audits, and the
working tree. This is a WORLDBUILDING project: the subject throughout is the
terrain generator for a simulated planet and the artifact it writes.*

`config/pipeline.yaml`'s `orogen` step carries a gate and refuses to compute it:

> the carve gate: nothing outstanding may still move an artifact the verdict is
> computed from. READ THE OPEN ISSUES AND DECIDE. Deliberately not computed, and
> not because a filter over `step:` labels performed badly -- because it answers
> a different question.

This is that reading, recorded so the next cycle re-takes it rather than
re-deriving it. It answers one question and not its neighbours: **does anything
outstanding change the artifact Orogen writes on the FIRST pass**, before any
loop input exists.

## The verdict

**No.** Nothing open changes the first export.

## What the first pass is a function of

The first pass through loop A cannot satisfy either of the `orogen` step's
`needs`. `carve_list` comes from a basin catalogue that does not exist yet and
`ice_mask` comes from a climatology that does not exist yet, which is why the
step declares that edge a PASS BOUNDARY rather than a cycle. So the first export
is determined by exactly four things:

- the planet code and the seed it decodes to,
- the three parameters the code does not carry and the recipe passes
  explicitly: radius, gravity and lithology strength, checked against
  `config/planet.yaml` per rule 2,
- the region count, and
- `vendor/orogen/js/` itself.

An outstanding item moves the first export only by moving one of those. Every
other row in the tracker reads the export, or reads something derived from it.

## How the sweep was taken

Three passes, because a tracker sweep alone misses a finding that never became a
row, and a label sweep answers a different question than the gate asks.

1. All 224 open rows read by title. The 17 whose title, description or notes
   name the generator, `vendor/orogen`, or one of its modules by filename read
   in full.
2. The three that survived as genuine candidates read end to end: `LITH-26`,
   `CONV-5`, `CLIM-63`.
3. `notes/audits/orogen-gravity.md`, `orogen-lithology.md` and
   `orogen-resolution.md` swept for findings carrying no row.

## Why each candidate does not reach the first pass

**`LITH-26`**, the shelf substrate classification, is the row the question is
usually about. Its replacement is driven by accepted sea-surface and column
state, carbonate chemistry, biological carbonate production and `ANUT-6`
terrigenous sediment supply. None of those exists before a climatology, so there
is nothing for a first pass to read even in principle, and the row's own scope
keeps the latitude map as a labelled bootstrap rather than removing it. Its
landed half, `SHELF_SUBSTRATE_BOOTSTRAP` at commit 35b6ae9e, changed no cell's
class and moved no terrain.

The row's block list is circular if read literally -- it requires a declared
canonical climatology, which requires a commissioned baseline, which requires a
build, which this row changes. Loop A resolves inputs of this class by
ITERATION, which is what the `orogen` step's `first_pass` note already records
for `ice_mask`. That circularity is a fact about which PASS the row lands in,
not about whether the first one is blocked.

**`CONV-5`** excludes the export in its own scope: a metadata edit under
`source/` is an overwrite rather than a new build, which rule 7 forbids.

**`CLIM-63`**, the SEMI surface mass balance candidate, feeds the ice mask,
which is a second-pass input by construction.

**`OCN-11`** builds its support and bathymetry contract FROM the accepted
export.

**`SPAT-5`** is the coastline fraction at the ATMOSPHERE boundary, which is the
`boundary_conditions` consumer. `notes/audits/coastline-threshold-cost.md`
leaves the fractional-boundary question open and points here; the 0.5 threshold
it settles is applied when the mesh is reduced to a Gaussian grid, after the
export exists.

**`SURF-4`, `SURF-7`, `MIN-6`** place surficial and supergene materials in
`pedology/` and `minerals/`.

**`BIO-28`, `GRAV-7`, `WET-7`, `WORLD-VEG9`** are downstream modules that need
to CONSUME this planet's gravity. None changes the declared value, which Orogen
is canonical for.

## The two audit findings with no row, and why neither is one

Both are in `notes/audits/orogen-lithology.md` and both resolve there as
intended behaviour rather than defects. They are recorded here because they are
what a sweep of the tracker alone would miss.

**The plutonic excess.** Granodiorite reaches 7.29 per cent of land against
GLiM's 0.4 per cent intermediate and 5.7 per cent acid plutonic, an order of
magnitude more exposed arc-root batholith than Earth carries. Two mechanisms
were proposed and MEASUREMENT refuted both: builds at 1.0, 2.0 and 3.0 km of arc
cover move granodiorite only from 7.29 to 7.10 per cent, because arc-belt
erosion is bimodal -- 94.1 per cent of granodiorite cells eroded their cover to
zero while cells that keep cover retain about half of it, and no cover thickness
sits inside that gap. The arc root is exposed because arc terrain erodes deeply,
and there is that much arc terrain because this world was given 100 plates. Both
are decisions already taken.

**The craton archetype fires on 0.75 per cent of land.** Not a regime mismatch:
granite runs 11.1 per cent at both 24 and 100 plates, so plate count does not
move the two classes the finding is about. The branch fires, responds to its
parameters and takes basement from the granite fallback as intended; what it
cannot do is reach the SURFACE, because the quiet ground it claims carries the
thickest cover, playa fill at 35.7 per cent and evaporite at 15.2 per cent.
That is the shield-versus-platform distinction, and Earth's shields are exposed
by glacial scour and epeirogenic uplift, which Orogen models neither of.

## What HAS changed, and is incorporated rather than outstanding

The generator has moved since `precarve-craton-10m`, so the next build differs
from it. Every one of these is a decision already landed:

| change | commit or row |
| --- | --- |
| `BASIN_MIN_AREA_KM2` 1000 to 850 | bb71cbf0, WORLD-BVL9 |
| the depth floor compared in physical km rather than model units | WORLD-YRIL |
| a relief-curve branch above the shape function's domain | 29963a28 |
| the export says the height curve does not stop at 1 | 35975c2b |
| `--glacial 0` rather than the code's slider | PHYS-13 |

## What this sweep does NOT say

It says nothing about whether the build survives its own loop. `LITH-26` still
moves terrain when it lands, and it lands on a later iteration; the same is true
of the ice mask and the carve list. A first export that nothing outstanding
changes is not a terminal one, and the gate does not ask for a terminal one.
