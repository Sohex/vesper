# What the binary coastline threshold costs, in both signs, on every rung

**Measured:** 2026-08-28, on `canonical-10m-carve2`, terrain hash
`f496ae9fd749`,
10,000,005 regions, at every rung of the T21/T42/T85/T127/T170 ladder.
`analysis/coastline_threshold_cost.py` reruns
`exoplasim/scripts/build_boundary_conditions.py:coastline_ledger` on every
rung. The machine-readable measurement is
`analysis/coastline_threshold_cost.json`; each ordinary boundary build also
reproduces its active rung in `boundary_conditions_report.json`.

This is worldbuilding. Vesper is an invented super-Earth; every quantity below
is a property of that planet's modelled surface or of the gridding this
pipeline applies to it. Nothing was simulated: no World Orogen generation, no
ExoPlaSim run.

SPAT-5's first clause asks for the binary coastline threshold's area, store and
flux cost on every ladder rung. This is the area and the store; the flux
follows from the area once a climatology exists and is not derivable here.

## The threshold, and why it is not a choice about the model

`config/planet.yaml:geography_land_threshold` is 0.5 and
`docs/src/reference/config-rationale.md` records it as DETERMINED, with its
cost stated as a coastline effect and not sized. The model's mask is binary
whatever this builder writes: `oceanmod.f90:256` hard-binarises `yls` at 0.5,
so there is no tiling and no partial-water cell to write into. The threshold is
therefore a rounding, and the only open questions are how large it is and
whether 0.5 is where it should sit.

## The net hides the rounding by a factor of seven

The report already carried the net, as `land_fraction_gauss_weighted` against
`mesh_land_fraction`. The rounding is two one-signed errors that partly cancel:
land in a sub-threshold cell is dropped into the slab ocean, and ocean in a
supra-threshold cell is promoted to dry land.

| rung | land dropped | water promoted | NET | rounding / net |
| --- | ---: | ---: | ---: | ---: |
| T21 | 11.58% | 10.56% | -1.14% | 10.1 |
| T42 | 7.74% | 6.82% | -0.99% | 7.8 |
| T85 | 5.38% | 4.43% | -1.00% | 5.4 |
| T127 | 4.37% | 3.44% | -0.96% | 4.5 |
| T170 | 3.67% | 2.93% | -0.76% | 4.8 |

Dropped as a share of the mesh's own land area; promoted as a share of the land
area the model receives.

**The net is nearly flat across the ladder and the rounding is not.** The net
moves by a factor 1.5 from T21 to T170 while the land actually misassigned
falls by 3.2. A convergence claim read off the net would report the coastline
as insensitive to support when what is happening is that two errors three times
larger are cancelling. That is the reading the two halves exist to prevent.

## The terrain the fork exists to preserve is what the threshold floods

`source/README.md`'s first rule rejects `land_mask` because it puts dry closed-
basin floor below sea level under water. The threshold does the same thing to
part of the same terrain, through a different door: a cell whose land share
falls below 0.5 goes to the slab ocean whatever its elevation.

| rung | share of below-datum land dropped to ocean | land volume closure |
| --- | ---: | ---: |
| T21 | 6.44% | +2.77% |
| T42 | 3.19% | +2.52% |
| T85 | 1.32% | +1.17% |
| T127 | 0.82% | +0.43% |
| T170 | 0.47% | +0.31% |

Land volume is the extensive closure: area times elevation about the datum,
summed over the mesh's land against the same sum over what the model receives.
It carries the sign the area does not, because the preserved terrain sits below
the datum.

This is the most rung-sensitive quantity in the ledger, falling by a factor of
14 across the ladder while the area rounding falls by 3.2. It is also the one
with a rule behind it rather than a preference, which is why it is tracked
separately from the coastline area.

## How much of the grid is partial at all

| rung | cells with a partial land share | share of planet area they hold | cells the mesh does not cover |
| --- | ---: | ---: | ---: |
| T21 | 1,053 | 55.5% | 0 |
| T42 | 2,715 | 36.0% | 0 |
| T85 | 6,833 | 22.6% | 0 |
| T127 | 11,578 | 17.1% | 32 |
| T170 | 16,610 | 13.8% | 365 |

This is the size of the problem a tile path would address, and it is the number
to weigh a tile implementation against: at T42 a third of the planet's area
sits in cells that are neither wholly land nor wholly ocean.

The last column is the opposite constraint and it grows the other way. Beyond
T85 the 10M mesh stops covering every Gaussian cell and `gridding.land_weighted`
falls back to the nearest region centre, which is a point sample of a
categorical field -- the exact reduction this builder exists to avoid. 365 cells
at T170 is small, but it means a finer rung is not a free route out of the
coastline rounding: past T85 the reduction starts substituting for the mesh.

## The third class is a latent hazard and reads zero today

`surface_class` partitions the mesh into ocean, subaerial land and inland
water, and the threshold is applied to the LAND share alone. Inland water
therefore competes with land for its own cell's classification, so a cell
holding no ocean at all goes to the slab ocean once its lakes take more than
half of it.

On this build the hazard is inert: `INLAND_WATER` holds 0 regions and 0.0% of
area, so the promoted inland-water area and the count of wholly terrestrial
cells sent to the ocean are both exactly zero. Both are in the ledger so that
the first build carrying inland water reports it rather than absorbing it into
the ocean share.

## What each candidate rule would cost instead

**Measured:** 2026-08-28, same build, same mesh, same ladder.
`build_boundary_conditions.py:coastline_ledger` now prices every candidate in
`candidate_rules` as part of each boundary build, so these are reproduced in
`boundary_conditions_report.json` rather than kept only here.

A threshold is a decision, and a decision taken against one rule's cost is not
a comparison. Three rules were proposed. They answer different questions and
only one of them has the fork's own rule behind it.

Below-datum land dropped as a share of the mesh's own; net land area against
the mesh's; water promoted as a share of the land the model receives; volume is
the extensive closure.

| rung | rule | threshold | below-datum dropped | net land | water promoted | volume closure |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| T21 | as configured | 0.5 | 6.44% | -1.14% | 10.56% | +2.77% |
| T21 | area-conserving | 0.491 | 6.05% | -0.01% | 11.01% | +3.67% |
| T21 | exempt any below-datum | 0.5 | 0.00% | +46.78% | 32.27% | +38.44% |
| T21 | exempt at 10% of the cell | 0.5 | 6.02% | -0.90% | 10.68% | +2.78% |
| T42 | as configured | 0.5 | 3.19% | -0.99% | 6.82% | +2.52% |
| T42 | area-conserving | 0.479 | 2.95% | -0.02% | 7.24% | +3.31% |
| T42 | exempt any below-datum | 0.5 | 0.00% | +21.53% | 19.12% | +13.82% |
| T85 | as configured | 0.5 | 1.32% | -1.00% | 4.43% | +1.17% |
| T85 | area-conserving | 0.474 | 1.13% | +0.01% | 4.90% | +1.69% |
| T85 | exempt any below-datum | 0.5 | 0.00% | +6.39% | 8.92% | +3.13% |
| T127 | as configured | 0.5 | 0.82% | -0.96% | 3.44% | +0.43% |
| T127 | area-conserving | 0.466 | 0.65% | +0.00% | 3.90% | +0.78% |
| T127 | exempt any below-datum | 0.5 | 0.00% | +2.47% | 5.53% | +1.16% |
| T170 | as configured | 0.5 | 0.47% | -0.76% | 2.93% | +0.31% |
| T170 | area-conserving | 0.465 | 0.39% | -0.00% | 3.30% | +0.56% |
| T170 | exempt any below-datum | 0.5 | 0.00% | +0.90% | 3.97% | +0.62% |

### The outright exemption trades a 6% error for a 32% one

Exempting any cell that holds below-datum land is the only rule that recovers
ALL of the terrain `source/README.md`'s first rule is about, and at T21 it does
so by inflating the land the model receives by 46.8%, promoting a THIRD of that
land out of open ocean, and taking the volume closure to +38.4%. That is the
same class of error -- a cell given the wrong surface -- at five times the
magnitude, and it reaches the volume the exemption exists to protect.

It is not absurd in itself: its cost falls by a factor 52 across the ladder and
at T170 it costs +0.90% of land area to recover the last 0.47% of below-datum
land, which is close to a fair trade. It is absurd at the rungs this world is
actually run at. **A rule cannot be adopted on the behaviour it would have at a
resolution the project does not use.**

### A share-based exemption is inert at any safe strength

Requiring the below-datum land to reach a share of the cell before the
exemption applies is the obvious way to blunt it. It does not work: at 25% of
the cell no cell qualifies at any rung, so the result is identical to 0.5
everywhere, and at 10% it recovers four tenths of a percentage point at T21.

The reason is the finding: **the below-datum land the threshold floods sits in
cells that are mostly OCEAN.** A sub-threshold cell is by definition less than
half land, and the below-datum part of that land is a fraction of a fraction.
So no share-based exemption strong enough to be safe is strong enough to
matter, and the only exemption that recovers the terrain is the one with no
share at all.

### Conserving area costs the closure that carries the sign

The area-conserving threshold is 0.491 at T21 falling to 0.465 at T170. It
conserves land area to within 0.02% at every rung, and it makes both of the
other numbers worse: it recovers almost none of the below-datum land, and it
takes the land volume closure from +2.77% to +3.67% at T21. It buys the
quantity that does not carry the sign of the preserved terrain by spending the
one that does.

### 0.5 is the only value that makes the rounding single

The other three rules are alternatives to a number. This is a reason for the
number itself. `oceanmod.f90:339-344` hard-binarises `yls` at 0.5 whatever this
builder writes. At a builder threshold of 0.5 the two agree and the cell is
rounded once; at any other value the mask written and the rounding the model
would apply to a fractional field disagree, and the cell is rounded twice by
two rules that do not know about each other. The area-conserving thresholds
above are all below 0.5, so every one of them is a second rounding as well as
a first.

### The verdict

**0.5 stays**, and what the measurement changes is the reason rather than the
value: the cost is sized, every alternative is priced beside it, and none of
them is better on the rule the fork exists for.

The route out of the residue is the RUNG. Below-datum land dropped falls by a
factor 14 from T21 to T170 while the coastline area rounding falls by 3.2 and
the net by 1.5. It is the only term in the ledger that support actually fixes,
and the constraint on going further is the mesh: past T85 the 10M reference
stops covering every Gaussian cell.

## What follows

- The two halves are now in every boundary-condition report. The net stays
  beside them; it is the right number for a mass budget and the wrong one for a
  convergence claim.
- Whether 0.5 is the right threshold is a separate question from whether the
  boundary should be fractional at all, and the section above answers the
  first: 0.5 stays, because every alternative is worse on the rule the fork
  exists for and only 0.5 rounds the cell once. The measured materiality test
  has answered the second in favour of a tile boundary. Code 1720 now carries
  the native-mesh subaerial share into every model start, separately named
  `dlf`, `ylf` and `xlf` in the land/atmosphere, ocean and sea-ice owners, while
  code 172 remains binary topology. A compiled endpoint-preserving combine
  exists so pure cells never multiply an absent tile by zero. Fourteen restart
  records now preserve both boundary publications, and land/glacier versus
  ocean/ice consume separate exchange bundles before a combined diagnostic is
  restored. Those bundles still mirror the one binary flux/radiation
  evaluation, so this is routed infrastructure and not yet tile physics.
- The flux half is now an executable, fail-closed measurement in
  `analysis/coastline_flux_bracket.py`. It requires the accepted baseline on
  this exact build, verifies its land mask, and has no bootstrap or stale-file
  fallback. Canonical-10m-carve2's accepted baseline now supplies that
  measurement: the global surface-energy bracket is -0.0719 to +0.1253 W m-2,
  against the same run's pre-existing 0.12 W m-2 state-storage tolerance, and
  the individual turbulent and radiative terms span 1.3 to 2.9 W m-2. The
  declared inequality therefore selects the conservative tile representation.
  `analysis/partial_surface_decision_gate.py` reproduces the choice, verifies
  the code-1720 carrier, and remains fail-closed until the model no longer
  hard-binarises the selected physical exchange.
