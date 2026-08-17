# What the retain fraction measures, and what it used to

Changed 2026-08-16, in `hydrography/scripts/export_carve_list.py`. Recorded here
because the verdict leaves this project and changes the terrain, so what the
number means has to be readable without reading the script.

## It used to measure our uncertainty, and hand that to a landscape

Retain is the fraction of a basin's rim Orogen leaves standing: 1 keeps it, 0
carves the outlet open, and a value between cuts a notch at the saddle and tapers
it over the divide band. It was computed as the distance from the basin's own
overflow threshold, in units of the evaporation uncertainty we happen to have,

    retain = min(1, ((E_penman - (P + critical * runoff)) / E_penman) / 0.25)

with every overflowing basin flattened to 0. So the notch Orogen cut encoded how
sure we were, not what the water could do, and a basin that trickled over its
sill was carved exactly as wide as one pouring a large catchment through it.

## It now measures what the overflow can cut, floored by what we do not know

Two quantities, computed separately and combined by taking the larger, because
either is a reason to leave a rim standing.

**The overflow itself**, as the water that has to leave at spill level:

    Q = runoff * (catchment - area_at_spill) - (E - P) * area_at_spill

**What that Q can cut**, through stream power, `dz/dt` going as `K Q^m S^n`:

    retain_incision = 1 - (Q / Q_full) ** 0.5,  clipped to [0, 1]

`Q_full` is 1 m3/s at land-mean rock, the perennial-stream scale: above it a
channel flows year round and works on its bed every year, below it the same
volume arrives as seasonal pulses that spend most of the relaxation window not
cutting. It is declared, not fitted, and it is recorded in the sidecar so a
verdict can be re-read against another choice. An order of magnitude either way
moves retain by about a factor of three at fixed discharge, so quote it with any
marginal count.

No absolute time-to-cut is attempted, and that is deliberate. What the mapping
carries is the ordering, which the old one did not have at all.

## The rock at the sill, and the number not to use

`Q_full` is divided by the erodibility of the basin's own sill, so a soft rim is
cut through by less water than a hard one. The export's `erodibility` field is
exactly the right quantity and says so in the manifest: a relative stream-power
multiplier from the exposed rock, mean-normalised to 1 over land. It enters where
K does.

**The contrast to apply is the expressed one, not the intact-rock one**, and the
two differ by orders of magnitude. Stock and Montgomery (1999) measure K between
1e-7 and 1e-2 across lithologies, which would make the rock the only term that
mattered and the discharge decorative. Zondervan et al. (2020) measure what a
real channel network expresses within one mountain belt at about a factor of 4,
against two orders of magnitude via intact strength, because channels adjust
their width and slope to the rock they are cutting. Both are read and indexed in
`references/INDEX.md`. Orogen's field spans 3.85x across the sills of this
build's basins, which is the right order, and it was grounded against Moosdorf
(2018), whose whole-Earth index spans 3.2x and was built for regional-to-global
erosion models.

The saddle is smaller than a mesh cell, so the regions either side of it bracket
the rock rather than naming it. The geometric mean of the two is used, that being
the right average for a multiplicative factor. They agree only loosely, at a
correlation of 0.34 on this build, so which one is picked would matter.

Measured on `precarve-craton`, 2026-08-16: all 3,621 basins resolve a sill on
both sides, erodibility there runs 0.465 to 1.791 with a median of 0.933, and the
resulting `Q_full` runs 0.56 to 2.15 m3/s. Sills sit in slightly harder rock than
the land as a whole, which is what a rim standing above a basin floor should do.
At an overflow of 0.5 m3/s the rock term alone spans retain 0.05 to 0.52, so it
is the same order as the discharge term rather than swamping it.

**The uncertainty is kept apart**, as the same margin expression as before, and
retain is the larger of the two. Where the two evaporation estimates disagree
about whether a basin overflows, the rim is kept.

**That second term has never once decided anything, and it cannot.** Taking the
larger of the two lets the margin raise a retain, never lower one, and the margin
is positive exactly when the basin does not overflow -- which is exactly when the
incision term is already 1. The two conditions are the same condition, read off
the sign of `Q` either way. Measured on the first real verdict: zero basins of
3,621 where `retain_margin` exceeds `retain_incision`, and `retain` equals
`retain_incision` everywhere. The expression is kept and reported in the sidecar
because it is the honest statement of how close a preserved basin sits to its
threshold, which is worth reading; it is not a second input to the cut.

## The same change fixes the dry pans, by construction

The old test divided by catchment runoff, so a basin with none was undefined and
was treated as never overflowing. It is not: a lake surface gaining more
precipitation than it evaporates grows, and since the imbalance does not close as
the surface expands, it grows until it spills. Those basins were being carved
wide open on one pass, when the margin expression lost its runoff term and
compared evaporation against precipitation alone, and then pinned shut on the
next by a guard that answered for the catchment and not for the lake.

`Q` needs no guard. With runoff zero it reduces to `(P - E) * area_at_spill`,
which is positive exactly when the lake surface gains, so those basins now go
through the marginal band on the same footing as every other. The run prints how
many came out that way.

## What is checked, and what is not

Checked on 200,000 synthetic basins: wherever catchment runoff is positive, the
discharge test and the ratio test it replaces agree on every basin, so the
verdict is unchanged where the old one was defined. Where runoff is zero the
ratio test carves nothing and the discharge test carves exactly the basins whose
lake surface gains water.

The sill lookup is checked on the real basin set rather than synthetically: the
join is by mesh region index and it resolves for every basin, with the terrain
hash on `basins.nc` verified against the export's before the lookup runs, since a
region index does not survive a terrain change.

What is still missing from the stream-power law is the slope term, `S^n`. The
saddle geometry that would give it is sub-grid, and unlike the rock it is not in
the export.

## The first real verdict, and what Earth says about Q_full

Measured 2026-08-17 on `precarve-craton` under the baseline climatology, 3,621
basins: 1,580 carved, 30 marginal, 2,011 preserved.

All 30 are marginal for want of discharge. None is marginal for uncertainty, for
the structural reason given above, and there are **no lake-fed spillers at all**:
every basin that overflows has catchment runoff doing it. So the case the
discharge test was rewritten to handle exists in the algebra and not on this
terrain, which is the right way round -- the old ratio test would have been
undefined for it either way.

### Thirty of 3,621 is not a small residual, it is a missing landform class

Marginal is not a hedge. It is the instruction that makes Orogen cut a notch and
taper it over the divide band, which is a through-flowing valley holding a
residual lake -- a distinct and conspicuous landform, and one of the few this
project decides rather than inherits. Thirty of them is a statement that Vesper
essentially does not have that landform. That is a result, and it should not
rest on a constant nobody has checked.

`Q_full` moves it by a factor of four either way within an order of magnitude:

| `Q_full` at land-mean rock | carve | marginal | preserve |
| ---: | ---: | ---: | ---: |
| 0.1 m3/s | 1605 | 5 | 2011 |
| 1.0 (declared) | 1580 | 30 | 2011 |
| 10 | 1490 | 120 | 2011 |
| 100 | 1127 | 483 | 2011 |

The preserved count does not move at any value, because it is set by the
overflow test. Everything in dispute is the split between a carved outlet and a
notched one.

### Earth falsifies Q_full = 1 m3/s outright

Measured 2026-08-17 from HydroLAKES v1.0 joined to HydroBASINS level 5, both
already in `hydrography/data/reference/` for HYD-4. The selection, declared
before the numbers were read: natural lakes (`Lake_type == 1`) with mean depth
above 20 m and area above 100 km2, whose pour point falls in a HydroBASINS basin
with `ENDO == 0`, at latitude under 35 degrees. Depth and area select an
impounded basin rather than a river widening; `ENDO == 0` selects one whose sill
a river is actually crossing, which is the situation `retain` describes; the
latitude cut keeps out sills that were under an ice sheet 20,000 years ago and
have had no time to be cut at all. That gives 34 lakes, and dropping five Amazon
floodplain water bodies, which have a river stage rather than a bedrock rim,
leaves 29. The extract is `hydrography/data/reference/exorheic_impounded_lakes.json`.

**Every one of the 29 has a through-flow above 1 m3/s.** The smallest is 1.3, the
first quartile is 23 and the median is 55. Under the mapping as declared, all 29
get `retain = 0` -- Tanganyika, Malawi, Albert, Toba and Titicaca's neighbours
included -- which is to say the rule as written erases every terrestrial example
of the landform it is supposed to produce.

### The deeper fault is that the rule never asks how deep the basin is

`Q_full` is a discharge, so it says the same thing about a 14 m pan and a 1,420 m
trough. Those are the fifth and ninety-fifth percentiles of depth at spill on
this build, and that field is in the export and currently does nothing. Earth is
unambiguous that it should: Tanganyika carries 1,491 m3/s across its sill and is
still a basin, because its floor is 577 m below the sill and the notch has taken
a small fraction of that. Depth is what decides whether a cut outlet leaves a
lake behind.

A mapping that asks the question would compare what the overflow can cut against
what has to go:

    retain = clip(1 - C * erodibility * Q**0.5 / depth_at_spill, 0, 1)

`C` absorbs the incision coefficient and the relaxation window into one constant
with units of metres per (m3/s)^0.5, which is a quantity Earth bounds directly:
every standing lake in the sample requires its own cut to be less than its own
depth, so `C` is bounded above by the minimum of `depth / sqrt(Q)` over the
sample. That minimum is 0.54, at Lake Albert.

| `C`, m per (m3/s)^0.5 | carve | marginal | preserve |
| ---: | ---: | ---: | ---: |
| 0.54, the Earth bound | 76 | 1534 | 2011 |
| 5 | 582 | 1028 | 2011 |
| 20 | 1093 | 517 | 2011 |
| 50 | 1335 | 275 | 2011 |
| 200 | 1534 | 76 | 2011 |

**The landform class is substantial at every value in that range**, which is the
finding. It was absent only because the rule could not see depth.

What is still open is `C` itself, and Earth bounds it from one side only. The 29
lakes bound the cut achieved over their own lifetimes, which for rift lakes is
short and is offset by subsidence; over the 1e5 to 1e7 years a landscape takes to
integrate, basins demonstrably do drain, so the true `C` is well above 0.54 and
nothing here says how far. That is HYD-8, and it is a declared parameter rather
than a measurement either way -- but declaring one that Earth's own basins
survive is a different act from declaring one that erases all of them.

