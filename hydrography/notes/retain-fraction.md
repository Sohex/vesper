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

## The first real verdict, and what it says about Q_full

Measured 2026-08-17 on `precarve-craton` under the baseline climatology, 3,621
basins: 1,580 carved, 30 marginal, 2,011 preserved.

All 30 are marginal for want of discharge. None is marginal for uncertainty, for
the structural reason given above, and there are **no lake-fed spillers at all**:
every basin that overflows has catchment runoff doing it. So the case the
discharge test was rewritten to handle exists in the algebra and not on this
terrain, which is the right way round -- the old ratio test would have been
undefined for it either way.

**`Q_full` is not worth calibrating**, and the discharge distribution is why. The
30 marginal basins overflow at 0.008 to 1.70 m3/s; the carved ones have a median
of 234 m3/s, a ninety-fifth percentile of 1,916 and a maximum of 38,650. The
boundary sits in an almost empty part of the distribution, so moving it costs
little:

| `Q_full` at land-mean rock | carve | marginal | preserve |
| ---: | ---: | ---: | ---: |
| 0.1 m3/s | 1605 | 5 | 2011 |
| 1.0 (declared) | 1580 | 30 | 2011 |
| 10 | 1490 | 120 | 2011 |
| 100 | 1127 | 483 | 2011 |

An order of magnitude either way takes the marginal count from 0.1% of basins to
3.3%. Two orders take it to 13%, which is where it would start to matter, and
100 m3/s is not a defensible reading of "perennial stream". The preserved count
does not move at all at any value, because it is set by the overflow test and not
by the incision mapping.

The number to keep an eye on instead is the 2,011: that is the overflow test, and
it is decided by the open-water evaporation, which is where the uncertainty
actually lives.
