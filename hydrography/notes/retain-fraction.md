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

## It now measures what the overflow can cut against what has to go

Two quantities, computed separately and combined by taking the larger, because
either was meant to be a reason to leave a rim standing. Only the first ever
decides anything; see below.

**The overflow itself**, as the water that has to leave at spill level:

    Q = runoff * (catchment - area_at_spill) - (E - P) * area_at_spill

**What that Q can cut, against how deep the basin is.** Stream power puts
`dz/dt` as `K Q^m S^n`, so what an overflow achieves over a relaxation window is
a LENGTH of incision. Whether that length empties the basin depends on the
distance from the spill point down to the floor:

    cut    = C * erodibility * Q ** 0.5                [metres]
    retain = clip(1 - cut / depth_at_spill, 0, 1)

Retain is then literally the fraction of the impoundment that survives, which is
what Orogen cuts with and what a reader of the map sees.

`C` absorbs the incision coefficient, the sub-grid slope term and the relaxation
window into one constant with units of metres per (m3/s)^0.5. **It is calibrated
against Earth rather than declared**, and the reason is that the relaxation
window is undefined rather than unmeasured: Orogen has no time axis, and
`notes/no-time-axis.md` carries the fact and the standing way round it. The
calibration is below and gives 161, with a bracket of 128 to 218.

A predecessor compared the cut against a declared discharge, `Q_full = 1 m3/s`,
and never asked about depth. It therefore said the same thing about a 14 m pan
and a 1,420 m trough -- the fifth and ninety-fifth percentiles of depth at spill
on this build -- and the marginal class came out as 30 basins of 3,621, all of
them trickles. Depth spans two orders of magnitude here and it is in the export
already.

No absolute time-to-cut is attempted even so, and none is available. What the
mapping carries is an ordering and a calibrated scale, where the version before
it had neither: a trickle and a torrent both went out at retain 0.

## The rock at the sill, and the number not to use

The cut is multiplied by the erodibility of the basin's own sill, so a soft rim
is cut through by less water than a hard one. The export's `erodibility` field is
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
resulting cut is scaled by 0.465 to 1.791. Sills sit in slightly harder rock than
the land as a whole, which is what a rim standing above a basin floor should do.
The rock term spans a factor of 3.85, against two orders of magnitude for depth
at spill and four for discharge, so it modulates the answer rather than setting
it.

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

## The first real verdict, and how it came to be calibrated

Measured 2026-08-17 on `precarve-craton` under the baseline climatology, 3,621
basins. The verdict as it stands is **1,512 carved, 98 marginal, 2,011
preserved**; the paragraphs below are how that number was arrived at, because the
first pass got it wrong in a way worth keeping.

There are **no lake-fed spillers at all**: every basin that overflows has
catchment runoff doing it. So the case the discharge test was rewritten to handle
exists in the algebra and not on this terrain, which is the right way round --
the old ratio test would have been undefined for it either way.

### The first pass gave 30, and 30 is not a small residual

It is a missing landform class. Marginal is not a hedge: it is the instruction
that makes Orogen cut a notch and taper it over the divide band, which is a
through-flowing valley holding a residual lake -- a distinct and conspicuous
landform, and one of the few this project decides rather than inherits. Thirty of
3,621 is a statement that Vesper essentially does not have that landform, and
that is a result rather than a rounding detail.

It was also unstable in exactly the way that should have prompted the check
rather than a shrug. Under the superseded discharge-only mapping, `Q_full` moved
it between 5 and 483 across two orders of magnitude, with the preserved count
never moving at all because that is set by the overflow test. A factor of four
either way on the size of a landform class is not a tolerance.

### The first check against Earth, and what it took to read it right

Measured 2026-08-17 from HydroLAKES v1.0 joined to HydroBASINS level 5, both
already in `hydrography/data/reference/` for HYD-4. The object being counted is
the terrestrial analogue of a marginal verdict: a lake that persists although a
river crosses its sill. Selection: natural lakes (`Lake_type == 1`) with mean
depth above 5 m whose pour point falls in a HydroBASINS basin with `ENDO == 0`,
at latitude under 35 degrees. `ENDO == 0` selects a sill a river actually
crosses; the latitude cut keeps out sills that were under an ice sheet 20,000
years ago and have had no time to be cut at all.

**Two readings of the same data give answers 300 times apart, and only one of
them is sound.** Both are recorded because the wrong one is the intuitive one.

**The survivor edge, which does not work.** Every standing lake requires its own
sill not to have been cut through its own depth, so `depth / sqrt(Q)` is bounded
below across the surviving population, and the bound looks like a measurement of
the incision constant. It gives 0.54, at Lake Albert. It is not usable: Earth's
large standing through-flowing lakes are almost all Tanganyika, Malawi, Albert,
Edward, Kivu, Toba, Towuti, Poso, Titicaca -- rift and volcano-tectonic basins,
maintained by active subsidence and young besides. What that edge measures is
how recently the basin floor last dropped, not how slowly the sill cuts.
Survivors are not a random sample of anything.

**The density, which does.** How many such basins stand per unit land is a
property of the population rather than of any member, and it is exactly the
expected-value argument `notes/no-time-axis.md` describes: Earth is one randomly
chosen moment, and so is this terrain.

**The size class decides the answer and nearly inverted it.** Earth has many
small through-flowing lakes and very few large ones, and this mesh cannot
represent the small ones: cells are 285 km2 and the overflowing basins on this
build have a median area at spill of 7,338 km2. Counting Earth's lakes above
10 km2 gives 495 within 35 degrees, a density that would put 1,136 standing
basins on Vesper against only 529 candidates -- which reads as "carve almost
nothing". Counting only those above 1,000 km2, which is what this mesh resolves,
gives **15**, a density of 0.19 per Mkm2, and 34 standing basins expected on
Vesper's 177 Mkm2 of land within the same latitude band. Against 529 candidates
that is 6.4 percent, and it reads as "carve almost everything". Match the size
class before reading a density.

### What that determines

Solving for the coefficient that leaves 34 of those 529 standing:

    cut    = C * erodibility * Q ** 0.5          [metres]
    retain = clip(1 - cut / depth_at_spill, 0, 1)

    C = 161 m per (m3/s)^0.5,  bracket 128 to 218

The bracket is the Poisson error on Earth's 15, which is the dominant
uncertainty and is stated rather than hidden. Across it the marginal count runs
64 to 121, so quote the bracket with the count.

| | carve | marginal | preserve |
| --- | ---: | ---: | ---: |
| C = 128 | 1489 | 121 | 2011 |
| C = 161 | 1512 | 98 | 2011 |
| C = 218 | 1546 | 64 | 2011 |
| superseded, Q_full = 1 m3/s | 1580 | 30 | 2011 |

**The counts barely moved and the composition changed completely**, which is the
result. Under the old rule the 98 were 30 basins overflowing at 0.008 to 1.70
m3/s: trickles, in whatever depth of basin they happened to sit. Under this one
they are deep basins with moderate discharge -- median depth 1,354 m against 115
m for overflowing basins as a whole, median overflow 38 m3/s against 228 -- and
their retain is spread across the range, median 0.35 with 53 of the 98 between
0.1 and 0.5. They are notched valleys holding substantial residual lakes, which
is the landform the marginal class exists to name, and membership is now decided
by the terrain rather than by an arbitrary discharge.

### Does the coefficient need a gravity term? No, and the reason is already in tree

Asked 2026-08-17 (PHYS-3). `C` was fixed by matching Earth's standing-basin
density, which is an Earth observable transferred to a 1.31 g world, and stream
power goes as `rho g Q S` -- so the objection was that this world cuts 31% faster
and should have FEWER standing basins than Earth, not the same number.

**It does not, because the 1/g relief scaling has already paid for it.**
`notes/audits/orogen-gravity.md` derives the compensation directly: at steady
state `K A^m S = U` gives `S ~ 1/K`, and stream power puts `K ~ rho g`, so
`S ~ 1/g`. Orogen scales every land height by `reliefScale = 0.7655464480874317`,
which is `1/(g_v/g_e)` to sixteen digits. The erosivity is up by g and the slopes
are down by g, and at grade the two cancel exactly -- as they must, because at
steady state the incision rate equals the uplift rate and gravity does not enter.

A basin sill is not at grade; it is a transient knickpoint. For a transient at
fixed local slope the residual is `g^(1-n)`, and that is where the whole question
lives:

| n | C factor | carve | marginal | preserve |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 | 1.1429 | 1618 | 107 | 1896 |
| **1.0, Orogen's own assumption** | **1.0000** | **1605** | **120** | **1896** |
| 2.0 | 0.7655 | 1553 | 172 | 1896 |

**So no correction is applied**, and the bracket is -52 to +13 basins rather than
the one-signed 30 to 40 the objection implied. `orogen-gravity.md` records that
"n = 1 is a choice, and nothing in the generator states it", which makes this the
same question as GRAV-4 and GRAV-5: **the gravity term and the slope exponent are
one question and cannot be settled separately.** Choosing a gravity factor here
while Orogen scales relief on n = 1 would double-count whatever n turns out to be.

**What is genuinely open, and it is not gravity.** This mapping treats the cut as
slope-INDEPENDENT: `S^n` is absorbed into `C`, so depth enters only as the amount
of rock to remove. Physically the sill's local slope is set by that same relief,
so `S` and `depth` are not independent, and a formulation carrying both would
have a different depth dependence -- at n = 1 with `S ~ depth`, the depth cancels
out of `retain` entirely. That is a structural question about the mapping rather
than a coefficient to scale, and it is worth more than the gravity factor was.

### What this model still cannot produce

It carves Tanganyika, Malawi and Albert. Under `C = 161` a sill carrying 1,491
m3/s is cut by far more than 577 m, so the model does not reproduce those
basins -- and it reproduces the right NUMBER of standing basins by a different
mechanism, depth and low discharge rather than active subsidence. That is honest
about what is in the model: there is no subsidence term, so a basin whose floor
keeps dropping is not representable. The density is matched; the identity of the
survivors is not, and should not be quoted as though it were.

