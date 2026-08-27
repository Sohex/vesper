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

    cut    = C * erodibility * (S / S_ref) ** n * Q ** 0.5      [metres]
    retain = clip(1 - cut / depth_at_spill, 0, 1)

Retain is then literally the fraction of the impoundment that survives, which is
what Orogen cuts with and what a reader of the map sees.

Three of the four terms are measured per basin. `erodibility` is the rock at the
sill; `S` is the gradient of the outflow channel below it, relative to the
population's, and the section on it below says why it cannot be folded into `C`;
`Q` is the overflow. `n = 1` is Orogen's own slope exponent, which GRAV-4 records
is baked into the Braun-Willett closed form the generator solves rather than
being a parameter of it.

`C` absorbs the incision coefficient and the relaxation window into one constant
with units of metres per (m3/s)^0.5 at land-mean rock and at the population's own
outlet gradient. **It is calibrated against Earth rather than declared**, and the
reason is that the relaxation window is undefined rather than unmeasured: Orogen
has no time axis, and `docs/src/reference/no-time-axis.md` carries the fact and the standing
way round it. What the calibration is solved against is a DENSITY, and a density
is a count over an area at a size class; the size class is a declared parameter
and is swept, under "The size floor is the largest lever" below.

**It is solved on every run rather than written down.** The calibration matches a
number of standing basins per unit land, and how many basins overflow at all is a
property of the climate, so the target moves whenever the verdict does. A literal
cannot follow it and did not: the 161 derived below was fitted against a verdict
taken before HYD-13, and against the verdict that replaced it that same 161
leaves 55 standing basins in the calibration band where Earth's density asks for
33.6. Recalibrating with no other change gives 235.9, outside the 128-to-218
bracket the 161 was published with.

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

## The slope term is in the export after all, and it is the basin's own depth

Settled 2026-08-17 as HYD-15, and it changed the mapping rather than a
coefficient. Recorded here because what stood in this place said the slope term
was unavailable, and that was true only of the place it was looked for.

**A saddle has no slope.** The cross-divide drop between `spill_region` and
`spill_exit_region` on this build has a median of 0.1 m with half the basins
negative, which is what a saddle is. But stream power is not written about the
divide; it is written about the channel the overflow runs DOWN, and that channel
is in the export. Following the drainage ten receiver steps from
`spill_exit_region` covers a median of 135 km, a dozen mesh regions. The estimate
is a property of the outlet rather than of the sampling: the same basins measured
at 5 steps and at 20 give r = 0.903 in the log with a median ratio of 1.00.

Measured on `precarve-craton` over the 2,937 basins whose outflow path runs at
least 50 km:

| | |
| --- | ---: |
| `r(log S, log depth_at_spill)` | **0.735** |
| `d log S / d log depth`, ordinary least squares | 1.04 |
| the same by reduced major axis | 1.41 |
| `r(log S, log Q)`, the confound check | 0.072 |

**The outlet's gradient and the basin's depth are the same relief counted
twice.** Under `S ~ depth^p` the effective depth exponent in `retain` is
`n*p - 1`, which at n = 1 and p = 1.04 is +0.02: the depth cancels. Absorbing
`S^n` into `C` while dividing by the depth was therefore not an approximation
with an unknown coefficient, it was a different mapping.

`S` is now carried explicitly and normalised by the population's geometric mean,
so `C` keeps the meaning its Earth calibration gives it and only the spread about
the median basin is new. The gradient earns its place on its own account too:
against discharge it measures r = 0.072, so it is new information rather than a
proxy for the water, and its spread is wider than depth's.

**The counts barely move and the composition changes completely**, which is the
same shape the previous revision of this mapping produced. Measured against the
verdict as it stood before the evaporation corrections of the same day, so that
the mapping is separated from the climate:

| | C | carve | marginal | preserve | standing < 35 deg |
| --- | ---: | ---: | ---: | ---: | ---: |
| as published, no slope term | 161.0 | 1598 | 123 | 1900 | 55 |
| recalibrated, still no slope term | 235.9 | 1648 | 73 | 1900 | 33 |
| slope term, n = 0.5 | 94.8 | 1628 | 93 | 1900 | 33 |
| **slope term, n = 1, adopted** | **102.9** | **1608** | **113** | **1900** | **34** |
| slope term, n = 2 | 330.2 | 1601 | 120 | 1900 | 33 |

The marginal class had a median depth at spill of 1,354 m against a population
median of 127; it now has 126, which is the population. The depth selection is
gone, exactly as the algebra says it should be. What selects a marginal basin now
is a flat outlet: median gradient 0.00015 against the population's 0.00068, four
and a half times flatter. 180 basins change side.

It does not collapse back to the 30 of the superseded discharge-only mapping,
because the outlet gradient has more spread than the depth does. A basin is now
marginal because its overflow leaves down a gentle channel, which is the right
reason for a notched valley to survive.

## The first real verdict, and how it came to be calibrated

Measured 2026-08-17 on `precarve-craton` under the baseline climatology, 3,621
basins. The verdict as it stands is **1,512 carved, 98 marginal, 2,011
preserved**; the paragraphs below are how that number was arrived at, because the
first pass got it wrong in a way worth keeping.

On that verdict there were **no lake-fed spillers at all**: every basin that
overflowed had catchment runoff doing it, so the case the discharge test was
rewritten to handle existed in the algebra and not on that terrain.

**On the corrected evaporation of 2026-08-17 there are 243 of them**, and they
are the reason the test was written that way. They sit at a median latitude of
70.3 degrees, where a lake evaporates 185 mm/yr against 269 mm/yr of
precipitation: cold basins fill from their own surface. They were invisible while
lake evaporation was floored at the model's land rate, which put open water at
the land's own evaporation and made `P > E` impossible almost everywhere. The
ratio test would have been undefined for every one of them.
`notes/audits/carve-criterion-terms.md` finding 1 carries the floor.

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
expected-value argument `docs/src/reference/no-time-axis.md` describes: Earth is one randomly
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

### What that determines, and why it is now solved rather than recorded

The density fixes the coefficient by solving for the value that leaves as many
standing basins in the band as Earth's density implies. On the verdict of
2026-08-17 that was 34 of 529 candidates and gave

    C = 161 m per (m3/s)^0.5,  bracket 128 to 218

The bracket is the Poisson error on Earth's 15, which is the dominant
uncertainty and is stated rather than hidden. Across it the marginal count ran
64 to 121, so quote the bracket with the count.

**Those numbers are the calibration of that day's verdict and are not the
coefficient in use.** The candidate count and the standing count both move with
the climate, so `export_carve_list.py` re-solves this on every run from the Earth
inputs -- 15 standing lakes over the land the HydroBASINS level 5 polygons cover
within 35 degrees, both measured by `--measure-earth-floors` -- against this
world's own land area in the same band. The denominator was 78.9 Mkm2 here and
in the script, which is not what those polygons cover; it is 77.05, and the
script now refuses when the two disagree. What is durable here is the Earth measurement and the
argument; the coefficient is a derived value and belongs in the sidecar with the
rest of them.

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

**What was genuinely open, and it was not gravity.** The mapping treated the cut
as slope-INDEPENDENT, absorbing `S^n` into `C` so that depth entered only as the
amount of rock to remove, while physically the sill's local slope is set by that
same relief. That was a structural question about the mapping rather than a
coefficient to scale, and it was worth more than the gravity factor. It is
settled above, under "The slope term is in the export after all": measured, the
two correlate at r = 0.735 with an exponent near 1, and at n = 1 the depth
cancels as predicted.

### What this model still cannot produce

It carves Tanganyika, Malawi and Albert. At any coefficient this calibration
gives, a sill carrying 1,491 m3/s is cut by far more than 577 m, so the model
does not reproduce those
basins -- and it reproduces the right NUMBER of standing basins by a different
mechanism, depth and low discharge rather than active subsidence. That is honest
about what is in the model: there is no subsidence term, so a basin whose floor
keeps dropping is not representable. The density is matched; the identity of the
survivors is not, and should not be quoted as though it were.


## The size floor is the largest lever, and it is now swept

Measured 2026-08-26 on `precarve-craton-10m` under the bootstrap climatology of
`run_2b20e3324bb0`, 9,419 basins, 5,565 of them overflowing, 177.52 Mkm2 of land
within 35 degrees. The bootstrap is the only climatology this build has, so what
follows is a statement about the SENSITIVITY and not a verdict; the ratios below
are what transfers, and the absolute coefficients move with the climate the way
the section above says they must.

**The floor on the Earth sample decides more of `C` than the calibration's own
error does.** The target is a density of standing through-flowing impounded
basins, and a density is a count over an area at a size class. Move the size
class and the count moves by two orders of magnitude while the area does not.
That was a single number defended in a comment as what the mesh resolves; it is
now `EARTH_SIZE_FLOOR_KM2` with `EARTH_STANDING_BY_FLOOR` beside it, and
`export_carve_list.py:sweep_size_floor` re-solves at every rung and writes the
result to the carve list's sidecar under
`method.calibration.size_floor_sensitivity`, beside the verdict.

### The span is derived from the build, not chosen

Both bounds are measured in the same run as the verdict, so they move when the
mesh or the terrain does, and neither can be picked to make the answer look
stable:

| | | on this build |
| --- | --- | ---: |
| lower | median land mesh cell area. Below it a counted Earth lake has no representable counterpart, because a depression smaller than one cell does not exist in the generator's output | 71.4 km2 |
| upper | median area at spill of the overflowing basins. Above it more than half the population being solved for sits outside the class the Earth sample stands for, so the Earth density is no longer a density of the same object | 2,687 km2 |

The floor in force, 1,000 km2, sits inside that span, which is the whole reason
the choice needed reporting rather than defending: it is a defensible value
among several defensible values, and nothing said what the others gave.

A rung is USABLE while the Earth sample keeps at least ten lakes. Below ten the
Poisson fractional error `1/sqrt(N)` exceeds 0.32, worse than the 0.258 the floor
in force already carries, and a density stops being a measurement. Thin rungs are
reported with their counts rather than dropped, because a reader has to see where
the sample runs out; a floor that leaves three lakes is a different failure from
one that leaves a hundred.

### The ladder

Earth counts re-measured from HydroLAKES v1.0 joined to HydroBASINS level 5, by
`export_carve_list.py --measure-earth-floors`, which is the query itself rather
than a record of it: natural lakes with mean depth above 5 m, pour point inside a
level-5 basin with `ENDO == 0`, within 35 degrees, counted strictly above each
floor on lake area. It reproduces the fifteen of 2026-08-17 exactly at 1,000 km2,
and gives 501 at 10 km2 against the 495 recorded then, a difference of six lakes
that moves nothing here. Eight pour points of 3,188 fall in no level-5 polygon
and are dropped.

Counts on the finished-depression basis, which is the basis `C` is solved on:

| floor km2 | Earth lakes | `C` | carve | marginal | preserve | |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 10 | 501 | 19.0 | 3354 | 2211 | 3854 | below the span |
| 20 | 311 | 34.9 | 4252 | 1313 | 3854 | below the span |
| 50 | 164 | 77.7 | 4897 | 668 | 3854 | below the span |
| 100 | 102 | 131.9 | 5143 | 422 | 3854 | |
| 200 | 57 | 230.9 | 5324 | 241 | 3854 | |
| 500 | 21 | 547.9 | 5477 | 88 | 3854 | |
| **1000** | **15** | **862.4** | **5504** | **61** | **3854** | **in force** |
| 2000 | 12 | 1076.9 | 5519 | 46 | 3854 | |
| 5000 | 5 | 1790.5 | 5543 | 22 | 3854 | too few lakes |
| 10000 | 2 | 4386.0 | 5556 | 9 | 3854 | too few lakes |

The preserved count never moves. That is the check that the sweep is measuring
what it claims to: preservation is set by the overflow test and `C` cannot reach
it, so a floor that moved it would mean the sweep was moving something else.

### The lever against the noise

The criterion was fixed before the sweep ran. `L` is the ratio of solved
coefficients over the rungs that are both inside the span and above the
ten-lake minimum; `P` is the ratio of the Poisson bracket on the Earth count at
the floor in force, solved on this same population. `L > P` is DECISIVE, the
floor moves the answer by more than the count's own error and has to be declared
and swept; `L <= P` is SUBORDINATE, it sits inside an uncertainty already
reported. Two classes, no gap and no overlap, and the same call is made a second
time on the marginal class.

| | across the span | Poisson at the floor in force | | |
| --- | ---: | ---: | ---: | --- |
| `C` | 131.9 to 1076.9, a factor of **8.17** | 610.8 to 1083.3, a factor of **1.77** | 4.6x | DECISIVE |
| marginal basins | 422 to 46, a factor of **9.17** | 81 to 45, a factor of **1.80** | 5.1x | DECISIVE |

**The floor is the larger lever by about five to one**, and the marginal class is
where that is felt: a landform class the project decides rather than inherits
runs from 422 basins to 46 across floors every one of which is defensible on the
resolution argument the single value was defended with. The Poisson factor of
1.80 on this build is the same instrument the earlier 1.9 measured on
`precarve-craton`, which is the point of quoting it: the two agree, so the
comparison is between a lever and a scatter that is known independently of it.

Neither call is close. `failure-modes.md` class 34 is the failure this comparison
exists to avoid, and it would have been a real risk had the ratio come out near
one: a sweep whose range never reaches the effect returns ordinary-looking
numbers.
