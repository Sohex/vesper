# The land support's mosaic crossing: where the tiles die, and what that costs

**Measured:** 2026-08-31 on `canonical-10m-carve2` at T21, terrain hash
`f496ae9f`. Nothing was simulated: no ExoPlaSim run, no hydrography solve, no
generation. `hydrography/scripts/lake_mosaic_cost.py` is the measurement and
`hydrography/analysis/lake_mosaic_cost.json` its output.

Worldbuilding. Vesper is an invented super-Earth. Everything below is about the
simulation of its land surface: the modelled lake set `surface_water.py`
solves, the bucket capacity `landmod.f90` runs on, and the two laws that read
that capacity.

LSHY-6 asks the native-mesh to climate-grid crossing to keep the upland-soil,
lake/playa and groundwater-fed fractions distinct, and to reject a cell-mean
available water capacity as an area proxy.
`notes/audits/ocean-support-nonlinear-reductions.md` measured this class of
loss on the ocean side; its method is reused here rather than re-derived, and
its central finding transfers: **a gate that asks only for conservation cannot
see this loss, because the affine quantity passes it exactly.** There the
affine quantity was ocean volume in depth. Here it is the bucket capacity in
the lake fraction.

## Result, in one table

| arm | verdict | in the bar's units |
| --- | --- | --- |
| the mosaic shares partition the cell | CONTROL, PASSES | float32 storage precision on the shares, exact on the areas |
| the crossing itself keeps the tiles distinct | PASSES, and this is not where the loss is | upland-soil, lake and barren shares are emitted separately |
| the whole lake tile, absent today | **NOT MATERIAL on the planet mean, and it is the larger of the two** | 0.63 to 0.80 times the bar, one-signed, and -6.9 W m-2 on the worst cell |
| the mixing term, if the blend were staged | NOT MATERIAL | 0.44 to 0.50 times the bar |
| the sign of the mixing term | **NOT ONE-SIGNED, and a pre-registered invariant said it was** | 13.1% of lake-area bins carry the opposite sign |
| the capacity the flux arms stand on | the bucket the climatology's own run integrated, and it is NOT the one staged now | zero cells above capacity, maximum fill fraction 0.99997 |
| the groundwater-fed tile | ABSENT, and refused by name rather than defaulted | no `water_table.nc` on the accepted build |

## 1. The crossing is not where the mosaic dies

`build_spatial_support.py` emits the cell's land as distinct conservative area
shares through `lib/gridding.py`'s categorical operator: `f_land`,
`f_solved_lake`, `f_barren` and `f_nonbarren`, with a `hydrologic_share` table
beside them. They partition the cell. Measured over the covered cells, the
surface-class shares sum to 1 within 3.0e-8, the barren pair within the same,
and the land, ocean and inland-water areas sum to the mesh area exactly. The
first two residuals are the float32 the shares are stored at; the areas are
float64 and close to zero.

**So the answer to LSHY-6's question as posed is that the crossing preserves
the mosaic.** The loss is one step further down, at the consumer, and this is
worth separating because a row that reads "preserve the mosaic through the
crossing" invites a fix to the crossing, which is already correct.

## 2. Where it dies: one scalar for a two-tile distribution

`build_surface_soil_water.py --lakes` writes

    C = (1 - f_lake) * C_soil + f_lake * lake_dwmax_m

into surface code 0229, and that single capacity is what the model's two land
laws read. The blend is literally the area-weighted mean of a two-point
capacity distribution, so it adds no curvature of its own and conserves the
mean capacity exactly. Both consumers are nonlinear in it:

- saturation-excess runoff, `max(0, w + F dt - C) / dt`, `landcolumn.f90:84-86`
- the evaporation wetness factor, `min(1, w / (0.4 C))`, `landcolumn.f90:370`

**The staged field carries no lake blend at all.** Its provenance record
carries `"lakes": null`, so what the world runs today is the upland-soil
capacity alone on every cell, and the lake tile is not mixed in badly -- it is
absent. That makes two different losses, and they are priced separately below,
because a bound on a term that is inert has to say which of the two it bounds.
The same is unrecorded for the field the run integrated: an inert blend is a
property of a generator invocation, and the record of that invocation stays
beside the staged field rather than travelling with the copy in the run
directory. `world-hl06` carries that gap.

## 3. What each loss is worth

The bar is `config/partial_surface.yaml`'s materiality test: the accepted
baseline run's own 0.12 W m-2 state-storage tolerance, which that file already
selected the tile operator against and which the land audit's albedo arm
reports in. A capacity error reaches it through the wetness factor and thence
the latent heat flux.

**The instrument is bounded rather than divided, and the check that forced this
is worth stating.** The latent heat flux is the wetness factor times a
potential flux, so recovering the potential flux as `hfls / beta` looks exact.
It is not usable: `beta` is under 0.05 on 48% of this planet's land bins and
the quotient reaches 1e17 W m-2. The potential flux is capped at the surface
energy actually available, `rss + rls`, and the result is a bracket whose upper
end gives every cell its whole net radiation to evaporation. A verdict that
holds at the upper end holds.

**The pair is one state.** The wetness factor is a function of the soil water
and the capacity together, so both come from the run whose climatology the soil
water is: `exoplasim/runs/run_67323a923013/N032_surf_0229.sra`, resolved through
`lib/provenance.py:run_surface_field`. `exoplasim/inputs/t21/` holds the field
the NEXT run will read and it is a later iteration -- the two differ on 1598 of
2048 cells. Section 5's arm is the one that asks about the next staging, and it
keeps the staged field.

| loss | planet mean, estimate | upper bound | multiples of the bar | worst cell |
| --- | ---: | ---: | ---: | ---: |
| the whole lake tile, absent today | -0.076 W m-2 | -0.096 W m-2 | 0.63 to 0.80 | -6.9 W m-2 |
| the mixing term, if the blend were staged | +0.052 W m-2 | +0.060 W m-2 | 0.44 to 0.50 | +5.3 W m-2 |

**Neither reaches the bar on the planet mean, at either end of the bracket**,
and the absent tile clears three quarters of it at the upper end. The absent
tile is about one and a half times the mixing term, which is the ordering that
matters for what to do: staging the blend buys back rather more than half of
what is missing, and refining the blend into tiles buys back the rest. The two
are close enough that neither step is the one that matters on its own.

The signs are not symmetric and the reason is the direction defect below. The
absent-tile loss is negative, meaning that staging the lake tile would LOWER
the wetness factor and lower the latent heat flux. That is the opposite of what
a lake is for.

### The conservation control, which is where the pairing was caught

The store is clipped at the capacity every timestep -- `bucket_step` at
`landcolumn.f90:75-89` sets it to `min(C, w + F dt)`, and the layered path
clips each layer at its share of the same capacity -- so a bin mean of clipped
values cannot exceed the clip and the right answer for "land cells whose soil
water exceeds the capacity" is exactly zero. Against the field the run
integrated it is zero, and the maximum fill fraction over land is 0.99997: the
bucket touches its capacity and does not pass it. That is the clip visible in
the data rather than argued from the source.

Against the STAGED field the same control returns 55 of 1,639 cells, worst
0.068 m on a 0.049 m bucket. That is not a defect in either artifact. It is the
signature of pairing a climatology with a bucket from a different iteration,
and it is the reason the instrument now refuses that pairing by name rather
than reporting around it.

## 4. The mixing term is not one-signed, and that changes the repair

`world-cyu3` pre-registered, before any measurement, that the mixing gap
"vanishes at `f_lake = 0` and at `f_lake = 1` and is one-signed in between,
because a two-point distribution has one sign of curvature". The endpoints hold
exactly. The one sign does not, and the selftest found it rather than the
measurement.

**The wetness factor is not convex in the capacity.** `min(1, w / (0.4 C))` is
flat for `C` below the saturation knee at `w / 0.4` and convex above it, so the
derivative jumps from zero down to `-1/C` at the knee and the function is
CONCAVE there. A two-point distribution is one-signed only when the law it is
pushed through has one sign of curvature, and this one does not. A cell whose
soil tile is saturated while its lake tile is not straddles the knee and
carries the opposite sign, and the sign reverses within a single cell as the
lake fraction is varied.

Measured over the solved lake area, bin by bin, on the capacity the run
integrated: 74.9% of it has both tiles on the unsaturated branch, where the
pre-registered sign holds, and **13.1% has tiles that straddle the knee**,
where it does not.

The runoff consumer is unaffected: `max(0, w + F dt - C) / dt` is a maximum of
zero and a decreasing affine function of the capacity, so it is convex
throughout and a two-point distribution through it is one-signed.

**The consequence is about the repair, not about the size.** A correction term
with one sign cannot track a gap with two, so if this crossing is ever taken
past the cell mean, what it must carry is the DISTRIBUTION and not a corrected
mean. That is the same disposition the ocean audit reached from a different
mechanism, and `lib/gridding.py:cell_quantiles` is the operator in both cases.

## 5. The blend's direction is inverted against both capacity fields

`config/planet.yaml` argues `lake_dwmax_m` shallower than 0.5 m so that a lake
cell saturates its wetness factor on less water and evaporates at the potential
rate. 0.5 m is ExoPlaSim's uniform default, and `soil_water_source` is
`pedology`, so 0.5 m is not what either field carries. `lake_dwmax_m` EXCEEDS
the soil capacity on cells holding **76.1% of lake-bearing land area on the
field the run integrated** and 55.8% on the staged field, so blending it in
RAISES the capacity there, which raises the water needed to reach 40% of it and
lowers both the wetness factor and the saturation-excess runoff. That is the
direction the comment argues against, and section 3's negative sign is that
inversion in W m-2.

The direction argument itself is sound. What is wrong is that it was written
against the uniform default and never re-read against the pedology field that
replaced it. `world-kvr` owns the repair and the value's derivation. The
defect is not an artifact of which iteration is read: it is worse on the field
that ran than on the one staged.

The perturbation the blend applies, as a fraction of the cell's own soil
capacity over lake-bearing land area, differs between the two fields where the
share above does not. On the field that ran it is 1.3% at the median, 21% at
p90 and reaches 50%, and no part of that area is perturbed by more than 100%.
On the staged field it is 0.5% at the median and 15% at p90, reaching 435%,
with 2.6% of the area perturbed by more than 100%. The staged field therefore
carries a thin tail of cells whose capacity the blend would multiply, and the
field that ran does not; both carry the same systematic direction. The term is
concentrated rather than absent, which is why the worst cell reaches -6.9 W m-2
while the planet mean stays under the bar.

## 6. The groundwater-fed tile, and what it would take

It does not exist and is refused by name in `hydrography/config/wetness.yaml`
rather than defaulted, which is the right disposition for a class with no
source. Two candidate routes remain and they fail differently:

- GW-26's saturated-area closure is withdrawn permanently. `world-n9sd`
  measured the instrument and found the gain the terrain half carries over the
  cell-mean depth at or below its own scatter on all eight arms of both Earth
  sets. An AREAL observation would license a revival; a narrower criterion
  would not.
- PLHY-5's `groundwater_access` needs no closure at all, and needs
  `water_table.nc`, which no build in `source/` has. **That is a checkpoint
  with a price rather than a run to start**: `config/pipeline.yaml` puts the
  `groundwater` step at `cost: hours`, and it was measured this week running
  alone under the host lock at 2600% CPU and 16 GB, passing 40 minutes without
  writing before it was killed. `world-wfge` carries the cost and the lever.

## What follows

- The crossing preserves the mosaic and the consumer collapses it. A fix aimed
  at `build_spatial_support.py` would be aimed at the half that is correct.
- Staging the lake tile is worth about one and a half times what refining the
  blend into tiles is worth, and both are under the bar on the planet mean.
  Neither is a reason to change the staging on materiality grounds; the reason
  to change it is section 5, where the term that IS staged has the wrong sign.
  The two terms are close enough that a repair which staged the blend and
  stopped would leave nearly as much on the table as it recovered.
- A tile operator for this crossing has to carry the distribution rather than a
  corrected mean, because the wetness factor's gap changes sign at the
  saturation knee and 13.1% of lake area sits across it.
- A capacity and a soil water have to come from one iteration. Both doors are
  in `lib/provenance.py`: `run_surface_field` for the field a run consumed,
  `staged_surface_field` for the field the next run will read, and the second
  refuses the cross-iteration pairing when a run is named.
- The solved lake set every share here is weighted by is forced by the
  bootstrap climatology while a baseline exists. `surface_water.py` already
  resolves the best available climatology, so this is an artifact older than
  the baseline rather than a step pinned to the wrong stage, and re-solving it
  costs minutes.
