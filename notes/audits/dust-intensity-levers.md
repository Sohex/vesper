# Which lever moves the simulated dust intensity

Measured 2026-08-24 on build `precarve-craton-10m` at T21, against the
bootstrap regular climatology
`exoplasim/analysis/climatology/bootstrap_regular_climatology.nc` and the
baseline dust artifact `aeolian/analysis/dust_baseline.json` generated
2026-08-24. Vesper is a simulated super-Earth; every quantity here is a
property of this project's offline dust component, not of Earth.

Evidence: `aeolian/analysis/dust_intensity_levers.json`, written by
`aeolian/scripts/dust_intensity_levers.py`. Every figure below is a RATIO of
emission totals. Transport is not re-run, and the licence for that is measured
in the artifact's `emission_to_optical_depth_transfer` block: on the baseline
artifact's own three roughness arms, the ratio of emissions and the ratio of
land-mean optical depths agree to 2.6% at worst, because the transport operator
is linear in the source at fixed winds. A 2.6% transfer error against effects of
15% to 1000% is an instrument that can see them.

world-03x. The question was whether the dust intensity this component reports is
an artefact of tabulating the aeolian roughness by lithology, and whether the
soil texture rather than the roughness is the lever that sets it. Both questions
come from Menut et al. (2013), read after world-c8m closed. Neither claim
transfers, and the reasons are structural rather than matters of degree.

## The criteria, fixed before the ratios were computed

Two numbers already on the record set the scale a tabulation bias has to reach
to be the explanation. Closing world-c8m moved the central land-mean optical
depth by a factor of 4.4, and the roughness bracket spans a factor of 10.6 on
emission. So a bias under 1.25 in either direction cannot be the explanation; a
bias of 4.4 in the direction that makes the tabulated answer too high would
account for the whole world-c8m jump; anything between is a partial account.
These are in the artifact's `tabulation_bias.criteria` block.

## 1. The tabulation biases the emission LOW, not high

Menut's caution is that a tabulated roughness gives dust fluxes that are "on
average, higher than the other model configurations", because tabulated values
"are more discrete and thus less variable and realistic" than a satellite
retrieval. Tabulating is exactly the collapse of a within-class distribution to
a point, and the emission is not linear in the roughness, so the collapse has a
sign and a size that can be computed.

The within-class geometric spread is measured rather than declared, from the one
source in this repository that reports several in-situ roughness values for the
SAME area: Prigent et al. (2005) Table 1 carries three such groups, all from
Greeley et al. (1997). Ten values over the Death Valley box span 0.0696 to
1.06 cm; two over the Namib box and two over the central Nevada box give the
extremes. Pooling the within-box variance of the logarithm over all three gives
a geometric standard deviation of 2.736 with eleven degrees of freedom, and the
smallest and largest single-box estimates, 1.297 and 5.095, give the bracket.
Both bracket ends rest on two points. `aeolian/config/dust.yaml` carries the
derivation under `drag_partition.within_class_z0`.

That figure is an UPPER bound on the within-lithology spread, deliberately. A
Death Valley box holds playa, alluvial fan, dune and lava, and nothing in the
table separates them, so its spread is across surface types as well as within
one.

Resolving the distribution instead of collapsing it, at the central roughness
arm, multiplies the emission by:

| within-class treatment | ratio to the tabulated point |
| --- | --- |
| lognormal at the pooled spread 2.736 | 1.249 |
| lognormal at the smallest measured spread 1.297 | 1.017 |
| lognormal at the largest measured spread 5.095 | 1.400 |
| uniform in the logarithm, as wide as the class bracket | 1.116 |

Every one of them is above 1. The tabulation makes the modelled emission too
LOW, not too high, and the caution being tested has the wrong sign for this
scheme. It cannot be why the intensity is what it is.

**Prefer the bounded shape where the two disagree.** A lognormal has no smooth
end, so its integral keeps creeping up as the quadrature reaches into a tail no
measurement supports: the lognormal arms above still drift with the node count
where the bounded one has converged. The bounded shape spans the width the
source actually measured.

**The correction is one-signed at the central and rough arms and not at the
smooth one**, and the asymmetry is the drag partition's clamp rather than
physics. The rough arm gives 1.058 to 2.979; the smooth arm gives 0.916 to
1.047, and the one value below 1 is the widest lognormal at a class centre of
1.0e-5 m, where a quarter of that distribution falls below the smooth-bed
roughness the partition clamps at and cannot get any smoother, leaving the rough
half of the distribution unopposed. Read the bounded shape there, which is 1.047.

The rough arm is the one that matters for the reopening test in `notes/dust.md`,
because it is the only arm of the baseline artifact below the 0.10 land-mean
optical depth threshold. Resolving the tabulation lifts it above that threshold
at every measured spread except the narrowest, so the crossing the test recorded
is a floor as well.

## 2. On this world the roughness is an INTENSITY lever, not a spatial one

Menut's decomposition puts the roughness on the spatial side: for a given soil
data set, two roughness data sets differ mainly in the spatial distribution of
the flux. That was measured between two roughness MAPS that differ in space.
This project has one map with a level bracket on it, and a level bracket cannot
behave that way.

Measured across the baseline artifact's three arms, on the emission field
normalised by its own total:

| pair | pattern correlation | fraction of emission displaced |
| --- | --- | --- |
| smooth against central | 0.9950 | 0.061 |
| central against rough | 0.9838 | 0.115 |
| smooth against rough | 0.9614 | 0.175 |

Across a factor of 10.6 in the total, the pattern barely moves. The roughness
bracket is an intensity bracket here, and the imported framing does not apply.

## 3. Soil texture cannot be the intensity lever here, and the reason is the clay cap

Menut puts the intensity on the texture side. On this world it cannot be, for a
structural reason: the erodible substrate is closed-basin fill and evaporite
crust, both clay-rich, so Kok's clay fraction sits at its declared cap over most
of the emission. Measured on the baseline artifact, 86.9% of the emission and
91.6% of the erodible area carry a clay fraction at or above the declared
`f_clay_max` of 0.2. Where the cap binds, the texture cannot move the flux
through that term at all, and only the Fecan residual is left, which runs the
other way because more clay holds more water before the threshold moves.

The sensitivity, taken against this project's own Earth validation of the
pedology texture model (`pedology/analysis/earth_validation.json`, 16 sites,
clay RMSE 0.0835, regression of observed on predicted with slope 1.1096 and
intercept -0.0382):

| texture arm | emission ratio |
| --- | --- |
| clay recalibrated by the Earth regression | 0.944 |
| clay raised by the validation RMSE | 1.326 |
| clay lowered by the validation RMSE | 0.608 |
| clay halved | 0.402 |
| clay doubled | 1.691 |

The recalibrated arm is the defensible one and it moves the intensity by 5.7%.
Even the two arms that perturb the clay by a factor of two span 4.2, against the
roughness bracket's 10.6. The texture is not carrying uncertainty that the
roughness bracket has been blamed for.

## 4. The roughness bracket is one lithology's, and that is where a narrowing has to come from

Sweeping one class across its bracket while the other is held at its centre:

| held at its centre | class swept | bracket factor on emission |
| --- | --- | --- |
| `evaporite` | `playa_clastic` | 10.415 |
| `playa_clastic` | `evaporite` | 1.006 |

The full bracket is 10.501. So `evaporite`'s whole declared range is worth 0.6%
of the emission and `playa_clastic`'s is worth all of it, which follows from the
erodible weights: `playa_clastic` carries 1.0 and `evaporite` 0.1.

Two consequences. Reinterpreting `evaporite`'s bracket does not narrow the
intensity, which matters because that bracket is MacKinnon's spread over one
measurement set and is arguably a within-class spatial range rather than a level
uncertainty. And any future narrowing of the dust intensity has to narrow
`playa_clastic`'s roughness specifically, whose bracket runs from a modelling
convention at the smooth end ("a low surface roughness on the order of 0.001 cm
is used to describe active dust sources") to the sand-desert boundary stepped
down one order at the rough end. Neither end is a measurement of a clastic playa.

### The criteria for re-sourcing that bracket, fixed before any measured value was read

world-ugkw, declared 2026-08-25. What follows was written and committed before
the roughness literature held in `references/` was opened, so that the bar is
not the answer.

The narrowing available is not a better single value for a clay plain. It is
that `playa_clastic` is not a clay plain. Its derivation is geometric and its
composition is measured in `notes/audits/orogen-lithology.md`: half the class
sits above 83 per cent of its basin's relief, and on steepest descent it splits
into a playa and mud-flat band, a sand-flat and distal-fan band, and an
alluvial-fan and bajada band. So the question the sources have to answer is not
"what is a clastic playa's z0" but "which measured surfaces is this class a
mixture of, and in what proportion".

1. **Every end of the replacement bracket is a z0 MEASURED over a named
   surface.** What is being replaced is a modelling convention at the smooth end
   and a sand-desert boundary stepped down one order at the rough end, so
   neither a convention nor an order-of-magnitude step may reappear in the
   replacement. A value that cannot be attributed to a site and a method does
   not qualify.

2. **The mixing weights are the steepest-descent band shares already measured**
   in `notes/audits/orogen-lithology.md`, and the two builds' columns are the
   weight bracket rather than a value to choose between. Those shares are a
   statement about the class's COMPOSITION and survive the carve. The class
   EXTENT does not survive it, and no step of this construction may use the
   extent.

3. **The combination is the erodible-area-weighted GEOMETRIC mean** already used
   between classes in `aeolian/config/dust.yaml`, for the reason recorded there:
   the drag goes as 1/ln(z/z0), so it is ln(z0) that averages. Within a class
   and between classes are the same operation.

4. **THE PREDICTION THAT CAN FAIL.** The gradient bands say a third to a half of
   the class lies on sand-flat and fan gradients, and a gravel apron is rougher
   than a sand sheet rather than smoother, which is the assumption the stepped-
   down rough end contradicts. So the mix-derived HIGH end must come out ABOVE
   the declared 1.0e-4 m. If it comes out at or below 1.0e-4 m, that prediction
   is wrong: the argument from the gradient bands does not reach the roughness,
   and it is recorded as having failed rather than reworded into a weaker claim.

5. **Whatever the mix does to the emission bracket is adopted.** A replacement
   that WIDENS the factor of 10.5 is the result and not a reason to keep the old
   bracket, because the old bracket's ends are not measurements and the new
   ones are. Physics is not a knob.

6. **If the held sources carry no measured z0 for a mud flat or clay pan, the
   bracket is not narrowed by argument.** The issue is reported blocked with the
   paper named, and the declared bracket stands unchanged.

### What the measurements said, and it is not a narrowing

Measured 2026-08-25 on the same build and climatology, after
`aeolian/config/dust.yaml` took the values `aeolian/scripts/playa_roughness_mix.py`
derives. Evidence: `aeolian/analysis/playa_roughness_mix.json` and the re-run
`aeolian/analysis/dust_intensity_levers.json`.

**The endmembers.** Greeley et al. (1997) Table 2 carries aerodynamic
wind-profile determinations for all three of the class's landform bands, which is
the same quantity `evaporite` already takes from MacKinnon et al. (2004). Lunar
Lake, a 14 km2 silty-clay playa in central Nevada, is the only aerodynamically
measured clastic playa in the sources this project holds. Death Valley supplies
ten alluvial-fan determinations and Gobabeb two interdune-flat ones. Greeley's
own Table 2 was used and not Prigent et al. (2005) Table 1, which reproduces
those rows and prints Golden Canyon NE 1989 a factor of about 9 high.

**The prediction held.** The clause fixed above said the mix-derived high end
must exceed 1.0e-4 m. It comes out at 8.0e-4 m. The margin is not marginal: the
SMOOTHEST playa anyone has aerodynamically measured, 1.26e-4 m, already sits
above the entire bracket the clay-plain argument produced. Both ends of the old
bracket were below every relevant measurement, and the central value was low by
a factor of 13.

**The class ordering flips.** `playa_clastic` now comes out rougher than
`evaporite` rather than five times smoother. Where the two surfaces have both
been measured aerodynamically they agree to within 30 per cent, so what
separates the classes here was never clay against salt: `evaporite` is the
lowest quarter of a basin's relief and is floor, `playa_clastic` is everything
above it and is mostly margin, and the margin is the rougher surface. The
Marticorena et al. (1997) smooth-clay-against-embossed-gypsum distinction is
about a clay plain, and this class is not one.

**AND THE INTENSITY BRACKET WIDENS BY TWO ORDERS OF MAGNITUDE.** The z0 bracket
narrowed, from a factor of 10 to a factor of 8, and the emission bracket it
implies went the other way:

| | z0 bracket | emission bracket factor | drag partition, with over without, at the centre |
| --- | --- | ---: | ---: |
| the clay-plain argument | 1.0e-5 to 1.0e-4 m | 10.50 | 0.1065 |
| the measured mixture | 1.0e-4 to 8.0e-4 m | 1142.6 | 3.33e-4 |

Sweeping one class while the other is held at its centre still puts almost all
of it on `playa_clastic`: 141.98 against `evaporite`'s 1.86, so the conclusion
that this is the one class whose roughness matters survives, in a stronger form.

The mechanism is the drag partition, and it is arithmetic rather than a surprise.
The efficiency falls from 0.527 at the new bracket's smooth end to 0.238 at its
rough end, against 0.846 to 0.527 across the old bracket. Because the threshold
friction velocity goes as 1/feff and the emission integrates a cubic above that
threshold over a Weibull, the same factor in z0 buys a far larger factor in
emission the further down feff has already fallen. The class did not become more
uncertain; it moved into a steeper part of a scheme it was previously sitting
near the flat end of.

Criterion 5 above governs the verdict and was fixed before any of this was
computed: a replacement that widens the emission bracket is the result and not a
reason to keep the old one. The old bracket's ends are not measurements.

**What the instrument can carry, checked before the number is quoted.**
`aeolian/config/dust.yaml` records that the Shao and Lu approximation the
partition uses "is within 20% of a rigorous integration for feff above 0.2". At
the new bracket's rough end feff is 0.238, so the whole bracket is inside that
range but its rough end is close to the edge, and a rougher value than 8.0e-4 m
would leave it. The second limit is the subgrid wind: at the new central
roughness the partition suppresses emission to 3.33e-4 of an unpartitioned bed
against 0.1065 before, so the surviving emission comes from further out in the
Weibull tail, and the declared shape parameter carries correspondingly more of
the answer than it did. Neither limit is a reason to prefer the old values, both
of which were further from any measurement; they say what the next narrowing has
to be about, and it is the wind tail rather than the roughness.

**One consequence for every consumer.** Absolute emission at the new central
roughness is smaller than at the old one by roughly the ratio of those two
partition factors. `aeolian/analysis/dust_baseline.json`, `dust_baseline.nc` and
everything derived from them are WORTHLESS rather than stale, in the sense rule
7 gives the word.

**What did NOT change, and it is the half that matters for sizing.** The class
EXTENT is a pre-carve limit and no step of this derivation used it. The band
shares are a composition statement and survive the carve. So the leverage this
issue rests on is still an upper bound on the carved world's, and this
re-sourcing does not change that: it changes what the class IS, not how much of
it there is.

## 5. What the export's geometry can and cannot supply

world-03x proposed building the within-class distribution from this project's
own geometry, on the grounds that a playa has a smooth interior and rougher
margins and the basin solution knows each basin's area, sink and outlet status.
That was tested and it does not work, for two measured reasons.

Connected components of the erodible substrate on the mesh are not playas
(artifact block `export_geometry_scale`). Taken on `substrate_class` restricted
to `surface_class == 1` and with the solved lake extent removed, the 836,954
erodible regions of the 10M-region export form 27,401 connected patches whose
area-weighted 99th percentile is 4.6e6 km2. Those are basin-fill provinces, not
depressions. The area-weighted median height above a patch's own floor is 81 m
and only 15% of erodible area lies within 10 m of it, so "height above the playa
floor" is not what the quantity measures.

The land mask there is `surface_class` and nothing else. `land_mask` would flood
the dry closed-basin floor this component emits from, and the naive
`land_mask | is_endorheic` reconstruction would miss the depressions too small
to have entered the basin catalogue, which are precisely the population a
playa-margin roughness question is about.

And the scale gap is six orders of magnitude. The mesh has a 7.60 km mean edge
and Orogen designs terrain down to about 20 km, while the drag partition
partitions stress between a bed and centimetre-scale roughness ELEMENTS. That is
the same category error `aeolian/config/dust.yaml` already records for the
grid-cell roughness field, restated at region scale: no geometric quantity in the
export can measure an aeolian roughness. The geometry can order sediment fineness
at basin scale, and the lithology split already carries the one distinction the
literature supports for that, which is the detrital-fill against evaporite
contrast in Marticorena et al. (1997).

This is also why the tabulation bias above depends on the SPREAD and not on how
the spread is arranged: a T21 cell is several hundred kilometres across and
already contains the whole within-class population, so the arrangement of
roughness inside it does not reach the emission integral.

## 6. Two collapses that were priced and kept

**The between-class collapse.** `source_fractions` replaces the two lithologies
in a cell with one erodible-area weighted geometric mean before the emission
integral, which is the same species of collapse as the tabulation and needs no
declared input to undo. Resolving it changes the emission by 0.3% at the smooth
arm, 0.6% at the central and 1.4% at the rough. The collapse stays; a correction
that size does not earn a second code path.

**The scalar in-model roughness.** Boundary field 1802 carries the drag partition
per cell, so the mosaic reaches the model there, but `dustsrc` divides by the
scalar namelist `dustz0` in the `ln(zref/z0)` term, and that cannot be folded
into 1802 because the model computes `zref` from each cell's own temperature.
Two questions follow with very different answers.

Whether the scalar was the RIGHT scalar: it was not. It carried
`aeolian_z0_m`, the class-blind value whose bracket ends predate world-c8m and
sit outside the per-lithology ones entirely, so the in-model arm emitted 0.646,
1.746 and 4.049 times the offline arm at the smooth, central and rough ends. That
is a comparison defeated before it starts, and it is fixed under world-h24h:
`DUSTZ0` is now the mosaic's own erodible-area weighted geometric mean, computed
from the fields being written rather than looked up in the config.

Whether a scalar can do at all: it can. With the mosaic's mean in it the in-model
emission sits at 0.997, 0.997 and 0.996 of the offline arm. That residual is the
whole of what a fourth boundary field carrying the roughness per cell would
recover, so the fourth field, its surface code, its gather, its staging and its
self-test arm are not worth adding.

## 7. Using a drag partition at all is worth far more here than on Earth

Menut's third caution is that using a drag partition lowers dust fluxes by a
factor of 2 to 3 against not using one, which matters when these numbers are set
beside a model that omits it. That is an Earth number over Earth's roughness
distribution, so it is measured here instead of borrowed, against a drag
efficiency of one, which is a bed whose only roughness is its own grains.

| roughness arm | emission with the partition, over emission without |
| --- | --- |
| smooth | 0.367 |
| central | 0.106 |
| rough | 0.0157 |

A factor of 2.7, 9.4 and 64. The partition is worth several times more here than
Menut's range, and steeply more at the rough arm, because this world's declared
playa roughnesses sit further above the smooth-bed value than the Earth
distribution he averaged over. This is not a bracket: the partition is the
physics, and the figure is what to quote when comparing against a model that
does not partition the stress.

## 8. What the in-model arm's missing vegetation bracket costs

`build_dust.py --variant arid_bare_ground` lets a cell drier than the declared
250 mm per Earth year emit regardless of lithology at the declared bare-bedrock
weight of 0.30. The test is per time bin, so a static boundary field cannot carry
it and the in-model arm has no counterpart; that is world-4qem. Measured, the
variant raises the emission by 1.165, 1.176 and 1.198 at the smooth, central and
rough ends. The gap is therefore one-signed and under a fifth, against a
roughness bracket of 10.6 and a subgrid wind shape worth a factor of 40. An
in-model total is a floor by that much rather than a quantity of unknown
standing, which is what makes the gap declarable rather than a defect.
