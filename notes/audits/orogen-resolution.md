# Orogen at higher region counts: what converges, what does not, and what breaks

Vesper is a generated world. This note asks whether World Orogen's terrain
generation is trustworthy above the region count the current build uses, and
what raising it would and would not buy. Measured on 2026-08-20 against the
`precarve-craton` build (2,500,001 regions, mesh spacing 15.19 km).

The occasion was a downstream question: several open rows (GW-6, GRAV-6,
DUST-16, SURF-7) are blocked on structure below the mesh scale, and it was
worth knowing whether that structure is simply waiting at a finer generation.

## What was run

Same planet code, radius, gravity and lithology strength as the build recipe in
`source/README.md`, varying only `--regions`. Explicit flags beat the decoded
planet code, so the region count is the only thing that moves.

| regions | mesh km | wall clock | peak RSS | raw output |
| --- | --- | --- | --- | --- |
| 100,000 | 75.95 | 4.8 s | 0.26 GB | no |
| 400,000 | 37.98 | 18.1 s | 0.54 GB | no |
| 5,000,003 | 10.74 | 399.5 s | 5.0 GB | no |
| 10,000,005 | 7.59 | 1092.8 s | 9.0 GB | yes, 4.7 GB |
| 25,000,001 | 4.80 | 3874.7 s | 21.1 GB | yes, 12 GB |

Mesh km is `avgEdgeKm`, `pi*R/sqrt(N)`, the generator's own definition and the
one the build's 15.19 km comes from.

Two control pairs at 100,000 and 400,000 with `--glacial 0.3` and `--glacial 0`
isolate the glacial contribution to the resolution shift.

The sanity, hypsometry and relief measurements are re-runnable as
`analysis/orogen_resolution.py`, registered in `config/pipeline.yaml` under
`one_offs`, which writes `analysis/orogen_resolution.json` with its provenance.
The tables below are that product. The generation runs themselves are not kept:
a build is disposable until a climate run has consumed one, and none has.

## Orogen's design ceiling is 2,560,000 regions, and the build sits at it

`js/detail-scale.js` maps the detail slider to a maximum of 2,560,000 regions.
The build's 2,500,001 is 98% of that and 25x the headless exporter's default of
100,000. Every result below is therefore outside the range the tool was tuned
in, which is why the mechanism audit matters more than the wall clock.

Cost scales close to linearly in both time and memory, about 1 KB of resident
memory per region. Nothing structural fails.

## Raising the region count refines the world, it does not regenerate it

Plate centres are chosen by farthest-point sampling from a single index-random
first seed, and continent seeds likewise over plate centroids. Only the first
seed's longitude moves with the region count; the geometry that follows is
determined. The consequence is that the same planet code gives the same world
at any region count, which the grid confirms.

| pair | cellwise r | RMS difference | land-mask agreement |
| --- | --- | --- | --- |
| 100,000 vs 400,000 | +0.9286 | 0.4273 km | not measured |
| 2,500,001 vs 5,000,003 | +0.9684 | 0.2883 km | 96.80% |
| 2,500,001 vs 10,000,005 | +0.9685 | 0.2861 km | 96.36% |

Gridded land fraction moves 0.5056, 0.5081, 0.5095 across the three.

The agreement does not improve between two and four times the region count,
which invites the reading that the build sits a fixed distance from a converged
limit. It does not mean that. A run at 2,600,001 regions, 4% away, agrees with
the build no better at +0.9622, so this correlation measures how far a PERTURBED
MESH moves the world, not how far coarse spacing does. See the noise floor
below before reading any number in this table as a resolution effect. What the
table does establish is the thing it was built for: the same planet code returns
the same world, and neither the region count nor the mesh perturbation that
comes with it regenerates the geography.

## The endorheic basin catalogue is resolution-limited, not physics-limited

`selectionCriteria` declares three floors: 0.05 km depth below spill, 1000 km2
floor area, and 12 mesh cells. The 12-cell floor is a signal-to-noise
criterion, asking whether the depression is a landform the mesh resolves, and
it is the one that binds at the build's region count.

| regions | cell km2 | 12-cell floor | binding floor | preserved | smallest preserved |
| --- | --- | --- | --- | --- | --- |
| 2,500,001 | 293.80 | 3526 km2 | cells | 3,621 | 2,646 km2 |
| 5,000,003 | 146.90 | 1763 km2 | cells | 6,345 | 1,227 km2 |
| 10,000,005 | 73.45 | 881 km2 | AREA | 9,419 | 1,000 km2 |

At the build's region count no depression below 2,646 km2 survives, against a
declared area floor of 1000 km2. The numerical floor is 3.5x the physical one,
and only near 8.8M regions does the declared criterion take over. Doubling the
region count nearly doubles the preserved set and halves the median preserved
area, while the largest basin agrees to 2.1% (3,059,665 km2 against
2,995,613 km2). The large features are converged; the small ones are being
admitted for the first time, correctly, as the mesh becomes able to hold them.

At four times the region count the smallest preserved depression is 1000 km2
exactly: the declared area floor has taken over from the mesh, which is the
crossover the arithmetic predicted near 8.8M regions.

| regions | mean edge | detected | preserved | smallest preserved |
| --- | --- | --- | --- | --- |
| 2,500,001 | 15.19 km | 81,871 | 3,621 | 2,646 km2 |
| 5,000,003 | 10.74 km | 154,865 | 6,345 | 1,227 km2 |
| 10,000,005 | 7.60 km | 286,534 | 9,419 | 1000 km2 |
| 25,000,001 | 4.80 km | 646,725 | 9,649 | 1000 km2 |

**And there it stops.** Ten times the build's region count adds 48% to the
preserved set over the step before it; twenty-five times adds 2.4%. Detected
depressions go on climbing, 286,534 to 646,725, because the mesh keeps resolving
smaller hollows, and every one of them is under the declared 1000 km2 floor and
none is kept. That is the crossover doing exactly what it was predicted to do:
while `minCells` binds, the catalogue grows with the mesh; once the physical
criterion takes over, the catalogue is a property of the world and refining
further buys nothing. The endorheic basin catalogue CONVERGES at about ten
million regions, and the 25M run cost 3,875 s and 21.1 GB, superlinear in time
against 10M's 1,093 s and 9.0 GB, to add 230 basins.

The large end does not settle, and that is the part that matters. The biggest
preserved basin is 3,059,665 km2 at the build, 2,995,613 at twice the regions
and 5,499,601 at four times. Its sink is at 80.5 degrees south on the build and
85.6 degrees south at four times the regions, so this is one south polar basin
whose mapped extent grows by 80%, not two basins merging. Total preserved area
goes 84.9, 93.0, 98.5 million km2.

Glacial erosion is not the cause. Repeating the 100,000 and 400,000 pair with
`--glacial 0` returns basin catalogues byte-identical to the `--glacial 0.3`
pair at both region counts, 107 and 574 preserved with the same three largest
areas, while the terrains themselves differ. The catalogue is computed on the
pre-conditioning surface, upstream of erosion, exactly as `catalogueStability`
in the manifest says. What is resolution-dependent is the depression
delineation itself.

This is the one finding that argues for a finer generation on its own merits,
because the carve list is drawn from this catalogue, and both ends of it move:
the small basins are admitted for the first time and the large ones are
remeasured.

## What in the generator is resolution-invariant, and what is not

The audit matters more than any single comparison, because the build already
sits at the design ceiling and nothing above it was calibrated.

| mechanism | carries a physical length? | verdict |
| --- | --- | --- |
| tectonic collision timestep | `dt` divided by `sqrt(N / 10000)` | scales, CFL-like |
| detail noise | fixed frequency in space, sampled finer | scales |
| hydraulic stream power | `flowSeed` is `cellArea / refCellKm2`, a real drainage area, and `distScale` carries the radius | scales, deliberately |
| lithology scarp band | `SCARP_EDGE_KM / edgeKm` hops, 150 km either way | scales |
| basin cell floor | 12 cells | relaxes correctly as the mesh resolves smaller depressions |
| `dist_*` export fields | BFS hops, declared `units='cells'` in the manifest | convertible, honestly labelled |
| erosion pass counts | fixed per slider, no region term | fixed erosion duration, which is what convergence wants |
| glacial ice accumulation | `iceFlow[target] += iceFlow[r]`, seeded from a per-cell index | DOES NOT SCALE, see below |
| scarp gradient thresholds | absolute gradient, 0.002 to 0.010, measured over ONE mesh edge | measurement depends on cell spacing; the thresholds are inherited, not anchored |

Two entries need their own argument.

**Glacial ice accumulation keeps the cell-count behaviour the hydraulic path
had already abandoned.** `flowSeed` exists precisely so accumulated flow is a
drainage area and not a cell tally, and the comment above it says so. The
glacial accumulator was not given the same treatment: it seeds from the
glaciation index and sums downstream, so the same valley collects roughly four
times the ice flow at four times the region count. Carving goes as
`iceFlow^0.6` and fjord incision as `iceFlow^0.5`, and the 0.1 and 0.5 flow
thresholds are crossed far more readily. The two accumulators disagree about
what a resolution-independent quantity is.

**Scarp gradient thresholds are absolute but the gradient a mesh can measure is
not.** Measured 2026-08-24 on the raw mesh of both builds, area-weighted over
`surface_class == 1`.

**This is NOT the compound topographic index's failure mode, and the difference
decides the remedy.** CTI shifts by exactly `ln 2` between these two builds
because its argument `a / tan(beta)` carries a LENGTH, so halving the cell edge
halves it deterministically and the shift can be subtracted. A gradient is a
length over a length and carries no length at all, so there is no dimensional
shift to remove. It fails to transport for a different reason: the surface is
self-affine, so the gradient measured over the sampling interval is itself a
function of that interval, and the shift is empirical rather than derivable. It
is also not one number. Land gradient on `elevation_km`, steepest descent to a
lower land neighbour:

| area quantile | p50 | p75 | p90 | p95 | p99 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2.5M, 15.19 km | 0.00132 | 0.00779 | 0.02878 | 0.04824 | 0.09768 |
| 10M, 7.60 km | 0.00186 | 0.01101 | 0.04424 | 0.07764 | 0.16557 |
| ratio | 1.408 | 1.412 | 1.537 | 1.610 | 1.695 |

The steep tail steepens faster than the median, so no single factor and no
single exponent removes the shift. The exported field moves with it: cells with
`scarp_potential` above zero go 0.2822 to 0.3410 of land, above 0.25 they go
0.0506 to 0.0949, and the mean where non-zero goes 0.1247 to 0.1924.

**A second defect was measured in the same pass and is unrelated to
resolution.** The relief term differenced `r_elevation`, the generator's
shaping parameter, against a denominator in kilometres. That is not a rise over
a run: `dh/de = 120 e^3 (1 - e)` is zero at sea level, peaks near 12.7 km per
unit at `e = 0.75`, and returns to zero at the ceiling, so one coded gradient
meant anything from flat to precipitous depending where on the hypsometric curve
the cell sat, and it was blind to this planet's 1/g relief scaling. The coded and
physical gradient distributions differ by a factor of seven at p99. The numerator
now converts through `scaledHeightKm`, the same converter `elevation_km` is built
with, and `reliefScale` is passed in.

**The neighbour filter was an elevation-sign test inside a function whose own
land test is `surface_class`.** `r_elevation[nb] <= 0` discarded every downhill
neighbour on a dry closed-basin floor below sea level, so a plateau margin
standing inside one could not be measured at all: 0.0528 of land at 2.5M and
0.0857 at 10M, of which 82 and 88 per cent respectively is endorheic. The
comment defended it as excluding bathymetry, which `isLand` does correctly and
without discarding the fork's own terrain. Fixed, and the three assertions in
`tools/test-lithology.mjs` that made the same substitution are fixed with it.

### What transports, measured against a bar fixed first

Bar declared before the runs: every land quantile from p50 to p99 must agree
between the two builds within 1.15x. Baseline length declared before the runs at
30 km, which is above Orogen's ~20 km terrain-information floor, one fifth of
`SCARP_EDGE_KM`, and a near-integer hop count on both meshes.

| relief estimator | p50 | p75 | p90 | p95 | p99 | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| drop to a lower neighbour over one mesh edge | 1.408 | 1.412 | 1.537 | 1.610 | 1.695 | fails |
| drop to the MINIMUM land elevation within 30 km, over 30 km | 1.285 | 1.053 | 1.110 | 1.206 | 1.247 | fails |
| drop below the MEAN land elevation within 30 km, over 30 km | n/a | 0.909 | 0.932 | 1.046 | 1.134 | passes where defined |

The minimum-over-a-ball form fails for a reason worth keeping: a minimum is an
extreme-value statistic and is biased by how many cells the ball holds, and the
fine mesh holds three times as many over the same ground. A mean is normalised
and carries no such bias. Its p50 is not a failure but a structural property: the
measure is one-sided, so half of land sits at or below its own 30 km mean and the
median is exactly zero, which leaves the bar undefined rather than missed.

**So a normalised form transports and a rank form is not needed.** A rank form
would also pin the areal extent, but it makes the field non-local, and this
measurement shows a local estimator is enough.

**What is NOT settled is the pair of numbers.** Adopting the 30 km form changes
the distribution the smoothstep reads by an order of magnitude, so 0.002 and
0.010 cannot be carried across, and choosing new ones by matching the areal
extent the present gate happens to produce would be calibrating against the
defect. What settles them is a measured escarpment: relief over 30 km at the foot
of a Great Escarpment-type plateau margin, which is an external number this
project does not hold. Until it does, the thresholds stay at their inherited
values and are labelled as inherited in the module, a scarp fraction is not
comparable between region counts, and `world-xgaj` carries the change.

## The realisation noise floor, and the two claims it killed

Changing `--regions` does not only change resolution. The first plate seed is
drawn as an INDEX, `randInt(numRegions)`, and the Fibonacci generator puts index
`f*N` at a fixed latitude but an N-dependent longitude, so every region count is
a slightly different realisation of the same world. Any global statistic
therefore differs between two region counts for two reasons, and only one of
them is resolution. Nothing in this note can be believed until the other one is
measured.

The control is three runs within 8% of each other in region count, where the
resolution is to all intents identical and anything that moves is scatter:

| statistic | 2,500,001 | 2,600,001 | 2,700,001 | spread | 2.5M to 10M |
| --- | --- | --- | --- | --- | --- |
| mean land elevation, km | 0.3957 | 0.4401 | 0.3879 | 0.0523 | -0.0355 |
| land area above 1 km | 0.12428 | 0.13666 | 0.12261 | 11.0% | -10.3% |
| gridded r against the build | -- | +0.9622 | +0.9737 | -- | +0.9685 |

**Quadrupling the region count moves mean land elevation and the area above
1 km by LESS than a 4% change in region count does.** The apparent hypsometric
bias is scatter. So is the gridded disagreement: 2,600,001 regions agrees with
the build no better than 10,000,005 does, +0.9622 against +0.9685, so the
correlation plateau near +0.97 measures how much a perturbed mesh moves this
world and not how far the build sits from a converged limit.

Two claims made earlier in this note's own life did not survive that control and
are withdrawn rather than softened:

- that land hypsometry is biased high and does not converge. The three-point
  sequence 0.3957, 0.4083, 0.3603 is not monotone, and the whole excursion fits
  inside the noise floor. There is no demonstrated resolution bias in
  hypsometry.
- that the build sits about 0.29 km RMS per T42 cell from converged terrain,
  with 3.6% of cells disagreeing about land against sea. That RMS is realisation
  scatter between mesh instances. A nearby region count reproduces it.

What DOES clear the floor, at more than twice the spread, is every statistic
whose change is mechanically predicted by cell size rather than statistically
inferred:

| statistic | 2,500,001 | 10,000,005 | noise | clears |
| --- | --- | --- | --- | --- |
| land area as channel, A above 1000 km2 | 33.32% | 15.54% | 1.53 | yes |
| preserved basins | 3,621 | 9,419 | 217 | yes |
| smallest preserved basin, km2 | 2,646 | 1,000 | 193 | yes |
| cells with scarp potential above zero | 12.46% | 15.00% | 0.52 | yes |
| mean scarp potential where non-zero | 0.1316 | 0.1986 | 0.005 | yes |
| land slope 99th percentile, deg | 5.015 | 8.427 | 0.35 | yes |
| upland concavity theta | 0.392 | 0.444 | 0.119 | NO |

Concavity is the instructive failure. At nearby region counts it reads 0.392,
0.273 and 0.354, so a spread of 0.119 swallows the 0.052 it gains at four times
the count. It was quoted in an earlier draft as the headline evidence that finer
terrain is more physical. It is not evidence of anything at this sample size.

The rule this leaves behind: on a generator whose mesh moves when its resolution
does, a resolution claim needs a same-resolution control before it is a claim at
all, and the cheap version is two runs a few percent apart.

## The glaciation index is Earth-calibrated and blind to this planet's climate

Separately from resolution, and larger: the ice that carves this terrain is not
this planet's ice.

    thresholdLat = PI/2 - glacialStrength * PI / GLACIAL_LAT_DIVISOR
    latFactor    = smoothstep(polarDist, thresholdLat, PI/2)
    elevFactor   = smoothstep(r_elevation[r], GLACIAL_ELEV_LOW, GLACIAL_ELEV_HIGH)
    glacIdx[r]   = max(latFactor, elevFactor * ...) * glacialStrength

The whole index is latitude and elevation. At the build's `glacialErosion` of
0.8 the ice line falls at 58 degrees latitude and ramps to full at the pole,
placed by the erosion slider and by nothing else. The altitude gate is worse
than it looks: `GLACIAL_ELEV_LOW` and `GLACIAL_ELEV_HIGH` are in the model's
dimensionless elevation parameter, so the snowline tracks maximum relief, which
scales as 1/g, rather than a freezing altitude, which is set by the lapse rate
and the surface temperature.

Nothing in it consults temperature. Not the host star's spectrum or
luminosity, not obliquity, not eccentricity, not the rotation period, and not
even Orogen's own temperature field, which cannot be an input because
`runPostProcessing` runs before `computeTemperature` in the same pass.

Every one of the declared parameters in `config/planet.yaml` that bears on an
ice line differs from the Earth these constants were tuned against, and the
three with the most leverage point the same way, toward less polar ice than a
58 degree threshold implies:

- obliquity 32 degrees against Earth's 23.4 sends more summer insolation to
  high latitudes and weakens the annual-mean pole-to-equator gradient
- the host is a K2.5V star, so the spectrum is shifted into the near infrared
  where snow and ice are markedly less reflective, weakening the ice-albedo
  feedback that sustains a low-latitude ice margin
- a 30 hour rotation period widens the overturning circulation and moves heat
  poleward more efficiently

So the build's mountains and fjords were cut by an ice mask that is probably
too extensive and reaches too far toward the equator, at a slider setting of
0.8 where the effect is strong rather than incidental.

The correction is available and does not need new physics. The ExoPlaSim
climatology on this build carries `glac`, `snd`, `sic`, `snm` and `alb`: the
model computes, under this planet's actual obliquity, orbit and stellar
spectrum, the ice distribution Orogen is guessing at. Driving `glacIdx` from
that field instead of from `polarDist` is a loop A feedback of exactly the kind
the loop exists for, and there is already a quasi-equilibrated run on
`precarve-craton` to take it from. The terrain-before-climate ordering that
blocks it inside one pass is not a constraint across iterations.

This subsumes the mechanism half of GRAV-6. That row parks a gravity term on
`iceFlow` as precision theatre on a heuristic accumulation, which is right, but
it understates the problem: the accumulation is also resolution-dependent, and
the index feeding it is calibrated to the wrong planet. A `g^3` Glen's law
velocity applied to an Earth ice mask would be precision on the wrong quantity.

## Relief is converged, so there is no sub-grid terrain waiting at a finer generation

This is the question the note was opened for, and it wants an estimator that
does not punish the coarse mesh twice. Averaging the standard deviation of
elevation inside a cap of given radius does exactly that: the coarse mesh has
fewer cells per cap, three against ten at a 15.19 km radius, and a standard
deviation from three samples is biased low, so the coarse mesh scores low partly
for having less relief and partly for having less data. Pooling every pair
instead removes the dependence on how many cells fall in a cap:

    gamma(r) = 0.5 * mean[(z_i - z_j)^2] over land pairs separated by at most r

with `sqrt(gamma)` in the units of the field. Measured over 6,000 land centres.

| lag km | 2,500,001 | 10,000,005 | ratio |
| --- | --- | --- | --- |
| 5.00 | 0.1571 | 0.1427 | 0.908 |
| 10.00 | 0.1820 | 0.1982 | 1.089 |
| 15.19 | 0.2270 | 0.2382 | 1.049 |
| 25.00 | 0.2771 | 0.3097 | 1.118 |
| 50.00 | 0.3807 | 0.4241 | 1.114 |
| 100.00 | 0.5075 | 0.5352 | 1.055 |
| 200.00 | 0.6292 | 0.6384 | 1.015 |

Quadrupling the region count adds at most 12% to relief at any separation from
5 to 200 km, and measured against the noise floor even that overstates it. The
same lags across 2,500,001, 2,600,001 and 2,700,001 regions scatter by 0.0237 at
15.19 km, 0.0193 at 25 km, 0.0277 at 100 km and 0.0377 at 200 km, against a
four-times-the-count difference of 0.0112, 0.0326, 0.0277 and 0.0092. The
resolution effect sits INSIDE the realisation scatter at every lag but 25 and
50 km, where it reaches 1.7 times it. The conclusion here is therefore stronger
than the 12% suggests, not weaker: most of what looked like extra relief is a
perturbed mesh, not a finer one. On `elevation_pre_erosion` the ratios run 0.95 to 1.00 from 15 km
outward, so the tectonic and noise substrate is converged outright and the
small residue is erosional texture the finer mesh cuts at its own scale.

The gridded within-cell measure agrees and is not even monotonic: `orog_std`
goes 0.4290, 0.4414, 0.4325 across the three region counts, and its 99th
percentile on land 1.8438, 1.8464, 1.7993. That is noise around a converged
value, not a trend.

Read this alongside the hypsometric result above, which does not converge. The
two are compatible and not redundant: relief BETWEEN neighbours at a fixed
separation has settled, while how much land AREA sits high has not. A finer
generation buys the second and not the first.

So the answer to the question that prompted this is no. A finer generation does
not unlock a reservoir of sub-grid relief, because the generator does not have
one to unlock: its noise sits at fixed physical wavelengths and its erosion
cuts at whatever the mesh scale is. Terrain-scale structure below the mesh has
to come from a parameterisation or from external data, not from more regions.

The arithmetic of the alternative is worth stating so it is not proposed again.
Mean edge is `pi*R/sqrt(N)`, the same `avgEdgeKm` the generator uses, which is
where the build's 15.19 km comes from. Resolving a glacial valley of 1 to 5 km
needs edges near 1.5 km, about 256 million regions, 26 times the run measured here and far outside anything the generator or the host will carry.
GRAV-6's judgement that the process is sub-grid holds at four times the region
count and would hold at forty.

## Scarps do strengthen with region count, and the thresholds are why

| measure | 2,500,001 | 10,000,005 |
| --- | --- | --- |
| cells with `scarp_potential` above zero | 12.456% | 14.999% |
| mean `scarp_potential` where non-zero | 0.1316 | 0.1986 |
| land fraction with gradient above 0.002 | 0.5036 | 0.5351 |
| land fraction with gradient above 0.010 | 0.2643 | 0.3008 |
| land slope 99th percentile | 5.015 deg | 8.427 deg |

Mean scarp potential rises 51% and the steep tail of the slope distribution
rises 68%, while relief at fixed separation rose at most 12%. The gap between
those numbers is the calibration: the landform is better captured, and on top
of that the absolute 0.002 to 0.010 smoothstep is crossed more readily by a
gradient measured over shorter runs. The two cannot be separated at the field
level, which is what LITH-25 records.

## Sanity of the four times run

Nothing structural fails at four times the design ceiling.

| check | 2,500,001 | 10,000,005 |
| --- | --- | --- |
| `cell_area` sum over `4*pi*R^2` | 1.0006809744 | 1.0006854371 |
| non-finite elevations | 0 | 0 |
| maximum elevation | 4.593 km | 4.593 km |
| minimum elevation | -8.894 km | -8.697 km |
| land fraction by region count | 0.4318 | 0.4329 |
| land fraction by area | 0.4317 | 0.4328 |
| deviation of positions from the unit sphere | 3.3e-16 | 3.3e-16 |

The area ratio is the centroidal dual's 1.00068 recorded in
`notes/audits/mesh-dual-area.md`, reproduced here to four more digits at four
times the region count, which is a check on the reader as much as on the mesh.
Maximum elevation is identical because it is the height curve's cap and not an
emergent quantity.

## Two hypotheses this note proposed and refuted

Both were mine, both were plausible from the code, and both cost one cheap
control to kill. Recording them because the pattern is the point.

Glacial erosion was the obvious cause of high ground lowering as the mesh
refines, since the glaciation index gates on altitude and `iceFlow` is
resolution-dependent. `--glacial 0` reproduces the profile to within 1%.

Glacial erosion was then the obvious cause of the largest basin growing 80%,
since that basin is polar and `latFactor` saturates at the pole. The catalogue
is byte-identical with ice off.

`iceFlow`'s cell-tally accumulation is still a real defect, recorded against
GRAV-6. It is simply not what either of these effects was. A mechanism that
predicts the sign, the gradient and the location of an effect can still not be
its cause, and the control that turns it off is minutes of compute.

## Do the landforms carry signal, or just more cells?

Converged relief could mean the extra cells add nothing. They do add something,
but not where the elevation statistics look: what improves is the DERIVED
structure, which depends on local gradient and connectivity rather than on how
much relief exists.

**Drainage concavity does not clear the noise floor, and is not evidence.**
Fluvial landscapes satisfy `S ~ A^-theta` with theta about 0.4 to 0.6, and on
upland terrain above 0.5 km the fit gives +0.392 at the build against +0.444 at
four times the region count, both inside the physical range and the second on a
tighter fit. That reads as the generator becoming more physical, and it was
quoted that way in an earlier draft. Nearby region counts give 0.392, 0.273 and
0.354. The measurement is too noisy at this sample size to carry the claim, and
it is withdrawn.

Worth keeping is what the attempt taught about conditioning. Over all land the
same fit gives +0.534 against +0.328 and appears to say the opposite, because
this world is mostly plains and a global statistic measures the flat regime
where the relation does not hold; conditioning on slope instead, `S > 0.01`,
destroys it from the other side by selecting on the dependent variable. Three
conditionings, three answers, and only one of them is a fluvial measurement.

**Channel articulation, which does clear it.** Cells with more than 1000 km2 draining through them
cover 33.3% of land area at the build's spacing and 15.5% at four times the
region count. A third of all land being channel is not a drainage network; the
finer mesh makes the threshold selective, because 1000 km2 is three cells
upstream at 15.19 km and fourteen at 7.59 km. The hillslope-to-channel break
moves from 759 km2 to 132 km2, toward the sub-km2 value real landscapes show.

So the honest summary of what a finer generation buys is: not more terrain, but
better-resolved drainage, basins and scarps. Every statistic that survives the
noise floor is one whose change cell size PREDICTS, and the one that failed was
the one inferred statistically. The elevation field has converged and the fields
derived from its connectivity have not, which is the same reason GW-6's water
table texture cannot be recovered by refining and HYD-20's basin catalogue can.

## T42 cannot see the ground that would glaciate, and that is why `glac` is zero

The ice mask PHYS-13 wants was to come from the commissioned climatology's
`glac`. Measured, that field is IDENTICALLY ZERO in all twelve output bins over
all 8,192 cells, while the run manifest records `"glaciers": true`. The module
ran and formed no permanent land ice. `snd` peaks at 0.29 m, which is seasonal
snow, and `snm` and `prsn` are not written in this run's output at all.

Zero is not a finding about the planet. It is a finding about the grid, and the
reason is the one thing a T42 cell cannot carry: the height of the ground inside
it. Against the mesh, the sub-grid peak excess `orog_max - orog_mean` on land is

| | km |
| --- | --- |
| mean | 1.065 |
| median | 0.499 |
| 90th percentile | 2.969 |
| maximum | 5.270 |

and the warm-season environmental lapse rate is 7.372 K/km, measured from the
climatology this project actually consumes rather than from a raw run file. So
the surface energy balance is evaluated about 7.8 K too warm for the high ground
in an average land cell and about 21.9 K too warm in the top tenth. Ice does not
survive that, and the model is right to report none on the orography it was
given.

These figures and the ones below are `analysis/ice_mask_freezing_height.py`,
which is registered and re-runnable. Earlier drafts of this section quoted
slightly different numbers, 7.271 K/km and 9.9 K, taken from the raw run rather
than the climatology product; where they disagree the script is what to believe,
because it is the one that can be checked.

Correcting the warmest-month surface temperature to each region's own elevation,
which is what `maps/build_basemap.py` already does for the basemap's ice
shading and what `lib/lapse.py` names the freezing-height criterion, moves the
answer by nearly three orders of magnitude:

| | fraction of land area below freezing in the warmest month |
| --- | --- |
| on the model's own orography | 0.002% |
| corrected to the 15.19 km mesh | 1.657% |

and it is a thermal criterion, so a margin is the honest sensitivity: 5 K colder
than freezing leaves 0.613% and 10 K leaves 0.181%. Repeating it on the 7.59 km
build with `--build` gives the finer-mesh figure, which an earlier draft put at
1.574% against 1.395% and which sits inside the realisation noise this audit
measures elsewhere, so it is not evidence that a finer mesh finds more glaciable
ground.

The consequence for PHYS-13 is a change of field, not of plan. `glac` is the
wrong input and would hand Orogen an ice-free world. The right input is the
warmest-month surface temperature, which is smooth, synoptic, and genuinely
resolved at T42, evaluated against the MESH's elevation through the measured
lapse rate. T42 is sufficient for the climate and insufficient for the mask,
and those are not the same statement.

Raising the model's resolution is not the alternative. T85 is about 150 km
against a mesh at 7.59, so nearly all of the peak excess above stays sub-grid,
and the correction is needed at any truncation the model can afford. This is the
same argument the downscaling machinery exists for.

Two honest limits on the 1.7%. The criterion is thermal and says where
ice can PERSIST, not where a glacier forms, which additionally needs
accumulation; `prsn` is absent from this run's output so accumulation cannot be
checked here. And the difference between the two meshes is inside the realisation noise floor
above, so it is not evidence that the finer mesh finds more glaciable ground.

## Orogen's information floor is about 20 km, and no region count buys past it

Sub-kilometre base topography is the natural next ask, and it fails twice: on
arithmetic, and then, more decisively, on there being nothing down there.

The arithmetic first. Mean edge is `pi*R/sqrt(N)`, so at this radius

| regions | mean edge km |
| --- | --- |
| 10,000,005 | 7.60 |
| 25,000,000 | 4.80 |
| 50,000,000 | 3.40 |
| 100,000,000 | 2.40 |
| 576,900,000 | 1.00 |

Fifty million regions is 3.40 km, not sub-kilometre. A one kilometre mean edge
needs about 577 million, and scaling the measured 10,000,005 run linearly that is
roughly 17 hours, 517 GB of resident memory and 271 GB of raw export. Fifty
million alone wants about 45 GB against a 61 GB machine. The first number that
matters here is that the target was an order of magnitude further away in N than
it looks, because edge length goes as the square root.

The second reason is the one that would hold on any hardware. Orogen's terrain
noise sits at FIXED physical wavelengths, and `terrain-config.js` states them:
`DETAIL_NOISE_FREQ` 5.0 is a 1280 km base and six octaves reach about 40 km, the
second pass doubles the frequency and reaches about 20 km, the regional term at
`DETAIL_NOISE_FREQ_MULT` 16 over four octaves reaches about 50 km, and the fine
term at `FINE_NOISE_FREQ_MULT` 20 over three octaves reaches about 80 km. The
finest feature the generator DESIGNS is therefore near 20 km, and the build's
15.19 km mesh already sits at its Nyquist.

The export confirms it. On `elevation_pre_erosion`, which is the designed
surface before any erosion, the semivariance gained per doubling of lag
collapses at short separation:

| lag | pre-erosion ratio | finished-terrain ratio |
| --- | --- | --- |
| 5 to 10 km | 1.081 | 1.389 |
| 10 to 15.19 km | 1.092 | 1.202 |
| 25 to 50 km | 1.257 | 1.369 |
| 50 to 100 km | 1.222 | 1.262 |

The designed terrain gains 8% per doubling below 10 km against 26% at 50 km: it
is smooth at those scales, which is what a 20 km noise floor looks like from
below. The finished terrain is not smooth there, and the difference between the
two columns is the whole of it. Every structure below about 15 km is EROSION
TEXTURE CUT AT THE MESH SCALE, carved by operators that work on the mesh graph
into a surface that carries no information at that wavelength.

So a finer mesh below the floor does not resolve terrain, it manufactures more
mesh-scale erosional texture whose placement is set by Voronoi jitter. That is a
property of the discretisation and not of the world, and no amount of it is
worth 17 hours.

Where to stop follows. At 7.60 km the mesh oversamples the 20 km floor by a
factor of 2.6 and the derived fields that motivated raising the count at all,
the basin catalogue and the drainage network, are resolved on it. 25,000,000
regions at 4.80 km is comfortable margin. Fifty million and beyond is spending
hours to sample a smooth surface more finely.

Genuine sub-kilometre relief has to be SYNTHESISED rather than generated: sub-grid
roughness conditioned on the resolved terrain and matched to the variogram
measured here, which is local, cheap, needs no global mesh, and is what GW-6 and
GRAV-6 have been asking for under the name downscaling. The floor measured in
this section is the specification such a synthesis has to meet at its upper end.
