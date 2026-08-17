# Hydrography

Turns the World Orogen geography into drainage: where water goes, which basins
it collects in, and how much each can hold. Everything here except the lake
solver's forcing is climate-independent and can be built before ExoPlaSim runs.

```bash
python hydrography/scripts/build_hydrography.py   # ~13 s, climate-independent
python hydrography/scripts/lake_balance.py        # solver smoke test and sweep
python hydrography/scripts/surface_water.py       # ~3 s, needs a climatology
```

## Why this component exists

Upstream World Orogen forced every land cell to drain to the sea. The fork stops
doing that, so closed basins survive into the export, and it deliberately does
not decide what happens next: which basins hold lakes is a precipitation-versus-
evaporation balance, and that belongs downstream.

What arrives is therefore unrouted. `drain_to` is raw steepest descent, and
`drainage_terminal` is -2 for anything whose chain ends in neither the ocean nor
a preserved basin sink. On this planet that is **most of the land**, draining
into unpreserved pits, most of them a single mesh cell of noise -- the exact
share is per build and is in `world_state.json`. Integrating
precipitation over catchments without resolving that would silently discard most
of the land's water.

## What it does

**Resolves drainage** with a priority flood over the 2.5M-region mesh, filling
the noise pits while keeping the surviving basins as genuine terminals. How many
survive is a property of the build and is in `world_state.json`; it named a
count and a build here, and both were several iterations out of date.
Every land region ends up assigned to the world ocean or to exactly one basin.
112,217 regions get filled, by a median of 8.4 m, which is the scale that
confirms these were noise rather than landforms.

**Recomputes hypsometry on the finished terrain.** The catalogue's
`hypsometry` is measured on the natural, pre-conditioning surface and overstates
what the basins can hold: the finished terrain holds 61% of the catalogue's
figure. Level, area and volume curves are rebuilt from the flood
and are exact at their sample points, verified against brute force to 4e-15
relative on volume and exactly on area.

**Finds spill levels and targets** on the finished terrain. Most basins overflow
into another basin rather than to the ocean, so filling has to be solved as a
cascade rather than basin by basin.

**Couples to the climate grid.** `coupling_*.nc` gives each basin's catchment
area per grid cell as a sparse matrix, which is what you integrate P-E over.
Built for T42 and T85.

## Products

| File | Contents |
| --- | --- |
| `data/regions.nc` | per mesh region: resolved terminal, depression-filled surface |
| `data/basins.nc` | per basin: hypsometric curves, spill level and target, catchment area, capacity |
| `data/coupling_<grid>.nc` | sparse basin-by-grid-cell catchment areas |
| `data/hydrography_report.json` | diagnostics, river mouths, marginal seas, provenance |
| `data/surface_water.nc` | per region: lake, lake depth, river discharge; per basin: solved area, level, volume, overflow |
| `analysis/lake_balance_sweep.json` | solver sensitivity under placeholder forcing |
| `analysis/surface_water_report.json` | the solved water balance and its forcing |

## The lake solver

`lake_balance.py` is machinery, not a result. At equilibrium

    r (C - A) + p A = e A      so      A = r C / (e - p + r)

for catchment area C, lake area A, runoff depth r over dry land, and
precipitation p and evaporation e over open water. Overflow cascades downstream
and the system iterates to a fixed point.

Running it directly sweeps uniform placeholder forcing, which exercises the
solver and shows the sensitivity. Under 10 to 200 mm/yr of runoff against 400 to
1600 mm/yr of lake evaporation, lake area lands between 0.21% and 8.7% of the
planet. That range is a property of the terrain, not a prediction.

## Surface water under a real climate

`surface_water.py` supplies the forcing the solver was waiting for, from the
`climatology_s096` baseline, and accumulates the same water down the drainage
network to get rivers. This is the first thing in the project to decide
`surface_class == 2`, which World Orogen deliberately leaves empty.

Current figures -- how many basins hold water, lake area, largest river -- are in
`world_state.json` and in the run's own report, not here. They change with every
carve iteration and every climate re-baseline, and this document is about how the
solver works.

Against Earth, which is the only calibration available, the ratios are the part
worth carrying: this world has roughly twice the land at somewhat lower runoff
depth, so total river discharge is of the same order, and the largest river is
comparable to the Amazon. Lakes and endorheic drainage are both several times
Earth's share, and by a similar factor -- which is the consistency check, since
they are two views of the same aridity.
Those last two are computed by different routes and agreeing is a check rather
than a coincidence.

Two choices in the forcing are worth stating, because neither is the obvious
field:

**Runoff is P-E, not `mrro`.** At steady state they are the same quantity:
whatever falls on land and does not evaporate has to leave. The baseline run's
global water budget closes to 0.03% of the mean and the land's surplus matches
the sea's deficit to three decimals, so the budget is trustworthy. Its `mrro`
diagnostic accounts for only 15% of that surplus, 0.069 against 0.469 mm/day
over land, so the diagnostic is the unreliable one and P-E is what a water
balance can be built on. Both are recorded in the report.

**Lake evaporation is Penman, shared with the carve verdict.** Over land the
model's `evap` is actual evapotranspiration, throttled by soil moisture, which
is the wrong quantity for open water: using it would make every lake as dry as
the desert it sits in. `carve_verdict.penman_open_water` answers the right
question instead, and using the same function here means the lakes and the
terrain they sit in are judged by one rule.

It is also the estimate that can be checked. Applied to ocean cells, which
already are open water:

| mm/day, area-weighted over ocean | | ratio to the model |
| --- | ---: | ---: |
| the model's own evaporation | 3.672 | |
| Penman | 3.763 | 1.025 |
| Priestley-Taylor at alpha 1.26 | 3.957 | 1.078 |

(Recomputed here with cos-latitude weights; `carve_verdict.py` quotes 1.017 for
Penman from its own weighting, so the two agree to under a percent.)

A Priestley-Taylor estimate stood here first and was reported as 3.13 mm/day
"just below" the model, which was wrong twice over: the figure was an unweighted
cell mean compared against an area-weighted one, and correcting that puts PT
7.8% *above* the model rather than below it. Penman is three times closer on the
only ground truth available, and it was already written. Switching moved basins at spill from 683 to 770, so the direction of that error
was the opposite of what the first reading suggested. (Both figures predate the
longitude fix below and are quoted only to compare the two estimates.)

Penman is floored at the model's land rate, the same floor the verdict uses: it
linearises around air temperature and can otherwise fall below the model's own
evaporation where the ground runs hotter, which is impossible for a saturated
surface under the same forcing. One caveat stands: the ocean validation is where
air is near-saturated and wind is well resolved, which is the opposite of an
inland arid basin, so 2.5% is an upper bound on its accuracy in the places that
decide the verdict.

### The coupling matrix was being read 180 degrees out

Worth stating plainly because it invalidated a result that had already been
applied. `basin_means` indexed a climatology field with the coupling matrix's
own cell numbering, and the two number their columns differently: the coupling
inherits the Orogen grid, which runs -180 to 180, while an ExoPlaSim
climatology runs 0 to 360. Every basin was therefore reading the runoff and
evaporation of its antipode. Nothing about it was visible from the outside: the
array shapes match, latitude is unaffected, and the resulting fields are
plausible everywhere.

It was caught by asking the coupling to average a field of known longitude and
comparing the answer to where the basins actually are. 99.2% of basins came back
within 20 degrees of the antipode; the same test on latitude was correct. That
test is cheap and should be rerun whenever either grid changes.

The fix makes the convention explicit rather than assumed. `coupling_*.nc` now
carries the `cell_lat` and `cell_lon` it was built on, `basin_means` takes the
longitude axis of the fields it is given and remaps columns before indexing
anything, and it refuses to run against a coupling file too old to state its own
convention.

**This reaches further than the lakes.** `carve_verdict.py` shares
`basin_means`, so the iteration-1 verdict, the 1,522 basins carved to produce
`carved-zoned`, was decided on climate read from the wrong side of the planet.
That verdict has since been regenerated: 1,089 carve, 170 marginal, 2,370
preserve, with only 58.3% of verdicts unchanged and retain fractions correlating
at 0.267. 850 of the original carves were unjustified and 417 were missed. It is
applied in `carved-zoned-v4`, and the corrected carve set is visibly more
physical -- carved basins carry 5.7x the median catchment runoff of preserved
ones, where under the old verdict the two were nearly indistinguishable.
The carve pattern in the current terrain does not correspond to the climate that
was supposed to justify it. Regenerating that verdict is no longer optional
tidying before iteration 2; it is a correction.

Lakes are painted by filling each basin to its solved *area*, in ascending order
of the flooded surface, rather than by thresholding on the solved level. A mesh
region is 230 km2, so thresholding gives every basin with a smaller lake a 15 km
blob, and a wet landscape of many small basins comes out speckled with lakes
that are mostly rounding error. Filling by area reproduces the solved total to
0.4% and simply does not draw a lake too small to reach a whole region.

Rivers come from the priority flood's own discovery pointer, now kept as
`receiver` in `regions.nc`. Steepest descent on the filled surface would not do:
a filled depression is flat and has no downhill neighbour, so whole tributaries
would be dropped. Accumulation peels leaves off that tree rather than sorting on
elevation, for the same reason.

**Overflow is routed onto the mesh**, which needs the saddle and not just its
height. `spill_levels` was taking the minimum over a basin's exit edges and
keeping only the value; it now keeps the argument too, so `basins.nc` carries
`spill_region` and `spill_exit_region` either side of the saddle. Without them a
basin pinned at its spill had a known outflow and nowhere to start it, and the
largest flows on the planet were missing from the network: routing the 770
overflowing basins took the biggest river from 63,800 to **170,300 m3/s**.

There is no double counting, because the solver has already resolved the
cascade: a basin's overflow is its final equilibrium value with everything
upstream included, and the paths are disjoint segments, one basin's saddle to
the next one's sink.

## The Earth comparator, and a correction

Earth's endorheic fraction is **about one-fifth of its land surface**, not the
13% this document and others here previously claimed. Wang et al. (2018),
"Recent global decline in endorheic basin water storages", Nature Geoscience 11,
926-932, states it in its opening sentence and measures 31.8 million km2 across
48,813 landlocked watersheds delineated from 15-arcsecond HydroSHEDS.

The units match ours, which is the part worth checking. Wang's figure is the
**catchment area that drains internally**, which is the same measure as this
component's endorheic share of land -- not the basin-floor area, which is a
different quantity about three and a half times smaller and which this project
has already confused once.

So this world is more endorheic than Earth by roughly a factor of two after
carving, and closer to four before it. Still a large difference and still the
fact this component exists to handle, but a smaller multiple than the 13%
comparator implied.

## Headline finding, and its limit

A far larger share of this planet's land drains to a closed basin than Earth's
about one-fifth, by roughly a factor of two once carved. That follows from
the fork preserving closed basins instead of carving drainage to them, and it is
the fact this component exists to handle.

The endorheic share falls with each carve iteration, which is the verdict doing
what it exists to do. It is also the number most likely to be quoted from a stale
product: these files are per-build, and a figure computed against one terrain
says nothing about another.


**That figure is the hyper-arid limit, not a property of the world.** It assumes
no basin ever overflows. A basin that overflows year on year incises its outlet,
and over the 1e4 to 1e6 years a landscape needs to relax, that drains the lake
and the depression stops existing. So a basin pinned at its spill is a transient,
not a landscape state, and the endorheic share depends on how many basins the
climate keeps overflowing:

| (E-P)/runoff over the catchment | basins that carve | endorheic land |
| ---: | ---: | ---: |
| 0.8 (humid) | 3593 | 9.0% |
| 1 | 3515 | 22.8% |
| 2 | 2583 | 47.6% |
| 3 | 1753 | 59.0% |
| 5 | 857 | 66.6% |
| 10 (arid) | 237 | 72.3% |
| 50 (hyper-arid) | 9 | 76.0% |

**That table is from the pre-carve build and has not been regenerated.** It
still describes the shape of the dependence, which is the point it is making,
but its counts are the 3,629-basin set. Regenerating it means running
`carve_verdict.py` against the carved terrain, and that is an iteration-2 carve
decision rather than a documentation chore, so it is left for whoever takes that
decision.

The decision variable is climate-free per basin and is stored as
`critical_aridity_index` in `basins.nc`. A basin overflows exactly when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

so the index is pure geometry and the climate supplies one number per basin.
Median here is 3.03 on the carved terrain, meaning the typical basin overflows
unless evaporative demand over open water exceeds about three times the runoff
depth. Earth's
surviving endorheic basins sit well above that, which is why Earth keeps only
13%. `carve_verdict()` in `lake_balance.py` applies it.

## What this component cannot fix

**We cannot carve the terrain here.** Fixing a rim means changing elevation,
which belongs upstream in World Orogen: incision depends on the `erodibility`
field, and a second copy of the terrain in this component would desynchronise
ExoPlaSim's orography from the maps and from everything else reading `source/`.

What this component can do is decide *which* basins should carve, which is what
`carve_verdict.py` and `export_carve_list.py` are for. The generator takes it:
`--preserve-basins FILE` accepts a retain fraction per basin, where 1 keeps the
rim, 0 carves it and anything between cuts a notch and tapers it over the divide
band. `--carve-basins FILE` is subtractive. Both are documented in the fork's
`tools/README.md`.

**An earlier version of this section said that hook did not exist, and that was
wrong.** It described the state before iteration 1, which used exactly this
interface: the count of retain-0 entries in `carve_list.json` matches
`carvedByRetainZero` in the resulting build's manifest, which is the check that
the list was consumed as written. Iteration 2 is a run, not a
generator change.

The terrain still is not in equilibrium: many basins fill to their spill under
this climate, which is the water balance saying their outlets should have been
cut. They hold a large share of both the land and the lake area on the map, so
what is drawn is largely a landscape that has not relaxed yet. That is the case
for another carve iteration, and the counts are in `world_state.json` and
`surface_water.nc` rather than here, because they move with every pass.

## Known approximations

- **A few overflow paths disagree with the cascade.** Most overflowing basins
  have a saddle opening into exactly the basin the solver routes them to. A
  handful are cycle-collapsed: their target was reassigned to the cycle's
  primary while their saddle still leads where it physically leads, so the
  routed river and the solver's bookkeeping name different destinations.
- **The two grids are offset half a cell in longitude.** Orogen's exported grid
  centres its columns at -180 + (k + 0.5) dlon; an ExoPlaSim climatology centres
  its own at k dlon. Remapping by nearest column leaves up to 1.4 degrees of
  slip, about 185 km at the equator. Whether the boundary conditions handed to
  ExoPlaSim carry the same offset is an `exoplasim/` question and has not been
  checked here.
- **Catchment forcing is binned, not integrated.** A basin's runoff is a
  catchment-area-weighted mean over whole grid cells, so a catchment that
  straddles a cell gets that cell's value entire. Comparing a basin's solved
  overflow against the discharge the mesh delivers to its sink, the ratio has a
  median of 0.99, exactly as it should for a lake passing on everything but its
  own evaporation, and a ninetieth percentile of 1.13. That spread is this
  binning.
- **Merged basins.** 205 basins were collapsed into a neighbour to break spill
  cycles, which arise where basins share a saddle and each names the other as
  its outlet. Physically they merge into one lake at that level. The survivor's
  hypsometry still describes only its own depression, so a merged group filled
  past the shared saddle holds more water than the curve says.
- Catchments are assigned to grid cells by region centre, matching the
  exporter's categorical resampling rule, so a cell straddling a divide is
  attributed whole.
- The solver finds equilibrium, not a seasonal cycle. Basins with large storage
  relative to annual throughput will lag; nothing here models that.
