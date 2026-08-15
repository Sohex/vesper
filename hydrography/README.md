# Hydrography

Turns the World Orogen geography into drainage: where water goes, which basins
it collects in, and how much each can hold. Everything here except the lake
solver's forcing is climate-independent and can be built before ExoPlaSim runs.

```bash
python hydrography/scripts/build_hydrography.py   # ~7 s
python hydrography/scripts/lake_balance.py        # solver smoke test and sweep
```

## Why this component exists

Upstream World Orogen forced every land cell to drain to the sea. The fork stops
doing that, so closed basins survive into the export, and it deliberately does
not decide what happens next: which basins hold lakes is a precipitation-versus-
evaporation balance, and that belongs downstream.

What arrives is therefore unrouted. `drain_to` is raw steepest descent, and
`drainage_terminal` is -2 for anything whose chain ends in neither the ocean nor
a preserved basin sink. On this planet that is **63% of the land**, draining into
220,649 unpreserved pits, most of them a single mesh cell of noise. Integrating
precipitation over catchments without resolving that would silently discard most
of the land's water.

## What it does

**Resolves drainage** with a priority flood over the 2.5M-region mesh, filling
the noise pits while keeping the 3,629 preserved basins as genuine terminals.
Every land region ends up assigned to the world ocean or to exactly one basin.
99,085 regions get filled, by a median of 8.5 m, which is the scale that
confirms these were noise rather than landforms.

**Recomputes hypsometry on the finished terrain.** The catalogue's
`hypsometry` is measured on the natural, pre-conditioning surface and overstates
what the basins can hold: total capacity there is about 1.5x what the finished
terrain actually holds. Level, area and volume curves are rebuilt from the flood
and are exact at their sample points, verified against brute force to 4e-15
relative on volume and exactly on area.

**Finds spill levels and targets** on the finished terrain. 2,597 basins overflow
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
| `analysis/lake_balance_sweep.json` | solver sensitivity under placeholder forcing |

## The lake solver

`lake_balance.py` is machinery, not a result. At equilibrium

    r (C - A) + p A = e A      so      A = r C / (e - p + r)

for catchment area C, lake area A, runoff depth r over dry land, and
precipitation p and evaporation e over open water. Overflow cascades downstream
and the system iterates to a fixed point.

Running it directly sweeps uniform placeholder forcing, which exercises the
solver and shows the sensitivity. Under 10 to 200 mm/yr of runoff against 400 to
1600 mm/yr of lake evaporation, lake area lands between 0.21% and 8.7% of the
planet. That range is a property of the terrain, not a prediction. The real
answer needs ExoPlaSim runoff and evaporation integrated through `coupling_*.nc`.

## Headline finding, and its limit

76% of the land drains to a closed basin rather than to the sea, against roughly
13% on Earth. That follows from the fork preserving closed basins instead of
carving drainage to them.

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

The decision variable is climate-free per basin and is stored as
`critical_aridity_index` in `basins.nc`. A basin overflows exactly when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

so the index is pure geometry and the climate supplies one number per basin.
Median here is 2.92, meaning the typical basin overflows unless evaporative
demand over open water exceeds about three times the runoff depth. Earth's
surviving endorheic basins sit well above that, which is why Earth keeps only
13%. `carve_verdict()` in `lake_balance.py` applies it.

## What this component cannot fix

**The terrain has not been carved and we cannot carve it here.** Hydraulic,
thermal and glacial erosion all acted on this surface, but drainage-enforcement
carving was disabled globally so that basins would survive to export. The result
is a landscape whose rims were never cut by the outflow that the water balance
says crosses them. For the basins that overflow, the terrain is not in
equilibrium with any climate.

Fixing that means changing elevation, which belongs upstream in World Orogen,
not here: incision depends on the `erodibility` field, and a second copy of the
terrain in this component would desynchronise ExoPlaSim's orography from the
maps and from everything else reading `source/`.

The generator does not currently expose the hook. `--preserve-basin ID` adds a
basin to the preserved set and `--no-basins` disables preservation globally, but
there is no way to say "carve these specific basins". Preservation retains about
80% of natural spill depth and carving retains about 3%, so the two settings are
all-or-nothing where what the water balance produces is a per-basin verdict.
Closing the loop needs something like `--carve-basin ID` or `--preserve-only`.

## Known approximations## Known approximations

- **Merged basins.** 323 basins were collapsed into a neighbour to break spill
  cycles, which arise where basins share a saddle and each names the other as
  its outlet. Physically they merge into one lake at that level. The survivor's
  hypsometry still describes only its own depression, so a merged group filled
  past the shared saddle holds more water than the curve says.
- Catchments are assigned to grid cells by region centre, matching the
  exporter's categorical resampling rule, so a cell straddling a divide is
  attributed whole.
- The solver finds equilibrium, not a seasonal cycle. Basins with large storage
  relative to annual throughput will lag; nothing here models that.
