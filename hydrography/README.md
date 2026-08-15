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

## Headline finding

**76% of this planet's land is endorheic**, against roughly 13% on Earth. That
is a direct consequence of the fork preserving closed basins instead of carving
drainage to the sea, and it is the single most important thing this component
says about the world. Rivers reaching the ocean are the exception here.

Two caveats on that number. It describes drainage on the dry terrain: once
basins fill and spill, their catchments join whatever is downstream, so the
effective endorheic share falls with a wetter climate. And it counts area
draining to a basin, not area under water; how much becomes lake is what the
solver decides.

## Known approximations

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
