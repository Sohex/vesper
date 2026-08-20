# The components, the data flow, and the register

## 1. The components

```
config/planet.yaml   Canonical planet, star, orbit, atmosphere. Every component reads it.
source/              World Orogen exports. Canonical, read-only.
lib/                 Shared readers: orogen.py (the export), gridding.py (mesh to grid).
hydrography/         Drainage, catchments, basin capacity, lake balance, carve verdict.
exoplasim/           Boundary conditions, climate integrations, climatology.
pedology/            Weathers lithology into soil, and into solute fluxes: CO2, silica, phosphorus.
biosphere/           LPJ-GUESS: vegetation, leaf area, carbon, PFT composition.
minerals/            Ore prospectivity per deposit type. Reads the export; never feeds climate.
aeolian/             Offline dust: emission, transport, deposition, optical depth.
                     Reads a climatology and the lake solution; feeds pedology's
                     loess and the radiation's dust forcing. NOT a leaf: it is
                     downstream of one climate and upstream of the next.
analysis/            Project-level products that belong to no single component:
                     the error budget, the dust optics and forcing.
maps/                Rendering. Terminal: nothing reads its output. Its generators
                     live in the component root rather than a scripts/ directory,
                     which is why an audit that globbed */scripts missed them.
```

## 2. The data flow

```
World Orogen (fork)
   |  seed + planet code -> terrain, lithology, closed basins, hydrology
   v
source/<build>/exoplasim-T21|T42|T85 + grid-512x256 + maps
   |  raw/ mesh lives in the T42 export and is identical for all of them
   v
lib/gridding.py            integrate mesh fields onto any model grid
   |                        |
   |                        v
   |                  exoplasim/build_boundary_conditions.py -> land mask, topography
   |                  exoplasim/build_surface_albedo.py      -> albedo, forest fraction
   |                  exoplasim/build_surface_soil_water.py  -> dwmax  (off by default)
   |                        |
   |                        v
   |                  ExoPlaSim spin-up
   |                        |
   |                  build_climatology.py -> averaged climatology
   |                                       -> per-orbit climatologies (--per-year)
   |                                       -> climate series, per bin per orbit
   |                        |
   |          +-------------+-------------+
   |          |                           |
   v          v                           v
hydrography/            pedology/build_soil.py        hydrography/carve_verdict.py
  drainage,               weathers lithology            which basins overflow
  catchments,             under the climate         <-- integrates climate over
  hypsometry,               |                           catchments
  coupling matrix           v
   |                  biosphere/build_lpj_driver.py -> one binary, N years
   |                        |
   |                        v
   |                  biosphere/run_lpj_guess.py  (LPJ-GUESS, MPI)
   |                        |
   |                        +--> cpool.out ---> back to build_soil.py
   |                        |                   (soil and biosphere iterate)
   |                        v
   |                  fpc.out -> build_surface_albedo.py --mode modelled
   |                        |     albedo and forest fraction from what grew
   |                        v
   |                  back to ExoPlaSim
   v
carve list -> back to World Orogen -> new terrain
```

Three loops close in that diagram and a fourth is cut across iterations;
section 4 says why each has to be what it is.

**Two branches hang off it that the diagram does not draw**, because they would
turn one picture into four. Both are in the register.

```
climatology + surface_water.nc -> aeolian/build_dust.py -> dust_baseline.{nc,json}
   |                                        |
   |   deposition -> pedology (loess, phosphorus)
   |   optical depth -> exoplasim/dust_optics.py -> dust_aerofile.py -> the model
   |                 -> exoplasim/dust_forcing.py -> carve_verdict.py --dust
   v
source/ + soil + drainage -> minerals/build_prospectivity.py            (terrain only)
                          -> minerals/build_downstream_prospectivity.py (carries a climate)
```

## 2b. The artifact register

**The register is `config/pipeline.yaml`, and it is not reproduced here** --
one graph in one place, referenced from the other.

    python scripts/pipeline.py --register     # artifact, step, and what reads it
    python scripts/pipeline.py --status       # what is present, and what is not
    python scripts/pipeline.py --plan <step>  # the ordered steps to reach a target
    python scripts/pipeline.py --purge <step> # what a change to that step makes worthless

Consumers are derived rather than declared: what reads an artifact is the set
of steps that `need` the step which writes it, plus the steps marked
`reads_export` for the export itself. Declaring both directions would be the
same duplication one level down.

`pipeline.py` plans and never RUNS A STEP. Several steps in the graph cost
hours, and a script that could start one by accident is worse than no script.
It also does not check whether artifacts AGREE -- that is
`check_consistency.py`, and asking either to do the other's job would give two
answers to one question.

`--purge` is the exception and it only ever deletes; it is a dry run until
`--execute`. Rule 7 makes "what is now worthless" the question after any
change, and that is a graph question -- doing it by hand produces a keep-list
rather than an answer. Two properties to know before using it. It never
crosses `orogen`: loop A is a cycle, so unrestricted reachability from any
climate step would come back round and offer to delete the terrain, which is
wrong rather than merely alarming -- a new climatology does not invalidate the
build it was computed on. And it does not touch climate RUNS, which are
UUID-named and have no path in the graph; it names them and points at
`scripts/archive_runs.py`, which extracts a run's identity and verifies the
archive before deleting anything.
