# The components, the data flow, and the register

## 1. The components

```
config/planet.yaml   Canonical planet, star, orbit, atmosphere. Every component reads it.
source/              World Orogen exports. Canonical, read-only.
lib/                 Shared readers: orogen.py, gridding.py, and the rest of the
                     list CLAUDE.md's Layout carries in full.
hydrography/         Drainage, catchments, basin capacity, lake balance, carve verdict,
                     and the steady-state water table, whose depth field is
                     UNVALIDATED: see hydrography/README.md before using it.
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
                     One UUID-named frame per state of the world under
                     data/<build>/, taken after every step the map is reachable
                     from, so a pass can be watched changing rather than only
                     its end state drawn.
```

## 2. The data flow

```
World Orogen (fork)
   |  seed + planet code -> terrain, lithology, closed basins, hydrology
   v
source/<build>/exoplasim-T21|T42|T85|T127|T170 + grid-512x256
   |  one native raw/ mesh is stored in the T42 export as a carrier and is
   |  identical beneath all five Gaussian grids; new support audits use the
   |  ~10M-region, 7.60 km reference export
   v
lib/gridding.py            integrate mesh fields onto any model grid
   |                        |
   |                        v
   |                  exoplasim/scripts/build_boundary_conditions.py -> land mask, topography
   |                  exoplasim/scripts/build_surface_albedo.py      -> albedo, forest fraction
   |                  exoplasim/scripts/build_surface_soil_water.py  -> dwmax  (off for the bootstrap)
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
hydrography/            pedology/scripts/build_soil.py        hydrography/scripts/carve_verdict.py
  drainage,               weathers lithology            which basins overflow
  catchments,             under the climate         <-- integrates climate over
  hypsometry,               |                           catchments
  coupling matrix           v
   |                  biosphere/scripts/build_lpj_driver.py -> one binary, N years
   |                        |
   |                        v
   |                  biosphere/scripts/run_lpj_guess.py  (LPJ-GUESS, MPI)
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
[section 4](loops.md) says why each has to be what it is.

**Several branches hang off it that the diagram does not draw** -- the offline
aerosols (dust, sea salt, volcanic sulfate) and the minerals overlay. All are
in the register; the two largest are sketched below.

```
climatology + surface_water.nc -> aeolian/scripts/build_dust.py -> dust_baseline.{nc,json}
   |                                        |
   |   deposition -> pedology (loess, phosphorus)
   |   optical depth -> exoplasim/scripts/dust_optics.py -> dust_aerofile.py -> build_surface_dust.py (1811) -> the model
   |                 -> exoplasim/scripts/dust_forcing.py -> carve_verdict.py --dust-forcing
   v
source/ + soil + drainage -> minerals/scripts/build_prospectivity.py            (terrain only)
                          -> minerals/scripts/build_downstream_prospectivity.py (carries a climate)
```

## 2b. The artifact register

**The register is `config/pipeline.yaml`.**

    python scripts/pipeline.py --register     # artifact, step, and what reads it
    python scripts/pipeline.py --status       # what is present, and what is not
    python scripts/pipeline.py --plan <step>  # the ordered steps to reach a target
    python scripts/pipeline.py --purge <step> # what a change to that step makes worthless

Consumers are derived rather than declared: what reads an artifact is the set
of steps that `need` the step which writes it, plus the steps marked
`reads_export` for the export itself. Declaring both directions would be the
same duplication one level down.

`--plan` carries a `MAP:` line after every step the map is reachable from,
which is most of loop A: a carve verdict changes the next terrain and the next
terrain is the picture. That set is derived from the same graph rather than
listed, and `maps/snapshot.py` keys a frame on the identity of the four inputs
a map is drawn from -- so a step that moves none of them is recorded against
the frame that already exists and costs no render. The result is one frame per
state of the world, which is what lets a pass be watched changing rather than
only its end state drawn. `--no-maps` drops the lines.

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
