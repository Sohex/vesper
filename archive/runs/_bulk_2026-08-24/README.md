# The 179 runs on disk when the project moved to the threaded build

`INDEX_AT_DELETION.json` is `exoplasim/runs/INDEX.json` as it stood at commit
`cf02ec01`, which is every row this project had: 179 runs, 224.8 GB.

**This archive is thinner than the ones beside it, and the reason is a mistake.**
The others carry each run's manifest, convergence assessment, climate series and
namelists, because those were copied out before the payload went. These carry
only the INDEX row, because the directories were removed by an `rm -rf` whose
argument came from an unguarded command substitution that expanded to nothing --
`rm -rf exoplasim/runs/$(...)` with the `$(...)` empty is `rm -rf
exoplasim/runs/`. The removal was going to happen; curating it first was not
optional and did not happen.

WHAT SURVIVES ELSEWHERE, and is the reason this is recoverable at all:
`exoplasim/analysis/` is a sibling of `runs/` and was untouched, so the
convergence assessments, the climatologies, the energy closures and the energy
term budgets are all still on disk. Those are the derived products anything
downstream ever read. What is gone is the raw NetCDF and the per-run manifests.

WHY IT IS NOT A LOSS OF RESULTS. `CLAUDE.md` rule 7: until the canonical
climatology lineage is declared, every build, run and climatology is disposable
whatever has consumed it. It is not declared. `config/planet.yaml` carried
`baseline_climatology: null`, so no component was reading a run. Every number
this project rests on is in the notes, and every one of these runs was on the
MPI transform path that `world-bdh` retired.

The lesson is in `docs/src/practice/failure-modes.md`; the operational form of
it is that `runs/` has no git history to fall back on and a destructive command
there takes a literal path or nothing.
