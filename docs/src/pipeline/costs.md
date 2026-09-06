# 0. What a re-commissioning costs

**The words are in [the vocabulary](../reference/vocabulary.md).** The rows
below are a RE-COMMISSIONING: the same build,
everything below the climatology redone. A bare re-run of the baseline is
cheaper still, and is only honest when nothing in the "yes" rows is reachable
from whatever changed.

| | redone? | why |
| --- | --- | --- |
| Orogen, `source/` | no | the terrain has not changed |
| drainage, basins, coupling | no | pure functions of terrain |
| land mask, topography, roughness | no | pure functions of terrain |
| model binaries | only if it was a patch | see CLAUDE.md rule 4 |
| the climate run | yes | that is the point |
| climatology | yes | follows the run |
| lakes, albedo 174-176, soil, code 229 | yes | all are built from a climatology |
| a second run on those rebuilt fields | usually | this is where the loop closes |
| carve verdict, dust, everything downstream | yes | below the climatology |

**Price an orbit in WALL CLOCK, not in model time, and always AT A RUNG.** What
a run waits on is the model plus postprocessing plus the ocean and ice stream
writes, and quoting the model half alone understates a run by about a third; a
per-orbit figure carrying no rung is worth nothing at all, since the rungs are a
factor of seven apart in `notes/audits/resolution-ladder-wall-clock.md`, which
is where the per-rung price lives. The wall figure for a run already on disk is
the gap between consecutive `MOST_REST.NNNNN`, or each segment's own
`started_utc` to `finished_utc` in the run manifest. Per-segment model-only
figures are `native_runtime_seconds` in the same manifest, and they are
CPU-SECONDS OVER THE WHOLE THREAD TEAM rather than wall seconds -- the model
reports `cpu_time` across the process, which was one MPI rank and is now the
whole team -- so divide by the thread count before comparing one with a wall
figure. On T21 at sixteen threads the two differ by about fourteen. The pyburn
fixes and what they did to the postprocessing share are measured in
`notes/audits/pyburn-postprocessing-cost.md`. A settling run off a
near-equilibrium restart plus ten clean orbits is hours, not minutes, and a
commissioning usually wants two of them. An
ITERATION is far more, because the generation in front of it regenerates the
terrain and every field that is a pure function of it.

**That affordability assumes the pyburn fixes are resident.** Unpatched,
postprocessing is roughly three quarters of the orbit rather than an overhead
on it, which triples a commissioning. The same audit note has the
measurements.
