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

**Price an orbit in WALL CLOCK, not in model time.** What a run waits on is the
model plus postprocessing plus the ocean and ice stream writes, and quoting the
model half alone understates a run by about a third. Per-segment model-only
figures are `native_runtime_seconds` in each run's manifest; the wall figure is
the gap between consecutive `MOST_REST.NNNNN`. The current per-orbit cost, and what the
pyburn fixes did to it, are measured in
`notes/audits/pyburn-postprocessing-cost.md`. A settling run off a
near-equilibrium restart plus ten clean orbits is hours, not minutes, and a
commissioning usually wants two of them. An
ITERATION is far more, because the generation in front of it regenerates the
terrain and every field that is a pure function of it.

**That affordability assumes the pyburn fixes are resident.** Unpatched,
postprocessing is roughly three quarters of the orbit rather than an overhead
on it, which triples a commissioning. The same audit note has the
measurements.
