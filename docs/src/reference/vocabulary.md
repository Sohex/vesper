# Vocabulary

These are the project's words. Use them and not synonyms.

**build** -- one World Orogen export, namespaced under `source/`. The ARTIFACT,
not the act that made it. Identified by its terrain hash, never by its name.
Superseded means wrong, not merely old.

**base / interim** -- a BASE is the geography loop A iterates FROM, and a base is
INTERIM while something outstanding would change the base itself. The test is not
"will this build be regenerated" -- every base is, repeatedly, and that is the
loop working. It is "would the thing outstanding change what the base IS".

The distinction is easy to get backwards and was, in this project, in the naming
of `interim-10m-base`. The two-pass recipe in `source/README.md` regenerates a
new base with an ice mask once it has a climatology, and every carve regenerates
it again -- but those are ITERATIONS. They consume the base and belong to the
lineage when they execute; they are steps in `config/pipeline.yaml` rather than
defects in their input. What makes a base interim is something like LITH-26,
which would change the lithology the base is generated WITH, so a build made
before it lands is a different starting geography rather than an earlier turn of
the same one.

A corollary worth keeping: a defect in a process the base does not run cannot
make the base interim. GRAV-6 is the case -- glacial erosion carries no gravity
term, and a base generated at `--glacial 0` runs no glacial erosion at all
(`terrain-post.js` gates it on `iter < gIters`), so the defect reaches the
ice-mask pass and not the base underneath it.

Do NOT call a build `canonical` on the strength of having no outstanding trigger
today. `canonical` is reserved below and names a LINEAGE that does not exist yet.

**generation** -- the act. One pass of World Orogen, consuming the planet code,
the seed and a carve list, producing a build. It is the only thing in this
project that changes the terrain, and it is the only way a build comes to exist.

**bootstrap run** -- the first climate run on a build, made with only the surface
fields that are pure functions of terrain. It exists to produce the climatology
that the remaining fields need. **Its numbers are not the baseline.**

**baseline run** -- the second climate run on a build, on the full surface
fields: lakes, the lake compositing in the albedo, and soil water capacity. On
the FIRST pass those are built from the bootstrap's climatology, because it is
the only one that exists. The climatology of a baseline run is what downstream
components read.

**commissioning** -- the whole process of taking a build to a settled baseline:
terrain-only surface fields, bootstrap run, every derived field rebuilt from the
BEST AVAILABLE climatology, baseline run, convergence. A build is UNCOMMISSIONED
until that finishes, and an uncommissioned build has no climatology that anything
downstream may read.

  Best available, not the bootstrap's, and the difference only shows on a
  second pass. A derived field whose answer depends on the climate STATE reads
  the baseline once one exists and the bootstrap before that, which on a first
  pass through a build is the bootstrap in every case. So the first pass reads
  exactly as it always did, and it is the pass after -- the one that re-runs
  the baseline on rebuilt fields -- where a field pinned to the bootstrap would
  be holding the loop at the earliest stage rather than the best determined
  one. `lib/paths.py:best_available_climatology` resolves it and stamps which
  stage each product was built at; `docs/src/pipeline/loops.md` argues it. A
  field that depends on the model calendar or the grid rather than the state
  reads the bootstrap at every pass, because the two climatologies carry the
  same answer.

  This is the word for "the work that extends from a build", and OROGEN IS NOT
  PART OF IT. If you are saying it and picturing a terrain being regenerated,
  you mean `iteration`.

**re-commissioning** -- commissioning the SAME build again, because something
that is not the terrain has invalidated its climatology: a model patch, a
configuration value, a physics correction. Same build, no generation, and the
bootstrap comes with it, because everything below a climatology is worthless
rather than stale once the climatology is, by `CLAUDE.md` rule 7.

**re-run the baseline** -- a new baseline run on the same build, and NOT a
re-commissioning. It covers two cases: derived fields that are still valid
because whatever changed does not reach them, and derived fields that have been
rebuilt on the standing baseline, which is one more turn of the loop rather
than a return to its start.

  The test is mechanical rather than a judgement, and it is asked of the
  CLIMATOLOGY and not of the fields. If what changed invalidates the
  climatology -- a model patch, a configuration value the run reads, a physics
  correction -- everything below it is worthless by `CLAUDE.md` rule 7, the
  bootstrap included, and it is a re-commissioning. If what changed invalidates
  only a derived field, that field is rebuilt from the standing baseline and
  the baseline is run again on it; no bootstrap is needed, because the
  climatology the field reads is still valid.

  Do not say "re-baseline". It reads as an iteration and means one of these two,
  and that gap has already cost one misunderstanding.

**iteration** -- one turn of loop A in `docs/src/pipeline/sequencing.md`, and the expensive
thing: a generation on the carve list the previous commissioning produced, then
the commissioning of the build that comes out of it. An iteration therefore
CONTAINS a generation and a commissioning, in that order, and it is what
"starting again from Orogen" means.

  Before writing `iteration`, check that you mean a generation plus its
  commissioning; if you mean only the second half, write `commissioning`. Most
  of a turn's work is the commissioning, which is why the check is worth
  making.

**segment** -- a contiguous block of orbits added to an existing run by
`continue_exoplasim.py`. Runs are made of segments; each records its I/O regime
and its purpose, and the purpose is DECLARED by the caller with `--purpose`
rather than inferred from the flags. `spinup` and `post_equilibrium_climatology`
are the run's own trajectory; `diagnostic` is orbits run to measure the model
rather than the planet, and those are kept out of convergence windows and
climatologies. See `exoplasim/scripts/segments.py`.

**canonical climatology lineage** -- the chain of runs the world's published
numbers will finally rest on. **It does not exist yet.** Until it is declared,
every build, every run and every climatology is disposable no matter what has
consumed it, and a defect found in the model is never weighed against the cost of
the output it invalidates. Declaring it is a decision, recorded in
`docs/src/reference/builds.md`, and it is the point from which output starts
being worth keeping.

  The word exists because "baseline" was doing this job badly. A baseline run is
  a POSITION IN ONE BUILD'S COMMISSIONING; the lineage is a statement about the
  whole project having stopped moving underneath it. A build can have a perfectly
  good baseline and still be disposable, and while the model is being corrected
  every one of them is.

**carve verdict** -- the finding: which basins overflow, per basin, with its
evidence. `hydrography/analysis/carve_verdict.json`.

**marginal** -- a LANDFORM, and only that. A basin whose retain lands strictly
between 0 and 1: the outlet is notched but not cut through to the basin floor, so
Orogen produces a through-flowing valley with a residual lake in it. It is a
statement about what the terrain looks like.

**bracketed** -- a basin the two bounding climates DISAGREE about, from the
intersection carve in `docs/src/pipeline/loops.md`. It is a statement about our
uncertainty, not about a landform.

  These two are not synonyms and were briefly the same word: a basin can be
  marginal in both arms, bracketed while marginal in neither, both, or
  neither, and the conflation once put a landform and an error bar under one
  name in the same file. **Say `bracketed` for the disagreement and reserve
  `marginal` for the notch.**

  A third neighbour, kept distinct: `disputed` in `export_carve_list.py` is where
  the two EVAPORATION ESTIMATORS disagree about one climate. Three concepts,
  three words: estimator disagreement, climate disagreement, landform.

**carve list** -- the artifact Orogen consumes, `carve_list.txt`, one retain
fraction per basin. The verdict is a conclusion; the list is an instruction.

**tuned / opaque / implicit-Earth** -- the three ways a constant can lack a
derivation a reader can use. They are three different defects with three
different repairs, and the words are not interchangeable.

**tuned** -- the number IS the residual of a fit. There is no derivation at all:
it was moved until a comparison came out, and what it records is the answer it
was fitted to rather than anything about the world. `CLAUDE.md`'s "No tuned
values" forbids these, and the disposition space is exactly four -- source it,
derive it, declare it with a bracket that gets swept, or record it as
irreducible with the argument. `notes/audits/tuned-values.md` enumerates them.

**opaque** -- a derivation EXISTS and sits somewhere a reader of the code cannot
reach: a vendored default, an upstream paper nothing cites, a value transcribed
through two secondaries. The number may be perfectly good. What is missing is the
route to it, so the repair is to walk the chain and record it rather than to
change the value. `notes/audits/opaque-constants.md`.

**implicit-Earth** -- the derivation is sound and is FOR THE WRONG PLANET. This
one is a DIAGNOSIS and never an endorsement, and it is the one most often read
backwards. Establishing that a constant is implicit-Earth is what proves it is
wrong HERE and names why; it is a reason to replace the value, not a reason to
keep it. `notes/audits/inherited-earth-constants.md` and
`notes/audits/ocean-tier-implicit-earth.md`.

  The classes are ordered by how much is known, not by how bad they are, and
  moving a constant between them is progress even when the number does not
  change. `EXOPLASIM_DZ0LAND_M = 2.0` is the worked example: it was carried as
  tuned, and WORLD-U8DS proved it is Earth's own area-averaged land roughness,
  recoverable by reducing the vendored Earth boundary dataset's code 173 over its
  own land. That moved it from tuned to implicit-Earth without moving the number
  by a millimetre -- and the value of doing so is that "arbitrary" and "provably
  Earth's" call for different work. The second is replaceable by derivation on
  this planet's own relief, lithology and land cover; the first gave nothing to
  derive from.

  **A number's class is not a verdict on the number.** Do not read "it is only
  implicit-Earth" as a defence, and do not read a reclassification as a repair.
  The repair is a value this planet's own physics produces, with Earth's figure
  kept as a COMPARISON to report the distance from rather than a target to solve
  onto. Solving onto it reports that distance as zero by construction, which is
  the failure the class exists to make visible.

`docs/src/pipeline/costs.md` prices what a re-commissioning costs, row by row, and is
where to look before assuming which of these you are in.

