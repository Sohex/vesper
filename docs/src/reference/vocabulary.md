# Vocabulary

These are the project's words. Use them and not synonyms.

**build** -- one World Orogen export, namespaced under `source/`. The ARTIFACT,
not the act that made it. Identified by its terrain hash, never by its name.
Superseded means wrong, not merely old.

**generation** -- the act. One pass of World Orogen, consuming the planet code,
the seed and a carve list, producing a build. It is the only thing in this
project that changes the terrain, and it is the only way a build comes to exist.

**bootstrap run** -- the first climate run on a build, made with only the surface
fields that are pure functions of terrain. It exists to produce the climatology
that the remaining fields need. **Its numbers are not the baseline.**

**baseline run** -- the second climate run on a build, on the full surface
fields: lakes, the lake compositing in the albedo, and soil water capacity, all
built from the bootstrap's climatology. The climatology of a baseline run is what
downstream components read.

**commissioning** -- the whole process of taking a build to a settled baseline:
terrain-only surface fields, bootstrap run, every derived field rebuilt from the
bootstrap's climatology, baseline run, convergence. A build is UNCOMMISSIONED
until that finishes, and an uncommissioned build has no climatology that anything
downstream may read.

  This is the word for "the work that extends from a build", and OROGEN IS NOT
  PART OF IT. If you are saying it and picturing a terrain being regenerated,
  you mean `iteration`.

**re-commissioning** -- commissioning the SAME build again, because something
that is not the terrain has invalidated its climatology: a model patch, a
configuration value, a physics correction. Same build, no generation, and the
bootstrap comes with it, because everything below a climatology is worthless
rather than stale once the climatology is, by rule 7 above.

**re-run the baseline** -- the narrow case, and NOT a re-commissioning: a new
baseline run on derived fields that are still valid, because whatever changed
does not reach them. The test is mechanical rather than a judgement: if lakes,
the albedo compositing, the soil or code 229 have to be rebuilt, their input is
a climatology, so it is a re-commissioning and the bootstrap is not optional.

  Do not say "re-baseline". It reads as an iteration and means one of these two,
  and that gap has already cost one misunderstanding.

**iteration** -- one turn of loop A in `docs/src/pipeline/sequencing.md`, and the expensive
thing: a generation on the carve list the previous commissioning produced, then
the commissioning of the build that comes out of it. An iteration therefore
CONTAINS a generation and a commissioning, in that order, and it is what
"starting again from Orogen" means.

  Say `generation` or `commissioning` when you mean one half. Saying `iteration`
  when you mean the second half is the collision this section exists for, and it
  is the easy one to make, because most of the work in a turn is the
  commissioning.

**segment** -- a contiguous block of orbits added to an existing run by
`continue_exoplasim.py`. Runs are made of segments; each records its I/O regime
and its purpose, and the purpose is DECLARED by the caller with `--purpose`
rather than inferred from the flags. `spinup` and `post_equilibrium_climatology`
are the run's own trajectory; `diagnostic` is orbits run to measure the model
rather than the planet, and those are kept out of convergence windows and
climatologies. See `exoplasim/scripts/segments.py`.

**carve verdict** -- the finding: which basins overflow, per basin, with its
evidence. `hydrography/analysis/carve_verdict.json`.

**marginal** -- a LANDFORM, and only that. A basin whose retain lands strictly
between 0 and 1: the outlet is notched but not cut through to the basin floor, so
Orogen produces a through-flowing valley with a residual lake in it. It is a
statement about what the terrain looks like.

**bracketed** -- a basin the two bounding climates DISAGREE about, from the
intersection carve in `docs/src/pipeline/loops.md`. It is a statement about our
uncertainty, not about a landform.

  **These two are not synonyms and were briefly the same word.** A basin can be
  marginal in both arms, bracketed while marginal in neither, both, or neither.
  "The marginal set" was used for the bracketed one in section 4 and leaked into
  the carve exporter, where `marginal` was already a verdict value written into
  the carve list Orogen reads -- so one word named a landform and an error bar in
  the same file. Say `bracketed` for the disagreement and reserve `marginal` for
  the notch.

  A third neighbour, kept distinct: `disputed` in `export_carve_list.py` is where
  the two EVAPORATION ESTIMATORS disagree about one climate. Three concepts,
  three words: estimator disagreement, climate disagreement, landform.

**carve list** -- the artifact Orogen consumes, `carve_list.txt`, one retain
fraction per basin. The verdict is a conclusion; the list is an instruction.

`docs/src/pipeline/costs.md` prices what a re-commissioning costs, row by row, and is
where to look before assuming which of these you are in.

