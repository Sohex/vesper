# CLAUDE.md

Worldbuilding project for **Vesper**, a super-Earth around a mid-K dwarf. The
geography comes from a fork of World Orogen; ExoPlaSim, hydrography, pedology and
the biosphere are consumers of it, and more components are expected alongside
them.

**This file is a map: identity, pointers, and the rules that bite. Nothing
else.** A lesson, an argument, an incident or a post-mortem does not land here;
it lands in the docs book or in `notes/`, and earns at most a one-line
directive or a pointer row below. An addition that is not a pointer or an
invariant is misfiled by its shape alone, wherever it came from.

| Read this | For |
| --- | --- |
| `docs/src/SUMMARY.md` | the docs book. Its pipeline chapters are CANONICAL for pipeline REASONING: what the components are, how they connect, why the order is what it is, and why it is a loop. **Read them before pipeline work** |
| `config/pipeline.yaml` | CANONICAL for the pipeline GRAPH: every step, what it writes, what must precede it, its cost, and each loop's exit predicate. The two do not overlap; the pipeline chapters reference step ids from here |
| `scripts/pipeline.py` | `--status` what exists, `--plan <step>` the ordered steps to a target, `--register` the artifact table, `--purge <step>` everything a change to that step makes worthless. Plans and never RUNS a step; `--purge` deletes, and is a dry run until `--execute` |
| `scripts/link_worktree.py` | a worktree sees the tracked tree and NOTHING else: no export payloads, no references, no run output, no `.venv`. Run it first in a new worktree and never hand-link one; `--check` says whether one is complete. What it deliberately does not link, and why, is in `docs/src/reference/environment.md` |
| `source/README.md` | how to read an export: field conventions, the land-mask rule, the traps |
| `vendor/orogen/tools/README.md` | the authoritative export format |
| `<component>/README.md` | what that component does and how to run it |
| `docs/src/practice/working-agreements.md` | how a session conducts itself, with the incidents that earned each agreement |
| `docs/src/practice/conventions.md` | how documents, numbers and results are kept, argued in full |
| `docs/src/practice/failure-modes.md` | how this project goes wrong, by class |
| `docs/src/reference/builds.md` | builds, identity, the registry, and the durable set |
| `docs/src/reference/vendored-upstreams.md` | the three subtrees: what each fork carries, how to pull upstream, the branch conventions |
| `docs/src/reference/design-intent.md` | the standing decisions that shape the world |
| `docs/src/reference/environment.md` | install commands, host packages, the model facts to know before touching it, and the command-line tools for reading an artifact, with the three traps that come of their assuming Earth |
| `docs/src/reference/no-time-axis.md` | Orogen has no time axis. Read before asking any component for a duration, an age, or a rate |
| `docs/src/reference/large-data.md` | batch, chunk, checkpoint, report. Required for any step whose input runs to GB |
| `docs/src/reference/external-data.md` | routes into data this project does not generate; check the AWS Registry of Open Data before an API |
| `notes/audits/` | findings: what is true, with its evidence |
| `bd` (beads) | what to do about a finding, tracked atomically. `bd ready` for what is unblocked, `bd show <id>` for one issue, `bd list --label clim` for an area. Every issue cites the document that justifies it |
| `world_state.json` | every current value |

## The rules that will bite you

These are the ones that have already cost real work. Each is enforced somewhere;
none of them is advice. The numbers are stable and cited from code and docs --
do not renumber.

1. **Take land from `surface_class`, never from `land_mask`.** They disagree over
   dry closed-basin floor below sea level, and preserving that terrain is the
   entire point of the fork. Do not reconstruct it as `land_mask | is_endorheic`
   either -- that misses the depressions too small to enter the basin catalogue.
   See `source/README.md`.
2. **`config/planet.yaml` is the single source of truth for the planet.** Do not
   copy its values into prose or into a script. Gravity is declared and Orogen's
   is canonical; mass follows from it and `derive()` raises if they disagree.
3. **Never match by longitude across the export/ExoPlaSim boundary.** The two
   label the same grid differently. Map by index, share one coordinate source,
   and never reconstruct one. This has silently matched zero cells on three
   separate scripts.
4. **After any change under `vendor/exoplasim`, rebuild every binary.**
   ExoPlaSim compiles one executable per (resolution, layers, ranks)
   configuration -- the grid and the thread count are `parameter`s compiled in,
   not settings a binary reads -- so a rebuild only refreshes the configuration
   you ran, and the rest go stale silently. `python exoplasim/scripts/rebuild_binaries.py`, then `--verify`,
   which compares every executable against the sha of the source it was built
   from.
5. **Pointing one component at another's output is deliberate, never
   defaulted.** Namespace per build and resolve with
   `builds.component_data(..., strict=True)`, or stamp `source_build` on the
   artifact and verify it with `lib/provenance.py:require_build` at the point of
   reading. That module says why; `check_consistency.py` enforces it.
6. **Generated runs get a UUID, never a derived name.** A parameter-built name
   separates runs only along the dimensions it encodes. ExoPlaSim and LPJ-GUESS
   both collided that way. Ask `exoplasim/runs/INDEX.json` what exists.
7. **`source/` is read-only: add a build, never overwrite one. What is durable
   is what the world is RECONSTRUCTED FROM; everything else is output.** An
   upstream change therefore makes everything below it worthless rather than
   stale, and the question after any change is "what is now worthless", never
   "what needs updating". The durable set is small and complete: the planet
   code and seed, the carve list, `config/planet.yaml`, the code, the build
   recipe in `source/README.md`, and the decisions and findings about
   MECHANISMS. A build is DISPOSABLE until a climate run has consumed it --
   and **until the canonical climatology lineage is declared, so is every run,
   climatology and build that has been consumed.** A defect in the model is
   therefore never weighed against the cost of the output it invalidates.
   `docs/src/reference/builds.md` has the argument.
8. **Before an expensive run**, and after changing `source_build`:

       python scripts/check_consistency.py     # do the artifacts agree?
       python scripts/smoke_test.py            # does the code that makes them?

9. **Read `docs/src/practice/failure-modes.md`** before quoting a geography
   number, adding a component, changing a quantity that more than one script
   consumes, or reaching for a measurement to justify a design decision. The one
   most likely to catch you first: a pre-carve build is a *limit*, not a state,
   and pre-carve numbers are what is physically sitting in `source/` at the
   start of every cycle.

## Conventions

The arguments and the incidents behind these are in
`docs/src/practice/conventions.md`; the directives are here.

- **Do not write current values into prose.** `world_state.json` is the only
  place they belong. A number earns a place in a document only as a decision, a
  threshold, an identity, or a magnitude an argument fails without. `notes/` is
  the exception and the opposite: dated records keep their numbers and gain a
  "measured on".
- **Rewrite superseded content; do not mark it.** A reader grepping for a number
  lands on the number, not on the warning above it.
- **A sentence earns its place by the future work it can inform.** Process
  narration, correction stories, dates that carry no identity, and commentary
  about the document itself are byproduct: delete them. Git carries how a
  document came to be right; the document carries what is right.
  `docs/src/practice/failure-modes.md` class 25.
- Prose in docs and reports uses ASCII punctuation and avoids em dashes; match it.
- **Write about the simulation in the simulation's terms.** Name the
  simulated subject explicitly -- the model's snow albedo, a gridcell's
  vegetation -- so no sentence parses as being about the real world or about
  people; keep established technical terms whose context already
  disambiguates them; and any document that can be read standalone states the
  worldbuilding frame up front, as `biosphere/notes/fire-model-audit.md` does.
  THIS RULE IS SUBJECT TO ITSELF: do not enumerate the phrases it exists to
  avoid, here or anywhere.
  The argument and the history are in `docs/src/practice/conventions.md`.
- Keep citations and table cells on one source line, even where that breaks
  column alignment. They get copied out.
- Every run and analysis product records its provenance (config hash, input
  hashes, software versions) in JSON. Keep that up when adding steps.
- **One host, many agents: claim it before timing anything.**
  `python scripts/machine.py --check` before starting heavy work,
  `--claim "<purpose>" --minutes N` before taking a wall-clock, throughput or
  per-step cost number, `--release` after. `--check` reads the PROCESS TABLE as
  well as the claim, because work that is insensitive to load still generates
  it and therefore never claims: a checker trusting the claim alone sees an
  empty file, calls the host quiet, and is wrong. A timing is a measurement of a
  MACHINE STATE as much as of a model, and this project runs many agents on one
  host, so a number taken while someone else is integrating is not a slow
  number, it is a number of a different experiment. Record the load beside any
  timing you keep; one without the machine state it was taken under cannot be
  compared against a later one. The claim file is deliberately outside the
  repository, because each fan-out agent works in its own worktree and an
  in-tree file would coordinate nothing. Where you can, price work in something
  the scheduler cannot move -- retired instructions under
  `OMP_WAIT_POLICY=passive` survive contention that wall clock does not.
- **A thread team's working set on one die targets 32 MB**, counting one copy
  per thread for anything threadprivate. Above it is a regression even when
  this machine gets faster: 32 MB is CCD1 here and is what a part without
  stacked cache is likely to have. `docs/src/reference/environment.md`.
- Thresholds are fixed before results are seen; a criterion chosen after the
  run it judges is not a criterion.
- Estimates that cannot be verified are bracketed rather than guessed, and the
  bracket is reported.
- Claims are checked against the artifact rather than the documentation.
- **Check the instrument against the size of the effect before believing a
  number.** Work out what the effect is worth in the units the instrument
  reports and compare that with its own scatter; if the effect is smaller, the
  number is noise however tidy it looks. A bed shorter than its startup, a sweep
  whose range never reaches the effect, and a refusal taken on a profile the
  model no longer has all return ordinary-looking numbers.
  `docs/src/practice/failure-modes.md` class 34.
- Convergence claims state their exact criteria and are labelled honestly when
  they miss. Preserve that standard rather than rounding results into passes.
- Scripts anchor their paths in a `_paths.py` and resolve from the file
  location, not the working directory, so they run from anywhere.
- **Physics is not a knob.** A process belongs in the model because it exists,
  not because including it improves a comparison. A correct term that worsens
  an agreement is information -- never a reason to remove the term.
  `docs/src/practice/failure-modes.md` class 16.
- **Test the implementation against something that can fail, not the outcome
  against something that can only differ.** A check needs a right answer: an
  identity, a conservation law, a quantity the other side already knows. If you
  cannot say in advance what result would mean "wrong", it is not a test.
  `docs/src/practice/failure-modes.md` class 17.
- **An undocumented component is not complete.** The test for done is "can
  someone who has never seen it find it and use it without reading the diff":
  a new module belongs in the `lib/` list, a new step in
  `config/pipeline.yaml`, a new component in the layout and in
  `docs/src/pipeline/components.md`, a new convention here. Record what a thing
  IS and where it lives, never what it currently SAYS.
- **`config/pipeline.yaml` is the pipeline graph, and an artifact that no step
  in it generates does not exist.** Consumers are derived from the `needs` of
  other steps, so each edge is stated once. When you add a generator, add its
  row in the same commit; rule 7 is only answerable from a graph that is
  complete.
- **Findings and tasks are kept apart.** A document under `notes/audits/` says
  what is true and carries its evidence; the `bd` issue tracker says what to do
  about it and cites the document.
- **No poison seeds.** Do not print what a rule forbids, announce tensions
  between rules, or leave untracked "worth checking" loops: state the
  positive form, the decision procedure, and an issue.
  `docs/src/practice/failure-modes.md` class 26.

## Working agreements

One line each; the argument and the recorded incidents are in
`docs/src/practice/working-agreements.md`.

- A task has TWO end states, completed or blocked; anything else is mid-task,
  however well it is reported.
- Every open thread gets its own verdict, named separately; one blocked item
  does not cover the one beside it.
- Naming the next step is not taking it: if you can say what would settle a
  question, you are not blocked.
- Ask whether the thing should exist BEFORE building or fixing it: deleting
  is a fix, so ask what breaks if the thing simply goes.
- Do not offer a defect as a decision: if one option is "leave the known-wrong
  thing as it is", the question is mis-framed.
- Stale derived artifacts are the resting state of the tree, not a work list;
  regenerate when a step you are running needs them, and otherwise leave them.

## Vocabulary

These are the project's words; use them and not synonyms. One line each here --
the definitions, the collisions and their cost history are in
`docs/src/reference/vocabulary.md`.

- **build** -- one Orogen export, namespaced under `source/`. The artifact,
  identified by terrain hash. Superseded means wrong, not merely old.
- **generation** -- one Orogen pass; the only thing that changes the terrain.
- **bootstrap run** -- the first climate run on a build, on terrain-only
  fields. Its numbers are NOT the baseline.
- **baseline run** -- the second, on the full surface fields; downstream
  components read its climatology.
- **commissioning** -- taking a build to a settled baseline. Orogen is NOT
  part of it.
- **re-commissioning** -- commissioning the SAME build again, bootstrap
  included.
- **re-run the baseline** -- a new baseline on still-valid derived fields.
  Do not say "re-baseline"; it is ambiguous between the two above.
- **iteration** -- one turn of loop A: a generation, then the commissioning of
  the build it produced.
- **segment** -- a contiguous block of orbits added to a run; its purpose is
  DECLARED by the caller, and diagnostics stay out of climatologies.
- **canonical climatology lineage** -- the run chain the world's numbers will
  finally rest on. It DOES NOT EXIST YET; until it is declared, every build, run
  and climatology is disposable whatever has consumed it.
- **carve verdict / carve list** -- the finding with its evidence / the
  instruction Orogen consumes.
- **marginal / bracketed / disputed** -- a landform / a disagreement between
  the two bounding climates / a disagreement between the two evaporation
  estimators. Three concepts, three words, never synonyms.

## Layout

```
config/planet.yaml     Canonical planet/star/orbit/atmosphere parameters. Project-level.
config/pipeline.yaml   Canonical pipeline GRAPH: steps, what each writes, what must
                       precede it, cost, and each loop's exit. Read by
                       scripts/pipeline.py and checked by smoke_test.py.
source/<build>/        World Orogen exports, one namespaced directory per build.
                       READ-ONLY. Only the ACTIVE build has a payload; superseded
                       ones are stubs in archive/builds/.
exoplasim/             The ExoPlaSim climate component (see exoplasim/README.md).
hydrography/           Drainage, catchments, basin capacity (see hydrography/README.md).
pedology/              Soil formation, texture, phosphorus (see pedology/README.md).
biosphere/             LPJ-GUESS port and the C-N-P fork (see biosphere/README.md).
minerals/              Ore prospectivity, a field not deposits (see minerals/README.md).
aeolian/               Offline dust: emission, transport, deposition (see aeolian/README.md).
maps/                  Rendering and cartography.
analysis/              Project-level analysis products (error budget, dust optics).
references/            Primary literature. PDFs untracked; INDEX.md tracked, and it
                       records which sources have actually been READ.
archive/               Identity of things whose payload has been deleted. Tracked.
docs/                  The docs book (mdBook; docs/src/SUMMARY.md is the index).
notes/                 Dated findings with their evidence; notes/audits/ for audits.
vendor/orogen/         World Orogen fork, git subtree. Generates the geography.
vendor/exoplasim/      ExoPlaSim fork, git subtree. THE model source: edited here,
                       compiled here, installed editable from here.
vendor/lpj-guess/      LPJ-GUESS CNP fork, git subtree. The Vesper input and
                       calendar port compile directly from this tree.
vendor/cgenie/         cGENIE.muffin, git subtree. The CANDIDATE offline ocean under
                       OCN-3, not an adopted component: it builds and runs here, and
                       nothing reads it. See vendored-upstreams.md.
vendor/lpjml/          LPJmL, git subtree. NOT the biosphere and NOT a replacement for
                       lpj-guess: a different model in the same family, vendored because
                       the intended work is fork-shaped. Nothing reads it yet.
.venv/                 Python 3.14. Untracked, and NEVER tracked: a symlink named
                       .venv was committed once and checked out in the main repo it
                       points at itself, which destroys the environment. See Environment.
```

New components get a sibling directory and read the same `config/planet.yaml`
and `source/`. `lib/` holds the shared readers, and the list is the whole of it
because "reuse them, do not reimplement" is unusable if it names half:
`orogen.py` (the export, and the build registry), `builds.py` (build to path),
`rungs.py` (the resolution ladder: rung to grid dimensions, which FFT module a
rung transforms on, the check that config's `resolution`, `latitudes` and
`longitudes` are one fact, the check that every restatement of the ladder
elsewhere in the tree agrees with it, and the THREE timestep quantities kept
apart -- each rung's measured stability ceiling, the escalation route's step per
rung, and what a commissioning-length run has shown about a pair),
`paths.py` (repo-relative paths, the climatology resolver, and its clean-I/O
and configured-grid guards), `gridding.py` (mesh-to-grid, the one grid
convention, and the reduction operators by field semantics: extensive,
intensive, categorical, moments, expectation, and the ledger of what a
reduction dropped), `orbit.py` (orbital period), `stellar.py` (spectrum, band split,
Rayleigh coefficient), `sensitivity.py` (the one flux-to-kelvin conversion),
`climatology.py` (time-bin weights),
`autocorrelation.py` (the integrated autocorrelation time, the effective
sample size, and the ONE standard error of a mean over a series whose samples
carry memory), `run_lengths.py` (how many orbits to BUY, derived from the
timescales: a commissioning span from the MEMORY time and a settling block from
the RELAXATION time, both bracketed because neither time is measured, and the
two never interchanged), `lapse.py` (lapse rates, and the
height of the lowest model level),
`surface_classes.py` (derived surface classes BY NAME),
`sea_water.py` (the four numbers salinity reaches the model through, read
from `icemod.f90` and the run's namelist rather than copied), `provenance.py` (build
stamping, config drift, the hashes of the derived files a generator read, and
`staged_surface_field`, the one door onto a staged `.sra` keyed by the rung). Rule 5 cites `builds.py` and `provenance.py` from this list.

Git does not track `exoplasim/runs/` (model output; but `runs/INDEX.json` IS
tracked, and since run ids are UUIDs it is the only record of what each run
was), the `source/` payloads, or `.venv/`. Those have no history to fall back
on, so be careful with destructive operations there; `.gitignore` says why
each is excluded.

## Environment

`.venv` is already active. ExoPlaSim is installed EDITABLE from
`vendor/exoplasim` -- the source you read is the source that compiles and runs
-- and is deliberately absent from `requirements.txt`. Model changes are
no-ops until a namelist key enables them, and the low-I/O change alters the
restart layout. `docs/src/reference/environment.md` has the install commands,
the host packages, and the argument behind each of those facts.


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:970c3bf2 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   bd dolt push
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->
