# CLAUDE.md

Worldbuilding project for **Vesper**, a super-Earth around a mid-K dwarf. The
geography comes from a fork of World Orogen; ExoPlaSim, hydrography, pedology and
the biosphere are consumers of it, and more components are expected alongside
them.

This file is directives and design intent. Detail lives with the thing it
describes, and the pointers below are the map:

| Read this | For |
| --- | --- |
| `WORKFLOW.md` | CANONICAL for pipeline REASONING: what the components are, how they connect, why the order is what it is, and why it is a loop. **Read it first.** |
| `config/pipeline.yaml` | CANONICAL for the pipeline GRAPH: every step, what it writes, what must precede it, its cost, and each loop's exit predicate. The two do not overlap; WORKFLOW references step ids from here |
| `scripts/pipeline.py` | `--status` what exists and what blocks the carve, `--plan <step>` the ordered steps to a target, `--register` the artifact table. Plans; never runs |
| `source/README.md` | how to read an export: field conventions, the land-mask rule, the traps |
| `vendor/orogen/tools/README.md` | the authoritative export format |
| `<component>/README.md` | what that component does and how to run it |
| `notes/failure-modes.md` | how this project goes wrong, by class |
| `notes/no-time-axis.md` | Orogen has no time axis. Read before asking any component for a duration, an age, or a rate |
| `notes/large-data.md` | batch, chunk, checkpoint, report. Required for any step whose input runs to GB |
| `notes/external-data.md` | routes into data this project does not generate; check the AWS Registry of Open Data before an API |
| `notes/audits/` | findings: what is true, with its evidence |
| `TASKS.md` | what to do about a finding, tracked atomically; closed ones move to `archive/tasks.md` |
| `world_state.json` | every current value |

## The rules that will bite you

These are the ones that have already cost real work. Each is enforced somewhere;
none of them is advice.

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
   ExoPlaSim compiles one executable per (resolution, layers, ranks) triple, so
   a rebuild only refreshes the configuration you ran, and the rest go stale
   silently. `python exoplasim/scripts/rebuild_binaries.py`, then `--verify`,
   which compares every executable against the sha of the source it was built
   from. This used to also warn that a `.venv` reinstall discarded the patch
   stack; it cannot any more, because the model is a subtree and a reinstall
   restores it rather than losing it.
5. **Pointing one component at another's output is deliberate, never
   defaulted.** Namespace per build and resolve with
   `builds.component_data(..., strict=True)`, or stamp `source_build` on the
   artifact and verify it with `lib/provenance.py:require_build` at the point of
   reading. That module says why; `check_consistency.py` enforces it.
6. **Generated runs get a UUID, never a derived name.** A parameter-built name
   separates runs only along the dimensions it encodes. ExoPlaSim and LPJ-GUESS
   both collided that way. Ask `exoplasim/runs/INDEX.json` what exists.
7. **`source/` is read-only. Add a build; never overwrite one.** But a build is
   DISPOSABLE until a climate run has consumed it: before that nothing depends
   on it, so the answer to a generator change is to regenerate rather than to
   migrate. After that it is not, because runs and verdicts start depending on
   it. What must survive either way is the recipe, in `source/README.md`.

   That is one case of the general rule, and the general rule is this. **What
   is durable is what the world is RECONSTRUCTED FROM; everything else is
   output.** The durable set is small and complete: the planet code and seed,
   the carve list, `config/planet.yaml`, the code, and the decisions and
   findings about MECHANISMS. Given those, every other artifact in this
   repository can be regenerated, which is why **an upstream change makes
   everything below it worthless rather than stale** -- there is nothing to
   reconcile, migrate or carry with a caveat, and the question after any change
   is "what is now worthless" rather than "what needs updating".

   The carve list is in the durable set because Orogen CONSUMES it, in the same
   way it consumes the seed, and not because incision is irreversible. Incision
   being irreversible is a fact about the world, not about the pipeline: a wrong
   verdict is recoverable, because a build is replaced wholesale rather than
   edited. Carving is an ordinary step that happens to be a convenient
   bottleneck, and treating it as sacred obscures the actual rule. Climatologies,
   run output, verdicts, soil, dust fields, prospectivity and `world_state.json`
   are all output. "The floor was binding on 73% of the overflowing set"
   survives every regeneration; "1,938 basins carve" does not.

8. **Before an expensive run**, and after changing `source_build`:

       python scripts/check_consistency.py     # do the artifacts agree?
       python scripts/smoke_test.py            # does the code that makes them?

9. **Read `notes/failure-modes.md`** before quoting a geography number, adding a
   component, or changing a quantity that more than one script consumes. The one
   most likely to catch you first: a pre-carve build is a *limit*, not a state,
   and pre-carve numbers are what is physically sitting in `source/` at the start
   of every cycle.

## Conventions

- **Do not write current values into prose.** `world_state.json` is generated
  from the artifacts and is the only place they belong. A number earns a place in
  a document only if it is a decision, a threshold, an identity, or if the
  magnitude carries an argument that fails without it: a convergence criterion, a
  terrain hash, "the two windows do not overlap", "exactly 180 degrees". Mean
  temperature, basin counts, land composition and the active build are none of
  those, and every one of them was wrong in these documents within a day of being
  written. `notes/` is the exception and the opposite: those are dated records of
  what was measured, so they keep their numbers and gain a "measured on".
- **Rewrite superseded content; do not mark it.** A reader grepping for a number
  lands on the number, not on the warning above it.
- Prose in docs and reports uses ASCII punctuation and avoids em dashes; match it.
- Keep citations and table cells on one source line, even where that breaks
  column alignment. They get copied out.
- Every run and analysis product records its provenance (config hash, input
  hashes, software versions) in JSON. Keep that up when adding steps.
- Claims about convergence and equilibration are stated with their exact criteria
  and are labelled honestly when they miss. One cold case is called
  "quasi-equilibrated" for missing its threshold by 0.004 W/m2. Preserve that
  standard rather than rounding results into passes.
- Scripts anchor their paths in a `_paths.py` and resolve from the file location,
  not the working directory, so they run from anywhere.
- **Physics is not a knob.** A process belongs in the model because it exists, not
  because including it improves a comparison. If adding a correct term makes an
  agreement worse, that is information about the implementation, the comparison,
  or a second error cancelling the first -- never a reason to remove the term.
  `notes/failure-modes.md` class 16.
- **Test the implementation against something that can fail, not the outcome
  against something that can only differ.** A check needs a right answer: an
  identity, a definition, a conservation law, or a quantity the other side
  already knows. Two valid formulations disagreeing is not evidence about
  either. If you cannot say in advance what result would mean "wrong", you
  have a number that will later be quoted as a validation, not a test.
  `notes/failure-modes.md` class 17.
- **An undocumented component is not complete.** Work here is picked up by
  someone with no memory of it -- assume a brick to the head between any two
  sessions, because a fresh session IS that. The test for done is not "does it
  work", it is "can someone who has never seen it find it and use it without
  reading the diff".

  What this does NOT mean is recording numbers. Anything that can be looked up
  or re-derived should be, and writing it into prose is the failure the first
  convention above exists to prevent. What it means is that the THING must be
  findable: a new module belongs in the `lib/` list, a new step in
  `config/pipeline.yaml`, a new component in the layout and in `WORKFLOW.md`, a
  model change as a commit under `vendor/exoplasim`, a new convention here. A
  component that works
  and is invisible will be reimplemented beside itself, which is how this project
  came to have four copies of a path resolver and three of a grid convention.

  The two rules pull in opposite directions and that is the point: **record what
  a thing IS and where it lives, never what it currently SAYS.**

- **Do not offer a defect as a decision.** A thing is a DECISION only if the
  project's declared truth does not already settle it. `config/planet.yaml`,
  the rules in this file, `WORKFLOW.md`'s ordering and the existing findings are
  declared truth; where they settle a question it is WORK, so do it and report
  it. It is a decision only where they conflict, or where it needs a preference
  or a threshold that nothing has fixed. The tell is unmistakable once you look
  for it: **if one of the options is "leave the known-wrong thing as it is",
  the question has been mis-framed.** Recorded 2026-08-17, after enabling a
  radiation weight derived for the wrong star and declaring the spectrum the
  config already named were both put to the user as choices. Neither was. By
  contrast the convergence criterion genuinely was one, because section 7's ban
  on retrofitting a criterion and conservation's verdict that the criterion read
  the wrong quantity pointed opposite ways and no document settled it.

- **`config/pipeline.yaml` is the pipeline graph, and an artifact that no step in
  it generates does not exist.** Every artifact this project keeps has to appear
  in `config/pipeline.yaml`, against the step that produces it. What reads it is
  DERIVED from the `needs` of other steps rather than declared, so there is one
  statement of each edge and not two. A script that writes something not in the register is either a
  missing step or a product nobody should be reading, and both are defects. This
  is not bookkeeping: it is what makes rule 7 usable, because "what is now
  worthless" can only be answered from a graph that is complete. When you add a
  generator, add its row in the same commit.

- **Findings and tasks are kept apart.** A document under `notes/audits/` says
  what is true and carries its evidence; `TASKS.md` says what to do about it and
  cites the document. That way a finding can be read without being re-litigated,
  and a task closed without editing the argument behind it.

## State, identity and history

`world_state.json` is what the project currently knows. It is generated by
`scripts/world_state.py`, never edited, and re-run after anything that changes a
build, a run or a verdict.

**It tracks the live thread only** -- the active build and the runs on it, not
the superseded builds or the runs belonging to them. A derived value whose
dependency has moved is *cleared*, not carried with a caveat: a stale derived
value is indistinguishable from a current one, so it emits an invalidation
record naming the dependency and what to run.

**Identity lives in the registry, state lives in world_state.** `lib/orogen.py`
keeps every build it has been checked against, so results already computed from a
superseded terrain stay readable and datable.

## Design intent

- **The orbit and the biosphere are one choice, not two.** The flux windows that
  put this world in its design temperature range do not overlap between a
  vegetated surface and a bare-rock one, and the endmember spread near the target
  is wider than the target band. No single flux is robust to the vegetation
  question. Chosen: vegetated. The band is narrow, so re-derive after anything
  that moves land albedo.
- **The pipeline is a loop, not a line.** Drainage depends on climate, climate on
  drainage, and both on the biosphere. `WORKFLOW.md` section 4 says why, and
  which loop is deliberately left open.
- **Do not reuse a sensitivity measured in one regime in another.** Bracket
  between two converged points that span the target rather than extrapolating
  from one. Doing the latter across the ice transition predicted 291.9 K for a
  run that converged at 287.47 K.
- **A correction's size depends on how much of the surface it acts on.** Estimate
  it against the state you are in, not the state it was first measured in.
- **A judgment made while a class is absent is not a judgment.** Values that
  nothing exposes go unchecked; three of them surfaced at once when a single
  tectonic bug was fixed. When a rule starts firing for the first time, audit
  everything it controls.

## Layout

```
config/planet.yaml     Canonical planet/star/orbit/atmosphere parameters. Project-level.
config/pipeline.yaml   Canonical pipeline GRAPH: steps, what each writes, what must
                       precede it, cost, and each loop's exit. Read by
                       scripts/pipeline.py and checked by smoke_test.py.
source/<build>/        World Orogen exports, one namespaced directory per build.
                       READ-ONLY. Add a build; never overwrite one. Only the
                       ACTIVE build has a payload; superseded ones are stubs in
                       archive/builds/.
exoplasim/             The ExoPlaSim climate component (see exoplasim/README.md).
hydrography/           Drainage, catchments, basin capacity (see hydrography/README.md).
pedology/              Soil formation, texture, phosphorus (see pedology/README.md).
biosphere/             LPJ-GUESS port and the C-N-P fork (see biosphere/README.md).
minerals/              Ore prospectivity, a field not deposits (see minerals/README.md).
aeolian/               Offline dust: emission, transport, deposition (see aeolian/README.md).
                       Reads a climatology and the lake solution; its deposition
                       feeds pedology, and its optical depth feeds the radiation.
maps/                  Rendering and cartography.
references/            Primary literature. PDFs untracked; INDEX.md tracked, and
                       it records which sources have actually been READ rather
                       than merely cited.
archive/               Identity of things whose payload has been deleted:
                       archive/runs/ for runs on superseded terrains,
                       archive/builds/ for the exports themselves. Tracked.
analysis/              Project-level analysis products (error budget, dust optics).
vendor/orogen/         World Orogen, a git subtree from the cf-fork branch of the
                       personal fork. Generates the geography.
vendor/exoplasim/      ExoPlaSim, a git subtree from the cf-fork branch of the
                       personal fork. THE model source: edited here, compiled
                       here, installed editable from here.
requirements.txt       Shared Python dependencies for .venv. ExoPlaSim is NOT here;
                       it is installed editable from vendor/exoplasim.
.venv/                 Python 3.12, already activated in this shell. Untracked,
                       and safe to reinstall: ExoPlaSim is installed editable
                       from vendor/exoplasim, so nothing is lost by rebuilding it.
```

`config/` and `source/` are project-level and shared. Component-specific work
lives under the component directory. New components get a sibling directory and
read the same `config/planet.yaml` and `source/`. `lib/` holds the shared
readers, and the list is the whole of it because "reuse them, do not reimplement"
is unusable if it names half: `orogen.py` for the export, `builds.py` for
resolving a build to a path, `paths.py` for repo-relative paths, `gridding.py`
for mesh-to-grid integration and the one grid convention, `orbit.py` for the
orbital period, `stellar.py` for the star's spectrum, band split and Rayleigh
coefficient, `sensitivity.py` for the one flux-to-kelvin conversion, and
`provenance.py` for build stamping and config drift. Reuse them; do not
reimplement. Rules 5 and 7 both cite modules from this list.

Git tracks the scripts, notes, configuration, and analysis products. It does
**not** track `exoplasim/runs/` (model output, 13 GB after triage; but
`exoplasim/runs/INDEX.json` **is** tracked, and since run ids are UUIDs it is the
only record of what each run was), the `source/` export
payloads, or `.venv/`. Those have no history to fall
back on, so be careful with destructive operations there; `.gitignore` says why
each is excluded.

## The vendored upstreams

Two of this project's components are other people's code, vendored as git
subtrees from personal forks so that a change to the component and the change to
whatever consumes it land in ONE commit, and so provenance is a commit in this
repository rather than the state of a directory outside it.

### World Orogen

The geography comes from a personal fork of World Orogen, vendored into this
repo at `vendor/orogen/` as a git subtree from the `cf-fork` branch of
`raguilar011095/planet_heightmap_generation`. It lives here so that a change to
the generator and the change to whatever consumes it land in ONE commit, and so
a build's provenance is a commit in this repository rather than the state of a
directory outside it. Pull upstream with `git subtree pull --prefix vendor/orogen
orogen-fork cf-fork --squash`. Its generated output stays untracked: the
subtree's own `.gitignore` excludes `out/`, which runs to 13 GB.

Its `tools/README.md` is the authoritative reference for the export format -- read it
before writing anything that consumes `source/`. The fork adds, over upstream:
lithology (rock class, erodibility, scarp potential), preserved endorheic basins,
a richer export manifest, non-Earth planet parameters, and direct emission onto
Gaussian (spectral) grids.

### ExoPlaSim

The climate model is a personal fork vendored at `vendor/exoplasim/`, a git
subtree from the `cf-fork` branch of `Sohex/ExoPlaSim`. Pull upstream with
`git subtree pull --prefix vendor/exoplasim exoplasim-fork cf-fork --squash`.

It is installed EDITABLE, so the source you read is the source that compiles and
the source that runs. Its build artifacts stay untracked: the subtree's own
`.gitignore` excludes everything `configure.sh` and `compile.sh` generate, which
is what keeps a five-binary rebuild from leaving the working tree dirty.

The fork carries, over upstream: the low-I/O restart and broadcast repairs, the
pyburn reader fix, the shortwave weights for a non-solar host, the dust and
aerosol stack, and the stellar cycle. `exoplasim/patches/README.md` says which
are upstream pull requests and which are ours to keep.

## Builds

`source/` is namespaced by build, because builds multiply faster than they can be
swapped in place and a result's provenance should depend on what it was computed
from, not on when. `config/planet.yaml` names the active one in `source_build`,
and `lib/builds.py` resolves it; nothing should hardcode a path under `source/`.

**No list of builds lives here, and no hash does either.** `lib/orogen.py` is the
registry: every build it has been checked against, keyed by
`manifest.hashes.finalElevation`, each with a note on what it is and what
superseded it. `world_state.json` records which one is active. The names are for
humans; the hash is the identity, and `lib/orogen.py` refuses a build it has not
been checked against.

A hash in prose is a current value wearing the costume of an identity. It reads
as durable, because a hash names a thing that never changes -- but *which* hash
is current changes every iteration, so the sentence around it goes quietly wrong
while still looking like a fact worth trusting. Cite the registry instead.

Basin ids are computed on the pre-conditioning surface, so they survive a carve
iteration and per-basin work generally resolves across builds. That is a
property of when the ids are assigned, not a promise -- `lib/orogen.py` keeps
the catalogue hashes it has verified, and that is the check.

Superseded means wrong, not merely old, and both are kept registered so results
already computed from them stay readable and datable. **Registered is not the
same as present**: only the active build has a payload under `source/`. The
others are identity stubs in `archive/builds/`, which is enough to recognise a
build and date a result. What each superseded build got wrong, and what replaced
it, is a note against its hash in `lib/orogen.py`.

A wrong verdict is recoverable and a wrong build is not the trap it sounds like:
Orogen regenerates terrain from the planet code plus a carve list in one pass,
so a build is replaced wholesale rather than edited.

## Environment

`.venv` is already active. Dependencies are pinned in `requirements.txt`; install
with `UV_CACHE_DIR=/tmp/world-uv-cache uv pip install --python .venv/bin/python -r requirements.txt`.
Building ExoPlaSim needs `gcc-fortran` and `openmpi` from the host (Arch).
Matplotlib is forced to `Agg` with its cache at `/tmp/world-matplotlib-cache`.

**ExoPlaSim is installed EDITABLE from its subtree, not from PyPI**, so a build
compiles in place and `requirements.txt` deliberately does not name it:

    uv pip install --python .venv/bin/python -e vendor/exoplasim

What the fork contains and where it came from is under "The vendored upstreams".

Two things about the model are worth knowing before you touch it. Several
changes are NO-OPS until a namelist key turns them on -- `h2osww` defaults to
1.0, `ndustrad` to 0, `nsolcycle` to 0 -- so a rebuilt binary reproduces the runs
that exist, and enabling one is a configuration decision that moves the mean.
And the low-I/O change ALTERS THE RESTART LAYOUT, so a run started before it
cannot be resumed by a binary built after it.

The stellar cycle is one of those switches rather than a separate build. There
is no cycle executable and no cycle tree: `nsolcycle` defaults to 0 and both
amplitudes to 0.0, which reduces the guard to `gsolinst = gsol0`, so the
ordinary binaries are bit-exact identical to an unpatched model until a cycle is
configured on.

