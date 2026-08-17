# CLAUDE.md

Worldbuilding project for **Vesper**, a super-Earth around a mid-K dwarf. The
geography comes from a fork of World Orogen; ExoPlaSim, hydrography, pedology and
the biosphere are consumers of it, and more components are expected alongside
them.

This file is directives and design intent. Detail lives with the thing it
describes, and the pointers below are the map:

| Read this | For |
| --- | --- |
| `WORKFLOW.md` | what the components are, how they connect, the order they run in, and why that order is a loop. **Read it first.** |
| `source/README.md` | how to read an export: field conventions, the land-mask rule, the traps |
| `vendor/orogen/tools/README.md` | the authoritative export format |
| `<component>/README.md` | what that component does and how to run it |
| `notes/failure-modes.md` | how this project goes wrong, by class |
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
4. **After any patch, and after any `.venv` reinstall, rebuild every binary.**
   ExoPlaSim compiles one executable per (resolution, layers, ranks) triple, so
   a rebuild only refreshes the configuration you ran. `.venv` is untracked, and
   reinstalling it silently discards every applied patch.
   `python exoplasim/scripts/rebuild_binaries.py`, then `--verify`.
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
source/<build>/        World Orogen exports, one namespaced directory per build.
                       READ-ONLY. Add a build; never overwrite one. Only the
                       ACTIVE build has a payload; superseded ones are stubs in
                       archive/builds/.
exoplasim/             The ExoPlaSim climate component (see exoplasim/README.md).
hydrography/           Drainage, catchments, basin capacity (see hydrography/README.md).
pedology/              Soil formation, texture, phosphorus (see pedology/README.md).
biosphere/             LPJ-GUESS port and the C-N-P fork (see biosphere/README.md).
minerals/              Ore prospectivity, a field not deposits (see minerals/README.md).
maps/                  Rendering and cartography.
references/            Primary literature. PDFs untracked; INDEX.md tracked, and
                       it records which sources have actually been READ rather
                       than merely cited.
archive/               Identity of things whose payload has been deleted:
                       archive/runs/ for runs on superseded terrains,
                       archive/builds/ for the exports themselves. Tracked.
analysis/              Project-level analysis products (error budget, dust optics).
requirements.txt       Shared Python dependencies for .venv.
.venv/                 Python 3.12, already activated in this shell. Untracked,
                       and a reinstall silently discards every applied patch --
                       see the ExoPlaSim section.
```

`config/` and `source/` are project-level and shared. Component-specific work
lives under the component directory. New components get a sibling directory and
read the same `config/planet.yaml` and `source/`. `lib/` holds the shared
readers: `orogen.py` for the export, `gridding.py` for mesh-to-grid integration,
`orbit.py` for the orbital period. Reuse them; do not reimplement.

Git tracks the scripts, notes, configuration, and analysis products. It does
**not** track `exoplasim/runs/` (model output, 13 GB after triage; but
`exoplasim/runs/INDEX.json` **is** tracked, and since run ids are UUIDs it is the
only record of what each run was), the `source/` export
payloads, or `.venv/`. Those have no history to fall
back on, so be careful with destructive operations there; `.gitignore` says why
each is excluded.

## The upstream generator

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

`exoplasim/patches/exoplasim-3.4.2-star-cycle.patch` adds a sinusoidal stellar-flux
cycle to `radmod.f90`. `build_star_cycle_exoplasim.sh` applies it, verifies the
pinned upstream SHA, rebuilds, copies the result to
`exoplasim/inputs/exoplasim_cycle_t42/`, and reverses the patch on exit; the
vendored ExoPlaSim tree in `.venv` is left clean.

