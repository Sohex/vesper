# The vendored upstreams

Four of this project's components are other people's code, vendored as git
subtrees so that a change to the component and the change to whatever consumes
it land together, and so provenance is a commit in this repository rather than
the state of a directory outside it. Each is a maintained fork rather than a
pinned dependency; `docs/src/reference/design-intent.md` says why that is the
expected end state for an external model here and how to weigh a candidate
against it.

## World Orogen

The geography comes from a personal fork of World Orogen, vendored into this
repo at `vendor/orogen/` as a git subtree from the `cf-fork` branch of
`raguilar011095/planet_heightmap_generation`. Pull upstream with `git subtree pull --prefix vendor/orogen
orogen-fork cf-fork --squash`. Its generated output stays untracked: the
subtree's own `.gitignore` excludes `out/`, which runs to GB.

Its `tools/README.md` is the authoritative reference for the export format -- read it
before writing anything that consumes `source/`. The fork adds, over upstream:
lithology (rock class, erodibility, scarp potential), preserved endorheic basins,
a richer export manifest, non-Earth planet parameters, and direct emission onto
Gaussian (spectral) grids.

## ExoPlaSim

The climate model is a personal fork vendored at `vendor/exoplasim/`, a git
subtree from the `master` branch of `Sohex/ExoPlaSim`. Pull upstream with
`git subtree pull --prefix vendor/exoplasim exoplasim-fork master --squash`.

The branch name differs from Orogen's deliberately. There, the fork's master
tracks upstream and `cf-fork` is a separate line of work; here the fork's master
IS the integrated line, so a second branch would be a copy with no owner. The
consequence to know: GitHub defaults the head of a new pull request from a fork
to that fork's default branch, so a PR against upstream must name its head
branch explicitly or it will offer the whole stack.

It is installed EDITABLE, so the source you read is the source that compiles and
the source that runs. Its build artifacts stay untracked: the subtree's own
`.gitignore` excludes everything `configure.sh` and `compile.sh` generate, which
is what keeps a five-binary rebuild from leaving the working tree dirty.

The fork carries, over upstream: the low-I/O restart and broadcast repairs, the
pyburn reader fix, the shortwave weights for a non-solar host, the dust and
aerosol stack, and the stellar cycle. `exoplasim/patches/README.md` says which
are upstream pull requests and which are ours to keep.

## LPJ-GUESS CNP

The terrestrial biosphere model is vendored at `vendor/lpj-guess/` from
`mateusdp/LPJ-GUESS-NTD`, pinned to tag `LPJ-GUESS-CNP_v1.0` and commit
`b368b893c4324840b43c56866e901c6916858afb`. Pull that exact upstream state
with:

    git subtree pull --prefix vendor/lpj-guess \
      https://github.com/mateusdp/LPJ-GUESS-NTD.git \
      LPJ-GUESS-CNP_v1.0 --squash

The tagged fork omits the MPL-2.0 licence files present in the LPJ-GUESS 4.1.1
distribution. The subtree restores those files verbatim from the md5-verified
Zenodo 8065737 tarball; an upstream pull must not delete them.

The subtree also carries the Vesper calendar and astronomy port and the
`vesperinput` module directly. There is no second source copy and no patch stack:
the source under `vendor/lpj-guess/` is what is compiled. `framework/vesper.h`
is the exception only because it is generated from `config/planet.yaml` before
each build and therefore ignored.

Vendoring the CNP source does not itself enable phosphorus limitation.
`data/ins/global.ins` and the run harness keep `ifplim 0` until the gridded
weathering, sorption, deposition, and replacement productivity prediction are
ready to land as one scientific change.

## cGENIE

The candidate offline ocean is vendored at `vendor/cgenie/`, a git subtree from
the `master` branch of `derpycode/cgenie.muffin`, MIT licensed. Pull upstream
with `git subtree pull --prefix vendor/cgenie cgenie-fork master --squash`.

The remote is `cgenie-fork`, pointing at `Sohex/cgenie.muffin`, which is a fork
of `derpycode/cgenie.muffin`. A fork is expected rather than optional here,
because OCN-19 and OCN-20 are fork-shaped by construction: one deletes dead
duplication and threads BIOGEM's tracer loops, the other replaces the barotropic
streamfunction solve. As with ExoPlaSim, the fork's `master` is the integrated
line, so upstream is pulled into the fork first and this subtree pulls from the
fork.

That is also why it is a subtree and not an extraction under `references/`.
External source this project reads and will not edit is held there instead, and
`references/INDEX.md` records which trees those are and why.

**It does not build where it stands.** `genie-main/user.mak` sets
`GENIE_ROOT = $(HOME)/cgenie.muffin` and `RUNTIME_ROOT = ../../cgenie.muffin`,
so the tree expects to sit at `~/cgenie.muffin`. Repointing those is part of
OCN-3's "whether it builds here" and is the first fork change anyone will make.
Run output goes to `OUT_DIR = $(HOME)/cgenie_output`, outside this repository,
which is why the ignore rules here cover only objects, archives and the
executable.

**What it carries that open rows already name.** `genie-goldstein` is the
frictional-geostrophic ocean OCN-19 and OCN-20 are about, and its `invert.f`
ordering is what OCN-20 quotes. `genie-knowngood` ships the four reference
configurations OCN-19 validates against. `genie-ecogem` is the ecosystem tier
OCN-4 names beside MARBL. And `genie-plasim` is not a coupling module but a
whole PlaSim: a strict SUBSET of `vendor/exoplasim`'s module set -- no `p_exo`,
no aerosol core, no MPI, no `hurricanemod` or `glaciermod` -- plus one file,
`geniemod.f90`, 71 lines of pure declaration. So this subtree carries the Earth
model ExoPlaSim's exoplanet capability was added to, and the entire
atmosphere-ocean coupling surface is one array-declaration module. It bears on
OCN-17's evaluation of the published ExoPlaSim-to-cGENIE path, and on OCN-10,
whose forcing contract can be read off it rather than derived. Note that the two
published paths are different arrangements: the offline regrid route moves wind
stress, winds and albedo, while this in-tree one hands over seventeen fields.
`notes/external-model-survey.md` section 30.

**IT IS NON-DIMENSIONALISED AGAINST EARTH'S RADIUS AND GRAVITY, and half of the
scales are configurable so the other half is easy to miss.**
`genie-goldstein/src/fortran/initialise_goldstein.F:376-391` hardcodes
`rsc = 6.37e6` and `gsc = 9.81` and exposes neither, while `sodaylen`,
`sidaylen`, `yearlen`, `nyear` and the depth scale `par_dsc` are all in
`ini_gold_nml`. Four derived scales then carry the hardcoded pair everywhere:
`tsc = rsc/usc` sets every non-dimensional time, `rhosc` carries both, `opsisc`
scales the reported overturning streamfunction and `rfluxsc` the heat flux.

At this planet's 1.20 Earth radii and 1.306 Earth gravity, three of those are
wrong by a factor of 1.20 and `rhosc` by about eight percent. **Nothing will
warn.** The model is entirely non-dimensional, so an overturning reported
through `opsisc` comes back twenty percent low, in the right units, looking
exactly like an answer.

Rotation is the one planetary constant the source does parameterise, and it
FAILS OPEN. `fsc` takes `4*pi/sidaylen` only when the solar and sidereal day
lengths differ by more than 0.001, and otherwise reverts to Earth's
`2*7.2921e-5` for backwards compatibility. A configuration that sets the two
equal gets Earth's Coriolis scaling silently. Set them to this planet's values
and check `fsc` took the intended branch before believing any run.

None of this is an argument against the model, and OCN-12 owns the full
inventory. It is here because it is the thing most likely to produce a
plausible wrong number for somebody who assumed a namelist covered the planet.
`notes/external-model-survey.md` section 10e has the derivation and the
precedent, including a scale that was already wrong once and corrected outside
the code.

**`genie-paleo`'s plots are deleted from this checkout, and the rest is kept.**
That directory arrived at 1.4 GB, of which 678 PostScript plots of Earth
paleogeographic reconstructions were 1290 MB. They are pictures, nothing here
will read them, and removing them took the subtree from 1.6 GB to 289 MB.

What was KEPT is the other 49 MB and the reason is specific: 486 base
configurations under `genie-main/configs` point into those directories through
`ea_1`, `go_1`, `gs_1` and `bg_par_pindir_name`, and the `.k1`, `.paths`,
`.psiles` and `.dat` files they resolve to are the format reference for what a
cGENIE bathymetry configuration IS. OCN-11 has to produce one from the Orogen
export, and OCN-18 asks specifically what muffingen's `.k1`, `.paths` and
`.psiles` generation does above 36 x 36. Deleting the whole directory would have
cost that and saved 49 MB.

The deletion is IN THE FORK, not only in this checkout, so a `git subtree pull`
will not reintroduce the plots. `git subtree split --prefix=vendor/cgenie`
reconstructs real upstream lineage rather than rooting at the squash commit, so
the split branch was a descendant of `008fd490` and pushed to the fork's master
as a fast-forward carrying only the deletion. That is the general recipe for
sending a change up: split, push the split branch to `cgenie-fork master`,
delete the split branch.

One consequence remains. The history here still carries the deleted blobs, so
`.git` did not shrink and will not; the saving is working-tree size, which is
what `grep` and `find` pay.

Unrelated and worth knowing before someone blames the deletion: four directories
that configs reference, `fkl_np10`, `fkl_pp01_DH`, `fm0450ab` and `wppcont1`,
are absent upstream and were never in the vendored commit.

## LPJmL

LPJmL is vendored at `vendor/lpjml/`, a git subtree squashed from the `master`
branch of `PIK-LPJmL/LPJmL` at commit `572e2b906ac2c55b2ee6661a93e4633b126254e4`.
Pull upstream with:

    git subtree pull --prefix vendor/lpjml \
      https://github.com/PIK-LPJmL/LPJmL.git master --squash

**It is a subtree and not a `references/` extraction for the usual reason**: the
work intended on it is fork-shaped and it is planned into the pipeline. It was
first pulled down as read-only comparison material for SPITFIRE, its
process-based fire scheme, because FIRE has ten open rows of ten issued and the
only scheme consulted was BLAZE, which arrived with LPJ-GUESS. That reading is
still worth doing and is now incidental to why it is here.

**It is NOT the biosphere, and it does not supersede LPJ-GUESS.** This project's
terrestrial biosphere is the LPJ-GUESS CNP fork above, and nothing in that
changes. LPJmL is a different model in the same LPJ family -- same lineage,
different code base, different scope -- and having both vendored is deliberate
rather than duplication. Read `vendor/lpj-guess/` when the question is about this
project's biosphere; read `vendor/lpjml/` when the question is about LPJmL.

**Nothing reads it yet**, exactly as with cGENIE: it is vendored source with no
pipeline step, so `config/pipeline.yaml` has no row for it and rule 7's "what is
now worthless" question does not reach it. When a step does consume it, that step
and its artifact go into the pipeline graph in the same commit.

**It ships its own `.gitignore`**, covering objects, the binaries `configure.sh`
generates into `bin/`, the `Makefile.inc` that script writes, and the `output/`
and `restart/` directories. Nested ignore files are honoured, so unlike cGENIE --
whose own ignore file carries only `.DS_Store` -- this subtree needs no build
rules added to the repository root.

**The licence is AGPL-3.0**, which is stronger copyleft than the other three
subtrees carry: MIT for cGENIE, MPL-2.0 for LPJ-GUESS, GPL for ExoPlaSim. That
is not a problem for reading or for local modification, and it is worth knowing
before any of it is copied into code that leaves this repository.
