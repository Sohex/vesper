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
OCN-4 names beside MARBL. And `genie-plasim` is a PlaSim coupling module, which
bears directly on OCN-17's evaluation of the published ExoPlaSim-to-cGENIE
forcing path.

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

Two consequences to know. The history still carries the deleted blobs, so `.git`
did not shrink and will not; the saving is working-tree size, which is what
`grep` and `find` pay. And a `git subtree pull` will reintroduce the plots until
the same deletion is made in the fork, which is a `git subtree push` to
`Sohex/cgenie.muffin` and has not been done.

Unrelated and worth knowing before someone blames the deletion: four directories
that configs reference, `fkl_np10`, `fkl_pp01_DH`, `fm0450ab` and `wppcont1`,
are absent upstream and were never in the vendored commit.
