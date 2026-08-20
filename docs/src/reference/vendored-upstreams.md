# The vendored upstreams

Two of this project's components are other people's code, vendored as git
subtrees from personal forks so that a change to the component and the change to
whatever consumes it land in ONE commit, and so provenance is a commit in this
repository rather than the state of a directory outside it.

## World Orogen

The geography comes from a personal fork of World Orogen, vendored into this
repo at `vendor/orogen/` as a git subtree from the `cf-fork` branch of
`raguilar011095/planet_heightmap_generation`. It lives here so that a change to
the generator and the change to whatever consumes it land in ONE commit, and so
a build's provenance is a commit in this repository rather than the state of a
directory outside it. Pull upstream with `git subtree pull --prefix vendor/orogen
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

