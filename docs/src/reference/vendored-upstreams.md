# The vendored upstreams

Two of this project's components are other people's code, vendored as git
subtrees from personal forks so that a change to the component and the change to
whatever consumes it land in ONE commit, and so provenance is a commit in this
repository rather than the state of a directory outside it.

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

### It is a HARD fork as of 2026-08-21

**Upstream compatibility is no longer a constraint on this component.** The
threaded build is the line this model develops along, and the changes it wants
next -- full-globe grid arrays for SHTns above all -- cannot be made while the
MPI decomposition also has to keep working. Declared deliberately, so that the
next person does not preserve a contract nobody is holding.

What that licenses: removing the MPI path, restructuring the parallelism,
changing array shapes and storage classes across the physics, and taking any
change that is faster without asking whether it could be contributed back.
Pulling upstream is not expected to be possible again, and a change worth
sending upstream now has to be written for upstream separately.

**What it does NOT license, and this is the part that bites.** The MPI build is
the reference every correctness check in this component compares against --
bit-identity at T21 on two, rounding scale at T170 on sixteen, every
`compare_restarts.py` arm. It caught the dv2uv planetary vorticity race, the
weight pre-scaling bug, and the mkdheat inertness result. **Do not delete it
until the serial build has replaced it as that reference**: threaded at N
threads against `mpimod_stub` at one, at rounding scale. The serial build is
arguably the better reference anyway, since it has no reduce-scatter ordering to
explain away.

And removing MPI is not itself a speedup. The threaded build does not link it;
`${MPIMOD}` selects a different file. The fork buys permission to BREAK the MPI
path when it blocks a restructure, which is a reason to drop it lazily at the
point it is in the way rather than as a task of its own.

### Before the fork

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

