# The patches, and what they are now

**These are no longer applied to anything.** The model source is the
`vendor/exoplasim` subtree, and every change here is a commit in it. Nothing
reads this directory at build time; `rebuild_binaries.py` compiles whatever is
in the subtree.

What they are is the AUTHORED RECORD: the header of each file argues the defect,
gives the evidence, and names the base it was written against. Nine of them
carry between thirty and two hundred lines of that, and it exists nowhere else
in this repository -- the subtree is imported with `--squash`, so the fork's
per-change commit messages are not visible from here.

So this directory answers "why is the model like that", and the subtree answers
"what is the model". Scripts and notes across the project cite these files by
path for exactly that reason, which is also why the paths have not moved.

**Do not edit them to change the model.** A change goes to the subtree, and to
the fork, as a commit. If a patch here and the subtree ever disagree, the
subtree is right and the file here is a historical document that should say so.

## What went upstream

Nine of these are defects in stock ExoPlaSim rather than choices for this world,
and are open as pull requests against `alphaparrot/ExoPlaSim`:

| patch | PR |
| --- | --- |
| `pyburn-quadratic-read` | #52 |
| `nlowio-broadcast` | #53 |
| `lowio-first-record` | #54, stacked on #53 |
| `arasc-output` | #55 |
| `rayleigh-reference-grid`, `makestellarspec` | #56 |
| `energy-diagnostics`, `denergy-accumulator` | #57 |
| `aerocore-defects` | #58 |

Plus #59, which is not a patch here: `compile.sh` runs `rm -rf *` after an
unguarded `cd` into a build directory that does not exist in a fresh checkout,
so it deletes the source tree. That was found by vendoring, and it deleted 1437
files before it was found.

The rest are capability rather than repair -- the shortwave weights, the dust
and aerosol stack, the stellar cycle -- and each defaults to reproducing stock
behaviour, so they are offerable upstream but have not been offered.
