# The patches, and what they are now

**These are no longer applied to anything.** The model source is the
`vendor/exoplasim` subtree, and every change here is a commit in it. Nothing
reads this directory at build time; `rebuild_binaries.py` compiles whatever is
in the subtree.

What they are is the AUTHORED RECORD: a patch's header argues the defect,
gives the evidence, and names the base it was written against. Twelve of the
eighteen carry between thirty and two hundred lines of that; the other six
open directly with the diff, and their argument lives only in the notes that
cite them. Either way it exists nowhere else in this repository -- the subtree is imported with `--squash`, so the fork's
per-change commit messages are not visible from here.

So this directory answers "why is the model like that", and the subtree answers
"what is the model". Scripts and notes across the project cite these files by
path for exactly that reason, which is also why the paths have not moved.

**Do not edit them to change the model.** A change goes to the subtree, and to
the fork, as a commit. If a patch here and the subtree ever disagree, the
subtree is right and the file here is a historical document that should say so.

## What went upstream

Most of these are open as pull requests against `alphaparrot/ExoPlaSim`:

| patch | PR | kind |
| --- | --- | --- |
| `pyburn-quadratic-read` | #52 | repair |
| `nlowio-broadcast` | #53 | repair |
| `lowio-first-record` | #54, stacked on #53 | repair |
| `arasc-output` | #55 | repair |
| `rayleigh-reference-grid`, `makestellarspec` | #56 | repair |
| `energy-diagnostics`, `denergy-accumulator` | #57 | repair |
| `aerocore-defects` | #58 | repair |
| `aerosol-apart` | #60 | repair |
| `aerosol-deposition` | #61 | capability, replacing a timestep-dependent sink |
| `h2o-shortwave-weight`, `ozone-band-weights`, `co2-shortwave` | #62 | capability |

Plus #59, which is not a patch here: `compile.sh` runs `rm -rf *` after an
unguarded `cd` into a build directory that does not exist in a fresh checkout,
so it deletes the source tree. That was found by vendoring, and it deleted 1437
files before it was found.

#58 and #61 are coupled and the coupling is stated on both: #58 fixes the
sedimentation sign while the 99%-per-step bottom-layer scrub is still in place,
so merging it alone gives two sinks where the paper describes one.
`notes/audits/aerosol-particle-radius.md` carries the table.

`aerosol-longwave`, `dust-emission`, `prescribed-dust` and `star-cycle` stay
local by decision: the dust pair while this project is its only user, the
other two until someone opens a task to offer them. Each defaults to
reproducing stock behaviour.
