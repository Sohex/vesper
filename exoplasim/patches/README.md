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

Plus #63, also not a patch here: two restart-write defects in `radmod.f90`.
`radstop` saved `zsolars`, a `real(2)`, through `mpputgp`, which gathers `NHOR`
per rank into a local `z(NUGP,klev)` and writes all of it, so each rank read 510
elements past the end of the array and the record was 8190 elements of buffer;
and `solarini`'s nine `put_restart_array` calls run before any restart file is
open on `nwriunit`, so they land in a stray `fort.34` and none of the nine
arrays reaches the restart. Both are upstream verbatim, checked against
`alphaparrot/ExoPlaSim` master before the PR was opened. The argument, the
evidence and the verification are in
`notes/audits/zsolars-restart-overread.md`; CLIM-37 and CLIM-38 closed them
here. Offered as ONE pull request with an offer to split, because the second is
the other half of the same incomplete edit and the commented-out read block is
what makes both self-evident.

Plus #64, also not a patch here: `legmod.f90` applies the wavenumber-dependent
physics filter inside the innermost loop of every spectral transform, where it
is invariant in the enclosing latitude loop. Every use of `qi`, `qj`, `qu` and
`qv` carries `skspgp` and every use of `qc`, `qe`, `qm` and `qq` carries
`skgpsp`, 94 call sites without exception, so the filter folds into the weight
matrices in `legini` and comes out of the loops. Offered as a SIMPLIFICATION
that is slightly faster rather than as an optimisation: it is worth 1.6% at
T127, 0.9% at T42 and nothing at T21, because those loops are bound by
streaming the weight matrices rather than by the multiplier. Bit-identical to
stock with the filter off and contraction disabled, which is what proves the
mode indexing; a reassociation under default flags, because deleting a multiply
changes what the compiler contracts to an FMA. The argument, the measurements
and the test are in `exoplasim/notes/legendre-filter-fold.md`.

Also not a patch here, and not yet offered: `plasim/src/make_plasim` declared
three of its own dependencies short. `glaciermod.o` did not depend on
`landmod.o`, `plasim.o` did not depend on `radmod.o`, and `carbonmod.o` did not
depend on `${RAINMOD}.o`, while each of those files USEs the module the other
defines. The first two link EARLIER in `OBJ` than the module they need, so a
build from an empty tree fails; upstream never sees it because `compile.sh`
empties `plasim/bld` only when switching between the MPI and serial builds, so
a `.mod` from some previous build is always lying there. The same three gaps
made `make -j` unsafe at any level. The argument and the proof that the fixed
makefile builds clean at `-j16` and produces a byte-identical executable are in
`notes/audits/aocl-and-model-build-flags.md`.

#58 and #61 are coupled and the coupling is stated on both: #58 fixes the
sedimentation sign while the 99%-per-step bottom-layer scrub is still in place,
so merging it alone gives two sinks where the paper describes one.
`notes/audits/aerosol-particle-radius.md` carries the table.

`aerosol-longwave`, `dust-emission`, `prescribed-dust`, `multi-species-aerosol`
and `star-cycle` stay local by decision: the dust set while this project is its
only user, the other two until someone opens a task to offer them. Each defaults
to reproducing stock behaviour.

`multi-species-aerosol` is the one most likely to be wanted upstream, because
the limit it removes is not specific to this world: stock ExoPlaSim carries one
aerosol with one global set of optical properties, so no two species with
different single-scattering albedos can be in the radiation at once. It is
offered when someone opens the task.
