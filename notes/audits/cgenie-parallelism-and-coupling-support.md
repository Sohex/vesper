# Can the EMIC ocean be parallelised, and what grid should it be coupled on

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean, wind stress, sea ice and tracer transport
name modelled quantities and numerical properties of a candidate component, not
observations of anything.

Measured 2026-08-25 on this machine. This is OCN-19's profile, the parallelism
half of OCN-20's precondition, and the support and regridding half of OCN-11 and
OCN-10. `analysis/cgenie_profile.py` is the driver and
`analysis/cgenie_profile.json` is what it wrote; the toolchain, the make
overrides and the host load are stamped in its provenance block rather than
repeated here.

It follows `notes/audits/cgenie-build-cost-and-grid-ceiling.md`, which prices a
configuration and says nothing about which routine spends it. Everything here
takes that document's costs as given and answers the question it left open.

**The vendored tree was not edited.** It was exported out of the repository with
`git archive` and built there, because a worktree's ignored build products are
symlinks into the shared checkout.

---

# 0. The decision rule, fixed before the profile was read

Written into `analysis/cgenie_profile.py` and committed before any profile
existed, because a threshold chosen after the run it judges is not a threshold.

A routine counts as **parallel over cells or tracers** when its loop nest runs
over the grid indices and any loop-carried dependence in it is a reuse of a face
flux the previous cell already computed -- recoverable by recomputing one row,
one column or one level per thread -- with no global reduction other than a sum
or a maximum. It counts as **serial** when it carries a dependence over the whole
domain that no such recomputation removes; the banded triangular solve in
`ubarsolv.f` is the case that motivates the definition. Anything that is neither
is **unclassified** and is counted against the parallel fraction, not for it.

The verdict on question one then follows from the SERIAL share at 36 x 36 x 16:

| serial share | verdict |
| --- | --- |
| below 5 per cent | parallelisation is possible and Amdahl is not what limits it at this grid |
| 5 to 20 per cent | possible, with a stated ceiling that must be quoted with every speedup claim |
| above 20 per cent | the algorithmic route of OCN-20 comes before the threading route of OCN-19 |

The bound at one grid is not the bound at another, because the serial part grows
faster with resolution than the parallel part does. Both grids are therefore
reported and the trend between them is part of the answer.

---

# 1. What is serial today, re-verified rather than inherited

`notes/audits/ocean-and-marine-biosphere.md` section 9h and the cost note's
section 1a both say the component is serial. Re-checked here over the whole of
`vendor/cgenie`, because the claim is load-bearing for everything below:

- **No live OpenMP sentinel exists anywhere in the tree.** Every `$omp` string in
  it sits in exactly three files -- `genie-goldstein/src/fortran/tstepo.F`,
  `genie-biogem/src/fortran/biogem.f90` and `genie-main/genie.F` -- and every one
  is commented out in column one, so none is a directive.
- **Four `use omp_lib` statements are live** and no call stands against any of
  them: `biogem.f90`, `tstepo.F`, `genie.F` and `genie_loop_wrappers.f90`. The
  earlier count of three missed the last.
- **`-fopenmp` is commented out** in `genie-main/makefile.arc`'s gfortran block,
  and enabling it would parallelise nothing.
- **MPI exists in one directory**, `genie-plasim/src/fortran`, which is the
  vendored Planet Simulator this project does not use.

So the question is not what to re-enable. It is what a decomposition would have
to be, and what it would buy.

---

# 2. Where the instructions actually go

Placeholder: filled from `analysis/cgenie_profile.json`.

---

# 3. What the parallelisation would cost to do

Placeholder.

---

# 4. The coupling resolution

Placeholder.

---

# 5. The regridding contract

Placeholder.

---

# 6. What this did NOT establish

Placeholder.
