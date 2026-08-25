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

The loop structure is not the obstacle. A build flag is.

## 3a. `-fno-automatic` makes every procedure-body local SHARED, silently

`makefile.arc`'s gfortran SHIP block compiles with
`-O2 -O3 -funroll-loops -msse -fno-automatic`. That last flag moves every
procedure-body local into static storage, and a variable in static storage is one
variable for the whole process. Inside an OpenMP region it is therefore shared
between the threads unless it is named in a `private` clause, and the compiler
will not tell you which ones you missed.

Demonstrated rather than argued, with a loop whose right answer is known:

    subroutine work(n, out)
      real :: scratch(512)
      !$omp parallel do private(j)          ! scratch NOT named
      do i = 1, n
         do j = 1, 512
            scratch(j) = real(i) + real(j)
         end do
         out(i) = scratch(1) + scratch(512) ! must be 2*i + 513
      end do

Built `gfortran -O2 -fno-automatic -fopenmp`, run over 200,000 iterations:

| threads | private clause | cells wrong |
| --- | --- | ---: |
| 1 | none | 0 of 200,000 |
| 8 | none | 177,007 of 200,000 |
| 8 | `private(scratch)` | 0 of 200,000 |

Three things follow, and each is a property of the work rather than of this
example. The wrong answer is a real answer -- finite, plausible, and dependent on
the thread count, which is the failure class `docs/src/practice/failure-modes.md`
already tracks. `private` on such a variable is ACCEPTED by gfortran rather than
rejected, so there is no compile-time gate to lean on. And the whole combination
draws exactly one diagnostic for the entire translation unit,
`Flag '-fno-automatic' overwrites '-frecursive' implied by '-fopenmp'`, which
says nothing about which variables are affected.

**So the cost of the threading route is an enumeration, not a directive.** Every
procedure-body local in every routine reachable from a parallel region has to be
classified and named. That is the same class as CLIM-51's 103 implicitly SAVEd
declarations in `plasim/src`, arriving here through a compiler flag instead of
through the language default.

## 3b. The flag looks removable, and the test for it already ships

The cheaper route is to stop compiling that way. `-frecursive` gives an unlimited
`-fmax-stack-var-size`, which puts the same locals on the stack where each thread
has its own. It costs stack rather than correctness, and `makefile.arc` already
carries the commented line `F90FLAGS += -frecursive` beside gfortran's own
warning about exactly this.

Whether the tree DEPENDS on the implicit SAVE that `-fno-automatic` supplies is
checkable, and the evidence says it very likely does not. Across
`genie-goldstein`, `genie-embm` and `genie-goldsteinseaice` there are ten `DATA`
statements and eight `save` statements, every one of them EXPLICIT and every one
in a driver, a diagnostic or a netCDF writer rather than in a per-cell kernel:
`goldstein.F:133,137` and `gold_seaice.F:65,69` (the first-call flags and the
initial energy and water inventories), `diag2.f:22`, `diagosc.F:56-58`,
`diagend.F:48` and `diagend_embm.F:47` (diagnostic accumulators), and the month
name tables in the three `netcdf*.F` writers. A routine that declares its
persistence keeps it under either flag.

**And the acceptance test needs no invention.** `genie-knowngood/` ships
per-component netCDF output for four configurations and `genie-main/src/c/compare.cpp`
builds `nccompare`, which takes a relative tolerance in units in the last place.
Dropping `-fno-automatic` is a change that should move nothing at all, so its
right answer is bit-for-bit or within a declared ULP count -- a test that can
fail, in the sense `CLAUDE.md` requires, rather than a comparison that can only
differ. That check is the first step of any parallelisation work here and it
needs no threads.

## 3c. What a decomposition would look like, read off the loop structure

Stated because the profile below only says where the time is, not whether it can
be divided. Each of these is a reading of the source, not a measurement.

- **`tstepo_flux.F`, the ocean tracer transport.** A `k, j, i, l` nest with a
  flux-recycling optimisation: `fw(l) = fe(l)`, `fs(l,i) = fn(l)` and
  `fb(l,i,j) = fa(l)` carry the face flux forward so it is computed once per
  face. That is a loop-carried dependence in all three space indices, and it is
  the recoverable kind -- each of those is a pure function of `ts1` and `u` at
  the face, so a thread taking a block of rows recomputes one row of meridional
  flux at its southern edge and owns the rest. The code already does exactly this
  at `i = 1` for the western doorway. The two genuine reductions, `limps` and
  `dmax`, are a count and a maximum.
- **`co.F` convection, `krausturner.F`, `jbar.f`, `velc.f`.** Column algorithms:
  a vertical recurrence inside each `(i,j)`, no horizontal dependence. Parallel
  over columns, which is thousands of them at any grid considered here.
- **`eos.f` is pure.** It declares no locals and is called from the innermost
  tracer loop, so it needs nothing.
- **EMBM's `tstipa.f` is a Jacobi iteration**, not a line solve: `tq2` holds the
  whole previous iterate and the update reads only `tq2` and `tq1`. Every one of
  its `nii` sweeps is data-parallel over cells with a barrier between them.
  `tstepa.f` is the same flux-recycling shape as `tstepo_flux.F` in two
  dimensions.
- **`surflux.F` is one `(i,j)` loop** running from line 683 to about 1600,
  including the sea-ice surface Newton iteration and the land column. Parallel
  over cells, and the routine with the most locals to classify by a wide margin.
- **`ubarsolv.f` is the exception.** Its forward elimination and back
  substitution over the banded factors are sequential along the band by
  construction, and no recomputation removes it. This is the routine the decision
  rule in section 0 is about.

## 3d. The working set is cache-resident at every grid considered

`CLAUDE.md`'s thread-team target is 32 MB on one die. The arrays a GOLDSTEIN
timestep actually touches are `ts`, `ts1`, `u` and `rho`: at 36 x 36 x 16 in
double precision that is roughly 2 MB, at 72 x 72 x 16 roughly 8 MB, and at
144 x 144 x 16 roughly 32 MB. So a thread team on this part would be dividing an
arithmetic-bound problem rather than a bandwidth-bound one, up to about the grid
where the horizontal doubling stops being affordable anyway. The gigabyte of
static COMMON that forces `-mcmodel=medium` at 72 x 72 is not this: it is
`embm.cmn`'s seasonal arrays and `ocean.cmn`'s `maxnyr`-dimensioned storage,
which the timestep does not read.

---

# 4. The coupling resolution

Placeholder.

---

# 5. The regridding contract

Two directions, and they are not one problem. Upscaling the atmosphere's fluxes
onto the ocean is a conservative reduction with an exact answer. Downscaling the
ocean's surface state onto the atmosphere invents structure the ocean never
resolved. They cost different things and the contract states them separately.

## 5a. The precedent chose matching, and the reason does not transfer intact

Holden et al. (2016) is the only existing PlaSim-to-GOLDSTEIN coupling and it
restricted itself to T21 against a matched 64 x 32 GOLDSTEIN grid, in its own
words "in order to avoid the need for interpolation", with 10 atmospheric levels
against 32 ocean levels. That coupling was SYNCHRONOUS: fluxes crossed every
coupling step, so an interpolation error entered the prognostic state, was
integrated, and could feed back through the atmosphere's own response to the sea
surface temperature it had just perturbed.

This project's coupling is offline and passes climatologies, so the feedback path
is gone. What is NOT gone, and is worse in one specific way, is the residual: a
climatology is a fixed field, so a non-conservative regrid of it is a permanent,
spatially structured spurious source term in an ocean that is then integrated for
of order 10,000 model years. A synchronous coupling's interpolation error is
transient and partly self-cancelling; an offline one's is a bias with unlimited
time to express itself.

**So the precedent's requirement changes shape rather than relaxing.** What the
offline architecture needs is not "avoid interpolation" but "the crossing
conserves exactly and the residual is reported", and that is a much more
tractable requirement than Holden's, because it is a property of a weight matrix
computed once rather than of an error budget accumulated every coupling step.
Grid matching is one way to get it and not the only one, and section 59c of
`notes/external-model-survey.md` already records that declining to match means
writing the regridder the offline architecture needs anyway.

## 5b. Upscaling the fluxes: exact, and cheaper than it looks

Both grids are separable in longitude and latitude. ExoPlaSim's is a Gaussian
grid: uniform in longitude, Gaussian latitudes. GOLDSTEIN's shipped `igrid = 0`
is uniform in longitude and uniform in the SINE of latitude, and
`initialise_goldstein.F:496` sets `asurf(j) = rsc*rsc*ds(j)*dphi` -- so every
ocean cell has the same area exactly, by construction.

Two consequences follow that make this direction cheap.

**The conservative weight matrix factorises.** Because neither grid's cell
boundaries depend on the other coordinate, the overlap area of an atmosphere cell
and an ocean cell is the product of a longitude overlap and a sine-of-latitude
overlap. The operator is therefore two one-dimensional sparse matrices rather
than a polygon intersection, it is exact in floating point, and the planetary
radius cancels out of it because every weight is a ratio of areas.

**The nonlinearity has already been resolved on the fine grid.** A surface heat
flux is not a linear function of the surface state, so the average of the flux is
not the flux of the average. In the architecture OCN-10 selects -- ExoPlaSim
computes the fluxes and there is no EMBM -- what crosses the boundary IS the
flux, already evaluated at 1.4 degrees or wherever the accepted rung sits, and
averaging it is a conservative linear operation. In the eliminated offline-regrid
architecture EMBM would have recomputed the fluxes from coarse state, putting the
nonlinearity back at the ocean's resolution. **That is an argument for the chosen
architecture that the record did not carry**, and it sits beside the tuning
argument in section 59b rather than replacing it.

What upscaling costs in fidelity is therefore exactly one thing: the subgrid
variance of the flux inside an ocean cell is discarded. `lib/gridding.py` already
owns the vocabulary for saying what that means field by field -- extensive,
intensive, categorical, moments, expectation -- and its `cell_moments` and
`cell_expectation` are the operators for the cases where a mean is not enough.
Which of them a field takes is a property of the field and belongs in the
contract, not in the call.

## 5c. Four rules the upscaling contract has to state

Read off `references/esmf/` (survey section 55), which is the only conservative
regridder in this project's reading with a stated conservation identity and a
test suite against it.

1. **Normalisation is a property of the field.** A destination value formed as
   intersection area over DESTINATION area is what makes an extensive total
   close; a destination value formed as intersection area over COVERED area is
   what makes a density or a rate come out right. A field of ones must remap to
   ones under the second and does not under the first. Net heat flux, freshwater
   flux and wind stress are rates per unit area and take the second; an integrated
   total takes the first.
2. **The coverage fraction travels with the result.** A fraction that is
   discarded is a conservation identity that can no longer be evaluated.
3. **An unmapped destination cell is an error, not a backfill.** `lib/gridding.py`
   currently backfills a grid cell with no mesh region from the nearest region
   centre and returns the count; the precedent refuses by default and makes the
   fallback the opt-in. For the ocean crossing the refusal is the right default,
   because an ocean cell with no atmosphere over it is a mask disagreement rather
   than a sparse-coverage artifact.
4. **The closure residual is registered in advance and stamped in the artifact.**
   `sum(F_src * A_src * frac_src)` against `sum(F_dst * A_dst)` has a right answer
   and the precedent holds itself to about 1e-9 relative on analytic fields.
   `lib/gridding.py:transfer_ledger` already emits the shape of this; what it has
   never been checked against is a field whose integral is known.

## 5d. The mask disagreement is the part that does not have a clean answer

The ocean's wet mask comes from the Orogen mesh through OCN-11, and the
atmosphere's comes from the same mesh at a different resolution. They will not
agree at the coast, and they cannot: a cell that is 40 per cent land at 1.4
degrees is inside a 5.6-degree ocean cell that is either wet or dry.

So there are atmosphere cells carrying an ocean flux that lie over ocean-grid
land, and ocean cells with no atmosphere-ocean area over them at all. Every joule
and every kilogram in the first class has to go somewhere named, or the global
heat and freshwater budgets do not close and the residual is not noise but a
coastline-shaped field. This is the same problem `references/climber-x/`'s
`coast_cells.f90` solves for river discharge, cited in OCN-11's notes, and the
contract's answer has to cover heat and salt as well as water.

**This is the part of the crossing that a finer ocean grid genuinely improves**,
and it is the argument for resolution that survives when the conservation
argument is settled by construction.

## 5e. Downscaling the ocean's surface state: what it invents

Three fields cross the other way, and none of them is a flux: the sea surface
temperature the next climate run is forced with, the ice state world-pt8 wants,
and the surface velocity `goldstein.F:780-781` exports as `ustar_ocn` and
`vstar_ocn`.

These are intensive state variables, so the crossing is an interpolation and
there is nothing to conserve. What it costs in fidelity is that the atmosphere
receives a field with no structure below the ocean's cell size, while its own
dynamics respond to gradients at its own. At an ocean with 5.6-degree longitudes
against a T85 atmosphere at 1.4, every sea surface temperature gradient the
atmosphere sees is a 5.6-degree gradient smeared over four cells: no fronts, no
western boundary currents, and no straits at all, because `notes/audits/ocean-and-marine-biosphere.md`
section 8b already records that straits and sills are not represented at these
resolutions in the first place.

Two rules follow, and both are about honesty rather than accuracy.

**Interpolate the state, never a flux, in this direction.** A bilinear or
higher-order interpolation of an intensive field is the right operator and a
conservative one is not, because there is no integral to preserve. The
temptation to route a heat transport this way -- an extensive quantity dressed as
a state -- is what the two words in `lib/gridding.py` exist to separate.

**The support mismatch is declared, not smoothed.** An adequacy claim about the
coupled system is a claim about a field carried across a factor of four in
longitude at T85 and a factor of eight at T170, and the honest form of it names
that factor. Section 5c of the cost note reaches the same place from the price
side: a finer ocean is affordable and a matched one is not, so the mismatch is a
fact to state rather than one to avoid.

---

# 6. What this did NOT establish

Placeholder.
