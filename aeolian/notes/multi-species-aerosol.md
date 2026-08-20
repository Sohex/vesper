# One aerosol slot, three species: what a second one takes

Worldbuilding. Vesper is an invented planet and this note is about the
simulation of it: a toy climate model's radiation scheme and the offline aerosol
components that feed it. Every quantity named here is a modelled field.

Written 2026-08-20 from source inspection of the `vendor/exoplasim` fork, in the
shape `in-model-dust.md` established for DUST-3. It sizes the work rather than
doing it, and the headline is that the cheap version is not cheaper.

## 1. The finding: the radiation carries exactly one aerosol

`aeolian/` produces three aerosol species and `analysis/error_budget.json` prices
all three. None of them is in the radiation, and the interface will not take more
than one.

`radmod.f90:173` declares `aeroqs(8,1)`: four optical constants per band, one
species. From it `radini` computes `ssa1`, `ssa2`, `bscat1`, `bscat2`, `qex1`,
`qex2` at lines 1106-1111, and every one is a SCALAR. In `swr` at 2328-2340 the
per-layer optical depth comes from either `ddustod` or `daerod` and never both,
and the longwave at 2875-2877 uses either `dustqlw` or `aeroqlw`. `radini`
refuses the two paths at once at 1162, with the fork's own reason recorded at
line 226: a prescribed column and a transported one are two aerosols, and adding
their optical depths would count one of them twice.

So only the optical depth varies in space. The optical PROPERTIES are one global
set of numbers for the whole planet.

## 2. Why the properties cannot be shared

Dust absorbs. Sea salt's single-scattering albedo is 1 to within 1e-5 in both
bands, measured in `aeolian/analysis/sea_salt_baseline.json`. Those are the two
extremes of the one quantity the code stores as a scalar, so no single value
describes both.

Geography does not rescue it. Dust is emitted over land and sea salt over ocean,
but dust transports over ocean -- that transport is the whole of the aeolian
deposition leg -- so the two are co-located in the same layers over a large part
of the planet, and a per-cell choice between two scalars would be wrong wherever
they overlap.

What is correct is an external mixture, per cell and per layer:

    tau_eff   = sum_s tau_s
    ssa_eff   = sum_s (ssa_s tau_s) / tau_eff
    bscat_eff = sum_s (bscat_s ssa_s tau_s) / sum_s (ssa_s tau_s)

and the same in band 2. That promotes `ssa1`, `ssa2`, `bscat1`, `bscat2` from
scalars to `(NHOR,NLEV)` fields, which moves the two-stream u-factors

    zaeru1 = SQRT((1.0-ssa1+2*bscat1*ssa1)/(1.0-ssa1))
    ztemp1 = SQRT((1.0-ssa1)*(1.0-ssa1+2*bscat1*ssa1))

from a once-per-call computation outside the layer loop to inside it.

**That restructure is the whole cost, and it does not depend on how many species
follow it.** It is required to put a second species in, and having done it, a
third costs a loop bound.

## 3. Sizing, against two patches that did comparable work

The fork carries 18 patches over 4158 lines. The comparables are
`exoplasim-3.4.2-prescribed-dust.patch` at 530 lines, which built the entire
single-species prescribed path including the `surfmod` code, the `radini`
checks, `dustprof` and its column-conservation self-check; and
`exoplasim-3.4.2-dust-emission.patch` at 827.

| work item | second column | arbitrary N |
| --- | --- | --- |
| per-layer mixing, u-factors into the loop | required | identical |
| longwave `exp(-1.66 sum_s qlw_s od_s)` | required | free: already additive in the exponent |
| `dustprof` per species with its own scale height | required | loop bound |
| `aeroqs` and the aerofile | widen to `(8,2)` | widen to `(8,N)` |
| namelist `dustsc`, `dusthsc`, `dustqlw` | a second copy of each | arrays |
| `surfmod` surfcode plus `.sra` codes | one line, 1812 | a block, 1811 onward |
| Python: aerofile writer, surface builder, `run_exoplasim` code set, pipeline rows | same | same |

A patch in the same 400 to 700 line band as `prescribed-dust`, plus small
`surfmod` and Python changes. The two-species form saves perhaps a fifth of that
and returns it the first time a third species lands, because a hand-duplicated
pair is what makes the third expensive. Two is the hardest number.

**The aerofile side is nearly free, because the species axis already exists.**
`readdat(filename,ndim,nitems,kdata)` in `utilities.f90:266` fills
`kdata(nitems,ndim)` and is called `readdat(aerofile,1,8,aeroqs)` against
`aeroqs(8,1)`. `ndim` is the species count and it is set to 1. Going to N is one
argument and one declaration, with the file gaining a column per species rather
than changing shape.

## 4. The trap: the conservative-scattering limit is untested code

`zaeru1` divides by `(1 - ssa1)`. Dust's single-scattering albedo is nowhere near
1, so nothing has ever driven that expression toward its limit. Sea salt sits at
`1 - 1e-5`, which would run it near-singular over the 57% of this planet that is
ocean, on every radiation step.

At the declared 8-byte precision this is a precision loss rather than an
overflow, which is the reason to treat it seriously: it fails by returning a
plausible number. A conservative-scattering branch is part of the work, and it
needs its own check rather than being assumed correct because the model still
runs.

## 5. The interactive path is not a separate question

The prescribed and transported paths are mutually exclusive today only because
they write into one slot. With a species array the transported tracer is one
species of it, and the exclusion at `radini:1162` goes away rather than being
preserved. That is the decision this work takes, and it is what keeps DUST-13
open rather than foreclosed: as the code stands, choosing interactive dust would
put sea salt permanently out of reach, because `ndustrad` and `iaerint` cannot
both be on.

## 6. What is worth routing in, and what is not

| species | forcing, W/m2 | route |
| --- | --- | --- |
| mineral dust | +0.34 to +0.61 | the field and the aerofile both exist |
| sea salt | -0.16 to -0.89 | the offline field exists; needs an `.sra` writer and an aerofile column |
| volcanic sulfate | -0.013 to -0.032 | one to two orders down; nearly free once the array exists, and it decides nothing |
| carbonaceous | not yet priced | CLIM-29, blocked on a biosphere run |

Sea salt is the one that changes an answer. It is also the cheaper half of the
pair to prescribe, because it is a pure scatterer over a dark surface that does
not change, so none of the interactive machinery DUST-13 argues about applies to
it.

`error_budget.json` already states why the pair cannot be reduced to its sum:
centrally the two are about -0.07 W/m2 together, which is a GLOBAL-MEAN
cancellation between an aerosol acting over land and one acting over ocean. It
does not hold by surface, and it does not hold in runoff, where dust is several
times sea salt because sea salt reaches land only through temperature. The
residual is a change in land-sea contrast, which is the circulation response
DUST-11 exists to measure, and runoff is the carve criterion's denominator.

## 7. The tests, declared before the work

**The reduction identity.** With every species after the first set to zero
optical depth, the run must be BIT-IDENTICAL to a single-species run of the same
configuration. Not close: identical. The mixing algebra collapses exactly when
one term is nonzero, and anything else means the restructure changed the
single-species answer.

**Conservative scattering.** An aerosol at ssa = 1 over a black surface absorbs
nothing. The column absorption must come out zero to round-off, which is the
check the near-singular branch needs and does not currently have.

**Column conservation per species.** `dustprof` already reports
`max |column - prescribed|` once per run and is written to fail if the
normalisation stops holding. Each species keeps its own, rather than one check
over the sum, or a species with the wrong scale height hides inside a correct
total.

**Mixing against the offline chain.** Dust and sea salt at their offline optical
depths, mixed in one layer, must reproduce the extinction-weighted single
scattering albedo the two components report separately. That is an identity
between the model and `analysis/` rather than an agreement between two estimates.

Rule 4 applies whichever form is built: every binary rebuilt, then `--verify`.

## 8. What the tests returned

Built and measured 2026-08-20 at T42 L10 on 8 ranks, 8-byte precision.
`exoplasim/patches/exoplasim-3.4.2-multi-species-aerosol.patch` carries the
full record; what belongs here is the one result that changes how any future
test in this project is written.

**The reduction identity was recorded as passing at one timestep**, the same
configuration through a binary built from the base and through the changed one
giving bit-identical gridpoint output.

**That claim was made at one timestep because a same-binary control appeared to
show the model diverging by sixteen. IT DOES NOT.** CLIM-44 re-took the control
as a declared factorial -- 8 and 16 ranks, output off and on, segments crossing
no write and two, three repeats each -- and every measurable cell is
bit-identical in `plasim_status`, `plasim_output`, `plasim_snapshot` and the
model half of `plasim_diag`. The same result holds on the pre-fix binaries this
worktree was built against, so it is not something the intervening fixes bought.

Two things were wrong with the control here. The rank count was not held fixed
across everything it was compared with, and 8 against 16 ranks really does
change the answer, at round-off and growing. And a gridpoint record is written
on `mod(nstep, nafter) == 0` over the ABSOLUTE step count, so on this project's
production restart neither a 1-step nor a 16-step segment writes one at all: the
comparison that reported agreement at one timestep could not have failed.

**The reduction identity is therefore UNVERIFIED rather than refuted.** It was
stated at the one horizon that turns out to carry no output, so it needs
re-running at a segment length chosen against the write cadence, at a fixed rank
count. `notes/audits/model-reproducibility.md` is the measurement and CLIM-45
is the re-test.

**Column conservation holds per species**, against a column maximum of 0.4896:
1.30e-18 for one species, and 8.67e-19 and exactly 0.0 for two, with the second
species deliberately given a different scale height so a shared one would have
shown.

**The conservative branch runs clean.** A species with `Qsca = Qext` exactly,
so `ssa` is 1 to the bit, produced no NaN anywhere and `T + R - 1 = 0` exactly.

**The single-precision case is worse than section 4 predicted**, and this is the
correction to it. At 4-byte precision `1.0 - ssa` rounds to exactly zero at
`1-ssa = 1e-8` and the exact expression returns **NaN**, not the plausible wrong
number the section anticipated. At 8-byte it is benign to `1-ssa = 1e-14`. The
threshold `SQRT(EPSILON(1.0))` puts sea salt on the exact branch in a double
build and the limit in a single one.

**The fast path is load-bearing, not cosmetic.** `(ssa*tau)/tau` returns
something other than `ssa` in 12,026 of 100,000 sampled cases at 8-byte
precision, so without it every single-species cell would move by an ulp and the
existing dust-only answer would not reproduce.

## 9. A prescribed species' aerofile column carries only RATIOS

Read out of the implemented path, and it changes what CLIM-40 has to derive.

For a PRESCRIBED species the absolute extinction efficiency never enters the
radiation. The optical depth comes from the column field through
`dustsc*ddustcol*zw/zsum`, which has no `Qext` in it, and every place the
optics are used afterwards is a ratio: `ssa` is `Qsca/Qext`, the backscatter
ratio is `Qback/Qsca`, the band split is `qex2/qex1`, and the longwave is
`dustqlw` per unit band-1 optical depth. Scale a prescribed species' whole
column of the aerofile by any constant and the answer does not move.

That is why `dust_aerofile.py`'s careful work is specific to the INTERACTIVE
path rather than general. There `aeroprof` builds the optical depth as
`nrho*PI*apart**2*qex1`, so the absolute `Qext` and the particle radius it is
declared against are both load-bearing, and the burden-matched rescale that
file documents exists for exactly that reason.

So a second PRESCRIBED species does not need a burden-matched radius, an
`apart`, or an absolute cross-section. It needs four ratios per band. That is a
much smaller derivation than the dust one, and it can be taken straight from
`analysis/sea_salt_optics.json` without touching the settling machinery.
