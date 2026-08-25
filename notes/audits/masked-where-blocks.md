# A `where` block protects the assignment and nothing else

Established 2026-08-25 by static analysis of
`vendor/exoplasim/exoplasim/plasim/src`, against the tree at the time. No model
run was involved and none could settle it: the question is which expressions the
compiler is licensed to evaluate, and that is a property of the source.

Worldbuilding frame: this is about how the Vesper climate model's Fortran is
compiled and which programming defects the build can and cannot detect. Nothing
here concerns the simulated planet.

## The mechanism

A Fortran `where` construct selects which elements of an array assignment are
STORED. It does not restrict which elements of the right-hand side are
EVALUATED. A scalar compilation walks the mask and skips the excluded elements
as a side effect of how it loops; a vectorising one computes the whole vector
and blends, because that is cheaper than a per-lane branch. Both are conforming.

The declared flag line carries `-ffpe-trap=invalid,zero,overflow`
(`config/planet.yaml`). So an excluded lane whose argument is out of domain does
not produce a discarded value: it raises SIGFPE at the instruction that
evaluated it. That is world-bhs, where a T42 run died eight orbits in and the
fault moved to a new site each time one was clamped -- the longwave, then kuo's
Tetens exponent, then the shortwave.

`-O2` rather than `-O3` removed the loop transforms that were firing at the
time. It is not a fix for the class: `-O2` vectorises too, so a different state
or a different rung can wake the same sites.
`notes/audits/model-build-flags.md` carries that measurement.

## Method, and the direction it errs in

`exoplasim/scripts/lint_masked_domains.py` is a lexical pass over the 38
translation units. It folds continuation lines, walks `where` / `elsewhere` /
`end where` nesting, and reports every call to `sqrt`, `log`, `log10`, `exp`,
`asin`, `acos` and their `alog` spellings that lies inside an open construct,
together with the mask that encloses it. Division and `merge` are reported as
separate classes, because both carry the identical hazard for the identical
reason.

**It over-reports, and that direction is what makes a clean run mean
anything.** Three places it errs towards reporting a site that is safe, and none
towards clearing one that is not:

1. The where-stack is popped only by an explicit `end where` or by a procedure
   boundary. A construct the parser fails to close keeps reporting sites after
   it, never fewer.
2. A guard is recognised only from a literal set of shapes -- an argument
   already wrapped in `max`, `min`, `amax1`, `amin1`, `abs` or `dim`. A correct
   guard written some other way is reported.
3. A division is reported wherever its divisor is parenthesised, whether or not
   the divisor can reach zero, which the pass cannot know.

Every reported site was then read in context. The population is 56
domain-sensitive intrinsics inside masked blocks: 20 already guarded by earlier
work, and 36 that this pass classified.

## The classification, and the three kinds it falls into

**The mask is not what makes the argument safe.** `swr`'s
`SQRT(273./dt(:,jlev))` sits under `where(losun(:))`, and the gridpoint
temperature has no bearing on whether the sun is up. A masked lane and a kept
lane are equally safe or equally unsafe, so the site is not of this class at
all. `seastep`'s `SQRT(dtaux**2+dtauy**2)` is the same shape for a different
reason: a sum of squares has no domain to leave. `lwr`'s `ALOG(ztau0)` is a
third: `ztau0` is clamped to `[zero, 1-zero]` two lines above the `where`,
unconditionally, so every lane is in domain before the mask is consulted.

**The mask WAS the guard, and the guard is now structural.** These are the
defects. The pattern is always the same -- a test written as a mask because the
author knew the expression was singular outside it, and a mask does not do that
job:

| Site | What the complement holds | Floor added |
| --- | --- | --- |
| `swr` aerosol u-factor, both bands | `lcons` excludes `1-ssa <= zepsc`; the complement includes exactly 0 and negative | `MAX(1.0, .../MAX(1.0-ssa, zepsc))` |
| `swr` aerosol `ztemp`, both bands | same complement, product of two factors one of which can be negative | `SQRT(MAX(0.0, ...))` |
| `swr` mixture single-scattering albedo | `knz == 0` lanes have `zext1` exactly 0, and the three quotients beside it were already floored | `MAX(zext1, TINY(1.0))` |
| `mkclouds` convective cloud cover | `dprc > 0.` excludes zero rain; `log(0)` raises the divide-by-zero | `log(max(1.E-30, ...))` |
| `mkclouds` overlap exponent | `1./real(icctot)` with `icctot` zero on the same lanes | `real(max(1, icctot))` |
| `mkclouds` cloud liquid water | `zzh > 0.` excludes a column with no water vapour, where `zzh` is exactly 0 | `max(1.E-30, zzh)` and the file's own `min(80.,max(-80.,·))` on the exponent |

The u-factor floor is at 1 rather than at 0 deliberately. The denominator below
it is `(u+1)^2 e^t - (u-1)^2 e^-t`, which for `u >= 1` and `t >= 0` is at least
`4u`; a floor at 0 would let that denominator reach zero and move the trap one
line down instead of removing it.

Each floor is a no-op on every lane the mask keeps, and that is the standard
each had to meet before it went in. Two of them meet it the way rainmod's Tetens
clamps already did -- the enclosing computation flattens the region the floor
can reach, `AMAX1(zccmin,·)` for the cloud cover and `MAX(dql,1.E-9)` for the
liquid water -- and the rest meet it because the mask's own predicate is the
floor's threshold.

**The mask was the guard and the value read was never written.** A `where` that
assigns a scratch array and reads it on the next line of the same block reads,
on an excluded lane, whatever the stack held. `-finit-real=zero` left the
production profile (`notes/audits/uninitialised-reads-and-implicit-save.md`), so
that is not a zero; it is an arbitrary word, and `EXP` of a non-finite one
raises the invalid the profile traps on.

`swr` had two chains of this shape. The cloud-optics chain runs
`zlwp -> ztau -> zlog -> zb2 -> zom0 -> zun -> zuz -> zu -> zexp -> zr`, each
written and then read inside one `where(losun .and. dcc > 0.)`; the aerosol
two-stream runs `zaeru -> ztemp -> zaertf -> zaerd`, each written and read
inside one `where(aod1 > 0. .and. .not. lcons)`. Both are now preset before
their loop, to the values the no-aerosol and no-cloud cases would produce. That
is the fix the same subroutine's magnification factor already carries, and the
comment above it already argued for, applied to the two chains it had not
reached.

Presetting cannot move a result the model uses, and the reason is worth stating
because it is what makes the change free: on every lane the mask keeps, each of
these is stored before it is read, in the same block, on the same pass. The
preset is visible only on lanes whose stores are discarded.

## What is not settled

The division class is reported and not resolved. 133 sites have a parenthesised
divisor inside a masked block, and most of them are sums that cannot vanish
(`1.+ztcon`, `zcap+zdiff`); the pass has no way to tell those from the handful
that can. The masked Tetens divisors within it are settled separately, by
`ra4d` in `plasimmod.f90`, which is world-6tp. The rest is world-d016,
which has to re-derive the population first: the pass reports a division only
where the divisor is parenthesised, and a division by a bare name is the same
hazard.
