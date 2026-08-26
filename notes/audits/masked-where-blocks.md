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

## What a source sweep cannot settle

The pass is a text parse and the guards are argued from the source, so what is
established here is that the sites this pass can see are floored -- not that the
class is gone. It cannot see an intrinsic reached through a call from inside a
`where`, or a mask whose complement is out of domain for a reason the physics
knows and the text does not. Only a run at the
optimisation level whose transforms fire says whether anything is left, and that
is the world-bhs probe at -O3 under `-ffpe-trap`, which is world-v9j9.

**Run 2026-08-26, and it did not fire.** A T42 arm built
`build_model.py --res T42 --ranks 16 --extra-flag=-O3 --no-publish` (sha
`401aa19d`, in its own build directory) integrated a full orbit from cold --
8774 steps at dt 30, 1531 s -- with `-ffpe-trap=invalid,zero,overflow` unmasked
and took **no SIGFPE**. `lint_masked_domains.py` reports zero remaining sites
over the same source. So a full orbit of this model at the optimisation level
that woke the class does not wake it.

**That is not yet the row's PASS**, and the reason is a different defect
entirely. The arm ended in SIGSEGV at `epilog_` -> `mpputgp_` -> `mpgagp_`, the
restart write, and so produced no finite restart to check; the `-O2` companion
from the same cold start ends at the identical frame. That is
`notes/audits/epilog-adenergy-use-after-free.md` -- `adenergy` freed and then
written -- and it is not optimisation-dependent, not resolution-dependent and
not thread-count-dependent. Until it is fixed, `compare_restarts.py` has nothing
to compare: both arms' `plasim_status` end one record marker short.

What the orbit does establish is the half that a trap can deliver on its own: a
SIGFPE would have named a surviving site and none appeared. What is still owed
is the finite restart, and it costs one rebuild plus one orbit.

## The division class, and what re-deriving the population found

The pass reported a division only where the divisor was PARENTHESISED, which
made the count of 133 a measurement of the parser rather than of the class. A
division by a bare name, and one by a name with a subscript or an argument
list, are the identical hazard: the mask does not stop the lane being
evaluated, and `-ffpe-trap=zero` does not care how the divisor was spelled.
Extending the pass to all three spellings took the population to 345 sites,
156 distinct divisor keys. The masked Tetens divisors within it are settled
separately, by `ra4d` in `plasimmod.f90`, which is world-6tp.

Every key has now been read in context; `radmod.f90`'s 45 are a section of
their own below, because they and the masked locals in the same file are one
reading. Most are what the first reading supposed: a scalar or a named constant the mask has no bearing
on, a field positive on every gridpoint, or a sum bounded away from zero. Each
carries a row in `lint_masked_domains.py`'s `CLASSIFIED` table with the
argument that says why it cannot vanish.

Nine were not, and they fall into two kinds.

### The divisor is the mask, negated

`initpm` computes `sidereal_day / (TWOPI * restim)` under `where (restim >
0.0)`, and the same for `tfrc`. Both arrays carry `= 0.0` in their `plasimmod`
declarations, so on a run that sets neither the divisor is zero on EVERY lane
and only `-O2`'s not having vectorised `initpm` has kept the trap from firing.
`mkicec` divides by `2.*picedo(:)` under `where (picedn(:) < picedo(:))`; a
lane the mask keeps has `picedo > picedn >= 0`, and a lane it discards is open
water with `picedo` exactly zero. `kuo` has five of the shape written out
plainly -- `where (zpt(:) > 0.) zat(:) = zi(:)*(1.-zbeta(:))/zpt(:)` and its
four siblings. `mkcflux`'s two conduction divisors, `zhsnow/CKAPSN +
xiced/CKAPI` and `zhsnow + xiced`, are zero on any discarded lane that is open
water with no snow, which is most of them.

Each is floored, and each floor is argued at the site to be a no-op on the
lanes the mask keeps. The floor is `1.0e-30` rather than `tiny()` throughout:
`-ffpe-trap=overflow` is declared beside the zero trap, and a floor at the
smallest normal number trades one trap for the other.

### The divisor is indeterminate memory

This is the kind the first reading did not anticipate, and it is the larger
one. `tands` declares `zcap`, `zdiff`, `zsoilz`, `zsnowz`, `zsoilz1`, `zdiff1`,
`zcap1`, `zsntop` and `zctop` as automatic locals with no initialiser, and
writes every one of them only inside `where (dls(:) > 0.0)`. So on a sea lane
the divisors below -- `zsoilz1(:)`, `zctop(:)`, `snowdiff*zsoilz(:,1) +
zdiff(:,1)*zsnowz(:)` -- were evaluating against whatever the stack held.
`mktsoil` has the same shape for `ztn`, `zebs`, `zcap` and `zdiff`, and
`mkdca` for `zsum1` and `zsumq1`.

The remedy is not a floor. `tands` and `mkdca` preset the locals on every lane,
and `mktsoil` drops the mask from its elimination and back-substitution
entirely and keeps it only on the two statements that write `dsoilt`. Both
leave the kept lanes bit-identical: every preset value is overwritten before it
is read on a lane the mask keeps, every consumer outside the mask is itself
masked on the same test, and `mktsoil`'s recursion stays positive on a sea lane
for the same reason it does on a land one, because the caller now presets
`pcap`, `pdiff` and `psoilz` positive everywhere.

The presets in `mkdca` were chosen to reproduce what the lane already holds:
`zsum1 = 1.`, `zsum2 = 0.` make the masked `ztht = zsum2/zsum1` evaluate to the
`ztht(:) = 0.` two lines above it.

## radmod.f90, where both remainders lived, and why they were one reading

The file was left unread twice, in the division class and in the masked-local
class, because it was being edited on another branch and a row keyed on
argument text from a version that no longer exists is a claim about nothing.
Read once the file settled: 45 divisor keys over 103 sites, and 56 locals.

They could not be read separately. Nine of the divisors ARE locals of the second
class. `swr` computes its whole shortwave two-stream -- the transmissivities,
the reflectivities, the band absorptions, the running column amounts and the
interface fluxes -- into 53 automatics whose only definitions are inside
`where (losun(:))`, and then divides by six of them, multiplies two pairs of
them together inside `1 - r*r`, and takes `EXP` and `LOG` of others. On a night
lane every one of those was an operation on a stack word.

The remedy is the one `tands` and `mkshallow` took: a preset on every lane,
valued at what the arithmetic produces on a lane the mask discards. Here that is
the transparent atmosphere over a black surface -- zero absorber amount, unit
transmissivity, zero reflectivity, zero flux -- and it settles the divisions at
the same time, because every masked quotient then divides by exactly 1. The
adding method's `1 - zra1s*zrb1s` becomes `1 - 0.*0.`; `zto3t` and its five
siblings are 1; and the three spelled-out absorptance denominators are
`1 - A(0)/zsolar`, which is 1 because the running column amounts the block above
already zeroes on every lane now have their per-level copies zeroed too.

Nothing stored moves. On every lane `losun` keeps, each of the 53 is assigned
inside the same masked block, on the same pass, before any read of it: the four
combined-layer pairs at the preset above the downward loop, the per-level R and
T at the head of each loop body, the column copies inside the absorber loop, and
`z1mrabr`, `zscf` and the four flux locals one statement before their use.

Of the remaining divisors, most are the shapes settled elsewhere here: a scalar
the mask has no lane bearing on (`ga`, `zmu00`, `zwfit`, `zsolar1`, `zsolar2`),
a sum of a positive literal and a non-negative term (`1.+ztcon`,
`1.+0.042*zo3+...`, `1.+bb*zmu0`), a divisor already floored at its site
(`zmu0+zero`, `zaerd1`, `1.+0.816*max(0.,zmu0)`), or a field positive on every
gridpoint (`dt(:,jlev)`, `qex1(jaer)`, which `radini` aborts on if it is not).
Three needed working out rather than matching:

* `zr`, the cloud two-stream denominator, is `(u+1)^2 e - (u-1)^2/e` with
  `e >= 1` and `u >= 0`, so it is at least `4u` and vanishes only where `u`
  does, which is where `zuz` reaches zero. `zuz` is `1 - zom0*(1-2*zb2)` with
  `zom0` capped at 0.9999 one line above and `zb2` non-negative, so it is at
  least 1e-4 whenever `2*zb2 <= 1`; `zb2` is at most `tswr2/ln 3`, which is
  0.059 at the 0.065 this project leaves `tswr2` at. The `max(0.,·)` guards
  around `zu` and `zexp` are stale-lane protections, not evidence that `zuz`
  can go negative on a real one.
* `lwr`'s `1.-ztau0` and `ALOG(ztau0)` are safe for the reason its `ALOG(ztau0)`
  already was: the unconditional `AMIN1(1.-zero, MAX(zero, ·))` two lines above
  the `where` holds `ztau0` in `[1e-6, 1-1e-6]`, so the complement is at least
  1e-6 and the logarithm at most about -1e-6.
* `sigma(NLEV)-sigma(NLEM)` is the gap between the two lowest full-level sigmas,
  which the vertical coordinate makes strictly positive and which no mask
  touches.

Ten of the rows rest on the beam cosine being non-negative, which was a property
of a namelist default rather than of the source: `solang` zeroes `gmu0` and
overwrites it only where the cosine exceeds `sin(dawn)`. `radini` now refuses a
negative `dawn`, so it is a property of the source. That changes no run, because
`dawn` is 0.0 by declaration and nothing in this project writes the key, and it
removes a latent trap of its own -- below zero `swr` forms `zmu0**1.7` for the
ECHAM6 ocean albedo, a real power of a negative base.

Both counts were checked by positive control on the tree before they were
believed, which is the discipline the 133 that became 345 earned. Dropping the
inner `ztau2` floor in `swr`'s cloud chain makes exactly one new division site
appear and nothing else move; dropping the `zto3t` preset makes exactly `zto3t`
reappear in the masked-local pass.

With radmod in, `lint_masked_domains.py` has nothing unclassified in any of its
three classes, so `--kind any` is its default and its gate, and
`lint_masked_locals.py`'s DEFERRED table is empty.

### What the reading found that is not a masked-domain defect

Thirteen divisors in `swr` are not bounded on a lane `losun` KEEPS, and the
mask never had anything to do with them: the expression is the same on a kept
lane as on a discarded one. They are two shapes, they were world-2223, and the
two shapes got different answers. Measured 2026-08-26 by
`exoplasim/scripts/swr_divisor_domain.py`, which evaluates each expression over
the model's own absorber ranges rather than arguing from the formulae. It
writes nothing: what it produces is a bound, and the bound belongs here and in
the CLASSIFIED rows of `exoplasim/scripts/lint_masked_domains.py` that cite it.

#### `1 - A(u)/zsolar_b`, nine sites: the family gets the clamp

The clear-sky band transmissivity, at `zto3t`, `zto3u`, `ztwvt`, `ztwvu`,
`ztco2t`, `ztco2u` and the three denominators written out. `A` is Lacis and
Hansen's absorptance as a fraction of TOTAL incident flux and `zsolar_b` is band
b's share of it, so the quotient is the fraction of the BAND the absorber removes
and the divisor vanishes when it removes all of it.

**All three absorptances can exceed their band's share, and one of them does so
at a finite asymptote this configuration puts 24 per cent above 1.**

| absorber | supremum of `A/zsolar_b` | crosses zero at | model's own maximum | margin |
| --- | --- | --- | --- | --- |
| ozone, band 1 | unbounded; the Huggins term goes as u^0.195 | 9.134e8 cm STP | 12.89 cm STP | 7.09e7x |
| water vapour, band 2 | 1.24053, the closed form `h2osww*h2oswl*2.9/5.925/zsolar2` | 2090.3 precipitable cm | 212.89 cm | 9.82x |
| CO2, band 2 | unbounded; two logarithms | 1.096e45 atmos-cm | 5403 atmos-cm | 2.03e41x |

The model ranges are derived and not assumed: the ozone column is `mko3`'s own
`a0o3+a1o3+aco3` times `ozone_scale`, the water and CO2 columns are `swr`'s own
`zwv` and `zco2` on the bootstrap climatology's `hus`, `ta` and `ps` with the
reconstructed sigma grid, and the magnification ceiling is `zm = 35`, which the
mask's own threshold lets a kept lane reach.

**The water vapour margin is a property of the host, not of the formula.** On
the Sun, with both weights at 1 and `radmod`'s declared `zsolar2 = 0.483`, the
supremum is 1.0134 and the crossing is at 5.75e6 cm. Here it is 1.2405 and the
crossing falls to 2090 cm: five orders of magnitude of margin removed by a K
dwarf's band split and by two re-weightings that are not confined to 1. Over the
declared `h2o_sw_level_bracket` the supremum runs 1.2043 to 1.2864 and the
crossing 3271 to 1296 cm, so the verdict survives its own bracket at both ends.

**So the divisor cannot reach zero on a lane this planet's configuration
produces -- the tightest of the nine floors at 0.2016 -- and nothing structural
keeps it there.** `swr` now applies the clamp `lwr` already applies to `ztaucs`,
`AMIN1(1.-zero, AMAX1(zero, .))` at each of the nine sites, through a `zcstr`
preset to 1 OUTSIDE the mask so a discarded lane still divides by exactly 1.
Being 0.2016 away from binding, the clamp is bit-identical to its control on
this configuration: it changes stored values only where the alternative was
dividing by a number approaching zero, which is the condition it exists for.

#### `1 - R_above*R_below`, four sites: the family gets a bound

The adding method's denominator, in the downward loop, the upward loop and the
flux loop. Both factors are layer reflectivities summed from a Rayleigh term, a
cloud term and an aerosol term without the sum being bounded, and the
direct-beam cloud reflectivity `1 - 1/(1 + zb1*ztau1/zmu0)` does approach 1 as
the beam cosine falls. The divisors use the SCATTERED pair, and that one is
bounded, by caps already in the code:

- **Rayleigh.** `zscf = rcoeff*ps/101100*9.80665/ga` is 0.5452 at this planet's
  gravity and surface pressure, giving a diffuse reflectance of 0.0813. Over
  every `rcoeff` candidate and every surface pressure the model produces, the
  ceiling is 0.1003.
- **Cloud.** `zlwp` is capped at 1000 g/m2, which caps `tau1` at 138.48 and
  `tau2` at 145.84, so `zrcl1s <= 0.9449` and `zrcl2s <= 0.7369`. The band-2
  maximum is interior, at 360 g/m2, not at the cap.
- **Aerosol.** The two-stream form is strictly below 1 at any optical depth:
  `(u-1)/(u+1)` ceilings of 0.5139 in band 1 and 0.5684 in band 2 on the optical
  constants the model reads. `naerosp` is 0 in the configured runs, so `iaeron`
  is 0 and the term is absent.
- **The layer sum.** Both terms are linear in `dcc`, so the maximum is exact and
  at an endpoint: 0.9449 at `dcc = 1` in band 1, 0.7369 in band 2. The sum
  cannot exceed 1.
- **The accumulation.** Every layer satisfies R + T <= 1 by construction, since
  `ztb1u` is 1 minus the ozone absorptance times the clear fraction, minus
  `zrb1s`, minus aerosol absorption. The adding recursion is therefore bounded
  by 1, with its fixed point at a conservative layer exactly 1; ten conservative
  layers at 0.9449 accumulate to 0.9942.

**The product's floor is therefore 0.005799 in band 1 and 0.034467 in band 2
against a perfect reflector, 0.1136 and 0.4730 at the model's own maximum
band-1 surface albedo, and 0.1298 and 0.4763 at the fields the model actually
carries.** It cannot reach zero, the bound is structural, and it is recorded
rather than clamped: a clamp here would be guarding against a state the caps
already forbid.

#### Which of the thirteen is tightest depends on the currency

By divisor VALUE it is the band-1 flux-loop denominator at 0.005799, but only
with every cap in the code saturated at once and a perfect reflector under it;
at the model's own fields it is 0.1298. By margin in the quantity that can
GROW it is the water vapour trio at 9.82x, and that is the one that matters:
the family-B margin sits behind hard caps, while the family-A margin sits on a
water column that grows with the modelled climate and on three namelist
quantities none of which is confined.

## The local whose only definition is inside a mask, as a class of its own

The nine locals in `tands`, the four in `mktsoil` and the four in `mkdca` were
found as DIVISORS, and that is an accident of which pass was being run. What
they share is not division: it is that a procedure declares an array local with
no initialiser, so its contents on entry are whatever the stack held, writes it
only inside `where` blocks, and reads it. The mask does not restrict which
elements of a right-hand side the compiler evaluates, so on every lane the mask
discarded the read is a read of indeterminate memory -- whether the read is a
divisor, an intrinsic argument, a factor in a product, or a term in a sum. And
where the read sits under a DIFFERENT mask than the write, or under none, the
indeterminate lane is not discarded at all: it is stored.

`exoplasim/scripts/lint_masked_locals.py` is that class enumerated. It reports a
local when every one of its definitions is inside a `where` construct and
something reads it, over-reporting in five directions it names and blind in two
it also names, and it is a separate pass from `lint_implicit_save.py` for a
reason argued in the module: the two are exact complements, one firing on the
presence of an initialiser and the other on its absence, and deleting an
initialiser is `lint_implicit_save`'s remedy and this class's precondition.

### Verification, because a count is a measurement of the pass until it is not

The 133 masked divisions above were a measurement of the parser. So this pass is
checked against a population whose answer is known before it is believed. Twenty
reduced fixtures, each right or wrong in a named way, run on every invocation
and gate the tree pass; they cover the mask, the preset, the dummy argument, the
initialised local, the element write, the one-line `where`, the one-line `if`,
the folded continuation and the `elsewhere` arm. The positive control is the
tree itself before world-d016: scanning `landmod.f90` and `rainmod.f90` at that
commit's parent returns exactly `tands`'s nine locals, `mktsoil`'s four and
`mkdca`'s four, and every one of them clears against the fixed source.

### The population, and what it holds

Outside `radmod.f90` the pass reported 24 locals. Seventeen in `rainmod.f90`,
five in `oceanmod.f90`, one in `icemod.f90`, one in `landmod.f90`.

Twenty-two were defects and are preset. `mkshallow` is the largest and is
`mktsoil`'s shape written out again: a tridiagonal elimination under
`where (kshallow(:) > 0)` in which each level reads the previous level's `zebs`
and `zqn`, so on a lane with no shallow convection the whole recursion ran on
stack words, and the back-substitution multiplied two of them together. `mklsp`
formed the Tetens exponent and `ra4d` from an indeterminate `ztn` on every lane
that was not supersaturated. `kuo`, `mkrain`, `icestep`, `mksst`, `addfc` and
`hdiffo` are each the shorter shape: a local written at the head of a masked
block and read two statements later in the same one.

The remedy is a preset on every lane, as in `tands` and `mkdca`, and the value
is chosen to be what the arithmetic produces on a lane the mask discards: the
unchanged state where the routine computes a change (`ztn = zt`, `zqn = zq`,
`zold = zsst`), zero where it computes a difference or a flux. `icestep`'s
`zsnowold` is the degenerate case -- its masked assignment is a bare copy of
`xsnow`, so removing the mask IS the preset.

None of them moves a stored value. On every lane its mask keeps, each local is
assigned before it is read, in the same block and on the same pass, and every
consumer of every one of them is itself masked on the same test.

Two were safe and carry a row in the pass's `CLASSIFIED` table instead. Both are
the `where` / `elsewhere` over-report: `mkradv`'s `zrop` is written in both arms
of a pair whose `elsewhere` is unconditional, so the two arms partition every
lane, and `hdiffo`'s `zdtx` is written in three index ranges that between them
cover `0:NLON`, two by such a pair and the third by an unmasked copy.

### radmod.f90's 56, and the three that were the pass looking at itself

`swr` held 53 and they are preset, for the reason and by the argument the
radmod section above gives. `lwr` held three -- `zaco2`, `zao3` and `zth2o` --
and all three are over-report 2, the same shape as `mkradv`'s `zrop`: each is
written in both arms of a `where` / `elsewhere` pair whose `elsewhere` is
unconditional, so the two arms partition every lane and the local is defined
everywhere before `ztaucs` reads it. Their sibling `zah2o` is not reported at
all, because the continuum term adds an unmasked definition below the pair,
which is the pass behaving exactly as its five declared directions say it will.

With those rows in, the DEFERRED table is empty and every file gates.
