# Every build flag, audited: what does nothing, what was refused on a profile that no longer exists

*Worldbuilding frame: a COMPUTE audit of the Vesper project's climate model on
this desktop. Nothing here is about the simulated planet. Measured 2026-08-22 at
645502a1, T170 on sixteen threads, on `bench/bed_t170cold`.*

`aocl-and-model-build-flags.md` measured the flag line on 2026-08-20 and closed
it. Everything it measured was measured **on the MPI build, with legmod as the
transform**. Since then SHTns has taken every per-timestep transform site, Stage
A and B moved the spectral state off per-thread stacks, `-finit-real=zero` has
left the production line, and the project runs the THREADED build. A refusal is
only as good as the profile it was taken on, and this profile has been replaced
twice over.

## The screen: does the flag change a single byte?

A flag whose removal leaves a byte-identical executable is doing NOTHING. There
is no trade-off to weigh and no benchmark to run. One build per flag, T21 omp on
four, baseline `97ca84b868bfb8e5`:

| flag | result |
| --- | --- |
| `-ffixed-line-length-132` | **identical** |
| `-cpp` | **identical** |
| `-ffpe-summary=none` | differs |
| `-ffpe-trap=invalid,zero,overflow` | differs |
| `-march=znver4` | differs |
| `-O3` | differs |

`-ffixed-line-length-132` **is removed**, and the reason is the language rather
than the build: every source is `.f90`, gfortran reads those as free form, and
the flag only sets the line length for FIXED form. It can never have applied.

`-cpp` **is kept**, byte-identical though it is. It is redundant only because
CMake supplies its own `-cpp` for the preprocessing pass its module scanner
needs -- a property of the build system, which can change, rather than of the
language. Stating the intent costs nothing.

The same screen kills three candidates without a benchmark: adding
`-mprefer-vector-width=512`, `-fno-math-errno` or `-fno-semantic-interposition`
each leaves the executable byte-identical. The first means GCC's znver4 tuning
already chooses its own vector width, the second that gfortran already assumes
no `errno` from maths functions, and the third that it does nothing for an
executable.

## The re-litigation, and two verdicts change

Paired and interleaved against the production line, four rounds, 300 steps. The
bar was fixed before the first arm ran: adopt only on a positive paired median,
the arm faster in 4 of 4 rounds, self-scatter under the 5% floor, and an
unchanged restart sha.

| arm | 2026-08-20 | now | verdict |
| --- | ---: | --- | --- |
| `-funroll-loops` | refused, -0.15% | **+11.10%** [+5.73, +11.49], 4/4, numerics unchanged | **ADOPTED** |
| `-fprefetch-loop-arrays` | refused, -4.75% | +2.53% [-0.11, +4.21], 3/4, numerics CHANGED | refused |
| `-fstack-arrays` | 0.2% | -0.15% [-0.66, +0.12], 2/4 | refused |
| `-flto` | refused, defect | **-1.60%** [-6.72, -0.52], 0/4, numerics CHANGED | refused, on a number at last |

**`-funroll-loops` is the reversal, and the mechanism is legible.** It was
refused on a build whose hot loops were legmod's Legendre sums. SHTns replaced
those, and what the flag now unrolls is the physics. 11% is four times what the
entire distance from scalar code to AVX-512 was worth on the old profile, which
is the measure of how much the hot code changed.

**`-flto` was never refused on a number.** It was refused for exposing the
`zsolars` restart over-read -- "a defect the flag exposed rather than caused" --
and CLIM-37 and CLIM-38 fixed that on 2026-08-20, so the objection had expired.
Re-measured properly it loses every round and changes the answer, so it is
refused again and this time for its own behaviour.

**`-fprefetch-loop-arrays` changes the restart**, where round three recorded it
numerics-unchanged. Its interval spans zero and it wins 3 of 4, so it fails the
bar on speed as well; but a flag that moves the answer is a numerics decision
and must not be folded into a timing comparison at all.

**`-fstack-arrays` was measured twice on purpose**, either side of CLIM-83,
because its value depends on what is on the stack and CLIM-83 took 503 MB of
per-team stack away. Before: -0.15% [-0.66, +0.12], 2 of 4. After: **-0.75%**
[-1.84, -0.25], 0 of 4, with scatter of 1.4% and 1.0%. The answer did not flip,
and it is now a clean negative rather than a null.

## `-g` is adopted, and not for speed

The bench reported +1.83% [+0.29, +6.43], which is NOISE and is not quoted as a
gain: the self-scatter was 7.2% and 7.5%, both above the floor, and `.text` is
**byte identical** with and without the flag, so no code changed and no speed
difference is possible. What it costs is 1.19 MB of executable. What it buys is
a line number in every backtrace. Twice in one session a fault had to be
re-created on a specially rebuilt binary to find out which line raised it --
once for the sNaN gate and once inside SHTns -- and that is the whole argument.

## `-ffpe-trap` stays unmasked, which closes CLIM-82

CLIM-82 asked whether the model should run SHTns under an unmasked FE_INVALID at
all, SHTns having been caught reading a lane of its own stack scratch that it
never wrote. The answer is yes, and it comes from measurements already taken
rather than from new ones.

**The lane never reaches a result.** The `snan_noinv` arm of
`the-zeroing-is-an-init-flag.md` initialised every local to a signalling NaN --
which is exactly what fills the stack region SHTns over-reads -- and ran with
FE_INVALID MASKED so a bad read would propagate a quiet NaN into the restart
instead of aborting. The restart was bit identical to the production build's.

**And in production that memory cannot hold a NaN.** `-ffpe-trap=invalid`
unmasks the invalid operation, so a NaN cannot be created by arithmetic without
trapping at the line that made it; `overflow` does the same for an infinity. The
stack under SHTns therefore holds ordinary numbers left by previous Fortran
frames, and an over-read of one is inert.

So the exposure is confined to the `poisoned` profile on the threaded build,
where `-finit-real=snan` deliberately writes the one bit pattern that trips it.
That is a documented limitation of that profile, not a reason to mask the trap
in production, and the trap's value is unchanged: in a spectral model one bad
gridpoint becomes every spectral coefficient within a few timesteps, and the
failure it prevents is a run that keeps going and reaches `assess_convergence`
as a plausible number.

Masking INVALID around the transform call was considered and refused. It would
mask real NaNs arising INSIDE the transform, which is exactly where a bad
gridpoint first shows up, so it would trade the whole value of the flag for a
theoretical exposure that is measured to be inert.

## `-O3` is replaced by `-O2`, which is the same speed and does not crash

*Measured 2026-08-24, T42 on sixteen threads, one orbit of 5850 steps from a
single restart, four interleaved rounds.*

| arm | rounds | mean |
| --- | --- | --- |
| `-O3` | 41.13, 41.98, 40.61, 41.19 | 41.23 s |
| `-O2` | 41.35, 41.42, 41.23, 41.38 | 41.35 s |

The gap is 0.3% and the `-O3` arm's own spread is 1.37 s, so this is not a
slowdown; it is two builds that run at the same speed. `-O3` was adopted on the
screen above, which only asked whether it changed a byte, and no benchmark was
ever run on it. One has now been.

`-O3` is not neutral, though, because of what its extra loop transforms do to
this particular source. Fortran evaluates BOTH sides of a `where` block and
masks only the assignment, and the physics modules lean on that: 66 calls to
`sqrt`, `log`, `log10` and `exp` sit inside masked blocks whose masked-out lanes
hold arguments the function is undefined on. `SQRT(zmu0)` under a `where(losun)`
is the plainest -- the cosine of the solar zenith angle is negative at night by
construction. Scalar code never reaches those lanes. Vectorised code does, and
with `-ffpe-trap` unmasked a discarded lane becomes SIGFPE.

That is what killed a T42 run eight orbits in. The evidence that it is the
optimisation and not the state: the same restart completes at `-O2`, at `-O1`
and with `-mno-avx512f`, and the state it starts from is ordinary -- gridpoint
temperature 244.7 to 327.2 K, humidity positive, nothing non-finite. It is not
one bad expression either. Clamping the site that faulted moved the fault to
`kuo`'s Tetens exponent; clamping all eighteen of those moved it to the
shortwave; and `-ffpe-trap=overflow`, `=invalid` and `=zero` each fire on their
own.

**`-O2` is a removal of the transforms that were firing, not a fix for the
class.** `-O2` vectorises too. A different state or a different rung can wake
the same 66 sites, and world-5a0 carries the work of making each argument safe
where it can be argued to be free -- as it can for the Tetens exponents, every
one of which is capped by `AMIN1(zqs,rdbrv)` on the next line.

Two alternatives were measured and refused. `-fno-tree-vectorize` on the whole
model costs 43% and buys the same thing `-O2` buys for nothing. The same flag on
`rainmod.f90` alone costs 2.8% and does not work at all: the faulting `exp` is a
scalar call into libm that the compiler speculates past the enclosing `if`
regardless.

## The line, after

    -O2 -cpp -ffpe-trap=invalid,zero,overflow -ffpe-summary=none
    -march=znver4 -funroll-loops -g

with `checked` adding `-fcheck=all`, `poisoned` adding
`-fcheck=all -finit-real=snan -Og`, and `build_model.py`
appending `-fdefault-real-8` for the declared precision, `-fopenmp -DOMPSHARED`
for the threaded build, and `-fno-omit-frame-pointer` for a profiling one.

## What this says about refusals generally

Three of the four re-measured arms kept their verdict and one reversed by more
than eleven points. The one that reversed is the one whose stated reason named a
part of the model that has since been replaced. **A refusal should record the
profile it was taken on**, because that is what tells a later reader whether it
still holds. The two refusals here that named a mechanism -- the Legendre
routines missing least, the prefetcher already doing the job -- are the ones
that turned out to be re-openable, and they were re-openable precisely because
they said why.
