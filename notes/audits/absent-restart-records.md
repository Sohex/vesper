# What happens when a restart does not carry a record the model asks for?

**Recorded:** 2026-08-25
**Scope:** the Vesper climate model's restart reader, `restartmod.f90`, and
every caller that reads a record which may not be in the file.

Worldbuilding frame: this is a correctness audit of the code that resumes a
simulation of the planet Vesper. Nothing here is about the real world.

## Result

There are two classes, and only one of them is silent.

**The silent class had exactly one member, and it is fixed.** `landmod`'s
`landini` lowered `nexcheck` and read four records that a restart written
before LSHY-3 does not carry, relying on `mpgetgp` to leave its array alone
when the record is absent. `mpgetgp` does not: it reads into a buffer of its
own, fills that buffer only on NROOT, and scatters it into the caller's array
whichever way the lookup went. The `-1.0` sentinel written into `dwatcl` was
destroyed by the call that was supposed to leave it, so the
`if (ALL(dwatcl(:,:) < 0.0))` rebuild test was decided on stack contents and
the rebuild was skipped unless the stack happened to be negative throughout.
`adrain`, the drainage accumulator, started its output window from the same
contents. Tracked and fixed as world-onw8.

**The loud class is a design question, and it is answered below.** A record read
under the default `nexcheck` that the file does not carry stops the model with
the record's name printed. That is `get_restart_array`'s stop at
`restartmod.f90`. Nothing is silently wrong there; what was open is which
records it should apply to, and "What a restart is allowed not to carry" settles
it per record.

## The mechanism

Three routines, and the defect lives in the seam between them.

`get_restart_array` scans `yresnam(1:nresnum)` for the name. On a hit it seeks
and reads; on a miss it stops the model when `nexcheck` is 1, and otherwise
falls through to a bare `return` WITHOUT TOUCHING its argument. That last
behaviour is what invites the sentinel idiom, and for a scalar or an array the
caller owns it is sound.

`mpgetgp` breaks it. Its buffer `z(NUGP,klev)` is an automatic array with no
initialiser, it is filled only under `if (mypid == NROOT)`, and `mpscgp(z,p,klev)`
is then called UNCONDITIONALLY. The caller's array is overwritten from `z`
whatever happened, so nothing the caller wrote into it before the call
survives. Pre-setting the array cannot help, and neither can zeroing `z`: that
makes the result deterministic and still destroys the sentinel, which moves
the rebuild test's error from one direction to the other.

`mpsurfgp`, in the same file, already had the shape the fix needs.
`get_surf_array` reports whether it read anything through an `iread` argument,
`mpbci` broadcasts that to the team, and the scatter is skipped when the
answer is no. The presence answer HAS to travel that way: `yresnam` and
`nresnum` are threadprivate and only NROOT has opened the restart, so a flag
that is correct on root and garbage elsewhere would be the same class of
defect, and `mpscgp` carries omp barriers, so a team that disagreed about
entering it would not simply compute the wrong answer.

## The population: every site that lowers nexcheck

Four, and the verdict for each is what it relies on when the record is absent.

| Site | Reads | Absent-record behaviour | Verdict |
| --- | --- | --- | --- |
| `plasim.f90` `read_atmos_restart`, the stamped geometry | `get_restart_integer` x4 into NROOT-only locals | the local keeps the `-1` set before the call, and the mismatch test is written to pass an unstamped restart | sound |
| `plasim.f90` the ecological stream, at the `ecovers` marker | `get_restart_real`, `get_restart_integer` into NROOT-only scalars, then `mpbcr`/`mpbci` | `zecovers` keeps its `-1.0` and the stream's interval starts clean | sound |
| `plasim.f90` the accumulator set, at the `accuvers` marker | `get_restart_array` DIRECTLY into NROOT-only spectral accumulators, `get_restart_real` for `denergyfix`, and `has_restart_array` for the two energy-diagnostic accumulators | the array reads sit behind `if (zaccuvers >= 2.0)`, and `epilog` writes the marker AFTER the whole set, so the marker is a guarantee that every member is present; `denergyfix` keeps its declared zero; `has_restart_array` reads nothing at all | sound |
| `landmod.f90` `landini`, the layered store and the drainage | `mpgetgp` x4 | the sentinel does not survive the scatter | THE DEFECT |

No other file lowers `nexcheck`, and no other call site passes a gridpoint
array to a reader that may not find its record.

The three sound sites are sound for the same reason: what they read under a
lowered `nexcheck` is a SCALAR on NROOT, or an array behind a version marker
that is written last. Neither shape goes through `mpgetgp`, so neither is
exposed to the scatter.

## What replaced it

`restartmod` gained `has_restart_array(yn,kfound)`, a name lookup that reads
nothing and is therefore safe at any `nexcheck`. `mpimod_omp` gained
`mpgetgp_found(yn,p,kdim,klev,kfound)`, which asks on NROOT, broadcasts the
answer with `mpbci`, skips `mpscgp` when the record is absent, and returns the
answer to the caller. `landmod` branches on that answer rather than on the
array's contents, and no longer lowers `nexcheck` at all, because
`mpgetgp_found` never asks for a name the file does not carry.

`ddrain` was the fourth record and the one with no pre-set value at all: the
block set `dwatcl`, `dsoili` and `adrain` before the reads and zeroed `ddrain`
only inside the rebuild branch. It is pre-set with the others now.

Neither routine takes an OPTIONAL argument. These are external subroutines
with no explicit interface, where OPTIONAL is not conforming and `PRESENT`
cannot be relied on; a separate entry point is the conforming shape.

`mpgetgp` itself now zeroes its buffer. That is not what makes it safe -- a
zeroed buffer scattered over a sentinel is still a destroyed sentinel -- and
it is worth having only because at the shipped flag line the alternative is
arbitrary stack that no instrument catches. What keeps `mpgetgp` from being
misused again is the pair: a record that MUST be present goes through
`mpgetgp`, where `get_restart_array` stops on a miss, and a record that may be
absent goes through `mpgetgp_found`.

## The check that can fail

`exoplasim/scripts/verify_absent_restart_record.sh` compiles the real
`restartmod.f90` against a checking program at `-Og -fdefault-real-8
-finit-real=snan`, writes a restart carrying one record, and asks for a second
one that is not in it. Four criteria, fixed before it was run: an absent name
leaves `get_restart_array`'s argument alone; `has_restart_array` answers 1 and
0 correctly and moves no read position; the PRE-FIX `mpgetgp` shape delivers a
NaN into the caller's array and makes `ALL(p < 0.0)` come back false; and the
`mpgetgp_found` shape leaves the array at its pre-set value and reports 0.

The third is the control and it is the reason the other three mean anything: a
gate nobody has seen fail is not a gate. It arms because the pre-fix buffer is
never written and `-finit-real=snan` fills it. At the shipped flag line the
same read returns ordinary stack, which is the same defect without an
instrument on it, and that is what a production run was doing.

The scatter is transcribed rather than linked: `mpscgp` needs the whole of
`pumamod` and a grid, so the two reader shapes in the checking program carry
the copy `mpscgp` performs on a one-thread team. What the transcription cannot
cover -- that every thread takes the same branch -- is structural rather than
measured, and rests on `mpbci` being the same broadcast `mpsurfgp` uses.

The whole model also compiles and links under the `poisoned` profile, which is
where the new external subroutines are shown to resolve.

## What a restart is allowed not to carry

The question is not which records could be made optional. It is which SHOULD be,
and the answer turns on what an absent record costs.

**An ACCUMULATOR over an output window costs a partial window.** The first
record after the resume is a mean over fewer timesteps than it declares, and
`outreset` makes even that go away by discarding the whole interval and zeroing
the counter with it. The loss is bounded by one output window and it shows up in
one record.

**A PROGNOSTIC costs state the cell integrates from**, and no reset makes it go
away. There is nothing bounded about it and nothing in the output says it
happened.

So the default is `required`: the model stops, loudly, with the record's name.
That is affordable because of rule 7 -- no canonical climatology lineage is
declared, so every run and every climatology is disposable and the cost of a
refused restart is a re-commissioning, not a loss. `optional` has to be earned
by the first class and is refused to the second.

### The records this fork adds, and the verdict for each

Derived by comparing the `put_restart_*`, `mpputgp` and `mpputsp` call sites in
`vendor/exoplasim/exoplasim/plasim/src` against the same sites in the vendored
subtree's base commit, which is the earliest state of the model this repository
has. Pristine upstream ExoPlaSim 3.4.2 is not reachable offline -- the subtree
arrived as one squashed commit with no history -- so the baseline is the fork as
it was vendored, and three records the fork already carried at that point
(`assolu`, `asthru`, `atsolu`) are listed with the rest because the question
asked of them is the same.

| record | class | verdict | what an absent one costs |
| --- | --- | --- | --- |
| `persistt` | prognostic clock | **required** | the snowpack persistence clock restarts at zero, so a cell one step from becoming glacier goes back to the start of `glacpersec`. Silent, and bounded only by the threshold |
| `dwatcl`, `dsoili` | prognostic reservoir | optional | `landini` rebuilds the layered store from the single-layer `dwatc` a pre-LSHY-3 restart carries. world-onw8 |
| `ddrain` | prognostic | optional | zero: no drainage is in flight at the resume. world-onw8 |
| `adrain` | accumulator | optional | zero: the drainage output window starts clean. world-onw8 |
| `adenergy`, `adener3d` | accumulator | optional | `outreset`: the whole output interval is discarded. world-5qy |
| `afdsw1`, `afdsw2` | accumulator | optional, behind `accuvers` | `outreset`, as above. WORLD-3QFZ |
| `aasosp`, `aaspsp`, `aastsp`, `aasqsp`, `aasdsp`, `aaszsp` | accumulator | optional, behind `accuvers` | `outreset`. CLIM-2, CLIM-5 |
| `aorbnu`, `alambm`, `azdecl`, `ardist`, `arasc` | accumulator | optional, behind `accuvers` | `outreset` |
| the nineteen `aeco*` | accumulator | optional, behind `ecovers` | `ecoreset`: the ecological interval starts clean and says so. EFOR-2 |
| `naccueco` | counter | optional, behind `ecovers` | zero, with the accumulators it counts |
| `denergyfix` | controller state | optional | 0.0, the declared value: that segment starts uncorrected and reconverges in one window. world-fsr |
| `accuvers`, `ecovers` | marker | optional | -1.0, which is what makes a marker a marker |
| `nlat`, `nlon`, `nlev`, `nrsp` | stamped geometry | optional | -1, and the mismatch test treats an unstamped restart as unstamped rather than as wrong. world-4yz |
| `assolu`, `asthru`, `atsolu` | accumulator | **required** | nothing, if they were made optional -- and nothing is what they would buy. They sit above the `accuvers` guard and are read unconditionally, but `outreset` runs after and zeroes them on the pre-marker path, so the partial window is already discarded where it has to be. Making them optional would add a branch that changes no outcome |

`persistt` is the whole of the argument. It is the one record here where the
absent-record behaviour is a state loss rather than a lost window, and the
cheaper answer is to re-commission rather than to resume: a run whose glacier
clock starts at zero looks exactly like a run whose glaciers form late.

### The check

`restart_schema.Policy.resume` carries the verdict per record and
`check_resume_policy` holds it against the source, which
`check_policy_covers_source` calls and `scripts/smoke_test.py` runs. It fails in
both directions, and that is the point:

- a record recorded as `optional` that no reader makes optional is a decision the
  code does not implement;
- a record recorded as `required` that the model already resumes without is a
  decision nobody took.

"Makes optional" is read out of the call sites: a read through `mpgetgp_found`,
a `get_restart_*` inside a lowered-`nexcheck` region, a name that anything asks
for with `has_restart_array`, or a read behind a version marker that is itself
optional. The alternative to a check here is a table that agrees with the model
on the day it is written, and the question it answers is only ever asked at the
start of the next segment.
