# Does SMT help a spectral GCM? The arms, declared before any of them ran

*This is a COMPUTE benchmark for the Vesper worldbuilding project: it measures
how fast a toy climate model integrates on this particular desktop, and nothing
in it is about the simulated planet. Written 2026-08-21, before any arm was
run.*

## Why reopen this

`docs/src/reference/config-rationale.md` rules out 32 ranks with an argument
rather than a measurement: *"those are SMT threads and MPI ranks on sibling
threads contend for the same FPU and cache, which does not help a compute-bound
spectral model."* Nobody has ever run it.

Two things now undercut the premise.

**The model is not purely compute-bound.** `spectral-transform-profile.md` puts
MPI at 11 to 27 percent of samples across ranks at T127, and most of it is
`opal_progress` -- ranks SPINNING at collectives, not moving bytes. At T42 it is
worse: 39% of samples against 14.5% on the busiest rank, which is roughly a
quarter of the machine idling on load imbalance. SMT is the mechanism that fills
stalls, and the argument against it is weakest exactly where the workload has
slack. This one has a lot.

**Concurrency already pays.** `rank-layout-benchmark.md` measured two concurrent
8-rank jobs beating two sequential 16-rank jobs by 33% on throughput, and
adopted it. Oversubscription is the untested continuation of a direction that
already won.

**What stays closed, and what is a genuinely different question.** That note
measured the two dies as indistinguishable, 93.4 s against 94.2 s, and said not
to reopen it without a reason that is not "the dies look different on paper".
A footprint argument looked like one -- 6 MB of Legendre weights per rank at
T127 against 4 MB of L3 per core on the frequency die -- and the cache-miss
profile refuted it: the Legendre routines miss LESS than anything else in the
model. **At one job per die, the question stays closed.**

It measured one job per die. Two jobs SHARING one die's L3 is a different cache
regime, not the same question asked again: CCD0 carries 96 MB and CCD1 32 MB,
and doubling the jobs on a die doubles the footprint competing for it. So the
matrix carries a 2x2 that attributes any contention penalty to a die rather than
confounding the two:

|  | one job | two jobs |
| --- | --- | --- |
| CCD0, 96 MB L3 | `T0_1x8` | `T2b_2x8_core` |
| CCD1, 32 MB L3 | `T0b_1x8_ccd1` | `T2c_2x8_ccd1` |

If `T2c/T0b` is worse than `T2b/T0`, the smaller L3 is paying for the contention
and the die matters once a die is shared. Expect it to be resolution-dependent
and to show at T127 if anywhere: the weights are 0.23 MB per rank at T42, where
sixteen ranks fit either die, and 6 MB per rank at T127, where eight want 48 MB
against CCD1's 32.

## The host

Physical core `i` owns hardware threads `i` and `i+16` (`lscpu -e=CPU,CORE`).
Cores 0-7 are the V-Cache die, 8-15 the frequency die.

## The arms

Each arm is N concurrent jobs with declared pinning. `core` binding gives a rank
both siblings of one core, which is what an unpinned launch does and therefore
what production runs; `hwthread` gives it one thread, which is the only way to
put two ranks on one core deliberately.

| arm | jobs x ranks | cpus | what it answers |
| --- | --- | --- | --- |
| `L1_1x16` | 1 x 16 | cores 0-15 | the incumbent. LATENCY control |
| `L2_1x32` | 1 x 32 | threads 0-31 | does SMT help LATENCY |
| `T0_1x8` | 1 x 8 | cores 0-7 | control for T2: one job on eight cores |
| `T1_2x8_dies` | 2 x 8 | cores 0-7 / 8-15 | the incumbent for THROUGHPUT |
| `T2_2x8_smt` | 2 x 8 | threads 0-7 / 16-23 | SMT isolated: two jobs on the SAME eight cores |
| `T3_2x16` | 2 x 16 | threads 0-15 / 16-31 | SMT for throughput, whole machine |
| `T4_4x8` | 4 x 8 | four thread groups | throughput at finer granularity |

T2 against T0 is the clean SMT question with the die held constant: same eight
physical cores, one job or two. T3 and T4 use the whole machine and are the
throughput candidates.

Resolutions: T42 first because it is cheap and validates the harness, then T127
because that is where the work is heading. T21 is excluded -- 32 ranks would give
it one latitude per rank.

## The gotcha that would have refused SMT for the wrong reason

`rank-layout-benchmark.md:135` is explicit, and the first cut of this matrix
ignored it:

> `--bind-to core` on this machine gives each rank the two SMT siblings of one
> core [...] That is deliberate and should be left alone: a rank that can move to
> its sibling can step around a daemon landing on its thread, where a rank
> pinned to one hardware thread simply waits -- and in a bulk-synchronous code
> that wait propagates to every other rank through the next collective.

It measured 19 to 383 migrations per rank over 80 s during a production run, so
stepping around background load is routine rather than hypothetical.

Every SMT arm above binds to `hwthread`, which is exactly the pinned-to-one-
thread case that describes. That handicaps all four of them against the
core-bound controls for a reason that has nothing to do with SMT, and the matrix
would have refused oversubscription on a confound.

So each oversubscribed arm gets a MIGRATION-TOLERANT twin: the same two ranks
per core, each still free to use either of that core's threads.

| arm | how | what it isolates |
| --- | --- | --- |
| `L2b_1x32_core` | one 32-rank job, `--map-by core --bind-to core:overload-allowed` | SMT latency, ranks can migrate |
| `T2b_2x8_core` | two 8-rank jobs, BOTH on cores 0-7, `--bind-to core` | SMT throughput on one die, ranks can migrate |
| `T3b_2x16_core` | two 16-rank jobs, BOTH on cores 0-15 | SMT throughput, whole machine |
| `T4b_4x8_core` | four 8-rank jobs, two per die | same, finer granularity |

The `a` and `b` pairs differ in ONE thing, so the gap between them prices the
migration freedom directly -- which is worth knowing on its own account, since
it is the first time that claim has been measured rather than argued.

Open MPI refuses to put two ranks on a core unless asked: the binding needs the
`overload-allowed` qualifier, which is one of the three that version accepts.

## The second dimension: does a waiting rank spin or step aside

`mpi_yield_when_idle` decides whether a rank blocked in a collective spins or
calls `sched_yield`. It defaults to FALSE, and Open MPI enables it by itself
only "when oversubscribing nodes" -- when ONE mpiexec has asked for more ranks
than slots.

**That auto-detection does not fire for the arms that need it most.** T2, T3 and
T4 are several INDEPENDENT mpiexec invocations of 8 or 16 ranks each, and every
one believes it has the machine to itself. So yield stays off and two aggressive
spinners share every physical core. Measuring SMT that way measures it at close
to its worst case, and would refuse it for the wrong reason.

It is therefore a declared dimension rather than a default, and the profile is
why it could decide the answer: MPI is 11 to 27 percent of samples at T127 and
most of it is `opal_progress` -- spinning, not moving bytes. A spinning rank on a
shared core is stealing issue slots from a sibling doing real work; a yielding
one is not.

Every arm runs twice, with `OMPI_MCA_mpi_yield_when_idle` at 0 and at 1. The
cheap resolution decides whether the setting matters at all; only the better of
the two is carried to the expensive one.

## Two Open MPI traps, both found by the pinning check before any arm was timed

Recorded because both fail SILENTLY in the direction of a plausible wrong
number, which is the same shape as the `:ordered` finding that motivated the
check in the first place.

**`pe-list` refuses ranges.** `pe-list=0-7` does not mean cpus 0 to 7; it fails
with "not enough CPUs in the specified PE-LIST". It has to be spelled
`0,1,2,3,4,5,6,7`.

**The PE index space changes meaning with the binding unit.** Under
`--bind-to core`, PE `p` is physical core `p` and the rank gets both its
threads. Under `--bind-to hwthread`, PE `p` is `(core p//2, thread p%2)`, so
cpu `p//2 + 16*(p%2)`. Verified directly: `pe-list=0,1 --bind-to hwthread` puts
two ranks on cpus 0 and 16, the two threads of ONE core, while `pe-list=0,16`
puts them on cpus 0 and 8, two different cores. An arm written in cpu numbers
lands somewhere else entirely and looks fine doing it -- which is exactly what
happened here, and what the check caught.

**A single job with more ranks than cores will not launch at all** without
`--oversubscribe`, because Open MPI counts one slot per physical core. It raises
the slot count and nothing else; 32 ranks still bind to 32 distinct cpus.

## The rules, declared before any arm ran

**Under 5% is no difference.** `rank-layout-benchmark.md`'s floor and
`aocl-and-model-build-flags.md`'s, not renegotiated here.

**Latency and throughput get SEPARATE verdicts and may disagree.** That is the
whole reason both are measured. Latency is seconds for one job; throughput is
jobs finished per wall hour across everything in flight.

- **`L2_1x32` is adopted for latency** if it beats `L1_1x16` by more than 5%.
- **A throughput arm is adopted** if it beats the incumbent `T1_2x8_dies` by
  more than 10%. Ten rather than five because adopting one means running several
  run directories at once, and that is operational complexity a marginal gain
  does not buy. This is the existing note's threshold and its reasoning.

**Pinning is verified before any arm is timed.** Every rank is asked for its own
`Cpus_allowed_list` and checked: the right number of CPUs for the binding unit,
inside the declared pe-list, and no two ranks sharing. A mismatch ABORTS. This
is not ceremony -- the existing note's central finding is that every spelling
without `:ordered` confines the job correctly and silently drops the per-rank
binding, so ranks stack and an arm reports a slow layout when the launch line
was at fault.

## What would mean the measurement is wrong rather than surprising

- Any arm whose wall time varies by more than the 5% floor across rounds.
  Something else was running.
- Concurrent jobs within one arm disagreeing with each other by more than 5%.
  They have identical work and identical core counts; a gap means the pinning
  did not take the way the check thought it did.
- A restart sha differing between arms **of the same rank count**. Across
  different rank counts the shas are EXPECTED to differ -- `NLPP` changes, so the
  latitude sum is grouped differently and the rounding with it. That is not a
  defect and must not be reported as one.

## What this decides beyond itself

The rank count sets `NLPP`, and `NLPP` is an input to the symmetry work: the
paired-latitude decomposition needs `NPRO` to divide `NLAT/2`. At 32 ranks that
holds for T42, T85, T127 and T170 and FAILS for T21, which would fall back to
the unpaired path. So a 32-rank production decision narrows where the Phase 2
and 3 work applies, and that is a reason to settle the layout before building
against it rather than after.

---

## T42, measured 2026-08-21

Cold bed, 600 steps, four rounds an arm, arm order rotated each round, quiet
machine, registered binaries. Every arm's self-scatter is inside the 5% floor
(0.2 to 3.6%), so all of it is a result.

| arm | jobs x ranks | latency s | beds/hour | vs incumbent |
| --- | --- | ---: | ---: | ---: |
| `L1_1x16` | 1 x 16 | **8.95** | 402 | -- |
| `L2_1x32` | 1 x 32 | 19.07 | 189 | |
| `L2b_1x32_core` | 1 x 32 | 19.23 | 187 | |
| `T0_1x8` | 1 x 8 | 9.33 | 386 | -40.9% |
| `T0b_1x8_ccd1` | 1 x 8 | 9.61 | 375 | -42.7% |
| `T1_2x8_dies` | 2 x 8 | 10.67 | **653** | incumbent |
| `T2_2x8_smt` | 2 x 8 | 15.34 | 468 | -28.3% |
| `T2b_2x8_core` | 2 x 8 | 15.18 | 473 | -27.6% |
| `T2c_2x8_ccd1` | 2 x 8 | 15.28 | 470 | -28.0% |
| `T3_2x16` | 2 x 16 | 16.15 | 439 | -32.8% |
| `T3b_2x16_core` | 2 x 16 | 16.12 | 438 | -33.0% |
| `T4_4x8` | 4 x 8 | 20.06 | 695 | **+6.4%** |
| `T4b_4x8_core` | 4 x 8 | 20.03 | 695 | +6.4% |

### SMT is refused for latency, and the number is confounded

32 ranks takes 19.07 s against 16 ranks' 8.95 s: **113% slower**, far outside any
floor. But it is NOT a clean SMT result and must not be quoted as one. At T42,
32 ranks is `NLPP = 2` latitudes a rank, which is over-decomposed whether or not
threads are shared -- the same effect that made 8 ranks beat 16 at T21 in
`rank-layout-benchmark.md`. The arm conflates SMT with too many ranks for the
resolution.

That confound weakens with resolution: 32 ranks is 6 latitudes a rank at T127 and
8 at T170. **The latency question is therefore reopened at the higher
resolutions rather than settled here.**

### Throughput: the incumbent holds

`T4_4x8` at +6.4% is the only arm to beat two-jobs-per-die, and the adoption
threshold is 10%. Keep the incumbent. Everything that puts two jobs on ONE die's
cores loses 28 to 33%, which is the plainest reading in the table: sharing cores
costs far more than sharing a machine.

### Both added dimensions are null, which is what the cheap resolution was for

**Migration freedom buys nothing here.** The `a`/`b` pairs differ only in whether
a rank is pinned to one thread or free to use both of its core's: 15.34 against
15.18, 16.15 against 16.12, 20.06 against 20.03. All inside the floor.

That does not refute `rank-layout-benchmark.md:135`, whose concern was
explicitly a daemon landing on a rank's thread under BACKGROUND LOAD. It says
the handicap does not appear on a quiet machine, which is the condition every
arm here runs under. The dimension is dropped from the expensive resolutions,
and the original note's advice stands for production.

**`mpi_yield_when_idle` buys nothing either.** Every arm is within about 2% of
its spin-pass twin, inside the floor: `L1` 8.95 against 9.07, `T4` 695 against
689, `T2b` 473 against 465. So the spinning that the profile shows -- 11 to 27
percent of samples in `opal_progress` -- is not recoverable by yielding, at
least at this resolution. Dropped from the expensive resolutions too.

### The dies are indistinguishable, as the footprint predicts

One job: 9.33 s on CCD0 against 9.61 s on CCD1, 3% and inside the floor. Two
jobs sharing one die: 473 against 470 beds/hour, half a percent. So the
contention penalty is the SAME on both dies and none of it is attributable to
L3.

That is the expected answer HERE and is not the test: T42's weights are 0.23 MB
a rank, so eight ranks want 1.8 MB and both dies hold it with room to spare. The
discriminating resolution is T127, where eight ranks want 50.7 MB -- inside
CCD0's 96 MB and outside CCD1's 32 MB. That arm is what decides the question.

### What runs at the expensive resolutions

With yield and migration both null, the cross-product collapses to one setting
and the arms that still ask something:

- **T127**, the die discriminator: `L1`, `L2`, `T0`, `T0b`, `T1`, `T2b`, `T2c`,
  `T4`. The 2x2 on the dies is the reason this resolution is run at all.
- **T170**, where the layout decision applies: `L1`, `L2`, `T1`, `T4`.

## T127 and T170, measured 2026-08-21

Cold beds, 300 steps, four rounds an arm, one setting: yield and
migration-freedom were both null at T42 and are not dimensions here.

### T127 -- the die discriminator

| arm | jobs x ranks | latency s | beds/hour | spread |
| --- | --- | ---: | ---: | ---: |
| `L1_1x16` | 1 x 16 | **39.67** | 90.7 | 0.9% |
| `L2_1x32` | 1 x 32 | 65.51 | 55.0 | 2.8% |
| `T0_1x8` CCD0 | 1 x 8 | 46.19 | 77.9 | 0.9% |
| `T0b_1x8_ccd1` | 1 x 8 | 51.59 | 69.8 | 1.3% |
| `T1_2x8_dies` | 2 x 8 | 63.73 | **106.8** | 0.5% |
| `T2b_2x8_core` CCD0 | 2 x 8 | 83.40 | 86.3 | 0.6% |
| `T2c_2x8_ccd1` | 2 x 8 | 105.87 | 68.0 | 1.5% |
| `T4_4x8` | 4 x 8 | 137.69 | 95.2 | 1.5% |

### T170

| arm | jobs x ranks | latency s | beds/hour | spread |
| --- | --- | ---: | ---: | ---: |
| `L1_1x16` | 1 x 16 | **85.71** | 42.0 | 0.6% |
| `L2_1x32` | 1 x 32 | 128.08 | 28.1 | 3.9% |
| `T1_2x8_dies` | 2 x 8 | 149.55 | **44.8** | 0.6% |
| `T4_4x8` | 4 x 8 | 336.40 | 37.8 | 1.1% |

The first `L2_1x32` at T170 came back with 36.1% self-scatter and was voided by
the guard rather than reported; the row above is a six-round re-run.

## THE DIES DIFFER AT T127, and the earlier closure was resolution-bound

| | one job | two jobs sharing it | gain from sharing |
| --- | ---: | ---: | ---: |
| CCD0, 96 MB L3 | 77.9 beds/hr | 86.3 | **+10.8%** |
| CCD1, 32 MB L3 | 69.8 | 68.0 | **-2.6%** |

Single job, CCD0 is **10.5% faster**: 46.19 s against 51.59 s, twice the floor.
And the contention behaviour diverges rather than merely differing in size --
sharing CCD0 gains eleven percent, sharing CCD1 LOSES. A thirteen point swing,
attributable to L3 and nothing else, since the 2x2 holds everything else fixed.

It lands where the footprint says it must. Eight ranks at T127 want 50.7 MB of
Legendre weights: inside CCD0's 96 MB and outside CCD1's 32 MB. T42 wanted
1.8 MB and both dies held it, which is why `rank-layout-benchmark.md` measured
them 0.9% apart and closed the question. **That closure was correct for the
resolution it was taken at and does not survive to T127.**

**The production consequence is sharper than the arm.** `L1_1x16` spans both
dies, so eight of its sixteen ranks sit on CCD1 with 50.7 MB against 32 MB. The
model is bulk-synchronous and the slowest rank paces every other one through
the next collective, so a T127 run is already being paced by ranks that are
thrashing. The table shows it: sixteen ranks beat eight-on-CCD0 by 14%
(39.67 against 46.19) where doubling the ranks should approach twice.

This also retires an argument made earlier in this project and got wrong. The
L3-thrashing hypothesis was raised, a cache-miss profile appeared to refute it --
the Legendre routines miss less than anything else in the model -- and it was
dropped. That profile was taken at T42-scale footprints. The hypothesis was
right and the refutation was measuring the wrong resolution.

## SMT is refused, and the confound is gone

| | 1x32 against 1x16 | latitudes a rank at 32 |
| --- | ---: | ---: |
| T42 | -113% | 2 |
| T127 | -65% | 6 |
| T170 | **-49.4%** | 8 |

The penalty shrinks exactly as the over-decomposition confound weakens, which is
what the T42 note predicted would happen and is the reason the question was
reopened at higher resolution rather than settled cheaply. It never approaches
parity. **So this is a genuine refusal of SMT for latency, not an artefact of
too few latitudes a rank.**

Two ranks per core is refused on a second, independent ground as well: a restart
written at sixteen ranks is not readable at thirty-two. `NESP = NSPP * NPRO`
with `NSPP = ceil(NRSP/NPRO)` gives 1904 at sixteen and 1920 at thirty-two, and
the reader runs past the end of the record. Eight ranks CAN resume a sixteen-rank
restart, so this is an asymmetry rather than a general limit -- but it means
adopting thirty-two ranks would strand every existing run.

## Throughput: the incumbent holds, and oversubscription gets worse with resolution

`T4_4x8` against two-jobs-per-die: **+6.4% at T42, -10.9% at T127, -15.5% at
T170.** It never cleared the 10% adoption bar and it turns negative exactly
where the cache pressure arrives. Consistent with the die result and with the
same mechanism.

## The verdicts

- **Latency: 16 ranks, spanning.** Unchanged. 32 is refused at every resolution
  and on the restart incompatibility besides.
- **Throughput: two concurrent 8-rank jobs, one per die.** Unchanged, and by a
  wider margin at high resolution than at T42.
- **The dies are NOT interchangeable at T127 and above.** New. A single job that
  must fit one die belongs on CCD0, and the 2x8 throughput layout should keep
  putting one job on each die rather than both on either.

## What this leaves open

The die result raises a question this matrix cannot answer: a 16-rank run gives
every rank the same number of latitudes, and the two dies do not have the same
cache to run them in. An uneven decomposition -- fewer latitudes on the CCD1
ranks -- would balance against the cache rather than against the core count. It
is a real option and it is not free: `mpimod` scatters equal `NHOR` blocks, so
uneven latitudes per rank is a deeper change than the paired-latitude work.
Worth a task rather than a paragraph.
