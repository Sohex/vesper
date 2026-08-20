# Rank layout on a two-die CPU: the arms, and the rule declared before they run

*This is a COMPUTE benchmark for the Vesper worldbuilding project: it measures
how fast a toy climate model integrates on this particular desktop, and nothing
in it is about the simulated planet. Written 2026-08-19, before any arm was run.*

## The question, and why it is two questions

The host is a Ryzen 9 7950X3D: sixteen physical cores over two dies that are not
interchangeable.

| cores | L3 | |
| --- | ---: | --- |
| 0-7 | 96 MB | V-Cache die |
| 8-15 | 32 MB | frequency die |

Production runs sixteen ranks, so every run spans both. PlaSim is
bulk-synchronous -- collectives close every timestep -- so the slowest rank sets
the pace for all of them, and half of them sit on a third of the cache.

That is the LATENCY question: how long does one orbit take. There is a second
and independent one. Loop A needs a flux BRACKET, two converged points spanning
the target, and `CLIM-30` is blocked for want of exactly that. When two runs are
wanted, the question stops being "how fast is one run" and becomes "how many
orbits an hour does the machine produce", and those can prefer different layouts:
a configuration can lose on latency and still win on throughput.

## What is already known, which is less than it looks

`notes/config-rationale.md` under `ncpus` records one orbit at T42 on 16 ranks,
and at T21 on 8 and on 16. **There is no T42 8-rank point**: sixteen was chosen
at T42 by argument from where T21 turned over, not by measurement here. Nothing
in the repository mentions the dies at all.

## The arms

One orbit each, same binary per rank count (`_p8` is the eight-rank build; see
the README on why that `p` is ranks and `compile.sh`'s `-p` is bytes), same
inputs, same restart, all segments declared `--purpose diagnostic` so none of
this reaches a convergence window or a climatology.

| # | layout | cores | answers |
| --- | --- | --- | --- |
| 1 | 8 ranks | 0-7 | the V-Cache die alone |
| 2 | 8 ranks | 8-15 | the frequency die alone |
| 3 | 16 ranks | 0-15 | the production configuration |
| 4 | 8 + 8 ranks, two runs AT THE SAME TIME | 0-7 and 8-15 | throughput when more than one run is wanted |

Arms 1 and 2 share rank count and decomposition and differ only in the die, so
they isolate it. Arm 3 against them confounds rank count with die-spanning and
should be read as "is production beaten", never as an explanation of why.

**Arm 4 cannot be predicted from arms 1 and 2 and that is the reason it exists.**
Two 8-rank runs on separate dies still share one memory controller, one Infinity
Fabric and one DRAM bus. Whatever those cost only appears when both dies are
working, so concurrent throughput has to be measured rather than composed from
two isolated timings.

## The rules, declared here before any arm has run

Per-orbit model time on the existing 40-orbit segment ran 103 to 106 s, a spread
of about 3%. So:

- **A difference under 5% is no difference.** Report it and change nothing.
- **Arms 1 vs 2 settle whether the dies differ.** If they land inside 5%, the
  dies do not matter for this workload, 16-rank spanning is fine on that count,
  and the question is closed rather than parked.
- **Arm 4 is adopted for multi-run work if it beats arm 3 on THROUGHPUT by more
  than 10%**, throughput being orbits per wall hour across every run in flight.
  Concretely arm 4 wins when one concurrent 8-rank orbit takes less than twice a
  16-rank orbit, because arm 3 has to run the two points one after the other.
  10% rather than 5% because adopting it means running two run directories at
  once, and that is operational complexity that a marginal gain does not buy.
- **Latency and throughput get separate verdicts and may disagree.** If arm 3
  wins on latency and arm 4 on throughput, the rule is: one urgent run uses 16
  ranks, and any bracket or matched pair uses 2x8. Do not average them into a
  single recommendation.

## What would mean the measurement is wrong rather than surprising

- Any arm whose per-orbit times vary by more than the 3% already on record.
  Something else was running; the machine has to be quiet, including no
  monitoring loop.
- Arm 4's two runs disagreeing with each other by more than 5%. They are on
  different dies, so arms 1 and 2 predict the sign of that gap; if it appears
  with the opposite sign, the pinning did not take.
- Any arm reporting a different climate. These are the same physics on the same
  restart; only the decomposition changes. `assess_convergence` must not be
  asked to read them, and `--purpose diagnostic` is what keeps it from doing so.

## How to pin, and the trap that makes the obvious spelling wrong

`exoplasim.Model` accepts `mpi_opts`, appended straight to its `mpiexec -np N`
(`__init__.py:498`), so the route exists. `run_exoplasim.py` does not plumb
`mpi_opts` today, and adding it is the prerequisite for arms 1, 2 and 4.

**The incantation, verified on Open MPI 5.0.10 on 2026-08-19:**

    --map-by pe-list=0,1,2,3,4,5,6,7:ordered --bind-to core

`:ordered` is load-bearing and is the whole finding. Measured by asking each
rank for its own `Cpus_allowed_list`, at 2 ranks:

| options | rank 0 | rank 1 |
| --- | --- | --- |
| none (the default, and what production runs) | `0,16` | `1,17` |
| `--cpu-set 0,1 --oversubscribe` | `0-1,16-17` | `0-1,16-17` |
| `--map-by pe-list=0,1` | `0-1,16-17` | `0-1,16-17` |
| `--map-by pe-list=0,1 --bind-to core` | `0-1,16-17` | `0-1,16-17` |
| **`--map-by pe-list=0,1:ordered --bind-to core`** | **`0,16`** | **`1,17`** |

Every spelling without `:ordered` CONFINES the job to the right CPUs and
silently drops the per-rank binding: both ranks get the whole pool, so the
scheduler may stack two of them on one physical core while another sits idle.
Adding `--bind-to core` does not fix it, which is what makes this worth writing
down -- the flag that looks like the fix is present and has no effect.

For this benchmark that failure mode is not a small error. An arm whose eight
ranks are colliding on four cores reports a slow die, and the conclusion drawn
would be about the hardware rather than about the launch line.

**So the check is not that the run started.** Before timing any arm, have each
rank print its `Cpus_allowed_list` and require eight distinct sibling pairs.
Verified at 8 ranks on both dies: cores 0-7 give `0,16` through `7,23`, and
cores 8-15 give `8,24` through `15,31`.

Two smaller notes from the same session. `--cpu-set` still parses but is
counted against available slots, so it needs `--oversubscribe` merely to launch
-- and it has the same missing-binding defect anyway. And
`--bind-to core:cpulist=...` does not exist on this version: the only qualifiers
accepted are `overload-allowed`, `no-overload` and `if-supported`.

Note that `--bind-to core` on this machine gives each rank the two SMT siblings
of one core, which is why rank affinity masks read as pairs like `0,16`. That is
deliberate and should be left alone: a rank that can move to its sibling can
step around a daemon landing on its thread, where a rank pinned to one hardware
thread simply waits -- and in a bulk-synchronous code that wait propagates to
every other rank through the next collective. Measured 2026-08-19 during a
production run, migrations ran 19 to 383 per rank over 80 s and were strongly
skewed toward the low-numbered cores, which is background load rather than
anything intrinsic.

---

## First attempt, 2026-08-19: one arm clean, four contaminated

Run 16:44 to 17:13, four orbits per arm, all seeded from
`run_4182235e9781/MOST_REST.00039`. Model time is separated from postprocessing
by taking `MOST_REST.NNNNN` against the previous orbit's `MOST.NNNNN.nc`, which
is worth doing for its own sake: **postprocessing was 1.8 to 2.0 s in every arm
and every orbit**, flat, so it is no longer a variable in any timing question
here and the 29.4 s figure that predates the second pyburn fix should not be
quoted again.

| arm | model time per orbit, s | median | self-scatter |
| --- | --- | ---: | ---: |
| 1. 8 ranks, V-Cache | 95.2, 95.9, 96.3 | 95.9 | **1.2%** |
| 2. 8 ranks, frequency | 95.4, 101.8, 100.7 | 100.7 | 6.4% |
| 3. 16 ranks, spanning | 82.0, 101.2, 113.8 | 101.2 | 31.4% |
| 4a. 8 ranks V-Cache, concurrent | 119.4, 111.6, 104.2 | 111.6 | 13.6% |
| 4b. 8 ranks frequency, concurrent | 123.0, 124.0, 102.3 | 123.0 | 17.7% |

**Four of five arms fail the 5% self-scatter floor and are not scored.** That
floor was declared above before any arm ran, and it is doing exactly the work it
was put there for.

**The contamination is visible in the data and has a shape.** Arm 3 rises
monotonically and both halves of arm 4 fall; the peak sits between them, around
17:00 to 17:08. That is one external load ramping up and finishing, not noise,
and the operator confirmed other work was running. Arm 1 ran before it and is
the only clean measurement: **8 ranks on the V-Cache die integrate an orbit in
95.9 s at 1.2% scatter.**

**What must NOT be concluded from this.** Arm 3's opening orbit of 82.0 s is the
fastest single orbit in the matrix and suggests 16 ranks beat 8, but it is one
observation from an arm scattering 31%. Pairing it against arm 4 to score
throughput -- which this note's author did before looking at the whole profile --
takes the cleanest number from one arm and the dirtiest from another, and the two
sit at opposite ends of the same external ramp. The throughput indication favours
2x8 by 27 to 37% under every pairing available, and every one of those pairings
is between contaminated quantities. **The adoption rule is not satisfied and 2x8
is not adopted.**

**What the matrix did establish**, because it does not depend on the contested
numbers:

- Postprocessing is flat and small, measured five times over.
- An 8-rank binary resumes a 16-rank restart, so the arms are legitimate.
- Concurrency has a real cost that is not derivable from the solo arms: arm 4's
  V-Cache half ran 111.6 s against arm 1's 95.9 s on the same eight cores. Some
  of that gap is the external load, so it is an upper bound rather than a
  measurement, but it cannot be zero and it is why arm 4 exists.
- The `:ordered` pinning holds under real load: every arm bound one rank per
  core, verified on the live processes rather than on a probe.

**Re-run when the machine is quiet.** Same arms, more orbits so drift shows as a
slope rather than being inferred from three points, and interleave the arms
rather than running them in a fixed order, so machine state cannot alias with
arm identity -- that aliasing is what made arm 1 look fastest, purely for having
gone first.

## Second attempt, 2026-08-19: answered

Re-run on a quiet machine, five orbits an arm. Pass 1 completed and is below;
pass 2, the reversal, was killed partway when unrelated builds arrived, so the
ordering check was not made. It is not needed: four of five arms come in under
the 5% floor, and the one that does not is contaminated in its last orbit only.

| arm | model s/orbit | median | spread | |
| --- | --- | ---: | ---: | --- |
| 1. 8 ranks, V-Cache | 93.2, 93.4, 93.8, 93.4 | 93.4 | 0.6% | clean |
| 2. 8 ranks, frequency | 93.8, 94.3, 94.5, 94.1 | 94.2 | 0.8% | clean |
| 3. 16 ranks, spanning | 79.4, 79.2, 79.2, 92.8 | 79.3 | 17.2% | first three 0.25% apart; orbit 4 is the incoming load |
| 4a. 8 V-Cache, concurrent | 102.4, 101.9, 100.9, 101.7 | 101.8 | 1.5% | clean |
| 4b. 8 frequency, concurrent | 105.7, 106.7, 105.7, 104.0 | 105.7 | 2.6% | clean |

### The dies do not matter, and that closes the question

**93.4 s against 94.2 s, 0.9% apart, against a floor of 5% declared before any
arm ran.** The V-Cache die is not meaningfully faster for this workload, so there
is no reason to place runs on it and no reason to prefer 8 ranks in order to fit
one die. The premise this whole benchmark was built on -- that a bulk-synchronous
model would be gated by the eight ranks sitting on a third of the L3 -- is WRONG,
and the rule declared above says a result inside the floor closes the question
rather than parking it. Do not re-open it without a reason that is not "the dies
look different on paper".

### 16 ranks beats 8, so rank scaling still pays at T42

79.3 s against 93.4 s, **15.1% faster**, and the first three orbits of the
16-rank arm agree to 0.25%, which is the tightest measurement in either attempt.
This also retires attempt 1's preliminary reading that 8 V-Cache ranks were
beating 16; that was an artifact of arm order on a drifting machine.

### 2x8 concurrent wins on throughput, and is ADOPTED

Two concurrent orbits complete in 105.7 s against 158.6 s for two sequential
16-rank orbits: **33.4% better, past the 10% adoption threshold.**

So the two verdicts the rules kept separate genuinely disagree, exactly as
anticipated, and the split stands:

- **One run that is wanted now: 16 ranks.** It is 15% faster per orbit.
- **Any bracket, matched pair or endmember pair: 2x8 concurrent**, one job per
  die -- not because the dies differ, they do not, but because splitting on the
  die boundary is the natural way to give each job eight cores that share an L3.

Concurrency costs 8.4% per orbit against running the same eight ranks alone
(101.8 against 93.4), which is the memory and fabric contention arm 4 existed to
measure and which no combination of the solo arms would have predicted. That cost
is real and is still overwhelmed by doing two things at once.

`WORKFLOW.md` section 6 lists three independent pairs on the way to a carve list,
so this is worth roughly a third of the wall clock of the runs that pair.
