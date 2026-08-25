# What a rung of the resolution ladder costs

*Worldbuilding frame: measurements of the Vesper project's climate model, an
ExoPlaSim fork, on this desktop. Nothing here is about the simulated planet.
Measured 2026-08-25 at 5dcdf64b, sixteen threads, `production` profile, T21/T42/
T85 at ten levels. SPAT-11.*

SPAT-8 defers cost to "record cost separately" and nothing said what a rung was
worth. The optimisation workstream -- SHTns, the thread port, single precision --
is aimed at the ladder rather than at any one rung, so its own justification
rested on a term nobody had written down.

## Which timestep is priced, and why the question has three answers

Cost at a fixed rung is proportional to steps per orbit and so to 1/dt, which
makes the timestep a factor of two in every projection. The tree carries three
per-rung timesteps and they are THREE DIFFERENT QUANTITIES sharing one name:

| where | T21 | T42 | T85 | what it is |
| --- | ---: | ---: | ---: | --- |
| `config/planet.yaml` `model.resolution_timestep_minutes` | 45 | 45 | 45 | a per-rung STABILITY CEILING, the coarsest step each rung was measured to carry |
| `docs/src/pipeline/sequencing.md` section D | 45 | 30 | 22.5 | the ESCALATION ROUTE: resolution and timestep alternate so every conversion happens at constant dt |
| `notes/audits/resolution-ladder.md` | 22.5 | 22.5 | 22.5 | one step across the ladder, so rungs differ by support and not by truncation error |

**This note prices the ESCALATION ROUTE**, because that is what will be run. The
other two are priced in the same table below, because once instructions per step
are known the rest is arithmetic and a reader should not have to redo it. The
route's timesteps are PROVISIONAL on world-37tn in the same way the ceilings
are. world-f997 carries the naming.

## The instrument, and why it is not wall clock

This machine is shared with a fan-out and with its user, so a wall-clock number
carries the machine state it was taken under and cannot be compared with a later
one. Three arms taken in three different windows earlier in the day were
discarded as contaminated for exactly that reason, and they are recorded as such
at the end of this note rather than deleted, because a number known to be
contaminated is evidence about the instrument and a deleted one is nothing.

The durable price is **retired instructions per timestep**, read with `perf stat
-e instructions` over a `make_profile_bed.py` bed at each rung -- output off, the
rung's own timestep, warm from that rung's own restart. Two lengths, 300 steps
and 900, so the startup subtracts out:

    instructions per step = (I_900 - I_300) / 600

**`OMP_WAIT_POLICY=passive` is part of the instrument and not a stray setting.**
A thread waiting at an OpenMP barrier spins by default, and a spinning thread
retires instructions in proportion to how long it waits -- which puts the machine
state straight back into the count the instrument exists to keep it out of.
Passive makes a waiting thread sleep, so the count is work. It costs wall clock,
which is why the bed's own elapsed time is not quoted anywhere below.

## The price

Six cells, each taken twice, at loads from 3.71 to 12.26 on a 32-thread host:

| rung | steps | repeats | relative range | loads |
| --- | ---: | --- | ---: | --- |
| T21 | 300 | 52861116713, 52881487960 | 0.0385% | 3.71, 3.71 |
| T21 | 900 | 157655984969, 157679073772 | 0.0146% | 3.81, 3.81 |
| T42 | 300 | 187074522921, 187099406105 | 0.0133% | 4.95, 4.79 |
| T42 | 900 | 558571140151, 558593615447 | 0.0040% | 5.13, 6.08 |
| T85 | 300 | 739210701649, 739258578830 | 0.0065% | 6.98, 8.64 |
| T85 | 900 | 2211023842924, 2211201667189 | 0.0080% | 9.93, 12.26 |

The effect being measured is a factor of fourteen between the end rungs. The
instrument's own scatter is four parts in ten thousand at worst, taken across a
three-fold range of machine load -- four orders below the effect, and it is the
pairing with the load range rather than the agreement figure alone that shows the
independence.

**Per timestep:**

| rung | instructions per step | vs T21 | vs the rung below | startup, instructions |
| --- | ---: | ---: | ---: | ---: |
| T21 | 1.7466e8 | 1.000 | -- | 4.73e8 |
| T42 | 6.1916e8 | 3.545 | 3.545 | 1.34e9 |
| T85 | 2.4531e9 | 14.045 | 3.962 | 3.30e9 |

Instructions per step are a property of the RUNG and not of the timestep, so one
measurement per rung prices every declared step. Steps per orbit are read off the
beds: 5850 at dt 45, 8774 at dt 30, 11699 at dt 22.5.

**Per orbit, by which of the three quantities is priced:**

| rung | route (45/30/22.5) | ceilings (45/45/45) | one step (22.5 throughout) |
| --- | ---: | ---: | ---: |
| T21 | 1.022e12 | 1.022e12 | 2.044e12 |
| T42 | 5.433e12 | 3.622e12 | 7.244e12 |
| T85 | 2.870e13 | 1.435e13 | 2.870e13 |
| T85 vs T21 | 28.1x | 14.0x | 14.0x |

## What that is in wall clock, and the machine state it was taken under

One window, 2026-08-25, load sampled every 15 s throughout: mean 17.15, maximum
18.46, 78 samples. The measuring job is itself 16 threads on 32 logical cores, so
that load is very largely its own. Each rung prepared once with `--run-years 0`
and again with `--run-years 2`; per-orbit is the difference over two.

| rung | dt | prepare | wall per orbit | model CPU per orbit | orbits per day |
| --- | ---: | ---: | ---: | ---: | ---: |
| T21 | 45.0 | 0.42 s | 13.45 s | 200.1 s | 6425 |
| T42 | 30.0 | 0.46 s | 75.74 s | 1155.9 s | 1141 |
| T85 | 22.5 | 0.46 s | 491.96 s | 7465.4 s | 176 |

One orbit is one local year, so the last column is local years a day.

Prepare does not scale with the rung and is not a term in any projection.

Postprocessing is not where the ladder gets expensive. Wall per orbit less model
CPU per orbit divided by sixteen leaves 0.94 s at T21, 3.50 s at T42 and 25.4 s
at T85: 7.0, 4.6 and 5.2 per cent of the orbit. The model is the cost.

## The two instruments disagree by a third at T85, and that is the finding

Per orbit on the route, instructions put T85 at 28.1 times T21. Wall puts it at
36.6 and model CPU at 37.3. The gap is a retire rate that falls with the rung:

| rung | instructions per model CPU-second, production run | same binary on the bed |
| --- | ---: | ---: |
| T21 | 5.11e9 | 5.38e9 |
| T42 | 4.70e9 | 6.56e9 |
| T85 | 3.85e9 | 5.12e9 |

On the bed the three rungs retire within five per cent of each other. In the
production run T85 retires a third fewer instructions per CPU-second than T21.
The bed differs from a production run in exactly two ways -- its output is off,
and its waiting threads sleep instead of spinning -- so the missing third is
CPU-seconds spent NOT retiring model instructions, and it grows with the rung.

**A projection composed from transform work ratios prices the instructions and
will under-predict T85's wall by about thirty per cent.** That is what SPAT-11's
own arithmetic does, and it is why the two instruments are both reported here
rather than one being reduced to the other. Which of the two terms it is --
barrier wait or output -- is not settled by these measurements;
`exoplasim/scripts/attribute_barrier_wait.sh` is the instrument for the first.

Against SPAT-11's naive composition of "roughly 240 local years a day at T85
against T42's thousand": T42's thousand holds at 1141. T85's 240 does not -- the
route gives 176. Priced at the ceilings instead, where T85 runs at dt 45, it
would be about 350.

## The configuration these numbers are for

The runs behind them declare `land_albedo_source: uniform`, `roughness_source:
uniform` and `soil_water_source: uniform`, because the staged surface families
for these rungs run through a full biosphere and baseline pass that does not
exist on this build yet, and `config/planet.yaml` is otherwise unmodified. Every
rung carries the same choice, so the comparison between them is clean; the
absolute per-orbit figures would move by whatever a field-valued land surface
costs, which is a startup read and no per-step work.

## The three contaminated arms, kept rather than deleted

Taken earlier the same day, in three different windows, before it was established
that the machine was shared. Load is known for the last of the three only.

| rung | dt | prepare + 2 orbits | load during | verdict |
| --- | ---: | ---: | --- | --- |
| T21 | 45.0 | 33.0 s | unrecorded | CONTAMINATED |
| T42 | 30.0 | 159.7 s | unrecorded | CONTAMINATED |
| T85 | 22.5 | 1069.9 s | mean 23.2, 29 samples | CONTAMINATED |

Compared with the clean window's 27.3, 151.9 and 984.4 s for the same three arms,
the contamination cost 21, 5 and 9 per cent. It does not scale evenly with the
rung, which is precisely why a rung-to-rung ratio taken across an inconsistent
set is worse than no ratio at all.

**A wall-clock instrument on this host has to record the machine state it was
taken under, or it cannot be compared against a later one.** Where a price can be
expressed in something the operating system is not free to reschedule, it should
be; the section above is that, and the contamination is what forced it.
