# The physics filter is load-bearing, and the model traps without it

*Worldbuilding frame: a numerical experiment on the Vesper project's climate
model. Nothing here is about the simulated planet. Measured 2026-08-23 at T42,
16 ranks, from a common restart.*

CLIM's unbooked-filter row proposed one experiment: measure `denergy26 -
denergy27` with `NFILTER = 2` and with `NFILTER = 0`, everything else identical,
and see whether the spectral core's apparent energy source survives the filter
being switched off. **The second arm does not exist.** The model takes a
floating-point trap on the first radiation call without the filter.

## The arms

Both from `MOST_REST.00039` of `run_4182235e9781`, same binary, same superseded
surface under `--superseded-surface-ok`, `NENERGY = 1`, `NLOWIO = 0`,
`MPSTEP = 45.0`. The configs are generated from `config/planet.yaml` at launch
and differ in `physics_filter` alone.

| arm | `physics_filter` | outcome |
| --- | --- | --- |
| A | `gp\|exp\|sp` | two orbits, clean |
| B | `""` | SIGFPE on the first radiation call |

Arm B, three times, identically:

    Program received signal SIGFPE: Floating-point exception
    #3  swr_      at radmod.f90:2486
    #4  radstep_  at radmod.f90:1486
    #5  gridpointd_ at plasim.f90:3575

`swr` is the shortwave. The build carries `-ffpe-trap=invalid,zero,overflow`, so
this is the trap firing at the instruction that made the bad value rather than a
NaN propagating silently into the state -- which is the whole argument for
keeping that flag, working.

**It is not a slow drift into instability.** No output record is written at all:
the crashed directory holds `MOST_REST.seed` and nothing else. The failure is
inside the first radiation step.

## What that settles, and what it does not

SETTLED: the filter is not inherited cargo. It came from ExoPlaSim's
tidally-locked presets rather than from a decision here, and the natural reading
of that provenance was that it might simply be wrong for a planet that is not
tidally locked. It is not wrong. At T42, on this terrain, at a 45-minute step,
the model does not run without it, so "turn it off and price the difference" is
not a question this configuration can be asked.

The mechanism is consistent with what the filter is for -- Gibbs ringing off
sharp gridpoint structure reaching a routine that cannot take it -- and this
world has deliberately sharp topography, since preserving carved closed-basin
floor below sea level is the point of the Orogen fork. Upstream names sharp
topography as a reason to filter even for Earth-like models. Not established
here: WHICH field's ringing reaches `swr`, and at which rungs the trap fires.

NOT SETTLED, and the row stays open for it: whether the adiabatic energy source
is the filter biasing its own diagnostic. `denergy26` carries the filter to the
first power through `sp2fl` and `denergy27` carries it squared through `dv2uv`,
so `26 - 27` cannot vanish under a nontrivial filter even for a core that
conserves exactly. An on/off arm would have separated that and cannot be run.

## The instrument reproduces on two orbits

Arm A, orbits 0-1, `close_term_energy.py`:

| quantity | this run, 2 orbits | the note's run, 10 orbits |
| --- | ---: | ---: |
| adiabatic non-conservation, `26 - 27` | **+0.3147 +/- 0.0150** | +0.3446 +/- 0.0129 |
| spectral diffusion of heat, `denergy24` | -0.0093 | -0.0099 |

Same sign, same size, on a different run at a twenty-to-one ratio of effect to
spread. So a two-orbit window resolves this quantity, and the finding is not an
artefact of the window the note happened to use.

## The experiment that replaces the refused one

Vary the filter's STRENGTH rather than its existence. `filterkappa` and
`nfilterexp` set `f(n) = exp(-kappa (n/NTRU)^gamma)`, and a weaker filter that
still stabilises gives the arms the on/off pair could not. If `26 - 27` scales
with filter strength, the diagnostic bias is the explanation; if it sits still
while the filter weakens, it is not, and the semi-implicit scheme keeps the
attribution. Both arms must be checked for the `swr` trap before their numbers
are read, since the weaker arm may simply reproduce arm B.

## The grid, and it refutes the bias hypothesis

*Measured 2026-08-23. Five arms, all from `MOST_REST.00039` of
`run_4182235e9781`, two orbits each, `NENERGY = 1`, `NLOWIO = 0`, same binary
and same superseded surface. Configs generated from `config/planet.yaml`,
differing only in the two knobs named.*

*AT gamma 8, AND ON DAMPING THAT REACHED ONE MODEL LEVEL. The model now runs
gamma 16, and these arms wrote `NDEL` and the four `TDISS*` as namelist scalars,
which in Fortran assign element one -- so levels 2 to 10 ran `plasimmod`'s
compiled defaults and, at T42, level one was read in seconds where days were
meant. Both are fixed in the writers. What the arms compare is still a
comparison at fixed damping, so the RELATIONS below stand; the absolute numbers
are for a configuration nothing runs.*

`denergy26 - denergy27`, which must be zero:

| arm | 26 - 27, W/m2 | spread | KE identity |
| --- | ---: | ---: | ---: |
| kappa 8, dt 45 | **0.3147** | 0.0150 | -0.0000 |
| kappa 2, dt 45 | 0.4020 | 0.0528 | 0.0371 |
| kappa 8, dt 22.5 | 0.1204 | 0.0114 | -0.0137 |
| kappa 2, dt 22.5 | 0.1637 | 0.0153 | 0.0064 |
| **filter off, dt 22.5** | **0.5575** | 0.0503 | 0.0312 |

The kinetic-energy identity closes in every arm against a conversion of about
2.1 W/m2, so the pairing that makes 26 - 27 readable holds throughout.

### The filter off RUNS at half the timestep

It traps at dt 45 and integrates two clean orbits at dt 22.5. **The prediction
stated before the arms ran -- that it would trap at either step, because the
failure is on the first radiation call and before time integration can amplify
anything -- is wrong.** The trap is a timestep-stability interaction, not
Gibbs ringing in the reconstruction alone.

So the filter is not required in the abstract. It is required AT A 45-MINUTE
STEP. Upstream pairs `gp|exp|sp` with 30 minutes; this configuration pairs it
with 45 and would need about 22.5 without it. **The filter buys timestep**, and
that is the finding that reaches SPAT-11: the largest scientifically acceptable
step is a function of filter strength, so the ladder has to sweep both or it
measures one arbitrary point on a trade.

### The adiabatic source is not the diagnostic biasing itself

The hypothesis was that `denergy26` carries the filter to the first power and
`denergy27` carries it squared, so `26 - 27` could not vanish under a
nontrivial filter even for a core that conserves exactly. **It is refuted, and
by the arm that makes the bias identically zero.** With the filter off, `f = 1`
everywhere, the two terms carry the same weight by construction and the bias
term is exactly nothing -- and the residual is the LARGEST of the five, 0.5575
against 0.1204 for kappa 8 at the same timestep. If bias were the cause,
removing it would drive the residual toward zero; it raises it by a factor of
4.6.

Both dependences point the same way instead:

- **It scales with the timestep.** Halving dt takes 0.3147 to 0.1204 at kappa 8
  and 0.4020 to 0.1637 at kappa 2 -- ratios of 0.38 and 0.41, between first and
  second order. A truncation error behaves like this; a spatial weighting
  artifact does not.
- **It grows as the filter weakens.** At dt 22.5 the sequence kappa 8, kappa 2,
  off is 0.1204, 0.1637, 0.5575, monotone. Less filtering leaves more
  small-scale structure for the semi-implicit step to mishandle.

So the closure note's attribution stands: the spectral core genuinely does not
conserve enthalpy plus kinetic energy, and the time scheme and the semi-implicit
conversion carry it. What is new is that the size is not a constant of the
model -- it is set by the timestep and by how much spectrum survives the filter,
and at the configuration this project runs it is 0.31 W/m2.

The unbooked-dissipation question is untouched by this: no term books what the
filter removes, and `mkdheat` still covers only Rayleigh friction and biharmonic
diffusion.

### The short-timestep arms are not contaminated by their restart

The restart was written at dt 45, and the roadmap warns that a timestep arm
needs its own compatible start because the file stores a step count and leapfrog
history. Checked rather than assumed, per orbit:

| arm | orbit 0 | orbit 1 |
| --- | ---: | ---: |
| kappa 8, dt 22.5 | 0.1284 | 0.1123 |
| kappa 8, dt 45 | 0.3253 | 0.3042 |

The short-step arm drifts between orbits by 0.016 and the control at its own
native step by 0.021. A leapfrog inconsistency would have made the first orbit
of the short-step arm stand out from its second, and it does not stand out any
more than the control does from its own.

## Two failure modes, and only one of them is the filter's

The 2026-08-23 sweep classified every SIGFPE as one event. They are two, and the
wall clock separates them cleanly.

**THE REFUSAL, at the first radiation call.** `swr_` through `radstep_` from
`gridpointd_`, in under a tenth of a minute, with no output record written at
all. This is the filter's failure: it is what the whole T42 boundary above is
made of, and it is a property of the CONFIGURATION -- the model will not take a
first step at that combination of step and filter.

**THE BLOW-UP, most of an orbit in.** `writegp_` at `outmod.f90:130`, reached
through `outgp_`, after twenty minutes and 4.3 GB of output. That line is
`zzf(:) = zf(:)`, and the declarations are the whole story: `zf` is `real`,
which is eight bytes under `-fdefault-real-8`, and `zzf` is `real(kind=4)`. So
every field this model writes passes through a DOUBLE-TO-SINGLE NARROWING, and
`-ffpe-trap=overflow` fires when a value exceeds about 3.4e38. The state had
already gone; the output writer is only where it became visible.

That makes the narrowing an accidental sanity check on the state, and a poor
one: it catches 1e39 and passes 1e10 K in silence. Worth knowing before anyone
reads a clean run as a checked one.

T127 AT dt 22.5 IS THE SECOND KIND, and it is the one cell in the matrix that
does not fit its pattern: T127 runs at dt 30 and blows up at 22.5, where a
shorter step should be safer. Three candidates, and the sweep cannot separate
them. It may be a genuine nonlinear instability, which is not obliged to be
monotone in the timestep. It may be the COLD START -- track B runs carry no
restart, so what is being integrated is a kick whose evolution depends on the
step, and that is a spin-up property rather than an equilibrium one. Or it may
be something else that moves with the step.

The first thing to establish is whether it lands in the same place twice, since
the model is deterministic and a cold start with a declared seed should
reproduce exactly. If it fails at the same step, it is a property of the
configuration; if it wanders, something in the run is not deterministic, and
archive CLIM-44 says how much that costs.

## The grid, re-measured on the damping the runs actually use

*Measured 2026-08-24 at a6d7e40c. 900-step probes, `--tau-scale 1`, all five
rungs, sixteen threads. The cost column was taken on a shared machine and is
provisional; the verdicts are not.*

The previous grid is replaced rather than corrected, because four separate
things were wrong with how it was taken and each alone invalidates a cell. The
replacement is itself superseded and is kept as a method rather than as an
answer: see "What it would take to re-take this grid".


- **The damping reached one level in ten.** `ndel` and the four `tdiss` are
  `(NLEV)` arrays and a namelist scalar assigns element one, so nine levels
  kept `readnl`'s presets -- and at T42, element one in seconds where days were
  meant.
- **The probe declared no damping at all.** It writes the `TDISS` keys only
  under `--tau-scale`, and the grid was taken without it, so every cell ran the
  model's own built-in values rather than the ones `config/planet.yaml`
  derives.
- **The probe launched under `mpiexec`.** It preferred the MPI binary and fell
  through to the threaded one when no MPI binary existed -- which is every rung
  now -- then ran it as `mpiexec -np 16`. That is sixteen independent models in
  one directory, not one model on sixteen threads.

A fourth thing was wrong with what could be measured at all: everything above
T42 died with SIGSEGV before writing a record, which reads as a refusal and is
a 16 MB process stack. T85, T127 and T170 appear here for the first time as
something other than impossible.

### What each rung will start

| rung | dt 60 | dt 45 | dt 30 | dt 22.5 | dt 15 | dt 10 |
| --- | --- | --- | --- | --- | --- | --- |
| T21 | runs | runs | runs | runs | runs | runs |
| T42 | runs | runs | runs | runs | runs | runs |
| T85 | **late** | runs | runs | runs | runs | runs |
| T127 | refuses | refuses | runs | runs | runs | runs |
| T170 | refuses | refuses | refuses | **late** | runs | runs |

`late` means the 300-step arm ran and the 900-step arm did not: a failure
inside the probe's own range, which is a different fact from a refusal and is
now reported as one.

### The filter buys exactly one step of headroom, at every rung that needs any

The same grid with the filter off:

| rung | dt 60 | dt 45 | dt 30 | dt 22.5 | dt 15 | dt 10 |
| --- | --- | --- | --- | --- | --- | --- |
| T21 | runs | runs | runs | runs | runs | runs |
| T42 | runs | runs | runs | runs | runs | runs |
| T85 | refuses | runs | runs | runs | runs | runs |
| T127 | refuses | refuses | refuses | runs | runs | runs |
| T170 | refuses | refuses | refuses | refuses | refuses | runs |

T85 goes from late to refusing at 60, T127 loses dt 30, and T170 loses dt 15.
One rung of the ladder, each time, and never more. That is a sharper statement
than the earlier "the filter buys timestep": it buys ONE step, and the amount
does not grow up the ladder even though the requirement does.

### The boundary does not follow 1/N all the way

Coarsest step each rung will start clean, at kappa 8: T85 45, T127 30, T170 15.
The middle of that is the 1/N rule -- 45 * 85/127 = 30.1 predicts T127 exactly.
T170 breaks it: 45 * 85/170 = 22.5, and 22.5 is the cell that starts and dies.
So the rule holds to T127 and over-predicts at T170 by a full step, which is
the rung it would have been used to plan.

### What this grid is, and what it is not

**It is a floor.** A cell marked `runs` has been shown to start and to survive
900 steps, which at dt 45 is about a seventh of an orbit. T42 at dt 45 is
`runs` here and dies in its sixty-seventh orbit. Nothing in this table
qualifies a step for a commissioning run, and the section below is what that
costs when the distinction is ignored.

### Cost, provisional

Minutes per orbit implied by the probe, kappa 8, where the cell runs. Taken on
a machine with other work on it, so these are an upper bound and the ratios are
sounder than the absolutes.

| rung | dt 45 | dt 30 | dt 22.5 | dt 15 |
| --- | ---: | ---: | ---: | ---: |
| T21 | 0.2 | 0.3 | 0.3 | 0.5 |
| T42 | 0.7 | 1.0 | 1.4 | 2.1 |
| T85 | 3.8 | 5.5 | 8.9 | 11.9 |
| T127 | -- | 15.6 | 20.5 | 30.2 |
| T170 | -- | -- | -- | 65.6 |

At its own coarsest clean step each rung costs 0.2, 0.7, 3.8, 15.6 and 65.6
minutes an orbit. **T170 is 94 times T42**, not the 16.5 the previous grid
reported, because it is both slower per step and pinned to a third of the step.


## T42 at dt 45 starts clean, runs for scores of orbits, and still dies

*Measured 2026-08-24, T42 on sixteen threads, `-O2`, kappa 8, `MPSTEP = 45.0`,
cold start, eighty-five orbits asked for. Two runs, differing only in whether
the damping reached every level.*

| run | damping | orbits before SIGFPE |
| --- | --- | ---: |
| `run_900548ae632e` | element one only | **47** |
| `run_3e1e116f99ee` | `10*` every level | **68** |

Correcting the damping bought twenty-one orbits and did not remove the failure.
So the per-level defect was real and expensive, and it was not the whole cause:
T42 at dt 45 is genuinely marginal, and the grid above marks it `runs` because
900 steps is a seventh of one orbit.

    #0  swr_       at radmod.f90:2530
    #1  radstep_   at radmod.f90:1523
    #2  gridpointd_ at plasim.f90:4048

The faulting instruction is `vsqrtsd`, and its operand is `273./dt(jhor,2)` --
the water-vapour amount's temperature scaling at the second level from the top.
Read out of the register at the trap, `dt` is **-12.81 K**. `dt` is
`dt(NHOR,NLEP)` in plasimmod and is temperature, not a tendency, so this is a
gridpoint at a negative absolute temperature. The instruction is reached down
the scalar `losun` branch, past a `cmpb`/`je` that skips the night lanes, so it
is a daylight cell and not a masked-lane artifact.

**Nothing was building towards it.** Per-level minimum temperature across the
run's own output, every sixth orbit and then every orbit to the end:

| orbit | 0 | 12 | 24 | 36 | 44 | 45 | 46 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| level 1 min, K | 181.1 | 185.2 | 186.4 | 184.9 | 186.7 | 184.5 | 184.3 |
| level 2 min, K | 183.3 | 190.9 | 191.0 | 190.1 | 191.8 | 190.1 | 190.9 |

Level 2's coldest cell sits at 190 K with no trend, and the four restarts
before the failure carry no non-positive temperature anywhere. So one cell went
from about 190 K to below zero inside a single orbit.

### What this changes about the grid

The grid's cells are 400-step probes and one full orbit, and this note already
says the probe cannot see a late blow-up. The correction is larger than that:
**an orbit cannot see one either.** T42 at dt 45 passes a two-orbit arm and
fails at the forty-seventh, so a cell marked as running has been shown to run
for as long as it was watched and no longer. That is not a criticism of the
measurement, which was taken for a two-orbit diagnostic and is sound for one.
It is a statement about what the grid can be used for: it qualifies a step for
a DIAGNOSTIC, and it does not qualify one for a commissioning run.

The ladder therefore runs T42 at dt 22.5, which this note measures as clean at
kappa 8 and which is also the step T170 needs -- so the rungs are compared at
one step rather than at each rung's own margin.

The 1/N law the grid fits does not predict this failure and is not contradicted
by it. `30 * 127/42` puts T42's stable step near 90 minutes, twice what failed.
Whatever ends the 47th orbit is not the linear stability boundary those probes
found, and world-td3 carries the question of what it is.

## The unbooked filter dissipation is bounded, by an identity that runs the wrong way

The filter preserves `n = 0` exactly -- `f(0) = exp(0) = 1` -- so it cannot move
the global mean of a linear field and the enthalpy budget is safe from it.
Kinetic energy is QUADRATIC in the winds, so damping high-n modes removes KE,
and nothing books that removal: `mkdheat` recomputes the wind field around
Rayleigh friction and biharmonic diffusion and returns the difference as heat,
`fluxmod` books surface and vertical-diffusion friction, and neither knows the
filter exists.

The instrument that would find it is the kinetic-energy steady-state identity,
`-denergy27 = denergy(21+22+23+25)`: in a settled run the adiabatic generation
must equal the frictional dissipation returned as heat, so an unbooked SINK
shows up as a gap. **A stronger filter removes more, so the gap must widen with
kappa if the filter is what it is missing.**

| kappa | KE identity residual, W/m2 |
| --- | ---: |
| 8 | -0.0002, -0.0114, -0.0109, -0.0113, -0.0134 |
| 4 | +0.0108, -0.0002, -0.0030 |
| 2 | +0.0074, +0.0092, +0.0052 |
| 1 | +0.1457, +0.0141, +0.0120 |
| off | +0.0429, +0.0387, +0.0325 |

**It narrows with kappa instead, and is widest with the filter OFF** -- where
the filter's dissipation is identically zero and there is nothing to book. So
the residual is not the missing term; it tracks filter WEAKNESS, which is what
the adiabatic residual does too and for the same reason: less filtering leaves
more small-scale structure for the time scheme to mishandle.

At the setting this project runs, kappa 8, the identity closes to within 0.013
W/m2 against a conversion of about 2.1 -- six tenths of a percent. So the
unbooked term is REAL IN PRINCIPLE AND BOUNDED BELOW WHAT THIS INSTRUMENT CAN
SEE at production settings. It is not worth a term in `mkdheat` on this
evidence, and the row closes on the bound rather than on an implementation.

What would reopen it: a configuration where the identity stops closing at the
adopted kappa, or a rung whose KE conversion is large enough that six tenths of
a percent is worth chasing.

## What it would take to re-take this grid

The grid above is the ladder's authority for which rung runs at which step, and
it is a boundary for a model this project no longer runs. Three things separate
it from the one that is wanted, and all three are now fixed in the writers
rather than in the reader:

1. **It declared no hyperdiffusion.** `stability_probe.py` wrote the `TDISS*`
   and `NDEL` keys only when `--tau-scale` was given and the grid was taken
   without it, so every cell ran `plasimmod`'s compiled defaults -- `ndel = 2`,
   `tdissz = 1.10 d` -- rather than the values `config/planet.yaml` derives.
   `--tau-scale` now defaults to 1.0 and `--inherited-damping` is the explicit
   arm for the old behaviour.
2. **The damping reached one model level.** A Fortran namelist scalar assigns
   element one of an `(NLEV)` array. Both writers now replicate over `NLEV`.
3. **It was taken at gamma 8.** The model runs `filter_power` 16, and T42 at dt
   90 is measured to refuse at the first and run at the second, so at least one
   ceiling in the table is wrong in the loose direction.
4. **And the probe cells cannot be attributed to a build.** The 22 entries in
   `analysis/stability_probe.json` name no executable, no sha256 and no
   template; the file carries one `generated` field that is rewritten on every
   write, so there are no per-entry timestamps; `find_template` selected its
   binary by modification time; and run directories are untracked, so git
   brackets when an entry was committed rather than what produced it. They are
   measurements of an unnamed model. The probe now checks its executable
   against `binary_manifest.json` and refuses on a mismatch, and every new entry
   records the name, sha256 and build profile -- so this is fixed forward and
   unfixable backward. world-qnue.

That fourth one is why the grid is re-taken from scratch rather than patched:
there is no cell in it that a new cell could be compared against.

And one thing about the instrument rather than the inputs: **a probe qualifies a
step against REFUSAL and against nothing else.** T42 at dt 45 passes a two-orbit
arm and dies in its forty-seventh, so refusal and endurance are two verdicts and
a re-taken grid has to carry both, with each endurance cell saying how many
orbits it actually survived.

### The command

Refusal first, because it is the half that does not need a quiet machine:

    python exoplasim/scripts/stability_probe.py --rung T170 \
        --sweep 45,30,22.5,15,10 --kappa 8,off --refusal-only --steps 600

once per rung in T21, T42, T85, T127, T170, then the cost half with the same
sweep and without `--refusal-only`. `--tau-scale` and `--gamma` default to
`config/planet.yaml`, which is the point: a cell now records what it was
measured on in its own `declared` block, the executable's name and sha beside
it, and the sha256 of the namelist and of every surface `.sra` in its bed. Pass
`--template` so the grid names the staging it was taken on rather than whichever
run directory was modified last.

### What has to be true before it is worth running

**A run directory to stage a bed from, and it is the binding constraint.** The
probe copies staging; it does not build it. It needs a namelist and the rung's
surface `.sra` family, and both come from a run directory that already exists.
With no run directory on disk there is no bed at any rung and the probe refuses
by name at `find_template` -- which is the state the tree is in after a run
purge. Getting one back is the whole staging chain: the surface family under
`exoplasim/inputs/<rung>/` has to be built, and then one arm has to be staged at
the rung, a failing one being enough. `python scripts/pipeline.py --status`
names the steps.

**Every binary rebuilt, and the manifest naming them.** The probe refuses an
executable `binary_manifest.json` does not match, so the staleness world-anl
found cannot recur -- but that guard only helps once the manifest describes the
tree. Note that a manifest is a description of source, not of a directory: a
worktree whose `vendor/exoplasim` files hash to what the manifest records still
has to build the executables locally, because build directories are deliberately
not linked between worktrees.

**A settled model source.** Rule 7 makes the grid worthless the moment the
source under `vendor/exoplasim` moves again, so it is re-taken after the model
settles and not before. A grid taken while a batch is still landing changes to
the surface builders is in the same position: the staging it was measured on
would be superseded on arrival.

**And the cost half needs the machine to itself.** It prices a step from wall
time by differencing two lengths, so a concurrent run does not add noise to it,
it adds a bias in one direction. `--refusal-only` exists so the half that does
not care can be taken anyway.
