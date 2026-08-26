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

## The grid, on the damping and the filter the runs declare

*Measured 2026-08-26 at model source `a041a1e9`. 600-step probes,
`--tau-scale 1`, `--gamma 16`, one process on an unlimited stack. Refusal cells
on eight threads from the registry (`most_plasim_t21_l10_p8.x` `877ab48b`,
`most_plasim_t42_l10_p8.x` `482ed4bd`, `most_plasim_t85_l10_p8.x` `6722f728`),
staged from `run_9947daa3f937`, `run_d03de74d6a72` and `run_2a50670d8f2b`, each
prepared from `config/planet.yaml` on the day so the namelist keys the probe
does not overwrite are the ones the runs carry. Every cell in
`exoplasim/analysis/stability_probe.json` names its executable and sha, its
build profile, the sha256 of its namelist and of each surface `.sra` in its bed,
and a `declared` block carrying gamma, `nhdiff`, `ndel` and the four
timescales.*

Every grid before this one is replaced rather than corrected. Four independent
things were wrong with how they were taken -- the damping reached one model
level in ten, the probe declared no hyperdiffusion at all, sixteen models ran in
one directory under `mpiexec`, and everything above T42 died on a 16 MB process
stack and was read as refusing -- and a fifth, that no cell named the executable
it was measured on, is why none of them could be compared against a new one.

**A sixth was wrong in the writer and is fixed here.** `NHDIFF` is the absolute
wavenumber the hyperdiffusion starts from and it is per-rung,
`cutoff_fraction * ntru`: 8 at T21, 16 at T42, 32 at T85. The probe wrote the
four timescales and `NDEL` over every level and left `NHDIFF` alone, so a bed
staged from a T21 run damped a T42 probe from T21's wavenumber -- 38 per cent of
T21's spectrum against 19 per cent of T42's. That is the contamination the
timescales had, one key over, and it reached every cell above T21 of every grid
this note has carried.

**And the orbit was a quarter too long.** `ORBIT_HOURS` multiplied
`lib/orbit.py`'s orbital period, which is in EARTH days of 24 hours, by the
planet's 30-hour rotation. The model settles it: `run_exoplasim.py` writes 5850
steps at dt 45 and 8774 at dt 30 for one orbit, both 4387 hours, where the
restated constant implied 7312 and 10968. Every per-orbit cost the probe
reported was 25 per cent high.

### Which rungs this grid covers

T21, T42 and **T85, which had never been probed on any source**: the rung the
route's most expensive block runs at had no surface family and so no bed. The
chain that gives it one -- `boundary_conditions`, `surface_roughness`,
`surface_albedo --mode vegetated`, then one prepared run -- took under three
minutes, and all four take `--config`, so the rung moves without touching
`config/planet.yaml`. T127 and T170 carry no cell, and the ceilings
`lib/rungs.py` declares for them are carried by no measurement this tree has
produced.

### What each rung will start, at kappa 8

| rung | dt 300 | dt 225 | dt 180 | dt 150 | dt 120 | dt 90 | dt 60 | dt 45 | dt 30 | dt 22.5 | dt 15 | dt 10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T21 | refuses | refuses | refuses | refuses | runs | runs | runs | runs | runs | runs | runs | runs |
| T42 | -- | -- | -- | -- | refuses | refuses | runs | runs | runs | runs | runs | runs |
| T85 | -- | -- | -- | -- | -- | refuses | refuses | runs | runs | runs | runs | -- |

Every swept step divides the 30-hour solar day into a whole number of model
steps, so a refusal is the model's and not a calendar's: `1800 / dt` is 6, 8,
10, 12, 15, 20, 30, 40, 60, 80, 120 and 180 across the row.

**The refusal boundaries: T21 between 120 and 150, T42 between 60 and 90, T85
between 45 and 60.** T85's ceiling is 45, which is the value `lib/rungs.py`
already declares and which now has a cell behind it. T21's is 120 against a
declared 60, and the declared 60 was never a boundary -- it was the coarsest
step any earlier grid had TESTED, and this one swept above it.

### The filter buys exactly one step, at two rungs of three

| rung | dt 300 | dt 225 | dt 180 | dt 150 | dt 120 | dt 90 | dt 60 | dt 45 | dt 30 | dt 22.5 | dt 15 | dt 10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T21 | refuses | refuses | refuses | refuses | **refuses** | runs | runs | runs | runs | runs | runs | runs |
| T42 | -- | -- | -- | -- | refuses | refuses | runs | runs | runs | runs | runs | runs |
| T85 | -- | -- | -- | -- | -- | refuses | refuses | **refuses** | runs | runs | runs | -- |

T21 loses dt 120 and T85 loses dt 45: one step each, at the two rungs where the
filter changes anything at all. T42's two columns are identical, and the honest
caveat is the sweep's resolution -- the T42 bracket 60 to 90 is a factor of 1.5
wide and one step of headroom could hide inside it.

### T42 with the filter off no longer refuses at dt 45, and two things changed

The 2026-08-23 arms above trap at T42, dt 45, filter off, and that is what "the
model does not run without it" rests on. This grid runs that cell for 600 steps.
**Two things differ between the two measurements and neither is isolated here.**
The damping is now derived, per level and in days, where those arms had element
one only and, at T42, that element in seconds. And the probe starts COLD where
those arms started from `MOST_REST.00039` of a developed run, which matters more
than it looks: the refusal is on the first radiation call, so what is in the
initial state is most of the question.

The decision procedure is one arm and it needs a thing this tree does not have:
a developed T42 restart. Run the filter-off, dt 45 cell from one, and the
disagreement is the damping if it starts and the initial state if it traps.
Until then the earlier claim is not superseded and this cell is not a
contradiction of it; they are two different tests.

### What a cell in this grid is, in the units the endurance question is asked in

A 600-step probe covers **0.10 of an orbit at dt 45**, 0.02 at dt 10 and 0.68 at
dt 300; every cell carries that number as `orbits_covered`. `run_900548ae632e`
ran 46 orbits of ordinary T42 climate at dt 45 and took a SIGFPE in the 47th,
and `run_3e1e116f99ee`, the same configuration with the damping corrected to
every level, reached 68. So this grid separates two verdicts and carries only
the first:

**REFUSAL.** Whether the configuration takes its first steps at all. Measured
here, and the failure it looks for happens at the first radiation call, so 600
steps is two orders of magnitude more than it needs.

**ENDURANCE.** Whether a commissioning-length run survives. Not measurable by
any probe, and not by an orbit either. It lives in
`lib/rungs.py:COMMISSIONING_EVIDENCE`, one row per (rung, step) an actual run
has reached. A pair with no row there has no endurance evidence, which is a
different thing from passing, and three of the six rows the escalation route
needs are missing. `exoplasim/notes/route-step-criteria.md` is where the route
is chosen from the two together.

**T42 at dt 45 is marked `runs` here and is the pair with both recorded
blow-ups.** Nothing in this grid qualifies a step for a commissioning run.

### Cost: the half that needs the machine, and this host did not have it

*Attempted 2026-08-26 under the host lock, one cell per rung, two passes, p16.*

| rung | dt | s/step | implied min/orbit |
| --- | ---: | ---: | ---: |
| T21 | 45 | 0.1480 | 14.4 |
| T42 | 30 | 0.1561 | 22.8 |
| T85 | 22.5 | 0.1962 | 38.2 |

**These numbers are refused, and the check that refuses them is inside the table
rather than beside it.** Cost per step is the model's own work, so it must grow
with the truncation: T85 solves sixteen times T21's gridpoints and its Legendre
transform grows faster still. The three cells span a factor of **1.33** where
the model's work spans something between one and two orders of magnitude. An
instrument that cannot separate T21 from T85 is not measuring the model.

**C-ROUTE-5, the criterion fixed before these cells were taken, PASSES them.**
It asked that the two passes agree to better than the saving being argued, and
they agree to 0.5 per cent at T21. That is not a rescue, it is the lesson: the
contention on this host was STEADY rather than bursty, and a steady bias
reproduces perfectly. Repeatability is not validity, and a scatter check between
two passes of the same instrument cannot see an error both passes share. The
check that catches these numbers is the one inside a single pass -- three rungs
that must differ and do not -- and it is the check C-ROUTE-5 should have been.

The mechanism is the host and it is legible. The lock is a mutex across AGENTS
and does not stop a process that outlives its holder's phase, and four other
sixteen-thread-equivalent models were resident throughout at a load average
between 42 and 50 on 32 cores. Worse, every one of them was PINNED: the probe
exports `OMP_PLACES=cores` with `OMP_PROC_BIND=close`, which binds from core
zero upward, so concurrent runs stack on the same low-numbered cores while the
upper half of the machine idles. `taskset` on three resident masters returned
cores 0, 8 and 16. Pinning without a disjoint placement is worse under sharing
than not pinning at all.

So the per-orbit costs the route wants are NOT MEASURED. What replaces them is
in `exoplasim/notes/route-step-criteria.md`: the route's saving is expressed in
MODEL STEPS, which is exact arithmetic and needs no clock at all. Steps cannot
give the weight of a T85 step against a T21 step, which is precisely the
quantity this table failed to deliver, so the rungs are reported separately and
never summed.

**What it would take**: the same six probes on a host with nothing else on it.
Under the load actually seen they cost about twenty minutes; on a quiet machine
the same six are a couple of minutes, and the earlier grid's T21 figure of about
0.2 minutes an orbit is the scale to expect rather than 14.4.

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

The grid's cells are 600-step probes, a tenth of an orbit at dt 45, and this note already
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

## What it takes to probe a rung, and what T127 and T170 still need

T21, T42 and T85 are measured on the current source and the declared damping.
T127 and T170 are not, and the obstacle is staging rather than the probe.

**A bed at the rung, and it is the binding constraint.** The probe copies
staging and does not build it: it needs a namelist and the rung's surface `.sra`
family, and both come from a run directory. Getting one is the surface chain --
`boundary_conditions`, `surface_roughness`, `surface_albedo --mode vegetated` --
run with `model.resolution` at that rung, and then one prepared run to stage
from. At T85 that whole chain took under three minutes, so it is not expensive;
it simply had never been done above T42, which is why the rung the route's most
expensive block runs at had never been probed at all.

The three builders and `run_exoplasim.py` all take `--config`, so the rung can
be moved without touching `config/planet.yaml`. The parsed configuration travels
into each artifact's `source_config` and reads `T85` there, which is what makes
the products attributable.

**The order, and the reason for it.** Refusal first with `--refusal-only`: it
carries no cost field at all, so it can be taken on a contended host, and it is
the half the route's necessary condition needs. The cost half differences two
wall times on one bed, which cancels startup but not contention, so it wants the
machine.

    python exoplasim/scripts/stability_probe.py --rung T127 \
        --sweep 60,45,30,22.5,15,10 --kappa 8,off --refusal-only \
        --steps 600 --template <a run directory at that rung>

`--tau-scale`, `--gamma` and `NHDIFF` all default to `config/planet.yaml`, and
every cell records what it was measured on: the executable's name, sha256 and
build profile, the sha256 of the namelist and of each surface `.sra` in the bed,
and a `declared` block carrying gamma, `nhdiff`, `ndel` and the four timescales.

**Sweep only steps that divide the solar day.** `1800 / dt` must be a whole
number of model steps or the refusal being measured is partly the calendar's.

**And the sweep has to reach a refusing cell.** A column of `runs` all the way
down is not a boundary, it is a range that never got there --
`docs/src/practice/failure-modes.md` class 34. The first pass of this grid swept
90 down to 10 at T21 and found nothing refusing at all; the boundary is between
120 and 150.

**Rule 7 sets when, not whether.** The grid rests on the model source, so a
change under `vendor/exoplasim` makes every cell worthless on arrival rather
than stale. Take it after the source settles, and re-take it after the next
change that reaches the dynamical core.

**What no probe at any rung can settle** is the endurance question, and that is
where the route's steps actually come from. `exoplasim/notes/route-step-criteria.md`
carries the rule and the evidence.

