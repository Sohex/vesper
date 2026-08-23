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

## T170, and the filter requirement gets STRONGER up the ladder

*Measured 2026-08-23. One full orbit at dt 22.5 through `run_exoplasim.py`; the
rest of the grid by `stability_probe.py`, 400 steps with output off.*

| dt | kappa 8 | kappa off |
| ---: | --- | --- |
| 45 | refuses | refuses |
| 30 | starts, blows up after 3.8 min | refuses |
| 22.5 | **completes an orbit, 51 min** | refuses |
| 15 | starts | refuses |
| 10 | starts | starts |

**The unfiltered model needs dt 10 at T170 where it needed 30 at T42.** That is
the resolution dependence upstream states as the reason for filtering at all --
a finer grid resolves sharper gradients off the same orography, so the ringing
worsens while a scale-free filter keeps cutting the same fraction -- and it is
now measured rather than quoted. It also settles the shape of the trade: the
filter is not a fixed overhead to be tuned once, it is buying more timestep the
further up the ladder the model goes.

The arithmetic pick landed. T127 runs at dt 30, the stable step scales roughly
as 1/N, and 30 * 127/170 = 22.4 predicted 22.5 -- which is exactly the coarsest
step that carries T170 through a whole orbit.

### The ladder priced, normalised to dt 45

| rung | s/orbit | x T42 | measured at |
| --- | ---: | ---: | --- |
| T21 | 37 | 0.40 | six timesteps |
| T42 | 93 | 1.00 | seven timesteps |
| T85 | 296 | 3.18 | four timesteps |
| T127 | 760 | 8.17 | dt 30 |
| T170 | 1530 | 16.5 | dt 22.5, one full orbit |

T170's declared bracket was 18.7 and it measures 16.5, so the brackets ran 13 to
30 percent high all the way up. **What a rung actually costs, at the step it can
actually run**, is the second column times the step ratio: T170 is not 16.5
times T42, it is 16.5 times T42 AND needs half the step, so an orbit costs 33
times what T42's does. That factor is the one SPAT-11 exists to state.

### What the probe is worth, and where it is not

Against the full-orbit measurement at dt 22.5 the probe reads 20% high --
60.8 minutes against 51.0 -- and the offset is consistent across the other
timesteps, since correcting by it reproduces the linear-in-steps law from the
one ground truth to within a percent. So the probe is a BRACKET for cost and a
verdict for refusal, and the two should not be quoted the same way.

Two lengths are differenced rather than one divided, because a single 200-step
probe read 69.2 minutes an orbit where the truth was 51: a bed shorter than its
own startup measuring its startup, which is class 34 and was caught here by
having a ground truth to check against.

**The probe cannot see a late blow-up**, and T170 at dt 30 is the demonstration:
it passes 400 steps and dies after 3.8 minutes of integration. A cell that
"runs" under the probe has only been shown not to REFUSE.

### Every ladder cell is ONE DRAW, and that is now fixed forward

The ladder was measured on COLD starts, and `initrandom` (`plasim.f90:2068`)
takes the namelist `SEED` when `seed(1)` is non-zero and the SYSTEM CLOCK
otherwise. Nothing wrote `SEED`. So each cell above is one draw of an initial
kick rather than a verdict that re-running would confirm.

The instant refusals are almost certainly robust anyway: they are monotone
across many cells and both filter settings, and a marginal result decided by a
random kick would not produce a clean staircase. The LATE failures are the
exposed ones -- T127 at dt 22.5 and T170 at dt 30 are exactly the shape of
result an initial condition can move.

They are kept as one draw deliberately. Both are loud, immediate to recognise
and understood well enough to act on if they turn up in other work, and a
re-measurement would cost orbits to confirm a boundary nothing is planning to
sit on. `model.cold_start_seed` is declared from here, so the next cold run is
reproducible even though these were not.

### What the T127 dt 22.5 blow-up actually wrote

Not a prognostic field. `outmod.f90:611` writes `aroff` -- ACCUMULATED RUNOFF,
code 160 -- on the line after `aroff(:) = aroff(:)/real(naccuout)`. So the first
quantity to exceed single precision is a diagnostic accumulator with a step
count in its denominator, which is not the same claim as the state having gone.
Whether the state was already bad is UNDETERMINED, and this note said otherwise
before the call site was read.

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
