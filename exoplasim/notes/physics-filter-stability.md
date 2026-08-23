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
