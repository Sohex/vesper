# What "tuned for coarse resolution" means in this model, concretely

*Worldbuilding frame: a source reading of the climate model the Vesper project
runs, checked against the model's own output. Nothing here is about the
simulated planet. Read 2026-08-23 at ExoPlaSim 3.4.2.*

The guidance around this model says, in effect, use it at coarse resolution and
be careful above that, without saying what moves. Four places in the source
branch on resolution, and they are not one phenomenon: **one is live and
silently wrong above T42, one is unreachable code at every resolution, and two
are narrow.**

## 1. Hyperdiffusion: LIVE, and every rung above T42 gets T21's

`plasim.f90:1388` is the whole of the resolution dependence:

    if(NTRU==42) then
     nhdiff=16
     ndel(:)=4
     tdissq(:)=0.1  * day_24hr
     tdisst(:)=0.76 * day_24hr
     tdissz(:)=0.3  * day_24hr
     tdissd(:)=0.06 * day_24hr
    endif

There is no branch for T85, T127 or T170, so they fall through to the module
defaults at `plasimmod.f90:383,816-819`, which are T21's:

| | T21 defaults, used at T85/T127/T170 | T42 branch |
| --- | --- | --- |
| `ndel`, operator order | 2, so grad^4 | 4, so grad^8 |
| `nhdiff`, critical wavenumber | 15 | 16 |
| `tdisst`, temperature | 5.60 d | 0.76 d |
| `tdissz`, vorticity | 1.10 d | 0.30 d |
| `tdissd`, divergence | 0.20 d | 0.06 d |
| `tdissq`, humidity | 0.10 d | 0.10 d |

**So the finest grids run the coarsest grid's damping**: a lower-order operator
and timescales three to seven times longer, on a mesh with four to sixteen times
the gridpoints. Confirmed from the model's own namelist echo rather than from
the source -- `NHDIFF=16` in a T42 run and `NHDIFF=15` in a T170 one, with
`NDEL= 10*4` against `NDEL= 10*2`.

The units are NOT a problem, which is worth recording because they look like
one. The T42 branch multiplies by `day_24hr` and the module defaults do not, so
the echo shows 65664 at T42 against 5.6 at T170. `dayseccheck`
(`plasim.f90:1588`) catches exactly that: below one timestep it assumes days and
converts, logging as it goes, and the T170 diag carries
`assuming [days] - converting to [sec]` for all three. The values differ; the
units do not.

WHAT THIS PREDICTS, and it is testable rather than asserted: under-damped small
scales at high resolution should show up as a shorter stable timestep and as
late blow-ups rather than clean refusals. Both are what the ladder measures --
T127 fails after twenty minutes at dt 22.5 and T170 after 3.8 at dt 30 -- so
the hyperdiffusion inherited from T21 is a live candidate for the ladder's
stability ceiling and has not been separated from the CFL limit.

## 2. Radiation: UNREACHABLE at every resolution, including T42

`radmod.f90:909-980` sets shortwave transmissivities and the water-vapour
continuum per resolution -- `tswr1`, `tswr2`, `tswr3`, `th2oc` -- through a
variable named `jtune`, with branches for T21/T1, T31 and T42. Every branch is
guarded:

    if(NDCYCLE==1) then
     jtune=0
    else
     ... tswr1=0.089 ... jtune=1

`ndcycle` defaults to 1 at `radmod.f90:186`, and **the block runs at line 909
while `read(11,radmod_nl)` is at line 988** -- seventy-nine lines later. So
`ndcycle` is always its compiled default when the test is made, whatever the
namelist says, and `jtune` is always 0.

The model announces it: *"No radiation setup for this resolution (NTRU,NLEV) /
using default setup. You may need to tune the radiation"*. **Every run this
project has ever made prints it** -- 88 run directories across T21, T42, T85,
T127 and T170, including all sixty-six T42 runs, where a branch exists and
cannot fire.

So the radiation half of "tuned for coarse resolution" does not apply to this
project at all, and the warning about it has always been correct and always been
noise, since no configuration can clear it.

## 3. Two narrow ones

`rainmod.f90:85` clears `nshallow` at T21 with 5 levels, and `rainmod.f90:88`
sets `gamma=0.007` at T42 with 10 levels -- which this project's configuration
does hit. `plasim.f90:1381` sets Rayleigh friction timescales when `NLEV==10`,
which is resolution-independent and applies here.

## What is NOT free-floating

The radiation's absorptances are fits -- Lacis-Hansen shortwave and Sasamori
longwave -- and this project has re-derived seven per-star corrections onto
them. Those are fitted to PHYSICS rather than to a resolution, and they are not
part of this. The distinction worth keeping: a coefficient fitted to a
line-by-line calculation is a model of an absorber, and a coefficient set by a
branch on `NTRU` is a knob compensating for a grid.
