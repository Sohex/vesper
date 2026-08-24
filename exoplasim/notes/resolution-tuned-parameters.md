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


## The derivation: what the damping has to absorb

*Added 2026-08-23, with the three sources in `references/`.*

**Neither source supplies a derivation, and both say so.** ECHAM5's model
description is explicit that horizontal diffusion "does not involve a physical
model of subgrid-scale processes, but rather a numerically convenient form of
scale selective diffusion with coefficients determined empirically to ensure a
realistic behavior of the resolved scales" (Report 349 section 4). Laursen and
Eliasen (1989), whose scheme PlaSim implements, chose theirs by running a series
of integrations at different diffusion strengths and comparing the resulting
kinetic-energy spectra against OBSERVED January means. So the numbers in the
PlaSim Reference Manual are two points from someone else's tuning exercise
against Earth observations, and nothing in the literature makes them anchors.

What the literature does supply is a form and, unexpectedly, a corroboration.

### ECHAM5's empirical table implies one physical rule

Report 349 Table 4.1 gives the e-folding damping time of the highest resolvable
wavenumber -- the same quantity PlaSim's `tau` is -- across seven truncations.
Read as an advective time at the smallest resolved scale, `tau_0 = dx / U` with
`dx = pi a / N`, it implies a wind:

| truncation | tau_0, hours | implied U, m/s |
| --- | ---: | ---: |
| T21 | 6 | 44.1 |
| T31 | 12 | 14.9 |
| T42 | 9 | 14.7 |
| T63 | 7 | 12.6 |
| T85 | 5 | 13.1 |
| T106 | 3 | 17.5 |
| T159 | 2 | 17.5 |

**From T31 to T159 the implied wind is 12.6 to 17.5 m/s** -- constant to within
the precision of values quoted only in whole hours, across a factor of five in
resolution, and squarely a real upper-tropospheric eddy speed for Earth. Their
empirical tuning is therefore consistent with a single physical statement:
**damp the smallest resolved scale on its own advective timescale.** T21 is the
outlier at 44 m/s, and its per-level order structure differs from every other
column in the table too.

That is a rule rather than a number, it contains the planet through `a` and the
flow through `U`, and it was not fitted to the values it reproduces.

### Applied to this planet

`U` is the EDDY wind -- the departure from the zonal mean, which is what
cascades to the truncation; the zonal-mean jet does not. Measured from
`run_4182235e9781` orbit 39, area-weighted over the column: **5.94 m/s**, against
14.7 implied for Earth. Vesper's eddies are slower and its grid boxes larger, so
the same rule asks for LESS damping at a given truncation, not more.

**The rule's `tau` belongs to VORTICITY.** The argument is about the enstrophy
cascade, and enstrophy is a vorticity quantity. That assignment is checkable on
Earth: the rule gives 9.0 h at T42 against PlaSim's Earth-tuned `tau_xi` of
7.2 h, agreeing to 25%. (An earlier version of this note compared the rule
against `tau_T` instead and reported errors of 10 to 20 times; that was the
wrong variable and the numbers below replace it.)

The other three timescales keep T42's ratios to vorticity -- 0.200 for
divergence, 2.533 for temperature, 0.333 for humidity. Those are INHERITED, not
derived: they are taken from the one rung whose base the rule corroborates, and
the T21 defaults disagree with them (5.09 for temperature against 2.533), which
is one more sign that the two documented columns are independent tunings rather
than a system.

| rung | tau_xi = dx/U | tau_D | tau_T | tau_q | what the model uses, vorticity |
| --- | ---: | ---: | ---: | ---: | ---: |
| T21 | 2.229 d | 0.446 | 5.646 | 0.743 | 0.49x -- **2.0x too strong** |
| T42 | 1.114 d | 0.223 | 2.823 | 0.371 | 0.27x -- **3.7x too strong** |
| T85 | 0.551 d | 0.110 | 1.395 | 0.184 | 2.00x -- 2x too weak |
| T127 | 0.368 d | 0.074 | 0.934 | 0.123 | 2.99x -- 3x too weak |
| T170 | 0.275 d | 0.055 | 0.697 | 0.092 | 4.00x -- 4x too weak |

So the error is not one-signed. The fixed defaults are too STRONG at the coarse
end and too WEAK at the fine end, crossing over between T42 and T85, which is
what a resolution-independent number does against a rule that scales as `1/N`.

### The shape is the larger defect, not the strength

Damping rate as a fraction of the LOCAL cascade rate, which is what says whether
the damping is confined to the scales it should be:

| setting | at n = N/2 | at 0.7N | at 0.9N |
| --- | ---: | ---: | ---: |
| T21 as documented, alpha 2, n*/N 0.714 | 0.000 | 0.000 | 0.469 |
| T42 as documented, alpha 4, n*/N 0.381 | 0.003 | 0.101 | 0.549 |
| **T21's settings at T170**, alpha 2, n*/N 0.088 | **0.408** | 0.643 | 0.881 |

**At T170 the inherited settings damp at 41% of the local cascade rate at HALF
the truncation**, against T42's 0.3%. That is the `n*` defect: an absolute cutoff
of 15 is 71% of the truncation at T21 and 8.8% of it at T170, so what is a
confined grid-scale filter at T21 becomes a broad drag across the resolved
spectrum at T170. Under-damped where the cascade arrives and over-damped
everywhere it should not reach, at once.

**Adopted: `alpha = 4` and `n*/N = 0.381` at every rung** -- T42's shape, applied
as a FRACTION so the confinement is the same at every truncation. This is a
choice within a family, not a unique derivation: the requirement is that damping
stay subdominant to the local cascade across the resolved range and become
comparable only at the truncation, and several `(alpha, n*/N)` pairs satisfy it.
This pair is the one attached to the corroborated base, and the spectral
diagnostic is what judges it.

### What the rule does not yet settle

- **One timescale against four.** The rule gives a single `tau_0`; PlaSim
  carries separate values for divergence, vorticity, temperature and humidity
  spanning a factor of twelve at T42. ECHAM5 quotes one `tau_0` for all three of
  its diffused variables. The per-variable ratios are a further undocumented
  choice and are not derived here.
- **`alpha` and `n*`.** The rule fixes the damping AT the truncation and says
  nothing about how far down the spectrum it should reach.
- **Whether `U` is resolution-dependent.** It was measured at T42; a finer grid
  resolves more eddy kinetic energy, so `U` may rise with truncation and the rule
  may need one iteration. Most eddy energy sits at large scales, so the effect is
  expected to be small, and it is a prediction to check rather than an assumption.
