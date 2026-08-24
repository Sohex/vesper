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

## The hyperdiffusion is not what damps this model

*Measured 2026-08-23, after the derived values were implemented and tested.*

The derived damping was applied at T42 and run two orbits from the settled
restart. **The kinetic-energy spectrum did not move**: inertial-range slope
-2.49 against the old configuration's -2.48, bite point m=28 in both, against a
change in `tau_vorticity` of 3.7x. A diagnostic insensitive to a factor of
nearly four in the quantity it is meant to judge is either broken or is watching
the wrong mechanism.

It is the wrong mechanism. The physics filter and the hyperdiffusion both damp
the top of the spectrum, and they are not close:

| n/N | hyperdiffusion, derived | physics filter | filter stronger by |
| ---: | ---: | ---: | ---: |
| 0.5 | 19579 h | 12.0 h | 1632x |
| 0.7 | 379 h | 0.81 h | 466x |
| 1.0 | 26.7 h | 0.047 h | 571x |

The filter is applied at BOTH transform directions every timestep, so a field is
multiplied by `f(n)^2` per step and the implied rate is `2 kappa (n/N)^gamma /
dt`. At the truncation that is an e-folding every 169 seconds, against the
hyperdiffusion's 26.7 hours. **Between two and three orders of magnitude, at
every scale.**

### Where the damping actually starts

Against the flow's own cascade rate at T42, the filter overtakes at
**n/N = 0.41**, and the spectrum's measured departure from its inertial range is
at **0.44**. Prediction and measurement agree to the resolution of the
diagnostic, which is what makes this an explanation rather than a coincidence.

So the model damps everything above about 0.44 of its truncation, and the
confinement criterion the derivation was built to -- damping subdominant to the
cascade until 0.6 or above -- is missed by the FILTER, not by the diffusion.

### What that overturns

**This row's premise.** The high-rung late failures were attributed to
hyperdiffusion inherited from T21 being too weak. It cannot be that: whatever
the hyperdiffusion was set to, the filter was providing several hundred times
more damping at the same scales. The inherited values were indefensible and are
now derived, which is worth having on its own, but they were never what was
holding those runs up or letting them go.

**The filter is the lever.** `filter_kappa` was chosen at 8 because it won on
stability and on the adiabatic residual, with no measurement of what it was
doing to the resolved spectrum. It is doing a great deal: reaching down to 0.41
of the truncation, where the design intent was 0.6 and above.

**And the rung-dependence still needs explaining.** The filter is scale-free in
`n/N`, so it damps the same FRACTION of the spectrum at every truncation; its
ratio to the cascade at the truncation moves only from 570 at T42 to 282 at
T170, through the timestep. Neither mechanism is strongly rung-dependent, so the
late failures at T85 and above are not explained by damping at all and the
ordinary candidate -- the advective CFL at the shorter grid spacing -- is back
in front.

## Corrected: the filter was not over-reaching, and the diagnostic was wrong

*2026-08-23. Three claims in the section above do not survive their own test, and
the cause was one defect in the instrument.*

**The diagnostic normalised by the wrong truncation.** `spectral_tail.py` took
the FFT's Nyquist as the truncation -- m=64 for a T42 run on 128 longitudes --
where the model represents nothing above m=42. Everything from m=43 to 64 sits
at 1e-15, which is roundoff. So a bite point at m=28 was reported as 0.44 of the
truncation when it is 0.67, and a filter comfortably inside its confinement
requirement was written up as damping away a third of the resolved spectrum.

What that retracts:

- **"The filter overtakes the cascade at 0.41 and the spectrum bites at 0.44,
  agreeing to the resolution of the diagnostic."** The agreement was between a
  correct calculation and a mis-normalised measurement. Against the model's
  actual truncation the bite at gamma 8 is at **0.60**, marginal against the
  0.60 floor rather than far below it.
- **"The model is damping away everything above 0.44 of its truncation."** It is
  not. Above 0.60 at gamma 8, and above 0.67 at gamma 16.
- **"kappa sets stability and gamma does not."** Measured false: T42 at dt 90
  refuses at gamma 8 and does not at gamma 16. `f(N) = exp(-kappa)` is indeed
  independent of gamma, so the grid-scale damping is unchanged, but the trap
  does not depend only on that.

### What the gamma change actually buys, and what it costs

Both arms from the settled T42 restart, two orbits, derived hyperdiffusion,
differing only in `filter_power`:

| | gamma 8 | gamma 16 |
| --- | ---: | ---: |
| spectral bite point | 0.60 of truncation | **0.67** |
| KE preserved at 0.6 of the FFT range | reference | **60x more** |
| T42 refusal at dt 90 | refuses | **runs** |
| adiabatic residual, 26-27 | -0.458 | **-1.443** |
| kinetic-energy identity | -0.000, closes | **+0.151, 7% open** |

So it is a TRADE and not a free improvement: sharper confinement and a longer
stable step, paid for with three times the adiabatic non-conservation and a
kinetic-energy identity that stops closing. The identity closing is what says
the dissipation is being booked, so 7% open is a real cost and not a cosmetic
one.

The derived hyperdiffusion carries its own share of that: at gamma 8 it moved
`26 - 27` from +0.315 under PlaSim's T42 branch to -0.458, a sign change. That
is the diffusion mattering to the energy budget while not mattering to the
spectrum, which is consistent -- the identity is sensitive to the spectral state
the diffusion shapes, not to the damping rate at the truncation alone.

### Where that leaves the filter row

Not where it was filed. The premise was that the filter reaches to 0.41 of the
truncation and throws away a third of every rung; it reaches to 0.60, which is
the floor rather than a violation of it. `gamma 16` improves confinement to 0.67
and buys timestep, and it is not obviously worth what it does to the energy
budget. That is a decision with two measured sides rather than a defect to fix,
and it should be taken against a longer window than two orbits.

## Scrutiny of the derivation, and the trade re-examined

*2026-08-23, after three corrections in one day made the rest worth re-testing.*

### The load-bearing claim holds

`tau_0 ~ 1/N` was the first rule tried, which is reason enough to test it against
others rather than accept it. Fitting `log tau_0 = c - p log N` to ECHAM5's
T31-T159 entries gives a free exponent of **p = 1.11**, and fixed exponents give
mean absolute errors of 60% at p=0.5, 26% at p=0.75, **10.6% at p=1.0**, 18% at
p=1.25 and 31% at p=1.5. The quoting precision alone -- whole hours, so plus or
minus half an hour -- is 8% of the mean. So `p = 1` fits about as well as values
quoted this coarsely can distinguish, and its neighbours are clearly worse. The
advective reading survives.

### The weak link is U, and it is weaker than the rule

`tau ~ 1/U`, and "the eddy wind" has several defensible definitions. At T42:
column eddy RMS 5.94 m/s gives 1.114 d, top-level eddy RMS 6.95 gives 0.952,
total RMS including the jet 19.77 gives 0.335. **A 3.3x spread in the definition
is a 3.3x spread in the answer.**

Worse, the constant is ASSUMED TO BE ONE. ECHAM's implied 14.7 m/s at Earth T42
is "one advective time" only if Earth's eddy wind, measured the same way as
Vesper's 5.94, is also about 14.7. Earth's column-mean transient eddy RMS is
nearer 8 to 10, which would make their constant 0.6 to 0.7 and this project's
derived `tau` correspondingly 1.4 to 1.7 times too long. **Unresolved**: settling
it needs Earth reanalysis reduced by the identical definition, which is not held
here. Until then the derived damping carries a factor of roughly 1.5 of
uncertainty on top of the definitional spread, and that is larger than most of
the differences this note has been arguing about.

### The other soft spots, named

- `alpha = 4` and `n*/N = 0.381` are a CHOICE inside a family, not a derivation.
- The confinement floor of 0.6 is a stated criterion, not a measured threshold.
- The bite-point diagnostic fits an "inertial range" over m=5-14 at T42 and
  measures a slope of about -2.0, where an enstrophy cascade would give -3. Ten
  wavenumbers at the large-scale end may not be an inertial range at all, in
  which case the bite point is measuring departure from something that is not a
  cascade.

### The trade is not an accounting artifact -- measured, not argued

If the energy cost of `gamma 16` were the filter's own unbooked dissipation
becoming visible, the kinetic-energy identity residual would EQUAL that
dissipation and would SHRINK as the filter sharpens. Computed from the measured
spectra, with the filter removing KE at `(1 - f^4)/dt` per mode:

| | filter's KE removal | KE identity residual |
| --- | ---: | ---: |
| gamma 8 | 0.0837 W/m2 | -0.0000 |
| gamma 16 | 0.0097 W/m2 | +0.1510 |

Sharpening the filter cuts its removal by 0.074 W/m2 and OPENS the identity by
0.151. Opposite directions, so the residual is not the filter's missing term.
That also re-confirms the bound this project closed the unbooked-filter row on.

**So the cost is real, and its mechanism is the time scheme.** `gamma 16`
preserves sixty times more kinetic energy near the truncation, and the
semi-implicit step does not conserve on those scales -- which is what
`water-and-energy-closure.md` measured independently and what the timestep
scaling of `26 - 27` says. The filter at `gamma 8` was not buying a good energy
budget; it was **destroying the scales on which the budget fails**, which is not
the same thing and looks identical in the diagnostics.

That reframes the question of whether speed can be bought without spending
physics. Between these two knobs it cannot: at fixed cost `gamma 8` has the
better budget and `gamma 16` the better spectrum, and `gamma 16`'s longer stable
step buys speed while making `26 - 27` worse in proportion. But the physics being
"spent" was never being preserved by the filter in the first place. The lever
that could give both is the TIME SCHEME, where the error actually lives, and
that is the semi-implicit attribution `water-and-energy-closure.md` names as
needing its own experiment.

## The Earth dependence is removable, and removing it answers the constant

*2026-08-23. The scrutiny above left the derived damping resting on a constant
assumed to be one, checkable only against Earth reanalysis this project does not
hold. That is the wrong shape for an Earth-agnostic system, and it is avoidable.*

### The two requirements are different and only one is physics

- **Absorb the cascade.** `tau <= dx/U`, which at T42 is 1.114 d. Derived from
  this planet's own radius and measured eddy wind, with no Earth in it.
- **Run at all.** The filter currently supplies an e-folding of 169 s at the
  truncation, which is **571 times more damping than the cascade requires**.

That gap is a property of the dynamical core, not of the planet, so the constant
was never "one advective time" in the first place -- it is whatever stability
demands. And what stability demands is measurable HERE, on this model, which is
what removes Earth from the derivation entirely. ECHAM's table keeps a role, but
only for the SCALING with truncation, `tau ~ 1/N`, which is a statement about
resolution rather than about a planet.

### Measured: the cascade-absorbing value is sufficient on its own

T42, **physics filter off entirely**, derived hyperdiffusion as the only damping,
dt 22.5, two orbits from the settled restart:

| | filter gamma 8, dt 45 | filter gamma 16, dt 45 | **no filter, dt 22.5** |
| --- | ---: | ---: | ---: |
| ran | yes | yes | **yes** |
| spectral bite point | 0.60 | 0.67 | **0.76** |
| pile-up (excess over the power law) | none | none | **none, 1.14x** |
| kinetic-energy identity | -0.000 | +0.151 | **+0.030** |
| adiabatic residual 26-27 | -0.458 | -1.443 | **-2.675** |

**The model runs on the derived damping alone**, and does so with the best
spectral confinement measured -- 0.76 of the truncation against the filter's
0.60 -- and no pile-up. So `C = 1`, the cascade-absorbing value, is sufficient.
It did not need calibrating against Earth, and the question of whether ECHAM's
constant is 1.0 or 0.65 does not arise: their constant is theirs, and this one
is measured here.

### What is still confounded, and what it would take

The three arms do not differ in one variable: the unfiltered one runs at dt 22.5
because dt 45 refuses without a filter. So the adiabatic residual of -2.675
cannot be read against -0.458 directly -- half the timestep should have REDUCED
it, and it is five times larger, which says the damping reduction dominates the
timestep improvement, but by how much is not separable from these three points.

A clean comparison needs the unfiltered arm at dt 45, which does not run, or the
filtered arms at dt 22.5, which do. The second is available and cheap and is the
next measurement rather than a conclusion.

What can be said without it: **the filter is not required at T42**, the derived
damping alone gives a better-confined spectrum than any filtered configuration
measured, and the kinetic-energy identity comes closer to closing without the
filter than with it sharpened. What it costs is the timestep, and whether that
is worth paying is a decision that needs the missing arm.

## The clean comparison: the trade is monotone and intrinsic to these knobs

*The confounded arms above are replaced by three that differ in ONE variable.
T42, dt 22.5 throughout, derived hyperdiffusion throughout, two orbits from the
same settled restart. Only the filter changes.*

| filter | spectral bite | adiabatic 26-27 | KE identity |
| --- | ---: | ---: | ---: |
| gamma 8 | 0.60 | -0.593 | -0.015 |
| gamma 16 | 0.67 | -1.683 | +0.019 |
| off | **0.76** | **-2.675** | +0.030 |

**Monotone in every column.** Weakening the filter improves confinement by 0.16
of the truncation and costs 2.08 W/m2 of adiabatic conservation, with no sweet
spot between. At fixed timestep the trade is intrinsic to these two knobs, and
the earlier reading -- that it might be an accounting artifact of the filter's
own unbooked dissipation -- is now excluded twice: once by the removal
calculation running the wrong way, and once by this monotone series.

The kinetic-energy identity moves the other way, -0.015 to +0.030, and stays
small throughout against a conversion of about 2 W/m2. So the dissipation is
being booked in every configuration; it is the adiabatic step that degrades.

### An unexplained reversal, recorded rather than smoothed

Under PlaSim's own T42 diffusion the adiabatic residual was POSITIVE and fell
with the timestep: +0.315 at dt 45, +0.189 at dt 30, +0.120 at dt 22.5. Under
the derived diffusion it is NEGATIVE and its magnitude GROWS as the timestep
shortens: -0.458 at dt 45 against -0.593 at dt 22.5.

A single term scaling with dt cannot do that. Two terms can -- a positive one
that shrinks with the timestep, which is what the semi-implicit truncation
should do, and a negative one that does not, which appeared when the diffusion
was weakened. Fitting the two points gives roughly +0.22 for the first at dt 45
and -0.68 for the second, but two points and a guessed functional form is not a
decomposition and it is not offered as one.

What it does say is that the adiabatic residual is not one thing, and that the
timestep scaling measured under the old diffusion does not transfer. Any
argument that leaned on "26 - 27 scales with dt" needs re-checking against the
configuration it is being applied to.

## Resolving the residual: it is paid for by the surface

*Three arms, T42, dt 22.5, derived hyperdiffusion, differing only in the filter.
Atmospheric budget, W/m2:*

| term | gamma 8 | gamma 16 | no filter | spread |
| --- | ---: | ---: | ---: | ---: |
| radiative convergence | -94.587 | -95.643 | -96.467 | 1.881 |
| sensible from the surface | 18.921 | 19.805 | 21.093 | 2.172 |
| latent from the surface | 75.850 | 76.714 | 77.326 | 1.476 |
| flux route total | +0.184 | +0.876 | +1.951 | 1.767 |
| latent asymmetry | -0.047 | +0.076 | +0.252 | 0.300 |
| **adiabatic non-conservation** | **-0.593** | **-1.683** | **-2.675** | **2.082** |
| spectral diffusion of heat | +0.269 | +0.541 | +0.271 | 0.272 |
| **closed total** | **-0.188** | **-0.190** | **-0.201** | **0.013** |

**The adiabatic term moves by 2.082 W/m2 and the closed total moves by 0.013.**
The compensation is the surface: sensible heat rises 2.17 and latent 1.48 as the
filter is removed, against 1.88 more radiative loss. So the energy the adiabatic
step destroys is replaced by drawing harder on the surface, and the model
settles at a different equilibrium rather than drifting.

That is the answer to where the residual goes, and it is not a benign one. The
column enthalpy barely moves and no budget fails to close, so nothing in the
ordinary diagnostics announces a problem -- but the surface exchange is inflated
by about two watts per square metre to feed a numerical sink. A surface flux
distorted to that degree is exactly the kind of thing downstream components read
as physics.

### What the residual is NOT, each eliminated rather than argued away

- **Not a time-level error in the diagnostic.** `adm` is copied from `sdm` at the
  head of `spectrala`, so it is the state at `t - dt`; `sdp` is a POINTER SLICE of
  `sd` rather than a "plus" array, so the update writes `sd` in place and the
  compared states are `t - dt` and `t + dt`. The step is `deltsec2` and the
  division is right.
- **Not the Robert-Asselin time filter.** `pnu` defaults to 0.0 and is absent
  from every namelist here, so `pnu21 = 1` and the filter is a no-op. The caveat
  in `water-and-energy-closure.md` about both terms carrying the time filter does
  not apply to this project's runs.
- **Not the physics filter's unbooked dissipation.** Sharpening the filter cuts
  its KE removal by 0.074 W/m2 and OPENS the identity by 0.151 -- opposite
  directions.
- **Not the kinetic-energy side.** The steady-state KE identity closes in every
  arm, -0.015 to +0.030 against a 2 W/m2 conversion, so the failure is on the
  enthalpy half of `26 - 27`.

### What is still open

The MECHANISM inside `spectrala`. What is established is that it acts on the
enthalpy side, scales with how much small-scale energy survives, and is
compensated by surface fluxes rather than by drift. What is not established is
which operation loses the energy -- the semi-implicit reference-state splitting
is the obvious candidate and has not been tested, and testing it wants a dry
adiabatic configuration where total energy must be exactly conserved and any
drift is unambiguously the numerics.

## The sink, found: the dynamical core loses energy and the filter conceals it

*Dry adiabatic runs, 2026-08-23. Radiation emptied (`nswr = nlwr = 0`, keeping
`nrad = 1` so the orbital geometry is still computed), evaporation, sensible
flux, surface stress and vertical diffusion off, and large-scale, convective,
shallow and dry-adjustment heating switched out of `rainstep`. Nothing heats or
cools the atmosphere, and there is no surface exchange left to absorb anything.*

**The atmosphere cools anyway.** T42, dt 22.5, one orbit from the settled
restart: mass-weighted column temperature falls 1.18 K, and the column enthalpy
trends at **-0.579 W/m2** -- the same number by two routes, since 1.18 K over an
orbit is 0.586 W/m2 at this column mass. Measured from the model's own
`denergy01`, not from the term decomposition that reports it.

So the sink is real. It is not the diagnostic, and with the surface removed it
has nowhere to hide.

### What it is not

| candidate | test | verdict |
| --- | --- | --- |
| mass drift | global-mean surface pressure over the orbit | +188 ppm, worth +0.02 W/m2. Not it |
| the physics filter | dry run with the filter off | sink grows to **-1.363**. The filter REDUCES it |
| time truncation | dt 22.5 against dt 11.25 | 0.85x per halving, where linear would be 0.50x |
| the kinetic-energy side | steady-state KE identity | closes at -0.035 |
| the Robert-Asselin filter | `pnu` | 0.0, a no-op |
| a diagnostic time-level error | `adm` is t-dt, `sd` after update is t+dt | the 2dt division is right |

### What it is

A loss that **grows as small-scale energy survives** and is **largely
independent of the timestep**. Both point away from the time scheme and toward
the spatial treatment of the nonlinear terms: quadratic products are dealiased
by the Gaussian grid at T42, but the primitive equations in sigma coordinates
carry non-polynomial terms through `exp(ln p_s)`, which the transform cannot
dealias.

That is a candidate class rather than a line, and identifying the line needs the
model instrumented rather than run.

### What it means for everything above

**The filter is not a damping cost. It is a mitigation of a dynamical-core
defect.** Removing it more than doubles the core's energy loss, which is why
every unfiltered configuration in this note showed a worse adiabatic residual.
"Recovering resolution" by weakening the filter is therefore not free and not
merely a trade against confinement: it exposes more of the scales the core
mishandles.

That also settles the shape of the earlier question. Speed cannot be bought by
spending physics here because the physics was already being spent -- at
0.6 W/m2 with the filter and 1.4 without it, in a configuration with no physics
in it at all.

### The resolution test was confounded, and says only what it says

Dry adiabatic, dt 22.5, one orbit, filter on:

| rung | sink | K per orbit | initial state |
| --- | ---: | ---: | --- |
| T21 | -0.041 | -0.076 | COLD |
| T42 | -0.579 | -1.182 | settled restart |
| T85 | -0.163 | -0.306 | COLD |

**Not a resolution scaling.** T21 and T85 were cold-started because no spun-up
state exists above T42, and a cold start carries almost no eddy energy for the
sink to act on. The ordering -- T42 worst, T85 next, T21 nearly clean -- tracks
how much small-scale energy each arm had, not its truncation.

It does reinforce the one dependence that has held throughout: **the sink tracks
resolved small-scale energy.** What it cannot do is separate spectral truncation
of the nonlinear terms from the alternatives, which was the point of running it.
Doing that properly needs a spun-up state at each rung, which is a commissioning
task rather than a test.

### Where instrumentation has to go, and why not elsewhere

The state changes in exactly two places per timestep. `spectrald` contributes
+0.27, matching the booked diffusion heating in `denergy24`. `26 - 27` is -0.87
against a total of -0.58, so the loss is in `spectrala`.

Inside `spectrala` the update is `stp = delt2*stt + atm` and its siblings --
arithmetic that cannot lose energy on its own. The tendencies are what carry it,
and they are formed in `gridpointa` by `mktend`, which computes the nonlinear
terms in gridpoint space and transforms them back with truncation.

**"Truncation discards energy" is not, by itself, a mechanism.** The discarded
part of a tendency lies entirely above the truncation, and the state it acts on
lies entirely below it; the two are orthogonal, so an exact projection removes
nothing the state could have received. Any energy argument that stops at "the
product exceeds the truncation" proves too much -- it would condemn every
spectral model ever written, including the ones that conserve energy to
round-off.

Two things break that orthogonality, and both live in the same routine:

  THE PAIRING IS MASS-WEIGHTED. Energy is not the plain inner product of the
  state with the tendency but the integral of one against the other weighted by
  the mass of the column, so the field the tendency is paired against is
  `ps * u`, not `u`. A product of two resolved fields is not resolved, so
  `ps * u` reaches above the truncation and is no longer orthogonal to what the
  projection throws away. Only the anomalous part of `ps` contributes, since a
  constant times a resolved field stays resolved.

  THE GRID DEALIASES PRODUCTS OF TWO FIELDS, AND `calcgp` FORMS PRODUCTS OF
  THREE. `NLON = 3*NTRU + 1` at every rung on the ladder, which is exactly the
  condition for a quadratic product to transform without aliasing. The
  sigma-coordinate tendencies are not quadratic: `zsdotp` is itself a sum of
  products, and it multiplies a vertical difference of temperature or wind to
  make `gtd`, `gud` and `gvd`. Aliasing is not a projection -- it folds content
  from above the truncation back onto resolved wavenumbers at the wrong phase --
  so nothing above protects against it.

Both are second order or higher in the flow, and every other candidate enters at
a different order. That is testable without touching the model: multiply the
state by a factor and re-measure, and each term shows itself by an exponent that
is predicted in advance. `exoplasim/scripts/scale_restart.py` builds the arms
and `exoplasim/scripts/dry_energy_order.py` fits the exponent, with the
attribution bands and the fitting window fixed in the scripts before any arm
ran. The spectral SHAPE is identical across arms, so the tail fraction near the
truncation is held while the amplitude moves -- which is the one control that
separates "how much energy is up there" from "how much flow there is".

That control was worth building because the obvious proxy failed. Across the
three dry arms the eddy kinetic energy above 0.6 of the truncation moves by 1.3x
while the sink moves by 2.4x, and the unfiltered arm carries 22 percent LESS
total eddy energy than the filtered one while sinking 2.4 times as fast. The
state's own tail is therefore not the controlling variable, and the range it
spans is too small to settle anything either way -- failure-modes class 34,
reached by an instrument whose dynamic range never covered the effect.

The semi-implicit splitting stays excluded by its own scaling: the loss per step
goes as `dt` because it is a tendency increment, so the loss RATE is
dt-independent -- 0.85 per halving against 0.50 for a first-order and 0.25 for a
second-order time error.

NOT AVAILABLE, and worth stating so it is not attempted: parking a decomposition
in unused `denergy` slots. All 28 are written somewhere, and radiation writes 9
and 10 in `gridpointd`, which runs after `spectrala` and would overwrite them
before output.

## The transform is exonerated, and the loss is in the conversion

The instrumentation went to `gridpointa` and came back with a refutation.

**The quadrature is exact.** The sink is a global-mean enthalpy loss, and the
global-mean temperature tendency comes through exactly one gridpoint field: in
`mktend` the (0,0) mode picks up neither the `fmm` term, which carries a factor
m, nor the `qmat` term, whose derivative vanishes there. So the transform can
only carry the sink if the quadrature that forms that mean is inexact. It is
not. `exoplasim/scripts/transform_exactness.py` builds each term twice from the
same T42 coefficients, once on the model's 64x128 grid and once on a 192x384
grid where the quadrature is exact to far higher degree, and every term agrees:
the quadratic products at 5e-15, `rcsq * U * dlnps/dlam` at 1e-15, and
`exp(ln ps)` at 4e-11. The two controls could have failed and did not -- a pure
quadratic, which the grid is built to dealias, and a cubic, which it is not,
but whose global mean needs only degree 126 against a 64-point rule exact to
127.

So the cubic aliasing this note reached for is real in the tendency FIELD and
absent from its global mean, which is the only part the sink can come from.

**The budget closes everywhere except one term.** In the dry adiabatic run
`spectrald` books +0.4177 to enthalpy, and that is exactly `denergy24` +
`denergy23` + `denergy25` -- spectral diffusion of heat plus the two friction
terms `mkdheat` converts. The kinetic side loses the same 0.1797 to friction.
Against the state, enthalpy trends at -0.6149, kinetic at -0.0430 and the
orographic term at -0.0006 for a total of -0.6585, and the terms reproduce that
to five percent by two independent routes.

**What is left is `denergy26 - denergy27`, and it is not a failed cancellation.**
Across `spectrala` the enthalpy falls at 1.0135 W/m2 while the kinetic energy
gains only 0.1472. The conversion between them, computed from the model's own
state rather than from either diagnostic -- omega from continuity on the refined
grid, against the specific volume -- comes to 0.45 +/- 0.22 W/m2. So neither
side matches it: the thermal side loses about 0.56 more than the conversion
accounts for and the kinetic side gains about 0.30 less, and the two shortfalls
sum to the 0.87 that goes missing. This is not one side of a large cancellation
being slightly wrong. It is a conversion whose two ends are discretised in
different variables and do not meet.

`wap` cannot be used for this and was not: the postprocessed vertical velocity
has an RMS of 2.6e-4 Pa/s where the divergence field implies 7e-2, a factor of
266, with a best-fit scale of 73 against a crude reference. That is a defect in
the diagnostic, not in the model, and it is tracked separately.

### The amplitude ladder, and why it does not run

Each candidate term enters the energy tendency at a different power of the flow
amplitude, so multiplying the state and re-measuring would name the term by an
exponent fixed in advance. `scale_restart.py` and `dry_energy_order.py` were
built for it and both are kept, because the machinery is sound and the finding
that killed the experiment is worth having.

**A scaled state is not a state this model will run.** At 0.5 and at 0.25 the
surface layer reached negative absolute temperatures within seconds and the run
took SIGFPE on `log(z/z0)`. Scaling one field group at a time located it exactly:
scaling the winds alone kills it, scaling temperature alone or humidity alone
does not. That is thermal-wind balance -- a quarter of the wind against the full
temperature gradient is not a state, and the adjustment radiates gravity waves
that blow the model up at dt 22.5.

Scaling everything together does not save it either, and finding out why
produced the one durable fact in this section: **ln(ps) is 96 percent
orographic.** Regressed on `so` across its coefficients it gives R^2 = 0.96,
because surface pressure is mostly a statement about how much atmosphere stands
above the terrain. Temperature and humidity carry no terrain signature at all,
below 0.005 at every level. So `scale_restart.py` projects the orographic
component out of ln(ps) and holds it -- but holding it while scaling the winds
is the same imbalance by another route, and the first attempt, which scaled
ln(ps) whole, lost one percent of the atmosphere's MASS. The mass drift is what
caught it.

The one arm that survived at 0.7 is therefore not evidence either: its eddy
kinetic energy came out more than double the unscaled run's, which is an
adjustment and not a scaled flow.

### What the decay says, and why it is not enough

A dry adiabatic atmosphere has no source. With the radiation emptied there is
nothing to maintain the temperature gradient the eddies feed on, so the flow
decays -- by a factor of four over one orbit -- and one run sweeps a range of
amplitudes on its own, in balance the whole way. The sink declines with it,
tracking total kinetic energy at r = 0.96 with an exponent of 1.84, which is
close to a fourth power in amplitude. But everything in a decaying flow decays
together: eddy kinetic energy falls to 0.489 of its value, the zonal mean to
0.627, the tail above 0.6 of the truncation to 0.313, and the sink to 0.483.
The measures are collinear and they do not fall in proportion, so no single
amplitude is being traced and the exponent from any one of them is an artefact
of which one was chosen. Recorded as suggestive, not as attribution.

The state's own spectral tail is separately excluded as the controlling
variable: across the dry arms it moves by 1.3x while the sink moves by 2.4x,
and the unfiltered arm carries 22 percent LESS total eddy energy than the
filtered one while sinking 2.4 times as fast.

## The term: the adiabatic conversion in the temperature equation

A control patch in `calcgp` split the global-mean temperature tendency into its
three pieces, using the model's own Gaussian quadrature. It can be split there
without a transform because `mktend`'s (0,0) coefficient takes nothing from
`gut` or `gvt` -- `fmm` carries a factor m and `qmat`'s derivative vanishes at
n=0 -- so the global-mean temperature tendency IS the global mean of `gtn`. The
vertical piece was taken by difference, so the three sum to `gtn` by
construction and the top level's separate form needed no special case.

    P1   gt*gd, the flux-form divergence term
    P2   akap*ztv2*(zvgpg - ztptb) + tkp*(zvgpg - ztpta), the conversion
    P3   the vertical advection, by difference

Late in the run, in K/s, the model's own numbers:

    nstep     P1          P2          P3          total
    597760   -3.0e-8     -1.58e-7    +1.19e-7    -6.9e-8
    600320   -3.1e-8     -9.4e-8     +6.6e-8     -5.9e-8
    600960   -4.2e-8     -8.8e-8     +9.2e-8     -3.8e-8

**P1 and P3 largely cancel, and P2 sets the total.** P2 is negative at every
print, from the first step to the last. The global-mean cooling that drives the
sink is the adiabatic conversion term, not the advection.

Alongside it, the same run printed the (0,0) coefficient of the temperature
tendency on either side of the semi-implicit correction. **They are equal to
every printed digit, and the divergence coefficient is exactly zero.** The
correction is `stt -= tau . sdt` and `sdt(1:2,:)` is forced to zero above it, so
the semi-implicit step cannot reach the global mean at all. It is exonerated,
from the model rather than by argument.

So, in the model's own arithmetic and with nothing rebuilt from output: across
`spectrala` the temperature equation removes about 1.0 W/m2 of enthalpy, led by
its conversion term, while the divergence equation delivers 0.147 W/m2 to
kinetic energy. The two ends of one conversion, and they do not meet.

### What it is not, each excluded by measurement rather than by argument

  the transform      the quadrature is exact for every term `calcgp` forms,
                     including `rcsq * U * dlnps/dlam` at 1e-15 and
                     `exp(ln ps)` at 4e-11, with a quadratic and a cubic
                     carried as controls that could have failed. Single-level,
                     so no layer thickness enters it
  its truncation     truncating the conversion costs 0.003 W/m2 net between the
                     two ends and truncating the advection 0.03, against 0.87.
                     Both are the same quantity measured with and without the
                     projection, so every error common to the two cancels and
                     only the truncation is left
  the semi-implicit
  scheme             it cannot reach the global mean, and separately the
                     matrices satisfy the adjoint identity conservation needs:
                     with `t0` isothermal, `t01s2` vanishes and the
                     construction gives tau(j,k) = tkp(k) g(k,j) dsigma(j) /
                     dsigma(k), so dsigma_k tau(j,k) = akap t0 dsigma_j g(k,j)
                     exactly
  the vertical
  discretisation     `c(a,b) = g(b,a) dsigma(a)/dsigma(b)` is built in, which
                     is the Simmons-Burridge relation the vertical scheme is
                     named for

### What is NOT established, and why the reconstruction cannot say it

Rebuilding the same three pieces from the model's gridpoint output puts P2 at
the opposite sign to the model's. That is not evidence against the model. **The
reconstruction has now produced two defects of its own and no independent
validation**, so its absolute values carry nothing:

  it rebuilt the half levels by the wrong recursion. PlaSim sets sigmah as the
  midpoint of adjacent SIGMAS (`plasim.f90:1644`); the inverse -- sigma as the
  midpoint of adjacent sigmah -- looks equally natural and is wrong by 28
  percent at the bottom layer and 16 percent at the top, where the mass is. It
  survived because the guard written against it, sigmah(NLEV) = 1, is satisfied
  by BOTH rules, as is a thickness sum of 1. A check that passes on the error it
  was written to catch is failure-modes class 17, and this one was mine.

  its one external check is inconclusive at the cadence available. The
  continuity equation ties `D + zvgpg` to the surface-pressure tendency, but the
  output's centred difference spans nearly eight days and averages away the
  field being predicted.

The conversion identity it does satisfy -- the two ends agreeing to a few
percent -- uses the same derivative operators on both sides, so it certifies
internal consistency and not correctness. The layer-thickness defect is fixed in
`dry_energy_order.py`, which now TAKES the thicknesses from `levp` rather than
rebuilding them, and says in the code why a check was not the answer.

The remedy is not a knob. Either the conversion's discretisation is corrected so
the two equations meet, or the model carries a global energy fixer of the kind
ECHAM, the IFS and CAM all carry for exactly this reason -- and a fixer is an
admission, not a fix: it restores the total without restoring where the energy
went.

### How far the reconstruction reaches, and where it stops

Two things were settled by pushing it to its limit, and both are worth keeping.

**`wap` is the right field divided by about a hundred.** Regressed level by
level against the same `omega/p` expression `calcgp` uses for `dw`, it
correlates at 0.95 to 0.98 everywhere below the top layer, and the fitted scale
times sigma is constant at 94 across every one of them. So `wap` is
`omega = (zvgpg - ztptb) * p`, uniformly scaled -- a units conversion, hPa/s
written as Pa/s, and not a physics error. `world-w51` carries that.

**`P2` reduces to one quantity.** Expanding it,

    P2 = akap * Tm * (omega/p)  +  akap * t0 * Sum_j c(j,k) D_j

and the second piece is linear in the divergence with constant coefficients, so
its global mean is zero if the divergence has no (0,0) mode. It does not: `sd`
and `sdm` carry exactly zero there in every restart checked, at every level. So
the model's global-mean temperature tendency is `akap * <Tm * omega/p>` and
nothing else.

That is also where the reconstruction stops. `<omega/p>` is a residual two parts
in a thousand of that field's own RMS, and the model's `wap` and the rebuilt
field -- correlated at 0.95 -- differ by more than that in the mean by a factor
of three. **The quantity is far below the level at which the two agree**, so
neither can settle its value, and the earlier attempt to read `P2` from outside
was measuring its own residual. Failure-modes class 34, reached honestly this
time: the instrument was pushed until it said so.

What survives is an exact constraint rather than a measurement. Since
`omega/p = D(ln p)/Dt`, mass conservation gives

    integral (omega/p) dm  =  d/dt integral (ln p) dm

and with sigma fixed that is the rate of change of the mass-weighted mean of
`ln ps`, which the run's own mass drift bounds at a few times 1e-11 per second.
The mass-weighted conversion therefore has to be the covariance term alone --
the ordinary baroclinic conversion, the 0.3 W/m2 the flow supports -- and cannot
be the 1.0 W/m2 the temperature equation removes.

Which END of the conversion is wrong is not decidable from outside the model,
and the route to it is on `world-0ov`: the kinetic side instrumented in the same
arithmetic as `P2`, in the same run.

## The budget closes: the conversion carries the sink

The mass-weighted enthalpy budget across `spectrala`, every term in the model's
own arithmetic and in the same units and by the same formula as `denergy02`, so
they are directly comparable. Averages over one orbit with the first three
prints dropped as transient, W/m2:

    P1   gt*gd                        -0.2989
    P2   the conversion               -1.0247
    P3   vertical advection           +0.6746
    ---------------------------------------
    gtn total, printed by the model   -0.6490
    denergy02, the full tendency      -0.9542
    => flux term <ps div(V T)>        -0.3051   by difference

The flux term is the piece that had never been looked at: it contributes nothing
to the global-mean TEMPERATURE tendency, because `mktend`'s (0,0) coefficient
takes nothing from `gut` or `gvt`, but it contributes to the mass-weighted
budget, which is the one the energy lives in. Recovering it closes both halves.

**THE ADVECTION CONSERVES ENTHALPY.** If it did not, the identity

    <ps (P1 + P3 + flux)>  +  <cp T dps/dt>  =  0

would fail. It comes to +0.0705 and +0.0196, summing to **+0.090 W/m2** against
a sink of 0.849 -- a tenth of it, and of the opposite sign.

**THE CONVERSION DOES NOT.** Its two ends must cancel:

    temperature equation, P2          -1.0247
    kinetic energy gained, -denergy27 +0.0860
    sum, which must be zero           -0.9387

And the two together account for the whole thing:

    -0.9387  (conversion)  +  0.0902  (advection)  =  -0.8486
    denergy26 - denergy27                          =  -0.8486

That is the sink, to the last digit, with nothing left over. **The temperature
equation's adiabatic conversion removes 1.02 W/m2 of enthalpy while the momentum
equations receive 0.09.**

### The mechanism candidate, measured

`omega/p` is `D(ln p)/Dt`, so mass conservation makes its mass-weighted integral
the rate of change of the mass-weighted mean of `ln p` -- with sigma fixed, of
`ln ps`, which this run's own 317 ppm mass drift bounds at a few times 1e-11 per
second. The model's own field, accumulated from the same `zvgpg - ztptb`
expression `calcgp` builds `dw` from, gives +9.6e-10 early and +1.6e-10 late:
**six to forty times the bound, positive at every print**, decaying with the flow
rather than fluctuating about zero.

The conversion's two ends are also split across the semi-implicit scheme in a
way that does not force them to match. The momentum equations carry
`ztv1 = T_v - t0`, the temperature ANOMALY only, and the reference part of the
pressure-gradient force is handled implicitly in `spectrala` through the
`z0 * spt` and `z0 * spm` terms of the divergence solve. The temperature
equation's reference conversion, `tkp * (zvgpg - ztpta)`, is explicit in
gridpoint space. One end in each place, and nothing between them.

### Two readings of the same decomposition, and only one is the budget

Both are true and they are not interchangeable.

  UNWEIGHTED and COMPLETE: the global-mean temperature tendency is exactly the
  global mean of `gtn`, because the flux term contributes nothing at (0,0).
  There P2 dominates and is negative.

  MASS-WEIGHTED and, until the flux term was recovered, INCOMPLETE: this is the
  enthalpy budget. There P2 is also dominant and negative, at -1.02.

An intermediate reading here put the mass-weighted P2 at small and positive.
That came from writing the reference term as `t0 * <omega/p>` instead of the
model's `t0 * (zvgpg - ztpta)`; the two differ by `<ps Sum_j c(j,k) D_j>`,
which vanishes unweighted and does not vanish weighted. The model's own
expression gives -1.02 and the table above is what stands.

## Which half of the conversion, and why that half

`P2` has two halves and the model prints them separately. Mass-weighted, in the
same units as `denergy02`, averaged over an orbit with the transient dropped:

    anomaly half      akap * ztv2 * (omega/p)      -0.0439
    reference half    tkp * (zvgpg - ztpta)        -0.9808
    P2                                             -1.0247

against a kinetic end of +0.086, so the miss is -0.939. **The anomaly half --
the physical conversion, the one the flow actually supports -- is 0.04. The
reference half carries all of it.**

That is the half whose counterpart is on the other side of the semi-implicit
split. `tkp = akap * t0` multiplies the reference temperature, and the reference
part of the momentum equations' pressure-gradient force is NOT in `calcgp` at
all: they carry `ztv1 = T_v - t0`, the anomaly only, and `R t0 grad(ln ps)` is
handled implicitly in `spectrala` through the `z0 * spt` and `z0 * spm` terms of
the divergence solve. One end explicit in gridpoint space, the other implicit in
spectral space, evaluated at different time levels and through different
operators, with nothing between them that forces them to agree.

So the defect is not the conversion in general. It is the REFERENCE conversion,
and the fix is to put its two ends on the same side of the split. That is a
change to the semi-implicit scheme -- `bm1` would need rederiving and the
stability revalidating -- which is why the correction below ships first.

## The energy fixer, which is a correction and not a fix

`model.energy_fixer`, default declared in `config/planet.yaml`, switched in the
model by `nenergyfix`, tracked as `world-mzy`. Each step it takes the imbalance
the model already reports, `denergy26 - denergy27`, and returns it as a uniform
temperature increment spread over the column heat capacity, so the energy it
gives back is proportional to the local mass.

It is an integral controller with unit gain rather than a one-shot correction:
the increment is already inside the imbalance measured after it, so subtracting
the residual leaves exactly minus the raw imbalance. It settles in one step and
then tracks. **Reading `denergy26 - denergy27` on a run with the fixer on
therefore reports the RESIDUAL, near zero, and the size of the defect is
`denergyfix` itself**, which the model prints and the manifest carries.

Why it is worth having even though it fixes nothing: without it the surface
supplies the shortfall silently, and the surface is what the rest of the project
reads. With the atmosphere no longer short those watts it stops drawing them
down, so the flux partitioning relaxes toward correct.

MEASURED, on the dry adiabatic arm, total energy trend in W/m2:

    fixer off                        -0.6595
    fixer targeting 26 - 27 only     +0.2203
    fixer targeting 26 - 27 + 24     -0.0063

with a settled correction of 0.53 to 0.60 W/m2 and a residual imbalance of
about 0.01. The target matters: aiming at the adiabatic step alone zeroes
`26 - 27` exactly, to +0.0004, and still leaves the total WORSE in the other
direction, because `denergy24` -- the temperature hyperdiffusion's heating --
is a source with no sink. The momentum diffusion's kinetic loss is booked back
as heat by `mkdheat` and cancels; damping the temperature anomaly changes the
mass-weighted mean with nothing to answer it. That first version passed its own
acceptance test perfectly while making the budget worse, and was caught only
because the criterion was the STATE trend rather than the term being aimed at.

Three further defects, each found by a run rather than by reading:

  IT DEADLOCKED. `nenergyfix` guards a block containing `mpsumbcr`, a
  collective, and namelists are read on NROOT only, so the other ranks skipped
  it and NROOT waited. `notes/audits/nlowio-collective-deadlock.md` records the
  same failure for `nlowio`, which is why that one is broadcast.

  IT WOUND UP WITHOUT BOUND. The correction was computed as a temperature
  INCREMENT and added to `stt`, a TENDENCY, so the leapfrog's `delt2` meant only
  a sixth of it landed and the controller never saw its own correction arrive.

  IT SWALLOWED A TRANSIENT. The first step out of a restart shows an imbalance
  of 252 W/m2 and unit gain takes the whole of it. Rate limited now, with a
  divergence guard that names the runaway instead of leaving it to appear as a
  blow-up in the dynamics -- and that guard is what produced the 252.

### The controller is windowed, and the reason it is was wrong

The correction is updated once per model day from the mean over that day, not
every step. The argument for it was that the PER-STEP imbalance swings by about
250 W/m2 either way -- the leapfrog's computational mode, undamped because `pnu`
is 0.0 -- so a step-by-step controller should be chasing noise, and averaging it
out should tighten the correction sharply.

**It did not.** Same configuration, same restart, second half of the orbit:

    per-step   mean +0.5594   range -0.200 to +1.830   span 2.030
    windowed   mean +0.5560   range +0.000 to +1.645   span 1.645

    total energy trend  -0.0063 per-step, -0.0079 windowed

A nineteen percent tightening, and the trend is the same within the scatter. So
the wander is NOT the computational mode: a 64-step average removes almost all
of an alternating signal and removed a fifth of this. What is left is genuine
day-to-day variation in the model's own energy imbalance, and no averaging
window makes that go away because it is the quantity being tracked.

The windowed form is kept anyway, on the grounds that survive: it targets the
systematic mean rather than a sample of it, and it changes the forcing once a
day instead of every step, which is less likely to interact with the dynamics.
Not on the grounds it was built for.

Only NROOT touches the accumulators, and only after the reduction. The obvious
cheaper design -- accumulate locally on every rank and reduce once per window --
races under the threaded build, where several threads run `spectrala` at once
against module-level accumulators. The per-step reduction is free anyway,
measured: 189 seconds an orbit with the fixer on against 191 with it off.

Why it is dangerous: it MASKS the defect it compensates. A fixer whose magnitude
nobody looks at turns a known 0.9 W/m2 into an unknown one that can grow. That is
the whole reason the applied increment is reported rather than absorbed, and the
reason `check_consistency.py` refuses a run whose manifest claims the fixer while
its namelist says otherwise.
