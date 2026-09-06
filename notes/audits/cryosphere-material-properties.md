# What the modelled sea ice and its snow are made of

**Derived:** 2026-08-25 and 2026-08-26, by `analysis/ice_properties.py` from
IAPWS R10-06(2009) and from the papers named at the foot of this note. No World
Orogen generation and no ExoPlaSim run: every number here is either a laboratory
standard evaluated, or a published relation evaluated at a declared density or
temperature, or a published relation integrated OFFLINE over an existing
climatology. The two places that read a climatology are the temperature span of
the modelled sea ice and the snow compaction pricing, and both are labelled
where they appear: a bootstrap's numbers are not the baseline.

This is worldbuilding. Vesper is an invented super-Earth around a mid-K dwarf,
and every quantity below is a material property of that planet's simulated sea
ice or of the simulated snow that lies on it. IAPWS-06 and the snow papers are
laboratory measurements of MATERIALS, and a material does not change when the
star does, which is the whole of why they transfer here.

WORLD-04OK and WORLD-A9S5 asked, of six compile-time constants in
`icemod.f90` and of one namelist key in `landmod.f90`, what each one's source
is and whether it should move for this world; WORLD-FG8W asked the same two
questions of the pair that describes the GLACIAL ice `glaciermod` grows. This
note is the answer, one row at a time.

## The finding, in one line each

| Constant | Was | Verdict |
| --- | --- | --- |
| `CRHOI`, density of the modelled sea ice | compile-time | **unchanged, DECLARED with its bound.** Not derivable in this model: sea-ice density is a function of brine volume and the model carries no ice salinity |
| `CPI`, specific heat of the modelled sea ice | compile-time | **unchanged, DECLARED with its bound.** Same reason; the effective heat capacity of sea ice carries the latent heat of its own internal melting |
| `CKAPI`, conductivity of the modelled sea ice | compile-time | **unchanged, DECLARED, and the bound is now a number.** Same reason; the conductivity of sea ice is pure ice's LESS a brine term. Yen (1981) supplies both the pure-ice arm and the band the brine term opens, and the declared value sits inside |
| `CLFSN`, melting enthalpy of the modelled snow | compile-time | **changed, DERIVED.** Snow has no brine, so IAPWS-06 gives it exactly. The compiled value was 0.08 per cent high |
| `CPSN`, specific heat of the modelled snow | compile-time | **deleted.** Declared and read nowhere in the vendored tree |
| `CKAPSN`, conductivity of the modelled snow | compile-time | **changed, HANDED.** A second statement of `landmod`'s `snowdiff`, which now follows the snow density |
| `snowdiff`, conductivity of the modelled snow on soil | `landmod_nl` key | **changed, DERIVED from `rhosnow`.** A namelist density beside a fixed conductivity is a broken relation |
| `Ksnow`, conductivity of the modelled snow in the vegetation model | Sturm et al. (1997), compiled | **changed, DERIVED from the same relation.** Below both arms of the bracket at every density the two components span; the relation is declared once in `lib/snow.py` and restated in both models under a check |
| `Csnow`, volumetric heat capacity of the modelled snow in the vegetation model | ice Ih's specific heat times a local `ice_density` of 917, compiled | **changed, DERIVED from `snowdens`.** A volumetric heat capacity is a density times a specific heat and the density is the snowpack's. Pinned at solid ice's, the modelled pack's thermal mass per unit water equivalent was 917/`snowdens` too high |
| `sicecap`, heat capacity of the modelled glacial ice | `landmod_nl` key | **changed, DERIVED from `rhoglac`.** It factorised as one thousand times ice's specific heat, which is LIQUID WATER's density; a bracket on the ice density moved the orography and left the ice's thermal mass behind |
| `sicediff`, conductivity of the modelled glacial ice | `landmod_nl` key | **changed, DERIVED from `rhoglac`.** Pure ice from Yen's Eq. (33), reduced by his Eq. (37) for the air the density implies. Glacial ice is bubbly and the reduction is the whole difference from pure ice |
| `rhosnow`, density of the modelled snow | `landmod_nl` key, no compaction | **unchanged, and now PRICED.** A prognostic density is worth a factor of six on the pack's conductive resistance and the gravity term inside it is worth 14 per cent, which is smaller than the vapour-kinetics bracket the conductivity already carries. On the modelled surface albedo it is worth exactly zero, because the snow-covered fraction is taken in water equivalent |

## The line between what is derivable here and what is not

The model carries its sea ice as **one number per cell: a thickness.** No
salinity, no brine volume, no porosity, no vertical structure. Every material
property of sea ice that differs from pure ice's differs *through the brine
volume*, and brine volume is a function of the ice's own salinity and
temperature. So the model cannot derive any of them, and a number that looks
derived would be a number derived from a salinity nobody declared.

`icemod.f90` already said exactly this about `CLFI`, the heat of fusion of sea
ice, which is DECLARED for that reason and reaches the model from
`config/planet.yaml`. **The same sentence covers `CRHOI`, `CPI` and `CKAPI`, and
recognising that is most of what these two rows are.** They are declarations
about what this world's sea ice is, and what `analysis/ice_properties.py`
supplies for them is the pure ice Ih value as the bound each sits against, so
that a declared number is at least placed against a computed one rather than
floating.

The snow is different in kind. Snow is ice Ih plus air, with no brine in it, and
the air carries a negligible share of its mass. So its melting enthalpy and its
specific heat ARE ice Ih's, and its conductivity is not a property of ice at all
but of a porous arrangement of it, which is a function of the density.

## The conductivity's bound is Yen's, not IAPWS-06's

Two of the three sea-ice declarations get their bound from IAPWS-06 and the
third cannot: the release is a Gibbs function, and a Gibbs function carries
density, specific heat and compressibility but no TRANSPORT property. Nothing in
it says how fast heat moves. So `CKAPI`'s bound had no arithmetic behind it at
all, only the sentence "pure ice's, less a brine term this model cannot
compute", and Yen (1981) is what turns that sentence into a number.

Yen fits `lambda = a exp(bT)` to the polycrystalline-ice measurements of seven
groups and reports three arms, because the data have a gap between 150 and 195
K. Two of them reach the temperatures the modelled sea ice occupies: the
whole-range arm he recommends for practical use, and the arm fitted only above
the gap. `analysis/ice_properties.py` carries both rather than averaging them,
because their disagreement at the melting point -- a few per cent, growing to
nothing at the cold end -- is the honest width of "pure ice's conductivity" and
is small against what is being bounded.

**The declared conductivity of the modelled sea ice sits below pure ice's over
the whole range, which is the side brine puts it on.** Yen's equation (71) is
explicit about the direction: sea ice is bubbly pure ice with the brine's much
lower conductivity mixed in, SUBTRACTED in proportion to the salinity over the
temperature, so sea ice conducts less than the fresh ice it is made of. His
Figure 23 computes that band over the warm end, and it is wide: across the
salinities and densities he plots, the sea-ice conductivity spans nearly a
factor of two, with the zero-salinity curves at the top.

**The declared value sits at the top of that band, against the zero-salinity
curves.** That is the correct place for it, and it is the answer to what the
declaration means: this model carries no salinity and no brine volume, so its
sea ice IS Yen's zero-salinity limit, and the declaration is that limit rather
than a value for salty ice. What Yen's figure adds is how far below it real
brine-bearing ice would sit, which is the size of the thing the model is not
carrying -- and it is much larger than the gap between the declared value and
pure ice's.

The temperature range this is checked over is the one
`analysis/ice_properties.py` already declares for the ice Ih profile, and it
contains the range the modelled sea ice occupies: in the bootstrap
climatology's monthly means, measured 2026-08-25, the surface over cells more
than half ice-covered spans 247 to 272.5 K. Those are a bootstrap's numbers and
monthly means clip the extremes at both ends, so they are used to confirm that
the declared range brackets the model's, not as the model's range.

## Gravity, checked rather than assumed

This project has found several gravity questions that cancel, and this is
another one. Surface gravity reaches a material property of ice through
overburden pressure and through nothing else. `analysis/ice_properties.py`
evaluates the isothermal compressibility of ice Ih against the overburden under
a column far thicker than this model carries, at this world's gravity and at
Earth's, and the density response is **parts per million in both cases.** The
difference between the two is a part in a hundred thousand of the answer. It
cancels.

Where the surface gravity does reach the modelled cryosphere is the SNOW's
density, which is set by compaction. That is GRAV-8's row and is not settled
here. What WORLD-A9S5 changes is that the route is now open: the snow
conductivity follows the density, so when GRAV-8 moves `rhosnow` the heat flux
through the pack moves with it instead of staying behind.

Atmospheric pressure would have been a second route, because the latent-heat
term in a snow pack's effective conductivity scales with the diffusion
coefficient of water vapour in air, which goes as the inverse of pressure. This
world's declared partial pressures sum to one bar to six figures, so that route
is closed too, and the snow relation below transfers at face value.

## The snow conductivity: which relation, and why not Sturm

WORLD-A9S5 named Sturm et al. (1997) and Yen (1981) as the two candidates and
observed that Sturm's quadratic at the shipped snow density returns about half
the value the model carried. That observation is correct and its conclusion --
that the two cannot both be right -- turns out to be the wrong conclusion, which
is why this row was worth checking rather than implementing.

**Sturm's regression is a needle-probe measurement set**, 488 samples at a mean
temperature near -14.6 C, and the paper itself states that its predicted values
carry an uncertainty of 0.1 W/m/K at 95 per cent confidence, which at the
shipped density is comparable to the value being predicted. That alone is the
project's own instrument-versus-effect test, and it very nearly fails it.

**The needle probe reads low on snow, and the two papers that show it disagree
about why.** Riche and Schneebeli (2013) applied three independent methods -- a
long-heating needle probe, a guarded heat flux plate, and direct numerical
simulation -- to IDENTICAL samples, and concluded that the numerical simulation
is the most reliable of the three, with a horizontally inserted needle probe
wrong by up to a quarter either way through the pack's anisotropy.

Calonne et al. (2011) reports the same direction from the other side and REFUSES
the mechanism. Its Figure 1 puts needle-probe measurements systematically and
significantly below every other method, and its section 4 then says the
microstructural-disturbance explanation was reached with an inadequate choice of
heating time and inference window, that modelling an air gap around the needle
shows no significant effect once the first thirty seconds are discarded, and
that thermal contact has been ruled out experimentally. What it offers instead
is that the three methods do not measure the same variable: the needle probe
assumes a homogeneous isotropic medium and works in transient mode, the
flux-gradient method sees only the vertical component in steady state, and the
numerical model counts only two of the processes that move heat through snow. On
that reading the gap is not an instrument defect at all but the non-conductive
processes showing up differently in each method. Fourteau et al. (2021) cites
both papers together for needle-probe bias, which is how the disagreement went
unnoticed.

Either way the conclusion for this row is the same and it does not rest on the
mechanism: Sturm's regression is not the relation to adopt. What changes is that
the reason is method-versus-method rather than a demonstrated instrument fault.

**And the physics the needle probe misses is a real term in this model's snow.**
A snow pack under a temperature gradient moves heat by conduction through the
ice matrix AND by water vapour diffusing through the pore space and releasing
latent heat where it deposits. Fourteau et al. (2021) computes the effective
conductivity WITH that term, on tomographic microstructures, and parameterises
the result as a quadratic in the ice volume fraction at each of five
temperatures. That is the relation adopted, evaluated at 263 K, the middle of
the five.

Three consequences worth recording:

1. **The shipped value was not wrong.** At the shipped snow density Fourteau's
   relation returns within a few per cent of what `landmod` and `icemod` already
   carried. Sturm's returns very nearly half of it. **A row that had implemented
   Sturm as named would have made the model worse while appearing to fix a
   defect**, and would
   have deepened a bias the error budget already carries in the same direction:
   `scripts/error_budget.py` holds "sea ice does not move" and the peer running
   this configuration reports excessive ice, and halving the conductivity of the
   snow on that ice insulates it further.
2. **The defect was real anyway, and it was the relation and not the value.**
   The conductivity was pinned while the density beside it was a namelist key,
   so a bracket on `rhosnow` moved the pack's thickness and its thermal mass and
   left the heat flux through it, which goes as `k/z`, on a value belonging to
   some other snow. It is the DEPENDENCE that was missing.
3. **The temperature is a declaration and it is second-order.** Across the whole
   223 to 273 K span Fourteau tabulates, the value at the declared density moves
   by about six per cent, against a factor of two between published relations,
   so the choice of temperature is not what this constant turns on.
   Making it a per-cell function would be a change to the snow scheme, not to
   this constant, and is not done here.

## The bracket, as a number

The unresolved question the literature itself leaves open is whether snow's
vapour deposition kinetics are fast or slow. Fourteau's section 2 poses both
limits and its section 4.1 says plainly that it is unclear which applies, citing
Krol and Lowe (2016) for isothermal metamorphism looking like the slow limit and
temperature-gradient metamorphism looking like the fast one. Its section 2.1
names Calonne et al. (2011) as the slow limit treated in detail, which is what
that paper is: it computes the effective conductivity on tomographic
microstructures with conduction through the ice and the interstitial air ONLY,
and says so in its abstract.

So the two arms are the same computation on the same kind of data, differing in
which processes they count, and the bracket between them is now a ratio rather
than a direction. `analysis/ice_properties.py` evaluates both across the density
range the modelled snow reaches, which is set by the two components between
them: ExoPlaSim declares one snow density and LPJ-GUESS ages its snow across a
range that contains it. The shape of the result:

- **The fast arm is above the slow arm everywhere**, by about a third at the
  bottom of the density range and about an eighth at the top. That narrowing
  with density is Fourteau's own Figure 6, which reports the relative difference
  as larger for low-density snow and for warmer snow, and the ratios here sit
  inside the band that figure gives near the melting point.
- **The bracket is wider than the noise of the fit that defines its lower arm.**
  Calonne states a residual standard deviation for his regression, and the gap
  between the two arms at the declared snow density is more than twice it. So
  this is a real disagreement between two published limits and not the scatter
  of one regression, which is the test that has to be passed before a bracket is
  worth reporting at all.
- **The adopted value is an ENDPOINT of the bracket, not a point inside it.**
  `landini` evaluates Fourteau's fast-kinetics row, so the question the bracket
  poses is one-sided: the modelled snow could conduct less than the model says,
  by up to the ratios above, and it cannot conduct more.
- **Yen (1981)'s own snow regression lands inside the bracket** at every density
  carried, between the two arms. That is not an accident of fitting: Yen states
  that a measured snow conductivity "includes vapor diffusion", so his curve is
  an APPARENT conductivity counting the very term that separates the two arms.
  Calonne's section 3.1 reports his purely conductive data agreeing closely with
  Yen's curve, and this note carries the third relation because two published
  claims that agree for different reasons are worth having visible.

The temperature choice inside the fast arm remains second-order against this,
which was already the finding: Fourteau's five temperature rows move the value
at the declared density by a few per cent, against the arm-to-arm gap above.

**What the bracket is worth in the model is not measured here.** It is a
one-signed change to the conductive resistance of every snow layer, on land and
on sea ice, and the run that would price it is named below.

## A cross-component finding: the two columns ran different snow

The bracket above is between two published limits of one quantity. The gap
inside this project was larger, and it was found while pricing the bracket:
`vendor/lpj-guess/modules/soil.cpp` computed its snow conductivity from Sturm et
al. (1997) -- the relation this row declined -- while `landmod` computed it from
Fourteau. Sturm's sits BELOW the slow arm at every density in the range, and at
the declared density it is close to half the value the climate column uses.

So the same snowfall, at the same density, insulated the ecology column's soil
about twice as well as the climate column's. That is not a bracket arm and it is
not a difference of opinion about snow: it is one material property with two
values inside one project, which is the class
`notes/audits/ocean-tier-implicit-earth.md` finding A1 was written about and the
class `CKAPSN` was below.

**It is one relation now, and one declaration rather than two agreeing copies.**
`lib/snow.py` states Fourteau's Eq. (18) once, and both compiled models restate
the adopted row as a literal because neither can import a Python module at
runtime. `snow.check_restatements()` holds each literal to that table, and
`scripts/smoke_test.py` and `biosphere/scripts/snow_thermal_gate.py` both run
it. `analysis/ice_properties.py` imports the relation rather than carrying a
third copy of it. WORLD-GJOV.

**Matching the numbers would not have been an answer.** The vegetation model's
snow density is prognostic, ramped from `snowdens_start` to `snowdens_end` by
its compaction scheme, where `landmod`'s `rhosnow` is a single namelist key
sitting inside that span. Equal conductivities at one density would have been a
coincidence at one point of two curves that diverge everywhere else, which is
what makes this separable from the density bracket GRAV-8 owns.

**What the ecology column's snow is worth after the change**, at the same water
equivalent, across the density range that model reaches: its conductivity rises
by a factor of 1.50 at the settled end of the compaction ramp and 2.30 at the
fresh end, so the insulating resistance of a modelled pack -- its thickness over
its conductivity -- falls to between 44 and 67 per cent of what it was, and the
temperature drop across the pack per unit ground heat flux falls with it. The
snow's thermal diffusivity rises by the same factor, so a modelled pack also
tracks the simulated air temperature faster. Which way the simulated soil
temperature moves is NOT settled by that, because it depends on the sign of the
gradient across the pack: the change couples the modelled soil more tightly to
the air above it in both directions. Decomposition reads the soil temperature,
and `biosphere/config/snow_thermal.yaml` names the run that would price it,
which cannot be made here because LPJ-GUESS does not build until a baseline
climatology exists.

The divergence is declared. `update_snow_properties` arrived byte-identical to
`guess_4.1/modules/soil.cpp`'s apart from a stripped licence header and is on by
default, so this is a change to default-on behaviour in a widely used community
model rather than a fork quirk; mainline's own lines stand verbatim beside the
changed one in the source, and the gate checks that the record is there, that it
is not what the model runs, and that the changed line still is.


## The third instance: the vegetation model's snow carried ice's heat capacity

`update_snow_properties` computed the modelled snowpack's volumetric heat
capacity as ice Ih's specific heat times a local `ice_density` of 917 kg/m3,
with the prognostic `snowdens` declared three lines above it and not read. A
volumetric heat capacity is a density times a specific heat, and the density in
question is the density of the substance occupying the volume.

**Only the density is wrong.** The SPECIFIC heat is ice Ih's and that is right
for this pack: snow is ice plus air, and at every density this model's snow
reaches the air carries under a thousandth of the mass, so a kilogram of pack
stores what a kilogram of ice stores. This is the same sentence that makes the
snow's melting enthalpy derivable and its conductivity not, and it is why this
repair needs no measurement: the argument is dimensional.

**What the wrong density did.** Snow enters that model as a water equivalent,
and the layer thickness its numerical solve gets is the water equivalent divided
by `snowdens`, so the thickness goes as 1/rho. A pack's thermal mass is its
thickness times its volumetric heat capacity, so at fixed water equivalent it
should not depend on the density at all: the two factors cancel. Pinning the
capacity at solid ice's broke the cancellation and left the modelled pack's
thermal mass per unit water equivalent a factor 917/`snowdens` too high, which
is 3.33 at `snowdens_start`, 1.83 at `snowdens_end` and 3.67 under
`snowdensityconstant`. `Ci[]` takes it directly for every active snow layer of
the multilayer solve, so that was the capacity of the whole modelled pack, and
`Dsnow = Ksnow/Csnow` carried the reciprocal into the diffusivity: at 273.15 K
the modelled pack's thermal diffusivity rises from 1.24e-7 to 4.12e-7 m2/s at
`snowdens_start` and from 3.50e-7 to 6.42e-7 at `snowdens_end`, so the depth a
temperature wave of a given period penetrates rises by 1.83 and 1.35. The
modelled pack was both too slow to warm and too slow to cool.

**It is one defect at a fourth site, and this is the one it reached in a second
component.** `snowcap` was a namelist key held fixed while `rhosnow` moved;
`snowdiff` was a fixed conductivity beside the same moving density; `sicecap`
factorised as ice's specific heat times LIQUID WATER's density; this one reads
solid ice's density beside a prognostic snow density three lines above it. All
four are now derived from the density of the substance whose volume they
describe, at the site where that density is set. The three in `landmod` are
`landini` and `glacierprep`; this one is `update_snow_properties`.

**A check that could have failed and did not.** The two columns' repaired lines
have the same CONSTRUCTION: the pack's own density times the specific heat of
the ice in it. What they did not then share was the specific heat, and the
section below is that residual.

**What is NOT settled by this.** Which way the simulated soil temperature moves,
for the same reason it is not settled for the conductivity: a faster pack
couples the modelled soil more tightly to the air above it in both directions,
and the sign depends on the gradient. Decomposition reads the soil temperature.
`biosphere/config/snow_thermal.yaml` names the arm that would price both
divergences together and says why one arm cannot separate them.

The divergence is declared beside the conductivity one, with mainline's own
lines standing verbatim in the source, and `biosphere/scripts/snow_thermal_gate.py`
checks that the record is there, that it is not what the model runs, and that
the changed line still is. WORLD-2AIJ.

## The third instance again: two specific heats of one ice

The two columns' constructions agreed and their SPECIFIC HEATS did not.
`landmod` declared `real, parameter :: CPSNOW = 2090.` with no citation;
`soil.cpp` evaluated Fukusako's linear relation in absolute temperature,
`(0.185 + 0.689 T * 0.01)` kJ/kg/K, which gives 2067.0 J/kg/K at 273.15 K and
1791.4 at 233.15 K. They agree to 1.1 per cent at the melting point and are 16.7
per cent apart at -40 degC, and only one of them varies with the simulated
temperature at all.

**They are not a disagreement to average, because one of them is the other
evaluated somewhere.** 2090 is Fukusako's line at 276.49 K and IAPWS-06's ice at
272.24 K. Both of those are at or above the melting point of the modelled snow,
so the fixed value is the relation evaluated where this snow only sits at the
moment it melts, applied at every temperature it reaches: 3.4 per cent high at
263 K and 15.8 per cent high at 233 K. That is the answer to "is the constant
the relation at a temperature this planet's modelled snow actually reaches" --
it is not.

**Neither number had to win.** IAPWS R10-06(2009), the equation of state 2006
for H2O ice Ih, gives the specific heat exactly, and this project already holds
it: `analysis/ice_properties.py` implements the Gibbs function and reproduces
every quantity at every state of the release's own Table 6 to 7e-10 relative,
and `glaciermod`'s `CPGLAC` already came from it. Fukusako (1990) is held and
read, and it carries the fact that decides the fixed value: its Eq. (2) is
stated for 273 K >= T >= 90 K, so the climate column's 2090 J/kg/K is that line
evaluated at 276.49 K -- 3.3 K past the correlation's own upper limit, and above
the melting point of the substance it describes. Against the standard,
Fukusako's line is 0.30 per cent low at 180 K and 1.42 per cent low at the
melting point.

**A Gibbs function in complex arithmetic is not something two compiled models
restate**, so `lib/snow.py` declares the closed form they do restate: a
least-squares quadratic in temperature, derived from IAPWS-06 over 170 K to the
triple point and within 0.554 J/kg/K of it there. 170 K is below anything the
simulated air over the modelled snow has reached -- the coldest monthly bin of
the two bounding climates is -75.8 degC, 197.4 K -- and the form degrades
gracefully rather than sharply outside it, staying within 0.9 J/kg/K of the
standard down to 140 K. The representation error is two orders below either
number it replaced.

**Each column evaluates the relation where it has a temperature.** `landmod`'s
`snowcap` is one compiled scalar, so `landini` evaluates the quadratic at
`TSNOWREF`, which is the same 263 K the conductivity row is taken at, so that
column's snow is one material stated at one temperature. `soil.cpp` has a
simulated daily air temperature in hand and evaluates it there. Sharing the
RELATION is what makes the two columns' snow one material; sharing a value would
have made them one material only on the days the ecology column's snow happened
to sit at 263 K.

**What it is worth.** The ecology column's specific heat rises by 0.30 per cent
at 180 K and 1.42 per cent at the melting point, and its thermal diffusivity
falls by the same fraction. The climate column's falls from 2090 to 2022.05
J/kg/K, so `snowcap` at the declared `rhosnow` of 330 falls from 689,700 to
667,275 J/m3/K, 3.25 per cent, and that column's modelled pack warms and cools
3.25 per cent faster per unit water equivalent. Neither is why this was done:
what it settles is that the two columns' snow is one material, and the
cross-column difference removed is the larger number, 16.7 per cent at -40 degC.

**The loop that keeps it a derivation.** `scripts/smoke_test.py` re-derives the
specific heat from `analysis/ice_properties.py`'s IAPWS-06 across the declared
domain and refuses if the closed form misses by more than the declared residual,
and `lib/snow.py`'s `check_restatements()` holds both compiled literals and both
of `landmod`'s placeholder defaults to the one declaration. Edit a coefficient,
widen the domain, or move `TSNOWREF` off the conductivity's row and a gate says
so. WORLD-A2LV.

## What `CKAPSN` was, which the earlier audit had missed

`icemod`'s `CKAPSN` stood at the same 0.31 `landmod` declared as `snowdiff`, and
both were used for the same thing: the conductive resistance of a snow layer,
its thickness over its conductivity. So the snow on the modelled sea ice and the
snow on the modelled soil were two independent statements of one material
property, with nothing comparing them. That is the class
`notes/audits/ocean-tier-implicit-earth.md` finding A1 was written about, and it
had been missed because only the DENSITY was declared twice in a way an earlier
pass could see. `iceini` takes it now, beside the snow density it is derived
from.

## `CPSN` reached nothing

`parameter(CPSN = 2090.)` was declared in `icemod.f90` and appears nowhere else
in the vendored tree -- not in an expression, not in a namelist, not in a
threadprivate list. Deleting it is the fix. Exposing it as WORLD-04OK's "Do:"
line asked would have put a key in `icemod_nl` that a bracket arm could move
with no effect whatever, which is a worse state than the one it was in.
`landmod`'s `landini` carries the live statement of that quantity -- IAPWS-06's
specific heat of ice Ih at `TSNOWREF`, restating `lib/snow.py` -- and is what
`snowcap` is derived from. 2090 was that relation a degree below the melting
point.

## The glacial-ice pair

`landmod.f90` carries two material constants of the modelled GLACIAL ice, for
the ice `glaciermod` grows and that weights the soil column's thermal properties
by `dglac`. Both were `landmod_nl` keys standing beside a density they did not
follow. Both are now derived from that density in `glacierprep`, which is where
the density is declared and which runs before `landini` in `surfini`, and
`analysis/ice_properties.py` computes both and holds the compiled literals to
what it computes.

### The heat capacity was computed at liquid water's density

`sicecap` factorised exactly as one thousand times the specific heat of ice, and
glacial ice is not water. **This was GRAV-8's defect, third instance.** GRAV-8
made `snowcap` follow `rhosnow` because a volumetric heat capacity is a density
times a specific heat, and WORLD-A9S5 made `snowdiff` follow it for the same
reason; `sicecap` followed nothing, so a bracket that moved `rhoglac` moved the
ice orography and left the thermal mass of the ice behind.

It is now `rhoglac` times the specific heat of ice Ih from IAPWS-06 at a
declared temperature, which is 17 per cent below the superseded value at the
shipped density.

### The conductivity is pure ice's, reduced for the air the density implies

`sicediff` stood at the same number as `icemod`'s `CKAPI` and was unsourced.
IAPWS-06 cannot settle it: a Gibbs function carries density, specific heat and
compressibility and no transport property at all. Yen (1981) does, in two steps
that are both his:

- **Pure ice, Eq. (33)**, `lambda = a exp(b T)`, the whole-range arm of his
  Table 3. He recommends it for practical use because it has the highest
  correlation coefficient of his three arms, 0.9313 against 0.5962 for the arm
  fitted above the 150-195 K data gap. At the declared temperature pure ice is
  about a twelfth ABOVE the superseded `sicediff`.
- **Bubbles, Eq. (37)**, `2 rho / (3 rhoice - rho)`, Schwerdtfeger's reduction of
  Maxwell's effective-medium result for randomly distributed spherical air
  inclusions once the conductivity of air is dropped against the ice's. Yen
  states the unreduced form twice -- Eq. (36) for dense snow and Eq. (70) for the
  bubbly ice inside sea ice -- and the two agree to two parts in a thousand at
  the shipped density, which the script checks rather than assumes. His Figure 22
  is the same equation drawn. Glacial ice IS bubbly, and at `rhoglac` against
  Yen's 917 kg/m3 pure ice the air fraction is about 7 per cent.

The product lands about 3 per cent BELOW the superseded value, so the correction
and the reduction pull opposite ways and nearly cancel. That near-cancellation is
a coincidence of the shipped density and is exactly why the two steps are carried
separately: a bracket on `rhoglac` moves the second and not the first.

### It is still not `CKAPI`, and Yen sharpens why

Sea ice is brine-bearing and glacial ice is fresh. Yen's Eq. (71) subtracts a
brine term from **exactly the bubbly ice computed above**, so the modelled sea
ice conducts less than modelled glacial ice of the same density, and the two
constants sit on opposite sides of pure ice's value for different reasons.
Deduplicating them would assert that this world's glaciers are salty.

### The declared temperature is the larger uncertainty, not the density

Both constants are evaluated at one declared temperature, and it is the
temperature `landmod` already evaluates the snow conductivity at, so the two
cryosphere materials are stated at one temperature rather than two. Over 233.15 K
to the melting point the conductivity spans a factor of about 1.26 and the heat
capacity about 1.16, against about 1.04 for the difference between Yen's two
regression arms at one temperature. So the temperature dominates, and the pair
is a **declaration with a bracket** rather than a derivation that leaves nothing
open.

Making the pair a function of the layer's own temperature is a change to the
soil heat solver's material model rather than to a constant -- the same shape of
change WORLD-JSFM makes for the soil half -- and it is named here rather than
made. Yen is himself inconsistent at the 1.3 per cent level about pure ice at the
melting point, his Eq. (33) and the 2.09 W/m/K his own sea-ice model uses
disagreeing by that much, which is a second reason the temperature has to be
stated rather than inherited.

### What it is worth today is nothing, and that is measurable

`dglac` is identically zero on every cell of the climatology this world has, so
neither constant reaches a result until `glaciermod` grows ice. The run that
would price them is therefore not a bracket arm at all: it is any run that
produces a glaciated cell.

## What snow compaction is worth, and where

Measured 2026-08-26 by `analysis/ice_properties.py --climatology`, offline
against an existing bootstrap climatology and before any change to the model.
GRAV-8's own sequencing, and the reason for it is that adding a prognostic
density imports three viscosity constants with no stated derivation, so the
question of what they buy has to be answered first. A bootstrap's numbers are
not the baseline, and this is labelled as one of the two places in this note
that reads one.

### The scheme, and what it costs to adopt

PALADYN's, Willeit and Ganopolski (2016) Eqs. (46) to (48): self-loading after
Kojima (1967) as implemented by Pitman et al. (1991), fresh-snow density after
Anderson (1976).

    d(rho)/dt = 0.5 g rho w / eta  +  P (rho_fresh - rho) / w
    eta       = eta_0 exp[k_T (T_0 - T_sn) + k_rho rho]
    rho_fresh = rho_min + 1.7 (T_a - T_0 + 15)**1.5

Gravity enters ONE term, linearly, so the compaction RATE runs 1.306 times
Earth's here. The three viscosity constants are tabulated with a description and
no source column, no range of validity and no sensitivity test, and reach the
paper from a 1991 technical report's implementation of a 1967 conference paper.
Adopting the scheme imports all three.

### The rate ratio is not the state ratio, and that is the first finding

Prognostic density over this world's modelled snow, at the median snowy land
cell and time bin: 161 kg/m3 at Earth's gravity, 170 at this world's. The
compaction rate is 1.306 times Earth's and the DENSITY it produces is 1.06
times, because the pack also relaxes toward the fresh-snow density at a rate set
by the snowfall and this world's modelled snow does not survive long enough for
the load term to dominate. Reading the rate ratio as a state ratio would
overstate the effect fivefold.

### The constant is worth six times what the gravity term is

The quantity is the pack's CONDUCTIVE RESISTANCE, thickness over conductivity,
which is what sets the temperature drop across the pack per unit ground heat
flux and is what `landmod` builds `zdiff1` from. Depth and conductivity both
move with the density and in opposite directions, so quoting either alone says
nothing.

| Case | Density | Physical depth | z/k |
| --- | --- | --- | --- |
| the compiled constant | 330 kg/m3 | 0.092 m | 0.291 m2K/W |
| prognostic at Earth's gravity | 161 kg/m3 | 0.204 m | 2.044 m2K/W |
| prognostic at this world's gravity | 170 kg/m3 | 0.189 m | 1.780 m2K/W |

Medians over 4030 snowy land cell-bins. The prognostic pack is about six times
more insulating than the constant makes it, and the gravity term inside that is
worth 14 per cent. **The model form is the large term and the gravity is the
small one**, which is the opposite of the order the row was filed in.

### The gravity term is smaller than a bracket the model already carries

At the prognostic density the gravity term moves the resistance by 16 per cent.
The vapour-kinetics bracket on the conductivity -- Fourteau's fast arm against
Calonne's slow arm at ONE density, which is an open question rather than an
error bar -- moves it by 53 per cent. So a pair of runs differing only in
gravity would not separate the gravity term from the arm the model happens to
evaluate. That is not an argument against representing compaction; it is an
argument that the gravity term is not what a run should be bought to measure.

### On albedo the answer is exactly zero, and it is a fact about the model

The row was filed expecting the density to reach surface albedo through snow
cover depth. It does not. `landmod` takes the snow-covered fraction as
`dsnow/(dsnow + snowcovz)` in WATER EQUIVALENT, so the density cancels out of it
entirely. The one path from density to albedo is the canopy burial depth handed
to `snowcanopymask`, and that argument is read only when `forhgt` is positive,
which it is not. So the density reaches the modelled surface albedo through
nothing at all, and a compaction term is worth exactly zero there until the
canopy gets a height under GRAV-7 and BIO-33.

Peer practice does not rescue that path either. ClimaLand multiplies snow albedo
by `min(1 - beta*(rho/rho_liq - x0), 1)`, and its `beta` is a free parameter with
no citation: 0.97 in its calibrated parameter set and 0 in its uncalibrated one.
Adopting it would import a tuned constant rather than a mechanism. PALADYN puts
no density in albedo at all and uses a snow AGE factor as the grain-size proxy
instead. So of the two peer schemes, one has no such term and the other's is
tuned.

### Where the density does still reach a result

The conductive resistance above, on land; the top soil layer's blended heat
capacity, but only where the physical pack is deeper than `dztop`, which is 16
per cent of snowy land cell-bins under the constant and more under a prognostic
density, because below that depth `snowcap * zsntop` is the water equivalent
times the specific heat of ice and the density cancels; and the modelled sea
ice, which takes the same density through `iceini` and uses it the same way.

### What this measurement cannot say

- The forcing is a climatology, so an intermittent pack appears as a persistent
  thin one. Compaction is linear in the load, so a time-mean water equivalent
  understates the compaction of a transient deep pack: those cells' densities
  here are a floor.
- The water equivalent is PRESCRIBED from a run that used the constant density,
  so this prices a density on that snow rather than on the snow a prognostic
  density would itself produce. The feedback is one-signed and named: a less
  dense, more insulating pack keeps the ground warmer under it.
- The snow layer temperature is proxied by the surface temperature capped at
  melting, which is the pack's warm end and therefore the soft end of the
  viscosity, so the compaction reported here is an upper bound.
- PALADYN itself neglects metamorphism and the effect of melting on density.

## What this did not establish

- **Nothing was run.** Every statement is against the source, the standard, or a
  published relation, with one exception that is labelled where it appears: the
  temperature span of the modelled sea ice is read off an EXISTING bootstrap
  climatology, and a bootstrap's numbers are not the baseline. The model was not
  compiled beyond a syntax check and no climate run was made, so **what these
  changes are worth in kelvin or in ice thickness is not measured here.**
- **The run that would measure it**, when the binaries are current: a pair of
  T21 baseline arms differing only in `landmod_nl`'s `rhosnow` -- the GRAV-8
  bracket that this row exists to make honest -- held at the same rung, the same
  timestep and the same orbit count, with everything else fixed. Before this
  change that pair moved the pack's thickness and thermal mass and held its
  conductivity; after it, all three move. The quantity to read is the ground and
  basal heat flux under snow-covered cells and the resulting ice thickness.
  **The pricing above says what that pair should be**, and it is not a gravity
  bracket: the two arms worth buying are the compiled constant against a
  prognostic density, which differ by a factor of six in the pack's conductive
  resistance, rather than the same prognostic density at two gravities, which
  differ by less than the kinetics bracket already open at one density.
- **What the snow bracket is worth is not measured either.** The run that would
  measure it is a second pair, on the same terms as the `rhosnow` pair above and
  not folded into it: two T21 baseline arms differing only in which arm of the
  kinetics bracket `landini` evaluates -- Fourteau's row against Calonne's
  equation (12) -- at one density, one rung, one timestep and one orbit count.
  The quantity to read is the same: ground and basal heat flux under
  snow-covered cells, and the ice thickness that follows. Running the two
  brackets as one four-arm sweep would confound them, because both act on the
  same conductive resistance.
- **The glacial pair reaches nothing yet, and that is checked rather than
  assumed.** `dglac` is identically zero on every cell of the climatology this
  world has, so `sicecap` and `sicediff` weight nothing in `tands`. What would
  price them is not an arm of a bracket but any run that produces a glaciated
  cell; until one exists, the temperature these two are declared at is the open
  question and no sweep of it can be read anywhere.
- **The sea-ice declarations are still not bracketed, and Yen says how wide the
  bracket would be.** Three of the six are stated positions, and the quantity
  that would set a range for them is a brine volume the model does not carry.
  Yen's Figure 23 shows what that costs for the conductivity: across the
  salinities and densities of ordinary sea ice the value spans nearly a factor
  of two, so the declaration is not a small approximation dressed as a choice.
  Giving it a bracket still needs an ice salinity, which is a change to the
  sea-ice scheme and not to a constant, and the size above is the argument for
  making it rather than a substitute for it.

## Sources

| Source | Status |
| --- | --- |
| `references/pdf/iapws_2009_revised-release-on-the-equation-of-state-2006-for-h2o-ice-ih.pdf` -- IAPWS R10-06(2009), *Revised Release on the Equation of State 2006 for H2O Ice Ih* | **read** -- Eq. (1), Table 2's coefficients and Table 6's numerical check values. `analysis/ice_properties.py` implements it and reproduces every entry of Table 6 to 7e-10 relative before reporting anything. Stated valid from 0 K to the melting curve and to 210 MPa, which every state used here is inside by a wide margin |
| `references/pdf/feistel_2006_a-new-equation-of-state-for-h2o-ice-ih.pdf` -- Feistel, Wagner (2006). J. Phys. Chem. Ref. Data 35, 1021-1047. `10.1063/1.2183324` | **held** -- the paper behind the release above. The release is self-contained for what is used here, so this is the derivation rather than the specification |
| `references/pdf/TEOS-10_Manual.pdf` -- IOC, SCOR, IAPSO (2010), Appendix I | **read in part** -- restates the same ice Ih Gibbs function and its coefficients, and was where they were found before the release itself was obtained. Its version is written in SEA pressure, the release's in absolute; the two agree once that is reconciled, and getting it wrong is a 1e-5 relative error that Table 6 catches |
| `references/pdf/sturm_1997_the-thermal-conductivity-of-seasonal-snow.pdf` -- Sturm, Holmgren, Konig, Morris (1997). J. Glaciol. 43(143), 26-41. `10.3189/S0022143000002781` | **read** -- the quadratic and its stated range, its R2 of 0.79, its 0.1 W/m/K uncertainty at 95 per cent confidence, and its statement that the regressions are strictly valid only near -14.6 C. NOT ADOPTED, for the reasons above |
| `references/pdf/riche_2013_thermal-conductivity-of-snow-measured-by-three-independent-methods-and.pdf` -- Riche, Schneebeli (2013). The Cryosphere 7, 217-227. `10.5194/tc-7-217-2013` | **read** -- three methods on identical samples, the conclusion that direct numerical simulation is the most reliable, and the up to plus or minus 25 per cent anisotropy error on a horizontally inserted needle probe. This is what decides against Sturm |
| `references/pdf/fourteau_2021_impact-of-water-vapor-diffusion-and-latent-heat-on-the-effective-therm.pdf` -- Fourteau, Domine, Hagenmuller (2021). The Cryosphere 15, 2739-2755. `10.5194/tc-15-2739-2021` | **read** -- Eq. (18), the vertical effective thermal conductivity under fast kinetics at five temperatures as a quadratic in the ice volume fraction, with 917 kg/m3 as the normalising ice density. ADOPTED. Also the statement that the fast and slow kinetics limits both remain plausible, which is the bracket recorded above |
| `references/pdf/calonne2011-effective-thermal-conductivity-of-snow.pdf` -- Calonne, Flin, Morin, Lesaffre, Rolland du Roscoat, Geindreau (2011). Geophys. Res. Lett. 38, L23501. `10.1029/2011GL049234` | **read** 2026-08-25, all six pages. SUPPLIED BY THE USER after every open-access route and Sci-Hub returned 403. Equation (12), its quadratic in density fitted so that the value goes to air's at zero density, with its correlation coefficient and the residual standard deviation that decides whether the bracket is bigger than the noise; the statement that only conduction through ice and interstitial air is counted, which is what makes it the slow arm; the table showing that dropping conduction through the pore AIR would lower the answer by a factor of two in dense snow and an order of magnitude in fresh snow; and section 3.1's agreement with Yen's snow curve. ITS NEEDLE-PROBE POSITION IS NOT WHAT THE EARLIER ROW EXPECTED, and section 4 above says what it is instead |
| `references/pdf/willeit_2016_paladyn-v1-0-a-comprehensive-land-surfacevegetationcarbon-cycle-model.pdf` -- Willeit, Ganopolski (2016). Geosci. Model Dev. 9, 3817-3857. `10.5194/gmd-9-3817-2016` | **read in part** 2026-08-26: Sect. 5.2's snow model with Eqs. (46) to (49), Table 2's surface model parameters, Sect. 3.1's snow albedo with Eqs. (17) to (20) and Appendix A's snow age factor, and Sect. 4's Eqs. (35) and (36). The compaction scheme priced above, and the source of the finding that PALADYN puts no density in its snow albedo at all. Its three viscosity constants carry no derivation, no range of validity and no sensitivity test, which is what makes adopting the scheme a cost rather than a free improvement |
| `references/pdf/yen1981-review-of-thermal-properties-of-snow-ice-and-sea-ice.pdf` -- Yen (1981). *Review of thermal properties of snow, ice and sea ice.* CRREL Report 81-10 | **read** 2026-08-25, the conductivity sections: "Thermal conductivity of ice" with equation (33) and Table 3's three regression arms, "Thermal conductivity of snow" with equation (34), and "Density and thermal conductivity of sea ice" and "Thermal conductivity model for sea ice" with equations (70) to (72) and Figures 22 and 23. This is `CKAPI`'s bound, which IAPWS-06 cannot give because a Gibbs function carries no transport property. Also the sea-ice conductivity model itself, which is the arithmetic behind the sentence that these constants are not derivable here: it needs a salinity and a temperature per cell, and shows the brine term SUBTRACTING |
