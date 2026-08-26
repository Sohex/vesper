# What the modelled sea ice and its snow are made of

**Derived:** 2026-08-25, by `analysis/ice_properties.py` from IAPWS R10-06(2009)
and from the papers named at the foot of this note. No World Orogen generation
and no ExoPlaSim run: every number here is either a laboratory standard
evaluated, or a published relation evaluated at a declared density or
temperature. One measured number is quoted from an existing bootstrap
climatology and is labelled where it appears.

This is worldbuilding. Vesper is an invented super-Earth around a mid-K dwarf,
and every quantity below is a material property of that planet's simulated sea
ice or of the simulated snow that lies on it. IAPWS-06 and the snow papers are
laboratory measurements of MATERIALS, and a material does not change when the
star does, which is the whole of why they transfer here.

WORLD-04OK and WORLD-A9S5 asked, of six compile-time constants in
`icemod.f90` and of one namelist key in `landmod.f90`, what each one's source
is and whether it should move for this world. This note is the answer, one row
at a time.

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

## A cross-component finding: the two columns run different snow

The bracket above is between two published limits of one quantity. The gap
inside this project is larger, and it was found while pricing the bracket:
`vendor/lpj-guess/modules/soil.cpp` computes its snow conductivity from Sturm et
al. (1997) -- the relation this row declined -- while `landmod` now computes it
from Fourteau. Sturm's sits BELOW the slow arm at every density in the range,
and at the declared density it is close to half the value the climate column
uses.

So the same snowfall, at the same density, insulates the ecology column's soil
about twice as well as the climate column's. That is not a bracket arm and it is
not a difference of opinion about snow: it is one material property with two
values inside one project, which is the class
`notes/audits/ocean-tier-implicit-earth.md` finding A1 was written about and the
class `CKAPSN` was below. It is `lshy-5`'s territory rather than this row's, and
it is recorded rather than resolved here because changing which relation the
biosphere runs is a change to `vendor/lpj-guess`, not to a constant.


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
`landmod`'s `CPSNOW`, at the same value, is the live declaration of that
quantity and is what `snowcap` is derived from.

## The glacial-ice pair, found beside these and not settled here

`landmod.f90` carries two more material constants of modelled ice, for the
GLACIAL ice that `glaciermod` grows and that weights the soil column's thermal
properties by `dglac`. Both were looked at while settling the six above and
neither is in scope for the ocean tier, so they are recorded here as a finding
and are not changed.

- **`sicecap`, the heat capacity per unit volume of the modelled glacial ice,
  is computed at LIQUID WATER'S DENSITY.** Its value factorises exactly as one
  thousand times the specific heat of ice, and glacial ice is not water: at
  `glaciermod`'s own declared `rhoglac` the product is about a sixth smaller,
  and at pure ice Ih's density about a twentieth smaller. **This is GRAV-8's
  defect, third instance.** GRAV-8 exposed `rhoglac` as a namelist key so the
  gravity bracket could be run, and made `snowcap` follow `rhosnow` because a
  volumetric heat capacity is a density times a specific heat; WORLD-A9S5 made
  `snowdiff` follow it for the same reason. `sicecap` follows nothing, so a
  bracket that moves `rhoglac` today moves the ice orography and leaves the
  thermal mass of the ice behind.
- **`sicediff`, the conductivity of the modelled glacial ice, stands at the
  same number as `icemod`'s `CKAPI`.** These are NOT the same quantity and must
  not be deduplicated into one: glacial ice is fresh and sea ice is
  brine-bearing, and their conductivities differ by a term in the ice's
  salinity over its temperature. What is true is that neither is sourced, and
  that the fresh one is the one a standard could settle, because it is a pure
  substance in the way the sea-ice constants are not.

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
- **What the snow bracket is worth is not measured either.** The run that would
  measure it is a second pair, on the same terms as the `rhosnow` pair above and
  not folded into it: two T21 baseline arms differing only in which arm of the
  kinetics bracket `landini` evaluates -- Fourteau's row against Calonne's
  equation (12) -- at one density, one rung, one timestep and one orbit count.
  The quantity to read is the same: ground and basal heat flux under
  snow-covered cells, and the ice thickness that follows. Running the two
  brackets as one four-arm sweep would confound them, because both act on the
  same conductive resistance.
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
| `references/iapws_2009_revised-release-on-the-equation-of-state-2006-for-h2o-ice-ih.pdf` -- IAPWS R10-06(2009), *Revised Release on the Equation of State 2006 for H2O Ice Ih* | **read** -- Eq. (1), Table 2's coefficients and Table 6's numerical check values. `analysis/ice_properties.py` implements it and reproduces every entry of Table 6 to 7e-10 relative before reporting anything. Stated valid from 0 K to the melting curve and to 210 MPa, which every state used here is inside by a wide margin |
| `references/feistel_2006_a-new-equation-of-state-for-h2o-ice-ih.pdf` -- Feistel, Wagner (2006). J. Phys. Chem. Ref. Data 35, 1021-1047. `10.1063/1.2183324` | **held** -- the paper behind the release above. The release is self-contained for what is used here, so this is the derivation rather than the specification |
| `references/TEOS-10_Manual.pdf` -- IOC, SCOR, IAPSO (2010), Appendix I | **read in part** -- restates the same ice Ih Gibbs function and its coefficients, and was where they were found before the release itself was obtained. Its version is written in SEA pressure, the release's in absolute; the two agree once that is reconciled, and getting it wrong is a 1e-5 relative error that Table 6 catches |
| `references/sturm_1997_the-thermal-conductivity-of-seasonal-snow.pdf` -- Sturm, Holmgren, Konig, Morris (1997). J. Glaciol. 43(143), 26-41. `10.3189/S0022143000002781` | **read** -- the quadratic and its stated range, its R2 of 0.79, its 0.1 W/m/K uncertainty at 95 per cent confidence, and its statement that the regressions are strictly valid only near -14.6 C. NOT ADOPTED, for the reasons above |
| `references/riche_2013_thermal-conductivity-of-snow-measured-by-three-independent-methods-and.pdf` -- Riche, Schneebeli (2013). The Cryosphere 7, 217-227. `10.5194/tc-7-217-2013` | **read** -- three methods on identical samples, the conclusion that direct numerical simulation is the most reliable, and the up to plus or minus 25 per cent anisotropy error on a horizontally inserted needle probe. This is what decides against Sturm |
| `references/fourteau_2021_impact-of-water-vapor-diffusion-and-latent-heat-on-the-effective-therm.pdf` -- Fourteau, Domine, Hagenmuller (2021). The Cryosphere 15, 2739-2755. `10.5194/tc-15-2739-2021` | **read** -- Eq. (18), the vertical effective thermal conductivity under fast kinetics at five temperatures as a quadratic in the ice volume fraction, with 917 kg/m3 as the normalising ice density. ADOPTED. Also the statement that the fast and slow kinetics limits both remain plausible, which is the bracket recorded above |
| `references/calonne2011-effective-thermal-conductivity-of-snow.pdf` -- Calonne, Flin, Morin, Lesaffre, Rolland du Roscoat, Geindreau (2011). Geophys. Res. Lett. 38, L23501. `10.1029/2011GL049234` | **read** 2026-08-25, all six pages. SUPPLIED BY THE USER after every open-access route and Sci-Hub returned 403. Equation (12), its quadratic in density fitted so that the value goes to air's at zero density, with its correlation coefficient and the residual standard deviation that decides whether the bracket is bigger than the noise; the statement that only conduction through ice and interstitial air is counted, which is what makes it the slow arm; the table showing that dropping conduction through the pore AIR would lower the answer by a factor of two in dense snow and an order of magnitude in fresh snow; and section 3.1's agreement with Yen's snow curve. ITS NEEDLE-PROBE POSITION IS NOT WHAT THE EARLIER ROW EXPECTED, and section 4 above says what it is instead |
| `references/yen1981-review-of-thermal-properties-of-snow-ice-and-sea-ice.pdf` -- Yen (1981). *Review of thermal properties of snow, ice and sea ice.* CRREL Report 81-10 | **read** 2026-08-25, the conductivity sections: "Thermal conductivity of ice" with equation (33) and Table 3's three regression arms, "Thermal conductivity of snow" with equation (34), and "Density and thermal conductivity of sea ice" and "Thermal conductivity model for sea ice" with equations (70) to (72) and Figures 22 and 23. This is `CKAPI`'s bound, which IAPWS-06 cannot give because a Gibbs function carries no transport property. Also the sea-ice conductivity model itself, which is the arithmetic behind the sentence that these constants are not derivable here: it needs a salinity and a temperature per cell, and shows the brine term SUBTRACTING |
