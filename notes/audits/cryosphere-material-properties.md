# What the modelled sea ice and its snow are made of

**Derived:** 2026-08-25, by `analysis/ice_properties.py` from IAPWS R10-06(2009),
and from the four snow-conductivity papers named at the foot of this note. No
World Orogen generation and no ExoPlaSim run: every number here is either a
laboratory standard evaluated, or a published relation evaluated at a declared
density.

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
| `CKAPI`, conductivity of the modelled sea ice | compile-time | **unchanged, DECLARED with its bound.** Same reason; the conductivity of sea ice is pure ice's plus a brine term |
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

**The needle probe has since been shown to be biased for snow.** Riche and
Schneebeli (2013) applied three independent methods -- a long-heating needle
probe, a guarded heat flux plate, and direct numerical simulation on the
microstructure -- to IDENTICAL samples, and concluded that the numerical
simulation is the most reliable of the three, with a horizontally inserted
needle probe wrong by up to a quarter either way through the pack's anisotropy.
Fourteau et al. (2021) restate that finding and cite Calonne et al. (2011) for
it as well.

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

The unresolved question the literature itself leaves open is whether snow's
vapour deposition kinetics are fast or slow; Fourteau studies both limits and
says plainly that it is unclear which applies. The slow limit is Calonne et al.
(2011)'s and is lower. That is a genuine bracket on this constant, and it is
recorded here rather than hidden inside a single number: the adopted value is
the fast-kinetics arm, and the slow arm sits below it.

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
  published relation. The model was not compiled beyond a syntax check and no
  climate run was made, so **what these changes are worth in kelvin or in ice
  thickness is not measured here.**
- **The run that would measure it**, when the binaries are current: a pair of
  T21 baseline arms differing only in `landmod_nl`'s `rhosnow` -- the GRAV-8
  bracket that this row exists to make honest -- held at the same rung, the same
  timestep and the same orbit count, with everything else fixed. Before this
  change that pair moved the pack's thickness and thermal mass and held its
  conductivity; after it, all three move. The quantity to read is the ground and
  basal heat flux under snow-covered cells and the resulting ice thickness.
- **The sea-ice declarations are not bracketed.** Three of the six are stated
  positions and no range is offered for what this world's sea ice might instead
  be made of, because the quantity that would set that range is a brine volume
  the model does not carry. Giving them a bracket needs an ice salinity, which
  is a change to the sea-ice scheme and not to a constant.

## Sources

| Source | Status |
| --- | --- |
| `references/iapws_2009_revised-release-on-the-equation-of-state-2006-for-h2o-ice-ih.pdf` -- IAPWS R10-06(2009), *Revised Release on the Equation of State 2006 for H2O Ice Ih* | **read** -- Eq. (1), Table 2's coefficients and Table 6's numerical check values. `analysis/ice_properties.py` implements it and reproduces every entry of Table 6 to 7e-10 relative before reporting anything. Stated valid from 0 K to the melting curve and to 210 MPa, which every state used here is inside by a wide margin |
| `references/feistel_2006_a-new-equation-of-state-for-h2o-ice-ih.pdf` -- Feistel, Wagner (2006). J. Phys. Chem. Ref. Data 35, 1021-1047. `10.1063/1.2183324` | **held** -- the paper behind the release above. The release is self-contained for what is used here, so this is the derivation rather than the specification |
| `references/TEOS-10_Manual.pdf` -- IOC, SCOR, IAPSO (2010), Appendix I | **read in part** -- restates the same ice Ih Gibbs function and its coefficients, and was where they were found before the release itself was obtained. Its version is written in SEA pressure, the release's in absolute; the two agree once that is reconciled, and getting it wrong is a 1e-5 relative error that Table 6 catches |
| `references/sturm_1997_the-thermal-conductivity-of-seasonal-snow.pdf` -- Sturm, Holmgren, Konig, Morris (1997). J. Glaciol. 43(143), 26-41. `10.3189/S0022143000002781` | **read** -- the quadratic and its stated range, its R2 of 0.79, its 0.1 W/m/K uncertainty at 95 per cent confidence, and its statement that the regressions are strictly valid only near -14.6 C. NOT ADOPTED, for the reasons above |
| `references/riche_2013_thermal-conductivity-of-snow-measured-by-three-independent-methods-and.pdf` -- Riche, Schneebeli (2013). The Cryosphere 7, 217-227. `10.5194/tc-7-217-2013` | **read** -- three methods on identical samples, the conclusion that direct numerical simulation is the most reliable, and the up to plus or minus 25 per cent anisotropy error on a horizontally inserted needle probe. This is what decides against Sturm |
| `references/fourteau_2021_impact-of-water-vapor-diffusion-and-latent-heat-on-the-effective-therm.pdf` -- Fourteau, Domine, Hagenmuller (2021). The Cryosphere 15, 2739-2755. `10.5194/tc-15-2739-2021` | **read** -- Eq. (18), the vertical effective thermal conductivity under fast kinetics at five temperatures as a quadratic in the ice volume fraction, with 917 kg/m3 as the normalising ice density. ADOPTED. Also the statement that the fast and slow kinetics limits both remain plausible, which is the bracket recorded above |
| Calonne et al. (2011). *Numerical and experimental investigations of the effective thermal conductivity of snow.* Geophys. Res. Lett. 38, L23501. `10.1029/2011GL049234` | **NOT OBTAINED.** Every open-access route and Sci-Hub returned 403. It is the slow-kinetics arm of the bracket named above, and having it would let that arm be stated as a number rather than as a direction. Nothing adopted here depends on it: its conclusion about the needle probe is carried by Riche and Schneebeli, which was obtained, and the adopted relation is Fourteau's |
| Yen (1981). *Review of thermal properties of snow, ice and sea ice.* CRREL Report 81-10 | **NOT OBTAINED**, and named by WORLD-A9S5 as the second candidate. It would have supplied a temperature-dependent conductivity for pure ice, which is the bound `CKAPI`'s declaration sits against; that bound is currently stated only as "pure ice's, plus a brine term this model cannot compute" |
