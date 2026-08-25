# Implicit Earth in the candidate ocean tiers

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from several real disciplines: **Vesper is a fictional planet and this is
engineering work on the simulation of it.** Ocean, salinity, plankton, sediment
and weathering name modelled quantities of an invented world throughout. Nothing
here is Earth science.

Audited 2026-08-24. This is OCN-12: the implicit-Earth audit applied to the
CANDIDATE ocean tiers before either can be chosen, and it is scoping rather than
adoption. `vendor/cgenie` is a candidate under OCN-3 and nothing reads it;
`references/marbl` is reference source. Nothing
below proposes to adopt, wire or fix anything.

**Method and vocabulary are `notes/audits/model-earth-centrism.md`'s and
`notes/audits/inherited-earth-constants.md`'s, deliberately unchanged.** The
question is the same one: where does an Earth number reach this world with
nothing saying it was chosen? The grouping is the same: live and load-bearing,
live but bounded or dormant, checked and clean, and what the audit did not
cover. The one thing this document adds is a column those two do not need,
because the candidates are external: **namelist-reachable or compile-time**, and
the finding of the audit as a whole is that the split falls in a place an evenly
sampled port would miss.

**The headline.** GOLDSTEIN and the geochemistry modules each compute the
physical size of every grid cell from a compile-time Earth radius, and they use
two different ones. Everything else on the circulation side is downstream of
that. On the ecosystem side the pattern inverts: ECOGEM's traits are entirely
namelist-reachable and its stoichiometry and photon conversion are buried in
source, while MARBL's per-population parameters are settings-reachable and its
light, burial and ballast constants are not reachable at all.

---

# Live and load-bearing

## 1. The grid's physical size comes from a compile-time Earth radius, stated twice

`vendor/cgenie/genie-goldstein/src/fortran/initialise_goldstein.F:377` sets
`rsc = 6.37e6` and `:496` builds every cell area from it as
`asurf(j) = rsc*rsc*ds(j)*dphi`. Independently,
`vendor/cgenie/genie-main/src/fortran/cmngem/gem_cmn.f90:802` declares
`const_rEarth = 6.37E+06` as a `PARAMETER`, and BIOGEM, SEDGEM and ATCHEM build
their own areas from that one:
`vendor/cgenie/genie-biogem/src/fortran/biogem_data.f90:1994`,
`vendor/cgenie/genie-sedgem/src/fortran/sedgem_data.f90:507`,
`vendor/cgenie/genie-atchem/src/fortran/atchem_data.f90:231`.

Neither is namelist-reachable. Cell area scales as radius squared, so at 1.20
radii the modelled surface is 0.694 of Vesper's: **ocean area, ocean volume and
every total tracer inventory are 30.6 per cent too small, and every per-area flux
is mis-scaled to match.** Silent, and the area is printed at
`initialise_goldstein.F:498` as an ordinary number.

**The porting trap is that fixing one does not fix the other.** `go_rsc` is
exported at `initialise_goldstein.F:2161` and consumed only by ENTS and the loop
wrapper; BIOGEM re-derives GOLDSTEIN's own scalings from `const_rEarth` instead,
at `biogem.f90:2620`, `:3140`, `:4277`,
`biogem_data_netCDF.f90:3837, 3916, 4036` and `biogem_data_ascii.f90:1432`. A
source edit that corrects `rsc` alone leaves the geochemistry on the other
radius, and nothing compares them.

## 2. The non-dimensionalisation carries Earth's radius and gravity into everything

`initialise_goldstein.F:370-404`. The model is entirely non-dimensional, so a
wrong scale does not produce a wrong-looking number; it produces a plausible one.

| symbol | line | literal or form | reach | wrong by |
| --- | --- | --- | --- | --- |
| `rsc` | `:377` | `6.37e6` | compile-time | 1.2002 |
| `gsc` | `:387` | `9.81` | compile-time | 1.306 |
| `dsc` | `:379` | `par_dsc`, 5000.0 | namelist | correct |
| `rhosc` | `:389` | `rh0sc*fsc*usc*rsc/gsc/dsc` | derived | 0.919, about 8 per cent |
| `tsc` | `:391` | `rsc/usc` | derived | 1.2002 |
| `opsisc` | `:394` | `dsc*usc*rsc*1e-6` | derived | 1.2002 |
| `rfluxsc` | `:398` | `rsc/(dsc*usc*rh0sc*cpsc)` | derived | 1.2002 |
| `rpmesco` | `:404` | `rsc*saln0/(dsc*usc)` | derived | 1.2002 |

`opsisc` scales the overturning streamfunction, so a reported overturning is
about twenty per cent low; `rfluxsc` scales the heat flux; `rpmesco` is the
freshwater-to-salinity conversion and was not previously in the inventory.
`initialise_goldstein.F:642` carries the source's own record of a scale having
been wrong before and corrected offline, which is the same failure arriving from
the other direction.

## 3. Rotation FAILS OPEN, in the one planetary constant the source parameterises

`initialise_goldstein.F:380-386`:

    c     CL (01/15/24) : used sidereal day length for Coriolis effect scaling factor
    c     AR (24/03/03) : check for sidereal and solar days being different, otherwise use original code for backwards compatabilty
          if (abs(sodaylen-sidaylen).gt.0.001)then
             fsc = 4*pi/sidaylen
          else
             fsc = 2*7.2921e-5
          endif

Both day lengths are in `ini_gold_nml` at `:207` and default to 86400.0 and
86164.0, so the default configuration takes the correct branch. **The failure
mode is an author who does not know the two are different and sets them equal**,
which reverts the Coriolis scaling to Earth's `2*7.2921e-5` silently. Vesper's
solar and sidereal days differ by 733.5 s, so a correctly specified
configuration takes the right branch and gets `fsc = 1.17151e-4`; the mistaken
one gets 1.45842e-4, **1.245 times too strong**, and the error does not scale
with how wrong the day length is. `fsc` is never printed except under
`debug_init`.

This is the pattern this project's conventions forbid, in the one planetary
constant the source does expose. `adrag` at `:1119` uses the sidereal day
correctly, so the drag timescale is rotation-aware and the Coriolis scaling is
the only fail-open.

## 4. Depth becomes pressure through Earth gravity, with no gravity on the line

NEW, and the largest item this audit adds.
`vendor/cgenie/genie-main/src/fortran/cmngem/gem_carbchem.f90:99`:

    loc_P        = dum_D/10.0

with the comment at `:56` reading "1 m depth approx = 1 dbar pressure". That is
seawater density times Earth gravity, compile-time, and it is the same shape as
the CLIMBER-X `p = 0.1*(-z)` the survey already records. It was not known to be
present in cGENIE as well.

`dum_D` is a cell or seafloor depth in metres, supplied by SEDGEM at
`sedgem.f90:154-155, 368, 377, 386, 722` and by BIOGEM at
`biogem_box.f90:2183`, `biogem.f90:829, 3310`, `biogem_data.f90:3493`. Every
pressure-corrected equilibrium constant runs through it: K1, K2, KB, KW, KSi,
KHF, KHSO4, KP1-3, calcite and aragonite solubility, KHS, KNH4.

Vesper's conversion is `1027.649 * 12.81 = 0.13156` bar/m against the model's
0.1. Computed through the model's own `fun_corr_p`
(`gem_util.f90:1238-1242`) at 2 C, calcite solubility is 1.066 times the
model's at 1000 m, 1.200 at 3000 m, **1.266 at 4000 m** and 1.331 at 5000 m.
**Calcite is about 27 per cent more soluble at 4 km than the model computes, so
the lysocline and the carbonate compensation depth sit substantially shallower
than SEDGEM would place them.** Silent: the saturation states look ordinary.
Bracketed: this uses the model's own reference density and ignores seawater
compression with depth, which adds a further couple of per cent in the same
direction by 5 km, so these are lower bounds.

**MARBL inverts on this axis and is clean.** Pressure is a driver-supplied
interior forcing, `references/marbl/src/marbl_interface_private_types.F90:1896-1898`,
consumed at `marbl_interior_tendency_mod.F90:306, 335, 372`, and MARBL hardcodes
no depth-to-pressure conversion anywhere. Gravity enters through the driver.

## 5. The equation of state is fitted to Earth seawater and divided by a wrong scale

`vendor/cgenie/genie-goldstein/src/fortran/eos.f:16, 19-20` is a linearised
Winton and Sarachik (1993) fit in temperature and salinity anomaly. Its
coefficients are at `initialise_goldstein.F:650-660` and **every one is divided
by `rhosc`**, so the equation of state carries finding 2's eight per cent error
into the buoyancy-to-Coriolis ratio. Compile-time. The thermobaricity term at
`:658` is gated by the namelist `ieos`, default 0, and its stated domain is
Earth's: the comment at `:654-656` reads "Optimised for `-1<deep T<6`,
`S=34.9`".

**The reference salinity is stated twice with different values, and one of them
shadows the other.** `saln0` is namelist-reachable (`ini_gold_nml` at `:228`,
default 34.9) and a shipped configuration sets 33.9
(`vendor/cgenie/genie-main/configs/cgenie.eb_go_gs_ac_bg.p0055c.BASES.config:141`),
while `vendor/cgenie/genie-goldstein/src/fortran/err_gold.F:234` hardcodes
`saln0 = 34.9` in a local of the same name. **That is a live one-PSU
inconsistency between the ocean and the error function in a standard
configuration, not only an Earth-centrism.** Silent.
`gem_cmn.f90:584` adds a third statement of the same substance,
`conv_m3_kg = 1027.649`, commented "from Winton and Sarachik [1993] @
34.7o/oo,0'C".

**The salinity range is a clamp, and it is in the chemistry rather than the
ocean.** GOLDSTEIN does not bound its salinity tracer at all. `gem_carbchem.f90:74-95`
substitutes at `par_carbchem_Smin`/`Smax`, defaults 26.0 and 43.0, and
`par_carbchem_Tmin`/`Tmax`, defaults 2.0 and 35.0 C, all namelist-reachable, with
no warning when the substitution happens. Vesper's declared salinity sits inside
the range, so the clamp is not binding at the mean, but the Mehrbach and Millero
fits underneath it are Earth seawater ionic composition and no clamp reaches
that.

## 6. Six statements of the year, none of them reachable together

Verified against the tree as it stands; the survey's inventory is correct and the
paths for two of them are not where a reader would look.

| statement | file:line | literal | reach |
| --- | --- | --- | --- |
| `conv_yr_d` | `gem_cmn.f90:585` | `365.25`, with `!360.00 !365.0` on the same line | compile-time `PARAMETER` |
| `conv_yr_hr` | `:586` | `24.0 * conv_yr_d`, the 24-hour day as a literal | compile-time |
| `conv_yr_s` | `:587` | 31,557,600 s | compile-time |
| sediment diagenesis | `sedgem_box_archer1991_sedflx.f90`, 12 occurrences from `:211` to `:838` | `3.15e7` | compile-time |
| `global_daysperyear` | `vendor/cgenie/genie-main/genie_control.f90:190` | `365.25` | compile-time `parameter` |
| the generated config | `vendor/cgenie/genie-main/runmuffin.sh:229` | `3600.0*24.0*365.25/...` as bash literals | shell |

`genie_control.f90` and `initialise_genie.F` are directly under `genie-main`, not
under `genie-main/src/fortran`. `initialise_genie.F:110` sets a one-Earth-hour
master timestep and `:115` an Earth solar constant of 1368.0 W/m2 five lines
below it, both namelist-overridable.

`conv_yr_s` is 31,557,600 s against this world's 15,794,006, so **every rate
carried through it is 2.00x out**, and the sediment module's `3.15e7` is 0.183
per cent shorter again, so a port that fixes `conv_yr_d` leaves diagenesis on an
Earth year and nothing fails.

**One place the Earth calendar fails LOUDLY, and it is unreachable by default.**
`vendor/cgenie/genie-goldstein/src/fortran/outm_netcdf.F:66-71` refuses unless
`mod(yearlen, 30.0)` is zero, which the model's own default 365.25 does not
satisfy; it survives only because `netout` defaults to off. The mirror is at
`vendor/cgenie/genie-embm/src/fortran/outm_netcdf_embm.F:61`. Vesper's 146.24
fails it too. It is the only member of this whole inventory that refuses rather
than proceeds.

A Vesper configuration can set `go_yearlen`, `go_sodaylen`, `go_sidaylen` and
`go_nyear` from the namelist. It cannot reach `conv_yr_s`, `conv_yr_hr`,
`global_daysperyear`, the sediment module's `3.15e7`, or `runmuffin.sh`, and
nothing checks the two halves agree.

## 7. The atmospheric mole inventory is Earth's by construction, and it is computable here

`vendor/cgenie/genie-atchem/src/fortran/atchem_lib.f90:133` carries
`par_atm_th = 7777.0` m as a compile-time `parameter` with no namelist path, and
`gem_cmn.f90:583` carries `conv_atm_mol = 1.7692e+020`. The same literal is
repeated inline at
`vendor/cgenie/genie-atchem/src/fortran/atchem_box.f90:328` where the named
constant was already in scope.

`gem_cmn.f90:569-582` states the derivation, which is what makes it portable:
7777 m is `R*T/(M*g)` at 273.15 K on Earth, 7994 m, scaled by 2.123/2.1839 so
that 1 ppm CO2 equals OCMIP's 2.123 PgC. **The Vesper equivalent of the 8000 m
is 6122 m and of the 7777 is 5952 m; `conv_atm_mol` scales as area times
thickness, 1.44 x 0.7658 = 1.103x, giving 1.951e20.** All three plus
`const_rEarth` must move together and none is namelist-reachable. Silent.

## 8. Latitude as a physical proxy, in four places, all in GOLDSTEIN

The project's rule is that latitude is valid for geometry, Coriolis and
diagnostics and is not a water-mass, productivity or regime classifier. Searched
GOLDSTEIN, BIOGEM, SEDGEM, ECOGEM and MARBL.

- **Equatorial drag enhancement.**
  `vendor/cgenie/genie-goldstein/src/fortran/drgset.f:31, 33` multiplies drag by
  `drgf*drgf` = 9 within `jeb` rows of the equator and by 3 in the next row,
  with `jeb = 1` and `drgf = 3.0` both compile-time at
  `initialise_goldstein.F:1121-1123`. The comment at `:18` says "Increase drag
  near equator, assuming domain is symmetric".
- **The initial hemisphere split.** `initialise_goldstein.F:1181` and `:1713`
  branch on `j.le.jmax/2` to apply `temp0` or `temp1`.
- **Polar mixed-layer substitution.**
  `vendor/cgenie/genie-goldstein/src/fortran/goldstein.F:248-249` replaces the
  local wind mixed-layer energy with its zonal mean in the two rows nearest each
  pole. Compile-time.
- **A freshwater hosing band at Earth's North Atlantic.**
  `initialise_goldstein.F:950-996` runs from `sin(50N)` to `sin(70N)` and is
  restricted to the "Atlantic" column range. Namelist-gated by `hosing`, default
  0.00, so dormant.

**And the basin identification is Earth's modern coastline written as grid
fractions.** `initialise_goldstein.F:879-919` carries three hardcoded index
hacks whose own comments read "This bit to get the Southern tip of Greenland
into the Atlantic", "This bit to get the Arctic all Atlantic" and "This bit to
get the Southern tip of Greenland out of the Pacific". They set the basin index
arrays that feed the basin overturning and freshwater diagnostics and the hosing
region. Compile-time, gated on `igrid.eq.0`. On any other coastline they are
meaningless.

**BIOGEM, SEDGEM, ECOGEM and MARBL are clean.** In the first three, latitude is
set from the grid and consumed only by netCDF and ASCII output. In MARBL a
case-insensitive search for latitude, hemisphere and equator over
`references/marbl/src/` returns nothing at all: latitude appears only in the test
driver and the documentation. MARBL is driver-agnostic on this axis by
construction.

## 9. PAR and the photon conversion, and BIOGEM has nowhere to put the correction

The named cites verify. ECOGEM's `PARfrac`, `k_w` and `k_chl` are all in
`ini_ecogem_nml` (`ecogem_lib.f90:106-117`, defaults 0.43, 0.04 and 0.03).
MARBL's `f_qsw_par = 0.45_r8` is a Fortran `parameter` at
`references/marbl/src/marbl_settings_mod.F90:160` with no YAML entry, and its
six attenuation literals are inline at
`marbl_interior_tendency_mod.F90:644-649`. ECOGEM's buried solar constant is at
`vendor/cgenie/genie-ecogem/src/fortran/ecogem_box.f90:297`,
`E0 = PARlocal/0.2174`, compile-time, and is the marine twin of the fixed 550 nm
coefficient BIO-25 found in LPJ photosynthesis.

Using the survey's own band and photon-energy ratios: `PARfrac` 0.43 becomes
**0.365**, `f_qsw_par` 0.45 becomes **0.382**, and the 0.2174 W per micromole
becomes **0.2140**, the last a 1.6 per cent correction because the tilt within
400 to 700 nm barely matters.

**The structural finding this audit adds: BIOGEM has no PAR fraction at all.**
`vendor/cgenie/genie-biogem/src/fortran/biogem_box.f90:701-704` attenuates the
FULL shortwave field by a single scalar e-folding depth,
`par_bio_I_eL` at 20.0 m from OCMIP-2, and compares the result against
`par_bio_c0_I` at 20.0 W/m2 from Doney et al. (2006). Both are namelist-reachable
and both are solar calibrations, but there is no band split to correct: 0.618 of
this star's shortwave sits above 0.75 um against the Sun's 0.483, and all of it
is deposited within about half a metre. ECOGEM's structure is sounder because
`PARfrac` removes the non-photosynthetic half before `k_w` is applied, so
`k_w = 0.04` is a within-band coefficient and barely moves. **The porting
difficulty on the light axis is not that a constant is wrong; it is that
BIOGEM's structure has nowhere to put the correction.**

## 10. ECOGEM buries in source what BIOGEM exposes in a namelist

BIOGEM's four Redfield ratios are in `ini_biogem_nml`
(`biogem_lib.f90:179-184`) with `par_bio_red_PC_flex` at `:184-187` switching
C:P to flexible entirely. ECOGEM's are literals:

| literal | file:line | what it is |
| --- | --- | --- |
| `138.0 / 106.0` | `vendor/cgenie/genie-ecogem/src/fortran/ecogem.f90:778` | O2 production from DIC uptake |
| `-16.0` | `.../ecogem.f90:783` | alkalinity drawdown from PO4 uptake |
| `16.0` | `.../ecogem_box.f90:311` | N uptake inferred from P uptake under `pquota` |
| `6.625` | `.../ecogem_data.f90:658` | chlorophyll initialisation |
| `40.0` | `.../ecogem_box.f90:308` | diazotroph N:P, with the source's own comment "Still need to check Moore et al (2002) + parameterise" |

The fifth was not in the inventory. All compile-time. **So the ecosystem tier
reverses the pattern the circulation tier sets**, and an audit that assumes one
rule for the whole tree spends its effort in the wrong place twice.

---

# Live, bounded or dormant

## 11. Closures are namelist-reachable, and the scale factors under them are not

`diff(1)` isopycnal 2000 m2/s and `diff(2)` diapycnal 1e-5 m2/s are namelist
(`definition.xml:2527-2538`) but are non-dimensionalised at
`initialise_goldstein.F:1125-1126` by `rsc` and `dsc`, so both carry finding 2.
`adrag` at 2.5 days is namelist and is correctly scaled by `sidaylen*fsc`. Three
things beside it are compile-time: `drgf = 3.0`, `kmxdrg = kmax/2` and `jeb = 1`
(`:1121-1123`). The isopycnal slope limits `ssmaxsurf` and `ssmaxdeep` are
namelist and their descriptions say they are "NOT a tunable parameter", while the
depth structure they are applied over, a 200 m e-folding centred at -300 m, is
compile-time at `:2152-2153`. Three mixed-layer coefficients are namelist and two
carry an annotation saying they were tuned in isolation by a named person.

**Convective adjustment is clean.** `vendor/cgenie/genie-goldstein/src/fortran/co.F:70`
mixes on density inversion alone, with no threshold, no timescale and no Earth
constant. It inherits the equation of state's problems and adds none.

## 12. There is no tidal scheme, and what stands in for it is fitted to Earth abyssal observations

No tidal forcing, tidal mixing or tidal dissipation term exists anywhere in
`genie-goldstein`. The bottom-intensified diapycnal profile that plays that role
is `vendor/cgenie/genie-goldstein/src/fortran/ediff.F`, and its comment at `:82`
names its own provenance as consistency with Earth observations from the internal
wave and abyssal mixing literature. A Levitus global-mean density profile with a
650 m e-folding is at `:75-76`, the mixing profile with a 2500 m reference depth
and a 700 m e-folding at `:84`, the cap at `:88`, and a Bryan and Lewis (1979)
alternative at `:97`. All compile-time; `iediff` defaults to 0, so the whole
scheme is dormant unless switched on.

**Geothermal is namelist-reachable, default off, and set by every shipped
configuration.** `par_Fgeothermal` defaults to 0.0
(`biogem_lib.f90:131-132`) and every user config sets 100 mW/m2, Earth's mean
ocean-floor heat flux, applied at `biogem.f90:1098-1107` where it is also
multiplied by `conv_yr_s` and so carries finding 6's 2.00x. **Vesper's is not
computable from anything this project determines**: a mass-over-area scaling
gives 131 mW/m2, so the honest statement is a bracket of 100 to 131 mW/m2 with
the specific radiogenic content undetermined, not a value.

## 13. Initial hydrography is two hemispheres of uniform water, and the salinity is not settable

`initialise_goldstein.F:1176-1191` sets temperature to `temp0` in the southern
half and `temp1` in the northern, both namelist and both defaulting to 5.0 C, and
sets the initial salinity ANOMALY to `0.0` with no way to change it, so the whole
ocean starts at `saln0`. There is no shipped three-dimensional initial state: the
Levitus `.silo` files in the tree are observational targets for the error
function, not initial conditions. The values are Earth-free; the Earth-shaped
parts are the hemisphere branch of finding 8 and the fact that a stratified
initial state is not expressible.

## 14. Sinking speeds are the gravity question the rotation answer did not cover

The survey settled the ROTATION half of this: `par_bio_remin_sinkingrate` and its
siblings are in m/d and `par_bio_remin_opal_K` and its siblings in /d, where the
day is 86400 s and nothing else, so a 30-hour rotation must not rescale them.
**Gravity is the axis a settling speed does depend on, and it was not covered.**

`par_bio_remin_sinkingrate` is namelist at 125.0 m/d
(`biogem_lib.f90:278-280`). MARBL states the same physics as a LENGTH:
`parm_POC_diss` 100 m, `parm_SiO2_diss` 650 m, `parm_CaCO3_diss` 500 m
(`references/marbl/src/marbl_settings_mod.F90:451-455`), all settings-reachable.
A length scale is a speed divided by a rate and MARBL exposes neither factor, so
**there is no defensible way to rescale MARBL's export profile for a different
gravity from inside the model, and cGENIE's tier does have the seam.**

**Bracketed, because the settling regime is not determined here.** Under Stokes
settling the speed is linear in gravity, giving 163 m/d and a 131 m
remineralisation length; marine aggregates are largely outside that regime, where
the exponent is nearer 0.5 to 2/3, giving 143 to 151 m/d. The bracket is 143 to
163 m/d and 114 to 131 m, and narrowing it needs a settling-regime determination
this tree does not carry. Silent either way: the model produces an export profile
that looks ordinary and is too shallow.

ECOGEM's own `biosink_a` and `biosink_b` default to 0.0 and the path is dead
code, which `notes/audits/ecosystem-tier-ecogem-marbl.md` section 9 establishes;
the gravity question in that tier is BIOGEM's.

## 15. A whole tuned physical parameter set, carried to sixteen significant figures

`vendor/cgenie/genie-main/configs/cgenie.eb_go_gs_ac_bg.p0055c.BASES.config:95-134`
sets the wind stress scaling, the ocean isopycnal and diapycnal diffusivities,
the inverse minimum drag, eight EMBM diffusion and advection parameters and the
sea-ice eddy diffusivity as sixteen-figure literals. `definition.xml:2516-2524`
names the provenance for the wind stress scaling: "Tuned values: 1.6674 (EnKF),
1.1841 (ACCPM), 1.3005 (NSGA-II)", three ensemble optimisations against Earth
observations.

OCN-17's verdict already treats the wind stress scaling as a declared bracket.
**The other six ocean and sea-ice values in that block are the same class and
have not been named.** Their observational targets are first-class build
products: `errfn_ea_go_gs.f90` and its siblings compile to binaries scoring the
model against Levitus fields through `err_gold.F`, and on Vesper every one is
inert.

Beside them, the shipped Earth fields a configuration reads: NCEP reanalysis
surface winds and four wind-stress interpolation files per palaeogeography in
`vendor/cgenie/genie-paleo/`, Albani Earth dust reconstructions in
`vendor/cgenie/genie-forcings/`, and prescribed windspeed and sea-ice fields
that BIOGEM's gas exchange reads. Gas transfer itself is
`biogem_box.f90:121`, whose `par_gastransfer_a` is namelist at Wanninkhof's
0.310 but whose `1.515E-3` is 1/660, CO2's Schmidt number in Earth seawater at
20 C, hardcoded with a comment saying so, and whose `conv_yr_hr` carries the
2.00x calendar.

## 16. The Earth nutrient inventory is configuration, and that is the good news

**Nothing in cGENIE source welds a nutrient inventory.** `ocn_init` is a
107-element namelist array (`biogem_lib.f90:23-24`) whose every element defaults
to zero. Earth's inventory arrives entirely through configuration files: DIC
2.244e-3 mol/kg, PO4 2.159e-6, O2 1.696e-4, alkalinity, sulfate, calcium and
chloride, all in the shipped base configs. So this row's "Earth nutrient
inventory" item is **a Vesper decision to be made rather than an Earth number to
be extracted**, and it is the one item in OCN-12's list that costs nothing to
port.

Two compile-time exceptions, both composition rather than inventory:
`gem_cmn.f90:809-813` carries `const_conc_Mg = 0.05282` and
`const_conc_MgtoCa = 5.155`, explicitly the modern mean ocean at S = 35.

**MARBL's optional abiotic-DIC module does weld three.**
`references/marbl/src/marbl_abio_dic_surface_flux_mod.F90:141, 143` set silicate
7.5 and phosphate 0.5 nmol/g from Orr et al. (2017), and `:165` sets surface
alkalinity as `2310 * rho_sw * sss/34.7`, Earth's mean surface alkalinity over
Earth's reference salinity.

---

# Checked and clean

Recorded so they are not re-derived.

- **GOLDSTEIN's freezing point is salinity-dependent, and ExoPlaSim's is not.**
  `vendor/cgenie/genie-goldstein/src/fortran/surf_ocn_sic.F:475-477` and
  `vendor/cgenie/genie-goldsteinseaice/src/fortran/surflux_goldstein_seaice.F:159`
  compute it from the local salinity. This is the class
  `model-earth-centrism.md` finding 10 opened on `TFREEZE`, and the candidate
  does it correctly. The pressure term is omitted, which is right at the
  surface.
- **The convective adjustment carries no Earth constant.**
- **The `.eco` plankton files are not Earth data.** They carry three columns,
  functional type, diameter and replicate count; every trait is generated from
  namelist allometric coefficient pairs, so the trait space is fully
  namelist-reachable. The Earth fitting is in those coefficients, which are the
  Ward et al. (2012) laboratory-culture fits, and being reachable they are
  declarable.
- **GOLDSTEIN's idealised wind profile is dead code**, commented out at
  `initialise_goldstein.F:1216-1229`; wind stress comes from the coupler.
- **MARBL ships no Earth field of its own** outside the abiotic-DIC module, and
  its restoring fields and timescales are driver-supplied.
- **Two constants that happen to be nearly right here.**
  `initialise_goldstein.F:1440` `rhoair = 1.25` and `:1443-1447` the Bolton
  (1980) saturation constants embed a 1000 hPa reference. Vesper's declared
  surface pressure is 1.0000 bar with an Earth-like composition, so both are
  within about a per cent. They would bite on any world with a different surface
  pressure.

---

# What would have to be true before a tier could be chosen

This is the scoping answer, and it is a list of preconditions rather than a
recommendation. Nothing here selects anything.

1. **The two radii would have to become one, and it would have to be
   Vesper's.** `rsc` and `const_rEarth` are independent statements of the same
   quantity in a model where cell area, volume and every inventory derive from
   them, and BIOGEM re-derives GOLDSTEIN's own scalings from the second. A port
   that changes one is worse than a port that changes neither, because the two
   halves then disagree silently. **The precondition is a check that they
   agree**, in the same shape as this project's own consistency gates, before
   any number is changed.
2. **The depth-to-pressure conversion would have to become a gravity.** Finding
   4 is a single division by ten standing for `rho*g`, and it reaches every
   carbonate equilibrium constant and therefore the carbonate compensation
   depth. This is the item that most changes what the tier would SAY about
   Vesper, as distinct from what it would compute.
3. **The frictional-geostrophic closure's Earth-fitted transport parameters
   would have to be declared as a bracket, all seven of them.** OCN-17 already
   does this for the wind stress scaling. Finding 15 says the other six in the
   same config block are the same class, and the observational targets that set
   them are inert here. The audit cannot say whether the closure survives 1.20
   radii; it can say that the answer is not readable from the source and that
   the model's own error function cannot be used to find it.
4. **The calendar would have to be settled in six places at once**, one of which
   is a shell script that generates the config every standard run uses. Four of
   the six are compile-time. Until then, any rate reported per year is 2.00x
   out and sediment diagenesis is out by a further 0.183 per cent in the other
   direction.
5. **The rotation fail-open would have to be closed, or the configuration
   audited for it.** Setting the solar and sidereal day lengths equal is a
   plausible mistake for an author who has not met the distinction, and it
   silently substitutes Earth's Coriolis scaling. This is cheap to guard and
   expensive to discover.
6. **A settling regime would have to be chosen** before either tier's export
   profile means anything, and only cGENIE's tier has the seam to apply it at.
   The bracket is finding 14's.
7. **A geothermal heat flux and a nutrient inventory would have to be
   DECLARED.** Neither is an Earth number to be corrected: the first is
   bracketed at 100 to 131 mW/m2 with its specific radiogenic content
   undetermined, and the second is configuration in cGENIE and therefore free.
   Both are Vesper design decisions that the models will not raise.

**And one instruction in the row is not executable from this tree.**
`ClimaParams.jl` is not present. It appears only as a dependency declaration in
`references/cloudmicrophysics/Project.toml:7, 28` and in this project's own
notes; there is no Julia depot, no vendored copy, and no "Planetary and Orbital
Parameters (Earth Defaults)" block anywhere. Acquiring it is a fetch rather than
a blocker, and this inventory was taken without it.

---

## What this audit did NOT cover

- **Nothing was compiled and nothing was run FOR THIS AUDIT.** Every claim here
  is from source and from shipped configuration.
  `notes/audits/cgenie-build-cost-and-grid-ceiling.md` is the separate document
  that builds and runs it.
- **Vesper's geothermal heat flux is bracketed, not derived.** The upper end is
  a mass-over-area scaling and the specific radiogenic content and thermal age
  are undetermined in this project.
- **The settling regime for marine aggregates is not determined**, which is what
  decides whether the sinking correction goes as gravity to the first power, the
  two-thirds power or the half. The bracket is reported rather than narrowed.
- **Whether the carbonate chemistry's salinity clamp binds in practice on Vesper
  is unknown.** The mean sits inside the range; the tails of a simulated
  distribution cannot be known without a run.
- **The island count Vesper's coastline generates is unknown**, so whether
  `GOLDSTEINMAXISLES` at 5 needs raising is unanswerable before a candidate
  configuration on an accepted terrain.
- **`genie-embm`, `genie-goldsteinseaice`, `genie-atchem` and `genie-ents` were
  not swept systematically**, only where they carry a constant the ocean or
  sediment tiers consume. Two known items are recorded and not pursued: the EMBM
  planetary albedo at
  `vendor/cgenie/genie-embm/src/fortran/initialise_embm.F:903-911`, which is a
  pure cosine series in latitude with the original Earth literals commented on
  `:906` and an alternative branch reading a one-dimensional per-latitude table,
  and the sea-ice constants at `initialise_goldstein.F:1462-1480`. The EMBM is
  the archetypal latitude-as-classifier, but the survey's architecture decision
  eliminated the path that would reach it, so it is out of scope and recorded
  here so it is not re-found.
- **MARBL's `defaults/` YAML settings were not read**, only the Fortran defaults
  they override. For at least one constant this matters: `parm_SiO2_gamma` is
  0.00 in Fortran and 0.1 in every YAML, and which is authoritative for a given
  run is decided outside `references/marbl`.
- **The silicon axis is not here.** It is
  `notes/audits/ecosystem-tier-ecogem-marbl.md` section 13, which closes it for
  OCN-4.

---

# The ADOPTED tier, on the same question

Audited 2026-08-25, and this section inverts the document's scope. Everything
above is about candidates. The ocean that RUNS is ExoPlaSim's slab --
`oceanmod.f90` integrates one mixed layer, `icemod.f90` runs the sea-ice
thermodynamics and owns the sea water, `seamod.f90` is the surface the
atmosphere sees -- and the same question asked of it returns the same headline
in a different place.

## A1. Two quantities were stated twice, and one of the splits was reachable from configuration

**The melting point.** `icemod.f90` declared `parameter(TMELT=273.16)` while
`tmelt` is a `planet_nl` key that `p_earth.f90` itself describes as the freezing
point every soil, snow, sea and ice routine tests against. Both statements were
live on the same modelled ice: `seamod`'s sea-ice albedo ramp anchors on
pumamod's and `icestep`'s skin-temperature melt anchored on icemod's. This is the
audit's finding-1 class with a worse property than the two radii have --
**it splits from configuration alone**, because one of the two is namelist-
reachable and the other could only be reached by a source edit, so a run that
moved the key would have left no line carrying the melting point at all.

**The sea-ice density.** `oceanmod.f90` and `icemod.f90` each declared
`parameter(CRHOI = 920.)`, neither reachable. `addfc`'s branch (c) builds its
melt flux `yiced*CRHOI*CLFI/dtmix` from oceanmod's while icemod weighs the same
ice, and its snow-ice flooding threshold `CRHOS-CRHOI`, from its own. This one
splits only under a source edit, which is what the audit says a port does first.

Both are now stated once and handed on, the melting point through `iceini`
beside the snow density GRAV-8 already routed that way, and the sea-ice density
through `oceanini` beside the sea-water quartet that is there for exactly this
reason. `exoplasim/scripts/ocean_tier_gate.py` is the check the audit's
precondition 1 asks for, in the shape it asks for it: the owner's declaration,
every receiver's assignment, the handoff, and the forms that must not come back.

## A2. A threshold was reachable in one of its two uses

`icemod.f90` tests the modelled ice compactness against two 0.5 thresholds.
`thicec` decides whether a cell is MASKED as iced, at `icestep` and `mkicec`, and
is an `icemod_nl` key. `cicemin` decides whether falling snow lands on the ice or
melts into the water, at `subsnow`'s two sites, and was reachable from nothing.
So a bracket that moved the mask left the snow partition on the compiled value.

Exposed rather than merged, and that is the substantive judgement here: the model
uses them for different decisions and nothing in the source says they are one
number, so making them one would be a physics change made for tidiness. This is
the remedy `hlead` got, which was a local variable credited to Hippler (1979) and
unreachable from any namelist.

## A3. The reachability column, for the tier that runs

Of the constants of the modelled sea and its ice that the declaration covers,
the split falls in a different place than it does on the candidate side. The
sea water quartet, the two sea-ice lengths, the lead scale, the albedo ramp and
the mixed-layer depth are all namelist keys and this project sets them from
`config/planet.yaml`. What is not reachable is the **material properties of the
modelled ice and snow**: density, the two specific heats, the two conductivities
and the heat of fusion of snow. Those are compile-time and each is an Earth
measurement standing where nothing says it was chosen. They are what
`ocean_tier_gate.py --strict` refuses on, together with the namelist keys that
are reachable but still carry an Earth number with no source.

The conductivity pair is the sharpest of them. `CKAPI` and `CKAPSN` set the
conductive flux through the modelled ice and the insulation the snow on it
provides, which is what fixes the equilibrium thickness at a given surface energy
balance -- and the peer that runs this configuration reports excessive ice.

## A4. What this section did NOT establish

- **Nothing was run.** Every statement is against the source. What the slab
  DELIVERS is verified by running the model, in section 11 and in
  `verify_ocean_flux_channel.py`; nothing here repeats that.
- **No number was changed.** Both duplicated quantities held the same value on
  both sides and `planet_nl`'s `tmelt` is what the compiled `TMELT` was, so
  every current configuration integrates identically. What changed is that the
  values can now move, and move in one place when they do.
- **The unreachable material properties are inventoried, not bracketed.** What
  ice and snow on this world are made of is a Vesper decision the models will not
  raise, and the gate names the six rather than guessing at them.
- **`icemod`'s freezing point is still salinity-INDEPENDENT**, which the
  "checked and clean" section above records from the other direction: the
  candidate computes it from the local salinity and the adopted tier takes it as
  a single declared number. That remains true and is not a defect of the
  handoffs above.

---

## Tasks

OCN-12 is the row this document answers, and it answers it as scoping: the
inventory above, and the preconditions. It does not adopt, wire or change
anything, and `vendor/cgenie` remains a candidate that nothing reads. Two
findings belong to other rows and are recorded there rather than acted on here:
the seven-parameter tuned block extends OCN-17's bracket from one parameter to
seven, and the depth-to-pressure conversion is the sharpest single item OCN-3's
viability question has to weigh.
